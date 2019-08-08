extern crate fnv;
#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate rand;
extern crate bit_vec;

use cpuprofiler::PROFILER;
use rand::Rng;
use std::f32;

use clap::{App};

use std::io::BufReader;
use std::io::BufRead;
use std::fs::File;

use hashbrown::{HashMap,HashSet};
use fnv::FnvHasher;
use std::hash::BuildHasherDefault;
type FnvHashMap<T,V> = HashMap<T,V, BuildHasherDefault<FnvHasher>>;
type FnvHashSet<T> = HashSet<T, BuildHasherDefault<FnvHasher>>;

use std::collections::BTreeSet;
use bit_vec::BitVec;

fn main() {
    let params = load_params();
    let whitelist = load_whitelist(&params); // map from barcode string to index in whitelist
    println!("whitelist loaded");
    let (variants, barcodes, info, phasing) = load_variants(&params, &whitelist); 
    //gather_info(variants, barcodes, info);
        // list from variant to list of molecule ids (negative means alt support, positive means ref support)
                                // map from molecules to list of variants (negative means alt support, positive means ref support)
    //println!("variant support loaded");
    //println!("testing, I have this many molecules {}", barcodes.len());
    let (variants, molecules) = extract_molecules(variants, barcodes);
    //experiment(variants, molecules);
    //println!("Hello, world!");
    phuzzy_phaser_master(molecules, variants, params.ploidy, phasing, info);
}

fn extract_molecules(variants: Vec<Vec<i32>>, barcodes: FnvHashMap<i32, Vec<i32>>) -> 
        (Vec<Vec<i32>>, FnvHashMap<i32, Vec<i32>>) { //may also want to make a Vec<Vec<i32>> mapping variants to other variants that have a > threshold link
    PROFILER.lock().unwrap().start("./my-prof.profile");
    let mut new_variants: Vec<Vec<i32>> = Vec::new();
    for _i in 0..(variants.len() + 1) { new_variants.push(Vec::new()); }
    let mut molecules: FnvHashMap<i32, Vec<i32>> = FnvHashMap::default();
    let mut molecule_number: i32 = 0;
    let mut barcode_molecules: FnvHashMap<i32, Vec<i32>> = FnvHashMap::default(); // mapping from barcode to list of molecules
    for (bcindex, (bc1, vars)) in barcodes.iter().enumerate() { //ok now for each barcode bc1
        if bcindex % 20000 == 0 { println!("{}", bcindex); }
        let mut disjoint_set: Vec<usize> = (0..vars.len()).collect(); // disjoint_set data structure starting with all separate sets 0..len 
        'outer: for (var1dex, var1) in vars.iter().enumerate() { // for each variant
            let var1 = var1.abs();
            let mut counts: Vec<usize> = vec![0; vars.len()]; // counts of barcodes connecting var1 to vardex 
            for bc2dex in 0..(variants[(var1 - 1) as usize].len()) { // for each barcode bc2 containing var1. So they are guaranteed to share var1
                let bc2 = variants[(var1 - 1) as usize][bc2dex].abs();
                let mut index1 = 0; let mut index2 = 0;
                let vars1 = vars; let vars2 = barcodes.get(&bc2).unwrap();
                if vars1.len() == 0 || vars2.len() == 0 { continue; }
                loop { // pairwise step through these sorted variant lists
                    let var1 = vars1[index1].abs(); let var2 = vars2[index2].abs();
                    if index1 >= var1dex { break; }
                    if var1 == var2 { counts[index1] += 1; 
                        if counts[index1] > 3 { disjoint_set[var1dex] = index1; continue 'outer; }
                        break; 
                    } // if they contain the same variant
                    else if var1 < var2 { index1 += 1; if index1 >= vars1.len() { break; } }
                    else { index2 += 1; if index2 >= vars2.len() { break; } }
                }
            }
        }
        for vardex in 0..disjoint_set.len() {
            let mut curdex = vardex;
            while disjoint_set[curdex] != curdex {
                curdex = disjoint_set[curdex];
            }
            disjoint_set[vardex] = curdex; // now everything points to their root node 
        }
        let mut touched: Vec<usize> = Vec::new();
        for _ in 0..disjoint_set.len() { touched.push(0); } //(0..disjoint_set.len()).collect(); //[i32; disjoint_set.len()] = [0; disjoint_set.len()];
        let mut mol_list: Vec<i32> = Vec::new(); // is the vector of molecules in this barcode
        for vardex in 0..disjoint_set.len() {
            if touched[disjoint_set[vardex]] == 0 {
                molecule_number += 1; // new molecule
                touched[disjoint_set[vardex]] = molecule_number as usize;
                mol_list.push(molecule_number);
                let mut newvars: Vec<i32> = Vec::new();
                newvars.push(vars[vardex]); // vars is the list of variants in this barcode
                molecules.insert(molecule_number, newvars);
                if vars[vardex] < 0 {
                    new_variants[vars[vardex].abs() as usize].push(-molecule_number);
                } else { new_variants[vars[vardex].abs() as usize].push(molecule_number); }
            } else {
                let molecule_name: i32 = touched[disjoint_set[vardex]] as i32;
                molecules.get_mut(&molecule_name).unwrap().push(vars[vardex]);
                if vars[vardex] < 0 {
                    new_variants[vars[vardex].abs() as usize].push(-molecule_name);
                } else { new_variants[vars[vardex].abs() as usize].push(molecule_name); }
            }
        }
        barcode_molecules.insert(*bc1, mol_list);
    }
    PROFILER.lock().unwrap().stop();
    (new_variants, molecules)
}

fn phuzzy_phaser_master(molecules: FnvHashMap<i32, Vec<i32>>, variants: Vec<Vec<i32>>, 
        num_clusters: usize, phasing: FnvHashMap<i32, String>, info: FnvHashMap<i32, (String, i32, String, String, String)>) {
    let mut cluster_support: Vec<f32> = Vec::new(); // to be reused
    for _cluster in 0..num_clusters { cluster_support.push(0.0); } //initialize
    let mut sums: Vec<Vec<f32>> = Vec::new(); let mut denoms: Vec<Vec<f32>> = Vec::new();
    let vars_per_step = 20;
    for cluster in 0..num_clusters {
        denoms.push(Vec::new()); sums.push(Vec::new());
        for _ in 0..vars_per_step { denoms[cluster].push(0.0); sums[cluster].push(0.0); }
    }
    let mut cluster_centers: Vec<Vec<f32>> = Vec::new();
    let mut rng = rand::thread_rng();
    for cluster in 0..num_clusters {
        cluster_centers.push(Vec::new());
        for _ in 0..variants.len() {
            let x: f32 = rng.gen();
            cluster_centers[cluster].push(x/2.0 + 0.25); // initialize from 0.25..0.75 uniformly 
        }
    }
    let mut touched_vars: Vec<bool> = Vec::new();
    for _ in 0..variants.len() {
        touched_vars.push(false);
    }


    let mut variant_window = VariantRecruiter::new(variants, molecules);
    //pick 10 starting variant
    let mut total_variants = 0;
    println!("chrom\t\t\tpos\tref\talt\tLR_ps\tLR_gt\tconnections\thap1\thap2");
    'phaseblockloop: loop {
        let mut new_variants: Vec<usize> = Vec::new();
        let mut current_molecules: Vec<Vec<i32>> = Vec::new();
        let mut used_mols: FnvHashSet<i32> = FnvHashSet::default();
        let mut mol_connections: Vec<u32> = Vec::new();
        for _i in 0..vars_per_step {
            total_variants += 1;
            let (variant, connections) = match variant_window.get_next_variant(50) {
                Some((x,y)) => (x,y),
                None => break 'phaseblockloop,
            };
            mol_connections.push(connections);
            let vardex = variant.abs() as usize;
            new_variants.push(vardex);
            for molecule in &variant_window.variants[vardex] {
                let moldex = molecule.abs();
                if !used_mols.contains(&moldex) {
                    let mut mol: Vec<i32> = Vec::new();
                    for var in variant_window.molecules.get(&moldex).unwrap() {
                        mol.push(*var);
                    }
                    current_molecules.push(mol);
                } 
            }     
            // add our new variant
            variant_window.add_variant(variant, connections);
        }
        new_variants.sort(); 
        for variant in &new_variants { touched_vars[*variant] = true; }
        next_cluster_step(current_molecules, &mut cluster_centers, &new_variants, &mut cluster_support, &touched_vars, &mut denoms, &mut sums);
        //println!("total variants so far {}",total_variants);
        for (index, vardex) in new_variants.iter().enumerate() {
            let (chrom, pos, reference, alt, ps) = info.get(&(*vardex as i32)).unwrap();
            println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}", chrom, pos, reference, alt, ps, 
                phasing.get(&(*vardex as i32)).unwrap(), mol_connections[index],cluster_centers[0][*vardex],cluster_centers[1][*vardex]);
        }
    }
    
    
    if false {
        // lets build a total test case 
        let mut mols: Vec<Vec<i32>> = Vec::new();
        for i in 1..8 {
            mols.push((i..(i+4)).collect());
            mols.push(((-i-3)..(-i+1)).rev().collect());
        }
        let new_variants: Vec<usize> = (1..11).collect();
        println!(" and off we go");
        next_cluster_step(mols, &mut cluster_centers, &new_variants, &mut cluster_support, &touched_vars, &mut denoms, &mut sums);
        for (i, cluster) in cluster_centers.iter().enumerate() {
            println!("cluster {} = {:?}",i,cluster);
        }
    }
    
    
}

fn log_sum_exp(p: &Vec<f32>) -> f32{
    let max_p: f32 = p.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let sum_rst: f32 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

fn next_cluster_step(molecules: Vec<Vec<i32>>, cluster_centers: &mut Vec<Vec<f32>>, 
    new_variants: &Vec<usize>, cluster_support: &mut Vec<f32>, touched: &Vec<bool>,
    denoms: &mut Vec<Vec<f32>>, sums: &mut Vec<Vec<f32>>) {
    let mut vardex = 0; let mut moldex = 0; let mut varval = 0; let mut molval = 0;
    let mut iteration = 0;
    let mut change = 1.0;
    while change > 0.0001 {
        iteration += 1;
        // reset reused memory to prior expectation, think of this like a pseudocount and a way to avoid div by 0
        for cluster in 0..denoms.len() {
            for vardex in 0..denoms[cluster].len() {
                denoms[cluster][vardex] = 0.001;
                sums[cluster][vardex] = 0.001*cluster_centers[cluster][new_variants[vardex]];
            }
        } 
        //println!("iteration {}",iteration);
        //for cluster in 0..2 {
            //println!("cluster center {} = {:?}",cluster,cluster_centers[cluster]);
        //}
        // Expectation
        for (mol_num, molecule_variants) in molecules.iter().enumerate() {
            let mut sum = 0.0;
            let mut count = 0;
            
            for variant in molecule_variants {
                let var = variant.abs() as usize;
                if touched[var] { count += 1; }
            }
            if count < 2 { continue; }
            for cluster in 0..cluster_centers.len() {
                cluster_support[cluster] = 0.0;
                for variant in molecule_variants {
                    let var = variant.abs() as usize;
                    if !touched[var] { continue; } // reach into old stuff and current stuff, not new new stuff.
                    if variant > &0 {
                        cluster_support[cluster] += cluster_centers[cluster][var].ln();
                    } else {
                        cluster_support[cluster] += (1.0 - cluster_centers[cluster][var]).ln();
                    }
                    assert!(!cluster_support[cluster].is_infinite(), "infinite value, cluster center {}", cluster_centers[cluster][var]);
                    if cluster_support[cluster].is_infinite() {
                        println!("infinite value, cluster_center {}",cluster_centers[cluster][var]);
                        assert!(false);
                    }
                }
                //sum += cluster_support[cluster].exp();
            }
            sum = log_sum_exp(cluster_support);
            for cluster in 0..cluster_centers.len() { 
                let tmp = sum;
                let tmp2 = cluster_support[cluster];
                cluster_support[cluster] -= sum; 
                cluster_support[cluster] = cluster_support[cluster].exp(); 
                if cluster_support[cluster].is_nan() || cluster_support[cluster].is_infinite() {
                    println!("PUKE0 {}\t{}\t{}",tmp,tmp2,tmp);
                }
            }
            
            //for cluster in 0..cluster_centers.len() {
            //    if iteration > 3 && cluster_support[cluster] < 0.95 && cluster_support[cluster] > 0.05 {
            //        println!("mol {}, cluster {}, support {}",mol_num, cluster, cluster_support[cluster]);
            //        println!("\t{:?}",molecule_variants);
            //    }
            //}
            vardex = 0; moldex = 0; varval = new_variants[vardex]; molval = molecule_variants[moldex].abs() as usize;
            //println!("mol {} count {}\tc1:\t{}\tc2:\t{}", mol_num, count, cluster_support[0], cluster_support[1]);
            //if cluster_support[0].is_nan() || cluster_support[1].is_nan() || cluster_support[0].is_infinite() || cluster_support[1].is_infinite() {
            //    println!("PUKE1 {}\t{}",cluster_support[0], cluster_support[1]);
            //}
            //print!("\t");
            loop {
                //println!("moldex {}, molval {}, vardex {}, varval {}", moldex, molval, vardex, varval);
                if varval > molval { 
                    moldex += 1; if moldex >= molecule_variants.len() { break; } molval = molecule_variants[moldex].abs() as usize; 
                } else if molval > varval {
                    vardex += 1; if vardex >= new_variants.len() { break; } varval = new_variants[vardex];
                } else {
                    // actually update things
                    for cluster in 0..cluster_centers.len() {
                        //print!("\t{}",molecule_variants[moldex]);
                        //println!("_{} denom = {} + {} = {}",cluster,denoms[cluster][vardex], cluster_support[cluster], denoms[cluster][vardex] + cluster_support[cluster]);
                        denoms[cluster][vardex] += cluster_support[cluster];
                        //if denoms[cluster][vardex].is_nan() || denoms[cluster][vardex].is_infinite() { println!("PUKE2 {}\t{}",denoms[cluster][vardex],cluster_support[cluster]);}            
                        //print!("\t{}",molecule_variants[moldex]);
                        if molecule_variants[moldex] > 0 {
                            //println!("_{} sum = {} + {} = {}",cluster,sums[cluster][vardex], cluster_support[cluster], sums[cluster][vardex] + cluster_support[cluster]);
                            sums[cluster][vardex] += cluster_support[cluster]; // 1*cluster_support[0]
                            //if sums[cluster][vardex].is_nan() || sums[cluster][vardex].is_infinite() { println!("PUKE3 {}\t{}",sums[cluster][vardex],cluster_support[cluster]);}            

                        } //else { 
                        //    println!("\tsum no change = {}",sums[cluster][vardex]);
                        //}// else add 0*cluster_support[cluster]
                    }//println!();
                    moldex += 1;  if moldex >= molecule_variants.len() { break; } molval = molecule_variants[moldex].abs() as usize;
                    vardex += 1;  if vardex >= new_variants.len() { break; }  varval = new_variants[vardex];
                }
            }
        }
        // Maximization
        change = 0.0;
        for (vardex, variant) in new_variants.iter().enumerate() {
            for cluster in 0..denoms.len() {
                let newval = (sums[cluster][vardex] / denoms[cluster][vardex]).max(0.00001).min(0.99999);
                //if newval.is_nan() {println!("PUKE nan works");}
                //println!("{}, {}, {}", sums[cluster][vardex], denoms[cluster][vardex], newval);
                change += (cluster_centers[cluster][*variant] - newval).abs();
                cluster_centers[cluster][*variant] = newval;
            }
        }
    }
    //(cluster_centers, sums, denoms, new_variants, cluster_support, touched) // return cluster centers and memory reuse things
}

fn gather_info(variants: Vec<Vec<i32>>, all_barcodes: FnvHashMap<i32, Vec<i32>>, info: Vec<(String, i32)>) {
    eprintln!("chrom1\tpos1\tchrom2\tchrom2\tbc1\tbc2\tshared");
    for var1dex in (0..variants.len()).step_by(5000) {
        let mut barcodes: FnvHashSet<i32> = FnvHashSet::default();
        for bc in &variants[var1dex] { barcodes.insert(bc.abs()); }
        let mut vars: FnvHashSet<usize> = FnvHashSet::default();
        for bc in &variants[var1dex] {
            let bc = bc.abs();
            for var2dex in all_barcodes.get(&bc).unwrap() {   //(var1dex+1)..variants.len() {
                let var2dex = var2dex.abs() as usize;
                if var2dex == var1dex { continue; }
                vars.insert(var2dex);
            }
        }
        for var2dex in vars {
            let mut count = 0;
            for bc2 in &variants[var2dex] {
                if barcodes.contains(&bc2.abs()) { count += 1; }
            }
            if count > 3 {
                let (chrom1, pos1) = &info[var1dex];
                let (chrom2, pos2) = &info[var2dex];
                eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{}",chrom1, pos1, chrom2, pos2, variants[var1dex].len(), variants[var2dex].len(), count);
            }
        }
    }
}


struct VariantRecruiter {
    variants_visited: BitVec,
    current_molecules_variants: FnvHashMap<i32, FnvHashSet<i32>>,
    current_variant_molecules: FnvHashMap<i32, FnvHashSet<i32>>,
    molecule_last_seen: FnvHashMap<i32, u32>,
    potential_variants: BTreeSet<(u32, i32)>,
    potential_variant_counts: FnvHashMap<i32, u32>,
    iteration: u32,
    variants: Vec<Vec<i32>>,
    molecules: FnvHashMap<i32, Vec<i32>>,
    rng: rand::rngs::ThreadRng,
}

impl VariantRecruiter { 
    fn new(variants: Vec<Vec<i32>>, molecules: FnvHashMap<i32, Vec<i32>>) -> VariantRecruiter {
        VariantRecruiter{
            variants_visited: BitVec::from_elem(variants.len(), false),
            current_molecules_variants: FnvHashMap::default(),
            current_variant_molecules: FnvHashMap::default(),
            molecule_last_seen: FnvHashMap::default(),
            potential_variants: BTreeSet::new(),
            potential_variant_counts: FnvHashMap::default(),
            iteration: 0,
            variants: variants,
            molecules: molecules,
            rng: rand::thread_rng(),
        }
    }

    fn add_variant(&mut self, variant: i32, connections: u32) {
        let variant = variant.abs();
        self.variants_visited.set(variant as usize, true);
        //println!("Adding variant {} with {} connections to current molecule set", variant, connections);
        let mut mols: FnvHashSet<i32> = FnvHashSet::default();
        for mol in &self.variants[variant as usize] {
            mols.insert(*mol);
            let mol = mol.abs();
            self.iteration += 1;
            self.molecule_last_seen.insert(mol, self.iteration);
            // for each new molecule we touch, go through their variants and add or update them to potential_variants and potential_variant_scores
            if !self.current_molecules_variants.contains_key(&mol) {
                for var in self.molecules.get(&mol).unwrap() {
                    let var = var.abs();
                    if self.variants_visited.get(var as usize).unwrap() { continue; }
                    let count = self.potential_variant_counts.entry(var).or_insert(0);
                    self.potential_variants.remove(&(*count, var));
                    *count += 1;
                    self.potential_variants.insert((*count, var));
                }
            }
            // add the variant to the current molecule sets
            let vars = self.current_molecules_variants.entry(mol).or_insert_with(|| FnvHashSet::default());
            vars.insert(variant);
        }
        // add molecules to variant
        self.current_variant_molecules.insert(variant, mols);
        //let mut molremove: Vec<i32> = Vec::new();
        //for (mol, last_seen) in self.molecule_last_seen.iter() {
        //    if self.iteration - last_seen > 200 {
        //        molremove.push(*mol);
        //    }
        //}
        //println!("{}\t{}\t{}\t{}",variant, connections, self.potential_variant_counts.len(),molremove.len());
        //for mol in molremove {
        //    self.remove_molecule(mol);
        //}
    }

    fn remove_molecule(&mut self, molecule: i32) {
        // so we will look at each variant in self.current_molecules_variants
        for var in self.current_molecules_variants.get(&molecule).unwrap().iter() {
            // and remove that molecule from self.current_variant_molecules
            self.current_variant_molecules.get_mut(var).unwrap().remove(&molecule); // what should I do when this goes to 0?
            if self.current_variant_molecules.get(var).unwrap().len() == 0 {
                self.current_variant_molecules.remove(var);
            }
        }
        // and we will go through the variants in self.molecules[molecule] and downgrade or remove them from self.potential_variants_counts and self.potential_variants
        for var in &self.molecules[&molecule] {
            let var = var.abs();
            if self.variants_visited.get(var as usize).unwrap() { continue; }
            let count = match self.potential_variant_counts.get(&var) {
                Some(x) => *x,
                None => 0,
            };
            if count > 0 {
                self.potential_variants.remove(&(count, var));
            }
            if count > 1 {
                self.potential_variants.insert((count-1, var));
                self.potential_variant_counts.insert(var, count-1);
            } else {
                self.potential_variant_counts.remove(&var);
            }
        }
        // then we will remove the molecule from self.current_molecules_variants
        self.current_molecules_variants.remove(&molecule);
        self.molecule_last_seen.remove(&molecule);
    }

    fn get_next_variant(&mut self, min: u32) -> Option<(i32, u32)> {
        //println!("choosing next variant");
        if self.potential_variant_counts.len() == 0 {
            return Some((self.rng.gen_range(0, self.variants.len() as i32), 0))
        } else {
            let (count, variant) = self.potential_variants.iter().next_back().unwrap();
            let (count, variant) = (*count, *variant);
            self.potential_variants.remove(&(count, variant));
            self.potential_variant_counts.remove(&variant);
            if count >= min {
                return Some((variant, count));
            }
        }
        println!("ending. why? I have visited {} variants out of {} variants", self.variants_visited.iter().filter(|x| *x).count(), self.variants.len());
        println!("{} variants in potential variants.", self.potential_variant_counts.len());
        let mut largest = 0;
        for (key, val) in self.potential_variant_counts.iter() {
            if *val > largest { largest = *val; }
        }
        println!("largest connectivity left in potential variants is {}",largest);
        None
    }
}

fn experiment(variants: Vec<Vec<i32>>, molecules: FnvHashMap<i32, Vec<i32>>) {
    let mut variant_window = VariantRecruiter::new(variants, molecules);
    //pick a starting variant
    loop {
        let (variant, connections) = match variant_window.get_next_variant(5) {
            Some((x,y)) => (x,y),
            None => break,
        };
        // add our new variant
        variant_window.add_variant(variant, connections);
    }
}

fn load_variants(params: &Params, whitelist: &FnvHashMap<String, i32>) -> 
        (Vec<Vec<i32>>, FnvHashMap<i32, Vec<i32>>, FnvHashMap<i32, (String, i32, String, String, String)>, FnvHashMap<i32,String>) {
    let reader = File::open(&params.variants).expect("cannot open file");
    let reader = BufReader::new(reader);
    let mut var_id: i32 = 1;
    let mut variants: Vec<Vec<i32>> = Vec::new();
    let mut molecules: FnvHashMap<i32, Vec<i32>> = FnvHashMap::default();
    let mut info: FnvHashMap<i32, (String, i32, String, String, String)> = FnvHashMap::default();
    let mut phasing: FnvHashMap<i32, String> = FnvHashMap::default();
    for line in reader.lines() {
        let line = line.expect("could not read line");
        if line.starts_with('#') { continue; }
        if var_id % 10000 == 0 { 
            //println!("{} variants processed",var_id); 
            if var_id >= 500000 { break; }
        }
        let tokens: Vec<&str> = line.trim().split_whitespace().collect();
        if tokens[6] != "PASS" { continue; }
        let chrom = tokens[0].to_string();
        let pos: i32 = tokens[1].to_string().parse::<i32>().unwrap();
        let reference = tokens[3].to_string();
        let alt = tokens[4].to_string();
        if !(reference.len() == 1 && alt.len() == 1) {
            continue;
        }

        let mut molsupport: Vec<i32> = Vec::new();
        let format_tokens: Vec<&str> = tokens[9].split(':').collect();
        let format_indexes: Vec<&str> = tokens[8].split(':').collect();
        let mut bx_index = 0;
        let mut gt_index = 0;
        let mut ps_index = 0;
        for (index, form) in format_indexes.iter().enumerate() {
            if form == &"BX" {
                bx_index = index;
            }
            if form == &"GT" {
                gt_index = index;
            }
            if form == &"PS" {
                ps_index = index;
            }
        }
        let ps = format_tokens[ps_index];
        let gt = format_tokens[gt_index];
        //println!("{}",gt);
        if gt == "0|0" || gt == "1|1" || gt == "0/0" || gt == "1/1" {
            //println!("continuing");
            continue;
        }
        phasing.insert(var_id, gt.to_string());
        let barcode_support: Vec<&str> = format_tokens[bx_index].split(',').collect();
        if barcode_support.len() < 2 {
            println!("assert going to fail on line: {}",line);
        }
        assert!(barcode_support.len() > 1);
        let reftokens: Vec<&str> = barcode_support[0].split(';').collect();
        let alttokens: Vec<&str> = barcode_support[1].split(';').collect();
        let mut refids: FnvHashSet<i32> = FnvHashSet::default();
        //let altids: FnvHashSet<MolId> = FnvHashSet::default();
        
        for bc in reftokens {
            let bc: Vec<&str> = bc.split('-').collect();
            let bc = bc[0];
            if bc.len() == 0 { continue; }
            let mol_id = whitelist.get(bc).unwrap();
            refids.insert(*mol_id);
            let varsupport = molecules.entry(*mol_id).or_insert(Vec::new());
            varsupport.push(var_id);
            molsupport.push(*mol_id);
        }
        for bc in alttokens {
            let bc: Vec<&str> = bc.split('-').collect();
            let bc = bc[0];
            if bc.len() == 0 { continue; }
            let mol_id = whitelist.get(bc).unwrap();
            //if refids.contains_key(id) {
            //    // do something
            //}
            let varsupport = molecules.entry(*mol_id).or_insert(Vec::new());
            varsupport.push(- var_id);
            molsupport.push(- *mol_id);
        }
        variants.push(molsupport);
        info.insert(var_id, (chrom, pos, reference, alt, ps.to_string()));
        var_id += 1;
    }
    println!("reducing to good molecules");
    molecules.retain(|_key, value| value.len() > 5 && value.len() < 5000); // remember to get rid of this, replace with > 1 or > 2 or 3
    for vardex in 0..variants.len() {
        let mut mol: Vec<i32> = Vec::new();
        for moldex in 0..variants[vardex].len() {
            if molecules.contains_key(&variants[vardex][moldex].abs()) {
                mol.push(variants[vardex][moldex]);
            }
        }
        variants[vardex] = mol;
    }
    (variants, molecules, info, phasing)
}

fn load_whitelist(params: &Params) -> FnvHashMap<String, i32> {
    let mut whitelist_map: FnvHashMap<String, i32> = FnvHashMap::default();
    let reader = File::open(&params.whitelist).expect("cannot open file");
    let reader = BufReader::new(reader);
    for (index, line) in reader.lines().enumerate() {
        let bcline = line.expect("could not read line").trim().to_string();
        whitelist_map.insert(bcline, (index + 1) as i32);
    }
    whitelist_map
}

#[derive(Clone)]
struct Params {
    variants: String,
    whitelist: String,
    ploidy: usize,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let variants = params.value_of("variants").unwrap();
    let ploidy = params.value_of("ploidy").unwrap_or("2");
    let ploidy: usize = ploidy.to_string().parse::<usize>().unwrap();
    let whitelist = params.value_of("whitelist").unwrap();
    Params {
        variants: variants.to_string(),
        ploidy: ploidy,
        whitelist: whitelist.to_string(),
    }
}
