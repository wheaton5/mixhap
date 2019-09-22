#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate rand;
extern crate bit_vec;
extern crate byteorder;

use byteorder::{ByteOrder, LittleEndian};

//use cpuprofiler::PROFILER;
use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::f32;

use clap::{App};

use std::io::BufReader;
use std::io::BufRead;
use std::io::Read;
use std::fs::File;

use hashbrown::{HashMap,HashSet};

use std::collections::BTreeSet;
use bit_vec::BitVec;

fn main() {
    let params = load_params();
    //let whitelist = load_whitelist(&params); // map from barcode string to index in whitelist
    //println!("whitelist loaded");
    //let (variants, barcodes, info, phasing) = load_variants(&params, &whitelist); 
    let (variants, barcodes) = load_molecule_kmers(&params);
    let phasing = None;
    let info = None;
    //gather_info(variants, barcodes, info);
        // list from variant to list of molecule ids (negative means alt support, positive means ref support)
                                // map from molecules to list of variants (negative means alt support, positive means ref support)
    //println!("variant support loaded");
    //println!("testing, I have this many molecules {}", barcodes.len());
    //let (variants, molecules) = extract_molecules(variants, barcodes);
    //go_ahead(&variants, &barcodes);
    phuzzy_phaser_master(barcodes, variants, params.ploidy, phasing, info);
}





fn phuzzy_phaser_master(molecules: HashMap<i32, Vec<i32>>, variants: Vec<Vec<i32>>, 
        num_clusters: usize, phasing: Option<HashMap<i32, String>>, 
        info: Option<HashMap<i32, (String, i32, String, String, String)>>) {
    let mut cluster_support: Vec<f32> = Vec::new(); // to be reused
    for _cluster in 0..num_clusters { cluster_support.push(0.0); } //initialize
    let mut sums: Vec<Vec<f32>> = Vec::new(); let mut denoms: Vec<Vec<f32>> = Vec::new();
    let vars_per_step = 20;
    for cluster in 0..num_clusters {
        denoms.push(Vec::new()); sums.push(Vec::new());
        for _ in 0..vars_per_step { denoms[cluster].push(0.0); sums[cluster].push(0.0); }
    }
    let mut cluster_centers: Vec<Vec<f32>> = Vec::new();
    let seed: [u8; 32] = [4; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
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

    let mut first_varset: Vec<i32> = Vec::new();
    let directions: [i32; 2] = [1, -1];
    let mut variant_window = VariantRecruiter::new(&variants, &molecules);
    println!("iteration\tvardex\tchrom\tpos\tref\talt\tLR_ps\tLR_gt\tconnections\thap1\thap2");
    for direction in &directions { 
        let mut total_variants: i32 = 0;
        if *direction == -1 && total_variants == 0 {
            variant_window.reset_except_visited();
            for vardex in &first_varset {
                variant_window.add_variant(&variants, &molecules, *vardex);
            }
        }
        //pick 20 starting variant
        'phaseblockloop: loop {
            let mut new_variants: Vec<usize> = Vec::new();
            let mut current_molecules: Vec<Vec<i32>> = Vec::new();
            let mut used_mols: HashSet<i32> = HashSet::default();
            let mut mol_connections: Vec<u32> = Vec::new();
            let mut new_iterations: Vec<i32> = Vec::new();
            for _i in 0..vars_per_step {
                new_iterations.push(total_variants);
                total_variants += direction;
                let (variant, connections) = match variant_window.get_next_variant(&variants, &molecules, 20) {
                    Some((x,y)) => (x,y),
                    None => { 
                        //println!("breaking phaseblock loop because i didnt get next variant"); 
                        break 'phaseblockloop },
                };
                // add our new variant
                //println!("{}\t{}",variant, connections);
                variant_window.add_variant(&variants, &molecules, variant);
                if *direction == 1 && total_variants < 10 { first_varset.push(variant); }
                mol_connections.push(connections);
                let vardex = variant.abs() as usize;
                new_variants.push(vardex);
                for molecule in &variants[vardex] {
                    let moldex = molecule.abs();
                    if !used_mols.contains(&moldex) {
                        let mut mol: Vec<i32> = Vec::new();
                        for var in molecules.get(&moldex).unwrap() {
                            let vardex = var.abs();
                            if variant_window.current_variant_molecules.contains_key(&vardex) { // TODO might want to remove
                                mol.push(*var);
                            }
                        }
                        current_molecules.push(mol);
                        used_mols.insert(moldex);
                    } 
                }     
            }
            new_variants.sort(); 
            for variant in &new_variants { touched_vars[*variant] = true; }
            next_cluster_step(current_molecules, &mut cluster_centers, &new_variants, &mut cluster_support, &touched_vars, &mut denoms, &mut sums);
            for (index, vardex) in new_variants.iter().enumerate() {
                let mut chrom: String = "none".to_string();
                let mut pos: i32 = 0; 
                let mut reference: String = "ref".to_string();
                let mut alt: String = "alt".to_string();
                let mut ps: String = "no_ps".to_string();
                if let Some(some_info) = &info {
                    let (chrom1, pos1, reference1, alt1, ps1) = some_info.get(&(*vardex as i32)).unwrap();
                    chrom = chrom1.to_string(); pos = *pos1; reference = reference1.to_string(); alt = alt1.to_string(); ps = ps1.to_string();
                }
                let mut phasing_str: String = "./.".to_string();
                if let Some(phasing_info) = &phasing {
                    phasing_str = phasing_info.get(&(*vardex as i32)).unwrap().to_string();
                }
                println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}", new_iterations[index], vardex, chrom, pos, reference, alt, ps, 
                    phasing_str, mol_connections[index],cluster_centers[0][*vardex],cluster_centers[1][*vardex]);
            }
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
    let mut vardex; let mut moldex; let mut varval; let mut molval;
    let mut iteration = 0;
    let mut change = 1.0;
    while change > 0.0001 {
        iteration += 1; // may want to limit iterations eventually, keeping for now
        // reset reused memory to prior expectation, think of this like a pseudocount and a way to avoid div by 0
        for cluster in 0..denoms.len() {
            for vardex in 0..denoms[cluster].len() {
                denoms[cluster][vardex] = 0.001;
                sums[cluster][vardex] = 0.001*cluster_centers[cluster][new_variants[vardex]];
            }
        } 
        // Expectation
        for (_mol_num, molecule_variants) in molecules.iter().enumerate() {
            //let mut sum = 0.0;
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
                    //assert!(!cluster_support[cluster].is_infinite(), "infinite value, cluster center {}", cluster_centers[cluster][var]);
                    //if cluster_support[cluster].is_infinite() {
                    //    println!("infinite value, cluster_center {}",cluster_centers[cluster][var]);
                    //    assert!(false);
                    //}
                }
            }
            let sum = log_sum_exp(cluster_support);
            for cluster in 0..cluster_centers.len() { 
                cluster_support[cluster] -= sum; 
                cluster_support[cluster] = cluster_support[cluster].exp(); 
            }
            vardex = 0; moldex = 0; varval = new_variants[vardex]; molval = molecule_variants[moldex].abs() as usize;
            loop {
                if varval > molval { 
                    moldex += 1; if moldex >= molecule_variants.len() { break; } molval = molecule_variants[moldex].abs() as usize; 
                } else if molval > varval {
                    vardex += 1; if vardex >= new_variants.len() { break; } varval = new_variants[vardex];
                } else {
                    // actually update things
                    for cluster in 0..cluster_centers.len() {
                        denoms[cluster][vardex] += cluster_support[cluster];
                        if molecule_variants[moldex] > 0 {
                            sums[cluster][vardex] += cluster_support[cluster]; // 1*cluster_support[0]
                        } // else add 0*cluster_support[cluster]
                    }
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
                change += (cluster_centers[cluster][*variant] - newval).abs();
                cluster_centers[cluster][*variant] = newval;
            }
        }
    }
}

#[allow(dead_code)]
struct PhaseBlock {
    variant_indices: Vec<usize>,
    cluster_centers: Vec<Vec<f32>>,
    visited: BitVec,
}

impl PhaseBlock {
    #[allow(dead_code)]
    fn new(centers: &Vec<Vec<f32>>, visited: &BitVec) -> PhaseBlock {
        let mut cluster_centers: Vec<Vec<f32>> = Vec::new();
        let mut variant_indices: Vec<usize> = Vec::new();
        let mut variants_visited: BitVec = BitVec::from_elem(visited.len(), false);
        for index in 0..visited.len() {
            if visited.get(index).unwrap() {
                variant_indices.push(index);
                variants_visited.set(index, true);
            }
        }
        for cluster in 0..centers.len() {
            let mut clust_center: Vec<f32> = Vec::new();
            for index in 0..visited.len() {
                if visited.get(index).unwrap() {
                    clust_center.push(centers[cluster][index]);
                }
            }
            cluster_centers.push(clust_center);
        }
        PhaseBlock{
            cluster_centers: cluster_centers,
            variant_indices: variant_indices,
            visited: variants_visited,
        }
    }
}

struct VariantRecruiter {
    variants_visited: BitVec,
    current_molecules_variants: HashMap<i32, HashSet<i32>>,
    current_variant_molecules: HashMap<i32, HashSet<i32>>,
    molecule_last_seen: HashMap<i32, u32>,
    potential_variants: BTreeSet<(u32, i32)>,
    potential_variant_counts: HashMap<i32, u32>,
    iteration: u32,
    rng: StdRng,
}

impl VariantRecruiter { 
    fn new(variants: &Vec<Vec<i32>>, _molecules: &HashMap<i32, Vec<i32>>) -> VariantRecruiter {
        let seed: [u8; 32] = [1; 32];
        let rng: StdRng = SeedableRng::from_seed(seed);
        VariantRecruiter{
            variants_visited: BitVec::from_elem(variants.len() + 1, false),
            current_molecules_variants: HashMap::default(),
            current_variant_molecules: HashMap::default(),
            molecule_last_seen: HashMap::default(),
            potential_variants: BTreeSet::new(),
            potential_variant_counts: HashMap::default(),
            iteration: 0,
            rng: rng,
        }
    }

    fn reset_except_visited(&mut self) {
        self.current_molecules_variants = HashMap::default();
        self.current_variant_molecules = HashMap::default();
        self.molecule_last_seen = HashMap::default();
        self.potential_variants = BTreeSet::new();
        self.potential_variant_counts = HashMap::default();
        self.iteration = 0;
    }

    fn add_variant(&mut self, variants: &Vec<Vec<i32>>, molecules: &HashMap<i32, Vec<i32>>, variant: i32) {
        let variant = variant.abs();
        self.variants_visited.set(variant as usize, true);
        //println!("Adding variant {} with {} connections to current molecule set", variant, connections);
        let mut mols: HashSet<i32> = HashSet::default();
        self.iteration += 1;
        for mol in &variants[variant as usize] {
            mols.insert(*mol);
            let mol = mol.abs();
            self.molecule_last_seen.insert(mol, self.iteration);
            // for each new molecule we touch, go through their variants and add or update them to potential_variants and potential_variant_scores
            if !self.current_molecules_variants.contains_key(&mol) {
                for var in molecules.get(&mol).unwrap() {
                    let var = var.abs();
                    assert!((var as usize) < self.variants_visited.len(), "{} < {}",var, self.variants_visited.len());
                    if self.variants_visited.get(var as usize).unwrap() { continue; }
                    let count = self.potential_variant_counts.entry(var).or_insert(0);
                    self.potential_variants.remove(&(*count, var));
                    //println!("{}", (1000.0/(variants[var as usize].len() as f32)).round() as u32);
                    *count += 1000/(variants[(var - 1) as usize].len() as u32);
                    self.potential_variants.insert((*count, var));
                }
            }
            // add the variant to the current molecule sets
            let vars = self.current_molecules_variants.entry(mol).or_insert_with(|| HashSet::default());
            vars.insert(variant);
        }
        // add molecules to variant
        self.current_variant_molecules.insert(variant, mols);
        if self.iteration % 20 == 0 {
            let mut molremove: Vec<i32> = Vec::new(); //TODO might want to remove
            for (mol, last_seen) in self.molecule_last_seen.iter() {
                if self.iteration - last_seen > 100 {
                    molremove.push(*mol);
                }
            }
            for mol in molremove.iter() {
                self.remove_molecule(variants, molecules, *mol);
            }
            // FOR DEBUG let (count, var) = self.potential_variants.iter().next_back().unwrap();
            //println!("removing {}, potential variants {}, new top connection {}", molremove.len(), self.potential_variants.len(), count);
        }
    }

    fn remove_molecule(&mut self, variants: &Vec<Vec<i32>>, molecules: &HashMap<i32, Vec<i32>>, molecule: i32) {
        // so we will look at each variant in self.current_molecules_variants
        for var in self.current_molecules_variants.get(&molecule).unwrap().iter() {
            // and remove that molecule from self.current_variant_molecules
            self.current_variant_molecules.get_mut(var).unwrap().remove(&molecule); // what should I do when this goes to 0?
            if self.current_variant_molecules.get(var).unwrap().len() == 0 {
                self.current_variant_molecules.remove(var);
            }
        }
        // and we will go through the variants in self.molecules[molecule] and downgrade or remove them from self.potential_variants_counts and self.potential_variants
        for var in &molecules[&molecule] {
            let var = var.abs();
            if self.variants_visited.get(var as usize).unwrap() { continue; }
            let count = match self.potential_variant_counts.get(&var) {
                Some(x) => *x,
                None => 0,
            };
            if count > 0 {
                self.potential_variants.remove(&(count, var));
            }
            let incr: u32 =  1000/(variants[(var - 1) as usize].len() as u32);
            if count > incr {
                self.potential_variants.insert((count-incr, var));
                self.potential_variant_counts.insert(var, count-incr);
            } else {
                self.potential_variant_counts.remove(&var);
            }
        }
        // then we will remove the molecule from self.current_molecules_variants
        self.current_molecules_variants.remove(&molecule);
        self.molecule_last_seen.remove(&molecule);
    }

    fn get_next_variant(&mut self, variants: &Vec<Vec<i32>>, _molecules: &HashMap<i32, Vec<i32>>, min: u32) -> Option<(i32, u32)> {
        //println!("choosing next variant");
        let min = min*10;
        if self.potential_variant_counts.len() == 0 {
            return Some((self.rng.gen_range(1, variants.len() as i32), 0))
        } else {    
            loop {
                if self.potential_variant_counts.len() == 0 { break; }
                let (count, variant) = self.potential_variants.iter().next_back().unwrap();
                let (count, variant) = (*count , *variant);
                self.potential_variants.remove(&(count, variant));
                self.potential_variant_counts.remove(&variant);
                if count >= min && (count as f32)/1000.0*(variants[(variant-1) as usize].len() as f32) > 10.0 {
                    return Some((variant, count));
                } //else {
                    //println!("rejecting {} because {} < {} or {} <= 10.0", variant, count, min, (count as f32)/1000.0*(variants[(variant-1) as usize].len() as f32));
                //}
            }
        }
        //println!("ending. why? I have visited {} variants out of {} variants", self.variants_visited.iter().filter(|x| *x).count(), variants.len());
        //println!("{} variants in potential variants.", self.potential_variant_counts.len());
        //let mut largest = 0;
        //for (key, val) in self.potential_variant_counts.iter() {
        //    if *val > largest { largest = *val; }
        //}
        //println!("largest connectivity left in potential variants is {}",largest);
        None
    }
}

fn eat_i32<R: Read>(reader: &mut BufReader<R>, buf: &mut [u8;4]) -> Option<i32> {
    let bytes = reader.read(buf).expect("could not read binary molecules file");
    if bytes < 4 { return None; }
    Some(LittleEndian::read_i32(buf))
}

fn load_molecule_kmers(params: &Params) -> (Vec<Vec<i32>>, HashMap<i32, Vec<i32>>){
    let f = File::open(params.variants.to_string()).expect("Unable to open file");
    let mut reader = BufReader::new(f);
    let mut variants: HashMap<i32, HashSet<i32>> = HashMap::new(); // map from variant_id to list of molecule_ids
    let mut molecules: HashMap<i32, HashSet<i32>> = HashMap::new(); //map from molecule_id to list of variant_ids
    let mut hic_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut bufi32 = [0u8; 4];
    let mut max_var = 0;
    let mut max_molid = 0;
    loop {
        if let Some(barcode_id) = eat_i32(&mut reader, &mut bufi32) {
            if let Some(kmer_id) = eat_i32(&mut reader, &mut bufi32) {
                if barcode_id == 0 && kmer_id == 0 { break; } // only dealing with 10x data so far
                let bc_vars = molecules.entry(barcode_id).or_insert(HashSet::new());
                bc_vars.insert(kmer_id);
                let var_bcs = variants.entry(kmer_id).or_insert(HashSet::new());
                if kmer_id < 0 { var_bcs.insert(-barcode_id); } else { var_bcs.insert(barcode_id); }
                if kmer_id.abs() > max_var { max_var = kmer_id.abs() }
                if barcode_id > max_molid { max_molid = barcode_id; }
            } else { break; }
        } else { break; }
    }
    let mut mol_id = max_molid + 1;
    loop { // now deal with hic data, format is 2 i32s, if i get to 2 0's we are done
        if let Some(kmer_id1) = eat_i32(&mut reader, &mut bufi32) {
            if let Some(kmer_id2) = eat_i32(&mut reader, &mut bufi32) {
                if kmer_id1 == 0 { break; }
                let mut blah: HashSet<i32> = HashSet::new();
                blah.insert(kmer_id1); blah.insert(kmer_id2);
                hic_molecules.insert(mol_id, blah);
                mol_id += 1;
            } else { break; }
        } else { break; }
    }
    mol_id += 1; //
    'outer: loop { // ok here we go again. Pacbio/longread data. Format is i32s until you hit a zero. when you hit two zeros you are done
        let mut vars: HashSet<i32> = HashSet::new();
        loop {
            if let Some(kmer_id) = eat_i32(&mut reader, &mut bufi32) {
                if kmer_id == 0 { break; }
                vars.insert(kmer_id);
                let var_bcs = variants.entry(kmer_id).or_insert(HashSet::new());
                if kmer_id < 0 { var_bcs.insert(-mol_id); } else { var_bcs.insert(mol_id); }
                if kmer_id.abs() > max_var { max_var = kmer_id.abs(); }
            } else { break 'outer; }
        }
        if vars.len() > 1 { molecules.insert(mol_id, vars); mol_id += 1; }
    }

    let mut vars: Vec<Vec<i32>> = Vec::new();
    for _i in 0..max_var {
        vars.push(Vec::new());
    }
    let mut mols: HashMap<i32, Vec<i32>> = HashMap::new();
    for (mol_id, kmer_ids) in molecules {
        let mut kmer_ids_sorted = kmer_ids.iter().cloned().collect::<Vec<i32>>();
        kmer_ids_sorted.sort_by(| a, b | a.abs().cmp(&b.abs()));
        for var in kmer_ids_sorted.iter() {
            if var < &0 {
                vars[(var.abs() - 1) as usize].push(-mol_id);
            } else { vars[(var - 1) as usize].push(mol_id) }
        }
        mols.insert(mol_id, kmer_ids_sorted);
    }
    for i in 0..max_var {
        vars[i as usize].sort_by(| a, b | a.abs().cmp(&b.abs()));
    }
    (vars, mols)
}

#[allow(dead_code)]
fn load_variants(params: &Params, whitelist: &HashMap<String, i32>) -> 
        (Vec<Vec<i32>>, HashMap<i32, Vec<i32>>, HashMap<i32, (String, i32, String, String, String)>, HashMap<i32,String>) {
    let reader = File::open(&params.variants).expect("cannot open file");
    let reader = BufReader::new(reader);
    let mut var_id: i32 = 1;
    let mut variants: Vec<Vec<i32>> = Vec::new();
    let mut molecules: HashMap<i32, Vec<i32>> = HashMap::default();
    let mut info: HashMap<i32, (String, i32, String, String, String)> = HashMap::default();
    let mut phasing: HashMap<i32, String> = HashMap::default();
    let mut index = 0;
    for line in reader.lines() {
        let line = line.expect("could not read line");
        if line.starts_with('#') { continue; }
        //if index < 50000 { index += 1; continue; }
        //if index % 10000 == 0 { 
            //println!("{} variants processed",var_id); 
            
            if index >= 500000 { break; }
        //}
        index += 1;
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
        let mut refids: HashSet<i32> = HashSet::default();
        //let altids: HashSet<MolId> = HashSet::default();
        
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
    println!("{} good molecules",molecules.len());
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

#[allow(dead_code)]
fn load_whitelist(params: &Params) -> HashMap<String, i32> {
    let mut whitelist_map: HashMap<String, i32> = HashMap::default();
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

// test code for sparse graph visualization
#[allow(dead_code)]
fn go_ahead(_variants: &Vec<Vec<i32>>, molecules: &HashMap<i32, Vec<i32>>) {
    let mut variant_connections: HashMap<(i32, i32), usize> = HashMap::default();
    for (_mol, vars) in molecules.iter() {
        for i in 0..vars.len() {
            for j in (i+1)..vars.len() {
                let v1 = vars[i].abs();
                let v2 = vars[j].abs();
                let count = variant_connections.entry((v1,v2)).or_insert(0);
                *count += 1;
            }
        }
    }
    println!("graph {{");
    for ((v1,v2), c) in variant_connections.iter() {
        if c > &20 && v1 % 5 == 0 && v2 % 5 == 0 {
            println!("{} -- {} [weight=\"{}\"];", v1, v2, c);
        }
    }
    println!("}}");
}

//Richard's 10x molecule detection algorithm has some problems.
#[allow(dead_code)]
fn extract_molecules(variants: Vec<Vec<i32>>, barcodes: HashMap<i32, Vec<i32>>) -> 
        (Vec<Vec<i32>>, HashMap<i32, Vec<i32>>) { //may also want to make a Vec<Vec<i32>> mapping variants to other variants that have a > threshold link
    //PROFILER.lock().unwrap().start("./my-prof.profile");
    let mut new_variants: Vec<Vec<i32>> = Vec::new();
    for _i in 0..(variants.len() + 1) { new_variants.push(Vec::new()); }
    let mut molecules: HashMap<i32, Vec<i32>> = HashMap::default();
    let mut molecule_number: i32 = 0;
    let mut barcode_molecules: HashMap<i32, Vec<i32>> = HashMap::default(); // mapping from barcode to list of molecules
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
    //PROFILER.lock().unwrap().stop();
    (new_variants, molecules)
}