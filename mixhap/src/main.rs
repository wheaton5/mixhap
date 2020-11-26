#![feature(map_first_last)]
#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate rand;
extern crate flate2; 
extern crate byteorder;
extern crate rayon;
extern crate itertools;

use flate2::read::GzDecoder;
use rayon::prelude::*;
use itertools::izip;
use std::fs;
use std::str;
use std::error::Error;
use std::fs::OpenOptions;
// Change the alias to `Box<error::Error>`.
type Result<T> = std::result::Result<T, Box<dyn Error>>;

use byteorder::{ByteOrder, LittleEndian};
use std::io::{BufWriter, Write};

//use cpuprofiler::PROFILER;
use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand::thread_rng;
use rand::seq::SliceRandom;
use std::f32;
use std::convert::TryInto;
use bit_vec::BitVec;

use clap::{App};

use std::io::BufReader;
use std::io::BufRead;
use std::io::Read;
use std::fs::File;

use hashbrown::{HashMap,HashSet};
use std::collections::VecDeque;

use std::collections::BTreeSet;
//use rayon::prelude::*;

fn main() {
    let params = load_params();
    //let whitelist = load_whitelist(&params); // map from barcode string to index in whitelist
    //println!("whitelist loaded");
    //let (variants, barcodes, info, phasing) = load_variants(&params, &whitelist); 
    //eprintln!("loading kmers");
    let kmers = load_kmers(&params.kmers);
    let crib = match &params.crib {
        Some(x) => {
            //eprintln!("loading crib");
            Some(load_crib(&x, &kmers))
        },
        None => None,
    };


    //eprintln!("loading molecules");
    let (variants, molecules) = load_molecule_kmers(&params, &kmers);


    /* 
    let data = "%\n%\n";
    let f1 = File::create("/Users/hh5/Documents/datasets/ribosomal_clusters/ref.mtx").expect("Unable to create file");
    let mut f1 = BufWriter::new(f1);
    f1.write_all(data.as_bytes()).expect("Unable to write data");

    let f2 = File::create("/Users/hh5/Documents/datasets/ribosomal_clusters/alt.mtx").expect("Unable to create file");
    let mut f2 = BufWriter::new(f2);
    f2.write_all(data.as_bytes()).expect("Unable to write data");

    
    
    let firstline = format!("{}\t{}\t{}\n",variants.long_read_variants.len(), molecules.long_read_molecules.len(), 100000);
    f1.write_all(firstline.as_bytes()).expect("boo");
    f2.write_all(firstline.as_bytes()).expect("boo");
    let mut variant_index: HashMap<i32, usize> = HashMap::new(); // variant space id to index
    let mut index_to_variants: HashMap<usize, (i32, i32)> = HashMap::new(); // variant index space to tuple of variant space ids 
    let mut variant_index_so_far: usize = 0; 
    let mut molecule_index: HashMap<i32, usize> = HashMap::new();
    let mut molecule_so_far: usize = 1;
    let mut bad_vars: HashSet<i32> = HashSet::new();
    let mut variant_molecules: HashMap<i32, Vec<i32>> = HashMap::new(); // variant space id to vector of molecules in molecule id space
    
    for (mol, vars) in molecules.long_read_molecules.iter() {
        let mol = mol.abs();
        if !molecule_index.contains_key(&mol) {
            molecule_so_far += 1;
            molecule_index.insert(mol, molecule_so_far);
        }
        let mut check: HashSet<i32> = HashSet::new();
        for var in vars.iter() {
            if check.contains(var) {
                //println!("mol {} has var {} twice?",mol, var);
                continue;
            } else if check.contains(&-var) {
                //println!("mol {} has var {} and its RC?", mol, var);
                continue;
            } else if var.abs() % 2 == 0 && (check.contains(&(var.abs() - 1)) || check.contains(&(-(var.abs() - 1)))) {
                //println!("mol {} has var {} and its pair???", mol, var);
                bad_vars.insert(var.abs());
                bad_vars.insert(var.abs()-1);
            } else if var.abs() % 2 == 1 && (check.contains(&(var.abs() + 1)) || check.contains(&(-(var.abs() + 1)))) {
                //println!("mol {} has var {} and its pair!!!", mol, var);
                bad_vars.insert(var.abs());
                bad_vars.insert(var.abs()+1);
            }
            check.insert(*var);
            let var = var.abs();
            
            
            let mut vardex: usize = 0;
            if variant_index.contains_key(&var) {
                vardex = *variant_index.get(&var).unwrap();
            } else {
                variant_index_so_far += 1;
                vardex = variant_index_so_far;
                variant_index.insert(var, vardex);

                if var % 2 == 0 {
                    variant_index.insert(var - 1, vardex);
                    index_to_variants.insert(vardex, (var - 1, var));
                } else {
                    variant_index.insert(var + 1, vardex);
                    index_to_variants.insert(vardex, (var, var + 1));
                }
            }
            let listme = variant_molecules.entry(var).or_insert(Vec::new());
            listme.push(mol);
        }
    }
    variant_index_so_far += 1;
    println!("number of bad vars is {}",bad_vars.len()/2);
    for vardex in 1..variant_index_so_far {
        //println!("vardex {}",vardex);
        let (var1, var2) = index_to_variants.get(&vardex).unwrap();
        if bad_vars.contains(var1) { continue; }
        let mut moldexes: Vec<i32> = Vec::new();
        if !variant_molecules.contains_key(var1) {
            continue;
        }
        for mol in variant_molecules.get(var1).unwrap().iter() {
            let moldex = molecule_index.get(mol).unwrap();
            moldexes.push(*moldex as i32);
        }
        if !variant_molecules.contains_key(var2) {
            continue;
        }
        for mol in variant_molecules.get(var2).unwrap().iter() { 
            let moldex = molecule_index.get(mol).unwrap();
            moldexes.push(-(*moldex as i32));
        }
        moldexes.sort_by(|a, b| a.abs().cmp(&b.abs()));
        for mol in moldexes {
            if mol < 0 {
                let line = format!("{}\t{}\t{}\n", vardex, mol.abs(), 0);
                f1.write_all(line.as_bytes()).expect("boo");
                let line = format!("{}\t{}\t{}\n", vardex, mol.abs(), 1);
                f2.write_all(line.as_bytes()).expect("boo");
            } else {
                let line = format!("{}\t{}\t{}\n", vardex, mol.abs(), 1);
                f1.write_all(line.as_bytes()).expect("boo");
                let line = format!("{}\t{}\t{}\n", vardex, mol.abs(), 0);
                f2.write_all(line.as_bytes()).expect("boo");
            }
        }
    }

    let mut adjacency_list: HashMap<i32, HashMap<i32,usize>> = HashMap::new();
    for (_mol, vars) in molecules.long_read_het_molecule_list.iter() {
        if vars.len() < 2 { continue; }
        for i in 0..(vars.len()-1) {
            //println!("its working?");
            let j = i+1;
            //println!("{}",i);

            //for j in (i+1)..vars.len() {
                let var1 = vars[i];
                let var2 = vars[j];
                {
                let thingy = adjacency_list.entry(var1).or_insert(HashMap::new());
                let count = thingy.entry(var2).or_insert(0);
                *count += 1;
                }
                {
                    let thingy = adjacency_list.entry(-var2).or_insert(HashMap::new());
                    let count = thingy.entry(-var1).or_insert(0);
                    *count += 1;
                }
            //}
        }
    }

    println!("graph {{");

    for (var, varcounts) in adjacency_list.iter() {
        let mut count_vec: Vec<(&i32, &usize)> = varcounts.iter().collect();
        count_vec.sort_by(|a, b| b.1.cmp(a.1));
        //if count_vec.len() > 1 {
        //    println!("var {} best {} count {} second best {} second best count {}", var, count_vec[0].0, count_vec[0].1, count_vec[1].0, count_vec[1].1);
        //} else {
        //    println!("var {} best {} count {}", var, count_vec[0].0, count_vec[0].1);
        //}
        if varcounts.len() == 1 || (*count_vec[0].1 as f32)/((count_vec[0].1 + count_vec[1].1) as f32) > 0.9 {
            println!("\t{} -- {};", kmers.kmers.get(&var.abs()).unwrap(), kmers.kmers.get(&count_vec[0].0.abs()).unwrap());
        }
        
    }
    println!("}}");
    */
   
    let (components, adjacency_list) = sparsembly(&molecules, &variants, &kmers, &crib);

    sparsembly2point0(&variants, &molecules, adjacency_list, crib, &kmers, &params);



    //let assignments = detect_loops(&adjacency_list, &varlist);
    //for (var, assignment) in assignments {
    //    println!("{} = {}",var, assignment);
    //}
    /*
    eprintln!("lets phase some blocks");
    let phasing = make_phase_blocks(&molecules, &variants);

    //eprintln!("creating consistency graph");
    //let consistency_graph = create_consistency_graph(&variants, &molecules);
    eprintln!("assigning long reads");


    //draw_block(1415, &consistency_graph);
    //draw_block(791, &consistency_graph);

    let hom_kmer_assignments = assign_homozygous(&phasing, &variants, &molecules);
    let long_read_assignment = assign_long_reads2(&phasing, &molecules, &hom_kmer_assignments);
    make_fastqs(&params, &long_read_assignment, &phasing.phase_block_sizes, &molecules);
    //let phasing = None;
    //let info = None;
    */

    //gather_info(variants, barcodes, info);
        // list from variant to list of molecule ids (negative means alt support, positive means ref support)
                                // map from molecules to list of variants (negative means alt support, positive means ref support)
    //println!("variant support loaded");
    //println!("testing, I have this many molecules {}", barcodes.len());
    //let (variants, molecules) = extract_molecules(variants, barcodes);
    
    //phuzzy_phaser_master(barcodes, variants, params.ploidy, None, None, kmer_ids_to_kmer);
}

#[derive(Debug)]
struct Status {
    r_h1: u16,
    r_h2: u16,
    r_none: u16,
    r_both: u16,
    a_h1: u16,
    a_h2: u16,
    a_none: u16,
    a_both: u16,
}

impl Status {
    fn new(var: i32, variants: &Variants) -> Status {
        let mut status = Status{
            r_h1: 0,
            r_h2: 0,
            r_none: variants.get_long_read_molecules(var.abs()).count() as u16,
            r_both: 0,
            a_h1: 0,
            a_h2: 0,
            a_none: variants.get_long_read_molecules(pair(var.abs())).count() as u16,
            a_both: 0,
        };
        status
    }
    fn empty() -> Status {
        Status{
            r_h1: 0,
            r_h2: 0,
            r_none: 0,
            r_both: 0,
            a_h1: 0,
            a_h2: 0,
            a_none: 0,
            a_both: 0,
        }
    }
    fn get_status(var: i32, variants: &Variants, molecules: &Molecules, phasing: &HashMap<i32, bool>, kmers: &Kmers) -> Status {
        let mut status = Status::empty();
        let mut refvar = var;
        let mut altvar = pair(var);
        if var % 2 != 0 {
            refvar = pair(var);
            altvar = var;
        }
        for mol in variants.get_long_read_molecules(refvar) {
            let mut hap1 = 0;
            let mut hap2 = 0;
            for var2 in molecules.get_long_read_variants(mol.abs()) {
                if let Some(phase) = phasing.get(&var2.abs()) {
                    if *phase {
                        hap1 += 1;
                    } else {
                        hap2 += 1;
                    }
                }
            }
            if hap1 > 0 && hap2 > 0 {
                //eprintln!("wtf kmers ref both so check kmers\n{}\n{}",
                //    kmers.kmers.get(&refvar).unwrap(), kmers.kmers.get(&altvar).unwrap());
                status.r_both += 1;
            } else if hap1 > 0 && hap2 == 0 {
                status.r_h1 += 1;
            } else if hap2 > 0 && hap1 == 0 {
                status.r_h2 += 1;
            } else {
                status.r_none += 1;
            }
        }
        for mol in variants.get_long_read_molecules(altvar) {
            let mut hap1 = 0;
            let mut hap2 = 0;
            for var2 in molecules.get_long_read_variants(mol.abs()) {
                if let Some(phase) = phasing.get(&var2.abs()) {
                    if *phase {
                        hap1 += 1;
                    } else {
                        hap2 += 1;
                    }
                }
            }
            if hap1 > 0 && hap2 > 0 {
                //eprintln!("wtf kmers alt both so check kmers\n{}\n{}",
                //    kmers.kmers.get(&refvar).unwrap(), kmers.kmers.get(&altvar).unwrap());
                status.a_both += 1;
            } else if hap1 > 0 && hap2 == 0 {
                status.a_h1 += 1;
            } else if hap2 > 0 && hap1 == 0 {
                status.a_h2 += 1;
            } else {
                status.a_none += 1;
            }
        }
        status
    }
}

struct Seeder {
    current_index: usize,
    seeds: Vec<i32>,
    visited: HashSet<i32>,
}

impl Seeder {
    fn next_seed(&mut self) -> Option<i32> {
        let mut seed: Option<i32> = None;
        loop {
            if self.current_index == self.seeds.len() { break; }
            let tmp = self.seeds[self.current_index];
            if !self.visited.contains(&tmp) { 
                seed = Some(tmp);
                self.visited.insert(tmp); self.visited.insert(-tmp); self.visited.insert(pair(tmp)); self.visited.insert(-pair(tmp)); 
                break;
            }
            self.current_index += 1;
        }
        seed
    }
}

fn sparsembly2point0(variants: &Variants, molecules: &Molecules, adjacency_list: HashMap<i32, HashMap<i32, (u16, u16)>>, 
        crib: Option<Crib>, kmers: &Kmers, params: &Params) {
    let mut reverse_edges: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut all_vars: Vec<i32> = Vec::new();
    for (var, varset) in adjacency_list.iter() {
        for (var2, _mols) in varset {
            let backset = reverse_edges.entry(*var2).or_insert(HashSet::new());
            backset.insert(*var);
            //all_vars.push(*var2);
        }
        all_vars.push(*var);
    }

    let mut both = 0;
    let mut one = 0;
    let mut neither = 0;
    for kmer in all_vars.iter() {
        let pair = pair(*kmer);
        match crib.as_ref() {
            Some(cheat) => {
                let mut first = false;
                let mut second = false;
                if let Some((contig, num, position)) = cheat.variants.get(&kmer.abs()) {
                    //eprintln!("STARTVAR ref chr{}\t{}", contig, position); any = true;
                    first = true;
                } 
                if let Some((contig, num, position)) = cheat.variants.get(&pair.abs()) {
                    //eprintln!("STARTVAR alt chr{}\t{}", contig, position); any = true;
                    second = true;
                } 
                if !(first || second) {
                    //eprintln!("STARTVAR  ref/alt not in crib");
                    neither +=1;
                } else if (first && !second) || (!first && second) {
                    one += 1;
                } else if first && second {
                    both += 1;
                }
            }
            None => (),
        }
    }
    eprintln!("one {}\tboth {}\tneither {}", one, both, neither);
    //all_vars.shuffle(&mut thread_rng());


    //let mut visited: HashSet<i32> = HashSet::new();
    let mut seeds: Seeder = Seeder{
        visited: HashSet::new(),
        seeds: all_vars,
        current_index: 0,
    };
    let mut deferred_seed: Option<i32> = None;
    let mut phase_block_number = 0;
    

    let mut is_real_block = false; 
    let mut crib_positions: HashMap<i32, Vec<usize>> = HashMap::new();
    let mut crib_chrom: i32 = 0;
    let mut direction = true;
    let mut phase_blocks: HashMap<usize, VecDeque<(i32, i32)>> = HashMap::new();
    let mut phasing: HashMap<i32, bool> = HashMap::new(); 
    let mut kmer_order: VecDeque<(i32, i32)> = VecDeque::new();
    loop {
        
        let mut startvar;
        //if visited.contains(&startvar) { continue; }
        if let Some(seed) = deferred_seed { 
            direction = false;
            if !is_real_block {
                deferred_seed = None;
                crib_positions = HashMap::new();
                continue; // if we didnt go forward dont attempt to go backwards
            }
            startvar = seed ; 
            deferred_seed = None ; 
            eprintln!("\nDeferred seed var {}", startvar);
        }
        else if let Some(seed) = seeds.next_seed() { 
            direction = true;
            startvar = seed; 
            deferred_seed = Some(-seed); 
            if startvar % 2 != 0 {
                startvar = pair(startvar); // ref
            }
            let startpair = pair(startvar); // alt
            

            if is_real_block {
                let mut order: VecDeque<(i32, i32)> = VecDeque::new();
                for (kmer1, kmer2) in kmer_order.iter() {
                    order.push_back((*kmer1, *kmer2));
                }
                phase_blocks.insert(phase_block_number, order);
                let mut blocks: Vec<(i32, usize, usize)> = Vec::new();
                for (chrom, positions) in crib_positions.iter() {
                    let mut sub_positions: Vec<usize> = Vec::new();
                    let mut last = 0;
                    if positions.len() == 0 {
                        blocks.push((-1,0,0));
                        continue;
                    }
                    
                    for pos in positions {
                        sub_positions.push(*pos);
                    }
                    sub_positions.sort();
                    let mut start = sub_positions[0];
                    for (index, pos) in sub_positions.iter().enumerate() {
                        if last != 0 && *pos - last > 100000 {
                            blocks.push((*chrom, start, last));
                            start = *pos;
                        }
                        last = *pos;
                    }
                    blocks.push((*chrom, start, last));
                }
                if blocks.len() > 1 {
                    eprintln!("MULTIPOSITION phaseblock\t{:?}", blocks);
                }
                for (chrom, start, end) in blocks {
                    eprintln!("PHASEBLOCK Complete. chr{}\t{}-{}\tlength\t{}", chrom, start, end, end-start);
                }
                kmer_order.clear();
                kmer_order.push_back((startvar, startpair));
                
            }
            eprintln!("\n\n\n\nNew phaseblock starting with var {}", startvar);
            is_real_block = false;  // must add at least one kmer pair to block to be a real block
            crib_positions = HashMap::new();
        }
        else { eprintln!("DONE"); break; }
        
        let startpair = pair(startvar); // alt
        

        match crib.as_ref() {
            Some(cheat) => {
                let mut any = false;
                if let Some((contig, num, position)) = cheat.variants.get(&startvar.abs()) {
                    eprintln!("STARTVAR ref{} chr{}\t{}\t{}\t{}", num, contig, position, 
                    kmers.kmer_counts.get(&startvar.abs()).unwrap(), kmers.kmer_counts.get(&startpair.abs()).unwrap()); any = true;
                    crib_chrom = *contig;
                    let poses = crib_positions.entry(*contig).or_insert(Vec::new());
                    poses.push(*position);
                } 
                if let Some((contig, num, position)) = cheat.variants.get(&startpair.abs()) {
                    
                    eprintln!("STARTVAR alt{} chr{}\t{}\t{}\t{}", num, contig, position, 
                    kmers.kmer_counts.get(&startvar.abs()).unwrap(), kmers.kmer_counts.get(&startpair.abs()).unwrap()); any = true;
                    let poses = crib_positions.entry(*contig).or_insert(Vec::new());
                    poses.push(*position);
                   
                } 
                if !any {
                    eprintln!("STARTVAR  ref/alt not in crib");
                }
            }
            None => { eprintln!("STARTVAR NoContig\t{}\t{}", 
                kmers.kmers.get(&startvar).unwrap(), kmers.kmers.get(&startpair).unwrap());
            }
        }
        
        //visited.insert(startvar); visited.insert(-startvar); visited.insert(startpair); visited.insert(-startpair);

        let mut forward_edges: HashMap<i32, [u8; 4]> = HashMap::new();
        
        phasing.insert(startvar.abs(), true);
        phasing.insert(startpair.abs(), false);
        let mut bfs_queue: VecDeque<i32> = VecDeque::new();
        for var1 in [startvar, startpair].iter() {
            //eprintln!("here");
            let var1 = *var1;
            if let Some(varset) = adjacency_list.get(&var1) {
                //eprintln!("varset len {}", varset.len());
                let var1mod2 = var1.abs() % 2;
                let mut var_vec: Vec<(i32, u16, u16)> = Vec::new();
                for (var2, (mols, dist)) in varset {
                    var_vec.push((*var2, *mols, *dist));
                }
                var_vec.sort_by(|a, b| a.2.cmp(&b.2));
                for (var2, _mols, _dist) in var_vec {
                    
                    let mut varhandle = var2.min(pair(var2));
                    if var2 > 0 {
                        varhandle = var2.max(pair(var2));
                    }
                    let var2mod2 = var2.abs() % 2;
                    //eprintln!("var2 {}, varhandle {}, var1mod2 {}, var2mod2 {}", var2, varhandle, var1mod2, var2mod2);
                    if !forward_edges.contains_key(&varhandle) {
                        bfs_queue.push_back(varhandle);
                    }
                    let counts = forward_edges.entry(varhandle).or_insert([0;4]);
                    if var1mod2 == 0 && var2mod2 == 0 {
                        counts[0] += 1;
                    } else if var1mod2 == 1 && var2mod2 == 1 {
                        counts[1] += 1;
                    } else if var1mod2 == 0 && var2mod2 == 1 {
                        counts[2] += 1;
                    } else if var1mod2 == 1 && var2mod2 == 0 {
                        counts[3] += 1;
                    } 
                }
            }
        }
        let mut bfs_iter = 0;
        while bfs_queue.len() != 0 {
            let mut refvar = bfs_queue.pop_front().unwrap(); 
            let mut altvar = pair(refvar);
            let counts = forward_edges.get(&refvar).unwrap();
            let cis = (counts[0] + counts[1]) as f32;
            let trans = (counts[2] + counts[3]) as f32;
            let total = cis + trans;
            let status = Status::get_status(refvar, variants, molecules, &phasing, &kmers);
            
            let mol_cis = (status.r_h1 + status.a_h2) as f32;
            let mol_trans = (status.r_h2 + status.a_h1) as f32;
            let mol_total = mol_cis + mol_trans;
            let mut add = false;
            let mut phase = true;
            if status.r_both + status.a_both == 0 {
                if (bfs_iter < 4 && total >= 2.0) || total >= 3.0 {
                    if mol_cis.max(mol_trans) / mol_total > 0.95 {
                        if cis > trans {
                            let minor = status.r_h1.min(status.a_h2) as f32;
                            if minor/mol_total > 0.15 {
                                add = true;
                            }
                        } else {
                            let minor = status.r_h2.min(status.a_h1) as f32;
                            if minor/mol_total > 0.15 {
                                add = true;
                                phase = false;
                                altvar = refvar;
                                refvar = pair(altvar);
                            }
                        }
                    }   
                }
            }
            if add {
                if !is_real_block {
                    phase_block_number += 1;
                    is_real_block = true;
                }
                if direction {
                    if phase {
                        kmer_order.push_back((refvar, altvar));
                    } else {
                         kmer_order.push_back((altvar, refvar));
                    }
                } else {
                    if phase {
                        kmer_order.push_front((refvar, altvar));
                    } else {
                        kmer_order.push_front((altvar, refvar));
                    }
                }

                match crib.as_ref() {
                    Some(cheat) => {
                        let mut any = false;
                        if let Some((contig, num, position)) = cheat.variants.get(&refvar.abs()) {
                            let poses = crib_positions.entry(*contig).or_insert(Vec::new());
                            poses.push(*position);
                            
                            eprintln!("ADD ref{} chr{}\t{}\t{:?}\t{:?}\t{}\t{}", num, contig, position, counts,
                                Status::get_status(refvar, variants, molecules, &phasing, &kmers),
                                kmers.kmer_counts.get(&refvar.abs()).unwrap(), kmers.kmer_counts.get(&altvar.abs()).unwrap()); any = true;
                                
                        } 
                        if let Some((contig, num, position)) = cheat.variants.get(&altvar.abs()) {
                            let poses = crib_positions.entry(*contig).or_insert(Vec::new());
                            poses.push(*position);
                            
                            eprintln!("ADD alt{} chr{}\t{}\t{:?}\t{:?}\t{}\t{}", num, contig, position, counts, 
                                Status::get_status(refvar, variants, molecules, &phasing, &kmers),
                                kmers.kmer_counts.get(&refvar.abs()).unwrap(), kmers.kmer_counts.get(&altvar.abs()).unwrap()); any = true;
                            
                        }   
                        if !any {
                            eprintln!("ADD alt/ref not in crib \t{:?}\t{:?}\t{}\t{}",  counts, 
                                Status::get_status(refvar, variants, molecules, &phasing, &kmers),
                                kmers.kmer_counts.get(&refvar.abs()).unwrap(), kmers.kmer_counts.get(&altvar.abs()).unwrap());
                        }
                    }
                    None => {}
                }
            } else {
                match crib.as_ref() {
                    Some(cheat) => {
                        let mut any = false;
                        if let Some((contig, num, position)) = cheat.variants.get(&refvar.abs()) {
                            eprintln!("    FAIL ref{} chr{}\t{}\t{:?}\t{:?}\t{}\t{}", num, contig, position, counts,
                                Status::get_status(refvar, variants, molecules, &phasing, &kmers),
                                kmers.kmer_counts.get(&refvar.abs()).unwrap(), kmers.kmer_counts.get(&altvar.abs()).unwrap()); any = true;
                        } 
                        if let Some((contig, num, position)) = cheat.variants.get(&altvar.abs()) {
                            eprintln!("    FAIL alt{} chr{}\t{}\t{:?}\t{:?}\t{}\t{}", num, contig, position, counts,
                                Status::get_status(refvar, variants, molecules, &phasing, &kmers),
                                kmers.kmer_counts.get(&refvar.abs()).unwrap(), kmers.kmer_counts.get(&altvar.abs()).unwrap()); any = true;
                        } 
                        if !any {
                            eprintln!("    FAIL alt/ref not in crib {:?}\t{:?}\t{}\t{}", counts, 
                                Status::get_status(refvar, variants, molecules, &phasing, &kmers),
                                kmers.kmer_counts.get(&refvar.abs()).unwrap(), kmers.kmer_counts.get(&altvar.abs()).unwrap());
                        } 
                    }
                    None => { 
                    }
                }
            }
            if add {
                phasing.insert(refvar.abs(), true);
                phasing.insert(altvar.abs(), false);
                seeds.visited.insert(altvar); seeds.visited.insert(-altvar); seeds.visited.insert(refvar); seeds.visited.insert(-refvar);
                for var1 in [refvar, altvar].iter() {
                    //eprintln!("here2");
                    let var1 = *var1;
                    if let Some(varset) = adjacency_list.get(&var1) {
                        //eprintln!("varset2 length {}", varset.len());
                        let mut var1mod2 = var1.abs() % 2;
                        if !phase { var1mod2 = (var1.abs() + 1) % 2 }
                        let mut var_vec: Vec<(i32, u16, u16)> = Vec::new();
                        for (var2, (mols, dist)) in varset {
                            var_vec.push((*var2, *mols, *dist));
                        }
                        var_vec.sort_by(|a, b| a.2.cmp(&b.2));
                        for (var2, _mols, _dist) in var_vec {
                            //eprintln!("var22 {} already visited? {}",var2, seeds.visited.contains(&var2));
                            if seeds.visited.contains(&var2) { continue; }
                            let var2mod2 = var2.abs() % 2;
                            let mut varhandle = var2.min(pair(var2));
                            if var2 > 0 {
                                varhandle = var2.max(pair(var2));
                            }
                            if !forward_edges.contains_key(&varhandle) {
                                bfs_queue.push_back(varhandle);
                            }
                            let counts = forward_edges.entry(varhandle).or_insert([0;4]);
                            if var1mod2 == 0 && var2mod2 == 0 {
                                counts[0] += 1;
                            } else if var1mod2 == 1 && var2mod2 == 1 {
                                counts[1] += 1;
                            } else if var1mod2 == 0 && var2mod2 == 1 {
                                counts[2] += 1;
                            } else  if var1mod2 == 1 && var2mod2 == 0 {
                                counts[3] += 1;
                            }
                            else { eprintln!("wtf");}
                        } 
                    }
                }
                
            }
            else {
                let hom1 = (counts[0] + counts[3]) as f32;
                let hom2 = (counts[1] + counts[2]) as f32;
                let total = hom1+hom2;
                if total > 4.0 && is_real_block {
                    seeds.visited.insert(altvar); seeds.visited.insert(-altvar); seeds.visited.insert(refvar); seeds.visited.insert(-refvar);
                }
                if hom1.max(hom2) / total > 0.95 {
                    if hom1 > hom2 {
                        let minor = status.r_h1.min(status.r_h2) as f32;
                        let major = status.r_h1.max(status.r_h2) as f32;
                        //if minor/total > 0.1 {
                        //    eprintln!("homozygous 1\t{:?}\n{:?}", counts, 
                        //        KmerPairStatus::get_status(refvar, variants, molecules, &phasing, &kmers));
                        //} else if major/total > 0.9 {
                        //    eprintln!("heterozygous indel 1\t{:?}\n{:?}", counts, 
                        //        KmerPairStatus::get_status(refvar, variants, molecules, &phasing, &kmers));
                        //} else {
                        //    eprintln!("no idea {:?}\n{:?}",counts, 
                        //        KmerPairStatus::get_status(refvar, variants, molecules, &phasing, &kmers));
                        //}
                    } else {
                        let minor = status.a_h1.min(status.a_h2) as f32;
                        let major = status.a_h1.max(status.a_h2) as f32;
                        //if minor/total > 0.1 {
                        //    eprintln!("homozygous 2\t{:?}\n{:?}", counts, 
                        //        KmerPairStatus::get_status(refvar, variants, molecules, &phasing, &kmers));
                        //} else if major/total > 0.9 {
                        //    eprintln!("heterozygous indel 2\t{:?}\n{:?}", counts, 
                        //        KmerPairStatus::get_status(refvar, variants, molecules, &phasing, &kmers));
                        //} else {
                        //    eprintln!("no idea {:?}\n{:?}",counts, 
                        //        KmerPairStatus::get_status(refvar, variants, molecules, &phasing, &kmers));
                        //}
                    }
                }
            }
            bfs_iter += 1;
        } 
        //if phasing.len()/2 > 4 {
        //    phase_block += 1;
        //}
        // ok lets see what we can do on linearizing
        //linearize(variants, molecules, phasing, kmers);
        
        //eprintln!("phase block contains {} variants",phasing.len()/2); 

    }


    // NOW WE are going to scaffold using linked reads and hic


    if params.txg_mols.len() > 0 {
        let mut scaffolding_phasing: HashMap<(i32, i32), [u32; 3]> = HashMap::new();
        //first we need to build Hashmap from kmer to phasing ends 
        let mut kmer_end_phasings: HashMap<i32, i32> = HashMap::new(); // kmer id to phase block id (phase block will be signed to indicate beginning or end)
        for (phase_block_id, kmer_ordering) in phase_blocks.iter() {
            for (index, (hap1mer, hap2mer)) in kmer_ordering.iter().enumerate() {
                if index < 200 {
                    kmer_end_phasings.insert(*hap1mer, *phase_block_id as i32);
                    kmer_end_phasings.insert(-hap1mer, *phase_block_id as i32);
                    kmer_end_phasings.insert(*hap2mer, *phase_block_id as i32);
                    kmer_end_phasings.insert(-hap2mer, *phase_block_id as i32);
                }
                if kmer_ordering.len() - index < 201 {
                    kmer_end_phasings.insert(*hap1mer, -(*phase_block_id as i32));
                    kmer_end_phasings.insert(-hap1mer, -(*phase_block_id as i32));
                    kmer_end_phasings.insert(*hap2mer, -(*phase_block_id as i32));
                    kmer_end_phasings.insert(-hap2mer, -(*phase_block_id as i32));
                }
            }
        }

        for mol in molecules.get_linked_read_molecules() {
            let mut counts: HashMap<(i32, i32), [u32; 3]> = HashMap::new();
            let mut vars: Vec<i32> = Vec::new();
            for var in molecules.get_linked_read_variants(*mol) {
                vars.push(*var);
            }
            if vars.len() < 2 { continue; }
            for vardex1 in 0..vars.len() {
                let var1 = vars[vardex1];
                if let Some(phase_block_end1) = kmer_end_phasings.get(&var1) {
                    for vardex2 in (vardex1+1)..vars.len() {
                        let var2 = vars[vardex2];
                        if let Some(phase_block_end2) = kmer_end_phasings.get(&var2) {
                            if phase_block_end1.abs() != phase_block_end2.abs() {
                                let min = phase_block_end1.min(phase_block_end2);
                                let max = phase_block_end1.max(phase_block_end2);
                                let count = counts.entry((*min, *max)).or_insert([0;3]);
                                
                                let phase1opt = phasing.get(&var1);
                                let phase2opt = phasing.get(&var2);
                                match phase1opt {
                                    Some(phase1) => match phase2opt {
                                        Some(phase2) => {
                                            if *phase1 == *phase2 && *phase1 {
                                                count[0] += 1;
                                            } else if *phase1 == *phase2 && !*phase1 {
                                                count[1] += 1;
                                            } else {
                                                count[2] += 1;
                                            }
                                        },
                                        None => {}
                                    },
                                    None => {}
                                }  
                            }
                        }
                    }
                }
            }
            for ((phaseend1, phaseend2), count) in counts.iter() {
                let consistent = count[0].max(count[1]) as f32;
                let total = count[0] as f32 + count[1] as f32 + count[2] as f32;
                if consistent > 5.0 && consistent / total > 0.9 {
                    eprintln!("phaseblocks {} and {} with counts {:?}", phaseend1, phaseend2, counts);
                    let consistency = scaffolding_phasing.entry((*phaseend1, *phaseend2)).or_insert([0;3]);
                    if count[0] > count[1] {
                        consistency[0] += 1;
                    } else {
                        consistency[1] += 1;
                    }
                } else if consistent > 5.0 {
                    let consistency = scaffolding_phasing.entry((*phaseend1, *phaseend2)).or_insert([0;3]);
                    consistency[2] += 1;
                    eprintln!("bad phaseblocks {} and {} with counts {:?}", phaseend1, phaseend2, counts);
                }

            }
            for ((pbe1, pbe2), counts) in scaffolding_phasing.iter() {
                let p1 = counts[0] as f32;
                let p2 = counts[1] as f32;
                let total = p1 + p2 + counts[2] as f32;
                if p1/total > 0.9 {
                    eprintln!("scaffolding link from {} -- {} with {:?}", pbe1, pbe2, counts);
                } else if p2/total > 0.9 {
                    eprintln!("scaffolding link from {} -- {} with {:?}", pbe1, pbe2, counts);
                }   else {
                    eprintln!("failure to scaffold from {} -- {} with {:?}", pbe1, pbe2, counts);
                }
            }
        }
    }



    // HIC scaffolding









    let mut single_jump: HashMap<(i32, i32), [u8; 4]> = HashMap::new();
    for mol in molecules.get_long_read_molecules() {
        let mut varlist: Vec<i32> = Vec::new();
        for var in molecules.get_long_read_variants_ordered(*mol) {
            varlist.push(*var);
        }
        if varlist.len() < 2 {continue;}
        for vardex in 0..(varlist.len()-1) {
            let var1 = varlist[vardex];
            for vardex2 in (vardex+1)..(varlist.len()) {
                //let vardex2 = vardex + 1;
                let var2 = varlist[vardex2];
                let min = var1.abs().min(var2.abs());
                let max = var1.abs().max(var2.abs());
                if vardex2 == vardex + 1 {
                    let counts = single_jump.entry((min, max)).or_insert([0;4]);
                    if var1.abs() < var2.abs() {
                        if (var1 > 0 && var2 > 0) || (var1 < 0 && var2 < 0) { counts[0] += 1; }
                        else if (var1 > 0 && var2 < 0) || (var1 < 0 && var2 > 0) { counts[1] += 1; }
                    } else {
                        if (var1 > 0 && var2 > 0) || (var1 < 0 && var2 < 0) { counts[2] += 1; }
                        else if (var1 > 0 && var2 < 0) || (var1 < 0 && var2 > 0) { counts[3] += 1; }
                    }
                }
            }
        }
    }
    let mut single_jump_adjacency_list: HashMap<i32, HashMap<i32, usize>> = HashMap::new();

    for ((var1, var2), counts) in single_jump.iter() {
        if seeds.visited.contains(var1) && seeds.visited.contains(var2) {
            let mut max = 0;
            let mut maxdex = 0;
            for (index, count) in counts.iter().enumerate() {
                if *count > max {
                    max = *count;
                    maxdex = index;
                }
            }
            if maxdex == 0 {
                let adjacency_map = single_jump_adjacency_list.entry(*var1).or_insert(HashMap::new());
                let count = adjacency_map.entry(*var2).or_insert(0);
                *count += 1;
                let adjacency_map = single_jump_adjacency_list.entry(-var2).or_insert(HashMap::new());
                let count = adjacency_map.entry(-var1).or_insert(0);
                *count += 1;
            } else if maxdex == 1 {
                let adjacency_map = single_jump_adjacency_list.entry(*var1).or_insert(HashMap::new());
                let count = adjacency_map.entry(-var2).or_insert(0);
                *count += 1;
                let adjacency_map = single_jump_adjacency_list.entry(*var2).or_insert(HashMap::new());
                let count = adjacency_map.entry(-var1).or_insert(0);
                *count += 1;
            } else if maxdex == 2 {
                let adjacency_map = single_jump_adjacency_list.entry(*var2).or_insert(HashMap::new());
                let count = adjacency_map.entry(*var1).or_insert(0);
                *count += 1;
                let adjacency_map = single_jump_adjacency_list.entry(-var1).or_insert(HashMap::new());
                let count = adjacency_map.entry(-var2).or_insert(0);
                *count += 1;
            } else {
                let adjacency_map = single_jump_adjacency_list.entry(-var2).or_insert(HashMap::new());
                let count = adjacency_map.entry(*var1).or_insert(0);
                *count += 1;
                let adjacency_map = single_jump_adjacency_list.entry(-var1).or_insert(HashMap::new());
                let count = adjacency_map.entry(*var2).or_insert(0);
                *count += 1;
            }
        }
    }

}

fn linearize(variants: &Variants, molecules: &Molecules, phasing: HashMap<i32, bool>, kmers: &Kmers) {

    let mut order_counts: HashMap<i32, HashMap<i32, usize>> = HashMap::new();

    for (var, _phase) in &phasing {
        let mut var = *var;
        for mol in variants.get_long_read_molecules(var) {
            let mut start = false;
            //if *mol < 0 { var = -var; }
            for var2 in molecules.get_long_read_variants_ordered(mol.abs()) {
                if var2.abs() == var {
                    start = true; 
                    if *var2 < 0 { var = -var; }
                    continue;
                }
                if start {
                    if phasing.contains_key(&var2.abs()) {
                        //eprintln!("from mol {} adding counts for {} -> {} and {} -> {}",mol, var, var2, -var2, -var);
                        {
                            let counts = order_counts.entry(var).or_insert(HashMap::new());
                            let count = counts.entry(*var2).or_insert(0);
                            *count += 1;
                        }
                        {
                            let counts = order_counts.entry(-var2).or_insert(HashMap::new());
                            let count = counts.entry(-var).or_insert(0);
                            *count += 1;
                        }
                        break;
                    }
                }
            }
        }

    }
    // calculate winners
    let mut order_map: HashMap<i32, i32> = HashMap::new();
    for (var, varcounts) in order_counts.iter() {
        let mut count_vec: Vec<(&i32, &usize)> = varcounts.iter().collect();
        count_vec.sort_by(|a, b| b.1.cmp(a.1));
        /*if count_vec.len() > 1 {
            eprintln!("var {} best {} count {} second best {} second best count {}", var, count_vec[0].0, count_vec[0].1, count_vec[1].0, count_vec[1].1);
        } else {
            eprintln!("var {} best {} count {}", var, count_vec[0].0, count_vec[0].1);
        }*/
        if varcounts.len() == 1 || (*count_vec[0].1 as f32)/((count_vec[0].1 + count_vec[1].1) as f32) > 0.9 {
            //eprintln!("\t{} -- {};", kmers.kmers.get(&var.abs()).unwrap(), kmers.kmers.get(&count_vec[0].0.abs()).unwrap());
            order_map.insert(*var, *count_vec[0].0);
        } //else {
        //    eprintln!("didn't have good support for best next {}",var);
        //}
    }
    
    
    let mut used: HashSet<i32> = HashSet::new();
    
    let mut orderings: Vec<Vec<i32>> = Vec::new();
    for (var, phase) in &phasing {
        if used.contains(var) { continue; }
        let mut cur_var = *var;
        used.insert(cur_var); 
        // starting at cur_var and walk till end
        loop {
            //eprintln!("cur var {}", cur_var);
            if order_map.contains_key(&cur_var) {
                let next_var = *order_map.get(&cur_var).unwrap();
                if !used.contains(&next_var) {
                    cur_var = next_var;
                    used.insert(cur_var);
                } else { break; }
            } else {
                //eprintln!("break");
                break;
            }
        }
        // start at end (cur_var) and walk till other end, this is a path
        let mut order: Vec<i32> = Vec::new();
        cur_var = -cur_var;
        used.insert(cur_var);
        order.push(cur_var);
        loop {
            //eprintln!("cur var2 {}",cur_var);
            if order_map.contains_key(&cur_var) {
                let next_var = *order_map.get(&cur_var).unwrap();
                if !used.contains(&next_var) {
                    cur_var = next_var;
                    used.insert(cur_var);
                    order.push(cur_var);
                } else {
                     orderings.push(order);
                    break;
                }
            } else {
                //eprintln!("break2");
                orderings.push(order);
                break;
            }
        }
       
    }
    /*eprintln!("for this phase block we have {} orderings",orderings.len());
    for (index, ordering) in orderings.iter().enumerate() {
        eprintln!("\tordering {}",index);
        for var in ordering {
            eprintln!("\t\t{}\t{}",var, phasing.get(&var.abs()).unwrap());
        }
    }*/

}

fn pair(v: i32) -> i32 {
    if v < 0 { 
        if v % 2 == 0 { v + 1 } 
        else { v - 1 } 
    }
    else { 
        if v % 2 == 0 { v - 1 } 
        else { v + 1 } 
    }
}

fn sparsembly(molecules: &Molecules, variants: &Variants, kmers: &Kmers, crib: &Option<Crib>) -> (HashMap<i32,HashSet<i32>>, HashMap<i32, HashMap<i32, (u16, u16)>>) {
    let mut graph: HashMap<(i32, i32), [u16; 8]> = HashMap::new();
    //eprintln!("sparsembly");
    for mol in molecules.get_long_read_molecules() {
        let mut varlist: Vec<(i32, i32)> = Vec::new();
        for (var, pos) in molecules.get_long_read_variants_and_positions(*mol) {
            varlist.push((*var, *pos));
        }
        if varlist.len() < 2 {continue;}
        for vardex in 0..(varlist.len()-1) {
            let var1 = varlist[vardex].0;
            let pos1 = varlist[vardex].1;
            for vardex2 in (vardex+1)..(varlist.len()) {
                let var2 = varlist[vardex2].0;
                let pos2 = varlist[vardex2].1;
                let min = var1.abs().min(var2.abs());
                let max = var1.abs().max(var2.abs());
                let counts = graph.entry((min, max)).or_insert([0;8]);
                // counts encodes the following
                // counts[0] = counts of for1 -> for2 + rev2 -> rev1
                // counts[1] = counts of for1 -> rev2 + for2 -> rev1
                // counts[2] = counts of for2 -> for1 + ref1 -> rev2
                // counts[3] = counts of rev2 -> for1 + for2 -> rev1
                // then counts[4-7] are the sum of distances between the kmers (so counts[4]/counts[0] is average distance)
                if var1.abs() < var2.abs() {
                    if var1 > 0 && var2 > 0 { counts[0] += 1; counts[4] += (pos2-pos1) as u16; }
                    else if var1 > 0 && var2 < 0 { counts[1] += 1; counts[5] += (pos2-pos1) as u16; }
                    else if var1 < 0 && var2 < 0 { counts[2] += 1; counts[6] += (pos2-pos1) as u16;}
                    else { counts[3] += 1; counts[7] += (pos2-pos1) as u16; }
                } else {
                    if var2 < 0 && var1 < 0 { counts[0] += 1; counts[4] += (pos2-pos1) as u16; }
                    else if var2 < 0 && var1 > 0 { counts[1] += 1; counts[5] += (pos2-pos1) as u16; }
                    else if var2 > 0 && var1 > 0 { counts[2] += 1; counts[6] += (pos2-pos1) as u16; }
                    else { counts[3] += 1; counts[7] += (pos2-pos1) as u16; }
                }
            }
        }
    }
    //eprintln!("graph size {}", graph.len());
    let mut adjacency_list: HashMap<i32, HashMap<i32, (u16, u16)>> = HashMap::new();
    

    let mut het_vars: Vec<i32> = Vec::new();
    let mut het_var_index: HashMap<i32, usize> = HashMap::new();
    //println!("graph longreads {{");
    for var in variants.get_variant_iter(KmerType::PairedHet) {
        //println!("\t{} -- {};", var, -var);
        het_var_index.insert(*var, het_vars.len());
        het_vars.push(*var);
        het_var_index.insert(-var, het_vars.len());
        het_vars.push(-var);
    }
    let mut union_find: Vec<usize> = (0..het_vars.len()).collect();
    //eprintln!("adjacency list has {} variants",adjacency_list.len());

    for ((var1, var2), counts) in graph.iter() {
        let mut max = 0;
        let mut maxdex = 0;
        let mut total = 0;
        for (index, count) in counts.iter().enumerate() {
            if index > 3 { break; }
            if *count > max {
                max = *count;
                maxdex = index;
            }
            total += *count;
        }
        if total < 3 { continue; }
        let var1dex = *het_var_index.get(var1).unwrap();
        let var1compdex = *het_var_index.get(&-var1).unwrap();
        let var2dex = *het_var_index.get(var2).unwrap();
        let var2compdex = *het_var_index.get(&-var2).unwrap();
        //eprintln!("duh {} {} {:?}", var1, var2, counts);
        let percent = (max as f32)/(total as f32);
        let distance = counts[maxdex + 4]/max;
        /*
        if *var1 == 2304 || *var2 == 2304 || *var1 == 1726 || *var2 ==1726 {
            eprintln!("var1 {}\tvar2 {}\tcounts {:?}", var1, var2, counts);
            match crib.as_ref() {
                Some(cheat) => {
                    //eprintln!("here, we have crib");
                    if let Some((contig, _order, position)) = cheat.variants.get(&pair(*var1)) {
                        eprintln!("cribsheet var1 {}\t{}", contig, position);
                    } 
                    if let Some((contig, _order, position)) = cheat.variants.get(&pair(*var2)) {
                        eprintln!("cribsheet var2 {}\t{}", contig, position);
                    }
                }
                None => { eprintln!("NoContig\t{}\t{}", 
                    kmers.kmers.get(var1).unwrap(), kmers.kmers.get(var2).unwrap());
                }
            }
        }
        */
        if percent > 0.9 {
            //merge(&mut union_find, var1dex, var2dex);
            if maxdex == 0 {
                merge(&mut union_find, var1dex, var2dex);
                merge(&mut union_find, var1compdex, var2compdex);
                //println!("n{} -> n{};", var1, var2);
                //println!("nC{} -> nC{}", var2, var1);
                let adjacency_set = adjacency_list.entry(*var1).or_insert(HashMap::new());
                adjacency_set.insert(*var2, (max, distance));
                let adjacency_set = adjacency_list.entry(-*var2).or_insert(HashMap::new());
                adjacency_set.insert(-*var1, (max, distance));

            } else if maxdex == 1 {
                merge(&mut union_find, var1dex, var2compdex);
                merge(&mut union_find, var1compdex, var2dex);
                //println!("n{} -> nC{};", var1, var2);
                //println!("n{} -> nC{}", var2, var1);
                let adjacency_set = adjacency_list.entry(*var1).or_insert(HashMap::new());
                adjacency_set.insert(-*var2, (max, distance));
                let adjacency_set = adjacency_list.entry(*var2).or_insert(HashMap::new());
                adjacency_set.insert(-*var1, (max, distance));
            } else if maxdex == 2 {
                merge(&mut union_find, var1dex, var2dex);
                merge(&mut union_find, var1compdex, var2compdex);
                //println!("n{} -> n{};", var2, var1);
                //println!("nC{} -> nC{};", var1, var2);
                let adjacency_set = adjacency_list.entry(*var2).or_insert(HashMap::new());
                adjacency_set.insert(*var1, (max, distance));
                let adjacency_set = adjacency_list.entry(-*var1).or_insert(HashMap::new());
                adjacency_set.insert(-*var2, (max, distance));
            } else if maxdex == 3 {
                merge(&mut union_find, var1dex, var2compdex);
                merge(&mut union_find, var1compdex, var2dex);
                //println!("nC{} -> n{};", var2, var1);
                //println!("nC{} -> n{}", var1, var2);
                let adjacency_set = adjacency_list.entry(-*var2).or_insert(HashMap::new());
                adjacency_set.insert(*var1, (max, distance));
                let adjacency_set = adjacency_list.entry(-*var1).or_insert(HashMap::new());
                adjacency_set.insert(*var2, (max, distance));
            }
        } 
    }
    //detect_loops(&mut adjacency_set);
    flatten_union_find(&mut union_find);
    let mut disjoint_sets: HashMap<usize, HashSet<usize>> = HashMap::new();
    for (index, parent) in union_find.iter().enumerate() {
        if index == *parent {
            let mut set: HashSet<usize> = HashSet::new();
            set.insert(index);
            disjoint_sets.insert(index, set);
        } else {
            match disjoint_sets.get_mut(parent) {
                Some(set) => { set.insert(index); },
                None => { assert!(false, "union find messed up"); },
            }
        }
    }
    let mut total = 0;
    let mut sorted_comps: Vec<(usize, HashSet<usize>)> = Vec::new();
    for (root, set) in disjoint_sets.iter() {
        sorted_comps.push((*root, set.clone()));
        total += set.len();
    }
    //eprintln!("total for N50 calculation {}",total);
    sorted_comps.sort_by(|a,b| b.1.len().cmp(&a.1.len()));
    let mut so_far = 0;
    let mut done = false;
    let mut N50 = 0;
    for (index, (root, set)) in sorted_comps.iter().enumerate() {
        so_far += set.len();
        if !done && (so_far as f32)/(total as f32) > 0.5 {
            done = true;
            N50 = set.len();
            //break;
        }
        //eprintln!("comp\t{}\tsize\t{}",root, set.len());
        if set.len() == 1 {
            let mut badvar = 0;
            for tmpvar in set { badvar = *tmpvar; }
            let badvar = het_vars[badvar];
            /*
            if let Some(kmer_string) = kmers.kmers.get(&badvar) {
                if let Some(kmer_count) = kmers.kmer_counts.get(&badvar) {
                    eprintln!("inspecting bad var {} = {} with {}",badvar, kmer_string, kmer_count);
                } else {
                    eprintln!("inspecting bad var {} = {} with {}",badvar, kmer_string, "no counts?");
                }
            } else {
                eprintln!("inspecting bad var {} = {} with {}",badvar, "i dont know what this kmer is?", "no counts?");
            }
            
            
            for mol in variants.get_long_read_molecules(badvar) {
                eprintln!("\ton long mol {}",mol);
                for var in molecules.get_long_read_variants_ordered(*mol) {
                    if let Some(counts) = graph.get(&(badvar.abs().min(var.abs()), badvar.abs().max(var.abs()))) {
                        eprintln!("\t\t{} -- {}\t{:?}",badvar, var, counts);
                    }
                }
                
            }
            for mol in variants.get_linked_read_molecules(badvar) {
                eprintln!("\t FR on linked-read mol {}",mol);
            }
            for mol in variants.get_linked_read_molecules(-badvar) {
                eprintln!("\t RC on linked-read mol {}",mol);
            }
            */
        }
    }
    //eprintln!("N50\t{}",N50);
    let biggest_block = sorted_comps[0].1.clone();

    let mut bfs: VecDeque<i32> = VecDeque::new();
    bfs.push_back(52);
    let mut vars_of_interest: HashSet<i32> = HashSet::new();
    vars_of_interest.insert(52);
    while vars_of_interest.len() < 4000 {
        if let Some(var) = bfs.pop_front() {
            if let Some(varlist) = adjacency_list.get(&var) {
                for (var2, mols) in varlist {
                    //eprintln!("node {}",var2);
                    if vars_of_interest.contains(var2) {
                        continue;
                    }
                    bfs.push_back(*var2);
                    vars_of_interest.insert(*var2);
                }   
            } else { continue; }
        } else {break;}
    }

    //let mut block_vars: HashSet<i32> = HashSet::new();
    //for index in biggest_block {
    //    block_vars.insert(het_vars[index]);
    //}
    let mut children: HashSet<i32> = HashSet::new();
    for node in bfs {
        let mut prefix = "";
        if node < 0 { prefix = "c"; }
        children.insert(node);
        //println!("\tn{}{} [color=red];",prefix, node.abs());
    }
    
    for var1 in adjacency_list.keys() {
        for (var2, mols) in adjacency_list.get(var1).unwrap() {
            if !vars_of_interest.contains(var1) || !vars_of_interest.contains(var2) {continue;}
            let mut prefix1 = "";
            let mut prefix2 = "";
            if *var1 < 0 {
                prefix1 = "c";
            }
            if *var2 < 0 { prefix2 = "c"; }
            let mut suffix = "";
            if children.contains(var1) && children.contains(var2) { suffix = "[color=red]"}
            //println!("\tn{}{} -> n{}{} {};", prefix1, var1.abs(), prefix2, var2.abs(), suffix);
        }
    }

    //println!("}}");
    let mut components: HashMap<i32, HashSet<i32>> = HashMap::new();
    for (vardex, vardexset) in disjoint_sets.iter() {
        let var = het_vars[*vardex];
        let mut varset: HashSet<i32> = HashSet::new();
        
        for sinkdex in vardexset.iter() {
            let var = het_vars[*sinkdex];
            varset.insert(var);
        }
        components.insert(var, varset);
    }
    (components, adjacency_list)
}

fn detect_loops(adjacency_list: &HashMap<i32, HashSet<i32>>, varlist: &Vec<i32>) -> HashMap<i32, usize> {
    let mut visited: HashSet<i32> = HashSet::new();
    let mut dfs_stack: VecDeque<i32> = VecDeque::new();
    let mut active_stack: VecDeque<usize> = VecDeque::new();
    let mut min_link: Vec<usize> = Vec::new();
    
    for (index, var) in varlist.iter().enumerate() {
        min_link.push(index);
    }
    let mut traverse_order_to_var: Vec<i32> = Vec::new();
    let mut var_to_index: HashMap<i32, usize> = HashMap::new();

    for (var, varset) in adjacency_list {
        if visited.contains(&var) { continue; }
        //println!("top level, var {}",var);
        dfs_stack.push_back(*var);
        let mut active: HashSet<usize> = HashSet::new();
        //println!("start dfs");
        while dfs_stack.len() != 0 {
            let var = dfs_stack.pop_back().unwrap();
            if visited.contains(&var) { 
                //println!("already visited var {} vardex {}, continuing",var, var_to_index.get(&var).unwrap()); 
                continue; 
            }
            let vardex;
            if !var_to_index.contains_key(&var) {
                vardex = traverse_order_to_var.len();
                traverse_order_to_var.push(var);
                //println!("var {} is traverse order {}",var, vardex);
                var_to_index.insert(var, vardex);
            } else {
                vardex = *var_to_index.get(&var).unwrap();
            }
            //let var_minlink = min_link[*vardex];
            //println!("popped var {} vardex {}",var, vardex);
            if active.contains(&vardex) {
                let mut minlink = min_link[vardex];
                //println!("found active var {} vardex {} with minlink {}", var, vardex, minlink);
                let mut vardex = active_stack.pop_back().unwrap();
                visited.insert(traverse_order_to_var[vardex]);
                active.remove(&vardex);
                //println!("popped var {} vardex {} from active stack, old minlink {} new minlink {}", traverse_order_to_var[vardex], vardex, min_link[vardex], min_link[vardex].min(minlink));
                min_link[vardex] = min_link[vardex].min(minlink);
                while min_link[vardex] != vardex {
                    minlink = min_link[vardex];
                    vardex = active_stack.pop_back().unwrap();
                    active.remove(&vardex);
                    visited.insert(traverse_order_to_var[vardex]);
                    //println!("in active rollback, popped var {} vardex {}, old minlink {} new minlink {}, stack length {}", traverse_order_to_var[vardex], vardex, minlink, min_link[vardex].min(minlink), active_stack.len());
                    min_link[vardex] = min_link[vardex].min(minlink);
                }
                continue;
            }
            active.insert(vardex);
            active_stack.push_back(vardex);
            if let Some(varset) = adjacency_list.get(&var) {
                for var2 in varset {
                    dfs_stack.push_back(*var2);
                }
            }
        }
    }

    let mut assignments: HashMap<i32, usize> = HashMap::new();
    for (index, link) in min_link.iter().enumerate() {
        assignments.insert(traverse_order_to_var[index], *link);
    }
    assignments
}

fn flatten_union_find(union_find: &mut Vec<usize>) {
    for index in (0..(union_find.len())).rev() {
        let mut root = index;
        while union_find[root] != root {
            root = union_find[root];
        }
        union_find[index] = root;
    }
}

fn merge(union_find: &mut Vec<usize>, var1dex: usize, var2dex: usize) {
    let mut root1 = var1dex;
    while union_find[root1] != root1 {
        root1 = union_find[root1];
    }
    let mut root2 = var2dex;
    while union_find[root2] != root2 {
        root2 = union_find[root2];
    }
    if root1 < root2 { union_find[root2] = root1; }
    else { union_find[root1] = root2; }
}

struct Phasing {
    variant_phasing: HashMap<i32, i32>,
    phase_block_sizes: HashMap<i32, usize>,
}

fn make_phase_blocks(molecules: &Molecules, variants: &Variants) -> Phasing {
    let mut order: Vec<i32> = Vec::new();
    for var in variants.get_variant_iter(KmerType::PairedHet) { order.push(*var); }
    let seed: [u8; 32] = [6; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    order.shuffle(&mut rng);
    let mut max = 0;
    for x in &order {
        if *x > max { max = *x; }
    }
    let mut used_variants = BitVec::from_elem((max+1) as usize, false);
    let mut removed_variants = BitVec::from_elem((max+1) as usize, false);
    let mut phasing: HashMap<i32, i32> = HashMap::new();
    let mut phase_block_sizes: HashMap<i32, usize> = HashMap::new();

    let mut phase_block = 1;
    for start_var in order {
        if used_variants.get(start_var as usize).unwrap() { continue; }
        println!("starting phase block with variant {}", start_var);
        let mut window: Window = Window::new();
        window.add_variant(start_var as i32, variants, molecules, &mut used_variants, &mut removed_variants);
        loop {
            if let Some(next_var) = window.get_next_variant(variants, molecules, 
                &mut used_variants, &phasing, &phase_block_sizes, &mut removed_variants) {
                //println!("adding variant {}", next_var);
                window.add_variant(next_var, variants, molecules, &mut used_variants, &mut removed_variants);
            } else { 
                if window.variants.len() == 1 { 
                    used_variants.set(start_var as usize, true);
                    removed_variants.set(start_var as usize, true);
                    continue; }
                println!("ending phase block {} with {} variants", phase_block, window.variants.len());
                for var in window.variants.iter() {
                    let mut block = phase_block;
                    if *var < 0 { block = -block; }
                    phasing.insert(var.abs(), block);
                }
                phase_block_sizes.insert(phase_block, window.variants.len());
                
                break; 
            }
        }
        println!("block has {} potential merges", window.merge_block_counts.len());
        for (block, counts) in window.merge_block_counts.iter() {
            println!("\tblock {} with counts {:?}", block, counts);
        }
        let mut legit_merges: Vec<i32> = Vec::new();
        for (block, counts) in window.merge_block_counts.iter() {
            let cis = counts[0] as f32;
            let trans = counts[1] as f32;
            let total = cis + trans;
            if cis.max(trans)/total > 0.9 && cis.max(trans) > 4.0 {
                legit_merges.push(*block);
            }
        }
        if legit_merges.len() == 1 {
            let block = legit_merges[0];
            let total_block_size;
            {
            let old_size = phase_block_sizes.get(&block).unwrap();
            let new_size = window.variants.len();
            println!("merging {} with {} with {} and {} variants, for a total of {}", phase_block, block, old_size, new_size, old_size+new_size);
            total_block_size = old_size + new_size;
            }
            let counts = window.merge_block_counts.get(&block).unwrap();
            let cis = counts[0];
            let trans = counts[1];
            let mut block_block_phase: i32= 1;
            if trans > cis { block_block_phase = -1; }
            for variant in window.variants.iter() {
                if *variant > 0 {
                    phasing.insert(variant.abs(), block * block_block_phase);
                } else {
                    phasing.insert(variant.abs(), -block * block_block_phase);
                }
            }
            phase_block_sizes.remove(&phase_block);
            
            phase_block_sizes.insert(block, total_block_size);
        }
        phase_block += 1;
    }
    Phasing { 
        variant_phasing: phasing,
        phase_block_sizes: phase_block_sizes,
    }
}

struct Window {
    current_molecules: HashMap<i32, [u32; 2]>, // set of current molecules and the haplotype counts of their kmers
    variants: HashSet<i32>, // set of variants with negative bit to indicate phase
    variants_tried: HashSet<i32>,
    molecule_last_seen: HashMap<i32, usize>,
    potential_variants: BTreeSet<(usize, i32)>, // sorted by number of molecules touching variant also in current_molecules
    potential_variants_counts: HashMap<i32, usize>, // mapping from var to count to be able to change potential_variants
    potential_num_variants_touching: HashMap<i32, usize>, // required to stipulate new variants must touch multiple variants in current phase block
    iteration: usize,
    merge_block_counts: HashMap<i32, [u32; 2]>,
}

impl Window {
    fn new() -> Window {
        Window{
            current_molecules: HashMap::new(),
            variants: HashSet::new(),
            variants_tried: HashSet::new(),
            molecule_last_seen: HashMap::new(),
            potential_variants: BTreeSet::new(),
            potential_variants_counts: HashMap::new(),
            potential_num_variants_touching: HashMap::new(),
            iteration: 0,
            merge_block_counts: HashMap::new(),
        }
    }

    fn remove_old_molecules(&mut self, molecules: &Molecules, variants: &Variants) {
        let mut to_remove: Vec<i32> = Vec::new();
        for (mol, last_seen) in self.molecule_last_seen.iter() {
            if self.iteration - last_seen > 1000 {
                to_remove.push(*mol);
            }
        }
        for mol in to_remove {
            for var in molecules.get_variants(&mol, KmerType::PairedHet) {
                let var = var.abs();
                if self.potential_variants_counts.contains_key(&var) {
                    let mut old_count = 0;
                    let mut old_percent_x_1000 = 0;
                    if let Some(percent_x_1000) = self.potential_variants_counts.get(&var) {
                        old_percent_x_1000 = *percent_x_1000;
                        old_count = ((*percent_x_1000 as f64) * (variants.get_num_molecules(&var, DataType::PBLR) as f64)/1000000.0).round() as usize;
                    }
                    if old_count <= 1 {
                        self.potential_variants.remove(&(old_percent_x_1000, var));
                        self.potential_variants_counts.remove(&var);
                        self.potential_num_variants_touching.remove(&var);
                    } else {
                        let new_count = old_count - 1;
                        let new_percent_x_1000 = (new_count as f64) / (variants.get_num_molecules(&var, DataType::PBLR) as f64) * 1000000.0;
                        let new_percent_x_1000: usize = new_percent_x_1000.round() as usize;
                        self.potential_variants.remove(&(old_percent_x_1000, var));
                        self.potential_variants.insert((new_percent_x_1000, var));
                        self.potential_variants_counts.insert(var, new_percent_x_1000);
                    }
                }
            }
            self.current_molecules.remove(&mol);
            self.molecule_last_seen.remove(&mol);
        }
    }

    fn add_variant(&mut self, var: i32, variants: &Variants, 
        molecules: &Molecules, used_variants: &mut BitVec, removed_variants: &mut BitVec) {
        self.variants.insert(var);
        self.remove_potential_variant(&var.abs());

        let mut touched_vars: HashSet<i32> = HashSet::new();
        for mol in variants.get_molecules(&var.abs(), DataType::PBLR) {
            // what if -phased_mol also in molecules?
            self.molecule_last_seen.insert(mol.abs(), self.iteration);
            for var2 in molecules.get_variants(&mol.abs(), KmerType::PairedHet) {
                let var2 = var2.abs();
                //let used = used_variants.get(var2 as usize).unwrap();
                //if used { continue; }
                if self.variants_tried.contains(&var2) { continue; }
                let removed = removed_variants.get(var2 as usize).unwrap();
                if removed { continue; }
                if self.variants.contains(&var2) || self.variants.contains(&(-var2)) { continue; } 
                touched_vars.insert(var2);
                if self.current_molecules.contains_key(&mol.abs()) { continue; }
                let old_count: f64;
                {
                    let percent_x_1000 = self.potential_variants_counts.entry(var2).or_insert(0);
                    self.potential_variants.remove(&(*percent_x_1000, var2));
                    old_count = ((*percent_x_1000 as f64) * variants.get_num_molecules(&var2, DataType::PBLR) as f64)/1000000.0;
                }
                let new_percent_x_1000 = (((old_count + 1.0) * 1000000.0) as f64) /
                    (variants.get_num_molecules(&var2, DataType::PBLR) as f64);
                let new_percent_x_1000: usize = new_percent_x_1000.round() as usize;
                self.potential_variants.insert((new_percent_x_1000, var2)); 
                self.potential_variants_counts.insert(var2, new_percent_x_1000);
            }

            let molphase = self.current_molecules.entry(mol.abs()).or_insert([0;2]);
            if (var > 0 && *mol > 0) || (var < 0 && *mol < 0) { molphase[0] += 1; } 
            else { molphase[1] += 1; }
        }
        for var2 in touched_vars {
            if self.potential_variants_counts.contains_key(&var2) {
                let count = self.potential_num_variants_touching.entry(var2).or_insert(0);
                *count += 1;
            }
        }
        used_variants.set(var.abs() as usize, true);
        self.iteration += 1;
        if self.iteration % 500 == 0 {  // TODO make this settable 
            self.remove_old_molecules(molecules, variants);
        }
    }

    fn get_next_variant(&mut self, variants: &Variants, _molecules: &Molecules, used_variants: &mut BitVec, 
        phase_blocks: &HashMap<i32, i32>, _block_sizes: &HashMap<i32, usize>, removed_variants: &mut BitVec) -> Option<i32> {
        let mut to_remove: Option<i32> = None;
        loop {
            if let Some(variant) = to_remove {
                self.remove_potential_variant(&variant);
                to_remove = None;
            }
            if self.potential_variants.len() == 0 { return None; }
            if let Some((percent_x_1000, variant)) = self.potential_variants.pop_last() {
                let count = (percent_x_1000 as f64) * (variants.get_num_molecules(&variant, DataType::PBLR) as f64) / 1000000.0;
                let count: usize = count.round() as usize;

                let mut phased_counts: [u32; 4] = [0;4];
                for molecule in variants.get_molecules(&variant, DataType::PBLR) {
                    if let Some(phase_counts) = self.current_molecules.get(&molecule.abs()) {
                        let hap1 = phase_counts[0] as f32;
                        let hap2 = phase_counts[1] as f32;
                        let total = hap1+hap2;
                        if hap1/total > 0.95 {
                            if *molecule > 0 { phased_counts[0] += 1; } else { phased_counts[2] += 1; }
                        } else if hap2/total > 0.95 {
                            if *molecule > 0 { phased_counts[3] += 1; } else { phased_counts[1] += 1; }
                        }
                    }
                }
                let percent = (percent_x_1000 as f32)/10000.0;
                if count < 10 { 
                    println!("iteration {}, bailing because best variant has count {}, double check {:?}, %{:.2}, mols {}", 
                        self.iteration, count, phased_counts, percent, variants.get_num_molecules(&variant, DataType::PBLR));
                    return None; 
                }
                let cis = (phased_counts[0] + phased_counts[1]) as f32;
                let trans = (phased_counts[2] + phased_counts[3]) as f32;
                let total = cis + trans;
                let variants_touched = self.potential_num_variants_touching.get(&variant).unwrap();
                let minor = if cis > trans { phased_counts[0].min(phased_counts[1]) as f32 } else { phased_counts[2].min(phased_counts[3]) as f32 };
                let major = if cis > trans { phased_counts[0].max(phased_counts[1]) as f32 } else { phased_counts[2].max(phased_counts[3]) as f32 };

                if used_variants.get(variant as usize).unwrap() {
                    if cis.max(trans)/total > 0.95 && minor/(minor+major) > 0.25 { // TODO make this settable
                        println!("wanting to jump to used variant. count {}, percent {:.2}, num mols {}, phasing {:?}, variants touched in current block {}, current block size {}",
                            count, percent, variants.get_num_molecules(&variant, DataType::PBLR), phased_counts, variants_touched, self.variants.len());
                        if let Some(block) = phase_blocks.get(&variant) {
                            let counts = self.merge_block_counts.entry(block.abs()).or_insert([0;2]);
                            if (cis > trans && *block > 0) || (trans > cis && *block < 0) {
                                counts[0] += 1;
                            } else {
                                counts[1] += 1;
                            }
                            //println!("should these merge? other phase block {}", block);

                            //if let Some(size) = block_sizes.get(&block.abs()) {
                             //   println!("with size {}", size);
                            //} else { println!("no block size?"); }
                        } // else { println!("used but no block to merge to?"); }
                    }
                    self.variants_tried.insert(variant);
                    continue;
                }
                if *variants_touched < self.variants.len().min(4) {
                    to_remove = Some(variant);
                    println!("continuing because variant only touched {} other variants {:?}, %{:.1}, #{}, mols {}", 
                        variants_touched, phased_counts, percent, count, variants.get_num_molecules(&variant, DataType::PBLR));
                    self.variants_tried.insert(variant);
                    continue;
                }


                if cis > trans && cis/total > 0.95 { // TODO make this settable
                    let minor = phased_counts[0].min(phased_counts[1]) as f32;
                    if minor/cis > 0.25 { // TODO make this settable
                        println!("adding variant {:?}, %{:.1}, #{}, mols {}", 
                            phased_counts, percent, count, variants.get_num_molecules(&variant, DataType::PBLR));
                        self.variants_tried.insert(variant);
                        return Some(variant);
                    } else {
                        to_remove = Some(variant);
                        //used_variants.set(variant as usize, true);
                        println!("continuing because variant not really phasing consistent {:?}, %{}, #{}, mols {}", 
                            phased_counts, percent_x_1000, count, variants.get_num_molecules(&variant, DataType::PBLR));
                        //if self.variants.len() > 4 { removed_variants.set(variant as usize, true); }
                        self.variants_tried.insert(variant);
                        
                        println!("continuing because variant not really phasing consistent {:?}, {}", phased_counts, percent_x_1000);
                        continue;
                    }
                } else if trans/total > 0.95 { // TODO make this settable
                    let minor = phased_counts[2].min(phased_counts[3]) as f32;
                    if minor/trans > 0.25 { // TODO make this settable
                        println!("adding -variant {:?}, %{:.1}, #{}, mols {}", 
                            phased_counts, percent, count, variants.get_num_molecules(&variant, DataType::PBLR));
                        //println!("adding -variant {:?}, {}", phased_counts, percent_x_1000);
                        self.variants_tried.insert(variant);
                        return Some(-variant);
                    } else {
                        to_remove = Some(variant);
                        //used_variants.set(variant as usize, true);
                        //if self.variants.len() > 4 { removed_variants.set(variant as usize, true); }
                        println!("continuing because variant not really phasing consistent {:?}, %{:.1}, #{}, mols {}", 
                            phased_counts, percent, count, variants.get_num_molecules(&variant, DataType::PBLR));
                        self.variants_tried.insert(variant);
                        continue;
                    }
                } else {
                    to_remove = Some(variant);
                    //used_variants.set(variant as usize, true);
                    println!("continuing because variant not really phasing consistent {:?}, %{:.2}, #{}, mols {}", 
                            phased_counts, percent, count, variants.get_num_molecules(&variant, DataType::PBLR));
                    
                    if self.variants.len() > 4 { removed_variants.set(variant as usize, true); }
                    self.variants_tried.insert(variant);
                    continue;
                    //println!("continuing because variant not really phasing consistent {:?}, {}", phased_counts, percent_x_1000);
                }
            }
        }
        //if let Some(variant) = to_remove {
       //     self.remove_potential_variant(&variant);
        //}
        //None
    }

    fn remove_potential_variant(&mut self, variant: &i32) {
        if let Some(count) = self.potential_variants_counts.get(&variant) {
            self.potential_variants.remove(&(*count, *variant));
        }
        self.potential_num_variants_touching.remove(variant);
        self.potential_variants_counts.remove(variant);
    }
}



#[allow(dead_code)]
fn draw_block(block: i32, consistency_graph: &PhasingConsistencyGraph) -> Result<()> {
    let writer = File::create(format!("block_{}.dot",block))
                .expect("Unable to create file");
    let mut writer = BufWriter::new(writer);
    let variants = consistency_graph.components.get(&block).unwrap();
    
    writer.write_all("graph {".as_bytes())?;
    for var in variants {
        for var2 in consistency_graph.adjacency_list.get(var).unwrap() {
            //if var.abs() < var2.abs() {
                
            writer.write_all(format!("\t{} -- {};\n",var.abs(), var2.abs()).as_bytes())?;
            //}
        }
    }
    writer.write_all("}".as_bytes())?;
    Ok(())
}

fn make_fastqs(params: &Params, long_read_assignment: &HashMap<i32, i32>, 
    phasing_block_sizes: &HashMap<i32, usize>, molecules: &Molecules) {
    eprintln!("making fastq");
    let mut molid_offsets: Vec<i32> = Vec::new();
    let mut last_offset = 0;
    for offset in molecules.longread_molid_offsets.iter() {
        last_offset = *offset;
        molid_offsets.push(*offset);
    }
    let reads = &params.longread_fqs;
    //write_2_zeros(); // delimiter, this file format isnt great

    let mut sizes: Vec<(&i32, &usize)> = phasing_block_sizes.iter().collect();
    sizes.sort_by(|a, b| b.1.cmp(a.1));


    let mut to_iterate: Vec<(usize, String, i32)> = Vec::new();
    for (molid_offset, (filenum, read)) in izip!(molid_offsets, reads.iter().enumerate()) {
        to_iterate.push((filenum, read.to_string(), molid_offset));
    }

    let mut block_molecules: HashMap<i32, Vec<i32>> = HashMap::new();
    for (mol, block) in long_read_assignment.iter() {
        let mols = block_molecules.entry(block.abs()).or_insert(Vec::new());
        mols.push(*mol);
    }

    let mut block_mols_list: Vec<(i32, Vec<i32>)> = Vec::new();
    let mut dropped = 0;
    for (block, mols) in block_molecules {
        //if mols.len() > 10 {
            if mols.len() < 20 {
                eprintln!("small block {}",block);
            }
            block_mols_list.push((block, mols));
        //} else {
        //    dropped += mols.len();
        //}
    }
    eprintln!("we dropped {} long reads to small phase blocks", dropped);
    eprintln!("we have {} phase blocks", block_mols_list.len());
    // now we are gonna have to do a 2 pass writing phase to not open too many files at once
    // we will first assign each phase block to a meta_block
    // we will then go through the fastqs and write 200 new fastqs of those meta blocks
    // we can then go through each meta block file and write each block out separately
    block_mols_list.sort_by(|a,b| b.1.len().cmp(&a.1.len()));
    
    let mut block_of_blocks: HashMap<i32, i32> = HashMap::new(); // mapping from phase block to meta block
    let num_block_of_blocks = 200;
    let mut meta_block_blocks: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut meta_block = 0;
    for (block, _) in block_mols_list.iter() {
        block_of_blocks.insert(*block, meta_block);
        let mols = meta_block_blocks.entry(meta_block).or_insert(Vec::new());
        mols.push(*block);
        meta_block += 1;
        if meta_block == num_block_of_blocks { meta_block = 0; }
    }
    let mut new_mol_ids: HashMap<i32, HashMap<i32, i32>> = HashMap::new(); // mapping from new molid to old molid
    let mut current_molid: HashMap<i32, i32> = HashMap::new();
    for meta_block in 0..num_block_of_blocks {
        current_molid.insert(meta_block, 0);
    }

    
    for (_filenum, read_file, molid_offset) in to_iterate.iter() {
        let mut writers: HashMap<i32, BufWriter<File>> = HashMap::new();
        for meta_block in 0..num_block_of_blocks {
            let writer1 = OpenOptions::new().write(true).create(true).open(
                format!("{}/longreads_metablock_{}.fq", params.output, meta_block))
                .expect("Unable to create file");
            let writer1 = BufWriter::new(writer1);
            writers.insert(meta_block, writer1);
        }
        
        let mut reader = get_reader(read_file.to_string());
        let mut read_num = *molid_offset;
        let mut last_read_name = "none".to_string();
        let mut buf1 = vec![]; let mut buf2 = vec![]; let mut buf3 = vec![]; let mut buf4 = vec![];
        loop {
            buf1.clear(); buf2.clear(); buf3.clear(); buf4.clear();
            let bytes = reader.read_until(b'\n', &mut buf1).expect("cannot read longread fastq");
            if bytes == 0 { break; }
            let qname = str::from_utf8(&buf1).unwrap();
            reader.read_until(b'\n', &mut buf2).expect("cannot read longread fastq");
            reader.read_until(b'\n', &mut buf3).expect("cannot read longread fastq");
            reader.read_until(b'\n', &mut buf4).expect("cannot read longread fastq");

            if qname == last_read_name { eprintln!("i do see this"); continue; }
            last_read_name = qname.to_string();
            if let Some(comp) = long_read_assignment.get(&(read_num+1)) {
                if let Some(meta_block) = block_of_blocks.get(&comp.abs()) {
                    let molids = new_mol_ids.entry(*meta_block).or_insert(HashMap::new());
                    let molid = *current_molid.get(meta_block).unwrap();
                    molids.insert(molid, read_num+1);
                    current_molid.insert(*meta_block, molid + 1);

                    if let Some(writer) = writers.get_mut(meta_block) {
                        writer.write_all(&buf1).expect("couldnt write");
                        writer.write_all(&buf2).expect("couldnt write");
                        writer.write_all(&buf3).expect("couldnt write");
                        writer.write_all(&buf4).expect("couldnt write");
                    }
                }
            }
            read_num += 1;
        }
        eprintln!("last longread {}",read_num-last_offset);
    }

    for meta_block in 0..num_block_of_blocks {
        let new_mol_ids = new_mol_ids.get(&meta_block).unwrap();
        let mut reader = get_reader(format!("{}/longreads_metablock_{}.fq", params.output, meta_block));
        let mut writers: HashMap<i32, BufWriter<File>> = HashMap::new();
        for phase_block in meta_block_blocks.get(&meta_block).unwrap() {
            let writer1 = OpenOptions::new().write(true).create(true).open(
                format!("{}/longreads_phase_block_{}_hap_{}.fq", params.output, phase_block, 0))
                .expect("Unable to create file");
            let writer1 = BufWriter::new(writer1);
            writers.insert(*phase_block, writer1);
            let writer2 = OpenOptions::new().write(true).create(true).open(
                format!("{}/longreads_phase_block_{}_hap_{}.fq", params.output, phase_block, 1))
                .expect("Unable to create file");
            let writer2 = BufWriter::new(writer2);
            writers.insert(-*phase_block, writer2);
            println!("making writers metablock {} phase block {}", meta_block, phase_block);
        }
        let mut read_num = 0;
        let mut buf1 = vec![];
        let mut buf2 = vec![];
        let mut buf3 = vec![];
        let mut buf4 = vec![];
        loop {
            buf1.clear(); buf2.clear(); buf3.clear(); buf4.clear();
            let bytes = reader.read_until(b'\n', &mut buf1).expect("cannot read longread fastq");
            if bytes == 0 { break; }
            reader.read_until(b'\n', &mut buf2).expect("cannot read longread fastq");
            reader.read_until(b'\n', &mut buf3).expect("cannot read longread fastq");
            reader.read_until(b'\n', &mut buf4).expect("cannot read longread fastq");
            if let Some(molid) = new_mol_ids.get(&read_num) {
                if let Some(comp) = long_read_assignment.get(molid) {
                    if let Some(writer) = writers.get_mut(comp) {
                        writer.write_all(&buf1).expect("couldn't write");
                        writer.write_all(&buf2).expect("couldn't write");
                        writer.write_all(&buf3).expect("couldn't write");
                        writer.write_all(&buf4).expect("couldn't write");
                    } else { println!("there was no writer for component {} in metablock {}", comp, meta_block); }
                } else { println!("there was no component for molid {} in metablock {}", molid, meta_block); }
            } else { println!("there was no molid for read num {} in metablock {}", read_num, meta_block); }
            read_num += 1;
        }
    }
    for meta_block in 0..num_block_of_blocks {
        fs::remove_file(format!("{}/longreads_metablock_{}.fq", params.output, meta_block));
    }
}

fn assign_long_reads2(phasing: &Phasing, molecules: &Molecules, hom_assignments: &HashMap<i32, i32>) -> HashMap<i32, i32> {
    let mut molecule_assignments: HashMap<i32, i32> = HashMap::new();
    let mut total_mols = 0;
    for mol in molecules.get_long_read_molecules() {
        total_mols += 1;
        let mut block_counts: HashMap<i32, [u32; 2]> = HashMap::new();
        let mut total = 0;
        for var in molecules.get_long_read_variants(*mol) {
            if let Some(block) = phasing.variant_phasing.get(&var.abs()) {
                let counts = block_counts.entry(block.abs()).or_insert([0;2]);
                if (*block < 0 && *var > 0) || (*block > 0 && *var < 0) {
                    counts[1] += 1;
                } else { counts[0] += 1; }
            }
            total += 1;
        }
        let mut best_block = -1;
        let mut best_block_count: u32 = 0;
        let mut second_best_block = -2;
        let mut second_best_block_count: u32 = 0;
        for (block, counts) in block_counts.iter() {
            let count = counts[0]+counts[1];
            if count > best_block_count {
                second_best_block = best_block;
                second_best_block_count = best_block_count;
                best_block = *block;
                best_block_count = count;
            } else if count > second_best_block_count {
                second_best_block = *block;
                second_best_block_count = count;
            }
        }
        let mut hom_best_block = -1;
        let mut hom_best_count = 0;
        let mut hom_block_counts: HashMap<i32, u32> = HashMap::new();
        for hom in molecules.get_hom_long_read_variants(mol) {
            if let Some(block) = hom_assignments.get(hom) {
                let count = hom_block_counts.entry(*block).or_insert(0);
                *count += 1;
            }
        }
        for (block, count) in hom_block_counts {
            if count > hom_best_count {
                hom_best_block = block;
                hom_best_count = count;
            }
        }
        if best_block != -1 && hom_best_block != -1 {
            let counts = block_counts.get(&best_block).unwrap();
            //let mut phase = 1;
            if counts[0] > counts[1] {
                molecule_assignments.insert(*mol, hom_best_block);
            } else {
                molecule_assignments.insert(*mol, -hom_best_block);
                //phase = -1;
            }
            //eprintln!("assigning {} het block {} == {} hom block? {}, counts {:?}, hom_count {}", 
            //    hom_best_block*phase, best_block, hom_best_block, best_block == hom_best_block, counts, hom_best_count);
        } else if best_block != -1 {
            let counts = block_counts.get(&best_block).unwrap();
            //let mut counts2 = [0;2];
            //if second_best_block != -1 {
            //    counts2 = *block_counts.get(&second_best_block).unwrap();
            //}
            //let mut phase = 0;
            //let mut best_count = counts[0];
            //if counts[1] > counts[0] { phase  = 1; best_count = counts[1]; }
            //eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
            //    mol, best_block, hom_best_block, counts[0], counts[1], best_count, phase, 
            //    second_best_block, counts2[0], counts2[1]);
            if counts[0] > counts[1] {
                molecule_assignments.insert(*mol, best_block);
            } else {
                molecule_assignments.insert(*mol, -best_block);
            }
        } else if hom_best_block != -1 {
            if mol % 2 == 0 {
                molecule_assignments.insert(*mol, hom_best_block);
            //    eprintln!("no het assingment, hom {}\t{}\t{}",mol, hom_best_block, hom_best_count);
            } else {
                molecule_assignments.insert(*mol, -hom_best_block);
            //    eprintln!("no het assingment, hom {}\t{}\t{}",mol, -hom_best_block, hom_best_count);
            }
        } else { eprintln!("\t\tNO BEST BLOCK with total variants {}", total); }
    }
    eprintln!("total long read molecules {}", total_mols);
    molecule_assignments
}

#[allow(dead_code)]
fn assign_long_reads(consistency_graph: &PhasingConsistencyGraph, molecules: &Molecules) -> HashMap<i32, i32> {
    let mut molecule_assignments: HashMap<i32, i32> = HashMap::new();
    eprintln!("molecule\tphase_block\thap1\thap2\tbest_hap\tphase\tsecond_best_block\tsecond_best_hap1\tsecond_best_hap2");
    for mol in molecules.get_long_read_molecules() {
        let mut block_counts: HashMap<i32, [u32;2]> = HashMap::new();
        //eprintln!("mol {}", mol);
        for var in molecules.get_long_read_variants(*mol) {
            let block = consistency_graph.variant_component.get(&var.abs()).unwrap();
            //let block_size = consistency_graph.components.get(&block).unwrap().len();
            let phase = consistency_graph.phasing_labels.get(&var.abs()).unwrap();
            //eprintln!("\tvar\t{}\tblock\t{}\tsize\t{}\tphase\t{}", var, block, block_size, phase);
            let counts = block_counts.entry(*block).or_insert([0;2]);
            if *phase {
                if *var > 0 {
                    counts[0] += 1;
                } else { counts[1] += 1; }
            } else {
                if *var > 0 {
                    counts[1] += 1;
                } else { counts[0] += 1; }
            }
        }
        let mut best_block = -1;
        let mut best_block_count: u32 = 0;
        let mut second_best_block = -2;
        let mut second_best_block_count: u32 = 0;
        for (block, counts) in block_counts.iter() {
            let count = counts[0]+counts[1];
            if count > best_block_count {
                second_best_block = best_block;
                second_best_block_count = best_block_count;
                best_block = *block;
                best_block_count = count;
            } else if count > second_best_block_count {
                second_best_block = *block;
                second_best_block_count = count;
            }
        }
        //eprintln!("biggest block has {} het variants", best_block_count);
        if best_block != -1 {
            let counts = block_counts.get(&best_block).unwrap();
            let mut counts2 = [0;2];
            if second_best_block != -1 {
                counts2 = *block_counts.get(&second_best_block).unwrap();
            }
            let mut phase = 0;
            let mut best_count = counts[0];
            if counts[1] > counts[0] { phase  = 1; best_count = counts[1]; }
            eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                mol, best_block, counts[0], counts[1], best_count, phase, 
                second_best_block, counts2[0], counts2[1]);
            if counts[0] > counts[1] {
                molecule_assignments.insert(*mol, best_block);
            } else {
                molecule_assignments.insert(*mol, -best_block);
            }
        } else { eprintln!("\t\tNO BEST BLOCK"); }
    }
    molecule_assignments
}

fn assign_homozygous(phasing: &Phasing, 
    _variants: &Variants, molecules: &Molecules) -> HashMap<i32, i32> {
    eprintln!("assigning homogyzous kmers to phase blocks");
    let mut homozygous_phase_blocks: HashMap<i32, i32> = HashMap::new();
    let mut homozygous_phase_counts: HashMap<i32, HashMap<i32, u32>> = HashMap::new();

    let mut counts: HashMap<i32, u32> = HashMap::new();
    for mol in molecules.get_long_read_molecules() {
        counts.clear();
        for var in molecules.get_long_read_variants(*mol) {
            let var = var.abs();
            if let Some(block) = phasing.variant_phasing.get(&var) {
                let block = block.abs();
                let count = counts.entry(block).or_insert(0);
                *count += 1;
            }
        }    
        for hom in molecules.get_hom_long_read_variants(mol) {
            for (block, count) in counts.iter() {
                let homcounts = homozygous_phase_counts.entry(*hom).or_insert(HashMap::new());
                let homcount = homcounts.entry(*block).or_insert(0);
                *homcount += *count;
            }
        }
    }

    for (hom_var, counts) in homozygous_phase_counts {
        let mut best_count = 0;
        let mut best_block = 0;
        let mut second_count = 0;
        for (block, count) in counts {
            if count > best_count {
                second_count = best_count;
                best_count  = count;
                best_block = block;
            }
        }
        best_count += 1;
        let portion = (best_count as f32) / ((best_count as f32) + (second_count as f32));
        //eprintln!("hom var likes block  {} with {}%", best_block, portion);
        if portion > 0.90 {
            homozygous_phase_blocks.insert(hom_var, best_block);
        }
    }
    
    homozygous_phase_blocks
}

#[allow(dead_code)]
fn create_consistency_graph(variants: &Variants, molecules: &Molecules) -> PhasingConsistencyGraph {
    let min_support = 9.0;
    eprintln!("{}",variants.len());
    //println!("graph tenx {{");
    eprintln!("cis1\tcis2\ttrans1\ttrans2\tcis_percent\ttrans_percent\tminor_allele_fraction\tphasing_consistent\tphasing_inconsistent");
    /*
    for var in 1..variants.len() {
        let mut linkedreadgraph: HashMap<(i32, i32), [u16; 4]> = HashMap::new();
        let var = var as i32;
        for mol in variants.get_linked_read_molecules(var) {
            for var2 in molecules.get_linked_read_variants(mol.abs()) {
                if var2.abs() > var {
                    if mol > &0 && var2 > &0 {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[0] += 1;
                    } else if mol < &0 && var2 < &0 {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[1] += 1;
                    } else if mol > &0 && var2 < &0 {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[2] += 1; 
                    } else {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[3] += 1; 
                    }
                }
            }
        }
        /*
        for mol in variants.get_hic_molecules(var) {
            for var2 in molecules.get_hic_variants(mol.abs()) {
                if var2.abs() > var {
                    if (mol > &0 && var2 > &0) {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[0] += 1;
                    } else if (mol < &0 && var2 < &0) {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[1] += 1;
                    } else if (mol > &0 && var2 < &0) {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[2] += 1; 
                    } else {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[3] += 1; 
                    }
                }
            }
        }
        */
        /*
        for mol in variants.get_long_read_molecules(var) {
            for var2 in molecules.get_long_read_variants(mol.abs()) {
                if var2.abs() > var {
                    if mol > &0 && var2 > &0 {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[0] += 1;
                    } else if mol < &0 && var2 < &0 {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[1] += 1;
                    } else if mol > &0 && var2 < &0 {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[2] += 1; 
                    } else {
                        let counts = linkedreadgraph.entry((var, var2.abs())).or_insert([0;4]);
                        counts[3] += 1; 
                    }
                }
            }
        }
        */
        for ((var1, var2), counts) in linkedreadgraph.iter() {
            let counts1 = (counts[0] + counts[1]) as f32;
            let counts2 = (counts[2] + counts[3]) as f32;
            let total = (counts1 + counts2) as f32;
            /*
            if total > 2.0 {
            let mut minor_allele_fraction = 0.0;
            let mut phasing_consistent = 0;
            let mut phasing_inconsistent = 0;
            if counts1 > counts2 {
                minor_allele_fraction = (counts[0].min(counts[1]) as f32)/counts1;
                phasing_consistent = counts1 as usize;
                phasing_inconsistent = counts2 as usize;
            } else if counts2 >= counts1 {
                minor_allele_fraction = (counts[2].min(counts[3]) as f32)/counts2;
                phasing_consistent = counts2 as usize;
                phasing_inconsistent = counts1 as usize;
            }
            //if var1 % 3 == 0 && var2 % 3 == 0 {
            //eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            //     counts[0],counts[1],counts[2],counts[3], counts1/total, counts2/total, 
            //    minor_allele_fraction, phasing_consistent, phasing_inconsistent);
            //}
            }
            */
            if (counts1.max(counts2) as f32)/((counts1+counts2) as f32) > 0.975 &&
                counts1.max(counts2) >= min_support {
                if counts1 >= min_support && ((counts[0] as f32)/(counts1 as f32) < 0.2 || (counts[1] as f32)/(counts1 as f32) < 0.2) {
                    //eprintln!("fake good variant pair {:?}", counts);
                    homozygous.insert(*var1);
                    homozygous.insert(*var2);
                    continue;
                } else if counts2 >= min_support && ((counts[2] as f32)/(counts2 as f32) < 0.2 || (counts[3] as f32)/(counts2 as f32) < 0.2) {
                    //eprintln!("fake good variant pair {:?}", counts);
                    homozygous.insert(*var1);
                    homozygous.insert(*var2);
                    continue;
                } else {
                    //println!("n{} -- n{} [weight={}];", var1, var2, counts1.max(counts2) - counts1.min(counts2));
                    {
                        let adjlist1 = adjacency_list.entry(*var1).or_insert(Vec::new());
                        if counts1 > counts2 {
                            adjlist1.push(*var2);
                        } else {
                            adjlist1.push(-var2);
                        }
                    }
                    {
                        let adjlist2 = adjacency_list.entry(*var2).or_insert(Vec::new());
                        if counts1 > counts2 {
                            adjlist2.push(*var1);
                        } else {
                            adjlist2.push(-var1);
                        }
                    }
                    //eprintln!("probably good pair {:?}",counts);
                }
            }
            if counts1 + counts2 > 5.0 {
                full_graph.insert((*var1, *var2), [counts[0], counts[1], counts[2], counts[3]]);
            }//else {
            //    eprintln!("bad variant pair {:?}", counts);
            //}
        }
    }
    //});
    
    eprintln!("homozygous variants {}",homozygous.len());

    let (components, var_comps) = graph_components(&full_graph, variants, min_support as u16);
    
    let labels = breadth_first_search_phase_label(&components, &adjacency_list);
    label_ends(&components, &adjacency_list);
    */
    let components: HashMap<i32, Vec<i32>> = HashMap::new();
    let thresholds = PhasingConsistencyThresholds{
        phasing_consistency: 0.975,
        min_support: min_support as f32,
        minor_allele_fraction: 0.2,
    };
    let phasing_consistency_graph = make_graph(&variants, &molecules, &thresholds);
    //check_labels(&labels, &adjacency_list, &var_comps, &components);
    /*
    let mut comp_mols: HashMap<usize, Vec<i32>> = HashMap::new();
    let mut mol_comps: HashMap<usize, Vec<i32>> = HashMap::new();
    for (comp1, vars1) in components.iter() {
        let mut molcounts1: HashMap<i32, [u8;2]> = HashMap::new();
        for var in vars1 {
            let varlabel = labels.get(var).unwrap();
            for mol in variants.get_linked_read_molecules(*var) {
                let mut mol = *mol;
                if !varlabel {
                    mol = -mol;
                }
                if mol > 0 {
                    let counts = molcounts1.entry(mol.abs()).or_insert([0;2]);
                    counts[0] += 1;
                } else {
                    let counts = molcounts1.entry(mol.abs()).or_insert([0;2]);
                    counts[1] += 1;
                }
            }
        }
        molcounts1.retain(|_key, value| 
            value[0]+value[1] > 0 && 
            (value[0].max(value[1]) as f32)/((value[0]+value[1]) as f32) > 0.9);  
        let mut mols: Vec<i32> = Vec::new();
        for (mol, counts) in molcounts1.iter() {
            if counts[0] > counts[1] { 
                mols.push(*mol); 
                let comps = mol_comps.entry(mol).or_insert(Vec::new());
                comps.push(*comp1);
            } 
            else { 
                mols.push(-mol); 
                let comps = mol_comps.entry(mol).or_insert(Vec::new());
                comps.push(-comp1);
            }
        }
        comp_mols.insert(*comp1, mols);
        
    }
    */
    
    eprintln!("got here");
    let thresholds = PhasingConsistencyThresholds{
        phasing_consistency: 0.95,
        min_support: 20.0,
        minor_allele_fraction: 0.2,
    };
    let phasing_consistency_graph = grow_graph(&phasing_consistency_graph, &variants, &molecules, &thresholds);

    let hypergraph: HashMap<(i32, i32), [u32; 4]> = HashMap::new();
    let hypergraph_edges: HashMap<(i32, i32), Vec<(i32, i32, bool)>> = HashMap::new();

    println!("graph tenx {{");
    for (e1, e2, weight) in phasing_consistency_graph.edges.iter() {
        println!("\tn{} -- n{} [weight={}];",e1,e2,weight);
    }
    /*
    for var in 1..variants.len() {
        let mut var = var as i32;
        let varlabel = labels.get(&var).unwrap();
        let mut varcomp = *var_comps.get(&var).unwrap();
        //for mol in variants.get_hic_molecules(var) {
        for mol in variants.get_linked_read_molecules(var) {
            let mut mol = *mol;
            //for var2 in molecules.get_hic_variants(mol.abs()) {
            for var2 in molecules.get_linked_read_variants(mol.abs()) {
                let mut var2 = *var2;
                let mut var2comp = *var_comps.get(&var2.abs()).unwrap();
                if var2comp > varcomp {//var2.abs() > var {
                    let var2label = labels.get(&var2.abs()).unwrap();
                    //if varcomp == var2comp { continue; }
                    //if var2comp < varcomp {
                    //    let tmp = varcomp;
                    //    varcomp = var2comp;
                    //    var2comp = tmp;
                    //}
                    //if mol < 0 {var = -var;}
                    var = var.abs();
                    if !varlabel {
                        mol = -mol;
                    }
                    if !var2label { 
                        var2 = -var2;
                    }
                    if mol > 0 && var2 > 0 {
                        let counts = hypergraph.get_mut(&(varcomp, var2comp)).unwrap();//.or_insert([0;4]);
                        counts[0] += 1;
                        //let edgelist = hypergraph_edges.get_mut(&(varcomp, var2comp)).unwrap();//.or_insert(Vec::new());
                        //edgelist.push((var, var2.abs(), true));
                    } else if mol < 0 && var2 < 0 {
                        let counts = hypergraph.get_mut(&(varcomp, var2comp)).unwrap();//.or_insert([0;4]);
                        counts[1] += 1;
                        //let edgelist = hypergraph_edges.get_mut(&(varcomp, var2comp)).unwrap();//.or_insert(Vec::new());
                        //edgelist.push((var, var2.abs(), true));
                    } else if mol > 0 && var2 < 0  {
                        let counts = hypergraph.get_mut(&(varcomp, var2comp)).unwrap();//.or_insert([0;4]);
                        counts[2] += 1; 
                        //let edgelist = hypergraph_edges.get_mut(&(varcomp, var2comp)).unwrap();//.or_insert(Vec::new());
                        //edgelist.push((var, var2.abs(), false));
                    } else {
                        let counts = hypergraph.get_mut(&(varcomp, var2comp)).unwrap();//.or_insert([0;4]);
                        counts[3] += 1; 
                        //let edgelist = hypergraph_edges.get_mut(&(varcomp, var2comp)).unwrap();//.or_insert(Vec::new());
                        //edgelist.push((var, var2.abs(), false));
                    }
                }
            }
        }
    }
    */
    eprintln!("cis1\tcis2\ttrans1\ttrans2\tcis_percent\ttrans_percent\tminor_allele_fraction\tphasing_consistent\tphasing_inconsistent\tsize1\tsize2");
    for ((comp1, comp2), counts) in hypergraph.iter() {
        let counts1 = (counts[0] + counts[1]) as f32;
        let counts2 = (counts[2] + counts[3]) as f32;
        let total = (counts1+counts2) as f32;

        let vars1 = components.get(comp1).unwrap();
        let vars2 = components.get(comp2).unwrap();
        if total > 8.0 {
            
            let mut minor_allele_fraction = 0.0;
            let mut phasing_consistent = 0;
            let mut phasing_inconsistent = 0;
            if counts1 > counts2 {
                minor_allele_fraction = (counts[0].min(counts[1]) as f32)/counts1;
                phasing_consistent = counts1 as usize;
                phasing_inconsistent = counts2 as usize;
            } else if counts2 > counts1 {
                minor_allele_fraction = (counts[2].min(counts[3]) as f32)/counts2;
                phasing_consistent = counts2 as usize;
                phasing_inconsistent = counts1 as usize;
            }
            eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                counts[0],counts[1],counts[2],counts[3], counts1/total, counts2/total, 
                minor_allele_fraction, phasing_consistent, phasing_inconsistent, vars1.len(), vars2.len());
        }
        if counts1 > 150.0 && counts1/total > 0.8 && (counts[0].min(counts[1]) as f32)/counts1 > 0.2 {
            //eprintln!("possible cis match, comps {} and {} with sizes {} and {} have counts {:?} {}",
            //    comp1, comp2, vars1.len(), vars2.len(), counts, counts1/total);
            let edgelist = hypergraph_edges.get(&(*comp1, *comp2)).unwrap();
            for (var1, var2, _thingy) in edgelist.iter() {
                println!("n{} -- n{} [weight=1,color=\"green\"];",var1, var2);
            }
        } else if counts2 > 150.0 && counts2/total > 0.8 && (counts[2].min(counts[3]) as f32)/counts2 > 0.2 {
            //eprintln!("possible trans match, comps {} and {} with sizes {} and {} have counts {:?} {}",
            //    comp1, comp2, vars1.len(), vars2.len(), counts, counts2/total);
            //let edgelist = hypergraph_edges.get(&(*comp1, *comp2)).unwrap();
            //for (var1, var2, _) in edgelist.iter() {
            //    println!("n{} -- n{} [weight=1, color=\"green\"];",var1, var2);
            //}
            let edgelist = hypergraph_edges.get(&(*comp1, *comp2)).unwrap();
            for (var1, var2, _thingy) in edgelist.iter() {
                println!("n{} -- n{} [weight=1,color=\"green\"];",var1, var2);
            }
        } //else if total > 50.0 {
            //eprintln!("bad hic, comps {} and {} with sizes {} and {} have counts {:?} {}",
            //    comp1, comp2, vars1.len(), vars2.len(), counts, counts1.max(counts2)/total);
        //}
    }


    println!("}}");

    phasing_consistency_graph
}

#[allow(dead_code)]
fn make_graph(variants: &Variants, _molecules: &Molecules, 
    phasing_consistency: &PhasingConsistencyThresholds) -> 
    PhasingConsistencyGraph {
    let mut var_mols: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut mol_vars: HashMap<i32, Vec<i32>> = HashMap::new();
    for var in 1..variants.len() {
        let mut mols: Vec<i32> = Vec::new();
        for mol in variants.get_long_read_molecules(var as i32) {
            mols.push(*mol);
            let vars = mol_vars.entry(mol.abs()).or_insert(Vec::new());
            if mol > &0 { vars.push(var as i32); }
            else {vars.push(-(var as i32)); }
        }
        var_mols.insert(var as i32, mols);
    }
    
    let (edges, adjacency_list) = graph_core(&var_mols, &mol_vars, phasing_consistency);
    let (new_components, var_comps) = graph_components2(&edges, variants);
    let labels = breadth_first_search_phase_label(&new_components, &adjacency_list);
    PhasingConsistencyGraph {
        edges: edges,
        adjacency_list: adjacency_list,
        components: new_components,
        variant_component: var_comps,
        phasing_labels: labels,
    }
}

struct PhasingConsistencyThresholds {
    min_support: f32,
    minor_allele_fraction: f32,
    phasing_consistency: f32,
}

#[allow(dead_code)]
fn graph_core(comp_mols: &HashMap<i32, Vec<i32>>, mol_comps: &HashMap<i32, Vec<i32>>,
    phasing_consistency: &PhasingConsistencyThresholds) -> 
    (HashSet<(i32, i32, u32)>, HashMap<i32, Vec<i32>>) {
    let mut edges: HashSet<(i32, i32, u32)> = HashSet::new();
    let mut adjacency_list: HashMap<i32, Vec<i32>> = HashMap::new();
    for comp1 in comp_mols.keys() {
        let mut graphcount: HashMap<i32, [u32;4]> = HashMap::new();
        for mol in comp_mols.get(comp1).unwrap() {
            for comp2 in mol_comps.get(&mol.abs()).unwrap() {
                if *comp1 >= comp2.abs() { continue; }
                if mol > &0 && comp2 > &0 {
                    let count = graphcount.entry(comp2.abs()).or_insert([0;4]);
                    count[0] += 1;
                } else if mol < &0 && comp2 < &0 {
                    let count = graphcount.entry(comp2.abs()).or_insert([0;4]);
                    count[1] += 1;
                } else if mol > &0 && comp2 < &0 {
                    let count = graphcount.entry(comp2.abs()).or_insert([0;4]);
                    count[2] += 1;
                } else if mol < &0 && comp2 > &0 {
                    let count = graphcount.entry(comp2.abs()).or_insert([0;4]);
                    count[3] += 1;
                }
            }
        }
        for (comp2, counts) in graphcount.iter() {
            let cis = (counts[0] + counts[1]) as f32;
            let trans = (counts[2] + counts[3]) as f32;
            let total = cis + trans;
            if cis.max(trans)/total >= phasing_consistency.phasing_consistency && 
                total >= phasing_consistency.min_support {
                if (cis > trans && (counts[0].min(counts[1]) as f32)/total >= phasing_consistency.minor_allele_fraction) ||
                    (trans > cis && (counts[2].min(counts[3]) as f32)/total >= phasing_consistency.minor_allele_fraction) {
                    //eprintln!("components {} and {} sizes {}, {} counts {:?}", comp1, comp2, 
                    //    components.get(comp1).unwrap().len(), components.get(comp2).unwrap().len(), counts);
                    if (*comp1 == 3951818) ||
                       (*comp2 == 3951818) {
                        eprintln!("BAD EDGE, {} -- {} with counts {:?}", comp1, comp2, counts);
                    } 
                    edges.insert((*comp1, *comp2, (cis.max(trans) - cis.min(trans)) as u32));
                    {
                        let varedges = adjacency_list.entry(*comp1).or_insert(Vec::new());
                        if cis > trans {
                            varedges.push(*comp2);
                        } else {
                            varedges.push(-comp2);
                        }
                    }
                    {
                        let varedges = adjacency_list.entry(*comp2).or_insert(Vec::new());
                        if cis > trans {
                            varedges.push(*comp1);
                        } else {
                            varedges.push(-comp1);
                        }
                    }
                }
            }
        }
    }
    (edges, adjacency_list)
}


#[allow(dead_code)]
fn grow_graph(putative_graph: &PhasingConsistencyGraph, variants: &Variants, 
    _molecules: &Molecules, phasing_consistency: &PhasingConsistencyThresholds) -> PhasingConsistencyGraph {
    let mut comp_mols: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut mol_comps: HashMap<i32, Vec<i32>> = HashMap::new();
    for (comp1, vars1) in putative_graph.components.iter() {
        let mut molcounts1: HashMap<i32, [u16;2]> = HashMap::new();
        for var in vars1 {
            assert!(putative_graph.phasing_labels.contains_key(var));
            let varlabel = putative_graph.phasing_labels.get(var).unwrap();
            for mol in variants.get_long_read_molecules(*var) {
                let mut mol = *mol;
                if !varlabel {
                    mol = -mol;
                }
                if mol > 0 {
                    let counts = molcounts1.entry(mol.abs()).or_insert([0;2]);
                    counts[0] += 1;
                } else {
                    let counts = molcounts1.entry(mol.abs()).or_insert([0;2]);
                    counts[1] += 1;
                }
            }
        }
        molcounts1.retain(|_key, value| value[0]+value[1] > 0 && 
            (value[0].max(value[1]) as f32)/((value[0]+value[1]) as f32) > phasing_consistency.phasing_consistency);  
        let mut mols: Vec<i32> = Vec::new();
        for (mol, counts) in molcounts1.iter() {
            if counts[0] > counts[1] { 
                mols.push(*mol); 
                let comps = mol_comps.entry(*mol).or_insert(Vec::new());
                comps.push(*comp1 as i32);
            } 
            else { 
                mols.push(-mol); 
                let comps = mol_comps.entry(*mol).or_insert(Vec::new());
                comps.push(-(*comp1 as i32));
            }
        }
        comp_mols.insert(*comp1, mols);
    }
    let (edges, adjacency_list) = graph_core(&comp_mols, &mol_comps, phasing_consistency);
    let mut full_edges: HashSet<(i32, i32, u32)> = HashSet::new();
    let mut full_adjacency_list: HashMap<i32, Vec<i32>> = HashMap::new();
    for (n1, n2, weight) in putative_graph.edges.iter() {
        full_edges.insert((*n1, *n2, *weight));
    }
    for (n1, n2, weight) in edges.iter() {
        full_edges.insert((*n1, *n2, *weight));
    }
    for (var1, vars) in putative_graph.adjacency_list.iter() {
        let mut vars2: Vec<i32> = Vec::new();
        for var in vars.iter() {
            vars2.push(*var);
        }
        full_adjacency_list.insert(*var1, vars2);
    }
    for (var, vars) in adjacency_list {
        for var2 in vars {
            let putative = full_adjacency_list.entry(var).or_insert(Vec::new());
            putative.push(var2);
        }
    }
    let (new_components, var_comps) = graph_components2(&full_edges, variants);
    let labels = breadth_first_search_phase_label(&new_components, &full_adjacency_list);
    PhasingConsistencyGraph {
        edges: full_edges,
        adjacency_list: full_adjacency_list,
        components: new_components,
        variant_component: var_comps,
        phasing_labels: labels,
    }
}
#[allow(dead_code)]
struct PhasingConsistencyGraph {
    edges: HashSet<(i32, i32, u32)>,
    adjacency_list: HashMap<i32, Vec<i32>>,
    components: HashMap<i32, Vec<i32>>,
    variant_component: HashMap<i32, i32>,
    phasing_labels: HashMap<i32, bool>,
}

#[allow(dead_code)]
fn graph_components2(graph: &HashSet<(i32,i32, u32)>, variants: &Variants) ->
    (HashMap<i32, Vec<i32>>, HashMap<i32, i32>) {
    let mut union_find: Vec<usize> = Vec::new();//(0..variants.len()).collect();
    for i in 0..variants.len() { union_find.push(i); }
    for (var1, var2, _) in graph.iter() {
        let mut cur = *var1 as usize;
        while cur != union_find[cur] {
            union_find[cur] = union_find[union_find[cur]];
            cur = union_find[cur];
        }
        let mut cur2 = *var2 as usize;
        while cur2 != union_find[cur2] {
            union_find[cur2] = union_find[union_find[cur2]];
            cur2 = union_find[cur2];
        }
        union_find[*var2 as usize] = cur.min(cur2);
        union_find[*var1 as usize] = cur.min(cur2);
        union_find[cur] = cur.min(cur2);
        union_find[cur2] = cur.min(cur2);
    }
    for i in (1..variants.len()).rev() {
        let mut place = i;
        let mut themall: Vec<usize> = Vec::new(); //it's "them all", a southernism 
        while place != union_find[place] {
            themall.push(place);
            place = union_find[place];
        }
        for j in themall {
            union_find[j] = place;
        }    
    }
    //let mut variant_comp: HashMap<i32, usize> = HashMap::new();
    let mut components: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut var_comps: HashMap<i32, i32> = HashMap::new();
    for i in 1..variants.len() {
        let component = components.entry(union_find[i].try_into().unwrap()).or_insert(Vec::new());
        component.push(i as i32);
        var_comps.insert(i as i32, union_find[i].try_into().unwrap());
    }
    (components, var_comps)
}

#[allow(dead_code)]
fn graph_components(full_graph: &HashMap<(i32, i32), [u16; 4]>, variants: &Variants, min_support: u16) -> 
    (HashMap<i32, Vec<i32>>, HashMap<i32, i32>) {
    let mut union_find: Vec<usize> = Vec::new();//(0..variants.len()).collect();
    for i in 0..variants.len() { union_find.push(i); }
    for ((var1, var2), counts) in full_graph.iter() {
        let counts1 = counts[0]+counts[1];
        let counts2 = counts[2]+counts[3];
        if (counts1.max(counts2) as f32)/((counts1+counts2) as f32) > 0.975 &&
                    counts1.max(counts2) >= min_support {
            if counts1 >= min_support && ((counts[0] as f32)/(counts1 as f32) < 0.2 || (counts[1] as f32)/(counts1 as f32) < 0.2) {
                continue;
            } else if counts2 >= min_support && ((counts[2] as f32)/(counts2 as f32) < 0.2 || (counts[3] as f32)/(counts2 as f32) < 0.2) {
                continue;
            } 
            let mut cur = *var1 as usize;
            while cur != union_find[cur] {
                union_find[cur] = union_find[union_find[cur]];
                cur = union_find[cur];
            }
            let mut cur2 = *var2 as usize;
            while cur2 != union_find[cur2] {
                union_find[cur2] = union_find[union_find[cur2]];
                cur2 = union_find[cur2];
            }
            union_find[*var2 as usize] = cur.min(cur2);
            union_find[*var1 as usize] = cur.min(cur2);
            union_find[cur] = cur.min(cur2);
            union_find[cur2] = cur.min(cur2);
        }
    }
    for i in (1..variants.len()).rev() {
        let mut place = i;
        let mut themall: Vec<usize> = Vec::new(); //it's "them all", but fine to interpret as "the mall" as well
        while place != union_find[place] {
            themall.push(place);
            place = union_find[place];
        }
        for j in themall {
            union_find[j] = place;
        }    
    }
    //let mut variant_comp: HashMap<i32, usize> = HashMap::new();
    let mut components: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut var_comps: HashMap<i32, i32> = HashMap::new();
    for i in 1..variants.len() {
        let component = components.entry(union_find[i].try_into().unwrap()).or_insert(Vec::new());
        component.push(i as i32);
        var_comps.insert(i as i32, union_find[i].try_into().unwrap());
    }
    (components, var_comps)
}

#[allow(dead_code)]
fn label_ends(components: &HashMap<i32, Vec<i32>>, adjacency_list: &HashMap<i32, Vec<i32>>) {
    //println!("\t{{");
    for (_comp, complist) in components.iter() {
        let mut visited: HashSet<i32> = HashSet::new();
        let startvar = complist[0];
        let mut priority_queue: VecDeque<i32> = VecDeque::new();
        let mut lastvar = startvar;
        if let Some(_varlist) = adjacency_list.get(&startvar) {
            for var in adjacency_list.get(&startvar).unwrap().iter() {
                priority_queue.push_back(var.abs());
            }
            while priority_queue.len() > 0 {
                let nextvar = priority_queue.pop_front().unwrap().abs();
                if visited.contains(&nextvar) { continue; }
                //eprintln!("nextvar {}",nextfar);
                visited.insert(nextvar);
                for var in adjacency_list.get(&nextvar).unwrap().iter() {
                    priority_queue.push_back(var.abs());
                }
                lastvar = nextvar;
            }
        }
        //println!("\t\tn{} [width=0.5 shape=circle style=filled fillcolor=blue size=100];", lastvar);
        let mut visited: HashSet<i32> = HashSet::new();
        let startvar = lastvar;
        let mut priority_queue: VecDeque<i32> = VecDeque::new();
        //let mut lastvar = startvar;
        if let Some(_varlist) = adjacency_list.get(&startvar) {
            for var in adjacency_list.get(&startvar).unwrap().iter() {
                priority_queue.push_back(var.abs());
            }
            while priority_queue.len() > 0 {
                let nextvar = priority_queue.pop_front().unwrap().abs();
                if visited.contains(&nextvar) { continue; }
                visited.insert(nextvar);
                for var in adjacency_list.get(&nextvar).unwrap().iter() {
                    priority_queue.push_back(var.abs());
                }
                //lastvar = nextvar;
            }
        }
        //println!("\t\tn{} [width=0.5 shape=circle style=filled fillcolor=red size=100];", lastvar);
    }
    //println!("\t}}");
}
#[allow(dead_code)]
fn breadth_first_search_phase_label(components: &HashMap<i32, Vec<i32>>, 
    adjacency_list: &HashMap<i32, Vec<i32>>) -> HashMap<i32, bool> {
    let mut labels: HashMap<i32, bool> = HashMap::new();
    let mut visited: HashSet<i32> = HashSet::new();
    for (_comp, complist) in components.iter() {
        let mut priority_queue: VecDeque<i32> = VecDeque::new();
        let startvar = complist[0];
        labels.insert(startvar, true);
        visited.insert(startvar);
        if let Some(_varlist) = adjacency_list.get(&startvar) {
            for var in adjacency_list.get(&startvar).unwrap().iter() {
                priority_queue.push_back(*var);
            }
            while priority_queue.len() > 0 {
                let nextvar = priority_queue.pop_front().unwrap();
                if visited.contains(&nextvar.abs()) {continue;} // maybe check label consistency
                visited.insert(nextvar.abs());
                if nextvar < 0 { labels.insert(nextvar.abs(), false); } 
                else { labels.insert(nextvar.abs(), true); }
                for var in adjacency_list.get(&nextvar.abs()).unwrap().iter() {
                    if nextvar < 0 {
                        priority_queue.push_back(-var);
                    } else {
                        priority_queue.push_back(*var);
                    }
                }
            }
        }
    }
    labels
}

#[allow(dead_code)]
fn check_labels(labels: &HashMap<i32, bool>, adjacency_list: &HashMap<i32, Vec<i32>>, 
    var_comps: &HashMap<i32, usize>, components: &HashMap<usize, Vec<i32>>) {
    let mut error_counts: HashMap<i32, [usize; 2]> = HashMap::new();
    for (var1, var2s) in adjacency_list.iter() {
        let var1phase = *labels.get(var1).unwrap();
        for var2 in var2s {
            let var2phase = *labels.get(&var2.abs()).unwrap();
            if var2 < &0 {
                if (var1phase || var2phase) && !(var1phase && var2phase) {
                    //eprintln!("good trans variants {} and {} phases {}/{}", var1, var2, var1phase, var2phase);
                    {let counts = error_counts.entry(*var1).or_insert([0;2]);
                    counts[0] += 1;}
                    {let counts = error_counts.entry(var2.abs()).or_insert([0;2]);
                    counts[0] += 1;}
                } else {
                    //eprintln!("bad trans variants {} and {} phases {}/{}", var1, var2, var1phase, var2phase);
                    {let counts = error_counts.entry(*var1).or_insert([0;2]);
                    counts[1] += 1;}
                    {let counts = error_counts.entry(var2.abs()).or_insert([0;2]);
                    counts[1] += 1;}
                }
            } else {
                if (var1phase && var2phase) || (!var1phase && !var2phase) {
                    //eprintln!("good cis variants {} and {} phases {}/{}", var1, var2, var1phase, var2phase);
                    {let counts = error_counts.entry(*var1).or_insert([0;2]);
                    counts[0] += 1;}
                    {let counts = error_counts.entry(*var2).or_insert([0;2]);
                    counts[0] += 1;}
                } else {
                    //eprintln!("bad cis variants {} and {} phases {}/{}", var1, var2, var1phase, var2phase);
                    {let counts = error_counts.entry(*var1).or_insert([0;2]);
                    counts[1] += 1;}
                    {let counts = error_counts.entry(*var2).or_insert([0;2]);
                    counts[1] += 1;}
                }
            }
        }
    }
    for (var, counts) in error_counts.iter() {
        if counts[1] != 0 {
            let comp = var_comps.get(var).unwrap();
            eprintln!("var {} from comp {} with size {} had {} bad counts and {} good counts", 
                var, comp, components.get(comp).unwrap().len(), counts[0], counts[1]);
        }
    }
}

fn get_seed(rng: &mut StdRng) -> [u8; 32] {
    let mut seed: [u8; 32] = [0; 32];
    for i in 0..32 {
        seed[i] = rng.gen();
    }
    seed
}

#[allow(dead_code)]
fn phuzzy_phaser_master(molecules: Molecules, variants: Variants, 
        num_clusters: usize, phasing: Option<HashMap<i32, String>>, 
        info: Option<HashMap<i32, (String, i32, String, String, String)>>,
        kmer_ids_to_kmer: HashMap<i32, String>) -> Vec<PhaseBlock> {
    let mut cluster_support: Vec<f32> = Vec::new(); // to be reused
    let mut phase_blocks: Vec<PhaseBlock> = Vec::new();
    for _cluster in 0..num_clusters { cluster_support.push(0.0); } //initialize
    let mut sums: Vec<Vec<f32>> = Vec::new(); let mut denoms: Vec<Vec<f32>> = Vec::new();
    let vars_per_step = 10;
    for cluster in 0..num_clusters {
        denoms.push(Vec::new()); sums.push(Vec::new());
        for _ in 0..vars_per_step { denoms[cluster].push(0.0); sums[cluster].push(0.0); }
    }
    let mut cluster_centers: Vec<Vec<f32>> = Vec::new();
    let seed: [u8; 32] = [6; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    for cluster in 0..num_clusters {
        cluster_centers.push(Vec::new());
        for _ in 0..variants.len() {
            let x: f32 = rng.gen();
            cluster_centers[cluster].push(x/2.0 + 0.25); // initialize from 0.25..0.75 uniformly 
        }
    }
    let mut touched_vars: Vec<bool> = Vec::new();
    let mut available_vars: HashSet<usize> = HashSet::new();
    for i in 1..variants.len() { 
        touched_vars.push(false);
        available_vars.insert(i);
    }


    let mut first_varset: Vec<i32> = Vec::new();
    let directions: [i32; 2] = [1, -1];
    let mut phase_block_id: usize = 0;
    println!("iteration\tvardex\tchrom\tpos\tref\talt\tLR_ps\tLR_gt\tconnections\tref_counts\talt_counts\tref_bc_counts\talt_bc_counts\thap1\thap2");
    'next_phase_block_loop: while available_vars.len() > 0 {
        println!("phase block {}", phase_block_id);
        phase_block_id += 1;
        let seed = get_seed(&mut rng);
        println!("starting new phase block with seed {:?}",seed);
        let mut variant_window = VariantRecruiter::new(&variants, &molecules, seed);
        println!("{} total variants left", available_vars.len());
        let mut variants_so_far = 0;
        for direction in &directions { 
            let mut total_variants: i32 = 0;
            if *direction == -1 && total_variants == 0 {
                variant_window.reset_except_visited();
                for vardex in &first_varset {
                    variant_window.add_variant(&variants, &molecules, *vardex, &available_vars);
                }
            }
            //pick n starting variants
            'phaseblockloop: loop {
                let mut new_variants: Vec<usize> = Vec::new();
                let mut mol_connections: Vec<u32> = Vec::new();
                let mut mol_percent_connections: Vec<f32> = Vec::new();
                let mut new_iterations: Vec<i32> = Vec::new();
                for _i in 0..vars_per_step {
                    new_iterations.push(total_variants);
                    total_variants += direction;
                    let mut min = 0.4;
                    if variants_so_far < vars_per_step { min = 0.1; }
                    let (variant, _percent_long_read, _count_long_read, percent_10x, count_10x, _count_hic) = 
                            match variant_window.get_next_variant(&variants, min, &available_vars) {
                        Some((x,y,z,a,b,c)) => (x,y,z,a,b,c),
                        None => { 
                            println!("breaking phaseblock");
                            if direction == &-1 {
                                println!("phaseblock size = {}", variants_so_far);
                                phase_blocks.push(variant_window.make_phaseblock(&cluster_centers));
                            }
                            for var in new_variants.iter() {
                                available_vars.insert(*var as usize);
                            }
                            if direction == &1 && variants_so_far < vars_per_step {
                                continue 'next_phase_block_loop;
                            }
                            break 'phaseblockloop },
                    };
                    variants_so_far += 1;
                    //println!("next variant {} of {} and {} connections {} {}", 
                        //variant, kmer_ids_to_kmer.get(&variant).unwrap(), kmer_ids_to_kmer.get(&-variant).unwrap(), 
                        //percent_long_read, count_long_read);
                    // add our new variant
                    //println!("{}\t{}",variant, connections);
                    variant_window.add_variant(&variants, &molecules, variant, &available_vars);
                    if *direction == 1 && total_variants < 10 { first_varset.push(variant); }
                    mol_connections.push(count_10x);
                    mol_percent_connections.push(percent_10x);
                    let vardex = variant.abs() as usize;
                    new_variants.push(vardex);
                }
                new_variants.sort();
                let mut current_molecules: Vec<Vec<i32>> = Vec::new();
                let mut used_mols: HashSet<i32> = HashSet::default();
                for vardex in new_variants.iter() {
                    for molecule in variants.get_linked_read_molecules(*vardex as i32) {//&variants[vardex - 1] {
                        let moldex = molecule.abs();
                        if !used_mols.contains(&moldex) { 
                            let mut mol: Vec<i32> = Vec::new();
                            for var in molecules.get_linked_read_variants(moldex) {
                                let vardex = var.abs();
                                if variant_window.current_variant_molecules.contains_key(&vardex) { 
                                    mol.push(*var);
                                }
                            }
                            current_molecules.push(mol);
                            used_mols.insert(moldex);
                        } 
                    }
                }
                for variant in &new_variants { 
                    touched_vars[*variant] = true; 
                    available_vars.remove(&(*variant as usize));
                }

                next_cluster_step(current_molecules, &mut cluster_centers, &new_variants, &mut cluster_support, &touched_vars, &mut denoms, &mut sums);
                for (index, vardex) in new_variants.iter().enumerate() {
                    let mut chrom: String = "n/a".to_string();
                    let mut pos: i32 = -1; 
                    let mut reference: String = kmer_ids_to_kmer.get(&(*vardex as i32)).unwrap().to_string();
                    let mut alt: String = kmer_ids_to_kmer.get(&-(*vardex as i32)).unwrap().to_string();
                    let mut ps: String = "n/a".to_string();
                    if let Some(some_info) = &info {
                        let (chrom1, pos1, reference1, alt1, ps1) = some_info.get(&(*vardex as i32)).unwrap();
                        chrom = chrom1.to_string(); pos = *pos1; reference = reference1.to_string(); alt = alt1.to_string(); ps = ps1.to_string();
                    }
                    let mut phasing_str: String = "./.".to_string();
                    if let Some(phasing_info) = &phasing {
                        phasing_str = phasing_info.get(&(*vardex as i32)).unwrap().to_string();
                    }
                    println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}", 
                        new_iterations[index], vardex, chrom, pos, reference, alt, ps, 
                        phasing_str, mol_connections[index], mol_percent_connections[index],
                        cluster_centers[0][*vardex], cluster_centers[1][*vardex]);
                }
            }
        }
    }
    phase_blocks
}
#[allow(dead_code)]
fn log_sum_exp(p: &Vec<f32>) -> f32{
    let max_p: f32 = p.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let sum_rst: f32 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

#[allow(dead_code)]
fn next_cluster_step(molecules: Vec<Vec<i32>>, cluster_centers: &mut Vec<Vec<f32>>, 
    new_variants: &Vec<usize>, cluster_support: &mut Vec<f32>, touched: &Vec<bool>,
    denoms: &mut Vec<Vec<f32>>, sums: &mut Vec<Vec<f32>>) {
    let mut vardex; let mut moldex; let mut varval; let mut molval;
    let mut change = 1.0;
    while change > 0.0001 {
        // reset reused memory to prior expectation, think of this like a pseudocount and a way to avoid div by 0
        for cluster in 0..denoms.len() {
            for vardex in 0..denoms[cluster].len() {
                denoms[cluster][vardex] = 0.001;
                sums[cluster][vardex] = 0.001*cluster_centers[cluster][new_variants[vardex]];
            }
        } 
        // Expectation
        for molecule_variants in molecules.iter() {
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
    potential_variants_linked_reads: BTreeSet<(u32, i32)>,
    potential_variant_counts_linked_reads: HashMap<i32, u32>,
    potential_variants_hic: BTreeSet<(u32, i32)>,
    potential_variant_counts_hic: HashMap<i32, u32>,
    potential_variants_long_reads: BTreeSet<(u32, i32)>,
    potential_variant_counts_long_reads: HashMap<i32, u32>,
    iteration: u32,
    rng: StdRng,
}

impl VariantRecruiter { 
    #[allow(dead_code)]
    fn make_phaseblock(&self, cluster_centers: &Vec<Vec<f32>>) -> PhaseBlock {
        PhaseBlock::new(cluster_centers, &self.variants_visited)
    }
    //fn new(variants: &Vec<Vec<i32>>, _molecules: &HashMap<i32, Vec<i32>>) -> VariantRecruiter {
    #[allow(dead_code)]
    fn new(variants: &Variants, _molecules: &Molecules, seed: [u8; 32]) -> VariantRecruiter {
        //let seed: [u8; 32] = [6; 32]; 
        let rng: StdRng = SeedableRng::from_seed(seed);
        VariantRecruiter{
            variants_visited: BitVec::from_elem(variants.len() + 1, false),
            current_molecules_variants: HashMap::default(),
            current_variant_molecules: HashMap::default(),
            molecule_last_seen: HashMap::default(),
            potential_variants_linked_reads: BTreeSet::new(),
            potential_variant_counts_linked_reads: HashMap::default(),
            potential_variants_long_reads: BTreeSet::new(),
            potential_variant_counts_long_reads: HashMap::default(),
            potential_variants_hic: BTreeSet::new(),
            potential_variant_counts_hic: HashMap::default(),
            iteration: 0,
            rng: rng,
        }
    }

    #[allow(dead_code)]
    fn reset_except_visited(&mut self) {
        self.current_molecules_variants = HashMap::default();
        self.current_variant_molecules = HashMap::default();
        self.molecule_last_seen = HashMap::default();
        self.potential_variants_linked_reads = BTreeSet::new();
        self.potential_variant_counts_linked_reads = HashMap::default();
        self.potential_variants_hic = BTreeSet::new();
        self.potential_variant_counts_hic = HashMap::default();
        self.potential_variants_long_reads = BTreeSet::new();
        self.potential_variant_counts_long_reads = HashMap::default();
        self.iteration = 0;
    }


    //fn add_variant(&mut self, variants: &Vec<Vec<i32>>, molecules: &HashMap<i32, Vec<i32>>, variant: i32) {
    #[allow(dead_code)]
    fn add_variant(&mut self, variants: &Variants, molecules: &Molecules, variant: i32, available_vars: &HashSet<usize>) {
        let variant = variant.abs();
        self.variants_visited.set(variant as usize, true); // variants_visited is the i32.abs() not the index
        //println!("Adding variant {}", variant);
        let count = self.potential_variant_counts_linked_reads.entry(variant).or_insert(0);
        self.potential_variants_linked_reads.remove(&(*count, variant));
        self.potential_variant_counts_linked_reads.remove(&variant);
        let count = self.potential_variant_counts_long_reads.entry(variant).or_insert(0);
        self.potential_variants_long_reads.remove(&(*count, variant));
        self.potential_variant_counts_long_reads.remove(&variant);
        let count = self.potential_variant_counts_hic.entry(variant).or_insert(0);
        self.potential_variants_hic.remove(&(*count, variant));
        self.potential_variant_counts_hic.remove(&variant);
        let mut mols: HashSet<i32> = HashSet::default();
        self.iteration += 1;
        for mol in variants.get_linked_read_molecules(variant) { //&variants[(variant - 1) as usize] {
            let mol = mol.abs();
            mols.insert(mol);
            self.molecule_last_seen.insert(mol, self.iteration);
            // for each new molecule we touch, go through their variants and add or update them to potential_variants and potential_variant_scores
            if !self.current_molecules_variants.contains_key(&mol) {
                //println!("adding mol {} to current molecule set", mol);
                for var in molecules.get_linked_read_variants(mol) {
                    let var = var.abs();
                    assert!((var as usize) < self.variants_visited.len(), "{} < {}",var, self.variants_visited.len());
                    if !available_vars.contains(&(var as usize)) { continue; }
                    if self.variants_visited.get(var as usize).unwrap() { continue; }
                    let count = self.potential_variant_counts_linked_reads.entry(var).or_insert(0);
                    self.potential_variants_linked_reads.remove(&(*count, var));
                    //println!("{}", (1000.0/(variants[var as usize].len() as f32)).round() as u32);
                    *count += 1000/variants.get_linked_read_molecules(var).len() as u32;//(variants[(var - 1) as usize].len() as u32);
                    //*count += 1000; // HIC TEST
                    //println!("potential variant {} adding with count {}",var, count);
                    self.potential_variants_linked_reads.insert((*count, var));
                }
            }
            // add the variant to the current molecule sets
            let vars = self.current_molecules_variants.entry(mol).or_insert_with(|| HashSet::default());
            vars.insert(variant);
        }

        for mol in variants.get_hic_molecules(variant) {
            let mol = mol.abs();
            mols.insert(mol);
            if !self.current_molecules_variants.contains_key(&mol) {
                for var in molecules.get_hic_variants(mol) {
                    let var = var.abs();
                    if !available_vars.contains(&(var as usize)) { continue; }
                    if self.variants_visited.get(var as usize).unwrap() { continue; }
                    let count = self.potential_variant_counts_hic.entry(var).or_insert(0);
                    self.potential_variants_hic.remove(&(*count, var));
                    *count += 1;
                    self.potential_variants_hic.insert((*count, var));
                }
            }
            let vars = self.current_molecules_variants.entry(mol).or_insert_with(|| HashSet::default());
            vars.insert(variant);
        }

        for mol in variants.get_long_read_molecules(variant) {
            let mol = mol.abs();
            mols.insert(mol);
            if !self.current_molecules_variants.contains_key(&mol) {
                for var in molecules.get_long_read_variants(mol) {
                    let var = var.abs();
                    if !available_vars.contains(&(var as usize)) { continue; }
                    if self.variants_visited.get(var as usize).unwrap() { continue; }
                    let count = self.potential_variant_counts_long_reads.entry(var).or_insert(0);
                    self.potential_variants_long_reads.remove(&(*count, var));
                    *count += 1000/variants.get_long_read_molecules(var).len() as u32;
                    self.potential_variants_long_reads.insert((*count, var));
                }
            }
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

    #[allow(dead_code)]
    fn remove_molecule(&mut self, variants: &Variants, molecules: &Molecules, molecule: i32) {//variants: &Vec<Vec<i32>>, molecules: &HashMap<i32, Vec<i32>>, molecule: i32) {
        // so we will look at each variant in self.current_molecules_variants
        for var in self.current_molecules_variants.get(&molecule).unwrap().iter() {
            // and remove that molecule from self.current_variant_molecules
            self.current_variant_molecules.get_mut(var).unwrap().remove(&molecule); // what should I do when this goes to 0?
            if self.current_variant_molecules.get(var).unwrap().len() == 0 {
                self.current_variant_molecules.remove(var);
            }
        }
        // and we will go through the variants in self.molecules[molecule] and downgrade or remove them from self.potential_variants_counts and self.potential_variants
        for var in molecules.get_linked_read_variants(molecule) {//&molecules[&molecule] {
            let var = var.abs();
            if self.variants_visited.get(var as usize).unwrap() { continue; }
            let count = match self.potential_variant_counts_linked_reads.get(&var) {
                Some(x) => *x,
                None => 0,
            };
            if count > 0 {
                self.potential_variants_linked_reads.remove(&(count, var));
            }
            let incr: u32 =  1000/(variants.get_linked_read_molecules(var).len() as u32);//(variants[(var - 1) as usize].len() as u32);
            //let incr = 1000; // HIC TEST
            if count > incr {
                self.potential_variants_linked_reads.insert((count-incr, var));
                self.potential_variant_counts_linked_reads.insert(var, count-incr);
            } else {
                self.potential_variant_counts_linked_reads.remove(&var);
            }
        }
        // then we will remove the molecule from self.current_molecules_variants
        self.current_molecules_variants.remove(&molecule); 
        self.molecule_last_seen.remove(&molecule);
    }

    #[allow(dead_code)]
    fn get_next_variant(&mut self, variants: &Variants, min: f32, available_vars: &HashSet<usize>) -> Option<(i32, f32, u32, f32, u32, u32)> {
        //let min = 0.50; // 0.4 
        if self.potential_variant_counts_long_reads.len() == 0 && 
            self.potential_variant_counts_linked_reads.len() == 0 &&
            self.potential_variant_counts_hic.len() == 0 {
            //let var = self.rng.gen_range(1, variants.len() as i32);
            //return Some((var, 0.0, 0, 0.0, 0, 0))
            println!("potential variants 0, getting random variant");
            let vardex = self.rng.gen_range(0, available_vars.len());
            for (index, variant) in available_vars.iter().enumerate() {
                if index == vardex {
                    return Some((*variant as i32, 0.0, 0, 0.0, 0, 0))
                }
            }
        } else { 
            //loop {
                
                if self.potential_variants_long_reads.len() == 0 && 
                    self.potential_variants_hic.len() == 0 &&
                    self.potential_variants_linked_reads.len() == 0 { return None; }
                let mut long_read_percent = 0.0;
                let mut long_read_count = 0;
                let mut long_read_count_full = 0;
                let mut long_read_variant = 0;
                let mut linked_read_percent = 0.0;
                let mut linked_read_count = 0;
                let mut linked_read_count_full = 0;
                let mut linked_read_variant = 0;
                let mut hic_variant = 0;
                let mut hic_count = 0;
                if self.potential_variants_long_reads.len() > 0 {
                    let (count, variant) = self.potential_variants_long_reads.iter().next_back().unwrap();
                    let (count, variant) = (*count , *variant);
                    let percent = (count as f32)/1000.0;
                    let count2: u32 = (percent * (variants.get_long_read_molecules(variant).len() as f32)).round() as u32;//(variants[(variant-1) as usize].len() as f32)).round() as i32; 
                    //println!("long read check, raw count {} molecules that hit this variant {} {}", count, count2, variants.get_long_read_molecules(variant).len());
                    long_read_percent = percent;
                    long_read_count = count2;
                    long_read_count_full = count;
                    long_read_variant = variant;
                }
                if self.potential_variants_linked_reads.len() > 0 {
                    let (count, variant) = self.potential_variants_linked_reads.iter().next_back().unwrap();
                    let (count, variant) = (*count, *variant);
                    let percent = (count as f32)/1000.0;
                    let count2 : u32 = (percent * (variants.get_linked_read_molecules(variant).len() as f32)).round() as u32;
                    //println!("best linked read variant {} {}",count2, percent);
                    linked_read_percent = percent;
                    linked_read_count = count2;
                    linked_read_count_full = count;
                    linked_read_variant = variant;
                }
                if self.potential_variants_hic.len() > 0 {
                    let (count, variant) = self.potential_variants_hic.iter().next_back().unwrap();
                    let (count, variant) = (*count, *variant);
                    hic_count = count;
                    hic_variant = variant;
                }
                println!("{} {} {} {} {}", long_read_percent, long_read_count, linked_read_percent, linked_read_count, hic_count);
                    //println!("{} {} {}", variant, percent, count2);
                    //let count = ((count as f32)/1000.0).round() as i32; // HIC TEST
                    //let percent = (count as f32)/(variants[(variant-1) as usize].len() as f32);

                    if long_read_percent >= min && long_read_count > 4 {
                        //println!("PB PB PB PB");
                        self.potential_variants_long_reads.remove(&(long_read_count_full, long_read_variant));
                        self.potential_variant_counts_long_reads.remove(&long_read_variant);
                        if let Some(count) = self.potential_variant_counts_linked_reads.get(&long_read_variant) {
                            self.potential_variants_linked_reads.remove(&(*count, long_read_variant));
                            self.potential_variant_counts_linked_reads.remove(&long_read_variant);
                        }
                        if let Some(count) = self.potential_variant_counts_hic.get(&long_read_variant) {
                            self.potential_variants_hic.remove(&(*count, long_read_variant));
                            self.potential_variant_counts_hic.remove(&long_read_variant);
                        }
                        assert!(available_vars.contains(&(long_read_variant as usize)));
                        return Some((long_read_variant, long_read_percent, long_read_count, linked_read_percent, linked_read_count, hic_count));
                    } else if linked_read_percent >= min {//&& linked_read_count > 10 {
                        //println!("10X 10X 10X 10X 10X 10X");
                        self.potential_variants_linked_reads.remove(&(linked_read_count_full, linked_read_variant));
                        self.potential_variant_counts_linked_reads.remove(&linked_read_variant);
                        if let Some(count) = self.potential_variant_counts_long_reads.get(&long_read_variant) {
                            self.potential_variants_long_reads.remove(&(*count, long_read_variant));
                            self.potential_variant_counts_long_reads.remove(&long_read_variant);
                        }
                        if let Some(count) = self.potential_variant_counts_hic.get(&long_read_variant) {
                            self.potential_variants_hic.remove(&(*count, long_read_variant));
                            self.potential_variant_counts_hic.remove(&long_read_variant);
                        }
                        assert!(available_vars.contains(&(linked_read_variant as usize)));
                        return Some((linked_read_variant, long_read_percent, long_read_count, linked_read_percent, linked_read_count, hic_count));
                    } else if hic_count > 50 {
                        println!("HIIIIIIIIIIIIIIIIC");
                        self.potential_variants_hic.remove(&(hic_count, hic_variant));
                        self.potential_variant_counts_hic.remove(&hic_variant);
                        if let Some(count) = self.potential_variant_counts_long_reads.get(&long_read_variant) {
                            self.potential_variants_long_reads.remove(&(*count, long_read_variant));
                            self.potential_variant_counts_long_reads.remove(&long_read_variant);
                        }
                        if let Some(count) = self.potential_variant_counts_linked_reads.get(&long_read_variant) {
                            self.potential_variants_linked_reads.remove(&(*count, long_read_variant));
                            self.potential_variant_counts_linked_reads.remove(&long_read_variant);
                        }
                        assert!(available_vars.contains(&(hic_variant as usize)));
                        return Some((hic_variant, long_read_percent, long_read_count, linked_read_percent, linked_read_count, hic_count));
                    } else {
                        return None
                    }
                } 
               
                    //println!("rejecting {} because {} < {} or {} <= 10.0", variant, count, min, (count as f32)/1000.0*(variants[(variant-1) as usize].len() as f32));
                
            //}
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

#[allow(dead_code)]
struct Molecules {
    //long_read_molecule_list: Vec<i32>,
    longread_molid_offsets: Vec<i32>,
    linked_read_molecules: HashMap<i32, Vec<i32>>,
    hic_molecules: HashMap<i32, Vec<i32>>,
    long_read_molecules: HashMap<i32, Vec<i32>>,
    long_read_molecule_list: HashMap<i32, Vec<i32>>,
    long_read_het_molecules: HashMap<i32, Vec<i32>>,
    long_read_het_molecule_list: HashMap<i32, Vec<i32>>,
    long_read_molecule_positions: HashMap<i32, Vec<i32>>,
    hom_linked_read_molecules: HashMap<i32, HashSet<i32>>,
    hom_long_read_molecules: HashMap<i32, HashSet<i32>>,
    hom_hic_molecules: HashMap<i32, HashSet<i32>>,
    het_linked_read_molecules: HashMap<i32, HashSet<i32>>,
    het_long_read_molecules: HashMap<i32, HashSet<i32>>,
    het_hic_molecules: HashMap<i32, HashSet<i32>>,
}

impl Molecules {
    fn get_linked_read_variants<'a>(&'a self, mol: i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.linked_read_molecules.get(&mol.abs()) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    //fn get_linked_read_variants(&self, mol: i32) -> std::slice::Iter<'_, i32> {
    //    return self.linked_read_molecules.get(&mol.abs()).unwrap().iter();
    //}
    fn get_long_read_molecules<'a>(&'a self) -> Box<dyn Iterator<Item=&i32>+'a> {
        Box::new(self.long_read_molecules.keys())
    }
    fn get_linked_read_molecules<'a>(&'a self) -> Box<dyn Iterator<Item=&i32>+'a> {
        Box::new(self.linked_read_molecules.keys())
    }

    fn get_long_read_variants<'a>(&'a self, mol: i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.long_read_molecules.get(&mol.abs()) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }

    fn get_long_read_variants_ordered<'a>(&'a self, mol: i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.long_read_molecule_list.get(&mol.abs()) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    fn get_long_read_variants_and_positions<'a>(&'a self, mol: i32) -> Box<dyn Iterator<Item=(&i32, &i32)>+'a> {
        match self.long_read_molecule_list.get(&mol.abs()) {

            Some(x) => {
                match self.long_read_molecule_positions.get(&mol.abs()) {
                    Some(y) => Box::new(x.iter().zip(y.iter())),
                    None => Box::new(std::iter::empty()),
                }
            },
            None => Box::new(std::iter::empty()),
        }
    }

    fn get_hic_variants<'a>(&'a self, mol: i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hic_molecules.get(&mol.abs()) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }

    fn get_hom_linked_read_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hom_linked_read_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    
    fn get_hom_long_read_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hom_long_read_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }

    fn get_hom_hic_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hom_hic_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }

    fn get_het_linked_read_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.het_linked_read_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    fn get_het_long_read_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.het_long_read_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }

    fn get_het_hic_variants<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.het_long_read_molecules.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }

    fn get_variants<'a>(&'a self, i: &i32, kmer_type: KmerType) -> Box<dyn Iterator<Item=&i32>+'a> {
        match kmer_type {
            KmerType::PairedHet => Box::new(
                    self.get_linked_read_variants(*i).chain(
                    self.get_long_read_variants(*i)).chain(
                    self.get_hic_variants(*i))),
            KmerType::UnpairedHet => Box::new(
                    self.get_het_linked_read_variants(i).chain(
                    self.get_het_long_read_variants(i)).chain(
                    self.get_het_hic_variants(i))),
            KmerType::Homozygous => Box::new(
                    self.get_hom_linked_read_variants(i).chain(
                    self.get_hom_long_read_variants(i)).chain(
                    self.get_hom_hic_variants(i))),
        }
    }
}

impl Variants {
    fn get_variant_iter<'a>(&'a self, kmer_type: KmerType) -> Box<dyn Iterator<Item=&i32>+'a> {
        match kmer_type {
            KmerType::PairedHet => Box::new(self.paired_het_variants.iter()),
            KmerType::UnpairedHet => Box::new(std::iter::empty()),
            KmerType::Homozygous => Box::new(std::iter::empty()),
        }
    }

    fn len(&self) -> usize {
        return self.linked_read_variants.len();
    }

    fn get_linked_read_molecules(&self, i: i32) -> std::slice::Iter<'_, i32> {
        return self.linked_read_variants[(i.abs() - 1) as usize].iter();
    }
    fn get_num_linked_read_molecules(&self, i: i32) -> usize {
        return self.linked_read_variants[(i.abs() - 1) as usize].len();
    }

    fn get_long_read_molecules(&self, i: i32) -> std::slice::Iter<'_, i32> {
        return self.long_read_variants[(i.abs() - 1) as usize].iter();
    }
    fn get_num_long_read_molecules(&self, i: i32) -> usize {
        return self.long_read_variants[(i.abs() - 1) as usize].len();
    }

    fn get_hic_molecules(&self, i: i32) -> std::slice::Iter<'_, i32> {
        return self.hic_variants[(i.abs() - 1) as usize].iter();
    }
    fn get_num_hic_molecules(&self, i: i32) -> usize {
        return self.hic_variants[(i.abs() - 1) as usize].len();
    }
    #[allow(dead_code)]
    fn get_hom_linked_read_molecules<'a>(&'a self, i: i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hom_linked_read_variants.get(&i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    #[allow(dead_code)]
    fn get_hom_long_read_molecules<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.hom_long_read_variants.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    #[allow(dead_code)]
    fn get_het_linked_read_molecules<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.het_linked_read_variants.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }
    #[allow(dead_code)]
    fn get_het_long_read_molecules<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
        match self.het_long_read_variants.get(i) {
            Some(x) => Box::new(x.iter()),
            None => Box::new(std::iter::empty()),
        }
    }

    //fn get_het_hic_molecules<'a>(&'a self, i: &i32) -> Box<dyn Iterator<Item=&i32>+'a> {
    //    match self.het_hic_variants.get(i) {
    //        Some(x) => Box::new(x.iter()),
    //        None => Box::new(std::iter::empty()),
    //    }
    //}

    fn get_molecules<'a>(&'a self, i: &i32, data_type: DataType) 
            -> Box<dyn Iterator<Item=&i32>+'a> {
        match data_type {
            DataType::HC => Box::new(self.get_hic_molecules(*i)),
            DataType::PB => Box::new(self.get_long_read_molecules(*i)),
            DataType::LR => Box::new(self.get_linked_read_molecules(*i)),
            DataType::PBLR => Box::new(self.get_long_read_molecules(*i).chain(self.get_linked_read_molecules(*i))),
        }
    }
    fn get_num_molecules(&self, i: &i32, data_type: DataType) -> usize {
        match data_type {
            DataType::HC => self.get_num_hic_molecules(*i),
            DataType::PB => self.get_num_long_read_molecules(*i),
            DataType::LR => self.get_num_linked_read_molecules(*i),
            DataType::PBLR => self.get_num_long_read_molecules(*i) + self.get_num_linked_read_molecules(*i),
        }
    }
}

#[allow(dead_code)]
enum DataType {
    HC, // hic
    PB, // pacbio
    LR, // linked reads
    PBLR, //pb and lr
}

struct Crib {
    variants: HashMap<i32, (i32, usize, usize)>, // map from kmer_id to crib molecule id and number seen and position
    molecules: HashMap<i32, HashMap<i32, (usize, usize)>>, // map from crib molecule id to a map from kmer_id to order index
}

fn load_crib(crib: &String, kmers: &Kmers) -> Crib {
    let mut mol_id = 1;
    let mut variants: HashMap<i32, (i32, usize, usize)> = HashMap::new();
    let mut molecules: HashMap<i32, HashMap<i32, (usize, usize)>> = HashMap::new();


    let mut paired_het_variants: HashSet<i32> = HashSet::new();
    let mut long_read_variants: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_molecule_list: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut hom_long_read_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut hom_long_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut bufi32 = [0u8; 4];
    let mut buf2 = [0u8; 4];
    let f = File::open(crib)
        .expect(&format!("Unable to open crib {}", crib));
    let mut reader = BufReader::new(f);
    let mut ordered_kmers: Vec<(i32, i32, usize)> = Vec::new();
    'outer: loop { // ok here we go again. Pacbio/longread data. Format is i32s until you hit a zero. when you hit two zeros you are done
        let mut vars: HashSet<i32> = HashSet::new();
        let mut varlist: Vec<i32> = Vec::new();
        let mut hom_kmers: HashSet<i32> = HashSet::new();
        //let mut het_kmers: HashSet<i32> = HashSet::new();
        let mut order = 0;
        loop {
            if let Some(kmer_id) = eat_i32(&mut reader, &mut bufi32) {
                if kmer_id == 0 { break; }
                if let Some(position) = eat_i32(&mut reader, &mut buf2){
                    let position = position + 10;
                    //eprintln!("position {}", position);
                    let mut molecule = 0;
                    let mut number = 0;
                    let mut position1 = 0;
                    let mut has = false;
                   

                    
                    let var_order = molecules.entry(mol_id).or_insert(HashMap::new());
                    var_order.insert(kmer_id.abs(), (order, position as usize));
                    match kmers.kmer_type.get(&kmer_id) {
                        Some(KmerType::PairedHet) => {
                            ordered_kmers.push((mol_id, kmer_id.abs(), position as usize));
                            //eprintln!("paired het");
                            vars.insert(kmer_id);
                            paired_het_variants.insert(kmer_id.abs());
                            varlist.push(kmer_id);
                             if let Some((mol, num, pos)) = variants.get(&kmer_id.abs()) {
                                molecule = *mol;
                                number = *num;
                                position1 = *pos;
                                has = true;
                            } else {
                                variants.insert(kmer_id.abs(), (mol_id, 1, position as usize));
                            }
                            if has {
                                variants.insert(kmer_id.abs(), (mol_id, number+1, position as usize));
                            }
                        },
                        Some(KmerType::UnpairedHet) => {
                            //het_kmers.insert(kmer_id);
                        },
                        Some(KmerType::Homozygous) => {
                            hom_kmers.insert(kmer_id);
                        },
                        None => { eprintln!("no kmer type? {}", kmer_id); }
                    }
                }
                
                
            } else { break 'outer; }
            order += 1;
        }
            for kmer_id in vars.iter() {
                let var_bcs = long_read_variants.entry(kmer_id.abs()).or_insert(HashSet::new());
                if kmer_id < &0 { var_bcs.insert(-mol_id); } else { var_bcs.insert(mol_id); }
            }
            long_read_molecules.insert(mol_id, vars); 
            long_read_molecule_list.insert(mol_id, varlist);
            
        if hom_kmers.len() > 0 {
            for kmer_id in hom_kmers.iter() {
                let var_bcs = hom_long_read_kmers.entry(*kmer_id).or_insert(HashSet::new());
                var_bcs.insert(mol_id);
            }
        }



        hom_long_read_molecules.insert(mol_id, hom_kmers);
        mol_id += 1;
        
    }
    println!("{}",ordered_kmers.len());
    println!("chrom\tpos\tkmer\thits\tkmer_count\thas_pair\tpair_chrom\tpair_pos\tpair_id\tpair_hits\tpair_count");
    for (mol, kmer, position) in ordered_kmers {
        let (_, num, _) = variants.get(&kmer).unwrap();
        if let Some((mol2, num2, position2)) = variants.get(&pair(kmer.abs())) {
            println!("{}\t{}\t{}\t{}\t{}\ttrue\t{}\t{}\t{}\t{}\t{}", mol, position, kmers.kmers.get(&kmer).unwrap(), num, kmers.kmer_counts.get(&kmer).unwrap(), 
                mol2, position2, kmers.kmers.get(&pair(kmer)).unwrap(), num2, kmers.kmer_counts.get(&pair(kmer.abs())).unwrap());
        } else {
            println!("{}\t{}\t{}\t{}\t{}\tfalse\tnone\tnone\tnone\tnone\tnone", mol, position, kmers.kmers.get(&kmer).unwrap(), num, kmers.kmer_counts.get(&kmer).unwrap());
        }
        
    } 

    Crib{
        variants: variants,
        molecules: molecules,
    }
}

#[allow(dead_code)]
struct Variants {
    paired_het_variants: HashSet<i32>,
    linked_read_variants: Vec<Vec<i32>>,
    hic_variants: Vec<Vec<i32>>,
    long_read_variants: Vec<Vec<i32>>,
    long_read_het_variants: Vec<Vec<i32>>,
    hom_linked_read_variants: HashMap<i32, HashSet<i32>>,
    hom_long_read_variants: HashMap<i32, HashSet<i32>>,
    hom_hic_variants: HashMap<i32, HashSet<i32>>,
    het_linked_read_variants: HashMap<i32, HashSet<i32>>,
    het_long_read_variants: HashMap<i32, HashSet<i32>>,
    het_hic_variants: HashMap<i32, HashSet<i32>>,
}

fn load_molecule_kmers(params: &Params, kmers: &Kmers) -> (Variants, Molecules){
    let mut linked_read_variants: HashMap<i32, HashSet<i32>> = HashMap::new(); // map from variant_id to list of molecule_ids
    let _hic_variants: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_variants: HashMap<i32, HashSet<i32>> = HashMap::new();
    
    let mut linked_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new(); //map from molecule_id to list of variant_ids
    let hic_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_molecule_list: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut long_read_molecule_positions: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut long_read_het_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut long_read_het_molecule_list: HashMap<i32, Vec<i32>> = HashMap::new();

    let mut paired_het_variants: HashSet<i32> = HashSet::new();

    let mut hom_linked_read_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut hom_long_read_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut hom_linked_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let mut hom_long_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let hom_hic_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let hom_hic_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();

    let het_linked_read_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let het_long_read_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let het_linked_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let het_long_read_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();
    let het_hic_kmers: HashMap<i32, HashSet<i32>> = HashMap::new();
    let het_hic_molecules: HashMap<i32, HashSet<i32>> = HashMap::new();

    let mut bufi32 = [0u8; 4];
    let mut buf2 = [0u8; 4];
    let mut max_var = 0;
    let mut max_molid = 0;
    for txg_file in params.txg_mols.iter() {
        println!("{}",txg_file);
        let f = File::open(txg_file.to_string())
            .expect(&format!("Unable to open txg file {}", txg_file));
        let mut reader = BufReader::new(f);
        loop {
            if let Some(barcode_id) = eat_i32(&mut reader, &mut bufi32) {
                if let Some(kmer_id) = eat_i32(&mut reader, &mut bufi32) {
                    //if barcode_id == 0 && kmer_id == 0 { break; } // only dealing with 10x data so far
                    if kmer_id.abs() > max_var { max_var = kmer_id.abs() }
                    if barcode_id > max_molid { max_molid = barcode_id; }
                    match kmers.kmer_type.get(&kmer_id).unwrap() {
                        KmerType::PairedHet => {
                            //eprintln!("paired het kmer");
                            let bc_vars = linked_read_molecules.entry(barcode_id).or_insert(HashSet::new());
                            bc_vars.insert(kmer_id);
                            paired_het_variants.insert(kmer_id.abs());
                        },
                        KmerType::UnpairedHet => {
                            //eprintln!("unpaired het kmer");
                            //let bc_vars = het_linked_read_molecules.entry(barcode_id).or_insert(HashSet::new());
                            //bc_vars.insert(kmer_id);
                            ()
                        },
                        KmerType::Homozygous => {
                            //eprintln!("homozygous kmer adding to hom_lingked_read_molecules");
                            let bc_vars = hom_linked_read_molecules.entry(barcode_id).or_insert(HashSet::new());
                            bc_vars.insert(kmer_id);
                            ()
                        },
                    }
                    
                    //println!("bc {} kmer id {}",barcode_id, kmer_id);
                    
                    
                } else { break; }
            } else { break; }
        }
    }
    
    eprintln!("reducing to good linked read molecules with > 5 variants, right now we have {}", linked_read_molecules.len());
    //linked_read_molecules.retain(|_key, value| value.len() > 10 && value.len() < 5000); 
    for (mol, varset) in linked_read_molecules.iter() {
        for var in varset.iter() {
            let var_bcs = linked_read_variants.entry(var.abs()).or_insert(HashSet::new());
            if *var < 0 { var_bcs.insert(-mol); } else { var_bcs.insert(*mol); }
        }
    }
    for (var, mols) in linked_read_variants.iter() {
        if mols.len() < 30 {
            for mol in mols.iter() {
                let vars = linked_read_molecules.get_mut(&mol.abs()).unwrap();
                vars.remove(&var);
                vars.remove(&-var);
            }
        }
    }
    eprintln!("{} good molecules", linked_read_molecules.len());
    /*
    for (mol, varset) in het_linked_read_molecules.iter() {
        for var in varset.iter() {
            let var_bcs = het_linked_read_kmers.entry(*var).or_insert(HashSet::new());
            var_bcs.insert(*mol);
        }
    }
    */
    for (mol, varset) in hom_linked_read_molecules.iter() {
        //eprintln!("mol {} with varset length {}", mol, varset.len());
        for var in varset.iter() {
            let var_bcs = hom_linked_read_kmers.entry(*var).or_insert(HashSet::new());
            var_bcs.insert(*mol);
        }
    }
    
    eprintln!("{} hom linked read molecules, {} hom linked read kmers",hom_linked_read_molecules.len(), hom_linked_read_kmers.len());
    //molecules.clear();

    /*
    let mut bad_vars = 0;
    let mut vars_to_remove: Vec<i32> = Vec::new();
    let mut bad_var_set: HashSet<i32> = HashSet::new();
    for (var, var_mols) in linked_read_variants.iter() {
        let mut mols: HashSet<i32> = HashSet::new();
        for mol in var_mols.iter() {
            mols.insert(mol.abs());
        }
        if mols.len() + 1 < var_mols.len() {
            bad_vars += 1;
            //println!("bad variant {}, {} mols with both versions, {} total mols with either.", kmer_id_to_kmer.get(var).unwrap(), var_mols.len() - mols.len(), mols.len());
            for mol in var_mols.iter() {
                if let Some(linked_read_molecule) = linked_read_molecules.get_mut(&mol.abs()) {
                    linked_read_molecules.remove(var);
                    linked_read_molecules.remove(&-var);
                }
            }
            vars_to_remove.push(var.abs());
            bad_var_set.insert(var.abs());
        }
    }
    */
    let mut mol_id = max_molid + 1;
    for hic_file in params.hic_mols.iter() {
        let f = File::open(hic_file.to_string())
            .expect(&format!("Unable to open hic file {}", hic_file));
        let mut reader = BufReader::new(f);
        'outerhic: loop { // now deal with hic data, format is i32s until you hit a 0 if i get to 2 0's we are done
            //break 'outerhic;
            let mut vars: HashSet<i32> = HashSet::new();
            let mut vars2: Vec<i32> = Vec::new();
            loop {
                if let Some(kmer_id) = eat_i32(&mut reader, &mut bufi32) {
                    //println!("kmer id {}", kmer_id);

                    if kmer_id == 0 { if vars.len() == 0 { break 'outerhic; } else { break; } }
                    match kmers.kmer_type.get(&kmer_id).unwrap() {
                        KmerType::PairedHet => {
                            //eprintln!("paired het kmer");
                                vars.insert(kmer_id); 
                                vars2.push(kmer_id);
                            ()
                        },
                        KmerType::UnpairedHet => {
                            //eprintln!("unpaired het kmer");
                            vars.insert(kmer_id); 
                            vars2.push(kmer_id);
                            //let bc_vars = het_linked_read_molecules.entry(barcode_id).or_insert(HashSet::new());
                            //bc_vars.insert(kmer_id);
                            ()
                        },
                        KmerType::Homozygous => {
                            //eprintln!("homozygous kmer");
                            vars.insert(kmer_id); 
                            vars2.push(kmer_id);
                            //let bc_vars = hom_linked_read_molecules.entry(barcode_id).or_insert(HashSet::new());
                            //bc_vars.insert(kmer_id);
                            ()
                        },
                    }
                    //if !bad_var_set.contains(&kmer_id.abs()) {
                        
                    
                    if kmer_id.abs() > max_var { max_var = kmer_id.abs() }
                    //}
                } else { break 'outerhic; }
            }
            mol_id += 1;
        }
    }
    
        /*
        if vars.len() > 1 {
            hic_molecules.insert(mol_id, vars);
            for kmer_id in vars2 {
                let var_bcs = hic_variants.entry(kmer_id.abs()).or_insert(HashSet::new());
                if kmer_id < 0 { var_bcs.insert(-mol_id); } else { var_bcs.insert(mol_id); }
            }
        }
        */
 
    eprintln!("num hic molecules is {}", hic_molecules.len());
    
    mol_id += 1; //
    let mut longread_mol_id_starts: Vec<i32> = Vec::new();
    for longread_file in params.longread_mols.iter() {
        longread_mol_id_starts.push(mol_id);
        let f = File::open(longread_file.to_string())
            .expect(&format!("Unable to open longread file {}", longread_file));
        let mut reader = BufReader::new(f);
        'outer: loop { // ok here we go again. Pacbio/longread data. Format is i32s until you hit a zero. when you hit two zeros you are done
            let mut vars: HashSet<i32> = HashSet::new();
            let mut varlist: Vec<i32> = Vec::new();
            let mut varlist_positions: Vec<i32> = Vec::new();
            let mut hom_kmers: HashSet<i32> = HashSet::new();
            let mut het_kmers: HashSet<i32> = HashSet::new();
            let mut het_kmers_list: Vec<i32> = Vec::new();
            loop {
                if let Some(kmer_id) = eat_i32(&mut reader, &mut bufi32) {
                    
                    if kmer_id == 0 { break; }
                    if let Some(position) = eat_i32(&mut reader, &mut buf2) {
                        let position = position + 10;
                        if kmer_id.abs() > max_var { max_var = kmer_id.abs(); }
                        match kmers.kmer_type.get(&kmer_id) {
                            Some(KmerType::PairedHet) => {
                                vars.insert(kmer_id);
                                paired_het_variants.insert(kmer_id.abs());
                                varlist.push(kmer_id);
                                varlist_positions.push(position);
                            },
                            Some(KmerType::UnpairedHet) => {
                                het_kmers.insert(kmer_id);
                                het_kmers_list.push(kmer_id);
                            },
                            Some(KmerType::Homozygous) => {
                                hom_kmers.insert(kmer_id);
                            },
                            None => { eprintln!("no kmer type? {}", kmer_id); }
                        }
                    }
                    //eat_i32(&mut reader, &mut bufi32); //TODO get rid of this old format
                    
                } else { break 'outer; }
            }
            //if vars.len() > 0 { 
                for kmer_id in vars.iter() {
                    let var_bcs = long_read_variants.entry(kmer_id.abs()).or_insert(HashSet::new());
                    if kmer_id < &0 { var_bcs.insert(-mol_id); } else { var_bcs.insert(mol_id); }
                }
                long_read_molecules.insert(mol_id, vars); 
                long_read_molecule_list.insert(mol_id, varlist);
                long_read_molecule_positions.insert(mol_id, varlist_positions);
                long_read_het_molecules.insert(mol_id, het_kmers);
                long_read_het_molecule_list.insert(mol_id, het_kmers_list);
                
                
            //}
            /*
            if het_kmers.len() > 0 {
                for kmer_id in het_kmers.iter() {
                    let var_bcs = het_long_read_kmers.entry(*kmer_id).or_insert(HashSet::new());
                    var_bcs.insert(mol_id);
                }
            }
            */
            if hom_kmers.len() > 0 {
                for kmer_id in hom_kmers.iter() {
                    let var_bcs = hom_long_read_kmers.entry(*kmer_id).or_insert(HashSet::new());
                    var_bcs.insert(mol_id);
                }
            }
            hom_long_read_molecules.insert(mol_id, hom_kmers);
            mol_id += 1;
            
        }
    }
    
    //long_read_molecules.retain(|_key, value| value.len() > 1 && value.len() < 5000); 
    eprintln!("{} good long read molecules", long_read_molecules.len());
    
    eprintln!("{} hom long read kmers, {} hom read long molecules", hom_long_read_kmers.len(), hom_long_read_molecules.len());

    
    //for var in &vars_to_remove {
    //    linked_read_variants.insert(*var, HashSet::new());
    //    hic_variants.insert(*var, HashSet::new());
    //    long_read_variants.insert(*var, HashSet::new());
    //}
    //println!("found {} bad vars", bad_vars);

    let mut txg_vars: Vec<Vec<i32>> = Vec::new();
    let mut hic_vars: Vec<Vec<i32>> = Vec::new();
    let mut pb_vars: Vec<Vec<i32>> = Vec::new();
    let mut pb_het_vars: Vec<Vec<i32>> = Vec::new();
    for _i in 0..(max_var+1) {
        txg_vars.push(Vec::new());
        hic_vars.push(Vec::new());
        pb_vars.push(Vec::new());
        pb_het_vars.push(Vec::new());
    }

    let mut txg_mols: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut hic_mols: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut pb_mols: HashMap<i32, Vec<i32>> = HashMap::new();
    let mut pb_het_mols: HashMap<i32, Vec<i32>> = HashMap::new();

    for (mol_id, kmer_ids) in linked_read_molecules {
        let mut kmer_ids_sorted = kmer_ids.iter().cloned().collect::<Vec<i32>>();
        kmer_ids_sorted.sort_by(| a, b | a.abs().cmp(&b.abs()));
        for var in kmer_ids_sorted.iter() {
            if var < &0 {
                txg_vars[(var.abs() - 1) as usize].push(-mol_id);
            } else { txg_vars[(var - 1) as usize].push(mol_id) }
        }
        txg_mols.insert(mol_id, kmer_ids_sorted);
    }

    for (mol_id, kmer_ids) in hic_molecules {
        let mut kmer_ids_sorted = kmer_ids.iter().cloned().collect::<Vec<i32>>();
        kmer_ids_sorted.sort_by(| a, b | a.abs().cmp(&b.abs()));
        for var in kmer_ids_sorted.iter() {
            if var < &0 {
                hic_vars[(var.abs() - 1) as usize].push(-mol_id);
            } else { hic_vars[(var - 1) as usize].push(mol_id) }
        }
        hic_mols.insert(mol_id, kmer_ids_sorted);
    }
    for (mol_id, kmer_ids) in long_read_molecules {
        let mut kmer_ids_sorted = kmer_ids.iter().cloned().collect::<Vec<i32>>();
        kmer_ids_sorted.sort_by(| a, b | a.abs().cmp(&b.abs()));
        for var in kmer_ids_sorted.iter() {
            if var < &0 {
                pb_vars[(var.abs() - 1) as usize].push(-mol_id);
            } else { pb_vars[(var - 1) as usize].push(mol_id) }

        }
        pb_mols.insert(mol_id, kmer_ids_sorted);
    }
    for (mol_id, kmer_ids) in long_read_het_molecules {
        let mut kmer_ids_sorted = kmer_ids.iter().cloned().collect::<Vec<i32>>();
        kmer_ids_sorted.sort_by(| a, b | a.abs().cmp(&b.abs()));
        for var in kmer_ids_sorted.iter() {
            if var < &0 {
                pb_het_vars[(var.abs() - 1) as usize].push(-mol_id);
            } else { pb_het_vars[(var - 1) as usize].push(mol_id) }

        }
        pb_het_mols.insert(mol_id, kmer_ids_sorted);
    }
    for i in 0..max_var {
        txg_vars[i as usize].sort_by(| a, b | a.abs().cmp(&b.abs()));
        hic_vars[i as usize].sort_by(| a, b | a.abs().cmp(&b.abs()));
        pb_vars[i as usize].sort_by(| a, b | a.abs().cmp(&b.abs()));
        pb_het_vars[i as usize].sort_by(| a, b | a.abs().cmp(&b.abs()));
    }
    let vars = Variants{
        paired_het_variants: paired_het_variants,
        linked_read_variants: txg_vars,
        hic_variants: hic_vars,
        long_read_variants: pb_vars,
        long_read_het_variants: pb_het_vars,
        hom_linked_read_variants: hom_linked_read_kmers,
        hom_long_read_variants: hom_long_read_kmers,
        hom_hic_variants: hom_hic_kmers,
        het_linked_read_variants: het_linked_read_kmers,
        het_long_read_variants: het_long_read_kmers,
        het_hic_variants: het_hic_kmers,
    };
    let mols = Molecules{
        longread_molid_offsets: longread_mol_id_starts,
        linked_read_molecules: txg_mols,
        hic_molecules: hic_mols,
        long_read_molecules: pb_mols,
        long_read_molecule_list: long_read_molecule_list,
        long_read_molecule_positions: long_read_molecule_positions,
        long_read_het_molecules: pb_het_mols,
        long_read_het_molecule_list: long_read_het_molecule_list,
        hom_linked_read_molecules: hom_linked_read_molecules,
        hom_long_read_molecules: hom_long_read_molecules,
        hom_hic_molecules: hom_hic_molecules,
        het_linked_read_molecules: het_linked_read_molecules,
        het_long_read_molecules: het_long_read_molecules,
        het_hic_molecules: het_hic_molecules,
    };
    (vars, mols)
}

#[derive(PartialEq)]
enum KmerType {
    PairedHet,
    UnpairedHet,
    Homozygous,
}

struct Kmers {
    kmers: HashMap<i32, String>,
    kmer_counts: HashMap<i32, i32>,
    kmer_type: HashMap<i32, KmerType>,
}

fn load_kmers(kmerfile: &String) -> Kmers {
    let mut kmers: HashMap<i32, String> = HashMap::new();
    let mut kmer_id: i32 = 1; // no 0 kmer id as we are using sign for which of the pair
    let mut kmer_type: HashMap<i32, KmerType> = HashMap::new();
    let mut kmer_counts: HashMap<i32, i32> = HashMap::new();
    let reader = File::open(kmerfile).expect("cannot open kmer file");
    let mut reader = BufReader::new(reader);
    let mut buf1 = vec![]; let mut buf2 = vec![];
    let mut buf3 = vec![]; let mut buf4 = vec![];
    loop {
        buf1.clear(); buf2.clear(); buf3.clear(); buf4.clear();

        let bytes1 = reader.read_until(b'\t', &mut buf1).expect("cannot read file");
        if bytes1 == 0 { break; } 
   
        let bytes2 = reader.read_until(b'\t', &mut buf2).expect("cannot read file");
        if bytes2 == 0 { break; } 

        let bytes3 = reader.read_until(b'\t', &mut buf3).expect("cannot read file");
        if bytes3 == 0 { break; } 

        let bytes4 = reader.read_until(b'\n', &mut buf4).expect("cannot read file");
        if bytes4 == 0 { break; } 

        let switchme = std::str::from_utf8(&buf3).unwrap();
        match switchme {
            "HET\t" => {
                kmer_type.insert(kmer_id, KmerType::UnpairedHet);
                kmer_type.insert(-kmer_id, KmerType::UnpairedHet);
                kmers.insert(kmer_id, std::str::from_utf8(&buf1[0..(bytes1-1)]).unwrap().to_string());
                //println!("")
            },
            "HOM\t" => {
                kmer_type.insert(kmer_id, KmerType::Homozygous);
                kmer_type.insert(-kmer_id, KmerType::Homozygous);
            },
            _ => {
                kmer_type.insert(kmer_id, KmerType::PairedHet);
                kmer_type.insert(-kmer_id, KmerType::PairedHet);
                kmer_type.insert(kmer_id+1, KmerType::PairedHet);
                kmer_type.insert(-(kmer_id+1), KmerType::PairedHet);
                kmers.insert(kmer_id, std::str::from_utf8(&buf1[0..(bytes1-1)]).unwrap().to_string());
                kmers.insert(-kmer_id, std::str::from_utf8(&buf1[0..(bytes1-1)]).unwrap().to_string());
                let count = std::str::from_utf8(&buf2[0..(bytes2-1)]).unwrap().to_string().parse::<i32>().unwrap();
                kmer_counts.insert(kmer_id, count);
                kmer_counts.insert(-kmer_id, count);
                kmers.insert(kmer_id+1, std::str::from_utf8(&buf3[0..(bytes3-1)]).unwrap().to_string());
                kmers.insert(-(kmer_id+1), std::str::from_utf8(&buf3[0..(bytes3-1)]).unwrap().to_string());
                let count = std::str::from_utf8(&buf4[0..(bytes4-1)]).unwrap().to_string().parse::<i32>().unwrap();
                kmer_counts.insert(kmer_id+1, count);
                kmer_counts.insert(-(kmer_id+1), count);
            },
        }
        kmer_id += 2;
    }
    //(kmers, kmer_counts)
    Kmers{
        kmers: kmers,
        kmer_counts: kmer_counts,
        kmer_type: kmer_type,
    }
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
    kmers: String,
    txg_mols: Vec<String>,
    longread_mols: Vec<String>,
    hic_mols: Vec<String>,
    longread_fqs: Vec<String>,
    output: String,
    crib: Option<String>,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let variants = params.value_of("variants").unwrap();
    let ploidy = params.value_of("ploidy").unwrap_or("2");
    let ploidy: usize = ploidy.to_string().parse::<usize>().unwrap();
    let whitelist = params.value_of("whitelist").unwrap();
    let kmers = params.value_of("kmers").unwrap();
    let output = params.value_of("output").unwrap();
    let txg_tmp = match params.values_of("txg_mols") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_mols: Vec<String> = Vec::new();
    for x in txg_tmp { txg_mols.push(x.to_string()); }
    
    let longread_tmp =  match params.values_of("longread_mols") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut longread_mols: Vec<String> = Vec::new();
    for x in longread_tmp { longread_mols.push(x.to_string()); }

    let longread_tmp =  match params.values_of("longread_fqs") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut longread_fqs: Vec<String> = Vec::new();
    for x in longread_tmp { longread_fqs.push(x.to_string()); }

    let hic_tmp =  match params.values_of("hic_mols") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut hic_mols: Vec<String> = Vec::new();
    for x in hic_tmp { hic_mols.push(x.to_string()); }
    let crib: Option<String> = match params.value_of("crib") {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    Params {
        variants: variants.to_string(),
        ploidy: ploidy,
        whitelist: whitelist.to_string(),
        kmers: kmers.to_string(),
        txg_mols: txg_mols,
        longread_mols: longread_mols,
        hic_mols: hic_mols,
        longread_fqs: longread_fqs,
        output: output.to_string(),
        crib: crib,
    }
}

fn get_reader(filename: String) -> BufReader<Box<dyn Read>> {
    let filetype: Vec<&str> = filename.split(".").collect();
    let filetype = filetype[filetype.len()-1];
    let file = match File::open(filename.clone()) {
        Ok(file) => file,
        Err(error) => panic!("There was a problem opening the file: {:?}", error),
    };
    let reader: Box<dyn Read> = match filetype { 
        "gz" => Box::new(GzDecoder::new(file)), 
        _ => Box::new(file),
    }; 
    BufReader::new(reader)
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
