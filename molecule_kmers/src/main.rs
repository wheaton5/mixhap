#[macro_use]
extern crate clap;
extern crate fnv;
extern crate hashbrown;
extern crate flate2; 
extern crate itertools;
extern crate needletail;
extern crate rayon;
extern crate byteorder;
use rayon::prelude::*;
use needletail::Sequence;
use flate2::read::GzDecoder;
use std::io::BufReader;
use std::io::BufRead;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::io;
use itertools::izip;
use byteorder::{WriteBytesExt, LittleEndian};

use hashbrown::{HashMap,HashSet};
//use fnv::FnvHasher;
//use std::hash::BuildHasherDefault;
//type FnvHashMap<T,V> = HashMap<T,V, BuildHasherDefault<FnvHasher>>;
//type FnvHashSet<T> = HashSet<T, BuildHasherDefault<FnvHasher>>;
//use hashbrown::hash_map::DefaultHashBuilder;

use clap::{App};

fn main() {
    let params = load_params();
    let kmers = load_kmers(&params);
    if let Some(barcodes) = &params.txg_barcodes {
        process_txg(&params, &kmers);
    }
}

fn load_kmers(params: &Params) -> HashMap<Vec<u8>, i32> {
    let mut kmers: HashMap<Vec<u8>, i32> = HashMap::default();
    let mut kmer_id: i32 = 1; // no 0 kmer id as we are using sign for which of the pair
    if let Some(het_kmers) = &params.paired_kmers {
        let reader = File::open(het_kmers).expect("cannot open barcode file");
        let mut reader = BufReader::new(reader);
        let mut buf = vec![];
        loop {
            let bytes = reader.read_until(b'\t', &mut buf).expect("cannot read file");
            if bytes == 0 { break; } 
            kmers.insert(buf[0..(bytes-1)].to_vec(), kmer_id);
            buf.clear();
            let bytes = reader.read_until(b'\t', &mut buf).expect("cannot read file");
            if bytes == 0 { break; } 
            buf.clear();
            let bytes = reader.read_until(b'\t', &mut buf).expect("cannot read file");
            if bytes == 0 { break; } 
            kmers.insert(buf[0..(bytes-1)].to_vec(), -kmer_id);           
            buf.clear();
            let bytes = reader.read_until(b'\n', &mut buf).expect("cannot read file");
            if bytes == 0 { break; } 
            buf.clear();
            kmer_id += 1;
        }
    }
    kmers
}

fn process_txg(params: &Params, kmer_ids: &HashMap<Vec<u8>, i32>) {
    let barcodes = load_whitelist(params.txg_barcodes.as_ref().unwrap());
    let read1s = &params.txg_r1s;
    let read2s = &params.txg_r2s;
    let read1_trims = &params.txg_trim_r1s;
    let read2_trims = &params.txg_trim_r2s;
    let mut to_iterate: Vec<(String, usize, String, usize)> = Vec::new();

    for (r1_file, r1_trim, r2_file, r2_trim) in izip!(read1s, read1_trims, read2s, read2_trims) {
        to_iterate.push((r1_file.to_string(), *r1_trim, r2_file.to_string(), *r2_trim));
    }
    to_iterate.par_iter().for_each(|(r1_file, r1_trim, r2_file, r2_trim)| {
        // get readers
        let mut r1_reader = get_reader(r1_file.to_string());
        let mut r2_reader = get_reader(r2_file.to_string());
        let mut buf1 = vec![];
        let mut buf2 = vec![];
        let mut linedex = 0;
        let mut result: Vec<u8> = Vec::new();
        loop {
            let bytes = r1_reader.read_until(b'\n', &mut buf1).expect("cannot read r1file");
            let bytes2 = r2_reader.read_until(b'\n', &mut buf2).expect("cannot read r2file");
            if bytes == 0 || bytes2 == 0 { break; }
            if linedex % 4 != 1 { 
                linedex += 1;
                buf1.clear();
                buf2.clear();
                continue; 
            }
            if let Some(barcode_id) = barcodes.get(&buf1[0..16]) {
                let r1_sequence = &buf1[(16+r1_trim)..].sequence();
                let r1_sequence = r1_sequence.normalize(false);
                let r1_rc = r1_sequence.reverse_complement();
                for (_, kmer, _) in r1_sequence.canonical_kmers(21, &r1_rc) {
                //r1_sequence.canonical_kmers(21, &r1_rc).collect::<Vec<(usize, &[u8], bool)>>().into_par_iter().for_each(|(_, kmer, _)| {
                    if let Some(kmer_id) = kmer_ids.get(kmer) {
                        //println!("{}\t{}", barcode_id, kmer_id);//std::str::from_utf8(&kmer).unwrap()); 
                        result.write_i32::<LittleEndian>(*barcode_id).expect("fail");
                        result.write_i32::<LittleEndian>(*kmer_id).expect("fail");
                        std::io::stdout().lock().write_all(&result).expect("fail");
                        result.clear();
                    }
                }
                let r2_sequence = &buf2[*r2_trim..].sequence();
                let r2_sequence = r2_sequence.normalize(false);
                let r2_rc = r2_sequence.reverse_complement();
                for (_, kmer, _) in r2_sequence.canonical_kmers(21, &r2_rc) {
                //r2_sequence.canonical_kmers(21, &r2_rc).collect::<Vec<(usize, &[u8], bool)>>().into_par_iter().for_each(|(_, kmer, _)| {
                    if let Some(kmer_id) = kmer_ids.get(kmer) {
                        //println!("{}\t{}", barcode_id, kmer_id);//std::str::from_utf8(&kmer).unwrap()); 
                        let mut result: Vec<u8> = Vec::new();
                        result.write_i32::<LittleEndian>(*barcode_id).expect("fail");
                        result.write_i32::<LittleEndian>(*kmer_id).expect("fail");
                        std::io::stdout().lock().write_all(&result).expect("fail");
                        result.clear();
                    }
                }
            }
            buf1.clear();
            buf2.clear();
            linedex += 1;
        }
    });
    std::io::stdout().flush().expect("flush failed");
}

fn load_whitelist(barcode_file: &String) -> HashMap<Vec<u8>, i32> {
    let mut buf: Vec<u8> = vec![];
    let mut barcodes: HashMap<Vec<u8>, i32> = HashMap::with_capacity(20000000);
    let mut barcode_index: i32 = 1;
    let reader = File::open(barcode_file).expect("cannot open barcode file");
    let mut reader = BufReader::new(reader);
    loop {
        let bytes = reader.read_until(b'\n', &mut buf).expect("cannot read file");
        if bytes == 0 { break; }
        barcodes.insert(buf[0..(bytes-1)].to_owned(), barcode_index);
        buf.clear();
        barcode_index += 1;
        //if barcode_index % 100000 == 0 {
        //    println!("whitelist {} {}", barcode_index, bytes);
        //}
    }
    barcodes 
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


#[derive(Clone)]
struct Params {
    paired_kmers: Option<String>,
    unpaired_kmers: Option<String>,
    txg_r1s: Vec<String>,
    txg_r2s: Vec<String>,
    txg_trim_r1s: Vec<usize>,
    txg_trim_r2s: Vec<usize>,
    txg_barcodes: Option<String>,
    long_reads: Vec<String>,
    hic: Vec<String>,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let paired_kmers = match params.value_of("paired_kmers") {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    let unpaired_kmers = match params.value_of("unpaired_kmers") {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    let txg_r1s_tmp = match params.values_of("txg_r1s") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_r1s: Vec<String> = Vec::new();
    for x in txg_r1s_tmp { txg_r1s.push(x.to_string()); }
    let txg_r2s_tmp = match params.values_of("txg_r2s") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_r2s: Vec<String> = Vec::new();
    for x in txg_r2s_tmp { txg_r2s.push(x.to_string()); }
    let txg_trim_tmp = match params.values_of("txg_trim_r1s") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_trim_r1s: Vec<usize> = Vec::new();
    for trim in txg_trim_tmp {
        txg_trim_r1s.push(trim.parse::<usize>().unwrap());
    }
    let txg_trim_tmp = match params.values_of("txg_trim_r2s") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_trim_r2s: Vec<usize> = Vec::new();
    for trim in txg_trim_tmp {
        txg_trim_r2s.push(trim.parse::<usize>().unwrap());
    }
    assert!(txg_trim_r1s.len() == txg_r1s.len(), "If you specify txg_r1s, must also specify trim length for each file (txg_trim_r1s)");
    assert!(txg_trim_r2s.len() == txg_r2s.len(), "If you specify txg_r2s, must also specify trim length for each file (txg_trim_r2s)");
    assert!(txg_r2s.len() == txg_r1s.len(), "If you specify txg_r1s or txg_r2s, must also specify the same number of the other.");

    let txg_barcodes = match params.value_of("txg_barcodes") {
        Some(x) => Some(x.to_string()),
        None => {assert!(txg_r1s.len() == 0, "if you specify txg_files, you must also specify txg_barcodes");  None},
    };
    let long_reads_tmp = match params.values_of("long_reads") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut long_reads: Vec<String> = Vec::new();
    for x in long_reads_tmp { long_reads.push(x.to_string()); }
    let hic_tmp: Vec<&str> = match params.values_of("hic") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut hic: Vec<String> = Vec::new();
    for x in hic_tmp { hic.push(x.to_string()); }
    Params{
        paired_kmers: paired_kmers,
        unpaired_kmers: unpaired_kmers,
        txg_r1s: txg_r1s,
        txg_r2s: txg_r2s,
        txg_trim_r1s: txg_trim_r1s, 
        txg_trim_r2s: txg_trim_r2s,
        txg_barcodes: txg_barcodes,
        long_reads: long_reads,
        hic: hic,
    }   
}
