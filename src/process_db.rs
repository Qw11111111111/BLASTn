use std::{
    collections::{BTreeMap, HashMap}, 
    fs::File, 
    io::BufReader,
    thread
};
use bio::io::fasta::Records;


pub fn get_kmers(k: usize, mut db: Records<BufReader<File>>) -> HashMap<usize, BTreeMap<u64, Vec<usize>>> {
    let mut handles = Vec::new();
    while let Some(record) = db.next() {
        handles.push(thread::spawn(move || {
            let mut kmers: HashMap<u64, Vec<usize>> = HashMap::new();
            let record = record.unwrap();
            for i in 0..record.seq().len() - k {
                let kmer = &record.seq()[i..i + k];
                let key: u64 = kmer.iter().enumerate().map(|(i, nt)| {
                    *nt as u64 * 4_u64.pow(i as u32)
                }).sum();
                let current = kmers.entry(key).or_insert_with(Vec::default);
                current.push(i);
            }
            kmers
        }));
    }

    let mut tree_map: HashMap<usize, BTreeMap<u64, Vec<usize>>> = HashMap::new();

    for (i, handle) in handles.into_iter().enumerate() {
        let map = handle.join().unwrap();
        let mut new_tree = BTreeMap::new();
        for (key, idx) in map.into_iter() {
            new_tree.insert(key, idx);
        }
        tree_map.insert(i, new_tree);
    }

    tree_map
}

pub fn save_processed_db() {
    //TODO
}

pub fn load_processed_db () {
    //TODO
}