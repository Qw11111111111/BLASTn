use crate::blastn::Searcher;

use std::{
    collections::HashMap, sync::{Arc, Mutex}, time::Instant
};

use bio::io::fasta::Reader;

pub fn benchmark(n_retries: usize, query_file: &str, databank: &str, threshhold: &usize, word_len: &usize) -> (f64, HashMap<(usize, usize), usize>) {
    let mut total_time = 0.0;
    let mut best_hits: HashMap<(usize, usize), usize> = HashMap::default();
    
    for _ in 0..n_retries {
        let query_reader = Reader::from_file(&query_file).unwrap();
        let db_reader = Reader::from_file(&databank).unwrap();
        
        let db = Arc::from(Mutex::from(db_reader.records()));
        let query = Arc::from(query_reader.records().next().unwrap().unwrap());
        
        let mut searcher = Searcher::new(query.clone(), db, threshhold.clone(), word_len.clone());
    
        let now = Instant::now();
        
        searcher.align();
        total_time += now.elapsed().as_secs_f64();
        
        let db_reader = Reader::from_file(&databank).unwrap();
        let summary = searcher.summary(&mut db_reader.records());

        best_hits.insert(summary.best_idx, summary.score);
    };  
    (total_time, best_hits)
}