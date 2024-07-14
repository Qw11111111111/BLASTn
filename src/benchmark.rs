use crate::blastn::{Searcher, Summary};

use std::{
    collections::HashMap, sync::{Arc, Mutex}, time::Instant, fmt
};

use bio::io::fasta::Reader;
use num::ToPrimitive;

pub fn benchmark(n_retries: usize, query_file: &str, databank: &str, threshhold: &usize, word_len: &usize) -> TopTen {
    let mut total_time = 0.0;
    let mut top = TopTen::default();
    
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
        top.insert(summary);
    };  
    top.time = total_time;
    top.retries = n_retries;
    top.sort();
    top
}

#[derive (Debug, Default)]
pub struct TopTen {
    pub time: f64, 
    pub hits: HashMap<usize, Summary>,
    pub keys: Vec<usize>,
    min: usize,
    min_idx: usize,
    retries: usize
}

impl TopTen {

    fn sort(&mut self) {
        for i in 1..self.keys.len() {
            for j in 0..self.keys.len() - i {
                let key = self.keys[j];
                let item = &self.hits[&key].score;
                let next_key = self.keys[j + 1];
                let next_item = &self.hits[&next_key].score;
                if next_item > item {
                    self.keys[i + 1] = key;
                    self.keys[i] = next_key;
                }
            }
        }
    }

    fn insert(&mut self, summary: Summary) {
        for (_k, s) in self.hits.iter() {
            if &summary == s {
                return;
            }
        }
        if self.keys.len() == 10 {
            if summary.score <= self.min {
                return;
            }
            self.hits.insert(self.min_idx, summary);
            let mut new_min = self.hits[&0].query.len();
            for (i, key) in self.keys.iter().enumerate() {
                let s = self.hits[key].score;
                if s < new_min {
                    new_min = s;
                    self.min_idx = i;
                }
            }
        }
        else {
            self.keys.push(self.keys.len());
            self.hits.insert(self.keys[self.keys.len() - 1], summary);
        }
    }

    pub fn print(&self, idx: usize) {
        if idx >= self.keys.len() {
            println!("\nIndex out of bounds\n");
            return;
        }
        println!("{}", self.hits[&self.keys[idx]]);
    }
}

impl fmt::Display for TopTen {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "\nTotal time used: {}s", self.time)?;
        writeln!(f, "time per run: {}s\n", self.time / self.retries.to_f64().unwrap())?;
        writeln!(f, "best hits: \n")?;

        for (i, key) in self.keys.iter().enumerate() {
            let s = &self.hits[key];
            writeln!(f, "{i}: Record: {0:>5} | Idx: {1:>8} | Similarity: {2}", s.id, s.best_idx.1, s.similarity)?;
        };
        write!(f, "\n")
    }
}