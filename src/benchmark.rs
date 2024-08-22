use crate::blastn::{Searcher, Summary};
use crate::dust::Dust;

use std::{
    collections::HashMap, 
    sync::{Arc, Mutex}, 
    time::Instant, 
    fmt
};

use bio::io::fasta::Reader;
use num::ToPrimitive;

pub fn benchmark(n_retries: usize, query_file: &str, databank: &str, threshhold: &usize, word_len: &usize, m_threshold: f64, mask: bool) -> Result<TopTen, String> {
    let mut total_time = 0.0;
    let mut top = TopTen::default();
    
    for _ in 0..n_retries {
        let query_reader = Reader::from_file(&query_file).unwrap();
        let db_reader = Reader::from_file(&databank).unwrap();
        
        let db = Arc::from(Mutex::from(db_reader.records()));
        let query = Arc::from(query_reader.records().next().unwrap().unwrap());
        let masked: Vec<u8>;
        if mask {
            let mut dust = Dust::new(78, m_threshold, query.seq().to_vec());
            masked = dust.mask_regions();
        }
        else {
            masked = query.seq().to_vec();
        }

        let mut searcher = Searcher::new(query.clone(), db, threshhold.clone(), word_len.clone(), masked);
    
        let now = Instant::now();
        
        searcher.align()?;
        total_time += now.elapsed().as_secs_f64();
        
        let db_reader = Reader::from_file(&databank).unwrap();
        let summary = searcher.summary(&mut db_reader.records());
        top.insert(summary);
    };  
    top.time = total_time;
    top.retries = n_retries;
    top.sort();
    top.show_time = true;
    Ok(top)
}

#[derive (Debug, Default)]
pub struct TopTen {
    pub time: f64, 
    pub hits: HashMap<usize, Summary>,
    pub keys: Vec<usize>,
    pub min: usize,
    min_idx: usize,
    pub retries: usize,
    pub show_time: bool,
    pub hits_found: Vec<usize>,
    pub no_result: usize
}

impl TopTen {

    pub fn sort(&mut self) {
        for i in 1..self.keys.len() {
            for j in 0..self.keys.len() - i {
                let item = &self.hits[&self.keys[j]].score;
                let next_item = &self.hits[&self.keys[j + 1]].score;
                if next_item > item {
                    (self.keys[j], self.keys[j + 1]) = (self.keys[j + 1], self.keys[j]);
                }
            }
        }
    }

    pub fn insert(&mut self, summary: Summary) {
        if summary.similarity == 0.0 {
            self.no_result += 1;
            return;
        }

        for key in self.keys.iter() {
            if self.hits[key] == summary {
                self.hits_found[*key] += 1;
                return;
            }
        }
        if self.keys.len() == 10 {
            if summary.score <= self.min {
                return;
            }
            self.hits.insert(self.min_idx, summary);
            self.hits_found[self.min_idx] = 1;
            let mut new_min = self.hits[&0].query.len(); // reference value
            for (i, key) in self.keys.iter().enumerate() {
                let s = self.hits[key].score;
                if s < new_min {
                    new_min = s;
                    self.min_idx = i;
                }
            }
        }
        else {
            if summary.score < self.min {
                self.min = summary.score;
                self.min_idx = self.keys.len();
            }
            self.hits.insert(self.keys.len(), summary);
            self.keys.push(self.keys.len());
            self.hits_found.push(1);
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
        if self.show_time {
            writeln!(f, "\nTotal time used: {}s", self.time)?;
            writeln!(f, "time per run: {}s\n", self.time / self.retries.to_f64().unwrap())?;
        }
        writeln!(f, "\nbest hits: \n")?;

        for (i, key) in self.keys.iter().enumerate() {
            let s = &self.hits[key];
            writeln!(f, "{i}: Record: {0:>15} | Idx: {1:>8} | Times found: {3:>5} | Similarity to masked query: {2} ", s.id, s.best_idx.1, s.similarity, self.hits_found[*key])?;
        };

        if self.keys.len() == 0 {
            writeln!(f, "No Match found")?;
        }

        writeln!(f, "\n{} search(es) produced no result", self.no_result)?;

        write!(f, "\n")
    }
}