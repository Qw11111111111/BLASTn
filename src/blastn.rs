use std::{
    collections::HashMap,
    fs::File,
    thread,
    io::BufReader,
    sync::{Arc, Mutex, RwLock}
};

use rand::{seq::index, thread_rng, Rng};


use bio::{data_structures, io::fasta::{Reader, Record, Records}};

pub struct Searcher {
    threshhold: u32,
    n_seeds: usize,
    word_start: usize,
    word_len: usize,
    best_hits: HashMap<usize, HashMap<usize, u32>>,
    query: Arc<Record>,
    db: Arc<Mutex<Records<BufReader<File>>>>,
    word: Arc<RwLock<Vec<u8>>>
}
fn get_score(seq1: Vec<u8>, seq2: &[u8]) -> u32 {
    seq2.iter().enumerate().map(|(i, &item)| {
        if seq1[i] == item {
            1
        }
        else {
            0
        }
    }).sum()
}

fn add_to_hits(hits: &mut HashMap<usize, u32>, score: u32, index: usize, max_size: usize) {
    hits.insert(index, score);
}

impl Searcher {

    pub fn new(query: Arc<Record>, db: Arc<Mutex<Records<BufReader<File>>>>, threshhold: u32, n_seeds: usize, word_len: usize) -> Self {
        Self {
            threshhold: threshhold,
            n_seeds: n_seeds,
            word_len: word_len,
            word_start: 0,
            best_hits: HashMap::new(),
            query: query,
            db: db,
            word: Arc::from(RwLock::from(vec![])),
        }
    }

    fn get_word(&mut self) {
        let mut rng = thread_rng();
        let length = self.query.seq().len();
        println!("l: {}", length);
        self.word_start = rng.gen_range(0..length - self.word_len);
        let mut word = self.word.write().unwrap();
        for i in self.word_start..self.word_len + self.word_start{
            word.push(self.query.seq()[i]);
        }
    }   

    pub fn align(&mut self) {
        self.get_word();
        self.seed();
        let val = &self.best_hits[&0];
        
        let val = &self.best_hits[&0];
        for (key, item) in val.iter() {
            println!("{:?}", (key, item));
        }
        
        println!("{}", val.len());
        /*
        self.set_query_as_word();

        let mut handles = vec![];
        let word_start = self.word_start;
        let mut database = self.db.try_lock().unwrap();
        let mut idx = 0;
        
        while let Some(record) = database.next() {
            let record = record.unwrap();
            let word = self.word.read().unwrap().clone();
            //let map = Arc::from(Mutex::from(self.best_hits[&idx]));
            let mut map = self.best_hits[&idx].clone();
            handles.push(thread::spawn(move || {
                //for (index, score) in map.try_lock().unwrap().iter_mut() {
                for (index, score) in map.iter_mut() {
                    *score = get_score(word.clone(), &record.seq()[index - word_start..index + word.len() - word_start])
                }
                map
            }));
            idx += 1;
        }
        let mut empty_map = HashMap::new();
        println!("{}", handles.len());
        for (i, handle) in handles.into_iter().enumerate() {
            let hits = handle.join().unwrap();
            empty_map.insert(i, hits);
        }
        self.best_hits = empty_map;
        let val = &self.best_hits[&0];
        for (key, item) in val.iter() {
            println!("{:?}", (key, item));
        }
        */

    }

    fn seed(&mut self) {
        let mut handles = vec![];
        let word_start = self.word_start;
        let word_len = self.word_len;
        let threshhold = self.threshhold;
        let word2 = self.set_query_as_word();
        //let database = self.db.clone();

        let mut database = self.db.try_lock().unwrap();

        while let Some(record) = database.next() {
            let record = record.unwrap();
            let word = self.word.read().unwrap().clone();
            let next_word = word2.clone();
            handles.push(thread::spawn(move || {
                let mut hits: HashMap<usize, u32> = HashMap::new();
                for i in word_start..record.seq().len() - word_len - word_start {
                    let score = get_score(word.clone(), &record.seq()[i..i + word_len]);
                    if score >= threshhold {
                        add_to_hits(&mut hits, score, i, 0);
                    }
                }
                println!("done");
                for (index, score) in hits.iter_mut() {
                    println!("{index}");
                    println!("{}", next_word.len());
                    println!("rec: {}", record.seq().len());
                    *score = get_score(next_word.clone(), &record.seq()[index - word_start..index + &next_word.len() - word_start])
                }
                println!("done2");
                hits
            }))
        }
        for (i, handle) in handles.into_iter().enumerate() {
            let hits = handle.join().unwrap();
            self.best_hits.insert(i, hits);
        }
    }

    fn set_query_as_word(&mut self) -> Vec<u8> {
        let mut word = vec![];
        for &val in self.query.seq().iter() { 
            word.push(val);
        }
        word
    }

    pub fn summary(&self, db: &mut Records<BufReader<File>>) -> Summary {
        let mut max = 0;
        let mut best_idx:(usize, usize) = (0, 0);
        for (&idx, map) in self.best_hits.iter() {
            for (&index, &score) in map.iter() {
                if score > max {
                    max = score;
                    best_idx = (idx, index);
                }
            }
        }
        println!("{}", best_idx.0);
        //let mut db = self.db.try_lock().unwrap();
        let binding = db.nth(best_idx.0).unwrap().unwrap();
        let seq = &binding.seq()[best_idx.1 - self.word_start..best_idx.1 - self.word_start + self.query.seq().len()];

        let mut word = vec![];
        for &val in seq.iter() { 
            word.push(val);
        }

        Summary {
            best_idx: best_idx,
            score: max,
            seq: word
        }
    }
}

pub struct Summary {
    pub best_idx: (usize, usize),
    pub score: u32,
    pub seq: Vec<u8>
}