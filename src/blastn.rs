use std::{
    collections::HashMap, fmt, fs::File, io::BufReader, sync::{Arc, Mutex, RwLock}, thread
};

use clap::builder::Str;
use num::ToPrimitive;

use rand::{thread_rng, Rng};

use bio::io::fasta::{Record, Records};

fn get_score(seq1: &Vec<u8>, seq2: &[u8]) -> u32 {
    seq2.iter().enumerate().map(|(i, &item)| {
        if seq1[i] == item {
            1
        }
        else {
            0
        }
    }).sum()
}

fn add_to_hits(hits: &mut HashMap<usize, u32>, score: u32, index: usize) {
    hits.insert(index, score);
}


pub struct Searcher {
    threshhold: u32,
    word_start: usize,
    word_len: usize,
    best_hits: HashMap<usize, HashMap<usize, u32>>,
    query: Arc<Record>,
    db: Arc<Mutex<Records<BufReader<File>>>>,
    word: Arc<RwLock<Vec<u8>>>
}

impl Searcher {

    pub fn new(query: Arc<Record>, db: Arc<Mutex<Records<BufReader<File>>>>, threshhold: u32, word_len: usize) -> Self {
        Self {
            threshhold: threshhold,
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
        self.word_start = rng.gen_range(0..length - self.word_len);
        let mut word = self.word.write().unwrap();
        for i in self.word_start..self.word_len + self.word_start{
            word.push(self.query.seq()[i]);
        }
    }   

    pub fn align(&mut self) {
        self.get_word();
        let mut handles = vec![];
        let word_start = self.word_start;
        let word_len = self.word_len;
        let threshhold = self.threshhold;
        let word2 = self.set_query_as_word();

        let mut database = self.db.try_lock().unwrap();
        let mut num = 0;

        while let Some(record) = database.next() {
            let rec = record.unwrap();
            let word = self.word.read().unwrap().clone();
            let next_word = word2.clone();
            handles.push(thread::spawn(move || {
                let num2 = num.clone();
                let mut hits: HashMap<usize, u32> = HashMap::new();
                for i in word_start..rec.seq().len() + word_len - word_start - next_word.len() {
                    let mut score = get_score(&word, &rec.seq()[i..i + word_len]);
                    if score >= threshhold {
                        score = get_score(&next_word, &rec.seq()[i - word_start..i + &next_word.len() - word_start]);
                        add_to_hits(&mut hits, score, i);
                    }
                }
                println!("thread {num2} done");
                hits
            }));
            num += 1;
        }
        for (i, handle) in handles.into_iter().enumerate() {
            let hits = handle.join().unwrap();
            self.best_hits.insert(i, hits);
        }
    }

    fn set_query_as_word(&self) -> Vec<u8> {
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
                    best_idx = (idx, index - self.word_start);
                }
            }
        }


        let query = self.set_query_as_word();
        let mut word: Vec<u8> = vec![];
        let similarity: f64;
        let name: String;
        if max == 0 {
            similarity = 0.0;
            name = "".to_string();
        }
        else {
            similarity = max.to_f64().unwrap() / query.len().to_f64().unwrap() *100.0;
            
            let binding = db.nth(best_idx.0).unwrap().unwrap();
            let seq = &binding.seq()[best_idx.1..best_idx.1 + self.query.seq().len()];
            name = binding.id().to_string();
            for &val in seq.iter() { 
                word.push(val);
        }
        }

        Summary {
            best_idx: best_idx,
            score: max,
            seq: word,
            query: query,
            similarity: similarity,
            id: name
        }
    }
}

pub struct Summary {
    pub best_idx: (usize, usize),
    pub score: u32,
    pub seq: Vec<u8>,
    pub similarity: f64,
    pub query: Vec<u8>,
    pub id: String
}

impl fmt::Display for Summary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {

        writeln!(f, "Best match found in record: {} at index: {}", self.best_idx.0, self.best_idx.1)?;
        writeln!(f, "Score: {}", self.score)?;
        writeln!(f, "\nQuery: \n")?;
        for char in self.query.iter() {
            write!(f, "{}", convert_to_ascii(char))?;
        }
        writeln!(f, "\n\nBest match: \n")?;
        for char in self.seq.iter() {
            write!(f, "{}", convert_to_ascii(char))?;
        }
        writeln!(f, "\n\nSimilarity: {}%", self.similarity)?;

        writeln!(f, "\nID: {}", self.id)
    }
}

fn convert_to_ascii(index: &u8) -> String {
    match index {
        65 => "A".into(),
        67 => "C".into(),
        71 => "G".into(),
        84 => "T".into(),
        _ => "?".into(),
    }
}