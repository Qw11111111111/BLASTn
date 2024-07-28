use std::{
    collections::{BTreeMap, HashMap}, fmt, fs::File, io::BufReader, sync::{Arc, Mutex}, thread  
};

use num::ToPrimitive;

use rand::{thread_rng, Rng};

use bio::io::fasta::{Record, Records, Reader};

use crate::benchmark::TopTen;


fn get_hsps(q: Arc<HashMap<u64, Vec<usize>>>, db: &BTreeMap<u64, Vec<usize>>) {

}

fn extend_match(q: &[u8], db: &[u8], t: usize, rev: bool) -> (usize, u64) {
    if rev {
        let mut s = 0;
        for (i, &base) in q.iter().enumerate().rev() {
            if base == db[i] {
                s += 1;
            }
            else {
                if s < t.to_u64().unwrap() {
                    return (i, 0);
                }
                s -= 1;
            }
        }
        (q.len(), s)
    }
    else {
        let mut s = 0;
        for (i, &base) in q.iter().enumerate() {
            if base == db[i] {
                s += 5;
            }
            else {
                if s < t.to_u64().unwrap() {
                    return (i, 0);
                }
                s -= 4
            }
        }
        (q.len(), s)
    }
}

#[derive (Clone, Debug)]
pub struct Match {
    wstart: usize,
    idx_in_db: usize,
    score: u64,
    similarity: f64,
    pub start_stop: (usize, usize)
}

pub struct BLASTn {
    masked_query: Arc<Vec<u8>>,
    hssp_threshold: usize,
    elongation_threshold: usize,
    k: usize,
    db: Records<BufReader<File>>
}

impl BLASTn {
    
    pub fn new(db: Records<BufReader<File>>, masked_query: Arc<Vec<u8>>, hssp_threshold: usize, elongation_threshold: usize, k: usize) -> Self {
        Self {
            masked_query: masked_query,
            hssp_threshold: hssp_threshold,
            elongation_threshold: elongation_threshold,
            k: k,
            db: db
        }
    }

    pub fn align(&mut self, db_kmers: Arc<HashMap<usize, BTreeMap<u64, Vec<usize>>>>) -> HashMap<usize, Vec<Match>> {
        let mut handles = Vec::new();
        let kmers = Arc::from(self.get_all_kmers());
        let mut idx = 0;
        let t = self.elongation_threshold;


        while let Some(record) = self.db.next() {
            let kmers = kmers.clone();
            let kmers_in_rec = Arc::from(db_kmers[&idx].clone());
            let rec = record.unwrap();
            let q = self.masked_query.clone();
            idx += 1;
            handles.push(thread::spawn(move || {
                let mut matches: Vec<Match> = Vec::new();
                for (key, idx) in kmers.iter() {
                    for idx_in_q in idx.iter() {
                        if let Some(index) = kmers_in_rec.get(key) {
                            for idx_in_db in index.iter() {
                                if idx_in_db < idx_in_q || (idx_in_db - idx_in_q + q.len()) > rec.seq().len() {
                                    continue;
                                }
                                let score1 = extend_match(&q[..*idx_in_q], &rec.seq()[idx_in_db - idx_in_q..*idx_in_db], t, true);
                                let score2 = extend_match(&q[*idx_in_q..], &rec.seq()[*idx_in_db..idx_in_db - idx_in_q + q.len()], t, false);
                                matches.push(Match { wstart: *idx_in_q, idx_in_db: *idx_in_db, score: score1.1 + score2.1, similarity: 0.0, start_stop: (idx_in_db - score1.0, score2.0 + idx_in_db) })
                            }
                        }
                    }
                }
                matches
            }));
        }
        let mut all_matches = HashMap::new();
        for (i, handle) in handles.into_iter().enumerate() {
            let matches = handle.join();
            let mut matches = matches.unwrap();
            matches.sort_by(|a, b| (a.start_stop.1 - a.start_stop.0).cmp(&(b.start_stop.1 - b.start_stop.0)));
            all_matches.insert(i, matches);
        }
        all_matches
    }

    fn get_all_kmers(&self) -> HashMap<u64, Vec<usize>> {
        let mut kmers: HashMap<u64, Vec<usize>> = HashMap::new();
        for i in 0..self.masked_query.len() - self.k {
            let kmer = &self.masked_query[i..i + self.k];
            let key: u64 = kmer.iter().enumerate().map(|(i, nt)| {
                nt.to_u64().unwrap() * 4_u64.pow(i.to_u32().unwrap())
            }).sum();
            if kmers.contains_key(&key) {
                let mut new_idx = kmers[&key].clone();
                new_idx.push(i);
                kmers.insert(key,  new_idx);
            }
            else {
                kmers.insert(key, vec![i]);
            }
       }
       kmers
    }

    fn get_hssps(&self) -> HashMap<u64, Vec<usize>> {
        let mut kmers = HashMap::new();


        kmers
    }   
}



fn get_score(seq1: &[u8], seq2: &[u8]) -> usize {
    let mut n = 0;
    for (i, &item) in seq1.iter().enumerate() {
        if item == seq2[i] {
            n += 1;
        }
    }
    n
}

pub struct Searcher {
    threshold: usize,
    word_start: usize,
    word_len: usize,
    best_hits: HashMap<usize, HashMap<usize, usize>>,
    query: Arc<Record>,
    db: Arc<Mutex<Records<BufReader<File>>>>,
    word: Arc<Vec<u8>>,
    masked: Vec<u8>
}

impl Searcher {

    pub fn new(query: Arc<Record>, db: Arc<Mutex<Records<BufReader<File>>>>, threshold: usize, word_len: usize, masked: Vec<u8>) -> Self {
        Self {
            threshold: threshold,
            word_len: word_len,
            word_start: 0,
            best_hits: HashMap::default(),
            query: query,
            db: db,
            word: Arc::default(),
            masked: masked
        }
    }

    fn get_word(&mut self) {
        let all_words = self.list_all_words();
        if all_words.len() == 0 {
            println!("No kmer of length {} could be found in the query sequence. Retry with lower word length, higher masking threshold or without masking", self.word_len);
            self.word = Arc::from(Vec::new());
            return;
        }
        let mut rng = thread_rng();
        loop {
            let idx = rng.gen_range(0..self.masked.len() - self.word_len);
            match all_words.get(&idx) {
                None => continue,
                Some(word) => {
                    self.word = Arc::from(word.to_vec());
                    self.word_start = idx;
                    break;
                }
            }
        } 
    }

    fn list_all_words(&self) -> HashMap<usize, &[u8]> {
        let mut all_words = HashMap::new();

        for i in 0..self.masked.len() - self.word_len {
            if self.masked[i..i + self.word_len].contains(&78) {
                continue;
            }
            all_words.insert(i, &self.masked[i..i + self.word_len]);
        }
        all_words
    }

    pub fn align(&mut self) {
        self.get_word();
        let mut handles = vec![];
        let word_start = self.word_start;
        let word_len = self.word_len;
        let threshold = self.threshold;
        let word2: Arc<[u8]> = Arc::from(self.query.seq().to_vec());
        //let word2: Arc<[u8]> = Arc::from(self.masked.clone());
        let mut database = self.db.try_lock().unwrap();

        while let Some(record) = database.next() {
            let rec = record.unwrap();
            let word = self.word.clone();
            let next_word = word2.clone();
            handles.push(thread::spawn(move || {
                let mut hits: HashMap<usize, usize> = HashMap::default();
                for i in word_start..rec.seq().len() + word_len - word_start - next_word.len() {
                    let mut score = get_score(&word, &rec.seq()[i..i + word_len]);
                    if score >= threshold {
                        score = get_score(&next_word, &rec.seq()[i - word_start..i + &next_word.len() - word_start]);
                        //score += get_score(&next_word[..word_start], &rec.seq()[i - word_start..i]);
                        //score += get_score(&next_word[word_start + word_len..], &rec.seq()[i + word_len..i - word_start + &next_word.len()]); //this is faster for very long word lengths or very low threshholds
                        hits.insert(i, score);
                    }
                }
                hits
            }));
        }
        for (i, handle) in handles.into_iter().enumerate() {
            let hits = handle.join().unwrap();
            self.best_hits.insert(i, hits);
        }
    }

    pub fn summary(&self, db: &mut Records<BufReader<File>>) -> Summary {
        //TODO: increase speed
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

        let mut word: Vec<u8> = vec![];
        let similarity: f64;
        let name: String;
        if max == 0 {
            similarity = 0.0;
            name = "".to_string();
        }
        else {
            similarity = max.to_f64().unwrap() / self.query.seq().len().to_f64().unwrap() * 100.0;
            let binding = db.nth(best_idx.0).unwrap().unwrap();
            word = binding.seq()[best_idx.1..best_idx.1 + self.query.seq().len()].to_vec();
            name = binding.id().to_string();
        }

        Summary {
            best_idx: best_idx,
            score: max,
            seq: word,
            query: self.query.seq().to_vec(),
            similarity: similarity,
            id: name
        }
    }
    
    pub fn sm(&self, db: &str) -> TopTen {
        //TODO: rework for increased efficiancy
        let mut tt = TopTen::default();
        let mut recs = vec![];
        for rec in Reader::from_file(db).unwrap().records() {
            recs.push(rec.unwrap());
        }
        for (i, map) in self.best_hits.iter() {
            for (j, item) in map.iter() {
                if tt.hits.len() == 10 {
                    if *item < tt.min {
                        continue;
                    }
                }
                let idx = (*i, j - self.word_start);
                let word = recs[idx.0].seq()[idx.1..idx.1 + self.query.seq().len()].to_vec();
                let name = recs[idx.0].id().to_string();
                let similarity = item.to_f64().unwrap() / self.query.seq().len().to_f64().unwrap() * 100.0;
                tt.insert(Summary { best_idx: idx, score: *item, seq: word, similarity: similarity, query: self.query.seq().to_vec(), id: name })
            }
        }
        tt.sort();
        if tt.hits.is_empty() {
            tt.no_result = 1;
        }
        tt
    }
    
}

#[derive (Debug, PartialEq)]
pub struct Summary {
    pub best_idx: (usize, usize),
    pub score: usize,
    pub seq: Vec<u8>,
    pub similarity: f64,
    pub query: Vec<u8>,
    pub id: String
}

impl fmt::Display for Summary {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {

        writeln!(f, "Best match found in record {} [{}] at index {}", self.id, self.best_idx.0, self.best_idx.1)?;
        writeln!(f, "Score: {}", self.score)?;
        writeln!(f, "\nQuery: \n")?;
        for char in self.query.iter() {
            write!(f, "{}", convert_to_ascii(char))?;
        }
        writeln!(f, "\n\nBest match: \n")?;
        for char in self.seq.iter() {
            write!(f, "{}", convert_to_ascii(char))?;
        }
        writeln!(f, "\n\nSimilarity: {}%", self.similarity)
    }
}

pub fn convert_to_ascii(index: &u8) -> String {
    match index {
        65 => "A".into(),
        67 => "C".into(),
        71 => "G".into(),
        84 => "T".into(),
        78 => "N".into(),
        _ => "?".into(),
    }
}