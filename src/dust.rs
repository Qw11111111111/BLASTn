use std::collections::{HashMap, VecDeque};

use num::ToPrimitive;

pub struct Dust {
    window_size: usize,
    threshold: f64,
    query: Vec<u8>,
    lut: HashMap<[u8; 3], usize>
}

impl Dust {
    //TODO: implement symmetric DUST
    pub fn new(window_size: usize, threshold: f64, query: Vec<u8>) -> Self {
        Self {
            window_size: window_size,
            threshold: threshold,
            query: query,
            lut: HashMap::default()
        }
    }

    pub fn mask_regions(&mut self) -> Vec<u8> {
        self.generate_lut();
        let mut current_word: VecDeque<usize> = VecDeque::new();
        let mut res : Vec<(usize, usize)> = Vec::new();

        for i in 2..self.query.len() + self.window_size - 3 {
            let wstart: usize;
            if i < self.window_size {
                wstart = 0;
            }
            else {
                wstart = i - self.window_size;
            }
            //let wend = std::cmp::min(self.query.len() - 1, i);
            if i < self.query.len() {
                let t = self.get_triplet(&self.query, i - 2);
                current_word.push_back(t);
            }
            if wstart > 0 {
                current_word.pop_front();
            }
            let (limit, mut max_score) = self.best_prefix(&current_word, 0, current_word.len() - 1);
            let mut start = 0;
            let mut finish = limit;
            if max_score > self.threshold {
                for s in 1..limit - 3 {
                    let (f, new_score) = self.best_prefix(&current_word, s, limit - 2);
                    if new_score > max_score {
                        max_score = new_score;
                        start = s;
                        finish = f;
                    }
                }
                res.push((start + wstart, finish + wstart));
            }
        }
        self.get_masked_seq(&res)

    }

    fn generate_lut(&mut self) {
        let bases: [u8; 4] = [65, 67, 71, 84];
        let mut count = 0;
        for &base1 in bases.iter() {
            for &base2 in bases.iter() {
                for &base3 in bases.iter() {
                    self.lut.insert([base1, base2, base3], count);
                    count += 1;
                }
            }
        }
    }

    fn get_triplet(&self, word: &[u8], idx: usize) -> usize {
        self.lut[&word[idx..idx + 3]]
    }

    fn best_prefix(&self, word: &VecDeque<usize>, istart: usize, ifinsh: usize) -> (usize, f64) {
        let mut scores = vec![0.0; 64];
        let mut running_score = 0.0;
        let mut finish = istart - 1;
        let mut max_score = 0.0;
        for i in istart..ifinsh {
            let t = word[i];
            running_score += scores[t];
            scores[t] += 1.0;
            if i > istart && running_score / (i - istart).to_f64().unwrap() > max_score {
                max_score = running_score / (i - istart).to_f64().unwrap();
                finish = i;
            }
        }
        (finish + 2, max_score)
    }

    fn get_masked_seq(&self, res: &Vec<(usize, usize)>) -> Vec<u8> {
        let mut masked = self.query.to_vec();

        for (start, finish) in res.iter() {
            for i in *start..*finish {
                masked[i] = 78;
            }
        }
        masked
    }
}