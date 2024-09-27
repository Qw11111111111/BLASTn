use std::{cell::RefCell, fmt, io, rc::Rc, sync::{mpsc, Arc}, thread::{self, JoinHandle}, time::Instant};

use crate::{make_db::{parse_fasta::{parse_small_fasta, parse_to_bytes}, read_db::{bytes_to_chars, parse_compressed_db_lazy, read_csv}, records::{Record, SimpleRecord}}, dust::Dust, blastn::convert_to_ascii};

//TODO: -precalculate all div_ceils, where possible
//      -implement support for word lengths k where k % 4 != 0
//      -implement support for CLI flags
//TODO: -No need to keep actual strings in HSP. Just keep track of idx in query, current length and gap indeces.

pub struct Params {
    pub k: usize,
    pub extension_threshhold: i16,
    pub scanning_threshhold: i16,
    pub extension_length: usize,
    pub verbose: bool,
    pub masking_threshhold: f64,
    pub masking: bool
}

#[derive (Default)]
struct ProcessedChunk {
    bytes: Vec<u8>,
    end_bit: usize,
    end_byte: usize,
    start_bit: usize, // what is the first relevant bit in the record
    id: String,
    start_in_rec: u128, // where is the first byte of the chunk releative to the record,
    start_byte: usize
}

struct Info {
    id: Vec<String>,
    length: Vec<usize>
}

impl Info {
    fn get_length_of(&self, id: &str) -> usize {
        for (i, l) in self.length.iter().enumerate() {
            if self.id[i] == id {
                return *l;
            }
        }
        0
    }
}

#[derive (Debug, Default, Clone)]
struct HSP {
    id: String,
    word: Vec<u8>,
    db_seq: Vec<u8>,
    idx_in_query: usize,
    idx_in_record: usize,
    padding_left: usize,
    padding_right: usize,
    score: i16,
    e_val: f64,
    is_extended_left: bool,
    is_extended_right: bool,
    is_joined: bool,
}

impl HSP {

    fn try_join(&mut self, other: &mut HSP, scheme: &Arc<ScoringScheme>, params: &Arc<Params>, query: &Arc<SimpleRecord>, max_d: usize) -> bool {
        //tries to joins other to the right side of self
        //TODO: expand logic to all cases
        if other.is_joined || self.id != other.id || self.idx_in_query > other.idx_in_query || self.idx_in_record + self.db_seq.len() < other.idx_in_record - 1 || self.idx_in_record > other.idx_in_record {
            return false;
        }
        if self.idx_in_query + self.word.len() >= other.idx_in_query && self.idx_in_record + self.padding_left + self.word.len() >= other.idx_in_record + other.padding_left {
            self.join_overlapping(other, scheme).expect("could not join overlapping");
        }
        else if (!other.is_extended_right || !self.is_extended_left) && self.idx_in_query + self.word.len() < other.idx_in_query - 1 {
            self.try_join_separate(other, scheme, params, query, max_d).expect("could not join separate");
        }
        else {
            self.join_overlapping(other, scheme).expect("could not join overlapping");
        }
        true
    }

    fn join_overlapping(&mut self, other: &mut HSP, scheme: &Arc<ScoringScheme>) -> Result<(), &str> {
        let overlap = (self.idx_in_query + self.word.len()).abs_diff(other.idx_in_query);
        if overlap > self.word.len() || overlap > other.word.len() {
            return Ok(());
        }
        let current_score = get_score(&self.word[self.word.len() - overlap..], &self.db_seq[self.db_seq.len() - self.padding_right - overlap..self.db_seq.len() - self.padding_right], scheme);
        let other_score = get_score(&other.word[..overlap], &other.db_seq[other.padding_left..other.padding_left + overlap], scheme);

        if other.idx_in_record + other.db_seq.len() > self.idx_in_record + self.db_seq.len() {
            self.db_seq.extend_from_slice(&other.db_seq[self.db_seq.len() - (other.idx_in_record - self.idx_in_record)..]);
            self.padding_right = other.padding_right;
        }
        else {
            self.padding_right = self.db_seq.len() + overlap - self.padding_left - self.word.len() - other.word.len();
        }
        
        if other_score > current_score {
            self.score = self.score - current_score + other_score + get_score(&other.word[overlap..], &other.db_seq[other.padding_left + overlap..other.padding_left + other.word.len()], scheme);
            self.word.splice(self.word.len() - overlap.., other.word.clone());
        }
        else {
            self.score += get_score(&other.word[overlap..], &other.db_seq[other.padding_left + overlap..other.padding_left + other.word.len()], scheme);
            self.word.extend_from_slice(&other.word[overlap..]);
        }
        other.is_joined = true;
        Ok(())
    }

    fn try_join_separate(&mut self, other: &mut HSP, scheme: &Arc<ScoringScheme>, params: &Arc<Params>, query: &Arc<SimpleRecord>, max_d: usize) -> Result<(), &str> {
        let distance = self.idx_in_query + self.word.len() - other.idx_in_query;
        let max_gaps = distance - (self.idx_in_record + self.padding_left + self.word.len() - (other.idx_in_record + other.padding_left));
        if !self.is_extended_left {
            self.extend_left(scheme, query, params, max_gaps, distance);
        }
        if !other.is_extended_right {
            other.extend_right(scheme, query, params, max_gaps, distance);
        }
        if self.idx_in_query + self.word.len() >= other.idx_in_query {
            self.join_overlapping(other, scheme)?;
        }
        else if other.idx_in_query - (self.idx_in_query + self.word.len()) < max_d {
        }
        Ok(())
    }

    fn _extend_gapped(&mut self, _max_gaps: usize) {
        //self.is_extended = true;
    }
    //TODO: bounds checking
    fn extend_left(&mut self, scheme: &Arc<ScoringScheme>, query: &Arc<SimpleRecord>, params: &Arc<Params>, mut max_gaps: usize, mut max_extension: usize) {
        let mut  max_score = self.score;
        let mut best_extension: usize = 0;
        let mut extended = 0;
        max_extension = max_extension.min(self.padding_left);
        for i in 0..max_extension {
            if self.idx_in_query == 0 || self.padding_left == 0 || self.score < params.extension_threshhold {
                continue;
            }
            if  self.extend_one_left_unchecked(scheme, query, max_gaps != 0) {
                max_gaps -= 1;
            }
            if self.score > max_score {
                max_score = self.score;
                best_extension = i;
            }
            extended += 1;
        }
        self.word = self.word.split_off(extended - best_extension);
        self.padding_left += extended - best_extension;
        self.score = max_score;
        self.is_extended_left = true;
        self.idx_in_query += extended - best_extension;
    }

    fn extend_right(&mut self, scheme: &Arc<ScoringScheme>, query: &Arc<SimpleRecord>, params: &Arc<Params>, mut max_gaps: usize, mut max_extension: usize) {
        let mut  max_score = self.score;
        let mut best_extension: usize = 0;
        let mut extended = 0;
        max_extension = max_extension.min(self.padding_right);
        for i in 0..max_extension {
            if self.idx_in_query + self.word.len() + 1 == query.seq.len() || self.padding_right == 0 || self.score < params.extension_threshhold {
                continue;
            }
            if self.extend_one_right_unchecked(scheme, query, max_gaps != 0) {
                max_gaps -= 1;
            }
            if self.score > max_score {
                max_score = self.score;
                best_extension = i;
            }
            extended += 1;
        }
        let _ = self.word.split_off(self.word.len() - (extended - best_extension));
        self.padding_right += extended - best_extension;
        self.score = max_score;
        self.is_extended_right = true;
    }

    fn extend_one_right_unchecked(&mut self, scheme: &Arc<ScoringScheme>, query: &Arc<SimpleRecord>, can_use_gaps: bool) -> bool {
        // extends the sequence by one nt to the left and returns a flag, which indicates a gap
         if query.seq[self.idx_in_query + self.word.len() + 1] == 'N' as u8 {
            self.word.push('N' as u8);
            // save the information of a masked region, as that influences stats.
            self.padding_right -= 1;
            return false;
        }
        if self.db_seq[self.padding_left + self.word.len()] == query.seq[self.idx_in_query + self.word.len() + 1] {
            self.score += scheme.hit;
            self.word.push(query.seq[self.idx_in_query + self.word.len() + 1]);
            self.padding_right -= 1;
            false
        }
        else {
            if !can_use_gaps || (self.padding_right < 2 || (query.seq.len() - self.idx_in_query - self.word.len()) < 3) || query.seq[self.idx_in_query + self.word.len() + 2] == self.db_seq[self.padding_left + self.word.len() + 1] {
                self.word.push(query.seq[self.idx_in_query + self.word.len() + 1]);
                self.score += scheme.miss;
                self.padding_right -= 1;
                return false;
            }
            self.word.push( '-' as u8);
            if self.word[self.word.len() - 2] == '-' as u8 {
                self.score += scheme.gap_extension;
            }
            else {
                self.score += scheme.gap_penalty;
            }
            self.padding_right -= 1;
            true
        }
    }

    fn extend_one_left_unchecked(&mut self, scheme: &Arc<ScoringScheme>, query: &Arc<SimpleRecord>, can_use_gaps: bool) -> bool {
        if query.seq[self.idx_in_query - 1] == 'N' as u8 {
            self.word.insert(0, 'N' as u8);
            // save the information of a masked region, as that influences stats.
            self.padding_left -= 1;
            self.idx_in_query -= 1;
            return false;
        }
        if self.db_seq[self.padding_left - 1] == query.seq[self.idx_in_query - 1] {
            self.score += scheme.hit;
            self.word.insert(0, query.seq[self.idx_in_query - 1]);
            self.padding_left -= 1;
            self.idx_in_query -= 1;
            false
        }
        else {
            if !can_use_gaps || (self.padding_left < 2 || self.idx_in_query < 2) || query.seq[self.idx_in_query - 2] == self.db_seq[self.padding_left - 2] {
                self.word.insert(0, query.seq[self.idx_in_query - 1]);
                self.score += scheme.miss;
                self.padding_left -= 1;
                self.idx_in_query -= 1;
                return false;
            }
            self.word.insert(0, '-' as u8);
            if self.word[1] == '-' as u8 {
                self.score += scheme.gap_extension;
            }
            else {
                self.score += scheme.gap_penalty;
            }
            self.padding_left -= 1;
            self.idx_in_query -= 1;
            true
        }
    }
}

impl fmt::Display for HSP {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Hit in db {}", self.id)?;
        writeln!(f, "score: {}", self.score)?;
        writeln!(f, "E value: {:e}", self.e_val)?;
        writeln!(f, "start_idx query: {}", self.idx_in_query)?;
        writeln!(f, "start_idx rec: {}", self.idx_in_record + self.padding_left)?;
        writeln!(f, "length: {}", self.word.len())?;
        writeln!(f, "joined: {}, extended: l{} r{}", self.is_joined, self.is_extended_left, self.is_extended_right)?;
        let chars_per_line = 50;
        let max_chars = 200;
        for j in (0..self.word.len().min(max_chars)).step_by(chars_per_line) {
            write!(f, "Query: ")?;
            for b in &self.word[j..self.word.len().min(j + chars_per_line)] {
                write!(f, "{}", convert_to_ascii(b))?;
            }
            writeln!(f, "")?;
            write!(f, "       ")?;
            for i in j..self.word.len().min(j + chars_per_line) {
                if self.word[i] == self.db_seq[self.padding_left + i] {
                    write!(f, "|")?;
                }
                else {
                    write!(f, " ")?;
                }
            }
            writeln!(f, "")?;
            write!(f, "DB:    ")?;
            for b in &self.db_seq[self.padding_left + j..(self.padding_left + self.word.len()).min(j + self.padding_left + chars_per_line)] {
                write!(f, "{}", convert_to_ascii(b))?;
            }
            writeln!(f, "")?;
        }
        Ok(())
    }
}

#[derive (Default)]
struct ScoringScheme {
    gap_penalty: i16,
    gap_extension: i16,
    hit: i16,
    miss: i16,
    lambda: f64,
    k: f64
}

pub fn align<'a>(path_to_db: &'a str, path_to_query: &str, num_workers: usize, params: Params) -> io::Result<()> {
    let total = Instant::now();
    let parser = Instant::now();
    let mut query = parse_small_fasta(path_to_query)?;

    query.seq = Dust::new(64, params.masking_threshhold, query.seq).mask_regions();

    let params = Arc::new(params);
    get_query_words(params.k, &mut query);

    let mut records = read_csv(&path_to_db)?;
    let ids = records.by_ref().map(|r| r.id.clone()).collect::<Vec<String>>();
    records.reset();
    let lengths = records.by_ref().map(|r| (r.end_byte as usize - 1) - (r.start_byte as usize - 1) + (4 - r.start_bit) + r.end_bit).collect::<Vec<usize>>();
    records.reset();
    let info = Info {
        id: ids,
        length: lengths
    };
    let (raw_chunk_sender, raw_chunk_receiver) = mpsc::channel::<Vec<u8>>();
    let p = path_to_db.to_string();
    let reader: JoinHandle<io::Result<()>> = thread::spawn(move || {
        parse_compressed_db_lazy(&p, 2048, raw_chunk_sender)?;
        Ok(())
    });

    let scheme = Arc::new(ScoringScheme {
        gap_penalty: -8,
        gap_extension: -6,
        hit: 5,
        miss: -4,
        lambda: 0.146,
        k: 0.039
    });
    let worker = Instant::now();
    let mut worker_senders = Vec::default();
    let mut workers = Vec::default();
    let (h_tx, h_rx) = mpsc::channel::<Vec<HSP>>();
    let query = Arc::new(query);
    for _ in 0..num_workers {
        let (tx, rx) = mpsc::channel::<ProcessedChunk>();
        let params_ = Arc::clone(&params);
        let query = Arc::clone(&query);
        let scheme = Arc::clone(&scheme);
        let h_tx = h_tx.clone();
        workers.push(thread::spawn(move || {
            while let Ok(chunk) = rx.recv() {
                let _ = scan(chunk, Arc::clone(&params_), Arc::clone(&query), Arc::clone(&scheme), &h_tx);
            }
        }));
        worker_senders.push(tx);
    }
    let distributor_t = Instant::now();
    let params_ = Arc::clone(&params);
    let distributor = thread::spawn(move || {
        distribute_chunks(raw_chunk_receiver, worker_senders, params_, records);
    });
    let params_ = Arc::clone(&params);
    let process = Instant::now();
    let hit_proccessor = thread::spawn(move || {
        process_hits(h_rx, query, scheme, params_, info);
    });

    let _ = reader.join();
    let _ = distributor.join();
    for handle in workers {
        let _ = handle.join();
    }
    drop(h_tx);
    let _ = hit_proccessor.join();
    if params.verbose {
        println!("parser: {:?}", parser.elapsed());
        println!("distributor: {:?}", distributor_t.elapsed());
        println!("workers: {:?}", worker.elapsed());
        println!("process: {:?}", process.elapsed());
        println!("total: {:?}", total.elapsed());
    }
    Ok(())
}

fn process_hits(rx: mpsc::Receiver<Vec<HSP>>, query: Arc<SimpleRecord>, scheme: Arc<ScoringScheme>, params: Arc<Params>, infos: Info) {
    // join close hits and extend
    let mut h = 0;
    let mut all_hits: Vec<Rc<RefCell<HSP>>> = Vec::default();
    while let Ok(hits) = rx.recv() {
        h += hits.len();
        all_hits.extend(hits.into_iter().map(|h| Rc::new(RefCell::new(h))));
        join_hits(&all_hits, &scheme, &params, &query, 5);
    }
    for h in all_hits.iter() {
        if h.borrow().is_joined {
            continue;
        }
        if !h.borrow().is_extended_left {
            h.borrow_mut().extend_left(&scheme, &query, &params, 20, 50);
        }
        if !h.borrow().is_extended_right {
            h.borrow_mut().extend_right(&scheme, &query, &params, 20, 50);
        }
    }
    let _ = all_hits.iter().map(|h|  {
        let s = h.borrow().score;
        let bit_s = (scheme.lambda * s as f64 - scheme.k.ln()) / 2.0_f64.ln();
        let e_val = (h.borrow().word.len() * infos.get_length_of(&h.borrow().id)) as f64 * 2.0_f64.powf(-bit_s);
        h.borrow_mut().e_val = e_val;
    }).collect::<Vec<()>>();
    all_hits.sort_by(|a, b| b.borrow().e_val.partial_cmp(&a.borrow().e_val).unwrap());
    println!("{} hits found, {} non joined hits", h, all_hits.iter().map(|i| if i.borrow().is_joined {0} else {1}).sum::<i32>());
    println!("\nBest hit: ");
    print!("\n{}\n", all_hits.last().unwrap_or(&Rc::new(RefCell::new(HSP::default()))).borrow());
}

fn join_hits(hits: &[Rc<RefCell<HSP>>], scheme: &Arc<ScoringScheme>, params: &Arc<Params>, query: &Arc<SimpleRecord>, max_distance: usize) {
    //using a greedy recursive approach
    //TODO: optimize this!! way too slow for many hits
    for i in 0..hits.len() - 1 {
        if hits[i].borrow().is_joined {
            continue;
        }
        let h = Rc::clone(&hits[i]);
        join_left(&hits[..i], h, scheme, params, query, max_distance);
        let h = Rc::clone(&hits[i]);
        join_right(&hits[i + 1..], h, scheme, params, query, max_distance);
    }
}

fn join_left(hits: &[Rc<RefCell<HSP>>], subject: Rc<RefCell<HSP>>, scheme: &Arc<ScoringScheme>, params: &Arc<Params>, query: &Arc<SimpleRecord>, max_distance: usize) {
    // recursively joins the given HSP with HSPs to its left, while possible
    if hits.len() == 0 {
        return;
    }
    let h: Rc<RefCell<HSP>>;
    if hits[hits.len() - 1].borrow_mut().try_join(&mut subject.borrow_mut(), scheme, params, query, max_distance) {
        h = Rc::clone(&hits[hits.len() - 1]);
    }
    else {
        h = Rc::clone(&subject);
    }
    join_left(&hits[..hits.len() - 1], h, scheme, params, query, max_distance);
}

fn join_right(hits: &[Rc<RefCell<HSP>>], subject: Rc<RefCell<HSP>>, scheme: &Arc<ScoringScheme>, params: &Arc<Params>, query: &Arc<SimpleRecord>, max_distance: usize) {
    // recursively joins the given HSP with HSPs to its right, while possible
    if hits.len() == 0 {
        return;
    }
    let h: Rc<RefCell<HSP>>;
    if subject.borrow_mut().try_join(&mut hits[0].borrow_mut(), scheme, params, query, max_distance) {
        h = Rc::clone(&hits[0]);
    }
    else {
        h = Rc::clone(&subject);
    }
    join_right(&hits[1..], h, scheme, params, query, max_distance);
}

fn distribute_chunks(chunk_receiver: mpsc::Receiver<Vec<u8>>, distributors: Vec<mpsc::Sender<ProcessedChunk>>, params: Arc<Params>, mut records: impl Iterator<Item = Arc<Record>>) {
    let mut i = 0;
    let mut total_bytes = 0;
    let mut overlap_buf: Vec<u8> = vec![];
    let mut current_record = records.next().unwrap();
    let mut bytes_in_current_record = 0;
    let mut not_finished = true;

    while not_finished {
        let received: Vec<u8>;
        match chunk_receiver.recv() {
            Ok(bytes) => received = bytes,
            Err(_e) => {
                received = vec![];
                not_finished = false;
            },
        }
        let bytes = received.len();
        total_bytes += bytes as u128;
        let seq: Vec<u8> = overlap_buf.iter().copied().chain(received.into_iter()).collect();
        let chunk = process_chunk(seq, &mut records, &mut current_record, &total_bytes, &params, &mut overlap_buf, bytes, &mut bytes_in_current_record);
        distributors[i % distributors.len()].send(chunk).expect(&format!("failed to send chunk {} to worker thread {}", i, i % distributors.len()));
        i += 1;
    }
}

fn process_chunk(seq: Vec<u8>, records: &mut impl Iterator<Item = Arc<Record>>, current_record: &mut Arc<Record>, total_bytes: &u128, params: &Params, overlap_buf: &mut Vec<u8>, received: usize, bytes_in_rec: &mut u128) -> ProcessedChunk {
    // generate a padded chunk and save the next padding in the buffer. chunk[extension_length..chunk.len() - extension_length], padded by at most extension_length in both directions, buffer = chunk[chunk.len() - word_size - 2 * extension_len..] if no new record else it starts at record boundary
    let chunk: ProcessedChunk;
    let start: usize;
    let start_in_rec: u128;
    let start_byte: usize;
    if seq.len() == 0  {
        return ProcessedChunk {
            bytes: vec![0],
            end_bit: 0,
            end_byte: 0, 
            start_bit: 0,
            id: "".to_string(),
            start_byte: 0,
            start_in_rec: 0
        };
    }
    if overlap_buf.len() != params.k.div_ceil(4) + 2 * params.extension_length.div_ceil(4) || overlap_buf.len() == 0 || seq.len() == overlap_buf.len() {
        start = current_record.start_bit;
        start_in_rec = 0;
        start_byte = 0;
    }
    else {
        start = 0;
        start_in_rec = *bytes_in_rec;
        start_byte = params.extension_length.div_ceil(4);
    }
    if current_record.end_byte < *total_bytes {
        let split_idx = seq.len() - (*total_bytes - current_record.end_byte) as usize;
        if current_record.end_bit == 3 {
            *overlap_buf = seq[split_idx + 1..].to_vec();
        }
        else {
            *overlap_buf = seq[split_idx..].to_vec();
        }
        chunk = ProcessedChunk {
            bytes: seq,
            id: current_record.id.clone(),
            end_bit: current_record.end_bit,
            end_byte: split_idx,
            start_in_rec: start_in_rec,
            start_bit: start,
            start_byte: split_idx.min(start_byte)
        };
        *bytes_in_rec = (received - split_idx) as u128;
        if let Some(rec) = records.next() {
            *current_record = rec;  
        }
    }
    else {
        *overlap_buf = seq[seq.len() - params.extension_length.div_ceil(4) * 2 - params.k.div_ceil(4)..].to_vec();
        chunk = ProcessedChunk {
            end_byte: seq.len() - params.extension_length.div_ceil(4) - 1,
            bytes: seq,
            id: current_record.id.clone(),
            end_bit: 3,
            start_in_rec: start_in_rec,
            start_bit: start as usize,
            start_byte: start_byte
        };
        *bytes_in_rec += received as u128;
    }
    chunk
}

fn _split_byte(byte: u8, split_idx: usize) -> (u8, u8) {
    // split the byte into two at the given index
    let (mut byte1, mut byte2) = (0b00u8, 0b00u8);
    for i in 0..split_idx {
        byte2 |= byte2 ^ (byte >> 2 * i & 0b11) << 6 - 2 * i;
    }
    for i in 0..4-split_idx {
        byte1 |= byte1 ^ (byte >> 2 * i & 0b11) << 6 - 2 * i;
    }
    (byte1, byte2)
}

fn scan(chunk: ProcessedChunk, params: Arc<Params>, query: Arc<SimpleRecord>, scheme: Arc<ScoringScheme>, hit_sender: &mpsc::Sender<Vec<HSP>>) -> Result<(), String> {
    let mut hits = Vec::default();
    for (i, slice1) in chunk.bytes[chunk.start_byte..=chunk.end_byte].windows(3).enumerate() {
        for (j, slice2) in query.words.iter().enumerate() {
            if !(slice1[0] == slice2[0] && slice1[1] == slice2[1]) {
                continue;
            }
            let mut score = scheme.hit * 8;
            score += get_bitwise_score(&slice1[2], &slice2[2], &scheme, 4);
            if score < params.scanning_threshhold {
                continue;
            } 
            //TODO: maybe calculate E-value and use that
            // send a hit   
            //TODO: check the range start in the following line
            let seq = bytes_to_chars(&chunk.bytes[0.max(i + chunk.start_byte - (params.extension_length.div_ceil(4)).min(i + chunk.start_byte))..chunk.bytes.len().min(i + params.extension_length.div_ceil(4) + params.k.div_ceil(4) + chunk.start_byte)], if chunk.bytes.len().min(i + params.extension_length.div_ceil(4) + params.k.div_ceil(4) + chunk.start_byte) == chunk.bytes.len() {chunk.end_bit} else {3}, chunk.start_bit);
            hits.push(HSP {
                id: chunk.id.clone(),
                score: score,
                padding_left: params.extension_length.min((chunk.start_byte + i) * 4),
                padding_right: params.extension_length.min((chunk.bytes.len() - i - params.k.div_ceil(4)) * 4 - 3 + chunk.end_bit),
                db_seq: seq,
                word: bytes_to_chars(&query.words[j], 3, 0),
                idx_in_query: j,
                idx_in_record: chunk.start_in_rec as usize + i * 4 + chunk.start_byte * 4 - params.extension_length.min((chunk.start_byte + i) * 4),
                is_extended_right: false,
                is_extended_left: false,
                is_joined: false,
                e_val: 0.0 //TODO
                });
        }
    }
    if hits.len() > 0 {
        if let Err(_e) = hit_sender.send(hits) {
            //eprintln!("could not send hit info, {:#?}", e);
        }
    }
    Ok(())
}

fn _get_score_extension(_s1: &[u8], _s2: &[u8], _scheme: &Arc<ScoringScheme>, _new_gap: bool) -> i16 {
    0
}

fn get_score(s1: &[u8], s2: &[u8], scheme: &Arc<ScoringScheme>) -> i16 {
    (0..s1.len()).map(|i| get_bitwise_score(&s1[i], &s2[i], scheme, 4)).sum()
}

fn get_bitwise_score(byte1: &u8, byte2: &u8, scheme: &Arc<ScoringScheme>, end_bit: usize) -> i16 {
    (0..4).map(|i| {
            if i >= end_bit {
                0
            }
            else if ((byte1 & 0b11 << i * 2) ^ (byte2 & 0b11 << i * 2)) == 0b00 {
                scheme.hit
            }
            else {
                scheme.miss
            }
        }).sum()
}

fn get_query_words(k: usize, query: &mut SimpleRecord) {
    //TODO: save the correct cutoffs of the bits per word/seq
    query.words = (0..query.seq.len() - k + 1)
        .into_iter()
        .map(|i| {
            if !query.seq[i..i+k].contains(&('N' as u8)) {
                parse_to_bytes(&query.seq[i..i+k])
            }
            else {
                Vec::default()
            }
        })
        .filter(|word| {
            word.len() == k.div_ceil(4)
        })
        .collect();
    query.k = k;
}

