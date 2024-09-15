use std::{io, sync::{mpsc, Arc}, thread::{self, JoinHandle}, time::Instant};

use crate::{make_db::{parse_fasta::{parse_small_fasta, parse_to_bytes}, read_db::{bytes_to_chars, parse_compressed_db_lazy, read_csv}, records::{Record, SimpleRecord}}, dust::Dust};

//TODO: -precalculate all div_ceils, where possible
//      -implement support for word lengths k where k % 4 != 0

pub struct Params {
    pub k: usize,
    pub extension_threshhold: i16,
    pub scanning_threshhold: i16,
    pub extension_length: usize,
    pub query_length: usize
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

#[derive (Debug, Default)]
struct JoinedHits {
    hits: Vec<Hit>,
    id: String
}

//TODO: rewrite
impl JoinedHits {
    fn join_overlapping(&mut self, other: JoinedHits, scheme: &Arc<ScoringScheme>, params: &Arc<Params>) {
        self.hits.last_mut().unwrap().join(&other.hits[0], scheme, params);
        self.hits.extend_from_slice(&other.hits[1..]);
    }

    fn join_separated(&mut self, mut other: JoinedHits) {
        self.hits.append(&mut other.hits);
    }

    fn get_end(&self) -> usize {
        self.hits.last().unwrap().right
    }

    fn get_start(&self) -> usize {
        self.hits[0].left
    }

    fn join_all(&mut self, scheme: &Arc<ScoringScheme>, params: &Arc<Params>) {
        let mut k = 0;
        for i in 0..self.hits.len() - 1 {
            if self.hits[i - k].word.len() >= params.query_length - params.k || self.hits[i - k + 1].word.len() >= params.query_length - params.k {
                continue;
            }
            if self.hits[i - k].right >= self.hits[i - k + 1].left {
                let hit = self.hits.remove(i - k + 1);
                self.hits[i - k].join(&hit, scheme, params);
                k += 1;
            }
        }
    }
}
#[derive (Default, Debug, Clone)]
struct Hit {
    score: i16,
    left: usize, // left position of the hit in the record
    right: usize, // right position of the hit in the record
    db_seq: Vec<u8>,
    word: Vec<u8>,
    is_extended: bool,
    position_in_query: usize,
    extension_length: usize
}

impl Hit {
    fn join(&mut self, other: &Hit, scheme: &Arc<ScoringScheme>, _params: &Arc<Params>) {
        if true {
            return;
        }
        let overlap = self.right.abs_diff(other.left);
        self.right = other.right;
        let current_score = _get_score_extension(&self.word[&self.word.len() - overlap..], &self.db_seq[self.extension_length + &self.db_seq.len() - overlap..], scheme, true);
        let new_score = _get_score_extension(&other.word[..overlap], &other.db_seq[other.extension_length..overlap], scheme, true);
        if new_score <= current_score {
            self.word.extend_from_slice(&other.word[overlap..]);
            self.score += _get_score_extension(&other.word[overlap..], &other.db_seq[overlap + other.extension_length..], scheme, true);
        }
        else {
            self.word.splice(self.word.len() - overlap.., other.word.clone());
            self.score += _get_score_extension(&other.word[overlap..], &other.db_seq[overlap + other.extension_length..], scheme, true) + new_score;
        }
    }
}

#[derive (Default)]
struct ScoringScheme {
    gap_penalty: i16,
    gap_extension: i16,
    hit: i16,
    miss: i16,
    _lambda: f64,
    _k: f64
}

pub fn align<'a>(path_to_db: &'a str, path_to_query: &str, num_workers: usize, mut params: Params) -> io::Result<()> {
    let parser = Instant::now();
    let mut query = parse_small_fasta(path_to_query)?;
    params.query_length = query.seq.len();

    query.seq = Dust::new(64, 1.0, query.seq).mask_regions();

    let params = Arc::new(params);
    get_query_words(params.k, &mut query);

    let records = read_csv(&path_to_db)?.into_iter();
    let (raw_chunk_sender, raw_chunk_receiver) = mpsc::channel::<Vec<u8>>();
    let p = path_to_db.to_string();
    let reader: JoinHandle<io::Result<()>> = thread::spawn(move || {
        parse_compressed_db_lazy(&p, 2048, raw_chunk_sender)?;
        Ok(())
    });

    let scheme = Arc::new(ScoringScheme {
        gap_penalty: -5,
        gap_extension: -1,
        hit: 5,
        miss: -4,
        _lambda: 1.0,
        _k: 1.0
    });
    let worker = Instant::now();
    let mut worker_senders = Vec::default();
    let mut workers = Vec::default();
    let (h_tx, h_rx) = mpsc::channel::<Vec<JoinedHits>>();
    let query = Arc::new(query);
    for _ in 0..num_workers {
        let (tx, rx) = mpsc::channel::<ProcessedChunk>();
        let params = Arc::clone(&params);
        let query = Arc::clone(&query);
        let scheme = Arc::clone(&scheme);
        let h_tx = h_tx.clone();
        workers.push(thread::spawn(move || {
            while let Ok(chunk) = rx.recv() {
                let _ = scan(chunk, Arc::clone(&params), Arc::clone(&query), Arc::clone(&scheme), &h_tx);
            }
        }));
        worker_senders.push(tx);
    }
    let distributor_t = Instant::now();
    let params_ = Arc::clone(&params);
    let distributor = thread::spawn(move || {
        distribute_chunks(raw_chunk_receiver, worker_senders, params_, records);
    });

    let hit_proccessor = thread::spawn(move || {
        process_hits(h_rx, query, scheme, params);
    });

    let _ = reader.join();
    println!("parser: {:?}", parser.elapsed());
    let _ = distributor.join();
    println!("distributor: {:?}", distributor_t.elapsed());
    for handle in workers {
        let _ = handle.join();
    }
    println!("workers: {:?}", worker.elapsed());
    drop(h_tx);
    let _ = hit_proccessor.join();
    Ok(())
}

fn process_hits(rx: mpsc::Receiver<Vec<JoinedHits>>, _query: Arc<SimpleRecord>, scheme: Arc<ScoringScheme>, params: Arc<Params>) {
    // join close hits and extend
    let mut h = 0;
    while let Ok(mut hits) = rx.recv() {
        h += hits.len();
        join_hits(&mut hits, 4, &scheme, &params);
        for _hit in hits {
            //println!("{}, {}, {}, {}", hit.hits[0].extension_length, hit.hits[0].left,  hit.hits[0].right,  hit.hits[0].position_in_query);
            //println!("{:?}", hit.hits[0].db_seq.iter().map(|i| convert_to_ascii(i)).collect::<String>());
            //println!("{:?}", hit.hits[0].word.iter().map(|i| convert_to_ascii(i)).collect::<String>());
        }
        /*
        h += hits.hits.len();
        hits.join_all(&scheme, &params);
        all_hits.push(hits);
        all_hits.sort_by(|a, b| a.hits[0].left.cmp(&b.hits[0].left)); // this can be optimized easily, as hits is already sorted and all_hits is sorted from last iteration
        join_hits(&mut all_hits, 5, &scheme, &params);
        extend_hits(&mut all_hits, &query, &scheme, &params);
        */
    }
    println!("{} hits found", h);
}

fn join_hits(hits: &mut Vec<JoinedHits>, max_distance: usize, scheme: &Arc<ScoringScheme>, params: &Arc<Params>) {
    let mut k = 0;
    hits[0].join_all(scheme, params);
    for i in 0..hits.len() - 1 {
        hits[i - k + 1].join_all(scheme, params);
        if hits[i - k].id != hits[i - k + 1].id || hits[i - k].hits.last().unwrap().word.len() >= params.query_length - params.k || hits[i - k + 1].hits[0].word.len() >= params.query_length - params.k {
            continue;
        }
        if hits[i - k].get_end() >= hits[i - k + 1].get_start(){
            let other = hits.remove(i - k + 1);
            hits[i - k].join_overlapping(other, scheme, params);
            k += 1;
        }
        else if hits[i - k].get_end().abs_diff(hits[i - k + 1].get_start()) <= max_distance {
            let other = hits.remove(i - k + 1);
            hits[i - k].join_separated(other);
            k += 1;
        }
    }
}

//TODO: rewrite extend
fn _extend_hits(hits: &mut Vec<JoinedHits>, query: &Arc<SimpleRecord>, scheme: &Arc<ScoringScheme>, params: &Arc<Params>) {
    for hit in hits {
        for i in 0..hit.hits.len() {
            if hit.hits[i].is_extended {
                continue;
            }
            let to_extend_right: usize;
            let contains_gaps_right: usize;
            let to_extend_left: usize;
            let contains_gaps_left: usize;
            if hit.hits.len() > i + 1 {
                to_extend_right = hit.hits[i + 1].left.abs_diff(hit.hits[i].right);
                contains_gaps_right = hit.hits[i + 1].left.abs_diff(hit.hits[i].right).abs_diff(hit.hits[i + 1].position_in_query.abs_diff(hit.hits[i].position_in_query));
            }
            else {
                to_extend_right = 1.max(hit.hits[i].extension_length) - 1;
                contains_gaps_right = hit.hits[i].extension_length - hit.hits[i].position_in_query + hit.hits[i].word.len();
            }
            if i >= 1 {
                to_extend_left = hit.hits[i].left.abs_diff(hit.hits[i - 1].right);
                contains_gaps_left = hit.hits[i].left.abs_diff(hit.hits[i - 1].right).abs_diff(hit.hits[i].position_in_query.abs_diff(hit.hits[i - 1].position_in_query));
            }
            else {
                to_extend_left = hit.hits[i].extension_length - 1;
                contains_gaps_left = hit.hits[i].extension_length - hit.hits[i].position_in_query;
            }
            _extend_left(&mut hit.hits[i], query, params, scheme, to_extend_left, contains_gaps_left);
            _extend_right(&mut hit.hits[i], query, params, scheme, to_extend_right, contains_gaps_right);
            hit.hits[i].is_extended = true;
        }
        hit.join_all(scheme, params);
    }
}

fn _extend_right(hit: &mut Hit, query: &Arc<SimpleRecord>, params: &Arc<Params>, scheme: &Arc<ScoringScheme>, to_extend: usize, mut n_gaps: usize) {
    let mut max = 0;
    let mut best_word = hit.word.clone();
    let mut buf: Vec<u8> = Vec::default();
    let mut extended = 0;
    let mut true_extended = 0;
    let mut in_gap: bool;
    while hit.score >= params.extension_threshhold && extended < to_extend {
        if hit.db_seq[hit.extension_length + hit.word.len() + 1] == query.seq[hit.position_in_query + hit.word.len() + 1] || query.seq[hit.position_in_query + hit.word.len() + 1] == 'N' as u8 {
            hit.word.push(query.seq[hit.position_in_query + hit.word.len() + 1]);
            buf.push(query.seq[hit.position_in_query + hit.word.len() + 1]);
            extended += 1;
            in_gap = false;
        }
        else if n_gaps > 0 {
            hit.word.push('-' as u8);
            buf.push('-' as u8);
            extended += 1;
            in_gap = true;
            n_gaps -= 1;
        }
        else {
            hit.word.push(query.seq[hit.position_in_query + hit.word.len() + extended + 1]);
            buf.push(query.seq[hit.position_in_query + hit.word.len() + extended + 1]);
            extended += 1;
            in_gap = false;
        }

        hit.score += _get_score_extension(&hit.word[&hit.word.len() - 1..], &hit.db_seq[&hit.word.len() - 1 + hit.extension_length..], scheme, in_gap);
        if hit.score >= max {
            max = hit.score;
            true_extended += buf.len();
            best_word.append(&mut buf);
        }
    }
    hit.word = best_word;
    hit.score = max;
    hit.right += true_extended;
}

fn _extend_left(hit: &mut Hit, query: &Arc<SimpleRecord>, params: &Arc<Params>, scheme: &Arc<ScoringScheme>, to_extend: usize, mut n_gaps: usize) {
    let mut max = 0;
    let mut best_word = hit.word.clone();
    let mut buf: Vec<u8> = Vec::default();
    let mut extended = 0;
    let mut true_extended = 0;
    let mut in_gap: bool;
    while hit.score >= params.extension_threshhold && extended <= to_extend {
        if hit.db_seq[(extended + 1).max(hit.extension_length) - extended - 1] == query.seq[hit.position_in_query - extended - 1] || query.seq[hit.position_in_query - extended - 1] == 'N' as u8 {
            hit.word.push(query.seq[hit.position_in_query - extended - 1]);
            buf.push(query.seq[hit.position_in_query - extended - 1]);
            extended += 1;
            in_gap = false;
        }
        else if n_gaps > 0 {
            hit.word.insert(0,'-' as u8);
            buf.insert(0, '-' as u8);
            extended += 1;
            in_gap = true;
            n_gaps -= 1;
        }
        else {
            hit.word.insert(0, query.seq[hit.position_in_query - extended - 1]);
            buf.insert(0, query.seq[hit.position_in_query - extended - 1]);
            extended += 1;
            in_gap = false;
        }

        hit.score += _get_score_extension(&[hit.word[0]], &[hit.db_seq[extended.max(hit.extension_length) - extended]], scheme, in_gap);
        if hit.score >= max {
            max = hit.score;
            true_extended += buf.len();
            best_word.insert(0, 0);
            best_word.splice(0..1, buf.clone());
        }
    }
    hit.word = best_word;
    hit.score = max;
    hit.left -= true_extended;
    hit.position_in_query -= true_extended;
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
    if overlap_buf.len() != params.k.div_ceil(4) + 2 * params.extension_length.div_ceil(4) || overlap_buf.len() == 0 {
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

fn scan(chunk: ProcessedChunk, params: Arc<Params>, query: Arc<SimpleRecord>, scheme: Arc<ScoringScheme>, hit_sender: &mpsc::Sender<Vec<JoinedHits>>) -> Result<(), String> {
    let mut hits = Vec::default();
    //println!("b: {:?}, {}, {}, {}", extract_str_from_bytes(&chunk.bytes), chunk.end_byte, (0..=chunk.end_byte).last().unwrap(), chunk.end_bit);
    //let mut  s = bytes_to_chars(&[chunk.bytes[0]], 3, chunk.start_bit);
    //s.extend(&bytes_to_chars(&chunk.bytes[chunk.start_byte + 1..=chunk.end_byte], chunk.end_bit, 0));
    //println!("s: {:?}", s.iter().map(|i| convert_to_ascii(i)).collect::<String>());
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
            let seq = bytes_to_chars(&chunk.bytes[0.max(i - params.extension_length.div_ceil(4).min(i))..chunk.bytes.len().min(i + params.extension_length.div_ceil(4) + params.k.div_ceil(4))], chunk.end_bit, chunk.start_bit);
            hits.push(JoinedHits {
                hits: vec![Hit {
                score: score,
                extension_length: params.extension_length.min((chunk.start_in_rec as usize).min(seq.len() - chunk.start_in_rec as usize - chunk.end_byte)),
                left: i * 4 + (chunk.start_in_rec * 4) as usize,
                right: i * 4 + params.k + (chunk.start_in_rec * 4) as usize,
                db_seq: seq,
                word: bytes_to_chars(&query.words[j], 3, 0),
                is_extended: false,
                position_in_query: j,
                }],
                id: chunk.id.clone(),
            });
        }
    }
    if hits.len() > 0 {
        if let Err(_e) = hit_sender.send(hits) {
            //println!("could not send hit info, {:#?}", e);
        }
    }
    Ok(())
}

fn _get_score_extension(_s1: &[u8], _s2: &[u8], _scheme: &Arc<ScoringScheme>, _new_gap: bool) -> i16 {
    0
}

fn _get_score(s1: &[u8], s2: &[u8], scheme: &Arc<ScoringScheme>) -> i16 {
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
    //query.seq = parse_to_bytes(&query.seq);
}

