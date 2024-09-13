use std::{io, ops::Div, sync::{mpsc, Arc}, thread::{self, JoinHandle}, time::Instant};

use crate::make_db::{parse_fasta::{parse_small_fasta, parse_to_bytes}, read_db::{extract_str_from_bytes, parse_compressed_db_lazy, read_csv}, records::{Record, SimpleRecord}};//, extract_str_from_bytes}};

//TODO: precalculated all div_ceils, where possible

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

#[derive (Default, Debug)]
struct Hit {
    id: String,
    score: Vec<i16>,
    left: Vec<usize>, // left position of the hit in the record
    right: Vec<usize>, // right position of the hit in the record
    db_seq: Vec<Vec<u8>>,
    word: Vec<Vec<u8>>
}

impl Hit {
    fn join_overlapping(&mut self, other: Hit, scheme: &Arc<ScoringScheme>) {
        let overlap = self.right.last().unwrap().abs_diff(other.left[0]);
        self.word.last_mut().unwrap().extend_from_slice(&other.word[0][overlap..]);
        *self.right.last_mut().unwrap() += other.word[0].len() - overlap;
        *self.score.last_mut().unwrap() += get_score(&other.word[0][overlap..], &other.db_seq[0][overlap..], &scheme);
        self.word.extend_from_slice(&other.word[1..]);
        self.db_seq.extend_from_slice(&other.db_seq[1..]);
        self.left.extend_from_slice(&other.left[1..]);
        self.right.extend_from_slice(&other.right[1..]);
        self.score.extend_from_slice(&other.score[1..]);
    }

    fn join_separated(&mut self, mut other: Hit) {
        self.score.append(&mut other.score);
        self.left.append(&mut other.left);
        self.right.append(&mut other.right);
        self.db_seq.append(&mut other.db_seq);
        self.word.append(&mut other.word);
    }
}

#[derive (Default)]
struct ScoringScheme {
    gap_penalty: i16,
    hit: i16,
    miss: i16,
    lambda: f64,
    k: f64
}

pub fn align<'a>(path_to_db: &'a str, path_to_query: &str, num_workers: usize, mut params: Params) -> io::Result<()> {
    let parser = Instant::now();
    let mut query = parse_small_fasta(path_to_query)?;
    params.query_length = query.seq.len() / 4;
    if query.seq.len() % 4 != 0 {
        params.query_length += 1;
    }
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
        gap_penalty: -1,
        hit: 5,
        miss: -4,
        lambda: 1.0,
        k: 1.0
    });
    let worker = Instant::now();
    let mut worker_senders = Vec::default();
    let mut workers = Vec::default();
    let (h_tx, h_rx) = mpsc::channel::<Vec<Hit>>();
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
    let params = Arc::clone(&params);
    let distributor = thread::spawn(move || {
        distribute_chunks(raw_chunk_receiver, worker_senders, params, records);
    });
    println!("distributor: {:?}", distributor_t.elapsed());
    let hit_proccessor = thread::spawn(move || {
        process_hits(h_rx, query, scheme);
    });

    let _ = reader.join();
    println!("parser: {:?}", parser.elapsed());
    let _ = distributor.join();
    for handle in workers {
        let _ = handle.join();
    }
    println!("workers: {:?}", worker.elapsed());
    drop(h_tx);
    let _ = hit_proccessor.join();
    Ok(())
}

fn process_hits(rx: mpsc::Receiver<Vec<Hit>>, query: Arc<SimpleRecord>, scheme: Arc<ScoringScheme>) {
    // join close hits and extend
    let mut h = 0;
    let mut all_hits = Vec::default();
    while let Ok(mut hits) = rx.recv() {
        h += hits.len();
        all_hits.append(&mut hits);
        all_hits.sort_by(|a, b| a.left.cmp(&b.left)); // this can be optimized easily, as hits is already sorted and all_hits is sorted from last iteration
        join_hits(&mut all_hits, 2, &scheme);
        extend_hits(&mut all_hits, &query, &scheme);
    }
    println!("{} hits found", h);
    println!("len: {}", all_hits.len());
}

fn join_hits(hits: &mut Vec<Hit>, max_distance: usize, scheme: &Arc<ScoringScheme>) {
    let mut k = 0;
    for i in 0..hits.len() - 1 {
        if hits[i - k].id != hits[i - k + 1].id {
            continue;
        }
        if *hits[i - k].right.last().unwrap() >= hits[i - k + 1].left[0] {
            let other = hits.remove(i - k + 1);
            hits[i - k].join_overlapping(other, &scheme);
            k += 1;
        }
        else if hits[i - k].right.last().unwrap().abs_diff(hits[i - k + 1].left[0]) <= max_distance {
            let other = hits.remove(i - k + 1);
            hits[i - k].join_separated(other);
            k += 1;
        }
    }
}

fn extend_hits(hits: &mut Vec<Hit>, _query: &Arc<SimpleRecord>, _scheme: &Arc<ScoringScheme>) {
    for _hit in hits {

    }
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
    // generate a padded chunk and save the next padding in the buffer. chunk[extension_lenght..chunk.len() - extension_length], buffer = chunk[chunk.len() - word_size - extension_len..] if no new record else it starts at record boundary
    let chunk: ProcessedChunk;
    let start: usize;
    let start_in_rec: u128;
    let start_byte: usize;
    if !overlap_buf.len() == params.k.div_ceil(4) + params.extension_length.div_ceil(4) || overlap_buf.len() == 0 {
        start = current_record.start_bit;
        start_in_rec = 0;
        start_byte = 0;
    }
    else {
        start = 0;
        start_in_rec = *bytes_in_rec;
        start_byte = params.extension_length.div_ceil(4);
    }
    if current_record.end_byte <= *total_bytes {
        let split_idx = seq.len() - (*total_bytes - current_record.end_byte) as usize;
        *overlap_buf = seq[split_idx..].to_vec();
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
        *overlap_buf = seq[0.max(seq.len().abs_diff(params.extension_length.div_ceil(4) * 2 + params.k.div_ceil(4)))..].to_vec();
        chunk = ProcessedChunk {
            end_byte: start_byte.max(seq.len() - params.extension_length.div_ceil(4) - 1),
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

fn scan(chunk: ProcessedChunk, params: Arc<Params>, query: Arc<SimpleRecord>, scheme: Arc<ScoringScheme>, hit_sender: &mpsc::Sender<Vec<Hit>>) -> Result<(), String> {
    let mut hits = Vec::default();
    for (i, slice1) in chunk.bytes[chunk.start_byte..=chunk.end_byte].windows(3).enumerate() {
        for (j, slice2) in query.words.iter().enumerate() {
            if !(slice1[0] == slice2[0] && slice1[1] == slice2[1]) {
                continue;
            }
            let mut score = scheme.hit * 8;
            score += get_bitwise_score(&slice1[2], &slice2[2], &scheme);
            if score < params.scanning_threshhold {
                continue;
            }
            // send a hit   
            hits.push(Hit {
                id: chunk.id.clone(),
                score: vec![score],
                left: vec![i + chunk.start_in_rec as usize],
                right: vec![i + params.k.div_ceil(4) + chunk.start_in_rec as usize],
                db_seq: vec![chunk.bytes[0.max(i - params.extension_length.div_ceil(4).min(i))..chunk.bytes.len().min(i + params.extension_length.div_ceil(4) + params.k.div_ceil(4))].to_vec()],
                word: vec![query.words[j].clone()]
            });
        }
    }
    if hits.len() > 0 {
        if let Err(e) = hit_sender.send(hits) {
            println!("could not send hit info, {:#?}", e);
        }
    }
    Ok(())
}

fn _extend(_db_seq: &[u8], _query_seq: &[u8]) -> i16 {
    0
}

fn get_score(s1: &[u8], s2: &[u8], scheme: &Arc<ScoringScheme>) -> i16 {
    (0..s1.len()).map(|i| get_bitwise_score(&s1[i], &s2[i], scheme)).sum()
}

fn get_bitwise_score(byte1: &u8, byte2: &u8, scheme: &Arc<ScoringScheme>) -> i16 {
    (0..4).map(|i| {
            if ((byte1 & 0b11 << i * 2) ^ (byte2 & 0b11 << i * 2)) == 0b00 {
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
            parse_to_bytes(&query.seq[i..i+k])
        })
        .collect();
    query.k = k;
    query.seq = parse_to_bytes(&query.seq);
}

