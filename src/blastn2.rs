use std::{io, sync::{mpsc, Arc}, thread::{self, JoinHandle}, time::Instant};

use crate::make_db::{records::{SimpleRecord, Record}, parse_fasta::{parse_small_fasta, parse_to_bytes}, read_db::{read_csv, parse_compressed_db_lazy, extract_str_from_bytes}};

pub struct Params {
    pub k: usize,
    pub extension_threshhold: i16,
    pub scanning_threshhold: i64,
    pub extension_length: usize,
    pub query_length: usize
}

#[derive (Default)]
struct ProcessedChunk {
    bytes: Vec<u8>,
    end_bit: usize,
    end_byte: usize,
    start_bit: usize, //what is the first relevant bit in the record
    id: String,
    start_in_rec: u128, // where is the first byte of the chunk releative to the record,
    start_byte: usize
}

#[derive (Default, Debug)]
struct Hit {
    id: String,
    score: i16,
    left: usize,
    right: usize
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
    let (h_tx, h_rx) = mpsc::channel::<Hit>();
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
    let distro = Instant::now();
    let params = Arc::clone(&params);
    let distributor = thread::spawn(move || {
        distribute_chunks(raw_chunk_receiver, worker_senders, params, records);
    });
    println!("distro: {:?}", distro.elapsed());
    let hit_proccessor = thread::spawn(move || {
        process_hits(h_rx);
    });

    let _ = reader.join();
    println!("parser: {:#?}", parser.elapsed());
    let _ = distributor.join();
    for handle in workers {
        let _ = handle.join();
    }
    println!("workers: {:?}", worker.elapsed());
    drop(h_tx);
    let _ = hit_proccessor.join();
    Ok(())
}

fn process_hits(rx: mpsc::Receiver<Hit>) {
    while let Ok(hit) = rx.recv() {
        println!("hit received: {:#?}", hit);
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
        //println!("received in distro: {:#?}", extract_str_from_bytes(&received));
        let bytes = received.len();
        total_bytes += bytes as u128;
        let seq: Vec<u8> = overlap_buf.iter().copied().chain(received.into_iter()).collect();
        let chunk = process_chunk(seq, &mut records, &mut current_record, &total_bytes, &params, &mut overlap_buf, bytes, &mut bytes_in_current_record);
        distributors[i % distributors.len()].send(chunk).expect(&format!("failed to send chunk {} to worker thread {}", i, i % distributors.len()));
        i += 1;
    }
}

fn process_chunk(seq: Vec<u8>, records: &mut impl Iterator<Item = Arc<Record>>, current_record: &mut Arc<Record>, total_bytes: &u128, params: &Params, overlap_buf: &mut Vec<u8>, received: usize, bytes_in_rec: &mut u128) -> ProcessedChunk {
    let chunk: ProcessedChunk;
    let start: usize;
    let start_in_rec: u128;
    if !overlap_buf.len() == params.k / 4 {
        start = current_record.start_bit;
        start_in_rec = 0;
    }
    else {
        start = 0;
        start_in_rec = *bytes_in_rec;
    }
    if current_record.end_byte <= *total_bytes {
        let split_idx = seq.len() - (*total_bytes - current_record.end_byte) as usize;
        *overlap_buf = seq[split_idx..].to_vec();
        //println!("processed in distro: {:#?}, buf: {:#?}, end_idx: {}, atcual: {}", extract_str_from_bytes(&seq), extract_str_from_bytes(&overlap_buf), split_idx, current_record.end_byte);
        chunk = ProcessedChunk {
            bytes: seq,
            id: current_record.id.clone(),
            end_bit: current_record.end_bit,
            end_byte: split_idx,
            start_in_rec: start_in_rec,
            start_bit: start,
            start_byte: 0
        };
        *bytes_in_rec = (received - split_idx) as u128;
        if let Some(rec) = records.next() {
            *current_record = rec;
        }
    }
    else {
        *overlap_buf = seq[0.max(seq.len().abs_diff(params.query_length))..].to_vec();
        chunk = ProcessedChunk {
            end_byte: seq.len() - 1,
            bytes: seq,
            id: current_record.id.clone(),
            end_bit: 3,
            start_in_rec: start_in_rec,
            start_bit: start as usize,
            start_byte: 0.max(params.query_length - params.k / 4)
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

fn scan(chunk: ProcessedChunk, params: Arc<Params>, query: Arc<SimpleRecord>, scheme: Arc<ScoringScheme>, hit_sender: &mpsc::Sender<Hit>) -> Result<(), String> {
    //println!("received: {:#?}, {}, {}, {}, {}", extract_str_from_bytes(&chunk.bytes[chunk.start_byte..=chunk.end_byte]), chunk.id, chunk.bytes.len(), query.seq.len(), query.words.len()/4);
    for (i, slice1) in chunk.bytes[chunk.start_byte..=chunk.end_byte].windows(3).enumerate() {
        for (j, slice2) in query.words.iter().enumerate() {
            //println!("{}, {}", i, j);
            //println!("word: {:#?}\nseq: {:#?}", extract_str_from_bytes(slice2), extract_str_from_bytes(slice1));
            if slice1[0] == slice2[0] && slice1[1] == slice2[1] && slice1[2] == slice2[2] { //TODO: rewrite this to extend hits with score T
                let mut score = scheme.hit * params.k as i16;
                //println!("{},j: {}", i, j);
                //println!("word: {:#?}", extract_str_from_bytes(slice2));
                score += extend_left(&chunk.bytes[0.max(i - (j/4).min(i))..i], &query.seq[..j/4]);
                score += extend_right(&chunk.bytes[i..chunk.bytes.len().min(i + params.query_length - j/4)], &query.seq[j/4..]);
                if score >= params.extension_threshhold {
                    if let Err(e) = hit_sender.send(Hit {
                        id: chunk.id.clone(),
                        score: score,
                        left: 0.max(i - (j/4).min(i)),
                        right: chunk.bytes.len().min(i + params.query_length - j/4)
                    }) 
                    {
                        println!("could not send hit info, {:#?}", e);
                        //return Err(format!("could not send hit: {:#?}", e));
                    }
                }
            }
        }
    }
    Ok(())
}

fn extend_left(db_seq: &[u8], query_seq: &[u8]) -> i16 {
    0
}

fn extend_right(db_seq: &[u8], query_seq: &[u8]) -> i16 {
    0
}

fn extend(db_seq: &[u8], query_seq: &[u8]) -> i8 {
    0
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

