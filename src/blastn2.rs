use std::{env::current_exe, fs, io, sync::{mpsc, Arc, RwLock}, thread::{self, JoinHandle}};

use crate::make_db::{records::{SimpleRecord, Record, VecRecord}, parse_fasta::{parse_small_fasta, parse_to_bytes}, read_db::{read_csv, parse_compressed_db_lazy}};


struct Params {
    k: usize,
    extension_threshhold: i64,
    scanning_threshhold: i64,
    extension_length: usize
}

struct ProcessedChunk {
    bytes: Vec<u8>,
    end_bit: usize,
    start_byte_not_full: bool,
    id: String
}

pub fn align<'a>(path_to_db: &'a str, path_to_query: &str, num_workers: usize, params: Arc<Params>) -> io::Result<()> {
    let mut query = parse_small_fasta(path_to_query)?;
    get_query_words(params.k, &mut query);

    let records = read_csv(&(path_to_db.to_string() + "records.csv"))?.into_iter();

    let (raw_chunk_sender, raw_chunk_receiver) = mpsc::channel::<Vec<u8>>();
    let p = path_to_db.to_string();
    let reader: JoinHandle<io::Result<()>> = thread::spawn(move || {
        parse_compressed_db_lazy(&p, 2048, raw_chunk_sender)?;
        Ok(())
    });

    let mut worker_senders = Vec::default();
    let mut workers = Vec::default();
    let query = Arc::new(query);
    for _ in 0..num_workers {
        let (tx, rx) = mpsc::channel::<ProcessedChunk>();
        let params = Arc::clone(&params);
        let query = Arc::clone(&query);
        workers.push(thread::spawn(move || {
            while let Ok(chunk) = rx.recv() {
                scan(chunk, Arc::clone(&params), Arc::clone(&query));
            }
        }));
        worker_senders.push(tx);
    }

    let params = Arc::clone(&params);
    let distributor = thread::spawn(move || {
        distribute_chunks(raw_chunk_receiver, worker_senders, 2048, params, records);
    });
    
    let _ = reader.join();
    let _ = distributor.join();
    for handle in workers {
        let _ = handle.join();
    }

    Ok(())
}

fn distribute_chunks(chunk_receiver: mpsc::Receiver<Vec<u8>>, distributors: Vec<mpsc::Sender<ProcessedChunk>>, chunk_size: usize, params: Arc<Params>, mut records: impl Iterator<Item = Arc<Record>>) {
    let mut i = 0;
    let mut total_bytes = 0;
    let mut overlap_buf: Vec<u8> = vec![];
    let mut current_record = records.next().unwrap();

    while let Ok(received) = chunk_receiver.recv() {
        let mut package = overlap_buf.clone();
        package.extend(received.iter());
        total_bytes += received.len() as u128;
        let (start_bit, end_bit) = process_chunk(&mut overlap_buf, &mut package, &mut total_bytes, &mut records, &mut current_record, &params);
        distributors[i % distributors.len()].send(ProcessedChunk {bytes: package, end_bit: end_bit, id: current_record.id.clone(), start_byte_not_full: start_bit}).expect(&format!("failed to send chunk {} to worker thread {}", i, i % distributors.len()));
        i += 1;
    }
}

fn process_chunk(buf: &mut Vec<u8>, package: &mut Vec<u8>, bytes: &mut u128, records: &mut impl Iterator<Item = Arc<Record>>, current_record: &mut Arc<Record>, params: &Arc<Params>) -> (bool, usize) {
    let (start, end): (bool, usize);
    if *bytes > current_record.end_byte {
        let byte_split = (*bytes - current_record.end_byte) as usize;
        let (byte1, byte2): (u8, u8);
        if current_record.end_bit == 3 {
            (byte1, byte2) = (package[byte_split], package[byte_split+1]);
            (start, end) = (false, 3);
        }
        else {
            (byte1, byte2) = split_byte(package[byte_split], current_record.end_bit);
            (start, end) = (true, current_record.end_bit);
        }
        *buf = vec![byte2];
        buf.append(&mut package.split_off(byte_split + 1));
        package[byte_split] = byte1;
        records.next();
    }
    else {
        *buf = vec![];
        buf.extend(package[package.len() - params.extension_length..].iter());
        (start, end) = (false, 3);
    }
    (start, end)
}

fn split_byte(byte: u8, split_idx: usize) -> (u8, u8) {
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

fn scan(chunk: ProcessedChunk, params: Arc<Params>, query: Arc<SimpleRecord>) {
    for slice1 in chunk.bytes.windows(3) {
        for word in &query.words {
            for slice2 in word.windows(3) {
                if slice1[0] == slice2[0] && slice1[1] == slice2[1] && slice1[2] == slice2[2] {
                    extend();
                }
            }
        }
    }
}

fn extend() {

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

