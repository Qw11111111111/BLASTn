use std::{fs, io::{self, Read}, sync, cmp::min, str};
use crate::make_db::{records::Record, save_db::*};

struct Id {
    id: Vec<u8>,
    start: usize
}

fn get_bit(ch: &u8) -> u8 {
    match ch {
        65 => 0b00, //A
        71 => 0b01, //G
        84 => 0b10, //T
        67 => 0b11, //C
        _ => panic!("Wrong nucleotide")
    }
}

pub fn seq_to_byte(seq: &[u8]) -> u8 {
    let mut byte = 0b00;
    for (i, ch) in seq.iter().enumerate() {
        if *ch == 10 {
            continue;
        }
        byte |= get_bit(ch) << (3 - i) * 2;
    }
    byte
}

fn parse_to_bytes(buf: &[u8]) -> Vec<u8> {
    (0..buf.len())
        .step_by(4)
        .map(|i| seq_to_byte(&buf[i..min(i+4, buf.len())]))
        .collect()
}

fn fill_byte(mut byte: u8, seq: &[u8]) -> Option<u8> {
    for (i, item) in seq.iter().enumerate() {
        if *item == 10 {
            return None;
        }
        byte |= get_bit(item) << (seq.len() - i) * 2;
    }
    Some(byte)
}

fn clean_up(last_byte: Option<u8>, tx: &sync::mpsc::Sender<Vec<u8>>, records: &mut Vec<Record>, total_bytes: u128, last_missing: usize) {
    //sends the last potentially not filled byte and sets the end of the last record to the end of the entire sequence
    //TODO: correct the eend bits and bytes values
    if let Some(byte) = last_byte {
        tx.send(vec![byte]).expect("failed tp send bytes");
        let length = records.len() - 1;
        records[length].end_byte = total_bytes;
        records[length].end_bit = last_missing;
    }
    else {
        let length = records.len() - 1;
        records[length].end_byte = total_bytes;
        records[length].end_bit = 3; 
    }
}

fn handle_ids(buffer: &mut Vec<u8>, bytes: &mut usize, ids: &mut Vec<Id>) {
    //appends all ids into ids and reduces the number of bytes correspondingly
    while let Some(id) = extract_ids(buffer) {
        *bytes -= id.id.len();
        ids.push(id);
    }
}

fn ensure_filled_bytes(filtered: Vec<u8>, last_byte: &mut Option<u8>, last_missing: &mut usize) -> (Vec<u8>, u128) {
    //makes sure that all bytes except for the last one in the database are filled with four nucleotides
    //TODO: this has a lot of code duplication and can be further abstracted
    let mut byte_seq: Vec<u8>;
    if filtered.len() % 4 != 0 {
        if let Some(byte) = *last_byte {
            let (filler, remainder) = filtered.split_at(*last_missing);
            let filled_byte_: u8;
            loop {
                if let Some(filled_byte) = fill_byte(byte, filler) {
                    filled_byte_ = filled_byte;
                    break;
                }
            }
            if remainder.len() % 4 != 0 {
                byte_seq = parse_to_bytes(&remainder);
                *last_byte = Some(byte_seq.remove(byte_seq.len() - 1));
                *last_missing = 4 - remainder.len() % 4;
                let mut next = vec![filled_byte_];
                next.append(&mut byte_seq);
                byte_seq = next;
            }
            else {
                byte_seq = parse_to_bytes(&remainder);
                let mut next = vec![filled_byte_];
                next.append(&mut byte_seq);
                byte_seq = next;
            }
        }
        else {
            byte_seq = parse_to_bytes(&filtered);
            *last_byte = Some(byte_seq.remove(byte_seq.len() - 1));
            *last_missing = 4 - filtered.len() % 4;
        }
    }
    else {
        if let Some(byte) = *last_byte {
            let (filler, remainder) = filtered.split_at(*last_missing);
            let filled_byte_: u8;
            loop {
                if let Some(filled_byte) = fill_byte(byte, filler) {
                    filled_byte_ = filled_byte;
                    break;
                }
            }
            if remainder.len() % 4 != 0 {
                byte_seq = parse_to_bytes(&remainder);
                *last_byte = Some(byte_seq.remove(byte_seq.len() - 1));
                *last_missing = 4 - remainder.len() % 4;
                let mut next = vec![filled_byte_];
                next.append(&mut byte_seq);
                byte_seq = next;
            }
            else {
                byte_seq = parse_to_bytes(&remainder);
                let mut next = vec![filled_byte_];
                next.append(&mut byte_seq);
                byte_seq = next;
            }
        }
        else {
            byte_seq = parse_to_bytes(&filtered);
        }
    }
    let next_bytes = byte_seq.len() as u128;
    (byte_seq, next_bytes)
}

fn append_to_records(ids: Vec<Id>, total_bytes: &mut u128, next_bytes: &mut u128, records: &mut Vec<Record>) {
    //appends and converts all ids into records, filling out their ranges and reducing the next_bytes accordingly
    //TODO: correct the calculations for start and end bytes and bits
    for id in ids {
        let start_byte = *total_bytes + id.start as u128 / 4;
        *next_bytes -= start_byte;
        *total_bytes = 0;
        if records.len() > 0 {
            let length = records.len() - 1;
            records[length].end_byte = start_byte - 1;
            records[length].end_bit = ((id.start as i32 - 1) % 4 - 3).abs() as usize;
        }
        records.push(
            Record {
                id: String::from_utf8(id.id).expect("failed to cast id to str"),
                start_byte: start_byte,
                start_bit: ((id.start as i32) % 4 - 3).abs() as usize,
                end_bit: 0,
                end_byte: 0
            }
        )
    }
}

fn extract_ids(buf: &mut Vec<u8>) -> Option<Id> {
    let mut start_of_id: Option<usize> = None;
    let mut end_of_id: Option<usize> = None;
    let mut n_new_lines_pre_id = 0;
    for (i, &item) in buf.iter().enumerate() {
        if item == 62 {
            start_of_id = Some(i);
        }
        if item == 10 {
            if let Some(_) = start_of_id {
                if end_of_id == None {
                    end_of_id = Some(i);
                }
            }
            else {
                n_new_lines_pre_id += 1;
            }
        }
        if let Some(start) = start_of_id {
            if let Some(end) = end_of_id {
                let id = Some(Id {id:buf[start..end].to_vec(), start: start - n_new_lines_pre_id});
                buf.splice(start..end, Vec::default());
                return id; 
            }
        }
    }
    None
}

pub fn parse_and_compress_fasta(path: &str, chunk_size: usize, tx: sync::mpsc::Sender<Vec<u8>>) -> io::Result<()> {
    let file = fs::File::open(path)?;
    let mut reader = io::BufReader::new(file);
    let mut buffer = vec![0; chunk_size];
    let mut records: Vec<Record> = Vec::default();
    let mut total_bytes: u128 = 0;
    let mut last_byte: Option<u8> = None;
    let mut last_missing: usize = 0;

    loop {
        let mut bytes = reader.read(&mut buffer)?;
        if bytes == 0 {
            clean_up(last_byte, &tx, &mut records, total_bytes, last_missing);
            break;
        }
        
        let mut ids = vec![];

        handle_ids(&mut buffer, &mut bytes, &mut ids);
        //filters all line feeds and fills all bytes except for the last one, which remains for the next iteration

        let filtered: Vec<u8> = buffer[..bytes]
            .iter()
            .copied()
            .filter(|&item| item != 10)
            .collect();
        
        let (byte_seq, mut next_bytes) = ensure_filled_bytes(filtered, &mut last_byte, &mut last_missing);
        append_to_records(ids, &mut total_bytes, &mut next_bytes, &mut records);

        total_bytes += next_bytes;
        tx.send(byte_seq).expect("failed tp send bytes");//.inspect_err(|e| eprintln!("failed to send bytes, {}", e));
    }
    save_to_csv(records, "genomes/records.csv")?;
    Ok(())
}