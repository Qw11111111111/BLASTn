use std::{cmp::{min, max}, fs, io::{self, Read}, str, sync, vec};
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
        _ => panic!("Wrong nucleotide: {}", ch)
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

fn fill_byte(byte: &mut u8, seq: &[u8]) {
    //fills the byte with a clean sequence of the number of missing bits
    for (i, item) in seq.iter().enumerate() {
        *byte |= get_bit(item) << (seq.len() - i) * 2;
    }
}

fn clean_up(last_byte: Option<u8>, tx: &sync::mpsc::Sender<Vec<u8>>, records: &mut Vec<Record>, total_bytes: u128, last_missing: usize) {
    //sends the last potentially not filled byte and sets the end of the last record to the end of the entire sequence
    if let Some(byte) = last_byte {
        tx.send(vec![byte]).expect("failed tp send bytes");
        let length = records.len() - 1;
        records[length].end_byte = total_bytes - 1;
        records[length].end_bit = last_missing;
    }
    else {
        let length = records.len() - 1;
        records[length].end_byte = total_bytes - 1;
        records[length].end_bit = 3; 
    }
}

fn handle_ids(buffer: &mut Vec<u8>, bytes: &mut usize, ids: &mut Vec<Id>) -> bool {
    //appends all ids into ids and reduces the number of bytes correspondingly
    loop {
        match extract_ids(buffer) {
            (Some(id), _is_valid_buffer) => {
                *bytes -= id.id.len();
                ids.push(id);
            },
            (None, is_valid_buffer) => return is_valid_buffer,
        }
    }
}

fn ensure_filled_bytes(filtered: &mut Vec<u8>, last_byte: &mut Option<u8>, last_missing: &mut usize) -> Option<Vec<u8>> {
    // generates a byte sequence and handles remaining bytes
    if let Some(byte) = last_byte {
        if *last_missing >= filtered.len() {
            fill_byte(byte, &filtered);
            *last_missing -= filtered.len();
            if *last_missing > 0 {
                return None;
            }
            else {
                let next_bytes = Some(vec![byte.clone()]);
                *last_byte = None;
                return next_bytes;
            }
        }
        else {
            let filler = filtered.split_off(*last_missing);
            fill_byte(byte, &filler);
            *last_missing = 0;
        }
    }
    Some(generate_next_seq(filtered, last_byte, last_missing))
}

fn generate_next_seq(filtered: &mut Vec<u8>, last_byte: &mut Option<u8>, last_missing: &mut usize) -> Vec<u8> {
    let mut next_bytes: Vec<u8>;

    if let Some(byte) = last_byte {
        next_bytes = vec![byte.clone()];
    }
    else {
        next_bytes = Vec::default();
    }

    if filtered.len() % 4 == 0 {
        let mut byte_seq = parse_to_bytes(&filtered);
        next_bytes.append(&mut byte_seq);
    }
    else {
        let mut byte_seq = parse_to_bytes(&filtered);
        *last_byte = byte_seq.pop();
        *last_missing = 4 - filtered.len() % 4;
        next_bytes.append(&mut byte_seq);
    }
    
    next_bytes
}

fn append_to_records(ids: Vec<Id>, total_bytes: &mut u128, records: &mut Vec<Record>, last_missing: usize) {
    //appends and converts all ids into records, filling out their ranges and reducing the next_bytes accordingly
    for id in ids {
        let start_byte = *total_bytes + id.start as u128 / 4;
        if records.len() > 0 {
            let length = records.len() - 1;
            let mut end_bit = ((id.start.abs_diff(last_missing) as i32) % 4) as usize;
            if end_bit == 0 {
                end_bit = 3;
            }
            else {
                end_bit -= 1;
            }
            let end_byte: u128;
            if end_bit == 3 {
                end_byte = max(start_byte, 1) - 1;
            }
            else {
                end_byte = start_byte;
            }
            records[length].end_bit = end_bit;
            records[length].end_byte = end_byte;
        }
        records.push(
            Record {
                id: String::from_utf8(id.id).expect("failed to cast id to str"),
                start_byte: start_byte,
                start_bit: ((id.start.abs_diff(last_missing) as i32) % 4) as usize,
                end_bit: 0,
                end_byte: 0
            }
        )
    }
}

fn extract_ids(buf: &mut Vec<u8>) -> (Option<Id>, bool) {
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
                return (id, true); 
            }
        }
    }
    (None, start_of_id == None)
}

pub fn parse_and_compress_fasta(path: &str, chunk_size: usize, tx: sync::mpsc::Sender<Vec<u8>>) -> io::Result<()> {
    let file = fs::File::open(path)?;
    let mut reader = io::BufReader::new(file);
    let mut records: Vec<Record> = Vec::default();
    let mut total_bytes: u128 = 0;
    let mut last_byte: Option<u8> = None;
    let mut last_missing: usize = 0;
    let mut bytes: usize = 0;
    let mut overflow: Vec<u8> = Vec::default();

    loop {
        let mut buffer = vec![0; chunk_size];
        bytes += reader.read(&mut buffer)?;
        overflow.append(&mut buffer);
        if bytes == 0 {
            clean_up(last_byte, &tx, &mut records, total_bytes, last_missing);
            break;
        }
        
        let mut ids = vec![];

        if !handle_ids(&mut overflow, &mut bytes, &mut ids) {
            continue;
        }
        //filters all line feeds and fills all bytes except for the last one, which remains for the next iteration

        let mut filtered: Vec<u8> = overflow[..bytes]
            .iter()
            .copied()
            .filter(|&item| item != 10)
            .collect();
        

        append_to_records(ids, &mut total_bytes, &mut records, last_missing);

        if let Some(next_bytes) = ensure_filled_bytes(&mut filtered, &mut last_byte, &mut last_missing) {
            total_bytes += next_bytes.len() as u128;
            tx.send(next_bytes).expect("failed to send bytes");
            overflow = Vec::default();
            bytes = 0;
        }
    }
    save_to_csv(records,  &(path.split('.').nth(0).unwrap().to_string() + ".csv"))?;
    Ok(())
}