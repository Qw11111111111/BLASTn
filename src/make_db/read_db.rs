use std::{fs, io::{self, Read}, path::PathBuf, str, sync::mpsc};
use super::records::VecRecord;

pub fn read_compressed_db(path: &str) -> io::Result<String> {
    let file = fs::File::open(path)?;
    let mut reader = io::BufReader::new(file);
    let mut buf = vec![0;1024];
    let mut all_nts = "".to_string();
    loop {
        let bytes = reader.read(&mut buf)?;
        if bytes == 0 {
            break;
        }
        //println!("bytes read: {}", bytes);
        let nts = extract_str_from_bytes(&buf[..bytes]);
        println!("{nts}");
        all_nts += &nts;
    }
    Ok(all_nts)
}

pub fn extract_str_from_bytes(bytes: &[u8]) -> String {
    let mut nts = "".to_string();
    for byte in bytes.iter() {
        let chars = (0..=6)
            .step_by(2)
            .map(|i| {
                let bit = (byte >> (6 - i)) & 0b11;
                match bit {
                    0b00 => 'A',
                    0b01 => 'G',
                    0b10 => 'T',
                    0b11 => 'C',
                    _ => 'N'
                }
            });
        nts.extend(chars);
    }
    nts
}

pub fn bytes_to_chars(bytes: &[u8], end_bit: usize, start_bit: usize) -> Vec<u8> {
    let mut nts = Vec::default();
    for (j, byte) in bytes.iter().enumerate() {
        let end: usize;
        let start: usize;
        if j == bytes.len() - 1 {
            end = end_bit * 2;
        }
        else {
            end = 6;
        }
        if j == 0 {
            start = start_bit * 2;
        }
        else {
            start = 0;
        }
        let chars = (start..=end)
            .step_by(2)
            .map(|i| {
                let bit = (byte >> (6 - i)) & 0b11;
                match bit {
                    0b00 => 'A' as u8,
                    0b01 => 'G' as u8,
                    0b10 => 'T' as u8,
                    0b11 => 'C' as u8,
                    _ => '-' as u8
                }
            });
        nts.extend(chars);
    }
    nts
}


pub fn parse_compressed_db_lazy(path: PathBuf, chunk_size: usize, sender: mpsc::Sender<Vec<u8>>) -> io::Result<()> {
    let file = fs::File::open(path.join("seq.bin"))?;
    let mut reader = io::BufReader::new(file);
    let mut buf = vec![0; chunk_size];

    loop {
        let bytes = reader.read(&mut buf)?;
        if bytes == 0 {
            break;
        }
        //println!("read: {:#?}", extract_str_from_bytes(&buf[..bytes]));
        sender.send(buf[..bytes].to_vec()).expect("couldnt send buf");
    }

    Ok(())
}

pub fn read_csv(path: PathBuf) -> io::Result<VecRecord> {
    let mut records = VecRecord::default();
    let mut reader = csv::Reader::from_path(path.join("records.csv"))?;
    for rec in reader.deserialize() {
        records.push(rec?);
    }

    Ok(records)
}