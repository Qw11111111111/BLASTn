use std::{fs, io::{self, Read}, str, sync::mpsc};
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
        println!("bytes read: {}", bytes);
        let nts = extract_str_from_bytes(&buf[..bytes]);
        println!("{nts}");
        all_nts += &nts;
    }
    Ok(all_nts)
}

fn extract_str_from_bytes(bytes: &[u8]) -> String {
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

pub fn parse_compressed_db_lazy(path: &str, chunk_size: usize, sender: mpsc::Sender<Vec<u8>>) -> io::Result<()> {
    let file = fs::File::open(path)?;
    let mut reader = io::BufReader::new(file);
    let mut buf = vec![0; chunk_size];

    loop {
        let bytes = reader.read(&mut buf)?;
        if bytes == 0 {
            break;
        }
        sender.send(buf[..bytes].to_vec()).expect("couldnt send buf");
    }

    Ok(())
}

pub fn read_csv(path: &str) -> io::Result<VecRecord> {
    let mut records = VecRecord::default();
    let mut reader = csv::Reader::from_path(path)?;
    for rec in reader.deserialize() {
        records.push(rec?);
    }

    Ok(records)
}