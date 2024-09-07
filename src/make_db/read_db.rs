use std::{fs, io::{self, Read}, str};
use crate::make_db::records::{Record, Records};

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

pub fn parse_compressed_db(path: &str) -> io::Result<Records> {
    
    let records = read_csv(&(path.to_string() + "records.csv"));

    

    Ok(Records::default())
}

fn read_csv(path: &str) -> Vec<Record> {
    let records = Vec::default();

    records
}