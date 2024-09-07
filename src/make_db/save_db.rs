use std::{fs, io::{self, Write}, sync, thread, str};
use crate::make_db::{records::Record, parse_fasta::parse_and_compress_fasta};

pub fn save_compressed_db(path: &str, rx: sync::mpsc::Receiver<Vec<u8>>) -> io::Result<()> {
    let file = fs::File::create(path)?;
    let mut writer = io::BufWriter::new(file);

    let worker: thread::JoinHandle<io::Result<()>> = thread::spawn(move || {
        while let Ok(bytes) = rx.recv() {
            writer.write_all(bytes.as_slice())?;
            writer.flush()?;
        }
        Ok(())
    });
    let _ = worker.join().expect("worker join failed");
    Ok(())
}

pub fn save_to_csv(recs: Vec<Record>, path: &str) -> io::Result<()> {
    let mut wrt = csv::Writer::from_path(path)?;
    for rec in recs {
        wrt.serialize(rec)?;
        wrt.flush()?;
    }
    Ok(())
}

pub fn example() -> io::Result<()> {
    let (tx, rx) = sync::mpsc::channel();
    parse_and_compress_fasta("genomes/seq3.fna", 12, tx)?;
    save_compressed_db("genomes/seq3.bin", rx)?;
    Ok(())
}