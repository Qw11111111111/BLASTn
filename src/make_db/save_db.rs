use std::{fs, io::{self, Write}, str, sync, thread, time::Instant};
use crate::make_db::{records::Record, parse_fasta::parse_and_compress_fasta};

pub fn save_compressed_db(path: &str, rx: sync::mpsc::Receiver<Vec<u8>>) -> io::Result<()> {
    let file = fs::File::create(path)?;
    let mut writer = io::BufWriter::new(file);

    while let Ok(bytes) = rx.recv() {
        writer.write_all(bytes.as_slice())?;
        writer.flush()?;
    }
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
    let now = Instant::now();
    let (tx, rx) = sync::mpsc::channel();
    let handle: thread::JoinHandle<io::Result<()>> = thread::spawn(move || {
        parse_and_compress_fasta("genomes/ecoli.fna", 2048, tx)?;
        Ok(())
    });
    println!("{:#?}", now.elapsed());
    let handle2: thread::JoinHandle<io::Result<()>> = thread::spawn(move || {
        save_compressed_db("genomes/ecoli.bin", rx)?;
        Ok(())
    });
    let _ = handle.join();
    let _ = handle2.join();
    Ok(())
}