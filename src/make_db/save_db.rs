use std::{fs, io::{self, Write}, str, sync, thread};
use super::{records::Record, parse_fasta::parse_and_compress_fasta};

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

pub fn generate_db(db_path: &str, out_path: &str) -> io::Result<()> {
    fs::create_dir_all(out_path)?;
    let db_path_ = db_path.to_string();
    let (tx, rx) = sync::mpsc::channel();
    let handle: thread::JoinHandle<io::Result<Vec<Record>>> = thread::spawn(move || {
        let records = parse_and_compress_fasta(&db_path_, 2048, tx)?;
        Ok(records)
    });
    let out_path_ = out_path.to_string();
    let handle2: thread::JoinHandle<io::Result<()>> = thread::spawn(move || {
        save_compressed_db(&(out_path_ + "seq.bin"), rx)?;
        Ok(())
    });
    let records = handle.join().unwrap().unwrap();
    save_to_csv(records,  &(out_path.to_string() + "records.csv"))?;
    let _ = handle2.join();
    Ok(())
}