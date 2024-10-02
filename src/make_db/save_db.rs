use std::{fs, io::{self, Write}, path::PathBuf, str, sync, thread};
use super::{records::Record, parse_fasta::parse_and_compress_fasta};

pub fn save_compressed_db(path: PathBuf, rx: sync::mpsc::Receiver<Vec<u8>>) -> io::Result<()> {
    let file = fs::File::create(path)?;
    let mut writer = io::BufWriter::new(file);

    while let Ok(bytes) = rx.recv() {
        writer.write_all(bytes.as_slice())?;
        writer.flush()?;
    }
    Ok(())
}

pub fn save_to_csv(recs: Vec<Record>, path: PathBuf) -> io::Result<()> {
    let mut wrt = csv::Writer::from_path(path)?;
    for rec in recs {
        wrt.serialize(rec)?;
        wrt.flush()?;
    }
    Ok(())
}

pub fn generate_db(db_path: &str, out_path: PathBuf) -> io::Result<()> {
    fs::create_dir_all(out_path.clone())?;
    let db_path_ = db_path.to_string();
    let (tx, rx) = sync::mpsc::channel();
    let handle: thread::JoinHandle<io::Result<Vec<Record>>> = thread::spawn(move || {
        let records = parse_and_compress_fasta(&db_path_, 2048, tx)?;
        Ok(records)
    });
    //let out_path_ = out_path.to_string();
    let out_path_ = out_path.clone();
    let handle2: thread::JoinHandle<io::Result<()>> = thread::spawn(move || {
        save_compressed_db(out_path_.join("seq.bin"), rx)?;
        Ok(())
    });
    let out_path_ = out_path.clone();
    let records = handle.join().unwrap().unwrap();
    save_to_csv(records,  out_path_.join("records.csv"))?;
    let _ = handle2.join();
    Ok(())
}