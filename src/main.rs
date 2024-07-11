pub mod blastn;
pub mod parser;

use clap::Parser;
use bio::io::fasta::Reader;
use std::{sync::{Arc, Mutex}, time::Instant};

use blastn::Searcher;
use parser::Args;

fn main() -> Result<(), String> {
    let args = Args::parse();
    
    let query_reader = Reader::from_file(args.query_file).unwrap();
    let db_reader = Reader::from_file(&args.db_file).unwrap();
    let db = Arc::from(Mutex::from(db_reader.records()));
    let query = Arc::from(query_reader.records().next().unwrap().unwrap());
    
    let mut searcher = Searcher::new(query.clone(), db, args.threshhold, args.length);
    let now = Instant::now();
    searcher.align();
    
    println!("Search finished after {:?}s", now.elapsed());

    let db_reader = Reader::from_file(args.db_file).unwrap();
    let summary = searcher.summary(&mut db_reader.records());
    println!("{}", summary);
    Ok(())
}
