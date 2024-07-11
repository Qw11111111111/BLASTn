pub mod blastn;
pub mod parser;

use clap::Parser;
use bio::io::fasta::Reader;
use std::{sync::{Arc, Mutex}, time::Instant};

use num::ToPrimitive;

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
    
    println!("\n Process finished \n");
    println!("Best match found in record: {} at index: {}", summary.best_idx.0, summary.best_idx.1);
    println!("Score: {}", summary.score);
    println!("\nQuery: \n");
    for char in query.clone().seq().iter() {
        print!("{}", convert_to_ascii(char));
    }
    println!("\n\nBestSeq: \n");
    for char in summary.seq.iter() {
        print!("{}", convert_to_ascii(char));
    }
    
    let similarity = summary.score.to_f64().unwrap() / query.clone().seq().len().to_f64().unwrap();
    println!("\nSimilarity: {}%", similarity * 100.0);
    
    Ok(())
}

fn convert_to_ascii(index: &u8) -> String {
    match index {
        65 => "A".into(),
        67 => "C".into(),
        71 => "G".into(),
        84 => "T".into(),
        _ => "?".into(),
    }
}
