pub mod blastn;
pub mod parser;
mod benchmark;

use clap::Parser;
use bio::io::fasta::Reader;
use num::ToPrimitive;
use std::{sync::{Arc, Mutex}, time::Instant};

use blastn::Searcher;
use parser::Args;
use crate::benchmark::benchmark;

fn main() -> Result<(), String> {
    //TODO: proper Error handling
    let args = Args::parse();

    if args.benchmark {
        let (t, best_hits) = benchmark(args.retries, &args.query_file, &args.db_file, &args.threshhold, &args.length);
        println!("\ntotal time: {:?}s", t);
        println!("time per run: {:?}s\n", t / args.retries.to_f64().unwrap());
        println!("best hits: ");
        for item in best_hits.iter() {
            println!("Record: [{0}] | Index: [{1}] | Score: {2}", item.0.0, item.0.1,item.1);

        }

        return Ok(());
    }

    let query_reader = Reader::from_file(&args.query_file).unwrap();
    let db_reader = Reader::from_file(&args.db_file).unwrap();
    
    let db = Arc::from(Mutex::from(db_reader.records()));
    let query = Arc::from(query_reader.records().next().unwrap().unwrap());
    
    let mut searcher = Searcher::new(query.clone(), db, args.threshhold, args.length);
    
    let now = Instant::now();
    
    searcher.align();
    
    if args.verbose {
        println!("Search finished after {:#?}\n", now.elapsed());
    }

    let db_reader = Reader::from_file(args.db_file).unwrap();
    let summary = searcher.summary(&mut db_reader.records());
    println!("{}", summary);

    Ok(())
}
