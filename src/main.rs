pub mod blastn;
pub mod parser;
mod benchmark;

use clap::Parser;
use bio::io::fasta::Reader;
use num::ToPrimitive;
use std::{io::{stdin, Read}, sync::{Arc, Mutex}, time::Instant};

use blastn::Searcher;
use parser::Args;
use benchmark::benchmark;

fn main() -> Result<(), String> {
    //TODO: proper Error handling
    let args = Args::parse();

    if args.benchmark {
        let best_hits = benchmark(args.retries, &args.query_file, &args.db_file, &args.threshhold, &args.length);
        
        let mut buf = [0, 0];
        loop {
            println!("{}", best_hits);
            println!("print a sequence [index]/[N]");
            let _  = stdin().read(&mut buf);

            if buf[0] < 58 && buf[0] > 47 {
                best_hits.print(get_idx_from_ascii(&buf[0]));
            }
            else {
                break;
            }
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

fn get_idx_from_ascii(num: &u8) -> usize {
    num.to_usize().unwrap() - 48
}
