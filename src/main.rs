pub mod blastn;
pub mod parser;
mod benchmark;

use clap::Parser;
use bio::io::fasta::Reader;
use num::ToPrimitive;
use std::{
    io::{stdin, Read}, 
    sync::{Arc, Mutex}, 
    time::Instant
};

use blastn::Searcher;
use parser::Args;
use benchmark::benchmark;

fn main() -> Result<(), String> {
    //TODO: proper Error handling
    let args = Args::parse();

    let mut t = args.threshhold;
        if t > args.length {
            t = args.length;
        }

    if args.recursive {
        let best_hits = benchmark(args.n_retries, &args.query_file, &args.db_file, &t, &args.length);
        
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
    
    let mut searcher = Searcher::new(query.clone(), db, t, args.length);
    
    let now = Instant::now();
    
    searcher.align();
    
    if args.verbose {
        println!("Search finished after {:#?}\n", now.elapsed());
    }
    if !args.extensive_result {
        let db_reader = Reader::from_file(args.db_file).unwrap();
        let summary = searcher.summary(&mut db_reader.records());
        println!("{}", summary);
    }
    else {
        let tt = searcher.sm(&args.db_file);
        let mut buf = [0, 0];
        loop {
            println!("{}", tt);
            println!("print a sequence [index]/[N]");
            let _  = stdin().read(&mut buf);

            if buf[0] < 58 && buf[0] > 47 {
                tt.print(get_idx_from_ascii(&buf[0]));
            }
            else {
                break;
            }
        }
    }

    Ok(())
}

fn get_idx_from_ascii(num: &u8) -> usize {
    num.to_usize().unwrap() - 48
}
