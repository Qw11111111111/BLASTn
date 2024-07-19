pub mod blastn;
pub mod parser;
pub mod dust;
mod benchmark;

use clap::Parser;
use bio::io::fasta::Reader;
use num::ToPrimitive;
use std::{
    io::{stdin, Read}, 
    sync::{Arc, Mutex}, 
    time::Instant
};

use blastn::{Searcher, convert_to_ascii};
use parser::Args;
use benchmark::benchmark;
use dust::Dust;

fn main() -> Result<(), String> {
    //TODO: proper Error handling
    //TODO: Rewrite the search to match the procedure outlined here: https://en.wikipedia.org/wiki/BLAST_(biotechnology), as the current implementation is rather naive.
    let args = Args::parse();

    let mut t = args.threshhold;
        if t > args.length {
            t = args.length;
        }

    if args.recursive > 1{
        let best_hits = benchmark(args.recursive, &args.query_file, &args.db_file, &t, &args.length);
        
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

    //TODO: CLI acces to window_size and threshold.
    let mut dust = Dust::new(64, 0.8, query.seq().to_vec());
    let res = dust.mask_regions();
    let mut searcher = Searcher::new(query.clone(), db, t, args.length, res.clone());
    let now = Instant::now();
    searcher.align();
    
    if args.verbose {
        println!("Search finished after {:#?}\n", now.elapsed());
        println!("Masked Query: ");
        for i in res.iter() {
            print!("{}", convert_to_ascii(i));
        }
        println!();
    }

    if !args.single_result {
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
