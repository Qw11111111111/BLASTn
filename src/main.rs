pub mod blastn;
pub mod parser;
pub mod dust;
pub mod process_db;
mod benchmark;

use clap::Parser;
use bio::io::fasta::Reader;
use num::ToPrimitive;
use process_db::get_kmers;
use std::{
    io::{stdin, Read}, 
    sync::{Arc, Mutex}, 
    time::Instant,
    collections::{HashMap, BTreeMap}
};

use blastn::{convert_to_ascii, BLASTn, Searcher};
use parser::Args;
use benchmark::benchmark;
use dust::Dust;

fn get_db(path: &str, k: usize) -> HashMap<usize, BTreeMap<u64, Vec<usize>>> {
    let db_reader = Reader::from_file(path).unwrap();
    let db = db_reader.records();
    get_kmers(k, db)
}

fn main() -> Result<(), String> {
    //TODO: proper Error handling
    //TODO: Rewrite the search to match the procedure outlined here: https://en.wikipedia.org/wiki/BLAST_(biotechnology), as the current implementation is rather naive.
    let args = Args::parse();

    let mut t = args.threshold;
        if t > args.length {
            t = args.length;
        }

    if args.recursive > 1{
        let best_hits = benchmark(args.recursive, &args.query_file, &args.db_file, &t, &args.length, args.masking_threshold, args.mask_low_complexity);
        
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
    let query = Arc::from(query_reader.records().next().unwrap().unwrap());

    if args.experimental {
        let mut dust = Dust::new(64, 0.8, query.seq().to_vec());
        let res = dust.mask_regions();
        
        let now = Instant::now();
        let db_kmers = get_db(&args.db_file, args.length);
        println!("\n{:#?} db built", now.elapsed());
        let db_reader = Reader::from_file(&args.db_file).unwrap();
        let mut blastn = BLASTn::new(db_reader.records(), Arc::from(res), args.threshold, args.threshold, args.length);

        let now = Instant::now();

        let result = blastn.align(Arc::from(db_kmers));

        println!("\n{:#?}, search finished", now.elapsed());
        println!("len of q: {}", query.seq().len());

        for res in result.iter() {
            println!("REC: {}",res.0);
            println!("len: {}", res.1.len());
            println!("best idx: {}, len: {}", res.1[res.1.len() - 1].start_stop.0, res.1[res.1.len() - 1].start_stop.1 - res.1[res.1.len() - 1].start_stop.0);
            println!();
        }

        return Ok(());
    }

    let db_reader = Reader::from_file(&args.db_file).unwrap();
    let db = Arc::from(Mutex::from(db_reader.records()));
    
    let res: Vec<u8>;
    if args.mask_low_complexity {
        let mut dust = Dust::new(64, args.masking_threshold, query.seq().to_vec());
        res = dust.mask_regions();
    }
    else {
        res = query.seq().to_vec();
    }

    let mut searcher = Searcher::new(query.clone(), db, t, args.length, res.clone());
    let now = Instant::now();
    searcher.align();
    
    if args.verbose {
        println!("\nSearch finished after {:#?}\n", now.elapsed());
        println!("Masked Query: ");
        for i in res.iter() {
            print!("{}", convert_to_ascii(i));
        }
        println!();
    }

    if args.single_result {
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
