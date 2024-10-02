pub mod blastn;
pub mod parser;
pub mod dust;
pub mod process_db;
pub mod make_db;
pub mod blastn2;
mod benchmark;

use clap::Parser;
use bio::io::fasta::Reader;
use process_db::get_kmers;
use std::{
    collections::{BTreeMap, HashMap}, env, fs, io::{stdin, Read}, path::PathBuf, sync::{Arc, Mutex}, time::Instant
};

use blastn::{convert_to_ascii, BLASTn, Searcher};
use blastn2::{align, Params};
use parser::Args;
use benchmark::benchmark;
use dust::Dust;
use make_db::save_db::generate_db;

fn get_db(path: &str, k: usize) -> HashMap<usize, BTreeMap<u64, Vec<usize>>> {
    let db_reader = Reader::from_file(path).unwrap();
    let db = db_reader.records();
    get_kmers(k, db)
}

fn main() -> Result<(), String> {
    let args = Args::parse();
    let test = !false;
    let exe_path = env::current_exe().expect("could not get path of executable");
    let p = &(args.db_file.split('.').nth(0).unwrap().to_string() + "/");
    println!("{:?}", exe_path);
    println!("{}", p);
    if test {
        let now = Instant::now();
        let out_p: PathBuf;
        if args.out_path == "." {
            if exe_path.ends_with("target/release/blast") {
                out_p = exe_path
                    .parent()
                    .and_then(|p| p.parent())
                    .and_then(|p| p.parent())
                    .map(|p_| p_.join(p))
                    .unwrap();
            }
            else {
                out_p = exe_path
                    .parent()
                    .map(|p_| p_.join(p))
                    .unwrap();
            }
        }
        else {
            out_p = args.out_path.into();
        }

        if args.db_file.starts_with('.') {
            eprintln!("DB Path should not start with a . as that case is currently not handled correctly");
            return Ok(());
        }

        if args.new_db || fs::exists(out_p.clone()).is_err() || !fs::exists(out_p.clone()).unwrap() {
            let _ = generate_db(&args.db_file, out_p.clone());
            println!("{:?}", now.elapsed());
        }

        let params = Params {
            k: if args.length % 4 == 0 {args.length} else {12},
            extension_threshold: 24,
            scanning_threshold: args.threshold as i16,
            extension_length: 64,
            verbose: args.verbose,
            masking_threshold: args.masking_threshold,
            masking: !args.no_masking
        };

        let _ = align(out_p, &args.query_file, args.num_workers, params);

        return Ok(());
    }
    let mut t = args.threshold;
        if t > args.length {
            t = args.length;
        }

    if args.recursive > 1{
        let best_hits = benchmark(args.recursive, &args.query_file, &args.db_file, &t, &args.length, args.masking_threshold, !args.no_masking)?;
        
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

    
    let res: Vec<u8>;
    if !args.no_masking {
        let mut dust = Dust::new(64, args.masking_threshold, query.seq().to_vec());
        res = dust.mask_regions();
    }
    else {
        res = query.seq().to_vec();
    }

    //let mut searcher = Searcher::new(query.clone(), db, t, args.length, res.clone());
    let mut n_retries = 0;
    loop {
        let db_reader = Reader::from_file(&args.db_file).unwrap();
        let db = Arc::from(Mutex::from(db_reader.records()));
        let mut searcher = Searcher::new(query.clone(), db, t, args.length, res.clone());
        let now = Instant::now();
        if searcher.align().is_err() {
            println!("Masked Query: ");
            for i in res.iter() {
                print!("{}", convert_to_ascii(i));
            }
            println!();
            return Ok(());
        }
        
        if args.verbose {
            println!("\nSearch finished after {:#?}\n", now.elapsed());
            println!("Masked Query: ");
            for i in res.iter() {
                print!("{}", convert_to_ascii(i));
            }
            println!();
        }

        if args.single_result {
            let db_reader = Reader::from_file(args.db_file.clone()).unwrap();
            let summary = searcher.summary(&mut db_reader.records());
            println!("{}", summary);
        }
        else {
            let tt = searcher.sm(&args.db_file);
            if tt.hits_found.len() == 0 {
                if args.verbose {
                    println!("\nNo hits found\n");
                }
                if n_retries > 5 {
                    println!("No hits found \naborting...");
                    return Ok(());
                }
                n_retries += 1;
                continue;
            }
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
        break;
    }

    Ok(())
}

fn get_idx_from_ascii(num: &u8) -> usize {
    *num as usize - 48
}
