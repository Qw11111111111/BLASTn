#![allow(
    clippy::too_many_arguments,
    clippy::upper_case_acronyms,
    unused_imports
)]
pub mod blastn2;
pub mod dust;
pub mod global;
pub mod make_db;
pub mod parser;

use clap::Parser;
use std::{env, fs, path::PathBuf, str::FromStr, time::Instant};

use blastn2::{align, Params};
use make_db::save_db::generate_db;
use parser::Args;

fn main() -> Result<(), String> {
    let args = Args::parse();
    let db_path: PathBuf = "genomes/".into();
    let exe_path = env::current_exe().expect("could not get path of executable");

    let exe_path = exe_path
        .parent()
        .and_then(|p| p.parent())
        .and_then(|p| p.parent());
    let now = Instant::now();

    // assuming the executable sits in dir/target/release/ or similar
    let out_p: PathBuf = if args.out_path == db_path {
        exe_path
            .unwrap()
            .join(args.out_path.join(args.db_file.parent().unwrap()))
    } else {
        args.out_path
    };

    if args.new_db || fs::exists(out_p.clone()).is_err() || !fs::exists(out_p.clone()).unwrap() {
        let _ = generate_db(args.db_file.clone(), out_p.clone());
        println!("db generated in {:?}", now.elapsed());
    }

    let params = Params {
        k: if args.length % 4 == 0 {
            args.length
        } else {
            12
        },
        extension_threshold: 24,
        scanning_threshold: args.threshold as i16,
        extension_length: 64,
        verbose: args.verbose,
        masking_threshold: args.masking_threshold,
        masking: !args.no_masking,
    };

    let _ = align(out_p, args.query_file, args.num_workers, params);
    Ok(())
}

pub fn convert_to_ascii(index: &u8) -> String {
    char::from_u32(*index as u32).unwrap().to_string()
}
