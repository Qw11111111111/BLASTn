pub mod parser;
pub mod dust;
pub mod make_db;
pub mod blastn2;

use clap::Parser;
use std::{
    env, fs, path::PathBuf, time::Instant
};

use blastn2::{align, Params};
use parser::Args;
use make_db::save_db::generate_db;

fn main() -> Result<(), String> {
    let args = Args::parse();
    let exe_path = env::current_exe().expect("could not get path of executable");
    let p = &(args.db_file.split('.').nth(0).unwrap().to_string() + "/");
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
    Ok(())
}

pub fn convert_to_ascii(index: &u8) -> String {
    char::from_u32(*index as u32).unwrap().to_string()
}
