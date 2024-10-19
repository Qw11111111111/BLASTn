use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// Path to query sequence
    #[arg(short, long)]
    pub query_file: PathBuf,

    /// Path to databank
    #[arg(short, long)]
    pub db_file: PathBuf,

    /// threshold for seeding
    #[arg(short, long, default_value = "60")]
    pub threshold: usize,

    /// length of word durign seeding
    #[arg(short, long, default_value = "12")]
    pub length: usize,

    /// Additional printout during run time
    #[arg(short, long, default_value = "false")]
    pub verbose: bool,

    /// do not mask low complexity regions with DUST
    #[arg(short, long, default_value = "false")]
    pub no_masking: bool,

    /// Threshold for a sequence to be masked
    #[arg(long, default_value = "1.0")]
    pub masking_threshold: f64,

    /// Force new database generation
    #[arg(long, default_value = "false")]
    pub new_db: bool,

    /// num workers for scanning
    #[arg(long, default_value = "10")]
    pub num_workers: usize,

    /// out path for db files
    #[arg(long, short, default_value = "genomes/")]
    pub out_path: PathBuf,
}
