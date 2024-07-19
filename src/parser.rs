use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// Path to query sequence
    #[arg(short, long)]
    pub query_file: String,

    /// Path to databank
    #[arg(short, long)]
    pub db_file: String,

    /// threshhold for seeding
    #[arg(short, long, default_value = "11")]
    pub threshhold: usize,

    /// length of word durign seeding
    #[arg(short, long, default_value = "11")]
    pub length: usize,

    /// Additional printout during run time
    #[arg(short, long, default_value = "false")]
    pub verbose: bool,

    /// perform the search multiple times
    #[arg(short, long, default_value = "1")]
    pub recursive: usize,

    /// retries for benchmark
    #[arg(short, long, default_value = "true")]
    pub single_result: bool
}