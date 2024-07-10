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

    /// number of seeds
    #[arg(short, long, default_value = "100")]
    pub n_saved: usize,

    /// threshhold for seeding
    #[arg(short, long, default_value = "7")]
    pub threshhold: u32,

    /// length of word durign seeding
    #[arg(short, long, default_value = "7")]
    pub length: usize,
}