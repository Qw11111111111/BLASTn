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

    /// threshold for seeding
    #[arg(short, long, default_value = "11")]
    pub threshold: usize,

    /// length of word durign seeding
    #[arg(short, long, default_value = "11")]
    pub length: usize,

    /// Additional printout during run time
    #[arg(short, long, default_value = "false")]
    pub verbose: bool,

    /// perform the search multiple times
    #[arg(short, long, default_value = "1")]
    pub recursive: usize,

    /// only show best match
    #[arg(short, long, default_value = "false")]
    pub single_result: bool,

    /// do not mask low complexity regions with DUST
    #[arg(short, long, default_value = "false")]
    pub mask_low_complexity: bool,

    /// Threshold for a sequence to be masked
    #[arg(long, default_value = "1.0")]
    pub masking_threshold: f64,

    /// Threshold for a sequence to be masked
    #[arg(long, default_value = "false")]
    pub experimental: bool,   
}