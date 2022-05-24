use clap::{AppSettings, Parser};
use colored::Colorize;
use std::process;
use pseqsid::{Matrix, SequenceLength};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
struct Cli {
    /// Multiple Sequence Alignment file
    msa: String,

    /// Calculate pairwise sequence identity
    #[clap(short, long)]
    identity: bool,

    /// Calculate pairwise sequence similarity
    #[clap(short, long)]
    similarity: bool,

    /// Calculate pairwise sequence Normalized Similarity Score
    #[clap(short, long)]
    nss: bool,

    /// Sequence length to be use for identity and similarity calculations
    #[clap(short, long, arg_enum, default_value_t = SequenceLength::Smallest)]
    length: SequenceLength,

    /// Similarity amino acid grouping definition file.
    /// A default one is created if required and not provided
    #[clap(short, long)]
    grouping: Option<String>,

    /// Type of matrix to be used for Normalized Similarity Score
    #[clap(short, long, arg_enum, default_value_t = Matrix::BLOSUM62)]
    matrix: Matrix,

    /// Number of threads to use.
    /// 0 use all available threads
    #[clap(short, long, default_value_t = 0)]
    threads: usize,


}

fn main() {
    let cli = Cli::parse();

    // eprintln!("{}\tmsa= {}","DEBUG".magenta().bold(), cli.msa);
    // eprintln!("{}\tidentity= {}","DEBUG".magenta().bold(), cli.identity);
    // eprintln!("{}\tsimilarity= {}","DEBUG".magenta().bold(), cli.similarity);
    // eprintln!("{}\tnss= {}","DEBUG".magenta().bold(), cli.nss);
    // eprintln!("{}\tlength= {:?}","DEBUG".magenta().bold(), cli.length);
    // eprintln!("{}\tgrouping= {:?}","DEBUG".magenta().bold(), cli.grouping);
    // eprintln!("{}\tmatrix= {:?}","DEBUG".magenta().bold(), cli.matrix);
    // eprintln!("{}\tthreads= {:?}","DEBUG".magenta().bold(), cli.threads);

    if !cli.identity && !cli.similarity && !cli.nss {
        println!("{}\tNothing to do", "INFO".green().bold());
        println!("{}\tPlease set -i, -s and/or -n options", "INFO".green().bold());
        println!("{}\tSee --help for more information", "INFO".green().bold());
        process::exit(0);
    }

    let aa_grouping_filepath = match cli.grouping {
        Some(s) => s.to_string(),
        None => match pseqsid::write_default_aa_sim_group() {
            Ok(s) => s,
            Err(e) => {
                eprintln!("{}\t{}", "ERROR".red().bold(), e);
                process::exit(0);
            },
        }
    };  

    if let Err(e) = pseqsid::run(&cli.msa, cli.identity, cli.similarity, cli.nss, cli.length, &aa_grouping_filepath, cli.matrix, cli.threads) {
        eprintln!("{}\t{}", "ERROR".red().bold(), e);
    }
}
