/// Multiple sequence alignment

use std::collections::{HashSet, HashMap};
use colored::Colorize;
use clap::ArgEnum;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use rayon::prelude::*;
use std::path::Path;
use std::ffi::OsStr;
mod data;
use data::StdAAnGap;
mod error;
use error::MSAError;

struct FastaSeq {
    identifier: String,
    description: String,
    sequence: Vec<char>,
}

impl FastaSeq {
    /// Creates a FastaSeq instance. 
    /// The parameters must be checked for errors before passing them.
    pub fn new(identifier: String, description: String, data: String) -> FastaSeq {
        assert!(identifier.len() > 0);
        assert!(data.len() > 0);

        FastaSeq {
            identifier,
            description,
            sequence: data.chars().collect(),
        }
    }

    /// Returns entry length, including gaps.
    pub fn whole_len(&self) -> usize {
        self.sequence.len()
    }

    /// Returns sequence length, excluding gaps.
    pub fn len(&self) -> usize {
        self.sequence.iter().filter(|&c| *c != '-').count()
    }

    pub fn sequence(&self) -> &Vec<char> {
        &self.sequence
    }

    pub fn identifier(&self) -> &str {
        &self.identifier
    }

}

fn is_seq_ok(sequence: &str) -> bool {
    let sequence_set : HashSet<char> = HashSet::from_iter(sequence.chars());
    let std_aa_set = StdAAnGap::create();

    sequence_set.is_subset(&std_aa_set) && sequence_set.len() > 0
}

/// Sequence length to be used
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
pub enum SequenceLength {
    Smallest,
    Mean,
    Largest,
    Alignment,
}

/// Type of matrix to be used for Normalized Similarity Score
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
pub enum Matrix {
    BLOSUM62,
    PAM250,
    GONNET,
}

enum OutputMode {
    Identity,
    Similarity,
    NSS,
}

/// A single pair of sequences.
/// Sequence comparison is done here.
struct SeqPair<'a>{
    first_index: usize, // Index of the first sequence, as it is in a MultipleSequenceAlignment instance.
    second_index: usize, // Index of the second sequence. first_index > second_index, always.
    first_sequence: &'a FastaSeq,
    second_sequence: &'a FastaSeq,
}

impl<'a> SeqPair<'a> {
    /// Creates and initialize a SeqPair instance.
    pub fn new(first_index: usize, second_index: usize, first_sequence: &'a FastaSeq, second_sequence: &'a FastaSeq) -> SeqPair<'a> {
        assert!(first_index > second_index);

        SeqPair {
            first_index,
            second_index,
            first_sequence,
            second_sequence,                
        }
    }

    pub fn identity(&self, length_mode: SequenceLength) -> f64 {
        let sequence_length = self.sequence_length(length_mode);
        assert!(self.first_sequence.whole_len() == self.second_sequence.whole_len());
        let mut identity_count: u32= 0;
        for i in 0..self.first_sequence.whole_len() {
            if self.first_sequence.sequence()[i] == self.second_sequence.sequence()[i] {
                if self.first_sequence.sequence()[i] != '-' {
                    identity_count += 1;
                } else if let SequenceLength::Alignment = length_mode  {
                    identity_count += 1;                        
                }
            }
        }
        (identity_count as f64 / sequence_length) * 100.0

    }

    fn sequence_length(&self, length_mode: SequenceLength) -> f64 {
        match length_mode {
            SequenceLength::Alignment => {
                assert!(self.first_sequence.whole_len() == self.second_sequence.whole_len());
                self.first_sequence.whole_len() as f64
            },
            SequenceLength::Largest => {
                if self.first_sequence.len() > self.second_sequence.len() {
                    self.first_sequence.len() as f64
                } else {
                    self.second_sequence.len() as f64
                }
            },
            SequenceLength::Mean => (self.first_sequence.len() as f64 + self.second_sequence.len() as f64)/2.0,
            SequenceLength::Smallest => {
                if self.first_sequence.len() < self.second_sequence.len() {
                    self.first_sequence.len() as f64
                } else {
                    self.second_sequence.len() as f64
                }
            }            
        }
    }

    fn index(&self) -> (usize, usize) {
        (self.first_index, self.second_index)
    }
    
}

struct MultipleSequenceAlignment {
    sequences: Vec<FastaSeq>,
}

impl MultipleSequenceAlignment {
    pub fn from_file(filepath: &str) -> Result<MultipleSequenceAlignment, Box<dyn Error>> {
        let file = File::open(filepath)?;
        let reader = BufReader::new(file);
        let mut seq_vec: Vec<FastaSeq> = Vec::new();
        let mut identifier = String::new();
        let mut description = String::new();
        let mut data = String::new();
        let mut open_seq = false;
        let mut some_seq = false;  

        for line in reader.lines() {
            let line_str = line?;

            if line_str.trim().starts_with('>') {
                if !some_seq {
                    some_seq = true;
                }

                if open_seq {
                    eprintln!("{} Ignoring line {}", "WARN".yellow().bold(), line_str.italic());
                    continue;
                } else {
                    if data.len() == 0 {
                        eprintln!("{} Sequence {} is empty. Ignoring entry.", "WARN".yellow().bold(), identifier.italic());
                        continue;
                    }

                    if !is_seq_ok(&data) {
                        eprintln!("{} Sequence {} contains non-standard amino acids. Ignoring entry.", "WARN".yellow().bold(), identifier.italic());
                        continue;
                    }

                    seq_vec.push(FastaSeq::new(identifier.clone(), description.clone(), data.clone()));

                }
                                
                description = String::from(line_str.trim().trim_start_matches('>'));
                identifier = match description.split_whitespace().next() {
                    Some(id) => String::from(id),
                    None => {
                        eprintln!("{} Sequence with no identifier detected. Ignoring line.", "WARN".yellow().bold());
                        continue;
                    },
                };
                open_seq = true;
            } else {
                    if !some_seq {
                        eprintln!("{} Ignoring line {}", "WARN".yellow().bold(), line_str.italic());
                        continue;
                    } else {
                        data += &line_str.trim();
                        open_seq = false;
                    }
                }
        }

        if data.len() == 0 {
            eprintln!("{} Sequence {} is empty. Ignoring entry.", "WARN".yellow().bold(), identifier.italic());
        } else if !is_seq_ok(&data) {
            eprintln!("{} Sequence {} contains non-standard amino acids. Ignoring entry.", "WARN".yellow().bold(), identifier.italic());
        } else {
            seq_vec.push(FastaSeq::new(identifier.clone(), description.clone(), data.clone()));
        }
        
        Ok(MultipleSequenceAlignment {
            sequences: seq_vec,
        })

    }

    pub fn get_seq_pairs(&self) -> Vec<SeqPair> {
        // At least two sequences are required.
        // This most be enforced at creation time.
        assert!(self.sequences.len() > 1);
        let mut seq_pair_vec: Vec<SeqPair> = Vec::new();
        for i in 0..self.sequences.len() {
            for j in 0..i {
                seq_pair_vec.push(SeqPair::new(i, j, &self.sequences[i], &self.sequences[j]));
            }
        }
        seq_pair_vec
    }

    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    fn get_sequences(&self) -> &Vec<FastaSeq> {
        &self.sequences
    }

    pub fn write_identity_matrix(&self, filepath: &str, identity_map: HashMap<(usize, usize), f64>) -> Result<(), Box<dyn Error>> {
        let mut output = File::create(filepath)?;
        let mut line = String::from("\"names\"");

        let sequences = self.get_sequences();

        // Write header
        for identifier in sequences.iter().map(|s| s.identifier()) {
            line += "\",\"";
            line += identifier;            
        }
        write!(output, "{}\n", line)?;

        // Write the matrix, with rownames
        for i in 0..sequences.len() {
            line = "\"".to_string() + sequences[i].identifier() + "\"";

            for j in 0..sequences.len() {
                if i > j {
                    let identity = match identity_map.get(&(i, j)) {
                        Some(x) => x,
                        None => {
                            return Err(Box::new(MSAError::new("missing identity value for sequence pair")));
                        },                    
                    };
                    line += &format!(",{:.2}", identity);

                } else {
                    line += ","
                }
            }
            write!(output, "{}\n", line)?;
        }

        Ok(())

    }

}

fn output_path(msa_filepath: &str, mode: OutputMode) -> String {
    let output_filename = match Path::new(msa_filepath).file_name() {
       Some(s) => s,
       None => OsStr::new(""),
    };

    let output_filename = match Path::new(output_filename).file_stem() {
        Some(s) => s,
        None => OsStr::new(""),
    };

    let output_filename = match output_filename.to_str() {
        Some(s) => s.to_string(),
        None => String::from(""),
    };    

    match mode {
        OutputMode::Identity => output_filename + "_identity.csv",
        OutputMode::Similarity => output_filename + "_similarity.csv",
        OutputMode::NSS => output_filename + "_nss.csv",
    }
}

pub fn run(msa_filepath: &str, identity: bool, similarity: bool, nss: bool, sim_def_filepath: &str, length_mode: SequenceLength, matrix: Matrix, threads: usize) -> Result<(), Box<dyn Error>> {
    // Initialize rayon.
    // This allows to control the number of threads to use.
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global()?;

    let msa = MultipleSequenceAlignment::from_file(msa_filepath)?;
    if msa.len() < 2 {
        return Err(Box::new(MSAError::new("less than two sequences in MSA file")));
    }

    let seqpair_vec: Vec<SeqPair> = msa.get_seq_pairs();

    if identity {
        let identity_vec: Vec<((usize, usize), f64)> = seqpair_vec.par_iter().map(|p| (p.index(), p.identity(length_mode))).collect();

        let mut identity_map: HashMap<(usize, usize), f64> = HashMap::new();
        for (index, identity) in identity_vec {
            identity_map.entry(index).or_insert(identity);
        }

        let output_filepath = output_path(msa_filepath, OutputMode::Identity);
        msa.write_identity_matrix(&output_filepath, identity_map)?;
    }

    if similarity {
        println!("Similarity related stuff goes here");
    }

    if nss {
        println!("Normalized similarity score related stuff goes here");

    }

    Ok(())

}

// #[cfg(test)]
// mod tests {
// //     use std::collections::HashSet;
// //    use super::*;

// //     #[test]
// //     fn seq_is_ok() {
// //         let seq = "ACDFRGSQWERTH";
// //         let std_aa_set: HashSet<char> = HashSet::from_iter(data::STD_AA_GAP.chars());
// //         assert!(msa::is_seq_ok(&seq, &std_aa_set));        
// //     }

// //     #[test]
// //     fn seq_is_not_ok() {
// //         let seq = "ZXCASQWERTGHBNMIKLOP";
// //         let std_aa_set: HashSet<char> = HashSet::from_iter(data::STD_AA_GAP.chars());
// //         assert!(!msa::is_seq_ok(&seq, &std_aa_set));        
// //     }
// }
