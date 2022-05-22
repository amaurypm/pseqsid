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
    //description: String, // Not interested in the description.
    sequence: Vec<char>,
}

impl FastaSeq {
    /// Creates a FastaSeq instance. 
    /// The parameters must be checked for errors before passing them.
    pub fn new(identifier: String, data: String) -> FastaSeq {
        assert!(identifier.len() > 0);
        assert!(data.len() > 0);

        FastaSeq {
            identifier,
            //description,
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


/// Sequence length to be used
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug)]
pub enum SequenceLength {
    Smallest,
    Mean,
    Largest,
    Alignment,
}

/// Type of matrix to be used for Normalized Similarity Score
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug)]
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

    pub fn similarity(&self, length_mode: SequenceLength, aa_sim_groups: &HashMap<String, HashSet<char>>) -> f64 {
        let sequence_length = self.sequence_length(length_mode);
        assert!(self.first_sequence.whole_len() == self.second_sequence.whole_len());
        let mut similarity_count: u32= 0;
        for i in 0..self.first_sequence.whole_len() {
            let mut aa_pair_set = HashSet::new();
            aa_pair_set.insert(self.first_sequence.sequence()[i]);
            aa_pair_set.insert(self.second_sequence.sequence()[i]);

            if self.first_sequence.sequence()[i] == self.second_sequence.sequence()[i] || aa_sim_groups.iter().any(|(k,v)| aa_pair_set.is_subset(&v)) {
                if self.first_sequence.sequence()[i] != '-' {
                    similarity_count += 1;
                } else if let SequenceLength::Alignment = length_mode  {
                    similarity_count += 1;                        
                }
            }
        }
        (similarity_count as f64 / sequence_length) * 100.0

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
    /// Parse Fasta file containing a multiple protein sequences alignment
    pub fn from_file(filepath: &str) -> Result<MultipleSequenceAlignment, Box<dyn Error>> {
        let file = File::open(filepath)?;
        let reader = BufReader::new(file);
        let mut seq_vec: Vec<FastaSeq> = Vec::new();
        let mut identifier = String::new();
        let mut data = String::new();
        let mut open_seq = false;
        let mut some_seq = false;  

        for line in reader.lines() {
            let line_str = line?;
            // eprintln!("{}\tline_str = {}", "DEBUG".magenta().bold(), &line_str);
            // eprintln!("{}\tline_str.trim() = {}", "DEBUG".magenta().bold(), &line_str.trim());

            if line_str.trim().starts_with('>') {
                if !some_seq {
                    some_seq = true;
                }

                if open_seq {
                    eprintln!("{}\tIgnoring line {}", "WARN".yellow().bold(), line_str.italic());
                    continue;
                } else {
                    if data.len() > 0 {
                        if !is_seq_ok(&data) {
                            eprintln!("{}\tSequence {} contains non-standard amino acids. Ignoring entry", "WARN".yellow().bold(), identifier.italic());
                            data = String::from("");
                        } else if identifier.len() == 0 {
                            eprintln!("{}\tData {} has an empty identifier. Ignoring data", "WARN".yellow().bold(), data.italic());
                            data = String::from("");
                        } else {
                            // eprintln!("{}\tidentifier = {}", "DEBUG".magenta().bold(), &identifier);
                            // eprintln!("{}\tdescription = {}", "DEBUG".magenta().bold(), &description);
                            // eprintln!("{}\tdata = {}", "DEBUG".magenta().bold(), &data);

                            seq_vec.push(FastaSeq::new(identifier.clone(), data.clone()));

                            data = String::from("");
                        }                       

                    }

                    let description = String::from(line_str.trim().trim_start_matches('>'));
                    identifier = match description.split_whitespace().next() {
                        Some(id) => String::from(id),
                        None => {
                            eprintln!("{}\tSequence with no identifier detected. Ignoring line", "WARN".yellow().bold());
                            continue;
                        },
                    };
                    open_seq = true;                   

                }
                                
                
            } else {
                    if !some_seq {
                        eprintln!("{}\tIgnoring line {}", "WARN".yellow().bold(), line_str.italic());
                        continue;
                    } else {
                        data += &line_str.trim();
                        open_seq = false;
                    }
                }
        }

        if data.len() == 0 {
            eprintln!("{}\tSequence {} is empty. Ignoring entry", "WARN".yellow().bold(), identifier.italic());
        } else if !is_seq_ok(&data) {
            eprintln!("{}\tSequence {} contains non-standard amino acids. Ignoring entry", "WARN".yellow().bold(), identifier.italic());
        } else if identifier.len() == 0 {
            eprintln!("{}\tData {} has an empty identifier. Ignoring data", "WARN".yellow().bold(), data.italic());
            data = String::from("");
        } else {
            // eprintln!("{}\tidentifier = {}", "DEBUG".magenta().bold(), &identifier);
            // eprintln!("{}\tdescription = {}", "DEBUG".magenta().bold(), &description);
            // eprintln!("{}\tdata = {}", "DEBUG".magenta().bold(), &data);

            seq_vec.push(FastaSeq::new(identifier.clone(), data.clone()));
        }
        
        if same_len(&seq_vec) {
            Ok(MultipleSequenceAlignment {
                sequences: seq_vec,
            })
        } else {
            Err(Box::new(MSAError::new("all entries in the MSA must have same length (including gaps)")))
        }

        

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

    pub fn write_matrix(&self, filepath: &str, value_map: HashMap<(usize, usize), f64>) -> Result<(), Box<dyn Error>> {
        let mut output = File::create(filepath)?;
        let mut line = String::from("\"names");

        let sequences = self.get_sequences();

        // Write header
        for identifier in sequences.iter().map(|s| s.identifier()) {
            line += "\",\"";
            line += identifier;            
        }
        line += "\"";
        write!(output, "{}\n", line)?;

        // Write the matrix, with rownames
        for i in 0..sequences.len() {
            line = "\"".to_string() + sequences[i].identifier() + "\"";

            for j in 0..sequences.len() {
                if i > j {
                    let identity = match value_map.get(&(i, j)) {
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

//Functions

/// `true` is sequence is composed by standard amino acids or gap only
fn is_seq_ok(sequence: &str) -> bool {
    let sequence_set : HashSet<char> = HashSet::from_iter(sequence.chars());
    let std_aa_set = StdAAnGap::create();

    sequence_set.is_subset(&std_aa_set) && sequence_set.len() > 0
}

/// `true` is char set is composed by standard amino acids (not gaps)
fn is_seq_set_ok(sequence_set: &HashSet<char>) -> bool {
    let mut std_aa_set = StdAAnGap::create();
    std_aa_set.remove(&'-');

    sequence_set.is_subset(&std_aa_set) && sequence_set.len() > 0
}

/// Check if all FastaSeq instances in the vector have the same entry length.
/// This count gaps too.
fn same_len(seq_vec: &Vec<FastaSeq>) -> bool {
    let length_set: HashSet<usize> = seq_vec.iter().map(|s| s.whole_len()).collect();
    length_set.len() == 1
}

/// Create output file path
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

pub fn write_default_aa_sim_group() -> Result<String, Box<dyn Error>> {
    let output_filepath = "default_aa_similarity_groups.txt";
    let mut output = File::create(output_filepath)?;

    write!(output, "# Default amino acid similarity groups definition file.\n")?;
    write!(output, "#\n")?;
    write!(output, "# File format is simple:\n")?;
    write!(output, "# aa group name: single_letter_aa_names\n")?;
    write!(output, "# '#' symbol comments out the rest of the line\n")?;
    write!(output, "#\n")?;
    write!(output, "# Each group name can be defined just once and\n")?;
    write!(output, "# any amino acid can belong to only one group (or none).\n")?;
    write!(output, "# Only standard amino acids will be accepted.\n")?;
    write!(output, "# Each group must have at least two amino acids.\n")?;
    write!(output, "# You can modify this file as you wish, as long as you\n")?;
    write!(output, "# comply with the previous definition rules.\n")?;
    write!(output, "#\n")?;
    write!(output, "# group: amino acids\n")?;
    write!(output, "aromatic: FWY\n")?;
    write!(output, "aliphatic: VIL\n")?;
    write!(output, "charged positive: RKH\n")?;
    write!(output, "charged negative: DE\n")?;
    write!(output, "small: ST\n")?;
    write!(output, "polar: NQ\n")?;
    write!(output, "\n")?;

    Ok(String::from(output_filepath))

}

/// Process the file containing amino acid group definitions for similarity calculation
fn process_aa_sim_group_file(filepath: &str) -> Result<HashMap<String, HashSet<char>>,  Box<dyn Error>> {
    let mut aa_sim_groups = HashMap::new();
    let file = File::open(filepath)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let mut line_str = line?;

        if line_str.contains("#") {
            if let Some((s1, s2)) = line_str.split_once('#') {
                line_str = s1.to_string();
            }
        }

        if line_str.len() > 0 {
            match line_str.trim().split_once(':') {
                Some((group_name, aa_str)) => {
                    if group_name.trim().len() == 0 {
                        return Err(Box::new(MSAError::new(&format!("empty group name in file {}", filepath.italic()))));
                    }

                    if aa_sim_groups.contains_key(group_name) {
                        return Err(Box::new(MSAError::new(&format!("group {} is defined more than once in file {}", group_name.bold(), filepath.italic()))));
                    }

                    let mut aa_set: HashSet<char> = HashSet::from_iter(aa_str.chars());
                    aa_set.remove(&' ');                    

                    if aa_set.len() < 2 {
                        return Err(Box::new(MSAError::new(&format!("less than two amino acids listed for group {} in file {}", group_name.bold(), filepath.italic()))));
                    }

                    if !is_seq_set_ok(&aa_set) {
                        return Err(Box::new(MSAError::new(&format!("group {} in file {} contains non standard amino acids", group_name.bold(), filepath.italic())))); 
                    }                    

                    if aa_sim_groups.iter().map(|(k,v)| aa_set.intersection(&v).collect()).any(|s: HashSet<&char>| s.len() > 0) {
                        return Err(Box::new(MSAError::new(&format!("amino acids belong to more than one group in file {}", filepath.italic()))));
                    }

                    aa_sim_groups.insert(group_name.to_string(), aa_set);
                },
                None => return Err(Box::new(MSAError::new(&format!("format error in file {}: ':' missing in declaration line", filepath.italic())))),
            }
        }
    }

    if aa_sim_groups.len() > 0 {
        Ok(aa_sim_groups)

    } else {
        Err(Box::new(MSAError::new(&format!("no group definition found in file {}: ':' missing in declaration line", filepath))))

    }

}

/// Function to be called from main
pub fn run(msa_filepath: &str, identity: bool, similarity: bool, nss: bool, length_mode: SequenceLength, aa_grouping_filepath: &str, matrix: Matrix, threads: usize) -> Result<(), Box<dyn Error>> {
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
        msa.write_matrix(&output_filepath, identity_map)?;
    }

    if similarity {
        let aa_sim_groups = process_aa_sim_group_file(aa_grouping_filepath)?;

        let similarity_vec: Vec<((usize, usize), f64)> = seqpair_vec.par_iter().map(|p| (p.index(), p.similarity(length_mode, &aa_sim_groups))).collect();

        let mut similarity_map: HashMap<(usize, usize), f64> = HashMap::new();
        for (index, similarity) in similarity_vec {
            similarity_map.entry(index).or_insert(similarity);
        }

        let output_filepath = output_path(msa_filepath, OutputMode::Similarity);
        msa.write_matrix(&output_filepath, similarity_map)?;
    }

    if nss {
        println!("{}\tNormalized similarity score not yet implemented", "INFO".yellow().bold());

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
