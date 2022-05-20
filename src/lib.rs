/// Multiple sequence alignment

use std::collections::HashSet;
use colored::Colorize;
use clap::ArgEnum;
mod data;
use data::StdAAnGap;

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

fn is_seq_ok(sequence: &str, std_aa_set: &HashSet<char>) -> bool {
    let sequence_set : HashSet<char> = HashSet::from_iter(sequence.chars());
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

/// A single pair of sequences.
/// Sequence comparison is done here.
struct SeqPair<'a>{
    first_index: usize, // Index of the first sequence, as it is in a MultipleSequenceAlignment instance.
    second_index: usize, // Index of the second sequence. first_index < second_index, always.
    first_sequence: &'a FastaSeq,
    second_sequence: &'a FastaSeq,
}

impl<'a> SeqPair<'a> {
    /// Creates and initialize a SeqPair instance.
    pub fn new(first_index: usize, second_index: usize, first_sequence: &'a FastaSeq, second_sequence: &'a FastaSeq) -> SeqPair<'a> {
        assert!(first_index < second_index);

        SeqPair {
            first_index,
            second_index,
            first_sequence,
            second_sequence,                
        }
    }

    pub fn identity(&mut self, length_mode: SequenceLength) -> f64 {
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
    
}

struct MultipleSequenceAlignment {
    sequences: Vec<FastaSeq>,
}




#[cfg(test)]
mod tests {
//     use std::collections::HashSet;
//    use super::*;

//     #[test]
//     fn seq_is_ok() {
//         let seq = "ACDFRGSQWERTH";
//         let std_aa_set: HashSet<char> = HashSet::from_iter(data::STD_AA_GAP.chars());
//         assert!(msa::is_seq_ok(&seq, &std_aa_set));        
//     }

//     #[test]
//     fn seq_is_not_ok() {
//         let seq = "ZXCASQWERTGHBNMIKLOP";
//         let std_aa_set: HashSet<char> = HashSet::from_iter(data::STD_AA_GAP.chars());
//         assert!(!msa::is_seq_ok(&seq, &std_aa_set));        
//     }
}
