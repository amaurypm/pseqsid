/// Multiple sequence alignment
mod msa {
    use std::collections::HashSet;
    use colored::Colorize;
    use crate::data::STD_AA_GAP;

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

}

/// Data associated to amino acids
mod data {
    //Complete
    pub const STD_AA_GAP: &str = "ACDEFGHIKLMNPQRSTVWY-";
}


#[cfg(test)]
mod tests {
//     use std::collections::HashSet;
    use super::*;
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
