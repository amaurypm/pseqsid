/// Data associated with amino acids

use std::collections::HashSet;

/// Standard amino acids and gap.
pub struct StdAAnGap;

impl StdAAnGap {
    /// Creates a HashSet of standard amino acids and gap as chars.
    pub fn create() -> HashSet<char> {
        let std_aa_gap = "ACDEFGHIKLMNPQRSTVWY-";
        let std_aa_set: HashSet<char> = HashSet::from_iter(std_aa_gap.chars());
        std_aa_set

    }
}

    