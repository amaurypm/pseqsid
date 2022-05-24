/// Data associated with amino acids

use std::collections::HashSet;
use clap::ArgEnum;

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

/// Type of matrix to be used for Normalized Similarity Score
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug)]
pub enum Matrix {
    BLOSUM62,
    PAM250,
    GONNET,
}
 
/// Substitution matrix
/// It will create a modified BLOSUM62, PAM250, or a GONNET
/// depending on the parameters of the new method.
/// This matrices will exclude values for ambiguous residues,
/// such as B, Z and X
pub struct SubstitutionMatrix {
    /// The actual matrix
    data: Vec<Vec<f64>>,

    /// Residues as index of each row and column. The order matters.
    names: Vec<char>,

}

impl SubstitutionMatrix {
    /// Create a new SubstitutionMatrix instance
    /// The instance created will contains the data and name order of
    /// BLOSUM62, PAM250 or GONNET depending on the Matrix value passed
    pub fn new(mat_type: Matrix) -> SubstitutionMatrix {
        match mat_type {
            Matrix::BLOSUM62 => {
                let names: Vec<char> = "ARNDCQEGHILKMFPSTWYV-".chars().collect();
                let data = Vec::from([
                    vec![4.0, -1.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0,  1.0,  0.0, -3.0, -2.0,  0.0, 0.0],
                    vec![-1.0, 5.0, 0.0, -2.0, -3.0, 1.0, 0.0, -2.0, 0.0, -3.0, -2.0, 2.0, -1.0, -3.0, -2.0, -1.0, -1.0, -3.0, -2.0, -3.0, 0.0],
                    vec![-2.0, 0.0, 6.0, 1.0, -3.0, 0.0, 0.0, 0.0, 1.0, -3.0, -3.0, 0.0, -2.0, -3.0, -2.0, 1.0, 0.0, -4.0, -2.0, -3.0, 0.0],
                    vec![-2.0, -2.0, 1.0, 6.0, -3.0, 0.0, 2.0, -1.0, -1.0, -3.0, -4.0, -1.0, -3.0, -3.0, -1.0, 0.0, -1.0, -4.0, -3.0, -3.0, 0.0],
                    vec![0.0, -3.0, -3.0, -3.0, 9.0, -3.0, -4.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0, -3.0, -1.0, -1.0, -2.0, -2.0, -1.0, 0.0],
                    vec![-1.0, 1.0, 0.0, 0.0, -3.0, 5.0, 2.0, -2.0, 0.0, -3.0, -2.0, 1.0, 0.0, -3.0, -1.0, 0.0, -1.0, -2.0, -1.0, -2.0, 0.0],
                    vec![-1.0, 0.0, 0.0, 2.0, -4.0, 2.0, 5.0, -2.0, 0.0, -3.0, -3.0, 1.0, -2.0, -3.0, -1.0, 0.0, -1.0, -3.0, -2.0, -2.0, 0.0],
                    vec![0.0, -2.0, 0.0, -1.0, -3.0, -2.0, -2.0, 6.0, -2.0, -4.0, -4.0, -2.0, -3.0, -3.0, -2.0, 0.0, -2.0, -2.0, -3.0, -3.0, 0.0],
                    vec![-2.0, 0.0, 1.0, -1.0, -3.0, 0.0, 0.0, -2.0, 8.0, -3.0, -3.0, -1.0, -2.0, -1.0, -2.0, -1.0, -2.0, -2.0, 2.0, -3.0, 0.0],
                    vec![-1.0, -3.0, -3.0, -3.0, -1.0, -3.0, -3.0, -4.0, -3.0, 4.0, 2.0, -3.0, 1.0, 0.0, -3.0, -2.0, -1.0, -3.0, -1.0, 3.0, 0.0],
                    vec![-1.0, -2.0, -3.0, -4.0, -1.0, -2.0, -3.0, -4.0, -3.0, 2.0, 4.0, -2.0, 2.0, 0.0, -3.0, -2.0, -1.0, -2.0, -1.0, 1.0, 0.0],
                    vec![-1.0, 2.0, 0.0, -1.0, -3.0, 1.0, 1.0, -2.0, -1.0, -3.0, -2.0, 5.0, -1.0, -3.0, -1.0, 0.0, -1.0, -3.0, -2.0, -2.0, 0.0],
                    vec![-1.0, -1.0, -2.0, -3.0, -1.0, 0.0, -2.0, -3.0, -2.0, 1.0, 2.0, -1.0, 5.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, 1.0, 0.0],
                    vec![-2.0, -3.0, -3.0, -3.0, -2.0, -3.0, -3.0, -3.0, -1.0, 0.0, 0.0, -3.0, 0.0, 6.0, -4.0, -2.0, -2.0, 1.0, 3.0, -1.0, 0.0],
                    vec![-1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -1.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -4.0, 7.0, -1.0, -1.0, -4.0, -3.0, -2.0, 0.0],
                    vec![1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -2.0, 0.0, -1.0, -2.0, -1.0, 4.0, 1.0, -3.0, -2.0, -2.0, 0.0],
                    vec![0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0, 1.0, 5.0, -2.0, -2.0, 0.0, 0.0],
                    vec![-3.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -2.0, -3.0, -1.0, 1.0, -4.0, -3.0, -2.0, 11.0, 2.0, -3.0, 0.0],
                    vec![-2.0, -2.0, -2.0, -3.0, -2.0, -1.0, -2.0, -3.0, 2.0, -1.0, -1.0, -2.0, -1.0, 3.0, -3.0, -2.0, -2.0, 2.0, 7.0, -1.0, 0.0],
                    vec![0.0, -3.0, -3.0, -3.0, -1.0, -2.0, -2.0, -3.0, -3.0, 3.0, 1.0, -2.0, 1.0, -1.0, -2.0, -2.0, 0.0, -3.0, -1.0, 4.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                ]);

                SubstitutionMatrix{
                    names: names,
                    data: data,
                }
            },
            Matrix::PAM250 => {
                let names: Vec<char> = "ARNDCQEGHILKMFPSTWYV-".chars().collect();
                let data = Vec::from([
                    vec![2.0, -2.0, 0.0, 0.0, -2.0, 0.0, 0.0, 1.0, -1.0, -1.0, -2.0, -1.0, -1.0, -3.0, 1.0, 1.0, 1.0, -6.0, -3.0, 0.0, 0.0],
                    vec![-2.0, 6.0, 0.0, -1.0, -4.0, 1.0, -1.0, -3.0, 2.0, -2.0, -3.0, 3.0, 0.0, -4.0, 0.0, 0.0, -1.0, 2.0, -4.0, -2.0, 0.0],
                    vec![0.0, 0.0, 2.0, 2.0, -4.0, 1.0, 1.0, 0.0, 2.0, -2.0, -3.0, 1.0, -2.0, -3.0, 0.0, 1.0, 0.0, -4.0, -2.0, -2.0, 0.0],
                    vec![0.0, -1.0, 2.0, 4.0, -5.0, 2.0, 3.0, 1.0, 1.0, -2.0, -4.0, 0.0, -3.0, -6.0, -1.0, 0.0, 0.0, -7.0, -4.0, -2.0, 0.0],
                    vec![-2.0, -4.0, -4.0, -5.0, 12.0, -5.0, -5.0, -3.0, -3.0, -2.0, -6.0, -5.0, -5.0, -4.0, -3.0, 0.0, -2.0, -8.0, 0.0, -2.0, 0.0],
                    vec![0.0, 1.0, 1.0, 2.0, -5.0, 4.0, 2.0, -1.0, 3.0, -2.0, -2.0, 1.0, -1.0, -5.0, 0.0, -1.0, -1.0, -5.0, -4.0, -2.0, 0.0],
                    vec![0.0, -1.0, 1.0, 3.0, -5.0, 2.0, 4.0, 0.0, 1.0, -2.0, -3.0, 0.0, -2.0, -5.0, -1.0, 0.0, 0.0, -7.0, -4.0, -2.0, 0.0],
                    vec![1.0, -3.0, 0.0, 1.0, -3.0, -1.0, 0.0, 5.0, -2.0, -3.0, -4.0, -2.0, -3.0, -5.0, 0.0, 1.0, 0.0, -7.0, -5.0, -1.0, 0.0],
                    vec![-1.0, 2.0, 2.0, 1.0, -3.0, 3.0, 1.0, -2.0, 6.0, -2.0, -2.0, 0.0, -2.0, -2.0, 0.0, -1.0, -1.0, -3.0, 0.0, -2.0, 0.0],
                    vec![-1.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -3.0, -2.0, 5.0, 2.0, -2.0, 2.0, 1.0, -2.0, -1.0, 0.0, -5.0, -1.0, 4.0, 0.0],
                    vec![-2.0, -3.0, -3.0, -4.0, -6.0, -2.0, -3.0, -4.0, -2.0, 2.0, 6.0, -3.0, 4.0, 2.0, -3.0, -3.0, -2.0, -2.0, -1.0, 2.0, 0.0],
                    vec![-1.0, 3.0, 1.0, 0.0, -5.0, 1.0, 0.0, -2.0, 0.0, -2.0, -3.0, 5.0, 0.0, -5.0, -1.0, 0.0, 0.0, -3.0, -4.0, -2.0, 0.0],
                    vec![-1.0, 0.0, -2.0, -3.0, -5.0, -1.0, -2.0, -3.0, -2.0, 2.0, 4.0, 0.0, 6.0, 0.0, -2.0, -2.0, -1.0, -4.0, -2.0, 2.0, 0.0],
                    vec![-3.0, -4.0, -3.0, -6.0, -4.0, -5.0, -5.0, -5.0, -2.0, 1.0, 2.0, -5.0, 0.0, 9.0, -5.0, -3.0, -3.0, 0.0, 7.0, -1.0, 0.0],
                    vec![1.0, 0.0, 0.0, -1.0, -3.0, 0.0, -1.0, 0.0, 0.0, -2.0, -3.0, -1.0, -2.0, -5.0, 6.0, 1.0, 0.0, -6.0, -5.0, -1.0, 0.0],
                    vec![1.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, -3.0, 0.0, -2.0, -3.0, 1.0, 2.0, 1.0, -2.0, -3.0, -1.0, 0.0],
                    vec![1.0, -1.0, 0.0, 0.0, -2.0, -1.0, 0.0, 0.0, -1.0, 0.0, -2.0, 0.0, -1.0, -3.0, 0.0, 1.0, 3.0, -5.0, -3.0, 0.0, 0.0],
                    vec![-6.0, 2.0, -4.0, -7.0, -8.0, -5.0, -7.0, -7.0, -3.0, -5.0, -2.0, -3.0, -4.0, 0.0, -6.0, -2.0, -5.0, 17.0, 0.0, -6.0, 0.0],
                    vec![-3.0, -4.0, -2.0, -4.0, 0.0, -4.0, -4.0, -5.0, 0.0, -1.0, -1.0, -4.0, -2.0, 7.0, -5.0, -3.0, -3.0, 0.0, 10.0, -2.0, 0.0],
                    vec![0.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -1.0, -2.0, 4.0, 2.0, -2.0, 2.0, -1.0, -1.0, -1.0, 0.0, -6.0, -2.0, 4.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                ]);

                SubstitutionMatrix{
                    names: names,
                    data: data,
                }
            },
            Matrix::GONNET => {
                let names: Vec<char> = "CSTPAGNDEQHRKMILVFYW-".chars().collect();
                let data = Vec::from([
                    vec![12.0, 0.0, 0.0, -3.0, 0.0, -2.0, -2.0, -3.0, -3.0, -2.0, -1.0, -2.0, -3.0, -1.0, -1.0, -2.0, 0.0, -1.0, 0.0, -1.0, 0.0],
                    vec![0.0, 2.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -2.0, -2.0, -1.0, -3.0, -2.0, -3.0, 0.0],
                    vec![0.0, 2.0, 2.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, 0.0, -2.0, -2.0, -4.0, 0.0],
                    vec![-3.0, 0.0, 0.0, 8.0, 0.0, -2.0, -1.0, -1.0, 0.0, 0.0, -1.0, -1.0, -1.0, -2.0, -3.0, -2.0, -2.0, -4.0, -3.0, -5.0, 0.0],
                    vec![0.0, 1.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, -1.0, -1.0, -1.0, 0.0, -2.0, -2.0, -4.0, 0.0],
                    vec![-2.0, 0.0, -1.0, -2.0, 0.0, 7.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -4.0, -4.0, -4.0, -3.0, -5.0, -4.0, -4.0, 0.0],
                    vec![-2.0, 1.0, 0.0, -1.0, 0.0, 0.0, 4.0, 2.0, 1.0, 1.0, 1.0, 0.0, 1.0, -2.0, -3.0, -3.0, -2.0, -3.0, -1.0, -4.0, 0.0],
                    vec![-3.0, 0.0, 0.0, -1.0, 0.0, 0.0, 2.0, 5.0, 3.0, 1.0, 0.0, 0.0, 0.0, -3.0, -4.0, -4.0, -3.0, -4.0, -3.0, -5.0, 0.0],
                    vec![-3.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 3.0, 4.0, 2.0, 0.0, 0.0, 1.0, -2.0, -3.0, -3.0, -2.0, -4.0, -3.0, -4.0, 0.0],
                    vec![-2.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 1.0, 2.0, 3.0, 1.0, 2.0, 2.0, -1.0, -2.0, -2.0, -2.0, -3.0, -2.0, -3.0, 0.0],
                    vec![-1.0, 0.0, 0.0, -1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 1.0, 6.0, 1.0, 1.0, -1.0, -2.0, -2.0, -2.0, 0.0, 2.0, -1.0, 0.0],
                    vec![-2.0, 0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 2.0, 1.0, 5.0, 3.0, -2.0, -2.0, -2.0, -2.0, -3.0, -2.0, -2.0, 0.0],
                    vec![-3.0, 0.0, 0.0, -1.0, 0.0, -1.0, 1.0, 0.0, 1.0, 2.0, 1.0, 3.0, 3.0, -1.0, -2.0, -2.0, -2.0, -3.0, -2.0, -4.0, 0.0],
                    vec![-1.0, -1.0, -1.0, -2.0, -1.0, -4.0, -2.0, -3.0, -2.0, -1.0, -1.0, -2.0, -1.0, 4.0, 2.0, 3.0, 2.0, 2.0, 0.0, -1.0, 0.0],
                    vec![-1.0, -2.0, -1.0, -3.0, -1.0, -4.0, -3.0, -4.0, -3.0, -2.0, -2.0, -2.0, -2.0, 2.0, 4.0, 3.0, 3.0, 1.0, -1.0, -2.0, 0.0],
                    vec![-2.0, -2.0, -1.0, -2.0, -1.0, -4.0, -3.0, -4.0, -3.0, -2.0, -2.0, -2.0, -2.0, 3.0, 3.0, 4.0, 2.0, 2.0, 0.0, -1.0, 0.0],
                    vec![0.0, -1.0, 0.0, -2.0, 0.0, -3.0, -2.0, -3.0, -2.0, -2.0, -2.0, -2.0, -2.0, 2.0, 3.0, 2.0, 3.0, 0.0, -1.0, -3.0, 0.0],
                    vec![-1.0, -3.0, -2.0, -4.0, -2.0, -5.0, -3.0, -4.0, -4.0, -3.0, 0.0, -3.0, -3.0, 2.0, 1.0, 2.0, 0.0, 7.0, 5.0, 4.0, 0.0],
                    vec![0.0, -2.0, -2.0, -3.0, -2.0, -4.0, -1.0, -3.0, -3.0, -2.0, 2.0, -2.0, -2.0, 0.0, -1.0, 0.0, -1.0, 5.0, 8.0, 4.0, 0.0],
                    vec![-1.0, -3.0, -4.0, -5.0, -4.0, -4.0, -4.0, -5.0, -4.0, -3.0, -1.0, -2.0, -4.0, -1.0, -2.0, -1.0, -3.0, 4.0, 4.0, 14.0, 0.0],
                    vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                ]);

                SubstitutionMatrix{
                    names: names,
                    data: data,
                }

            },
        }

    }

    /// Return the value r,c of the matrix, if index in range
    /// If not return None
    fn get_by_index(&self, r: usize, c: usize) -> Option<f64> {
        if r < self.names.len() && c < self.names.len() {
            Some(self.data[r][c])

        } else {
            None
        }
    }

    fn char_index(&self, aa: char) -> Option<usize> {
        for i in 0..self.names.len() {
            if aa == self.names[i] {
                return Some(i);
            }
        }

        None
    }

    /// Retrieve the values in the matrix converting the amino
    /// acid names into appropriated indices.
    fn get_by_char(&self, aar: char, aac:char) -> Option<f64> {        

        let r = match self.char_index(aar) {
            Some(x) => x,
            None => return None
        };

        let c = match self.char_index(aac) {
            Some(x) => x,
            None => return None
        };

        self.get_by_index(r, c)
    }

    /// Use get_by_char to retrieve matrix values.
    pub fn get(&self, aar: char, aac:char) -> Option<f64> {
        self.get_by_char(aar, aac)
    }
}