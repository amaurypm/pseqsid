# pseqsid
Calculates pairwise sequence identity, similarity and normalized similarity score of proteins in a multiple sequence alignment.

## Usage
```
USAGE:
    pseqsid [OPTIONS] <MSA>

ARGS:
    <MSA>    Multiple Sequence Alignment file

OPTIONS:
    -i, --identity               Calculate pairwise sequence identity
    -s, --similarity             Calculate pairwise sequence similarity
    -n, --nss                    Calculate pairwise sequence Normalized Similarity Score
    -l, --length <LENGTH>        Sequence length to be use for identity and similarity calculations
                                 [default: smallest] [possible values: smallest, mean, largest,
                                 alignment]
    -g, --grouping <GROUPING>    Similarity amino acid grouping definition file. A default one is
                                 created if required and not provided
    -m, --matrix <MATRIX>        Type of matrix to be used for Normalized Similarity Score [default:
                                 blosum62] [possible values: blosum62, pam250, gonnet]
    -p, --po <PO>                Gap opening penalty (Po) [default: 10.0]
    -e, --pe <PE>                Gap extending penalty (Pe) [default: 0.5]
    -t, --threads <THREADS>      Number of threads to use. 0 use all available threads [default: 0]
    -h, --help                   Print help information
    -V, --version                Print version information
```

## About the argument and options
**`<MSA>`**&emsp;is the filepath of a protein multiple sequence alignment in FASTA format. Only standard amino acids and gap (`-`) are accepted. An alignment of DNA sequences will be interpreted as proteins with a bunch of Ala, Gly, Thr and Cys residues, giving meaningless similarity and Normalized Similarity Score (NSS) values.

**`-i`**&emsp;**`--identity`**&emsp;activates the calculation of pairwise identity between the entries in the alignment. Sequence identity is defined as the percentage of identical residues in equivalent positions with respect to the sequence length. Which sequence length type is used is defined by **`-l`** option. Gaps-only columns are ignored for all sequence length types, except for `alignment`.

**`-l`**&emsp;**`--length`**&emsp;sequence length to be used. If the given sequence pair has gaps there are four possible different sequence length values to be used: the smallest, the average (mean), largest and the alignment length (this last considers gap-only columns, that are ignored for the previous three). Which type of length to select will depend on the intended use of your analysis. For example, if you have a hypothetical protein `prot_short` and you are using a template protein `template_long` to model `prot_short` structure, and `prot_short` has shorter sequence than `template_long`, then using `-l` `smallest` will tell you the identity value you need to asses how feasible you homology model could be. In other situation, if you want to use the identity value to know how similar any given pair of sequence is in your alignment, you should probably use the `-l` `mean` or `alignment`.

**`-s`**&emsp;**`--similarity`**&emsp;activates the calculation of pairwise similarity between the entries in the alignment. Sequence similarity is defined as the percentage of identical or similar residues in equivalent positions with respect to the sequence length. Most of what was explained above for identity applies to similarity too. Amino acid similarity groups can be defined in a file and provided with the `-g` option.

**`-g`**&emsp;**`--grouping`**&emsp;filepath to the similarity amino acid grouping file. If none is provided and `-s` option is given, then a default one, named *default_aa_similarity_groups.txt* is created and used. This file format is simple: `aa group name`: `single_letter_aa_names`. `#` symbol comments out the rest of the line.  Each group name can be defined just once and any amino acid can belong to only one group (or none). Only standard amino acids are accepted. Each group must have at least two amino acids. Use default file *default_aa_similarity_groups.txt* as a template for custom grouping definitions if needed.

**`-n`**&emsp;**`--nss`**&emsp;activates the calculation of pairwise Normalized Similarity Score, using the formula:

$S = \frac{(\sum M_{ij} - oP_o -eP_e)(\sum M_{ii} + \sum M_{jj})}{2\sum M_{ii}\sum M_{jj}}$

Where:
$M_{ij}$ are log-odds values for each pair of aligned amino acids *ij* obtained from substitution matrices. Which substitution matrix to use is defined by the `-m` option.
$M_{ii}$ and $M_{jj}$ are the log-odds values corresponding to conserving the given residue for the first and the second sequence in the pair, respectively.
$o$ is the number of gap openings in the sequence pair.
$P_o$ is the gap opening penalty, as defined by `-p` option.
$e$ is the total number of gaps in the sequence pair. Gap-only columns are ignored.
$P_e$ gap extension penalty, as defined by `-e` option. 

**`-m`**&emsp;**`--matrix`**&emsp;substitution matrix to be used for normalized similarity score calculation. Values corresponding to ambiguous amino acid definitions (such as `B`, `Z`, `X`) are ignored and a new column and row is added for gaps `-` with value zero. All three matrices are public, and available from multiple sources, in case you want to inspect them.

**`-p`**&emsp;**`--po`**&emsp;gap opening penalty to be used for normalized similarity score calculation.

**`-e`**&emsp;**`--pe`**&emsp;gap extension penalty to be used for normalized similarity score calculation.

## Output
Identity, similarity and NSS matrices are saved as CSV files, which filenames are derived from `<MSA>`, appending a proper substring. As these are symmetric matrices only the lower half is saved. Diagonals are ignored as identity and similarity of each protein sequence with itself is 100%, and NSS is 1.
All results are also saved in a pairwise tabular format (a TSV file named by the input FASTA file plus _table.tsv) to facilitate integration with potential downstream tasks.

## Installation

### Using [snap](https://snapcraft.io/):

`snap install pseqsid`

### Using [cargo](https://www.rust-lang.org/tools/install):

`cargo install pseqsid`

Or you can download and build the crate yourself.
