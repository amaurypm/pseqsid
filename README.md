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
    -t, --threads <THREADS>      Number of threads to use. 0 use all available threads [default: 0]
    -h, --help                   Print help information
    -V, --version                Print version information
```
At least one option from `-i`, `-s`, and/or `-n` must be supplied for the program to do something.

## Notes
*At this moment just the  Normalized Similarity Score calculation is not yet implemented.*
More information will be supplied when the implementation of all features progress.

## Installation
`cargo install pseqsid`

Or you can download and build the crate yourself.

## Examples
pseqsid -i `msa_protein.fasta`
