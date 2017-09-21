# reaper
The REAlign error PronE Reads variant caller.

Reaper is a prototype SNV caller for long error prone reads such as Pacific Biosciences SMRT.

Reaper requires Rust 1.2 which comes with the cargo package manager. You can install rust with the
following command:
```
curl https://sh.rustup.rs -sSf | sh
```
To build reaper, type:
```
cargo build --release
```
This will compile the code with optimizations in release mode, and the binary will be
in ```target/release/reaper```. It will automatically install dependencies for reaper,
such rust-bio and rust-htslib.

To run unit tests, type:
```
cargo test
```

Usage:
```
$ ./target/release/reaper --help

Reaper (REAlign error PronE Reads) 0.1
Peter Edge <edge.peterj@gmail.com>
SNV caller for Third-Generation Sequencing reads

USAGE:
    reaper [FLAGS] [OPTIONS] --bam <BAM> --ref <FASTA> --out <VCF>

FLAGS:
    -f, --fast_alignment    Use non-numerically stable alignment algorithm. Is significantly faster but may be less accurate or have unexpected behaviour.
    -h, --help              Prints help information
    -V, --version           Prints version information

OPTIONS:
    -b, --bam <BAM>                  sorted, indexed BAM file with error-prone reads
    -r, --ref <FASTA>                indexed fasta reference that BAM file is aligned to
    -o, --out <VCF>                  output VCF file with called variants.
    -a, --min_alt_count <int>        Minimum number of occurrences of a variant base in pileup to consider it as a potential SNV. [default: 2]
    -A, --min_alt_frac <float>       Minimum fraction of variant base in pileup to consider it as a potential SNV. [default: 0.125]
    -c, --min_cov <int>              Minimum coverage to consider position as a potential SNV. [default: 5]
    -C, --max_cov <int>              Maximum coverage to consider position as a potential SNV.
    -q, --min_mapq <int>             Minimum mapping quality to use a read. [default: 30]
    -k, --anchor_k <int>             A filter for low-complexity anchor sequences. A valid anchor must have no duplicates in the kmers that overlap it. [default: 5]
    -d, --snv_distance <int>         SNVs separated by distance less than d will be considered together and all their possible haplotypes will be aligned against. [default: 20]
    -m, --max_snvs <int>             Cut off short haplotypes after this many SNVs. 2^m haplotypes must be aligned against per read for a variant cluster of size m. [default: 5]
    -w, --min_window <int>           Ignore a variant/short haplotype if the realignment window for read or reference is smaller than w bases. [default: 15]
    -W, --max_window <int>           Ignore a variant/short haplotype if the realignment window for read or reference is larger than w bases. [default: 200]
    -B, --band_width <Band width>    Minimum width of alignment band. Band will increase in size if sequences are different lengths. [default: 20]
    -l, --anchor_length <int>        Length of indel-free anchor sequence on the left and right side of read realignment window. [default: 10]
```
