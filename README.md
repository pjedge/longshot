# longshot

Longshot is a variant calling tool for diploid genomes using long error prone reads such as Pacific Biosciences (PacBio) SMRT and Oxford Nanopore Technologies (ONT). It takes as input an aligned BAM/CRAM file and outputs a phased VCF file with variants and haplotype information. It can also genotype and phase input VCF files. It can output haplotype-separated BAM files that can be used for downstream analysis. Currently, it only calls single nucleotide variants (SNVs), but it can genotype indels if they are given in an input VCF.

## citation
If you use Longshot, please cite the publication:

[Edge, P. and Bansal, V., 2019. Longshot enables accurate variant calling in diploid genomes from single-molecule long read sequencing. Nature communications, 10(1), pp.1-10.](https://www.nature.com/articles/s41467-019-12493-y)

## supported operating systems
Longshot has been tested using Ubuntu 16.04 and 18.04, CentOS 6.6, Manjaro Linux 17.1.11, and Mac OS 10.14.2 Mojave.
It should work on any linux-based system that has Rust and Cargo installed.

## dependencies

* rust >= 1.40.0
* zlib >= 1.2.11
* xz >= 5.2.3
* clangdev >= 7.0.1
* gcc >= 7.3.0
* libc-dev
* make
* various rust dependencies (automatically managed by cargo)

(older versions may work but have not been tested)
## installation

### installation using Bioconda

It is recommended to install Longshot using [Bioconda](https://bioconda.github.io/):
```
conda install longshot
```
This method supports Linux and Mac.
If you do not have Bioconda, you can install it with these steps:
First, install Miniconda (or Anaconda). Miniconda can be installed using the
scripts [here](https://docs.conda.io/en/latest/miniconda.html). 

The Bioconda channel can then be added using these commands:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
### manual installation using apt for dependencies (Ubuntu 18.04)
If you are using Ubuntu 18.04, you can install the dependencies using apt. Then, the Rust cargo package manager is used to compile Longshot. 
```
sudo apt-get install cargo zlib1g-dev xz-utils \
         libclang-dev clang cmake build-essential curl git  # install dependencies 
git clone https://github.com/pjedge/longshot          # clone the Longshot repository
cd longshot                                           # change directory
cargo install --path .                                # install Longshot
export PATH=$PATH:/home/$USER/.cargo/bin              # add cargo binaries to path
```
Installation should take around 4 minutes on a typical desktop machine and will use between 400 MB (counting cargo) and 1.2 GB (counting all dependencies) of disk space.
It is recommended to add the line ```export PATH=$PATH:/home/$USER/.cargo/bin``` to the end of your ```~/.bashrc``` file so that the longshot binary is in the PATH for future shell sessions.

## usage:
After installation, execute the longshot binary as so:
```
$ longshot [FLAGS] [OPTIONS] --bam <BAM/CRAM> --ref <FASTA> --out <VCF>
```

## execution on an example dataset
The directory ```example_data``` contains a simulated toy dataset that can be used to test out Longshot:
- Reference genome containing 3 contigs each with length 200 kb (```example_data/genome.fa```)
- 30x coverage simulated pacbio reads generated using [SimLoRD](https://bitbucket.org/genomeinformatics/simlord/) (```example_data/pacbio_reads_30x.bam```)
- The 714 "true" variants for validation (```example_data/ground_truth_variants.vcf```)

Run Longshot on the example data as so:
```
longshot --bam example_data/pacbio_reads_30x.bam --ref example_data/genome.fa --out example_data/longshot_output.vcf
```

Execution should take around 5 to 10 seconds on a typical desktop machine. The output can be compared to ```ground_truth_variants.vcf``` for accuracy.

## command line options
```
$ longshot --help

Longshot: variant caller (SNVs) for long-read sequencing data 

USAGE:
    longshot [FLAGS] [OPTIONS] --bam <BAM/CRAM> --ref <FASTA> --out <VCF>

FLAGS:
    -A, --auto_max_cov        Automatically calculate mean coverage for region and set max coverage to mean_coverage +
                              5*sqrt(mean_coverage). (SLOWER)
    -S, --stable_alignment    Use numerically-stable (logspace) pair HMM forward algorithm. Is significantly slower but
                              may be more accurate. Tests have shown this not to be necessary for highly error prone
                              reads (PacBio CLR).
    -F, --force_overwrite     If output files (VCF or variant debug directory) exist, delete and overwrite them.
    -x, --max_alignment       Use max scoring alignment algorithm rather than pair HMM forward algorithm.
    -n, --no_haps             Don't call HapCUT2 to phase variants.
        --output-ref          print reference genotypes (non-variant), use this option only in combination with -v
                              option.
    -h, --help                Prints help information
    -V, --version             Prints version information

OPTIONS:
    -b, --bam <BAM>                            sorted, indexed BAM file with error-prone reads (CRAM files also supported)
    -f, --ref <FASTA>                          indexed FASTA reference that BAM file is aligned to
    -o, --out <VCF>                            output VCF file with called variants.
    -r, --region <string>                      Region in format <chrom> or <chrom:start-stop> in which to call variants
                                               (1-based, inclusive).
    -v, --potential_variants <VCF>             Genotype and phase the variants in this VCF instead of using pileup
                                               method to find variants. NOTES: VCF must be gzipped and tabix indexed or
                                               contain contig information. Use with caution because excessive false
                                               potential variants can lead to inaccurate results. Every variant is used
                                               and only the allele fields are considered -- Genotypes, filters,
                                               qualities etc are ignored. Indel variants will be genotyped but not
                                               phased. Structural variants (length > 50 bp) are currently not supported.
    -O, --out_bam <BAM>                        Write new bam file with haplotype tags (HP:i:1 and HP:i:2) for reads
                                               assigned to each haplotype, any existing HP and PS tags are removed
    -c, --min_cov <int>                        Minimum coverage (of reads passing filters) to consider position as a
                                               potential SNV. [default: 6]
    -C, --max_cov <int>                        Maximum coverage (of reads passing filters) to consider position as a
                                               potential SNV. [default: 8000]
    -q, --min_mapq <int>                       Minimum mapping quality to use a read. [default: 20]
    -a, --min_allele_qual <float>              Minimum estimated quality (Phred-scaled) of allele observation on read to
                                               use for genotyping/haplotyping. [default: 7.0]
    -y, --hap_assignment_qual <float>          Minimum quality (Phred-scaled) of read->haplotype assignment (for read
                                               separation). [default: 20.0]
    -Q, --potential_snv_cutoff <float>         Consider a site as a potential SNV if the original PHRED-scaled QUAL
                                               score for 0/0 genotype is below this amount (a larger value considers
                                               more potential SNV sites). [default: 20.0]
    -e, --min_alt_count <int>                  Require a potential SNV to have at least this many alternate allele
                                               observations. [default: 3]
    -E, --min_alt_frac <float>                 Require a potential SNV to have at least this fraction of alternate
                                               allele observations. [default: 0.125]
    -L, --hap_converge_delta <float>           Terminate the haplotype/genotype iteration when the relative change in
                                               log-likelihood falls below this amount. Setting a larger value results in
                                               faster termination but potentially less accurate results. [default:
                                               0.0001]
    -l, --anchor_length <int>                  Length of indel-free anchor sequence on the left and right side of read
                                               realignment window. [default: 6]
    -m, --max_snvs <int>                       Cut off variant clusters after this many variants. 2^m haplotypes must be
                                               aligned against per read for a variant cluster of size m. [default: 3]
    -W, --max_window <int>                     Maximum "padding" bases on either side of variant realignment window
                                               [default: 50]
    -I, --max_cigar_indel <int>                Throw away a read-variant during allelotyping if there is a CIGAR indel
                                               (I/D/N) longer than this amount in its window. [default: 20]
    -B, --band_width <Band width>              Minimum width of alignment band. Band will increase in size if sequences
                                               are different lengths. [default: 20]
    -D, --density_params <string>              Parameters to flag a variant as part of a "dense cluster". Format
                                               <n>:<l>:<gq>. If there are at least n variants within l base pairs with
                                               genotype quality >=gq, then these variants are flagged as "dn" [default:
                                               10:500:50]
    -s, --sample_id <string>                   Specify a sample ID to write to the output VCF [default: SAMPLE]
        --hom_snv_rate <float>                 Specify the homozygous SNV Rate for genotype prior estimation [default:
                                               0.0005]
        --het_snv_rate <float>                 Specify the heterozygous SNV Rate for genotype prior estimation [default:
                                               0.001]
        --ts_tv_ratio <float>                  Specify the transition/transversion rate for genotype grior estimation
                                               [default: 0.5]
    -P, --strand_bias_pvalue_cutoff <float>    Remove a variant if the allele observations are biased toward one strand
                                               (forward or reverse) according to Fisher's exact test. Use this cutoff
                                               for the two-tailed P-value. [default: 0.01]
    -d, --variant_debug_dir <path>             write out current information about variants at each step of algorithm to
                                               files in this directory
```

## usage examples
Call variants with default parameters:
```
longshot --bam pacbio.bam --ref ref.fa --out output.vcf
```
Call variants for chromosome 1 only using the automatic max coverage cutoff:
```
longshot -A -r chr1 --bam pacbio.bam --ref ref.fa --out output.vcf
```
Call variants in a 500 kb region and then output the reads into ```reads.bam``` using a haplotype assignment threshold of 30:
```
longshot -r chr1:1000000-1500000 -y 30 -O reads.bam --bam pacbio.bam --ref ref.fa --out output.vcf
```
If a read has an assigned haplotype, it will get a tag `HP:i:1` or `HP:i:2` and tag `PS:i:x` where `x` is a phase set number of the variants it covers.

## important considerations
- It is highly recommended to use reads with at least 30x coverage.
- It is recommended to process chromosomes separately using the ```--region``` option.
- Longshot has only been tested using data from humans. Results may vary with organisms with significantly higher or lower SNV rate.
- It is important to set a reasonable max read coverage cutoff (```-C``` option) to filter out sites coinciding with genomic features such as CNVs which can be problematic for variant calling. If the ```-A``` option is used, Longshot will estimate the mean read coverage and set the max coverage to ```mean_cov+5*sqrt(mean_cov)```, which we have found to be a reasonable filter in practice for humans.
- CNVs and mapping issues can result in dense clusters of false positive SNVs. Longshot will attempt to find clusters like this and mark them as "dn" in the FILTER field. The ```--density_params``` option is used to control which variants are flagged as "dn". The default parameters have been found to be effective for human sequencing data, but this option may need to be tweaked for other organisms with SNV rates significantly different from human.
- Oxford Nanopore Technology (ONT) SMS reads are now officially supported. It is recommended to use the default ```--strand_bias_pvalue_cutoff``` of 0.01 for ONT reads, since this option filters out false SNV sites prior to variant calling.

## installation troubleshooting

### older version of Rust
Check that the Rust version is 1.30.0 or higher:
```
rustc --version
```
If not, update Rust using this command:
```
rustup update
```


### linker errors
For example:
```
error: linking with `cc` failed: exit code: 1
...
...
...
= note: Non-UTF-8 output: /usr/bin/ld: /home/pedge/temp/longshot/target/release/build/longshot-347f3774e75b380c/out/libhapcut2.a(common.o)(.text.fprintf_time+0x81): unresolvable H\x89\\$\xe8H\x89l$\xf0H\x89\xf3L\x89d$\xf8H\x83\xec\x18H\x8bG\x10H\x89\xfdI\x89\xd4H\x89\xd6H\x8b;\xffPxH\x8bE\x10I\x8dt$\x08H\x8b{\x08\xffPxH\x8bE\x10H\x8b{\x10I\x8dt$\x10H\x8b\x1c$H\x8bl$\x08L\x8bd$\x10H\x8b@xH\x83\xc4\x18\xff\xe0f\x90H\x89\\$\xe8H\x89l$\xf0H\x89\xfbL\x89d$\xf8H\x83\xec\x18H\x8bG\x10I\x89\xd4H\x89\xf5H\x89\xf7\xffPhI\x89\x04$H\x8bC\x10H\x8d}\x08\xffPhH\x8b\x1c$I\x89D$\x08H\x8bl$\x08L\x8bd$\x10H\x83\xc4\x18\xc3\x0f\x1f relocation against symbol `time@@GLIBC_2.2.5\'\n/usr/bin/ld: BFD version 2.20.51.0.2-5.42.el6 20100205 internal error, aborting at reloc.c line 443 in bfd_get_reloc_size\n\n/usr/bin/ld: Please report this bug.\n\ncollect2: ld returned 1 exit status\n
...
...
...
```
Your system may have multiple versions of your linker that are causing a conflict. Rustc may be calling to a different or old version of the linker. In this case, specify the linker (in linux, gcc) as follows:
```
rustc -vV
```
Note the build target after "host: ", i.e. "x86_64-unknown-linux-gnu".
```
mkdir .cargo
nano .cargo/config
```
edit the config file to have these contents:
```
[target.<target-name>]
linker = "</path/to/linker>"
```
for example,
```
[target.x86_64-unknown-linux-gnu]
linker = "/opt/gnu/gcc/bin/gcc"
```
then,
```
cargo clean
cargo build --release
```
