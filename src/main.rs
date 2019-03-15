//! Longshot
//!
//! Longshot is a diploid SNV caller for single molecule sequencing (SMS) reads such as PacBio
//!
//! Author: Peter Edge
//!
//! Contact: edge.peterj@gmail.com

#![allow(dead_code)]
// `error_chain!` can recurse deeply
#![recursion_limit = "1024"]

// external crates
extern crate bio;
extern crate chrono;
extern crate clap;
extern crate core;
extern crate rand;
extern crate rust_htslib;
#[macro_use]
extern crate error_chain;
extern crate fishers_exact;
extern crate hashbrown;

// import modules
mod call_genotypes;
mod call_potential_snvs;
mod errors;
mod estimate_alignment_parameters;
mod estimate_read_coverage;
mod extract_fragments; //mod extract_fragments_debug;
mod genotype_probs;
mod haplotype_assembly;
mod print_output;
mod realignment;
mod spoa;
mod util;
mod variants_and_fragments;

// use declarations
use bio::stats::{LogProb, PHREDProb, Prob};
use call_genotypes::*;
use clap::{App, Arg};
use errors::*;
use estimate_alignment_parameters::estimate_alignment_parameters;
use estimate_read_coverage::calculate_mean_coverage;
use extract_fragments::ExtractFragmentParameters;
use fishers_exact::fishers_exact;
use genotype_probs::{Genotype, GenotypePriors};
use haplotype_assembly::*;
use print_output::{print_variant_debug, print_vcf};
use realignment::{AlignmentType, AlignmentParameters, EmissionProbs, TransitionProbs, baum_welch_non_numerically_stable};
use std::fs::create_dir;
use std::fs::remove_dir_all;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use util::*;
use util::{
    parse_flag, parse_positive_f64, parse_prob_into_logprob, parse_u32, parse_u8, parse_usize,
};
use variants_and_fragments::{VarList,VarFilter};

//use variants_and_fragments::parse_VCF_potential_variants;
//use haplotype_assembly::separate_reads_by_haplotype;
//use realignment::{AlignmentParameters, TransitionProbs, EmissionProbs};

/// The main function
///
/// The main function follows the [error-chain recommended practice](http://brson.github.io/2016/11/30/starting-with-error-chain)
/// It execute the run() function (which contains the entirety of the program logic)
/// If the program fails, grab the chain of errors incurred and print them with nonzero exit code (and any backtrace)
fn main() {
    if let Err(ref e) = run() {
        println!("error: {}", e);
        for e in e.iter().skip(1) {
            println!("caused by: {}", e);
        }
        if let Some(backtrace) = e.backtrace() {
            println!("backtrace: {:?}", backtrace);
        }
        std::process::exit(1);
    }
}

/// The run function
///
/// The run function contains the whole program logic as per [error-chain recommended practice](http://brson.github.io/2016/11/30/starting-with-error-chain)
///
/// # Outline of steps performed
/// 1. Retrieve and process command line arguments
/// 2. If --auto_max_cov is set, calculate mean read coverage and use it to estimate max read coverage
/// 3. Traverse the CIGARs in the bam file to estimate the HMM alignment parameters
/// 4. Use pileup-based approach to call potential SNVs
/// 5. Extract haplotype fragment information (alleles/qualities per site per read) from reads using Pair-HMM realignment
/// 6. Calculate genotypes for each site using the estimated alleles/qualities without using haplotype information
/// 7. Iteratively assemble haplotypes and refine genotypes
/// 8. Print output VCF
///
/// # Errors
/// - errors returned by any of the function calls deeper in the program
/// - if various command-line arguments are specified multiple times or with invalid values
/// - if a file/directory already exists and -F option isn't set (e.g. vcf output or vcf debug directory)
/// - input bam file isn't indexed
fn run() -> Result<()> {
    /***********************************************************************************************/
    // READ COMMAND LINE ARGUMENTS
    /***********************************************************************************************/

    {

        let v = vec!['A','C','T','G'];
        let w = vec!['A','C','C','G'];
        let alignment_parameters = AlignmentParameters {
            emission_probs: EmissionProbs {equal: 0.982, not_equal: 0.006, insertion: 1.0, deletion:1.0},
            transition_probs: TransitionProbs {match_from_match: 0.879,
                insertion_from_match: 0.08,
                deletion_from_match: 0.04,
                match_from_deletion: 0.9,
                deletion_from_deletion: 0.1,
                match_from_insertion: 0.8,
                insertion_from_insertion: 0.2}
        };

        baum_welch_non_numerically_stable(&v,&w,alignment_parameters,20);
        assert!(true);
    }

    eprintln!("");

    let input_args = App::new("Longshot")
        .version("0.3.2")
        .author("Peter Edge <edge.peterj@gmail.com>")
        .about("SNV caller for Third-Generation Sequencing reads")
        .arg(Arg::with_name("Input BAM")
                .short("b")
                .long("bam")
                .value_name("BAM")
                .help("sorted, indexed BAM file with error-prone reads")
                .display_order(10)
                .required(true)
                .takes_value(true))
        .arg(Arg::with_name("Input FASTA")
                .short("f")
                .long("ref")
                .value_name("FASTA")
                .help("indexed FASTA reference that BAM file is aligned to")
                .display_order(20)
                .required(true)
                .takes_value(true))
        .arg(Arg::with_name("Output VCF")
                .short("o")
                .long("out")
                .value_name("VCF")
                .help("output VCF file with called variants.")
                .display_order(30)
                .required(true)
                .takes_value(true))
        .arg(Arg::with_name("Region")
                .short("r")
                .long("region")
                .value_name("string")
                .help("Region in format <chrom> or <chrom:start-stop> in which to call variants (1-based, inclusive).")
                .display_order(40)
                //.required(true)
                .takes_value(true))
        /*.arg(Arg::with_name("Potential Variants VCF")
            .short("v")
            .long("potential_variants")
            .value_name("VCF")
            .help("Use the variants in this VCF as the potential variants instead of using pileup method. NOTE: every variant is used and only the allele fields are considered! Genotypes, filters, qualities etc are ignored!")
            .display_order(45)
            .takes_value(true))*/
        .arg(Arg::with_name("Haplotype Bam Prefix")
            .short("p")
            .long("hap_bam_prefix")
            .value_name("BAM")
            .help("Write haplotype-separated reads to 3 bam files using this prefix: <prefix>.hap1.bam, <prefix>.hap2.bam, <prefix>.unassigned.bam")
            .display_order(50))
        .arg(Arg::with_name("Auto max coverage")
            .short("A")
            .long("auto_max_cov")
            .help("Automatically calculate mean coverage for region and set max coverage to mean_coverage + 5*sqrt(mean_coverage). (SLOWER)")
            .display_order(75))
        .arg(Arg::with_name("Min coverage")
            .short("c")
            .long("min_cov")
            .value_name("int")
            .help("Minimum coverage (of reads passing filters) to consider position as a potential SNV.")
            .display_order(78)
            .default_value("6"))
        .arg(Arg::with_name("Max coverage")
                .short("C")
                .long("max_cov")
                .value_name("int")
                .help("Maximum coverage (of reads passing filters) to consider position as a potential SNV.")
                .display_order(80)
                .default_value("8000"))
        .arg(Arg::with_name("Min mapq")
                .short("q")
                .long("min_mapq")
                .value_name("int")
                .help("Minimum mapping quality to use a read.")
                .display_order(90)
                .default_value("30"))
        .arg(Arg::with_name("Min allele quality")
            .short("a")
            .long("min_allele_qual")
            .value_name("float")
            .help("Minimum estimated quality (Phred-scaled) of allele observation on read to use for genotyping/haplotyping.")
            .display_order(92)
            .default_value("7.0"))
        .arg(Arg::with_name("Haplotype assignment quality")
            .short("y")
            .long("hap_assignment_qual")
            .value_name("float")
            .help("Minimum quality (Phred-scaled) of read->haplotype assignment (for read separation).")
            .display_order(94)
            .default_value("20.0"))
        .arg(Arg::with_name("Potential SNV Cutoff")
            .long("potential_snv_cutoff")
            .short("Q")
            .value_name("float")
            .help("Consider a site as a potential SNV if the original PHRED-scaled QUAL score for 0/0 genotype is below this amount (a larger value considers more potential SNV sites).")
            .display_order(96)
            .default_value("20.0"))
        .arg(Arg::with_name("Potential SNV Min Alt Count")
            .long("min_alt_count")
            .short("e")
            .value_name("int")
            .help("Require a potential SNV to have at least this many alternate allele observations.")
            .display_order(97)
            .default_value("3"))
        .arg(Arg::with_name("Potential SNV Min Alt Fraction")
            .long("min_alt_frac")
            .short("E")
            .value_name("float")
            .help("Require a potential SNV to have at least this fraction of alternate allele observations.")
            .display_order(98)
            .default_value("0.125"))
        .arg(Arg::with_name("Haplotype Convergence Delta")
            .long("hap_converge_delta")
            .short("L")
            .value_name("float")
            .help("Terminate the haplotype/genotype iteration when the relative change in log-likelihood falls below this amount. Setting a larger value results in faster termination but potentially less accurate results.")
            .display_order(99)
            .default_value(&"0.0001"))
        .arg(Arg::with_name("Anchor length")
                .short("l")
                .long("anchor_length")
                .value_name("int")
                .help("Length of indel-free anchor sequence on the left and right side of read realignment window.")
                .display_order(100)
                .default_value("6"))
        .arg(Arg::with_name("Variant cluster max size")
                .short("m")
                .long("max_snvs")
                .value_name("int")
                .help("Cut off variant clusters after this many variants. 2^m haplotypes must be aligned against per read for a variant cluster of size m.")
                .display_order(130)
                .default_value("3"))
        .arg(Arg::with_name("Use POA")
            .short("M")
            .long("poa")
            .help("EXPERIMENTAL: Run the algorithm twice, using Partial-Order-Alignment on phased reads to find new candidate SNVs and Indels the second time.")
            .display_order(130))
        .arg(Arg::with_name("Max window padding")
                .short("W")
                .long("max_window")
                .value_name("int")
                .help("Maximum \"padding\" bases on either side of variant realignment window")
                .display_order(150)
                .default_value("50"))
        .arg(Arg::with_name("Max CIGAR indel")
                .short("I")
                .long("max_cigar_indel")
                .value_name("int")
                .default_value("20")
                .help("Throw away a read-variant during allelotyping if there is a CIGAR indel (I/D/N) longer than this amount in its window.")
                .display_order(151))
        .arg(Arg::with_name("Numerically stable alignment")
            .short("S")
            .long("stable_alignment")
            .help("Use numerically-stable (logspace) pair HMM forward algorithm. Is significantly slower but may be more accurate. Tests have shown this not to be necessary for highly error prone reads (PacBio CLR).")
            .display_order(161))
        .arg(Arg::with_name("Force overwrite")
            .short("F")
            .long("force_overwrite")
            .help("If output files (VCF or variant debug directory) exist, delete and overwrite them.")
            .display_order(164))
        .arg(Arg::with_name("Max alignment")
            .short("x")
            .long("max_alignment")
            .help("Use max scoring alignment algorithm rather than pair HMM forward algorithm.")
            .display_order(166))
        .arg(Arg::with_name("Band width")
                .short("B")
                .long("band_width")
                .help("Minimum width of alignment band. Band will increase in size if sequences are different lengths.")
                .display_order(170)
                .default_value("20"))
        .arg(Arg::with_name("Density parameters")
            .short("D")
            .long("density_params")
            .value_name("string")
            .help("Parameters to flag a variant as part of a \"dense cluster\". Format <n>:<l>:<gq>. \
                     If there are at least n variants within l base pairs with genotype quality >=gq, \
                     then these variants are flagged as \"dn\"")
            .display_order(172)
            .default_value("10:500:50"))
        .arg(Arg::with_name("Sample ID")
            .short("s")
            .long("sample_id")
            .value_name("string")
            .help("Specify a sample ID to write to the output VCF")
            .display_order(174)
            .default_value(&"SAMPLE"))
        .arg(Arg::with_name("Homozygous SNV Rate")
            .long("hom_snv_rate")
            .value_name("float")
            .help("Specify the homozygous SNV Rate for genotype prior estimation")
            .display_order(176)
            .default_value(&"0.0005"))
        .arg(Arg::with_name("Heterozygous SNV Rate")
            .long("het_snv_rate")
            .value_name("float")
            .help("Specify the heterozygous SNV Rate for genotype prior estimation")
            .display_order(178)
            .default_value(&"0.001"))
        .arg(Arg::with_name("Homozygous Indel Rate")
            .long("hom_indel_rate")
            .value_name("float")
            .help("Specify the homozygous Indel Rate for genotype prior estimation")
            .display_order(180)
            .hidden(true)
            .default_value(&"0.00005"))
        .arg(Arg::with_name("Heterozygous Indel Rate")
            .long("het_indel_rate")
            .value_name("float")
            .help("Specify the heterozygous Indel Rate for genotype prior estimation")
            .display_order(182)
            .hidden(true)
            .default_value(&"0.0001"))
        .arg(Arg::with_name("ts/tv Ratio")
            .long("ts_tv_ratio")
            .value_name("float")
            .help("Specify the transition/transversion rate for genotype grior estimation")
            .display_order(184)
            .default_value(&"0.5"))
        .arg(Arg::with_name("Strand Bias P-value cutoff")
            .short("P")
            .long("strand_bias_pvalue_cutoff")
            .value_name("float")
            .help("Remove a variant if the allele observations are biased toward one strand (forward or reverse) according to Fisher's exact test. Use this cutoff for the two-tailed P-value.")
            .display_order(185)
            .default_value(&"0.01"))
        .arg(Arg::with_name("Min Genotype Quality")
            .short("g")
            .long("min_gq")
            .value_name("float")
            .help("Minimum Genotype Quality for final genotypes")
            .display_order(187)
            .default_value(&"0.0"))
        .arg(Arg::with_name("No haplotypes")
                .short("n")
                .long("no_haps")
                .help("Don't call HapCUT2 to phase variants.")
                .display_order(190))
        .arg(Arg::with_name("Variant debug directory")
            .short("d")
            .long("variant_debug_dir")
            .value_name("path")
            .help("write out current information about variants at each step of algorithm to files in this directory")
            .display_order(210))
        .get_matches();

    // parse the input arguments and throw errors if inputs are invalid
    let bamfile_name = input_args
        .value_of("Input BAM")
        .chain_err(|| "Input BAM file not defined.")?
        .to_string();
    let fasta_file = input_args
        .value_of("Input FASTA")
        .chain_err(|| "Input FASTA file not defined.")?
        .to_string();
    let output_vcf_file = input_args
        .value_of("Output VCF")
        .chain_err(|| "Output VCF file not defined.")?
        .to_string();
    let interval: Option<GenomicInterval> =
        parse_region_string(input_args.value_of("Region"), &bamfile_name)?;
    let hap_bam_prefix: Option<&str> = input_args.value_of("Haplotype Bam Prefix");
    let force = parse_flag(&input_args, "Force overwrite")?;
    let use_poa = parse_flag(&input_args, "Use POA")?;
    let no_haps = parse_flag(&input_args, "No haplotypes")?;
    let min_mapq: u8 = parse_u8(&input_args, "Min mapq")?;
    let anchor_length: usize = parse_usize(&input_args, "Anchor length")?;
    let variant_cluster_max_size: usize = parse_usize(&input_args, "Variant cluster max size")?;
    let max_window_padding: usize = parse_usize(&input_args, "Max window padding")?;
    let max_cigar_indel: usize = parse_usize(&input_args, "Max CIGAR indel")?;
    let min_allele_qual: f64 = parse_nonnegative_f64(&input_args, "Min allele quality")?;
    let min_gq: f64 = parse_nonnegative_f64(&input_args, "Min Genotype Quality")?;

    let strand_bias_pvalue_cutoff: f64 =
        parse_nonnegative_f64(&input_args, "Strand Bias P-value cutoff")?;
    let hap_assignment_qual: f64 =
        parse_nonnegative_f64(&input_args, "Haplotype assignment quality")?;
    let ll_delta: f64 = parse_positive_f64(&input_args, "Haplotype Convergence Delta")?;
    let potential_snv_cutoff_phred = parse_positive_f64(&input_args, "Potential SNV Cutoff")?;
    let potential_snv_min_alt_count: usize =
        parse_usize(&input_args, "Potential SNV Min Alt Count")?;
    let potential_snv_min_alt_frac: f64 =
        parse_nonnegative_f64(&input_args, "Potential SNV Min Alt Fraction")?;
    let hom_snv_rate: LogProb = parse_prob_into_logprob(&input_args, "Homozygous SNV Rate")?;
    let het_snv_rate: LogProb = parse_prob_into_logprob(&input_args, "Heterozygous SNV Rate")?;
    let hom_indel_rate: LogProb = parse_prob_into_logprob(&input_args, "Homozygous Indel Rate")?;
    let het_indel_rate: LogProb = parse_prob_into_logprob(&input_args, "Heterozygous Indel Rate")?;
    let sample_name: String = input_args
        .value_of(&"Sample ID")
        .chain_err(|| "Sample ID not defined.")?
        .to_string();
    //let potential_variants_file: Option<&str> = input_args.value_of("Potential Variants VCF");

    // sanity checks on values that aren't covered by parsing functions
    ensure!(
        ll_delta < 1.0,
        format!("Haplotype Convergence Delta must be less than 1.0!")
    );

    // manipulations to get some of the option values into forms we want
    let max_p_miscall: f64 = *Prob::from(PHREDProb(min_allele_qual));
    let hap_max_p_misassign: f64 = *Prob::from(PHREDProb(hap_assignment_qual));
    let potential_snv_cutoff: LogProb = LogProb::from(PHREDProb(potential_snv_cutoff_phred));

    // if VCF file exists, throw error unless --force_overwrite option is set
    let vcf = Path::new(&output_vcf_file);
    ensure!(
        !vcf.is_file() || force,
        "Variant output file already exists. Rerun with -F option to force overwrite."
    );

    // ensure that BAM file is indexed
    let bai_str = bamfile_name.clone() + ".bai";
    ensure!(Path::new(&bai_str).is_file(), "BAM file must be indexed with samtools index. Index file should have same name as BAM file with .bai appended.");

    // ensure that FASTA file is indexed
    let fai_str = fasta_file.clone() + ".fai";
    ensure!(Path::new(&fai_str).is_file(), "FASTA reference file must be indexed with samtools faidx. Index file should have same name as FASTA file with .fai appended.");

    // check if variant debug directory exists
    // if it does, delete the directory if --force_overwrite option is set or throw an error
    let variant_debug_directory: Option<String> = match input_args
        .value_of("Variant debug directory")
    {
        Some(dir) => {
            let p = Path::new(&dir);
            if p.exists() {
                if force {
                    remove_dir_all(p).chain_err(|| "Error removing variant debug directory.")?;
                } else {
                    bail!("Variant debug directory already exists. Rerun with -F option to force overwrite.");
                }
            }
            create_dir(&dir).chain_err(|| "Error creating variant debug directory.")?;
            Some(dir.to_string())
        }
        None => None,
    };

    // multiply by 2.0 because internally we use this value as the probability of transition to
    // a single transversion base, not the combined probability of transversion to either one
    let ts_tv_ratio = 2.0 * parse_positive_f64(&input_args, "ts/tv Ratio")?;

    let dn_params = input_args
        .value_of("Density parameters")
        .chain_err(|| "Density parameters not defined.")?
        .split(":")
        .collect::<Vec<&str>>();

    if dn_params.len() != 3 {
        bail!(
            "Format for density params should be <n>:<l>:<gq>, with all 3 values being integers."
        );
    }

    let dn_count = dn_params[0].parse::<usize>().chain_err(|| {
        "Format for density params should be <n>:<l>:<gq>, with all 3 values being integers."
    })?;
    let dn_len = dn_params[1].parse::<usize>().chain_err(|| {
        "Format for density params should be <n>:<l>:<gq>, with all 3 values being integers."
    })?;
    let dn_gq = dn_params[2].parse::<usize>().chain_err(|| {
        "Format for density params should be <n>:<l>:<gq>, with all 3 values being integers."
    })?;

    let density_params = DensityParameters {
        n: dn_count,
        len: dn_len,
        gq: dn_gq as f64,
    };

    let alignment_type = match (
        parse_flag(&input_args, "Numerically stable alignment")?,
        parse_flag(&input_args, "Max alignment")?,
    ) {
        (false, false) => AlignmentType::ForwardAlgorithmNonNumericallyStable,
        (true, false) => AlignmentType::ForwardAlgorithmNumericallyStable,
        (false, true) => AlignmentType::ViterbiMaxScoringAlignment,
        (true, true) => {
            bail!(
                "Numerically stable alignment option and max alignment options are incompatible."
            );
        }
    };

    let band_width: usize = parse_usize(&input_args, "Band width")?;
    //let use_poa = parse_flag(&input_args, "Use POA");
    let min_cov: u32 = parse_u32(&input_args, "Min coverage")?;

    let max_cov: u32 = match parse_flag(&input_args, "Auto max coverage")? {
        false => {
            // manually assigned coverage cutoff from user
            parse_u32(&input_args, "Max coverage")?
        }
        true => {
            eprintln!(
                "{} Automatically determining max read coverage.",
                print_time()
            );
            eprintln!("{} Estimating mean read coverage...", print_time());
            let mean_coverage: f64 = calculate_mean_coverage(&bamfile_name, &interval)
                .chain_err(|| "Error calculating mean coverage for BAM file.")?;
            let calculated_max_cov =
                (mean_coverage as f64 + 5.0 * (mean_coverage as f64).sqrt()) as u32;
            eprintln!("{} Mean read coverage: {:.2}", print_time(), mean_coverage);

            calculated_max_cov
        }
    };

    if max_cov > 0 {
        eprintln!("{} Min read coverage set to {}.", print_time(), min_cov);
        eprintln!("{} Max read coverage set to {}.", print_time(), max_cov);
    } else {
        bail!("{} ERROR: Max read coverage set to 0.");
    }

    // we store the read IDs if we will be separating the reads by haplotype
    // we will compute sets holding the separated read IDs and then refer back to the original BAM
    // and write to separate files based on set membership
    let store_read_id = hap_bam_prefix != None;

    let extract_fragment_parameters = ExtractFragmentParameters {
        min_mapq,
        alignment_type,
        band_width,
        anchor_length,
        variant_cluster_max_size: variant_cluster_max_size,
        max_window_padding,
        max_cigar_indel,
        store_read_id,
    };


    eprintln!("{} Estimating alignment parameters...", print_time());
    let alignment_parameters = estimate_alignment_parameters(
        &bamfile_name,
        &fasta_file,
        &interval,
        min_mapq,
        max_cigar_indel as u32,
    )
    .chain_err(|| "Error estimating alignment parameters.")?;

    /*
    let alignment_parameters = AlignmentParameters {
        emission_probs: EmissionProbs {equal: 0.982, not_equal: 0.006, insertion: 1.0, deletion:1.0},
        transition_probs: TransitionProbs {match_from_match: 0.879,
            insertion_from_match: 0.06,
            deletion_from_match: 0.06,
            match_from_deletion: 0.85,
            deletion_from_deletion: 0.15,
            match_from_insertion: 0.85,
            insertion_from_insertion: 0.15}
    };
    */

    /***********************************************************************************************/
    // GET GENOTYPE PRIORS
    /***********************************************************************************************/

    let genotype_priors = GenotypePriors::new(
        hom_snv_rate,
        het_snv_rate,
        hom_indel_rate,
        het_indel_rate,
        ts_tv_ratio,
    )
    .chain_err(|| "Error estimating genotype priors.")?;

    /***********************************************************************************************/
    // FIND INITIAL SNVS WITH READ PILEUP
    /***********************************************************************************************/

    //let bam_file: String = "test_data/test.bam".to_string();
    eprintln!("{} Calling potential SNVs using pileup...", print_time());
    /*let mut varlist = match potential_variants_file {
        Some(file) => { parse_VCF_potential_variants(&file.to_string(), &bamfile_name) }
        None => { call_potential_snvs::call_potential_snvs(&bamfile_name,
                                                 &fasta_file,
                                                 &interval,
                                                 &genotype_priors,
                                                 min_cov,
                                                 max_cov,
                                                 min_mapq,
                                                 max_p_miscall,
                                                 alignment_parameters.ln()) }
    };*/
    let mut varlist = call_potential_snvs::call_potential_snvs(
        &bamfile_name,
        &fasta_file,
        &interval,
        &genotype_priors,
        min_cov,
        max_cov,
        potential_snv_min_alt_count,
        potential_snv_min_alt_frac,
        min_mapq,
        alignment_parameters.ln(),
        potential_snv_cutoff,
    )
    .chain_err(|| "Error calling potential SNVs.")?;

    print_variant_debug(
        &mut varlist,
        &interval,
        &variant_debug_directory,
        &"1.0.potential_SNVs.vcf",
        max_cov,
        &density_params,
        &sample_name,
    )?;

    eprintln!(
        "{} {} potential SNVs identified.",
        print_time(),
        varlist.lst.len()
    );

    if varlist.lst.len() == 0 {
        return Ok(());
    }

    /***********************************************************************************************/
    // EXTRACT FRAGMENT INFORMATION FROM READS
    /***********************************************************************************************/

    eprintln!(
        "{} Generating haplotype fragments from reads...",
        print_time()
    );
    let mut flist = extract_fragments::extract_fragments(
        &bamfile_name,
        &fasta_file,
        &mut varlist,
        &interval,
        extract_fragment_parameters,
        alignment_parameters,
    )
    .chain_err(|| "Error generating haplotype fragments from BAM reads.")?;

    // if we're printing out variant "debug" information, print out a fragment file to that debug directory
    match &variant_debug_directory {
        &Some(ref debug_dir) => {
            let ffn = match Path::new(&debug_dir).join(&"fragments.txt").to_str() {
                Some(s) => s.to_owned(),
                None => {
                    bail!("Invalid unicode provided for variant debug directory");
                }
            };
            // normally phase_variant is used to select which variants are heterozygous, so that
            // we only pass to HapCUT2 heterozygous variants
            // in this case, we set them all to 1 so we generate fragments for all variants
            let phase_variant: Vec<bool> = vec![true; varlist.lst.len()];
            // generate_flist_buffer generates a Vec<Vec<u8>> where each inner vector is a file line
            // together the lines represent the contents of a fragment file in HapCUT-like format
            let mut fragment_buffer =
                generate_flist_buffer(&flist, &phase_variant, max_p_miscall, true)
                    .chain_err(|| "Error generating fragment list buffer.")?;

            // convert the buffer of u8s into strings and print them to the fragment file
            let fragment_file_path = Path::new(&ffn);
            let mut fragment_file = File::create(&fragment_file_path)
                .chain_err(|| "Could not open fragment file for writing.")?;
            for mut line_u8 in fragment_buffer {
                line_u8.pop();
                writeln!(fragment_file, "{}", u8_to_string(&line_u8)?)
                    .chain_err(|| "Error writing to fragment file.")?;
            }
        }
        &None => {}
    }

    /***********************************************************************************************/
    // CALL GENOTYPES USING REFINED QUALITY SCORES
    /***********************************************************************************************/

    eprintln!(
        "{} Calling initial genotypes using pair-HMM realignment...",
        print_time()
    );
    call_genotypes_no_haplotypes(&flist, &mut varlist, &genotype_priors, max_p_miscall)
        .chain_err(|| "Error calling initial genotypes with estimated allele qualities.")?;

    // use Fishers exact test to check if allele observations are biased toward one strand or the other
    for mut var in &mut varlist.lst {
        if !var.alleles.len() == 2 {
            continue;
        }
        let counts: [u32; 4] = [
            var.allele_counts_forward[0] as u32,
            var.allele_counts_reverse[0] as u32,
            var.allele_counts_forward[1] as u32,
            var.allele_counts_reverse[1] as u32,
        ];
        let fishers_exact_pvalues = fishers_exact(&counts)
            .chain_err(|| "Error calculating Fisher's exact test for strand bias.")?;;

        //println!("{:?} {:?} {:?}  {:?}",&counts, fishers_exact_pvalues.two_tail_pvalue, fishers_exact_pvalues.less_pvalue, fishers_exact_pvalues.greater_pvalue);
        var.strand_bias_pvalue = if fishers_exact_pvalues.two_tail_pvalue <= 500.0 {
            *PHREDProb::from(Prob(fishers_exact_pvalues.two_tail_pvalue))
        } else {
            500.0
        };

        if fishers_exact_pvalues.two_tail_pvalue < strand_bias_pvalue_cutoff {
            var.filter.add_filter(VarFilter::StrandBias);
            var.genotype = Genotype(0, 0);
            var.gq = 0.0;
        }
    }

    for f in 0..flist.len() {
        &flist[f].calls.retain(|&c| {
            !varlist.lst[c.var_ix as usize]
                .filter
                .has_filter(VarFilter::StrandBias)
        });
    }

    print_variant_debug(
        &mut varlist,
        &interval,
        &variant_debug_directory,
        &"2.0.realigned_genotypes.vcf",
        max_cov,
        &density_params,
        &sample_name,
    )?;

    // if haplotype information usage is turned off, immediately print VCF and terminate.
    if no_haps {
        print_vcf(
            &mut varlist,
            &interval,
            &None,
            &output_vcf_file,
            false,
            max_cov,
            &density_params,
            &sample_name,
            false,
            min_gq
        )
        .chain_err(|| "Error printing VCF output.")?;
        return Ok(());
    }
    /***********************************************************************************************/
    // ITERATIVELY ASSEMBLE HAPLOTYPES AND CALL GENOTYPES
    /***********************************************************************************************/

    eprintln!(
        "{} Iteratively assembling haplotypes and refining genotypes...",
        print_time()
    );
    call_genotypes_with_haplotypes(
        &mut flist,
        &mut varlist,
        &interval,
        &genotype_priors,
        &variant_debug_directory,
        3,
        max_cov,
        &density_params,
        max_p_miscall,
        &sample_name,
        ll_delta,
    )
    .chain_err(|| "Error during haplotype/genotype iteration procedure.")?;

    // if haplotype-based read separation is turned on,
    // write BAM files for h1,h2, and unassigned
    match hap_bam_prefix {
        Some(p) => {
            eprintln!(
                "{} Calculating fraction of reads assigned to either haplotype...",
                print_time()
            );
            // h1 and h2 are hash-sets containing the qnames of the reads assigned to haplotype 1 and 2 respectively.
            let (h1, h2) = separate_fragments_by_haplotype(
                &flist,
                LogProb::from(Prob(1.0 - hap_max_p_misassign)),
            )?;

            eprintln!(
                "{} Writing haplotype-assigned reads to bam files...",
                print_time()
            );
            separate_bam_reads_by_haplotype(
                &bamfile_name,
                &interval,
                p.to_string(),
                &h1,
                &h2,
                min_mapq,
            )
            .chain_err(|| "Error separating BAM reads by haplotype.")?;
        }
        None => {}
    }

    if use_poa {
        /***********************************************************************************************/
        // PERFORM PARTIAL ORDER ALIGNMENT TO FIND NEW VARIANTS
        /***********************************************************************************************/
        let (h1, h2) = separate_fragments_by_haplotype(
            &flist,
            LogProb::from(Prob(0.51)),
        )?;

        eprintln!(
            "{} Using Partial Order Alignment (POA) to find new variants...",
            print_time()
        );

        let mut varlist_poa = call_potential_snvs::call_potential_variants_poa(
            &bamfile_name,
            &fasta_file,
            &interval,
            &h1,
            &h2,
            min_mapq,
            alignment_parameters.ln(),
        )
        .chain_err(|| "Error calling new potential variants using POA.")?;

        // want to remove duplicates from POA variant list
        let mut varlist_empty: VarList = VarList::new(vec![], varlist.target_names.clone())?;
        varlist_poa.combine(&mut varlist_empty)?;

        print_variant_debug(
            &mut varlist_poa,
            &interval,
            &variant_debug_directory,
            &"4.0.potential_SNVs_identified_with_POA.vcf",
            max_cov,
            &density_params,
            &sample_name,
        )?;

        eprintln!("{} Merging POA variants with pileup SNVs...", print_time());

        varlist.combine(&mut varlist_poa)?;

        print_variant_debug(
            &mut varlist,
            &interval,
            &variant_debug_directory,
            &"5.0.merged_potential_SNVs_pileup_and_POA.vcf",
            max_cov,
            &density_params,
            &sample_name,
        )?;

        eprintln!(
            "{} {} potential variants after POA.",
            print_time(),
            varlist.lst.len()
        );

        /***********************************************************************************************/
        // PRODUCE FRAGMENT DATA FOR NEW VARIANTS
        /***********************************************************************************************/
        eprintln!(
            "{} Producing condensed read data for POA variants...",
            print_time()
        );
        let mut flist2 = extract_fragments::extract_fragments(
            &bamfile_name,
            &fasta_file,
            &mut varlist,
            &interval,
            extract_fragment_parameters,
            alignment_parameters,
        )
        .chain_err(|| "Error extracting fragments for new POA variants.")?;

        call_genotypes_no_haplotypes(&flist2, &mut varlist, &genotype_priors, max_p_miscall)
            .chain_err(|| "Error calling initial genotypes for new POA variants.")?;

        print_variant_debug(
            &mut varlist,
            &interval,
            &variant_debug_directory,
            &"6.0.realigned_genotypes_after_POA.vcf",
            max_cov,
            &density_params,
            &sample_name,
        )?;

        eprintln!(
            "{} Iteratively assembling haplotypes and refining genotypes (with POA variants)...",
            print_time()
        );
        call_genotypes_with_haplotypes(
            &mut flist2,
            &mut varlist,
            &interval,
            &genotype_priors,
            &variant_debug_directory,
            7,
            max_cov,
            &density_params,
            max_p_miscall,
            &sample_name,
            ll_delta,
        )
        .chain_err(|| {
            "Error during haplotype/genotype iteration procedure with new POA variants."
        })?;;

        /***********************************************************************************************/
        // PERFORM FINAL FILTERING STEPS AND PRINT OUTPUT VCF
        /***********************************************************************************************/
        //calculate_mec(&flist2, &mut varlist);
    }

    // calculate MEC-based statistics for variants and blocks
    calculate_mec(&flist, &mut varlist, max_p_miscall)
        .chain_err(|| "Error calculating MEC for haplotype blocks.")?;

    let debug_filename = if use_poa {
        "8.0.final_genotypes.vcf"
    } else {
        "4.0.final_genotypes.vcf"
    };

    // Print the final VCF output
    eprintln!("{} Printing VCF file...", print_time());
    print_variant_debug(
        &mut varlist,
        &interval,
        &variant_debug_directory,
        debug_filename,
        max_cov,
        &density_params,
        &sample_name,
    )?;
    print_vcf(
        &mut varlist,
        &interval,
        &Some(fasta_file),
        &output_vcf_file,
        false,
        max_cov,
        &density_params,
        &sample_name,
        false,
        min_gq
    )
    .chain_err(|| "Error printing VCF output.")?;

    Ok(())
}
