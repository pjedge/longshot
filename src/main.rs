#![allow(dead_code)]

extern crate bio;
extern crate clap;
extern crate rust_htslib;
#[macro_use]
extern crate quick_error;
extern crate core;
extern crate chrono;
extern crate rand;

mod haplotype_assembly;
mod call_potential_snvs;
mod extract_fragments;
mod call_genotypes;
mod realignment;
mod util;
//mod de_bruijn;
use clap::{Arg, App};
use chrono::prelude::*;

use call_genotypes::{call_genotypes, call_realigned_genotypes_no_haplotypes, print_vcf, var_filter};
use util::{GenomicInterval, ExtractFragmentParameters, AlignmentParameters, parse_region_string, AlignmentType};

static PACBIO_ALIGNMENT_PARAMETERS: AlignmentParameters = AlignmentParameters {
    // non-homopolymer probabilities
    match_from_match: 0.95,
    mismatch_from_match: 0.01,
    insertion_from_match: 0.02,
    deletion_from_match: 0.02,
    extend_from_insertion: 0.13,
    match_from_insertion: 0.8609375,
    mismatch_from_insertion: 0.0090625,
    extend_from_deletion: 0.06,
    match_from_deletion: 0.9290697,
    mismatch_from_deletion: 0.010930,
    // homopolymer probabilities
    match_from_match_homopolymer: 0.83,
    mismatch_from_match_homopolymer: 0.01,
    insertion_from_match_homopolymer: 0.12,
    deletion_from_match_homopolymer: 0.04,
    extend_from_insertion_homopolymer: 0.5,
    match_from_insertion_homopolymer: 0.49397590,
    mismatch_from_insertion_homopolymer: 0.00602409,
    extend_from_deletion_homopolymer: 0.25,
    match_from_deletion_homopolymer: 0.7409638,
    mismatch_from_deletion_homopolymer: 0.00903614,
};

static ONT_ALIGNMENT_PARAMETERS: AlignmentParameters = AlignmentParameters {
    // non-homopolymer probabilities
    match_from_match: 0.82,
    mismatch_from_match: 0.05,
    insertion_from_match: 0.05,
    deletion_from_match: 0.08,
    extend_from_insertion: 0.25,
    match_from_insertion: 0.7069,
    mismatch_from_insertion: 0.0431,
    extend_from_deletion: 0.35,
    match_from_deletion: 0.6126,
    mismatch_from_deletion: 0.0373,
    // homopolymer probabilities
    match_from_match_homopolymer: 0.82,
    mismatch_from_match_homopolymer: 0.05,
    insertion_from_match_homopolymer: 0.05,
    deletion_from_match_homopolymer: 0.08,
    extend_from_insertion_homopolymer: 0.25,
    match_from_insertion_homopolymer: 0.7069,
    mismatch_from_insertion_homopolymer: 0.0431,
    extend_from_deletion_homopolymer: 0.35,
    match_from_deletion_homopolymer: 0.6126,
    mismatch_from_deletion_homopolymer: 0.0373,
};

fn main() {

    let input_args = App::new("Reaper (REAlign error PronE Reads)")
        .version("0.1")
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
                .help("indexed fasta reference that BAM file is aligned to")
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
                .value_name("VCF")
                .help("Region in format <chrom> or <chrom:start-stop> in which to call variants.")
                .display_order(40)
                .takes_value(true))
        .arg(Arg::with_name("Max coverage")
                .short("C")
                .long("max_cov")
                .value_name("int")
                .help("Maximum coverage to consider position as a potential SNV.")
                .display_order(80))
        .arg(Arg::with_name("Min mapq")
                .short("q")
                .long("min_mapq")
                .value_name("int")
                .help("Minimum mapping quality to use a read.")
                .display_order(90)
                .default_value("30"))
        .arg(Arg::with_name("Anchor length")
                .short("l")
                .long("anchor_length")
                .value_name("int")
                .help("Length of indel-free anchor sequence on the left and right side of read realignment window.")
                .display_order(100)
                .default_value("10"))
        .arg(Arg::with_name("Anchor k")
                .short("k")
                .long("anchor_k")
                .value_name("int")
                .help("A filter for low-complexity anchor sequences. A valid anchor must have no duplicates in the kmers that overlap it.")
                .display_order(110)
                .default_value("5"))
        .arg(Arg::with_name("Short haplotype max SNVs")
                .short("m")
                .long("max_snvs")
                .value_name("int")
                .help("Cut off short haplotypes after this many SNVs. 2^m haplotypes must be aligned against per read for a variant cluster of size m.")
                .display_order(130)
                .default_value("5"))
        .arg(Arg::with_name("Max window length")
                .short("W")
                .long("max_window")
                .value_name("int")
                .help("Ignore a variant/short haplotype if the realignment window for read or reference is larger than w bases.")
                .display_order(150)
                .default_value("200"))
        .arg(Arg::with_name("Fast alignment")
                .short("z")
                .long("fast_alignment")
                .help("Use non-numerically stable alignment algorithm. Is significantly faster but may be less accurate or have unexpected behaviour.")
                .display_order(160))
        .arg(Arg::with_name("Band width")
                .short("B")
                .long("band_width")
                .help("Minimum width of alignment band. Band will increase in size if sequences are different lengths.")
                .display_order(170)
                .default_value("20"))
        .arg(Arg::with_name("Read technology")
                .short("t")
                .long("read_technology")
                .help("Which read technology is being used (\"ont\" or \"pacbio\").")
                .default_value("pacbio")
                .display_order(180))
        .arg(Arg::with_name("No haplotypes")
                .short("n")
                .long("no_haps")
                .help("Don't call HapCUT2 to phase variants.")
                .display_order(190))
        .arg(Arg::with_name("Indels")
            .short("i")
            .long("indels")
            .help("Report short indels -- called indels are wildly inaccurate but are used internally to avoid false SNVs.")
            .display_order(200))
        .arg(Arg::with_name("Max alignment")
            .short("x")
            .long("max_alignment")
            .help("Use max alignment algorithm rather than all-alignments algorithm.")
            .display_order(165))
        .get_matches();

    // should be safe just to unwrap these because they're required options for clap
    let bamfile_name = input_args.value_of("Input BAM").unwrap().to_string();
    let fasta_file = input_args.value_of("Input FASTA").unwrap().to_string();
    let output_vcf_file = input_args.value_of("Output VCF").unwrap().to_string();

    let interval: Option<GenomicInterval> = parse_region_string(input_args.value_of("Region"),
                                                                &bamfile_name);
    // TODO check that min_alt_frac is between 0 and 1

    let max_cov: Option<u32> = match input_args.value_of("Max coverage") {
        Some(cov) => {
            Some(cov.parse::<u32>().expect("Argument max_cov must be a positive integer!"))
        }
        None => None,
    };

    let min_mapq: u8 = input_args.value_of("Min mapq")
        .unwrap()
        .parse::<u8>()
        .expect("Argument min_mapq must be an integer between 0 and 255!");

    let anchor_length: usize = input_args.value_of("Anchor length")
        .unwrap()
        .parse::<usize>()
        .expect("Argument anchormust be a positive integer!");

    let anchor_k: usize = input_args.value_of("Anchor k")
        .unwrap()
        .parse::<usize>()
        .expect("Argument anchor_k must be a positive integer!");
    // TODO check that anchor_k is less than or equal to anchor length

    let short_hap_max_snvs: usize = input_args.value_of("Short haplotype max SNVs")
        .unwrap()
        .parse::<usize>()
        .expect("Argument max_snvs must be a positive integer!");

    let max_window_length: usize = input_args.value_of("Max window length")
        .unwrap()
        .parse::<usize>()
        .expect("Argument max_window must be a positive integer!");

    let mut alignment_type = AlignmentType::NumericallyStableAllAlignment;

    match input_args.occurrences_of("Fast alignment") {
        0 => {},
        1 => {alignment_type = AlignmentType::FastAllAlignment;},
        _ => {
            panic!("fast_alignment specified multiple times");
        }
    };

    match input_args.occurrences_of("Max alignment") {
        0 => {},
        1 => {
            if alignment_type == AlignmentType::FastAllAlignment {
                panic!("Fast alignment and max alignment options are incompatible. The max alignment uses sums rather than pows/logs so it is already fast compared to default.");
            }
            alignment_type = AlignmentType::MaxAlignment;
        },
        _ => {
            panic!("max_alignment specified multiple times");
        }
    };

    let mut band_width: usize = input_args.value_of("Band width")
        .unwrap()
        .parse::<usize>()
        .expect("Argument band_width must be a positive integer!");

    let alignment_parameters;
    match input_args.value_of("Read technology").unwrap() {
        "pacbio" | "Pacbio" | "PacBio" | "PACBIO" => {
            alignment_parameters = PACBIO_ALIGNMENT_PARAMETERS.clone();
        }
        "ont" | "ONT" => {
            alignment_parameters = ONT_ALIGNMENT_PARAMETERS.clone();
            band_width = 50
        }
        _ => {
            panic!("Invalid read technology argument.");
        }
    }


    let assemble_haps = match input_args.occurrences_of("No haplotypes") {
        0 => true,
        1 => false,
        _ => {
            panic!("no_haps specified multiple times");
        }
    };


    let indels = match input_args.occurrences_of("Indels") {
        0 => false,
        1 => true,
        _ => {
            panic!("Indels specified multiple times");
        }
    };

    let print_time: fn() -> String = || Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

    //let bam_file: String = "test_data/test.bam".to_string();
    eprintln!("{} Calling potential SNVs using pileup...",print_time());
    let mut varlist = call_potential_snvs::call_potential_snvs(&bamfile_name,
                                                           &fasta_file,
                                                           &interval,
                                                               max_cov,
                                                               min_mapq,
                                                           alignment_parameters.ln());

    eprintln!("{} {} potential SNVs identified.", print_time(),varlist.lst.len());

    let extract_fragment_parameters = ExtractFragmentParameters {
        min_mapq: min_mapq,
        alignment_type: alignment_type,
        band_width: band_width,
        anchor_length: anchor_length,
        anchor_k: anchor_k,
        short_hap_max_snvs: short_hap_max_snvs,
        max_window_length: max_window_length,
    };

    eprintln!("{} Generating condensed read data for SNVs...",print_time());
    let flist = extract_fragments::extract_fragments(&bamfile_name,
                                                         &fasta_file,
                                                         &varlist,
                                                         &interval,
                                                         extract_fragment_parameters,
                                                         alignment_parameters);

    eprintln!("{} Calling initial genotypes and assembling haplotypes...", print_time());
    match assemble_haps {
        true => {
            call_realigned_genotypes_no_haplotypes(&flist, &mut varlist);
        },
        false => {panic!("Calling genotypes without haplotypes not currently supported.")},
    };

    eprintln!("{} Refining genotypes with haplotype information...",print_time());
    call_genotypes(&flist, &mut varlist, &interval);

    eprintln!("{} Adding filter flags based on depth and variant density...",print_time());
    var_filter(&mut varlist, 50.0, 500, 10, max_cov);

    eprintln!("{} Printing VCF file...",print_time());
    print_vcf(&varlist, &interval, indels, output_vcf_file);
}
