#![allow(dead_code)]

extern crate bio;
extern crate clap;
extern crate rust_htslib;
#[macro_use]
extern crate quick_error;
extern crate core;

mod haplotype_assembly;
mod call_potential_snvs;
mod extract_fragments;
mod call_genotypes;
mod realignment;
mod util;
//mod de_bruijn;
use clap::{Arg, App};

use call_genotypes::{call_genotypes, call_haplotypes};
use util::{GenomicInterval, ExtractFragmentParameters, AlignmentParameters, parse_region_string};

static PACBIO_ALIGNMENT_PARAMETERS: AlignmentParameters = AlignmentParameters {
    match_from_match: 0.89,
    mismatch_from_match: 0.01,
    insertion_from_match: 0.07,
    deletion_from_match: 0.03,
    extend_from_insertion: 0.26,
    match_from_insertion: 0.7317777,
    mismatch_from_insertion: 0.0082222,
    extend_from_deletion: 0.12,
    match_from_deletion: 0.8702222,
    mismatch_from_deletion: 0.0097777,
};

static ONT_ALIGNMENT_PARAMETERS: AlignmentParameters = AlignmentParameters {
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
                .display_order(1)
                .required(true)
                .takes_value(true))
        .arg(Arg::with_name("Input FASTA")
                .short("f")
                .long("ref")
                .value_name("FASTA")
                .help("indexed fasta reference that BAM file is aligned to")
                .display_order(2)
                .required(true)
                .takes_value(true))
        .arg(Arg::with_name("Output VCF")
                .short("o")
                .long("out")
                .value_name("VCF")
                .help("output VCF file with called variants.")
                .display_order(3)
                .required(true)
                .takes_value(true))
        .arg(Arg::with_name("Region")
                .short("r")
                .long("region")
                .value_name("VCF")
                .help("Region in format <chrom> or <chrom:start-stop> in which to call variants.")
                .display_order(4)
                .takes_value(true))
        .arg(Arg::with_name("Min alt count")
                .short("a")
                .long("min_alt_count")
                .value_name("int")
                .help("Minimum number of occurrences of a variant base in pileup to consider it as a potential SNV.")
                .display_order(5)
                .default_value("2"))
        .arg(Arg::with_name("Min alt frac")
                .short("A")
                .long("min_alt_frac")
                .value_name("float")
                .help("Minimum fraction of variant base in pileup to consider it as a potential SNV.")
                .display_order(6)
                .default_value("0.125"))
        .arg(Arg::with_name("Min coverage")
                .short("c")
                .long("min_cov")
                .value_name("int")
                .help("Minimum coverage to consider position as a potential SNV.")
                .display_order(7)
                .default_value("5"))
        .arg(Arg::with_name("Max coverage")
                .short("C")
                .long("max_cov")
                .value_name("int")
                .help("Maximum coverage to consider position as a potential SNV.")
                .display_order(8))
        .arg(Arg::with_name("Min mapq")
                .short("q")
                .long("min_mapq")
                .value_name("int")
                .help("Minimum mapping quality to use a read.")
                .display_order(9)
                .default_value("30"))
        .arg(Arg::with_name("Anchor length")
                .short("l")
                .long("anchor_length")
                .value_name("int")
                .help("Length of indel-free anchor sequence on the left and right side of read realignment window.")
                .display_order(10)
                .default_value("10"))
        .arg(Arg::with_name("Anchor k")
                .short("k")
                .long("anchor_k")
                .value_name("int")
                .help("A filter for low-complexity anchor sequences. A valid anchor must have no duplicates in the kmers that overlap it.")
                .display_order(11)
                .default_value("5"))
        .arg(Arg::with_name("Short haplotype SNV distance")
                .short("d")
                .long("snv_distance")
                .value_name("int")
                .help("SNVs separated by distance less than d will be considered together and all their possible haplotypes will be aligned against.")
                .display_order(12)
                .default_value("20"))
        .arg(Arg::with_name("Short haplotype max SNVs")
                .short("m")
                .long("max_snvs")
                .value_name("int")
                .help("Cut off short haplotypes after this many SNVs. 2^m haplotypes must be aligned against per read for a variant cluster of size m.")
                .display_order(13)
                .default_value("5"))
        .arg(Arg::with_name("Min window length")
                .short("w")
                .long("min_window")
                .value_name("int")
                .help("Ignore a variant/short haplotype if the realignment window for read or reference is smaller than w bases.")
                .display_order(14)
                .default_value("15"))
        .arg(Arg::with_name("Max window length")
                .short("W")
                .long("max_window")
                .value_name("int")
                .help("Ignore a variant/short haplotype if the realignment window for read or reference is larger than w bases.")
                .display_order(15)
                .default_value("200"))
        .arg(Arg::with_name("Fast alignment")
                .short("z")
                .long("fast_alignment")
                .help("Use non-numerically stable alignment algorithm. Is significantly faster but may be less accurate or have unexpected behaviour.")
                .display_order(16))
        .arg(Arg::with_name("Band width")
                .short("B")
                .long("band_width")
                .help("Minimum width of alignment band. Band will increase in size if sequences are different lengths.")
                .display_order(17)
                .default_value("20"))
        .arg(Arg::with_name("Read technology")
                .short("t")
                .long("read_technology")
                .help("Which read technology is being used (\"ont\" or \"pacbio\").")
                .default_value("pacbio")
                .display_order(18))
        .arg(Arg::with_name("No haplotypes")
                .short("n")
                .long("no_haps")
                .help("Don't call HapCUT2 to phase variants.")
                .display_order(19))
        .get_matches();

    // should be safe just to unwrap these because they're required options for clap
    let bamfile_name = input_args.value_of("Input BAM").unwrap().to_string();
    let fasta_file = input_args.value_of("Input FASTA").unwrap().to_string();
    let output_vcf_file = input_args.value_of("Output VCF").unwrap().to_string();

    let interval: Option<GenomicInterval> = parse_region_string(input_args.value_of("Region"),
                                                                &bamfile_name);

    let min_alt_count: u32 = input_args.value_of("Min alt count")
        .unwrap()
        .parse::<u32>()
        .expect("Argument min_alt_count must be a positive integer!");

    let min_alt_frac: f32 = input_args.value_of("Min alt frac")
        .unwrap()
        .parse::<f32>()
        .expect("Argument min_alt_frac must be a float!");

    // TODO check that min_alt_frac is between 0 and 1


    let min_cov: u32 = input_args.value_of("Min coverage")
        .unwrap()
        .parse::<u32>()
        .expect("Argument min_cov must be a positive integer!");

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

    let short_hap_snv_distance: usize =
        input_args.value_of("Short haplotype SNV distance")
            .unwrap()
            .parse::<usize>()
            .expect("Argument snv_distance must be a positive integer!");

    let short_hap_max_snvs: usize = input_args.value_of("Short haplotype max SNVs")
        .unwrap()
        .parse::<usize>()
        .expect("Argument max_snvs must be a positive integer!");

    let min_window_length: usize = input_args.value_of("Min window length")
        .unwrap()
        .parse::<usize>()
        .expect("Argument min_window must be a positive integer!");

    let max_window_length: usize = input_args.value_of("Max window length")
        .unwrap()
        .parse::<usize>()
        .expect("Argument max_window must be a positive integer!");

    let numerically_stable_alignment = match input_args.occurrences_of("Fast alignment") {
        0 => true,
        1 => false,
        _ => {
            panic!("Fast_alignment specified multiple times");
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

    //let bam_file: String = "test_data/test.bam".to_string();
    eprintln!("Calling potential SNVs using pileup...");
    let varlist = call_potential_snvs::call_potential_snvs(&bamfile_name,
                                                           &fasta_file,
                                                           &interval,
                                                           min_alt_count,
                                                           min_alt_frac,
                                                           min_cov,
                                                           max_cov,
                                                           min_mapq);
    eprintln!("{} potential SNVs identified.", varlist.lst.len());

    let extract_fragment_parameters = ExtractFragmentParameters {
        numerically_stable_alignment: numerically_stable_alignment,
        band_width: band_width,
        anchor_length: anchor_length,
        anchor_k: anchor_k,
        short_hap_snv_distance: short_hap_snv_distance,
        short_hap_max_snvs: short_hap_max_snvs,
        min_window_length: min_window_length,
        max_window_length: max_window_length,
    };

    eprintln!("Generating condensed read data for SNVs...");
    let mut flist = extract_fragments::extract_fragments(&bamfile_name,
                                                         &fasta_file,
                                                         &varlist,
                                                         &interval,
                                                         extract_fragment_parameters,
                                                         alignment_parameters);

    eprintln!("Calling genotypes/haplotypes...");
    let hap: Option<Vec<char>> = match assemble_haps {
        true => Some(call_haplotypes(&flist, &varlist)),
        false => None,
    };

    match hap {
        Some(ref h) => {
            for i in 0..flist.len() {
                flist[i].assign_haps(&h);
            }
        }
        None => {}
    }

    //eprintln!("Calling genotypes...");
    call_genotypes(&flist, &varlist, &interval, &hap, output_vcf_file);
}
