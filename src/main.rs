#![allow(dead_code)]

extern crate bio;
extern crate clap;
extern crate rust_htslib;
#[macro_use]
extern crate quick_error;
extern crate core;
extern crate chrono;
extern crate rand;
extern crate petgraph;

mod haplotype_assembly;
mod call_potential_snvs;
mod extract_fragments;
mod call_genotypes;
mod realignment;
mod util;
mod estimate_read_coverage;
mod estimate_alignment_parameters;
mod spoa;
mod variants_and_fragments;
mod print_output;
mod genotype_probs;

use clap::{Arg, App};
use std::fs::create_dir;
use std::path::Path;
use std::fs::remove_dir_all;

use call_genotypes::*;
use util::*;
use estimate_read_coverage::calculate_mean_coverage;
use estimate_alignment_parameters::estimate_alignment_parameters;
use bio::stats::{LogProb,Prob, PHREDProb};
use haplotype_assembly::separate_reads_by_haplotype;
use print_output::{print_variant_debug, print_vcf};
use realignment::{AlignmentType};
//use realignment::{AlignmentParameters, TransitionProbs, EmissionProbs};
use genotype_probs::GenotypePriors;
use extract_fragments::ExtractFragmentParameters;

fn main() {

    /***********************************************************************************************/
    // READ STANDARD INPUT
    /***********************************************************************************************/

    eprintln!("");

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
                .value_name("string")
                .help("Region in format <chrom> or <chrom:start-stop> in which to call variants.")
                .display_order(40)
                .takes_value(true))
        .arg(Arg::with_name("Max coverage fraction")
            .short("c")
            .long("max_cov_frac")
            .value_name("float")
            .help("Maximum coverage (of reads passing filters) to consider position as a potential SNV, as a fraction of the mean coverage. For example, \"-c 2.0\" throws out positions with more than twice the mean coverage.")
            .display_order(79)
            .default_value("2.0"))
        .arg(Arg::with_name("Max coverage")
                .short("C")
                .long("max_cov")
                .value_name("int")
                .help("Maximum coverage (of reads passing filters) to consider position as a potential SNV. Overrides --max_cov_frac parameter.")
                .display_order(80))
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
            .display_order(90)
            .default_value("7.0"))
        .arg(Arg::with_name("Min genotype quality for haplotype scaffold")
            .short("Q")
            .long("min_hap_gq")
            .value_name("float")
            .help("Minimum genotype quality (Phred-scaled) of a variant to use its haplotype phase to genotype other variants.")
            .display_order(90)
            .default_value("50.0"))
        .arg(Arg::with_name("Anchor length")
                .short("l")
                .long("anchor_length")
                .value_name("int")
                .help("Length of indel-free anchor sequence on the left and right side of read realignment window.")
                .display_order(100)
                .default_value("6"))
        .arg(Arg::with_name("Short haplotype max SNVs")
                .short("m")
                .long("max_snvs")
                .value_name("int")
                .help("Cut off short haplotypes after this many SNVs. 2^m haplotypes must be aligned against per read for a variant cluster of size m.")
                .display_order(130)
                .default_value("3"))
        .arg(Arg::with_name("Use POA")
            .short("p")
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
        .arg(Arg::with_name("Fast alignment")
                .short("z")
                .long("fast_alignment")
                .help("Use non-numerically stable pair HMM algorithm. Is significantly faster but may be less accurate or have unexpected behaviour.")
                .display_order(160))
        .arg(Arg::with_name("Force overwrite")
            .short("F")
            .long("force_overwrite")
            .help("If output files (VCF or variant debug directory) exist, delete and overwrite them.")
            .display_order(164))
        .arg(Arg::with_name("Band width")
                .short("B")
                .long("band_width")
                .help("Minimum width of alignment band. Band will increase in size if sequences are different lengths.")
                .display_order(170)
                .default_value("20"))
        .arg(Arg::with_name("Max iters since likelihood improvement")
            .short("I")
            .long("max_iters")
            .help("Maximum rounds of haplotype assembly without likelihood improvement before termination.")
            .display_order(170)
            .default_value("5"))
        .arg(Arg::with_name("Sample ID")
            .short("s")
            .long("sample_id")
            .value_name("string")
            .help("Specify a sample ID to write to the output VCF")
            .display_order(177)
            .default_value(&"SAMPLE"))
        /*.arg(Arg::with_name("No haplotypes")
                .short("n")
                .long("no_haps")
                .help("Don't call HapCUT2 to phase variants.")
                .display_order(190))*/
        .arg(Arg::with_name("Max alignment")
            .short("x")
            .long("max_alignment")
            .help("Use max scoring alignment algorithm rather than pair HMM forward algorithm.")
            .display_order(165))
        .arg(Arg::with_name("Variant debug directory")
            .short("d")
            .long("variant_debug_dir")
            .value_name("path")
            .help("write out current information about variants at each step of algorithm to files in this directory")
            .display_order(210))
        .get_matches();

    // should be safe just to unwrap these because they're required options for clap
    let bamfile_name = input_args.value_of("Input BAM").unwrap().to_string();
    let fasta_file = input_args.value_of("Input FASTA").unwrap().to_string();
    let output_vcf_file = input_args.value_of("Output VCF").unwrap().to_string();

    let interval: Option<GenomicInterval> = parse_region_string(input_args.value_of("Region"),
                                                                &bamfile_name);

    let force = match input_args.occurrences_of("Force overwrite") {
        0 => false,
        1 => true,
        _ => {
            panic!("force_overwrite specified multiple times");
        }
    };

    let vcf = Path::new(&output_vcf_file);
    if vcf.is_file() {
        if !force {
            eprintln!("{} ERROR: Variant output file already exists. Rerun with -F option to force overwrite.", print_time());
            return;
        }
    }

    let bai_str = bamfile_name.clone() + ".bai";
    if !Path::new(&bai_str).is_file() {
        eprintln!("{} ERROR: BAM file must be indexed with samtools index. Index file should have same name with .bai appended.",print_time());
        return;
    }

    let sample_name: String = input_args.value_of(&"Sample ID").unwrap().to_string();

    let variant_debug_directory: Option<String> = match input_args.value_of("Variant debug directory") {
        Some(dir) => {
            let p = Path::new(&dir);
            if p.exists() {
                if force {
                    match remove_dir_all(p) {
                        Ok(_) => {},
                        Err(_) => {panic!("{} Error removing variant debug directory.", print_time())}
                    }
                } else {
                    eprintln!("{} ERROR: Variant debug directory already exists. Rerun with -F option to force overwrite.",print_time());
                    return;
                }
            }
            match create_dir(&dir) {
                Ok(_) => {},
                Err(_) => {panic!("{} Error creating variant debug directory.",print_time());}
            }
            Some(dir.to_string())
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
        .expect("Argument anchor must be a positive integer!");

    let short_hap_max_snvs: usize = input_args.value_of("Short haplotype max SNVs")
        .unwrap()
        .parse::<usize>()
        .expect("Argument max_snvs must be a positive integer!");

    let max_window_padding: usize = input_args.value_of("Max window padding")
        .unwrap()
        .parse::<usize>()
        .expect("Argument max_window must be a positive integer!");

    let min_allele_qual: f64 = input_args.value_of("Min allele quality")
        .unwrap()
        .parse::<f64>()
        .expect("Argument max_mec_frac must be a positive float!");

    if min_allele_qual <= 0.0 {
        panic!("Min allele quality must be a positive float.");
    }

    let min_hap_gq: f64 = input_args.value_of("Min genotype quality for haplotype scaffold")
        .unwrap()
        .parse::<f64>()
        .expect("Argument min_hap_gq must be a positive float!");

    if min_hap_gq <= 0.0 {
        panic!("min_hap_gq must be a positive float.");
    }

    let max_p_miscall: f64 = *Prob::from(PHREDProb(min_allele_qual));

    /*
    let max_mec_frac = match input_args.occurrences_of("Max MEC Fraction") {
        0 => None,
        1 => {
            let frac:f64 = input_args.value_of("Max MEC Fraction")
                .unwrap()
                .parse::<f64>()
                .expect("Argument max_mec_frac must be a positive float!");

            if !(frac >= 0.0 && frac <= 1.0) {
                panic!("Max MEC Fraction must be a float between 0.0 and 1.0.");
            }

            Some(frac)
        },
        _ => {
            panic!("max_mec_frac specified multiple times");
        }
    };*/

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

    let band_width: usize = input_args.value_of("Band width")
        .unwrap()
        .parse::<usize>()
        .expect("Argument band_width must be a positive integer!");

    let max_iters_since_improvement: usize = input_args.value_of("Max iters since likelihood improvement")
        .unwrap()
        .parse::<usize>()
        .expect("Argument max_iters must be a positive integer!");

    let use_poa = match input_args.occurrences_of("Use POA") {
        0 => false,
        1 => true,
        _ => {
            panic!("Use POA specified multiple times");
        }
    };

    if input_args.occurrences_of("Max coverage") > 0 && input_args.occurrences_of("Max coverage fraction") > 0 {
        eprintln!("{} ERROR: -c and -C options cannot be used at the same time.",print_time());
        return;
    }

    let max_cov: Option<u32> = match input_args.value_of("Max coverage") {
        Some(cov) => {
            Some(cov.parse::<u32>().expect("Argument max_cov must be a positive integer!"))
        }
        None => match input_args.value_of("Max coverage fraction") {
            Some(f) => {
                let frac = f.parse::<f64>().expect("Argument max_cov_frac must be a floating point value!");
                if frac <= 1.0 {
                    eprintln!("{} ERROR: Argument max_cov_frac should not be less than 1.0. This would filter out the majority of positions!", print_time());
                    return;
                }
                eprintln!("{} Estimating mean read coverage...",print_time());
                let mean_coverage: f64 = calculate_mean_coverage(&bamfile_name, &interval, min_mapq);
                let calculated_cov = (mean_coverage as f64 * frac) as u32;

                eprintln!("{} Mean read coverage: {:.2}",print_time(), mean_coverage);
                if calculated_cov > 0 {
                    eprintln!("{} Max read coverage set to {}.",print_time(), calculated_cov);
                } else {
                    eprintln!("{} ERROR: Max read coverage set to 0.",print_time());
                    return;
                }
                Some(calculated_cov)
            }
            None => None,
        },
    };

    let extract_fragment_parameters = ExtractFragmentParameters {
        min_mapq: min_mapq,
        alignment_type: alignment_type,
        band_width: band_width,
        anchor_length: anchor_length,
        short_hap_max_snvs: short_hap_max_snvs,
        max_window_padding: max_window_padding
    };

    eprintln!("{} Estimating alignment parameters...",print_time());
    let alignment_parameters = estimate_alignment_parameters(&bamfile_name, &fasta_file, &interval);
    /*
    let alignment_parameters = AlignmentParameters {
        transition_probs: TransitionProbs {
                            match_from_match: 0.9,
                            insertion_from_match: 0.08,
                            deletion_from_match: 0.02,
                            insertion_from_insertion: 0.25,
                            match_from_insertion: 0.75,
                            deletion_from_deletion: 0.1,
                            match_from_deletion: 0.9,
        },
        emission_probs: EmissionProbs {
                            equal: 0.99,
                            not_equal: 0.01 / 3.0,
                            deletion: 1.0,
                            insertion: 1.0
        }
    };*/


    /***********************************************************************************************/
    // GET HUMAN GENOTYPE PRIORS
    /***********************************************************************************************/

    let genotype_priors = GenotypePriors::new();

    /***********************************************************************************************/
    // FIND INITIAL SNVS WITH READ PILEUP
    /***********************************************************************************************/

    //let bam_file: String = "test_data/test.bam".to_string();
    eprintln!("{} Calling potential SNVs using pileup...",print_time());
    let mut varlist = call_potential_snvs::call_potential_snvs(&bamfile_name,
                                                           &fasta_file,
                                                           &interval,
                                                               &genotype_priors,
                                                               max_cov,
                                                               min_mapq,
                                                               max_p_miscall,
                                                           alignment_parameters.ln());

    // back up the variant indices
    // they will be needed later when we try to re-use fragment alleles that don't change
    // as the variant list expands
    varlist.backup_indices();

    print_variant_debug(&mut varlist, &interval, &variant_debug_directory,&"1.0.potential_SNVs.vcf", max_cov, &sample_name);

    eprintln!("{} {} potential SNVs identified.", print_time(),varlist.lst.len());

    if varlist.lst.len() == 0 {
        return;
    }
    /***********************************************************************************************/
    // EXTRACT FRAGMENT INFORMATION FROM READS
    /***********************************************************************************************/

    eprintln!("{} Generating condensed read data for SNVs...",print_time());
    let mut flist = extract_fragments::extract_fragments(&bamfile_name,
                                                         &fasta_file,
                                                         &varlist,
                                                         &interval,
                                                         extract_fragment_parameters,
                                                         alignment_parameters,
                                                          None);

    /***********************************************************************************************/
    // CALL GENOTYPES USING REFINED QUALITY SCORES
    /***********************************************************************************************/

    eprintln!("{} Calling initial genotypes using pair-HMM realignment...", print_time());
    call_genotypes_no_haplotypes(&flist, &mut varlist, &genotype_priors, max_p_miscall);

    print_variant_debug(&mut varlist, &interval, &variant_debug_directory,&"2.0.realigned_genotypes.vcf", max_cov, &sample_name);


    /***********************************************************************************************/
    // ITERATIVELY ASSEMBLE HAPLOTYPES AND CALL GENOTYPES
    /***********************************************************************************************/

    eprintln!("{} Iteratively assembling haplotypes and refining genotypes...",print_time());
    call_genotypes_with_haplotypes(&mut flist, &mut varlist, &interval, &genotype_priors,
                                   &variant_debug_directory, 3, max_cov, max_p_miscall,
                                   min_hap_gq, max_iters_since_improvement, &sample_name);


    if use_poa {
        /***********************************************************************************************/
        // PERFORM PARTIAL ORDER ALIGNMENT TO FIND NEW VARIANTS
        /***********************************************************************************************/

        let (h1,h2) = separate_reads_by_haplotype(&flist, LogProb::from(Prob(0.99)));

        eprintln!("{} Using Partial Order Alignment (POA) to find new variants...", print_time());

        let mut varlist_poa = call_potential_snvs::call_potential_variants_poa(&bamfile_name,
                                                                   &fasta_file,
                                                                   &interval,
                                                                   &h1,
                                                                   &h2,
                                                                   max_cov,
                                                                   min_mapq,
                                                                   alignment_parameters.ln());

        eprintln!("{} Merging POA variants with pileup SNVs...",print_time());


        varlist.combine(&mut varlist_poa);

        print_variant_debug(&mut varlist, &interval, &variant_debug_directory,&"4.0.new_potential_SNVs_after_POA.vcf", max_cov, &sample_name);

        eprintln!("{} {} potential variants after POA.", print_time(),varlist.lst.len());

        /***********************************************************************************************/
        // PRODUCE FRAGMENT DATA FOR NEW VARIANTS
        /***********************************************************************************************/

        eprintln!("{} Producing condensed read data for POA variants...",print_time());
        let mut flist2 = extract_fragments::extract_fragments(&bamfile_name,
                                                             &fasta_file,
                                                             &varlist,
                                                             &interval,
                                                             extract_fragment_parameters,
                                                             alignment_parameters,
                                                              None);  // Some(flist)


        call_genotypes_no_haplotypes(&flist2, &mut varlist, &genotype_priors, max_p_miscall); // temporary
        print_variant_debug(&mut varlist, &interval, &variant_debug_directory,&"5.0.realigned_genotypes_after_POA.vcf", max_cov, &sample_name);

        eprintln!("{} Iteratively assembling haplotypes and refining genotypes (with POA variants)...",print_time());
        call_genotypes_with_haplotypes(&mut flist2, &mut varlist, &interval, &genotype_priors,
            &variant_debug_directory, 6, max_cov, max_p_miscall, min_hap_gq, max_iters_since_improvement,
            &sample_name);

        /***********************************************************************************************/
        // PERFORM FINAL FILTERING STEPS AND PRINT OUTPUT VCF
        /***********************************************************************************************/

        //calculate_mec(&flist2, &mut varlist);
    }

    let debug_filename = if use_poa {
        "7.0.final_genotypes.vcf"
    } else {
        "4.0.final_genotypes.vcf"
    };

    eprintln!("{} Printing VCF file...",print_time());
    print_variant_debug(&mut varlist, &interval,&variant_debug_directory, &debug_filename, max_cov, &sample_name);
    print_vcf(&mut varlist, &interval, &output_vcf_file, false, max_cov, &sample_name);
}
