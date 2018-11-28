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
mod extract_fragments; mod extract_fragments_debug;
mod call_genotypes;
mod realignment;
mod util;
mod estimate_read_coverage;
mod estimate_alignment_parameters;
//mod spoa;
mod variants_and_fragments;
mod print_output;
mod genotype_probs;

use clap::{Arg, App};
use std::fs::create_dir;
use std::fs::remove_dir_all;
use std::error::Error;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;

use call_genotypes::*;
use util::*;
use estimate_read_coverage::calculate_mean_coverage;
use estimate_alignment_parameters::estimate_alignment_parameters;
use bio::stats::{LogProb,Prob, PHREDProb};
//use haplotype_assembly::separate_reads_by_haplotype;
use print_output::{print_variant_debug, print_vcf};
use realignment::{AlignmentType};
//use realignment::{AlignmentParameters, TransitionProbs, EmissionProbs};
use genotype_probs::GenotypePriors;
use extract_fragments::ExtractFragmentParameters;
use haplotype_assembly::*;
//use variants_and_fragments::parse_VCF_potential_variants;


fn main() {

    /***********************************************************************************************/
    // READ STANDARD INPUT
    /***********************************************************************************************/

    eprintln!("");

    let input_args = App::new("Longshot")
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
            .display_order(77))
        .arg(Arg::with_name("Min coverage")
            .short("c")
            .long("min_cov")
            .value_name("int")
            .help("Minimum coverage (of reads passing filters) to consider position as a potential SNV.")
            .display_order(78)
            .default_value("0"))
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
            .display_order(90)
            .default_value("7.0"))
        .arg(Arg::with_name("Haplotype Assignment Quality")
            .short("y")
            .long("hap_assignment_qual")
            .value_name("float")
            .help("Minimum quality (Phred-scaled) of read->haplotype assignment (for read separation).")
            .display_order(90)
            .default_value("20.0"))
        .arg(Arg::with_name("Potential SNV Cutoff")
            .long("potential_snv_cutoff")
            .short("Q")
            .value_name("float")
            .help("Consider a site as a potential SNV if the original PHRED-scaled QUAL score for 0/0 genotype is below this amount (a larger value considers more potential SNV sites).")
            .display_order(91)
            .default_value("30.0"))
        .arg(Arg::with_name("Haplotype Convergence Delta")
            .long("hap_converge_delta")
            .short("L")
            .value_name("float")
            .help("Terminate the haplotype/genotype iteration when the relative change in log-likelihood falls below this amount. Setting a larger value results in faster termination but potentially less accurate results.")
            .display_order(92)
            .default_value(&"0.0001"))
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
        /*.arg(Arg::with_name("Use POA")
            .short("p")
            .long("poa")
            .help("EXPERIMENTAL: Run the algorithm twice, using Partial-Order-Alignment on phased reads to find new candidate SNVs and Indels the second time.")
            .display_order(130))*/
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
            .display_order(175)
            .default_value("10:500:50"))
        .arg(Arg::with_name("Sample ID")
            .short("s")
            .long("sample_id")
            .value_name("string")
            .help("Specify a sample ID to write to the output VCF")
            .display_order(177)
            .default_value(&"SAMPLE"))
        .arg(Arg::with_name("Homozygous SNV Rate")
            .long("hom_snv_rate")
            .value_name("float")
            .help("Specify the homozygous SNV Rate for genotype prior estimation")
            .display_order(178)
            .default_value(&"0.0005"))
        .arg(Arg::with_name("Heterozygous SNV Rate")
            .long("het_snv_rate")
            .value_name("float")
            .help("Specify the heterozygous SNV Rate for genotype prior estimation")
            .display_order(179)
            .default_value(&"0.001"))
        .arg(Arg::with_name("Homozygous Indel Rate")
            .long("hom_indel_rate")
            .value_name("float")
            .help("Specify the homozygous Indel Rate for genotype prior estimation")
            .display_order(180)
            .hidden(true)
            .default_value(&"0.0"))
        .arg(Arg::with_name("Heterozygous Indel Rate")
            .long("het_indel_rate")
            .value_name("float")
            .help("Specify the heterozygous Indel Rate for genotype prior estimation")
            .display_order(181)
            .hidden(true)
            .default_value(&"0.0"))
        .arg(Arg::with_name("ts/tv Ratio")
            .long("ts_tv_ratio")
            .value_name("float")
            .help("Specify the transition/transversion rate for genotype grior estimation")
            .display_order(182)
            .default_value(&"0.5"))
        .arg(Arg::with_name("No haplotypes")
                .short("n")
                .long("no_haps")
                .help("Don't call HapCUT2 to phase variants.")
                .display_order(190))
        .arg(Arg::with_name("Max alignment")
            .short("x")
            .long("max_alignment")
            .help("Use max scoring alignment algorithm rather than pair HMM forward algorithm.")
            .display_order(165))
        .arg(Arg::with_name("Debug Allele Realignment")
            .short("u")
            .long("debug_allele_realignment")
            .hidden(true)
            .help("Do NOT call variants, only print out meta-stats (number of expected operations etc) for allele realignment.")
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

    //let potential_variants_file: Option<&str> = input_args.value_of("Potential Variants VCF");
    let hap_bam_prefix: Option<&str> = input_args.value_of("Haplotype Bam Prefix");

    let force = match input_args.occurrences_of("Force overwrite") {
        0 => false,
        1 => true,
        _ => {
            panic!("force_overwrite specified multiple times");
        }
    };

    let no_haps = match input_args.occurrences_of("No haplotypes") {
        0 => false,
        1 => true,
        _ => {
            panic!("no_haps specified multiple times");
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

    let max_cigar_indel: usize = input_args.value_of("Max CIGAR indel")
        .unwrap()
        .parse::<usize>()
        .expect("Argument max_cigar_indel must be a positive integer!");

    let min_allele_qual: f64 = input_args.value_of("Min allele quality")
        .unwrap()
        .parse::<f64>()
        .expect("Argument max_mec_frac must be a positive float!");

    if min_allele_qual <= 0.0 {
        panic!("Min allele quality must be a positive float.");
    }

    let hap_assignment_qual: f64 = input_args.value_of("Haplotype Assignment Quality")
        .unwrap()
        .parse::<f64>()
        .expect("Argument max_mec_frac must be a positive float!");

    if hap_assignment_qual <= 0.0 {
        panic!("Haplotype assignment quality must be a positive float.");
    }

    let max_p_miscall: f64 = *Prob::from(PHREDProb(min_allele_qual));

    let hap_max_p_misassign: f64 = *Prob::from(PHREDProb(min_allele_qual));

    let ll_delta: f64 = input_args.value_of("Haplotype Convergence Delta")
        .unwrap()
        .parse::<f64>()
        .expect("Argument hap_converge_delta must be a positive float!");

    let potential_snv_cutoff: LogProb = LogProb::from(PHREDProb(input_args.value_of("Potential SNV Cutoff")
        .unwrap()
        .parse::<f64>()
        .expect("Argument potential_snv_cutoff must be a positive float!")));

    let hom_snv_rate: LogProb = LogProb::from(Prob(input_args.value_of("Homozygous SNV Rate")
        .unwrap()
        .parse::<f64>()
        .expect("Argument hom_snv_rate must be a positive float!")));

    let het_snv_rate: LogProb = LogProb::from(Prob(input_args.value_of("Heterozygous SNV Rate")
        .unwrap()
        .parse::<f64>()
        .expect("Argument het_snv_rate must be a positive float!")));

    let hom_indel_rate: LogProb = LogProb::from(Prob(input_args.value_of("Homozygous Indel Rate")
        .unwrap()
        .parse::<f64>()
        .expect("Argument hom_indel_rate must be a positive float!")));

    let het_indel_rate: LogProb = LogProb::from(Prob(input_args.value_of("Heterozygous Indel Rate")
        .unwrap()
        .parse::<f64>()
        .expect("Argument het_indel_rate must be a positive float!")));

    // multiply by 2.0 because internally we use this value as the probability of transition to
    // a single transversion base, not the combined probability of transversion to either one
    let ts_tv_ratio = 2.0 * input_args.value_of("ts/tv Ratio")
        .unwrap()
        .parse::<f64>()
        .expect("Argument ts_tv_ratio must be a positive float!");

    let dn_params = input_args.value_of("Density parameters").unwrap().split(":").collect::<Vec<&str>>();

    if dn_params.len() != 3 {
        panic!("Format for density params should be <n>:<l>:<gq>, with all 3 values being integers.");
    }

    let dn_count = dn_params[0].parse::<usize>().expect("Format for density params should be <n>:<l>:<gq>, with all 3 values being integers.");
    let dn_len = dn_params[1].parse::<usize>().expect("Format for density params should be <n>:<l>:<gq>, with all 3 values being integers.");
    let dn_gq = dn_params[2].parse::<usize>().expect("Format for density params should be <n>:<l>:<gq>, with all 3 values being integers.");

    let density_params = DensityParameters{n: dn_count, len: dn_len, gq: dn_gq as f64};
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

    let mut alignment_type = AlignmentType::FastAllAlignment;

    match input_args.occurrences_of("Numerically stable alignment") {
        0 => {},
        1 => {alignment_type = AlignmentType::NumericallyStableAllAlignment;},
        _ => {
            panic!("stable_alignment specified multiple times");
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

    let debug_allele_realignment: bool = match input_args.occurrences_of("Debug Allele Realignment") {
        0 => {false},
        1 => {true},
        _ => {panic!("debug_allele_realignment specified multiple times");}
    };

    let band_width: usize = input_args.value_of("Band width")
        .unwrap()
        .parse::<usize>()
        .expect("Argument band_width must be a positive integer!");

    /*
    let use_poa = match input_args.occurrences_of("Use POA") {
        0 => false,
        1 => true,
        _ => {
            panic!("Use POA specified multiple times");
        }
    };*/

    let min_cov: u32 = input_args.value_of("Min coverage").unwrap().parse::<u32>().expect("Argument min_cov must be a positive integer!");

    let max_cov: u32 = match input_args.occurrences_of("Auto max coverage") {
        0 => {
            // manually assigned coverage cutoff from user
             input_args.value_of("Max coverage").unwrap().parse::<u32>().expect("Argument max_cov must be a positive integer!")
        },
        1 => {
            eprintln!("{} Automatically determining max read coverage.",print_time());
            eprintln!("{} Estimating mean read coverage...",print_time());
            let mean_coverage: f64 = calculate_mean_coverage(&bamfile_name, &interval);
            let calculated_max_cov = (mean_coverage as f64 + 5.0 * (mean_coverage as f64).sqrt()) as u32;

            eprintln!("{} Mean read coverage: {:.2}",print_time(), mean_coverage);

            calculated_max_cov
        },
        _ => {
            panic!("auto_min_max_cov specified multiple times");
        }
    };

    if max_cov > 0 {
        eprintln!("{} Min read coverage set to {}.", print_time(), min_cov);
        eprintln!("{} Max read coverage set to {}.", print_time(), max_cov);
    } else {
        eprintln!("{} ERROR: Max read coverage set to 0.",print_time());
        return;
    }

    let extract_fragment_parameters = ExtractFragmentParameters {
        min_mapq: min_mapq,
        alignment_type: alignment_type,
        band_width: band_width,
        anchor_length: anchor_length,
        short_hap_max_snvs: short_hap_max_snvs,
        max_window_padding: max_window_padding,
        max_cigar_indel: max_cigar_indel
    };

    eprintln!("{} Estimating alignment parameters...",print_time());
    let alignment_parameters = estimate_alignment_parameters(&bamfile_name, &fasta_file, &interval, min_mapq, max_cigar_indel as u32);
    /*let alignment_parameters = AlignmentParameters{
        transition_probs: TransitionProbs {
            match_from_match: 0.879,
            insertion_from_match: 0.080,
            deletion_from_match: 0.041,
            insertion_from_insertion: 0.240,
            match_from_insertion: 0.760,
            deletion_from_deletion: 0.113,
            match_from_deletion: 0.887,
        },
        emission_probs: EmissionProbs {
            equal: 0.982,
            not_equal: 0.006,
            insertion: 1.0,
            deletion: 1.0
        }
    };*/

    /***********************************************************************************************/
    // GET HUMAN GENOTYPE PRIORS
    /***********************************************************************************************/

    let genotype_priors = GenotypePriors::new(hom_snv_rate, het_snv_rate,
                                                            hom_indel_rate, het_indel_rate,
                                                            ts_tv_ratio);

    /***********************************************************************************************/
    // FIND INITIAL SNVS WITH READ PILEUP
    /***********************************************************************************************/

    //let bam_file: String = "test_data/test.bam".to_string();
    eprintln!("{} Calling potential SNVs using pileup...",print_time());
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
    let mut varlist = call_potential_snvs::call_potential_snvs(&bamfile_name,
                                             &fasta_file,
                                             &interval,
                                             &genotype_priors,
                                             min_cov,
                                             max_cov,
                                             min_mapq,
                                             max_p_miscall,
                                             alignment_parameters.ln(),
                                             potential_snv_cutoff);

    // back up the variant indices
    // they will be needed later when we try to re-use fragment alleles that don't change
    // as the variant list expands
    varlist.backup_indices();

    print_variant_debug(&mut varlist, &interval, &variant_debug_directory,&"1.0.potential_SNVs.vcf", max_cov, &density_params, &sample_name);

    eprintln!("{} {} potential SNVs identified.", print_time(),varlist.lst.len());

    if varlist.lst.len() == 0 {
        return;
    }

    if debug_allele_realignment {

        extract_fragments_debug::extract_fragments_debug(&bamfile_name,
                                             &fasta_file,
                                             &varlist,
                                             &interval,
                                             extract_fragment_parameters,
                                             alignment_parameters,
                                             None);

        eprintln!("{} Allele realignment debugging complete. Exiting...",print_time());
        return;
    }
    /***********************************************************************************************/
    // EXTRACT FRAGMENT INFORMATION FROM READS
    /***********************************************************************************************/

    eprintln!("{} Generating condensed read data for SNVs...",print_time());
    let mut flist = extract_fragments::extract_fragments(&bamfile_name,
                                                         &fasta_file,
                                                         &mut varlist,
                                                         &interval,
                                                         extract_fragment_parameters,
                                                         alignment_parameters,
                                                          None);


    //for f in 0..flist.len() {
    //    for call in &flist[f].calls {
    //        if varlist.lst[call.var_ix].chrom == "chr1".to_string() && varlist.lst[call.var_ix].pos0 == 12067043 {
    //            println!("{}\t{}\t{:.2}", &flist[f].id, varlist.lst[call.var_ix].alleles[call.allele as usize], *PHREDProb::from(call.qual));
    //        }
    //    }
    //}
    //return;

    // if we're printing out variant "debug" information, print out a fragment file to that debug directory
    match &variant_debug_directory {
        &Some(ref debug_dir) => {
            // normally phase_variant is used to select which variants are heterozygous, so that
            // we only pass to HapCUT2 heterozygous variants
            // in this case, we set them all to 1 so we generate fragments for all variants
            let ffn = match Path::new(&debug_dir).join(&"fragments.txt").to_str() {
                Some(s) => {s.to_owned()},
                None => {panic!("Invalid unicode provided for variant debug directory");}
            };
            let phase_variant: Vec<bool> = vec![true; varlist.lst.len()];
            let mut fragment_buffer = generate_flist_buffer(&flist, &phase_variant, max_p_miscall);

            let fragment_file_path = Path::new(&ffn);
            let fragment_file_display = fragment_file_path.display();
            let mut fragment_file = match File::create(&fragment_file_path) {
                Err(why) => panic!("couldn't create {}: {}", fragment_file_display, why.description()),
                Ok(file) => file,
            };
            for mut line_u8 in fragment_buffer {
                line_u8.pop();
                match writeln!(fragment_file, "{}", u8_to_string(&line_u8)) {
                    Err(why) => panic!("couldn't write to {}: {}", fragment_file_display, why.description()),
                    Ok(_) => {}
                }
            }
        },
        &None => {}
    }

    /***********************************************************************************************/
    // CALL GENOTYPES USING REFINED QUALITY SCORES
    /***********************************************************************************************/

    eprintln!("{} Calling initial genotypes using pair-HMM realignment...", print_time());
    call_genotypes_no_haplotypes(&flist, &mut varlist, &genotype_priors, max_p_miscall);

    print_variant_debug(&mut varlist, &interval, &variant_debug_directory,&"2.0.realigned_genotypes.vcf", max_cov, &density_params, &sample_name);


    /***********************************************************************************************/
    // ITERATIVELY ASSEMBLE HAPLOTYPES AND CALL GENOTYPES
    /***********************************************************************************************/
    if no_haps {
        print_vcf(&mut varlist, &interval, &output_vcf_file, false, max_cov, &density_params, &sample_name, false);
        return;
    }
    eprintln!("{} Iteratively assembling haplotypes and refining genotypes...",print_time());
    call_genotypes_with_haplotypes(&mut flist, &mut varlist, &interval, &genotype_priors,
                                   &variant_debug_directory, 3, max_cov, &density_params, max_p_miscall,
                                   &sample_name, ll_delta);

    /*
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

        print_variant_debug(&mut varlist, &interval, &variant_debug_directory,&"4.0.new_potential_SNVs_after_POA.vcf", max_cov, &density_params, &sample_name);

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
        print_variant_debug(&mut varlist, &interval, &variant_debug_directory,&"5.0.realigned_genotypes_after_POA.vcf", max_cov, &density_params, &sample_name);

        eprintln!("{} Iteratively assembling haplotypes and refining genotypes (with POA variants)...",print_time());
        call_genotypes_with_haplotypes(&mut flist2, &mut varlist, &interval, &genotype_priors,
            &variant_debug_directory, 6, max_cov, max_p_miscall, &sample_name, ll_delta);

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
    */

    calculate_mec(&flist, &mut varlist, max_p_miscall);

    let debug_filename = "4.0.final_genotypes.vcf";

    eprintln!("{} Calculating fraction of reads assigned to either haplotype...",print_time());
    let (h1,h2) = separate_fragments_by_haplotype(&flist,
                                                  LogProb::from(Prob(1.0 - hap_max_p_misassign)));
    match hap_bam_prefix {
        Some(p) => {
            eprintln!("{} Writing haplotype-assigned reads to bam files...",print_time());
            separate_bam_reads_by_haplotype(&bamfile_name, &interval, p.to_string(), &h1, &h2, min_mapq);
        },
        None => {}
    }


    eprintln!("{} Printing VCF file...",print_time());
    print_variant_debug(&mut varlist, &interval,&variant_debug_directory, &debug_filename, max_cov, &density_params, &sample_name);
    print_vcf(&mut varlist, &interval, &output_vcf_file, false, max_cov, &density_params, &sample_name, false);
}
