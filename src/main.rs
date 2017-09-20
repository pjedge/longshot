#![allow(dead_code)]

extern crate bio;
extern crate rust_htslib;
#[macro_use]
extern crate quick_error;

mod call_potential_snvs;
mod extract_fragments;
mod call_genotypes;
mod realignment;
mod util;
use bio::stats::{LogProb, Prob, PHREDProb};

use call_genotypes::call_genotypes;
use util::AlignmentParameters;
use std::env;

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

fn main() {

    let mut args: Vec<String> = vec![];
    for a in env::args() {
        args.push(a);
    }
    if args.len() < 3 {
        println!("usage: ./reaper [bamfile] [fasta_file] > [output.vcf]");
        return;
    }

    let bamfile_name = args[1].clone();
    let fasta_file = args[2].clone();

    //let bam_file: String = "test_data/test.bam".to_string();
    eprintln!("Calling potential SNVs using pileup...");
    let varlist = call_potential_snvs::call_potential_snvs(&bamfile_name,
                                                           &fasta_file,
                                                           2,
                                                           1.0 / 8.0,
                                                           5,
                                                           63,
                                                           30);

    eprintln!("Generating condensed read data for SNVs...");
    let flist = extract_fragments::extract_fragments(&bamfile_name,
                                                     &fasta_file,
                                                     &varlist,
                                                     PACBIO_ALIGNMENT_PARAMETERS);
    eprintln!("Calling Genotypes...");
    call_genotypes(flist, &varlist);

}
