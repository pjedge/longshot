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

use call_genotypes::call_genotypes;
use util::{u8_to_string, dna_vec};

fn main() {

    let bamfile_name = "test_data/test.bam".to_string();
    let fasta_file = "/home/peter/data/hg19_pacbio_chr20.fa".to_string();

    //let bam_file: String = "test_data/test.bam".to_string();
    let varlist =
        call_potential_snvs::call_potential_snvs(&bamfile_name, &fasta_file, 2, 1.0 / 8.0, 63, 30);

    let flist = extract_fragments::extract_fragments(&bamfile_name, &fasta_file, &varlist);

    call_genotypes(flist, &varlist);

}
