extern crate bio;
extern crate rust_htslib;

mod call_potential_snvs;

fn main() {
    let bam_file: String = "test_data/test.bam".to_string();
    let fasta_file: String = "/home/peter/data/hg19_pacbio_chr20.fa".to_string();
    call_potential_snvs::call_potential_snvs(bam_file, fasta_file);
}
