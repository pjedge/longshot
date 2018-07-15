
use bio::stats::{PHREDProb};
use util::*; //{MAX_VCF_QUAL, ln_sum_matrix, GenotypePriors, VarList, Fragment, FragCall, GenomicInterval};
use variants_and_fragments::{VarList, var_filter};
use genotype_probs::Genotype;
use std::error::Error;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;


pub fn print_vcf(varlist: &mut VarList, interval: &Option<GenomicInterval>, output_vcf_file: &String, print_reference_genotype: bool, max_cov: Option<u32>, sample_name: &String) {

    // first, add filter flags for variant density
    var_filter(varlist, 50.0, 500, 10, max_cov);

    let vcf_path = Path::new(output_vcf_file);
    let vcf_display = vcf_path.display();
    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match File::create(&vcf_path) {
        Err(why) => panic!("couldn't create {}: {}", vcf_display, why.description()),
        Ok(file) => file,
    };

    let headerstr = format!("##fileformat=VCFv4.2
##source=ReaperV0.1
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AC,Number=R,Type=Integer,Description=\"Number of Observations of Each Allele\">
##INFO=<ID=AM,Number=1,Type=Integer,Description=\"Number of Ambiguous Allele Observations\">
##INFO=<ID=MC,Number=1,Type=Integer,Description=\"Minimum Error Correction (MEC) for this single variant\">
##INFO=<ID=MF,Number=1,Type=Float,Description=\"Minimum Error Correction (MEC) Fraction for this variant.\">
##INFO=<ID=MF,Number=1,Type=Float,Description=\"Minimum Error Correction (MEC) Fraction for this variant's haplotype block.\">
##INFO=<ID=AQ,Number=1,Type=Float,Description=\"Mean Allele Quality value (PHRED-scaled).\">
##INFO=<ID=PH,Number=G,Type=Integer,Description=\"PHRED-scaled Probabilities of Phased Genotypes\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", sample_name);

    //"##INFO=<ID=MEC,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Score for Variant\">
    //##INFO=<ID=MF,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Fraction for Variant\">
    //##INFO=<ID=PSMF,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Fraction for Phase Set\">"
    match writeln!(file, "{}", headerstr) {
        Err(why) => panic!("couldn't write to {}: {}", vcf_display, why.description()),
        Ok(_) => {}
    }


    for var in &varlist.lst {

        assert!(var.alleles.len() >= 2);
        assert!(var.allele_counts.len() == var.alleles.len());
        assert!(var.genotype_post.n_alleles() == var.alleles.len());

        match interval {
            &Some(ref iv) => {
                if var.chrom != iv.chrom ||
                    var.pos0 < iv.start_pos as usize ||
                    var.pos0 > iv.end_pos as usize {
                    continue;
                }
            }
            &None => {}
        }

        if !print_reference_genotype {
            if var.genotype == Genotype(0,0) {
                continue;
            }
        }

        let ps = match var.phase_set {
            Some(ps) => format!("{}", ps),
            None => ".".to_string()
        };

        let allele_counts_str = var.allele_counts.clone().iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",");
        let mut post_vec: Vec<String> = vec![];
        for g in var.possible_genotypes() {
            post_vec.push(format!("{:.2}", *PHREDProb::from(var.genotype_post.get(g))));
        }

        let post_str = post_vec.join(",");

        let mut var_alleles: Vec<String> = vec![];
        for allele in var.alleles[1..].iter() {
            var_alleles.push(allele.clone());
        }

        assert!(var_alleles.len() == var.alleles.len() - 1);

        let sep = match var.phase_set {
            Some(_) => "|",
            None => "/"
        };

        let genotype_str = vec![var.genotype.0.to_string(), var.genotype.1.to_string()].join(sep);


        match writeln!(file,
                       "{}\t{}\t.\t{}\t{}\t{:.2}\t{}\tDP={};AC={};AM={};MC={};MF={:.3};MB={:.3};AQ={:.2};PH={};\tGT:PS:GQ\t{}:{}:{:.2}",
                       var.chrom,
                       var.pos0 + 1,
                       var.alleles[0],
                       var_alleles.join(","),
                       var.qual,
                       var.filter,
                       var.dp,
                       allele_counts_str,
                       var.ambiguous_count,
                       var.mec,
                       var.mec_frac_variant,
                       var.mec_frac_block,
                       var.mean_allele_qual,
                       post_str,
                       genotype_str,
                       ps,
                       var.gq) {
            Err(why) => panic!("couldn't write to {}: {}", vcf_display, why.description()),
            Ok(_) => {}
        }
    }
}

pub fn print_variant_debug(varlist: &mut VarList, interval: &Option<GenomicInterval>, variant_debug_directory: &Option<String>, debug_filename: &str, max_cov: Option<u32>, sample_name: &String){
    match variant_debug_directory {
        &Some(ref dir) => {
            let outfile = match Path::new(&dir).join(&debug_filename).to_str() {
                Some(s) => {s.to_owned()},
                None => {panic!("Invalid unicode provided for variant debug directory");}
            };
            print_vcf(varlist, &interval,&outfile, true, max_cov, sample_name);
        }
        &None => {}
    };
}
