
use bio::stats::{PHREDProb};
use util::*; //{MAX_VCF_QUAL, ln_sum_matrix, GenotypePriors, VarList, Fragment, FragCall, GenomicInterval};
use variants_and_fragments::VarList;
use std::error::Error;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;

pub fn print_vcf(varlist: &VarList, interval: &Option<GenomicInterval>, indels: bool, output_vcf_file: &String, print_whole_varlist: bool) {
    let vcf_path = Path::new(output_vcf_file);
    let vcf_display = vcf_path.display();
    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match File::create(&vcf_path) {
        Err(why) => panic!("couldn't create {}: {}", vcf_display, why.description()),
        Ok(file) => file,
    };

    let headerstr = "##fileformat=VCFv4.2
##source=ReaperV0.1
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Number of Observations of Each Allele\">
##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Number of Ambiguous Alleles\">
##INFO=<ID=PH,Number=1,Type=Integer,Description=\"Phred-scaled Probabilities of Phased Genotypes\">
##INFO=<ID=MEC,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Score for Variant\">
##INFO=<ID=MF,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Fraction for Variant\">
##INFO=<ID=PSMF,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Fraction for Phase Set\">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">
##FORMAT=<ID=GQ,Number=2,Type=Float,Description=\"Genotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        .to_string();

    match writeln!(file, "{}", headerstr) {
        Err(why) => panic!("couldn't write to {}: {}", vcf_display, why.description()),
        Ok(_) => {}
    }


    for var in &varlist.lst {

        assert!(var.alleles.len() >= 2);
        assert!(var.genotype.len() >= 2);
        assert!(var.allele_counts.len() == var.alleles.len());
        assert!(var.genotype_post.size() == var.alleles.len());

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

        if !print_whole_varlist {
            if var.genotype == [0u8, 0u8] {
                continue;
            }

            if !indels {
                let mut found_indel = false;

                for ref allele in &var.alleles {
                    if allele.len() > 1 {
                        found_indel = true;
                        break;
                    }
                }

                if found_indel {
                    continue;
                }
            }
        }

        let ps = match var.phase_set {
            Some(ps) => format!("{}", ps),
            None => ".".to_string()
        };

        let allele_counts_str = var.allele_counts.clone().iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",");
        let mut post_vec: Vec<String> = vec![];
        for g1 in 0..var.alleles.len() {
            for g2 in 0..var.alleles.len(){
                post_vec.push(format!("{:.2}", *PHREDProb::from(var.genotype_post.get(g1,g2))));
            }
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

        let genotype_str = vec![var.genotype[0].to_string(), var.genotype[1].to_string()].join(sep);

        match writeln!(file,
                       "{}\t{}\t.\t{}\t{}\t{:.2}\t{}\tDP={};AC={};NA={};PH={:.2};MEC={};MF={};PSMF={:.5}\tGT:PS:GQ\t{}:{}:{:.2}",
                       var.chrom,
                       var.pos0 + 1,
                       var.alleles[0],
                       var_alleles.join(","),
                       var.qual,
                       var.filter,
                       var.dp,
                       allele_counts_str,
                       var.ambiguous_count,
                       post_str,
                       0,
                       0,
                       0,
                       genotype_str,
                       ps,
                       var.gq) {
            Err(why) => panic!("couldn't write to {}: {}", vcf_display, why.description()),
            Ok(_) => {}
        }
    }
}

pub fn print_variant_debug(varlist: &VarList, interval: &Option<GenomicInterval>, variant_debug_directory: &Option<String>, debug_filename: &str){
    match variant_debug_directory {
        &Some(ref dir) => {
            let outfile = match Path::new(&dir).join(&debug_filename).to_str() {
                Some(s) => {s.to_owned()},
                None => {panic!("Invalid unicode provided for variant debug directory");}
            };
            print_vcf(&varlist, &interval, true, &outfile, true);
        }
        &None => {}
    };
}
