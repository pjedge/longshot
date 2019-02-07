//! Print Longshot output in VCF format

use bio::stats::PHREDProb;
use errors::*;
use genotype_probs::Genotype;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use util::*; //{MAX_VCF_QUAL, ln_sum_matrix, GenotypePriors, VarList, Fragment, FragCall, GenomicInterval};
use variants_and_fragments::{var_filter, VarList};

pub fn print_vcf(
    varlist: &mut VarList,
    interval: &Option<GenomicInterval>,
    output_vcf_file: &String,
    print_reference_genotype: bool,
    max_cov: u32,
    density_params: &DensityParameters,
    sample_name: &String,
    print_outside_region: bool,
) -> Result<()> {
    // first, add filter flags for variant density
    var_filter(
        varlist,
        density_params.gq,
        density_params.len,
        density_params.n,
        max_cov,
    );

    let vcf_path = Path::new(output_vcf_file);
    let vcf_display = vcf_path.display();
    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = File::create(&vcf_path)
        .chain_err(|| ErrorKind::CreateFileError(vcf_display.to_string()))?;

    let headerstr = format!("##fileformat=VCFv4.2
##source=Longshot v0.2.1
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth of reads passing MAPQ filter\">
##INFO=<ID=AC,Number=R,Type=Integer,Description=\"Number of Observations of Each Allele\">
##INFO=<ID=AM,Number=1,Type=Integer,Description=\"Number of Ambiguous Allele Observations\">
##INFO=<ID=FC,Number=R,Type=Integer,Description=\"Number of Observations of Each Allele on the Forward Strand\">
##INFO=<ID=RC,Number=R,Type=Integer,Description=\"Number of Observations of Each Allele on the Reverse Strand\">
##INFO=<ID=MC,Number=1,Type=Integer,Description=\"Minimum Error Correction (MEC) for this single variant\">
##INFO=<ID=MF,Number=1,Type=Float,Description=\"Minimum Error Correction (MEC) Fraction for this variant.\">
##INFO=<ID=MB,Number=1,Type=Float,Description=\"Minimum Error Correction (MEC) Fraction for this variant's haplotype block.\">
##INFO=<ID=AQ,Number=1,Type=Float,Description=\"Mean Allele Quality value (PHRED-scaled).\">
##INFO=<ID=GM,Number=1,Type=Integer,Description=\"Phased genotype matches unphased genotype (boolean).\">
##INFO=<ID=DA,Number=1,Type=Integer,Description=\"Total Depth of reads at any MAPQ (but passing samtools filter 0xF00).\">
##INFO=<ID=MQ10,Number=1,Type=Float,Description=\"Fraction of reads (passing 0xF00) with MAPQ>=10.\">
##INFO=<ID=MQ20,Number=1,Type=Float,Description=\"Fraction of reads (passing 0xF00) with MAPQ>=20.\">
##INFO=<ID=MQ30,Number=1,Type=Float,Description=\"Fraction of reads (passing 0xF00) with MAPQ>=30.\">
##INFO=<ID=MQ40,Number=1,Type=Float,Description=\"Fraction of reads (passing 0xF00) with MAPQ>=40.\">
##INFO=<ID=MQ50,Number=1,Type=Float,Description=\"Fraction of reads (passing 0xF00) with MAPQ>=50.\">
##INFO=<ID=PH,Number=G,Type=Integer,Description=\"PHRED-scaled Probabilities of Phased Genotypes\">
##INFO=<ID=SC,Number=1,Type=String,Description=\"Reference Sequence in 21-bp window around variant.\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">
##FORMAT=<ID=UG,Number=1,Type=String,Description=\"Unphased Genotype (pre-haplotype-assembly)\">
##FORMAT=<ID=UQ,Number=1,Type=Float,Description=\"Unphased Genotype Quality (pre-haplotype-assembly)\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", sample_name);

    //"##INFO=<ID=MEC,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Score for Variant\">
    //##INFO=<ID=MF,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Fraction for Variant\">
    //##INFO=<ID=PSMF,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Fraction for Phase Set\">"
    writeln!(file, "{}", headerstr)
        .chain_err(|| ErrorKind::FileWriteError(vcf_display.to_string()))?;

    for var in &varlist.lst {
        assert!(var.alleles.len() >= 2);
        assert!(var.allele_counts.len() == var.alleles.len());
        assert!(var.genotype_post.n_alleles() == var.alleles.len());

        match interval {
            &Some(ref iv) => {
                if !print_outside_region
                    && (var.chrom != iv.chrom
                        || var.pos0 < iv.start_pos as usize
                        || var.pos0 > iv.end_pos as usize)
                {
                    continue;
                }
            }
            &None => {}
        }

        if !print_reference_genotype {
            if var.genotype == Genotype(0, 0) {
                continue;
            }
        }

        let ps = match var.phase_set {
            Some(ps) => format!("{}", ps),
            None => ".".to_string(),
        };

        let allele_counts_str = var
            .allele_counts
            .clone()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(",");

        let forward_counts_str = var
            .allele_counts_forward
            .clone()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(",");

        let reverse_counts_str = var
            .allele_counts_reverse
            .clone()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(",");

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
            None => "/",
        };

        let genotype_str = vec![var.genotype.0.to_string(), var.genotype.1.to_string()].join(sep);
        let unphased_genotype_str = vec![
            var.unphased_genotype.0.to_string(),
            var.unphased_genotype.1.to_string(),
        ].join("/");

        let genotypes_match: usize = (var.genotype == var.unphased_genotype
            || Genotype(var.genotype.1, var.genotype.0) == var.unphased_genotype)
            as usize;

        writeln!(file,
                       "{}\t{}\t.\t{}\t{}\t{:.2}\t{}\tDP={};AC={};AM={};FC={};RC={};MC={};MF={:.3};MB={:.3};AQ={:.2};GM={};DA={};MQ10={:.2};MQ20={:.2};MQ30={:.2};MQ40={:.2};MQ50={:.2};PH={};SC={};\tGT:GQ:PS:UG:UQ\t{}:{:.2}:{}:{}:{:.2}",
                       var.chrom,
                       var.pos0 + 1,
                       var.alleles[0],
                       var_alleles.join(","),
                       var.qual,
                       var.filter,
                       var.dp,
                       allele_counts_str,
                       var.ambiguous_count,
                       forward_counts_str,
                       reverse_counts_str,
                       var.mec,
                       var.mec_frac_variant,
                       var.mec_frac_block,
                       var.mean_allele_qual,
                       genotypes_match,
                       var.dp_any_mq,
                       var.mq10_frac,
                       var.mq20_frac,
                       var.mq30_frac,
                       var.mq40_frac,
                       var.mq50_frac,
                       post_str,
                       var.sequence_context,
                       genotype_str,
                       var.gq,
                       ps,
                       unphased_genotype_str,
                       var.unphased_gq).chain_err(|| ErrorKind::FileWriteError(vcf_display.to_string()))?;
    }
    Ok(())
}

pub fn print_variant_debug(
    varlist: &mut VarList,
    interval: &Option<GenomicInterval>,
    variant_debug_directory: &Option<String>,
    debug_filename: &str,
    max_cov: u32,
    density_params: &DensityParameters,
    sample_name: &String,
) -> Result<()> {
    match variant_debug_directory {
        &Some(ref dir) => {
            let outfile = match Path::new(&dir).join(&debug_filename).to_str() {
                Some(s) => s.to_owned(),
                None => {
                    bail!("Invalid unicode provided for variant debug directory");
                }
            };
            print_vcf(
                varlist,
                &interval,
                &outfile,
                true,
                max_cov,
                density_params,
                sample_name,
                true,
            ).chain_err(|| "Error printing debug VCF file.")?;
        }
        &None => {}
    };
    Ok(())
}
