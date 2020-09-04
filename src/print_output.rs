//! Print Longshot output in VCF format

use bio::io::fasta;
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
    fasta_file: &Option<String>,
    output_vcf_file: &String,
    print_reference_genotype: bool,
    max_cov: u32,
    density_params: &DensityParameters,
    sample_name: &String,
    print_outside_region: bool,
    used_potential_variants_vcf: bool
) -> Result<()> {
    // first, add filter flags for variant density
    var_filter(
        varlist,
        density_params.gq,
        density_params.len,
        density_params.n,
        max_cov,
    );

    let mut fasta = match fasta_file {
        &Some(ref ff) => Some(
            fasta::IndexedReader::from_file(&ff).chain_err(|| ErrorKind::IndexedFastaOpenError)?,
        ),
        None => None,
    };

    let mut ref_seq: Vec<char> = vec![]; // this vector will be used to hold the reference sequence
    let mut prev_tid = 4294967295;

    let vcf_path = Path::new(output_vcf_file);
    let vcf_display = vcf_path.display();
    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = File::create(&vcf_path)
        .chain_err(|| ErrorKind::CreateFileError(vcf_display.to_string()))?;

    // first part of the header
    let headerstr1 = "##fileformat=VCFv4.2
##source=Longshot v0.4.0
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth of reads passing MAPQ filter\">
##INFO=<ID=AC,Number=R,Type=Integer,Description=\"Number of Observations of Each Allele\">
##INFO=<ID=AM,Number=1,Type=Integer,Description=\"Number of Ambiguous Allele Observations\">
##INFO=<ID=MC,Number=1,Type=Integer,Description=\"Minimum Error Correction (MEC) for this single variant\">
##INFO=<ID=MF,Number=1,Type=Float,Description=\"Minimum Error Correction (MEC) Fraction for this variant.\">
##INFO=<ID=MB,Number=1,Type=Float,Description=\"Minimum Error Correction (MEC) Fraction for this variant's haplotype block.\">
##INFO=<ID=AQ,Number=1,Type=Float,Description=\"Mean Allele Quality value (PHRED-scaled).\">
##INFO=<ID=GM,Number=1,Type=Integer,Description=\"Phased genotype matches unphased genotype (boolean).\">".to_string();

    writeln!(file, "{}", headerstr1)
        .chain_err(|| ErrorKind::FileWriteError(vcf_display.to_string()))?;

    // these header lines are not printed if a potential variants input VCF is used
    if !used_potential_variants_vcf {
        let headerstr2 = "##INFO=<ID=DA,Number=1,Type=Integer,Description=\"Total Depth of reads at any MAPQ (but passing samtools filter 0xF00).\">
##INFO=<ID=MQ10,Number=1,Type=Float,Description=\"Fraction of reads (passing 0xF00) with MAPQ>=10.\">
##INFO=<ID=MQ20,Number=1,Type=Float,Description=\"Fraction of reads (passing 0xF00) with MAPQ>=20.\">
##INFO=<ID=MQ30,Number=1,Type=Float,Description=\"Fraction of reads (passing 0xF00) with MAPQ>=30.\">
##INFO=<ID=MQ40,Number=1,Type=Float,Description=\"Fraction of reads (passing 0xF00) with MAPQ>=40.\">
##INFO=<ID=MQ50,Number=1,Type=Float,Description=\"Fraction of reads (passing 0xF00) with MAPQ>=50.\">".to_string();

        writeln!(file, "{}", headerstr2)
            .chain_err(|| ErrorKind::FileWriteError(vcf_display.to_string()))?;
    }
    // last part of the header
    let headerstr3 = format!("##INFO=<ID=PH,Number=G,Type=Integer,Description=\"PHRED-scaled Probabilities of Phased Genotypes\">
##INFO=<ID=SC,Number=1,Type=String,Description=\"Reference Sequence in 21-bp window around variant.\">
##FILTER=<ID=dn,Description=\"In a dense cluster of variants\">
##FILTER=<ID=dp,Description=\"Exceeds maximum depth\">
##FILTER=<ID=sb,Description=\"Allelic strand bias\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">
##FORMAT=<ID=UG,Number=1,Type=String,Description=\"Unphased Genotype (pre-haplotype-assembly)\">
##FORMAT=<ID=UQ,Number=1,Type=Float,Description=\"Unphased Genotype Quality (pre-haplotype-assembly)\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", sample_name);


    writeln!(file, "{}", headerstr3)
        .chain_err(|| ErrorKind::FileWriteError(vcf_display.to_string()))?;

    for var in &varlist.lst {
        assert!(var.alleles.len() >= 2);
        assert!(var.allele_counts.len() == var.alleles.len());
        assert!(var.genotype_post.n_alleles() == var.alleles.len());

        match interval {
            &Some(ref iv) => {
                if !print_outside_region
                    && (var.tid != iv.tid
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

        match fasta {
            Some(ref mut fa) => {
                // if we're on a different contig/chrom, we need to read in the sequence for that
                // contig/chrom from the FASTA into the ref_seq vector
                if var.tid != prev_tid {
                    let mut ref_seq_u8: Vec<u8> = vec![];
                    fa.fetch_all(&varlist.target_names[var.tid as usize])
                        .chain_err(|| ErrorKind::IndexedFastaReadError)?;
                    fa.read(&mut ref_seq_u8)
                        .chain_err(|| ErrorKind::IndexedFastaReadError)?;
                    ref_seq = dna_vec(&ref_seq_u8);
                }
                prev_tid = var.tid;
            }
            None => {}
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
        let mut post_vec: Vec<String> = vec![];
        for g in var.possible_genotypes() {
            let mut p = *PHREDProb::from(var.genotype_post.get(g));
            if p > 500.0 {
                p = 500.0;
            }
            post_vec.push(format!("{:.2}", p));
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
        ]
        .join("/");

        let genotypes_match: usize = (var.genotype == var.unphased_genotype
            || Genotype(var.genotype.1, var.genotype.0) == var.unphased_genotype)
            as usize;

        // we want to save the sequence context (21 bp window around variant on reference)
        // this will be printed to the VCF later and may help diagnose variant calling
        // issues e.g. if the variant occurs inside a large homopolymer or etc.
        let sequence_context: String = match fasta {
            Some(_) => {
                // get the position 10 bases to the left
                let l_window = if var.pos0 >= 10 {
                    var.pos0 as usize - 10
                } else {
                    0
                };
                // get the position 11 bases to the right
                let mut r_window = var.pos0 as usize + 11;
                if r_window >= ref_seq.len() {
                    r_window = ref_seq.len();
                }

                (ref_seq[l_window..r_window]).iter().collect::<String>()
            }
            None => "None".to_string(),
        };

        write!(file,
                       "{}\t{}\t.\t{}\t{}\t{:.0}\t{}\tDP={};AC={};AM={};MC={};MF={:.3};MB={:.3};AQ={:.2};GM={};",
                       varlist.target_names[var.tid as usize],
                       var.pos0 + 1,
                       var.alleles[0],
                       var_alleles.join(","),
                       var.qual+0.5, // round off to integer, 09/04/2020
                       var.filter,
                       var.dp,
                       allele_counts_str,
                       var.ambiguous_count,
                       var.mec,
                       var.mec_frac_variant,
                       var.mec_frac_block,
                       var.mean_allele_qual,
                       genotypes_match).chain_err(|| ErrorKind::FileWriteError(vcf_display.to_string()))?;

        if !used_potential_variants_vcf {
            write!(file,
                     "DA={};MQ10={:.2};MQ20={:.2};MQ30={:.2};MQ40={:.2};MQ50={:.2};",
                     var.dp_any_mq,
                     var.mq10_frac,
                     var.mq20_frac,
                     var.mq30_frac,
                     var.mq40_frac,
                     var.mq50_frac).chain_err(|| ErrorKind::FileWriteError(vcf_display.to_string()))?;
        }
        writeln!(file,
                 "PH={};SC={};\tGT:GQ:PS:UG:UQ\t{}:{:.2}:{}:{}:{:.2}",
                 post_str,
                 sequence_context,
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
                &None,
                &outfile,
                true,
                max_cov,
                density_params,
                sample_name,
                true,
                true // don't print MQ statistics in VCF because they may or may not be present
            )
            .chain_err(|| "Error printing debug VCF file.")?;
        }
        &None => {}
    };
    Ok(())
}
