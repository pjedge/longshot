//! This module contains the function used for identifying potential SNV sites.

extern crate rust_htslib;

use std::char;
use std::collections::HashMap;

use bio::io::fasta;
//,HashSet};
use bio::stats::LogProb;
use rust_htslib::bam;
use rust_htslib::bam::pileup::Indel;
use rust_htslib::prelude::*;
use rand::{Rng, SeedableRng, StdRng};

use call_genotypes::calculate_genotype_posteriors_no_haplotypes;
use errors::*;
//use std::str;
//use bio::alignment::Alignment;
//use bio::alignment::pairwise::banded::*;
//use bio::alignment::AlignmentOperation::*;
use genotype_probs::*;
//use spoa::poa_multiple_sequence_alignment;
use realignment::LnAlignmentParameters;
use util::*;
// {FragCall, GenotypePriors, LnAlignmentParameters, GenomicInterval, Var, VarList, parse_target_names, u8_to_string};
use variants_and_fragments::*;

pub static VARLIST_CAPACITY: usize = 1000000;
static VERBOSE: bool = false; //true;

/// Calls potential SNV sites using a pileup-based genotyping calculation
///
/// Potential SNVs are identified by performing a relatively standard pileup-based
/// genotyping calculation. Despite the fact that SMS read alignments around SNVs are inaccurate,
/// some of reads should be correctly aligned and indicate the presence of a SNV. So for every site
/// in the genome, we perform a genotyping calculation and consider every site meeting a low threshold
/// for variant evidence (by default, probability of non-reference genotype > 0.001) as a potential SNV
///
/// # Arguments
/// -```bam_file```: the input BAM file name as a string
/// -```fasta_file```: the input FASTA file name as a string
/// -```interval```: optional struct holding the genomic interval to call variants in
/// -```genotype_priors```: struct holding the genotype priors
/// -```min_coverage```: the minimum read coverage to consider a site as a potential variant
/// -```max_coverage```: the maximum read coverage to consider a site as a potential variant
/// -```min_mapq```: the minimum mapping quality to use a read in variant calling
/// -```max_p_miscall```: the maximum probability of an allele miscall to count the allele (equivalent
/////                     to the minimum allowed allele quality, but represented as a normal probability
/////                     rather than PHRED-scaled)
/// -```ln_align_params```: natural-log-scaled parameters for read alignment (Pair-HMM). Note that
///                         no alignment is performed in this function; these parameters are only
///                         required so that the base substitution rate can be used as interim
///                         "allele quality values" for the genotyping calculation.
/// -```potential_snv_cutoff```: natural-log-scaled cutoff for the probability of non-reference
///                              genotype (by default 0.001). Any site with probability of
///                              non-reference genotype greater than this amount will be kept and
///                              considered as a potential SNV site.
///
/// # Returns
/// Returns a result that wraps a VarList struct, representing the list of potential variants.
/// The data structure has fields corresponding to data that will be written to the output VCF file.
/// This variant list will be passed mutably to other functions that perform genotyping and
/// haplotype assembly steps.
///
/// #Errors
/// - ```IndexedFastaOpenError```: error opening the indexed FASTA file
/// - ```IndexedBamOpenError```: error opening the indexed BAM file
/// - ```IndexedBamFetchError```: error fetching a region from the indexed BAM file
/// - ```IndexedBamPileupReadError```: error retrieving pileup from indexed BAM file
/// - ```IndexedFastaReadError```: error reading an entry from the indexed FASTA file
/// - ```IndexedBamPileupQueryPositionError```: error accessing the query position from the BAM pileup
/// - Error calculating genotype posteriors for reference genotype qual calculation (usually having to
///          do with accessing invalid genotypes)
pub fn call_potential_snvs(
    bam_file: &String,
    fasta_file: &String,
    interval: &Option<GenomicInterval>,
    genotype_priors: &GenotypePriors,
    min_coverage: u32,
    max_coverage: u32,
    min_mapq: u8,
    max_p_miscall: f64,
    ln_align_params: LnAlignmentParameters,
    potential_snv_cutoff: LogProb,
) -> Result<VarList> {

    // the list of target (contig) names from the bam file
    let target_names = parse_target_names(&bam_file)?;

    let mut fasta = fasta::IndexedReader::from_file(&fasta_file)
        .chain_err(|| ErrorKind::IndexedFastaOpenError)?;

    let mut varlist: Vec<Var> = Vec::with_capacity(VARLIST_CAPACITY);

    // pileup over all covered sites
    let mut ref_seq: Vec<char> = vec![]; // this vector will be used to hold the reference sequence
    let mut prev_tid = 4294967295;

    // the 4 bases as String type
    let a_str = "A".to_string();
    let c_str = "C".to_string();
    let g_str = "G".to_string();
    let t_str = "T".to_string();

    // the strategy used for iterating over the BAM entries is as follows:
    // if a genomic region was specified, then ```get_interval_lst``` puts that genomic region into a vector (interval_lst)
    // containing only that one region. Then we iterate over the region list and seek every region,
    // which effectively means we just seek that single genomic interval.
    //
    // if a genomic region was not specified, then we want to iterate over the whole BAM file.
    // so get_interval_lst returns a list of genomic intervals that contain each whole chromosome
    // as described by the BAM header SQ lines. Then we iterate over the interval lst and seek each
    // interval separately, which effectively just iterates over all the BAM entries.
    //
    // the reason for this strange design (instead of either fetching a region beforehand or not and
    // then just iterating over all of ```bam_ix.pileup()```) is the following:
    // if an indexed reader is used, and fetch is never called, pileup() hangs.

    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bam_file, interval)?;
    let mut bam_ix =
        bam::IndexedReader::from_path(bam_file).chain_err(|| ErrorKind::IndexedBamOpenError)?;

    // interval_lst has either a single genomic interval (if --region was specified) or a list of
    // genomic intervals for each chromosome covering the entire genome
    for iv in interval_lst {
        bam_ix
            .fetch(iv.tid as u32, iv.start_pos as u32, iv.end_pos as u32 + 1)
            .chain_err(|| ErrorKind::IndexedBamFetchError)?;
        let bam_pileup = bam_ix.pileup();

        // this variable is used to avoid having a variant inside a previous variant's deletion.
        let mut next_valid_pos = 0;

        // iterate over base pileup
        for p in bam_pileup {
            let pileup = p.chain_err(|| ErrorKind::IndexedBamPileupReadError)?;

            let tid: usize = pileup.tid() as usize;
            let chrom: String = target_names[tid].clone();

            // if we're on a different contig/chrom, we need to read in the sequence for that
            // contig/chrom from the FASTA into the ref_seq vector
            if tid != prev_tid {
                let mut ref_seq_u8: Vec<u8> = vec![];
                fasta
                    .read_all(&chrom, &mut ref_seq_u8)
                    .chain_err(|| ErrorKind::IndexedFastaReadError)?;
                ref_seq = dna_vec(&ref_seq_u8);
                next_valid_pos = 0;
            }

            prev_tid = tid;

            // this is specifically to avoid having a variant inside a previous variant's deletion.
            if pileup.pos() < next_valid_pos {
                continue;
            }

            // we want to save the sequence context (21 bp window around variant on reference)
            // this will be printed to the VCF later and may help diagnose variant calling
            // issues e.g. if the variant occurs inside a large homopolymer or etc.
            // get the position 10 bases to the left
            let l_window = if pileup.pos() >= 10 {
                pileup.pos() as usize - 10
            } else {
                0
            };
            // get the position 11 bases to the right
            let mut r_window = pileup.pos() as usize + 11;
            if r_window >= ref_seq.len() {
                r_window = ref_seq.len();
            }
            let sequence_context: String = (ref_seq[l_window..r_window]).iter().collect::<String>();

            let ref_base_str = (ref_seq[pileup.pos() as usize]).to_string();

            if ref_base_str.contains("N") {
                continue;
            }

            // the dna_vec conversion function should remove any non-ACGT bases
            assert!(
                ref_base_str == a_str
                    || ref_base_str == c_str
                    || ref_base_str == g_str
                    || ref_base_str == t_str
            );

            //let mut counts = [0; 5]; // A,C,G,T,N
            let mut counts: HashMap<(String, String), usize> = HashMap::new();
            // use a counter instead of pileup.depth() since that would include qc_fail bases, low mapq, etc.
            let mut depth: usize = 0;
            let mut pileup_alleles: Vec<(String, String)> = vec![];
            let pos: usize = pileup.pos() as usize;

            let mut passing_reads: usize = 0; // total number of reads in this pileup at any mapq
            let mut mq10: f64 = 0.0; // total number of reads in this pileup with mapq >= 10
            let mut mq20: f64 = 0.0; // total number of reads in this pileup with mapq >= 20
            let mut mq30: f64 = 0.0; // total number of reads in this pileup with mapq >= 30
            let mut mq40: f64 = 0.0; // total number of reads in this pileup with mapq >= 40
            let mut mq50: f64 = 0.0; // total number of reads in this pileup with mapq >= 50

            // pileup the bases for a single position and count number of each base
            for alignment in pileup.alignments() {
                let record = alignment.record();

                // check that the read doesn't fail any standard filters
                if record.is_secondary()
                    || record.is_quality_check_failed()
                    || record.is_duplicate()
                    || record.is_supplementary()
                {
                    continue;
                }

                // we're keeping track of how many reads meeting different MAPQ cutoff we see
                // this information will be printed to the VCF and may help identify poor-mappability sites
                // iterate counters if read meets a MAPQ cutoff
                passing_reads += 1;
                if record.mapq() >= 10 {
                    mq10 += 1.0
                };
                if record.mapq() >= 20 {
                    mq20 += 1.0
                };
                if record.mapq() >= 30 {
                    mq30 += 1.0
                };
                if record.mapq() >= 40 {
                    mq40 += 1.0
                };
                if record.mapq() >= 50 {
                    mq50 += 1.0
                };

                if record.is_unmapped() || record.mapq() < min_mapq {
                    continue;
                }

                depth += 1; // depth counter does not consider unmapped low-mapq reads

                // handle the base/indel observed on the read
                if !alignment.is_del() && !alignment.is_refskip() {
                    let ref_allele;
                    let var_allele;

                    match alignment.indel() {
                        Indel::None => {
                            // read is NOT an indel (a base is observed)
                            ref_allele = (ref_seq[pos] as char).to_string().to_uppercase();

                            let base: char =
                                alignment.record().seq()
                                    [alignment.qpos().chain_err(|| {
                                        ErrorKind::IndexedBamPileupQueryPositionError
                                    })?] as char;

                            var_allele = base.to_string().to_uppercase();
                        }
                        Indel::Ins(l) => {
                            // read is an insertion
                            let start = alignment
                                .qpos()
                                .chain_err(|| ErrorKind::IndexedBamPileupQueryPositionError)?;
                            let end = start + l as usize + 1;
                            // don't want to convert whole seq to bytes...
                            let mut var_char: Vec<char> = vec![];

                            let record_len = alignment.record().seq().len();

                            for i in start..end {
                                if i < record_len {
                                    var_char.push(alignment.record().seq()[i] as char);
                                } else {
                                    var_char.push('N');
                                }
                            }

                            ref_allele = match ref_seq.get(pos) {
                                Some(&r) => (r as char).to_string().to_uppercase(),
                                None => "N".to_string(),
                            };

                            var_allele = var_char.into_iter().collect::<String>().to_uppercase();
                        }
                        Indel::Del(l) => {
                            // read is a deletion
                            let start = pos;
                            let end: usize = pos + l as usize + 1;

                            ref_allele = ref_seq[start..end].iter().collect();
                            var_allele = (ref_seq[pos] as char).to_string().to_uppercase();
                        }
                    }

                    // add the ref and var alleles to a vector representing the pileup
                    pileup_alleles.push((ref_allele.clone(), var_allele.clone()));
                    // iterate the counts for this observation in the counts hashmap
                    *counts
                        .entry((ref_allele.clone(), var_allele.clone()))
                        .or_insert(0) += 1;
                }
            }

            // these values representing read mappability will be saved and written to the VCF later
            let mq10_frac = mq10 / (passing_reads as f64);
            let mq20_frac = mq20 / (passing_reads as f64);
            let mq30_frac = mq30 / (passing_reads as f64);
            let mq40_frac = mq40 / (passing_reads as f64);
            let mq50_frac = mq50 / (passing_reads as f64);

            if depth < min_coverage as usize {
                continue;
            }

            if depth > max_coverage as usize {
                continue;
            }

            let mut snv_max_count: u32 = 0;
            let mut snv_ref_allele = 'N'.to_string();
            let mut snv_var_allele = 'N'.to_string();

            // iterate over all the counts in the hashmap.
            for (&(ref r, ref v), &count) in &counts {
                if has_non_acgt(&r) || has_non_acgt(&v) {
                    continue;
                }

                if r.len() == 1 && v.len() == 1 {
                    // potential SNV
                    if count > snv_max_count as usize && (r, v) != (&ref_base_str, &ref_base_str) {
                        snv_max_count = count as u32;
                        snv_ref_allele = r.clone();
                        snv_var_allele = v.clone();
                    }
                }
            }

            // vectors contain entries with (call, qual)
            // call is '0' or '1' (or potentially '2'...)
            // qual is a LogProb probability of miscall
            let mut snv_pileup_calls: Vec<FragCall> = vec![];
            //let mut indel_pileup_calls: Vec<(char, LogProb)> = vec![];

            // iterate over the observed alleles and count which ones match the most common one
            // (snv_var_allele)
            for (ref_allele, var_allele) in pileup_alleles {
                let a = if (ref_allele.clone(), var_allele.clone())
                    == (snv_ref_allele.clone(), snv_var_allele.clone())
                {
                    1u8
                } else {
                    0u8
                };

                let qual = ln_align_params.emission_probs.not_equal;

                let call = FragCall {
                    frag_ix: None,                                    // index into fragment list
                    var_ix: 0,                                        // index into variant list
                    allele: a,                                        // allele call
                    qual: qual, // LogProb probability the call is an error
                    one_minus_qual: LogProb::ln_one_minus_exp(&qual), // LogProb probability the call is correct
                };
                snv_pileup_calls.push(call);
            }

            // use a basic genotype likelihood calculation to call SNVs
            // snv_qual is the LogProb probability of a non-reference base observation
            let alleles = vec![snv_ref_allele.clone(), snv_var_allele.clone()];
            let snv_qual = if !snv_ref_allele.contains("N") && !snv_var_allele.contains("N") {
                let snv_post = calculate_genotype_posteriors_no_haplotypes(
                    &snv_pileup_calls,
                    &genotype_priors,
                    &alleles,
                    max_p_miscall,
                ).chain_err(|| "Error getting genotype posteriors for calling potential SNVs.")?;
                LogProb::ln_one_minus_exp(&snv_post.get(Genotype(0, 0)))
            } else {
                LogProb::ln_zero()
            };

            let (ref_allele, var_allele, qual) = (snv_ref_allele, snv_var_allele, snv_qual);

            next_valid_pos = (pos + 1) as u32;

            // check if SNV meets our quality criteria for a potential SNV
            // if it does, make a new variant and add it to the list of potential SNVs.
            if qual > potential_snv_cutoff
                && !ref_allele.contains("N")
                && !var_allele.contains("N")
                && (ref_allele.clone(), var_allele.clone())
                    != (ref_base_str.clone(), ref_base_str.clone())
            {
                let tid: usize = pileup.tid() as usize;
                let new_var = Var {
                    ix: 0,
                    old_ix: None,
                    // these will be set automatically,
                    tid: tid,
                    chrom: target_names[tid].clone(),
                    pos0: pos,
                    alleles: vec![ref_allele.clone(), var_allele.clone()],
                    dp: depth,
                    allele_counts: vec![0, 0],
                    ambiguous_count: 0,
                    qual: 0.0,
                    filter: ".".to_string(),
                    genotype: Genotype(0, 0),
                    gq: 0.0,
                    unphased_genotype: Genotype(0, 0),
                    unphased_gq: 0.0,
                    genotype_post: GenotypeProbs::uniform(2),
                    phase_set: None,
                    mec: 0,
                    mec_frac_variant: 0.0, // mec fraction for this variant
                    mec_frac_block: 0.0,   // mec fraction for this haplotype block
                    mean_allele_qual: 0.0,
                    dp_any_mq: passing_reads,
                    mq10_frac,
                    mq20_frac,
                    mq30_frac,
                    mq40_frac,
                    mq50_frac,
                    sequence_context,
                    called: false,
                };

                // we don't want potential SNVs that are inside a deletion, for instance.
                next_valid_pos = pileup.pos() + ref_allele.len() as u32;

                varlist.push(new_var);
            }
        }
    }
    // return the vector of Vars as a VarList struct
    Ok(VarList::new(varlist)?)
}


/// Calls potential SNV sites using a pileup-based genotyping calculation
///
/// Potential SNVs are identified by performing a relatively standard pileup-based
/// genotyping calculation. Despite the fact that SMS read alignments around SNVs are inaccurate,
/// some of reads should be correctly aligned and indicate the presence of a SNV. So for every site
/// in the genome, we perform a genotyping calculation and consider every site meeting a low threshold
/// for variant evidence (by default, probability of non-reference genotype > 0.001) as a potential SNV
///
/// # Arguments
/// -```bam_file```: the input BAM file name as a string
/// -```fasta_file```: the input FASTA file name as a string
/// -```interval```: optional struct holding the genomic interval to call variants in
/// -```genotype_priors```: struct holding the genotype priors
/// -```min_coverage```: the minimum read coverage to consider a site as a potential variant
/// -```max_coverage```: the maximum read coverage to consider a site as a potential variant
/// -```min_mapq```: the minimum mapping quality to use a read in variant calling
/// -```max_p_miscall```: the maximum probability of an allele miscall to count the allele (equivalent
/////                     to the minimum allowed allele quality, but represented as a normal probability
/////                     rather than PHRED-scaled)
/// -```ln_align_params```: natural-log-scaled parameters for read alignment (Pair-HMM). Note that
///                         no alignment is performed in this function; these parameters are only
///                         required so that the base substitution rate can be used as interim
///                         "allele quality values" for the genotyping calculation.
/// -```potential_snv_cutoff```: natural-log-scaled cutoff for the probability of non-reference
///                              genotype (by default 0.001). Any site with probability of
///                              non-reference genotype greater than this amount will be kept and
///                              considered as a potential SNV site.
///
/// # Returns
/// Returns a result that wraps a VarList struct, representing the list of potential variants.
/// The data structure has fields corresponding to data that will be written to the output VCF file.
/// This variant list will be passed mutably to other functions that perform genotyping and
/// haplotype assembly steps.
///
/// #Errors
/// - ```IndexedFastaOpenError```: error opening the indexed FASTA file
/// - ```IndexedBamOpenError```: error opening the indexed BAM file
/// - ```IndexedBamFetchError```: error fetching a region from the indexed BAM file
/// - ```IndexedBamPileupReadError```: error retrieving pileup from indexed BAM file
/// - ```IndexedFastaReadError```: error reading an entry from the indexed FASTA file
/// - ```IndexedBamPileupQueryPositionError```: error accessing the query position from the BAM pileup
/// - Error calculating genotype posteriors for reference genotype qual calculation (usually having to
///          do with accessing invalid genotypes)
pub fn select_random_potential_snvs(
    bam_file: &String,
    fasta_file: &String,
    interval: &Option<GenomicInterval>,
    sample_frequency: f64,
    min_coverage: u32,
    max_coverage: u32,
    min_mapq: u8,
) -> Result<VarList> {

    // the list of target (contig) names from the bam file
    let target_names = parse_target_names(&bam_file)?;

    let mut fasta = fasta::IndexedReader::from_file(&fasta_file)
        .chain_err(|| ErrorKind::IndexedFastaOpenError)?;

    let mut varlist: Vec<Var> = Vec::with_capacity(VARLIST_CAPACITY);

    // pileup over all covered sites
    let mut ref_seq: Vec<char> = vec![]; // this vector will be used to hold the reference sequence
    let mut prev_tid = 4294967295;
    let mut rng: StdRng = StdRng::from_seed(&[0]);

    // the 4 bases as String type
    let a_str = "A".to_string();
    let c_str = "C".to_string();
    let g_str = "G".to_string();
    let t_str = "T".to_string();

    // the strategy used for iterating over the BAM entries is as follows:
    // if a genomic region was specified, then ```get_interval_lst``` puts that genomic region into a vector (interval_lst)
    // containing only that one region. Then we iterate over the region list and seek every region,
    // which effectively means we just seek that single genomic interval.
    //
    // if a genomic region was not specified, then we want to iterate over the whole BAM file.
    // so get_interval_lst returns a list of genomic intervals that contain each whole chromosome
    // as described by the BAM header SQ lines. Then we iterate over the interval lst and seek each
    // interval separately, which effectively just iterates over all the BAM entries.
    //
    // the reason for this strange design (instead of either fetching a region beforehand or not and
    // then just iterating over all of ```bam_ix.pileup()```) is the following:
    // if an indexed reader is used, and fetch is never called, pileup() hangs.

    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bam_file, interval)?;
    let mut bam_ix =
        bam::IndexedReader::from_path(bam_file).chain_err(|| ErrorKind::IndexedBamOpenError)?;

    // interval_lst has either a single genomic interval (if --region was specified) or a list of
    // genomic intervals for each chromosome covering the entire genome
    for iv in interval_lst {
        bam_ix
            .fetch(iv.tid as u32, iv.start_pos as u32, iv.end_pos as u32 + 1)
            .chain_err(|| ErrorKind::IndexedBamFetchError)?;
        let bam_pileup = bam_ix.pileup();

        // this variable is used to avoid having a variant inside a previous variant's deletion.
        let mut next_valid_pos = 0;

        // iterate over base pileup
        for p in bam_pileup {

            if rng.next_f64() > sample_frequency {
                continue;
            }

            let pileup = p.chain_err(|| ErrorKind::IndexedBamPileupReadError)?;

            let tid: usize = pileup.tid() as usize;
            let chrom: String = target_names[tid].clone();

            // if we're on a different contig/chrom, we need to read in the sequence for that
            // contig/chrom from the FASTA into the ref_seq vector
            if tid != prev_tid {
                let mut ref_seq_u8: Vec<u8> = vec![];
                fasta
                    .read_all(&chrom, &mut ref_seq_u8)
                    .chain_err(|| ErrorKind::IndexedFastaReadError)?;
                ref_seq = dna_vec(&ref_seq_u8);
                next_valid_pos = 0;
            }

            prev_tid = tid;

            // this is specifically to avoid having a variant inside a previous variant's deletion.
            if pileup.pos() < next_valid_pos {
                continue;
            }

            // we want to save the sequence context (21 bp window around variant on reference)
            // this will be printed to the VCF later and may help diagnose variant calling
            // issues e.g. if the variant occurs inside a large homopolymer or etc.
            // get the position 10 bases to the left
            let l_window = if pileup.pos() >= 10 {
                pileup.pos() as usize - 10
            } else {
                0
            };
            // get the position 11 bases to the right
            let mut r_window = pileup.pos() as usize + 11;
            if r_window >= ref_seq.len() {
                r_window = ref_seq.len();
            }
            let sequence_context: String = (ref_seq[l_window..r_window]).iter().collect::<String>();

            let ref_base = ref_seq[pileup.pos() as usize];
            let ref_base_str = (ref_base).to_string();

            if ref_base_str.contains("N") {
                continue;
            }

            // the dna_vec conversion function should remove any non-ACGT bases
            assert!(
                ref_base_str == a_str
                    || ref_base_str == c_str
                    || ref_base_str == g_str
                    || ref_base_str == t_str
            );

            let mut depth: usize = 0;
            let pos: usize = pileup.pos() as usize;


            // pileup the bases for a single position and count number of each base
            for alignment in pileup.alignments() {
                let record = alignment.record();

                // check that the read doesn't fail any standard filters
                if record.is_secondary()
                    || record.is_quality_check_failed()
                    || record.is_duplicate()
                    || record.is_supplementary()
                    {
                        continue;
                    }

                if record.is_unmapped() || record.mapq() < min_mapq {
                    continue;
                }

                depth += 1; // depth counter does not consider unmapped low-mapq reads

            }

            if depth < min_coverage as usize {
                continue;
            }

            if depth > max_coverage as usize {
                continue;
            }

            let mut ref_allele = ref_base_str.clone();

            let var_allele = match ref_base {
                'A' => {rng.choose(&['C','G','T']).unwrap()},
                'C' => {rng.choose(&['A','G','T']).unwrap()},
                'G' => {rng.choose(&['A','C','T']).unwrap()},
                'T' => {rng.choose(&['A','C','G']).unwrap()},
                _ => {bail!("Invalid reference base while choosing random variant positions");}
            }.to_string();

            next_valid_pos = (pos + 1) as u32;

            // check if SNV meets our quality criteria for a potential SNV
            // if it does, make a new variant and add it to the list of potential SNVs.
            if !ref_allele.contains("N")
                && !var_allele.contains("N")
                {
                    let tid: usize = pileup.tid() as usize;
                    let new_var = Var {
                        ix: 0,
                        old_ix: None,
                        // these will be set automatically,
                        tid: tid,
                        chrom: target_names[tid].clone(),
                        pos0: pos,
                        alleles: vec![ref_allele.clone(), var_allele.clone()],
                        dp: depth,
                        allele_counts: vec![0, 0],
                        ambiguous_count: 0,
                        qual: 0.0,
                        filter: ".".to_string(),
                        genotype: Genotype(0, 0),
                        gq: 0.0,
                        unphased_genotype: Genotype(0, 0),
                        unphased_gq: 0.0,
                        genotype_post: GenotypeProbs::uniform(2),
                        phase_set: None,
                        mec: 0,
                        mec_frac_variant: 0.0, // mec fraction for this variant
                        mec_frac_block: 0.0,   // mec fraction for this haplotype block
                        mean_allele_qual: 0.0,
                        dp_any_mq: depth,
                        mq10_frac: 0.0,
                        mq20_frac: 0.0,
                        mq30_frac: 0.0,
                        mq40_frac: 0.0,
                        mq50_frac: 0.0,
                        sequence_context,
                        called: false,
                    };

                    // we don't want potential SNVs that are inside a deletion, for instance.
                    next_valid_pos = pileup.pos() + ref_allele.len() as u32;

                    varlist.push(new_var);
                }
        }
    }
    // return the vector of Vars as a VarList struct
    Ok(VarList::new(varlist)?)
}

// alignment: a rust-bio alignment object where x is a read consensus window, and y is the window from the reference
//
// l_ref: the 0-indexed position on the reference of the start of the reference window
/*
fn extract_variants_from_alignment(alignment: &Alignment,
                                   consensus: &Vec<u8>,
                                   ref_window: &Vec<u8>,
                                   l_ref: usize,
                                   tid: usize,
                                   chrom: String,
                                   depth: usize,
                                   boundary: usize) -> Result<Vec<Var>> {

    let mut ref_pos = alignment.ystart;
    let mut read_pos = alignment.xstart;

    let mut new_vars: Vec<Var> = vec![];
    for o in 0..alignment.operations.len() {
        //println!("read: {} ref: {}", ref_window[ref_pos] as char, consensus_all[read_pos] as char);
        let op = alignment.operations[o];
        let next_op = if o + 1 < alignment.operations.len() {
            alignment.operations[o + 1]
        } else {
            Match
        };

        if ref_pos >= boundary && ref_pos <= ref_window.len() - boundary {
            if op == Subst {
                let new_var = Var {
                    ix: 0,
                    old_ix: None,
                    // this will be set automatically,
                    tid: tid,
                    chrom: chrom.clone(),
                    pos0: l_ref + ref_pos,
                    alleles: vec![(ref_window[ref_pos] as char).to_string(), (consensus[read_pos] as char).to_string()],
                    dp: depth,
                    allele_counts: vec![0,0],
                    ambiguous_count: 0,
                    qual: 0.0,
                    filter: ".".to_string(),
                    genotype: Genotype(0,0), // this will be refined later
                    gq: 0.0,
                    genotype_post: GenotypeProbs::uniform(2),
                    phase_set: None,
                    mec: 0,                 // mec for variant
                    mec_frac_variant: 0.0,  // mec fraction for this variant
                    mec_frac_block: 0.0,    // mec fraction for this haplotype block
                    mean_allele_qual: 0.0,
                    dp_any_mq: depth,
                    mq10_frac: //TODO,
                    mq20_frac: //TODO,
                    called: false
                };
                new_vars.push(new_var);
            }

            if (op == Match || op == Subst) && next_op == Ins {
                let mut ref_allele = vec![ref_window[ref_pos]];
                let mut var_allele = vec![consensus[read_pos]];
                let mut offset = 1;
                for o2 in o + 1..alignment.operations.len() {
                    if alignment.operations[o2] == Ins {
                        var_allele.push(consensus[read_pos + offset]);
                    } else {
                        break;
                    }
                    offset += 1;
                }

                let new_var = Var {
                    ix: 0,
                    old_ix: None,
                    // this will be set automatically,
                    tid: tid,
                    chrom: chrom.clone(),
                    pos0: l_ref + ref_pos,
                    alleles: vec![u8_to_string(&ref_allele), u8_to_string(&var_allele)],
                    dp: depth,
                    allele_counts: vec![0,0],
                    ambiguous_count: 0,
                    qual: 0.0,
                    filter: ".".to_string(),
                    genotype: Genotype(0,0), // this will be refined later
                    gq: 0.0,
                    genotype_post: GenotypeProbs::uniform(2),
                    phase_set: None,
                    mec: 0,                 // mec for variant
                    mec_frac_variant: 0.0,  // mec fraction for this variant
                    mec_frac_block: 0.0,    // mec fraction for this haplotype block
                    mean_allele_qual: 0.0,
                    called: false
                };
                new_vars.push(new_var);
            }

            if (op == Match || op == Subst) && next_op == Del {
                let mut ref_allele = vec![ref_window[ref_pos]];
                let mut var_allele = vec![consensus[read_pos]];
                let mut offset = 1;

                for o2 in o + 1..alignment.operations.len() {
                    if alignment.operations[o2] == Del {
                        ref_allele.push(ref_window[ref_pos + offset]);
                    } else {
                        break;
                    }
                    offset += 1;
                }
                let new_var = Var {
                    ix: 0,
                    old_ix: None,
                    // this will be set automatically
                    tid: tid,
                    chrom: chrom.clone(),
                    pos0: l_ref + ref_pos,
                    alleles: vec![u8_to_string(&ref_allele), u8_to_string(&var_allele)],
                    dp: depth,
                    allele_counts: vec![0,0],
                    ambiguous_count: 0,
                    qual: 0.0,
                    filter: ".".to_string(),
                    genotype: Genotype(0,0), // this will be refined later
                    gq: 0.0,
                    genotype_post: GenotypeProbs::uniform(2),
                    phase_set: None,
                    mec: 0,                 // mec for variant
                    mec_frac_variant: 0.0,  // mec fraction for this variant
                    mec_frac_block: 0.0,    // mec fraction for this haplotype block
                    mean_allele_qual: 0.0,
                    called: false
                };
                new_vars.push(new_var);
            }
        }


        match op {
            Match => {
                ref_pos += 1;
                read_pos += 1;
            },
            Subst => {
                ref_pos += 1;
                read_pos += 1;
            },
            Ins => {
                read_pos += 1;
            },
            Del => {
                ref_pos += 1;
            },
            _ => {panic!("Unexpected alignment operation.");}
        }
    }

    new_vars
}
*/
/*
pub fn call_potential_variants_poa(bam_file: &String,
                                   fasta_file: &String,
                                   interval: &Option<GenomicInterval>,
                                   h1: &HashSet<String>,
                                   h2: &HashSet<String>,
                                   _max_coverage: Option<u32>,
                                   min_mapq: u8,
                                   _ln_align_params: LnAlignmentParameters)
                                   -> VarList {

    //let potential_snv_qual = LogProb::from(Prob(0.5));
    let target_names = parse_target_names(&bam_file);

    //let genotype_priors = estimate_genotype_priors();
    let mut fasta = fasta::IndexedReader::from_file(&fasta_file).chain_err(|| ErrorKind::IndexedFastaOpenError)?;

    let mut varlist: Vec<Var> = Vec::with_capacity(VARLIST_CAPACITY);

    // pileup over all covered sites
    let mut ref_seq: Vec<char> = vec![];
    let mut prev_tid = 4294967295;

    // there is a really weird bug going on here,
    // hence the duplicate file handles to the bam file.
    // if an indexed reader is used, and fetch is never called, pileup() hangs.
    // so we need to iterate over the fetched indexed pileup if there's a region,
    // or a totally separate pileup from the unindexed file if not.
    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bam_file, interval);

    let mut bam_ix = bam::IndexedReader::from_path(bam_file).chain_err(|| ErrorKind::IndexedBamOpenError)?;

    let alignment_type: i32 = 0;
    let match_score: i32 = 5;
    let mismatch_score: i32 = -4;
    let gap_score: i32 = -2;

    let d: usize = 50;

    for iv in interval_lst {
        bam_ix.fetch(iv.tid as u32, iv.start_pos as u32, iv.end_pos as u32 + 1).chain_err(|| ErrorKind::IndexedBamFetchError)?;
        let bam_pileup = bam_ix.pileup();

        for p in bam_pileup {
            let pileup = p.chain_err(||ErrorKind::IndexedBamPileupReadError)?;

            let tid: usize = pileup.tid() as usize;
            let chrom: String = target_names[tid].clone();

            if tid != prev_tid {
                let mut ref_seq_u8: Vec<u8> = vec![];
                fasta.read_all(&chrom, &mut ref_seq_u8).chain_err(|| ErrorKind::IndexedFastaReadError)?;
                ref_seq = dna_vec(&ref_seq_u8);
            }

            if pileup.pos() % d as u32 != 0 {
                prev_tid = tid;
                continue;
            }

            let pos_ref = pileup.pos() as usize;
            let l_ref: usize = pos_ref - d;
            let r_ref: usize = pos_ref + d;
            let mut ref_window: Vec<u8> = vec![];
            for i in l_ref..r_ref + 1 {
                ref_window.push(ref_seq[i]);
            }

            let mut ref_window_nullterm = ref_window.clone();
            ref_window_nullterm.push('\0' as u8);

            let mut all_read_seqs: Vec<Vec<u8>> = vec![];
            let mut h1_read_seqs: Vec<Vec<u8>> = vec![];
            let mut h2_read_seqs: Vec<Vec<u8>> = vec![];

            let mut all_seq_count: usize = 0;
            let mut h1_seq_count: usize = 0;
            let mut h2_seq_count: usize = 0;

            // pileup the bases for a single position and count number of each base
            for alignment in pileup.alignments() {
                let record = alignment.record();

                // may be faster to implement this as bitwise operation on raw flag in the future?
                if record.mapq() < min_mapq || record.is_unmapped() || record.is_secondary() ||
                    record.is_quality_check_failed() ||
                    record.is_duplicate() || record.is_supplementary() {
                    continue;
                }

                let pos_read = match alignment.qpos() {
                    Some(t) => t as usize,
                    None => { continue; }
                };

                let l_read = if pos_read >= d { pos_read - d } else { 0 };
                let mut r_read = pos_read + d;

                let len = alignment.record().seq().len();
                if r_read >= len {
                    r_read = len - 1;
                }

                let mut read_seq: Vec<u8> = vec![];
                for i in l_read..r_read {
                    let c = alignment.record().seq()[i].to_ascii_uppercase();
                    read_seq.push(c)
                }
                read_seq.push('\0' as u8);

                all_read_seqs.push(read_seq.clone());
                all_seq_count += 1;

                let read_id = u8_to_string(alignment.record().qname());

                if h1.contains(&read_id) {
                    h1_read_seqs.push(read_seq.clone());
                    h1_seq_count += 1;
                }

                if h2.contains(&read_id) {
                    h2_read_seqs.push(read_seq.clone());
                    h2_seq_count += 1;
                }
            }

            if VERBOSE {
                println!("{}:{}-{}", target_names[tid].clone(), l_ref, r_ref);
                println!("-------");

                println!("All read seqs:", );
                for (i, seq) in all_read_seqs.iter().enumerate() {
                    println!(">seq{}", i);
                    println!("{}", str::from_utf8(&seq.clone()).chain_err(|| "Error converting read sequences to valid UTF8")?);
                }
                println!("-------");
                println!("H1 read seqs:", );
                for (i, seq) in h1_read_seqs.iter().enumerate() {
                    println!(">seq{}", i);
                    println!("{}", str::from_utf8(&seq.clone()).chain_err(|| "Error converting read sequences to valid UTF8")?);
                }

                println!("-------");
                println!("H2 read seqs:", );
                for (i, seq) in h2_read_seqs.iter().enumerate() {
                    println!(">seq{}", i);
                    println!("{}", str::from_utf8(&seq.clone()).chain_err(|| "Error converting read sequences to valid UTF8")?);
                }
                println!("-------");
            }


            let score = |a: u8, b: u8| if a == b { 5i32 } else { -4i32 };
            let k = 6;  // kmer match length
            let w = 20;  // Window size for creating the band
            let mut aligner = Aligner::new(-8, -2, score, k, w);

            let consensus_max_len = 200;
            let min_reads = 5;

            let mut h1_vars: Option<Vec<Var>> = if h1_seq_count >= min_reads {
                let mut consensus_h1: Vec<u8> = vec![0u8; consensus_max_len];
                poa_multiple_sequence_alignment(&h1_read_seqs, ref_window_nullterm.clone(),
                                                1i32, &mut consensus_h1, alignment_type, //  (h1_seq_count/2) as i32
                                                match_score, mismatch_score, gap_score);
                let h1_alignment = aligner.local(&consensus_h1, &ref_window);
                //println!("{}\n", h1_alignment.pretty(&consensus_h1, &ref_window));
                Some(extract_variants_from_alignment(&h1_alignment,
                                                     &consensus_h1,
                                                     &ref_window,
                                                     l_ref,
                                                     tid,
                                                     target_names[tid].clone(),
                                                     all_seq_count,
                                                     25))
            } else {
                None
            };

            let mut h2_vars: Option<Vec<Var>> = if h2_seq_count >= min_reads {
                let mut consensus_h2: Vec<u8> = vec![0u8; consensus_max_len];
                poa_multiple_sequence_alignment(&h2_read_seqs, ref_window_nullterm.clone(),
                                                1i32, &mut consensus_h2, alignment_type, // (h1_seq_count/2) as i32
                                                match_score, mismatch_score, gap_score);
                let h2_alignment = aligner.local(&consensus_h2, &ref_window);
                //println!("{}\n", h2_alignment.pretty(&consensus_h2, &ref_window));
                Some(extract_variants_from_alignment(&h2_alignment,
                                                     &consensus_h2,
                                                     &ref_window,
                                                     l_ref,
                                                     tid,
                                                     target_names[tid].clone(),
                                                     all_seq_count,
                                                     25))
            } else {
                None
            };

            let mut all_vars: Option<Vec<Var>> = if h1_vars == None
                && h2_vars == None
                && all_seq_count >= min_reads {
                let mut consensus_all: Vec<u8> = vec![0u8; consensus_max_len];
                poa_multiple_sequence_alignment(&all_read_seqs, ref_window_nullterm.clone(),
                                                1i32, // (all_seq_count/2) as i32
                                                &mut consensus_all, alignment_type,
                                                match_score, mismatch_score, gap_score);
                let all_alignment = aligner.local(&consensus_all, &ref_window);
                //println!("{}\n", all_alignment.pretty(&consensus_all, &ref_window));
                Some(extract_variants_from_alignment(&all_alignment,
                                                     &consensus_all,
                                                     &ref_window,
                                                     l_ref,
                                                     tid,
                                                     target_names[tid].clone(),
                                                     all_seq_count,
                                                     25))
            } else {
                None
            };

            if VERBOSE {
                /*
                match &all_vars {
                    &Some(ref vars) => {
                        println!("ALL READS VARS:");
                        for var in vars {
                            println!("{}\t{}\t{}\t{}", var.chrom, var.pos0+1, var.ref_allele, var.var_allele);
                        }
                    },
                    &None => {}
                }*/

                match &h1_vars {
                    &Some(ref vars) => {
                        println!("H1 VARS:");
                        for var in vars {
                            println!("{}\t{}\t{}\t{}", var.chrom, var.pos0 + 1, var.alleles[0], var.alleles[1]);
                        }
                    },
                    &None => {}
                }

                match &h2_vars {
                    &Some(ref vars) => {
                        println!("H2 VARS:");
                        for var in vars {
                            println!("{}\t{}\t{}\t{}", var.chrom, var.pos0 + 1, var.alleles[0], var.alleles[1]);
                        }
                    },
                    &None => {}
                }
                println!("----------------------------------------------------------------------------");
            }

            match all_vars {
                Some(mut vars) => { varlist.append(&mut vars); }
                None => {}
            }

            match h1_vars {
                Some(mut vars) => { varlist.append(&mut vars); }
                None => {}
            }

            match h2_vars {
                Some(mut vars) => { varlist.append(&mut vars); }
                None => {}
            }

            prev_tid = tid;
        }
    }
    Ok(VarList::new(varlist))
}
*/
