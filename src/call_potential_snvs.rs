extern crate rust_htslib;

use rust_htslib::bam;
use rust_htslib::prelude::*;
use bio::io::fasta;
use std::char;
use util::{LnAlignmentParameters, GenomicInterval, Var, VarList, parse_target_names, u8_to_string};
use rust_htslib::bam::pileup::Indel;
use std::collections::HashMap;
use bio::stats::{LogProb, Prob};
use call_genotypes::{estimate_genotype_priors, calculate_genotypes_without_haplotypes};

static VARLIST_CAPACITY: usize = 1000000;

pub fn call_potential_snvs(bam_file: &String,
                           fasta_file: &String,
                           interval: &Option<GenomicInterval>,
                           max_coverage: Option<u32>,
                           min_mapq: u8,
                           ln_align_params: LnAlignmentParameters)
                           -> VarList {

    let potential_snv_qual = LogProb::from(Prob(0.5));
    let target_names = parse_target_names(&bam_file);

    let genotype_priors = estimate_genotype_priors();
    let mut fasta = fasta::IndexedReader::from_file(&fasta_file).unwrap();

    let mut varlist: Vec<Var> = Vec::with_capacity(VARLIST_CAPACITY);

    // pileup over all covered sites
    let mut ref_seq: Vec<u8> = vec![];
    let mut prev_tid = 4294967295;

    // there is a really weird bug going on here,
    // hence the duplicate file handles to the bam file.
    // if an indexed reader is used, and fetch is never called, pileup() hangs.
    // so we need to iterate over the fetched indexed pileup if there's a region,
    // or a totally separate pileup from the unindexed file if not.
    // TODO: try to reproduce as a minimal example and possibly raise issue on Rust-htslib repo
    let bam = bam::Reader::from_path(bam_file).unwrap();
    let mut bam_ix = bam::IndexedReader::from_path(bam_file).unwrap();
    let bam_pileup = match interval {
        &Some(ref iv) => {
            let iv_tid = bam_ix.header().tid(iv.chrom.as_bytes()).unwrap();
            bam_ix.fetch(iv_tid, iv.start_pos, iv.end_pos + 1).ok().expect("Error seeking BAM file while extracting fragments.");
            bam_ix.pileup()
        }
        &None => bam.pileup(),
    };

    let mut next_valid_pos = 0;

    for p in bam_pileup {

        let pileup = p.unwrap();

        let tid: usize = pileup.tid() as usize;

        if tid != prev_tid {
            fasta.read_all(&target_names[tid], &mut ref_seq).expect("Failed to read fasta sequence record.");
            next_valid_pos = 0;
        }

        // this is specifically to avoid having a variant inside a previous variant's deletion.
        if pileup.pos() < next_valid_pos {
            continue;
        }


        let ref_base_str = (ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();

        if ref_base_str.contains("N") {
            continue;
        }

        //let mut counts = [0; 5]; // A,C,G,T,N
        let mut counts: HashMap<(String, String), usize> = HashMap::new();
        // use a counter instead of pileup.depth() since that would include qc_fail bases, low mapq, etc.
        let mut depth: usize = 0;
        let mut pileup_alleles: Vec<(String,String, bool)> = vec![];

        // pileup the bases for a single position and count number of each base
        for alignment in pileup.alignments() {
            let record = alignment.record();

            // may be faster to implement this as bitwise operation on raw flag in the future?
            if record.mapq() < min_mapq || record.is_unmapped() || record.is_secondary() ||
                record.is_quality_check_failed() ||
                record.is_duplicate() {
                continue;
            }

            depth += 1;

            if !alignment.is_del() && !alignment.is_refskip() {

                let ref_allele;
                let var_allele;
                let mut homopolymer;

                match alignment.indel() {
                    Indel::None => {
                        // unwrapping a None value here
                        ref_allele =
                            (ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();

                        let ref_allele_prev =
                            (ref_seq[pileup.pos() as usize - 1] as char).to_string().to_uppercase();
                        let ref_allele_next =
                            (ref_seq[pileup.pos() as usize + 1] as char).to_string().to_uppercase();
                        let base: char = alignment.record().seq()[alignment.qpos().unwrap()] as char;
                        var_allele = base.to_string();
                        homopolymer = ref_allele == ref_allele_prev || ref_allele == ref_allele_next;
                    },
                    Indel::Ins(l) => {

                        let start = alignment.qpos().unwrap();
                        let end = start + l as usize + 1;
                        // don't want to convert whole seq to bytes...
                        let mut var_char: Vec<char> = vec![];
                        for i in start..end {
                            var_char.push(alignment.record().seq()[i] as char);
                        }
                        ref_allele =
                            (ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();
                        var_allele = var_char.into_iter().collect::<String>();


                        let ref_allele_prev =
                        (ref_seq[start - 1] as char).to_uppercase().nth(0).unwrap();
                        let ref_allele_next =
                            (ref_seq[end] as char).to_uppercase().nth(0).unwrap();
                        homopolymer = true;
                        for c in var_allele.chars() {
                            if ! (c == ref_allele_prev && c == ref_allele_next) {
                                homopolymer = false;
                            }
                        }

                    },
                    Indel::Del(l) => {
                        let start = pileup.pos() as usize;
                        let end = (pileup.pos() + l + 1) as usize;
                        ref_allele = u8_to_string(&ref_seq[start..end]).to_uppercase();
                        var_allele = (ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();

                        let ref_allele_prev =
                        (ref_seq[start - 1] as char).to_uppercase().nth(0).unwrap();
                        let ref_allele_next =
                            (ref_seq[end] as char).to_uppercase().nth(0).unwrap();
                        //(ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();

                        homopolymer = true;
                        for c in ref_allele.chars() {
                            if ! (c == ref_allele_prev && c == ref_allele_next) {
                                homopolymer = false;
                            }
                        }
                    },
                }

                pileup_alleles.push((ref_allele.clone(), var_allele.clone(), homopolymer));
                *counts.entry((ref_allele.clone(), var_allele.clone())).or_insert(0) += 1;
            }
        }

        //if depth < min_coverage as usize {
        //    continue
        //}

        match max_coverage {
            Some(cov) if depth > cov as usize => {continue;}
            _ => {}
        }

        //let mut var_allele = "N".to_string();
        let mut snv_max_count: u32 = 0;
        let mut snv_ref_allele = 'N'.to_string();
        let mut snv_var_allele = 'N'.to_string();
        //let mut var_allele = "N".to_string();
        let mut indel_max_count: u32 = 0;
        let mut indel_ref_allele = 'N'.to_string();
        let mut indel_var_allele = 'N'.to_string();

        // iterate over everything.

        for (&(ref r, ref v), &count) in &counts {

            if r.contains("N") || v.contains("N") {
                continue;
            }

            if r.len() == 1 && v.len() == 1 {
                // potential SNV
                if count > snv_max_count as usize && (r, v) != (&ref_base_str, &ref_base_str) {
                    snv_max_count = count as u32;
                    snv_ref_allele = r.clone();
                    snv_var_allele = v.clone();
                }
            } else {
                // potential indel
                if count > indel_max_count as usize {
                    indel_max_count = count as u32;
                    indel_ref_allele = r.clone();
                    indel_var_allele = v.clone();
                }
            }
        }

        // vectors contain entries with (call, qual)
        // call is '0' or '1' (or potentially '2'...)
        // qual is a LogProb probability of miscall
        let mut snv_pileup_calls: Vec<(char, LogProb)> = vec![];
        let mut indel_pileup_calls: Vec<(char, LogProb)> = vec![];

        for (ref_allele, var_allele, homopolymer) in pileup_alleles {

            if (ref_allele.clone(), var_allele.clone()) == (snv_ref_allele.clone(), snv_var_allele.clone()) {
                let qual = if homopolymer {ln_align_params.mismatch_from_match_homopolymer} else {ln_align_params.mismatch_from_match};
                snv_pileup_calls.push(('1', qual));
            } else if (ref_allele.clone(), var_allele.clone()) == (indel_ref_allele.clone(), indel_var_allele.clone()) {
                let qual = if indel_ref_allele.len() > indel_var_allele.len() {
                    if homopolymer {ln_align_params.deletion_from_match_homopolymer} else {ln_align_params.deletion_from_match}
                } else {
                    if homopolymer {ln_align_params.insertion_from_match_homopolymer} else {ln_align_params.insertion_from_match}
                };

                indel_pileup_calls.push(('1', qual));
            }

            if (ref_allele.clone(), var_allele.clone()) == (ref_allele.clone(), ref_allele.clone()){
                let qual = if homopolymer {ln_align_params.mismatch_from_match_homopolymer} else {ln_align_params.mismatch_from_match};

                snv_pileup_calls.push(('0', qual));
                indel_pileup_calls.push(('0', qual));
            }
        }

        // use a basic genotype likelihood calculation to call SNVs
        let snv_qual = if !snv_ref_allele.contains("N") && !snv_var_allele.contains("N") {
            let (_snv_post00, snv_post01, snv_post11) = calculate_genotypes_without_haplotypes(&snv_pileup_calls, &genotype_priors, &snv_ref_allele, &snv_var_allele);
            LogProb::ln_add_exp(snv_post01, snv_post11)
        } else {
            LogProb::ln_zero()
        };

        let indel_qual = if !indel_ref_allele.contains("N") && !indel_var_allele.contains("N") {
            let (_indel_post00, indel_post01, indel_post11) = calculate_genotypes_without_haplotypes(&indel_pileup_calls, &genotype_priors, &indel_ref_allele, &indel_var_allele);
            LogProb::ln_add_exp(indel_post01, indel_post11)
        } else {
            LogProb::ln_zero()
        };

        //println!("{} {}",*Prob::from(snv_qual),*Prob::from(indel_qual));

        let (ref_allele, var_allele, qual) = if snv_qual > indel_qual {
            (snv_ref_allele, snv_var_allele, snv_qual)
        } else {
            (indel_ref_allele, indel_var_allele, indel_qual)
        };

        next_valid_pos = pileup.pos()+1;

        if qual > potential_snv_qual &&
            !ref_allele.contains("N") && !var_allele.contains("N") &&
            (ref_allele.clone(), var_allele.clone()) !=
                (ref_base_str.clone(), ref_base_str.clone()){

            let tid: usize = pileup.tid() as usize;
            let new_var = Var {
                ix: 0,
                // this will be set automatically
                chrom: target_names[tid].clone(),
                pos0: pileup.pos() as usize,
                ref_allele: ref_allele.clone(),
                var_allele: var_allele.clone(),
                dp: depth,
                ra: 0,
                aa: 0,
                na: 0,
                qual: 0.0,
                filter: ".".to_string(),
                genotype: "./.".to_string(),
                gq: 0.0,
                genotype_post: [LogProb::from(Prob(0.25)); 4],
                phase_set: None
            };

            // we don't want potential SNVs that are inside a deletion, for instance.
            next_valid_pos = pileup.pos() + ref_allele.len() as u32;

            varlist.push(new_var);
        }

        prev_tid = tid;
    }
    VarList::new(varlist)
}
