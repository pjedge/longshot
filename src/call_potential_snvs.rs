extern crate rust_htslib;

use rust_htslib::bam;
use rust_htslib::prelude::*;
use bio::io::fasta;
use std::char;
use util::{GenomicInterval, PotentialVar, VarList, parse_target_names, u8_to_string};
use rust_htslib::bam::pileup::Indel;
use std::collections::HashMap;

static VARLIST_CAPACITY: usize = 1000000;
static INDEL_MIN_COUNT: u32 = 6;
static INDEL_MIN_FRAC: f32 = 0.3;

pub fn call_potential_snvs(bam_file: &String,
                           fasta_file: &String,
                           interval: &Option<GenomicInterval>,
                           min_alt_count: u32,
                           min_alt_frac: f32,
                           min_coverage: u32,
                           max_coverage: Option<u32>,
                           min_mapq: u8,
                           call_indels: bool)
                           -> VarList {
    let target_names = parse_target_names(&bam_file);

    let mut fasta = fasta::IndexedReader::from_file(&fasta_file).unwrap();

    let mut varlist: Vec<PotentialVar> = Vec::with_capacity(VARLIST_CAPACITY);

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

    for p in bam_pileup {
        let pileup = p.unwrap();

        let tid: usize = pileup.tid() as usize;

        match max_coverage {
            Some(cov) if pileup.depth() > cov => {
                continue;
            }
            _ => {}
        }

        if pileup.depth() < min_coverage {
            continue;
        }

        if tid != prev_tid {
            fasta.read_all(&target_names[tid], &mut ref_seq).expect("Failed to read fasta sequence record.");
        }

        //let mut counts = [0; 5]; // A,C,G,T,N
        let mut counts: HashMap<(String, String), usize> = HashMap::new();

        // pileup the bases for a single position and count number of each base
        for alignment in pileup.alignments() {
            let record = alignment.record();

            // may be faster to implement this as bitwise operation on raw flag in the future?
            if record.mapq() < min_mapq || record.is_unmapped() || record.is_secondary() ||
                record.is_quality_check_failed() ||
                record.is_duplicate() || record.is_supplementary() {
                continue;
            }
            if !alignment.is_del() && !alignment.is_refskip() {
                let base: char = alignment.record().seq()[alignment.qpos().unwrap()] as char;
                //let b = match base {
                //    'A' | 'a' => 0,
                //    'C' | 'c' => 1,
                //    'G' | 'g' => 2,
                //    'T' | 't' => 3,
                //    'N' | 'n' => 4,
                //    _ => panic!("Invalid base read from BAM file."),
                //};

                let ref_allele;
                let var_allele;

                if call_indels {
                    match alignment.indel() {
                        Indel::None => {
                            ref_allele =
                                (ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();
                            var_allele = base.to_string();
                        }
                        Indel::Ins(l) => {
                            ref_allele =
                                (ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();

                            let start = alignment.qpos().unwrap();
                            let end = start + l as usize + 1;
                            // don't want to convert whole seq to bytes...
                            let mut var_char: Vec<char> = vec![];
                            for i in start..end {
                                var_char.push(alignment.record().seq()[i] as char);
                            }
                            var_allele = var_char.into_iter().collect::<String>();
                        }
                        Indel::Del(l) => {
                            let start = pileup.pos() as usize;
                            let end = (pileup.pos() + l + 1) as usize;
                            ref_allele = u8_to_string(&ref_seq[start..end]).to_uppercase();
                            var_allele = base.to_string();
                            //(ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();
                        }
                    }
                } else {
                    ref_allele =
                        (ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();
                    var_allele = base.to_string();
                }

                //counts[(ref)] += 1;
                *counts.entry((ref_allele, var_allele)).or_insert(0) += 1;
            }
        }

        //let mut var_allele = "N".to_string();
        let mut max_count: u32 = 0;
        let mut ref_allele = 'N'.to_string();
        let mut var_allele = 'N'.to_string();

        let ref_base_str = (ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();

        // iterate over everything.

        for (&(ref r, ref v), &count) in &counts {
            //eprintln!("{}, {}: {}", r, v, count);
            let frac = (count as f32) / (pileup.depth() as f32);
            if count > max_count as usize && (r, v) != (&ref_base_str, &ref_base_str) &&
                (r.len() == 1 && v.len() == 1 ||
                    count >= INDEL_MIN_COUNT as usize && frac >= INDEL_MIN_FRAC) {
                max_count = count as u32;
                ref_allele = r.clone();
                var_allele = v.clone();
            }
        }

        //eprintln!("--------------");

        let alt_frac: f32 = (max_count as f32) / (pileup.depth() as f32); //(max_count as f32) / (base_cov as f32);

        if !ref_allele.contains("N") && !var_allele.contains("N") &&
            (ref_allele.clone(), var_allele.clone()) !=
                (ref_base_str.clone(), ref_base_str.clone()) &&
            (ref_allele.len() == 1 && var_allele.len() == 1 && max_count >= min_alt_count &&
                alt_frac >= min_alt_frac ||
                max_count >= INDEL_MIN_COUNT && alt_frac >= INDEL_MIN_FRAC) {
            let tid: usize = pileup.tid() as usize;
            let new_var = PotentialVar {
                ix: 0,
                // this will be set automatically
                chrom: target_names[tid].clone(),
                pos0: pileup.pos() as usize,
                ref_allele: ref_allele,
                var_allele: var_allele,
                dp: pileup.depth() as usize
            };

            varlist.push(new_var);
        }

        prev_tid = tid;
    }
    VarList::new(varlist)
}
