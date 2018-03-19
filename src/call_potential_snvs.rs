extern crate rust_htslib;

use rust_htslib::bam;
use rust_htslib::prelude::*;
use bio::io::fasta;
use std::char;
use util::{LnAlignmentParameters, GenomicInterval, Var, VarList, parse_target_names, u8_to_string};
use rust_htslib::bam::pileup::Indel;
use std::collections::{HashMap,HashSet};
use bio::stats::{LogProb, Prob};
use call_genotypes::{estimate_genotype_priors, calculate_genotypes_without_haplotypes};
use spoa::poa_multiple_sequence_alignment;
use std::ascii::AsciiExt;

//use std::str;
use bio::alignment::Alignment;
use bio::alignment::pairwise::banded::*;
use bio::alignment::AlignmentOperation::*;
//use bio::alignment::sparse::hash_kmers;
//use bio::alignment::pairwise::{MIN_SCORE, Scoring};
//use bio::alignment::AlignmentOperation::*;

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
    /*
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
    */
    let mut bam_ix = bam::IndexedReader::from_path(bam_file).unwrap();

    match interval {
        &Some(ref iv) => {
            let iv_tid = bam_ix.header().tid(iv.chrom.as_bytes()).unwrap();
            bam_ix.fetch(iv_tid, iv.start_pos, iv.end_pos + 1).ok().expect("Error seeking BAM file while extracting fragments.");
        }
        &None => {},
    };

    let bam_pileup = bam_ix.pileup();

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
        let mut pileup_alleles: Vec<(String,String)> = vec![];
        let pos: usize = pileup.pos() as usize;

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

                match alignment.indel() {
                    Indel::None => {
                        // unwrapping a None value here
                        ref_allele =
                            (ref_seq[pos] as char).to_string().to_uppercase();

                        let base: char = alignment.record().seq()[alignment.qpos().unwrap()] as char;

                        var_allele = base.to_string();
                    },
                    Indel::Ins(l) => {

                        let start = alignment.qpos().unwrap();
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
                            None => "N".to_string()
                        };

                        var_allele = var_char.into_iter().collect::<String>();

                    },
                    Indel::Del(l) => {
                        let start = pos;
                        let end: usize = pos + l as usize + 1;

                        ref_allele = u8_to_string(&ref_seq[start..end]).to_uppercase();
                        var_allele = (ref_seq[pos] as char).to_string().to_uppercase();

                    },
                }

                pileup_alleles.push((ref_allele.clone(), var_allele.clone()));
                *counts.entry((ref_allele.clone(), var_allele.clone())).or_insert(0) += 1;
            }
        }

        match max_coverage {
            Some(cov) if depth > cov as usize => {continue;}
            _ => {}
        }

        let mut snv_max_count: u32 = 0;
        let mut snv_ref_allele = 'N'.to_string();
        let mut snv_var_allele = 'N'.to_string();

        //let mut indel_max_count: u32 = 0;
        //let mut indel_ref_allele = 'N'.to_string();
        //let mut indel_var_allele = 'N'.to_string();

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
            } /*else {
                // potential indel
                if count > indel_max_count as usize {
                    indel_max_count = count as u32;
                    indel_ref_allele = r.clone();
                    indel_var_allele = v.clone();
                }
            }*/
        }

        // vectors contain entries with (call, qual)
        // call is '0' or '1' (or potentially '2'...)
        // qual is a LogProb probability of miscall
        let mut snv_pileup_calls: Vec<(char, LogProb)> = vec![];
        //let mut indel_pileup_calls: Vec<(char, LogProb)> = vec![];

        for (ref_allele, var_allele) in pileup_alleles {

            if (ref_allele.clone(), var_allele.clone()) == (snv_ref_allele.clone(), snv_var_allele.clone()) {
                let qual = ln_align_params.emission_probs.not_equal;
                snv_pileup_calls.push(('1', qual));
            } /*else if (ref_allele.clone(), var_allele.clone()) == (indel_ref_allele.clone(), indel_var_allele.clone()) {
                let qual = if indel_ref_allele.len() > indel_var_allele.len() {
                    ln_align_params.emission_probs.deletion
                } else {
                    ln_align_params.emission_probs.insertion
                };

                indel_pileup_calls.push(('1', qual));
            } */

            if (ref_allele.clone(), var_allele.clone()) == (ref_allele.clone(), ref_allele.clone()){
                let qual = ln_align_params.emission_probs.not_equal;

                snv_pileup_calls.push(('0', qual));
                //indel_pileup_calls.push(('0', qual));
            }
        }

        // use a basic genotype likelihood calculation to call SNVs
        let snv_qual = if !snv_ref_allele.contains("N") && !snv_var_allele.contains("N") {
            let (_snv_post00, snv_post01, snv_post11) = calculate_genotypes_without_haplotypes(&snv_pileup_calls, &genotype_priors, &snv_ref_allele, &snv_var_allele);
            LogProb::ln_add_exp(snv_post01, snv_post11)
        } else {
            LogProb::ln_zero()
        };

        /*
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
        */

        let (ref_allele, var_allele, qual) = (snv_ref_allele, snv_var_allele, snv_qual);

        next_valid_pos = (pos+1) as u32;

        if qual > potential_snv_qual &&
            !ref_allele.contains("N") && !var_allele.contains("N") &&
            (ref_allele.clone(), var_allele.clone()) !=
                (ref_base_str.clone(), ref_base_str.clone()){

            let tid: usize = pileup.tid() as usize;
            let new_var = Var {
                ix: 0,
                // this will be set automatically
                chrom: target_names[tid].clone(),
                pos0: pos,
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
                phase_set: None,
                mec: 0,
                mec_frac: 0.0
            };

            // we don't want potential SNVs that are inside a deletion, for instance.
            next_valid_pos = pileup.pos() + ref_allele.len() as u32;

            varlist.push(new_var);
        }

        prev_tid = tid;
    }
    VarList::new(varlist)
}

// alignment: a rust-bio alignment object where x is a read consensus window, and y is the window from the reference
//
// l_ref: the 0-indexed position on the reference of the start of the reference window
fn extract_variants_from_alignment(alignment: &Alignment,
                                   consensus: &Vec<u8>,
                                   ref_window: &Vec<u8>,
                                   l_ref: usize,
                                   chrom: String,
                                   depth: usize,
                                   boundary: usize) -> Vec<Var>{

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
                    // this will be set automatically
                    chrom: chrom.clone(),
                    pos0: l_ref + ref_pos,
                    ref_allele: (ref_window[ref_pos] as char).to_string(),
                    var_allele: (consensus[read_pos] as char).to_string(),
                    dp: depth,
                    ra: 0,
                    aa: 0,
                    na: 0,
                    qual: 0.0,
                    filter: ".".to_string(),
                    genotype: "./.".to_string(),
                    gq: 0.0,
                    genotype_post: [LogProb::from(Prob(0.25)); 4],
                    phase_set: None,
                    mec: 0,
                    mec_frac: 0.0
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
                    // this will be set automatically
                    chrom: chrom.clone(),
                    pos0: l_ref + ref_pos,
                    ref_allele: u8_to_string(&ref_allele),
                    var_allele: u8_to_string(&var_allele),
                    dp: depth,
                    ra: 0,
                    aa: 0,
                    na: 0,
                    qual: 0.0,
                    filter: ".".to_string(),
                    genotype: "./.".to_string(),
                    gq: 0.0,
                    genotype_post: [LogProb::from(Prob(0.25)); 4],
                    phase_set: None,
                    mec: 0,
                    mec_frac: 0.0
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
                    // this will be set automatically
                    chrom: chrom.clone(),
                    pos0: l_ref + ref_pos,
                    ref_allele: u8_to_string(&ref_allele),
                    var_allele: u8_to_string(&var_allele),
                    dp: depth,
                    ra: 0,
                    aa: 0,
                    na: 0,
                    qual: 0.0,
                    filter: ".".to_string(),
                    genotype: "./.".to_string(),
                    gq: 0.0,
                    genotype_post: [LogProb::from(Prob(0.25)); 4],
                    phase_set: None,
                    mec: 0,
                    mec_frac: 0.0
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
    let mut fasta = fasta::IndexedReader::from_file(&fasta_file).unwrap();

    let varlist: Vec<Var> = Vec::with_capacity(VARLIST_CAPACITY);

    // pileup over all covered sites
    let mut ref_seq: Vec<u8> = vec![];
    let mut prev_tid = 4294967295;

    // there is a really weird bug going on here,
    // hence the duplicate file handles to the bam file.
    // if an indexed reader is used, and fetch is never called, pileup() hangs.
    // so we need to iterate over the fetched indexed pileup if there's a region,
    // or a totally separate pileup from the unindexed file if not.
    // TODO: try to reproduce as a minimal example and possibly raise issue on Rust-htslib repo

    let mut bam_ix = bam::IndexedReader::from_path(bam_file).unwrap();

    match interval {
        &Some(ref iv) => {
            let iv_tid = bam_ix.header().tid(iv.chrom.as_bytes()).unwrap();
            bam_ix.fetch(iv_tid, iv.start_pos, iv.end_pos + 1).ok().expect("Error seeking BAM file while extracting fragments.");
        }
        &None => {},
    };

    let bam_pileup = bam_ix.pileup();

    //let mut next_valid_pos = 0;

    let d:usize = 50;

    for p in bam_pileup {

        let pileup = p.unwrap();

        let tid: usize = pileup.tid() as usize;

        if tid != prev_tid {
            fasta.read_all(&target_names[tid], &mut ref_seq).expect("Failed to read fasta sequence record.");
            //next_valid_pos = 0;
        }

        if pileup.pos() % d as u32 != 0 {
            prev_tid = tid;
            continue;
        }

        let mut all_read_seqs: Vec<Vec<u8>> = vec![];
        let mut h1_read_seqs: Vec<Vec<u8>> = vec![];
        let mut h2_read_seqs: Vec<Vec<u8>> = vec![];

        let pos_ref = pileup.pos() as usize;
        let l_ref: usize = pos_ref - d;
        let r_ref: usize = pos_ref + d;
        let mut ref_window: Vec<u8> = vec![];
        for i in l_ref..r_ref + 1 {
            ref_window.push(AsciiExt::to_ascii_uppercase(&ref_seq[i]));
        }

        // pileup the bases for a single position and count number of each base
        for alignment in pileup.alignments() {
            let record = alignment.record();

            // may be faster to implement this as bitwise operation on raw flag in the future?
            if record.mapq() < min_mapq || record.is_unmapped() || record.is_secondary() ||
                record.is_quality_check_failed() ||
                record.is_duplicate() {
                continue;
            }

            let pos_read = match alignment.qpos() {
                Some(t) => t as usize,
                None => {continue;}
            };

            let l_read = if pos_read >= d {pos_read - d} else {0};
            let mut r_read = pos_read + d;

            let len = alignment.record().seq().len();
            if r_read >= len {
                r_read = len - 1;
            }

            let mut read_seq: Vec<u8> = vec![];
            for i in l_read..r_read {
                let c = AsciiExt::to_ascii_uppercase(&alignment.record().seq()[i]);
                read_seq.push(c)
            }
            read_seq.push('\0' as u8);
            //let read_seq = dna_vec(&alignment.record().seq()[l_read..r_read]);

            all_read_seqs.push(read_seq.clone());

            let read_id = u8_to_string(alignment.record().qname());

            if h1.contains(&read_id) {
                h1_read_seqs.push(read_seq.clone());
            }

            if h2.contains(&read_id) {
                h2_read_seqs.push(read_seq.clone());
            }
        }
        /*
        println!("{}:{}-{}",target_names[tid].clone(),l_ref, r_ref);
        println!("-------");

        println!("All read seqs:", );
        for (i,seq) in all_read_seqs.iter().enumerate() {
            println!(">seq{}",i);
            println!("{}", str::from_utf8(&seq.clone()).unwrap());
        }
        println!("-------");
        println!("H1 read seqs:", );
        for (i,seq) in h1_read_seqs.iter().enumerate() {
            println!(">seq{}",i);
            println!("{}", str::from_utf8(&seq.clone()).unwrap());
        }

        println!("-------");
        println!("H2 read seqs:", );
        for (i,seq) in h2_read_seqs.iter().enumerate() {
            println!(">seq{}",i);
            println!("{}", str::from_utf8(&seq.clone()).unwrap());
        }
        println!("-------");
        */

        let score = |a: u8, b: u8| if a == b {5i32} else {-4i32};
        let k = 6;  // kmer match length
        let w = 20;  // Window size for creating the band
        let mut aligner = Aligner::new(-8, -2, score, k, w);

        let consensus_max_len = 200;
        let min_reads = 5;

        let all_vars: Option<Vec<Var>> = if all_read_seqs.len() >= min_reads {
            let mut consensus_all: Vec<u8> = vec![0u8; consensus_max_len];
            poa_multiple_sequence_alignment(&all_read_seqs, &mut consensus_all);
            let all_alignment = aligner.local(&consensus_all, &ref_window);
            //println!("{}\n", all_alignment.pretty(&consensus_all, &ref_window));
            Some(extract_variants_from_alignment(&all_alignment,
                                            &consensus_all,
                                            &ref_window, l_ref,
                                            target_names[tid].clone(),
                                            all_read_seqs.len(),
                                             25))

        } else {
            None
        };

        let h1_vars: Option<Vec<Var>> = if h1_read_seqs.len() >= min_reads {
            let mut consensus_h1: Vec<u8> = vec![0u8; consensus_max_len];
            poa_multiple_sequence_alignment(&h1_read_seqs, &mut consensus_h1);
            let h1_alignment = aligner.local(&consensus_h1, &ref_window);
            //println!("{}\n", h1_alignment.pretty(&consensus_h1, &ref_window));
            Some(extract_variants_from_alignment(&h1_alignment,
                                              &consensus_h1,
                                              &ref_window, l_ref,
                                              target_names[tid].clone(),
                                              all_read_seqs.len(),
                                               25))

        } else {
            None
        };

        let h2_vars: Option<Vec<Var>> = if h2_read_seqs.len() >= min_reads {
            let mut consensus_h2: Vec<u8> = vec![0u8; consensus_max_len];
            poa_multiple_sequence_alignment(&h2_read_seqs, &mut consensus_h2);
            let h2_alignment = aligner.local(&consensus_h2, &ref_window);
            //println!("{}\n", h2_alignment.pretty(&consensus_h2, &ref_window));
            Some(extract_variants_from_alignment(&h2_alignment,
                                            &consensus_h2,
                                            &ref_window, l_ref,
                                            target_names[tid].clone(),
                                            all_read_seqs.len(),
                                            25))

        } else {
            None
        };


        match all_vars {
            Some(vars) => {
                for var in vars {
                    println!("{}\t{}\t{}\t{}", var.chrom, var.pos0+1, var.ref_allele, var.var_allele);
                }
            },
            None => {}
        }

        match h1_vars {
            Some(vars) => {
                for var in vars {
                    println!("{}\t{}\t{}\t{}", var.chrom, var.pos0+1, var.ref_allele, var.var_allele);
                }
            },
            None => {}
        }

        match h2_vars {
            Some(vars) => {
                for var in vars {
                    println!("{}\t{}\t{}\t{}", var.chrom, var.pos0+1, var.ref_allele, var.var_allele);
                }
            },
            None => {}
        }


        //println!("------------------------------------------------------");

        /*
        println!("{:?}",u8_to_string(&consensus_all));
        println!("{:?}",u8_to_string(&consensus_h1));
        println!("{:?}",u8_to_string(&consensus_h2));
        println!("------------------------------------------------------");
        */

        /*
        if qual > potential_snv_qual &&
            !ref_allele.contains("N") && !var_allele.contains("N") &&
            (ref_allele.clone(), var_allele.clone()) !=
                (ref_base_str.clone(), ref_base_str.clone()){

            let tid: usize = pileup.tid() as usize;
            let new_var = Var {
                ix: 0,
                // this will be set automatically
                chrom: target_names[tid].clone(),
                pos0: pos,
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
                phase_set: None,
                mec: 0,
                mec_frac: 0.0
            };

            varlist.push(new_var);

        }
        */
        prev_tid = tid;
    }
    VarList::new(varlist)
}
