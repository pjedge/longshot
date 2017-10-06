use std::collections::HashSet;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam::record::Cigar;
use std::error::Error;
use util::*;
use bio::stats::{LogProb, Prob, PHREDProb};
use bio::io::fasta;
use realignment;

static VERBOSE: bool = false;

// returns true if kmers in seq are unique, else false
fn check_kmers_unique(seq: &Vec<char>, k: usize) -> (bool) {
    //seq_str = seq.iter().cloned().collect();
    let mut kmers: HashSet<&[char]> = HashSet::new();

    for i in 0..(seq.len() - k + 1) {
        let kmer = &seq[i..(i + k)]; // (&seq[i..(i+k)]).iter().cloned().collect();
        //kmer.clone_from_slice(&seq[i..(i + k)]);
        if kmers.contains(kmer) {
            return false;
        }

        kmers.insert(kmer);
    }

    true
}

// an extension of the rust-htslib cigar representation
// has the cigar operation as well as the position on the reference of the operation start,
// and the position on the read of the operation start
pub struct CigarPos {
    cig: Cigar,
    ref_pos: u32,
    read_pos: u32,
}

//************************************************************************************************
// BEGINNING OF RUST-HTSLIB BASED CODE *****************************************************************
//************************************************************************************************

// This function (create_augmented_cigarlist) is a modification of Rust-Htslib function rust_htslib::bam::record::CigarStringView::read_pos
// https://github.com/rust-bio/rust-htslib/blob/master/src/bam/record.rs
// Rust-Htslib license is copied here as per its terms:
//The MIT License (MIT)
//
//Copyright (c) 2016 Johannes Köster, the Rust-Htslib team.

//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

quick_error! {
    #[derive(Debug)]
    pub enum CigarOrAnchorError {
        UnsupportedOperation(msg: String) {
            description("Unsupported CIGAR operation")
            display(x) -> ("{}: {}", x.description(), msg)
        }
        UnexpectedOperation(msg: String) {
            description("CIGAR operation not allowed at this point")
            display(x) -> ("{}: {}", x.description(), msg)
        }
        AnchorRangeOutsideRead(msg: String) {
            description("Attempted to find sequence anchors for a range completely outside of the sequence.")
            display(x) -> ("{}: {}", x.description(), msg)
        }
    }
}

/// Creates a vector of CigarPos, which contain each cigar operation and the reference position
/// and read position at the start of the operation
///
/// # Arguments
///
/// * `cigar_string_view` - cigar string to create cigar position list for
///
pub fn create_augmented_cigarlist(refpos: u32,
                                  cigar_string_view: &CigarStringView)
                                  -> Result<Vec<CigarPos>, CigarOrAnchorError> {
    let mut rpos = refpos;
    let mut qpos = 0u32; // position within read
    let mut j = 0; // index into cigar operation vector

    let mut cigar_list: Vec<CigarPos> = Vec::with_capacity(cigar_string_view.len());
    // find first cigar operation referring to qpos = 0 (and thus bases in record.seq()),
    // because all augmentations of qpos and rpos before that are invalid
    for (i, c) in cigar_string_view.iter().enumerate() {
        match c {
            &Cigar::Match(_) |
            &Cigar::Diff(_) |
            &Cigar::Equal(_) |
            // this is unexpected, but bwa + GATK indel realignment can produce insertions
            // before matching positions
            &Cigar::Ins(_) => {
                j = i;
                break;
            }
            &Cigar::SoftClip(_) => {
                j = i;
                break;
            }
            &Cigar::Del(_) => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "'deletion' (D) found before any operation describing read sequence".to_owned()
                ));
            }
            &Cigar::Back(_) => {
                return Err(CigarOrAnchorError::UnsupportedOperation(
                    "'back' (B) operation is deprecated according to htslib/bam_plcmd.c and is not in SAMv1 spec".to_owned()
                ));
            }
            &Cigar::RefSkip(_) => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "'reference skip' (N) found before any operation describing read sequence".to_owned()
                ));
            }
            &Cigar::HardClip(_) if i > 0 && i < cigar_string_view.len() - 1 => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                ));
            }
            // if we have reached the end of the CigarString with only pads and hard clips, we have no read position matching the variant
            &Cigar::Pad(_) | &Cigar::HardClip(_) if i == cigar_string_view.len() - 1 => return Ok(cigar_list),
            // skip leading HardClips and Pads, as they consume neither read sequence nor reference sequence
            &Cigar::Pad(_) | &Cigar::HardClip(_) => ()
        }
    }

    while j < cigar_string_view.len() {
        // append CigarPos to list (contains cigar op, ref position, read position)
        match &cigar_string_view[j] {
            &Cigar::Match(_) |
            &Cigar::Diff(_) |
            &Cigar::Equal(_) |
            &Cigar::Ins(_) |
            &Cigar::Del(_) |
            &Cigar::RefSkip(_) => {
                // push the current cigar operation, ref position, and query position onto the list
                cigar_list.push(CigarPos {
                    cig: cigar_string_view[j].clone(),
                    ref_pos: rpos,
                    read_pos: qpos,
                });
            }
            &Cigar::SoftClip(_) |
            &Cigar::Pad(_) |
            &Cigar::HardClip(_) |
            &Cigar::Back(_) => {}
        }

        // move the reference and query positions forward
        match &cigar_string_view[j] {
            &Cigar::Match(l) | &Cigar::Diff(l) | &Cigar::Equal(l) => {
                rpos += l;
                qpos += l;
                j += 1;
            }
            &Cigar::SoftClip(l) => {
                qpos += l;
                j += 1;
            }
            &Cigar::Ins(l) => {
                qpos += l;
                j += 1;
            }
            &Cigar::RefSkip(l) |
            &Cigar::Del(l) => {
                rpos += l;
                j += 1;
            }
            &Cigar::Pad(_) => {
                j += 1;
            }
            &Cigar::HardClip(_) if j < cigar_string_view.len() - 1 => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                ));
            }
            &Cigar::HardClip(_) => {}
            &Cigar::Back(_) => {
                return Err(CigarOrAnchorError::UnsupportedOperation(
                    "'back' (B) operation is deprecated according to htslib/bam_plcmd.c and is not in SAMv1 spec".to_owned()
                ));
            }
        }
    }

    Ok(cigar_list)
}

pub struct AnchorPositions {
    pub left_anchor_ref: u32,
    pub right_anchor_ref: u32,
    pub left_anchor_read: u32,
    pub right_anchor_read: u32,
}

pub fn find_anchors(bam_record: &Record,
                    var_interval: GenomicInterval,
                    anchor_length: u32,
                    ref_seq: &Vec<char>,
                    target_names: &Vec<String>,
                    extract_params: ExtractFragmentParameters)
                    -> Result<Option<AnchorPositions>, CigarOrAnchorError> {
    let anchor_k = extract_params.anchor_k;
    let min_window_length = extract_params.min_window_length;
    let max_window_length = extract_params.max_window_length;

    if var_interval.chrom != target_names[bam_record.tid() as usize] ||
        (var_interval.start_pos as i32) >=
            bam_record.cigar().end_pos().expect("Error while accessing CIGAR end position") ||
        (var_interval.end_pos as i32) < bam_record.pos() {
        return Err(CigarOrAnchorError::AnchorRangeOutsideRead(
            "Attempted to find sequence anchors for a range completely outside of the sequence.".to_owned()
        ));
    }

    let bam_cig: CigarStringView = bam_record.cigar();
    let rpos: u32 = bam_record.pos() as u32;
    let cigarpos_list: Vec<CigarPos> =
        create_augmented_cigarlist(rpos, &bam_cig).expect("Error creating augmented cigarlist.");

    let mut left_anchor_ref: u32 = 0;
    let mut right_anchor_ref: u32 = 0;
    let mut left_anchor_read: u32 = 0;
    let mut right_anchor_read: u32 = 0;

    let mut left_ix = 0;
    let mut right_ix = 0;

    let contains_pos = |pos: u32, cigar_op_start: u32, cigar_op_length: u32| {
        cigar_op_start <= pos && cigar_op_start + cigar_op_length > pos
    };

    // find left_ix and right_ix.
    // these are the indexes into the CIGAR of the CIGAR operation holding the
    // left side of the range and right side of the range
    for (i, cigarpos) in cigarpos_list.iter().enumerate() {
        match cigarpos.cig {
            Cigar::Match(l) |
            Cigar::Diff(l) |
            Cigar::Equal(l) |
            Cigar::Ins(l) |
            Cigar::Del(l) |
            Cigar::RefSkip(l) => {
                if contains_pos(var_interval.start_pos, cigarpos.ref_pos, l) {
                    left_ix = i;
                }
                if contains_pos(var_interval.end_pos, cigarpos.ref_pos, l) {
                    right_ix = i;
                    break;
                }
            }
            Cigar::Pad(_) |
            Cigar::Back(_) |
            Cigar::SoftClip(_) |
            Cigar::HardClip(_) => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "CIGAR operation found in cigarpos_list that should have been removed already.".to_owned()
                ));
            }
        }
    }

    // step backwards to find left anchor
    let mut match_len_left = 0;
    let mut seen_indel_left = false;
    let mut found_anchor_left = false;
    for i in (0..left_ix + 1).rev() {
        match cigarpos_list[i].cig {
            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                match_len_left += l;
                if seen_indel_left && match_len_left > anchor_length {
                    // we have found a sufficiently long match and it's NOT the cigar op containing var_interval.start_pos
                    left_anchor_ref = cigarpos_list[i].ref_pos + match_len_left - anchor_length;
                    left_anchor_read = cigarpos_list[i].read_pos + match_len_left - anchor_length;
                    let l: usize = left_anchor_ref as usize - anchor_k;
                    let r: usize = (left_anchor_ref + anchor_length) as usize + anchor_k;

                    if l > 0 && r < ref_seq.len() &&
                        check_kmers_unique(&ref_seq[l..r + 1].to_vec(), anchor_k) {
                        found_anchor_left = true;
                        break;
                    }
                } else if !seen_indel_left &&
                    cigarpos_list[i].ref_pos < var_interval.start_pos - anchor_length {
                    // we have found a sufficiently long match but it's the cigar ops surrounding var_interval.start_pos
                    // the var_interval.start_pos we are trying to anchor is just inside one huge match
                    left_anchor_ref = var_interval.start_pos - anchor_length;
                    left_anchor_read = cigarpos_list[i].read_pos +
                        (var_interval.start_pos - anchor_length -
                            cigarpos_list[i].ref_pos);
                    let l: usize = left_anchor_ref as usize - anchor_k;
                    let r: usize = (left_anchor_ref + anchor_length) as usize + anchor_k;

                    if l > 0 && r < ref_seq.len() &&
                        check_kmers_unique(&ref_seq[l..r + 1].to_vec(), anchor_k) {
                        found_anchor_left = true;
                        break;
                    }
                }
            }
            Cigar::Ins(_) |
            Cigar::Del(_) |
            Cigar::RefSkip(_) => {
                match_len_left = 0;
                seen_indel_left = true;
            }
            Cigar::Pad(_) |
            Cigar::Back(_) |
            Cigar::SoftClip(_) |
            Cigar::HardClip(_) => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "CIGAR operation found in cigarpos_list that should have been removed already.".to_owned()
                ));
            }
        }
    }

    if !found_anchor_left {
        return Ok(None); // failed to find a left anchor
    }

    // step forwards to find right anchor
    let mut match_len_right = 0;
    let mut seen_indel_right = false;
    let mut found_anchor_right = false;
    for i in right_ix..cigarpos_list.len() {
        match cigarpos_list[i].cig {
            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                match_len_right += l;
                if seen_indel_right && match_len_right > anchor_length {
                    // we have found a sufficiently long match and it's NOT the cigar op containing var_interval.end_pos
                    right_anchor_ref = cigarpos_list[i].ref_pos + l - match_len_right +
                        anchor_length;
                    right_anchor_read = cigarpos_list[i].read_pos + l - match_len_right +
                        anchor_length;
                    let l: usize = (right_anchor_ref - anchor_length) as usize - anchor_k;
                    let r: usize = right_anchor_ref as usize + anchor_k;

                    if l > 0 && r < ref_seq.len() &&
                        check_kmers_unique(&ref_seq[l..r + 1].to_vec(), anchor_k) {
                        found_anchor_right = true;
                        break;
                    }
                } else if !seen_indel_right &&
                    cigarpos_list[i].ref_pos + l > var_interval.end_pos + anchor_length {
                    // we have found a sufficiently long match but it's the cigar ops surrounding var_interval.end_pos
                    // the var_interval.end_pos we are trying to anchor is just inside one huge match
                    right_anchor_ref = var_interval.end_pos + anchor_length;
                    right_anchor_read = cigarpos_list[i].read_pos +
                        (var_interval.end_pos + anchor_length -
                            cigarpos_list[i].ref_pos);
                    let l: usize = (right_anchor_ref - anchor_length) as usize - anchor_k;
                    let r: usize = right_anchor_ref as usize + anchor_k;

                    if l > 0 && r < ref_seq.len() &&
                        check_kmers_unique(&ref_seq[l..r + 1].to_vec(), anchor_k) {
                        found_anchor_right = true;
                        break;
                    }
                }
            }
            Cigar::Ins(_) |
            Cigar::Del(_) |
            Cigar::RefSkip(_) => {
                match_len_right = 0;
                seen_indel_right = true;
            }
            Cigar::Pad(_) |
            Cigar::Back(_) |
            Cigar::SoftClip(_) |
            Cigar::HardClip(_) => {
                return Err(CigarOrAnchorError::UnexpectedOperation(
                    "CIGAR operation found in cigarpos_list that should have been removed already.".to_owned()
                ));
            }
        }
    }

    if !found_anchor_right {
        return Ok(None); // failed to find a right anchor
    }

    // return none if read window or ref window is larger or smaller than the allowed lengths
    let ref_window_len = (right_anchor_ref - left_anchor_ref) as usize;
    let read_window_len = (right_anchor_read - left_anchor_read) as usize;
    if ref_window_len < min_window_length || read_window_len < min_window_length ||
        ref_window_len > max_window_length || read_window_len > max_window_length {
        return Ok(None);
    }

    // return none if any of the anchors are out of bounds
    if right_anchor_ref as usize >= ref_seq.len() ||
        right_anchor_read as usize >= bam_record.seq().len() {
        return Ok(None);
    }

    // return anchor pos
    Ok(Some(AnchorPositions {
        left_anchor_ref: left_anchor_ref,
        right_anchor_ref: right_anchor_ref,
        left_anchor_read: left_anchor_read,
        right_anchor_read: right_anchor_read,
    }))
}

fn extract_var_cluster(read_seq: &Vec<char>,
                       ref_seq: &Vec<char>,
                       var_cluster: Vec<Var>,
                       anchors: AnchorPositions,
                       extract_params: ExtractFragmentParameters,
                       align_params: AlignmentParameters)
                       -> Vec<FragCall> {
    let mut calls: Vec<FragCall> = vec![];

    //let ref_window = ref_seq[(anchors.left_anchor_ref as usize)..
    //(anchors.right_anchor_ref as usize) + 1]
    //        .to_vec();
    let window_capacity = (anchors.right_anchor_ref - anchors.left_anchor_ref + 10) as usize;

    let read_window: Vec<char> = read_seq[(anchors.left_anchor_read as usize)..
        (anchors.right_anchor_read as usize) + 1]
        .to_vec();

    let mut max_score: LogProb = LogProb::ln_zero();
    let mut max_hap: usize = 0;
    let n_vars: usize = var_cluster.len() as usize; // number of variants in cluster
    let n_haps: usize = 2usize.pow(n_vars as u32); // number of possible haplotypes for cluster

    let in_hap = |var: usize, hap: usize| (hap & 2usize.pow(var as u32)) > 0;

    // allele_scores0[i] contains the Log sum of the probabilities of all short haplotypes
    // that had a '0' at the ith variant of the cluster.
    // similarly for allele_scores1, for '1' variants.
    let mut allele_scores0: Vec<LogProb> = vec![LogProb::ln_zero(); n_vars];
    let mut allele_scores1: Vec<LogProb> = vec![LogProb::ln_zero(); n_vars];
    let mut score_total: LogProb = LogProb::ln_zero();

    if VERBOSE {
        for var in var_cluster.clone() {
            println!("{} {} {} {}",
                     var.chrom,
                     var.pos0,
                     var.ref_allele,
                     var.var_allele);
        }
        let read_seq_str: String = read_window.clone().into_iter().collect();
        println!("read: {}", read_seq_str);
    }


    for hap in 0..n_haps {
        let mut hap_window: Vec<char> = Vec::with_capacity(window_capacity);
        let mut i: usize = anchors.left_anchor_ref as usize;
        for var in 0..n_vars {
            while i < var_cluster[var].pos0 {
                hap_window.push(ref_seq[i]);
                i += 1;
            }
            if in_hap(var, hap) {
                for c in var_cluster[var].var_allele.chars() {
                    hap_window.push(c);
                }
            } else {
                for c in var_cluster[var].ref_allele.chars() {
                    hap_window.push(c);
                }
            }
            i += var_cluster[var].ref_allele.len();
        }

        while i <= anchors.right_anchor_ref as usize {
            hap_window.push(ref_seq[i]);
            i += 1;
        }

        // we now want to score hap_window
        let score: LogProb = match extract_params.alignment_type{
            AlignmentType::NumericallyStableAllAlignment => {
                realignment::sum_all_alignments_numerically_stable(&read_window,
                                                                   &hap_window,
                                                                   align_params.ln(),
                                                                   extract_params.band_width)
            }
            AlignmentType::FastAllAlignment => {
                realignment::sum_all_alignments(&read_window,
                                                &hap_window,
                                                align_params,
                                                extract_params.band_width)
            }
            AlignmentType::MaxAlignment => {
                realignment::max_alignment(&read_window,
                                           &hap_window,
                                           align_params.ln(),
                                           extract_params.band_width)
            }
        };

        for var in 0..n_vars {
            if in_hap(var, hap) {
                allele_scores1[var] = LogProb::ln_add_exp(allele_scores1[var], score);
            } else {
                allele_scores0[var] = LogProb::ln_add_exp(allele_scores0[var], score);
            }
        }
        if VERBOSE {
            let hap_seq_str: String = hap_window.into_iter().collect();
            println!("hap:{} {} {}", hap, hap_seq_str, *PHREDProb::from(score));
        }
        // add current alignment score to the total score sum
        score_total = LogProb::ln_add_exp(score_total, score);

        if score > max_score {
            max_score = score;
            max_hap = hap;
        }
    }
    if VERBOSE {
        println!("--------------------------------------");
    }

    for v in 0..n_vars {
        if in_hap(v, max_hap) {
            // the best haplotype has a '1' variant allele at this position

            // quality score is the probability call is wrong, in other words the ratio of the
            // sum of haplotype scores that had a '0' here over the total sum of scores
            let qual = allele_scores0[v] - score_total;

            calls.push(FragCall {
                var_ix: var_cluster[v].ix,
                allele: '1',
                qual: qual,
                p_hap1: LogProb::from(Prob(0.5)),
                // haplotype unknown at this time
                p_hap2: LogProb::from(Prob(0.5)),
                // haplotype unknown at this time
            });
        } else {
            // the best haplotype has a '0' ref allele at this position

            // quality score is the probability call is wrong, in other words the ratio of the
            // sum of haplotype scores that had a '1' here over the total sum of scores
            let qual = allele_scores1[v] - score_total;

            calls.push(FragCall {
                var_ix: var_cluster[v].ix,
                allele: '0',
                qual: qual,
                p_hap1: LogProb::from(Prob(0.5)),
                // haplotype unknown at this time
                p_hap2: LogProb::from(Prob(0.5)),
                // haplotype unknown at this time
            });
        }
    }

    calls
}

pub fn extract_fragment(bam_record: &Record,
                        vars: Vec<Var>,
                        ref_seq: &Vec<char>,
                        target_names: &Vec<String>,
                        extract_params: ExtractFragmentParameters,
                        align_params: AlignmentParameters)
                        -> Option<Fragment> {

    if bam_record.is_quality_check_failed() || bam_record.is_duplicate() ||
        bam_record.is_secondary() || bam_record.is_unmapped() || bam_record.mapq() < extract_params.min_mapq
        {
            return None;
        }

    // TODO assert that every single variant in vars is on the same chromosome
    let id: String = u8_to_string(bam_record.qname());
    let mut fragment = Fragment {
        id: id,
        calls: vec![],
    };
    let read_seq: Vec<char> = dna_vec(&bam_record.seq().as_bytes());
    let mut cluster_lst: Vec<Vec<Var>> = vec![];

    {
        // populate cluster_lst with variant clusters
        let mut var_cluster: Vec<Var> = vec![];

        // generate clusters of SNVs that should be considered
        for var in vars {
            if var_cluster.len() == 0 ||
                (var.pos0 - var_cluster[var_cluster.len() - 1].pos0 <
                    extract_params.short_hap_snv_distance &&
                    var_cluster.len() <= extract_params.short_hap_max_snvs) {
                var_cluster.push(var);
            } else {
                cluster_lst.push(var_cluster.clone());
                var_cluster.clear();
                var_cluster.push(var);
            }
        }
        if var_cluster.len() > 0 {
            cluster_lst.push(var_cluster.clone());
            var_cluster.clear();
        }
    }

    for var_cluster in cluster_lst {
        let var_interval = GenomicInterval {
            chrom: target_names[bam_record.tid() as usize].clone(),
            start_pos: var_cluster[0].pos0 as u32,
            end_pos: var_cluster[var_cluster.len() - 1].pos0 as u32,
        };

        match find_anchors(bam_record,
                           var_interval,
                           10,
                           &ref_seq,
                           target_names,
                           extract_params)
            .expect("CIGAR or Anchor Error while finding anchor sequences.") {
            Some(anchors) => {
                for call in extract_var_cluster(&read_seq,
                                                ref_seq,
                                                var_cluster,
                                                anchors,
                                                extract_params,
                                                align_params) {
                    fragment.calls.push(call);
                }
            }
            None => {
                continue;
            }
        }
    }

    if fragment.calls.len() > 0 {
        Some(fragment)
    } else {
        None
    }
}

pub fn extract_fragments(bamfile_name: &String,
                         fastafile_name: &String,
                         varlist: &VarList,
                         interval: &Option<GenomicInterval>,
                         extract_params: ExtractFragmentParameters,
                         align_params: AlignmentParameters)
                         -> Vec<Fragment> {
    let t_names = parse_target_names(&bamfile_name);

    let mut prev_tid = 4294967295; // huge value so that tid != prev_tid on first iter
    let mut fasta = fasta::IndexedReader::from_file(fastafile_name).unwrap();
    let mut ref_seq: Vec<char> = vec![];

    let mut flist: Vec<Fragment> = vec![];

    // TODO: this uses a lot of duplicate code, need to figure out a better solution.
    match interval {
        &Some(ref iv) => {
            let mut bam = bam::IndexedReader::from_path(bamfile_name).unwrap();
            let iv_tid = bam.header().tid(iv.chrom.as_bytes()).unwrap();
            bam.fetch(iv_tid, iv.start_pos, iv.end_pos + 1).ok().expect("Error seeking BAM file while extracting fragments.");
            for r in bam.records() {
                let record = r.unwrap();

                let tid: usize = record.tid() as usize;
                let chrom: String = t_names[record.tid() as usize].clone();

                if tid != prev_tid {
                    let mut ref_seq_u8: Vec<u8> = vec![];
                    fasta.read_all(&chrom, &mut ref_seq_u8).expect("Failed to read fasta sequence record.");
                    ref_seq = dna_vec(&ref_seq_u8);
                }

                let start_pos = record.pos();
                let end_pos =
                    record.cigar().end_pos().expect("Error while accessing CIGAR end position") - 1;

                let interval = GenomicInterval {
                    chrom: chrom,
                    start_pos: start_pos as u32,
                    end_pos: end_pos as u32,
                };

                // get the list of variants that overlap this read
                let read_vars = varlist.get_variants_range(interval);

                let frag = extract_fragment(&record,
                                            read_vars,
                                            &ref_seq,
                                            &t_names,
                                            extract_params,
                                            align_params);

                match frag {
                    Some(f) => {flist.push(f);}
                    None => {}
                }

                prev_tid = tid;
            }
        }
        &None => {
            let bam = bam::Reader::from_path(bamfile_name).unwrap();
            for r in bam.records() {
                let record = r.unwrap();

                let tid: usize = record.tid() as usize;
                let chrom: String = t_names[record.tid() as usize].clone();

                if tid != prev_tid {
                    let mut ref_seq_u8: Vec<u8> = vec![];
                    fasta.read_all(&chrom, &mut ref_seq_u8).expect("Failed to read fasta sequence record.");
                    ref_seq = dna_vec(&ref_seq_u8);
                }

                let start_pos = record.pos();
                let end_pos =
                    record.cigar().end_pos().expect("Error while accessing CIGAR end position") - 1;

                let interval = GenomicInterval {
                    chrom: chrom,
                    start_pos: start_pos as u32,
                    end_pos: end_pos as u32,
                };

                // get the list of variants that overlap this read
                let read_vars = varlist.get_variants_range(interval);

                let frag = extract_fragment(&record,
                                            read_vars,
                                            &ref_seq,
                                            &t_names,
                                            extract_params,
                                            align_params);

                match frag {
                    Some(f) => {flist.push(f);}
                    None => {}
                }

                prev_tid = tid;
            }
        }
    };


    flist
}

//************************************************************************************************
// END OF RUST-HTSLIB BASED CODE *****************************************************************
//************************************************************************************************


#[cfg(test)]
mod tests {
    use super::*;

    // helper function for more efficient testing of check_kmers_unique function
    fn test_check_kmers_helper(test_str: String, k: usize, exp_result: bool) {
        let mut test_vec: Vec<char> = vec![];
        let mut test_vec_backup: Vec<char> = vec![];

        for c in test_str.chars() {
            let u: char = c as char;
            test_vec.push(u);
            test_vec_backup.push(u);
        }
        // test that function computes the expected result
        assert_eq!(exp_result, check_kmers_unique(&test_vec, k));
        // check that the original vector is not consumed or mutated by comparing to a backup
        assert_eq!(test_vec, test_vec_backup);
    }

    #[test]
    fn test_check_kmers_unique_expect_pass_1() {
        // expect true. basic test, all 3-mers are unique
        test_check_kmers_helper("AAATGCGCCA".to_string(), 3, true);
    }

    #[test]
    fn test_check_kmers_unique_expect_pass_2() {
        // expect true. AAAT, a 4-mer, is repeated, but no 5-mers are repeated
        test_check_kmers_helper("AAATGCGAAATCCAAGAATGA".to_string(), 5, true);
    }

    #[test]
    fn test_check_kmers_unique_expect_fail_1() {
        //expect false because AAAT, a 4-mer, is repeated.
        test_check_kmers_helper("AAATGCGAAATCCA".to_string(), 4, false);
    }

    #[test]
    fn test_check_kmers_unique_expect_fail_end() {
        // check that result is false if last kmer ('GAATGA') is repeated.
        test_check_kmers_helper("AAATGAATGACAAGAATGA".to_string(), 6, false);
    }

    #[test]
    fn test_check_kmers_unique_expect_fail_beginning() {
        // check that result is false if first kmer ('GAATGA') is repeated.
        test_check_kmers_helper("GAATGAAAATGAATGACAAG".to_string(), 6, false);
    }

    #[test]
    fn test_check_kmers_unique_expect_fail_overlap() {
        // check that result is false if overlapping kmers ('TGATG' -> 'TGATGATG') are repeated
        test_check_kmers_helper("CCCCTGATGATGACCCC".to_string(), 5, false);
    }

    #[test]
    fn test_check_kmers_unique_expect_pass_3() {
        // check that overlap test returns true if we insert a 'C'
        test_check_kmers_helper("CCCCTGATCATGACCCC".to_string(), 5, true);
    }
}
