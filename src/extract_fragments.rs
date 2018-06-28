use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam::record::Cigar;
use std::error::Error;
use util::*;
use variants_and_fragments::*;
use std::collections::HashMap;
use bio::stats::{LogProb, Prob, PHREDProb};
use bio::io::fasta;
use bio::pattern_matching::bndm;
use realignment::*;

static VERBOSE: bool = false;
static IGNORE_INDEL_ONLY_CLUSTERS: bool = false;

#[derive(Clone, Copy)]
pub struct ExtractFragmentParameters {
    pub min_mapq: u8,
    pub alignment_type: AlignmentType,
    pub band_width: usize,
    pub anchor_length: usize,
    pub short_hap_max_snvs: usize,
    pub max_window_padding: usize,
    pub max_cigar_indel: usize
}

// an extension of the rust-htslib cigar representation
// has the cigar operation as well as the position on the reference of the operation start,
// and the position on the read of the operation start
pub struct CigarPos {
    pub cig: Cigar,
    pub ref_pos: u32,
    pub read_pos: u32,
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
            &Cigar::HardClip(_) => {
                j += 1;
            }
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
                    cigarpos_list: &Vec<CigarPos>,
                    var_interval: GenomicInterval,
                    ref_seq: &Vec<char>,
                    read_seq: &Vec<char>,
                    target_names: &Vec<String>,
                    extract_params: ExtractFragmentParameters)
                    -> Result<Option<AnchorPositions>, CigarOrAnchorError> {

    //let min_window_length = extract_params.min_window_length;
    let max_window_padding = extract_params.max_window_padding;

    let anchor_length = extract_params.anchor_length as u32;

    let l_max = if var_interval.start_pos as usize >= max_window_padding {
        var_interval.start_pos as usize - max_window_padding
    } else {
        0
    };

    let r_max = if var_interval.end_pos as usize + max_window_padding < ref_seq.len() {
        var_interval.end_pos as usize + max_window_padding
    } else {
        ref_seq.len() - 1
    };
    let mut ref_seq_max_window: Vec<u8> = vec![];


    for c in ref_seq[l_max..r_max+1].iter() {
        ref_seq_max_window.push(*c as u8);
    }

    if VERBOSE {
        eprintln!("Finding anchors for variant at {} {} {}:",var_interval.chrom, var_interval.start_pos, var_interval.end_pos);
    }

    if var_interval.chrom != target_names[bam_record.tid() as usize] ||
        (var_interval.start_pos as i32) >=
            bam_record.cigar().end_pos().expect("Error while accessing CIGAR end position") ||
        (var_interval.end_pos as i32) < bam_record.pos() {
        eprintln!("var_interval: {}\t{}\t{}",var_interval.chrom, var_interval.start_pos, var_interval.end_pos);
        eprintln!("bam_record:   {}\t{}\t{}",target_names[bam_record.tid() as usize], bam_record.pos(), bam_record.cigar().end_pos().expect("Error while accessing CIGAR end position"));

        return Err(CigarOrAnchorError::AnchorRangeOutsideRead(
            "Attempted to find sequence anchors for a range completely outside of the sequence.".to_owned()
        ));
    }

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

    // how slow is this? could be sped up with binary search in the future.
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
    //println!("**************************************************");
    //println!("searching for left anchor...");

    if VERBOSE {
        eprint!("Finding left anchor: ");
    }

    for i in (0..left_ix + 1).rev() {

        if cigarpos_list[i].ref_pos <= anchor_length
            || cigarpos_list[i].ref_pos >= ref_seq.len() as u32 - anchor_length {
            return Ok(None);
        }

        match cigarpos_list[i].cig {

            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {

                match_len_left += l;
                let mut potential_anchor = false;
                if seen_indel_left && match_len_left > anchor_length {
                    // we have found a sufficiently long match and it's NOT the cigar op containing var_interval.start_pos
                    left_anchor_ref = cigarpos_list[i].ref_pos + match_len_left - anchor_length;
                    left_anchor_read = cigarpos_list[i].read_pos + match_len_left - anchor_length;
                    potential_anchor = true;
                } else if !seen_indel_left &&
                    cigarpos_list[i].ref_pos < var_interval.start_pos - anchor_length {
                    // we have found a sufficiently long match but it's the cigar ops surrounding var_interval.start_pos
                    // the var_interval.start_pos we are trying to anchor is just inside one huge match
                    left_anchor_ref = var_interval.start_pos - anchor_length;
                    left_anchor_read = cigarpos_list[i].read_pos +
                        (var_interval.start_pos - anchor_length -
                            cigarpos_list[i].ref_pos);
                    potential_anchor = true;
                }

                if VERBOSE {
                    eprint!("M/=/X {},", l);
                }

                if potential_anchor {
                    while left_anchor_ref >= cigarpos_list[i].ref_pos {

                        let l_anc: usize = left_anchor_ref as usize;
                        let r_anc: usize = (left_anchor_ref + anchor_length) as usize;
                        if l_anc <= 0 || r_anc >= ref_seq.len() {
                            break;
                        }

                        // check if there is an exact match between anchor sequence on read and ref
                        let anchor_on_read: Vec<char> = read_seq[(left_anchor_read as usize)..
                            (left_anchor_read + anchor_length) as usize]
                            .to_vec();
                        assert_eq!(anchor_on_read.len(),anchor_length as usize);
                        let anchor_on_ref: Vec<char> = ref_seq[l_anc..r_anc].to_vec();
                        let anchor_match = anchor_on_read == anchor_on_ref;

                        // check that the anchor sequence is unique in the region
                        let mut pattern: Vec<u8>= vec![];
                        for c in anchor_on_ref.iter() {
                            pattern.push(*c as u8);
                        };
                        assert_eq!(pattern.len(), anchor_length as usize);
                        let bndm = bndm::BNDM::new(&pattern);
                        let occ: Vec<usize> = bndm.find_all(&ref_seq_max_window).collect();

                        //println!("ref_seq_max_window: {}", String::from_utf8(ref_seq_max_window.clone()).unwrap());
                        //println!("anchor_on_ref: {}", String::from_utf8(pattern.clone()).unwrap());
                        //let s: String = anchor_on_read.into_iter().collect();
                        //println!("anchor_on_read: {}", s);
                        //println!("pattern: {}", String::from_utf8(pattern.clone()).unwrap());
                        //println!("anchor_match: {} occ.len(): {} l_anc <= l_max: {}", anchor_match, occ.len(), l_anc <= l_max);
                        //println!("***************");

                        if (anchor_match && occ.len() == 1) || l_anc <= l_max {
                            found_anchor_left = true;
                            break;
                        }

                        left_anchor_ref -= 1;
                        left_anchor_read -= 1;
                    }
                    if found_anchor_left {
                        break;
                    }
                }
            }
            Cigar::Ins(l) |
            Cigar::Del(l) |
            Cigar::RefSkip(l) => {
                if l > extract_params.max_cigar_indel as u32 {
                    return Ok(None); // cigar indel too long
                }
                if VERBOSE {
                    eprint!("I/D/N {},", l);
                }
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

    if VERBOSE {
        eprintln!("");
    }

    if !found_anchor_left {
        return Ok(None); // failed to find a left anchor
    }

    //println!("**************************************************");
    if VERBOSE {
        eprint!("Finding right anchor: ");
    }    // step forwards to find right anchor
    let mut match_len_right = 0;
    let mut seen_indel_right = false;
    let mut found_anchor_right = false;
    for i in right_ix..cigarpos_list.len() {

        if cigarpos_list[i].ref_pos <= anchor_length
            || cigarpos_list[i].ref_pos >= ref_seq.len() as u32 - anchor_length {
            return Ok(None);
        }

        match cigarpos_list[i].cig {
            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                match_len_right += l;
                let mut potential_anchor = false;
                if seen_indel_right && match_len_right > anchor_length {
                    // we have found a sufficiently long match and it's NOT the cigar op containing var_interval.end_pos
                    right_anchor_ref = cigarpos_list[i].ref_pos + l - match_len_right +
                        anchor_length;
                    right_anchor_read = cigarpos_list[i].read_pos + l - match_len_right +
                        anchor_length;
                    potential_anchor = true;
                } else if !seen_indel_right &&
                    cigarpos_list[i].ref_pos + l > var_interval.end_pos + anchor_length {
                    // we have found a sufficiently long match but it's the cigar ops surrounding var_interval.end_pos
                    // the var_interval.end_pos we are trying to anchor is just inside one huge match
                    right_anchor_ref = var_interval.end_pos + anchor_length;
                    right_anchor_read = cigarpos_list[i].read_pos +
                        (var_interval.end_pos + anchor_length -
                            cigarpos_list[i].ref_pos);
                    potential_anchor = true;
                }

                if VERBOSE {
                    eprint!("M/=/X {}, ", l);
                }

                if potential_anchor {
                    while right_anchor_ref < cigarpos_list[i].ref_pos + l {

                        let l_anc: usize = (right_anchor_ref - anchor_length) as usize;
                        let r_anc: usize = right_anchor_ref as usize;
                        if l_anc <= 0 || r_anc >= ref_seq.len() {
                            break;
                        }

                        // check if there is an exact match between anchor sequence on read and ref
                        let anchor_on_read: Vec<char> = read_seq[(right_anchor_read - anchor_length) as usize..
                            right_anchor_read as usize]
                            .to_vec();
                        assert_eq!(anchor_on_read.len(),anchor_length as usize);
                        let anchor_on_ref: Vec<char> = ref_seq[l_anc..r_anc].to_vec();
                        let anchor_match = anchor_on_read == anchor_on_ref;

                        // check that the anchor sequence is unique in the region
                        let mut pattern: Vec<u8> = vec![];
                        for c in anchor_on_ref.iter() {
                            pattern.push(*c as u8);
                        };
                        assert_eq!(pattern.len(), anchor_length as usize);
                        let bndm = bndm::BNDM::new(&pattern);
                        let occ: Vec<usize> = bndm.find_all(&ref_seq_max_window).collect();

                        //println!("ref_seq_max_window: {}", String::from_utf8(ref_seq_max_window.clone()).unwrap());
                        //println!("anchor_on_ref: {}", String::from_utf8(pattern.clone()).unwrap());
                        //let s: String = anchor_on_read.into_iter().collect();
                        //println!("anchor_on_read: {}", s);
                        //println!("pattern: {}", String::from_utf8(pattern.clone()).unwrap());
                        //println!("anchor_match: {} occ.len(): {} r_anc >= r_max: {}", anchor_match, occ.len(), r_anc >= r_max);
                        //println!("***************");

                        if (anchor_match && occ.len() == 1) || r_anc >= r_max {
                            found_anchor_right = true;
                            break;
                        }

                        right_anchor_ref += 1;
                        right_anchor_read += 1;
                    }
                    if found_anchor_right {
                        break;
                    }
                }
            }
            Cigar::Ins(l) |
            Cigar::Del(l) |
            Cigar::RefSkip(l) => {
                if l > extract_params.max_cigar_indel as u32 {
                    return Ok(None); // cigar indel too long
                }
                if VERBOSE {
                    eprint!("I/D/N {}, ", l);
                }
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

    if VERBOSE {
        eprintln!("");
    }

    if !found_anchor_right {
        return Ok(None); // failed to find a right anchor
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

// input: a reference to a vector of variants length n, and a number k
// output: a vector of all possible haplotypes from the k-th variant onward,
// where each haplotype is a vector of u8s representing sequence of alleles
fn generate_haps_k_onward(var_cluster: &Vec<Var>, k: usize) -> Vec<Vec<u8>> {

    if k >= var_cluster.len() {
        return vec![vec![]];
    }

    let hap_suffixes = generate_haps_k_onward(var_cluster, k + 1);
    let c = (hap_suffixes.len() * var_cluster[k].alleles.len()) as usize;
    let mut hap_list: Vec<Vec<u8>> = Vec::with_capacity(c);

    for a in 0..var_cluster[k].alleles.len() {
        for mut hap_suffix in hap_suffixes.clone() {
            let mut hap = Vec::with_capacity(var_cluster.len() - k);
            hap.push(a as u8);
            hap.append(&mut hap_suffix);
            hap_list.push(hap);
        }
    }

    hap_list
}

// input: a reference to a vector of variants length n
// output: a vector of all possible haplotypes
// where each haplotype is a vector of u8s length n representing sequence of alleles for the variants
fn generate_haps(var_cluster: &Vec<Var>) -> Vec<Vec<u8>> {
    generate_haps_k_onward(var_cluster, 0)
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
    let mut max_hap: Vec<u8> = vec![0u8; var_cluster.len()];
    let n_vars: usize = var_cluster.len() as usize; // number of variants in cluster
    assert!(n_vars <= extract_params.short_hap_max_snvs);

    // allele_scores[i][j] contains the Log sum of the probabilities of all short haplotypes
    // that had the jth allele at the ith variant of the cluster.
    let mut allele_scores: Vec<Vec<LogProb>> = Vec::with_capacity(n_vars);
    for v in 0..n_vars {
        let n_alleles = var_cluster[v].alleles.len();
        allele_scores.push(vec![LogProb::ln_zero(); n_alleles]);
    }

    let mut score_total: LogProb = LogProb::ln_zero();

    if VERBOSE {
        for var in var_cluster.clone() {
            eprint!("{} {}",
                     var.chrom,
                     var.pos0);
            for allele in var.alleles {
                eprint!(" {}", allele);
            }
            eprintln!("");
        }
        let read_seq_str: String = read_window.clone().into_iter().collect();
        eprintln!("read: {}", read_seq_str);
    }

    let haps = generate_haps(&var_cluster);

    for ref hap in haps {

        assert!(hap.len() > 0);
        let mut hap_window: Vec<char> = Vec::with_capacity(window_capacity);
        let mut i: usize = anchors.left_anchor_ref as usize;
        for var in 0..n_vars {
            while i < var_cluster[var].pos0 {
                hap_window.push(ref_seq[i]);
                i += 1;
            }

            for c in var_cluster[var].alleles[hap[var] as usize].chars() {
                hap_window.push(c);
            }

            i += var_cluster[var].alleles[0].len();
        }

        while i <= anchors.right_anchor_ref as usize {
            hap_window.push(ref_seq[i]);
            i += 1;
        }

        // we now want to score hap_window
        let score: LogProb = match extract_params.alignment_type{
            AlignmentType::NumericallyStableAllAlignment => {
                sum_all_alignments_numerically_stable(&read_window,
                                                                   &hap_window,
                                                                   align_params.ln(),
                                                                   extract_params.band_width)
            }
            AlignmentType::FastAllAlignment => {
                sum_all_alignments(&read_window,
                                                &hap_window,
                                                align_params,
                                                extract_params.band_width)
            }
            AlignmentType::MaxAlignment => {
                max_alignment(&read_window,
                                           &hap_window,
                                           align_params.ln(),
                                           extract_params.band_width)
            }
        };

        assert!(score > LogProb::ln_zero());

        for var in 0..n_vars {

            allele_scores[var][hap[var] as usize] = LogProb::ln_add_exp(allele_scores[var][hap[var] as usize], score);

        }
        if VERBOSE {
            let hap_seq_str: String = hap_window.into_iter().collect();
            eprintln!("hap:{:?} {} PHRED: {}", hap, hap_seq_str, *PHREDProb::from(score));
        }
        // add current alignment score to the total score sum
        score_total = LogProb::ln_add_exp(score_total, score);

        if score > max_score {
            max_score = score;
            max_hap = hap.clone();
        }
    }


    for v in 0..n_vars {

        let best_allele = max_hap[v];
        assert_ne!(allele_scores[v][best_allele as usize], LogProb::ln_zero());
        let mut qual = LogProb::ln_one_minus_exp(&(allele_scores[v][best_allele as usize] - score_total));

        //assert_ne!(qual, LogProb::ln_zero());

        // TODO: BUG: qual should never ever be 0.
        // need to investigate why this happens
        if qual == LogProb::ln_zero() {
            //eprintln!("WARNING: Qual being set to ln(0.00001) due to zero-probability value.");
            qual = LogProb::from(Prob(0.00001));
        }

        if VERBOSE {
            eprint!("adding call: {} {}", var_cluster[v].chrom, var_cluster[v].pos0);
            for allele in &var_cluster[v].alleles {
                eprint!(" {}", allele);
            }

            eprint!("; allele = {};", best_allele);
            eprintln!(" qual = {};", *Prob::from(qual));

        }

        calls.push(FragCall {
            frag_ix: None,
            var_ix: var_cluster[v].ix,
            allele: best_allele,
            qual: qual,
            one_minus_qual: LogProb::ln_one_minus_exp(&qual)
        });
    }

    if VERBOSE {
        eprintln!("--------------------------------------");
    }

    calls
}


pub fn extract_fragment(bam_record: &Record,
                        cigarpos_list: &Vec<CigarPos>,
                        vars: Vec<Var>,
                        ref_seq: &Vec<char>,
                        target_names: &Vec<String>,
                        extract_params: ExtractFragmentParameters,
                        align_params: AlignmentParameters,
                        old_frag: Option<Fragment>)
                        -> Fragment {

    // TODO assert that every single variant in vars is on the same chromosome
    let id: String = u8_to_string(bam_record.qname());

    if VERBOSE {eprintln!("Extracting fragment for read {}...", id);}

    let mut fragment = Fragment {
        id: id,
        calls: vec![],
        p_read_hap: [LogProb::from(Prob(0.5)),LogProb::from(Prob(0.5))]
    };

    if bam_record.is_quality_check_failed() || bam_record.is_duplicate() ||
        bam_record.is_secondary() || bam_record.is_unmapped() || bam_record.mapq() < extract_params.min_mapq
        || bam_record.is_supplementary(){
        return fragment;
    }

    let read_seq: Vec<char> = dna_vec(&bam_record.seq().as_bytes());
    let mut cluster_lst: Vec<(AnchorPositions, Vec<Var>)> = vec![];
    let mut var_anchor_lst: Vec<(Var, AnchorPositions)> = vec![];

    // populate a list with tuples of each variant, and anchor sequences for its alignment
    for ref var in vars {
        let var_interval = GenomicInterval{
            tid: var.tid as u32,
            chrom: var.chrom.clone(),
            start_pos: var.pos0 as u32,
            end_pos: var.pos0 as u32
        };
        match find_anchors(&bam_record,
                                   &cigarpos_list,
                                   var_interval,
                                   &ref_seq,
                                   &read_seq,
                                   &target_names,
                                   extract_params).expect("CIGAR or Anchor Error while finding anchor sequences.") {
            Some(anchors) => {var_anchor_lst.push((var.clone(), anchors));},
            _ => {}
        };
    }

    // now that we have anchors for each var the read covers,
    // group the variants into clusters to align together if adjacent anchors overlap
    // we generate anchors for the whole cluster by taking the first-left and last-right anchor pos of the cluster
    {
        // populate cluster_lst with variant clusters
        let mut var_cluster: Vec<Var> = vec![];
        let mut var_anchors: Vec<AnchorPositions> = vec![];
        let mut l;

        // generate clusters of SNVs that should be considered
        for (var, anc) in var_anchor_lst {
            l = var_cluster.len();
            if l == 0 ||
                (anc.left_anchor_ref < var_anchors[l - 1].right_anchor_ref
                    && l < extract_params.short_hap_max_snvs) {
                var_cluster.push(var);
                var_anchors.push(anc);
            } else {

                // sequence anchor that covers the whole cluster of variants
                let combined_anchor = AnchorPositions{
                    left_anchor_ref: var_anchors[0].left_anchor_ref,
                    right_anchor_ref: var_anchors[l - 1].right_anchor_ref,
                    left_anchor_read: var_anchors[0].left_anchor_read,
                    right_anchor_read: var_anchors[l - 1].right_anchor_read
                };

                cluster_lst.push((combined_anchor, var_cluster.clone()));

                var_cluster.clear();
                var_anchors.clear();
                var_cluster.push(var);
                var_anchors.push(anc);
            }
        }

        l = var_cluster.len();

        if l > 0 {
            // sequence anchor that covers the whole cluster of variants
            let combined_anchor = AnchorPositions{
                left_anchor_ref: var_anchors[0].left_anchor_ref,
                right_anchor_ref: var_anchors[l - 1].right_anchor_ref,
                left_anchor_read: var_anchors[0].left_anchor_read,
                right_anchor_read: var_anchors[l - 1].right_anchor_read
            };

            cluster_lst.push((combined_anchor, var_cluster.clone()));
        }
    }

    // now extract alleles for the variant cluster

    let mut old_frag_ix = 0; // hold onto index in old fragment
    // save reallocation of this hashmap for every variant
    let mut new_var_ix: HashMap<usize, usize> = HashMap::with_capacity(extract_params.short_hap_max_snvs);

    for (anchors, var_cluster) in cluster_lst {

        if IGNORE_INDEL_ONLY_CLUSTERS {
            let mut seen_snv = false;
            for var in var_cluster.iter() {
                if seen_snv {
                    break;
                }
                for allele in var.alleles[1..].iter() {
                    if allele.len() == var.alleles[0].len() {
                        let mut diff = 0;
                        for (c1, c2) in allele.chars().zip(var.alleles[0].chars()) {
                            if c1 != c2 {
                                diff += 1;
                            }
                        }
                        if diff == 1 {
                            seen_snv = true;
                            break;
                        }
                    }
                }
            }

            if !seen_snv {
                continue;
            }
        }

        // some new variants may have been added after performing POA multiple sequence alignment.
        // if these new variants are too close to existing variant cluster, they will all need to be
        // realigned. otherwise, if this "variant cluster" didn't change, just reuse the old fragment calls
        match &old_frag {
            &Some(ref old_f) => {
                // if var.called is true, then this variant existed before POA
                let mut changed = false;
                for ref var in &var_cluster {
                    if !var.called { // this is a new variant from POA
                        assert_eq!(var.old_ix, None);
                        changed = true;
                        break;
                    }
                }

                // the fragment cluster didn't change, so convert the old varlist indices to new ones
                // and just re-use the calls without realigning anything.
                if !changed {

                    // map old varlist indices to new ones
                    assert!(new_var_ix.is_empty()); // shouldn't ever have to clear this hashmap

                    for ref var in &var_cluster {
                        match var.old_ix {
                            Some(old_ix) => {
                                new_var_ix.insert(old_ix, var.ix);
                            },
                            None => {panic!("Attempted to reuse fragment call in an unchanged variant \
                                 cluster, but var.old_ix was None.")}
                        }
                    }

                    // convert the varlist indices in the old calls and push them to the fragment
                    while !new_var_ix.is_empty() {
                        let call = old_f.calls[old_frag_ix];

                        match new_var_ix.remove(&call.var_ix){
                            Some(new_ix) => {
                                let mut new_call = call.clone();
                                new_call.var_ix = new_ix;
                                fragment.calls.push(new_call);
                            },
                            None => {}
                        }

                        old_frag_ix += 1;
                    }

                    // skip to the next var_cluster or else we would double-add the fragment calls!
                    continue;
                }
            }
            &None => {}
        }

        // extract the calls for the fragment
        for call in extract_var_cluster(&read_seq,
                                        ref_seq,
                                        var_cluster,
                                        anchors,
                                        extract_params,
                                        align_params) {
            fragment.calls.push(call);
        }
    }

    fragment
}

pub fn extract_fragments(bam_file: &String,
                         fastafile_name: &String,
                         varlist: &VarList,
                         interval: &Option<GenomicInterval>,
                         extract_params: ExtractFragmentParameters,
                         align_params: AlignmentParameters,
                         old_flist: Option<Vec<Fragment>>)
                         -> Vec<Fragment> {

    let t_names = parse_target_names(&bam_file);

    let mut prev_tid = 4294967295; // huge value so that tid != prev_tid on first iter
    let mut fasta = fasta::IndexedReader::from_file(fastafile_name).unwrap();
    let mut ref_seq: Vec<char> = vec![];

    let mut flist: Vec<Fragment> = vec![];

    // TODO: this uses a lot of duplicate code, need to figure out a better solution.
    let mut complete = 0;

    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bam_file, interval);
    let mut bam_ix = bam::IndexedReader::from_path(bam_file).unwrap();

    for iv in interval_lst {
        bam_ix.fetch(iv.tid, iv.start_pos, iv.end_pos + 1).ok().expect("Error seeking BAM file while extracting fragments.");

        for (i, r) in bam_ix.records().enumerate() {
            let record = r.unwrap();

            if record.is_quality_check_failed() || record.is_duplicate() ||
                record.is_secondary() || record.is_unmapped() || record.mapq() < extract_params.min_mapq
                || record.is_supplementary() {
                continue;
            }

            let old_frag: Option<Fragment> = match &old_flist {
                &Some(ref fl) => {
                    assert_eq!(fl[i].id, u8_to_string(record.qname()));
                    Some(fl[i].clone())
                }
                &None => { None }
            };

            let tid: usize = record.tid() as usize;
            let chrom: String = t_names[tid].clone();

            if tid != prev_tid {
                let mut ref_seq_u8: Vec<u8> = vec![];
                fasta.read_all(&chrom, &mut ref_seq_u8).expect("Failed to read fasta sequence record.");
                ref_seq = dna_vec(&ref_seq_u8);
            }

            let start_pos = record.pos();
            let end_pos =
                record.cigar().end_pos().expect("Error while accessing CIGAR end position") - 1;

            let bam_cig: CigarStringView = record.cigar();
            let cigarpos_list: Vec<CigarPos> =
                create_augmented_cigarlist(start_pos as u32, &bam_cig).expect("Error creating augmented cigarlist.");


            let interval = GenomicInterval {
                tid: tid as u32,
                chrom: chrom,
                start_pos: start_pos as u32,
                end_pos: end_pos as u32,
            };

            // get the list of variants that overlap this read
            let read_vars = varlist.get_variants_range(interval);

            // print the percentage of variants processed every 10%
            if read_vars.len() > 0 && ((read_vars[0].ix as f64 / varlist.lst.len() as f64) * 10.0) as usize > complete {
                complete = ((read_vars[0].ix as f64 / varlist.lst.len() as f64) * 10.0) as usize;
                if complete < 10 {
                    eprintln!("{}    {}% of variants processed...", print_time(), complete * 10);
                }
            }

            let frag = extract_fragment(&record,
                                        &cigarpos_list,
                                        read_vars,
                                        &ref_seq,
                                        &t_names,
                                        extract_params,
                                        align_params,
                                        old_frag);

            flist.push(frag);

            prev_tid = tid;
        }
    }
    eprintln!("{}    100% of variants processed.", print_time());


    // label every fragment call with its index in the fragment list.
    for i in 0..flist.len() {

        for j in 0..flist[i].calls.len() {
            flist[i].calls[j].frag_ix = Some(i);
        }
    }

    flist
}

//************************************************************************************************
// END OF RUST-HTSLIB BASED CODE *****************************************************************
//************************************************************************************************


#[cfg(test)]
mod tests {
    use super::*;
    use genotype_probs::*;

    /*
    #[test]
    fn test_extended_cigar() {
        let mut fasta = fasta::IndexedReader::from_file(&"/home/peter/git/reaper/study/data/genomes/hs37d5.fa").unwrap();

        let mut ref_seq_u8: Vec<u8> = vec![];
        fasta.read_all(&"20", &mut ref_seq_u8).expect("Failed to read fasta sequence record.");
        let ref_seq: Vec<char> = dna_vec(&ref_seq_u8);

        let mut bam = bam::IndexedReader::from_path(&"/home/peter/git/reaper/test_data/test_extended_cigar.bam").unwrap();

        for r in bam.records() {
            let record = r.unwrap();

            let read_seq: Vec<char> = dna_vec(&record.seq().as_bytes());

            let start_pos = record.pos();
            let bam_cig: CigarStringView = record.cigar();
            let cigarpos_list: Vec<CigarPos> =
                create_augmented_cigarlist(start_pos as u32, &bam_cig).expect("Error creating augmented cigarlist.");

            for cigarpos in cigarpos_list.iter() {
                match cigarpos.cig {
                    Cigar::Diff(l) => {
                        for i in 0..l {
                            assert_ne!(ref_seq[cigarpos.ref_pos as usize + i as usize], read_seq[cigarpos.read_pos as usize + i as usize]);
                        }
                    }
                    Cigar::Equal(l) => {
                        for i in 0..l {
                            assert_eq!(ref_seq[cigarpos.ref_pos as usize + i as usize], read_seq[cigarpos.read_pos as usize + i as usize]);
                        }
                    }
                    Cigar::Match(_) |
                    Cigar::Ins(_) |
                    Cigar::Del(_) |
                    Cigar::RefSkip(_) => {}
                    Cigar::Pad(_) |
                    Cigar::Back(_) |
                    Cigar::SoftClip(_) |
                    Cigar::HardClip(_) => {
                        panic!("CIGAR operation found in cigarpos_list that should have been removed already.".to_owned());
                    }
                }
            }
        }
    }
    */
    fn generate_var2(ix: usize, tid: usize, chrom: String, pos0: usize, alleles: Vec<String>) -> Var {
        Var {
            ix: ix,
            old_ix: None,
            tid: tid,
            chrom: chrom,
            pos0: pos0,
            alleles: alleles,
            dp: 40,
            allele_counts: vec![20,20],
            ambiguous_count: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: Genotype(0,1),
            gq: 0.0,
            genotype_post: GenotypeProbs::uniform(2),
            phase_set: None,
            called: true
        }
    }

    #[test]
    fn test_generate_haplotypes_basic() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 1, vec!["A".to_string(), "G".to_string()]));
        let mut haps = generate_haps(&lst1);

        haps.sort();
        let mut exp = vec![vec![0u8],
                           vec![1u8]];

        exp.sort();
        assert_eq!(haps, exp);
    }

    #[test]
    fn test_generate_haplotypes_multivariant() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 1, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["A".to_string(), "T".to_string()]));
        lst1.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut haps = generate_haps(&lst1);

        haps.sort();
        let mut exp = vec![vec![0u8,0u8,0u8],
                       vec![0u8,0u8,1u8],
                       vec![0u8,1u8,0u8],
                       vec![1u8,0u8,0u8],
                       vec![1u8,1u8,0u8],
                       vec![0u8,1u8,1u8],
                       vec![1u8,0u8,1u8],
                       vec![1u8,1u8,1u8]];

        exp.sort();
        assert_eq!(haps, exp);
    }

    #[test]
    fn test_generate_haplotypes_multiallelic() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 1, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["A".to_string(), "T".to_string(), "C".to_string()]));
        let mut haps = generate_haps(&lst1);

        haps.sort();
        let mut exp = vec![vec![0u8,0u8],
                           vec![0u8,1u8],
                           vec![1u8,0u8],
                           vec![1u8,1u8],
                           vec![0u8,2u8],
                           vec![1u8,2u8]];

        exp.sort();
        assert_eq!(haps, exp);
    }

    #[test]
    fn test_generate_haplotypes_multiallelic2() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 1, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["A".to_string(), "T".to_string(), "C".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 200, vec!["A".to_string(), "T".to_string(), "C".to_string()]));

        let mut haps = generate_haps(&lst1);

        haps.sort();
        let mut exp = vec![vec![0u8,0u8,0u8],
                           vec![0u8,0u8,1u8],
                           vec![0u8,0u8,2u8],
                           vec![0u8,1u8,0u8],
                           vec![0u8,1u8,1u8],
                           vec![0u8,1u8,2u8],
                           vec![0u8,2u8,0u8],
                           vec![0u8,2u8,1u8],
                           vec![0u8,2u8,2u8],
                           vec![1u8,0u8,0u8],
                           vec![1u8,0u8,1u8],
                           vec![1u8,0u8,2u8],
                           vec![1u8,1u8,0u8],
                           vec![1u8,1u8,1u8],
                           vec![1u8,1u8,2u8],
                           vec![1u8,2u8,0u8],
                           vec![1u8,2u8,1u8],
                           vec![1u8,2u8,2u8],];

        exp.sort();
        assert_eq!(haps, exp);
    }
}
