//! This module contains functions for estimating the alignment parameters for the sequence alignment
//! Pair Hidden Markov Model (Pair-HMM).
//!
//! For more information about sequence alignment Pair-HMM:
//! Durbin, R., Eddy, S.R., Krogh, A. and Mitchison, G., 1998. Biological sequence analysis:
//! probabilistic models of proteins and nucleic acids. Cambridge university press.

use bio::io::fasta;
use errors::*;
use extract_fragments::{create_augmented_cigarlist, CigarPos};
use realignment::*;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam::Read;
use util::*;
use std::collections::HashMap;

/// represents counts of different transitions from state-to-state in the Pair-HMM
///
/// the counts in this struct are used to keep track of the different alignment transitions that
/// are observed in the BAM file. These counts can be used to help us directly estimate alignment parameters
#[derive(Clone)]
pub struct TransitionCounts {
    pub p: HashMap<(String, String), usize>
}

/// represents counts of different emission events (aligned bases being equal or not equal) in the BAM
#[derive(Clone)]
pub struct EmissionCounts {
    pub p: HashMap<(String, (String, String)), usize>
}

/// a struct to hold counts of transition events and emission events from the BAM file
///
/// together these counts represent all of the counts necessary to estimate the HMM parameters
#[derive(Clone)]
struct AlignmentCounts {
    num_states: usize,
    transition_counts: TransitionCounts,
    emission_counts: EmissionCounts,
    state_lst: Vec<String>,
    kmer_len: usize,
    valid_prev_states: Vec<Vec<usize>>
}

impl TransitionCounts {
    /// converts TransitionCounts to TransitionProbs
    ///
    /// Using this function, HMM transition counts that are empirically observed and counted in the
    /// BAM file can be converted into probabilities that can be used for sequence alignment.
    /// The transition counts out of a single state (e.g. insertion -> other state) are converted
    /// into a fraction of the total transition counts out of that state.
    fn to_probs(&self) -> TransitionProbs {
        // make a hashtable mapping string to usize
        // this will hold the total outgoing for each state
        let mut start_state_counts: HashMap<String, usize> = HashMap::new();
        let mut trans_probs: HashMap<(String,String), f64> = HashMap::new();

        // iterate over all of the transition counts
        // add the count of the transitions to the hash for the starting state
        for (&(ref state1, ref _state2), &count) in &self.p {
            *start_state_counts.entry(state1.clone()).or_insert(0) += count;
        }
        // this means that for each starting state we have the total count of outgoing transitions

        // then iterate over the transition counts again
        // for each transition, divide the count by the total for that starting state
        for (&(ref state1, ref state2), &count) in &self.p {
            let frac = count as f64 / *start_state_counts.get(&state1.clone()).unwrap() as f64;
            trans_probs.insert((state1.clone(),state2.clone()), frac);

        }
        // return as a TransitionProbs
        TransitionProbs {
            p: trans_probs
        }
    }

    /// add the corresponding counts inside two ```TransitionCount```s together
    ///
    /// this function allows us to create separate ```TransitionCount``` structs for each newly
    /// observed read, and then just add those counts onto a grand total count for the BAM file.
    fn add(&mut self, other: TransitionCounts) {
        for (&(ref state1, ref state2), &count) in &other.p {
            *self.p.entry((state1.clone(), state2.clone()))
                .or_insert(0) += count;
        }
    }
}

impl EmissionCounts {
    /// converts EmissionCounts to EmissionProbs
    ///
    /// Using this function, HMM emission counts that are empirically observed and counted in the
    /// BAM file can be converted into probabilities that can be used for sequence alignment.
    /// The emission counts for equal and not equal bases are turned into a fraction of the total.
    /// The insertion and deletion emission probabilities are just fixed to 1.0.
    /// This probably isn't ideal, but it works alright in practice.
    fn to_probs(&self) -> EmissionProbs {
        // make a hashtable mapping string to usize
        // this will hold the total counts of symbol pairs observed at a state
        let mut state_counts: HashMap<String, usize> = HashMap::new();
        let mut emit_probs: HashMap<(String,String), f64> = HashMap::new();

        // iterate over all of the emission counts
        // add the count of the pair of emitted symbols to the count for the state their emitted from
        for (&(ref state, (ref _symbol1, ref _symbol2)), &count) in &self.p {
            *state_counts.entry(state.clone()).or_insert(0) += count;
        }
        // this means that for each state we have the total count of emitted symbol pairs

        // then iterate over the emission counts again
        // for each emission, divide the count by the total emitted for that state
        for (&(ref state, (ref symbol1, ref symbol2)), &count) in &self.p {
            let frac = (count as f64) / (*state_counts.get(&state.clone()).unwrap() as f64);
            emit_probs.insert((symbol1.clone(),symbol2.clone()), frac);
        }

        // return as a TransitionProbs
        EmissionProbs {
            probs: emit_probs
        }
    }

    /// add the corresponding counts inside two ```EmissionCount```s together
    ///
    /// this function allows us to create separate ```EmissionCount``` structs for each newly
    /// observed read, and then just add those counts onto a grand total count for the BAM file.
    fn add(&mut self, other: EmissionCounts) {
        for (&(ref state, (ref symbol1, ref symbol2)), &count) in &other.p {
            *self.p.entry((state.clone(),(symbol1.clone(), symbol2.clone())))
                .or_insert(0) += count;
        }
    }
}

impl AlignmentCounts {
    /// convert ```AlignmentCounts``` into ```AlignmentParameters``` by converting the transition
    /// and emission counts into probabilities
    fn to_parameters(&self) -> AlignmentParameters {
        AlignmentParameters {
            num_states: self.num_states,
            transition_probs: self.transition_counts.to_probs(),
            emission_probs: self.emission_counts.to_probs(),
            state_lst: self.state_lst.clone(),
            kmer_len: self.kmer_len,
            valid_prev_states: self.valid_prev_states.clone()
        }
    }
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

/// Count the number of alignment events (state transitions and emissions) observed for a single
/// read from a BAM file
///
/// #Arguments
/// -```cigarpos_list```: a vector of ```CigarPos``` which hold the CIGAR operations for the read,
///                       along with the reference and read positions where they occur
/// -```ref_seq```: the reference sequence for the chromosome that the read is aligned to
/// -```read_seq```: the sequence of the read from the BAM file
/// -```max_cigar_indel```: the maximum length of a CIGAR operation in order to count it.
///                         this is meant to filter out large indels observed in the BAM alignment
///                         that are due to misalignment or structural variations instead of
///                         random sequencing error.
///
/// #Returns
/// Returns a Result containing a TransitionCounts and an EmissionCounts.
/// The TransitionCounts represents the total count of observed HMM transitions described by the
/// BAM alignment. The EmissionCounts represents the observed HMM emissions (base match / base
/// substitution) in the BAM alignment.
///
/// #Errors
/// Can throw an error if ```cigarpos_list``` contains a cigar operation that should already have
///   been filtered out (Pad,Back,Softclip,Hardclip).
pub fn count_alignment_events(
    cigarpos_list: &Vec<CigarPos>,
    ref_seq: &Vec<char>,
    read_seq: &Vec<char>,
    kmer_len: usize,
    max_cigar_indel: u32,
) -> Result<(TransitionCounts, EmissionCounts)> {

    let mut aln_seq = vec![];
    let mut aln_ref = vec![];

    // for each cigar operation in the BAM record
    for cigarpos in cigarpos_list {
        match cigarpos.cig {
            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                // the current operation is a match

                let mut ref_pos = cigarpos.ref_pos as usize;
                let mut read_pos = cigarpos.read_pos as usize;

                // traverse over the length of the match operation
                for _ in 0..l {
                    if ref_pos >= ref_seq.len() - 1 || read_pos >= read_seq.len() - 1 {
                        break; // break if we've reached the end of the read or reference
                    }

                    // check if the aligned bases match or mismatch and use these to iterate the
                    // emission counts (whether bases match or mismatch)
                    if ref_seq[ref_pos] != 'N' && read_seq[read_pos] != 'N' {
                        aln_seq.push(read_seq[read_pos]);
                        aln_ref.push(ref_seq[ref_pos]);
                    }

                    // it's a match operation so both read and reference move forward one position
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Cigar::Ins(l) => {
                // the current operation is insertion

                if l > max_cigar_indel {
                    continue; // skip this operation if the insertion is too long
                }

                let mut ref_pos = cigarpos.ref_pos as usize;
                let mut read_pos = cigarpos.read_pos as usize;

                // traverse over the length of the insertion operation
                for _ in 0..l {
                    if ref_pos >= ref_seq.len() - 1 || read_pos >= read_seq.len() - 1 {
                        break; // break if we've reached the end of the read or reference
                    }

                    if read_seq[read_pos] != 'N' {
                        aln_seq.push(read_seq[read_pos]);
                        aln_ref.push('-');
                    }

                    // this is an insertion so only the read position is moved forward
                    read_pos += 1;
                }
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => {
                // the current operation is deletion

                if l > max_cigar_indel {
                    continue; // skip this operation if the insertion is too long
                }

                let mut ref_pos = cigarpos.ref_pos as usize;
                let mut read_pos = cigarpos.read_pos as usize;

                // traverse over the length of the deletion operation
                for _ in 0..l {
                    if ref_pos >= ref_seq.len() - 1 || read_pos >= read_seq.len() - 1 {
                        break; // break if we've reached the end of the read or reference
                    }

                    if ref_seq[ref_pos] != 'N' {
                        aln_seq.push('-');
                        aln_ref.push(ref_seq[ref_pos]);
                    }
                    // this is a deletion so only the reference position is moved forward
                    ref_pos += 1;
                }
            }
            Cigar::Pad(_) | Cigar::Back(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) => {
                bail!(
                    "CIGAR operation found in cigarpos_list that should have been removed already."
                );
            }
        }
    }

    assert_eq!(aln_seq.len(), aln_ref.len());

    // iterate over the BAM alignment as a pair of strings
    // e.g.
    // AA--ACTTCC
    // AATGGCT-CC

    let mut prev_state = "-".to_string();
    let mut tcounts: HashMap<(String,String),usize> = HashMap::new();
    let mut ecounts: HashMap<(String,(String,String)),usize> = HashMap::new();

    for i in 0..aln_seq.len()-kmer_len+1 {
        // iterate over the aligned kmers in the alignment
        // e.g. the first pair of aligned 3-mers in the above alignment is
        // (AA-,AAT) and the second is (A--, ATG)
        let kmer1 = aln_seq[i..i+kmer_len].to_vec();
        let kmer2 = aln_ref[i..i+kmer_len].to_vec();
        let mut state = vec![];

        // figure out what the HMM state is for this pair of aligned kmers
        // e.g. (AA-, AAT) -> MMD (match match deletion)
        //      (--A, TGG) -> DDM (deletion deletion match)
        for j in 0..kmer_len {
            if kmer1[j] != '-' && kmer2[j] != '-' {
                state.push('M');
            }else if kmer1[j] == '-' && kmer2[j] != '-' {
                state.push('D');
            }else if kmer1[j] != '-' && kmer2[j] == '-' {
                state.push('I');
            } else {
                panic!("error building emitted kmers from BAM alignment.")
            }
        }

        // convert emitted kmers and state to strings
        let kmer1_str: String = kmer1.into_iter().collect();
        let kmer2_str: String = kmer2.into_iter().collect();
        let state_str: String = state.into_iter().collect();

        // add the emission of these two aligned kmers, in the current state, to the emission counts
        *ecounts
            .entry((state_str.clone(),(kmer1_str.clone(), kmer2_str.clone())))
            .or_insert(0) += 1;

        if i > 0 {
        // add the transition from the previous state, to the current state, to the transition counts
        *tcounts
            .entry((prev_state.clone(), state_str.clone()))
                       .or_insert(0) += 1;
        }

        // assign the current state to the previous state for next iteration
        prev_state = state_str.clone();
    }

    return Ok((TransitionCounts{p:tcounts}, EmissionCounts{p:ecounts}));
}

fn get_all_kmers_vec(len: usize) -> Vec<Vec<char>> {
    if len == 0 {
        return vec![vec![]];
    }

    // we're defining the set of all kmers recursively. The set of all kmers of length k
    // is equal to the set of kmers length k-1, with the 5 ending symbols added to each one
    let kmers_prefixes = get_all_kmers_vec(len-1);

    let mut kmer_lst = vec![];

    // for each of the possible kmer prefixes, add the 5 possible suffix symbols to each one.
    for prefix in kmers_prefixes {
        let mut prefix_a = prefix.clone(); prefix_a.push('A');
        let mut prefix_c = prefix.clone(); prefix_c.push('C');
        let mut prefix_g = prefix.clone(); prefix_g.push('G');
        let mut prefix_t = prefix.clone(); prefix_t.push('T');
        let mut prefix_gap = prefix.clone(); prefix_gap.push('-');

        kmer_lst.push(prefix_a);
        kmer_lst.push(prefix_c);
        kmer_lst.push(prefix_g);
        kmer_lst.push(prefix_t);
        kmer_lst.push(prefix_gap);
    }

    kmer_lst
}


fn get_all_kmers(len: usize) -> Vec<String> {
    let mut kmer_lst: Vec<String> = vec![];

    for kmer_vec in get_all_kmers_vec(len).iter() {
        kmer_lst.push(kmer_vec.into_iter().collect::<String>());
    }

    kmer_lst
}

fn get_all_states_vec(len: usize) -> Vec<Vec<char>> {
    if len == 0 {
        return vec![vec![]];
    }

    // we're defining the set of all states recursively. The set of all states of length k
    // is equal to the set of states length k-1, with the 5 ending symbols added to each one
    let states_prefixes = get_all_states_vec(len-1);

    let mut state_lst = vec![];

    // for each of the possible state prefixes, add the 5 possible suffix symbols to each one.
    for prefix in states_prefixes {
        let mut prefix_m = prefix.clone(); prefix_m.push('M');
        let mut prefix_d = prefix.clone(); prefix_d.push('D');
        let mut prefix_i = prefix.clone(); prefix_i.push('I');

        state_lst.push(prefix_m);
        state_lst.push(prefix_d);
        state_lst.push(prefix_i);
    }

    state_lst
}


fn get_all_states(len: usize) -> Vec<String> {
    let mut state_lst: Vec<String> = vec![];

    for state_vec in get_all_states_vec(len).iter() {
        state_lst.push(state_vec.into_iter().collect::<String>());
    }

    state_lst
}

//************************************************************************************************
// END OF RUST-HTSLIB BASED CODE *****************************************************************
//************************************************************************************************

/// Estimates the alignment parameters needed for the Pair Hidden Markov Model (Pair-HMM) directly
/// from the alignments in a BAM file
///
/// #Arguments
/// -```bam_file```: the input BAM file name
/// -```fasta_file```: the input FASTA file name
/// -```interval```: the (optional) GenomicInterval within which variants should be called
///                  the reads that are used for estimating the alignment parameters are also
///                  limited to this region.
/// -```min_mapq```: the minimium mapping quality to use a read
/// -```max_cigar_indel```: the maximum length of a CIGAR operation in order to count it.
///                         this is meant to filter out large indels observed in the BAM alignment
///                         that are due to misalignment or structural variations instead of
///                         random sequencing error.
///
/// #Returns
/// Returns a result contain an ```AlignmentParameters``` struct. This struct contains the alignment
/// parameters estimated from the BAM file.
///
/// #Errors
/// - ```IndexedFastaOpenError```: error opening the indexed FASTA file
/// - ```IndexedBamOpenError```: error opening the indexed BAM file
/// - ```IndexedBamFetchError```: error fetching region from the indexed BAM file
/// - ```IndexedBamRecordReadError```: error reading a record from the BAM
/// - ```IndexedFastaReadError```: error reading a record from the FASTA
/// - Any errors incurred while creating the augmented cigar list or counting alignment events.
pub fn estimate_alignment_parameters(
    bam_file: &String,
    fasta_file: &String,
    interval: &Option<GenomicInterval>,
    kmer_len: usize,
    min_mapq: u8,
    max_cigar_indel: u32,
) -> Result<AlignmentParameters> {
    let t_names = parse_target_names(&bam_file)?;

    let mut prev_tid = 4294967295; // huge value so that tid != prev_tid on first iter
    let mut fasta = fasta::IndexedReader::from_file(fasta_file)
        .chain_err(|| ErrorKind::IndexedFastaOpenError)?;
    let mut ref_seq: Vec<char> = vec![];

    // initial transition and emission counts
    // set everything to 1 so that it's impossible to have e.g. divide by 0 errors
    let mut init_trans: HashMap<(String,String),usize> = HashMap::new();
    let mut init_emit: HashMap<(String, (String,String)), usize> = HashMap::new();

    let possible_states: Vec<String> = get_all_states(kmer_len);
    let possible_kmers: Vec<String> = get_all_kmers(kmer_len);

    for kmer1 in &possible_kmers {
        for kmer2 in &possible_kmers {

            let mut invalid = false;
            let mut state: Vec<char> = vec![];

            for (c1,c2) in kmer1.chars().zip(kmer2.chars()) {
                if c1 != '-' && c2 != '-' {
                    state.push('M');
                } else if c1 != '-' && c2 == '-' {
                    state.push('I');
                } else if c1 == '-' && c2 != '-' {
                    state.push('D');
                } else {
                    invalid = true;
                }
            }

            if invalid {
                continue;
            }

            let state_str: String = state.into_iter().collect();

            init_emit.insert((state_str.clone(),(kmer1.clone(),kmer2.clone())),1);
        }
    }


    let mut emission_counts = EmissionCounts {p: init_emit.clone()};

    for state1 in &possible_states {
        for state2 in &possible_states {
            init_trans.insert((state1.clone(),state2.clone()),1);
        }
    }

    let mut transition_counts = TransitionCounts {p: init_trans.clone()};

    // interval_lst has either the single specified genomic region, or list of regions covering all chromosomes
    // for more information about this design decision, see get_interval_lst implementation in util.rs
    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bam_file, interval)?;
    let mut bam_ix =
        bam::IndexedReader::from_path(bam_file).chain_err(|| ErrorKind::IndexedBamOpenError)?;

    for iv in interval_lst {
        bam_ix
            .fetch(iv.tid, iv.start_pos, iv.end_pos + 1)
            .chain_err(|| ErrorKind::IndexedBamFetchError)?;

        for r in bam_ix.records() {
            let record =
                r.chain_err(|| ErrorKind::IndexedBamRecordReadError)?;

            // check that the read doesn't fail any standard filters
            if record.mapq() < min_mapq
                || record.is_unmapped()
                || record.is_secondary()
                || record.is_quality_check_failed()
                || record.is_duplicate()
                || record.is_supplementary()
            {
                continue;
            }

            // if we're on a different contig/chrom than the previous BAM record, we need to read
            // in the sequence for that contig/chrom from the FASTA into the ref_seq vector
            let tid: usize = record.tid() as usize;
            let chrom: String = t_names[record.tid() as usize].clone();
            if tid != prev_tid {
                let mut ref_seq_u8: Vec<u8> = vec![];
                fasta
                    .read_all(&chrom, &mut ref_seq_u8)
                    .chain_err(|| ErrorKind::IndexedFastaReadError)?;
                ref_seq = dna_vec(&ref_seq_u8);
            }

            // create a vector of CigarPos from the BAM record
            // these contain the CIGAR operation as well as the reference and read positions
            // where it occurs in the BAM alignment
            let read_seq: Vec<char> = dna_vec(&record.seq().as_bytes());
            let bam_cig: CigarStringView = record.cigar();
            let cigarpos_list: Vec<CigarPos> =
                create_augmented_cigarlist(record.pos() as u32, &bam_cig)
                    .chain_err(|| "Error creating augmented cigarlist.")?;

            // count the emission and transition events directly from the record's CIGAR and sequences
            let (read_transition_counts, read_emission_counts) =
                count_alignment_events(&cigarpos_list, &ref_seq, &read_seq, kmer_len, max_cigar_indel)
                    .chain_err(|| "Error counting cigar alignment events.")?;

            // add emission and transition counts to the running total
            transition_counts.add(read_transition_counts);
            emission_counts.add(read_emission_counts);

            prev_tid = tid;
        }
    }

    let mut valid_prev_states: Vec<Vec<usize>> = vec![vec![]; possible_states.len()];

    for i in 0..possible_states.len() {
        let mut state_vec: Vec<char> = possible_states[i].chars().collect();
        state_vec.pop();

        let mut state_vec_m = state_vec.clone();
        state_vec_m.insert(0, 'M');
        let state_vec_m_str: String = state_vec_m.iter().collect();

        let mut state_vec_d = state_vec.clone();
        state_vec_d.insert(0, 'D');
        let state_vec_d_str: String = state_vec_d.iter().collect();

        let mut state_vec_i = state_vec.clone();
        state_vec_i.insert(0, 'I');
        let state_vec_i_str: String = state_vec_i.iter().collect();

        for j in 0..possible_states.len() {
            if possible_states[j] == state_vec_m_str {
                valid_prev_states[i].push(j)
            } else if possible_states[j] == state_vec_d_str {
                valid_prev_states[i].push(j)
            } else if possible_states[j] == state_vec_i_str {
                valid_prev_states[i].push(j)
            }
        }
        assert!(valid_prev_states[i].len() == 3)
    }

    // place the transition and emission counts together in an AlignmentCounts struct
    let alignment_counts = AlignmentCounts {
        num_states: possible_states.len(),
        transition_counts: transition_counts,
        emission_counts: emission_counts,
        state_lst: possible_states.clone(),
        kmer_len: kmer_len,
        valid_prev_states: valid_prev_states.clone()
    };

    // convert the alignment counts from the BAM into probabilities
    let params = alignment_counts.to_parameters();

    // print the estimated alignment parameters to STDERR
    eprintln!("{} Done estimating alignment parameters.", print_time());
    eprintln!("");

    eprintln!("{} Transition Probabilities:", SPACER);
    for (&(ref state1, ref state2), &prob)in &params.transition_probs.p {
        eprintln!(
            "{} {} -> {}:          {:.10}",
            SPACER, state1, state2, prob
        );
    }

    eprintln!("");

    eprintln!("{} Emission Probabilities:", SPACER);
    for (&(ref kmer1, ref kmer2), &prob) in &params.emission_probs.probs {
        eprintln!(
            "{} ({},{}):          {:.10}",
            SPACER, kmer1, kmer2, prob
        );
    }

    eprintln!("");

    Ok(params)
}

mod tests {
    //use super::*;

    #[test]
    fn test_get_all_states0() {
        let states: Vec<String> = get_all_states(0);
        let expected = vec![
            "".to_string()
        ];
        assert_eq!(states, expected);
    }

    #[test]
    fn test_get_all_states1() {
        let states: Vec<String> = get_all_states(1);
        let expected = vec![
            "M".to_string(),
            "D".to_string(),
            "I".to_string()
        ];
        assert_eq!(states, expected);
    }

    #[test]
    fn test_get_all_states2() {
        let states: Vec<String> = get_all_states(2);
        let expected = vec![
            "MM".to_string(),
            "MD".to_string(),
            "MI".to_string(),
            "DM".to_string(),
            "DD".to_string(),
            "DI".to_string(),
            "IM".to_string(),
            "ID".to_string(),
            "II".to_string(),
        ];
        assert_eq!(states, expected);
    }

    #[test]
    fn test_get_all_states3() {
        let states: Vec<String> = get_all_states(3);
        let expected = vec![
            "MMM".to_string(),
            "MMD".to_string(),
            "MMI".to_string(),
            "MDM".to_string(),
            "MDD".to_string(),
            "MDI".to_string(),
            "MIM".to_string(),
            "MID".to_string(),
            "MII".to_string(),
            "DMM".to_string(),
            "DMD".to_string(),
            "DMI".to_string(),
            "DDM".to_string(),
            "DDD".to_string(),
            "DDI".to_string(),
            "DIM".to_string(),
            "DID".to_string(),
            "DII".to_string(),
            "IMM".to_string(),
            "IMD".to_string(),
            "IMI".to_string(),
            "IDM".to_string(),
            "IDD".to_string(),
            "IDI".to_string(),
            "IIM".to_string(),
            "IID".to_string(),
            "III".to_string()
        ];
        assert_eq!(states, expected);
    }

    #[test]
    fn test_get_all_kmers0() {
        let kmers: Vec<String> = get_all_kmers(0);
        let expected = vec![
            "".to_string()
        ];
        assert_eq!(kmers, expected);
    }

    #[test]
    fn test_get_all_kmers1() {
        let kmers: Vec<String> = get_all_kmers(1);
        let expected = vec![
            "A".to_string(),
            "C".to_string(),
            "G".to_string(),
            "T".to_string(),
            "-".to_string(),
        ];
        assert_eq!(kmers, expected);
    }


    #[test]
    fn test_get_all_kmers2() {
        let kmers: Vec<String> = get_all_kmers(2);
        let expected = vec![
            "AA".to_string(),
            "AC".to_string(),
            "AG".to_string(),
            "AT".to_string(),
            "A-".to_string(),
            "CA".to_string(),
            "CC".to_string(),
            "CG".to_string(),
            "CT".to_string(),
            "C-".to_string(),
            "GA".to_string(),
            "GC".to_string(),
            "GG".to_string(),
            "GT".to_string(),
            "G-".to_string(),
            "TA".to_string(),
            "TC".to_string(),
            "TG".to_string(),
            "TT".to_string(),
            "T-".to_string(),
            "-A".to_string(),
            "-C".to_string(),
            "-G".to_string(),
            "-T".to_string(),
            "--".to_string(),
        ];
        assert_eq!(kmers, expected);
    }
}
