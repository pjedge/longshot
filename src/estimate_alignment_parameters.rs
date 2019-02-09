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

/// represents the 3 states for the sequence alignment Pair-HMM
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub enum AlignmentState {
    Match,
    Insertion,
    Deletion,
}

/// represents a transition between two states in the sequence alignment Pair-HMM
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct StateTransition {
    pub prev_state: AlignmentState,
    pub current_state: AlignmentState,
}

/// represents counts of different transitions from state-to-state in the Pair-HMM
///
/// the counts in this struct are used to keep track of the different alignment transitions that
/// are observed in the BAM file. These counts can be used to help us directly estimate alignment parameters
#[derive(Clone, Copy)]
pub struct TransitionCounts {
    match_from_match: usize,
    insertion_from_match: usize,
    deletion_from_match: usize,
    insertion_from_insertion: usize,
    match_from_insertion: usize,
    deletion_from_deletion: usize,
    match_from_deletion: usize,
}

/// represents counts of different emission events (aligned bases being equal or not equal) in the BAM
#[derive(Clone, Copy)]
pub struct EmissionCounts {
    equal: usize,
    not_equal: usize,
}

/// a struct to hold counts of transition events and emission events from the BAM file
///
/// together these counts represent all of the counts necessary to estimate the HMM parameters
#[derive(Clone, Copy)]
struct AlignmentCounts {
    transition_counts: TransitionCounts,
    emission_counts: EmissionCounts,
}

impl TransitionCounts {
    /// converts TransitionCounts to TransitionProbs
    ///
    /// Using this function, HMM transition counts that are empirically observed and counted in the
    /// BAM file can be converted into probabilities that can be used for sequence alignment.
    /// The transition counts out of a single state (e.g. insertion -> other state) are converted
    /// into a fraction of the total transition counts out of that state.
    fn to_probs(&self) -> TransitionProbs {
        let total_from_match: f64 = self.match_from_match as f64
            + self.insertion_from_match as f64
            + self.deletion_from_match as f64;
        let total_from_insertion: f64 =
            self.insertion_from_insertion as f64 + self.match_from_insertion as f64;
        let total_from_deletion: f64 =
            self.deletion_from_deletion as f64 + self.match_from_deletion as f64;

        TransitionProbs {
            match_from_match: self.match_from_match as f64 / total_from_match,
            insertion_from_match: self.insertion_from_match as f64 / total_from_match,
            deletion_from_match: self.deletion_from_match as f64 / total_from_match,
            insertion_from_insertion: self.insertion_from_insertion as f64 / total_from_insertion,
            match_from_insertion: self.match_from_insertion as f64 / total_from_insertion,
            deletion_from_deletion: self.deletion_from_deletion as f64 / total_from_deletion,
            match_from_deletion: self.match_from_deletion as f64 / total_from_deletion,
        }
    }

    /// add the corresponding counts inside two ```TransitionCount```s together
    ///
    /// this function allows us to create separate ```TransitionCount``` structs for each newly
    /// observed read, and then just add those counts onto a grand total count for the BAM file.
    fn add(&mut self, other: TransitionCounts) {
        self.match_from_match += other.match_from_match;
        self.insertion_from_match += other.insertion_from_match;
        self.deletion_from_match += other.deletion_from_match;
        self.insertion_from_insertion += other.insertion_from_insertion;
        self.match_from_insertion += other.match_from_insertion;
        self.deletion_from_deletion += other.deletion_from_deletion;
        self.match_from_deletion += other.match_from_deletion;
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
        let total: f64 = self.equal as f64 + self.not_equal as f64;

        EmissionProbs {
            equal: self.equal as f64 / total,
            not_equal: self.not_equal as f64 / total / 3.0, // 3 possible bases to mismatch to
            insertion: 1.0,
            deletion: 1.0,
        }
    }

    /// add the corresponding counts inside two ```EmissionCount```s together
    ///
    /// this function allows us to create separate ```EmissionCount``` structs for each newly
    /// observed read, and then just add those counts onto a grand total count for the BAM file.
    fn add(&mut self, other: EmissionCounts) {
        self.equal += other.equal;
        self.not_equal += other.not_equal;
    }
}

impl AlignmentCounts {
    /// convert ```AlignmentCounts``` into ```AlignmentParameters``` by converting the transition
    /// and emission counts into probabilities
    fn to_parameters(&self) -> AlignmentParameters {
        AlignmentParameters {
            transition_probs: self.transition_counts.to_probs(),
            emission_probs: self.emission_counts.to_probs(),
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
    max_cigar_indel: u32,
) -> Result<(TransitionCounts, EmissionCounts)> {

    // initialize TransitionCounts and EmissionCounts with all counts set to 0
    let mut transition_counts = TransitionCounts {
        match_from_match: 0,
        insertion_from_match: 0,
        deletion_from_match: 0,
        insertion_from_insertion: 0,
        match_from_insertion: 0,
        deletion_from_deletion: 0,
        match_from_deletion: 0,
    };

    let mut emission_counts = EmissionCounts {
        equal: 0,
        not_equal: 0,
    };

    // assume the initial state is MATCH
    let mut state: AlignmentState = AlignmentState::Match;

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

                    // we add the transition from the current state
                    // to the new state (which is match)
                    match state {
                        AlignmentState::Match => {
                            transition_counts.match_from_match += 1;
                        }
                        AlignmentState::Deletion => {
                            transition_counts.match_from_deletion += 1;
                        }
                        AlignmentState::Insertion => {
                            transition_counts.match_from_insertion += 1;
                        }
                    }

                    // we have transitioned to a match so set the current state to match
                    state = AlignmentState::Match;

                    // check if the aligned bases match or mismatch and use these to iterate the
                    // emission counts (whether bases match or mismatch)
                    if ref_seq[ref_pos] != 'N' && read_seq[read_pos] != 'N' {
                        if ref_seq[ref_pos] == read_seq[read_pos] {
                            emission_counts.equal += 1;
                        } else {
                            emission_counts.not_equal += 1;
                        }
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

                    // we add the transition from the current state
                    // to the new state (which is insertion)
                    match state {
                        AlignmentState::Insertion => {
                            transition_counts.insertion_from_insertion += 1;
                        }
                        AlignmentState::Match => {
                            transition_counts.insertion_from_match += 1;
                        }
                        AlignmentState::Deletion => {
                            // MINIMAP2 sometimes goes directly from insertion <-> deletion
                            // we will just add an implicit deletion -> match -> insertion
                            transition_counts.match_from_deletion += 1;
                            transition_counts.insertion_from_match += 1;
                        }
                    }

                    // we have transitioned to an insertion so set the current state as insertion
                    state = AlignmentState::Insertion;
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
                    // we add the transition from the current state
                    // to the new state (which is deletion)
                    match state {
                        AlignmentState::Deletion => {
                            transition_counts.deletion_from_deletion += 1;
                        }
                        AlignmentState::Match => {
                            transition_counts.deletion_from_match += 1;
                        }
                        AlignmentState::Insertion => {
                            // MINIMAP2 sometimes goes directly from insertion <-> deletion
                            // we will just add an implicit Insertion -> match -> deletion
                            transition_counts.match_from_insertion += 1;
                            transition_counts.deletion_from_match += 1;
                        }
                    }

                    // we have transitioned to deletion so set the current state to deletion
                    state = AlignmentState::Deletion;
                    // this is a deletion so only the reference position is moved forward
                    ref_pos += 1;
                }
            }
            Cigar::Pad(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) => {
                bail!(
                    "CIGAR operation found in cigarpos_list that should have been removed already."
                );
            }
        }
    }

    return Ok((transition_counts, emission_counts));
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
    let mut transition_counts = TransitionCounts {
        match_from_match: 1,
        insertion_from_match: 1,
        deletion_from_match: 1,
        insertion_from_insertion: 1,
        match_from_insertion: 1,
        deletion_from_deletion: 1,
        match_from_deletion: 1,
    };

    let mut emission_counts = EmissionCounts {
        equal: 1,
        not_equal: 1,
    };

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
                    .fetch_all(&chrom)
                    .chain_err(|| ErrorKind::IndexedFastaReadError)?;
                fasta
                    .read(&mut ref_seq_u8)
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
                count_alignment_events(&cigarpos_list, &ref_seq, &read_seq, max_cigar_indel)
                    .chain_err(|| "Error counting cigar alignment events.")?;

            // add emission and transition counts to the running total
            transition_counts.add(read_transition_counts);
            emission_counts.add(read_emission_counts);

            prev_tid = tid;
        }
    }

    // place the transition and emission counts together in an AlignmentCounts struct
    let alignment_counts = AlignmentCounts {
        transition_counts: transition_counts,
        emission_counts: emission_counts,
    };

    // convert the alignment counts from the BAM into probabilities
    let params = alignment_counts.to_parameters();

    // print the estimated alignment parameters to STDERR
    eprintln!("{} Done estimating alignment parameters.", print_time());
    eprintln!("");

    eprintln!("{} Transition Probabilities:", SPACER);
    eprintln!(
        "{} match -> match:          {:.3}",
        SPACER, params.transition_probs.match_from_match
    );
    eprintln!(
        "{} match -> insertion:      {:.3}",
        SPACER, params.transition_probs.insertion_from_match
    );
    eprintln!(
        "{} match -> deletion:       {:.3}",
        SPACER, params.transition_probs.deletion_from_match
    );
    eprintln!(
        "{} deletion -> match:       {:.3}",
        SPACER, params.transition_probs.match_from_deletion
    );
    eprintln!(
        "{} deletion -> deletion:    {:.3}",
        SPACER, params.transition_probs.deletion_from_deletion
    );
    eprintln!(
        "{} insertion -> match:      {:.3}",
        SPACER, params.transition_probs.match_from_insertion
    );
    eprintln!(
        "{} insertion -> insertion:  {:.3}",
        SPACER, params.transition_probs.insertion_from_insertion
    );
    eprintln!("");

    eprintln!("{} Emission Probabilities:", SPACER);
    eprintln!(
        "{} match (equal):           {:.3}",
        SPACER, params.emission_probs.equal
    );
    eprintln!(
        "{} match (not equal):       {:.3}",
        SPACER, params.emission_probs.not_equal
    );
    eprintln!(
        "{} insertion:               {:.3}",
        SPACER, params.emission_probs.insertion
    );
    eprintln!(
        "{} deletion:                {:.3}",
        SPACER, params.emission_probs.deletion
    );
    eprintln!("");

    Ok(params)
}
