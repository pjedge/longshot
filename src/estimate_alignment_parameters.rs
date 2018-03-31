use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam::record::Cigar;
use std::error::Error;
use util::*;
use realignment::*;
use bio::io::fasta;
use extract_fragments::{CigarPos, create_augmented_cigarlist};

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub enum AlignmentState {
    Match,
    Insertion,
    Deletion
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct StateTransition {
    pub prev_state: AlignmentState,
    pub current_state: AlignmentState
}

// these counts can be used to help us estimate alignment parameters
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

// these counts can be used to help us estimate alignment parameters
#[derive(Clone, Copy)]
pub struct EmissionCounts {
    equal: usize,
    not_equal: usize
}

#[derive(Clone, Copy)]
struct AlignmentCounts {
    transition_counts: TransitionCounts,
    emission_counts: EmissionCounts
}

impl TransitionCounts {
    fn to_probs(&self) -> TransitionProbs {
        let total_from_match: f64 = self.match_from_match as f64 + self.insertion_from_match as f64 + self.deletion_from_match as f64;
        let total_from_insertion: f64 = self.insertion_from_insertion as f64 + self.match_from_insertion as f64;
        let total_from_deletion: f64 = self.deletion_from_deletion as f64 + self.match_from_deletion as f64;

        TransitionProbs {
            match_from_match: self.match_from_match as f64 / total_from_match,
            insertion_from_match: self.insertion_from_match as f64 / total_from_match,
            deletion_from_match: self.deletion_from_match as f64 / total_from_match,
            insertion_from_insertion: self.insertion_from_insertion as f64 / total_from_insertion,
            match_from_insertion: self.match_from_insertion as f64 / total_from_insertion,
            deletion_from_deletion: self.deletion_from_deletion as f64 / total_from_deletion,
            match_from_deletion: self.match_from_deletion as f64 / total_from_deletion
        }
    }

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
    fn to_probs(&self) -> EmissionProbs {
        let total: f64 = self.equal as f64 + self.not_equal as f64;

        EmissionProbs {
            equal: self.equal as f64 / total,
            not_equal: self.not_equal as f64 / total / 3.0, // 3 possible bases to mismatch to
        }
    }

    fn add(&mut self, other: EmissionCounts) {
        self.equal += other.equal;
        self.not_equal += other.not_equal;
    }
}


impl AlignmentCounts {
    fn to_parameters(&self) -> AlignmentParameters {
        AlignmentParameters {
            transition_probs: self.transition_counts.to_probs(),
            emission_probs: self.emission_counts.to_probs()
        }
    }
}


// for a given read

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
    pub enum CigarError {
        UnsupportedOperation(msg: String) {
            description("Unsupported CIGAR operation")
            display(x) -> ("{}: {}", x.description(), msg)
        }
        UnexpectedOperation(msg: String) {
            description("CIGAR operation not allowed at this point")
            display(x) -> ("{}: {}", x.description(), msg)
        }
    }
}

pub fn count_alignment_events(cigarpos_list: &Vec<CigarPos>,
                              ref_seq: &Vec<char>,
                              read_seq: &Vec<char>)
                              -> Result<(TransitionCounts, EmissionCounts), CigarError>{

    // key is a state transition
    // value is the number of times that transition was observed
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
        not_equal: 0
    };

    let mut state: AlignmentState = AlignmentState::Match;

    for cigarpos in cigarpos_list {

        match cigarpos.cig {
            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {

                let mut ref_pos = cigarpos.ref_pos as usize;
                let mut read_pos = cigarpos.read_pos as usize;

                for _ in 0..l {

                    if ref_pos >= ref_seq.len() - 1 || read_pos >= read_seq.len() - 1 {
                        break;
                    }

                    match state {
                        AlignmentState::Match => {transition_counts.match_from_match += 1;},
                        AlignmentState::Deletion => {transition_counts.match_from_deletion += 1;},
                        AlignmentState::Insertion => {transition_counts.match_from_insertion += 1;},
                    }

                    state = AlignmentState::Match;

                    if ref_seq[ref_pos] != 'N' && read_seq[read_pos] != 'N' {
                        if ref_seq[ref_pos] == read_seq[read_pos] {
                            emission_counts.equal += 1;
                        } else {
                            emission_counts.not_equal += 1;
                        }
                    }

                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Cigar::Ins(l) => {

                let mut ref_pos = cigarpos.ref_pos as usize;
                let mut read_pos = cigarpos.read_pos as usize;

                for _ in 0..l {

                    if ref_pos >= ref_seq.len() - 1 || read_pos >= read_seq.len() - 1 {
                        break;
                    }

                    match state {
                        AlignmentState::Insertion => {transition_counts.insertion_from_insertion += 1;},
                        AlignmentState::Match => {transition_counts.insertion_from_match += 1;},
                        AlignmentState::Deletion => {panic!("Unexpected DELETION alignment state transition observed during INSERTION in bam file.");}
                    }

                    state = AlignmentState::Insertion;

                    read_pos += 1;
                }
            }
            Cigar::Del(l) |
            Cigar::RefSkip(l) => {

                let mut ref_pos = cigarpos.ref_pos as usize;
                let mut read_pos = cigarpos.read_pos as usize;

                for _ in 0..l {

                    if ref_pos >= ref_seq.len() - 1 || read_pos >= read_seq.len() - 1 {
                        break;
                    }

                    match state {
                        AlignmentState::Deletion => {transition_counts.deletion_from_deletion += 1;},
                        AlignmentState::Match => {transition_counts.deletion_from_match += 1;},
                        AlignmentState::Insertion => {panic!("Unexpected INSERTION alignment state transition observed during DELETION in bam file.");}
                    }

                    state = AlignmentState::Deletion;

                    ref_pos += 1;
                }
            }
            Cigar::Pad(_) |
            Cigar::Back(_) |
            Cigar::SoftClip(_) |
            Cigar::HardClip(_) => {
                return Err(CigarError::UnexpectedOperation(
                    "CIGAR operation found in cigarpos_list that should have been removed already.".to_owned()
                ));
            }
        }
    }

    return Ok((transition_counts, emission_counts))
}

//************************************************************************************************
// END OF RUST-HTSLIB BASED CODE *****************************************************************
//************************************************************************************************


pub fn estimate_alignment_parameters(bamfile_name: &String,
                                     fastafile_name: &String,
                                     interval: &Option<GenomicInterval>)
                                     -> AlignmentParameters {

    let t_names = parse_target_names(&bamfile_name);

    let mut prev_tid = 4294967295; // huge value so that tid != prev_tid on first iter
    let mut fasta = fasta::IndexedReader::from_file(fastafile_name).unwrap();
    let mut ref_seq: Vec<char> = vec![];

    // TODO: this uses a lot of duplicate code, need to figure out a better solution.
    // key is a state transition
    // value is the number of times that transition was observed
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
        not_equal: 1
    };

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

                let read_seq: Vec<char> = dna_vec(&record.seq().as_bytes());
                let bam_cig: CigarStringView = record.cigar();
                let cigarpos_list: Vec<CigarPos> =
                    create_augmented_cigarlist(record.pos() as u32, &bam_cig).expect("Error creating augmented cigarlist.");

                let (read_transition_counts, read_emission_counts) =
                    count_alignment_events(&cigarpos_list, &ref_seq, &read_seq).expect("Error counting cigar alignment events.");

                transition_counts.add(read_transition_counts);
                emission_counts.add(read_emission_counts);

                prev_tid = tid;
            }
        }
        &None => {
            let mut bam = bam::Reader::from_path(bamfile_name).unwrap();
            for r in bam.records() {
                let record = r.unwrap();

                let tid: usize = record.tid() as usize;
                let chrom: String = t_names[record.tid() as usize].clone();

                if tid != prev_tid {
                    let mut ref_seq_u8: Vec<u8> = vec![];
                    fasta.read_all(&chrom, &mut ref_seq_u8).expect("Failed to read fasta sequence record.");
                    ref_seq = dna_vec(&ref_seq_u8);
                }

                let read_seq: Vec<char> = dna_vec(&record.seq().as_bytes());
                let bam_cig: CigarStringView = record.cigar();
                let cigarpos_list: Vec<CigarPos> =
                    create_augmented_cigarlist(record.pos() as u32, &bam_cig).expect("Error creating augmented cigarlist.");

                let (read_transition_counts, read_emission_counts) =
                    count_alignment_events(&cigarpos_list, &ref_seq, &read_seq).expect("Error counting cigar alignment events.");

                transition_counts.add(read_transition_counts);
                emission_counts.add(read_emission_counts);

                prev_tid = tid;
            }
        }
    };

    let alignment_counts = AlignmentCounts {
        transition_counts: transition_counts,
        emission_counts: emission_counts
    };

    let params = alignment_counts.to_parameters();

    eprintln!("{} Done estimating alignment parameters.",print_time());
    eprintln!("");

    eprintln!("{} Transition Probabilities:",SPACER);
    eprintln!("{} match -> match:          {:.3}", SPACER, params.transition_probs.match_from_match);
    eprintln!("{} match -> insertion:      {:.3}", SPACER, params.transition_probs.insertion_from_match);
    eprintln!("{} match -> deletion:       {:.3}", SPACER, params.transition_probs.deletion_from_match);
    eprintln!("{} deletion -> match:       {:.3}", SPACER, params.transition_probs.match_from_deletion);
    eprintln!("{} deletion -> deletion:    {:.3}", SPACER, params.transition_probs.deletion_from_deletion);
    eprintln!("{} insertion -> match:      {:.3}", SPACER, params.transition_probs.match_from_insertion);
    eprintln!("{} insertion -> insertion:  {:.3}", SPACER, params.transition_probs.insertion_from_insertion);
    eprintln!("");

    eprintln!("{} Emission Probabilities:", SPACER);
    eprintln!("{} match (equal):           {:.3}", SPACER, params.emission_probs.equal);
    eprintln!("{} match (not equal):       {:.3}", SPACER, params.emission_probs.not_equal);
    eprintln!("");

    params
}
