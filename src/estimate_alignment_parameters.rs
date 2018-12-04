use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::CigarStringView;
use rust_htslib::bam::record::Cigar;
use errors::*;
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
            insertion: 1.0,
            deletion: 1.0
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

pub fn count_alignment_events(cigarpos_list: &Vec<CigarPos>,
                              ref_seq: &Vec<char>,
                              read_seq: &Vec<char>,
                              max_cigar_indel: u32)
                              -> Result<(TransitionCounts, EmissionCounts)> {

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

                if l > max_cigar_indel {
                    continue;
                }

                let mut ref_pos = cigarpos.ref_pos as usize;
                let mut read_pos = cigarpos.read_pos as usize;

                for _ in 0..l {

                    if ref_pos >= ref_seq.len() - 1 || read_pos >= read_seq.len() - 1 {
                        break;
                    }

                    match state {
                        AlignmentState::Insertion => {transition_counts.insertion_from_insertion += 1;},
                        AlignmentState::Match => {transition_counts.insertion_from_match += 1;},
                        AlignmentState::Deletion => {
                            // MINIMAP2 sometimes goes directly from insertion <-> deletion
                            // we will just add an implicit deletion -> match -> insertion
                            transition_counts.match_from_deletion += 1;
                            transition_counts.insertion_from_match += 1;
                        }
                    }

                    state = AlignmentState::Insertion;

                    read_pos += 1;
                }
            }
            Cigar::Del(l) |
            Cigar::RefSkip(l) => {

                if l > max_cigar_indel {
                    continue;
                }

                let mut ref_pos = cigarpos.ref_pos as usize;
                let mut read_pos = cigarpos.read_pos as usize;

                for _ in 0..l {

                    if ref_pos >= ref_seq.len() - 1 || read_pos >= read_seq.len() - 1 {
                        break;
                    }

                    match state {
                        AlignmentState::Deletion => {transition_counts.deletion_from_deletion += 1;},
                        AlignmentState::Match => {transition_counts.deletion_from_match += 1;},
                        AlignmentState::Insertion => {
                            // MINIMAP2 sometimes goes directly from insertion <-> deletion
                            // we will just add an implicit Insertion -> match -> deletion
                            transition_counts.match_from_insertion += 1;
                            transition_counts.deletion_from_match += 1;
                        }
                    }

                    state = AlignmentState::Deletion;

                    ref_pos += 1;
                }
            }
            Cigar::Pad(_) |
            Cigar::Back(_) |
            Cigar::SoftClip(_) |
            Cigar::HardClip(_) => {
                bail!("CIGAR operation found in cigarpos_list that should have been removed already.");
            }
        }
    }

    return Ok((transition_counts, emission_counts))
}

//************************************************************************************************
// END OF RUST-HTSLIB BASED CODE *****************************************************************
//************************************************************************************************


pub fn estimate_alignment_parameters(bam_file: &String,
                                     fasta_file: &String,
                                     interval: &Option<GenomicInterval>,
                                     min_mapq: u8,
                                     max_cigar_indel: u32)
                                     -> Result<AlignmentParameters> {

    let t_names = parse_target_names(&bam_file)?;

    let mut prev_tid = 4294967295; // huge value so that tid != prev_tid on first iter
    let mut fasta = fasta::IndexedReader::from_file(fasta_file).chain_err(|| ErrorKind::IndexedFastaOpenError)?;
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

    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bam_file, interval)?;
    let mut bam_ix = bam::IndexedReader::from_path(bam_file).chain_err(|| ErrorKind::IndexedBamOpenError)?;

    for iv in interval_lst {
        bam_ix.fetch(iv.tid, iv.start_pos, iv.end_pos + 1).chain_err(|| "Error seeking BAM file while estimating alignment parameters.")?;

        for r in bam_ix.records() {
            let record = r.chain_err(|| "Error reading BAM record while estimating alignment parameters.")?;

            // may be faster to implement this as bitwise operation on raw flag in the future?
            if record.mapq() < min_mapq || record.is_unmapped() || record.is_secondary() ||
                record.is_quality_check_failed() ||
                record.is_duplicate() || record.is_supplementary() {
                continue;
            }

            let tid: usize = record.tid() as usize;
            let chrom: String = t_names[record.tid() as usize].clone();
            //println!("{}",u8_to_string(record.qname()));
            if tid != prev_tid {
                let mut ref_seq_u8: Vec<u8> = vec![];
                fasta.read_all(&chrom, &mut ref_seq_u8).chain_err(|| "Failed to read fasta sequence record.")?;
                ref_seq = dna_vec(&ref_seq_u8);
            }

            let read_seq: Vec<char> = dna_vec(&record.seq().as_bytes());
            let bam_cig: CigarStringView = record.cigar();
            let cigarpos_list: Vec<CigarPos> =
                create_augmented_cigarlist(record.pos() as u32, &bam_cig).chain_err(|| "Error creating augmented cigarlist.")?;

            let (read_transition_counts, read_emission_counts) =
                count_alignment_events(&cigarpos_list, &ref_seq, &read_seq, max_cigar_indel).chain_err(|| "Error counting cigar alignment events.")?;

            transition_counts.add(read_transition_counts);
            emission_counts.add(read_emission_counts);

            prev_tid = tid;
        }
    }

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
    eprintln!("{} insertion:               {:.3}", SPACER, params.emission_probs.insertion);
    eprintln!("{} deletion:                {:.3}", SPACER, params.emission_probs.deletion);
    eprintln!("");

    Ok(params)
}
