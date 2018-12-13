//! This module contains a function for estimating the mean read coverage from a BAM file.
//!
//! This is used in Longshot if the user specifies the ```-A``` option, which calculates the mean read
//! coverage and then estimates a maximum read coverage based on this mean coverage.

// extern crates
extern crate rust_htslib;

// use declarations
use errors::*;
use rust_htslib::bam;
use rust_htslib::prelude::*;
use util::{get_interval_lst, print_time, GenomicInterval};

///
pub fn calculate_mean_coverage(
    bam_file: &String,
    interval: &Option<GenomicInterval>,
) -> Result<f64> {
    let bam = bam::Reader::from_path(bam_file).chain_err(|| ErrorKind::BamOpenError)?;

    let mut prev_tid = 4294967295;
    let mut bam_covered_positions = 0;
    let mut total_read_bases = 0;
    let mut total_bam_ref_positions = 0;

    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bam_file, interval)?;

    let mut bam_ix =
        bam::IndexedReader::from_path(bam_file).chain_err(|| ErrorKind::IndexedBamOpenError)?;

    for iv in interval_lst {
        bam_ix
            .fetch(iv.tid, iv.start_pos, iv.end_pos + 1)
            .chain_err(|| ErrorKind::IndexedBamFetchError)?;

        for p in bam_ix.pileup() {
            let pileup = p.chain_err(|| ErrorKind::IndexedBamPileupReadError)?;

            let tid: u32 = pileup.tid();

            if tid != prev_tid {
                total_bam_ref_positions += bam
                    .header()
                    .target_len(tid)
                    .chain_err(|| ErrorKind::BamHeaderTargetLenAccessError)?;
            }

            if tid != iv.tid || pileup.pos() < iv.start_pos || pileup.pos() > iv.end_pos {
                prev_tid = tid;
                continue;
            }

            let mut depth: usize = 0;

            // pileup the bases for a single position and count number of each base
            for alignment in pileup.alignments() {
                let record = alignment.record();

                // may be faster to implement this as bitwise operation on raw flag in the future?
                if record.is_unmapped()
                    || record.is_secondary()
                    || record.is_quality_check_failed()
                    || record.is_duplicate()
                    || record.is_supplementary()
                {
                    continue;
                }

                depth += 1;
            }

            bam_covered_positions += 1;
            total_read_bases += depth;
            prev_tid = tid;
        }
    }

    let total_ref_positions: usize = match interval {
        &Some(ref iv) => (iv.end_pos - iv.start_pos + 1) as usize,
        &None => {
            // output a warning if the number of covered bases is significantly less than the ref positions
            if total_bam_ref_positions / 2 > bam_covered_positions {
                eprintln!("{} WARNING: Max coverage calculation is highly likely to be incorrect. The number of reference \
                              bases covered by the bam file ({}) differs significantly from the expected number of positions in the \
                              reference ({}). If you are using a bam file that only covers part of the genome, please specify \
                              this region exactly with the --region argument so the number of reference bases is known. \
                              Alternatively, disable maximum coverage filtering by setting -C to a large number.",
                          print_time(), bam_covered_positions, total_bam_ref_positions);
            }

            total_bam_ref_positions as usize
        }
    };

    eprintln!(
        "{} Total reference positions: {}",
        print_time(),
        total_ref_positions
    );
    eprintln!("{} Total bases in bam: {}", print_time(), total_read_bases);

    Ok(total_read_bases as f64 / total_ref_positions as f64)
}
