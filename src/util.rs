use bio::stats::{LogProb, Prob};
use chrono::prelude::*;
use clap::ArgMatches;
use errors::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;

pub static INDEX_FREQ: usize = 1000;
pub static MAX_VCF_QUAL: f64 = 500.0;

pub fn print_time() -> String {
    Local::now().format("%Y-%m-%d %H:%M:%S").to_string()
}

// use this spacer instead of calling print_time() to have the spaces match up with
// lines that document the time
pub static SPACER: &str = "                   ";

pub fn parse_u8(argmatch: &ArgMatches, arg_name: &str) -> Result<u8> {
    let parse_result: u8 = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<u8>()
        .chain_err(|| format!("{} must be an integer between 0 and 255!", arg_name))?;
    Ok(parse_result)
}

pub fn parse_u32(argmatch: &ArgMatches, arg_name: &str) -> Result<u32> {
    let parse_result: u32 = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<u32>()
        .chain_err(|| format!("{} must be a positive integer!", arg_name))?;
    Ok(parse_result)
}

pub fn parse_usize(argmatch: &ArgMatches, arg_name: &str) -> Result<usize> {
    let parse_result: usize = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<usize>()
        .chain_err(|| format!("{} must be a positive integer!", arg_name))?;
    Ok(parse_result)
}

pub fn parse_positive_f64(argmatch: &ArgMatches, arg_name: &str) -> Result<f64> {
    let parse_result: f64 = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<f64>()
        .chain_err(|| format!("{} must be a positive float!", arg_name))?;
    ensure!(
        parse_result > 0.0,
        format!("{} must be a positive float!", arg_name)
    );
    Ok(parse_result)
}

pub fn parse_prob_into_logprob(argmatch: &ArgMatches, arg_name: &str) -> Result<LogProb> {
    let parse_result: f64 = argmatch
        .value_of(arg_name)
        .chain_err(|| format!("{} not defined.", arg_name))?
        .parse::<f64>()
        .chain_err(|| format!("{} must be a float between 0.0 and 1.0!", arg_name))?;
    ensure!(
        parse_result >= 0.0 && parse_result <= 1.0,
        format!("{} must be a float between 0.0 and 1.0!", arg_name)
    );
    Ok(LogProb::from(Prob(parse_result)))
}

pub fn parse_flag(argmatch: &ArgMatches, arg_name: &str) -> Result<bool> {
    let is_flag_set = match argmatch.occurrences_of(arg_name) {
        0 => false,
        1 => true,
        _ => {
            bail!(format!("{} cannot be specified multiple times.", arg_name));
        }
    };
    Ok(is_flag_set)
}

// this is really ugly. TODO a less verbose implementation
pub fn parse_region_string(
    region_string: Option<&str>,
    bamfile_name: &String,
) -> Result<Option<GenomicInterval>> {
    let bam = bam::Reader::from_path(bamfile_name).chain_err(|| ErrorKind::BamOpenError)?;

    match region_string {
        Some(r) if r.contains(":") && r.contains("-") => {
            let split1: Vec<&str> = r.split(":").collect();
            if split1.len() != 2 {
                bail!("Invalid format for region. Please use <chrom> or <chrom:start-stop>");
            }
            let split2: Vec<&str> = split1[1].split("-").collect();
            if split2.len() != 2 {
                bail!("Invalid format for region. Please use <chrom> or <chrom:start-stop>");
            }
            let iv_chrom = split1[0].to_string();
            let iv_start = split2[0]
                .parse::<u32>()
                .chain_err(|| "Invalid position value specified in region string.")?; // read in as 1-based inclusive range
            let iv_end = split2[1]
                .parse::<u32>()
                .chain_err(|| "Invalid position value specified in region string.")?; // read in as 1-based inclusive range

            let mut tid: u32 = 0;
            for name in bam.header().target_names() {
                if u8_to_string(name)? == iv_chrom {
                    break;
                }
                tid += 1;
            }
            if tid as usize == bam.header().target_names().len() {
                bail!("Chromosome name for region is not in BAM file.");
            }

            let tlen = bam
                .header()
                .target_len(tid)
                .chain_err(|| ErrorKind::BamHeaderTargetLenAccessError)?;

            ensure!(
                iv_start > 0,
                "--region start position must be greater than 0."
            );
            ensure!(
                iv_start < tlen,
                "--region start position exceeds the length of the contig."
            );
            ensure!(
                iv_end <= tlen,
                "--region end position exceeds the length of the contig."
            );

            Ok(Some(GenomicInterval {
                tid: tid,
                chrom: iv_chrom,
                start_pos: iv_start - 1, // convert to 0-based inclusive range
                end_pos: iv_end - 1,     // convert to 0-based inclusive range
            }))
        }
        Some(r) => {
            let r_str = r.to_string();
            let mut tid: u32 = 0;
            for name in bam.header().target_names() {
                if u8_to_string(name)? == r_str {
                    break;
                }
                tid += 1;
            }
            if tid as usize == bam.header().target_names().len() {
                bail!("Chromosome name for region is not in BAM file.");
            }

            let tlen = bam
                .header()
                .target_len(tid)
                .chain_err(|| ErrorKind::BamHeaderTargetLenAccessError)?;
            Ok(Some(GenomicInterval {
                tid: tid,
                chrom: r_str,
                start_pos: 0,
                end_pos: tlen - 1,
            }))
        }
        None => Ok(None),
    }
}

#[derive(Clone)]
pub struct GenomicInterval {
    pub tid: u32,
    pub chrom: String,
    // chromosome name
    pub start_pos: u32,
    // start of interval (0-indexed)
    pub end_pos: u32,
    // end of interval (0-indexed, inclusive)
}

#[derive(Clone)]
pub struct DensityParameters {
    pub n: usize,
    pub len: usize,
    pub gq: f64,
}

pub fn u8_to_string(u: &[u8]) -> Result<String> {
    Ok(String::from_utf8(u.to_vec()).chain_err(|| "Error converting u8 to String.")?)
}

//
pub fn dna_vec(u: &[u8]) -> (Vec<char>) {
    let mut v: Vec<char> = Vec::with_capacity(u.len());
    for cu in u.to_ascii_uppercase() {
        let c = cu as char;
        //assert!(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
        if c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' {
            v.push(c);
        } else {
            eprintln!(
                "Warning: Unexpected base \"{}\" encountered. Replaced with \"N\".",
                c
            );
            v.push('N');
        }
    }
    v
}

pub fn has_non_acgt(s: &String) -> bool {
    for c in s.chars() {
        if !(c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            return true;
        }
    }
    false
}

pub fn parse_target_names(bam_file: &String) -> Result<Vec<String>> {
    let bam = bam::Reader::from_path(bam_file).chain_err(|| ErrorKind::BamOpenError)?;
    let header_view = bam.header();
    let target_names_dec: Vec<&[u8]> = header_view.target_names();
    let mut target_names: Vec<String> = vec![];

    for t_name_dec in target_names_dec {
        let mut name_vec: Vec<char> = vec![];
        for decr in t_name_dec {
            let dec: u8 = *decr;
            name_vec.push(dec as char);
        }
        let name_string: String = name_vec.into_iter().collect();
        target_names.push(name_string);
    }

    Ok(target_names)
}

pub fn get_whole_genome_intervals(bam_file: &String) -> Result<Vec<GenomicInterval>> {
    let bam = bam::Reader::from_path(bam_file).chain_err(|| ErrorKind::BamOpenError)?;
    let header_view = bam.header();
    let target_names_dec: Vec<&[u8]> = header_view.target_names();
    let mut intervals: Vec<GenomicInterval> = vec![];

    for (tid, t_name_dec) in target_names_dec.iter().enumerate() {
        let mut name_vec: Vec<char> = vec![];
        for decr in t_name_dec.iter() {
            let dec: u8 = *decr;
            name_vec.push(dec as char);
        }
        let name_string: String = name_vec.into_iter().collect();
        intervals.push(GenomicInterval {
            tid: tid as u32,
            chrom: name_string,
            start_pos: 0,
            end_pos: header_view
                .target_len(tid as u32)
                .chain_err(|| format!("Error accessing target len for tid {}", tid))?
                - 1,
        });
    }

    Ok(intervals)
}

// given a bam file name and a possible genomic interval,
// if the interval exists then just return a vector holding that lone interval
// otherwise, if the interval is None,
// return a vector holding GenomicIntervals representing the whole genome.
pub fn get_interval_lst(
    bam_file: &String,
    interval: &Option<GenomicInterval>,
) -> Result<Vec<GenomicInterval>> {
    match interval {
        &Some(ref iv) => Ok(vec![iv.clone()]),
        &None => get_whole_genome_intervals(bam_file),
    }
}
