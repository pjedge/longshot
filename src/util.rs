use std::ascii::AsciiExt;
use bio::stats::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use chrono::prelude::*;

pub static INDEX_FREQ: usize = 1000;
pub static MAX_VCF_QUAL: f64 = 500.0;
// THESE SHOULD BE COMMAND LINE PARAMS
pub static MAX_P_MISCALL_F64: f64 = 0.2;
pub static MIN_GQ_FOR_PHASING: f64 = 50.0;

pub fn print_time() -> String {
    Local::now().format("%Y-%m-%d %H:%M:%S").to_string()
}

// use this spacer instead of calling print_time() to have the spaces match up with
// lines that document the time
pub static SPACER: &str = "                   ";

#[derive(Debug, Clone)]
pub struct LogProbTable {
    tab: Vec<Vec<LogProb>>
}

impl LogProbTable {
    pub fn zeros(size: usize) -> LogProbTable {
        LogProbTable { tab: vec![vec![LogProb::ln_zero(); size]; size] }
    }

    pub fn ones(size: usize) -> LogProbTable {
        LogProbTable { tab: vec![vec![LogProb::ln_one(); size]; size] }
    }

    pub fn uniform(size: usize) -> LogProbTable {
        LogProbTable { tab: vec![vec![LogProb::from(Prob(1.0 / ((size*size) as f64))); size]; size] }
    }

    pub fn size(&self) -> usize {
        self.tab.len()
    }

    pub fn max_ix(&self) -> (usize, usize) {
        let mut max_i = 0;
        let mut max_j = 0;

        for i in 0..self.size() {
            for j in 0..self.size() {
                if self.tab[i][j] > self.tab[max_i][max_j] {
                    max_i = i;
                    max_j = j;
                }
            }
        }

        (max_i, max_j)
    }

    pub fn ln_sum_table(&self) -> LogProb {
        let mut total: LogProb = LogProb::ln_zero();
        //for row in table {
        //    total = LogProb::ln_add_exp(total, LogProb::ln_sum_exp(row));
        //}
        for i in 0..self.size() {
            for j in 0..self.size() {
                total = LogProb::ln_add_exp(total, self.tab[i][j]);
            }
        }

        total
    }

    pub fn ln_table_normalize(&mut self) -> LogProbTable {

        let total: LogProb = self.ln_sum_table();
        let mut norm: LogProbTable = LogProbTable::zeros(self.size());

        for i in 0..self.size() {
            for j in 0..self.size() {

                norm.tab[i][j] = self.tab[i][j] - total;
                assert!(norm.tab[i][j] <= LogProb::ln_one());
                assert!(norm.tab[i][j] >= LogProb::ln_zero());

                if norm.tab[i][j].is_nan() {
                    return LogProbTable::uniform(self.size());
                }
            }
        }

        //self.print_phred();

        norm.assert_approx_normalized();
        norm
    }

    pub fn assert_approx_normalized(&self) {

        let total: LogProb = self.ln_sum_table();
        let err: LogProb = if total > LogProb::ln_one() {
            LogProb::ln_sub_exp(total, LogProb::ln_one())
        } else {
            LogProb::ln_sub_exp(LogProb::ln_one(), total)
        };

        let margin = 0.00001;
        if err > LogProb::from(Prob(margin)) {
            println!("ERROR");
            println!("{:?}", self.tab);
            println!("{}",*Prob::from(total))
        }
        assert!(err < LogProb::from(Prob(margin)));
    }

    pub fn plus_equals(&mut self, i: usize, j: usize, p: LogProb) {
        self.tab[i][j] = self.tab[i][j] + p;
    }

    pub fn set(&mut self, i: usize, j: usize, p: LogProb) {
        self.tab[i][j] = p;
    }

    pub fn get(&self, i: usize, j: usize) -> LogProb {
        self.tab[i][j]
    }

    pub fn print_phred(&self) {
        for i in 0..self.size() {
            for j in 0..self.size() {
                print!("{} ", *PHREDProb::from(self.tab[i][j]));
            }
            println!("");
        }
        println!("------------------------");
    }

    pub fn print_prob(&self) {
        for i in 0..self.size() {
            for j in 0..self.size() {
                print!("{} ", *Prob::from(self.tab[i][j]));
            }
            println!("");
        }
        println!("------------------------");
    }
}

// this is really ugly. TODO a less verbose implementation
pub fn parse_region_string(region_string: Option<&str>,
                           bamfile_name: &String)
                           -> Option<GenomicInterval> {
    let bam = bam::Reader::from_path(bamfile_name).unwrap();

    match region_string {
        Some(r) if r.contains(":") && r.contains("-") => {
            let split1: Vec<&str> = r.split(":").collect();
            if split1.len() != 2 {
                panic!("Invalid format for region. Please use <chrom> or <chrom:start-stop>");
            }
            let split2: Vec<&str> = split1[1].split("-").collect();
            if split2.len() != 2 {
                panic!("Invalid format for region. Please use <chrom> or <chrom:start-stop>");
            }
            let iv_chrom = split1[0].to_string();
            let iv_start = split2[0].parse::<u32>().expect("Invalid position value specified in region string.");
            let iv_end = split2[1].parse::<u32>().expect("Invalid position value specified in region string.");

            let mut tid: u32 = 0;
            for name in bam.header().target_names() {
                if u8_to_string(name) == iv_chrom {
                    break;
                }
                tid += 1;
            }
            if tid as usize == bam.header().target_names().len() {
                panic!("Chromosome name for region is not in BAM file.");
            }

            Some(GenomicInterval {
                chrom: iv_chrom,
                start_pos: iv_start,
                end_pos: iv_end - 1,
            })
        }
        Some(r) => {
            let r_str = r.to_string();
            let mut tid: u32 = 0;
            for name in bam.header().target_names() {
                if u8_to_string(name) == r_str {
                    break;
                }
                tid += 1;
            }
            if tid as usize == bam.header().target_names().len() {
                panic!("Chromosome name for region is not in BAM file.");
            }

            let tlen = bam.header().target_len(tid).unwrap();
            Some(GenomicInterval {
                chrom: r_str,
                start_pos: 0,
                end_pos: tlen - 1,
            })
        }
        None => None,
    }
}

#[derive(Clone)]
pub struct GenomicInterval {
    pub chrom: String,
    // chromosome name
    pub start_pos: u32,
    // start of interval
    pub end_pos: u32,
    // end of interval (inclusive)
}


pub fn u8_to_string(u: &[u8]) -> String {
    String::from_utf8(u.to_vec()).unwrap()
}

//
pub fn dna_vec(u: &[u8]) -> (Vec<char>) {
    let mut v: Vec<char> = Vec::with_capacity(u.len());
    for cu in AsciiExt::to_ascii_uppercase(u) {
        let c = cu as char;
        //assert!(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
        if c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' {
            v.push(c);
        } else {
            eprintln!("Warning: Unexpected base \"{}\" encountered. Replaced with \"N\".",
                      c);
            v.push('N');
        }
    }
    v
}

pub fn parse_target_names(bam_file: &String) -> Vec<String> {
    let bam = bam::Reader::from_path(bam_file).unwrap();
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

    target_names
}


