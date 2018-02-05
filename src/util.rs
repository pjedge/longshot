use std::ascii::AsciiExt;
use std::collections::HashMap;
use bio::stats::{LogProb, Prob};
use rust_htslib::bam;
use rust_htslib::bam::Read;

static INDEX_FREQ: usize = 1000;

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

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum AlignmentType {
    FastAllAlignment,
    NumericallyStableAllAlignment,
    MaxAlignment
}

#[derive(Clone, Copy)]
pub struct ExtractFragmentParameters {
    pub min_mapq: u8,
    pub alignment_type: AlignmentType,
    pub band_width: usize,
    pub simple_anchors: bool,
    pub anchor_length: usize,
    pub anchor_k: usize,
    pub short_hap_snv_distance: usize,
    pub short_hap_max_snvs: usize,
    pub max_window_length: usize,
}

#[derive(Clone, Copy)]
pub struct AlignmentParameters {
    pub match_from_match: f64,
    pub mismatch_from_match: f64,
    pub insertion_from_match: f64,
    pub deletion_from_match: f64,
    pub extend_from_insertion: f64,
    pub match_from_insertion: f64,
    pub mismatch_from_insertion: f64,
    pub extend_from_deletion: f64,
    pub match_from_deletion: f64,
    pub mismatch_from_deletion: f64,
}

#[derive(Clone, Copy)]
pub struct LnAlignmentParameters {
    pub match_from_match: LogProb,
    pub mismatch_from_match: LogProb,
    pub insertion_from_match: LogProb,
    pub deletion_from_match: LogProb,
    pub extend_from_insertion: LogProb,
    pub match_from_insertion: LogProb,
    pub mismatch_from_insertion: LogProb,
    pub extend_from_deletion: LogProb,
    pub match_from_deletion: LogProb,
    pub mismatch_from_deletion: LogProb,
}

impl AlignmentParameters {
    pub fn ln(&self) -> LnAlignmentParameters {
        LnAlignmentParameters {
            match_from_match: LogProb::from(Prob(self.match_from_match)),
            mismatch_from_match: LogProb::from(Prob(self.mismatch_from_match)),
            insertion_from_match: LogProb::from(Prob(self.insertion_from_match)),
            deletion_from_match: LogProb::from(Prob(self.deletion_from_match)),
            extend_from_insertion: LogProb::from(Prob(self.extend_from_insertion)),
            match_from_insertion: LogProb::from(Prob(self.match_from_insertion)),
            mismatch_from_insertion: LogProb::from(Prob(self.mismatch_from_insertion)),
            extend_from_deletion: LogProb::from(Prob(self.extend_from_deletion)),
            match_from_deletion: LogProb::from(Prob(self.match_from_deletion)),
            mismatch_from_deletion: LogProb::from(Prob(self.mismatch_from_deletion)),
        }
    }
}

#[derive(Clone, Copy)]
pub struct AlignmentCounts {
    pub match_from_match: usize,
    pub mismatch_from_match: usize,
    pub insertion_from_match: usize,
    pub deletion_from_match: usize,
    pub extend_from_insertion: usize,
    pub match_from_insertion: usize,
    pub mismatch_from_insertion: usize,
    pub extend_from_deletion: usize,
    pub match_from_deletion: usize,
    pub mismatch_from_deletion: usize,
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

#[derive(Clone, Copy)]
pub struct FragCall {
    pub frag_ix: Option<usize>, // index into variant list
    pub var_ix: usize, // index into variant list
    pub allele: char, // allele call
    pub qual: LogProb, // LogProb probability the call is an error
    pub one_minus_qual: LogProb, // LogProb probability the call is correct
}

#[derive(Clone)]
pub struct Fragment {
    pub id: String,
    pub calls: Vec<FragCall>,
}


#[derive(Debug, Clone)]
pub struct Var {
    pub ix: usize,
    // index of this variant in the global var list
    pub chrom: String,
    pub pos0: usize,
    pub ref_allele: String,
    pub var_allele: String,
    pub dp: usize,
    // depth of coverage
    pub ra: usize,
    pub aa: usize,
    pub qual: f64,
    pub filter: String,
    pub genotype: String,
    pub gq: f64,
    pub genotype_post: [LogProb; 4],  // genotype posteriors... [p00, p01, p10, p11]
    //pub pileup: Option(Vec<PileupElement>),
}

#[derive(Debug, Clone)]
pub struct VarList {
    pub lst: Vec<Var>,
    ix: HashMap<String, Vec<usize>>,
}

impl VarList {
    pub fn new(lst: Vec<Var>) -> VarList {
        let mut v = VarList {
            lst: lst,
            ix: HashMap::new(),
        };
        for i in 0..v.lst.len() {
            v.lst[i].ix = i;
        }
        v.index_lst();
        v
    }
    fn index_lst(&mut self) {
        // need to throw error if list isn't sorted

        // for every chromosome, get the position of the last variant
        // and the varlist index of the first variant on that chromosome
        let mut max_positions: HashMap<String, usize> = HashMap::new();
        let mut first_ix: HashMap<String, usize> = HashMap::new();
        for (i, var) in self.lst.iter().enumerate() {
            let mpos: usize = *(max_positions.entry(var.chrom.clone()).or_insert(0));
            if var.pos0 > mpos {
                max_positions.insert(var.chrom.clone(), var.pos0);
            }
            if !first_ix.contains_key(&var.chrom) {
                first_ix.insert(var.chrom.clone(), i); // insert varlist index of first variant on chrom
            }
        }

        // for every chrom, iterate up to its last potential variant and create an index that
        // returns the lst index of the next SNV for some position mod 1000
        for (chrom, max_pos) in &max_positions {
            let mut v: Vec<usize> = vec![];
            // the first variant after position 0 is the first variant on the chromosome
            let e = match first_ix.get(chrom) {
                Some(x) => *x,
                None => {
                    panic!("This dictionary is missing a chromosome!");
                }
            };
            v.push(e);
            // for every subsequent position (1000,2000,3000...) we start at the last position
            // and step forward until we get the first variant at or after that position
            for i in 1..(max_pos / INDEX_FREQ + 1) {
                let index_pos = INDEX_FREQ * i;
                let start = v[i - 1];

                let mut found_var = false;
                for j in start..self.lst.len() {
                    if self.lst[j].pos0 >= index_pos {
                        v.push(j);
                        found_var = true;
                        break;
                    }
                }
                if !found_var {
                    break;
                }
            }
            self.ix.insert(chrom.clone(), v); // insert the index vector into the ix dictionary
        }
    }

    pub fn get_variants_range(&self, interval: GenomicInterval) -> (Vec<Var>) {
        // vector of variants to fill and return
        let mut vlst: Vec<Var> = vec![];
        // get the varlist index of a nearby position on the left

        let index_pos = (interval.start_pos as usize) / INDEX_FREQ;

        if index_pos >=
            self.ix
                .get(&interval.chrom)
                .unwrap()
                .len() {
            return vlst;
        }

        let mut i = self.ix.get(&interval.chrom).unwrap()[index_pos];

        while i < self.lst.len() && self.lst[i].pos0 <= interval.end_pos as usize {
            if self.lst[i].pos0 >= interval.start_pos as usize {
                vlst.push(self.lst[i].clone());
            }
            i += 1;
        }
        vlst
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pos_alleles_eq(vlst1: Vec<Var>, vlst2: Vec<Var>) -> bool {
        for (v1, v2) in vlst1.iter().zip(vlst2.iter()){
            if v1.ix != v2.ix || v1.chrom != v2.chrom || v1.pos0 != v2.pos0 ||
                v1.ref_allele != v2.ref_allele || v1.var_allele != v2.var_allele {
                return false;
            }
        }
        true
    }

    fn generate_test_lst1() -> VarList {
        let mut lst: Vec<Var> = vec![];
        lst.push(Var {
            ix: 0,
            chrom: "chr1".to_string(),
            pos0: 5,
            ref_allele: "A".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 1,
            chrom: "chr1".to_string(),
            pos0: 1000,
            ref_allele: "T".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 2,
            chrom: "chr1".to_string(),
            pos0: 2005,
            ref_allele: "T".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 3,
            chrom: "chr1".to_string(),
            pos0: 2900,
            ref_allele: "C".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 4,
            chrom: "chr1".to_string(),
            pos0: 6000,
            ref_allele: "C".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 5,
            chrom: "chr1".to_string(),
            pos0: 10000,
            ref_allele: "C".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 6,
            chrom: "chr2".to_string(),
            pos0: 5,
            ref_allele: "A".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 7,
            chrom: "chr2".to_string(),
            pos0: 1000,
            ref_allele: "T".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 8,
            chrom: "chr2".to_string(),
            pos0: 2005,
            ref_allele: "T".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 9,
            chrom: "chr2".to_string(),
            pos0: 2900,
            ref_allele: "C".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 10,
            chrom: "chr2".to_string(),
            pos0: 6000,
            ref_allele: "C".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 11,
            chrom: "chr2".to_string(),
            pos0: 10000,
            ref_allele: "C".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 12,
            chrom: "chr3".to_string(),
            pos0: 20200,
            ref_allele: "C".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 13,
            chrom: "chr3".to_string(),
            pos0: 25100,
            ref_allele: "A".to_string(),
            var_allele: "C".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        lst.push(Var {
            ix: 14,
            chrom: "chr3".to_string(),
            pos0: 30400,
            ref_allele: "C".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        VarList::new(lst)
    }

    #[test]
    fn test_varlist_get_variants_range1() {
        let vlst = generate_test_lst1();
        let c = "chr1".to_string();
        let p1 = 2500;
        let p2 = 8000;
        let interval = GenomicInterval {
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<Var> = vec![];
        exp.push(Var {
            ix: 3,
            chrom: "chr1".to_string(),
            pos0: 2900,
            ref_allele: "C".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        exp.push(Var {
            ix: 4,
            chrom: "chr1".to_string(),
            pos0: 6000,
            ref_allele: "C".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range2() {
        let vlst = generate_test_lst1();
        let c = "chr2".to_string();
        let p1 = 0;
        let p2 = 3000;
        let interval = GenomicInterval {
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<Var> = vec![];
        exp.push(Var {
            ix: 6,
            chrom: "chr2".to_string(),
            pos0: 5,
            ref_allele: "A".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        exp.push(Var {
            ix: 7,
            chrom: "chr2".to_string(),
            pos0: 1000,
            ref_allele: "T".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        exp.push(Var {
            ix: 8,
            chrom: "chr2".to_string(),
            pos0: 2005,
            ref_allele: "T".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        exp.push(Var {
            ix: 9,
            chrom: "chr2".to_string(),
            pos0: 2900,
            ref_allele: "C".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range3() {
        let vlst = generate_test_lst1();
        let c = "chr2".to_string();
        let p1 = 6000;
        let p2 = 10000;
        let interval = GenomicInterval {
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<Var> = vec![];
        exp.push(Var {
            ix: 10,
            chrom: "chr2".to_string(),
            pos0: 6000,
            ref_allele: "C".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        exp.push(Var {
            ix: 11,
            chrom: "chr2".to_string(),
            pos0: 10000,
            ref_allele: "C".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range4() {
        let vlst = generate_test_lst1();
        let c = "chr3".to_string();
        let p1 = 20100;
        let p2 = 20200;
        let interval = GenomicInterval {
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<Var> = vec![];
        exp.push(Var {
            ix: 12,
            chrom: "chr3".to_string(),
            pos0: 20200,
            ref_allele: "C".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range5() {
        let vlst = generate_test_lst1();
        let c = "chr3".to_string();
        let p1 = 20200;
        let p2 = 20200;
        let interval = GenomicInterval {
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<Var> = vec![];
        exp.push(Var {
            ix: 12,
            chrom: "chr3".to_string(),
            pos0: 20200,
            ref_allele: "C".to_string(),
            var_allele: "G".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range6() {
        let vlst = generate_test_lst1();
        let c = "chr3".to_string();
        let p1 = 25000;
        let p2 = 30500;
        let interval = GenomicInterval {
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<Var> = vec![];
        exp.push(Var {
            ix: 13,
            chrom: "chr3".to_string(),
            pos0: 25100,
            ref_allele: "A".to_string(),
            var_allele: "C".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });
        exp.push(Var {
            ix: 14,
            chrom: "chr3".to_string(),
            pos0: 30400,
            ref_allele: "C".to_string(),
            var_allele: "A".to_string(),
            dp: 40,
            ra: 0,
            aa: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: "./.".to_string(),
            gq: 0.0
        });

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range7() {
        let vlst = generate_test_lst1();
        let c = "chr3".to_string();
        let p1 = 100000;
        let p2 = 200000;
        let interval = GenomicInterval {
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let exp: Vec<Var> = vec![];
        assert!(pos_alleles_eq(vars, exp));
    }
}
