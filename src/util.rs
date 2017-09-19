
use std::ascii::AsciiExt;
use std::collections::HashMap;
use bio::stats::{LogProb, Prob};
use rust_htslib::bam;
use rust_htslib::bam::Read;

static INDEX_FREQ: usize = 1000;

pub struct GenomicInterval {
    pub tid: Option<u32>, // target ID for chromosome in bam
    pub chrom: Option<String>, // chromosome name corresponding to tid
    pub start_pos: u32, // start of interval
    pub end_pos: u32, // end of interval (inclusive)
}

pub fn u8_to_string(u: &[u8]) -> String {
    String::from_utf8(u.to_vec()).unwrap()
}

//
pub fn dna_vec(u: &[u8]) -> (Vec<char>) {
    let mut v: Vec<char> = Vec::with_capacity(u.len());
    for cu in AsciiExt::to_ascii_uppercase(u) {
        let c = cu as char;
        assert!(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
        v.push(c);
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
    pub var_ix: usize, // index into variant list
    pub allele: char, // allele call
    pub qual: LogProb, // LogProb probability the call is an error
}

pub struct Fragment {
    pub calls: Vec<FragCall>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PotentialVar {
    pub ix: usize, // index of this variant in the global var list
    pub chrom: String,
    pub pos0: usize,
    pub ref_allele: String,
    pub var_allele: String, 
    //pub pileup: Option(Vec<PileupElement>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VarList {
    pub lst: Vec<PotentialVar>,
    ix: HashMap<String, Vec<usize>>,
}

impl VarList {
    pub fn new(lst: Vec<PotentialVar>) -> VarList {
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

    pub fn get_variants_range(&self, interval: GenomicInterval) -> (Vec<PotentialVar>) {

        // vector of variants to fill and return
        let mut vlst: Vec<PotentialVar> = vec![];
        // get the varlist index of a nearby position on the left
        let chrom = match interval.chrom {
            Some(str) => str,
            None => panic!("Called get_variants_range with a range without chrom field"),
        };
        let mut i = self.ix.get(&chrom).unwrap()[(interval.start_pos as usize) / INDEX_FREQ];

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

    fn generate_test_lst1() -> VarList {
        let mut lst: Vec<PotentialVar> = vec![];
        lst.push(PotentialVar {
                     chrom: "chr1".to_string(),
                     pos0: 5,
                     ref_allele: "A".to_string(),
                     var_allele: "G".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr1".to_string(),
                     pos0: 1000,
                     ref_allele: "T".to_string(),
                     var_allele: "A".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr1".to_string(),
                     pos0: 2005,
                     ref_allele: "T".to_string(),
                     var_allele: "G".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr1".to_string(),
                     pos0: 2900,
                     ref_allele: "C".to_string(),
                     var_allele: "G".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr1".to_string(),
                     pos0: 6000,
                     ref_allele: "C".to_string(),
                     var_allele: "A".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr1".to_string(),
                     pos0: 10000,
                     ref_allele: "C".to_string(),
                     var_allele: "A".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 5,
                     ref_allele: "A".to_string(),
                     var_allele: "G".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 1000,
                     ref_allele: "T".to_string(),
                     var_allele: "A".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 2005,
                     ref_allele: "T".to_string(),
                     var_allele: "G".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 2900,
                     ref_allele: "C".to_string(),
                     var_allele: "G".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 6000,
                     ref_allele: "C".to_string(),
                     var_allele: "A".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 10000,
                     ref_allele: "C".to_string(),
                     var_allele: "A".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr3".to_string(),
                     pos0: 20200,
                     ref_allele: "C".to_string(),
                     var_allele: "G".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr3".to_string(),
                     pos0: 25100,
                     ref_allele: "A".to_string(),
                     var_allele: "C".to_string(),
                 });
        lst.push(PotentialVar {
                     chrom: "chr3".to_string(),
                     pos0: 30400,
                     ref_allele: "C".to_string(),
                     var_allele: "A".to_string(),
                 });
        VarList::new(lst)
    }

    #[test]
    fn test_varlist_get_variants_range1() {

        let mut vlst = generate_test_lst1();
        let c = "chr1".to_string();
        let p1 = 2500;
        let p2 = 8000;
        let interval = GenomicInterval {
            tid: 2,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<PotentialVar> = vec![];
        exp.push(PotentialVar {
                     chrom: "chr1".to_string(),
                     pos0: 2900,
                     ref_allele: "C".to_string(),
                     var_allele: "G".to_string(),
                 });
        exp.push(PotentialVar {
                     chrom: "chr1".to_string(),
                     pos0: 6000,
                     ref_allele: "C".to_string(),
                     var_allele: "A".to_string(),
                 });

        assert_eq!(vars, exp);
    }

    #[test]
    fn test_varlist_get_variants_range2() {

        let mut vlst = generate_test_lst1();
        let c = "chr2".to_string();
        let p1 = 0;
        let p2 = 3000;
        let interval = GenomicInterval {
            tid: 3,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<PotentialVar> = vec![];
        exp.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 5,
                     ref_allele: "A".to_string(),
                     var_allele: "G".to_string(),
                 });
        exp.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 1000,
                     ref_allele: "T".to_string(),
                     var_allele: "A".to_string(),
                 });
        exp.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 2005,
                     ref_allele: "T".to_string(),
                     var_allele: "G".to_string(),
                 });
        exp.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 2900,
                     ref_allele: "C".to_string(),
                     var_allele: "G".to_string(),
                 });

        assert_eq!(vars, exp);
    }

    #[test]
    fn test_varlist_get_variants_range3() {

        let mut vlst = generate_test_lst1();
        let c = "chr2".to_string();
        let p1 = 6000;
        let p2 = 10000;
        let interval = GenomicInterval {
            tid: 3,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<PotentialVar> = vec![];
        exp.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 6000,
                     ref_allele: "C".to_string(),
                     var_allele: "A".to_string(),
                 });
        exp.push(PotentialVar {
                     chrom: "chr2".to_string(),
                     pos0: 10000,
                     ref_allele: "C".to_string(),
                     var_allele: "A".to_string(),
                 });

        assert_eq!(vars, exp);
    }

    #[test]
    fn test_varlist_get_variants_range4() {

        let mut vlst = generate_test_lst1();
        let c = "chr3".to_string();
        let p1 = 20100;
        let p2 = 20200;
        let interval = GenomicInterval {
            tid: 3,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<PotentialVar> = vec![];
        exp.push(PotentialVar {
                     chrom: "chr3".to_string(),
                     pos0: 20200,
                     ref_allele: "C".to_string(),
                     var_allele: "G".to_string(),
                 });

        assert_eq!(vars, exp);
    }

    #[test]
    fn test_varlist_get_variants_range5() {

        let mut vlst = generate_test_lst1();
        let c = "chr3".to_string();
        let p1 = 20200;
        let p2 = 20200;
        let interval = GenomicInterval {
            tid: 3,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<PotentialVar> = vec![];
        exp.push(PotentialVar {
                     chrom: "chr3".to_string(),
                     pos0: 20200,
                     ref_allele: "C".to_string(),
                     var_allele: "G".to_string(),
                 });

        assert_eq!(vars, exp);
    }

    #[test]
    fn test_varlist_get_variants_range6() {

        let mut vlst = generate_test_lst1();
        let c = "chr3".to_string();
        let p1 = 25000;
        let p2 = 30500;
        let interval = GenomicInterval {
            tid: 3,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval);

        let mut exp: Vec<PotentialVar> = vec![];
        exp.push(PotentialVar {
                     chrom: "chr3".to_string(),
                     pos0: 25100,
                     ref_allele: "A".to_string(),
                     var_allele: "C".to_string(),
                 });
        exp.push(PotentialVar {
                     chrom: "chr3".to_string(),
                     pos0: 30400,
                     ref_allele: "C".to_string(),
                     var_allele: "A".to_string(),
                 });

        assert_eq!(vars, exp);
    }

}
