use std::ascii::AsciiExt;
use std::collections::HashMap;
use bio::stats::{LogProb, Prob};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use chrono::prelude::*;
use std::cmp::Ordering;

pub static INDEX_FREQ: usize = 1000;

pub fn print_time() -> String {
    Local::now().format("%Y-%m-%d %H:%M:%S").to_string()
}

// use this spacer instead of calling print_time() to have the spaces match up with
// lines that document the time
pub static SPACER: &str = "                   ";

#[derive(Clone)]
pub struct GenotypePriors {
    priors_dict: HashMap<(char, (char, char)), LogProb>, // (ref_allele, (allele1, allele2)) -> P(G)
}

impl GenotypePriors {

    // create a struct representing human variant priors
    // estimate prior probability of genotypes using strategy described here:
    // http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694485/
    // "prior probability of each genotype"
    // we've modified it in this case to have two extra options: 'I' and 'D'
    // these represent short insertions and short deletions.
    pub fn new() -> GenotypePriors {

        let hom_snp_rate = LogProb::from(Prob(0.0005));
        let het_snp_rate = LogProb::from(Prob(0.001));
        let hom_indel_rate = LogProb::from(Prob(0.00005));
        let het_indel_rate = LogProb::from(Prob(0.0001));

        // key of diploid_genotype_priors is (char,(char,char)) (ref_allele, G=(allele1,allele2))
        // key of haploid priors is (char, char) which is (ref_allele, allele1)
        let mut diploid_genotype_priors: HashMap<(char, (char, char)), LogProb> = HashMap::new();
        let mut haploid_genotype_priors: HashMap<(char, char), LogProb> = HashMap::new();

        // define base transitions
        let mut transition: HashMap<char, char> = HashMap::new();
        transition.insert('A', 'G');
        transition.insert('G', 'A');
        transition.insert('T', 'C');
        transition.insert('C', 'T');

        let alleles: Vec<char> = vec!['A', 'C', 'G', 'T'];
        let genotypes: Vec<(char, char)> = vec![ // all combinations of DNA bases
                                                ('A','A'),('A','C'),('A','G'),('A','T'),
                                                ('C', 'C'),('C', 'G'),('C', 'T'),
                                                ('G', 'G'),('G', 'T'),
                                                ('T', 'T'),
                                                 // all combinations of a DNA base and an Insertion or Deletion
                                                ('A','D'),('A','I'),
                                                ('C','D'),('C','I'),
                                                ('G','D'),('G','I'),
                                                ('T','D'),('T','I'),
                                                 // all combinations of multiple Insertion or Deletion
                                                ('D','D'),('I','I'),('D','I')];


        for aref in &alleles {
            let allele = *aref;
            // priors on haploid alleles
            haploid_genotype_priors.insert((allele, allele), LogProb::ln_one_minus_exp(&(het_snp_rate + het_indel_rate)));
            haploid_genotype_priors.insert((allele, *transition.get(&allele).unwrap()), het_snp_rate + LogProb::from(Prob(4.0 / 6.0)));

            for tref in &alleles {
                let transversion = *tref;
                if haploid_genotype_priors.contains_key(&(allele, transversion)) {
                    continue;
                }
                haploid_genotype_priors.insert((allele, transversion), het_snp_rate + LogProb::from(Prob(1.0 / 6.0)));
            }

            // assume indel has a 0.5 chance of being insertion, 0.5 chance of deletion
            // scale the prior by 1/3 since the SNP events are divided into 3 types... (transition to each base)
            haploid_genotype_priors.insert((allele, 'D'), het_indel_rate + LogProb::from(Prob(0.5)) + LogProb::from(Prob(1.0 / 4.0)));
            haploid_genotype_priors.insert((allele, 'I'), het_indel_rate + LogProb::from(Prob(0.5)) + LogProb::from(Prob(1.0 / 4.0)));

            for gt in &genotypes {
                let (g1, g2) = *gt;
                // probability of homozygous reference is the probability of neither het or hom SNP
                if g1 == g2 && g1 == allele {
                    // g1 and g2 are the reference bases
                    let var_rate = LogProb::ln_sum_exp(&[hom_snp_rate, het_snp_rate, hom_indel_rate, het_indel_rate]);
                    let one_minus_snp_rate = LogProb::ln_one_minus_exp(&var_rate);
                    diploid_genotype_priors.insert((allele, *gt), one_minus_snp_rate);
                } else if g1 == g2 && g1 != allele {
                    // this could be multiple indels, must handle indel case first
                    if g2 == 'D' || g2 == 'I' {
                        diploid_genotype_priors.insert((allele, *gt), hom_indel_rate + LogProb::from(Prob(0.5)) + LogProb::from(Prob(1.0 / 4.0)));
                    } else if g1 == *transition.get(&allele).unwrap() { // otherwise it is a homozygous SNV
                        // transitions are 4 times as likely as transversions
                        diploid_genotype_priors.insert((allele, *gt), hom_snp_rate + LogProb::from(Prob(4.0 / 6.0)));
                    } else {
                        // transversion
                        diploid_genotype_priors.insert((allele, *gt), hom_snp_rate + LogProb::from(Prob(1.0 / 6.0)));
                    }
                } else { // else it's the product of the haploid priors
                    diploid_genotype_priors.insert((allele, *gt), *haploid_genotype_priors.get(&(allele, g1)).unwrap() +
                        *(haploid_genotype_priors.get(&(allele, g2))).unwrap());
                }
            }
        }

        eprintln!("{} GENOTYPE PRIORS:", SPACER);
        eprintln!("{} REF G1/G2 PROB", SPACER);
        //let mut priors_vec: Vec<_> = diploid_genotype_priors.iter()
        //    .map(|(&(ra,(g1,g2)),&p)| (&(ra,(g1,g2)),&*Prob::from(p)))
        //    .collect::<Vec<_>>();
        //priors_vec.sort_by(|a, b| a.partial_cmp(b).unwrap());

        for (&(ref ra, (ref g1, ref g2)), &p) in &diploid_genotype_priors {
            eprintln!("{} {} {}/{} {}", SPACER, ra, g1, g2, *Prob::from(p));
        }

        GenotypePriors{priors_dict: diploid_genotype_priors}
    }

    // takes a vector of strings representing alleles (i.e. from Var.alleles), with the 0-th allele being reference
    // and a phased genotype
    // represented as indices into the alleles vector
    // returns the prior probability of that genotype

    // TODO: currently MNPs are going to be assigned the prior of the first SNV in the MNP
    // not ideal, but should suffice for the time being
    pub fn get(&self, alleles: &Vec<String>, genotype: [u8; 2]) -> LogProb {

        let mut ra = alleles[0].chars().nth(0).unwrap();

        let mut g1 = if alleles[genotype[0] as usize].len() == alleles[0].len() {
            alleles[genotype[0] as usize].chars().nth(0).unwrap()
        } else if alleles[genotype[0] as usize].len() > alleles[0].len() {
            'I'
        } else {
            'D'
        };

        if ra == g1 {
            for i in 0..alleles[genotype[0] as usize].len() {
                ra = alleles[0].chars().nth(i).unwrap();
                g1 = alleles[genotype[0] as usize].chars().nth(i).unwrap();
                if ra != g1 {
                    break;
                }
            }
        }

        let g2 = if alleles[genotype[1] as usize].len() == alleles[0].len() {
            alleles[genotype[1] as usize].chars().nth(0).unwrap()
        } else if alleles[genotype[1] as usize].len() > alleles[0].len() {
            'I'
        } else {
            'D'
        };

        if ra == g2 {
            for i in 0..alleles[genotype[1] as usize].len() {
                ra = alleles[0].chars().nth(i).unwrap();
                g1 = alleles[genotype[1] as usize].chars().nth(i).unwrap();
                if ra != g2 {
                    break;
                }
            }
        }

        let ln_half = LogProb::from(Prob(0.5));

        match self.priors_dict.get(&(ra, (g1, g2))) {
            Some(p) => {
                if g1 != g2 {return ln_half+*p;} else {return *p;}
            },
            None => {
                match self.priors_dict.get(&(ra, (g2, g1))) {
                    Some(p) => {
                        if g1 != g2 {return ln_half+*p;} else {return *p;}
                    },
                    None => { println!("{} ({},{})",ra,g2,g1); panic!("Genotype not in genotype priors."); }
                };
            }
        }
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
    pub anchor_length: usize,
    pub short_hap_max_snvs: usize,
    pub max_window_padding: usize,
}

// these parameters describe state transition probabilities for a pair HMM
// there are two kinds: "eq" transition probs and "neq" transition_probs
// the correct kind to use depends on sequence context.
// the "eq" probabilities are for the case where accepting a match results in equal bases
// the "neq" probabilities are for the case where accepting a match results in different bases

#[derive(Clone, Copy)]
pub struct TransitionProbs {
    pub match_from_match: f64,
    pub insertion_from_match: f64,
    pub deletion_from_match: f64,
    pub insertion_from_insertion: f64,
    pub match_from_insertion: f64,
    pub deletion_from_deletion: f64,
    pub match_from_deletion: f64,
}

#[derive(Clone, Copy)]
pub struct LnTransitionProbs {
    pub match_from_match: LogProb,
    pub insertion_from_match: LogProb,
    pub deletion_from_match: LogProb,
    pub insertion_from_insertion: LogProb,
    pub match_from_insertion: LogProb,
    pub deletion_from_deletion: LogProb,
    pub match_from_deletion: LogProb,
}

impl TransitionProbs {
    pub fn ln(&self) -> LnTransitionProbs {
        LnTransitionProbs {
            match_from_match: LogProb::from(Prob(self.match_from_match)),
            insertion_from_match: LogProb::from(Prob(self.insertion_from_match)),
            deletion_from_match: LogProb::from(Prob(self.deletion_from_match)),
            insertion_from_insertion: LogProb::from(Prob(self.insertion_from_insertion)),
            match_from_insertion: LogProb::from(Prob(self.match_from_insertion)),
            deletion_from_deletion: LogProb::from(Prob(self.deletion_from_deletion)),
            match_from_deletion: LogProb::from(Prob(self.match_from_deletion)),
        }
    }
}

#[derive(Clone, Copy)]
pub struct EmissionProbs {
    pub equal: f64,
    pub not_equal: f64
}

#[derive(Clone, Copy)]
pub struct LnEmissionProbs {
    pub equal: LogProb,
    pub not_equal: LogProb
}

impl EmissionProbs {
    pub fn ln(&self) -> LnEmissionProbs {
        LnEmissionProbs {
            equal: LogProb::from(Prob(self.equal)),
            not_equal: LogProb::from(Prob(self.not_equal))
        }
    }
}

#[derive(Clone, Copy)]
pub struct AlignmentParameters {
    pub transition_probs: TransitionProbs,
    pub emission_probs: EmissionProbs
}

#[derive(Clone, Copy)]
pub struct LnAlignmentParameters {
    pub transition_probs: LnTransitionProbs,
    pub emission_probs: LnEmissionProbs
}

impl AlignmentParameters {
    pub fn ln(&self) -> LnAlignmentParameters {
        LnAlignmentParameters {
            transition_probs: self.transition_probs.ln(),
            emission_probs: self.emission_probs.ln()
        }
    }
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
    pub frag_ix: Option<usize>, // index into fragment list
    pub var_ix: usize, // index into variant list
    pub allele: u8, // allele call
    pub qual: LogProb, // LogProb probability the call is an error
    pub one_minus_qual: LogProb, // LogProb probability the call is correct
}

#[derive(Clone)]
pub struct Fragment {
    pub id: String,
    pub calls: Vec<FragCall>,
    pub p_read_hap: [LogProb; 2]
}


#[derive(Debug, Clone)]
pub struct Var {
    pub ix: usize,
    pub old_ix: Option<usize>,
    // index of this variant in the global var list
    pub tid: usize,
    pub chrom: String,
    pub pos0: usize,
    pub alleles: Vec<String>, // ref allele is alleles[0] and each that follows is a variant allele
    pub dp: usize,
    // depth of coverage
    pub allele_counts: Vec<usize>, // indices match up with those of Var.alleles
    pub ambiguous_count: usize,
    pub qual: f64,
    pub filter: String,
    pub genotype: [u8; 2],
    pub gq: f64,
    pub genotype_post: Vec<Vec<LogProb>>,  // genotype posteriors[a1][a2] is log posterior of phased a1|a2 haplotype
    // e.g. genotype_posteriors[2][0] is the log posterior probability of 2|0 haplotype
    pub phase_set: Option<usize>,
    pub mec: usize,
    pub mec_frac: f64,
    pub called: bool
    //pub pileup: Option(Vec<PileupElement>),
}

impl Var {
    fn longest_allele_len(&self) -> usize {
        self.alleles.iter().map(|x| x.len()).max().unwrap()
    }
}
impl Ord for Var {
    fn cmp(&self, other: &Var) -> Ordering {
        if self.tid == other.tid {
            self.pos0.cmp(&other.pos0)
        } else {
            self.tid.cmp(&other.tid)
        }
    }
}

impl PartialOrd for Var {
    fn partial_cmp(&self, other: &Var) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Var {
    fn eq(&self, other: &Var) -> bool {
        self.tid == other.tid && self.pos0 == other.pos0
    }
}

impl Eq for Var {}

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
    pub fn index_lst(&mut self) {
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

        while i < self.lst.len() &&
            self.lst[i].pos0 + self.lst[i].longest_allele_len() <= interval.end_pos as usize {
            if self.lst[i].pos0 >= interval.start_pos as usize {
                vlst.push(self.lst[i].clone());
            }
            i += 1;
        }
        vlst
    }

    pub fn backup_indices(&mut self) {
        for ref mut var in &mut self.lst {
            var.old_ix = Some(var.ix);
        }
    }

    pub fn combine (&mut self, other: &mut VarList) {
        other.lst.append(&mut self.lst);
        other.lst.sort();

        let mut new_vlst: Vec<Var> = vec![];
        // any variants that start at or after area_start, but before area_end, are added to the group
        let mut var_group: Vec<Var> = vec![];
        let mut area_tid = other.lst[0].tid;
        let mut area_start = other.lst[0].pos0;
        let mut area_end = other.lst[0].pos0 + other.lst[0].alleles[0].len();

        for var in &other.lst {
            if var.tid == area_tid && var.pos0 >= area_start && var.pos0 < area_end {
                var_group.push(var.clone());
                // the new var might overlap, but also extend the variant area
                if var.pos0 + var.alleles[0].len() > area_end {
                    area_end = var.pos0 + var.alleles[0].len();
                }
            } else {

                ///////////////////////////////////////////////////////////////////////////////////
                // PROCEDURE TO MERGE MULTIPLE VCF VARIANTS, EACH WITH MULTIPLE ALLELES
                // 1. find the earliest variant position in the variant list.
                // 2. for every variant with a later position, pad the variant's alleles
                //        at the beginning with enough ref chars to match this position, and
                //        update their position so all the positions are the same.
                // 3. find the longest ref allele out of all the vars. This will be the new ref allele.
                // 4. for every other variant (not longest ref), pad EACH allele at the end with ref chars
                //        until the ref allele is the same length as (step 3)
                // 5. add the unified ref_allele as 0-th allele of new allele list
                // 6. combine the remaining alleles from every variant (non-ref)
                //        into a single vector, sort, and remove duplicates. (vec.dedup())
                // 7. append these unique variant alleles to the new allele list.

                // 1. find the earliest variant position in the variant list.
                let mut min_pos0 = var_group[0].pos0;
                let mut min_pos0_refseq = var_group[0].alleles[0].clone();
                for v in &var_group {
                    if v.pos0 < min_pos0
                        || (v.pos0 == min_pos0 && v.alleles[0].len() > min_pos0_refseq.len()) {
                        min_pos0 = v.pos0;
                        min_pos0_refseq = v.alleles[0].clone();
                    }
                }

                // 2. for every variant with a later position, pad the variant's alleles
                //        at the beginning with enough ref chars to match this position, and
                //        update their position so all the positions are the same.

                for mut v in &mut var_group {
                    if v.pos0 > min_pos0 {
                        let diff = v.pos0 - min_pos0;
                        // we need to steal the first diff bases from min_pos0_refseq
                        let prefix_seq: String = min_pos0_refseq.chars().take(diff).collect();
                        let mut new_alleles = vec![];

                        for ref allele in &v.alleles {
                            let mut a = prefix_seq.clone();
                            a.push_str(&allele);
                            new_alleles.push(a);
                        }

                        v.alleles = new_alleles;
                        v.pos0 = min_pos0
                    }
                }

                for v in &var_group {
                    assert!(v.pos0 == min_pos0);
                }

                // 3. find the longest ref allele out of all the vars
                let mut longest_ref = var_group[0].alleles[0].clone();
                for  v in var_group.iter() {
                    if v.alleles[0].len() > longest_ref.len() {
                        longest_ref = v.alleles[0].clone();
                    }
                }

                // 4. for every other variant (not longest ref), pad EACH allele at the end with ref chars
                //        until the ref allele is the same length as (step 3)

                for ref mut v in &mut var_group {
                    if v.alleles[0].len() < longest_ref.len() {
                        let diff = longest_ref.len() - v.alleles[0].len();
                        // we need to steal the last diff bases from longest_ref
                        let suffix_seq: String = longest_ref.chars().skip(v.alleles[0].len()).take(diff).collect();
                        let mut new_alleles = vec![];

                        for ref allele in &v.alleles {
                            let mut a = (*allele).clone();
                            a.push_str(&suffix_seq.clone());
                            new_alleles.push(a);
                        }

                        v.alleles = new_alleles;
                    }
                }

                for v in &var_group {
                    assert!(v.alleles[0] == longest_ref);
                }

                // 5. add the unified ref_allele as 0-th allele of new allele list
                let mut new_allele_lst = vec![var_group[0].alleles[0].clone()];

                // 6. combine the remaining alleles from every variant (non-ref)
                //        into a single vector, sort, and remove duplicates. (vec.dedup())
                let mut var_alleles: Vec<String> = vec![];
                for v in &var_group {
                    for a in v.alleles[1..].iter() {
                        var_alleles.push(a.clone());
                    }
                }
                var_alleles.sort();
                var_alleles.dedup();

                // 7. append these unique variant alleles to the new allele list.

                new_allele_lst.append(&mut var_alleles);
                ///////////////////////////////////////////////////////////////////////////////////

                let mut new_v = var_group[0].clone();
                new_v.alleles = new_allele_lst;
                new_v.genotype = [0u8,0u8];
                new_v.gq = 0.0;
                new_v.genotype_post = vec![vec![]];
                new_v.phase_set = None;

                new_vlst.push(new_v);
                // clear out the variant group and add the current variant
                var_group.clear();
                var_group.push(var.clone());
                // set the new "variant area" within which new variants will be said to overlap with this one
                area_tid = var.tid;
                area_start = var.pos0;
                area_end = var.pos0 + var.alleles[0].len();
            }
        }

        // set the list to the new merged list and re-index it
        self.lst = new_vlst;
        for i in 0..self.lst.len() {
            self.lst[i].ix = i;
        }
        self.ix = HashMap::new();
        self.index_lst();

        // clear out the other VarList since we've mutated it beyond saving
        other.lst.clear();
        other.ix.clear();
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
