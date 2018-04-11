
use bio::stats::{LogProb};
use util::*;
use std::collections::HashMap;
use std::cmp::Ordering;
use genotype_probs::*;

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
    pub genotype: Genotype,
    pub gq: f64,
    pub genotype_post: GenotypeProbs,  // genotype posteriors[a1][a2] is log posterior of phased a1|a2 haplotype
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

    pub fn possible_genotypes(&self) -> Vec<Genotype> {
        possible_genotypes(&self.alleles)
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
        v.sort();
        v
    }

    pub fn sort(&mut self) {
        self.lst.sort();
        self.add_ix();
        self.ix.clear();
        self.index_lst();
    }

    pub fn assert_sorted(&self) {
        if self.lst.len() == 0 {
            return;
        }
        for i in 0..self.lst.len()-1{
            assert!((self.lst[i].tid < self.lst[i+1].tid) ||
                    (self.lst[i].tid == self.lst[i+1].tid && self.lst[i].pos0 <= self.lst[i+1].pos0));
            assert_eq!(self.lst[i].ix, i);
            assert_eq!(self.lst[i+1].ix, i+1);
        }
    }

    fn add_ix(&mut self) {
        for i in 0..self.lst.len() {
            self.lst[i].ix = i;
        }
    }

    fn index_lst(&mut self) {
        // need to throw error if list isn't sorted
        self.assert_sorted();
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

    fn combine_variant_group(var_group: &mut Vec<Var>) -> Var {
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
        for v in var_group.iter() {
            if v.pos0 < min_pos0
                || (v.pos0 == min_pos0 && v.alleles[0].len() > min_pos0_refseq.len()) {
                min_pos0 = v.pos0;
                min_pos0_refseq = v.alleles[0].clone();
            }
        }

        // find the last variant end position
        let mut max_end = var_group[0].pos0;
        for v in var_group.iter() {
            if v.pos0 + v.alleles[0].len() > max_end {
            max_end = v.pos0 + v.alleles[0].len();
            }

        }

        // we are going to combine the ref alleles from each variant to infer the longer
        // representation of the reference sequence for the variant region.
        let mut ref_seq = vec!['N'; max_end - min_pos0];
        for v in var_group.iter() {
            let ref_allele_vec: Vec<char> = v.alleles[0].chars().collect::<Vec<char>>();
            let vs: usize = v.pos0 - min_pos0; // variant position in ref_seq vector
            let ve: usize = vs + ref_allele_vec.len();
            for (a,r) in (vs..ve).enumerate() {
                if ref_seq[r] == 'N' {
                    ref_seq[r] = ref_allele_vec[a];
                } else {
                    assert_eq!(ref_seq[r], ref_allele_vec[a]);
                }
            }
        }

        for c in &ref_seq {
            assert_ne!(c, &'N');
        }

        // 2. for every variant with a later position, pad the variant's alleles
        //        at the beginning with enough ref chars to match this position, and
        //        update their position so all the positions are the same.

        for mut v in var_group.iter_mut() {
            if v.pos0 > min_pos0 {
                let diff = v.pos0 - min_pos0;
                // we need to steal the first diff bases from refseq
                let prefix_seq: String = ref_seq[0..diff].iter().collect();
                let mut new_alleles = vec![];

                for allele in v.alleles.iter_mut() {
                    let mut a = prefix_seq.clone();
                    a.push_str(allele);
                    new_alleles.push(a);
                }

                v.alleles = new_alleles;
                v.pos0 = min_pos0
            }
        }

        for v in var_group.iter() {
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

        for v in var_group.iter_mut() {
            if v.alleles[0].len() < longest_ref.len() {
                let diff = longest_ref.len() - v.alleles[0].len();
                // we need to steal the last diff bases from longest_ref
                let suffix_seq: String = longest_ref.chars().skip(v.alleles[0].len()).take(diff).collect();
                //let mut new_alleles = vec![];

                for mut allele in v.alleles.iter_mut() {
                    allele.push_str(&suffix_seq.clone());
                    //new_alleles.push(allele);
                }

                //v.alleles = new_alleles;
            }
        }

        for v in var_group.iter() {
            assert!(v.alleles[0] == longest_ref);
        }

        // 5. add the unified ref_allele as 0-th allele of new allele list
        let mut new_allele_lst = vec![var_group[0].alleles[0].clone()];

        // 6. combine the remaining alleles from every variant (non-ref)
        //        into a single vector, sort, and remove duplicates. (vec.dedup())
        let mut var_alleles: Vec<String> = vec![];
        for v in var_group.iter() {
            for a in v.alleles[1..].iter() {
                var_alleles.push(a.clone());
            }
        }
        var_alleles.sort();
        var_alleles.dedup();

        // 7. append these unique variant alleles to the new allele list.

        new_allele_lst.append(&mut var_alleles);
        ///////////////////////////////////////////////////////////////////////////////////

        assert!(new_allele_lst.len() >= 2);

        let mut new_v = var_group[0].clone();
        new_v.allele_counts = vec![0; new_allele_lst.len()];
        new_v.alleles = new_allele_lst.clone();
        new_v.genotype = Genotype(0,0);
        new_v.gq = 0.0;
        new_v.genotype_post = GenotypeProbs::uniform(new_v.alleles.len());
        new_v.phase_set = None;

        new_v
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

                let new_v = VarList::combine_variant_group(&mut var_group);

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

        // the var_group might have one last variant left in it
        if var_group.len() > 0 {
            let new_v = VarList::combine_variant_group(&mut var_group);
            new_vlst.push(new_v);
        }

        // set the list to the new merged list and re-index it
        self.lst = new_vlst;

        self.sort();
        self.assert_sorted();

        // clear out the other VarList since we've mutated it beyond saving
        other.lst.clear();
        other.ix.clear();
    }
}

pub fn var_filter(varlist: &mut VarList, density_qual: f64, density_dist: usize, density_count: usize, max_depth: Option<u32>) {

    for i in 0..varlist.lst.len() {
        if varlist.lst[i].qual < density_qual { continue; }

        let mut count = 0;
        for j in i + 1..varlist.lst.len() {
            if varlist.lst[j].pos0 - varlist.lst[i].pos0 > density_dist {
                break;
            }
            if varlist.lst[j].qual < density_qual { continue; }
            count += 1;
            if count > density_count {
                for k in i..j + 1 {
                    varlist.lst[k].filter = "dn".to_string();
                }
            }
        }
    }
    match max_depth {
        Some(dp) => {
            for i in 0..varlist.lst.len() {
                if varlist.lst[i].dp > dp as usize {
                    if varlist.lst[i].filter == ".".to_string() || varlist.lst[i].filter == "PASS".to_string() {
                        varlist.lst[i].filter = "dp".to_string();
                    } else {
                        varlist.lst[i].filter.push_str(";dp");
                    }
                }
            }
        }
        None => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /**********************************************************************************************/
    // TEST VARIANT RANGE LOOKUP
    /**********************************************************************************************/

    fn pos_alleles_eq(vlst1: Vec<Var>, vlst2: Vec<Var>) -> bool {
        for (v1, v2) in vlst1.iter().zip(vlst2.iter()){
            if v1 != v2 || v1.chrom != v2.chrom || v1.alleles != v2.alleles || v1.ix != v2.ix  {
                return false;
            }
        }
        true
    }

    fn varlist_pos_alleles_eq(vlst1: VarList, vlst2: VarList) -> bool {
        pos_alleles_eq(vlst1.lst, vlst2.lst)
    }

    fn generate_var1(ix: usize, tid: usize, chrom: String, pos0: usize, ra: String, aa: String) -> Var {
        Var {
            ix: ix,
            old_ix: None,
            tid: tid,
            chrom: chrom,
            pos0: pos0,
            alleles: vec![ra,aa],
            dp: 40,
            allele_counts: vec![20,20],
            ambiguous_count: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: Genotype(0,1),
            gq: 0.0,
            genotype_post: GenotypeProbs::uniform(2),
            phase_set: None,
            mec: 0,
            mec_frac: 0.0,
            called: true
        }
    }

    fn generate_test_lst1() -> VarList {
        let mut lst: Vec<Var> = vec![];

        lst.push( generate_var1(0, 0, "chr1".to_string(), 5, "A".to_string(), "G".to_string()));
        lst.push( generate_var1(1, 0, "chr1".to_string(), 1000, "T".to_string(), "A".to_string()));
        lst.push( generate_var1(2, 0, "chr1".to_string(), 2005, "T".to_string(), "G".to_string()));
        lst.push( generate_var1(3, 0, "chr1".to_string(), 2900, "C".to_string(), "G".to_string()));
        lst.push( generate_var1(4, 0, "chr1".to_string(), 6000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(5, 0, "chr1".to_string(), 10000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(6, 1, "chr2".to_string(), 5, "A".to_string(), "G".to_string()));
        lst.push( generate_var1(7, 1, "chr2".to_string(), 1000, "T".to_string(), "A".to_string()));
        lst.push( generate_var1(8, 1, "chr2".to_string(), 2005, "T".to_string(), "G".to_string()));
        lst.push( generate_var1(9, 1, "chr2".to_string(), 2900, "C".to_string(), "G".to_string()));
        lst.push( generate_var1(10, 1, "chr2".to_string(), 6000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(11, 1, "chr2".to_string(), 10000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(12, 2, "chr3".to_string(), 20200, "C".to_string(), "G".to_string()));
        lst.push( generate_var1(13, 2, "chr3".to_string(), 25100, "A".to_string(), "C".to_string()));
        lst.push( generate_var1(14, 2, "chr3".to_string(), 30400, "C".to_string(), "A".to_string()));

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
        exp.push( generate_var1(3, 0, "chr1".to_string(), 2900, "C".to_string(), "G".to_string()));
        exp.push( generate_var1(4, 0, "chr1".to_string(), 6000, "C".to_string(), "A".to_string()));

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
        exp.push( generate_var1(6, 1, "chr2".to_string(), 5, "A".to_string(), "G".to_string()));
        exp.push( generate_var1(7, 1, "chr2".to_string(), 1000, "T".to_string(), "A".to_string()));
        exp.push( generate_var1(8, 1, "chr2".to_string(), 2005, "T".to_string(), "G".to_string()));
        exp.push( generate_var1(9, 1, "chr2".to_string(), 2900, "C".to_string(), "G".to_string()));

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
        exp.push( generate_var1(10, 1, "chr2".to_string(), 6000, "C".to_string(), "A".to_string()));
        exp.push( generate_var1(11, 1, "chr2".to_string(), 10000, "C".to_string(), "A".to_string()));

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
        exp.push( generate_var1(12, 2, "chr3".to_string(), 20200, "C".to_string(), "G".to_string()));

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
        exp.push( generate_var1(12, 2, "chr3".to_string(), 20200, "C".to_string(), "G".to_string()));


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

        exp.push( generate_var1(13, 2, "chr3".to_string(), 25100, "A".to_string(), "C".to_string()));
        exp.push( generate_var1(14, 2, "chr3".to_string(), 30400, "C".to_string(), "A".to_string()));

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

    /**********************************************************************************************/
    // TEST VARIANT SORTING
    /**********************************************************************************************/

    // generate an unsorted version of test_lst1, with meaningless indices. for testing sorting.
    fn generate_test_lst1_unsorted1() -> VarList {
        let mut lst: Vec<Var> = vec![];

        lst.push( generate_var1(0, 1, "chr2".to_string(), 2005, "T".to_string(), "G".to_string()));
        lst.push( generate_var1(1000, 1, "chr2".to_string(), 1000, "T".to_string(), "A".to_string()));
        lst.push( generate_var1(0, 0, "chr1".to_string(), 2900, "C".to_string(), "G".to_string()));
        lst.push( generate_var1(100, 0, "chr1".to_string(), 5, "A".to_string(), "G".to_string()));
        lst.push( generate_var1(4, 0, "chr1".to_string(), 6000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(14, 2, "chr3".to_string(), 30400, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(13, 2, "chr3".to_string(), 25100, "A".to_string(), "C".to_string()));
        lst.push( generate_var1(10, 0, "chr1".to_string(), 10000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(1, 1, "chr2".to_string(), 5, "A".to_string(), "G".to_string()));
        lst.push( generate_var1(2, 0, "chr1".to_string(), 2005, "T".to_string(), "G".to_string()));
        lst.push( generate_var1(10, 1, "chr2".to_string(), 6000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(10, 1, "chr2".to_string(), 2900, "C".to_string(), "G".to_string()));
        lst.push( generate_var1(11, 1, "chr2".to_string(), 10000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(2, 0, "chr1".to_string(), 1000, "T".to_string(), "A".to_string()));
        lst.push( generate_var1(12, 2, "chr3".to_string(), 20200, "C".to_string(), "G".to_string()));

        VarList::new(lst)
    }

    // generate another unsorted version of test_lst1, with meaningless indices. for testing sorting.
    fn generate_test_lst1_unsorted2() -> VarList {
        let mut lst: Vec<Var> = vec![];

        lst.push( generate_var1(10, 1, "chr2".to_string(), 5, "A".to_string(), "G".to_string()));
        lst.push( generate_var1(1, 0, "chr1".to_string(), 6000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(0, 1, "chr2".to_string(), 2005, "T".to_string(), "G".to_string()));
        lst.push( generate_var1(0, 0, "chr1".to_string(), 5, "A".to_string(), "G".to_string()));
        lst.push( generate_var1(1000, 0, "chr1".to_string(), 1000, "T".to_string(), "A".to_string()));
        lst.push( generate_var1(0, 1, "chr2".to_string(), 1000, "T".to_string(), "A".to_string()));
        lst.push( generate_var1(12, 2, "chr3".to_string(), 20200, "C".to_string(), "G".to_string()));
        lst.push( generate_var1(2, 2, "chr3".to_string(), 25100, "A".to_string(), "C".to_string()));
        lst.push( generate_var1(100, 0, "chr1".to_string(), 2900, "C".to_string(), "G".to_string()));
        lst.push( generate_var1(9, 0, "chr1".to_string(), 2005, "T".to_string(), "G".to_string()));
        lst.push( generate_var1(5, 0, "chr1".to_string(), 10000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(3, 1, "chr2".to_string(), 2900, "C".to_string(), "G".to_string()));
        lst.push( generate_var1(10, 1, "chr2".to_string(), 6000, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(50, 2, "chr3".to_string(), 30400, "C".to_string(), "A".to_string()));
        lst.push( generate_var1(3, 1, "chr2".to_string(), 10000, "C".to_string(), "A".to_string()));

        VarList::new(lst)
    }

    #[test]
    fn test_varlist_sort1() {
        let mut vlst_unsorted = generate_test_lst1_unsorted1();

        vlst_unsorted.sort();
        vlst_unsorted.assert_sorted();
    }

    #[test]
    fn test_varlist_sort2() {
        let mut vlst_unsorted = generate_test_lst1_unsorted2();

        vlst_unsorted.sort();
        vlst_unsorted.assert_sorted();
    }

    /**********************************************************************************************/
    // TEST COMBINING VARLISTS
    /**********************************************************************************************/

    fn generate_var2(ix: usize, tid: usize, chrom: String, pos0: usize, alleles: Vec<String>) -> Var {
        Var {
            ix: ix,
            old_ix: None,
            tid: tid,
            chrom: chrom,
            pos0: pos0,
            alleles: alleles,
            dp: 40,
            allele_counts: vec![20,20],
            ambiguous_count: 0,
            qual: 0.0,
            filter: ".".to_string(),
            genotype: Genotype(0,1),
            gq: 0.0,
            genotype_post: GenotypeProbs::uniform(2),
            phase_set: None,
            mec: 0,
            mec_frac: 0.0,
            called: true
        }
    }

    // combining two identical lists should result in the same list
    #[test]
    fn test_varlist_combine_same_lists() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "A".to_string()]));
        lst1.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let vlst1_bak = vlst1.clone();

        let mut lst2: Vec<Var> = vec![];
        lst2.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst2.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "A".to_string()]));
        lst2.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst2 = VarList::new(lst2);

        vlst1.combine(&mut vlst2);
        assert!(varlist_pos_alleles_eq(vlst1, vlst1_bak));
    }

    // combine varlists with different variants into essentially the union of the variant set
    #[test]
    fn test_varlist_combine_not_same() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "A".to_string()]));
        lst1.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let mut lst2: Vec<Var> = vec![];
        lst2.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst2.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "A".to_string()]));
        let mut vlst2 = VarList::new(lst2);

        let mut lst3: Vec<Var> = vec![];
        lst3.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst3.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "A".to_string()]));
        lst3.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let exp = VarList::new(lst3);

        vlst1.combine(&mut vlst2);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // combine a variant where a SNVs overlaps a deletion
    #[test]
    fn test_varlist_combine_deletion() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["TTT".to_string(), "T".to_string()]));
        lst1.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let mut lst2: Vec<Var> = vec![];
        lst2.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst2.push( generate_var2(1, 0, "chr1".to_string(), 101, vec!["T".to_string(), "C".to_string()]));
        lst2.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst2 = VarList::new(lst2);

        let mut lst3: Vec<Var> = vec![];
        lst3.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string(),]));
        lst3.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["TTT".to_string(), "T".to_string(), "TCT".to_string(),]));
        lst3.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let exp = VarList::new(lst3);

        vlst1.combine(&mut vlst2);
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // combine a variant where a SNVs overlaps an insertion
    #[test]
    fn test_varlist_combine_insertion() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "C".to_string()]));
        lst1.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let mut lst2: Vec<Var> = vec![];
        lst2.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst2.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "TAAAA".to_string()]));
        lst2.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst2 = VarList::new(lst2);

        let mut lst3: Vec<Var> = vec![];
        lst3.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string(),]));
        lst3.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "C".to_string(), "TAAAA".to_string(),]));
        lst3.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let exp = VarList::new(lst3);

        vlst1.combine(&mut vlst2);
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // make sure that if the insertion is only adjacent, it is left separate from SNV.
    #[test]
    fn test_varlist_combine_nearby_insertion() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 101, vec!["T".to_string(), "C".to_string()]));
        lst1.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let mut lst2: Vec<Var> = vec![];
        lst2.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst2.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "TAAAA".to_string()]));
        lst2.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst2 = VarList::new(lst2);

        let mut lst3: Vec<Var> = vec![];
        lst3.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string(),]));
        lst3.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "TAAAA".to_string(),]));
        lst3.push( generate_var2(2, 0, "chr1".to_string(), 101, vec!["T".to_string(), "C".to_string()]));
        lst3.push( generate_var2(3, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let exp = VarList::new(lst3);

        vlst1.combine(&mut vlst2);
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert_eq!(vlst1.lst[2].alleles, exp.lst[2].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // combine a variant where multiple SNVs overlap the same deletion
    #[test]
    fn test_varlist_combine_deletion_multiple_snp() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["TTT".to_string(), "T".to_string()]));
        lst1.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let mut lst2: Vec<Var> = vec![];
        lst2.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst2.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "C".to_string()]));
        lst2.push( generate_var2(2, 0, "chr1".to_string(), 102, vec!["T".to_string(), "C".to_string()]));
        lst2.push( generate_var2(3, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst2 = VarList::new(lst2);

        let mut lst3: Vec<Var> = vec![];
        lst3.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string(),]));
        lst3.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["TTT".to_string(), "CTT".to_string(), "T".to_string(), "TTC".to_string()]));
        lst3.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let exp = VarList::new(lst3);

        vlst1.combine(&mut vlst2);
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // combine a variant where multiple SNVs overlap the same deletion
    #[test]
    fn test_varlist_combine_overlapping_deletions() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["TTT".to_string(), "T".to_string()]));
        lst1.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let mut lst2: Vec<Var> = vec![];
        lst2.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst2.push( generate_var2(1, 0, "chr1".to_string(), 102, vec!["TCC".to_string(), "T".to_string()]));
        lst2.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst2 = VarList::new(lst2);

        let mut lst3: Vec<Var> = vec![];
        lst3.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string(),]));
        lst3.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["TTTCC".to_string(), "TCC".to_string(),"TTT".to_string()]));
        lst3.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let exp = VarList::new(lst3);

        vlst1.combine(&mut vlst2);
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // combine a variant where multiple SNVs overlap the same deletion
    #[test]
    fn test_varlist_combine_overlapping_insertion_deletion() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["TTT".to_string(), "T".to_string()]));
        lst1.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let mut lst2: Vec<Var> = vec![];
        lst2.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst2.push( generate_var2(1, 0, "chr1".to_string(), 101, vec!["T".to_string(), "TCCC".to_string()]));
        lst2.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst2 = VarList::new(lst2);

        let mut lst3: Vec<Var> = vec![];
        lst3.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string(),]));
        lst3.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["TTT".to_string(), "T".to_string(), "TTCCCT".to_string()]));
        lst3.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let exp = VarList::new(lst3);

        vlst1.combine(&mut vlst2);
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    #[test]
    fn test_varlist_combine_overlapping_insertion_deletion_snv() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["TTT".to_string(), "T".to_string()]));
        lst1.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let mut lst2: Vec<Var> = vec![];
        lst2.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string()]));
        lst2.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["T".to_string(), "G".to_string(), "C".to_string()]));
        lst2.push( generate_var2(2, 0, "chr1".to_string(), 101, vec!["T".to_string(), "TCCC".to_string()]));
        lst2.push( generate_var2(3, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let mut vlst2 = VarList::new(lst2);

        let mut lst3: Vec<Var> = vec![];
        lst3.push( generate_var2(0, 0, "chr1".to_string(), 5, vec!["A".to_string(), "G".to_string(),]));
        lst3.push( generate_var2(1, 0, "chr1".to_string(), 100, vec!["TTT".to_string(), "CTT".to_string(), "GTT".to_string(), "T".to_string(), "TTCCCT".to_string()]));
        lst3.push( generate_var2(2, 0, "chr1".to_string(), 200, vec!["T".to_string(), "G".to_string()]));
        let exp = VarList::new(lst3);

        vlst1.combine(&mut vlst2);
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    #[test]
    fn test_varlist_combine_chr20_bug_case1() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 19, "chr20".to_string(), 5898747, vec!["GAA".to_string(), "G".to_string()]));
        lst1.push( generate_var2(1, 19, "chr20".to_string(), 5898749, vec!["ACT".to_string(), "A".to_string()]));
        lst1.push( generate_var2(2, 19, "chr20".to_string(), 5898751, vec!["T".to_string(), "A".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let mut lst2: Vec<Var> = vec![];
        lst2.push( generate_var2(0, 19, "chr20".to_string(), 5898748, vec!["A".to_string(), "C".to_string()]));
        lst2.push( generate_var2(1, 19, "chr20".to_string(), 5898751, vec!["T".to_string(), "A".to_string()]));
        let mut vlst2 = VarList::new(lst2);

        let mut lst3: Vec<Var> = vec![];
        lst3.push( generate_var2(0, 19, "chr20".to_string(), 5898747, vec!["GAACT".to_string(), "GAA".to_string(), "GAACA".to_string(), "GCACT".to_string(), "GCT".to_string()]));
        let exp = VarList::new(lst3);

        vlst1.combine(&mut vlst2);


        assert_ne!(vlst1.lst.len(), 0);
        assert_eq!(vlst1.lst[0].alleles, exp.lst[0].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    #[test]
    fn test_varlist_combine_chr20_bug_case2() {

        let mut lst1: Vec<Var> = vec![];
        lst1.push( generate_var2(0, 19, "chr20".to_string(), 42449920, vec!["AAAGCTT".to_string(), "A".to_string()]));
        lst1.push( generate_var2(1, 19, "chr20".to_string(), 42449926, vec!["T".to_string(), "TAA".to_string()]));
        lst1.push( generate_var2(2, 19, "chr20".to_string(), 42449926, vec!["T".to_string(), "TAA".to_string()]));
        let mut vlst1 = VarList::new(lst1);

        let mut vlst2 = VarList::new(vec![]);

        let mut lst3: Vec<Var> = vec![];
        lst3.push( generate_var2(0, 19, "chr20".to_string(), 42449920, vec!["AAAGCTT".to_string(), "A".to_string(), "AAAGCTTAA".to_string()]));
        let exp = VarList::new(lst3);

        vlst1.combine(&mut vlst2);
        assert_ne!(vlst1.lst.len(), 0);
        assert_eq!(vlst1.lst[0].alleles, exp.lst[0].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }
}
