//! Data structures to represent variants and haplotype fragments defined over those variants.

use bio::stats::LogProb;
use call_potential_snvs::VARLIST_CAPACITY;
use errors::*;
use genotype_probs::*;
use half::f16;
use hashbrown::HashMap;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bcf;
use rust_htslib::bcf::Read as bcfread;
use std::cmp::Ordering;
use std::convert::From;
use std::fmt;
use util::*;

#[derive(Clone, Copy)]
pub struct FragCall {
    pub frag_ix: u32,  // index into fragment list
    pub var_ix: u32,           // index into variant list
    pub allele: u8,              // allele call
    pub qual: LogProb,           // LogProb probability the call is an error
}

#[derive(Clone)]
pub struct Fragment {
    pub id: Option<String>,
    pub calls: Vec<FragCall>,
    pub p_read_hap: [f16; 2],
    pub reverse_strand: bool
}

#[repr(u8)]
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum VarFilter {
    Pass = 0,
    Density = 1,
    Depth = 2,
    DensityAndDepth = 3,
    StrandBias = 4,
    DensityAndStrandBias = 5,
    DepthAndStrandBias = 6,
    DensityAndDepthAndStrandBias = 7,
}

impl fmt::Display for VarFilter {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            VarFilter::Pass => write!(f, "PASS"),
            VarFilter::Density => write!(f, "dn"),
            VarFilter::Depth => write!(f, "dp"),
            VarFilter::DensityAndDepth => write!(f, "dn;dp"),
            VarFilter::StrandBias => write!(f, "sb"),
            VarFilter::DensityAndStrandBias => write!(f, "dn;sb"),
            VarFilter::DepthAndStrandBias => write!(f, "dp;sb"),
            VarFilter::DensityAndDepthAndStrandBias => write!(f, "dn;dp;sb"),
        }
    }
}

impl From<usize> for VarFilter {
    fn from(item: usize) -> Self {
        match item {
            0 => VarFilter::Pass,
            1 => VarFilter::Density,
            2 => VarFilter::Depth,
            3 => VarFilter::DensityAndDepth,
            4 => VarFilter::StrandBias,
            5 => VarFilter::DensityAndStrandBias,
            6 => VarFilter::DepthAndStrandBias,
            7 => VarFilter::DensityAndDepthAndStrandBias,
            _ => {
                panic!("Invalid value while combining variant filters");
            }
        }
    }
}

impl VarFilter {
    // rhs is filter to add
    pub fn add_filter(&mut self, filter: VarFilter) {
        *self = VarFilter::from(*self as usize | filter as usize)
    }
    pub fn has_filter(&self, filter: VarFilter) -> bool {
        (*self as usize & filter as usize) != 0
    }
}

#[derive(Debug, Clone)]
pub struct Var {
    pub ix: usize,
    pub tid: u32,
    pub pos0: usize,
    pub alleles: Vec<String>, // ref allele is alleles[0] and each that follows is a variant allele
    pub dp: u16,
    // depth of coverage
    pub allele_counts: Vec<u16>, // indices match up with those of Var.alleles
    pub allele_counts_forward: Vec<u16>, // indices match up with those of Var.alleles
    pub allele_counts_reverse: Vec<u16>, // indices match up with those of Var.alleles
    pub ambiguous_count: u16,
    pub qual: f16,
    pub filter: VarFilter, // bitwise flag holding filter info: 0 == PASS, 1 == dp, 2 == dn, 3 == dp && dn
    pub genotype: Genotype,
    //pub unphased: bool, // whether the variant has been explicity flagged as unphased
    pub gq: f16,
    pub unphased_genotype: Genotype,
    pub unphased_gq: f16,
    pub genotype_post: GenotypeProbs, // genotype posteriors[a1][a2] is log posterior of phased a1|a2 haplotype
    // e.g. genotype_posteriors[2][0] is the log posterior probability of 2|0 haplotype
    pub phase_set: Option<usize>,
    pub strand_bias_pvalue: f16, // fisher's exact test strand bias Pvalue
    pub mec: u16,                // mec for variant
    pub mec_frac_variant: f16,   // mec fraction for this variant
    pub mec_frac_block: f16,     // mec fraction for this haplotype block
    pub mean_allele_qual: f16,
    pub dp_any_mq: u16,
    pub mq10_frac: f16,
    pub mq20_frac: f16,
    pub mq30_frac: f16,
    pub mq40_frac: f16,
    pub mq50_frac: f16,
}

impl Var {
    fn longest_allele_len(&self) -> Result<usize> {
        Ok(self
            .alleles
            .iter()
            .map(|x| x.len())
            .max()
            .chain_err(|| "Error obtaining max allele length.")?)
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
    ix: HashMap<u32, Vec<usize>>, // key is TID
    pub target_names: Vec<String>,
}

pub fn parse_vcf_potential_variants(
    vcffile_name: &String,
    bamfile_name: &String,
) -> Result<VarList> {
    // must assert that the VCF file is sorted correctly
    // can we just read it in and then check that it's sorted using the check_sorted function vs the bam's tlist?

    let mut vcf = bcf::Reader::from_path(vcffile_name).chain_err(|| ErrorKind::BCFOpenError)?;
    let vcfh = bcf::Reader::from_path(vcffile_name).chain_err(|| ErrorKind::BCFOpenError)?;

    let target_names = parse_target_names(&bamfile_name)?;
    let bam = bam::Reader::from_path(bamfile_name).chain_err(|| ErrorKind::BamOpenError)?;
    let mut chrom2tid: HashMap<String, usize> = HashMap::new();

    for (t, name) in bam.header().target_names().iter().enumerate() {
        let s: String = u8_to_string(name)?;
        chrom2tid.insert(s, t);
    }

    let mut varlist: Vec<Var> = Vec::with_capacity(VARLIST_CAPACITY);

    for r in vcf.records() {
        let record = r.chain_err(|| ErrorKind::BCFReadError)?;
        // map the VCF rid to chrom name
        // use BAM header to map this to TID
        // fill in the tid, chrom, and pos
        // get the alleles from the vcf record as well

        let rid = record.rid().chain_err(|| "Error accessing vcf RID")?;
        let chrom: String = u8_to_string(vcfh.header().rid2name(rid))?;

        if !chrom2tid.contains_key(&chrom) {
            eprintln!(
                "WARNING: Potential variant VCF contains contig {} not found in BAM contigs.",
                chrom
            );
        }

        let mut alleles: Vec<String> = vec![];
        for a in record.alleles().iter() {
            let s = u8_to_string(a)?;
            alleles.push(s);
        }

        let new_var = Var {
            ix: 0,
            tid: *chrom2tid
                .get(&chrom)
                .chain_err(|| "Error accessing tid from chrom2tid data structure")?
                as u32,
            pos0: record.pos() as usize,
            alleles: alleles.clone(),
            dp: 0,
            allele_counts: vec![0; alleles.len()],
            allele_counts_forward: vec![0; alleles.len()],
            allele_counts_reverse: vec![0; alleles.len()],
            ambiguous_count: 0,
            qual: f16::from_f64(0.0),
            filter: VarFilter::Pass,
            genotype: Genotype(0, 0),
            //unphased: false,
            gq: f16::from_f64(0.0),
            unphased_genotype: Genotype(0, 0),
            unphased_gq: f16::from_f64(0.0),
            genotype_post: GenotypeProbs::uniform(alleles.len()),
            phase_set: None,
            strand_bias_pvalue: f16::from_f64(0.0),
            mec: 0,
            mec_frac_variant: f16::from_f64(0.0), // mec fraction for this variant
            mec_frac_block: f16::from_f64(0.0),   // mec fraction for this haplotype block
            mean_allele_qual: f16::from_f64(0.0),
            dp_any_mq: 0,
            mq10_frac: f16::from_f64(0.0),
            mq20_frac: f16::from_f64(0.0),
            mq30_frac: f16::from_f64(0.0),
            mq40_frac: f16::from_f64(0.0),
            mq50_frac: f16::from_f64(0.0),
        };
        varlist.push(new_var);
    }

    let vlst = VarList::new(varlist, target_names.clone())?;
    vlst.assert_sorted();

    Ok(vlst)
}

/*
impl Index<VarIx> for VarList {
    type Output = Var;

    fn index(&self, var_ix: VarIx) -> &Var {
        &self.lst[var_ix.ix]
    }
}

impl Index<Range<VarIx>> for VarList {
    type Output = VarList;

    fn index(&self, r: Range<VarIx>) -> &Var {
        &self.lst[r.start.ix..r.end.ix]
    }
}*/

impl VarList {
    pub fn new(lst: Vec<Var>, target_names: Vec<String>) -> Result<VarList> {
        let mut v = VarList {
            lst: lst,
            ix: HashMap::new(),
            target_names: target_names,
        };
        v.sort()?;
        Ok(v)
    }

    pub fn len(&self) -> usize {
        self.lst.len()
    }

    pub fn sort(&mut self) -> Result<()> {
        self.lst.sort();
        self.add_ix();
        self.ix.clear();
        self.index_lst()?;
        Ok(())
    }

    pub fn assert_sorted(&self) {
        if self.lst.len() == 0 {
            return;
        }
        for i in 0..self.lst.len() - 1 {
            assert!(
                (self.lst[i].tid < self.lst[i + 1].tid)
                    || (self.lst[i].tid == self.lst[i + 1].tid
                        && self.lst[i].pos0 <= self.lst[i + 1].pos0)
            );
            assert_eq!(self.lst[i].ix, i);
            assert_eq!(self.lst[i + 1].ix, i + 1);
        }
    }

    fn add_ix(&mut self) {
        for i in 0..self.lst.len() {
            self.lst[i].ix = i;
        }
    }

    fn index_lst(&mut self) -> Result<()> {
        // need to throw error if list isn't sorted
        self.assert_sorted();
        // for every chromosome, get the position of the last variant
        // and the varlist index of the first variant on that chromosome
        let mut max_positions: HashMap<u32, usize> = HashMap::new();
        let mut first_ix: HashMap<u32, usize> = HashMap::new();
        for (i, var) in self.lst.iter().enumerate() {
            let mpos: usize = *(max_positions.entry(var.tid).or_insert(0));
            if var.pos0 > mpos {
                max_positions.insert(var.tid, var.pos0);
            }
            if !first_ix.contains_key(&var.tid) {
                first_ix.insert(var.tid, i); // insert varlist index of first variant on chrom
            }
        }

        // for every chrom, iterate up to its last potential variant and create an index that
        // returns the lst index of the next SNV for some position mod 1000
        for (tid, max_pos) in &max_positions {
            let mut v: Vec<usize> = vec![];
            // the first variant after position 0 is the first variant on the chromosome
            let e = match first_ix.get(tid) {
                Some(x) => *x,
                None => {
                    bail!("This dictionary is missing a chromosome!");
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
            self.ix.insert(*tid, v); // insert the index vector into the ix dictionary
        }
        Ok(())
    }

    pub fn get_variants_range(&self, interval: GenomicInterval) -> Result<Vec<Var>> {
        // vector of variants to fill and return
        let mut vlst: Vec<Var> = vec![];
        // get the varlist index of a nearby position on the left

        let index_pos = (interval.start_pos as usize) / INDEX_FREQ;

        if index_pos
            >= self
                .ix
                .get(&interval.tid)
                .chain_err(|| {
                    format!(
                        "Error accessing chromosome {} from variant list index.",
                        &interval.chrom
                    )
                })?
                .len()
        {
            return Ok(vlst);
        }

        let mut i = self.ix.get(&interval.tid).chain_err(|| {
            format!(
                "Error accessing chromosome {} from variant list index.",
                &interval.chrom
            )
        })?[index_pos];

        while i < self.lst.len()
            && self.lst[i].tid == interval.tid
            && self.lst[i].pos0 + self.lst[i].longest_allele_len()? <= interval.end_pos as usize
        {
            if self.lst[i].pos0 >= interval.start_pos as usize {
                vlst.push(self.lst[i].clone());
            }
            i += 1;
        }

        for var in &vlst {
            assert!(var.tid == interval.tid);
            assert!(var.pos0 >= interval.start_pos as usize);
            assert!(var.pos0 <= interval.end_pos as usize);
        }

        Ok(vlst)
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
                || (v.pos0 == min_pos0 && v.alleles[0].len() > min_pos0_refseq.len())
            {
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
            for (a, r) in (vs..ve).enumerate() {
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
        for v in var_group.iter() {
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
                let suffix_seq: String = longest_ref
                    .chars()
                    .skip(v.alleles[0].len())
                    .take(diff)
                    .collect();
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
        new_v.genotype = Genotype(0, 0);
        new_v.gq = f16::from_f64(0.0);
        new_v.genotype_post = GenotypeProbs::uniform(new_v.alleles.len());
        new_v.phase_set = None;

        new_v
    }

    pub fn combine(&mut self, other: &mut VarList) -> Result<()> {
        if self.target_names != other.target_names {
            bail!("Target names of variant lists that are being combined are not the same.");
        }

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

        self.sort()?;
        self.assert_sorted();

        // clear out the other VarList since we've mutated it beyond saving
        other.lst.clear();
        other.ix.clear();
        Ok(())
    }
}

pub fn var_filter(
    varlist: &mut VarList,
    density_qual: f64,
    density_dist: usize,
    density_count: usize,
    max_depth: u32,
) {
    for i in 0..varlist.lst.len() {
        if varlist.lst[i].qual < f16::from_f64(density_qual) {
            continue;
        }

        let mut count = 0;
        for j in i + 1..varlist.lst.len() {
            if varlist.lst[j].pos0 - varlist.lst[i].pos0 > density_dist {
                break;
            }
            if varlist.lst[j].qual < f16::from_f64(density_qual) {
                continue;
            }
            count += 1;
            if count > density_count {
                for k in i..j + 1 {
                    varlist.lst[k].filter.add_filter(VarFilter::Density);
                }
            }
        }
    }

    for i in 0..varlist.lst.len() {
        if varlist.lst[i].dp > max_depth as u16 {
            varlist.lst[i].filter.add_filter(VarFilter::Depth);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_varfilter_cast() {
        assert_eq!(VarFilter::Pass as usize, 0);
        assert_eq!(VarFilter::Density as usize, 1);
        assert_eq!(VarFilter::Depth as usize, 2);
        assert_eq!(VarFilter::DensityAndDepth as usize, 3);
        assert_eq!(VarFilter::StrandBias as usize, 4);
    }

    #[test]
    fn test_varfilter_enum() {
        let pass = VarFilter::Pass;
        let dp = VarFilter::Depth;
        let dn = VarFilter::Density;
        let sb = VarFilter::StrandBias;

        let mut f1 = VarFilter::Pass;
        assert!(!f1.has_filter(dp));
        f1.add_filter(pass);
        assert_eq!(f1, pass);
        f1.add_filter(dp);
        assert!(f1.has_filter(dp));
        assert!(!f1.has_filter(dn));
        assert_eq!(f1, dp);
        f1.add_filter(dn);
        assert!(f1.has_filter(dp));
        assert!(f1.has_filter(dn));
        assert!(!f1.has_filter(sb));
        assert_eq!(f1, VarFilter::DensityAndDepth);
        f1.add_filter(sb);
        assert_eq!(f1, VarFilter::DensityAndDepthAndStrandBias);
        assert!(f1.has_filter(dp));
        assert!(f1.has_filter(dn));
        assert!(f1.has_filter(sb));
    }

    /**********************************************************************************************/
    // TEST VARIANT RANGE LOOKUP
    /**********************************************************************************************/

    fn pos_alleles_eq(vlst1: Vec<Var>, vlst2: Vec<Var>) -> bool {
        for (v1, v2) in vlst1.iter().zip(vlst2.iter()) {
            if v1 != v2 || v1.tid != v2.tid || v1.alleles != v2.alleles || v1.ix != v2.ix {
                return false;
            }
        }
        true
    }

    fn varlist_pos_alleles_eq(vlst1: VarList, vlst2: VarList) -> bool {
        pos_alleles_eq(vlst1.lst, vlst2.lst)
    }

    fn generate_var1(ix: usize, tid: usize, pos0: usize, ra: String, aa: String) -> Var {
        Var {
            ix: ix,
            tid: tid as u32,
            pos0: pos0,
            alleles: vec![ra, aa],
            dp: 40,
            allele_counts: vec![20, 20],
            allele_counts_forward: vec![10, 10],
            allele_counts_reverse: vec![10, 10],
            ambiguous_count: 0,
            qual: f16::from_f64(0.0),
            filter: VarFilter::Pass,
            genotype: Genotype(0, 1),
            gq: f16::from_f64(0.0),
            mean_allele_qual: f16::from_f64(0.0),
            mec: 0 as u16,
            strand_bias_pvalue: f16::from_f64(0.0),
            mec_frac_block: f16::from_f64(0.0),
            mec_frac_variant: f16::from_f64(0.0),
            dp_any_mq: 40 as u16,
            mq10_frac: f16::from_f64(1.0),
            mq20_frac: f16::from_f64(1.0),
            mq30_frac: f16::from_f64(1.0),
            mq40_frac: f16::from_f64(1.0),
            mq50_frac: f16::from_f64(1.0),
            unphased_genotype: Genotype(0, 1),
            unphased_gq: f16::from_f64(0.0),
            genotype_post: GenotypeProbs::uniform(2),
            phase_set: None,
        }
    }

    fn generate_test_lst1() -> VarList {
        let mut lst: Vec<Var> = vec![];

        lst.push(generate_var1(0, 0, 5, "A".to_string(), "G".to_string()));
        lst.push(generate_var1(1, 0, 1000, "T".to_string(), "A".to_string()));
        lst.push(generate_var1(2, 0, 2005, "T".to_string(), "G".to_string()));
        lst.push(generate_var1(3, 0, 2900, "C".to_string(), "G".to_string()));
        lst.push(generate_var1(4, 0, 6000, "C".to_string(), "A".to_string()));
        lst.push(generate_var1(5, 0, 10000, "C".to_string(), "A".to_string()));
        lst.push(generate_var1(6, 1, 5, "A".to_string(), "G".to_string()));
        lst.push(generate_var1(7, 1, 1000, "T".to_string(), "A".to_string()));
        lst.push(generate_var1(8, 1, 2005, "T".to_string(), "G".to_string()));
        lst.push(generate_var1(9, 1, 2900, "C".to_string(), "G".to_string()));
        lst.push(generate_var1(10, 1, 6000, "C".to_string(), "A".to_string()));
        lst.push(generate_var1(
            11,
            1,
            10000,
            "C".to_string(),
            "A".to_string(),
        ));
        lst.push(generate_var1(
            12,
            2,
            20200,
            "C".to_string(),
            "G".to_string(),
        ));
        lst.push(generate_var1(
            13,
            2,
            25100,
            "A".to_string(),
            "C".to_string(),
        ));
        lst.push(generate_var1(
            14,
            2,
            30400,
            "C".to_string(),
            "A".to_string(),
        ));

        VarList::new(
            lst,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap()
    }

    #[test]
    fn test_varlist_get_variants_range1() {
        let vlst = generate_test_lst1();
        let t = 0;
        let c = "chr1".to_string();
        let p1 = 2500;
        let p2 = 8000;
        let interval = GenomicInterval {
            tid: t,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval).unwrap();

        let mut exp: Vec<Var> = vec![];
        exp.push(generate_var1(3, 0, 2900, "C".to_string(), "G".to_string()));
        exp.push(generate_var1(4, 0, 6000, "C".to_string(), "A".to_string()));

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range2() {
        let vlst = generate_test_lst1();
        let t = 1;
        let c: String = "chr1".to_string();
        let p1 = 0;
        let p2 = 3000;
        let interval = GenomicInterval {
            tid: t,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval).unwrap();

        let mut exp: Vec<Var> = vec![];
        exp.push(generate_var1(6, 1, 5, "A".to_string(), "G".to_string()));
        exp.push(generate_var1(7, 1, 1000, "T".to_string(), "A".to_string()));
        exp.push(generate_var1(8, 1, 2005, "T".to_string(), "G".to_string()));
        exp.push(generate_var1(9, 1, 2900, "C".to_string(), "G".to_string()));

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range3() {
        let vlst = generate_test_lst1();
        let t = 1;
        let c = "chr2".to_string();
        let p1 = 6000;
        let p2 = 10000;
        let interval = GenomicInterval {
            tid: t,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval).unwrap();

        let mut exp: Vec<Var> = vec![];
        exp.push(generate_var1(10, 1, 6000, "C".to_string(), "A".to_string()));
        exp.push(generate_var1(
            11,
            1,
            10000,
            "C".to_string(),
            "A".to_string(),
        ));

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range4() {
        let vlst = generate_test_lst1();
        let t = 2;
        let c = "chr3".to_string();
        let p1 = 20100;
        let p2 = 20200;
        let interval = GenomicInterval {
            tid: t,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval).unwrap();

        let mut exp: Vec<Var> = vec![];
        exp.push(generate_var1(
            12,
            2,
            20200,
            "C".to_string(),
            "G".to_string(),
        ));

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range5() {
        let vlst = generate_test_lst1();
        let t = 2;
        let c = "chr3".to_string();
        let p1 = 20200;
        let p2 = 20200;
        let interval = GenomicInterval {
            tid: t,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval).unwrap();

        let mut exp: Vec<Var> = vec![];
        exp.push(generate_var1(
            12,
            2,
            20200,
            "C".to_string(),
            "G".to_string(),
        ));

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range6() {
        let vlst = generate_test_lst1();
        let t = 2;
        let c = "chr3".to_string();
        let p1 = 25000;
        let p2 = 30500;
        let interval = GenomicInterval {
            tid: t,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval).unwrap();

        let mut exp: Vec<Var> = vec![];

        exp.push(generate_var1(
            13,
            2,
            25100,
            "A".to_string(),
            "C".to_string(),
        ));
        exp.push(generate_var1(
            14,
            2,
            30400,
            "C".to_string(),
            "A".to_string(),
        ));

        assert!(pos_alleles_eq(vars, exp));
    }

    #[test]
    fn test_varlist_get_variants_range7() {
        let vlst = generate_test_lst1();
        let t = 2;
        let c = "chr3".to_string();
        let p1 = 100000;
        let p2 = 200000;
        let interval = GenomicInterval {
            tid: t,
            chrom: c,
            start_pos: p1,
            end_pos: p2,
        };
        let vars = vlst.get_variants_range(interval).unwrap();

        let exp: Vec<Var> = vec![];
        assert!(pos_alleles_eq(vars, exp));
    }

    /**********************************************************************************************/
    // TEST VARIANT SORTING
    /**********************************************************************************************/

    // generate an unsorted version of test_lst1, with meaningless indices. for testing sorting.
    fn generate_test_lst1_unsorted1() -> VarList {
        let mut lst: Vec<Var> = vec![];

        lst.push(generate_var1(0, 1, 2005, "T".to_string(), "G".to_string()));
        lst.push(generate_var1(
            1000,
            1,
            1000,
            "T".to_string(),
            "A".to_string(),
        ));
        lst.push(generate_var1(0, 0, 2900, "C".to_string(), "G".to_string()));
        lst.push(generate_var1(100, 0, 5, "A".to_string(), "G".to_string()));
        lst.push(generate_var1(4, 0, 6000, "C".to_string(), "A".to_string()));
        lst.push(generate_var1(
            14,
            2,
            30400,
            "C".to_string(),
            "A".to_string(),
        ));
        lst.push(generate_var1(
            13,
            2,
            25100,
            "A".to_string(),
            "C".to_string(),
        ));
        lst.push(generate_var1(
            10,
            0,
            10000,
            "C".to_string(),
            "A".to_string(),
        ));
        lst.push(generate_var1(1, 1, 5, "A".to_string(), "G".to_string()));
        lst.push(generate_var1(2, 0, 2005, "T".to_string(), "G".to_string()));
        lst.push(generate_var1(10, 1, 6000, "C".to_string(), "A".to_string()));
        lst.push(generate_var1(10, 1, 2900, "C".to_string(), "G".to_string()));
        lst.push(generate_var1(
            11,
            1,
            10000,
            "C".to_string(),
            "A".to_string(),
        ));
        lst.push(generate_var1(2, 0, 1000, "T".to_string(), "A".to_string()));
        lst.push(generate_var1(
            12,
            2,
            20200,
            "C".to_string(),
            "G".to_string(),
        ));

        VarList::new(
            lst,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap()
    }

    // generate another unsorted version of test_lst1, with meaningless indices. for testing sorting.
    fn generate_test_lst1_unsorted2() -> VarList {
        let mut lst: Vec<Var> = vec![];

        lst.push(generate_var1(10, 1, 5, "A".to_string(), "G".to_string()));
        lst.push(generate_var1(1, 0, 6000, "C".to_string(), "A".to_string()));
        lst.push(generate_var1(0, 1, 2005, "T".to_string(), "G".to_string()));
        lst.push(generate_var1(0, 0, 5, "A".to_string(), "G".to_string()));
        lst.push(generate_var1(
            1000,
            0,
            1000,
            "T".to_string(),
            "A".to_string(),
        ));
        lst.push(generate_var1(0, 1, 1000, "T".to_string(), "A".to_string()));
        lst.push(generate_var1(
            12,
            2,
            20200,
            "C".to_string(),
            "G".to_string(),
        ));
        lst.push(generate_var1(2, 2, 25100, "A".to_string(), "C".to_string()));
        lst.push(generate_var1(
            100,
            0,
            2900,
            "C".to_string(),
            "G".to_string(),
        ));
        lst.push(generate_var1(9, 0, 2005, "T".to_string(), "G".to_string()));
        lst.push(generate_var1(5, 0, 10000, "C".to_string(), "A".to_string()));
        lst.push(generate_var1(3, 1, 2900, "C".to_string(), "G".to_string()));
        lst.push(generate_var1(10, 1, 6000, "C".to_string(), "A".to_string()));
        lst.push(generate_var1(
            50,
            2,
            30400,
            "C".to_string(),
            "A".to_string(),
        ));
        lst.push(generate_var1(3, 1, 10000, "C".to_string(), "A".to_string()));

        VarList::new(
            lst,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap()
    }

    #[test]
    fn test_varlist_sort1() {
        let mut vlst_unsorted = generate_test_lst1_unsorted1();

        vlst_unsorted.sort().unwrap();
        vlst_unsorted.assert_sorted();
    }

    #[test]
    fn test_varlist_sort2() {
        let mut vlst_unsorted = generate_test_lst1_unsorted2();

        vlst_unsorted.sort().unwrap();
        vlst_unsorted.assert_sorted();
    }

    /**********************************************************************************************/
    // TEST COMBINING VARLISTS
    /**********************************************************************************************/

    fn generate_var2(ix: usize, tid: usize, pos0: usize, alleles: Vec<String>) -> Var {
        Var {
            ix: ix,
            tid: tid as u32,
            pos0: pos0,
            alleles: alleles,
            dp: 40,
            allele_counts: vec![20, 20],
            allele_counts_forward: vec![10, 10],
            allele_counts_reverse: vec![10, 10],
            ambiguous_count: 0,
            qual: f16::from_f64(0.0),
            filter: VarFilter::Pass,
            genotype: Genotype(0, 1),
            gq: f16::from_f64(0.0),
            mean_allele_qual: f16::from_f64(0.0),
            strand_bias_pvalue: f16::from_f64(0.0),
            mec: 0,
            mec_frac_block: f16::from_f64(0.0),
            mec_frac_variant: f16::from_f64(0.0),
            dp_any_mq: 40,
            mq10_frac: f16::from_f64(1.0),
            mq20_frac: f16::from_f64(1.0),
            mq30_frac: f16::from_f64(1.0),
            mq40_frac: f16::from_f64(1.0),
            mq50_frac: f16::from_f64(1.0),
            unphased_genotype: Genotype(0, 1),
            unphased_gq: f16::from_f64(0.0),
            genotype_post: GenotypeProbs::uniform(2),
            phase_set: None,
        }
    }

    // combining two identical lists should result in the same list
    #[test]
    fn test_varlist_combine_same_lists() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst1.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "A".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let vlst1_bak = vlst1.clone();

        let mut lst2: Vec<Var> = vec![];
        lst2.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst2.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "A".to_string()],
        ));
        lst2.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst2 = VarList::new(
            lst2,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();
        assert!(varlist_pos_alleles_eq(vlst1, vlst1_bak));
    }

    // combine varlists with different variants into essentially the union of the variant set
    #[test]
    fn test_varlist_combine_not_same() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "A".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst2: Vec<Var> = vec![];
        lst2.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst2.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "A".to_string()],
        ));
        let mut vlst2 = VarList::new(
            lst2,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst3: Vec<Var> = vec![];
        lst3.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst3.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "A".to_string()],
        ));
        lst3.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let exp = VarList::new(
            lst3,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // combine a variant where a SNVs overlaps a deletion
    #[test]
    fn test_varlist_combine_deletion() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst1.push(generate_var2(
            1,
            0,
            100,
            vec!["TTT".to_string(), "T".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst2: Vec<Var> = vec![];
        lst2.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst2.push(generate_var2(
            1,
            0,
            101,
            vec!["T".to_string(), "C".to_string()],
        ));
        lst2.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst2 = VarList::new(
            lst2,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst3: Vec<Var> = vec![];
        lst3.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst3.push(generate_var2(
            1,
            0,
            100,
            vec!["TTT".to_string(), "T".to_string(), "TCT".to_string()],
        ));
        lst3.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let exp = VarList::new(
            lst3,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // combine a variant where a SNVs overlaps an insertion
    #[test]
    fn test_varlist_combine_insertion() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst1.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "C".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst2: Vec<Var> = vec![];
        lst2.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst2.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "TAAAA".to_string()],
        ));
        lst2.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst2 = VarList::new(
            lst2,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst3: Vec<Var> = vec![];
        lst3.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst3.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "C".to_string(), "TAAAA".to_string()],
        ));
        lst3.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let exp = VarList::new(
            lst3,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // make sure that if the insertion is only adjacent, it is left separate from SNV.
    #[test]
    fn test_varlist_combine_nearby_insertion() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst1.push(generate_var2(
            1,
            0,
            101,
            vec!["T".to_string(), "C".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst2: Vec<Var> = vec![];
        lst2.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst2.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "TAAAA".to_string()],
        ));
        lst2.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst2 = VarList::new(
            lst2,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst3: Vec<Var> = vec![];
        lst3.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst3.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "TAAAA".to_string()],
        ));
        lst3.push(generate_var2(
            2,
            0,
            101,
            vec!["T".to_string(), "C".to_string()],
        ));
        lst3.push(generate_var2(
            3,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let exp = VarList::new(
            lst3,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert_eq!(vlst1.lst[2].alleles, exp.lst[2].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // combine a variant where multiple SNVs overlap the same deletion
    #[test]
    fn test_varlist_combine_deletion_multiple_snp() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst1.push(generate_var2(
            1,
            0,
            100,
            vec!["TTT".to_string(), "T".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst2: Vec<Var> = vec![];
        lst2.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst2.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "C".to_string()],
        ));
        lst2.push(generate_var2(
            2,
            0,
            102,
            vec!["T".to_string(), "C".to_string()],
        ));
        lst2.push(generate_var2(
            3,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst2 = VarList::new(
            lst2,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst3: Vec<Var> = vec![];
        lst3.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst3.push(generate_var2(
            1,
            0,
            100,
            vec![
                "TTT".to_string(),
                "CTT".to_string(),
                "T".to_string(),
                "TTC".to_string(),
            ],
        ));
        lst3.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let exp = VarList::new(
            lst3,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // combine a variant where multiple SNVs overlap the same deletion
    #[test]
    fn test_varlist_combine_overlapping_deletions() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst1.push(generate_var2(
            1,
            0,
            100,
            vec!["TTT".to_string(), "T".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst2: Vec<Var> = vec![];
        lst2.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst2.push(generate_var2(
            1,
            0,
            102,
            vec!["TCC".to_string(), "T".to_string()],
        ));
        lst2.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst2 = VarList::new(
            lst2,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst3: Vec<Var> = vec![];
        lst3.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst3.push(generate_var2(
            1,
            0,
            100,
            vec!["TTTCC".to_string(), "TCC".to_string(), "TTT".to_string()],
        ));
        lst3.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let exp = VarList::new(
            lst3,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    // combine a variant where multiple SNVs overlap the same deletion
    #[test]
    fn test_varlist_combine_overlapping_insertion_deletion() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst1.push(generate_var2(
            1,
            0,
            100,
            vec!["TTT".to_string(), "T".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst2: Vec<Var> = vec![];
        lst2.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst2.push(generate_var2(
            1,
            0,
            101,
            vec!["T".to_string(), "TCCC".to_string()],
        ));
        lst2.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst2 = VarList::new(
            lst2,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst3: Vec<Var> = vec![];
        lst3.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst3.push(generate_var2(
            1,
            0,
            100,
            vec!["TTT".to_string(), "T".to_string(), "TTCCCT".to_string()],
        ));
        lst3.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let exp = VarList::new(
            lst3,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    #[test]
    fn test_varlist_combine_overlapping_insertion_deletion_snv() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst1.push(generate_var2(
            1,
            0,
            100,
            vec!["TTT".to_string(), "T".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst2: Vec<Var> = vec![];
        lst2.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst2.push(generate_var2(
            1,
            0,
            100,
            vec!["T".to_string(), "G".to_string(), "C".to_string()],
        ));
        lst2.push(generate_var2(
            2,
            0,
            101,
            vec!["T".to_string(), "TCCC".to_string()],
        ));
        lst2.push(generate_var2(
            3,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let mut vlst2 = VarList::new(
            lst2,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst3: Vec<Var> = vec![];
        lst3.push(generate_var2(
            0,
            0,
            5,
            vec!["A".to_string(), "G".to_string()],
        ));
        lst3.push(generate_var2(
            1,
            0,
            100,
            vec![
                "TTT".to_string(),
                "CTT".to_string(),
                "GTT".to_string(),
                "T".to_string(),
                "TTCCCT".to_string(),
            ],
        ));
        lst3.push(generate_var2(
            2,
            0,
            200,
            vec!["T".to_string(), "G".to_string()],
        ));
        let exp = VarList::new(
            lst3,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();
        assert_eq!(vlst1.lst[1].alleles, exp.lst[1].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    #[test]
    fn test_varlist_combine_chr20_bug_case1() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            0,
            19,
            5898747,
            vec!["GAA".to_string(), "G".to_string()],
        ));
        lst1.push(generate_var2(
            1,
            19,
            5898749,
            vec!["ACT".to_string(), "A".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            19,
            5898751,
            vec!["T".to_string(), "A".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst2: Vec<Var> = vec![];
        lst2.push(generate_var2(
            0,
            19,
            5898748,
            vec!["A".to_string(), "C".to_string()],
        ));
        lst2.push(generate_var2(
            1,
            19,
            5898751,
            vec!["T".to_string(), "A".to_string()],
        ));
        let mut vlst2 = VarList::new(
            lst2,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst3: Vec<Var> = vec![];
        lst3.push(generate_var2(
            0,
            19,
            5898747,
            vec![
                "GAACT".to_string(),
                "GAA".to_string(),
                "GAACA".to_string(),
                "GCACT".to_string(),
                "GCT".to_string(),
            ],
        ));
        let exp = VarList::new(
            lst3,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();

        assert_ne!(vlst1.lst.len(), 0);
        assert_eq!(vlst1.lst[0].alleles, exp.lst[0].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }

    #[test]
    fn test_varlist_combine_chr20_bug_case2() {
        let mut lst1: Vec<Var> = vec![];
        lst1.push(generate_var2(
            0,
            19,
            42449920,
            vec!["AAAGCTT".to_string(), "A".to_string()],
        ));
        lst1.push(generate_var2(
            1,
            19,
            42449926,
            vec!["T".to_string(), "TAA".to_string()],
        ));
        lst1.push(generate_var2(
            2,
            19,
            42449926,
            vec!["T".to_string(), "TAA".to_string()],
        ));
        let mut vlst1 = VarList::new(
            lst1,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut vlst2 = VarList::new(
            vec![],
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        let mut lst3: Vec<Var> = vec![];
        lst3.push(generate_var2(
            0,
            19,
            42449920,
            vec![
                "AAAGCTT".to_string(),
                "A".to_string(),
                "AAAGCTTAA".to_string(),
            ],
        ));
        let exp = VarList::new(
            lst3,
            vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()],
        )
        .unwrap();

        vlst1.combine(&mut vlst2).unwrap();
        assert_ne!(vlst1.lst.len(), 0);
        assert_eq!(vlst1.lst[0].alleles, exp.lst[0].alleles);
        assert!(varlist_pos_alleles_eq(vlst1, exp));
    }
}
