use bio::stats::{PHREDProb, LogProb, Prob};
use util::{VarList, Fragment, FragCall, GenomicInterval};
use haplotype_assembly::{generate_flist_buffer, call_hapcut2};
use std::collections::HashMap;
use std::error::Error;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use chrono::prelude::*;

//use std::ops::range;
use rand::{Rng,StdRng,SeedableRng};

static MAX_P_MISCALL_F64: f64 = 0.2;
static MIN_GQ_FOR_PHASING: f64 = 50.0;

fn generate_simple_pileup(flist: &Vec<Fragment>, n_var: usize) -> Vec<Vec<(char,LogProb)>> {
    let mut pileup_lst: Vec<Vec<(char,LogProb)>> = vec![vec![]; n_var];

    for fragment in flist {
        for call in fragment.clone().calls {
            //if call.qual < LogProb::from(Prob(MAX_P_MISCALL_F64)) {
            pileup_lst[call.var_ix].push((call.allele, call.qual));
            //}
        }
    }
    pileup_lst
}


fn generate_fragcall_pileup(flist: &Vec<Fragment>, n_var: usize) -> Vec<Vec<FragCall>> {
    let mut pileup_lst: Vec<Vec<FragCall>> = vec![vec![]; n_var];

    for fragment in flist {
        for call in fragment.clone().calls {
            //if call.qual < LogProb::from(Prob(MAX_P_MISCALL_F64)) {
            pileup_lst[call.var_ix].push(call);
           // }
        }
    }
    pileup_lst
}

pub fn estimate_genotype_priors() -> HashMap<(char, (char, char)), LogProb>{
    // estimate prior probability of genotypes using strategy described here:
    // http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694485/
    // "prior probability of each genotype"
    let het_snp_rate = LogProb::from(Prob(0.0005));
    let hom_snp_rate = LogProb::from(Prob(0.001));

    // key of diploid_genotype_priors is (char,(char,char)) (ref_allele, G=(allele1,allele2))
    // key of haploid priors is (char, char) which is (ref_allele, allele1)
    let mut diploid_genotype_priors: HashMap<(char,(char,char)), LogProb> = HashMap::new();
    let mut haploid_genotype_priors: HashMap<(char,char), LogProb> = HashMap::new();

    // define base transitions
    let mut transition: HashMap<char, char> = HashMap::new();
    transition.insert('A','G');
    transition.insert('G','A');
    transition.insert('T','C');
    transition.insert('C','T');

    let alleles: Vec<char> = vec!['A','C','G','T'];
    let genotypes: Vec<(char,char)> = vec![('A', 'A'), ('A', 'C'), ('A', 'G'), ('A', 'T'),
                                       ('C', 'C'), ('C', 'G'), ('C', 'T'),
                                       ('G', 'G'), ('G', 'T'),
                                       ('T', 'T')];


    for aref in &alleles {
        let allele = *aref;
        // priors on haploid alleles
        haploid_genotype_priors.insert((allele,allele), LogProb::ln_one_minus_exp(&hom_snp_rate));
        haploid_genotype_priors.insert((allele,*transition.get(&allele).unwrap()), hom_snp_rate + LogProb::from(Prob(4.0/6.0)));

        for tref in &alleles{
            let transversion = *tref;
            if haploid_genotype_priors.contains_key(&(allele,transversion)){
                continue;
            }
            haploid_genotype_priors.insert((allele,transversion), hom_snp_rate + LogProb::from(Prob(1.0/6.0)));

        }

        for gt in &genotypes {
            let (g1, g2) = *gt;
            // probability of homozygous reference is the probability of neither het or hom SNP
            if g1 == g2 && g1 == allele {
                let snp_rate = LogProb::ln_add_exp(het_snp_rate, hom_snp_rate);
                let one_minus_snp_rate = LogProb::ln_one_minus_exp(&snp_rate);
                diploid_genotype_priors.insert((allele,*gt), one_minus_snp_rate);
            } else if g1 == g2 && g1 != allele {
                // transitions are 4 times as likely as transversions
                if g1 == *transition.get(&allele).unwrap() {
                    diploid_genotype_priors.insert((allele,*gt), het_snp_rate + LogProb::from(Prob(4.0/6.0)));
                }else {
                    diploid_genotype_priors.insert((allele,*gt), het_snp_rate + LogProb::from(Prob(1.0/6.0)));
                }
            } else { // else it's the product of the haploid priors
                diploid_genotype_priors.insert((allele,*gt), *haploid_genotype_priors.get(&(allele,g1)).unwrap() +
                                                     *(haploid_genotype_priors.get(&(allele,g2))).unwrap());
            }
        }
    }

    diploid_genotype_priors
}

fn count_alleles(pileup: &Vec<FragCall>) -> (usize, usize, usize) {
    // compute probabilities of data given each genotype
    let mut count0 = 0;
    let mut count1 = 0;
    let mut count_amb = 0; // ambiguous
    let ln_max_p_miscall = LogProb::from(Prob(MAX_P_MISCALL_F64));

    for call in pileup {
        match call.allele {
            '0' => {
                if call.qual < ln_max_p_miscall {
                    //println!("{}",*Prob::from(call.qual));
                    count0 += 1;
                } else {
                    count_amb += 1;
                }
            }
            '1' => {
                if call.qual < ln_max_p_miscall {
                    //println!("{}",*Prob::from(call.qual));
                    count1 += 1;
                } else {
                    count_amb += 1;
                }
            }
            _ => {
                panic!("Unexpected allele observed in pileup.");
            }
        }
    }

    (count0, count1, count_amb)
}

struct HapPost {
    post00: LogProb,
    post01: LogProb,
    post10: LogProb,
    post11: LogProb,
}

pub fn calculate_genotypes_without_haplotypes(pileup: &Vec<(char, LogProb)>,
                                          genotype_priors: &HashMap<(char, (char, char)),LogProb>,
                                          ref_allele: &String,
                                          var_allele: &String) -> (LogProb, LogProb, LogProb) {

    let ln_max_p_miscall = LogProb::from(Prob(MAX_P_MISCALL_F64));
    let ln_half = LogProb::from(Prob(0.5));
    let mut priors = vec![LogProb::from(Prob(0.25)); 4];

    for g in 0..4 {
        if ref_allele.len() == 1 && var_allele.len() == 1 {

            let ra = ref_allele.chars().nth(0).unwrap();
            let g1 = if g == 0 || g == 1 { ref_allele.chars().nth(0).unwrap() } else { var_allele.chars().nth(0).unwrap() };
            let g2 = if g == 0 || g == 2 { ref_allele.chars().nth(0).unwrap() } else { var_allele.chars().nth(0).unwrap() };

            match genotype_priors.get(&(ra, (g1, g2))) {
                Some(p) => {
                    priors[g] = if g == 1 || g == 2 {ln_half+*p} else {*p};
                },
                None => {
                    match genotype_priors.get(&(ra, (g2, g1))) {
                        Some(p) => {
                            priors[g] = if g == 1 || g == 2 {ln_half+*p} else {*p};
                        },
                        None => { println!("{} ({},{})",ra,g2,g1); panic!("Genotype not in genotype priors."); }
                    };
                }
            }
        } else {
            // we just assume that the rate of indels is 1e-4
            priors[0] = LogProb::from(Prob(0.9999));
            priors[1] = LogProb::from(Prob(0.0001 / 3.0));
            priors[2] = LogProb::from(Prob(0.0001 / 3.0));
            priors[3] = LogProb::from(Prob(0.0001 / 3.0));
        }
    }

    // compute probabilities of data given each genotype
    let mut prob00 = priors[0];
    let mut prob01 = LogProb::ln_add_exp(priors[1], priors[2]);
    let mut prob11 = priors[3];

    for &(allele, qual) in pileup {

        if qual >= ln_max_p_miscall {
            continue;
        }

        let p_call = LogProb::ln_sub_exp(LogProb::ln_one(), qual);
        let p_miscall = qual;

        match allele {
            '0' => {
                prob00 = prob00 + p_call;
                prob01 = prob01 + LogProb::ln_add_exp(ln_half + p_call, ln_half + p_miscall);
                prob11 = prob11 + p_miscall;
            }
            '1' => {
                prob00 = prob00 + p_miscall;
                prob01 = prob01 + LogProb::ln_add_exp(ln_half + p_call, ln_half + p_miscall);
                prob11 = prob11 + p_call;
            }
            _ => {
                panic!("Unexpected allele observed in pileup.");
            }
        }
    }

    // sum of data probabilities
    let posts: Vec<LogProb> = vec![prob00, prob01, prob11];
    let total = LogProb::ln_sum_exp(&posts);

    // genotype posterior probabilities
    let post00 = prob00 - total;
    let post01 = prob01 - total;
    let post11 = prob11 - total;

    (post00, post01, post11)
}

pub fn call_realigned_genotypes_no_haplotypes(flist: &Vec<Fragment>, varlist: &mut VarList) {
    let pileup_lst = generate_simple_pileup(&flist, varlist.lst.len());

    assert_eq!(pileup_lst.len(), varlist.lst.len());

    let genotype_priors = estimate_genotype_priors();

    for i in 0..varlist.lst.len() {
        let pileup = &pileup_lst[i];
        let var = &mut varlist.lst[i];


        let (post00, post01, post11): (LogProb, LogProb, LogProb) = calculate_genotypes_without_haplotypes(&pileup,
                                                                                                           &genotype_priors,
                                                                                                           &var.ref_allele,
                                                                                                           &var.var_allele);

        let genotype: String;
        let genotype_qual: f64;

        if post00 > post01 && post00 > post11 {
            let p_call_wrong = LogProb::ln_add_exp(post01, post11);
            genotype_qual = *PHREDProb::from(p_call_wrong);
            genotype = "0/0".to_string();
        } else if post01 > post00 && post01 > post11 {
            let p_call_wrong = LogProb::ln_add_exp(post00, post11);
            genotype_qual = *PHREDProb::from(p_call_wrong);
            genotype = "0/1".to_string();
        } else {
            let p_call_wrong = LogProb::ln_add_exp(post00, post01);
            genotype_qual = *PHREDProb::from(p_call_wrong);
            genotype = "1/1".to_string();
        }

        //let (count_ref, count_var): (usize, usize) = count_alleles(&pileup);

        var.qual = *PHREDProb::from(post00);
        var.ra = 0;
        var.aa = 0;
        var.genotype = genotype.clone();
        var.gq = genotype_qual;
    }
}

pub fn call_genotypes(flist: &Vec<Fragment>,
                      varlist: &mut VarList,
                      interval: &Option<GenomicInterval>,
                      variant_debug_directory: &Option<String>) {

    let genotype_priors = estimate_genotype_priors();
    let n_var = varlist.lst.len();
    let pileup_lst = generate_fragcall_pileup(&flist, varlist.lst.len());
    assert_eq!(pileup_lst.len(), varlist.lst.len());

    let max_iterations: usize = 1000000;
    //let sample_every: usize = 1;
    let ln_half = LogProb::from(Prob(0.5));
    let mut rng: StdRng = StdRng::from_seed(&[1]);
    let print_time: fn() -> String = || Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

    let hap_ixs = vec![0, 1];

    // obtain var_frags. var_frags[i] contains a vector with the indices of fragments overlapping
    // the i-th variant
    let mut var_frags: Vec<Vec<usize>> = vec![vec![]; varlist.lst.len()];

    for i in 0..flist.len() {
        for call in &flist[i].calls {
            var_frags[call.var_ix].push(i);
        }
    }

    let ln_max_p_miscall = LogProb::from(Prob(MAX_P_MISCALL_F64));
    let mut haps: Vec<Vec<char>> = vec![vec!['0'; n_var]; 2];
    //let mut prev_total_likelihood = LogProb::ln_zero();
    let mut prev_num_phased = 0;

    for hapcut2_iter in 0..max_iterations {

        //let varlist_bak = (*varlist).clone();

        // print the haplotype assembly iteration
        eprintln!("{}    Round {} of haplotype assembly...",print_time(), hapcut2_iter+1);

        let mut var_phased: Vec<bool> = vec![false; varlist.lst.len()];
        let mut vcf_buffer: Vec<Vec<u8>> = Vec::with_capacity(varlist.lst.len());
        let mut hap1: Vec<u8> = vec![0u8; varlist.lst.len()];

        for i in 0..varlist.lst.len() {

            let var = &varlist.lst[i];
            if (var.genotype == "0/1" || var.genotype == "0|1" || var.genotype == "1|0")
                && var.ref_allele.len() == 1 && var.var_allele.len() == 1
                && var.gq > MIN_GQ_FOR_PHASING {
                var_phased[i] = true;

                if var.genotype == "0|1" {
                    hap1[i] = '0' as u8;
                } else if var.genotype == "1|0" {
                    hap1[i] = '1' as u8;
                } else {
                    if rng.next_f64() < 0.5 {
                        hap1[i] = '0' as u8;
                    } else {
                        hap1[i] = '1' as u8;
                    }
                }
            } else {
                hap1[i] = '-' as u8;
            }

            let line: String = format!("{}\t{}\t.\t{}\t{}\t.\t.\tRA={};AA={}\tGT:GQ\t{}:{}",
                                       var.chrom,
                                       var.pos0 + 1,
                                       var.ref_allele,
                                       var.var_allele,
                                       0,
                                       0,
                                       var.genotype,
                                       var.gq);

            let mut vcf_line: Vec<u8> = vec![];
            for u in line.into_bytes() {
                vcf_line.push(u as u8);
            }
            vcf_line.push('\n' as u8);
            vcf_line.push('\0' as u8);
            vcf_buffer.push(vcf_line);
        }

        let frag_buffer = generate_flist_buffer(&flist, &var_phased);
        let mut phase_sets: Vec<i32> = vec![-1i32; varlist.lst.len()];

        call_hapcut2(&frag_buffer,
                     &vcf_buffer,
                     frag_buffer.len(),
                     vcf_buffer.len(),
                     &mut hap1,
                     &mut phase_sets);

        let m = <usize>::max_value();
        let mut min_pos_ps: Vec<usize> = vec![m; varlist.lst.len()];

        for (i,p) in phase_sets.iter().enumerate() {
            if p < &0 {
                continue;
            }
            if varlist.lst[i].pos0 < min_pos_ps[*p as usize] {
                min_pos_ps[*p as usize] = varlist.lst[i].pos0 + 1;
            }
        }

        for i in 0..hap1.len() {
            match hap1[i] as char {
                '0' => {
                    varlist.lst[i].genotype = "0|1".to_string();
                    if phase_sets[i] >= 0 {
                        varlist.lst[i].phase_set = Some(min_pos_ps[phase_sets[i] as usize]);
                    } else {
                        varlist.lst[i].phase_set = None;
                    }
                }
                '1' => {
                    varlist.lst[i].genotype = "1|0".to_string();
                    if phase_sets[i] >= 0 {
                        varlist.lst[i].phase_set = Some(min_pos_ps[phase_sets[i] as usize]);
                    } else {
                        varlist.lst[i].phase_set = None;
                    }
                }
                _ => {}
            }
        }

        for v in 0..varlist.lst.len() {
            let var = &mut varlist.lst[v];


            if var.genotype == "0|0".to_string() || var.genotype == "0/0".to_string() {
                haps[0][v] = '0';
                haps[1][v] = '0';
            } else if var.genotype == "0|1".to_string() {
                haps[0][v] = '0';
                haps[1][v] = '1';
            } else if var.genotype == "1|0".to_string() {
                haps[0][v] = '1';
                haps[1][v] = '0';
            } else if var.genotype == "1|1".to_string() || var.genotype == "1/1".to_string() {
                haps[0][v] = '1';
                haps[1][v] = '1';
            } else if var.genotype == "0/1".to_string() || var.genotype == "1/0".to_string() {
                if rng.next_f64() < 0.5 {
                    haps[0][v] = '0';
                    haps[1][v] = '1';
                } else {
                    haps[0][v] = '1';
                    haps[1][v] = '0';
                }
            } else {
                panic!("Invalid genotype string \"{}\" in Gibbs Sampler genotyping.", var.genotype);
            }
        }

        // p_read_hap[i][j] contains P(R_j | H_i)
        let mut p_read_hap: Vec<Vec<LogProb>> = vec![vec![LogProb::ln_one(); flist.len()]; 2];

        for hap_ix in &hap_ixs {
            for f in 0..flist.len() {
                for call in &flist[f].calls {
                    if var_phased[call.var_ix] && call.qual < ln_max_p_miscall {
                        // read allele matches haplotype allele
                        if call.allele == haps[*hap_ix][call.var_ix] {
                            p_read_hap[*hap_ix][f] = p_read_hap[*hap_ix][f] + &call.one_minus_qual;
                        } else { // read allele does not match haplotype allele
                            p_read_hap[*hap_ix][f] = p_read_hap[*hap_ix][f] + call.qual;
                        }
                    }
                }
            }
        }

        eprintln!("{}    Refining genotypes using round {} haplotypes...",print_time(), hapcut2_iter+1);

        for _ in 0..max_iterations {
            let mut changed = false;

            // take the variants v in random order
            let mut ixvec: Vec<usize> = (0..varlist.lst.len()).collect();
            let ixslice: &mut [usize] = ixvec.as_mut_slice();
            rng.shuffle(ixslice);

            for v_r in ixslice {
                let v = *v_r;

                let var = &mut varlist.lst[v];

                /*
                match interval {
                    &Some(ref iv) => {
                        if var.chrom != iv.chrom ||
                            var.pos0 < iv.start_pos as usize ||
                            var.pos0 > iv.end_pos as usize {
                            continue;
                        }
                    }
                    &None => {}
                }
                */
                let mut p_reads: [LogProb; 4] = [LogProb::ln_one(); 4];

                // let g be the current genotype being considered to switch to
                // then p_read_g[g] contains a vector of tuples (frag_ix, p_read_h0, p_read_h0)
                // that have the probability of each fragment under the new genotypes
                let mut p_read_g: Vec<Vec<(usize, LogProb, LogProb)>> = vec![vec![]; 4];

                for g in 0..4 {
                    for call in &pileup_lst[v] {
                        if call.qual >= ln_max_p_miscall {
                            continue;
                        }

                        // get the value of the read likelihood given each haplotype
                        let (mut p_read_h0, mut p_read_h1) = match call.frag_ix {
                            Some(frag_ix) => (p_read_hap[0][frag_ix], p_read_hap[1][frag_ix]),
                            None => panic!("ERROR: Fragment index is missing in pileup iteration.")
                        };

                        // for each haplotype allele
                        // if that allele on the haplotype changes in g,
                        // then we divide out the old value and multiply in the new value

                        // haplotype 0 at site j was a '1', but now it will be '0'
                        // if call on read was '1', we divide out (1 - q) and multiply in q
                        // if call on read was '0', we divide out q and multiply in (1 - q)

                        if (g == 0 || g == 1) && (haps[0][v] == '1') {
                            if call.allele == '1' {
                                p_read_h0 = p_read_h0 - call.one_minus_qual + call.qual;
                            } else if call.allele == '0' {
                                p_read_h0 = p_read_h0 - call.qual + call.one_minus_qual;
                            } else {
                                panic!("Unsupported allele");
                            }
                        }

                        // haplotype 0 at site j was a '0', but now it will be '1'
                        // if call on read was '0', we divide out (1 - q) and multiply in q
                        // if call on read was '1', we divide out q and multiply in (1 - q)
                        if (g == 2 || g == 3) && (haps[0][v] == '0') {
                            if call.allele == '0' {
                                p_read_h0 = p_read_h0 - call.one_minus_qual + call.qual;
                            } else if call.allele == '1' {
                                p_read_h0 = p_read_h0 - call.qual + call.one_minus_qual;
                            } else {
                                panic!("Unsupported allele");
                            }
                        }

                        // haplotype 1 at site j was a '1', but now it will be '0'
                        // if call on read was '1', we divide out (1 - q) and multiply in q
                        // if call on read was '0', we divide out q and multiply in (1 - q)

                        if (g == 0 || g == 2) && (haps[1][v] == '1') {
                            if call.allele == '1' {
                                p_read_h1 = p_read_h1 - call.one_minus_qual + call.qual;
                            } else if call.allele == '0' {
                                p_read_h1 = p_read_h1 - call.qual + call.one_minus_qual;
                            } else {
                                panic!("Unsupported allele");
                            }
                        }

                        // haplotype 1 at site j was a '0', but now it will be '1'
                        // if call on read was '0', we divide out (1 - q) and multiply in q
                        // if call on read was '1', we divide out q and multiply in (1 - q)
                        if (g == 1 || g == 3) && (haps[1][v] == '0') {
                            if call.allele == '0' {
                                p_read_h1 = p_read_h1 - call.one_minus_qual + call.qual;
                            } else if call.allele == '1' {
                                p_read_h1 = p_read_h1 - call.qual + call.one_minus_qual;
                            } else {
                                panic!("Unsupported allele");
                            }
                        }

                        match call.frag_ix {
                            Some(frag_ix) => { p_read_g[g].push((frag_ix, p_read_h0, p_read_h1)); },
                            None => { panic!("Fragment index in pileup call is None.") },
                        }

                        let p_read = LogProb::ln_add_exp(ln_half + p_read_h0, ln_half + p_read_h1);
                        p_reads[g] = p_reads[g] + p_read;
                    }

                    let mut prior;
                    if var.ref_allele.len() == 1 && var.var_allele.len() == 1 {
                        let ra = var.ref_allele.chars().nth(0).unwrap();
                        let g1 = if g == 0 || g == 1 { var.ref_allele.chars().nth(0).unwrap() } else { var.var_allele.chars().nth(0).unwrap() };
                        let g2 = if g == 0 || g == 2 { var.ref_allele.chars().nth(0).unwrap() } else { var.var_allele.chars().nth(0).unwrap() };

                        match genotype_priors.get(&(ra, (g1, g2))) {
                            Some(p) => {
                                prior = if g == 1 || g == 2 { ln_half + *p } else { *p };
                            },
                            None => {
                                match genotype_priors.get(&(ra, (g2, g1))) {
                                    Some(p) => {
                                        prior = if g == 1 || g == 2 { ln_half + *p } else { *p };
                                    },
                                    None => { panic!("Genotype not in genotype priors."); }
                                };
                            }
                        }
                    } else {
                        // we just assume that the rate of indels is 1e-4
                        prior = if g == 0 {
                            LogProb::from(Prob(0.9999))
                        } else {
                            LogProb::from(Prob(0.0001 / 3.0))
                        };
                    }

                    //println!("{}",*prior);
                    p_reads[g] = p_reads[g] + prior;
                }

                let p_total = LogProb::ln_sum_exp(&p_reads);

                var.genotype_post[0] = p_reads[0] - p_total;
                var.genotype_post[1] = p_reads[1] - p_total;
                var.genotype_post[2] = p_reads[2] - p_total;
                var.genotype_post[3] = p_reads[3] - p_total;

                let mut max_ix = 5;
                let mut max_val = LogProb::from(Prob(0.0));

                for k in 0..4 {
                    if var.genotype_post[k] > max_val {
                        max_ix = k;
                        max_val = var.genotype_post[k];
                    }
                }

                match max_ix {
                    0 => {
                        if !(haps[0][v] == '0' && haps[1][v] == '0') { changed = true };
                        haps[0][v] = '0';
                        haps[1][v] = '0';
                    },
                    1 => {
                        if !(haps[0][v] == '0' && haps[1][v] == '1') { changed = true };
                        haps[0][v] = '0';
                        haps[1][v] = '1';
                    },
                    2 => {
                        if !(haps[0][v] == '1' && haps[1][v] == '0') { changed = true };
                        haps[0][v] = '1';
                        haps[1][v] = '0';
                    },
                    3 => {
                        if !(haps[0][v] == '1' && haps[1][v] == '1') { changed = true };
                        haps[0][v] = '1';
                        haps[1][v] = '1';
                    },
                    _ => { panic!("Invalid genotype  index.") }
                }

                for &(frag_ix, p_read_h0, p_read_h1) in &p_read_g[max_ix] {
                    p_read_hap[0][frag_ix] = p_read_h0;
                    p_read_hap[1][frag_ix] = p_read_h1;
                }

                let new_qual: f64 = *PHREDProb::from(max_val);

                // the quality score of this variant has increased above the limit
                // we add it from the pool of phased variants
                if new_qual >= MIN_GQ_FOR_PHASING && var_phased[v] == false {

                    // need to visit each fragment overlapping vth variant and multiply (add) in the
                    // amount the call at this site contributes to fragment likelihoods
                    for call in &pileup_lst[v] {
                        if call.qual >= ln_max_p_miscall {
                            continue;
                        }

                        let frag_ix = match call.frag_ix {
                            Some(f_ix) => f_ix,
                            None => panic!("ERROR: Fragment index is missing in pileup iteration.")
                        };

                        match call.allele {
                            '0' => {
                                match haps[0][v] {
                                    '0' => { p_read_hap[0][frag_ix] = p_read_hap[0][frag_ix] + call.one_minus_qual; },
                                    '1' => { p_read_hap[0][frag_ix] = p_read_hap[0][frag_ix] + call.qual; },
                                    _ => { panic!("Invalid allele in haplotype."); }
                                }
                                match haps[1][v] {
                                    '0' => { p_read_hap[1][frag_ix] = p_read_hap[1][frag_ix] + call.one_minus_qual; },
                                    '1' => { p_read_hap[1][frag_ix] = p_read_hap[1][frag_ix] + call.qual; },
                                    _ => { panic!("Invalid allele in haplotype."); }
                                }
                            },
                            '1' => {
                                match haps[0][v] {
                                    '0' => { p_read_hap[0][frag_ix] = p_read_hap[0][frag_ix] + call.qual; },
                                    '1' => { p_read_hap[0][frag_ix] = p_read_hap[0][frag_ix] + call.one_minus_qual; },
                                    _ => { panic!("Invalid allele in haplotype."); }
                                }
                                match haps[1][v] {
                                    '0' => { p_read_hap[1][frag_ix] = p_read_hap[1][frag_ix] + call.qual; },
                                    '1' => { p_read_hap[1][frag_ix] = p_read_hap[1][frag_ix] + call.one_minus_qual; },
                                    _ => { panic!("Invalid allele in haplotype."); }
                                }
                            },
                            _ => { panic!("Invalid allele in fragment.") }
                        }
                    }

                    var_phased[v] = true;
                }

                // the quality score of this variant has dipped below the limit
                // we remove it from the pool of phased variants
                if new_qual < MIN_GQ_FOR_PHASING && var_phased[v] == true {

                    // need to visit each fragment overlapping vth variant and divide (subtract) out the
                    // amount the call at this site contributes to fragment likelihoods

                    for call in &pileup_lst[v] {
                        if call.qual >= ln_max_p_miscall {
                            continue;
                        }

                        let frag_ix = match call.frag_ix {
                            Some(f_ix) => f_ix,
                            None => panic!("ERROR: Fragment index is missing in pileup iteration.")
                        };

                        match call.allele {
                            '0' => {
                                match haps[0][v] {
                                    '0' => { p_read_hap[0][frag_ix] = p_read_hap[0][frag_ix] - call.one_minus_qual; },
                                    '1' => { p_read_hap[0][frag_ix] = p_read_hap[0][frag_ix] - call.qual; },
                                    _ => { panic!("Invalid allele in haplotype."); }
                                }
                                match haps[1][v] {
                                    '0' => { p_read_hap[1][frag_ix] = p_read_hap[1][frag_ix] - call.one_minus_qual; },
                                    '1' => { p_read_hap[1][frag_ix] = p_read_hap[1][frag_ix] - call.qual; },
                                    _ => { panic!("Invalid allele in haplotype."); }
                                }
                            },
                            '1' => {
                                match haps[0][v] {
                                    '0' => { p_read_hap[0][frag_ix] = p_read_hap[0][frag_ix] - call.qual; },
                                    '1' => { p_read_hap[0][frag_ix] = p_read_hap[0][frag_ix] - call.one_minus_qual; },
                                    _ => { panic!("Invalid allele in haplotype."); }
                                }
                                match haps[1][v] {
                                    '0' => { p_read_hap[1][frag_ix] = p_read_hap[1][frag_ix] - call.qual; },
                                    '1' => { p_read_hap[1][frag_ix] = p_read_hap[1][frag_ix] - call.one_minus_qual; },
                                    _ => { panic!("Invalid allele in haplotype."); }
                                }
                            },
                            _ => { panic!("Invalid allele in fragment.") }
                        }
                    }

                    var_phased[v] = false;
                }
            }


            // if the haplotypes have not changed in this iteration, then we break
            if !changed {
                break;
            }
        }

        let mut total_likelihood: LogProb = LogProb::ln_one();
        let mut num_phased = 0;

        for var in varlist.lst.iter() {
            if (var.genotype == "0/1" || var.genotype == "0|1" || var.genotype == "1|0")
                && var.ref_allele.len() == 1 && var.var_allele.len() == 1
                && var.gq > MIN_GQ_FOR_PHASING {
                num_phased += 1;
            }
        }

        for f in 0..flist.len() {
            let mut pr: Vec<LogProb> = vec![LogProb::ln_one(); 2];
            for hap_ix in &hap_ixs {
                for call in &flist[f].calls {
                    if call.qual < ln_max_p_miscall {
                        // read allele matches haplotype allele
                        if call.allele == haps[*hap_ix][call.var_ix] {
                            pr[*hap_ix] = pr[*hap_ix] + &call.one_minus_qual;
                        } else { // read allele does not match haplotype allele
                            pr[*hap_ix] = pr[*hap_ix] + call.qual;
                        }
                    }
                }
            }
            total_likelihood = total_likelihood + LogProb::ln_add_exp(ln_half + pr[0], ln_half + pr[1]);
        }


        // update the various fields for the variant.

        for i in 0..varlist.lst.len() {
            let pileup = &pileup_lst[i];
            let var = &mut varlist.lst[i];

            //let total: f64 = var.genotype_counts[0] +  var.genotype_counts[1] +
            //                var.genotype_counts[2] +  var.genotype_counts[3];

            let post00: LogProb = var.genotype_post[0];
            let post01: LogProb = var.genotype_post[1];
            let post10: LogProb = var.genotype_post[2];
            let post11: LogProb = var.genotype_post[3];

            let genotype: String;
            let mut genotype_qual: f64;


            let post00_unphased = post00;
            let post01_unphased = LogProb::ln_add_exp(post01, post10);
            let post11_unphased = post11;

            if post00_unphased > post01_unphased && post00_unphased > post11_unphased {
                let p_call_wrong = LogProb::ln_add_exp(post01_unphased, post11_unphased);
                genotype_qual = *PHREDProb::from(p_call_wrong);
                genotype = match var.phase_set {
                    Some(_) => {"0|0".to_string()}
                    None => {"0/0".to_string()}
                };
            } else if post01_unphased > post00_unphased && post01_unphased > post11_unphased {
                let p_call_wrong = LogProb::ln_add_exp(post00_unphased, post11_unphased);
                genotype_qual = *PHREDProb::from(p_call_wrong);

                // we only show phase if variant was phased in the first round
                // we do consider "flipping", if posteriors shifted when genotyping in round 2.

                genotype = match var.phase_set {
                    Some(_) => {
                        if post01 > post10 {
                            "0|1".to_string()
                        } else {
                            "1|0".to_string()
                        }
                    }
                    None => {"0/1".to_string()}
                };

            } else {
                let p_call_wrong = LogProb::ln_add_exp(post00_unphased, post01_unphased);
                genotype_qual = *PHREDProb::from(p_call_wrong);
                genotype = match var.phase_set {
                    Some(_) => {"1|1".to_string()}
                    None => {"1/1".to_string()}
                };
            }

            let (count_ref, count_var, count_amb): (usize, usize, usize) = count_alleles(&pileup);

            var.qual = *PHREDProb::from(post00);
            var.ra = count_ref;
            var.aa = count_var;
            var.na = count_amb;
            var.genotype = genotype;
            var.gq = genotype_qual;
            var.filter = "PASS".to_string();

        }

        let debug_vcf_str = format!("3.{}.haplotype_genotype_iteration.vcf", hapcut2_iter).to_owned();
        print_variant_debug(&varlist, &interval, &variant_debug_directory,&debug_vcf_str);

        eprintln!("{}    Total phased heterozygous SNVs: {}  Total likelihood (phred): {:.2}",print_time(), num_phased, *PHREDProb::from(total_likelihood));

        if num_phased <= prev_num_phased { //if total_likelihood <= prev_total_likelihood {

            // restore the previous varlist
            //for i in 0..varlist.lst.len() {
            //    varlist.lst[i] = varlist_bak.lst[i].clone();
            //}

            break;
        }

        prev_num_phased = num_phased;
        //prev_total_likelihood = total_likelihood;

    }
}

pub fn calculate_mec(flist: &Vec<Fragment>, varlist: &mut VarList) {

    let hap_ixs = vec![0, 1];
    let ln_max_p_miscall = LogProb::from(Prob(MAX_P_MISCALL_F64));

    for mut var in &mut varlist.lst {
        var.mec = 0;
        var.mec_frac = 0.0;
    }

    for f in 0..flist.len() {

        let mut mismatched_vars: Vec<Vec<usize>> = vec![vec![], vec![]];

        for &hap_ix in &hap_ixs {
            for call in &flist[f].calls {
                if call.qual < ln_max_p_miscall {
                    let mut chs = varlist.lst[call.var_ix].genotype.chars();
                    let g1 = chs.next();
                    let sep = chs.next();
                    let g2 = chs.next();
                    let g = vec![g1.unwrap(),g2.unwrap()];
                    if sep.unwrap() != '|' {
                        continue; // only care about phased variants.
                    }
                    // read allele matches haplotype allele
                    if call.allele != g[hap_ix] {
                        mismatched_vars[hap_ix].push(call.var_ix);
                    }
                }
            }
        }

        let min_error_hap = if mismatched_vars[0].len() < mismatched_vars[1].len() { 0 } else { 1 };

        for &ix in &mismatched_vars[min_error_hap] {
            varlist.lst[ix].mec += 1;
        }
    }

    let mut block_mec: HashMap<usize, usize> = HashMap::new();
    let mut block_total: HashMap<usize, usize> = HashMap::new();

    for mut var in &mut varlist.lst {
        match var.phase_set {
            Some(ps) => {
                *block_mec.entry(ps).or_insert(0) += var.mec;
                *block_total.entry(ps).or_insert(0) += var.ra + var.aa;
            }
            None => {}
        }
    }

    for mut var in &mut varlist.lst {
        match var.phase_set {
            Some(ps) => {
                var.mec_frac = *block_mec.get(&ps).unwrap() as f64 / *block_total.get(&ps).unwrap() as f64;
            }
            None => {}
        }
    }
}

pub fn var_filter(varlist: &mut VarList, density_qual: f64, density_dist: usize, density_count: usize, max_depth: Option<u32>, max_mec_frac: Option<f64>) {

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
    match max_mec_frac {
        Some(frac) => {
            for i in 0..varlist.lst.len() {
                if varlist.lst[i].mec_frac >= frac {
                    if varlist.lst[i].filter == ".".to_string() || varlist.lst[i].filter == "PASS".to_string() {
                        varlist.lst[i].filter = "psmf".to_string();
                    } else {
                        varlist.lst[i].filter.push_str(";psmf");
                    }
                }
            }
        },
        None => {}
    }
}

pub fn print_vcf(varlist: &VarList, interval: &Option<GenomicInterval>, indels: bool, output_vcf_file: &String, print_whole_varlist: bool) {
    let vcf_path = Path::new(output_vcf_file);
    let vcf_display = vcf_path.display();
    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match File::create(&vcf_path) {
        Err(why) => panic!("couldn't create {}: {}", vcf_display, why.description()),
        Ok(file) => file,
    };

    let headerstr = "##fileformat=VCFv4.2
##source=ReaperV0.1
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=RA,Number=1,Type=Integer,Description=\"Number of Observed Reference Alleles\">
##INFO=<ID=AA,Number=1,Type=Integer,Description=\"Number of Observed Variant Alleles\">
##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Number of Observed Ambiguous Alleles\">
##INFO=<ID=P00,Number=1,Type=Integer,Description=\"Phred-scaled Probability of 0|0 Phased Genotype\">
##INFO=<ID=P01,Number=1,Type=Integer,Description=\"Phred-scaled Probability of 0|1 Phased Genotype\">
##INFO=<ID=P10,Number=1,Type=Integer,Description=\"Phred-scaled Probability of 1|0 Phased Genotype\">
##INFO=<ID=P11,Number=1,Type=Integer,Description=\"Phred-scaled Probability of 1|1 Phased Genotype\">
##INFO=<ID=MEC,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Score for Variant\">
##INFO=<ID=MF,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Fraction for Variant\">
##INFO=<ID=PSMF,Number=1,Type=Integer,Description=\"Minimum Error Criterion (MEC) Fraction for Phase Set\">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">
##FORMAT=<ID=GQ,Number=2,Type=Float,Description=\"Genotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        .to_string();

    match writeln!(file, "{}", headerstr) {
        Err(why) => panic!("couldn't write to {}: {}", vcf_display, why.description()),
        Ok(_) => {}
    }


    for var in &varlist.lst {
        match interval {
            &Some(ref iv) => {
                if var.chrom != iv.chrom ||
                    var.pos0 < iv.start_pos as usize ||
                    var.pos0 > iv.end_pos as usize {
                    continue;
                }
            }
            &None => {}
        }

        if !print_whole_varlist {
            if var.genotype == "0/0".to_string() ||
                var.genotype == "0|0".to_string(){
                continue;
            }

            if !indels && (var.ref_allele.len() != 1 || var.var_allele.len() != 1) {
                continue;
            }
        }
        let ps = match var.phase_set {
            Some(ps) => format!("{}", ps),
            None => ".".to_string()
        };

        match writeln!(file,
                       "{}\t{}\t.\t{}\t{}\t{:.2}\t{}\tDP={};RA={};AA={};NA={};P00={:.2};P01={:.2};P10={:.2};P11={:.2};MEC={};MF={};PSMF={:.5}\tGT:PS:GQ\t{}:{}:{:.2}",
                       var.chrom,
                       var.pos0 + 1,
                       var.ref_allele,
                       var.var_allele,
                       var.qual,
                       var.filter,
                       var.dp,
                       var.ra,
                       var.aa,
                       var.na,
                       *PHREDProb::from(var.genotype_post[0]),
                       *PHREDProb::from(var.genotype_post[1]),
                       *PHREDProb::from(var.genotype_post[2]),
                       *PHREDProb::from(var.genotype_post[3]),
                       var.mec,
                       (var.mec as f64 / (var.ra + var.aa) as f64),
                       var.mec_frac,
                       var.genotype,
                       ps,
                       var.gq) {
            Err(why) => panic!("couldn't write to {}: {}", vcf_display, why.description()),
            Ok(_) => {}
        }
    }
}

pub fn print_variant_debug(varlist: &VarList, interval: &Option<GenomicInterval>, variant_debug_directory: &Option<String>, debug_filename: &str){
    match variant_debug_directory {
        &Some(ref dir) => {
            let outfile = match Path::new(&dir).join(&debug_filename).to_str() {
            Some(s) => {s.to_owned()},
            None => {panic!("Invalid unicode provided for variant debug directory");}
            };
            print_vcf(&varlist, &interval, true, &outfile, true);
        }
        &None => {}
    };
}
