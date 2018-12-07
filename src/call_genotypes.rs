//! This module contains functions for calling diploid genotypes.
//!
//! It has functions for calling genotypes without haplotype information (basic pileup-based
//! calculation similar to samtools), as well as a function for refining genotypes by iteratively
//! assembling haplotypes.

// use declarations
use bio::stats::{PHREDProb, LogProb, Prob};
use util::{MAX_VCF_QUAL, GenomicInterval, DensityParameters};
use variants_and_fragments::*;
use print_output::*;
use genotype_probs::*;
use haplotype_assembly::{generate_flist_buffer, call_hapcut2};
use chrono::prelude::*;
use errors::*;
use rand::{Rng,StdRng,SeedableRng};

/// Takes a vector of fragments and returns a vector of "allele pileups"
///
/// The allele pileup of a variant site is the list of all allele observations (FragCalls) that occur at that site.
///
/// # Arguments
/// - flist: a vector of Fragments
/// - n_var: the total number of variants (e.g. the length of VarList)
///
/// # Returns
/// Returns the allele pileups for all variant sites as a vector of vectors.
/// The outer vector is indexed by variant index (same index as VarList) and contains the allele
/// pileup for that variant, represented as a vector of of FragCalls.
///
/// # Example
/// see ```call_genotypes::tests::test_generate_fragcall_pileup()```
fn generate_fragcall_pileup(flist: &Vec<Fragment>, n_var: usize) -> Vec<Vec<FragCall>> {
    let mut pileup_lst: Vec<Vec<FragCall>> = vec![vec![]; n_var];
    for fragment in flist {
        for call in fragment.clone().calls {
            pileup_lst[call.var_ix].push(call);
        }
    }
    pileup_lst
}

/// Counts the number of each allele in an allele pileup
///
/// # Arguments
/// - pileup: an allele pileup for some variant site (represented as a vector of ```FragCalls```)
/// - num_alleles: how many alleles does this variant site have (2 for biallelic, 3 for triallelic...)
/// - max_p_miscall: the maximum probability of an allele miscall to count the allele (equivalent
///                  to the minimum allowed allele quality, but represented as a normal probability
///                  rather than PHRED-scaled)
///
/// # Returns
/// Returns a tuple containing ```(counts, count_amb)``` where ```counts``` has length ```num_alleles```
/// and contains the count for each allele (```counts[0]``` is the reference allele, etc).
/// ```counts_amb``` is the number of ambiquous alleles that fell beneath the quality cutoff.
fn count_alleles(pileup: &Vec<FragCall>, num_alleles: usize, max_p_miscall: f64) -> (Vec<usize>, usize) {
    // compute probabilities of data given each genotype
    let mut counts: Vec<usize> = vec![0; num_alleles]; // counts for each allele
    let mut count_amb = 0; // ambiguous allele count
    let ln_max_p_miscall = LogProb::from(Prob(max_p_miscall));

    for call in pileup {
        if call.qual < ln_max_p_miscall {
            counts[call.allele as usize] += 1;
        } else {
            count_amb += 1;
        }
    }

    (counts, count_amb)
}

/// Calculates the posterior probabilities for a pileup-based genotyping calculation (without using
/// haplotype information)
///
/// # Arguments
/// - pileup: an allele pileup for some variant site (represented as a vector of ```FragCalls```)
/// - genotype_priors: a struct holding the genotype prior probabilities
/// - alleles: a vector of the alleles (as Strings) for this variant site. ```alleles[0]``` should
///            be the ref allele, and alleles must be in same order as the allele indices held in the
///            pileup ```FragCalls```.
///
/// # Returns
/// Returns a Result holding a ```GenotypeProbs``` struct.
/// These ```GenotypeProbs``` hold the posterior genotype probabilities for each genotype
/// calculated from the allele pileup
///
/// # Errors
/// - Can throw an error if attempts to query ```genotype_priors``` using an invalid genotype
pub fn calculate_genotype_posteriors_no_haplotypes(pileup: &Vec<FragCall>,
                                                   genotype_priors: &GenotypePriors,
                                                   alleles: &Vec<String>,
                                                   max_p_miscall: f64) -> Result<GenotypeProbs> {

    let ln_max_p_miscall: LogProb = LogProb::from(Prob(max_p_miscall));
    let ln_half: LogProb = LogProb::from(Prob(0.5));

    // this matrix holds P(data | g) * p(g)
    let mut probs: GenotypeProbs = genotype_priors.get_all_priors(alleles).chain_err(|| "Error getting all genotype priors while calculating genotypes.")?;

    for &call in pileup {

        let allele = call.allele;
        let p_miscall = call.qual;
        let p_call = call.one_minus_qual;

        if p_miscall >= ln_max_p_miscall {
            continue;
        }

        // for each possible genotype, update the P(data | g) * p(g)
        for g in possible_genotypes(alleles) {
            if g.0 == allele && g.1 == allele {
                probs.ln_times_equals(g, p_call);
            } else if g.0 != allele && g.1 != allele {
                probs.ln_times_equals(g, p_miscall);
            } else {
                let p_het = LogProb::ln_add_exp(ln_half + p_call, ln_half + p_miscall);
                probs.ln_times_equals(g, p_het);
            }
        }
    }

    // posterior probabilities
    let posts = probs.normalize();

    Ok(posts)
}

/// Calls genotypes for each variant in the ```VarList``` using a pileup-based genotyping calculation
/// (without using haplotype information)
///
/// #Arguments
/// -flist: a vector of Fragment structs
/// - 
pub fn call_genotypes_no_haplotypes(flist: &Vec<Fragment>, varlist: &mut VarList, genotype_priors: &GenotypePriors, max_p_miscall: f64) -> Result<()> {
    let pileup_lst = generate_fragcall_pileup(&flist, varlist.lst.len());

    assert_eq!(pileup_lst.len(), varlist.lst.len());

    for i in 0..varlist.lst.len() {
        let pileup = &pileup_lst[i];
        let var = &mut varlist.lst[i];

        let posts: GenotypeProbs = calculate_genotype_posteriors_no_haplotypes(&pileup,
                                                                               &genotype_priors,
                                                                               &var.alleles,
                                                                               max_p_miscall).chain_err(|| "Error calculating genotype posteriors for haplotype-free genotyping")?;

        let (max_g, max_post) = posts.max_genotype(false, false);

        let genotype_qual:f64 = *PHREDProb::from(LogProb::ln_one_minus_exp(&max_post));

        let (allele_counts, ambig_count) = count_alleles(&pileup, var.alleles.len(), max_p_miscall);
        let allele_total: usize = allele_counts.iter().sum::<usize>() + ambig_count;

        if var.dp < allele_total {
            var.dp = allele_total;
        }

        var.qual = *PHREDProb::from(posts.get(Genotype(0,0)));
        if var.qual > MAX_VCF_QUAL {
            var.qual = MAX_VCF_QUAL;
        }

        var.genotype = max_g;
        var.allele_counts = allele_counts;
        var.ambiguous_count = ambig_count;
        var.unphased_genotype = max_g;
        var.gq = genotype_qual;
        var.unphased_gq = genotype_qual;
        if var.gq > MAX_VCF_QUAL {
            var.gq = MAX_VCF_QUAL;
        }

        if var.unphased_gq > MAX_VCF_QUAL {
            var.unphased_gq = MAX_VCF_QUAL;
        }

        var.phase_set = None;
    }
    Ok(())
}

pub fn call_genotypes_with_haplotypes(flist: &mut Vec<Fragment>,
                                      varlist: &mut VarList,
                                      interval: &Option<GenomicInterval>,
                                      genotype_priors: &GenotypePriors,
                                      variant_debug_directory: &Option<String>,
                                      program_step: usize,
                                      max_cov: u32,
                                      density_params: &DensityParameters,
                                      max_p_miscall: f64,
                                      sample_name: &String,
                                      ll_delta: f64) -> Result<()> {

    let min_hap_gq: f64 = 1.0;

    let n_var = varlist.lst.len();
    let pileup_lst = generate_fragcall_pileup(&flist, varlist.lst.len());
    assert_eq!(pileup_lst.len(), varlist.lst.len());

    let max_iterations: usize = 1000000;
    //let sample_every: usize = 1;
    let ln_half = LogProb::from(Prob(0.5));
    let mut rng: StdRng = StdRng::from_seed(&[0]);
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

    let ln_max_p_miscall = LogProb::from(Prob(max_p_miscall));
    let mut haps: Vec<Vec<u8>> = vec![vec![0u8; n_var]; 2];
    let mut prev_likelihood = LogProb::ln_zero();

    // for all basic biallelic heterozygous variants
    // randomly shuffle the phase of the variant
    for i in 0..varlist.lst.len() {
        let var = &mut varlist.lst[i];
        if var.alleles.len() == 2 && (var.genotype == Genotype(0,1) || var.genotype == Genotype(1,0))
            && var.alleles[0].len() == 1 && var.alleles[1].len() == 1 {
            if rng.next_f64() < 0.5 {
                var.genotype = Genotype(0,1);
            } else {
                var.genotype = Genotype(1,0);
            }
        }
    }

    for hapcut2_iter in 0..max_iterations {

        // print the haplotype assembly iteration
        eprintln!("{}    Round {} of haplotype assembly...",print_time(), hapcut2_iter+1);

        ////////////////////////////////////////////////////////////////////////////////////////////

        // count how many variants meet the criteria for "phased"
        let mut num_phased = 0;
        for var in varlist.lst.iter() {
            if var.alleles.len() == 2 && (var.genotype == Genotype(0,1) || var.genotype == Genotype(1,0))
                && var.alleles[0].len() == 1 && var.alleles[1].len() == 1
                && var.gq > min_hap_gq {
                num_phased += 1;
            }
        }

        let mut total_likelihood: LogProb = LogProb::ln_one();

        for v in 0..varlist.lst.len() {
            let g = Genotype(haps[0][v], haps[1][v]);
            total_likelihood = total_likelihood + genotype_priors.get_prior(&varlist.lst[v].alleles, g)?;
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

        eprintln!("{}    (Before HapCUT2) Total phased heterozygous SNVs: {}  Total likelihood (phred): {:.2}",print_time(), num_phased, *PHREDProb::from(total_likelihood));

        ////////////////////////////////////////////////////////////////////////////////////////////

        let mut var_phased: Vec<bool> = vec![false; varlist.lst.len()];
        let mut vcf_buffer: Vec<Vec<u8>> = Vec::with_capacity(varlist.lst.len());
        let mut hap1: Vec<u8> = vec!['-' as u8; varlist.lst.len()];

        for i in 0..varlist.lst.len() {

            let var = &varlist.lst[i];
            if var.alleles.len() == 2 && (var.genotype == Genotype(0,1) || var.genotype == Genotype(1,0))
                && var.alleles[0].len() == 1 && var.alleles[1].len() == 1
                && var.gq > min_hap_gq {
                var_phased[i] = true;

                if var.genotype == Genotype(0,1) {
                    hap1[i] = '0' as u8;
                } else if var.genotype == Genotype(1,0) {
                    hap1[i] = '1' as u8;
                }
            } else {
                hap1[i] = '-' as u8;
            }

            let genotype_str = vec![var.genotype.0.to_string(), var.genotype.1.to_string()].join("/");
            let line: String = format!("{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT:GQ\t{}:{}",
                                       var.chrom,
                                       var.pos0 + 1,
                                       var.alleles[0],
                                       var.alleles[1],
                                       genotype_str,
                                       var.gq);

            //println!("{}", line);

            let mut vcf_line: Vec<u8> = vec![];
            for u in line.into_bytes() {
                vcf_line.push(u as u8);
            }
            vcf_line.push('\n' as u8);
            vcf_line.push('\0' as u8);
            vcf_buffer.push(vcf_line);
        }

        let frag_buffer = generate_flist_buffer(&flist, &var_phased, max_p_miscall, false).chain_err(|| "Error generating fragment list buffer.")?;
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
                    varlist.lst[i].genotype = Genotype(0,1);

                    if phase_sets[i] >= 0 {
                        varlist.lst[i].phase_set = Some(min_pos_ps[phase_sets[i] as usize]);
                    } else {
                        varlist.lst[i].phase_set = None;
                    }
                }
                '1' => {
                    varlist.lst[i].genotype = Genotype(1,0);

                    if phase_sets[i] >= 0 {
                        varlist.lst[i].phase_set = Some(min_pos_ps[phase_sets[i] as usize]);
                    } else {
                        varlist.lst[i].phase_set = None;
                    }
                }
                _ => {
                    var_phased[i] = false;
                }
            }

            haps[0][i] = varlist.lst[i].genotype.0;
            haps[1][i] = varlist.lst[i].genotype.1;

        }


        ////////////////////////////////////////////////////////////////////////////////////////////

        // count how many variants meet the criteria for "phased"
        num_phased = 0;
        for var in varlist.lst.iter() {
            if var.alleles.len() == 2 && (var.genotype == Genotype(0,1) || var.genotype == Genotype(1,0))
                && var.alleles[0].len() == 1 && var.alleles[1].len() == 1
                && var.gq > min_hap_gq {
                num_phased += 1;
            }
        }

        total_likelihood = LogProb::ln_one();

        for v in 0..varlist.lst.len() {
            let g = Genotype(haps[0][v], haps[1][v]);
            total_likelihood = total_likelihood + genotype_priors.get_prior(&varlist.lst[v].alleles, g)?;
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

        eprintln!("{}    (After HapCUT2)  Total phased heterozygous SNVs: {}  Total likelihood (phred): {:.2}",print_time(), num_phased, *PHREDProb::from(total_likelihood));

        ////////////////////////////////////////////////////////////////////////////////////////////

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

        for _ in 0..max_iterations {
            let mut changed = false;

            // take the variants v in random order
            let mut ixvec: Vec<usize> = (0..varlist.lst.len()).collect();
            let ixslice: &mut [usize] = ixvec.as_mut_slice();
            rng.shuffle(ixslice);

            for v_r in ixslice {
                let v = *v_r;

                let var = &mut varlist.lst[v];

                assert_eq!(v, var.ix);

                let mut p_reads: GenotypeProbs = genotype_priors.get_all_priors(&var.alleles).chain_err(|| "Error getting all genotype priors while calculating haplotype-informed genotypes")?;

                // let (g1,g2) be the current genotype being considered to switch to
                // then p_read_g[g1][g2] contains a vector of tuples (frag_ix, p_read_h0, p_read_h1
                // that have the probability of each fragment under the new genotypes
                let mut p_reads_g: Vec<Vec<Vec<(usize, LogProb, LogProb)>>> = vec![vec![vec![]; var.alleles.len()]; var.alleles.len()];


                // THOUGHT -- p_read_hap contains the probability of reads given haplotypes,
                // ONLY for variants in variant_phased.
                // so to correctly calculate p_reads for a variant, if the variant isn't in var_phased
                // we should simply multiply in the values
                // if it is in var phased we do the whole divide out, multiply in stuff
                for g in var.possible_genotypes() {

                    for call in &pileup_lst[v] {
                        if call.qual >= ln_max_p_miscall {
                            continue;
                        }

                        // get the value of the read likelihood given each haplotype
                        let (mut p_read_h0, mut p_read_h1) = match call.frag_ix {
                            Some(frag_ix) => (p_read_hap[0][frag_ix], p_read_hap[1][frag_ix]),
                            None => panic!("ERROR: Fragment index is missing in pileup iteration.")
                        };

                        if var_phased[v] {
                            // for each haplotype allele
                            // if that allele on the haplotype changes in g = [g1,g2],
                            // then we divide out the old value and multiply in the new value

                            // haplotype 0 at site j will change under this genotype
                            // therefore p_read_h0 needs to change
                            if haps[0][v] != g.0 as u8 {
                                if haps[0][v] == call.allele && g.0 != call.allele {
                                    // fragment call matched old h0 but doesn't match new h0
                                    // divide out the p(call), and multiply in p(miscall)
                                    p_read_h0 = p_read_h0 - call.one_minus_qual + call.qual;
                                } else if haps[0][v] != call.allele && g.0 == call.allele {
                                    // fragment call didn't match old h0 but matches new h0
                                    // divide out the p(miscall), and multiply in p(call)
                                    p_read_h0 = p_read_h0 - call.qual + call.one_minus_qual;
                                }
                            }

                            // haplotype 1 at site j will change under this genotype
                            // therefore p_read_h1 needs to change
                            if haps[1][v] != g.1 as u8 {
                                if haps[1][v] == call.allele && g.1 != call.allele {
                                    // fragment call matched old h1 but doesn't match new h1
                                    // divide out the p(call), and multiply in p(miscall)
                                    p_read_h1 = p_read_h1 - call.one_minus_qual + call.qual;
                                } else if haps[1][v] != call.allele && g.1 == call.allele {
                                    // fragment call didn't match old h1 but matches new h1
                                    // divide out the p(miscall), and multiply in p(call)
                                    p_read_h1 = p_read_h1 - call.qual + call.one_minus_qual;
                                }
                            }

                        } else {
                            if g.0 == call.allele  {
                                p_read_h0 = p_read_h0 + call.one_minus_qual;
                            } else {
                                p_read_h0 = p_read_h0 + call.qual;
                            }

                            if g.1 == call.allele  {
                                p_read_h1 = p_read_h1 + call.one_minus_qual;
                            } else {
                                p_read_h1 = p_read_h1 + call.qual;
                            }
                        }

                        if var_phased[v]{
                            match call.frag_ix {
                                Some(frag_ix) => { p_reads_g[g.0 as usize][g.1 as usize].push((frag_ix, p_read_h0, p_read_h1)); },
                                None => { panic!("Fragment index in pileup call is None.") },
                            }
                        }

                        let p_read = LogProb::ln_add_exp(ln_half + p_read_h0, ln_half + p_read_h1);

                        p_reads.ln_times_equals(g, p_read);
                    }
                }

                let posts: GenotypeProbs = p_reads.normalize();

                let (max_g, _) = posts.max_genotype(true, false);

                var.genotype_post = posts.clone();
                // TODO: should we reassign var.gq here?

                // we need to track if any changes occured for termination
                if haps[0][v] != max_g.0 || haps[1][v] != max_g.1 {
                    changed = true;

                    // if this variant was phased with HapCUT2 and used in calculating P(read | h),
                    // then we need to update the P(read | h1) and P(read | h2) values that changed
                    // when we changed h1 and h2
                    if var_phased[v] {
                        for &(frag_ix, p_read_h0, p_read_h1) in &p_reads_g[max_g.0 as usize][max_g.1 as usize] {
                            p_read_hap[0][frag_ix] = p_read_h0;
                            p_read_hap[1][frag_ix] = p_read_h1;
                        }
                    }
                }

                // update the haplotype vectors with the max scoring phased genotype
                haps[0][v] = max_g.0;
                haps[1][v] = max_g.1;
            }

            // if the haplotypes have not changed in this iteration, then we break
            if !changed {
                break;
            }
        }

        // count how many variants meet the criteria for "phased"
        num_phased = 0;
        for var in varlist.lst.iter() {
            if var.alleles.len() == 2 && (var.genotype == Genotype(0,1) || var.genotype == Genotype(1,0))
                && var.alleles[0].len() == 1 && var.alleles[1].len() == 1
                && var.gq > min_hap_gq {
                num_phased += 1;
            }
        }

        total_likelihood = LogProb::ln_one();

        for v in 0..varlist.lst.len() {
            let g = Genotype(haps[0][v], haps[1][v]);
            total_likelihood = total_likelihood + genotype_priors.get_prior(&varlist.lst[v].alleles, g)?;
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

            let (max_g, _) = var.genotype_post.max_genotype(true, false);

            // we computed the max phased genotype but we want the unphased genotype quality
            // sum all of the genotypes that aren't max_g, or the flipped phase version of max_g
            let mut non_max_post: Vec<LogProb> = vec![];
            for g in var.possible_genotypes() {
                if g != max_g && g != Genotype(max_g.1, max_g.0) {
                    non_max_post.push(var.genotype_post.get(g));
                }
            }
            let p_call_wrong: LogProb = LogProb::ln_sum_exp(&non_max_post);
            //let p_call_wrong: LogProb = LogProb::ln_one_minus_exp(&max_post);
            let genotype_qual: f64 = *PHREDProb::from(p_call_wrong);

            let (allele_counts, ambig_count) = count_alleles(&pileup, var.alleles.len(), max_p_miscall);
            let allele_total: usize = allele_counts.iter().sum::<usize>() + ambig_count;

            if var.dp < allele_total {
                var.dp = allele_total;
            }

            var.qual = *PHREDProb::from(var.genotype_post.get(Genotype(0,0)));
            var.allele_counts = allele_counts;
            var.ambiguous_count = ambig_count;
            var.genotype = max_g;
            var.gq = genotype_qual;
            var.filter = "PASS".to_string();
            var.called = true;

            if var.qual > MAX_VCF_QUAL {
                var.qual = MAX_VCF_QUAL;
            }

            if var.gq > MAX_VCF_QUAL {
                var.gq = MAX_VCF_QUAL;
            }
        }

        for i in 0..flist.len() {
            flist[i].p_read_hap = [p_read_hap[0][i], p_read_hap[1][i]];
        }

        let debug_vcf_str = format!("{}.{}.haplotype_genotype_iteration.vcf", program_step, hapcut2_iter).to_owned();
        print_variant_debug(varlist, &interval, &variant_debug_directory,&debug_vcf_str, max_cov, density_params, sample_name)?;

        eprintln!("{}    (After Greedy)   Total phased heterozygous SNVs: {}  Total likelihood (phred): {:.2}",print_time(), num_phased, *PHREDProb::from(total_likelihood));

        // convert logprob value to base 10
        let b10 = |x: LogProb| (*PHREDProb::from(x) / -10.0) as f64;

        if ((b10(total_likelihood) - b10(prev_likelihood)) / b10(prev_likelihood)).abs() < ll_delta {
            break;
        }

        prev_likelihood = total_likelihood;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    //use rand::Rng;

    #[test]
    fn test_generate_fragcall_pileup() {

        let fcall = |f_ix,v_ix,a| {
            FragCall {
                frag_ix: Some(f_ix), // index into fragment list
                var_ix: v_ix, // index into variant list
                allele: a, // allele call
                qual: LogProb::from(Prob(0.01)), // LogProb probability the call is an error
                one_minus_qual: LogProb::from(Prob(0.99)), // LogProb probability the call is correct
            }
        };
        let p50 = LogProb::from(Prob(0.5));
        // in this example assume the haplotype pair is (0000,1111)

        // first fragment
        let f0v0 = fcall(0,0,1);
        let f0v1 = fcall(0,1,1);
        let f0v2 = fcall(0,2,1);
        let f0v3 = fcall(0,3,1);
        let f0 = Fragment{id:"f0".to_string(),
                          calls: vec![f0v0,f0v1,f0v2,f0v3],
                          p_read_hap: [p50,p50]};
        // second fragment
        let f1v0 = fcall(1,0,0);
        let f1v1 = fcall(1,1,0);
        let f1v2 = fcall(1,2,0);
        let f1 = Fragment{id:"f1".to_string(),
            calls: vec![f1v0,f1v1,f1v2],
            p_read_hap: [p50,p50]};
        // third fragment
        let f2v1 = fcall(2,1,1);
        let f2v2 = fcall(2,2,1);
        let f2v3 = fcall(2,3,1);
        let f2 = Fragment{id:"f2".to_string(),
            calls: vec![f2v1,f2v2,f2v3],
            p_read_hap: [p50,p50]};

        // the fragment list looks like this (rows are fragments and columns are variant sites)

        // f0: 1111
        // f1: 000-
        // f2: -111

        // so the pileups should look like this
        // rows are variant sites containing alleles over that site and columns are not significant
        // p0: 10
        // p1: 101
        // p2: 101
        // p3: 11

        let expected_pileups = vec![vec![f0v0,f1v0],
                                                     vec![f0v1,f1v1,f2v1],
                                                     vec![f0v2,f1v2,f2v2],
                                                     vec![f0v3,f2v3]];

        let pileups = generate_fragcall_pileup(&vec![f0,f1,f2], 4);

        for i in 0..4 {
            assert_eq!(pileups[i].len(), expected_pileups[i].len());
            for j in 0..pileups[i].len() {
                assert_eq!(pileups[i][j].frag_ix, expected_pileups[i][j].frag_ix);
                assert_eq!(pileups[i][j].var_ix, expected_pileups[i][j].var_ix);
                assert_eq!(pileups[i][j].allele, expected_pileups[i][j].allele);
            }
        }
    }
}