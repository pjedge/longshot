
use bio::stats::{PHREDProb, LogProb, Prob};
use util::{PotentialVar, VarList, Fragment, FragCall, dna_vec, parse_target_names};


static MAX_P_MISCALL_F64: f64 = 0.25;

fn generate_realigned_pileup(flist: Vec<Fragment>, n_var: usize) -> Vec<Vec<FragCall>> {
    let mut pileup_lst: Vec<Vec<FragCall>> = vec![vec![]; n_var];

    for fragment in flist {
        for call in fragment.calls {
            if call.qual < LogProb::from(Prob(MAX_P_MISCALL_F64)) {
                pileup_lst[call.var_ix].push(call);
            }
        }
    }
    pileup_lst
}

fn count_alleles(pileup: &Vec<FragCall>) -> (usize, usize) {

    // compute probabilities of data given each genotype
    let mut count0 = 0;
    let mut count1 = 0;

    for call in pileup {

        match call.allele {
            '0' => {
                count0 += 1;
            }
            '1' => {
                count1 += 1;
            }
            _ => {
                panic!("Unexpected allele observed in pileup.");
            }
        }
    }

    (count0, count1)
}


fn compute_posteriors(pileup: &Vec<FragCall>) -> (LogProb, LogProb, LogProb) {

    // compute probabilities of data given each genotype
    let mut prob00 = LogProb::ln_one();
    let mut prob01 = LogProb::ln_one();
    let mut prob11 = LogProb::ln_one();
    let ln_half = LogProb::from(Prob(0.5));

    for call in pileup {

        let p_call = LogProb::ln_sub_exp(LogProb::ln_one(), call.qual);
        let p_miscall = call.qual;

        match call.allele {
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

pub fn call_genotype(pileup: Vec<FragCall>, var: PotentialVar) {

    let (post00, post01, post11): (LogProb, LogProb, LogProb) = compute_posteriors(&pileup);

    let mut genotype: String = "./.".to_string();
    let mut genotype_qual: f64 = 0.0;

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

    let (count_ref, count_var): (usize, usize) = count_alleles(&pileup);

    println!("{}\t{}\t.\t{}\t{}\t.\t.\t00={};01={};11={};RA={};AA={}\tGT:GQ\t{}:{}",
             var.chrom,
             var.pos0 + 1,
             var.ref_allele,
             var.var_allele,
             *PHREDProb::from(post00),
             *PHREDProb::from(post01),
             *PHREDProb::from(post11),
             count_ref,
             count_var,
             genotype,
             genotype_qual);

}

pub fn call_genotypes(flist: Vec<Fragment>, varlist: &VarList) {

    let pileup_lst = generate_realigned_pileup(flist, varlist.lst.len());

    assert_eq!(pileup_lst.len(), varlist.lst.len());

    for i in 0..varlist.lst.len() {

        call_genotype(pileup_lst[i].clone(), varlist.lst[i].clone());

    }
}
