
use bio::stats::{PHREDProb, LogProb, Prob};
use util::{VarList, Fragment, FragCall, GenomicInterval};
use std::error::Error;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;

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

pub fn call_genotypes(flist: Vec<Fragment>,
                      varlist: &VarList,
                      interval: &Option<GenomicInterval>,
                      output_vcf_file: String) {

    let vcf_path = Path::new(&output_vcf_file);
    let vcf_display = vcf_path.display();

    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match File::create(&vcf_path) {
        Err(why) => panic!("couldn't create {}: {}", vcf_display, why.description()),
        Ok(file) => file,
    };

    let pileup_lst = generate_realigned_pileup(flist, varlist.lst.len());

    assert_eq!(pileup_lst.len(), varlist.lst.len());

    for i in 0..varlist.lst.len() {

        let pileup = &pileup_lst[i];
        let var = &varlist.lst[i];

        match interval {
            &Some(ref iv) => {
                if var.pos0 < iv.start_pos as usize || var.pos0 > iv.end_pos as usize {
                    continue;
                }
            }
            &None => {}
        }

        let (post00, post01, post11): (LogProb, LogProb, LogProb) = compute_posteriors(&pileup);

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

        let (count_ref, count_var): (usize, usize) = count_alleles(&pileup);

        // Write the `LOREM_IPSUM` string to `file`, returns `io::Result<()>`
        match writeln!(file,
                       "{}\t{}\t.\t{}\t{}\t.\t.\tRA={};AA={}\tGT:GQ\t{}:{}",
                       var.chrom,
                       var.pos0 + 1,
                       var.ref_allele,
                       var.var_allele,
                       count_ref,
                       count_var,
                       genotype,
                       genotype_qual) {
            Err(why) => panic!("couldn't write to {}: {}", vcf_display, why.description()),
            Ok(_) => {}
        }

    }
}
