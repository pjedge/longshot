// original author: Ben Pullman
// modified: Peter Edge, September 2017

use std::f64::MIN;
use std::cmp::max;
use bio::stats::{LogProb, Prob};
use util::{AlignmentParameters, LnAlignmentParameters};

pub fn sum_all_alignments(v: &Vec<char>,
                          w: &Vec<char>,
                          params: AlignmentParameters,
                          min_band_width: usize)
                          -> (LogProb) {

    let len_diff = ((v.len() as i32) - (w.len() as i32)).abs() as usize;
    let band_width = min_band_width + len_diff;

    let mut lower_prev: Vec<f64> = vec![0.0; w.len() + 1];
    let mut middle_prev: Vec<f64> = vec![0.0; w.len() + 1];
    let mut upper_prev: Vec<f64> = vec![0.0; w.len() + 1];
    let mut lower_curr: Vec<f64> = vec![0.0; w.len() + 1];
    let mut middle_curr: Vec<f64> = vec![0.0; w.len() + 1];
    let mut upper_curr: Vec<f64> = vec![0.0; w.len() + 1];

    middle_prev[0] = 1.0;
    upper_prev[1] = params.deletion_from_match;
    //lower_prev[1] = params.deletion_from_match;

    for j in 2..(w.len() + 1) {
        upper_prev[j] = upper_prev[j - 1] * params.extend_from_deletion;
        middle_prev[j] = 0.0;
        //lower_prev[j] = 0.0;
    }

    for i in 1..(v.len() + 1) {

        let band_middle = (w.len() * i) / v.len();
        let band_start = if band_middle >= band_width / 2 + 1 {
            band_middle - band_width / 2
        } else {
            1
        };
        let band_end = if band_middle + band_width / 2 <= w.len() {
            band_middle + band_width / 2
        } else {
            w.len()
        };

        if band_start == 1 {
            upper_curr[0] = 0.0;;
            middle_curr[0] = 0.0;
            lower_curr[0] = lower_prev[0] + params.extend_from_insertion;
        }

        for j in band_start..(band_end + 1) {
            //for j in 1..(w.len() + 1) {
            let score = match v.get(i - 1) {
                Some(v_i) => {
                    match w.get(j - 1) {
                        Some(w_j) => {
                            if v_i == w_j {
                                params.match_from_match
                            } else {
                                params.mismatch_from_match
                            }
                        }
                        None => panic!("w_j shouldn't be empty"),
                    }
                }
                None => panic!("v_i shouldn't be empty"),
            };

            let lower_continue = lower_prev[j] * params.extend_from_insertion;
            let lower_from_middle = middle_prev[j] * params.insertion_from_match;
            lower_curr[j] = lower_continue + lower_from_middle;

            let upper_continue = upper_curr[j - 1] * params.extend_from_deletion;
            let upper_from_middle = middle_curr[j - 1] * params.deletion_from_match;
            upper_curr[j] = upper_continue + upper_from_middle;

            let middle_from_lower = lower_curr[j];
            let middle_continue = middle_prev[j - 1] * score;
            let middle_from_upper = upper_curr[j];
            middle_curr[j] = middle_from_lower + middle_continue + middle_from_upper;
        }

        for j in band_start..(band_end + 1) {
            upper_prev[j] = upper_curr[j];
            middle_prev[j] = middle_curr[j];
            lower_prev[j] = lower_curr[j];
        }

        upper_curr[band_start] = 0.0;
        middle_curr[band_start] = 0.0;
        lower_curr[band_start] = 0.0;

    }

    if middle_prev[w.len()] != 0.0 {
        LogProb::from(Prob(middle_prev[w.len()]))
    } else {
        sum_all_alignments_numerically_stable(v, w, params.ln(), band_width)
    }
}

pub fn sum_all_alignments_numerically_stable(v: &Vec<char>,
                                             w: &Vec<char>,
                                             params: LnAlignmentParameters,
                                             min_band_width: usize)
                                             -> (LogProb) {

    let len_diff = ((v.len() as i32) - (w.len() as i32)).abs() as usize;
    let band_width = min_band_width + len_diff;

    let mut lower_prev: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut middle_prev: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut upper_prev: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut lower_curr: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut middle_curr: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut upper_curr: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];

    middle_prev[0] = LogProb::ln_one();
    upper_prev[1] = params.deletion_from_match;

    for j in 2..(w.len() + 1) {
        upper_prev[j] = upper_prev[j - 1] + params.extend_from_deletion;
        //middle_prev[j] = LogProb::ln_zero();
        //lower_prev[j] = lower_prev[j - 1] + gap_extend_penalty;
    }

    for i in 1..(v.len() + 1) {

        let band_middle = (w.len() * i) / v.len();
        let band_start = if band_middle >= band_width / 2 + 1 {
            band_middle - band_width / 2
        } else {
            1
        };
        let band_end = if band_middle + band_width / 2 <= w.len() {
            band_middle + band_width / 2
        } else {
            w.len()
        };

        if band_start == 1 {
            //upper_curr[0] = upper_prev[0] + gap_extend_penalty;
            middle_curr[0] = LogProb::ln_zero();
            lower_curr[0] = lower_prev[0] + params.extend_from_insertion;
        }

        for j in band_start..(band_end + 1) {
            //for j in 1..(w.len() + 1) {
            let score = match v.get(i - 1) {
                Some(v_i) => {
                    match w.get(j - 1) {
                        Some(w_j) => {
                            if v_i == w_j {
                                params.match_from_match
                            } else {
                                params.mismatch_from_match
                            }
                        }
                        None => panic!("w_j shouldn't be empty"),
                    }
                }
                None => panic!("v_i shouldn't be empty"),
            };

            let lower_continue = lower_prev[j] + params.extend_from_insertion;
            let lower_from_middle = middle_prev[j] + params.insertion_from_match;
            lower_curr[j] = LogProb::ln_add_exp(lower_continue, lower_from_middle);

            let upper_continue = upper_curr[j - 1] + params.extend_from_deletion;
            let upper_from_middle = middle_curr[j - 1] + params.deletion_from_match;
            upper_curr[j] = LogProb::ln_add_exp(upper_continue, upper_from_middle);

            let middle_from_lower = lower_curr[j];
            let middle_continue = middle_prev[j - 1] + score;
            let middle_from_upper = upper_curr[j];
            let options3 = [middle_from_lower, middle_continue, middle_from_upper];
            middle_curr[j] = LogProb::ln_sum_exp(&options3);
        }

        for j in band_start..(band_end + 1) {
            upper_prev[j] = upper_curr[j];
            middle_prev[j] = middle_curr[j];
            lower_prev[j] = lower_curr[j];
        }

        upper_curr[band_start] = LogProb::ln_zero();
        middle_curr[band_start] = LogProb::ln_zero();
        lower_curr[band_start] = LogProb::ln_zero();

    }

    middle_prev[w.len()]
}
