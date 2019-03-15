//! Structs and functions for sequence alignment (Smith-waterman and Pair-HMM forward algorithm)
// original author: Ben Pullman
// modified: Peter Edge, September 2017

use bio::stats::{LogProb, Prob};
use std::f64;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum AlignmentType {
    ForwardAlgorithmNonNumericallyStable,
    ForwardAlgorithmNumericallyStable,
    ViterbiMaxScoringAlignment,
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
    pub not_equal: f64,
    pub insertion: f64,
    pub deletion: f64,
}

#[derive(Clone, Copy)]
pub struct LnEmissionProbs {
    pub equal: LogProb,
    pub not_equal: LogProb,
    pub insertion: LogProb,
    pub deletion: LogProb,
}

impl EmissionProbs {
    pub fn ln(&self) -> LnEmissionProbs {
        LnEmissionProbs {
            equal: LogProb::from(Prob(self.equal)),
            not_equal: LogProb::from(Prob(self.not_equal)),
            insertion: LogProb::from(Prob(self.insertion)),
            deletion: LogProb::from(Prob(self.deletion)),
        }
    }
}

#[derive(Clone, Copy)]
pub struct AlignmentParameters {
    pub transition_probs: TransitionProbs,
    pub emission_probs: EmissionProbs,
}

#[derive(Clone, Copy)]
pub struct LnAlignmentParameters {
    pub transition_probs: LnTransitionProbs,
    pub emission_probs: LnEmissionProbs,
}

impl AlignmentParameters {
    pub fn ln(&self) -> LnAlignmentParameters {
        LnAlignmentParameters {
            transition_probs: self.transition_probs.ln(),
            emission_probs: self.emission_probs.ln(),
        }
    }
}

pub fn forward_algorithm_non_numerically_stable(
    v: &Vec<char>,
    w: &Vec<char>,
    params: AlignmentParameters,
    min_band_width: usize,
) -> (LogProb) {
    let len_diff = ((v.len() as i32) - (w.len() as i32)).abs() as usize;
    let band_width = min_band_width + len_diff;

    let mut lower_prev: Vec<f64> = vec![0.0; w.len() + 1];
    let mut middle_prev: Vec<f64> = vec![0.0; w.len() + 1];
    let mut upper_prev: Vec<f64> = vec![0.0; w.len() + 1];
    let mut lower_curr: Vec<f64> = vec![0.0; w.len() + 1];
    let mut middle_curr: Vec<f64> = vec![0.0; w.len() + 1];
    let mut upper_curr: Vec<f64> = vec![0.0; w.len() + 1];

    middle_prev[0] = 1.0;

    upper_prev[1] = params.transition_probs.deletion_from_match;
    for j in 2..(w.len() + 1) {
        upper_prev[j] = upper_prev[j - 1] * params.transition_probs.deletion_from_deletion;
        middle_prev[j] = 0.0;
    }


    let t = params.transition_probs;
    let e = params.emission_probs;

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
            upper_curr[0] = 0.0;
            middle_curr[0] = 0.0;
            if i == 1 {
                lower_curr[0] = params.transition_probs.insertion_from_match
            } else {
                lower_curr[0] =
                    lower_prev[0] * params.transition_probs.insertion_from_insertion;
            }
        }


        for j in band_start..(band_end + 1) {
            let lower_continue = lower_prev[j] * t.insertion_from_insertion;
            let lower_from_middle = middle_prev[j] * t.insertion_from_match;
            lower_curr[j] = e.insertion * (lower_continue + lower_from_middle);

            let upper_continue = upper_curr[j - 1] * t.deletion_from_deletion;
            let upper_from_middle = middle_curr[j - 1] * t.deletion_from_match;
            upper_curr[j] = e.deletion * (upper_continue + upper_from_middle);

            let middle_from_lower = lower_prev[j - 1] * t.match_from_insertion;
            let middle_continue = middle_prev[j - 1] * t.match_from_match;
            let middle_from_upper = upper_prev[j - 1] * t.match_from_deletion;
            //eprintln!("Aligning v[{}]=={} to w[{}]=={}",i-1,v[i-1],j-1,w[j-1]);

            let match_emission: f64 = if v[i - 1] == w[j - 1] {
                e.equal
            } else {
                e.not_equal
            };
            middle_curr[j] =
                match_emission * (middle_from_lower + middle_continue + middle_from_upper);
        }

        for j in (band_start-1)..(band_end + 1) {
            upper_prev[j] = upper_curr[j];
            middle_prev[j] = middle_curr[j];
            lower_prev[j] = lower_curr[j];
        }
        // we previously had a bug at the left boundary of the band... set these to NaN to make sure they
        // aren't used again
        if band_start >= 2 {
            upper_prev[band_start-2] = f64::NAN;
            middle_prev[band_start-2] = f64::NAN;
            lower_prev[band_start-2] = f64::NAN;
        }

        upper_curr[band_start] = 0.0;
        middle_curr[band_start] = 0.0;
        lower_curr[band_start] = 0.0;
    }
    //eprintln!("{:?}", middle_prev);

    if middle_prev[w.len()] != 0.0 {
        LogProb::from(Prob(middle_prev[w.len()]))
    } else {
        forward_algorithm_numerically_stable(v, w, params.ln(), band_width)
    }
}

pub fn forward_algorithm_numerically_stable(
    v: &Vec<char>,
    w: &Vec<char>,
    params: LnAlignmentParameters,
    min_band_width: usize,
) -> (LogProb) {
    let len_diff = ((v.len() as i32) - (w.len() as i32)).abs() as usize;
    let band_width = min_band_width + len_diff;

    let mut lower_prev: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut middle_prev: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut upper_prev: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut lower_curr: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut middle_curr: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut upper_curr: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];

    middle_prev[0] = LogProb::ln_one();
    let t = params.transition_probs;
    let e = params.emission_probs;

    upper_prev[1] = params.transition_probs.deletion_from_match;
    for j in 2..(w.len() + 1) {
        upper_prev[j] = upper_prev[j - 1] + params.transition_probs.deletion_from_deletion;
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
            middle_curr[0] = LogProb::ln_zero();
            if i == 1 {
                lower_curr[0] = params.transition_probs.insertion_from_match
            } else {
                lower_curr[0] =
                    lower_prev[0] + params.transition_probs.insertion_from_insertion;
            }
        }

        for j in band_start..(band_end + 1) {
            let lower_continue = lower_prev[j] + t.insertion_from_insertion;
            let lower_from_middle = middle_prev[j] + t.insertion_from_match;
            lower_curr[j] = e.insertion + LogProb::ln_add_exp(lower_continue, lower_from_middle);

            let upper_continue = upper_curr[j - 1] + t.deletion_from_deletion;
            let upper_from_middle = middle_curr[j - 1] + t.deletion_from_match;
            upper_curr[j] = e.deletion + LogProb::ln_add_exp(upper_continue, upper_from_middle);

            let middle_from_lower = lower_prev[j - 1] + t.match_from_insertion;
            let middle_continue = middle_prev[j - 1] + t.match_from_match;
            let middle_from_upper = upper_prev[j - 1] + t.match_from_deletion;
            let options3 = [middle_from_lower, middle_continue, middle_from_upper];
            let match_emission: LogProb = if v[i - 1] == w[j - 1] {
                e.equal
            } else {
                e.not_equal
            };
            middle_curr[j] = match_emission + LogProb::ln_sum_exp(&options3);
        }

        for j in (band_start-1)..(band_end + 1) {
            upper_prev[j] = upper_curr[j];
            middle_prev[j] = middle_curr[j];
            lower_prev[j] = lower_curr[j];
        }
        // we previously had a bug at the left boundary of the band... set these to NaN to make sure they
        // aren't used again
        if band_start >= 2 {
            upper_prev[band_start-2] = LogProb(f64::NAN);
            middle_prev[band_start-2] = LogProb(f64::NAN);
            lower_prev[band_start-2] = LogProb(f64::NAN);
        }

        upper_curr[band_start] = LogProb::ln_zero();
        middle_curr[band_start] = LogProb::ln_zero();
        lower_curr[band_start] = LogProb::ln_zero();
    }

    middle_prev[w.len()]
}

pub fn viterbi_max_scoring_alignment(
    v: &Vec<char>,
    w: &Vec<char>,
    params: LnAlignmentParameters,
    min_band_width: usize,
) -> LogProb {
    let len_diff = ((v.len() as i32) - (w.len() as i32)).abs() as usize;
    let band_width = min_band_width + len_diff;

    let mut lower_prev: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut middle_prev: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut upper_prev: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut lower_curr: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut middle_curr: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];
    let mut upper_curr: Vec<LogProb> = vec![LogProb::ln_zero(); w.len() + 1];

    middle_prev[0] = LogProb::ln_one();
    let t = params.transition_probs;
    let e = params.emission_probs;

    upper_prev[1] = params.transition_probs.deletion_from_match;
    for j in 2..(w.len() + 1) {
        upper_prev[j] = upper_prev[j - 1] + params.transition_probs.deletion_from_deletion;
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
            middle_curr[0] = LogProb::ln_zero();
            if i == 1 {
                lower_curr[0] = params.transition_probs.insertion_from_match
            } else {
                lower_curr[0] =
                    lower_prev[0] + params.transition_probs.insertion_from_insertion;
            }
        }


        for j in band_start..(band_end + 1) {
            let lower_continue = lower_prev[j] + t.insertion_from_insertion;
            let lower_from_middle = middle_prev[j] + t.insertion_from_match;
            lower_curr[j] = if lower_continue > lower_from_middle {
                e.insertion + lower_continue
            } else {
                e.insertion + lower_from_middle
            };

            let upper_continue = upper_curr[j - 1] + t.deletion_from_deletion;
            let upper_from_middle = middle_curr[j - 1] + t.deletion_from_match;
            upper_curr[j] = if upper_continue > upper_from_middle {
                e.deletion + upper_continue
            } else {
                e.deletion + upper_from_middle
            };

            let middle_from_lower = lower_prev[j - 1] + t.match_from_insertion;
            let middle_continue = middle_prev[j - 1] + t.match_from_match;
            let middle_from_upper = upper_prev[j - 1] + t.match_from_deletion;
            let mut max_option: LogProb = LogProb::ln_zero();
            let options: Vec<LogProb> = vec![middle_from_lower, middle_continue, middle_from_upper];
            for option in options {
                if option > max_option {
                    max_option = option;
                }
            }
            let match_emission: LogProb = if v[i - 1] == w[j - 1] {
                e.equal
            } else {
                e.not_equal
            };
            middle_curr[j] = match_emission + max_option;
        }

        for j in (band_start-1)..(band_end + 1) {
            upper_prev[j] = upper_curr[j];
            middle_prev[j] = middle_curr[j];
            lower_prev[j] = lower_curr[j];
        }
        // we previously had a bug at the left boundary of the band... set these to NaN to make sure they
        // aren't used again
        if band_start >= 2 {
            upper_prev[band_start-2] = LogProb(f64::NAN);
            middle_prev[band_start-2] = LogProb(f64::NAN);
            lower_prev[band_start-2] = LogProb(f64::NAN);
        }

        upper_curr[band_start] = LogProb::ln_zero();
        middle_curr[band_start] = LogProb::ln_zero();
        lower_curr[band_start] = LogProb::ln_zero();
    }

    middle_prev[w.len()]
}

pub fn baum_welch_non_numerically_stable(
    v: &Vec<char>,
    w: &Vec<char>,
    params: AlignmentParameters,
    min_band_width: usize,
) -> AlignmentParameters {
    let len_diff = ((v.len() as i32) - (w.len() as i32)).abs() as usize;
    let band_width = min_band_width + len_diff;

    // calculate forward probabilities

    let mut forward_lower: Vec<Vec<f64>> = vec![vec![0.0; w.len() + 1]; v.len()+1];
    let mut forward_middle: Vec<Vec<f64>> = vec![vec![0.0; w.len() + 1]; v.len()+1];
    let mut forward_upper: Vec<Vec<f64>> = vec![vec![0.0; w.len() + 1]; v.len()+1];

    forward_middle[0][0] = 1.0;

    forward_upper[0][1] = params.transition_probs.deletion_from_match;
    for j in 2..(w.len() + 1) {
        forward_upper[0][j] = forward_upper[0][j - 1] * params.transition_probs.deletion_from_deletion;
    }

    forward_lower[1][0] = params.transition_probs.insertion_from_match;
    for i in 2..(v.len() + 1) {
        forward_lower[i][0] = forward_lower[i-1][0] * params.transition_probs.insertion_from_insertion;
    }

    let t = params.transition_probs;
    let e = params.emission_probs;

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

        for j in band_start..(band_end + 1) {
            let lower_continue = forward_lower[i-1][j] * t.insertion_from_insertion;
            let lower_from_middle = forward_middle[i-1][j] * t.insertion_from_match;
            forward_lower[i][j] = e.insertion * (lower_continue + lower_from_middle);

            let upper_continue = forward_upper[i][j - 1] * t.deletion_from_deletion;
            let upper_from_middle = forward_middle[i][j - 1] * t.deletion_from_match;
            forward_upper[i][j] = e.deletion * (upper_continue + upper_from_middle);

            let middle_from_lower = forward_lower[i-1][j - 1] * t.match_from_insertion;
            let middle_continue = forward_middle[i-1][j - 1] * t.match_from_match;
            let middle_from_upper = forward_upper[i-1][j - 1] * t.match_from_deletion;

            let match_emission: f64 = if v[i - 1] == w[j - 1] {
                e.equal
            } else {
                e.not_equal
            };
            forward_middle[i][j] =
                match_emission * (middle_from_lower + middle_continue + middle_from_upper);
        }
    }


    // calculate backward probabilities

    let mut backward_lower: Vec<Vec<f64>> = vec![vec![0.0; w.len() + 1]; v.len()+1];
    let mut backward_middle: Vec<Vec<f64>> = vec![vec![0.0; w.len() + 1]; v.len()+1];
    let mut backward_upper: Vec<Vec<f64>> = vec![vec![0.0; w.len() + 1]; v.len()+1];

    backward_middle[v.len()][w.len()] = 1.0;

    backward_upper[v.len()][w.len()] = params.transition_probs.deletion_from_deletion;
    for j in (0..(w.len())).rev() {
        backward_upper[v.len()][j] = backward_upper[v.len()][j+1] * params.transition_probs.deletion_from_deletion;
    }

    backward_lower[v.len()][w.len()] = params.transition_probs.insertion_from_insertion;
    for i in (0..(v.len())).rev() {
        backward_lower[i][w.len()] = backward_lower[i+1][w.len()] * params.transition_probs.insertion_from_insertion;
    }

    let t = params.transition_probs;
    let e = params.emission_probs;

    for i in (0..v.len()).rev() {
        let band_middle = (w.len() * i) / v.len();
        let band_start = if band_middle >= band_width / 2 {
            band_middle - band_width / 2
        } else {
            0
        };
        let band_end = if band_middle + band_width / 2 < w.len() {
            band_middle + band_width / 2
        } else {
            w.len()-1
        };

        for j in (band_start..(band_end + 1)).rev() {
            let lower_continue = backward_lower[i+1][j] * t.insertion_from_insertion;
            let lower_from_middle = backward_middle[i+1][j] * t.match_from_insertion;
            backward_lower[i][j] = e.insertion * (lower_continue + lower_from_middle);

            let upper_continue = backward_upper[i][j+1] * t.deletion_from_deletion;
            let upper_from_middle = backward_middle[i][j+1] * t.match_from_deletion;
            backward_upper[i][j] = e.deletion * (upper_continue + upper_from_middle);

            let middle_from_lower = backward_lower[i+1][j+1] * t.insertion_from_match;
            let middle_continue = backward_middle[i+1][j+1] * t.match_from_match;
            let middle_from_upper = backward_upper[i+1][j+1] * t.deletion_from_match;

            let match_emission: f64 = if v[i] == w[j] {
                e.equal
            } else {
                e.not_equal
            };
            backward_middle[i][j] =
                match_emission * (middle_from_lower + middle_continue + middle_from_upper);
        }
    }

    eprintln!("forward_lower:");
    for i in 0..(v.len()+1) {
        eprintln!("{:?}",forward_lower[i]);
    }

    eprintln!("forward_middle:");
    for i in 0..(v.len()+1) {
        eprintln!("{:?}",forward_middle[i]);
    }
    eprintln!("forward_upper:");
    for i in 0..(v.len()+1) {
        eprintln!("{:?}",forward_upper[i]);
    }

    eprintln!("\n\nbackward_lower:");
    for i in 0..(v.len()+1) {
        eprintln!("{:?}",backward_lower[i]);
    }

    eprintln!("backward_middle:");
    for i in 0..(v.len()+1) {
        eprintln!("{:?}",backward_middle[i]);
    }
    eprintln!("backward_upper:");
    for i in 0..(v.len()+1) {
        eprintln!("{:?}",backward_upper[i]);
    }

    let x = forward_lower[v.len()][w.len()]*params.transition_probs.match_from_insertion+
    forward_middle[v.len()][w.len()]*params.transition_probs.match_from_match+
        forward_upper[v.len()][w.len()]*params.transition_probs.deletion_from_insertion;

    let y = backward_lower[0][0]*params.transition_probs.insertion_from_match+
        backward_middle[0][0]*params.transition_probs.match_from_match+
        backward_upper[0][0]*params.transition_probs.deletion_from_match;

    println!("{} {}", x, y);

    //(forward_middle[v.len()][w.len()],backward_middle[0][0])
    params.clone()
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_baumwelch() {

        let v = vec!['A','C','T','G'];
        let w = vec!['A','C','C','G'];
        let alignment_parameters = AlignmentParameters {
            emission_probs: EmissionProbs {equal: 0.982, not_equal: 0.006, insertion: 1.0, deletion:1.0},
            transition_probs: TransitionProbs {match_from_match: 0.879,
                insertion_from_match: 0.06,
                deletion_from_match: 0.06,
                match_from_deletion: 0.85,
                deletion_from_deletion: 0.15,
                match_from_insertion: 0.85,
                insertion_from_insertion: 0.15}
        };

        baum_welch_non_numerically_stable(&v,&w,alignment_parameters,20);
        assert!(true);
    }
}