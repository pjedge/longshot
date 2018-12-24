//! Structs and functions for sequence alignment (Smith-waterman and Pair-HMM forward algorithm)
// original author: Ben Pullman
// modified: Peter Edge, September 2017

use bio::stats::{PHREDProb, LogProb, Prob};
use std::collections::HashMap;

static ALLOW_END_GAPS: bool = false;

// these parameters describe state transition probabilities for a pair HMM
// there are two kinds: "eq" transition probs and "neq" transition_probs
// the correct kind to use depends on sequence context.
// the "eq" probabilities are for the case where accepting a match results in equal bases
// the "neq" probabilities are for the case where accepting a match results in different bases

#[derive(Clone)]
pub struct TransitionProbs {
    pub p: HashMap<(String, String), f64>
}

#[derive(Clone)]
pub struct LnTransitionProbs {
    pub p: HashMap<(String, String), LogProb>
}

impl TransitionProbs {
    pub fn ln(&self) -> LnTransitionProbs {
        let mut trans_probs: HashMap<(String,String), LogProb> = HashMap::new();

        for (&(ref state1, ref state2), &prob) in &self.p {
            trans_probs.insert((state1.clone(),state2.clone()), LogProb::from(Prob(prob)));
        }

        LnTransitionProbs { p: trans_probs }
    }
}

#[derive(Clone)]
pub struct EmissionProbs {
    pub probs: HashMap<(String, String), f64>
}

#[derive(Clone)]
pub struct LnEmissionProbs {
    pub probs: HashMap<(String, String), LogProb>
}

impl EmissionProbs {
    pub fn ln(&self) -> LnEmissionProbs {
        let mut probs: HashMap<(String, String), LogProb> = HashMap::new();

        for (&(ref s1, ref s2), &p) in &self.probs {
            probs.insert((s1.to_string(),s2.to_string()),LogProb::from(Prob(p)));
        }

        LnEmissionProbs { probs }
    }
}

#[derive(Clone)]
pub struct AlignmentParameters {
    pub num_states: usize,
    pub transition_probs: TransitionProbs,
    pub emission_probs: EmissionProbs,
    pub state_lst: Vec<String>,
    pub kmer_len: usize,
    pub valid_prev_states: Vec<Vec<usize>>
}

#[derive(Clone)]
pub struct LnAlignmentParameters {
    pub num_states: usize,
    pub transition_probs: LnTransitionProbs,
    pub emission_probs: LnEmissionProbs,
    pub state_lst: Vec<String>,
    pub kmer_len: usize,
    pub valid_prev_states: Vec<Vec<usize>>
}

impl AlignmentParameters {
    pub fn ln(&self) -> LnAlignmentParameters {
        LnAlignmentParameters {
            num_states: self.num_states,
            transition_probs: self.transition_probs.ln(),
            emission_probs: self.emission_probs.ln(),
            state_lst: self.state_lst.clone(),
            kmer_len: self.kmer_len,
            valid_prev_states: self.valid_prev_states.clone()
        }
    }
}

pub fn forward_algorithm_higher_order_numerically_stable(
    v: &Vec<char>,
    w: &Vec<char>,
    params: &LnAlignmentParameters,
    min_band_width: usize,
) -> (LogProb) {
    let kmer_len = params.kmer_len;
    let num_states = params.num_states;
    let starting_state = 0;

    let len_diff = ((v.len() as i32) - (w.len() as i32)).abs() as usize;
    let band_width = min_band_width + len_diff;

    let mut prev: Vec<Vec<LogProb>> = vec![];
    for _i in 0..num_states {
        prev.push(vec![LogProb::ln_zero(); w.len() + 1]);
    }

    let mut curr: Vec<Vec<LogProb>> = vec![];
    for _i in 0..num_states {
        curr.push(vec![LogProb::ln_zero(); w.len() + 1]);
    }

    curr[starting_state][kmer_len] = LogProb::ln_one();
    prev[starting_state][kmer_len] = LogProb::ln_one();

    let t = &params.transition_probs;
    let e = &params.emission_probs;

    for i in kmer_len..(v.len() + 1) {
        let band_middle = (w.len() * i) / v.len();
        let band_start = if band_middle >= band_width / 2 + kmer_len {
            band_middle - band_width / 2
        } else {
            kmer_len
        };
        let band_end = if band_middle + band_width / 2 <= w.len() {
            band_middle + band_width / 2
        } else {
            w.len()
        };

        for j in band_start..(band_end + 1) {
            for c_s in 0..num_states {
                let curr_state = params.state_lst[c_s].clone();
                let curr_state_vec: Vec<char> = curr_state.chars().collect::<Vec<char>>();
                // consider the k-th state as previous state
                let mut prob_from_prev_states = LogProb::ln_zero();

                let mut last_char = 'M';
                for char in curr_state.chars() {
                    last_char = char;
                }

                for p_s in params.valid_prev_states[c_s].iter() {
                    let prev_state = params.state_lst[*p_s].clone();

                    // this is all messed up, we actually need to trace back by the total amount of
                    // operations that were performed
                    match last_char {
                        'M' => {
                            prob_from_prev_states = LogProb::ln_add_exp(prob_from_prev_states,
                                                                        prev[*p_s][j-1]+t.p.get(&(prev_state, curr_state.clone())).unwrap());
                        },
                        'I' => {
                            prob_from_prev_states = LogProb::ln_add_exp(prob_from_prev_states,
                                                                        prev[*p_s][j]+t.p.get(&(prev_state, curr_state.clone())).unwrap());
                        },
                        'D' => {
                            prob_from_prev_states = LogProb::ln_add_exp(prob_from_prev_states,
                                                                        curr[*p_s][j-1]+t.p.get(&(prev_state, curr_state.clone())).unwrap());
                        }
                        _ => { panic!("unexpected state") }
                    }
                }

                // need to consider the current sequence indices and the current state
                // "look ahead" by the operations specified in the state{
                // build the two aligned kmers at this position for the indices and state
                // then get the emission probability for those kmers

                let mut i_prime = i;
                let mut j_prime = j;
                let mut kmer1 = vec![];
                let mut kmer2 = vec![];

                for c in curr_state_vec.iter().rev() {
                    match c {
                        'M' => {
                            kmer1.push(v[i_prime - 1]);
                            kmer2.push(w[j_prime - 1]);
                            i_prime -= 1;
                            j_prime -= 1;
                        },
                        'I' => {
                            kmer1.push(v[i_prime - 1]);
                            kmer2.push('-');
                            i_prime -= 1;
                        },
                        'D' => {
                            kmer1.push('-');
                            kmer2.push(w[j_prime - 1]);
                            j_prime -= 1;
                        }
                        _ => { panic!("invalid state in state string"); }
                    }
                }
                let kmer1_str = kmer1.iter().rev().collect::<String>();
                let kmer2_str = kmer2.iter().rev().collect::<String>();

                curr[c_s][j] = prob_from_prev_states + e.probs.get(&(kmer1_str.clone(), kmer2_str.clone())).unwrap();
            }
        }

        for j in band_start..(band_end + 1) {
            for k in 0..num_states {
                //eprintln!("{} {} {} {}",k,i,j,*PHREDProb::from(curr[k][j]));
                prev[k][j] = curr[k][j];
            }
        }


        for k in 0..num_states {
            curr[k][band_start] = LogProb::ln_zero();
        }
    }

    prev[starting_state][w.len()]

}

pub fn forward_algorithm_higher_order_non_numerically_stable(
    v: &Vec<char>,
    w: &Vec<char>,
    params: &AlignmentParameters,
    min_band_width: usize,
) -> (LogProb) {
    let kmer_len = params.kmer_len;
    let num_states = params.num_states;
    let starting_state = 0;

    let len_diff = ((v.len() as i32) - (w.len() as i32)).abs() as usize;
    let band_width = min_band_width + len_diff;

    let mut prev: Vec<Vec<f64>> = vec![];
    for _i in 0..num_states {
        prev.push(vec![0.0; w.len() + 1]);
    }

    let mut curr: Vec<Vec<f64>> = vec![];
    for _i in 0..num_states {
        curr.push(vec![0.0; w.len() + 1]);
    }

    //curr[starting_state][0] = 1.0;
    curr[starting_state][kmer_len] = 1.0;
    prev[starting_state][kmer_len] = 1.0;

    let t = &params.transition_probs;
    let e = &params.emission_probs;

    for i in kmer_len..(v.len() + 1) {
        let band_middle = (w.len() * i) / v.len();
        let band_start = if band_middle >= band_width / 2 + kmer_len {
            band_middle - band_width / 2
        } else {
            kmer_len
        };
        let band_end = if band_middle + band_width / 2 <= w.len() {
            band_middle + band_width / 2
        } else {
            w.len()
        };

        for j in band_start..(band_end + 1) {

            for c_s in 0..num_states {
                let curr_state = params.state_lst[c_s].clone();
                let curr_state_vec: Vec<char> = curr_state.chars().collect::<Vec<char>>();
                // consider the k-th state as previous state
                let mut prob_from_prev_states = 0.0;

                let mut last_char = 'M';
                for char in curr_state.chars() {
                    last_char = char;
                }

                for p_s in params.valid_prev_states[c_s].iter() {
                    let prev_state = params.state_lst[*p_s].clone();

                    // this is all messed up, we actually need to trace back by the total amount of
                    // operations that were performed
                    match last_char {
                        'M' => {
                            prob_from_prev_states += prev[*p_s][j-1] * t.p.get(&(prev_state, curr_state.clone())).unwrap();
                        },
                        'I' => {
                            prob_from_prev_states += prev[*p_s][j] * t.p.get(&(prev_state, curr_state.clone())).unwrap();
                        },
                        'D' => {
                            prob_from_prev_states += curr[*p_s][j-1] * t.p.get(&(prev_state, curr_state.clone())).unwrap();
                        }
                        _ => {panic!("unexpected state")}
                    }
                }

                // need to consider the current sequence indices and the current state
                // "look ahead" by the operations specified in the state{
                // build the two aligned kmers at this position for the indices and state
                // then get the emission probability for those kmers

                let mut i_prime = i;
                let mut j_prime = j;
                let mut kmer1 = vec![];
                let mut kmer2 = vec![];

                for c in curr_state_vec.iter().rev() {
                    match c {
                        'M' => {
                            kmer1.push(v[i_prime-1]);
                            kmer2.push(w[j_prime-1]);
                            i_prime -= 1;
                            j_prime -= 1;
                        },
                        'I' => {
                            kmer1.push(v[i_prime-1]);
                            kmer2.push('-');
                            i_prime -= 1;
                        },
                        'D' => {
                            kmer1.push('-');
                            kmer2.push(w[j_prime-1]);
                            j_prime -= 1;
                        }
                        _ => {panic!("invalid state in state string");}
                    }
                }
                let kmer1_str = kmer1.iter().rev().collect::<String>();
                let kmer2_str = kmer2.iter().rev().collect::<String>();

                curr[c_s][j] = prob_from_prev_states * e.probs.get(&(kmer1_str.clone(), kmer2_str.clone())).unwrap();
            }
        }

        for j in band_start..(band_end + 1) {
            for k in 0..num_states{
                //eprintln!("{} {} {} {}",k,i,j,*PHREDProb::from(Prob(curr[k][j])));
                prev[k][j] = curr[k][j];
            }
        }


        for k in 0..num_states{
            curr[k][band_start] = 0.0;
        }
    }

    if prev[starting_state][w.len()] != 0.0 {
        LogProb::from(Prob(prev[starting_state][w.len()]))
    } else {
        LogProb::from(Prob(1e-100))
    }
}
/*

#[cfg(test)]
mod tests {
    use super::*;

    //use rand::Rng;

    #[test]
    fn test_higher() {

    }
}
*/

