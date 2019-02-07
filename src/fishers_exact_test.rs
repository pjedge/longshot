// Copyright 2018 fishers_exact Developers
//
// Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
// http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
// http://opensource.org/licenses/MIT>, at your option. This file may not be
// copied, modified, or distributed except according to those terms.

//! Fisher's exact test.
//!
//! Implements a 2×2 Fishers exact test. Use this to test the independence of two
//! categorical variables when the sample sizes are small.
//!
//! For an approachable explanation of Fisher's exact test, see
//! [Fisher's exact test of independence](http://www.biostathandbook.com/fishers.html) by
//! John H. McDonald in the [Handbook of Biological Statistics](http://www.biostathandbook.com/).
//!
//! The test is computed using code ported from Øyvind Langsrud's JavaScript
//! implementation at [http://www.langsrud.com/fisher.htm](http://www.langsrud.com/fisher.htm),
//! used with permission.

use std::fmt;

fn lngamm(z: i32) -> f64
// Reference: "Lanczos, C. 'A precision approximation
// of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
// Translation of  Alan Miller's FORTRAN-implementation
// See http://lib.stat.cmu.edu/apstat/245
{
    let z = z as f64;
    let mut x = 0.0;
    x += 0.1659470187408462e-06 / (z + 7.0);
    x += 0.9934937113930748e-05 / (z + 6.0);
    x -= 0.1385710331296526 / (z + 5.0);
    x += 12.50734324009056 / (z + 4.0);
    x -= 176.6150291498386 / (z + 3.0);
    x += 771.3234287757674 / (z + 2.0);
    x -= 1259.139216722289 / (z + 1.0);
    x += 676.5203681218835 / (z);
    x += 0.9999999999995183;
    x.ln() - 5.58106146679532777 - z + (z - 0.5) * (z + 6.5).ln()
}

fn lnfact(n: i32) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    lngamm(n + 1)
}

fn lnbico(n: i32, k: i32) -> f64 {
    lnfact(n) - lnfact(k) - lnfact(n - k)
}

fn hyper_323(n11: i32, n1_: i32, n_1: i32, n: i32) -> f64 {
    (lnbico(n1_, n11) + lnbico(n - n1_, n_1 - n11) - lnbico(n, n_1)).exp()
}

fn hyper(s: &mut HyperState, n11: i32) -> f64 {
    hyper0(s, n11, 0, 0, 0)
}

struct HyperState {
    n11: i32,
    n1_: i32,
    n_1: i32,
    n: i32,
    prob: f64,
    valid: bool,
}

impl HyperState {
    fn new() -> HyperState {
        HyperState {
            n11: 0,
            n1_: 0,
            n_1: 0,
            n: 0,
            prob: 0.0,
            valid: false,
        }
    }
}

fn hyper0(s: &mut HyperState, n11i: i32, n1_i: i32, n_1i: i32, ni: i32) -> f64 {
    if s.valid && (n1_i | n_1i | ni) == 0 {
        if !(n11i % 10 == 0) {
            if n11i == s.n11 + 1 {
                s.prob *= ((s.n1_ - s.n11) as f64 / n11i as f64)
                    * ((s.n_1 - s.n11) as f64 / (n11i + s.n - s.n1_ - s.n_1) as f64);
                s.n11 = n11i;
                return s.prob;
            }
            if n11i == s.n11 - 1 {
                s.prob *= ((s.n11 as f64) / (s.n1_ - n11i) as f64)
                    * ((s.n11 + s.n - s.n1_ - s.n_1) as f64 / (s.n_1 - n11i) as f64);
                s.n11 = n11i;
                return s.prob;
            }
        }
        s.n11 = n11i;
    } else {
        s.n11 = n11i;
        s.n1_ = n1_i;
        s.n_1 = n_1i;
        s.n = ni;
        s.valid = true
    }
    s.prob = hyper_323(s.n11, s.n1_, s.n_1, s.n);
    return s.prob;
}

// Returns prob,sleft,sright,sless,slarg
fn exact(n11: i32, n1_: i32, n_1: i32, n: i32) -> (f64, f64, f64, f64, f64) {
    let mut sleft: f64;
    let mut sright: f64;
    let sless: f64;
    let slarg: f64;
    let mut p: f64;
    let mut i;
    let mut j;
    let prob: f64;
    let mut max = n1_;
    if n_1 < max {
        max = n_1;
    }
    let mut min = n1_ + n_1 - n;
    if min < 0 {
        min = 0;
    }
    if min == max {
        return (1.0, 1.0, 1.0, 1.0, 1.0);
    }
    let mut s = HyperState::new();
    prob = hyper0(&mut s, n11, n1_, n_1, n);
    sleft = 0.0;
    p = hyper(&mut s, min);
    i = min + 1;
    while p <= 0.99999999 * prob {
        sleft += p;
        p = hyper(&mut s, i);
        i += 1;
    }
    i -= 1;
    if p <= 1.00000001 * prob {
        sleft += p;
    } else {
        i += 1;
    }
    sright = 0.0;
    p = hyper(&mut s, max);
    j = max - 1;
    while p <= 0.99999999 * prob {
        sright += p;
        p = hyper(&mut s, j);
        j -= 1;
    }
    j += 1;
    if p <= 1.00000001 * prob {
        sright += p;
    } else {
        j += 1;
    }
    if (i - n11).abs() < (j - n11).abs() {
        sless = sleft;
        slarg = 1.0 - sleft + prob;
    } else {
        sless = 1.0 - sright + prob;
        slarg = sright;
    }
    return (prob, sleft, sright, sless, slarg);
}

/// `FishersExactPvalues` holds the pvalues calculated by the `fishers_exact` function.
#[derive(Clone, Copy, Debug)]
pub struct FishersExactPvalues {
    /// pvalue for the two-tailed test. Use this when there is no prior alternative.
    pub two_tail_pvalue: f64,
    /// pvalue for the "left" or "lesser" tail. Use this when the alternative to
    /// independence is that there is negative association between the variables.
    /// That is, the observations tend to lie in lower left and upper right.
    pub less_pvalue: f64,
    /// Use this when the alternative to independence is that there is positive
    /// association between the variables. That is, the observations tend to lie
    /// in upper left and lower right.
    pub greater_pvalue: f64,
}

impl fmt::Display for FishersExactPvalues {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "(less_pvalue={}, greater_pvalue={}, two_tail_pvalue={})",
            self.less_pvalue, self.greater_pvalue, self.two_tail_pvalue
        )
    }
}


fn exact22(n11: i32, n12: i32, n21: i32, n22: i32) -> FishersExactPvalues {
    let (left, right, mut twotail);

    let n1_ = n11 + n12;
    let n_1 = n11 + n21;
    let n = n11 + n12 + n21 + n22;
    let (_, sleft, sright, sless, slarg) = exact(n11, n1_, n_1, n);
    left = sless;
    right = slarg;
    twotail = sleft + sright;
    if twotail > 1.0 {
        twotail = 1.0;
    }

    return FishersExactPvalues {
        two_tail_pvalue: twotail,
        less_pvalue: left,
        greater_pvalue: right,
    };
}

/// Computes the Fisher's exact pvales to determine if there are nonrandom associations between two
/// categorical variables, in a two by two contingency table.
///
/// The test is computed using code ported from Øyvind Langsrud's JavaScript
/// implementation at [http://www.langsrud.com/fisher.htm](http://www.langsrud.com/fisher.htm).
///
/// Use this when sample sizes are small. For large samples, other statistical tests of independence
/// are more appropriate.
///
/// # Examples
/// ```
/// use fishers_exact::fishers_exact;
///
/// let p = fishers_exact(&[1,9,11,3]).unwrap();
///
/// assert!((p.less_pvalue - 0.001346).abs() < 0.0001);
/// assert!((p.greater_pvalue - 0.9999663).abs() < 0.0001);
/// assert!((p.two_tail_pvalue - 0.0027594).abs() < 0.0001);
/// ```
///
/// # Errors
/// Returns `TooLargeValueError` if any member in the `table` array is too large. Currently
/// "too large" is defined as greater than `std::i32::MAX`.
pub fn fishers_exact(table: &[u32; 4]) -> FishersExactPvalues {
    exact22(table[0] as i32, table[1] as i32, table[2] as i32, table[3] as i32)
}

#[cfg(test)]
mod tests {
    use super::fishers_exact;

    fn fuzzy_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < 0.000001
    }

    #[test]
    fn test_fishers_exact() {
        // 20 cases randomly generated via scipy.
        // ([a,b,c,d], less, greater, two-tail)
        let cases = [
            (
                [61, 118, 2, 1],
                0.27535061623455315,
                0.9598172545684959,
                0.27535061623455315,
            ),
            (
                [172, 46, 90, 127],
                1.0,
                6.662405187351769e-16,
                9.041009036528785e-16,
            ),
            (
                [127, 38, 112, 43],
                0.8637599357870167,
                0.20040942958644145,
                0.3687862842650179,
            ),
            (
                [186, 177, 111, 154],
                0.9918518696328176,
                0.012550663906725129,
                0.023439141644624434,
            ),
            (
                [137, 49, 135, 183],
                0.999999999998533,
                5.6517533666400615e-12,
                8.870999836202932e-12,
            ),
            (
                [37, 115, 37, 152],
                0.8834621182590621,
                0.17638403366123565,
                0.29400927608021704,
            ),
            (
                [124, 117, 119, 175],
                0.9956704915461392,
                0.007134712391455461,
                0.011588218284387445,
            ),
            (
                [70, 114, 41, 118],
                0.9945558498544903,
                0.010384865876586255,
                0.020438291037108678,
            ),
            (
                [173, 21, 89, 7],
                0.2303739114068352,
                0.8808002774812677,
                0.4027047267306024,
            ),
            (
                [18, 147, 123, 58],
                4.077820702304103e-29,
                0.9999999999999817,
                7.686224774594537e-29,
            ),
            (
                [116, 20, 92, 186],
                0.9999999999998267,
                6.598118571034892e-25,
                8.164831402188242e-25,
            ),
            (
                [9, 22, 44, 38],
                0.01584272038710196,
                0.9951463496539362,
                0.021581786662999272,
            ),
            (
                [9, 101, 135, 7],
                3.3336213533847776e-50,
                1.0,
                3.3336213533847776e-50,
            ),
            (
                [153, 27, 191, 144],
                0.9999999999950817,
                2.473736787266208e-11,
                3.185816623300107e-11,
            ),
            (
                [111, 195, 189, 69],
                6.665245982898848e-19,
                0.9999999999994574,
                1.0735744915712542e-18,
            ),
            (
                [125, 21, 31, 131],
                0.99999999999974,
                9.720661317939016e-34,
                1.0352129312860277e-33,
            ),
            (
                [201, 192, 69, 179],
                0.9999999988714893,
                3.1477232259550017e-09,
                4.761075937088169e-09,
            ),
            (
                [167, 184, 141, 28],
                7.045789653297585e-16,
                1.0,
                9.362858503272341e-16,
            ),
            (
                [194, 74, 141, 182],
                0.9999999999999848,
                1.2268868025030845e-12,
                1.8076995960009742e-12,
            ),
            (
                [124, 138, 159, 160],
                0.30153826772785475,
                0.7538974235759873,
                0.5601766196310243,
            ),
        ];

        for &(table, expected_left, expected_right, expected_two_tails) in cases.iter() {
            let p = fishers_exact(&table).unwrap();
            println!(
                "{:?} expect={},{},{} observed={},{},{}",
                table,
                expected_left,
                expected_right,
                expected_two_tails,
                p.less_pvalue,
                p.greater_pvalue,
                p.two_tail_pvalue
            );
            assert!(fuzzy_eq(p.less_pvalue, expected_left));
            assert!(fuzzy_eq(p.greater_pvalue, expected_right));
            assert!(fuzzy_eq(p.two_tail_pvalue, expected_two_tails));
        }

        match fishers_exact(&[std::i32::MAX as u32 + 1, 1, 1, 1]) {
            Err(e) => println!("Error: {}", e),
            _ => assert!(false), // Should have errored.
        }
    }
}
