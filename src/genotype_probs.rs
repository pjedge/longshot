//! Contains structs for representing tables of genotype probabilities, as well as estimation of
//!  genotype priors.

use bio::stats::*;
use errors::*;
use hashbrown::HashMap;
use util::*;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Genotype(pub u8, pub u8);

#[derive(Debug, Clone)]
pub struct GenotypeProbs {
    pub tab: Vec<Vec<LogProb>>,
}

impl GenotypeProbs {
    pub fn zeros(n_alleles: usize) -> GenotypeProbs {
        GenotypeProbs {
            tab: vec![vec![LogProb::ln_zero(); n_alleles]; n_alleles],
        }
    }

    pub fn ones(n_alleles: usize) -> GenotypeProbs {
        GenotypeProbs {
            tab: vec![vec![LogProb::ln_one(); n_alleles]; n_alleles],
        }
    }

    pub fn uniform(n_alleles: usize) -> GenotypeProbs {
        GenotypeProbs {
            tab: vec![
                vec![LogProb::from(Prob(1.0 / ((n_alleles * n_alleles) as f64))); n_alleles];
                n_alleles
            ],
        }
    }

    pub fn n_alleles(&self) -> usize {
        self.tab.len()
    }

    pub fn max_prob(&self) -> (Genotype, LogProb) {
        let mut max_post: LogProb = LogProb::ln_zero();
        let mut max_i = 0;
        let mut max_j = 0;

        for i in 0..self.n_alleles() {
            for j in 0..self.n_alleles() {
                let post: LogProb = self.tab[i][j];

                if post > max_post {
                    max_post = post;
                    max_i = i;
                    max_j = j;
                }
            }
        }

        (Genotype(max_i as u8, max_j as u8), max_post)
    }

    pub fn max_genotype_post(&self, phased: bool, force_nonreference: bool) -> (Genotype, LogProb) {
        let mut max_post: LogProb = LogProb::ln_zero();
        let mut max_i = 0;
        let mut max_j = 0;

        for i in 0..self.n_alleles() {
            for j in 0..self.n_alleles() {
                if i == 0 && j == 0 && force_nonreference {
                    continue;
                }

                let post: LogProb = if i == j || phased {
                    self.tab[i][j]
                } else {
                    LogProb::ln_add_exp(self.tab[i][j], self.tab[j][i])
                };

                if post > max_post {
                    max_post = post;
                    max_i = i;
                    max_j = j;
                }
            }
        }

        if max_post > LogProb::ln_one() {
            let err: LogProb = LogProb::ln_sub_exp(max_post, LogProb::ln_one());
            assert!(err < LogProb::from(Prob(0.00001)));
            max_post = LogProb::ln_one();
        }

        assert!(max_post <= LogProb::ln_one());
        assert!(max_post > LogProb::ln_zero());

        (Genotype(max_i as u8, max_j as u8), max_post)
    }

    pub fn sum(&self) -> LogProb {
        let mut total: LogProb = LogProb::ln_zero();
        for row in &self.tab {
            total = LogProb::ln_add_exp(total, LogProb::ln_sum_exp(row));
        }
        //let mut flattened = vec![];
        /*
        for i in 0..self.n_alleles() {
            for j in 0..self.n_alleles() {
                //flattened.push(self.tab[i][j]);
                total = LogProb::ln_add_exp(total, self.tab[i][j]);
            }
        }
        */

        //LogProb::ln_sum_exp(&flattened)
        total
    }

    pub fn sum_genotypes_with_allele(&self, a: u8) -> LogProb {
        let i: usize = a as usize;
        let mut total: LogProb = self.tab[i][i];

        for j in 0..self.n_alleles() {
            if i == j {
                continue;
            }
            total = LogProb::ln_add_exp(total, self.tab[i][j]);
            total = LogProb::ln_add_exp(total, self.tab[j][i]);
        }

        total
    }

    pub fn normalize(&mut self) -> GenotypeProbs {
        let total: LogProb = self.sum();
        let mut norm: GenotypeProbs = GenotypeProbs::zeros(self.n_alleles());

        for i in 0..self.n_alleles() {
            for j in 0..self.n_alleles() {
                norm.tab[i][j] = self.tab[i][j] - total;
                assert!(norm.tab[i][j] <= LogProb::ln_one());
                assert!(norm.tab[i][j] >= LogProb::ln_zero());
            }
        }

        norm.assert_approx_normalized();
        norm
    }

    pub fn assert_approx_normalized(&self) {
        let total: LogProb = self.sum();
        let err: LogProb = if total > LogProb::ln_one() {
            LogProb::ln_sub_exp(total, LogProb::ln_one())
        } else {
            LogProb::ln_sub_exp(LogProb::ln_one(), total)
        };

        let margin = 0.00001;
        if err > LogProb::from(Prob(margin)) {
            println!("ERROR");
            println!("{:?}", self.tab);
            println!("{}", *Prob::from(total))
        }
        assert!(err < LogProb::from(Prob(margin)));
    }

    pub fn ln_times_equals(&mut self, g: Genotype, p: LogProb) {
        let g0 = g.0 as usize;
        let g1 = g.1 as usize;
        self.tab[g0][g1] = self.tab[g0][g1] + p;
    }

    pub fn set(&mut self, g: Genotype, p: LogProb) {
        self.tab[g.0 as usize][g.1 as usize] = p;
    }

    pub fn get(&self, g: Genotype) -> LogProb {
        self.tab[g.0 as usize][g.1 as usize]
    }

    pub fn print_phred(&self) {
        for i in 0..self.n_alleles() {
            for j in 0..self.n_alleles() {
                print!("{} ", *PHREDProb::from(self.tab[i][j]));
            }
            println!("");
        }
        println!("------------------------");
    }

    pub fn print_prob(&self) {
        for i in 0..self.n_alleles() {
            for j in 0..self.n_alleles() {
                print!("{} ", *Prob::from(self.tab[i][j]));
            }
            println!("");
        }
        println!("------------------------");
    }
}

pub fn possible_genotypes(alleles: &Vec<String>) -> Vec<Genotype> {
    let n_alleles = alleles.len();
    let mut genotypes = Vec::with_capacity(n_alleles * n_alleles);
    for i in 0..n_alleles {
        for j in 0..n_alleles {
            genotypes.push(Genotype(i as u8, j as u8));
        }
    }
    return genotypes;
}

#[derive(Clone)]
pub struct GenotypePriors {
    priors_dict: HashMap<(char, (char, char)), LogProb>, // (ref_allele, (allele1, allele2)) -> P(G)
}

impl GenotypePriors {
    // create a struct representing human variant priors
    // estimate prior probability of genotypes using strategy described here:
    // http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694485/
    // "prior probability of each genotype"
    // we've modified it in this case to have two extra options: 'I' and 'D'
    // these represent short insertions and short deletions.
    pub fn new(
        hom_snv_rate: LogProb,
        het_snv_rate: LogProb,
        hom_indel_rate: LogProb,
        het_indel_rate: LogProb,
        ts_tv_ratio: f64,
    ) -> Result<GenotypePriors> {
        //let hom_snv_rate = LogProb::from(Prob(0.0005));
        //let het_snv_rate = LogProb::from(Prob(0.001));
        //let hom_indel_rate = LogProb::from(Prob(0.00005));
        //let het_indel_rate = LogProb::from(Prob(0.0001));

        let ts: LogProb = LogProb::from(Prob(ts_tv_ratio / (ts_tv_ratio + 2.0)));
        let tv: LogProb = LogProb::from(Prob(1.0 / (ts_tv_ratio + 2.0)));
        // key of diploid_genotype_priors is (char,(char,char)) (ref_allele, G=(allele1,allele2))
        // key of haploid priors is (char, char) which is (ref_allele, allele1)
        let mut diploid_genotype_priors: HashMap<(char, (char, char)), LogProb> = HashMap::new();
        let mut haploid_genotype_priors: HashMap<(char, char), LogProb> = HashMap::new();

        // define base transitions
        let mut transition: HashMap<char, char> = HashMap::new();
        transition.insert('A', 'G');
        transition.insert('G', 'A');
        transition.insert('T', 'C');
        transition.insert('C', 'T');

        let alleles: Vec<char> = vec!['A', 'C', 'G', 'T'];
        let genotypes: Vec<(char, char)> = vec![
            // all combinations of DNA bases
            ('A', 'A'),
            ('A', 'C'),
            ('A', 'G'),
            ('A', 'T'),
            ('C', 'C'),
            ('C', 'G'),
            ('C', 'T'),
            ('G', 'G'),
            ('G', 'T'),
            ('T', 'T'),
            // all combinations of a DNA base and an Insertion or Deletion
            ('A', 'D'),
            ('A', 'I'),
            ('C', 'D'),
            ('C', 'I'),
            ('G', 'D'),
            ('G', 'I'),
            ('T', 'D'),
            ('T', 'I'),
            // all combinations of multiple Insertion or Deletion
            ('D', 'D'),
            ('I', 'I'),
            ('D', 'I'),
        ];

        for aref in &alleles {
            let allele = *aref;
            // priors on haploid alleles

            haploid_genotype_priors.insert(
                (allele, allele),
                LogProb::ln_one_minus_exp(&(het_snv_rate + het_indel_rate)),
            );
            haploid_genotype_priors.insert(
                (
                    allele,
                    *transition
                        .get(&allele)
                        .chain_err(|| ErrorKind::InvalidTransitionBase(allele.to_string()))?,
                ),
                het_snv_rate + ts,
            );

            for tref in &alleles {
                let transversion = *tref;
                if haploid_genotype_priors.contains_key(&(allele, transversion)) {
                    continue;
                }
                haploid_genotype_priors.insert((allele, transversion), het_snv_rate + tv);
            }

            // assume indel has a 0.5 chance of being insertion, 0.5 chance of deletion
            // scale the prior by 1/4 since the SNP events are divided into 3 types... (transition to each base)
            haploid_genotype_priors.insert(
                (allele, 'D'),
                het_indel_rate + LogProb::from(Prob(0.5)) + LogProb::from(Prob(1.0 / 4.0)),
            );
            haploid_genotype_priors.insert(
                (allele, 'I'),
                het_indel_rate + LogProb::from(Prob(0.5)) + LogProb::from(Prob(1.0 / 4.0)),
            );

            for gt in &genotypes {
                let (g1, g2) = *gt;
                // probability of homozygous reference is the probability of neither het or hom SNP
                if g1 == g2 && g1 == allele {
                    // g1 and g2 are the reference bases
                    let var_rate = LogProb::ln_sum_exp(&[
                        hom_snv_rate,
                        het_snv_rate,
                        hom_indel_rate,
                        het_indel_rate,
                    ]);
                    let one_minus_snp_rate = LogProb::ln_one_minus_exp(&var_rate);
                    diploid_genotype_priors.insert((allele, *gt), one_minus_snp_rate);
                } else if g1 == g2 && g1 != allele {
                    // this could be multiple indels, must handle indel case first
                    if g2 == 'D' || g2 == 'I' {
                        diploid_genotype_priors.insert(
                            (allele, *gt),
                            hom_indel_rate
                                + LogProb::from(Prob(0.5))
                                + LogProb::from(Prob(1.0 / 4.0)),
                        );
                    } else if g1
                        == *transition
                            .get(&allele)
                            .chain_err(|| ErrorKind::InvalidTransitionBase(allele.to_string()))?
                    {
                        // otherwise it is a homozygous SNV
                        // transitions are 4 times as likely as transversions
                        diploid_genotype_priors.insert((allele, *gt), hom_snv_rate + ts);
                    } else {
                        // transversion
                        diploid_genotype_priors.insert((allele, *gt), hom_snv_rate + tv);
                    }
                } else {
                    // else it's the product of the haploid priors
                    diploid_genotype_priors.insert(
                        (allele, *gt),
                        *haploid_genotype_priors.get(&(allele, g1)).chain_err(|| {
                            "Invalid haploid genotype accessed from haploid genotype priors."
                        })? + *(haploid_genotype_priors.get(&(allele, g2))).chain_err(|| {
                            "Invalid haploid genotype accessed from haploid genotype priors."
                        })?,
                    );
                }
            }
        }

        eprintln!("{} GENOTYPE PRIORS:", SPACER);
        eprintln!("{} REF G1/G2 PROB", SPACER);

        for (&(ref ra, (ref g1, ref g2)), &p) in &diploid_genotype_priors {
            eprintln!("{} {} {}/{} {}", SPACER, ra, g1, g2, *Prob::from(p));
        }

        Ok(GenotypePriors {
            priors_dict: diploid_genotype_priors,
        })
    }

    // takes a vector of strings representing alleles (i.e. from Var.alleles), with the 0-th allele being reference
    // and a phased genotype
    // represented as indices into the alleles vector
    // returns the prior probability of that genotype

    // TODO: currently MNPs are going to be assigned the prior of the first SNV in the MNP
    // not ideal, but should suffice for the time being
    pub fn get_prior(&self, alleles: &Vec<String>, genotype: Genotype) -> Result<LogProb> {
        let nth_char: fn(&String, usize) -> Result<char> = |s: &String, i: usize| {
            s.chars()
                .nth(i)
                .chain_err(|| "Error accessing nth character of String")
        };
        //let first_char: fn(&String) -> Result<char> = |s: &String| nth_char(s, 0);
        let mut ra = nth_char(&alleles[0], 0)?;

        let mut g0 = if alleles[genotype.0 as usize].len() == alleles[0].len() {
            nth_char(&alleles[genotype.0 as usize], 0)?
        } else if alleles[genotype.0 as usize].len() > alleles[0].len() {
            'I'
        } else {
            'D'
        };

        if ra == g0 {
            for i in 0..alleles[genotype.0 as usize].len() {
                ra = nth_char(&alleles[0], i)?;
                g0 = nth_char(&alleles[genotype.0 as usize], i)?;
                if ra != g0 {
                    break;
                }
            }
        }

        let mut g1 = if alleles[genotype.1 as usize].len() == alleles[0].len() {
            nth_char(&alleles[genotype.1 as usize], 0)?
        } else if alleles[genotype.1 as usize].len() > alleles[0].len() {
            'I'
        } else {
            'D'
        };

        if ra == g1 {
            for i in 0..alleles[genotype.1 as usize].len() {
                ra = nth_char(&alleles[0], i)?;
                g1 = nth_char(&alleles[genotype.1 as usize], i)?;
                if ra != g1 {
                    break;
                }
            }
        }

        let ln_half = LogProb::from(Prob(0.5));

        match self.priors_dict.get(&(ra, (g0, g1))) {
            Some(p) => {
                if g0 != g1 {
                    return Ok(ln_half + *p);
                } else {
                    return Ok(*p);
                }
            }
            None => match self.priors_dict.get(&(ra, (g1, g0))) {
                Some(p) => {
                    if g0 != g1 {
                        return Ok(ln_half + *p);
                    } else {
                        return Ok(*p);
                    }
                }
                None => bail!(ErrorKind::GenotypeNotInGenotypePriorsError(ra, g0, g1)),
            },
        }
    }

    pub fn get_all_priors(&self, alleles: &Vec<String>) -> Result<GenotypeProbs> {
        let mut priors = GenotypeProbs::zeros(alleles.len());

        for g0 in 0..alleles.len() {
            for g1 in 0..alleles.len() {
                priors.tab[g0][g1] = self
                    .get_prior(alleles, Genotype(g0 as u8, g1 as u8))
                    .chain_err(|| "Error while accessing genotype prior in get_all_priors()")?;
            }
        }

        Ok(priors)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /**********************************************************************************************/
    // TEST GENOTYPE TABLE
    /**********************************************************************************************/
    fn lp(x: f64) -> LogProb {
        LogProb::from(Prob(x))
    }

    fn generate_genotype_probs1() -> GenotypeProbs {
        GenotypeProbs {
            tab: vec![vec![lp(0.001), lp(0.997)], vec![lp(0.001), lp(0.001)]],
        }
    }

    //#[test]
    //fn test_max_prob (){
    //
    //}

}
