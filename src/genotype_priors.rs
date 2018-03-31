
use std::collections::HashMap;
use bio::stats::*;

use util::*;

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
    pub fn new() -> GenotypePriors {

        let hom_snp_rate = LogProb::from(Prob(0.0005));
        let het_snp_rate = LogProb::from(Prob(0.001));
        let hom_indel_rate = LogProb::from(Prob(0.00005));
        let het_indel_rate = LogProb::from(Prob(0.0001));

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
        let genotypes: Vec<(char, char)> = vec![ // all combinations of DNA bases
                                                 ('A','A'),('A','C'),('A','G'),('A','T'),
                                                 ('C', 'C'),('C', 'G'),('C', 'T'),
                                                 ('G', 'G'),('G', 'T'),
                                                 ('T', 'T'),
                                                 // all combinations of a DNA base and an Insertion or Deletion
                                                 ('A','D'),('A','I'),
                                                 ('C','D'),('C','I'),
                                                 ('G','D'),('G','I'),
                                                 ('T','D'),('T','I'),
                                                 // all combinations of multiple Insertion or Deletion
                                                 ('D','D'),('I','I'),('D','I')];


        for aref in &alleles {
            let allele = *aref;
            // priors on haploid alleles
            haploid_genotype_priors.insert((allele, allele), LogProb::ln_one_minus_exp(&(het_snp_rate + het_indel_rate)));
            haploid_genotype_priors.insert((allele, *transition.get(&allele).unwrap()), het_snp_rate + LogProb::from(Prob(4.0 / 6.0)));

            for tref in &alleles {
                let transversion = *tref;
                if haploid_genotype_priors.contains_key(&(allele, transversion)) {
                    continue;
                }
                haploid_genotype_priors.insert((allele, transversion), het_snp_rate + LogProb::from(Prob(1.0 / 6.0)));
            }

            // assume indel has a 0.5 chance of being insertion, 0.5 chance of deletion
            // scale the prior by 1/3 since the SNP events are divided into 3 types... (transition to each base)
            haploid_genotype_priors.insert((allele, 'D'), het_indel_rate + LogProb::from(Prob(0.5)) + LogProb::from(Prob(1.0 / 4.0)));
            haploid_genotype_priors.insert((allele, 'I'), het_indel_rate + LogProb::from(Prob(0.5)) + LogProb::from(Prob(1.0 / 4.0)));

            for gt in &genotypes {
                let (g1, g2) = *gt;
                // probability of homozygous reference is the probability of neither het or hom SNP
                if g1 == g2 && g1 == allele {
                    // g1 and g2 are the reference bases
                    let var_rate = LogProb::ln_sum_exp(&[hom_snp_rate, het_snp_rate, hom_indel_rate, het_indel_rate]);
                    let one_minus_snp_rate = LogProb::ln_one_minus_exp(&var_rate);
                    diploid_genotype_priors.insert((allele, *gt), one_minus_snp_rate);
                } else if g1 == g2 && g1 != allele {
                    // this could be multiple indels, must handle indel case first
                    if g2 == 'D' || g2 == 'I' {
                        diploid_genotype_priors.insert((allele, *gt), hom_indel_rate + LogProb::from(Prob(0.5)) + LogProb::from(Prob(1.0 / 4.0)));
                    } else if g1 == *transition.get(&allele).unwrap() { // otherwise it is a homozygous SNV
                        // transitions are 4 times as likely as transversions
                        diploid_genotype_priors.insert((allele, *gt), hom_snp_rate + LogProb::from(Prob(4.0 / 6.0)));
                    } else {
                        // transversion
                        diploid_genotype_priors.insert((allele, *gt), hom_snp_rate + LogProb::from(Prob(1.0 / 6.0)));
                    }
                } else { // else it's the product of the haploid priors
                    diploid_genotype_priors.insert((allele, *gt), *haploid_genotype_priors.get(&(allele, g1)).unwrap() +
                        *(haploid_genotype_priors.get(&(allele, g2))).unwrap());
                }
            }
        }

        eprintln!("{} GENOTYPE PRIORS:", SPACER);
        eprintln!("{} REF G1/G2 PROB", SPACER);
        //let mut priors_vec: Vec<_> = diploid_genotype_priors.iter()
        //    .map(|(&(ra,(g1,g2)),&p)| (&(ra,(g1,g2)),&*Prob::from(p)))
        //    .collect::<Vec<_>>();
        //priors_vec.sort_by(|a, b| a.partial_cmp(b).unwrap());

        for (&(ref ra, (ref g1, ref g2)), &p) in &diploid_genotype_priors {
            eprintln!("{} {} {}/{} {}", SPACER, ra, g1, g2, *Prob::from(p));
        }

        GenotypePriors{priors_dict: diploid_genotype_priors}
    }

    // takes a vector of strings representing alleles (i.e. from Var.alleles), with the 0-th allele being reference
    // and a phased genotype
    // represented as indices into the alleles vector
    // returns the prior probability of that genotype

    // TODO: currently MNPs are going to be assigned the prior of the first SNV in the MNP
    // not ideal, but should suffice for the time being
    pub fn get(&self, alleles: &Vec<String>, genotype: [u8; 2]) -> LogProb {

        let mut ra = alleles[0].chars().nth(0).unwrap();

        let mut g1 = if alleles[genotype[0] as usize].len() == alleles[0].len() {
            alleles[genotype[0] as usize].chars().nth(0).unwrap()
        } else if alleles[genotype[0] as usize].len() > alleles[0].len() {
            'I'
        } else {
            'D'
        };

        if ra == g1 {
            for i in 0..alleles[genotype[0] as usize].len() {
                ra = alleles[0].chars().nth(i).unwrap();
                g1 = alleles[genotype[0] as usize].chars().nth(i).unwrap();
                if ra != g1 {
                    break;
                }
            }
        }

        let mut g2 = if alleles[genotype[1] as usize].len() == alleles[0].len() {
            alleles[genotype[1] as usize].chars().nth(0).unwrap()
        } else if alleles[genotype[1] as usize].len() > alleles[0].len() {
            'I'
        } else {
            'D'
        };

        if ra == g2 {
            for i in 0..alleles[genotype[1] as usize].len() {
                ra = alleles[0].chars().nth(i).unwrap();
                g2 = alleles[genotype[1] as usize].chars().nth(i).unwrap();
                if ra != g2 {
                    break;
                }
            }
        }

        let ln_half = LogProb::from(Prob(0.5));

        match self.priors_dict.get(&(ra, (g1, g2))) {
            Some(p) => {
                if g1 != g2 {return ln_half+*p;} else {return *p;}
            },
            None => {
                match self.priors_dict.get(&(ra, (g2, g1))) {
                    Some(p) => {
                        if g1 != g2 {return ln_half+*p;} else {return *p;}
                    },
                    None => { println!("{} ({},{})",ra,g2,g1); panic!("Genotype not in genotype priors."); }
                };
            }
        }
    }
}
