extern crate rust_htslib;

use rust_htslib::bam;
use rust_htslib::prelude::*;
use bio::io::fasta;
use std::char;

static VARLIST_CAPACITY: usize = 1000000;

pub struct PotentialVar {
    chrom: String,
    pos0: usize,
    ref_allele: String,
    var_allele: String,
}

pub fn parse_target_names(bam_file: String) -> Vec<String> {
    let bam = bam::Reader::from_path(bam_file).unwrap();
    let header_view = bam.header();
    let target_names_dec: Vec<&[u8]> = header_view.target_names();
    let mut target_names: Vec<String> = vec![];

    for t_name_dec in target_names_dec {
        let mut name_vec: Vec<char> = vec![];
        for decr in t_name_dec {
            let dec: u8 = *decr;
            name_vec.push(dec as char);
        }
        let name_string: String = name_vec.into_iter().collect();
        target_names.push(name_string);
    }

    target_names
}

pub fn call_potential_snvs(bam_file: String, fasta_file: String) -> Vec<PotentialVar> {

    let target_names = parse_target_names(bam_file.clone());

    let bam = bam::Reader::from_path(bam_file).unwrap();
    let mut fasta = fasta::IndexedReader::from_file(&fasta_file).unwrap();

    let mut varlist: Vec<PotentialVar> = Vec::with_capacity(VARLIST_CAPACITY);

    // pileup over all covered sites
    let bases = ['A', 'C', 'G', 'T', 'N'];
    let mut ref_seq: Vec<u8> = vec![];
    let mut prev_tid = 4294967295;

    for p in bam.pileup() {
        let pileup = p.unwrap();
        let tid: usize = pileup.tid() as usize;

        if tid != prev_tid {
            fasta.read_all(&target_names[tid], &mut ref_seq).expect("Failed to read fasta sequence record.");
        }
        //println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());

        let mut counts = [0; 5]; // A,C,G,T,N

        // pileup the bases for a single position and count number of each base
        for alignment in pileup.alignments() {
            if !alignment.is_del() && !alignment.is_refskip() {

                let base: char = alignment.record().seq()[alignment.qpos().unwrap()] as char;
                let b = match base {
                    'A' | 'a' => 0,
                    'C' | 'c' => 1,
                    'G' | 'g' => 2,
                    'T' | 't' => 3,
                    'N' | 'n' => 4,
                    _ => panic!("Invalid base read from BAM file."),
                };

                counts[b] += 1;
            }
        }


        //println!("counts[A]: {}", counts[0]);
        //println!("counts[C]: {}", counts[1]);
        //println!("counts[G]: {}", counts[2]);
        //println!("counts[T]: {}", counts[3]);
        //println!("counts[N]: {}", counts[4]);
        //println!("---------------------\n");

        let mut var_allele = "N".to_string();
        let mut max_count = -1;
        let mut max_base = 'N';

        for i in 0..5 {
            if counts[i] > max_count {
                max_count = counts[i];
                max_base = bases[i];
            }
        }
        //println!("{}", max_base);

        let ref_allele: String =
            (ref_seq[pileup.pos() as usize] as char).to_string().to_uppercase();
        //let pos0: usize = pileup.pos() as usize;

        if max_count > 3 && !(max_base.to_string() == ref_allele) {
            var_allele = max_base.to_string();

            println!("{}\t{}\t{}\t{}",
                     target_names[tid].clone(),
                     pileup.pos(),
                     ref_allele,
                     var_allele);

            let tid: usize = pileup.tid() as usize;
            let new_var = PotentialVar {
                chrom: target_names[tid].clone(),
                pos0: pileup.pos() as usize,
                ref_allele: ref_allele,
                var_allele: var_allele,
            };

            varlist.push(new_var);
        }

        prev_tid = tid;
    }
    varlist
}
