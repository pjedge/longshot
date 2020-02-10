//! Contains functions related to haplotype assembly, including a FFI wrapper function for HapCUT2,
//! MEC criteria, haplotype read separation, etc.
use bio::stats::{LogProb, PHREDProb, Prob};
use errors::*;
use hashbrown::{HashMap, HashSet};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::char::from_digit;
use util::*;
use variants_and_fragments::*;

pub fn separate_fragments_by_haplotype(
    flist: &Vec<Fragment>,
    threshold: LogProb,
) -> Result<(HashSet<String>, HashSet<String>)> {
    let mut h1 = HashSet::new();
    let mut h2 = HashSet::new();

    let mut h1_count = 0;
    let mut h2_count = 0;
    let mut unassigned_count = 0;

    for ref f in flist {

        // we store p_read_hap as ln-scaled f16s to save space. need to convert back.
        let p_read_hap0 = LogProb(f64::from(f.p_read_hap[0]));
        let p_read_hap1 = LogProb(f64::from(f.p_read_hap[1]));

        let total: LogProb = LogProb::ln_add_exp(p_read_hap0, p_read_hap1);
        let p_read_hap0: LogProb = p_read_hap0 - total;
        let p_read_hap1: LogProb = p_read_hap1 - total;

        match &f.id {
            &Some(ref fid) => {
                if p_read_hap0 > threshold {
                    h1_count += 1;
                    h1.insert(fid.clone());
                } else if p_read_hap1 > threshold {
                    h2_count += 1;
                    h2.insert(fid.clone());
                } else {
                    unassigned_count += 1;
                }
            }
            &None => {bail!("Fragment without read ID found while separating reads by haplotype.");}
        }
    }

    // count the number assigned to either haplotype
    let total: f64 = (h1_count + h2_count + unassigned_count) as f64;
    let h1_percent: f64 = 100.0 * h1_count as f64 / total;
    let h2_percent: f64 = 100.0 * h2_count as f64 / total;
    let unassigned_percent: f64 = 100.0 * unassigned_count as f64 / total;

    eprintln!(
        "{}     {} reads ({:.2}%) assigned to haplotype 1",
        print_time(),
        h1_count,
        h1_percent
    );
    eprintln!(
        "{}     {} reads ({:.2}%) assigned to haplotype 2",
        print_time(),
        h2_count,
        h2_percent
    );
    eprintln!(
        "{}     {} reads ({:.2}%) unassigned.",
        print_time(),
        unassigned_count,
        unassigned_percent
    );

    Ok((h1, h2))
}

pub fn separate_bam_reads_by_haplotype(
    bamfile_name: &String,
    interval: &Option<GenomicInterval>,
    hap_bam_prefix: String,
    h1: &HashSet<String>,
    h2: &HashSet<String>,
    min_mapq: u8,
) -> Result<()> {
    let interval_lst: Vec<GenomicInterval> = get_interval_lst(bamfile_name, interval)
        .chain_err(|| "Error getting genomic interval list.")?;

    let mut bam_ix =
        bam::IndexedReader::from_path(bamfile_name).chain_err(|| ErrorKind::IndexedBamOpenError)?;

    let h1_bam_file = format!("{}.hap1.bam", &hap_bam_prefix);
    let h2_bam_file = format!("{}.hap2.bam", &hap_bam_prefix);
    let unassigned_bam_file = format!("{}.unassigned.bam", &hap_bam_prefix);

    let header = bam::Header::from_template(&bam_ix.header());
    let mut h1_bam = bam::Writer::from_path(&h1_bam_file, &header,bam::Format::BAM)
        .chain_err(|| ErrorKind::BamWriterOpenError(h1_bam_file))?;
    let mut h2_bam = bam::Writer::from_path(&h2_bam_file, &header,bam::Format::BAM)
        .chain_err(|| ErrorKind::BamWriterOpenError(h2_bam_file))?;
    let mut unassigned_bam = bam::Writer::from_path(&unassigned_bam_file, &header,bam::Format::BAM)
        .chain_err(|| ErrorKind::BamWriterOpenError(unassigned_bam_file))?;

    for iv in interval_lst {
        bam_ix
            .fetch(iv.tid, iv.start_pos, iv.end_pos + 1)
            .chain_err(|| ErrorKind::IndexedBamFetchError)?;

        for r in bam_ix.records() {
            let record = r.chain_err(|| ErrorKind::IndexedBamRecordReadError)?;

            if record.is_quality_check_failed()
                || record.is_duplicate()
                || record.is_secondary()
                || record.is_unmapped()
                || record.mapq() < min_mapq
                || record.is_supplementary()
            {
                continue;
            }

            let qname = u8_to_string(record.qname())?;

            if h1.contains(&qname) {
                h1_bam
                    .write(&record)
                    .chain_err(|| ErrorKind::BamRecordWriteError(qname))?;
            } else if h2.contains(&qname) {
                h2_bam
                    .write(&record)
                    .chain_err(|| ErrorKind::BamRecordWriteError(qname))?;
            } else {
                unassigned_bam
                    .write(&record)
                    .chain_err(|| ErrorKind::BamRecordWriteError(qname))?;
            }
        }
    }

    Ok(())
}

pub fn generate_flist_buffer(
    flist: &Vec<Fragment>,
    phase_variant: &Vec<bool>,
    max_p_miscall: f64,
    single_reads: bool,
) -> Result<Vec<Vec<u8>>> {
    let mut buffer: Vec<Vec<u8>> = vec![];
    let mut frag_num = 0;
    for frag in flist {
        let mut prev_call = phase_variant.len() + 1;
        let mut quals: Vec<u8> = vec![];
        let mut blocks: usize = 0;
        let mut n_calls: usize = 0;

        for c in frag.clone().calls {
            if phase_variant[c.var_ix as usize] && c.qual < LogProb::from(Prob(max_p_miscall)) {
                n_calls += 1;
                if prev_call > phase_variant.len() || c.var_ix as usize - prev_call != 1 {
                    blocks += 1;
                }
                prev_call = c.var_ix as usize;
            }
        }
        if !single_reads && n_calls == 1 {
            continue;
        }
        if n_calls == 0 {
            continue;
        }

        let mut line: Vec<u8> = vec![];
        for u in blocks.to_string().into_bytes() {
            line.push(u as u8);
        }
        line.push(' ' as u8);

        let fid = match &frag.id {
            Some(ref fid) => {fid.clone()}
            None => {frag_num.to_string()}
        };

        for u in fid.clone().into_bytes() {
            line.push(u as u8);
        }
        //line.push(' ' as u8);

        let mut prev_call = phase_variant.len() + 1;

        for c in frag.clone().calls {
            if phase_variant[c.var_ix as usize] && c.qual < LogProb::from(Prob(max_p_miscall)) {
                if prev_call < c.var_ix && c.var_ix - prev_call == 1 {
                    ensure!(
                        c.allele == 0 as u8 || c.allele == 1 as u8,
                        "Allele is not valid for incorporation into fragment file."
                    );
                    line.push(
                        from_digit(c.allele as u32, 10)
                            .chain_err(|| "Error converting allele digit to char.")?
                            as u8,
                    )
                } else {
                    line.push(' ' as u8);
                    for u in (c.var_ix + 1).to_string().into_bytes() {
                        line.push(u as u8);
                    }
                    line.push(' ' as u8);
                    ensure!(
                        c.allele == 0 as u8 || c.allele == 1 as u8,
                        "Allele is not valid for incorporation into fragment file."
                    );
                    line.push(
                        from_digit(c.allele as u32, 10)
                            .chain_err(|| "Error converting allele digit to char.")?
                            as u8,
                    )
                }
                let mut qint = *PHREDProb::from(c.qual) as u32 + 33;
                if qint > 126 {
                    qint = 126;
                }
                quals.push(qint as u8);
                prev_call = c.var_ix
            }
        }
        line.push(' ' as u8);
        line.append(&mut quals);
        //line.push('\n' as u8);
        line.push('\0' as u8);

        let mut charline: Vec<char> = vec![];
        for u in line.clone() {
            charline.push(u as char)
        }

        //println!("{}", charline.iter().collect::<String>());

        buffer.push(line);
        frag_num += 1;
    }
    Ok(buffer)
}

extern "C" {
    fn hapcut2(
        fragmentbuffer: *const *const u8,
        fragments: usize,
        snps: usize,
        hap1: *mut u8,
        phase_sets: *mut i32,
    );
}


pub fn call_hapcut2(
    frag_buffer: &Vec<Vec<u8>>,
    fragments: usize,
    snps: usize,
    hap1: &mut Vec<u8>,
    phase_sets: &mut Vec<i32>,
) {
    unsafe {
        let mut frag_ptrs: Vec<*const u8> = Vec::with_capacity(frag_buffer.len());

        for line in frag_buffer {
            frag_ptrs.push(line.as_ptr());
        }

        hapcut2(
            frag_ptrs.as_ptr(),
            fragments,
            snps,
            hap1.as_mut_ptr(),
            phase_sets.as_mut_ptr(),
        );
    }
}

pub fn calculate_mec(
    flist: &Vec<Fragment>,
    varlist: &mut VarList,
    max_p_miscall: f64,
) -> Result<()> {
    let hap_ixs = vec![0, 1];
    let ln_max_p_miscall = LogProb::from(Prob(max_p_miscall));

    for mut var in &mut varlist.lst {
        var.mec = 0;
        var.mec_frac_variant = 0.0;
        var.mec_frac_block = 0.0;
    }

    for f in 0..flist.len() {
        let mut mismatched_vars: Vec<Vec<usize>> = vec![vec![], vec![]];

        for &hap_ix in &hap_ixs {
            for call in &flist[f].calls {
                if call.qual < ln_max_p_miscall {
                    let var = &varlist.lst[call.var_ix as usize]; // the variant that the fragment call covers

                    if var.phase_set == None {
                        continue; // only care about phased variants.
                    }

                    let hap_allele = if hap_ix == 0 {
                        var.genotype.0
                    } else {
                        var.genotype.1
                    };

                    // read allele matches haplotype allele
                    if call.allele != hap_allele {
                        mismatched_vars[hap_ix].push(call.var_ix);
                    }
                }
            }
        }

        let min_error_hap = if mismatched_vars[0].len() < mismatched_vars[1].len() {
            0
        } else {
            1
        };

        for &ix in &mismatched_vars[min_error_hap] {
            varlist.lst[ix as usize].mec += 1;
        }
    }

    let mut block_mec: HashMap<usize, usize> = HashMap::new();
    let mut block_total: HashMap<usize, usize> = HashMap::new();

    for var in &mut varlist.lst {
        match var.phase_set {
            Some(ps) => {
                *block_mec.entry(ps).or_insert(0) += var.mec;
                *block_total.entry(ps).or_insert(0) += var.allele_counts.iter().sum::<u16>() as usize;
            }
            None => {}
        }
    }

    for mut var in &mut varlist.lst {
        match var.phase_set {
            Some(ps) => {
                var.mec_frac_block = *block_mec
                        .get(&ps)
                        .chain_err(|| "Error retrieving MEC for phase set.")?
                        as f64
                        / *block_total
                            .get(&ps)
                            .chain_err(|| "Error retrieving MEC for phase set.")?
                            as f64;
                var.mec_frac_variant =
                    var.mec as f64 / var.allele_counts.iter().sum::<u16>() as f64;
            }
            None => {}
        }
    }

    Ok(())
}
