// calls HapCUT2 as a static library
use util::Fragment;
use bio::stats::{LogProb, Prob, PHREDProb};

static MAX_QUAL_F64: f64 = 0.2;

pub fn generate_flist_buffer(flist: &Vec<Fragment>, phase_variant: &Vec<bool>) -> Vec<Vec<u8>> {
    let mut buffer: Vec<Vec<u8>> = vec![];
    for frag in flist {
        let prev_call = phase_variant.len() + 1;
        let mut quals: Vec<u8> = vec![];
        let mut blocks = 0;
        let mut n_calls = 0;

        for c in frag.clone().calls {
            if phase_variant[c.var_ix] && c.qual < LogProb::from(Prob(MAX_QUAL_F64)) {
                n_calls += 1;
                if c.var_ix - prev_call != 1 {
                    blocks += 1;
                }
            }
        }
        if n_calls < 2 {
            continue;
        }

        let mut line: Vec<u8> = vec![];
        for u in blocks.to_string().into_bytes() {
            line.push(u as u8);
        }
        line.push(' ' as u8);

        for u in frag.id.clone().into_bytes() {
            line.push(u as u8);
        }
        line.push(' ' as u8);

        for c in frag.clone().calls {
            if phase_variant[c.var_ix] {
                if c.var_ix - prev_call == 1 {
                    line.push(c.allele as u8)
                } else {
                    line.push(' ' as u8);
                    for u in (c.var_ix + 1).to_string().into_bytes() {
                        line.push(u as u8);
                    }
                    line.push(' ' as u8);
                    line.push(c.allele as u8)
                }
                let mut qint = *PHREDProb::from(c.qual) as u32 + 33;
                if qint > 126 {
                    qint = 126;
                }
                quals.push(qint as u8);
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

        //print!("{}", charline.iter().collect::<String>());

        buffer.push(line);
    }
    buffer
}

extern "C" {
    fn hapcut2(fragmentbuffer: *const *const u8,
               variantbuffer: *const *const u8,
               fragments: usize,
               snps: usize,
               hap1: *mut u8,
               phase_sets: *mut i32,);
}

pub fn call_hapcut2(frag_buffer: &Vec<Vec<u8>>,
                    vcf_buffer: &Vec<Vec<u8>>,
                    fragments: usize,
                    snps: usize,
                    hap1: &mut Vec<u8>,
                    phase_sets: &mut Vec<i32>) {

    unsafe {
        let mut frag_ptrs: Vec<*const u8> = Vec::with_capacity(frag_buffer.len());
        let mut vcf_ptrs: Vec<*const u8> = Vec::with_capacity(vcf_buffer.len());

        for line in frag_buffer {
            frag_ptrs.push(line.as_ptr());
        }

        for line in vcf_buffer {
            vcf_ptrs.push(line.as_ptr());
        }

        hapcut2(frag_ptrs.as_ptr(),
                vcf_ptrs.as_ptr(),
                fragments,
                snps,
                hap1.as_mut_ptr(),
                phase_sets.as_mut_ptr());
    }
}
