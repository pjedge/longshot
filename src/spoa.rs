
extern "C" {
    fn poa_func(seqs: *const *const u8,
               num_seqs: usize,
               consensus: *const u8,
               consensus_len: usize) -> u32;
}

pub fn poa_multiple_sequence_alignment(seqs: &Vec<Vec<u8>>, consensus: &mut Vec<u8>) {

    unsafe {

        let num_seqs: usize = seqs.len() as usize;
        let consensus_len: usize = consensus.len() as usize;

        let mut seq_ptrs: Vec<*const u8> = Vec::with_capacity(seqs.len());

        for seq in seqs {
            seq_ptrs.push(seq.as_ptr());
        }

        let len = poa_func(seq_ptrs.as_ptr(),
                num_seqs,
                consensus.as_ptr(),
                consensus_len);

        consensus.truncate(len as usize);
    }
}
