//! Wrapper functions for the SPOA (simd-accelerated partial order alignment) library

extern "C" {
    fn poa_func(
        seqs: *const *const u8,
        num_seqs: i32,
        seed_seq: *const u8,
        num_seeds: i32,
        consensus: *const u8,
        consensus_len: i32,
        alignment_type: i32, // 0 = local, 1 = global, 2 = gapped
        match_score: i32,
        mismatch_score: i32,
        gap_score: i32,
    ) -> u32;
}

pub fn poa_multiple_sequence_alignment(
    seqs: &Vec<Vec<u8>>,
    seed_seq: Vec<u8>,
    num_seeds: i32,
    consensus: &mut Vec<u8>,
    alignment_type: i32,
    match_score: i32,
    mismatch_score: i32,
    gap_score: i32,
) {
    unsafe {
        let num_seqs = seqs.len() as i32;
        let consensus_len = consensus.len() as i32;

        let mut seq_ptrs: Vec<*const u8> = Vec::with_capacity(seqs.len());

        for seq in seqs {
            seq_ptrs.push(seq.as_ptr());
        }

        let len = poa_func(
            seq_ptrs.as_ptr(),
            num_seqs,
            seed_seq.as_ptr(),
            num_seeds,
            consensus.as_ptr(),
            consensus_len,
            alignment_type,
            match_score,
            mismatch_score,
            gap_score,
        );

        consensus.truncate(len as usize);
    }
}
