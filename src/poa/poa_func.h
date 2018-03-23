#ifndef POA_FUNC_H
#define POA_FUNC_H

#ifdef __cplusplus
extern "C" {
#endif
unsigned poa_func(char** seqs,        // the sequences (null-terminated) to perform multiple-sequence-alignment with.
                  int num_seqs,       // the number of sequences being multiply aligned
                  char* seed_seq,     // a seed sequence to start the multiple alignment. for example, a reference sequence.
                  int num_seeds,      // the weight or multiplicity (in number of sequences) to give the seed sequence.
                                      // a higher number biases the multiple alignment toward the seed more!
                  char* consensus,    // this chunk of memory will hold the return value, the consensus of the multiple alignment
                  int consensus_len,  // the amount of memory allocated to consensus, i.e. use the MAXIMUM length of the consensus
                  int alignment_type, // the alignment type: 0 = local align, 1 = global align, 2 = semi-global
                  int match_score,    // the score to give a sequence match in alignment, e.g. 5
                  int mismatch_score, // the score to give a sequence mismatch in alignment, e.g. -4
                  int gap_score);     // the score to give a sequence gap in alignment, e.g. -8

#ifdef __cplusplus
}
#endif

#endif // POA_FUNC_H
