#include "poa_func.h"
#include "spoa/spoa.hpp"

extern "C" {

    // see the C header file (poa_func.h) for detailed descriptions of each argument
    unsigned poa_func(char** seqs, int num_seqs,
                      char* seed_seq, int num_seeds,
                      char* consensus, int consensus_len,
                      int alignment_type, int match_score, int mismatch_score, int gap_open, int gap_extend) {

        if (num_seqs == 0) {
            return (unsigned) 0;
        }

        // populate the list of sequences
        std::vector<std::string> sequences;
        for (int i = 0; i < num_seqs; i++){
            sequences.push_back((std::string) seqs[i]);
        }

        std::string seed_sequence = (std::string) seed_seq;

        auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(alignment_type),
                                                            (int8_t) match_score,
                                                            (int8_t) mismatch_score,
                                                            (int8_t) gap_open,
                                                            (int8_t) gap_extend);
        auto graph = spoa::createGraph();

        // seed the POA with a seed sequence (e.g. reference sequence)
        // this is useful e.g. to reduce noise from noisy PacBio sequence reads.
        // we set num_seeds to num_seqs / 2, i.e. the reference sequence has half the "weight" or multiplicity of the noisy reads
        // the first alignment is added to initialize nodes in the graph
        // but each subsequent alignment is going to be the same, and is only to increase the seed sequence's weight.
        // so for each subsequent addition we save computation by just adding the same alignment onto the graph over and over.
        if (num_seeds > 0) {
            auto first_seed_aln = (*alignment_engine)(seed_sequence, graph);
            graph->add_alignment(first_seed_aln, seed_sequence);

            auto next_seed_aln = (*alignment_engine)(seed_sequence, graph);
            for (int i = 1; i < num_seeds; i++){
                graph->add_alignment(next_seed_aln, seed_sequence);
            }
        }

        // add each of the real sequences (e.g. noisy sequence reads) to the graph
        for (const auto& it: sequences) {
            auto alignment = (*alignment_engine)(it, graph);
            graph->add_alignment(alignment, it);
        }

        // generate the consensus sequence, assign it to the allocated memory block, and return the consensus length.
        std::string cns = graph->generate_consensus();

        int l = cns.length();
        if (l > consensus_len) {
            l = consensus_len;
        }

        for (int i = 0; i < l; i++){
            consensus[i] = cns[i];
        }

        sequences.clear();

        return (unsigned) l;
    }
}
