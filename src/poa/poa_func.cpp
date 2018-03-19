#include "poa_func.h"
#include "spoa/spoa.hpp"

extern "C" {

    unsigned poa_func(char** seqs, int num_seqs, char* consensus, int consensus_len) {

        if (num_seqs == 0) {
            return (unsigned) 0;
        }

        std::vector<std::string> sequences;

        for (int i = 0; i < num_seqs; i++){
            sequences.push_back((std::string) seqs[i]);
        }

        auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(0), 5, -4, -8);
        auto graph = spoa::createGraph();

        for (const auto& it: sequences) {
            auto alignment = alignment_engine->align_sequence_with_graph(it, graph);
            graph->add_alignment(alignment, it);
        }

        std::string cns = graph->generate_consensus();

        int l = cns.length();
        if (l > consensus_len) {
            l = consensus_len;
        }

        for (int i = 0; i < l; i++){
            consensus[i] = cns[i];
        }

        sequences.clear();
        //delete &sequences;
        //delete &cns;
        //delete &graph;
        //delete &alignment_engine;

        return (unsigned) l;
    }
}
