
/*
CODE for likelihood based cut calculation that preserves long reads and operates on fragments rather than edges

code first implemented 03/02/15 | likelihood based model for initializing two shores of the cut (startnode... secondnode)

##for each read store 4 floats  P(R | H), P(R | complement(H)) and P(R | H_new), P(R | complement(H_new)) where likelihood is only for nodes added to the cut till now
                                scores[0] scores[1]                scores[2]     scores[3]

where H_new is new haplotype formed by flipping the phase of vertices in shore2 relative to shore1

score of a variant snpfrag[node].score = log10( P(R|H) + P(R | complement(H)) - log10( P(R|H_new) + P(R | complement(H_new) )
 */
#include "common.h"
#include "math.h"
// simple calculation for difference between likelihood of old and new solution based on the Flist[f].scores used for building max-cut

float cut_score(struct fragment* Flist, struct BLOCK* component, char* hap) {
    int i = 0, t = 0, f = 0;
    float oldLL = 0, newLL = 0;
    float scores[4];   // normal scores
    float Ln, fL;
    // use the 4 values fragment.scores[] to calculate the likelihoods
    for (i = 0; i < component->frags; i++) {
        f = component->flist[i];
        for (t = 0; t < 4; t++) scores[t] = Flist[f].scores[t];

        Ln = addlogs(scores[2], scores[3]);

        calculate_fragscore(Flist, f, hap, &fL);
        oldLL += fL; // ready to go with h-trans calculations
        newLL += Ln;

    }
    if (newLL > oldLL && DEBUG) fprintf_time(stderr, "component %d old %f new %f %f \n", component->phased, oldLL, newLL, component->SCORE);
    if (DEBUG) fprintf(stdout, "FRAG component %d old %f new %f %f \n\n", component->phased, oldLL, newLL, component->SCORE);
    return oldLL - newLL;
}

// initialize the scores for initial two nodes
void init_fragment_scores(struct SNPfrags* snpfrag, struct fragment* Flist, char* hap, int startnode, int secondnode) {
    int i=0, j = 0, f = 0, t = 0, n = 0, k = 0, node = 0, node1;
    float prob, prob2;   // prob = P(miscall). prob2 == P(correct)
    float scores[4];     // a temp variable
    float Lo, Ln;        // old likelihood, new likelihood

    for (i = 0; i < 2; i++) {

        node1 = i ? secondnode:startnode;
        if (i && node1 < 0) break; // if secondnode is not valid, we should exit, otherwise segfault

        for (n = 0; n < snpfrag[node1].frags; n++) {
            f = snpfrag[node1].flist[n]; // index into Flist global
            //if (Flist[f].init == '0') continue;
            j = snpfrag[node1].jlist[n];
            k = snpfrag[node1].klist[n];

            node = node1;
            if (hap[node] == '-' || (int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;

            prob = (QVoffset - (int) Flist[f].list[j].qv[k]); prob /= 10;
            prob2 = Flist[f].list[j].p1[k];

            if (node == startnode) {
                if (hap[node] == Flist[f].list[j].hap[k]){
                    Flist[f].scores[0] += prob2;
                    Flist[f].scores[1] += prob;
                    Flist[f].scores[2] += prob2;
                    Flist[f].scores[3] += prob;
                } else {
                    Flist[f].scores[0] += prob;
                    Flist[f].scores[1] += prob2;
                    Flist[f].scores[2] += prob;
                    Flist[f].scores[3] += prob2;
                }
            } else if (node == secondnode){ // probabilities are flipped compared to start node

                if (hap[node] == Flist[f].list[j].hap[k]) {
                    Flist[f].scores[0] += prob;
                    Flist[f].scores[1] += prob2;
                    Flist[f].scores[2] += prob;
                    Flist[f].scores[3] += prob2;
                } else {
                    Flist[f].scores[0] += prob2;
                    Flist[f].scores[1] += prob;
                    Flist[f].scores[2] += prob2;
                    Flist[f].scores[3] += prob;
                }
            }

            // update score of every node outside 2 shores covered by 'f'
            for (j = 0; j < Flist[f].blocks; j++){

                for (k = 0; k < Flist[f].list[j].len; k++) {
                    node = Flist[f].list[j].offset + k;
                    if (hap[node] == '-' || (int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
                    if (node == startnode || node == secondnode) continue;

                    prob = (QVoffset - (int) Flist[f].list[j].qv[k]); prob /= 10;
                    prob2 = Flist[f].list[j].p1[k];

                    for (t = 0; t < 4; t++) scores[t] = Flist[f].scores[t];
                    if (hap[node] == Flist[f].list[j].hap[k]) {
                        scores[0] += prob2;
                        scores[1] += prob; // add node to startnode side of cut
                        scores[2] += prob;
                        scores[3] += prob2; // add to other side of cut
                    } else {
                        scores[0] += prob;
                        scores[1] += prob2; // allele mismatch so flip probability
                        scores[2] += prob2;
                        scores[3] += prob;
                    }

                    // Lo is likelihood of fragment if SNP is added to 'startnode' and Ln if SNP is added to 'secondnode'
                    Lo = addlogs(scores[0], scores[1]);
                    Ln = addlogs(scores[2], scores[3]);

                    snpfrag[node].score += Lo - Ln;

                }
            }
        }
    }
}

// updates fragment scores and variant scores as a result of adding the node 'node_added' to the growing max-cut

void update_fragment_scores(struct SNPfrags* snpfrag, struct fragment* Flist, char* hap, int startnode, int secondnode, int node_added, struct PHEAP* pheap, int* slist) {
    int j = 0, f = 0, t = 0, n = 0, k = 0, node = 0;
    float prob, prob2;
    float scores[4];
    float Lo, Ln;
    float f_scores[4];
    float oldscore;

    for (n = 0; n < snpfrag[node_added].frags; n++){
        f = snpfrag[node_added].flist[n]; // index into Flist global
        for (t = 0; t < 4; t++) f_scores[t] = Flist[f].scores[t]; // store previous fragment scores before updating

        j = snpfrag[node_added].jlist[n];
        k = snpfrag[node_added].klist[n];

        if (hap[Flist[f].list[j].offset + k] == '-' || (int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
        node = Flist[f].list[j].offset + k;

        if (node != node_added) continue;

        prob = QVoffset - (int) Flist[f].list[j].qv[k];
        prob /= 10; // log10(e)
        prob2 = Flist[f].list[j].p1[k];

        if (snpfrag[node].parent == startnode){ // if node is added to 'startnode', original and new likelihoods are updated identically

            if (hap[node] == Flist[f].list[j].hap[k]){
                Flist[f].scores[0] += prob2;
                Flist[f].scores[1] += prob;
                Flist[f].scores[2] += prob2;
                Flist[f].scores[3] += prob;
            } else {
                Flist[f].scores[0] += prob;
                Flist[f].scores[1] += prob2;
                Flist[f].scores[2] += prob;
                Flist[f].scores[3] += prob2;
            }
        } else if (snpfrag[node].parent == secondnode){
            if (hap[node] == Flist[f].list[j].hap[k]){
                Flist[f].scores[0] += prob;
                Flist[f].scores[1] += prob2;
                Flist[f].scores[2] += prob;
                Flist[f].scores[3] += prob2;
            } else {
                Flist[f].scores[0] += prob2;
                Flist[f].scores[1] += prob;
                Flist[f].scores[2] += prob2;
                Flist[f].scores[3] += prob;
            }
        }

        for (j = 0; j < Flist[f].blocks; j++){ // update score of every node outside 2 shores covered by 'f'

            for (k = 0; k < Flist[f].list[j].len; k++) {
                if (hap[Flist[f].list[j].offset + k] == '-' || (int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
                node = Flist[f].list[j].offset + k;
                if (snpfrag[node].parent != startnode && snpfrag[node].parent != secondnode && node != node_added) {
                    oldscore = snpfrag[node].score; // store old score
                    prob = QVoffset - (int) Flist[f].list[j].qv[k];
                    prob /= 10; // log10(e)
                    //prob1 = 1.0 - pow(10,prob); prob2 = log10(prob1);
                    prob2 = Flist[f].list[j].p1[k];
                    for (t = 0; t < 4; t++) scores[t] = f_scores[t];
                    if (hap[node] == Flist[f].list[j].hap[k]) {
                        scores[0] += prob2;
                        scores[1] += prob;
                        scores[2] += prob;
                        scores[3] += prob2;
                    } else {
                        scores[0] += prob;
                        scores[1] += prob2;
                        scores[2] += prob2;
                        scores[3] += prob;
                    }

                    Lo = addlogs(scores[0], scores[1]);
                    Ln = addlogs(scores[2], scores[3]);

                    snpfrag[node].score -= Lo - Ln;; // subtract old score for variant

                    for (t = 0; t < 4; t++) scores[t] = Flist[f].scores[t];
                    if (hap[node] == Flist[f].list[j].hap[k]) {
                        scores[0] += prob2;
                        scores[1] += prob;
                        scores[2] += prob;
                        scores[3] += prob2;
                    } else {
                        scores[0] += prob;
                        scores[1] += prob2;
                        scores[2] += prob2;
                        scores[3] += prob;
                    }

                    Lo = addlogs(scores[0], scores[1]);
                    Ln = addlogs(scores[2], scores[3]);

                    snpfrag[node].score += Lo - Ln; // add new delta LL for variant

                    if (fabsf(oldscore) > fabsf(snpfrag[node].score)){ // score decreased
                        pmaxHeapify(pheap, snpfrag[node].heaploc, snpfrag, slist);

                    } else pbubbleUp(pheap, snpfrag[node].heaploc, snpfrag, slist);
                }
            }
        }
    }
}
