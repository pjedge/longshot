/* functions to calculate likelihoods P(read| haplotype) for sequencing errors and chimeric fragments */

#include "common.h"
#include <assert.h>     /* assert */

void calculate_fragscore(struct fragment* Flist, int f, char* h, float* ll) {
    int j = 0, k = 0;
    float p0 = 0, p1 = 0, prob = 0, prob2 = 0;

    for (j = 0; j < Flist[f].blocks; j++) {
        for (k = 0; k < Flist[f].list[j].len; k++) {

            if (h[Flist[f].list[j].offset + k] == '-' || (int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
            prob = (QVoffset - (int) Flist[f].list[j].qv[k]); prob /= 10;
            prob2 = Flist[f].list[j].p1[k];

            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) {
                p0 += prob2;
                p1 += prob;
            } else {
                p0 += prob;
                p1 += prob2;
            }
        }
    }

    // normal LL calculation
    *ll = addlogs(p0,p1);
}

void update_fragscore(struct fragment* Flist, int f, char* h) {
    int j = 0, k = 0;
    float p0 = 0, p1 = 0, prob = 0, prob1 = 0, prob2 = 0;
    Flist[f].calls = 0;
    float good = 0, bad = 0;
    int switches = 0;
    int m = 0;
    if (h[Flist[f].list[0].offset] == Flist[f].list[0].hap[0]) m = 1;
    else m = -1; // initialize
    for (j = 0; j < Flist[f].blocks; j++) {
        Flist[f].calls += Flist[f].list[j].len;
        for (k = 0; k < Flist[f].list[j].len; k++) {
            if (h[Flist[f].list[j].offset + k] == '-') continue; // { fprintf(stdout,"fragment error"); continue;}
            if ((int) Flist[f].list[j].qv[k] - QVoffset < MINQ) continue;
            prob = QVoffset - (int) Flist[f].list[j].qv[k];
            prob /= 10;
            prob1 = 1.0 - pow(10, prob); //prob2 = log10(prob1);
            prob2 = Flist[f].list[j].p1[k];

            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) good += prob1;
            else bad += prob1;
            //if (h[Flist[f].list[j].offset+k] == Flist[f].list[j].hap[k]) good++; else bad++;
            // this is likelihood based calculation
            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k]) {
                p0 += prob2;
                p1 += prob;
            } else {
                p0 += prob;
                p1 += prob2;
            }


            if (h[Flist[f].list[j].offset + k] == Flist[f].list[j].hap[k] && m == -1) {
                m = 1;
                switches++;
            } else if (h[Flist[f].list[j].offset + k] != Flist[f].list[j].hap[k] && m == 1) {
                m = -1;
                switches++;
            }
        }
    }

    // normal LL calculation
    Flist[f].ll = addlogs(p0,p1);
    Flist[f].currscore = -1 * Flist[f].ll;
}

// compute score of fragment
// don't mutate Flist or other structures.
// return score as return value
// homozygous: 0-based index of a homozygous position. -1 if no homozygous pos

// switch_ix1: 0-based index of the first switch error being tested, -1 if none
// switch_ix2: 0-based index of the second switch error being tested, -1 if none

float fragment_ll(struct fragment* Flist, int f, char* h, int homozygous, int switch_ix1, int switch_ix2) {
    int j = 0, k = 0;
    float p0 = 0, p1 = 0, prob = 0, prob2 = 0;
    float ll = 0;
    int snp_ix, switched;

    // normal LL calculation, no Hi-C h-trans
    for (j = 0; j < Flist[f].blocks; j++) {
        for (k = 0; k < Flist[f].list[j].len; k++) {
            snp_ix = Flist[f].list[j].offset + k; // index of current position with respect to all SNPs

            // conditions to skip this base
            if (h[snp_ix] == '-') continue;

            prob = QVoffset - (int) Flist[f].list[j].qv[k];
            prob /= 10;
            prob2 = Flist[f].list[j].p1[k];

            // this is likelihood based calculation
            assert(switch_ix2 == -1);
            switched = (switch_ix1 != -1 && snp_ix >= switch_ix1);
            if ((h[snp_ix] == Flist[f].list[j].hap[k]) != switched) { // true if match, or not match but switched
                p0 += prob2;
                if (snp_ix != homozygous) p1 += prob;
                else p1 += prob2;
            } else {
                p0 += prob;
                if (snp_ix != homozygous) p1 += prob2;
                else p1 += prob;
            }
        }
    }

    ll = addlogs(p0, p1);

    return ll;
}
