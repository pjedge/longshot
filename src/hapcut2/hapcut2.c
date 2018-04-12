// CODE STARTED SEPT 10 2007 4pm //  april 8 2008 this code used for producing results in ECCB 2008 paper //
// author: VIKAS BANSAL (vbansal@scripps.edu) last modified December 23, 2010

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include "common.h"
#include "fragmatrix.h"
#include "pointerheap.h"
#include "readinputbuffers.h"

// Printing related
int VERBOSE = 0;
int PRINT_FRAGMENT_SCORES = 0; // output the MEC/switch error score of erroneous reads/fragments to a file for evaluation

// Quality-score related parameters
int QVoffset = 33;
int MINQ = 6; // additional base quality filter in hapcut added april 18 2012

// Number of iterations
int MAXITER = 10000;     // maximum number of global iterations
int MAXCUT_ITER = 10000; // maximum number of iterations for max-cut algorithm, if this is proportional to 'N' -> complexity is 'N^2', added march 13 2013
int CONVERGE = 5; // stop iterations on a given block if exceed this many iterations since improvement

// Post-processing related variables
float THRESHOLD = 0.8;
float SPLIT_THRESHOLD = 0.8;
float HOMOZYGOUS_PRIOR = -80; // in log form. assumed to be really unlikely
int CALL_HOMOZYGOUS = 0;
int SPLIT_BLOCKS = 0;
int DISCRETE_PRUNING = 0;
int ERROR_ANALYSIS_MODE = 0;
int SKIP_PRUNE = 0;

int AUTODETECT_LONGREADS = 1;
int LONG_READS = 0; // if this variable is 1, the data contains long read data

#include "find_maxcut.c"   // function compute_good_cut
#include "post_processing.c"  // post-processing functions

int hapcut2(char** fragmentbuffer, char** variantbuffer, int fragments, int snps, char* HAP1, int* PS) {
    // IMP NOTE: all SNPs start from 1 instead of 0 and all offsets are 1+

    srand(1);
    if (VERBOSE) fprintf_time(stderr, "Calling Max-Likelihood-Cut based haplotype assembly algorithm\n");

    int iter = 0, components = 0;
    int i = 0, k = 0;
    int* slist;
    int flag = 0;
    float bestscore = 0, miscalls = 0;
    struct SNPfrags* snpfrag = NULL;
    struct BLOCK* clist;
    int converged_count=0, component;

    // READ FRAGMENT MATRIX
    struct fragment* Flist;
    Flist     = (struct fragment*) malloc(sizeof (struct fragment)* fragments);
    flag = read_fragment_buffer(fragmentbuffer, Flist, fragments);

    if (flag < 0) {
        fprintf_time(stderr, "unable to read fragment info buffer\n");
        return -1;
    }

    float mean_snps_per_read = 0;
    if (AUTODETECT_LONGREADS){
        for (i = 0; i < fragments; i++){
            mean_snps_per_read += Flist[i].calls;
        }
        mean_snps_per_read /= fragments;
        if (mean_snps_per_read >= 3){
            LONG_READS = 1;
        }else{
            LONG_READS = 0;
        }
    }

    snpfrag = (struct SNPfrags*) malloc(sizeof (struct SNPfrags)*snps);
    update_snpfrags(Flist, fragments, snpfrag, snps);

    // 10/25/2014, edges are only added between adjacent nodes in each fragment and used for determining connected components...
    for (i = 0; i < snps; i++) snpfrag[i].elist = (struct edge*) malloc(sizeof (struct edge)*(snpfrag[i].edges+1));
    if (LONG_READS ==0){
        add_edges(Flist,fragments,snpfrag,snps,&components);
    }else if (LONG_READS >=1){
        add_edges_fosmids(Flist,fragments,snpfrag,snps,&components);
    }

    for (i = 0; i < snps; i++) snpfrag[i].telist = (struct edge*) malloc(sizeof (struct edge)*(snpfrag[i].edges+1));

    // this considers only components with at least two nodes
    if (VERBOSE) fprintf_time(stderr, "fragments %d snps %d component(blocks) %d\n", fragments, snps, components);

    // BUILD COMPONENT LIST
    clist = (struct BLOCK*) malloc(sizeof (struct BLOCK)*components);
    generate_clist_structure(Flist, fragments, snpfrag, snps, components, clist);

    // READ VCF FILE
    read_vcf_buffer(variantbuffer, snpfrag, snps);

    // INITIALIZE RANDOM HAPLOTYPES
    //char* HAP1 = (char*)malloc(snps+1);

    //for (i = 0; i < snps; i++) {
    //    if (snpfrag[i].frags == 0) {
    //        HAP1[i] = '-';
    //    } else if (drand48() < 0.5) {
    //        HAP1[i] = '0';
    //    } else {
    //        HAP1[i] = '1';
    //    }
    //}

    // for each block, we maintain best haplotype solution under MFR criterion
    // compute the component-wise score for 'initHAP' haplotype
    miscalls = 0;
    bestscore = 0;
    for (k = 0; k < components; k++) {
        clist[k].SCORE = 0;
        clist[k].bestSCORE = 0;
        for (i = 0; i < clist[k].frags; i++) {
            update_fragscore(Flist, clist[k].flist[i], HAP1);
            clist[k].SCORE += Flist[clist[k].flist[i]].currscore;
        }
        clist[k].bestSCORE = clist[k].SCORE;
        bestscore += clist[k].bestSCORE;
        miscalls += clist[k].SCORE;
    }

    if (VERBOSE)fprintf_time(stderr, "processed fragment file and variant file: fragments %d variants %d\n", fragments, snps);

    slist = (int*) malloc(sizeof (int)*snps);

    for (k = 0; k < components; k++){
        clist[k].iters_since_improvement = 0;
    }
    for (i=0; i<snps; i++){
        snpfrag[i].post_hap = 0;
    }
    // RUN THE MAX_CUT ALGORITHM ITERATIVELY TO IMPROVE LIKELIHOOD

    for (iter = 0; iter < MAXITER; iter++) {
        if (VERBOSE)
            fprintf_time(stdout, "PHASING ITER %d\n", iter);
        converged_count = 0;
        for (k = 0; k < components; k++){
            if(VERBOSE && iter == 0)
                fprintf_time(stdout, "component %d length %d phased %d %d...%d\n", k, clist[k].length, clist[k].phased, clist[k].offset, clist[k].lastvar);
            if (clist[k].SCORE > 0)
                converged_count += evaluate_cut_component(Flist, snpfrag, clist, k, slist, HAP1);
            else converged_count++;
        }

        if (converged_count == components) {
            //fprintf(stdout, "Haplotype assembly terminated early because no improvement seen in blocks after %d iterations\n", CONVERGE);
            break;
        }
    }

    // BLOCK SPLITTING
    /*
    new_components = components;
    if (SPLIT_BLOCKS){
        split_count = 0;
        for (k=0; k<components; k++){
            // attempt to split block
            split_count += split_block(HAP1, clist, k, Flist, snpfrag, &new_components);
        }
        if (split_count > 0){
            // regenerate clist if necessary
            free(clist);
            clist = (struct BLOCK*) malloc(sizeof (struct BLOCK)*new_components);
            generate_clist_structure(Flist, fragments, snpfrag, snps, new_components, clist);
        }
        components = new_components;
    }else if(ERROR_ANALYSIS_MODE){
        for (k=0; k<components; k++){
            // run split_block but don't actually split, just get posterior probabilities
            split_block(HAP1, clist, k, Flist, snpfrag, &new_components);
        }
    }

    // PRUNE SNPS
    if (!SKIP_PRUNE){
        likelihood_pruning(snps, Flist, snpfrag, HAP1, CALL_HOMOZYGOUS);
    }
    */

    for (i = 0; i < snps; i++){

        if (HAP1[i] == '-' || snpfrag[i].genotypes[0] == snpfrag[i].genotypes[2]) {
            HAP1[i] = '-';
        } else {
            PS[i] = snpfrag[i].component;
        }
    }

    // FREE UP MEMORY
    //free(HAP1);
    for (i = 0; i < snps; i++) free(snpfrag[i].elist);
    for (i = 0; i < snps; i++) free(snpfrag[i].telist);
    component = 0;
    for (i = 0; i < snps; i++) {
        free(snpfrag[i].flist);
        free(snpfrag[i].alist);
        free(snpfrag[i].jlist);
        free(snpfrag[i].klist);

        if (snpfrag[i].component == i && snpfrag[i].csize > 1) // root node of component
        {
            free(clist[component].slist);
            component++;
        }
    }

    for (i = 0; i < components; i++) free(clist[i].flist);
    free(snpfrag);
    free(clist);
    free(Flist);

    return 0;
}

//int main(){
//}
