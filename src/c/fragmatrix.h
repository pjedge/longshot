#ifndef _FRAGMATRIX_H
#define _FRAGMATRIX_H

#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "common.h"

float edge_weight(char* hap, int i, int j, char* p, struct fragment* Flist, int f);

void update_fragscore(struct fragment* Flist, int f, char* h);
void calculate_fragscore(struct fragment* Flist, int f, char* h, float* mec_ll);
float fragment_ll(struct fragment* Flist, int f, char* h, int homozygous, int switch_ix1, int switch_ix2);

//void label_node(struct SNPfrags* snpfrag, int node, int comp);
void add_edges(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components);
void add_edges_fosmids(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int* components);

void update_snpfrags(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps);

int edge_compare(const void *a, const void *b);
int fragment_compare(const void *a, const void *b);

void generate_clist_structure(struct fragment* Flist, int fragments, struct SNPfrags* snpfrag, int snps, int components, struct BLOCK* clist);

#endif
