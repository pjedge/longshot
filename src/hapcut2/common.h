
#ifndef _COMMON_H
#define _COMMON_H
#include <stdint.h>

//Tue May 29 23:13:29 PDT 2007
extern int QVoffset;
extern int MINQ;

#define MAXBUF 100000

// given a=log10(x) and b=log10(y), returns log10(x+y)
#define addlogs(a, b) (((a) > (b)) ? ((a) + log10(1.0 + pow(10.0, (b) - (a)))) : ((b) + log10(1.0 + pow(10.0, (a) - (b)))))
// given a=log10(x) and b=log10(y), returns log10(x-y)
#define subtractlogs(a, b) (((a) > (b)) ? ((a) + log10(1.0 - pow(10, (b) - (a)))) : ((b) + log10(1.0 - pow(10.0, (a) - (b)))))

#define flip(allele) if (allele == '1') allele = '0'; else if (allele == '0') allele = '1'

struct block {
    int offset;
    char* hap;
    short len;
    float* pv;
    char* qv;
    float* p1;
};

struct fragment {
    char* id;
    short blocks;
    struct block* list;
    int component;
    float currscore;
    int calls;
    float ll;
    float scores[4]; // added 03/02/15
};

// haplotype block

struct BLOCK {
    int offset;
    int length;
    int phased;
    char* haplotype;
    int* flist;
    int frags;
    float SCORE, bestSCORE, lastSCORE;
    int* slist; // ordered list of variants in this connected component
    int lastvar; // index of the first and last variants in this connected component
    int iters_since_improvement;
};

struct edge {
    int snp;
    int frag;
    char p[2];
    float w;
};

typedef struct EDGE {
    int s, t;
    float w;
} EDGE;

struct SNPfrags {
    int* flist;
    int* jlist; // list of j indexes used to index into Flist[f].list
    int* klist; // list of k indexes used to index into Flist[f].list[j]
    int frags;
    char* alist; // alist is the alleles corresponding to each fragment that covers this SNP
    int component;
    int edges; // those that span the interval between this snp and the next snp
    int csize;
    struct edge* elist;
    int bcomp; // index of clist to which this snp belongs: reverse mapping
    struct edge* telist;
    int tedges; // temporary edge list and number of edges for MIN CUT computation
    int parent;
    float score;
    int heaploc;
};

int fprintf_time(FILE *stream, const char *format, ...);
float phred(float x);
float unphred(float x);

#endif
