#ifndef READINPUTBUFFERS_H
#define READINPUTBUFFERS_H

#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include "common.h"
#include "fragmatrix.h"

int read_fragment_buffer(char** fragmentbuffer, struct fragment* Flist, int fragments);

int read_vcf_buffer(char** vcfbuffer, struct SNPfrags* snpfrag, int snps);

#endif
