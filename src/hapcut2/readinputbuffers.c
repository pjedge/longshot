#include "readinputbuffers.h"
#include "common.h"
#include <assert.h>

extern int HOMOZYGOUS_PRIOR;

int fragment_compare(const void *a, const void *b) {
    const struct fragment *ia = (const struct fragment*) a;
    const struct fragment *ib = (const struct fragment*) b;
    if (ia->list[0].offset == ib->list[0].offset) {
        return ia->list[ia->blocks - 1].offset + ia->list[ia->blocks - 1].len - ib->list[ib->blocks - 1].offset - ib->list[ib->blocks - 1].len;
        //return ia->blocks - ib->blocks;
    } else return ia->list[0].offset - ib->list[0].offset;
}


int read_fragment_buffer(char** fragmentbuffer, struct fragment* Flist, int fragments) {
    int i = 0, j = 0, k = 0, t = 0, t1 = 0;
    int blocks = 0, type = 0, l = 0, biter = 0, offset = 0;
    char* buffer = fragmentbuffer[0];
    char blockseq[100000];
    //for (i=0;i<MAXBUF;i++) buffer[i] = 0;
    for (i=0;i<100000;i++) blockseq[i] = 0;
    //int num_fields;
    //int expected_num_fields;

    for (i = 0; i < fragments; i++) {
        //		fprintf(stdout,"%s \n",buffer);
        buffer = fragmentbuffer[i];
        j = 0;
        while (buffer[j] != '\0') {
            j++;
        };
        k = 0;
        t = 0;
        type = 0;
        while (k < j) {
            while (k < j && buffer[k] != ' ' && buffer[k] != '\t' && buffer[k] != '\0') {
                blockseq[t] = buffer[k];
                t++;
                k++;
            }
            k++;
            while (k < j && (buffer[k] == ' ' || buffer[k] == '\t')) k++;
            blockseq[t] = '\0';

            if (type == 0) // read the number of blocks in fragment
            {
                blocks = 0;
                for (l = 0; l < t; l++) {
                    blocks = 10 * blocks + (int) (blockseq[l] - 48);
                }
                type = 1;
                Flist[i].blocks = blocks;
                Flist[i].list = (struct block*) malloc(sizeof (struct block)*(blocks));
                biter = 0;

            } else if (type == 1) // read the fragment id, changed to allow dynamic length feb202011
            {
                type = 2; // old format, skip right to reading alleles
            } else if (type == 2 && biter < blocks) {
                offset = 0;
                for (l = 0; l < t; l++) {
                    offset = 10 * offset + (int) (blockseq[l] - 48);
                }
                type = 3;
                Flist[i].list[biter].offset = offset - 1;
                //printf("block %d %d ",biter,offset-1);
            } else if (type == 2 && biter == blocks) {
                offset = 0;
                Flist[i].calls = 0;
                for (l = 0; l < blocks; l++) {
                    for (t1 = 0; t1 < Flist[i].list[l].len; t1++) Flist[i].list[l].pv[t1] = pow(0.1, (float) (blockseq[offset + t1] - QVoffset) / 10);
                    //for (t1=0;t1<Flist[i].list[l].len;t1++) printf("qv %f %d ",pow(0.1,(float)(blockseq[offset+t1]-QVoffset)/10),blockseq[offset+t1]-33);
                    for (t1 = 0; t1 < Flist[i].list[l].len; t1++) Flist[i].list[l].qv[t1] = blockseq[offset + t1];
                    for (t1 = 0; t1 < Flist[i].list[l].len; t1++) Flist[i].list[l].p1[t1] = log10(1.0 - Flist[i].list[l].pv[t1]); // added 03/03/15
                    offset += Flist[i].list[l].len;
                    Flist[i].calls += Flist[i].list[l].len;
                }
            } else if (type == 3) {
                Flist[i].list[biter].hap = (char*) malloc(t + 1);
                Flist[i].list[biter].qv = (char*) malloc(t + 1);
                Flist[i].list[biter].len = t;
                Flist[i].list[biter].pv = (float*) malloc(sizeof (float)*Flist[i].list[biter].len);
                Flist[i].list[biter].p1 = (float*) malloc(sizeof (float)*Flist[i].list[biter].len);

                for (l = 0; l < t; l++) Flist[i].list[biter].hap[l] = blockseq[l];

                type = 2;
                biter++;
            }
            t = 0;
        }
    }

    qsort(Flist, fragments, sizeof (struct fragment), fragment_compare);

    return fragments;
}
