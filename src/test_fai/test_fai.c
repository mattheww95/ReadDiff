#include <stdio.h>
#include <stdlib.h>
#include "htslib/faidx.h"



int main(int argc, char** argv){

    faidx_t *fasta_file = fai_load(argv[1]);
    printf("Successfully loaded fasta file.\n");
    // get sequence from the index
    hts_pos_t seq_len;
    //seqlen is 64 bit signed integer
    char* sequence = fai_fetch64(fasta_file, "MN908947.3", &seq_len); // can get full sequence from index or just parts
    //printf("%s\n", sequence);
    printf("%d\n", seq_len);
    free(sequence);
    fai_destroy(fasta_file);
    return 0;
}