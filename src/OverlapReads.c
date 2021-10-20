#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

//-----HTSLIB HEADERS----
#include "htslib/sam.h"

int main(int argc, char** argv){

    // extend reads with matepair if over lap
    // then look for another read with the maximum exact mathc

    samFile* sp = hts_open(argv[1], 'r');
    bam_hdr_t *head = sam_hdr_read(sp);
    bam1_t *bam = bam_init1();
    fprintf(stderr, "Reading Multi-Sequence Alignment\n");
    while(sam_read1(sp, head, bam) >= 0){
        int32_t pos = bam->core.pos + 1;
        // ignoring the contif for now

        uint32_t q2 = bam->core.qual;
        if(q2 > 40){
            uint32_t seq_len = bam->core.l_qseq;
            uint8_t *q = bam_get_seq(bam);

            // Apparently the proper way to cast malloc
            char *sequence = malloc(sizeof * sequence * (seq_len + 1));
            for(size_t i = 0; i < seq_len; i++){
                sequence[i] = seq_nt16_str[bam_seqi(q, i)];
            }
            sequence[seq_len] = "\0"; // terminate with null byte
            // check for overlap regardless of read pair
            


        }

       

    }

    free(sp);
    bam_hdr_destroy(head);
    bam_destroy1(bam);


}