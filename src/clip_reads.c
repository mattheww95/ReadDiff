#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/sam.h"
#include "clip_reads.h"

#define ARRAY_SIZE(arr)     (sizeof(arr) / sizeof((arr)[0]))

void clip_reads(bam1_t *b_, char* sequence_unclipped){
    // No return type as operate on sequence itself
    /// b is an initialized bam file to a bam1_t pointer type
    /// sequence unclipped is a pointer to the sequence to nuke nucleotides of
    /// Make alter sequence array into on that contains mismatches

    // TODO: Alter read apply cigar characterisitcs

    uint32_t tracker = 0;



    for(uint32_t i = 0; i < b_->core.n_cigar; i++){
        // sum of the lengths in cigar is sequence length
        // S soft clipped, H hard clipped, P padding, = sequence match, X mismatch
        uint32_t *cigar = bam_get_cigar(b_);
        printf("%i%c", bam_cigar_oplen(cigar[i]), BAM_CIGAR_STR[bam_cigar_op(cigar[i])]);
        int cigar_l = bam_cigar_oplen(cigar[i])- 1; // number of objects 1 based 
        char cigar_char = BAM_CIGAR_STR[bam_cigar_op(cigar[i])];
        switch (cigar_char)
        {
        case '=': // dont compare values if theyre equal
            //TODO
            break;
        case 'X': // Dont bother comparing vals sequence mismatch
            //TODO
            break;
        case 'S': // nuke the the values in string
            for(uint32_t soft = 0; soft <= cigar_l; soft++){
                sequence_unclipped[tracker] = 'X';
                tracker++;
            }
            break;
        case 'H':
            for(uint32_t hard = 0; hard <= cigar_l; hard++){
                sequence_unclipped[tracker] = 'X';
                tracker++;
            }
            break;
        case 'P':
            //for(uint32_t pad = 0; pad <= cigar_l; pad++){
            //    sequence_unclipped[tracker] = 'X';
            //    tracker++;
            //}
            break;
        case 'M':
            tracker = tracker + cigar_l + 1; // starting from zero
            break;
        case 'I':
            tracker = tracker + cigar_l + 1;
            break;
        case 'D':
            //tracker = tracker + cigar_l; //dont increment deletion as no nucleotides present
            break;
        case 'N':
            tracker = tracker + cigar_l + 1;
            break;
        default:
            break;
        }  
    }
    printf("\n");
    //return sequence_unclipped;
}


    //for(size_t i = 0; i<b->core.n_cigar; i++){
    //    int op = bam_cigar_op(cigar[i]);
    //    int op_len = bam_cigar_oplen(cigar[i]);
    //    
    //    //printf("\n");
    //    //printf("%d\n", b->core.n_cigar);
    //}


uint32_t calculate_start_pos(unsigned char* clipped_sequence, int32_t start_pos){
    /// clipped sequence is the ouptut of the above function
    /// Need to move start positon of the read foward to deal with clipping
    // start pos is the curretn starting position
    int32_t new_start = start_pos;
    for(uint32_t i = 0; i < ARRAY_SIZE(clipped_sequence); i++){
        if(clipped_sequence[i] == 'X'){
            new_start++;
        }else{
            break;
        }
    }
    return new_start;

}



uint32_t cigar_length(bam1_t *b_){
    // Make the read aligned visually appear aligned as if in a vcf 
    // get read length with no cigar ops
    uint32_t new_len = 0;
    for(uint32_t cig = 0; cig < b_->core.n_cigar; ++cig){
        uint32_t *cigar = bam_get_cigar(b_);
        uint32_t cigar_l = bam_cigar_oplen(cigar[cig]);
        char cigar_char = BAM_CIGAR_STR[bam_cigar_op(cigar[cig])];
        if(cigar_char == 'D' || cigar_char == 'M' || cigar_char == 'I'){ // only a deletion stretches the read and adds in dels
            new_len += cigar_l;
        }
    }
    return new_len;
}