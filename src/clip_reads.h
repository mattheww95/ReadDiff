#ifndef __CLIP_READS__
#define __CLIP_READS__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/sam.h"



void clip_reads(bam1_t *b, char* sequence_unclipped);
uint32_t calculate_start_pos(unsigned char* clipped_sequence, int32_t start_pos);
uint32_t cigar_length(bam1_t *b_);


#endif