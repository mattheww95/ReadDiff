#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

//-----HTSLIB HEADERS-----
#include "htslib/sam.h"
#include "htslib/khash.h"
#include "htslib/faidx.h"

//-----USER DEFINED HEADERS-----
#include "clip_reads.h"


#define ARRAY_SIZE(arr)     (sizeof(arr) / sizeof((arr)[0]))
KHASH_MAP_INIT_STR(str, int)
#define MAX_STR_LENGTH 4096

/*
TODO: Apparently I should not be type casting malloc, nor am i using sizeof properly

Fix in later iterations

*/


uint32_t count_combos(const char* combos, uint32_t min){
    // combos is the string of combinations
    // min is the number of mutations required to be put out
    uint32_t tracker_ = 0;
    int muts_ = 0; 
    while(combos[tracker_] != '\0' && muts_ != min){
        if('|' == combos[tracker_]){ // Can probably refactor this to clean it up by defining if not ==
            muts_++;
            tracker_++;
        }else{
            tracker_++;
        }

    }
    return muts_;
}


void usage(){
    fprintf(stderr, "Provided an indexed bam file and fasta with index,\n tabulate all mutations on reads compared to a reference.\n");
    fprintf(stderr, "Arguments are expected in a in the order of:\n");
    fprintf(stderr, "\t1) Path to bamfile. The bam index file is expected to be found in the directory.\n");
    fprintf(stderr, "\t2) Path to reference fasta. The fasta index is expected to be located in the same directory\n");
    fprintf(stderr, "\t3) Minimum number of mutations required to be reported.\n");
}


int main(int argc, char **argv){

    if(argc < 4){
        fprintf(stderr, "Not enough arguments.\n");
        usage();
        return -1;
    }

    // Make a cmd line arg, but define combo lengths to keep
    uint32_t combos_;
    sscanf(argv[3], "%d", &combos_); // make cmd line arg
    
    samFile *fp = hts_open(argv[1], "r");

    bam_hdr_t *h = sam_hdr_read(fp);
    bam1_t *b = bam_init1();

    // read in the fasta information
    faidx_t *ref_fasta = fai_load(argv[2]);




    // ----KAHSH DECLERATIONS----
    // Set up hash table
    char s[MAX_STR_LENGTH]; // max string length might be able to use the malloc'd sequence later
    // Making this a malloc'd will be a todo
    khash_t(str) *hash; 
    khint_t k;
    hash = kh_init(str);
    

    // Fasta inofrmation
    hts_pos_t seq_len; // the length of the buffer the fasta will be read into
    char* prev_chr = h->target_name[b->core.tid]; // get the first sequence
    char* fasta_sequence = fai_fetch64(ref_fasta, prev_chr, &seq_len);
    
    //fprintf(stderr, "Read fasta sequence %s\n", fasta_sequence);
    //printf("Read\tREF\tALT\tGenome_Position\n");
    fprintf(stderr, "Reading Multi-Sequence Alignment\n");
    while(sam_read1(fp, h, b) >= 0){
        int32_t pos = b->core.pos +1; //left most position of alignment in zero based coordianate (+1)
        char *chr = h->target_name[b->core.tid]; //contig name (chromosome), Use this in the FAI program
        
        
        // make sure were reading from the right sequence eachtime, could be multiple chromosomes
        if(strcmp(prev_chr, chr) != 0){
            char* fasta_sequence = fai_fetch64(ref_fasta, chr, &seq_len);
            prev_chr = chr;
        }

        uint32_t len = b->core.l_qseq; //length of the read.
        

        uint8_t *q = bam_get_seq(b); //quality string
        uint32_t q2 = b->core.qual ; //mapping quality

        unsigned char *qseq = (unsigned char *)malloc((len + 1) * sizeof(unsigned char));

        for(size_t i=0; i < len ; i++){
            qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
        }
        qseq[len] = '\0'; // allocate null byte
        // printf("%d\n",b->core.flag);
        // check if is revers from flags  b->core.flag & BAM_FREVERSE
        // check if each segment is mapped accordingly

        //uint16_t reverse =  b->core.flag & BAM_FREVERSE;
        
        if(q2 > 40){ // need to check bit flags to get if segment is unmapped && reverse != 16
            //clip_reads(b, qseq);
            
            uint32_t new_len = cigar_length(b);
            unsigned char* aln_seq = malloc((new_len + 1) * sizeof(unsigned char)); // for the nullbyte!
            unsigned char* ref_seq = malloc((new_len + 1) * sizeof(unsigned char));

            
            uint32_t tracker = 0;
            uint32_t seqeunce_tracker = 0;
            uint32_t ref_seq_tracker = pos - 1; // track in genome position
            uint32_t new_alignment = 0;
            for(uint32_t i = 0; i < b->core.n_cigar; i++){
                uint32_t *cigar = bam_get_cigar(b);
                //printf("%i%c", bam_cigar_oplen(cigar[i]), BAM_CIGAR_STR[bam_cigar_op(cigar[i])]);

                int cigar_l = bam_cigar_oplen(cigar[i]); // number of objects 1 based 
                char cigar_char = BAM_CIGAR_STR[bam_cigar_op(cigar[i])];
                // will need to check for reallocation of memory later
                //printf("%c%d", cigar_char, cigar_l);
                if(cigar_char == 'D'){
                    for(uint32_t del_ = 0; del_ < cigar_l; del_++){
                        aln_seq[tracker] = '-';
                        //seqeunce_tracker++;
                        tracker++;   
                    }
                    for(uint32_t nuc = 0; nuc < cigar_l; nuc++){ // append reference nucleotides
                        ref_seq[new_alignment] = fasta_sequence[ref_seq_tracker];
                        ref_seq_tracker++;
                        new_alignment++;
                    }
                }else if(cigar_char == 'M'){
                    for(uint32_t nuc = 0; nuc < cigar_l; nuc++){ // append read nucleotides
                        aln_seq[tracker] = qseq[seqeunce_tracker];
                        tracker++;
                        seqeunce_tracker++;
                    }
                    for(uint32_t nuc = 0; nuc < cigar_l; nuc++){ // append reference nucleotides
                        ref_seq[new_alignment] = fasta_sequence[ref_seq_tracker];
                        ref_seq_tracker++;
                        new_alignment++;
                    }

                }else if(cigar_char =='I'){
                    // nothing to be done to the read here, need to add insertions to the ref_sequence
                    for(uint32_t nuc = 0; nuc < cigar_l; nuc++){ // append reference nucleotides

                        ref_seq[new_alignment] = '-'; // add in insertions into the references

                        //ref_seq_tracker++; do not increment the reference sequence
                        new_alignment++;
                    }
                    for(uint32_t nuc = 0; nuc < cigar_l; nuc++){ // as still nucleotides present in sequence
                        aln_seq[tracker] = qseq[seqeunce_tracker];
                        tracker++;
                        seqeunce_tracker++;
                    }

                }else{
                    seqeunce_tracker = seqeunce_tracker + cigar_l;
                }
            }
           
            aln_seq[new_len] = '\0';
            ref_seq[new_len] = '\0'; // null terminate the new sequences for security :D

            uint32_t mut_check = 0;
            uint32_t safety_margin = new_len * 4;


            unsigned char mutation_combos[MAX_STR_LENGTH];
            uint16_u combo_count = 0;
            bool mut_present = false;
            while(mut_check <= new_len){
                // deltions are one off compared to ivar
                if(aln_seq[mut_check] != ref_seq[mut_check]){
                    //printf("%c%d%c\t", ref_seq[mut_check], pos + mut_check, aln_seq[mut_check]);
                    mut_present = true;
                    mutation_combos[combo_count] = '|';
                    combo_count++;
                    
                    mutation_combos[combo_count] = ref_seq[mut_check];
                    combo_count++;
                    unsigned char position_buf[50]; // should cover the worlds longest read...
                    int mod_pos = pos + mut_check;
                    
                    int num_chars = sprintf(position_buf, "%d", mod_pos);
                    //uint8_t string_count = 0;
                    for(int v_cars = 0; v_cars < num_chars; v_cars++){
                        //printf("%d\n", mod_pos);
                        mutation_combos[combo_count] = position_buf[v_cars];
                        combo_count++;
                    }
                    //while(position_buf[string_count] != '\0'){
                    //    //printf("Corruption\n");
                    //    mutation_combos[combo_count] = position_buf[string_count];
                    //    combo_count++;
                    //    string_count++;
                    //}
                    mutation_combos[combo_count] = aln_seq[mut_check];
                    combo_count++;
                    //memset(position_buf, 0, 50); // clear the array
                }
                mut_check++;
                //printf("%s\n", mutation_combos);
            }
           
            mutation_combos[combo_count] = '\0';
            if(mut_present){
                //printf("%s\n", mutation_combos);
                int absent;
                int count = 1;
                k = kh_put(str, hash, mutation_combos, &absent);
                
                //printf("Corruption During Hashing %s %d\n", mutation_combos, mut_check);
                if(absent){ 
                    kh_key(hash, k) = strdup(mutation_combos);
                    kh_val(hash, k) = count;
                } // MAKE SURE THIS IS FREED  
                else{
                    kh_val(hash, k)++;
                }
            }
           
            //mutation_combos[combo_count] = '\0';
        }
        // Check for mutations first 
        //forward_read_count++;
    }
    fprintf(stderr, "Number of distinct SNP Combinations: %d\n", kh_size(hash));
    printf("Frequency\tCombination\n");
    for (k = 0; k < kh_end(hash); ++k){

        if (kh_exist(hash, k)){
            const char* mut_term = kh_key(hash, k);
            uint32_t num_muts = count_combos(mut_term, combos_);
            if(num_muts >= combos_){
                printf("%d\t%s\n", kh_val(hash, k), kh_key(hash, k));
                //free((char*)kh_key(hash, k));
            }
            free((char*)kh_key(hash, k));
        }     
    }
    kh_destroy(str,hash);
    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fp);
    free(fasta_sequence);
    fai_destroy(ref_fasta);
   

   return 0;
}





//if(unmapped){
//    printf("Segment unmappeds: %d\n", unmapped);
//
//}else{
//
//    //printf("Does not fail qc: %d\n", unmapped);
//}


//printf("%s\t%d\t%d\n", h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
//printf("%d\t%s\n", b->core.qual, bam_get_seq(b)); //, bam_seqi(bam_get_seq(b), 20));

// trying to get cigar inf
//printf("%s\t%d\t%d\t%s\n",chr,pos,len,qseq);
// qseq is an array of the sequence 
// printf("%d\t%d\n", pos, len);
// printf("%s\n", qseq);
//for(size_t position = 0; position < len; position++){
//    uint32_t new_pos = (pos - 1) + position;
//    printf("%c", fasta_sequence[new_pos]);
//}
// printf("\n");

            //printf("%d\t%d\n", pos, len); // mpos gives map position in genome seems to be working
            //printf("%s\n", qseq);
            //printf("%d\n", b->id);
                    //unit32_t len = b->core.isize;
        //int32_t mpos = b->core.mpos; mpos is mate position


//uint32_t adjusted_length = 0;
            //for(size_t i = 0; i < len; i++){
            //    if(qseq[i] != 'X'){
            //        printf("%c", qseq[i]);
            //        adjusted_length++;
            //    }
            //}
            //printf("\n");
            //int32_t new_pos = pos - 1;
            ////int32_t new_start = (calculate_start_pos(qseq, pos) - 1);
            //for(size_t position = 0; position < adjusted_length; position++){        
            //   printf("%c", fasta_sequence[new_pos]);
            //   new_pos++;
            //}
            //printf("\n\n");

// From converting integer to char array...

                    //while(ref_pos > 10){
                    //    printf("%d\n", ref_pos);
                    //    uint8_t pos_val = ref_pos % 10;
                    //    mutation_combos[combo_count] = pos_val + '0';
                    //    combo_count++;
                    //    ref_pos/=10;
                    //    
                    //}
                    //printf("%d\n", ref_pos);
                    //mutation_combos[combo_count] = ref_pos + '0';
                    //combo_count++;

                    /*
    // Read in reference fasta seqeunce
    FILE *ref_fasta = fopen(argv[2], "r"); // returs a const char to the file
    // to read the file
    // start with some arbitrarily long buffer to read the file into
    //char fasta_sequence[40000]
    //TODO: Fasta header is trimmed currently fix that later
    char nucleotide;
    size_t tracking = 0
    if(ref_fasta == NULL){
        fprintf(stderr, "Could not open file %s\n", argv[2]);
        return 1;
    }


    unsigned char *fasta_sequence = (unsigned char*)malloc(40000 * sizeof(char));


    while ( (nucleotide = fgetc(ref_fasta)) != EOF ){
        fasta_sequence[tracking++] = nucleotide; // TODO: need to realloc if this is getting full
    }
    fasta_sequence [tracking] = '\0'; // add the null byte
    free(ref_fasta);
    */
    
    //uint32_t forward_read_count = 0;