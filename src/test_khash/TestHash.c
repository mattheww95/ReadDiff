// To run this program: `echo a bc a cd bc|./this_prog`
#include <stdio.h>
#include <string.h>
#include "htslib/khash.h"


// How to create a dictionary eith khash

KHASH_MAP_INIT_STR(str, int)


#define ARRAY_SIZE(arr)     (sizeof(arr) / sizeof((arr)[0]))



int main(int argc, char *argv[])
{
    char s[4096]; // max string length: 4095 characters
    khash_t(str) *h;
    khint_t k;

    h = kh_init(str);
    while (scanf("%s", s) > 0) {
        int absent;
        //combo test = {some_arr: s, some: 0};
        //combo *new_test = &test;
        int count = 1;
        k = kh_put(str, h, s, &absent);
        kh_val(h, k) = count;
        if (absent) kh_key(h, k) = strdup(s);
        else{
            kh_val(h, k)++;
        }
        // else, the key is not touched; we do nothing
    }
    printf("# of distinct words: %d\n", kh_size(h));
    // IMPORTANT: free memory allocated by strdup() above
    for (k = 0; k < kh_end(h); ++k){
        //printf("%d\n", k);
        if (kh_exist(h, k)){
            printf("%s\t%d\n", kh_key(h, k), kh_val(h, k));
            free((char*)kh_key(h, k));
        }     
    }


;
    kh_destroy(str, h);
    return 0;
}