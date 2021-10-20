#CC=gcc
#CFLAGS=-I ~/ReadDiff/src ./clip_reads.c -I ~/htslib/ -L ~/htslib/ -lhts

CC = gcc

all: 
	$(CC) ~/ReadDiff/src/main_htslib.c -I ~/htslib -I ~/ReadDiff/src ./clip_reads.c -L ~/htslib/ -lhts -o ReadDiff

clean:
	rm ~/ReadDiff/src/*.o ~/ReadDiff/src/ReadDiff

#all: main_htslib.o clip_reads.o
#	$(CC) main_htslib.o clip_reads.o -o ReadDiff

#clip_reads.o: clip_reads.c clip_reads.h
#	$(CC) -c clip_reads.c
#
#main_htslib.o: 