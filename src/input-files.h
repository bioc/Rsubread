#ifndef __INPUT_FILES_H_
#define __INPUT_FILES_H_

#include "subread.h"

#define GENE_SPACE_BASE 1
#define GENE_SPACE_COLOR 2

#define GENE_INPUT_PLAIN 0
#define GENE_INPUT_FASTQ 1
#define GENE_INPUT_FASTA 2

#include <stdlib.h>
#include <stdio.h>




int chars2color(char c1, char c2);

int genekey2color(char last_base,char * key);

// Convert a string key into an integer key
int genekey2int(char key [], int space_type);

// Open a read file. This function automatically decides its type.
// Return 0 if successfully opened, or 1 if error occurred
int geinput_open(const char * filename, gene_input_t * input);

// Get the next read from the input file
// Return the length of this read or -1 if EOF. 
// The memory space for read_string must be at least 512 bytes.
int geinput_next_read(gene_input_t * input, char * read_name, char * read_string, char * quality_string);

// Close the input file
void geinput_close(gene_input_t * input);

// return the next ATGC char from a input file.
// return 0 if this read segment reaches the end; return -1 if EOF; return -2 if error
int geinput_next_char(gene_input_t * input);

// line buffer has to be at least 300 bytes
// it returns the length of reading
int geinput_readline(gene_input_t * input, char * linebuffer, int conv_to_upper) ;

// read a line into the buff,
// the line should not be longer than 300 chars or the remaining part will be discarded.
// therefore the buff has to be at least 300 chars.
int read_line(FILE * fp, char * buff, int conv_to_upper);

// count the number of reads in a flie
double guess_reads_density(char * fname) ;

// guess the size of the chromosome lib
// return the number of bases, or (-index-1) if the file at the index is not found.
long long int guess_gene_bases(char ** files, int file_number);


void reverse_read(char * ReadString, int Length, int space_type);

void reverse_quality(char * QualtyString, int Length);

#define reverse_char(c)	((c)=='A'?'T':((c)=='G'?'C':((c)=='C'?'G':'A')))

int find_subread_end(int len, int  TOTAL_SUBREADS,int subread) ;
#endif
