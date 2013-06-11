/***************************************************************

   The Subread and Rsubread software packages are free
   software packages:
 
   you can redistribute it and/or modify it under the terms
   of the GNU General Public License as published by the 
   Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Subread is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
   See the GNU General Public License for more details.

   Authors: Drs Yang Liao and Wei Shi

  ***************************************************************/
  
  
#ifndef __INPUT_FILES_H_
#define __INPUT_FILES_H_

#include "subread.h"
#include "hashtable.h"

#define GENE_SPACE_BASE 1
#define GENE_SPACE_COLOR 2

#define GENE_INPUT_PLAIN 0
#define GENE_INPUT_FASTQ 1
#define GENE_INPUT_FASTA 2

#define GENE_INPUT_SAM_SINGLE   93
#define GENE_INPUT_SAM_PAIR_1   94
#define GENE_INPUT_SAM_PAIR_2   95

#include <stdlib.h>
#include <stdio.h>




void fastq_64_to_33(char * qs);

int chars2color(char c1, char c2);

int genekey2color(char last_base,char * key);

// Convert a string key into an integer key
int genekey2int(char key [], int space_type);

// Open a read file. This function automatically decides its type.
// Return 0 if successfully opened, or 1 if error occurred
int geinput_open(const char * filename, gene_input_t * input);

// Open a sam file. Parameter half_no indicates which read of the pair is concerned.
// half_no = 1: the first read; half_no = 2: the last read; half_no = 0: single-end
// Return 0 if successfully opened, or 1 if error occurred
int geinput_open_sam(const char * filename, gene_input_t * input, int half_no);

// Read a line from the input and reward the pointer to the last position.
int geinput_readline_back(gene_input_t * input, char * linebuffer) ;

// Get the next read from the input file
// Return the length of this read or -1 if EOF. 
// The memory space for read_string must be at least 512 bytes.
int geinput_next_read(gene_input_t * input, char * read_name, char * read_string, char * quality_string);
int geinput_next_read_sam(gene_input_t * input, char * read_name, char * read_string, char * quality_string, gene_offset_t* offsets, unsigned int * pos, int * mapping_quality, int * mapping_flags, int need_reversed);

void geinput_jump_read(gene_input_t * input);

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
int read_line(int max_len, FILE * fp, char * buff, int conv_to_upper);

// count the number of reads in a flie
double guess_reads_density(char * fname, int is_sam) ;

// guess the size of the chromosome lib
// return the number of bases, or (-index-1) if the file at the index is not found.
long long int guess_gene_bases(char ** files, int file_number);


void reverse_read(char * ReadString, int Length, int space_type);

void reverse_quality(char * QualtyString, int Length);

unsigned int read_numbers(gene_input_t * input);

//This function returns 0 if the line is a mapped read; -1 if the line is in a wrong format and 1 if the read is unmapped.
int parse_SAM_line(char * sam_line, char * read_name, int * flags, char * chro, unsigned int * pos, char * cigar, int * mapping_quality, char * sequence , char * quality_string, int * rl);

#define reverse_char(c)	((c)=='A'?'T':((c)=='G'?'C':((c)=='C'?'G':'A')))

int find_subread_end(int len, int  TOTAL_SUBREADS,int subread) ;

int break_SAM_file(char * in_SAM_file, char * temp_location, unsigned int * real_read_count, chromosome_t * known_chromosomes, int is_sequence_needed, int base_ignored_head_tail);

int load_exon_annotation(char * annotation_file_name, gene_t ** output_genes, gene_offset_t* offsets );

int is_in_exon_annotations(gene_t *output_genes, unsigned int offset, int is_start);

void colorread2base(char * read_buffer, int read_len);

#endif
