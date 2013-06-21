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
  
  
#ifndef _EXON_ALGORITHMS_H_
#define _EXON_ALGORITHMS_H_

#include "gene-algorithms.h" 

#define EXON_ARENA_SIZE 20

struct gene_exon_arena
{
	unsigned int forward_positions [EXON_ARENA_SIZE];
	unsigned int backward_positions [EXON_ARENA_SIZE];

	char forward_isreversed [EXON_ARENA_SIZE];
	char backward_isreversed [EXON_ARENA_SIZE];

	char forward_isgood[EXON_ARENA_SIZE];
	char backward_isgood[EXON_ARENA_SIZE];

	short forward_goodbases [EXON_ARENA_SIZE];
	short backward_goodbases [EXON_ARENA_SIZE];

};

typedef struct
{
	int max_len;
	short * read_lens;
	struct gene_exon_arena * arenas;
} gene_exon_allrecords_t;



int match_read_exon(const char read_str[], int read_len, unsigned int potential_position,  gene_value_index_t * my_array_index, int indel_tolerance, int * split_point, char from_end) ;

void init_exon_arena(gene_exon_allrecords_t * arena, int all_queries);
void clear_exon_arena(gene_exon_allrecords_t * arena);

void add_exon_arena(gene_exon_allrecords_t * arena, int queries, unsigned int potential_position, int my_good_bases, char from_end, char is_reversed, int read_len, int is_good_bound);

int find_best_edges(short read_len, char * read_str, struct gene_exon_arena * arena, unsigned int * ex1p, unsigned int * ex2p, short * ex1len, short * ex2len, char * ex1_reversed, char * ex2_reversed, char * ex1_isgood, char * ex2_isgood) ;

int test_read_conjunction(char * read_str,  int read_len, gene_vote_t * vote, gene_value_index_t * my_array_index,  int indel_tolerance, gene_exon_allrecords_t * arena, int queries, char is_reversed, gene_value_index_t * index) ;

int is_good_boundary(unsigned int test_boundary, gene_value_index_t * index, unsigned int * real_bound, int is_reversed);
#endif
