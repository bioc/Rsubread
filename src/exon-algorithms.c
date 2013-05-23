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
  
  

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/timeb.h>
#include "input-files.h"
#include "exon-algorithms.h"
#include "sorted-hashtable.h"


#define front2(str, bias)	(*((str)+(bias))+*((str)+1+(bias)))
#define front4(str, bias)	(front2(str, bias)+front2(str, bias+2))
#define front8(str, bias)	(front4(str, bias)+front4(str, bias+4))
#define front16(str, bias)	(front8(str, bias)+front8(str, bias+8))



int match_read_exon(const char read_str[], int read_len, unsigned int potential_position,  gene_value_index_t * my_array_index, int indel_tolerance, int * split_point, char from_end)
{
	int i, bias;
	char read_matchingness [7][1250];

	if(indel_tolerance>3) indel_tolerance = 3;

	for(bias=-indel_tolerance; bias<=indel_tolerance;bias++)
	{
		for(i=0;i<read_len; i++)
		{
			char base_int = base2int(read_str[i]);
			int is_matched_base =  gvindex_match_base(my_array_index, potential_position+i+bias, base_int);
			read_matchingness[bias+indel_tolerance][i] = is_matched_base; 
		}
	}

	int sum_match = 0, test_pos = 0;
	for(i=0; i<read_len ; i+=8)
	{
		int max_matchness = -1;
		int j;
		test_pos =  i;
		for (j=-indel_tolerance; j<=indel_tolerance; j++)
		{
			int m = front8(read_matchingness[j],test_pos);
			if(m > max_matchness)
			{
				max_matchness = m;
			}
		}
		sum_match += max_matchness;
	}
	if (sum_match *1. / i > 0.75)
		return 0;

	sum_match = 0;
	test_pos = 0;
	for(i=0; i<read_len/2+25 ; i+=1)
	{
		int j;
		int max_matchness = -1;
		test_pos = from_end?(read_len - i - 9):i;
		for (j=-indel_tolerance; j<=indel_tolerance; j++)
		{
			int m = front8(read_matchingness[j],test_pos);
			if(m > max_matchness)
			{
				max_matchness = m;
			}
		}

		if (max_matchness <= 5)
		{
			if(sum_match > 28*8)
			{
				if (from_end)
					*split_point = test_pos + 4;
				else
					*split_point = test_pos + 5;
				return 1;
			}
			sum_match = 0;
		}
		else
			sum_match += max_matchness;
	}
	if(sum_match > 28*8)
	{
		if (from_end)
			*split_point = test_pos + 4;
		else
			*split_point = test_pos + 5;
		return 0;
	}
	return 0;
}


void init_exon_arena(gene_exon_allrecords_t * arena, int all_queries)
{
	arena -> max_len = all_queries;
	arena -> read_lens = malloc(sizeof(short) * all_queries);
	arena -> arenas = malloc(sizeof(struct gene_exon_arena) * all_queries);
 
	clear_exon_arena(arena);
}

void add_exon_arena(gene_exon_allrecords_t * arena, int queries, unsigned int potential_position, int my_good_bases, char from_end, char im_reversed, int read_len, int is_good_bound)
{
	struct gene_exon_arena * my_arena = &(arena -> arenas[queries]);
	if (arena -> read_lens[queries] == 0)
	{
		my_arena -> forward_goodbases[0] = 0; 
		my_arena -> backward_goodbases[0] = 0; 
		arena -> read_lens[queries] = read_len ;
	}

	short * goodbases;
	unsigned int * positions;
	char * is_reversed, *is_good;
	if (from_end)
	{
		goodbases = my_arena -> backward_goodbases;
		positions = my_arena -> backward_positions;
		is_reversed = my_arena -> backward_isreversed;
		is_good = my_arena -> backward_isgood;
	}
	else
	{
		goodbases = my_arena -> forward_goodbases;
		positions = my_arena -> forward_positions;
		is_reversed = my_arena -> forward_isreversed;
		is_good = my_arena -> forward_isgood;
	}

	int i;
	short min_good_base = 29999;
	int min_index = -1;
	for (i=0; i<EXON_ARENA_SIZE; i++)
	{
		if(goodbases[i]<1)
			break;
		if(min_good_base > goodbases[i])
		{
			min_good_base = goodbases[i];
			min_index = i;
		}
	}

	if (i == EXON_ARENA_SIZE)
	{
		if (min_good_base < my_good_bases)
		{
			goodbases [min_index] = my_good_bases;
			positions [min_index] = potential_position;
			is_reversed [min_index] = im_reversed;
			is_good [min_index] = is_good_bound;
		}
	}
	else
	{
		goodbases[i] = my_good_bases;
		positions[i] = potential_position;
		is_reversed [i] = im_reversed;
		is_good[i] = is_good_bound;
	}
}

void clear_exon_arena(gene_exon_allrecords_t * arena)
{
	bzero (arena -> read_lens, sizeof(short) * arena->max_len);
}

int find_best_edges(short read_len, char * read_str, struct gene_exon_arena * arena, unsigned int * ex1p, unsigned int * ex2p, short * ex1len, short * ex2len, char * ex1_reversed, char * ex2_reversed, char * ex1_isgood, char * ex2_isgood)
{
	int best_total_len = -1;
	int i, j;

	for (i =0; i < EXON_ARENA_SIZE; i++)
	{
		if(!arena -> forward_goodbases [i]) break;
		for (j =0; j < EXON_ARENA_SIZE; j++)
		{
			if(!arena -> backward_goodbases [j]) break;
			if (arena -> forward_goodbases [i] + arena -> backward_goodbases [j] < read_len + 20) 
			{
				long long int value = 1000000;
				if (arena->forward_positions [i] > arena->backward_positions [j])
					value -=  arena->forward_positions [i] - arena->backward_positions [j];
				else
					value -=  arena->backward_positions [j] - arena->forward_positions [i];
				if (value >  999989) value = 1;
				if (value < 0) value = 2;
				if (value > best_total_len)
				//if(arena -> forward_goodbases [i] + arena -> backward_goodbases [j] > best_total_len)
				{
					//best_total_len = arena -> forward_goodbases [i] + arena -> backward_goodbases [j] ;
					best_total_len = value;
					*ex1p = arena -> forward_positions [i];
					*ex2p = arena -> backward_positions [j];

					*ex1len = arena -> forward_goodbases [i];
					*ex2len = arena -> backward_goodbases [j];

					*ex1_reversed = arena -> forward_isreversed[i];
					*ex2_reversed = arena -> backward_isreversed[j];

					* ex1_isgood = arena -> forward_isgood[i];
					* ex2_isgood = arena -> backward_isgood[j];
				}
			}
		}
	
	}

	return best_total_len >0;

}


int test_read_conjunction(char * read_str,  int read_len, gene_vote_t * vote, gene_value_index_t * my_array_index,  int indel_tolerance, gene_exon_allrecords_t * arena, int queries, char is_reversed, gene_value_index_t * index)
{
	int i, j, from_end;
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			for(from_end = 0; from_end <2; from_end++)
			{
				unsigned int potential_position = vote -> pos[i][j];
				int split_point;
				int is_exoned = match_read_exon(read_str, read_len, potential_position, my_array_index, indel_tolerance, & split_point, from_end);


				if (is_exoned)
				{
					int good_bases;
					if (from_end)
						good_bases = read_len - split_point;
					else
						good_bases = split_point;

					int is_good_b = is_good_boundary(potential_position + split_point,index, &potential_position, is_reversed);
					add_exon_arena(arena, queries, potential_position, good_bases, (from_end  && !is_reversed) || ((!from_end) && is_reversed), is_reversed, read_len, is_good_b);
//					printf("EXON end=%d, rev=%d, pos=%u, rl=%d, spl=%d, ppos=%u, GOOD=%d\n", from_end, is_reversed, potential_position, read_len, split_point, (is_reversed?(read_len-split_point):split_point), is_good_b);
//					printf("EXON: + %d, %c  @%u\n", split_point, (from_end+is_reversed==1)?'B':'F', potential_position);
				}
//				else
//					printf("NOEXON end=%d, rev=%d\n", from_end, is_reversed);
			}
		}
	return 0;
}


int is_good_boundary(unsigned int test_boundary, gene_value_index_t * index, unsigned int * real_boundary, int is_reversed)
{
	int i;
	for (i=-3; i<2; i++)
	{
		char c1 = gvindex_get(index, test_boundary+i);
		char c2 = gvindex_get(index, test_boundary+i+1);
//		printf ("%c%c\n", c1,c2);
		if ((c1=='C' && c2=='T') || (c1=='A' && c2=='C') ||(c1=='A' && c2=='G') || (c1=='G' && c2=='T'))
		{
//			*real_boundary = test_boundary + i;
			return 1;
		}
	}

	return 0;
}
