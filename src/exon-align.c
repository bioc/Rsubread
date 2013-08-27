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
#include <assert.h>
#include <getopt.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "hashtable.h"
#include "exon-algorithms.h"
#include "gene-algorithms.h"
#include "gene-value-index.h"
#include "input-files.h"
#include "sorted-hashtable.h"

gene_offset_t _global_offsets;
float accepted_support_rate = 0.3;
int EXON_ALL_THREADS=1;
int EXON_EXTENDING_SCAN=0;
int REPORT_SAM_FILE = 1;
int EXON_JUNCTION_READS_ONLY = 0;

int TOTAL_SUBREADS;
float EXON_MAJOR_HALF_VOTES_RATE = 0.1;
float EXON_MIN_HALF_VOTES_RATE = 0.15;


int MIN_VOTE2_TMP_VAR = 1;
//int EXON_SUBREAD_GAP = 6;
int EXON_LARGE_WINDOW = 60;
int EXON_LONG_READ_LENGTH = 120;
int EXON_MAX_METHYLATION_C_NUMBER = 0;


int ACCEPT_MINOR_SUBREADS;
int INDEX_THRESHOLD;
int EXON_MAX_PAIRED_DISTANCE = 600;
int EXON_MIN_PAIRED_DISTANCE = 50;
int EXON_INDEL_TOLERANCE = 6;
int EXON_QUALITY_SCALE = QUALITY_SCALE_NONE;
int EXON_USE_VALUE_ARRAY_INDEX = 1;
int EXON_FIRST_READ_REVERSE = 0;
int EXON_SECOND_READ_REVERSE = 1;
int EXON_NUMBER_OF_ANCHORS_PAIRED = 50;
int EXON_DIRECT_READS = 0;
int EXON_NO_TOLERABLE_SCAN = 1;
//int EXON_MAX_CIGAR_LEN = 48;

int EXON_IS_STEP1_RUN = 1;
int EXON_IS_STEP2_RUN = 1;

int EXON_FASTQ_FORMAT = FASTQ_PHRED33;

int IS_SAM_INPUT = 0;
int EXON_MIN_HALF_LENGTH = 0;
float EXON_HALF_MATCH_PERCENTAGE =.7f;
double reads_density;
long long int thread_block_locations[200];
unsigned int thread_read_start[200];
unsigned int thread_read_end[200];
unsigned int total_input_reads;

int EXON_FUSION_DETECTION = 0;

#define EXON_DONOR_TEST_WINDOW 17
#define EXON_GROUPING_SIZE 1
#define EXON_MAX_BIGMARGIN_OVERLAPPING  5
//#define DEBUG

#define is_donar_chars_full(cc) (((cc)[0]=='G' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='G') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') || \
			    ((cc)[0]=='C' && (cc)[1]=='T') || \
			    ((cc)[0]=='G' && (cc)[1]=='C') || \
			    ((cc)[0]=='A' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') ) 


#define is_donar_chars_part(cc) (((cc)[0]=='G' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='G') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') || \
			    ((cc)[0]=='C' && (cc)[1]=='T')) 

#define is_donar_chars is_donar_chars_part


typedef struct{
	unsigned int supporting_reads;
	unsigned int feed_supporting_reads;
	char strand;
	unsigned short left_extend;
	unsigned short right_extend;
} exon_junction_t;

typedef struct {
	char i_am_reversed;
	unsigned int connect_to [MAX_EXON_CONNECTIONS];
	char is_opposite_reversed [MAX_EXON_CONNECTIONS];
} connect_to_t;

typedef struct {
	int max_len;

	unsigned int * best_pos1_list;
	unsigned int * best_pos2_list;
	unsigned char * best_vote1_list;
	unsigned char * best_vote2_list;
	char * is_abnormal_list;
	char * is_reversed_list;
	short * half_marks_list;
	short * splice_point_list;
	short * best1_read_start_pos;
	short * best1_read_end_pos;
	float * final_quality;
	char * cigar_string_buffer;
	short * read_coverage_start;
	short * read_coverage_end;
	char * splice_point_offset_1;
	char * splice_point_offset_2;
	char * indel_in_piece1;
	char * indel_in_piece2;


}halves_record_t;

int init_halves_record(halves_record_t* halves_record, int items)
{

	halves_record -> max_len = items;
	halves_record -> best_pos1_list = (unsigned int * )malloc(sizeof(unsigned int)*items);
	halves_record -> best_pos2_list = (unsigned int * )malloc(sizeof(unsigned int)*items);
	halves_record -> best_vote1_list = (unsigned char * )malloc(sizeof(unsigned char)*items);
	halves_record -> best_vote2_list = (unsigned char * )malloc(sizeof(unsigned char)*items);
	halves_record -> is_abnormal_list = (char * )malloc(sizeof(char)*items);
	halves_record -> is_reversed_list = (char * )malloc(sizeof(char)*items);
	halves_record -> half_marks_list = (short * )malloc(sizeof(short)*items);
	halves_record -> splice_point_list = (short * )malloc(sizeof(short)*items);
	halves_record -> splice_point_offset_1 = (char * )malloc(sizeof(char)*items);
	halves_record -> splice_point_offset_2 = (char * )malloc(sizeof(char)*items);
	halves_record -> indel_in_piece1 = (char * )malloc(sizeof(char)*items);
	halves_record -> indel_in_piece2 = (char * )malloc(sizeof(char)*items);

	halves_record -> best1_read_start_pos = (short * )malloc(sizeof(short)*items);
	halves_record -> best1_read_end_pos = (short * )malloc(sizeof(short)*items);
	halves_record -> final_quality = (float *)malloc(sizeof(float)*items);
	halves_record -> read_coverage_start = (short *)malloc(sizeof(short)*items);
	halves_record -> read_coverage_end = (short *)malloc(sizeof(short)*items);



	#define CIGAR_STRING_
	halves_record -> cigar_string_buffer = (char * )malloc(sizeof(char)*items*EXON_MAX_CIGAR_LEN);

	if(!halves_record -> cigar_string_buffer)
	{
		SUBREADputs(MESSAGE_OUT_OF_MEMORY);
		return 1;
	}


	memset(halves_record -> best_vote1_list, 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> best_vote2_list, 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> cigar_string_buffer , 0 , sizeof(char)*items*EXON_MAX_CIGAR_LEN);
	memset(halves_record -> half_marks_list, 0, sizeof(short)*halves_record -> max_len);
	memset(halves_record -> indel_in_piece1 , 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> indel_in_piece2 , 0, sizeof(char)*halves_record -> max_len);
	return 0;
}


void destory_halves_record(halves_record_t* halves_record)
{
	free(halves_record -> best_pos1_list);
	free(halves_record -> best_pos2_list);
	free(halves_record -> best_vote1_list);
	free(halves_record -> best_vote2_list);
	free(halves_record -> is_abnormal_list);
	free(halves_record -> is_reversed_list);
	free(halves_record -> half_marks_list);
	free(halves_record -> splice_point_list);
	free(halves_record -> splice_point_offset_1);
	free(halves_record -> splice_point_offset_2);
	free(halves_record -> indel_in_piece1);
	free(halves_record -> indel_in_piece2);
	free(halves_record -> best1_read_start_pos);
	free(halves_record -> best1_read_end_pos);
	free(halves_record -> final_quality);
	free(halves_record -> read_coverage_start);
	free(halves_record -> read_coverage_end);
	free(halves_record -> cigar_string_buffer);
}
void clear_processed_marks(halves_record_t* halves_record)
{	int i;
	for (i = 0; i < halves_record -> max_len; i++)
	{
		halves_record -> half_marks_list[i] &= ~(IS_PROCESSED_READ | IS_PROCESSED_READ_R2);
		if (!(halves_record -> half_marks_list[i] & IS_FINALISED_PROCESSING))
		{
			halves_record -> best_vote1_list[i] =0;
			halves_record -> best_vote2_list[i]=0;
//			if (i<10)printf ("\nmCLR MARK %d\n", i);
		}
	}
}

void clear_halve_record(halves_record_t* halves_record)
{
	memset(halves_record -> best_vote1_list, 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> best_vote2_list, 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> half_marks_list, 0, sizeof(short)*halves_record -> max_len);
	memset(halves_record -> indel_in_piece1 , 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> indel_in_piece2 , 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> cigar_string_buffer , 0 , sizeof(char)*halves_record -> max_len*EXON_MAX_CIGAR_LEN);
}

void add_best_matching_halves(halves_record_t * halves_record, unsigned int best_pos1, unsigned int best_pos2, unsigned char best_vote1 , unsigned char best_vote2, char is_abnormal, char is_reversed, short splice_point, short half_marks, int read_number, int best_pos_start, int best_pos_end, short read_coverage_start, short read_coverage_end, char indel_in_p1, char indel_in_p2)
{
	if ((halves_record -> best_vote2_list[read_number] == 0 && best_vote2>0) || (halves_record -> best_vote1_list[read_number] < best_vote1))
	{
		halves_record -> best_pos1_list[read_number] = best_pos1;
		halves_record -> best_pos2_list[read_number] = best_pos2;
		halves_record -> best_vote1_list[read_number] = best_vote1;
		halves_record -> best_vote2_list[read_number] = best_vote2;
		halves_record -> is_abnormal_list[read_number] = is_abnormal;
		halves_record -> is_reversed_list[read_number] = is_reversed;
		halves_record -> half_marks_list[read_number] = half_marks;
		halves_record -> splice_point_list[read_number] = splice_point;
		halves_record -> best1_read_start_pos[read_number] = best_pos_start;
		halves_record -> best1_read_end_pos[read_number] = best_pos_end;
		halves_record -> read_coverage_start[read_number] = read_coverage_start;
		halves_record -> read_coverage_end[read_number] = read_coverage_end;
		halves_record -> indel_in_piece1[read_number] = indel_in_p1;
		halves_record -> indel_in_piece2[read_number] = indel_in_p2;

		//SUBREADprintf("ADD_BEST=%d, v1=%d, v2=%d\n", splice_point, best_vote1, best_vote2);
	}
}

int findOverlappingGap(unsigned int matching_position, char * read, int read_len, int breakpoint, int search_to_end,  gene_value_index_t * my_value_array_index)
	// search-to-end : from the breakpoint - 16, scan toward the right hand direction, finding for the ovserlapped area
	// !search-to-end : from the breakpoint + 16, scan toward the left hand direction, finding for the ovserlapped area
{
	#define OVERLAPPING_TEST_LENGTH 10

	int test_start = breakpoint + (search_to_end?-16:16);

	int pos;
	int overlap_buffer = 0;
	
	int overlap_mask = (2<<OVERLAPPING_TEST_LENGTH)-1;
	int overlap_start = -1;
	int overlap_end = -1;

	for (pos = test_start; pos >=OVERLAPPING_TEST_LENGTH && pos < read_len-OVERLAPPING_TEST_LENGTH; pos+= (search_to_end?-1:1))
	{
		overlap_buffer <<=1;
		char base = gvindex_get( my_value_array_index, pos + matching_position);
		if (base == base2int(read[pos]))
			overlap_buffer |= 1;
		else
			overlap_buffer &= 0xfffffffe;
		overlap_buffer &= overlap_mask;

		int i, sum_1=0;
		for (i = 0; i<OVERLAPPING_TEST_LENGTH; i++)
			sum_1 += ((overlap_buffer << i)&1) ;
		if(sum_1*1./OVERLAPPING_TEST_LENGTH>0.65)
		{
			if (search_to_end)
			{
				if(base == base2int(read[pos]) && overlap_start<0)
					overlap_start = pos;

				base = gvindex_get( my_value_array_index, pos + matching_position + OVERLAPPING_TEST_LENGTH);
				if(base == base2int(read[pos + OVERLAPPING_TEST_LENGTH]))
					overlap_end = pos + OVERLAPPING_TEST_LENGTH;

			}else
			{
				if(base == base2int(read[pos]) && overlap_end<0)
					overlap_end = pos;

				base = gvindex_get( my_value_array_index, pos + matching_position - OVERLAPPING_TEST_LENGTH);
				if(base == base2int(read[pos-OVERLAPPING_TEST_LENGTH]))
					overlap_start = pos - OVERLAPPING_TEST_LENGTH;
			}

		}
		else if(overlap_end >0 || overlap_start>0) break;
	}

	return search_to_end?overlap_start:overlap_end;
}

int match_last_bases(int * masks, int is_correct)
{
	*masks =0xff & ((* masks << 1 )| is_correct);
	int tmask = *masks;
	int i, ret=0;
	for(i=0;i<8;i++)
	{
		ret += tmask &1;
		tmask = tmask >>1;
	}
	return ret; 
}



int select_best_matching_halves_maxone(gene_vote_t * vote, unsigned int * best_pos1, unsigned int * best_pos2, int * best_vote1, int * best_vote2, char * is_abnormal, short * half_marks, int * is_reversed_halves, float accept_rate, int read_len, long long int hint_pos, int tolerable_bases, short * read_coverage_start, short * read_coverage_end, char * indel_in_p1, char * indel_in_p2, gehash_data_t max_pos, gene_vote_number_t max_votes, short max_start, short max_end, short max_mask, char * max_indel_recorder, int* best_select_max_votes, int rl)
{
	int best_splicing_point = -1, i,j;
	char * best_chro_name, is_reversed;
	unsigned int best_chro_pos;
	int selected_max_votes = -1;


	is_reversed = (max_mask & IS_NEGATIVE_STRAND)?1:0;
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			char * chro_name;
			char is_partner_reversed;
			unsigned int chro_pos;

			int overlapped_len, overlap_start, overlap_end;
			// All logical conditions

			//if( (vote->votes[i][j] < vote-> coverage_start[i][j]) < 12 && (vote-> coverage_end[i][j] > rl - 12 )) continue;

			is_partner_reversed = (vote->masks [i][j] & IS_NEGATIVE_STRAND) ? 1:0;
			overlap_start = max(max_start , vote->coverage_start[i][j]);
			overlap_end   = min(max_end , vote->coverage_end[i][j]);
			overlapped_len =overlap_end - overlap_start;

			int coverage_len = max_end - max_start + vote->coverage_end[i][j] - vote->coverage_start[i][j];
			if (overlapped_len >0)coverage_len -= overlapped_len;
			//SUBREADprintf("MAX: %d-%d   OTHER %d-%d    COV=%d   OVLP=%d\n", max_start, max_end, vote->coverage_start[i][j], vote->coverage_end[i][j], coverage_len, overlapped_len);



			if(overlapped_len >=14)
				continue;

			long long int dist = vote->pos[i][j];
			dist -= max_pos;

			//SUBREADprintf ("D=%lld\n", abs(dist));
			if (abs(dist)<6)
				continue;

			int support_r1 = (int) (TOTAL_SUBREADS * EXON_MAJOR_HALF_VOTES_RATE); 
			int support_r2 = 1;

			if (max_votes < support_r1 || vote->votes[i][j]<support_r2)
				continue;

			// Same chromosome
			if ((vote->coverage_start[i][j] < max_start) + is_reversed == 1)
			{
				locate_gene_position(max_pos + read_len, &_global_offsets, &best_chro_name, &best_chro_pos);
				locate_gene_position(vote->pos[i][j] , &_global_offsets, &chro_name, &chro_pos);
			}else
			{
				locate_gene_position(max_pos , &_global_offsets, &best_chro_name, &best_chro_pos);
				locate_gene_position(vote->pos[i][j] +read_len, &_global_offsets, &chro_name, &chro_pos);
			}

			if (!(EXON_FUSION_DETECTION || chro_name == best_chro_name))	// The pointers can be compared because they can be the same.
				continue;

			int is_fusion = 0;

			if(is_reversed != is_partner_reversed) is_fusion = 1; 

			if( is_reversed && ((max_pos > vote->pos[i][j]) + (vote->coverage_start[i][j] < max_start) != 1))is_fusion = 1;
			if((! is_reversed) && ((max_pos > vote->pos[i][j]) + (vote->coverage_start[i][j] > max_start) != 1)) is_fusion = 1;

			if(abs(dist) > 500000 || chro_name != best_chro_name)
				is_fusion = 1;

			if (is_fusion && !EXON_FUSION_DETECTION) continue;

			int test_vote_value ;
			test_vote_value = 8888888 +  vote->votes[i][j]* 1000000 - abs(dist);
			if (hint_pos>=0)
			{
				long long int hint_dist = hint_pos;
				hint_dist -= vote->pos[i][j];
				if (abs (hint_dist) < 100000)
					test_vote_value += 100;
				if (abs (hint_dist) < 5000)
					test_vote_value += 100;
			}

			if (test_vote_value<selected_max_votes)continue;
			// Conditions of order of R3 and R5
			*half_marks &= ~IS_REVERSED_HALVES;
			if (vote->coverage_start[i][j] < max_start && (((max_pos < vote->pos[i][j]) && !is_reversed) || ((max_pos > vote->pos[i][j]) && is_reversed) ) )
				*half_marks |= IS_REVERSED_HALVES;
			if (vote->coverage_start[i][j] >= max_end  &&  (((max_pos > vote->pos[i][j]) && !is_reversed) || ((max_pos < vote->pos[i][j]) && is_reversed) ) )
				*half_marks |= IS_REVERSED_HALVES;

			if (vote->coverage_start[i][j] < max_start)
			{
				(*half_marks) = (*half_marks) & ~IS_R1_CLOSE_TO_5;
			}
			else
			{
				(*half_marks) |= IS_R1_CLOSE_TO_5;
			}

			if(max_mask & IS_NEGATIVE_STRAND)
				*half_marks = (*half_marks) |   IS_NEGATIVE_STRAND_R1;
			else
				*half_marks = (*half_marks) &  ~IS_NEGATIVE_STRAND_R1;

			if(vote->masks[i][j] & IS_NEGATIVE_STRAND)
				*half_marks = (*half_marks) |   IS_NEGATIVE_STRAND_R2;
			else
				*half_marks = (*half_marks) &  ~IS_NEGATIVE_STRAND_R2;
	

			
			best_splicing_point = ((vote->coverage_start[i][j] < max_start)? (vote->coverage_end[i][j]):(max_end)) + ((vote->coverage_start[i][j] < max_start)? (max_start):(vote->coverage_start[i][j]));


			best_splicing_point /=2;

			* best_pos1 = max_pos ;
			* best_pos2 = vote->pos[i][j] ;
			* best_vote1 = max_votes ;
			* best_vote2 = vote->votes[i][j] ;
			* read_coverage_start = min(vote->coverage_start[i][j] , max_start);
			* read_coverage_end = max(vote->coverage_end[i][j] , max_end);

			* read_coverage_start = max_start;
			* read_coverage_end = max_end;
			
			int k;
			for(k=0; k<MAX_INDEL_TOLERANCE ; k+=3)
				if(!max_indel_recorder[k+3])break;
			* indel_in_p1 = max_indel_recorder[k+2];

			for(k=0; k<MAX_INDEL_TOLERANCE ; k+=3)
				if(!vote->indel_recorder[i][j][k+3])break;
			* indel_in_p2 = vote->indel_recorder[i][j][k+2];


			* is_reversed_halves = is_reversed;

			if (test_vote_value >=100)
				*half_marks = (*half_marks) | IS_PAIRED_HINTED;
			else
				*half_marks = (*half_marks) & ~(IS_PAIRED_HINTED);

			if (is_fusion)
				*half_marks = (*half_marks)    | IS_FUSION;
			else
				*half_marks = (*half_marks) & ~( IS_FUSION);
	

			selected_max_votes = test_vote_value; 

		}
	*best_select_max_votes = selected_max_votes ;
	return best_splicing_point;
}



int select_best_matching_halves(gene_vote_t * vote, unsigned int * best_pos1, unsigned int * best_pos2, int * best_vote1, int * best_vote2, char * is_abnormal, short * half_marks, int * is_reversed_halves, float accept_rate, int read_len, long long int hint_pos, int tolerable_bases, short * read_coverage_start, short * read_coverage_end, char * indel_in_p1, char * indel_in_p2 , int * max_cover_start, int * max_cover_end, int rl, int repeated_pos_base, int is_negative, char * repeat_record, unsigned int index_valid_range)
{
	unsigned int tmp_best_pos1=0, tmp_best_pos2=0;
	int tmp_best_vote1=0, tmp_best_vote2=0, tmp_is_reversed_halves=0;
	char tmp_is_abnormal=0, tmp_indel_in_p1=0, tmp_indel_in_p2=0;
	short tmp_half_marks=0, tmp_read_coverage_start=0, tmp_read_coverage_end=0;
	int ret = 0, best_ret = 0;	

	int i,j;
	int test_select_votes=-1, best_select_votes = 1000000;
	//int max_minor = 0;

	/*
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if(vote->votes[i][j] < vote->max_vote)continue;
			int ii,jj;
			for (ii=0; ii<GENE_VOTE_TABLE_SIZE;ii++)
				for(jj=0; jj< vote->items[ii]; jj++)
				{
					if(max_minor >= vote->votes[ii][jj]) continue;
					if(ii==i && jj==j)continue;
					long long int dist =  vote->pos[ii][jj];
					dist =abs(dist - vote->pos[i][j]);
					if(dist > 500000)
						continue;
					max_minor = vote->votes[ii][jj];
				}

		}

	int encountered = 0;


	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if(vote->votes[i][j] < vote->max_vote)continue;
			int ii,jj;
			for (ii=0; ii<GENE_VOTE_TABLE_SIZE;ii++)
				for(jj=0; jj< vote->items[ii]; jj++)
				{
					if(max_minor != vote->votes[ii][jj]) continue;
					if(ii==i && jj==j)continue;
					long long int dist =  vote->pos[ii][jj];
					dist =abs(dist - vote->pos[i][j]);
					if(dist > 500000)
						continue;
					encountered++;
				}

		}
	*/

	int repeated_pos = repeated_pos_base;
	int offset_shifting = (rl > 220)?4:0;
	//int encounter = 0;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			/*if((vote->votes[i][j] >=  vote->max_vote -1) && (vote->max_coverage_start >= vote-> coverage_start[i][j] - EXON_MAX_BIGMARGIN_OVERLAPPING ) &&  (vote->max_coverage_end <= vote-> coverage_end[i][j] + EXON_MAX_BIGMARGIN_OVERLAPPING))
				encounter++;*/
			if(repeated_pos_base>=0 && vote->pos[i][j]<=index_valid_range)
				if(vote->votes[i][j] >=  vote->max_vote && repeated_pos < repeated_pos_base+12)
				{
					repeat_record[repeated_pos] = (vote-> coverage_start[i][j] >> offset_shifting);
					repeat_record[repeated_pos+1] = (vote-> coverage_end[i][j] >> offset_shifting);
					repeat_record[repeated_pos+2] = (is_negative?0x80:0) | (vote->votes[i][j]&0x7f);
					repeated_pos+=3;
				}
		}
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if(repeated_pos_base>=0 && vote->pos[i][j]<=index_valid_range)
				if(vote->votes[i][j] ==  vote->max_vote -1 && repeated_pos < repeated_pos_base+12)
				{
					repeat_record[repeated_pos] = (vote-> coverage_start[i][j] >> offset_shifting);
					repeat_record[repeated_pos+1] = (vote-> coverage_end[i][j] >> offset_shifting);
					repeat_record[repeated_pos+2] = (is_negative?0x80:0) | (vote->votes[i][j]&0x7f);
					repeated_pos+=3;
				}
		}


	/*
	if(encounter>=2)
		return 0;
	*/

	ret = select_best_matching_halves_maxone(vote, &tmp_best_pos1, &tmp_best_pos2, &tmp_best_vote1, &tmp_best_vote2,  &tmp_is_abnormal,&tmp_half_marks, &tmp_is_reversed_halves, accept_rate, read_len, hint_pos,  tolerable_bases, &tmp_read_coverage_start, &tmp_read_coverage_end, &tmp_indel_in_p1, &tmp_indel_in_p2, vote -> max_position,  vote->max_vote, vote-> max_coverage_start, vote-> max_coverage_end,  vote-> max_mask, vote->max_indel_recorder, &test_select_votes, rl);
	test_select_votes += vote->max_vote*1000000;
			//SUBREADprintf("TSV=%d\n",test_select_votes);

	if(test_select_votes > best_select_votes)
	{
		best_select_votes = test_select_votes;
		*best_pos1 = tmp_best_pos1;
		*best_pos2 = tmp_best_pos2;
		*is_reversed_halves= tmp_is_reversed_halves;
		
		*best_vote1 = tmp_best_vote1;
		*best_vote2 = tmp_best_vote2;
		*is_abnormal = tmp_is_abnormal;
		*indel_in_p1 = tmp_indel_in_p1;
		*indel_in_p2 = tmp_indel_in_p2;
				
		*half_marks = tmp_half_marks;
		*read_coverage_start = tmp_read_coverage_start;
		*read_coverage_end = tmp_read_coverage_end;

		* max_cover_start = vote-> max_coverage_start;
		* max_cover_end = vote-> max_coverage_end;
		best_ret = ret;
	}		
	return best_ret;
}



int pointercmp_forbed(const void *pointer1, const void *pointer2)
{
	paired_exon_key *p1 = (paired_exon_key *)pointer1;
	paired_exon_key *p2 = (paired_exon_key *)pointer2;
	return !((p1-> big_key == p2 -> big_key) && (p2-> small_key == p1-> small_key));
}

unsigned long pointerHashFunction_forbed(const void *pointer)
{
	paired_exon_key *p  = (paired_exon_key *)pointer;
	return p-> big_key ^ p-> small_key;
}

int pointercmp_forpos(const void *pointer1, const void *pointer2)
{
	return pointer1 != pointer2;
}

unsigned long pointerHashFunction_forpos(const void *pointer)
{
	return (unsigned long) pointer & 0xffffffff;
}




#define ceq(c,t) ((c)[0]==(t)[0] && (c)[1]==(t)[1])
#define c2eq(ch1, ch2, tg1, tg2) ((ceq(ch1, tg1) && ceq(ch2, tg2)) || (ceq(ch1, tg2) && ceq(ch2, tg1)) )

int paired_chars_full(char * ch1, char * ch2, int is_reverse)
{
	if (c2eq(ch1, ch2, "GT", "AG") || c2eq(ch1, ch2, "CT", "AC"))
	{
		if (is_reverse) if (ceq(ch1, "AG") || ceq(ch1, "AC")) return 2;
		if (!is_reverse) if (ceq(ch1, "CT") || ceq(ch1, "GT")) return 2;
	}
	else if ( c2eq(ch1, ch2,"GC","AG") || c2eq(ch1, ch2,"GC","CT") || c2eq(ch1, ch2,"AT","AC") || c2eq(ch1, ch2,"GT","AT"))
	{
		if (is_reverse) if (ceq(ch1, "GC") || ceq(ch1, "AT")  || ceq(ch1, "AG") || ceq(ch1, "AC")) return 1;
		if (!is_reverse) if (ceq(ch1, "GC") || ceq(ch1, "AT") ||ceq(ch1, "GT") || ceq(ch1, "CT")) return 1;
	}
	return 0;
}

int paired_chars_part(char * ch1, char * ch2, int is_reverse)
{
	if (c2eq(ch1, ch2, "GT", "AG") || c2eq(ch1, ch2, "CT", "AC"))
	{
		if (is_reverse) if (ceq(ch1, "AG") || ceq(ch1, "AC")) return 1;
		if (!is_reverse) if (ceq(ch1, "CT") || ceq(ch1, "GT")) return 1;
	}
	return 0;
}

#define  paired_chars paired_chars_part


void get_chro_2base(char *buf, gene_value_index_t * index, unsigned int pos, int is_negative_strand)
{
	gvindex_get_string (buf, index, pos, 2, is_negative_strand);
}


// pos1 must be small than pos2.
int test_donor(char *read, int read_len, unsigned int pos1, unsigned int pos2, int guess_break_point, char negative_strand, int test_range, char is_soft_condition, int EXON_INDEL_TOLERANCE, int* real_break_point, gene_value_index_t * my_value_array_index, int indel_offset1, int indel_offset2, int is_reversed, int space_type, int confidence, int * best_donor_score, int * is_GTAG)
{
	int bps_pos_x;
	int search_start = guess_break_point - test_range ;
	int search_end   = guess_break_point + test_range ;
	char h1_2ch[3], h2_2ch[3];

	h1_2ch[2] = h2_2ch[2]=0;
	search_start=max(10, search_start);
	search_end = min(read_len-10, search_end);
	int best_break = -1;
	int min_x = -9099;

	for (bps_pos_x = search_start; bps_pos_x < search_end ; bps_pos_x ++)
	{
		int paired_score = 0;
		get_chro_2base(h1_2ch, my_value_array_index, pos1 - indel_offset1+ bps_pos_x , is_reversed);
		get_chro_2base(h2_2ch, my_value_array_index, pos2 - 2 - indel_offset2 + bps_pos_x, is_reversed);


		//if(!is_reversed)
		//SUBREADprintf("C1=%s @%u, C2=%s @%u\n",h1_2ch, pos1 + bps_pos_x, h2_2ch,pos2 - 2 + indel_offset + bps_pos_x);
		if(h1_2ch[0]==h2_2ch[0] && h1_2ch[1]==h2_2ch[1]) continue;

		if(is_donar_chars(h1_2ch) && is_donar_chars(h2_2ch))
		{

			paired_score = paired_chars(h1_2ch, h2_2ch, is_reversed);

			if(paired_score)
			{
				int m1, m2, x1, x2;
				int break_point_half = is_reversed?(read_len - bps_pos_x):bps_pos_x;
				int first_exon_end,second_half_start;
				int donar_conf_len = 0;

				donar_conf_len = min(break_point_half , EXON_DONOR_TEST_WINDOW);
				donar_conf_len = min(read_len - break_point_half, donar_conf_len);
				//SUBREADprintf("DONOR_CONF_LEN=%d\n", donar_conf_len);

				if (is_reversed)
				{
					first_exon_end = pos2 + bps_pos_x - indel_offset2;
					second_half_start = pos1 + bps_pos_x- indel_offset1;

					m1 = match_chro(read + break_point_half - donar_conf_len , my_value_array_index, first_exon_end, donar_conf_len, is_reversed, space_type);
					m2 = match_chro(read + break_point_half , my_value_array_index, second_half_start-donar_conf_len , donar_conf_len, is_reversed, space_type);

					x1 = match_chro(read + break_point_half ,  my_value_array_index, first_exon_end - donar_conf_len, donar_conf_len , is_reversed, space_type);
					x2 = match_chro(read + break_point_half - donar_conf_len ,  my_value_array_index, second_half_start , donar_conf_len, is_reversed, space_type);
				}
				else
				{
					first_exon_end = pos1 + bps_pos_x - indel_offset1;
					second_half_start = pos2 + bps_pos_x - indel_offset2;

					m1 = match_chro(read + break_point_half - donar_conf_len, my_value_array_index, first_exon_end-donar_conf_len , donar_conf_len, is_reversed, space_type);
					m2 = match_chro(read + break_point_half , my_value_array_index, second_half_start, donar_conf_len, is_reversed, space_type);

					x1 = match_chro(read + break_point_half ,  my_value_array_index, first_exon_end, donar_conf_len , is_reversed,space_type);
					x2 = match_chro(read + break_point_half - donar_conf_len,  my_value_array_index, second_half_start - donar_conf_len, donar_conf_len , is_reversed,space_type);
				}

				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
				{
					SUBREADprintf("DONOR TEST STR=%s, %s ; pos=%d    %d %d ; M=%d %d ; X=%d %d\n", h1_2ch, h2_2ch, bps_pos_x, indel_offset1, indel_offset2, m1, m2, x1, x2);
				}
				#endif
	
				int threshold = 3;
				if (paired_score == 1)
					threshold = 3;

				#ifdef QUALITY_KILL
				if (m1 >= donar_conf_len-1    && m2>=donar_conf_len-1 )
					if(x1<donar_conf_len - threshold  && x2<donar_conf_len- threshold )
				#else
				if (m1 >= donar_conf_len-1    && m2>=donar_conf_len -1)
					if(x1<donar_conf_len - threshold  && x2<donar_conf_len - threshold)
				#endif
					{
						int score =  3000-(x1 + x2) + (m1+ m2) ;
						if (min_x < score)
						{
							min_x = score;
							best_break = bps_pos_x;
							*is_GTAG = 1==((is_reversed) + (h1_2ch[0]=='G' || h1_2ch[1]=='G'));	//"GT" or "AG"
							*best_donor_score = score;
						}
					}
			}
		}
	}

	if (best_break>0)
	{
				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
					SUBREADprintf("SELECRED!!!_BREAKPOINT=%d, RAW POS=%u,%u, R=%s\n",  best_break, pos1 , pos2, read);
				#endif
		//SUBREADprintf ("FINAL BREAK: %d   ; REV = %d\n ", best_break, is_reversed);
		*real_break_point = best_break;
		return 1;
	}
	else
	{
				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
					SUBREADprintf("KILLED!!!_BREAKPOINT=%d, R=%s\n",  best_break+ pos1, read);
				#endif
	}
	return 0;
}




unsigned int get_grouped_position(HashTable * pos_table, unsigned int pos)
{

	int delta_pos;
	unsigned int group_anchor = pos / EXON_GROUPING_SIZE;
	unsigned int grouped_pos = 0;
	for (delta_pos = 0 ; delta_pos < 1; delta_pos ++)
	{
		grouped_pos = (unsigned int) (HashTableGet(pos_table, (NULL + group_anchor)) - NULL);
		if(grouped_pos) break;
		group_anchor ++;
	}
	if(!grouped_pos)
	{
		grouped_pos = pos;
		HashTablePut(pos_table, (NULL+(pos / EXON_GROUPING_SIZE)), (NULL +grouped_pos));
	}
//	else	SUBREADprintf("REUSE:%d\n", grouped_pos);
	return grouped_pos;
}

#define EXON_EXPLAIN_DEPTH 5
//#define DEBUG

void junction_tree_b_explorer(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table, char * read , int rl, int full_rl ,unsigned int read_tail_pos, int number_of_piece, unsigned int * cigar_recorder, int total_matched_bases, int * result_total_matched_bases,unsigned int * result_cigar_recorder, int *result_number_pieces, gene_value_index_t * my_value_array_index, gene_input_t * ginp,int subread_votes_cover_start, char * quality_str, int quality_format, float match_score, unsigned int total_jump_length, unsigned int * max_jump_length)
{
	unsigned int delta_pos;
	unsigned int iii,group_anchor = read_tail_pos / EXON_GROUPING_SIZE-1;
	float max_matched_bases_rate = -1;
	int max_matched_bases = -1;
	int indels=0, max_indels=0, indel_point = 0, max_indel_point=0;
	connect_to_t * next_jump = NULL;
	int max_piece_len = rl, max_grouped_pos=read_tail_pos;
	max_matched_bases = -1.;//*match_chro(read, my_value_array_index, read_tail_pos - rl, rl, 0, ginp->space_type)/rl;
	int max_piece_quality_good = 0;


	int matched_bases =  match_indel_chro_to_back(read, my_value_array_index, read_tail_pos, rl,&indels, &indel_point, EXON_INDEL_TOLERANCE, 0);
	int effect_tested_len = rl -  max(0, indels);

	max_matched_bases_rate =  (matched_bases)*1./effect_tested_len ;
	max_indels = indels;
	max_indel_point = indel_point;
	max_piece_quality_good = 1;
	

	max_matched_bases_rate = max(-1965, max_matched_bases_rate);

	//if (max_matched_bases *1. / rl < EXON_HALF_MATCH_PERCENTAGE)max_matched_bases =-1;

	unsigned int low_border = my_value_array_index -> start_base_offset;
	unsigned int high_border = my_value_array_index -> start_base_offset + my_value_array_index -> length; 


	if(rl > 7)
	for (delta_pos = 0 ; delta_pos < rl ; delta_pos ++)
	{
		unsigned int grouped_pos = (unsigned int) (HashTableGet(pos_table, (NULL + read_tail_pos - delta_pos)) - NULL);

		if (grouped_pos && grouped_pos <= read_tail_pos)
		{
			int test_piece_len = read_tail_pos - grouped_pos;
			if(test_piece_len < 8) continue;
			
			matched_bases = match_indel_chro_to_back(read+rl-test_piece_len, my_value_array_index, read_tail_pos - test_piece_len, test_piece_len,&indels, &indel_point, EXON_INDEL_TOLERANCE,  test_piece_len - rl);

			#ifdef DEBUG
			SUBREADprintf("BSEARCH: INDEL=%d\n",indels);
			#endif

			effect_tested_len = (test_piece_len- max(0, indels));

			float test_matched_bases_rate = (matched_bases)*1./effect_tested_len ;
			if ( test_matched_bases_rate > max_matched_bases_rate)
			{
				connect_to_t * connect_to = (connect_to_t *)HashTableGet(connection_table, NULL+grouped_pos);
				int accepted = 0;
				for (iii=0; iii<MAX_EXON_CONNECTIONS;iii++)
				{
					unsigned int connect_to_iii = connect_to -> connect_to[iii];
					if (! connect_to_iii)break;
					if (connect_to_iii < grouped_pos)
					{
						if(connect_to_iii - (rl-test_piece_len) >= low_border && connect_to_iii < high_border)
						{
							accepted=1;
							break;
						}
					}
				}
				if (accepted && (grouped_pos > read_tail_pos - rl +  EXON_MIN_HALF_LENGTH/* -1*/)){
					max_matched_bases_rate = test_matched_bases_rate;
					max_matched_bases = matched_bases;
					next_jump = connect_to;
					max_piece_len = test_piece_len;
					max_grouped_pos = grouped_pos;
					max_indels = indels;
					max_indel_point = indel_point;
					max_piece_quality_good = 1;
				}
			}

		}

		group_anchor--;
	}

	int digged = 0;
	if(next_jump && EXON_EXPLAIN_DEPTH > number_of_piece)
	{
		cigar_recorder [number_of_piece * 4+1] = rl - max_piece_len + max_indels;
		cigar_recorder [number_of_piece * 4+3] = ((max_indel_point + cigar_recorder [number_of_piece*4+1]) << 16) + (0xffff&( max_indels << 4)) + (max_piece_quality_good << 3);
		#ifdef DEBUG
		SUBREADprintf("KKK1: %d, %d, %u\n", max_indel_point + cigar_recorder [number_of_piece*4+1], max_indels, cigar_recorder [number_of_piece * 4+3] );
		#endif
		for (iii=0; iii<MAX_EXON_CONNECTIONS;iii++)
		{
			unsigned int connect_to_iii = next_jump -> connect_to[iii];
			if (! connect_to_iii)break;

			if (connect_to_iii > max_grouped_pos)
				continue;
			if (max_grouped_pos - connect_to_iii < 10)
				continue;
			if(!(connect_to_iii - (rl-max_piece_len) >= low_border && connect_to_iii < high_border))
				continue;

			cigar_recorder [number_of_piece * 4+4] = connect_to_iii;
			cigar_recorder [number_of_piece * 4+6] = rl + max_indels - max_piece_len;

			unsigned int this_jump = max_grouped_pos - connect_to_iii;

			#ifdef DEBUG
			SUBREADprintf ("DIG-IN B-SEARCH REMAIN_LEN=%d, POS_TAIL=%u\n", rl + max_indels- max_piece_len, connect_to_iii);
			#endif
			junction_tree_b_explorer(bed_table, pos_table, connection_table, read, rl + max_indels - max_piece_len, full_rl, connect_to_iii , 1+number_of_piece , cigar_recorder, total_matched_bases + max_matched_bases, result_total_matched_bases , result_cigar_recorder, result_number_pieces,  my_value_array_index, ginp, subread_votes_cover_start, quality_str , quality_format, match_score, total_jump_length + this_jump, max_jump_length);
			digged = 1;
		}
	}
		// the "read_tail_pos" is the exon that matches through this read; it has to be well matched to accept the read.
	//if (!digged)
	{
		indels = 0;
		indel_point =0;


		int max_piece_quality_good = 1;

		int matched_bases = match_indel_chro_to_back(read, my_value_array_index, read_tail_pos - rl, rl, &indels, &indel_point, EXON_INDEL_TOLERANCE, 0);
		if ( matched_bases > -1965 * rl )
		{
			//put this explaination to the result, if the matched bases is the max.
			if ((matched_bases +total_matched_bases > (*result_total_matched_bases)) || (matched_bases +total_matched_bases ==  (*result_total_matched_bases) && total_jump_length < (*max_jump_length)))
			{
				cigar_recorder [number_of_piece * 4 + 1] = 0;
				cigar_recorder [number_of_piece * 4 + 3] = ((indel_point + cigar_recorder [number_of_piece*4+1]) << 16) + (0xffff&(indels << 4)) + (max_piece_quality_good << 3);

				memcpy(result_cigar_recorder, cigar_recorder, 138);
				(*result_total_matched_bases) = matched_bases + total_matched_bases;
				(*result_number_pieces) = number_of_piece+1;
				(*max_jump_length) = total_jump_length;
			}
		}
	}
}


void junction_tree_f_explorer(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table, char * read , int rl, int full_rl ,unsigned int read_head_pos, int number_of_piece, unsigned int * cigar_recorder, int total_matched_bases, int * result_total_matched_bases,unsigned int * result_cigar_recorder, int *result_number_pieces, gene_value_index_t * my_value_array_index, gene_input_t * ginp, int next_head_modify, char * quality_str, int quality_format, float match_score, unsigned int total_jump_length, unsigned int * max_jump_length)
{
	unsigned int delta_pos;
	unsigned int iii,group_anchor = read_head_pos / EXON_GROUPING_SIZE-1;
	float max_matched_bases_rate = -1;
	int max_matched_bases = -1;
	connect_to_t * next_jump = NULL;
	int max_piece_len = 0;
	int max_grouped_pos = read_head_pos;
	max_matched_bases = -1.;//match_chro(read, my_value_array_index, read_head_pos, rl, 0, ginp->space_type)*1./rl;
	int indels = 0, max_indels=0, indel_point = 0, max_indel_point = 0;
	int max_piece_quality_good = 0;



	int matched_bases =  match_indel_chro_to_front(read, my_value_array_index, read_head_pos, rl,&indels, &indel_point, EXON_INDEL_TOLERANCE, rl);
	int effect_tested_len = rl -  max(0, indels);

	max_matched_bases_rate =  (matched_bases)*1./effect_tested_len;
	max_indels = indels;
	max_indel_point = indel_point;
	max_piece_quality_good = 1;
	
	max_matched_bases_rate = max(-1965, max_matched_bases_rate);
	//if (max_matched_bases *1. / rl < EXON_HALF_MATCH_PERCENTAGE)max_matched_bases =-1;

	unsigned int low_border = my_value_array_index -> start_base_offset;
	unsigned int high_border = my_value_array_index -> start_base_offset + my_value_array_index -> length; 

	if (rl>7)
	for (delta_pos = 0 ; delta_pos < rl ; delta_pos ++)
	{
		unsigned int grouped_pos = (unsigned int) (HashTableGet(pos_table, (NULL + delta_pos + read_head_pos)) - NULL);

		//SUBREADprintf("TESTPOS=%u > %u\n", grouped_pos, read_head_pos);

		if (grouped_pos && grouped_pos >= read_head_pos)
		{
			int test_piece_len = grouped_pos - read_head_pos;
			if(test_piece_len < 8) continue;
			int matched_bases;

			matched_bases = match_indel_chro_to_front(read, my_value_array_index, read_head_pos, test_piece_len, &indels, &indel_point, EXON_INDEL_TOLERANCE, rl);


			#ifdef DEBUG
			SUBREADprintf("TEST match = %d ; piece_len = %d; GRP_POS=%u; \n", matched_bases, test_piece_len, grouped_pos);
			SUBREADprintf("FSEARCH: INDEL=%d\n",indels);
			#endif

			int effect_tested_len = (test_piece_len- max(0, indels));

			float test_matched_bases_rate = (matched_bases)*1./effect_tested_len;

			if (test_matched_bases_rate > max_matched_bases_rate)
			{
				connect_to_t * connect_to = (connect_to_t *)HashTableGet(connection_table, NULL+grouped_pos);
				int accepted = 0;
				for (iii=0; iii<MAX_EXON_CONNECTIONS;iii++)
				{
					unsigned int connect_to_iii = connect_to -> connect_to[iii];
					if (! connect_to_iii)break;
					if (connect_to_iii > grouped_pos)
					{
						if(connect_to_iii >= low_border && connect_to_iii + rl < high_border)
						{
							accepted=1; 
							break;
						}
					}
				}
				if (accepted && (grouped_pos < read_head_pos + rl - EXON_MIN_HALF_LENGTH/* +1*/)){
					max_matched_bases_rate = test_matched_bases_rate;
					max_matched_bases = matched_bases ;// + test_piece_len;
					next_jump = connect_to;
					max_piece_len = test_piece_len;
					max_grouped_pos = grouped_pos;
					max_indels = indels;
					max_indel_point = indel_point;
					max_piece_quality_good = 1;
				}
			}

		}

		group_anchor++;
	}
	
	int digged = 0;
	if(next_jump && EXON_EXPLAIN_DEPTH > number_of_piece)
	{
		cigar_recorder [number_of_piece * 4+2] = full_rl - rl + max_piece_len - max_indels;
		cigar_recorder [number_of_piece * 4+3] = ((max_indel_point + cigar_recorder [number_of_piece*4+1]) << 16) + (0xffff&( max_indels << 4))  + (max_piece_quality_good << 3);

#ifdef DEBUG
		SUBREADprintf("KKK3: %d, %d, %u\n", max_indel_point + cigar_recorder [number_of_piece*4+1], max_indels, cigar_recorder [number_of_piece * 4+3] );
#endif

		for (iii=0; iii<MAX_EXON_CONNECTIONS;iii++)
		{
			unsigned int connect_to_iii = next_jump -> connect_to[iii];
			if (! connect_to_iii)break;

			if (connect_to_iii < max_grouped_pos)
				continue;
			if (connect_to_iii - max_grouped_pos < 10)
				continue;
			if(!(connect_to_iii >= low_border && connect_to_iii + rl < high_border))
				continue;

			cigar_recorder [number_of_piece * 4+4] = connect_to_iii;
			cigar_recorder [number_of_piece * 4+5] = full_rl - rl - max_indels + max_piece_len;
			next_head_modify = 0;
#ifdef TEST_TARGET
			if(read_head_pos >2216991666 - 10000 && read_head_pos < 2216991666 + 10000)
			SUBREADprintf ("DIG-IN:match/piece_len = %d/%d, rl=%d head_pos=%u  READ=%s\n", max_matched_bases, max_piece_len, rl - max_piece_len, connect_to_iii, read+max_piece_len - max_indels, rl);
#endif
			//SUBREADprintf("DIGIN: TOTAL=%d NEWADD=%d\n", total_matched_bases, max_matched_bases);
			//SUBREADprintf("F-MAXINDEL=%d ; NEXT_HEAD_MODIFY=%d ; #PIECES=%d\n", max_indels, next_head_modify, number_of_piece );
			unsigned int this_jump_size = connect_to_iii - max_grouped_pos;
			junction_tree_f_explorer(bed_table, pos_table, connection_table, read+max_piece_len-max_indels, rl - max_piece_len + max_indels, full_rl , connect_to_iii , 1+number_of_piece , cigar_recorder, total_matched_bases + max_matched_bases, result_total_matched_bases , result_cigar_recorder,result_number_pieces, my_value_array_index, ginp, 0 , quality_str + max_piece_len-max_indels , quality_format, match_score, total_jump_length + this_jump_size, max_jump_length);
			digged = 1;
		}
	}
		// the "read_head_pos" is the exon that matches through this read; it has to be well matched to accept the read.
	//if(!digged)
	{
		indels = 0;
		indel_point =0;



		int max_piece_quality_good = 1;

		int matched_bases = match_indel_chro_to_front(read, my_value_array_index, read_head_pos, rl, &indels, &indel_point, EXON_INDEL_TOLERANCE, rl);

		if ( matched_bases > -1965 * rl )
		{
			//put this explaination to the result, if the matched bases is the max.
			if ((matched_bases + total_matched_bases > (*result_total_matched_bases)) || (matched_bases + total_matched_bases ==  (*result_total_matched_bases) && total_jump_length < (*max_jump_length) ) )
			{
				cigar_recorder [number_of_piece * 4 + 2] = full_rl;
				cigar_recorder [number_of_piece * 4 + 3] = ((indel_point + cigar_recorder [number_of_piece*4+1]) << 16) +  ((0xfff&indels) << 4)  + (max_piece_quality_good << 3);

#ifdef TEST_TARGET
				if(read_head_pos >2216991666 - 10000 && read_head_pos < 2216991666 + 10000)
				SUBREADprintf("KKK4: %d+%d, %d, %u ; Tests: %d %d %d %u %s %d TOL=%d\n", indel_point,cigar_recorder [number_of_piece*4+1], indels, cigar_recorder [number_of_piece * 4+3] , matched_bases + total_matched_bases, matched_bases + total_matched_bases,  total_jump_length, read_head_pos, read , rl, EXON_INDEL_TOLERANCE);
#endif

				memcpy(result_cigar_recorder, cigar_recorder, 138);
				(*result_total_matched_bases) = matched_bases + total_matched_bases;
				(*result_number_pieces) = number_of_piece+1;
				(*max_jump_length) = total_jump_length;
			}
		}
	}
}

struct feed_exonbed_parameters
{
	int thread_id;

	HashTable * bed_table;
	HashTable * pos_table;
	HashTable * connection_table;
	halves_record_t * halves_record;
	gene_input_t* ginp;
	gene_input_t * ginp2;
	FILE * out_fp;
	char * index_prefix;
	unsigned int processed_reads;
	unsigned long long int all_processed_reads;
	unsigned long long int *succeed_reads;
	gene_value_index_t ** my_value_array_index_set;
	int all_tables;
	int tolerable_scan;
	long long int * ginp1_end_pos;
	long long int * ginp2_end_pos;


	pthread_spinlock_t * init_lock;
};




void explorer_junc_exonbed(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table, halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads, gene_value_index_t ** my_value_array_index_set, int all_tables, int tolerable_scan, int thread_id, long long int * ginp1_end_pos, long long int * ginp2_end_pos)
{
	int i, ic=0, j, result_head_modify=-1;
	//return;
	gene_value_index_t * my_value_array_index;

	int my_range_start , my_range_end;

	if(thread_id==0)
		SUBREADprintf("Realigning junction reads:\n");
	int mapping_quality=0, mapping_flags=0;
	double reads_per_thread = processed_reads*(1+(ginp2!=NULL)) * 1.0 / EXON_ALL_THREADS;
	my_range_start = max(0,(int)(reads_per_thread * thread_id - 0.5));
	my_range_end = min((int)(reads_per_thread*(1+thread_id) + 0.5 -1 ), processed_reads*(1+(ginp2!=NULL)));
	

	if(thread_block_locations[thread_id]>=0)
	{

		fseeko(ginp->input_fp, thread_block_locations[thread_id], SEEK_SET);
		i = my_range_start;
	}
	else
	{
		for (i=0; i<my_range_start; i++)
		{
			geinput_jump_read(ginp);
	
		}

		thread_block_locations[thread_id] = ftello(ginp->input_fp);
	}


	while (1)
	{
		char nameb[1201], inb[1201], qualityb[1201];
		int rl = -1;
		int is_settle = 1;


		char negative_strand = halves_record -> is_reversed_list[i];
		if (i < my_range_end)
		{

			rl = geinput_next_read_sam(ginp, nameb, inb, qualityb, NULL, NULL, &mapping_quality, &mapping_flags, (ginp2 && (i % 2)) + negative_strand == 1);
		}

		/*
		if(memcmp(nameb, "HWI-ST667_0137:6:2308:15204:7849#CCGTCC", strlen("HWI-ST667_0137:6:2308:15204:7849#CCGTCC"))==0)
			SUBREADprintf("\n\n PROCESSED: HWI-ST667_0137:6:2308:15204:7849#CCGTCC [%d] by %d \n", i, thread_id);
		*/

		if(rl<0)
		{
			if(ginp1_end_pos)
			{
				*ginp1_end_pos = ftello(ginp->input_fp);
			}
	
			break;
		}



		if(i % (processed_reads/14) ==0 && i>1 && thread_id==0)
			print_text_scrolling_bar("", i*EXON_ALL_THREADS*1./(1+(ginp2!=NULL))/processed_reads, 80, &ic);

		if(halves_record->best_vote1_list[i] >=MAX_QUALITY_TO_EXPLORER_JUNCTION) 
		{
			i++;
			continue;
		}


		if(halves_record -> best_vote1_list[i]< TOTAL_SUBREADS * EXON_MAJOR_HALF_VOTES_RATE)
		{
			i+=1;
			continue;
		}


		unsigned int pos = halves_record -> best_pos1_list[i];

		my_value_array_index = NULL;
		for(j=0; j<all_tables; j++)
		{
			unsigned int this_index_last_base = my_value_array_index_set[j] -> start_base_offset + my_value_array_index_set[j] -> length;
			if( pos  > my_value_array_index_set[j] -> start_base_offset + (j?900000:0) && pos  < this_index_last_base - ((j==all_tables-1)? 0:900000))
			{
				my_value_array_index = my_value_array_index_set[j];
				break;
			}

		}

		if(!my_value_array_index)
		{
			i++;
			continue;
		}

		if(0 == EXON_IS_STEP1_RUN)
			halves_record -> best_vote2_list[i] = 0;


		unsigned int explain_buff[138], explain_result[138];
		int result_number_pieces = 0,  result_total_matched_bases = 0, total_pieces = 0, total_total_matched_bases;
		float match_score = EXON_HALF_MATCH_PERCENTAGE;

		//int before_best_offset = (pos >   halves_record -> best_pos2_list[i])? halves_record -> splice_point_offset_2[i] : halves_record ->splice_point_offset_1[i];

		int b_search_tail = halves_record -> best1_read_end_pos[i] - 5 ;
		int f_search_head = halves_record -> best1_read_start_pos[i] + 5 ;


		if (((ginp2 && (i % 2)) + negative_strand) == 1){
			int tmp_btail = b_search_tail;
			b_search_tail = rl -  f_search_head ;
			f_search_head = rl -  tmp_btail ;
		}

		explain_buff [0] = pos + b_search_tail;
		explain_buff [2] = b_search_tail;

#ifdef TEST_TARGET
		if(memcmp(inb, TEST_TARGET, 15)==0)
			SUBREADprintf ("\n%s %s\nB-SEARCH RANGE [%d - %d] ; POS=%u ; P1_Indel=%d ; P2_Indel = %d\n",nameb, inb ,f_search_head , b_search_tail,pos + b_search_tail ,  halves_record ->indel_in_piece1[i],  halves_record ->indel_in_piece2[i]);
#endif

		unsigned max_jump_length = 999999999;

		junction_tree_b_explorer(bed_table, pos_table, connection_table, inb , b_search_tail, rl , pos + b_search_tail + halves_record ->indel_in_piece1[i] , 0, explain_buff, 0, &result_total_matched_bases , explain_result,  &result_number_pieces, my_value_array_index, ginp, f_search_head, qualityb, EXON_FASTQ_FORMAT, match_score, 0, &max_jump_length);
	//	SUBREADprintf(" POS BBB RL=%d, MATCHED_TOTAL=%d LEN=%d\n", rl, result_total_matched_bases, b_search_tail + halves_record ->indel_in_piece1[i]);

#ifdef DEBUG

		int xx;
		SUBREADprintf("BSEARCH MATCH BASES=%d/%d; PIECES=%d\n", result_total_matched_bases,b_search_tail, result_number_pieces );
		for (xx=0; xx<result_number_pieces;xx++)
		{

			char * chro_name;
			unsigned int chro_pos;

			unsigned int piece_end_pos = explain_result[xx*4];
			int piece_start = explain_result[xx*4+1];
			int piece_end = explain_result[xx*4+2];

			int indels_in_section = (explain_result[4*xx+3] & 0xffff)>>4;
			if(indels_in_section > 0x800)indels_in_section -= 0x1000;
			int indels_point = explain_result[4*xx+3]>>16;

			if(indels_point > f_search_head && indels_point < b_search_tail) indels_in_section=0;


			locate_gene_position(piece_end_pos, &_global_offsets, &chro_name, &chro_pos);
			chro_pos -= piece_end+indels_in_section;

			SUBREADprintf ("BS-RES: [%d ~ %d  %c %s,%d] \n", piece_start, piece_end, ((ginp2 && (i % 2)) + negative_strand == 1)?'~':'@', chro_name , chro_pos);
		}

#endif







		if(result_number_pieces>0)
		{
			total_pieces = result_number_pieces;
			if (result_number_pieces<1)is_settle=0;
			for (j = 0; j < total_pieces/2; j++)
			{
				unsigned int tmp;
				tmp = explain_result[4*(total_pieces - j-1)];
				explain_result[4*(total_pieces - j-1)] = explain_result[4*j];
				explain_result[4*j] = tmp;

				tmp = explain_result[4*(total_pieces - j-1)+1];
				explain_result[4*(total_pieces - j-1)+1] = explain_result[4*j+1];
				explain_result[4*j+1] = tmp;

				tmp = explain_result[4*(total_pieces - j-1)+2];
				explain_result[4*(total_pieces - j-1)+2] = explain_result[4*j+2];
				explain_result[4*j+2] = tmp;

				tmp = explain_result[4*(total_pieces - j-1)+3];
				explain_result[4*(total_pieces - j-1)+3] = explain_result[4*j+3];
				explain_result[4*j+3] = tmp;
			}
			for (j = 0; j < total_pieces; j++)
			{
				int indels_in_section = (explain_result[4*j+3] & 0xffff)>>4;
				if(indels_in_section > 0x800)indels_in_section -= 0x1000;
				int indels_point = explain_result[4*j+3]>>16;
				//			SUBREADprintf("indels_in_section=%d\n",indels_in_section);
				int good_quality = explain_result[4*j+3] & 0x8;
				if(indels_point > f_search_head && indels_point < b_search_tail) indels_in_section=0;
				if(j >0 || good_quality)
					explain_result[4*j] -= (explain_result[4*j+2]-explain_result[4*j+1] + indels_in_section);
				else
					explain_result[4*j] -= (explain_result[4*j+2]-explain_result[4*j+1]);
			}

			if (total_pieces >0)f_search_head = explain_result[(total_pieces-1)*4 + 1];

			explain_buff [0] = explain_result[4*(total_pieces-1)];
			explain_buff [1] = f_search_head;


#ifdef TEST_TARGET
		if(memcmp(inb, TEST_TARGET, 15)==0)
			SUBREADprintf ("\n%s %s\nF-SEARCH RANGE [%d - END] ; POS=%u ; P1_Indel=%d ; P2_Indel = %d\n",nameb, inb ,f_search_head ,explain_result[(total_pieces-1)*4],  halves_record ->indel_in_piece1[i],  halves_record ->indel_in_piece2[i]);
#endif
		}
		{
			result_number_pieces = 0,  result_total_matched_bases = 0;
			max_jump_length = 999999999;
			junction_tree_f_explorer(bed_table, pos_table, connection_table, inb + f_search_head , rl - f_search_head, rl , explain_result[(total_pieces-1)*4], 0, explain_buff, 0, &result_total_matched_bases , explain_result+ 4 * (total_pieces>0?(total_pieces-1):0),  &result_number_pieces, my_value_array_index, ginp, result_head_modify, qualityb+ f_search_head,  EXON_FASTQ_FORMAT, match_score, 0 , &max_jump_length);

			total_total_matched_bases = result_total_matched_bases;

			total_pieces --;
		}

		if (result_number_pieces >0)
			total_pieces += result_number_pieces ;
		else{	
			total_pieces=0;
			is_settle=0;
		}


	//	SUBREADprintf(" POS AAA RL=%d, MATCHED_TOTAL=%d TEST_LEB=%d\n", rl, total_total_matched_bases, rl-f_search_head);

#ifdef DEBUG



		SUBREADprintf ("\nFSEARCH FROM %d MATCH BASES=%d/%d; PIECES=%d\n" , f_search_head, result_total_matched_bases, rl - f_search_head, result_number_pieces );
		for (xx=0; xx<total_pieces;xx++)
		{
			char * chro_name;
			unsigned int chro_pos;

			unsigned int piece_start_pos = explain_result[xx*4];
			int piece_start = explain_result[xx*4+1];
			int piece_end = explain_result[xx*4+2];

			locate_gene_position(piece_start_pos, &_global_offsets, &chro_name, &chro_pos);
			chro_pos -= piece_start ;

			SUBREADprintf ("FS-RES(ALL): [%d ~ %d  %c %s,%d] \n", piece_start, piece_end, ((ginp2 && (i % 2)) + negative_strand == 1)?'~':'@', chro_name , chro_pos);
		}

#endif


		if (total_pieces>4)
			is_settle = 0;

		char cigar_buf[100];

		if(halves_record -> best_vote1_list[i]< TOTAL_SUBREADS * EXON_MAJOR_HALF_VOTES_RATE ){ total_pieces=0 ; is_settle = 0;}

		if (total_pieces>=2)
		{

			if(halves_record -> best_vote2_list[i]<1 || halves_record -> best_vote1_list[i]<2) 
			{
				if (!((halves_record -> half_marks_list[i] & IS_RECOVERED_JUNCTION_READ_STEP4) || (halves_record -> half_marks_list[i] & IS_RECOVERED_JUNCTION_READ)))
				{
					if (tolerable_scan)
						halves_record -> half_marks_list[i] |= IS_RECOVERED_JUNCTION_READ_STEP4;
					else
						halves_record -> half_marks_list[i] |= IS_RECOVERED_JUNCTION_READ;
				}
				halves_record -> best_vote2_list[i]= max(MIN_VOTE2_TMP_VAR, halves_record -> best_vote2_list[i]);
				halves_record -> best_vote1_list[i]= max(9, halves_record -> best_vote1_list[i]);

				halves_record -> best_pos1_list[i] = explain_result[0];
				halves_record -> best_pos2_list[i] = explain_result[0];
				halves_record -> splice_point_list[i] = 0;
			}

			cigar_buf[0]=0;
			if (is_settle)
			{
				int last_indel = 0;
				for (j = 0; j < total_pieces; j++)
				{
					char cigar_piece[100];

					int indel_pos, indel_length;
					int is_good_piece = explain_result[ 4*j+3] & 0x8;
					indel_length = (explain_result[ 4*j+3] & 0xffff)>>4;
					if(indel_length > 0x800)indel_length -= 0x1000; 					
					indel_pos = explain_result[ 4*j+3]>>16;

					if(indel_length && is_good_piece)
					{
						int len_p1 =  indel_pos - explain_result[ 4*j+1] , len_p2 = explain_result[ 4*j+2] - indel_pos +min(0, indel_length);

						if( ( len_p1>0)  && (len_p2 > 0))
							sprintf (cigar_piece, "%dM%d%c%dM", len_p1,  abs(indel_length), indel_length>0?'D':'I', len_p2);
						else
						{
							halves_record -> best_vote1_list[i]=0;
							halves_record -> best_pos1_list[i]=0;
							halves_record -> best_vote2_list[i]=0;
							//SUBREADprintf("READ KILLED: %s\n", nameb);
						}
						//if(j == total_pieces -1 && halves_record -> best_pos2_list[i] == 0xffffffff)
						//	halves_record -> best_pos2_list[i] = explain_result[4*j] + len_p1 + len_p2 + indel_length;
					}
					else{
						if(explain_result[ 4*j+2]- explain_result[ 4*j+1] >0)
							sprintf (cigar_piece, "%dM", explain_result[ 4*j+2]- explain_result[ 4*j+1] );
						else
						{
							halves_record -> best_vote1_list[i]=0;
							halves_record -> best_vote2_list[i]=0;
#ifdef TEST_TARGET
						if(memcmp(inb, TEST_TARGET, 15)==0)
							SUBREADprintf("READ KILLED: %s\n", nameb);
#endif

						}
					}

					strcat(cigar_buf, cigar_piece);
					if (j<total_pieces-1)
					{
#ifdef TEST_TARGET
						if(memcmp(inb, TEST_TARGET, 15)==0)
							SUBREADprintf("CIGAR- %dN p1=%u p2=%u pos_def=%d, half1=%d, indel=%d\n", explain_result[4*j+4]- (explain_result[ 4*j] + explain_result[ 4*j+2] - explain_result[4*j+1] + indel_length), explain_result[4*j+4]  , explain_result[4*j]  , explain_result[4*j+4] - explain_result[4*j],  explain_result[ 4*j+2] - explain_result[4*j+1] + indel_length,  indel_length);
#endif
						sprintf (cigar_piece, "%dN",explain_result[4*j+4]- (explain_result[ 4*j] + explain_result[ 4*j+2] - explain_result[4*j+1] + indel_length));
						strcat(cigar_buf, cigar_piece);
						last_indel = indel_length;

					}
				}
			}

			int mismatch = 0;
			int is_safeguarded= 0;
			compress_cigar(cigar_buf, rl, inb, NULL, NULL, &is_safeguarded);

			halves_record -> final_quality [i] = final_mapping_quality(my_value_array_index, explain_result[0], inb, qualityb[0]?qualityb:NULL, cigar_buf, EXON_FASTQ_FORMAT, &mismatch, rl, 0);


			#ifdef QUALITY_KILL
			//SUBREADprintf("QUAL=%.5f; MM=%d\n", halves_record -> final_quality [i], mismatch);
			if(mismatch > (200-QUALITY_KILL)/2*rl/100)
			{
				halves_record -> best_vote1_list[i] = 0;
				//halves_record -> best_pos1_list[i]=0;
				cigar_buf[0]=0;
				
			}
			#endif

			strncpy(halves_record -> cigar_string_buffer + i * EXON_MAX_CIGAR_LEN, cigar_buf, EXON_MAX_CIGAR_LEN-1);
			if (tolerable_scan
				#ifdef QUALITY_KILL
				 && mismatch <= (200-QUALITY_KILL)/2*rl/100
				#endif
			)
			{
				for(j=0 ; j < total_pieces -1; j++)
				{
					paired_exon_key search_key;

					int this_step_indels = (explain_result[4*j+3]&0xffff)>>4;
					if(this_step_indels > 0x800) this_step_indels -= 0x1000;
					int is_good_piece = (explain_result[4*j+3]&0xf)>>3;

					if(!is_good_piece)continue;

					unsigned int pos_small = min(explain_result[4*j] + explain_result[4*j+2] - explain_result[4*j+1] + this_step_indels, explain_result[4*j+4]);
					unsigned int pos_big   = max(explain_result[4*j] + explain_result[4*j+2] - explain_result[4*j+1] + this_step_indels, explain_result[4*j+4]);
					search_key.small_key = pos_small;
					search_key.big_key = pos_big;


					exon_junction_t *search_res = (exon_junction_t*) HashTableGet(bed_table, &search_key);
					//assert(search_res);
					if(search_res)
					{
 
						if(halves_record -> best_vote1_list[i])
						{
							short left_len = explain_result[4*j+2] - explain_result[4*j+1];
							short right_len = explain_result[4*j+6] - explain_result[4*j+5];

							assert(left_len<MAX_READ_LENGTH);
							assert(right_len<MAX_READ_LENGTH);

							search_res -> supporting_reads++;
							search_res -> left_extend = max(left_len, search_res -> left_extend);
							search_res -> right_extend = max(right_len, search_res -> right_extend);

							assert(search_res -> left_extend <MAX_READ_LENGTH);
							assert(search_res -> right_extend <MAX_READ_LENGTH);
						}

						long long int d1=search_key.small_key, d2=search_key.big_key;
						d1 -=  min(halves_record -> best_pos1_list[i],halves_record -> best_pos2_list[i]);
						d2 -=  max(halves_record -> best_pos1_list[i],halves_record -> best_pos2_list[i]);
						//if(abs(d1)< 2000 && abs(d2)<2000)
						//	(*search_res)|= 0x80000000;
	
#ifdef DEBUG
						SUBREADprintf("FINAL RECOVERED-IS-FINISH\n");
#endif
					}

				}
				
				if(is_settle)
				{
					halves_record -> best_pos1_list[i] = explain_result[0];
					halves_record -> best_pos2_list[i] = explain_result[0];
				}
			}
		}
		else if(is_settle)
		{
			halves_record -> best_vote1_list[i] = max(9, halves_record -> best_vote1_list[i]);
			halves_record -> best_vote2_list[i] = 0;
			halves_record -> best_pos1_list[i] = explain_result[0];

#ifdef TEST_TARGET
			if(memcmp(inb, TEST_TARGET, 15)==0)
				SUBREADprintf("READ KILLED2: %s\n", nameb);
#endif

			if (!((halves_record -> half_marks_list[i] & IS_RECOVERED_JUNCTION_READ_STEP4) || (halves_record -> half_marks_list[i] & IS_RECOVERED_JUNCTION_READ)))
			{
				if (tolerable_scan)
					halves_record -> half_marks_list[i] |= IS_RECOVERED_JUNCTION_READ_STEP4;
				else
					halves_record -> half_marks_list[i] |= IS_RECOVERED_JUNCTION_READ;
			}
		}
		else
		{
#ifdef TEST_TARGET
			if(memcmp(inb, TEST_TARGET, 15)==0)
				SUBREADprintf("READ KILLED3: %s\n", nameb);
#endif
			halves_record -> best_vote2_list[i]=0;
		}

		if (is_settle)
			halves_record -> half_marks_list[i] |= IS_FINALISED_PROCESSING;
		i+=1;
	}

	if(thread_id==0)
		SUBREADprintf("\n");
}



void * explorer_junc_exonbed_wrapper(void * parameters)
{
	struct feed_exonbed_parameters * fep = parameters;
	int thread_id = fep -> thread_id;
	HashTable * bed_table2 = fep->bed_table;

	gene_input_t * ginp = fep->ginp, *ginp2 = fep->ginp2;
	long long int * ginp1_end_pos = fep->ginp1_end_pos , *ginp2_end_pos = fep->ginp2_end_pos; 

	pthread_spin_unlock(fep -> init_lock);
	explorer_junc_exonbed(bed_table2 , fep -> pos_table, fep -> connection_table, fep->halves_record, ginp, ginp2, fep->out_fp, fep->index_prefix, fep->processed_reads, fep->all_processed_reads, fep-> succeed_reads, fep->my_value_array_index_set, fep->all_tables, fep->tolerable_scan, thread_id, ginp1_end_pos, ginp2_end_pos);


	return NULL;

}


// This function adds items that are in NEW but not in OLD into OLD.
// If an item is in both NEW and OLD, its copy in NEW will be deallocated.
// It compares pos_small and pos_big to decide if an item is in NEW or OLD.
void merge_bed_table(HashTable *old, HashTable* new)
{
	int bucket;
	KeyValuePair * cursor;

	for(bucket=0; bucket<new -> numOfBuckets; bucket++)
	{
		cursor = new -> bucketArray[bucket];
		while (1)
		{
			if(!cursor) break;
			paired_exon_key * p = (paired_exon_key * ) cursor -> key;
			exon_junction_t *counter = (exon_junction_t*) cursor ->value;
			exon_junction_t *oldcounter = HashTableGet(old, p);

			if(oldcounter)
			{
				oldcounter->feed_supporting_reads += counter->feed_supporting_reads;
				free(p);
				free(counter);
			}
			else
			{
				HashTablePut(old, p , counter);
			}
			cursor = cursor->next;
		}
	}
}

// This function easily add non-existing positions in NEW.
void merge_pos_table(HashTable *old, HashTable* new)
{
	int bucket;
	KeyValuePair * cursor;

	for(bucket=0; bucket<new -> numOfBuckets; bucket++)
	{
		cursor = new -> bucketArray[bucket];
		while (1)
		{
			if(!cursor) break;
			unsigned int pos = (unsigned int)(cursor -> key - NULL);

			if(!HashTableContainsKey(old, NULL+pos))
				HashTablePut(old, NULL+pos , NULL+pos);

			cursor = cursor->next;
		}
	}
}


// This function first add non-existing postions into NEW
// If a position in OLD is also in NEW, it adds non-existing items of that postion into NEW
// It then deallocate 
void merge_connection_table(HashTable *old, HashTable* new)
{
	int bucket;
	KeyValuePair * cursor;
	for(bucket=0; bucket<new -> numOfBuckets; bucket++)
	{
		cursor = new -> bucketArray[bucket];
		while (1)
		{
			if(!cursor) break;
			unsigned int pos = (unsigned int)(cursor -> key - NULL);

			connect_to_t * connect_old = HashTableGet(old, NULL+pos);

			if(connect_old)
			{
				connect_to_t * connect_new = cursor -> value;
				int i;

				for(i=0; i< MAX_EXON_CONNECTIONS; i++)
				{
					unsigned int connect_to_pos = connect_new -> connect_to[i];
					if(!connect_to_pos) break;
					int j;
					int found = 0;
					for(j=0; j< MAX_EXON_CONNECTIONS; j++)
					{
						if(!connect_old -> connect_to[j]) break;
						if(connect_old -> connect_to[j] == connect_to_pos)
						{
							found =1;
							break;
						}
					}

					if(!found && j<MAX_EXON_CONNECTIONS)
					{
						connect_old -> connect_to[j] = connect_to_pos;
						connect_old -> is_opposite_reversed[j] = connect_new -> is_opposite_reversed[i];
	
						if(j+1<MAX_EXON_CONNECTIONS)
							connect_old -> connect_to[j+1] = 0;
					}
				}
				free(connect_new);
			}
			else
			{
				HashTablePut(old, NULL+pos, cursor-> value);
			}

			cursor = cursor->next;
		}
	}
}



// This function will add supporting counts to ALL, and deallocate the memory space for counters in NEW, then destory NEW.
void add_bed_table(HashTable * all , HashTable * new)
{
	
	int bucket;
	KeyValuePair * cursor;
	for(bucket=0; bucket<new -> numOfBuckets; bucket++)
	{
		cursor = new -> bucketArray[bucket];
		while (1)
		{
			if(!cursor) break;
			paired_exon_key * p = (paired_exon_key * ) cursor -> key;
			exon_junction_t * c = (exon_junction_t *) cursor -> value;

			exon_junction_t * counter = (exon_junction_t*) HashTableGet(all, p);
			assert(counter);


			counter -> supporting_reads += c->supporting_reads;
			counter -> left_extend = max(counter -> left_extend,c -> left_extend);
			counter -> right_extend = max(counter -> right_extend,c -> right_extend);

			

			free(cursor -> value);
			cursor = cursor->next;
		}
	}
}

#define MAX_TABLE_SIZE_JUNCTIONS 5000000

void remove_neighbours(HashTable * bed_table, HashTable * connection_table, HashTable * pos_table)
{

	paired_exon_key ** to_remove_keys;
	int removed_keys=0, bucket;
	KeyValuePair * cursor;

	to_remove_keys = (paired_exon_key **)malloc(sizeof(paired_exon_key *) * MAX_TABLE_SIZE_JUNCTIONS); 

	for(bucket=0; bucket<bed_table -> numOfBuckets; bucket++)
	{
		cursor = bed_table -> bucketArray[bucket];
		while (1)
		{
			if (!cursor) break;
			paired_exon_key * p = (paired_exon_key * ) cursor -> key;
			exon_junction_t *counter = (exon_junction_t *) cursor ->value;

			if(1)
			{
				int delta;
				int indels;
				int is_delete=0;
				for(indels=-6; indels<=6; indels++)
				{
					for(delta=-11; delta < 12; delta++)
					{
						if(delta==0 && indels==0)continue;
						paired_exon_key nbkey;
						nbkey.big_key = p->big_key + indels + delta;
						nbkey.small_key = p->small_key + delta;
						exon_junction_t * nb = HashTableGet(bed_table, &nbkey);
						if(nb && nb -> feed_supporting_reads > counter -> feed_supporting_reads)
						{
							is_delete = 1;
							break;
						}

					}
					if(is_delete)break;
				}
				if(is_delete)to_remove_keys[removed_keys++] = p;
			}

			cursor = cursor->next;
		}
	}

	SUBREADprintf("\n There are %d low-confidence junction table items.\n", removed_keys);

	int i,j;
	for(i=0; i<removed_keys; i++)
	{
		exon_junction_t * rmd = HashTableGet(bed_table, to_remove_keys[i]);
		HashTableRemove(bed_table, to_remove_keys[i]);

		for(j=0;j<2;j++)
		{
			unsigned int c0 = j? to_remove_keys[i] -> big_key: to_remove_keys[i] -> small_key;
			unsigned int c1 = (!j)? to_remove_keys[i] -> big_key: to_remove_keys[i] -> small_key;

			connect_to_t * next_jump = HashTableGet(connection_table, (NULL+c0));
			
			int k;
			for(k=0; k<MAX_EXON_CONNECTIONS; k++)
			{
				if(!next_jump->connect_to[k])break;
				if(next_jump->connect_to[k] == c1)
				{
					int l;
					next_jump->connect_to[k] =0;
				
					for(l=k+1; l<MAX_EXON_CONNECTIONS; l++)
					{
						if(!next_jump->connect_to[l])break;
						next_jump->connect_to[l-1]=next_jump->connect_to[l];
						next_jump->is_opposite_reversed[l-1]=next_jump->is_opposite_reversed[l];
						next_jump->connect_to[l]=0;
					}
					if(k==0 && l==k+1)
					{
						HashTableRemove(pos_table, (NULL+c0));
						HashTableRemove(connection_table, (NULL+c0));
						free(next_jump);
					}
				}
			}
		}

		free(rmd);
		free(to_remove_keys[i]);
	}
	free(to_remove_keys);
	
}

// This function copies all items in old into the returnned HashTable.
// It does not allocate new memory spaces for keys; it only allocates new memory spaces for the number values.
HashTable * bed_table_copy(HashTable * old)
{

	HashTable * ret = HashTableCreate(399997);

	HashTableSetKeyComparisonFunction(ret, pointercmp_forbed);
	HashTableSetHashFunction(ret, pointerHashFunction_forbed);



	int bucket;
	KeyValuePair * cursor;

	for(bucket=0; bucket<old -> numOfBuckets; bucket++)
	{
		cursor = old -> bucketArray[bucket];
		while (1)
		{
			if(!cursor) break;
			paired_exon_key * p = (paired_exon_key * ) cursor -> key;
			exon_junction_t *counter = (exon_junction_t *) cursor ->value;
			exon_junction_t * new_counter = (exon_junction_t*) malloc(sizeof(exon_junction_t));
			new_counter -> supporting_reads = 0;// counter-> supporting_reads;
			new_counter -> feed_supporting_reads = 0;
			new_counter -> strand = counter -> strand;
			new_counter -> left_extend = 0;
			new_counter -> right_extend = 0;

			HashTablePut(ret, p, new_counter);

			cursor = cursor->next;
		}
	}
	
	return ret;
}

void explorer_junc_exonbed_maybe_threads(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table, halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads, gene_value_index_t ** my_value_array_index_set, int all_tables, int tolerable_scan)
{

	memset(halves_record->cigar_string_buffer,0, halves_record->max_len*EXON_MAX_CIGAR_LEN);
	if(EXON_ALL_THREADS < 2)
		explorer_junc_exonbed(bed_table, pos_table, connection_table, halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, all_processed_reads,  succeed_reads, my_value_array_index_set, all_tables, tolerable_scan, 0,  NULL, NULL);
	else
	{
		struct feed_exonbed_parameters fep;
		pthread_spinlock_t  init_lock;
		pthread_t runners [EXON_ALL_THREADS];
		int i;
		long long int ginp1_end_pos[EXON_ALL_THREADS];
		long long int ginp2_end_pos[EXON_ALL_THREADS];
		gene_input_t ginp1s[EXON_ALL_THREADS];
		gene_input_t ginp2s[EXON_ALL_THREADS];
		

		fep.bed_table=bed_table;
		fep.pos_table=pos_table;
		fep.connection_table=connection_table;
		fep.halves_record=halves_record;
		fep.out_fp=out_fp;
		fep.index_prefix=index_prefix;
		fep.processed_reads=processed_reads;
		fep.all_processed_reads=all_processed_reads;
		fep.succeed_reads=succeed_reads;
		fep.my_value_array_index_set=my_value_array_index_set;
		fep.all_tables=all_tables;
		fep.tolerable_scan=tolerable_scan;
		fep.init_lock = &init_lock;
		fep.ginp2 = NULL;

		pthread_spin_init(&init_lock, PTHREAD_PROCESS_PRIVATE);
		pthread_spin_lock(&init_lock);

		HashTable * hashtables [500];

		for (i=0; i< EXON_ALL_THREADS; i++)
		{
			long long int fposition;
			fep.thread_id=i;
			fep.ginp1_end_pos = ginp1_end_pos + i;
			fep.ginp2_end_pos = ginp2_end_pos + i;

			ginp1s[i].space_type = ginp->space_type;
			ginp1s[i].file_type = ginp->file_type;
			ginp1s[i].input_fp = fopen(ginp->filename,"r");
			fposition = ftello(ginp -> input_fp );
			fseeko(ginp1s[i].input_fp, fposition, SEEK_SET); 
			fep.ginp = ginp1s+i;


			if(ginp2)
				fep.ginp2 = ginp2s+i;

			if(i>0)
			{
				HashTable * bed_table_2 = bed_table_copy(bed_table);	
				fep.bed_table = bed_table_2;
				hashtables[i] = bed_table_2;
			}

			pthread_create(runners+i, NULL,  explorer_junc_exonbed_wrapper, &fep);
			pthread_spin_lock(&init_lock);
		}


		long long int max_fpos1 = -1;
		for (i=0; i< EXON_ALL_THREADS; i++)
		{
			pthread_join(*(runners+i), NULL);


			if(!i) SUBREADprintf("\rWaiting for other threads. This can take several minutes...       \r");

			max_fpos1 = max(ginp1_end_pos[i], max_fpos1);

			fclose(ginp1s[i].input_fp);

			if(i>0)
			{
				add_bed_table(bed_table, hashtables[i]);
				HashTableDestroy(hashtables[i]);
			}
		}
	
		fseeko(ginp->input_fp, max_fpos1, SEEK_SET);

		pthread_spin_destroy(&init_lock);
	}
}



void put_connection_table(HashTable *connection_table, unsigned int p1, unsigned int p2, int is_p1_reversed, int is_p2_reversed)
{
	connect_to_t * connect_to ;
	int i;

	connect_to = (connect_to_t *)HashTableGet(connection_table, NULL + p1);
	if (!connect_to)
	{
		connect_to = (connect_to_t * ) malloc(sizeof(connect_to_t));
		connect_to-> connect_to[0] = 0;
		connect_to-> i_am_reversed = is_p1_reversed;
		HashTablePut(connection_table, NULL + p1, connect_to);
	}

	for (i=0; i<MAX_EXON_CONNECTIONS;)
	{
		if (connect_to ->connect_to[i] == p2)
		{
			i = -1;
			break;
		}
		if (!connect_to-> connect_to[i])break;
		i++;
	}

	if (i<0)return;	// exists
	if (i == MAX_EXON_CONNECTIONS){
		//SUBREADprintf("WARNING: TOO MANY CONNECTIONS.\n");
		return;	// too many connections
	}

	connect_to ->connect_to[i] = p2;
	connect_to ->is_opposite_reversed[i] = is_p2_reversed;

	if (i<MAX_EXON_CONNECTIONS-1) connect_to -> connect_to[1+i]=0;

}


#define SHORT_EXON_MIN_LENGTH 18
#define SHORT_EXON_WINDOW 6 
#define SHORT_EXON_EXTEND 5000



void search_short_exons(char * inb0, char * qualityb0, int rl, HashTable * connection_table, HashTable * bed_table, HashTable * pos_table, unsigned int P1_Pos, unsigned int P2_Pos, short read_coverage_start, short read_coverage_end,  gene_value_index_t *base_index, int space_type, int is_negative, int tolerable_scan)
{
	char inb[1201], qualityb[1201];
	if ( (rl <= EXON_LONG_READ_LENGTH ) && (!EXON_EXTENDING_SCAN)) return;

	strcpy(inb, inb0);
	strcpy(qualityb, qualityb0);
	if (is_negative)
	{
		short tmps;
		tmps = read_coverage_end;
		read_coverage_end = read_coverage_start;
		read_coverage_start = tmps;

		reverse_read(inb, rl, space_type);
		if(qualityb[0])reverse_quality(qualityb, rl);
	}
	unsigned int pos_small=min(P1_Pos, P2_Pos), pos_big = max(P1_Pos, P2_Pos);

	int max_score , test_score;
	unsigned int best_j1_edge=0 , best_j2_edge=0;
	int need_to_test = 0;

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
// SCAN TO THE HEAD  /////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////



	if (read_coverage_start  > SHORT_EXON_MIN_LENGTH)
	{
		max_score = -1;

		int need_check2 = 1;
		if(qualityb[0])
		{
			float head_quality = read_quality_score(qualityb , SHORT_EXON_MIN_LENGTH , EXON_FASTQ_FORMAT); 
			if(head_quality < 6 )
				need_check2 = 0;
		}


		if(need_check2)
			if(SHORT_EXON_MIN_LENGTH *0.6 < match_chro(inb, base_index, pos_small, SHORT_EXON_MIN_LENGTH , 0, space_type))
				need_check2 = 0; 


		if(need_check2)
		{

			int delta_pos, is_indel = 0;
			for(delta_pos=-3; delta_pos <=3; delta_pos ++)
			{
				if(match_chro(inb, base_index, pos_small + delta_pos, SHORT_EXON_MIN_LENGTH , 0, space_type) >= SHORT_EXON_MIN_LENGTH*.7)
				{
					is_indel = 1;
					break;
				}
			}
			// The head of the read is incorrect. Do we need to search a long way?
			// See if there is a donor in the head area.
			int test_donor_pos;
			char cc[3];
			cc[2]=0;

			if(!is_indel)
				for(test_donor_pos = SHORT_EXON_MIN_LENGTH ; test_donor_pos < read_coverage_start ; test_donor_pos ++)
				{
					get_chro_2base(cc, base_index, pos_small + test_donor_pos, 0);
					if(is_donar_chars_part(cc))
					{
						need_to_test = 1;
						break;
					}
				}
		}
	}

	max_score = -999;
	int max_is_GTAG = 0;

	if(need_to_test)
	{
		unsigned int test_end = pos_small - SHORT_EXON_EXTEND;
		if(SHORT_EXON_EXTEND > pos_small) test_end = 0;

		unsigned int new_pos = pos_small-16;
		while(1)
		{
			new_pos = match_chro_range(inb,  base_index, new_pos, 7 , new_pos - test_end , SEARCH_BACK);
			if(new_pos==0xffffffff) break;
			// There is an exact match. See if the donor/receptors are matched.
			// new_pos is the new head position of the read.
			int splice_point;
			for(splice_point = SHORT_EXON_MIN_LENGTH; splice_point < read_coverage_start ; splice_point ++)
			{
				char cc[3];
				cc[2]=0;
				char cc2[3];
				cc2[2]=0;

				get_chro_2base(cc, base_index, pos_small + splice_point -2, 0);
				if(is_donar_chars_part(cc))
				{
					// EXON---|CC2---INTRON---CC|---EXON
					get_chro_2base(cc2, base_index, new_pos + splice_point, 0);
					if(is_donar_chars_part(cc2) && paired_chars_part(cc2 , cc, 0)) 
					{
						int matched_in_exon_old = match_chro(inb + splice_point, base_index, pos_small + splice_point , SHORT_EXON_WINDOW , 0, space_type);
						int matched_in_exon_new = match_chro(inb, base_index, new_pos , splice_point, 0, space_type);

						
						test_score = 1000000+ (matched_in_exon_new )*10000  + matched_in_exon_old * 1000 + new_pos - test_end;
						if(test_score <= max_score) continue;
						max_score = test_score + 39999 ;

						if(matched_in_exon_new < splice_point || matched_in_exon_old < SHORT_EXON_WINDOW ) 
							continue;

						max_is_GTAG = (cc2[0]=='G');
						best_j1_edge = new_pos + splice_point;
						best_j2_edge = pos_small + splice_point;
					}
				}
			}
		}
	}


	if(best_j1_edge>0)
	{

		unsigned int pos_small_x, pos_big_x;
		pos_small_x = get_grouped_position(pos_table, best_j1_edge);
		pos_big_x = get_grouped_position(pos_table, best_j2_edge);

		paired_exon_key search_key; 
		search_key.small_key = pos_small_x;
		search_key.big_key = pos_big_x;

		put_connection_table(connection_table, pos_small_x, pos_big_x, is_negative, is_negative);
		put_connection_table(connection_table, pos_big_x, pos_small_x, is_negative, is_negative);

		//if (tolerable_scan && !(halves_record ->half_marks_list[i] & IS_FINALISED_PROCESSING))
		//	SUBREADprintf("FOUND IN STEP 4: %u %u\n", pos_R1, pos_R2);
		exon_junction_t *search_res = (exon_junction_t*) HashTableGet(bed_table, &search_key);
		if(search_res)
			search_res -> feed_supporting_reads ++;
		else
		{
			search_res = (exon_junction_t *)malloc(sizeof(exon_junction_t));
			search_res -> supporting_reads = 0;
			search_res -> left_extend = 0;
			search_res -> right_extend = 0;
			search_res -> strand = !max_is_GTAG;
			paired_exon_key * new_key = (paired_exon_key *) malloc(sizeof(paired_exon_key));
			new_key->small_key = pos_small_x;
			new_key->big_key = pos_big_x;
			HashTablePut(bed_table, new_key, search_res);
		}
	}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
// SCAN TO THE TAIL  /////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

	need_to_test = 0;
	max_score = -999;


	if (read_coverage_end< rl - SHORT_EXON_MIN_LENGTH)
	{
		int need_check2 = 1;
		if(qualityb[0])
		{
			float head_quality = read_quality_score(qualityb + rl - SHORT_EXON_MIN_LENGTH , SHORT_EXON_MIN_LENGTH , EXON_FASTQ_FORMAT); 
			if(head_quality < 6 )
				need_check2 = 0;
		}


		if(SHORT_EXON_MIN_LENGTH *0.6 < match_chro(inb + rl - SHORT_EXON_MIN_LENGTH, base_index, pos_big + rl - SHORT_EXON_MIN_LENGTH , SHORT_EXON_MIN_LENGTH , 0, space_type))
			need_check2 = 0; 
		if(need_check2)
		{
			int delta_pos, is_indel = 0;
			for(delta_pos=-3; delta_pos <=3; delta_pos ++)
			{
				if(match_chro(inb + rl - SHORT_EXON_MIN_LENGTH, base_index, pos_big + rl - SHORT_EXON_MIN_LENGTH + delta_pos, SHORT_EXON_MIN_LENGTH , 0, space_type) >= SHORT_EXON_MIN_LENGTH*.7)
				{
					is_indel = 1;
					break;
				}
			}
			// The head of the read is incorrect. Do we need to search a long way?
			// See if there is a donor in the head area.
			int test_donor_pos;
			char cc[3];
			cc[2]=0;

			if(!is_indel)
				for(test_donor_pos = read_coverage_end  ; test_donor_pos < rl ; test_donor_pos ++)
				{
					get_chro_2base(cc, base_index, pos_big + test_donor_pos, 0);
					if(is_donar_chars_part(cc))
					{
						need_to_test = 1;
						break;
					}
				}
		}
	}

	best_j1_edge = 0;

	if(need_to_test)
	{
		unsigned int test_end = pos_big + SHORT_EXON_EXTEND;
		if(test_end > base_index -> length + base_index -> start_point) test_end = base_index -> length + base_index -> start_point;

		unsigned int new_pos = pos_big +rl - SHORT_EXON_MIN_LENGTH +16;
		int max_is_GTAG = 0;

		while(1)
		{
			new_pos = match_chro_range(inb + rl - SHORT_EXON_MIN_LENGTH,  base_index, new_pos, 7 , test_end - new_pos , SEARCH_FRONT);
			if(new_pos==0xffffffff) break;
			// There is an exact match. See if the donor/receptors are matched.
			// (new_pos + SHORT_EXON_MIN_LENGTH -rl + splice_point) is the new exon start.

			int splice_point;
			for(splice_point = read_coverage_end ; splice_point < rl -  SHORT_EXON_MIN_LENGTH; splice_point ++)
			{
				char cc[3];
				cc[2]=0;
				char cc2[3];
				cc2[2]=0;

				unsigned int new_pos_tail = (new_pos + SHORT_EXON_MIN_LENGTH -rl + splice_point);

				get_chro_2base(cc, base_index, pos_big + splice_point, 0);
				if(is_donar_chars_part(cc))
				{
					get_chro_2base(cc2, base_index, new_pos_tail -2, 0);
					if(is_donar_chars_part(cc2) && paired_chars_part(cc , cc2, 0)) 
					{
						int matched_in_exon_new = match_chro(inb + splice_point, base_index, new_pos_tail , rl - splice_point , 0, space_type);
						int matched_in_exon_old = match_chro(inb + splice_point - SHORT_EXON_WINDOW , base_index, pos_big + splice_point - SHORT_EXON_WINDOW , SHORT_EXON_WINDOW, 0, space_type);

						test_score = 1000000+ (matched_in_exon_new)*10000 + matched_in_exon_old * 1000  + test_end - new_pos;
						if(test_score <= max_score) continue;
						max_score = test_score + 39999;

						if(matched_in_exon_new < (rl - splice_point) || matched_in_exon_old < SHORT_EXON_WINDOW)
							continue;

						// EXON ---|CC---INTRON---CC2|--- EXON 
						max_is_GTAG = 1==(cc[0]=='G');
						best_j1_edge = pos_big + splice_point;
						best_j2_edge = new_pos_tail;
					}
				}
			}

		}
	}


	if(best_j1_edge>0)
	{

		unsigned int pos_small_x, pos_big_x;
		pos_small_x = get_grouped_position(pos_table, best_j1_edge);
		pos_big_x = get_grouped_position(pos_table, best_j2_edge);

		paired_exon_key search_key; 
		search_key.small_key = pos_small_x;
		search_key.big_key = pos_big_x;

		put_connection_table(connection_table, pos_small_x, pos_big_x, is_negative, is_negative);
		put_connection_table(connection_table, pos_big_x, pos_small_x, is_negative, is_negative);

		//if (tolerable_scan && !(halves_record ->half_marks_list[i] & IS_FINALISED_PROCESSING))
		//	SUBREADprintf("FOUND IN STEP 4: %u %u\n", pos_R1, pos_R2);
		exon_junction_t *search_res = (exon_junction_t*) HashTableGet(bed_table, &search_key);
		if(search_res)
			search_res -> feed_supporting_reads ++;
		else
		{
			search_res = (exon_junction_t *)malloc(sizeof(exon_junction_t));
			search_res -> supporting_reads = 0;
			search_res -> feed_supporting_reads = 0;
			search_res -> left_extend = 0;
			search_res -> right_extend = 0;
			search_res -> strand = !max_is_GTAG;

			paired_exon_key * new_key = (paired_exon_key *) malloc(sizeof(paired_exon_key));
			new_key->small_key = pos_small_x;
			new_key->big_key = pos_big_x;
			HashTablePut(bed_table, new_key, search_res);
		}
	}
}

int is_repeated_region(unsigned char * repeat_recorder, int is_negative, int votes, int coverage_start, int coverage_end, int rl)
{
	int i;
	int offset_shifting = (rl > 220)?4:0;

	int encounter = 0;
	if(votes<1) return 1;

	coverage_start = ((coverage_start  - EXON_MAX_BIGMARGIN_OVERLAPPING) >> offset_shifting) -1;
	coverage_end = ((coverage_end  + EXON_MAX_BIGMARGIN_OVERLAPPING) >> offset_shifting) +1;

	for(i=0; i<48; i+=3)
	{
		if((!repeat_recorder[i+2]) || (repeat_recorder[i+2] &0x80) != 0x80*is_negative) continue;
		if((repeat_recorder[i+2] & 0x7f)  > votes-2)
			if(repeat_recorder[i] >= coverage_start && repeat_recorder[i+1] <= coverage_end)
				encounter ++;
		if(encounter>=2)return 1;
	}
	return 0;
}


void feed_exonbed(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table,  halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads, gene_value_index_t ** my_value_array_index_set, int all_tables, int tolerable_scan, int thread_id, long long int * ginp1_end_pos,  long long int * ginp2_end_pos)
{
	int i, ic=0, j;

	unsigned int my_range_start , my_range_end;
	gene_value_index_t * my_value_array_index;

	if(thread_id==0)
		SUBREADprintf("Detecting junctions:\n");

	double reads_per_thread = processed_reads*(1+(ginp2!=NULL)) * 1.0 / EXON_ALL_THREADS;
	my_range_start = max((int)(reads_per_thread * thread_id - 0.5),0);
	my_range_end = min((int)(reads_per_thread*(1+thread_id) + 0.5), processed_reads*(1+(ginp2!=NULL)));
	
	if(thread_block_locations[thread_id]>=0)
	{
		fseeko(ginp->input_fp, thread_block_locations[thread_id], SEEK_SET);
		i = my_range_start;
	}
	else
	{
		for (i=0; i<my_range_start; i++)
		{
			geinput_jump_read(ginp);
		}

		thread_block_locations[thread_id] = ftello(ginp->input_fp);
	}

	while (1)
	{
		char nameb[1201], inb[1201], qualityb[1201];
		int rl = -1, mapping_quality=0, mapping_flags;

		if (i < my_range_end)
			rl = geinput_next_read_sam(ginp, nameb, inb, qualityb, NULL, NULL, &mapping_quality, &mapping_flags, ginp2 && (i % 2));


		if (rl<0){
			if(ginp1_end_pos)
				*ginp1_end_pos = ftello(ginp->input_fp);
			break;
		}

		/*if (ginp2 && (i % 2))
			reverse_read(inb, rl, ginp->space_type);*/

		if(i % (processed_reads/14) ==0 && i>1 && !thread_id)
			print_text_scrolling_bar("", i* EXON_ALL_THREADS*1./(1+(ginp2!=NULL))/processed_reads, 80, &ic);

		unsigned int pos = ( IS_R1_CLOSE_TO_5 & halves_record -> half_marks_list[i] ) ?halves_record -> best_pos1_list[i]:halves_record -> best_pos2_list[i];
		unsigned int pos2 = ( IS_R1_CLOSE_TO_5 & halves_record -> half_marks_list[i] ) ?halves_record -> best_pos2_list[i]:halves_record -> best_pos1_list[i];



#ifdef TEST_TARGET
		if(memcmp(TEST_TARGET, inb,15)==0)
		{
			SUBREADprintf("RAW i=%d PRE TEST DONOR = %u, %u\n", i, halves_record -> best_pos1_list[i] ,halves_record -> best_pos2_list[i]);
		}
#endif

		char negative_strand = halves_record -> is_reversed_list[i];
		int is_repeated = 0;
		is_repeated = is_repeated_region((unsigned char *)halves_record -> cigar_string_buffer + i *EXON_MAX_CIGAR_LEN, negative_strand?1:0, halves_record -> best_vote1_list[i], halves_record -> read_coverage_start[i], halves_record -> read_coverage_end[i], rl);
		if(is_repeated || mapping_quality>=MAX_QUALITY_TO_CALL_JUNCTION)
		{
			i++;
			continue;
		}


		my_value_array_index = NULL;
		for(j=0; j<all_tables; j++)
		{
			unsigned int this_index_last_base = my_value_array_index_set[j] -> start_base_offset + my_value_array_index_set[j] -> length;
			if(halves_record -> best_vote2_list[i]>=1)
			{
				if(min(pos, pos2) > my_value_array_index_set[j] -> start_base_offset + (j?900000:0) && max(pos, pos2) < this_index_last_base - ((j==all_tables-1)? 0:900000))
				{
					my_value_array_index = my_value_array_index_set[j];
					break;
				}
			}
			else if(halves_record -> best_vote2_list[i]==0)
			{
				if( halves_record -> best_pos1_list[i]  > my_value_array_index_set[j] -> start_base_offset + (j?900000:0) && halves_record ->best_pos1_list[i]  < this_index_last_base - ((j==all_tables-1)? 0:900000))
				{
					my_value_array_index = my_value_array_index_set[j];
					break;
				}

			}
		}

		if(!my_value_array_index)
		{
			i++;
			continue;
		}

		//SUBREADprintf("\n P1=%u  P2=%u \n",  halves_record -> best_pos1_list[i],  halves_record -> best_pos2_list[i]);

		//SUBREADprintf("FEED-P2 v1=%d, v2=%d\n", halves_record -> best_vote1_list[i] ,halves_record -> best_vote2_list[i] );
		long long int dist = halves_record -> best_pos1_list[i] ;
		dist-= halves_record -> best_pos2_list[i];

		if(halves_record -> best_vote1_list[i]>= (int)(TOTAL_SUBREADS* EXON_MAJOR_HALF_VOTES_RATE) && halves_record -> best_vote2_list[i]>=MIN_VOTE2_TMP_VAR )
		{
			//SUBREADprintf("FEED-P3\n");
			//char * chro_name;
			//unsigned int chro_pos, chro_pos_small, chro_pos_large;
			//char is_strong = (min(halves_record ->best_vote1_list[i],halves_record ->best_vote2_list[i])>5);

			search_short_exons(inb, qualityb, rl, connection_table, bed_table, pos_table, halves_record -> best_pos1_list[i], halves_record -> best_pos2_list[i], halves_record -> read_coverage_start[i], halves_record -> read_coverage_end[i], my_value_array_index, ginp->space_type, halves_record -> is_reversed_list[i], tolerable_scan);
			int real_break_point;
			char is_soft_condition = 0;
			int test_range = rl / 4;
			char is_accepted = 0;
			int indel_offset_pos1=0, indel_offset_pos2=0;
			int guess_break_point = (halves_record -> is_reversed_list[i]) ? (rl - halves_record ->splice_point_list[i]) : halves_record ->splice_point_list[i];
			int best_donor_score =-1;
			int best_GTAG= 0;

			if (min(halves_record ->best_vote1_list[i],halves_record ->best_vote2_list[i]) >= MIN_VOTE2_TMP_VAR)// && short_overlap)
			{
				int best_indel_p1 = -9999, best_indel_p2 = -9999;
				int indel_x1;
				for(indel_x1 = 0; indel_x1 < 2*EXON_INDEL_TOLERANCE +1 ; indel_x1 ++)
				{
					int indel_offset1 = ((indel_x1 %2)?1:-1) * ((indel_x1+1)>>1);

					int indel_x2;
					for(indel_x2 = 0; indel_x2 < 2*EXON_INDEL_TOLERANCE +1 ; indel_x2 ++)
					{

						//if(indel_x1>2 || indel_x2>2)continue;

						int indel_offset2 = ((indel_x2 %2)?1:-1) * ((indel_x2+1)>>1);
						int test_donor_score=-1;
						int test_real_break_point;
						int is_GTAG=0;

						is_accepted = test_donor(inb, rl, min(pos, pos2), max(pos,pos2), guess_break_point, negative_strand, test_range, is_soft_condition, EXON_INDEL_TOLERANCE, &test_real_break_point, my_value_array_index, indel_offset1, indel_offset2, negative_strand, ginp->space_type, halves_record ->best_vote2_list[i], &test_donor_score, &is_GTAG);

						test_donor_score -= abs(indel_offset2)*1 + abs(indel_offset1)*1;

						if (is_accepted  && (test_donor_score > best_donor_score)){
							//if(best_donor_score >0)
							//	SUBREADprintf("TEST SCORE=%d MAX=%d O1=%d O2=%d\n", test_donor_score , best_donor_score , indel_offset1, indel_offset2);
							best_indel_p1 = indel_offset1;
							best_indel_p2 = indel_offset2;
							best_donor_score = test_donor_score;
							real_break_point = test_real_break_point;
							best_GTAG = is_GTAG;
						}
					}
				}

				if(best_donor_score >0)
				{
					//SUBREADprintf("FINAL BEST SCORE=%d\n", best_donor_score);
					indel_offset_pos1 = best_indel_p1;
					indel_offset_pos2 = best_indel_p2;
					halves_record ->splice_point_offset_1[i] = best_indel_p1;
					halves_record ->splice_point_offset_2[i] = best_indel_p2;
					is_accepted =1;
				}
			}
			else
			{
				is_accepted = test_donor(inb, rl, min(pos, pos2), max(pos,pos2), guess_break_point, negative_strand, test_range, is_soft_condition, EXON_INDEL_TOLERANCE, &real_break_point, my_value_array_index,0,0, negative_strand, ginp->space_type, halves_record ->best_vote2_list[i], &best_donor_score, & best_GTAG);


				halves_record ->splice_point_offset_2[i]=0;
				halves_record ->splice_point_offset_1[i]=0;
			}

			#ifdef QUALITY_KILL

			/*
			if(is_accepted)
			{
				char cigar[100];

				if( negative_strand )
					reverse_read(inb, rl, ginp->space_type);

				sprintf(cigar, "%dM%dN%dM", real_break_point, max(pos,pos2) -  min(pos, pos2) , rl-real_break_point);
				int qq = final_mapping_quality(my_value_array_index, min(pos, pos2), inb, NULL, cigar,0);
				if (qq<QUALITY_KILL) is_accepted=0;
			}*/
			#endif

			//if ((i % 2) && is_accepted)printf ("22 SECOND!\n");

			//SUBREADprintf ("\nPRE-ACCEPT=%d\n", is_accepted);

			if((! is_accepted) && EXON_FUSION_DETECTION && (halves_record -> half_marks_list[i] & IS_FUSION))
			{
				int p1_offset = (halves_record ->half_marks_list[i] & IS_R1_CLOSE_TO_5)? 20:0; 
				int p2_offset = (halves_record ->half_marks_list[i] & IS_R1_CLOSE_TO_5)? 0:20;
				int c1_offset = (!(halves_record ->half_marks_list[i] & IS_R1_CLOSE_TO_5) ^ !(halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1))? 20:0; 
				int c2_offset = (!(halves_record ->half_marks_list[i] & IS_R1_CLOSE_TO_5) ^ !(halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R2))? 0:20; 
				int m1 = match_chro(inb+p1_offset, my_value_array_index, halves_record -> best_pos1_list[i]+c1_offset , rl -20, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1)?1:0, ginp->space_type);
				int m2 = match_chro(inb+p2_offset, my_value_array_index, halves_record -> best_pos2_list[i]+c2_offset , rl -20, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R2)?1:0, ginp->space_type);
				//SUBREADprintf ("M1=%d   M2=%d    H=%f\n", m1, m2, rl * .85);
				if ((m1 < (rl-20) *.85) && (m2 < (rl-20) *.85))
				{
					is_accepted=1;
					real_break_point = guess_break_point;
				}
			}


			if (is_accepted)
			{
				// real_break_point is on chromosome: pos+real_break_point / pos2 + real_break_point are end/start of exons
				unsigned int pos_R1 = real_break_point+ halves_record -> best_pos1_list[i] - ((  halves_record -> best_pos1_list[i] < halves_record -> best_pos2_list[i]  )?indel_offset_pos1:indel_offset_pos2);
				unsigned int pos_R2 = real_break_point+ halves_record -> best_pos2_list[i] - ((  halves_record -> best_pos1_list[i] < halves_record -> best_pos2_list[i]  )?indel_offset_pos2:indel_offset_pos1);

				//SUBREADprintf("ACCP: %u, %u / indel = %d, %d?\n", pos_R1, pos_R2, indel_offset_pos1, indel_offset_pos2);

				halves_record -> splice_point_list[i] = real_break_point;

				pos_R1 = get_grouped_position(pos_table, pos_R1);
				pos_R2 = get_grouped_position(pos_table, pos_R2);

				unsigned int pos_small = min(pos_R1, pos_R2);
				unsigned int pos_big   = max(pos_R1, pos_R2);

#ifdef TEST_TARGET
				if(memcmp(TEST_TARGET, inb, 15)==0) 
					SUBREADprintf("ACCEPTED: %u, %u\nRAW: %u, %u\nBP:%d\tOFFSETS=%d %d\n", pos_small, pos_big, halves_record -> best_pos1_list[i] , halves_record -> best_pos2_list[i], real_break_point, indel_offset_pos1, indel_offset_pos2);
#endif

				paired_exon_key search_key; 
				search_key.small_key = pos_small;
				search_key.big_key = pos_big;
				put_connection_table(connection_table, pos_R1, pos_R2, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1)?1:0, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R2)?1:0);
				put_connection_table(connection_table, pos_R2, pos_R1, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R2)?1:0, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1)?1:0);

				//if (tolerable_scan && !(halves_record ->half_marks_list[i] & IS_FINALISED_PROCESSING))
				//	SUBREADprintf("FOUND IN STEP 4: %u %u\n", pos_R1, pos_R2);

				exon_junction_t *search_res = (exon_junction_t*) HashTableGet(bed_table, &search_key);
				if(search_res)
					search_res -> feed_supporting_reads ++;
				else
				{
					search_res = (exon_junction_t *)malloc(sizeof(exon_junction_t));
					search_res -> supporting_reads = 0;// (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1)? 0x80000000 : 0;
					search_res -> feed_supporting_reads =0;
					search_res -> strand = !best_GTAG;// 0 : positive; 1: negative
					search_res -> left_extend = 0;
					search_res -> right_extend = 0;
					paired_exon_key * new_key = (paired_exon_key *) malloc(sizeof(paired_exon_key));
					new_key->small_key = pos_small;
					new_key->big_key = pos_big;
					//SUBREADprintf("INSERT0: POS=%u , %u\n",  pos_small, pos_big);
					HashTablePut(bed_table, new_key, search_res);
				}

				//(*search_res)++;

			}
			else
			{
				halves_record -> best_vote2_list [i]=0;
			}
		}
		else if (halves_record -> best_vote1_list[i]>= (int)(EXON_MAJOR_HALF_VOTES_RATE * TOTAL_SUBREADS/* (rl -15)/EXON_SUBREAD_GAP*/))
		{
			search_short_exons(inb, qualityb, rl, connection_table, bed_table, pos_table, halves_record -> best_pos1_list[i] , halves_record -> best_pos1_list[i], halves_record -> read_coverage_start[i], halves_record -> read_coverage_end[i], my_value_array_index, ginp->space_type, halves_record -> is_reversed_list[i], tolerable_scan);
		}

		i++;
	}
	if(thread_id==0)
		SUBREADprintf("\n");
}




void * feed_exonbed_wrapper(void * parameters)
{
	struct feed_exonbed_parameters * fep = parameters;
	int thread_id = fep -> thread_id;

	HashTable * bed_table2 = fep-> bed_table;
	HashTable * pos_table2 = fep-> pos_table;
	HashTable * connection_table2 = fep-> connection_table;


	gene_input_t * ginp = fep->ginp, *ginp2 = fep->ginp2;
	long long int * ginp1_end_pos = fep->ginp1_end_pos , *ginp2_end_pos = fep->ginp2_end_pos; 


	pthread_spin_unlock(fep -> init_lock);
	feed_exonbed(bed_table2, pos_table2, connection_table2, fep->halves_record, ginp,ginp2, fep->out_fp, fep->index_prefix, fep->processed_reads, fep->all_processed_reads, fep-> succeed_reads, fep->my_value_array_index_set,  fep->all_tables, fep->tolerable_scan, thread_id, ginp1_end_pos, ginp2_end_pos);


	return NULL;
}


void feed_exonbed_maybe_threads(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table,  halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads, gene_value_index_t ** my_value_array_index_set, int all_tables, int tolerable_scan)
{
	if(EXON_ALL_THREADS < 2)
		feed_exonbed(bed_table, pos_table, connection_table, halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, all_processed_reads,  succeed_reads, my_value_array_index_set, all_tables, tolerable_scan, 0, NULL , NULL);
	else
	{
		struct feed_exonbed_parameters fep;
		pthread_spinlock_t  init_lock;
		pthread_t runners [EXON_ALL_THREADS];
		int i;
		long long int ginp1_end_pos[EXON_ALL_THREADS];
		long long int ginp2_end_pos[EXON_ALL_THREADS];
		gene_input_t ginp1s[EXON_ALL_THREADS];
		gene_input_t ginp2s[EXON_ALL_THREADS];
	

		fep.bed_table=bed_table;
		fep.pos_table=pos_table;
		fep.connection_table=connection_table;
		fep.halves_record=halves_record;
		fep.out_fp=out_fp;
		fep.index_prefix=index_prefix;
		fep.processed_reads=processed_reads;
		fep.all_processed_reads=all_processed_reads;
		fep.succeed_reads=succeed_reads;
		fep.my_value_array_index_set=my_value_array_index_set;
		fep.all_tables=all_tables;
		fep.tolerable_scan=tolerable_scan;
		fep.init_lock = &init_lock;
		fep.ginp2=NULL;

		HashTable * hashtables [3][500];

		pthread_spin_init(&init_lock, PTHREAD_PROCESS_PRIVATE);
		pthread_spin_lock(&init_lock);
		for (i=0; i< EXON_ALL_THREADS; i++)
		{
			long long int fposition;
			fep.thread_id=i;
			fep.ginp1_end_pos = ginp1_end_pos + i;
			fep.ginp2_end_pos = ginp2_end_pos + i;

			ginp1s[i].space_type = ginp->space_type;
			ginp1s[i].file_type = ginp->file_type;
			ginp1s[i].input_fp = fopen(ginp->filename,"r");
			fposition = ftello(ginp -> input_fp );
			fseeko(ginp1s[i].input_fp, fposition, SEEK_SET); 
			fep.ginp = ginp1s+i;

			if(ginp2)
			{
				fep.ginp2 = ginp2s+i;
			}

			if(i>0)
			{
				HashTable * bed_table_2 = HashTableCreate(399997);
				HashTable * pos_table_2 = HashTableCreate(399997);
				HashTable * connection_table_2 = HashTableCreate(399997);


				HashTableSetKeyComparisonFunction(connection_table_2, pointercmp_forpos);
				HashTableSetHashFunction(connection_table_2, pointerHashFunction_forpos);

				HashTableSetKeyComparisonFunction(pos_table_2, pointercmp_forpos);
				HashTableSetHashFunction(pos_table_2, pointerHashFunction_forpos);

				
				HashTableSetKeyComparisonFunction(bed_table_2, pointercmp_forbed);
				HashTableSetHashFunction(bed_table_2, pointerHashFunction_forbed);


				fep.bed_table=bed_table_2;
				fep.pos_table=pos_table_2;
				fep.connection_table=connection_table_2;

				hashtables[0][i] = bed_table_2;
				hashtables[1][i] = pos_table_2;
				hashtables[2][i] = connection_table_2;
			}

			pthread_create(runners+i, NULL, feed_exonbed_wrapper, &fep);
			pthread_spin_lock(&init_lock);
		}
		long long int max_fpos1 = -1;

		for (i=0; i< EXON_ALL_THREADS; i++)
		{
			pthread_join(*(runners+i), NULL);


			max_fpos1 = max(ginp1_end_pos[i], max_fpos1);

			fclose(ginp1s[i].input_fp);


			if(i>0)
			{
				merge_bed_table(bed_table , hashtables[0][i]);
				merge_pos_table(pos_table , hashtables[1][i]);
				merge_connection_table(connection_table , hashtables[2][i]);
				HashTableDestroy(hashtables[0][i]);
				HashTableDestroy(hashtables[1][i]);
				HashTableDestroy(hashtables[2][i]);
			}
		}
	
		fseeko(ginp->input_fp, max_fpos1, SEEK_SET);

		pthread_spin_destroy(&init_lock);
	}

}



void print_exon_res_single(gene_value_index_t *array_index , halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads)
{
	int i=0, ic=0 ;
	int mapping_quality=0, mapping_flags=0;


	SUBREADprintf("%u reads were processed. Saving the results for them:\n", processed_reads);
	if(ftello(out_fp)<1)
	{
		unsigned int last_offset = 0;
		i=0;
		while(_global_offsets.read_offsets[i])
		{
			fprintf(out_fp, "@SQ\tSN:%s\tLN:%u\n", _global_offsets.read_names + i*MAX_READ_NAME_LEN, _global_offsets.read_offsets[i] - last_offset+16);
			last_offset = _global_offsets.read_offsets[i];
			i++;
		}
	}

	i=0;

	while (1)
	{
		char nameb[1201], inb[1201], qualityb[1201];
		char old_line[3001];
		int rl;
		if(i >= processed_reads*(1+(ginp2!=NULL)))break;


		if(i % (processed_reads/14) ==0 && i>1)
			print_text_scrolling_bar("", i*1./(1+(ginp2!=NULL))/processed_reads, 80, &ic);

		old_line[0]=0;

		if (ginp2 && (i % 2))
		{
			if(ginp->file_type >= GENE_INPUT_SAM_SINGLE)
				geinput_readline_back(ginp, old_line);

			rl = geinput_next_read_sam(ginp, nameb, inb, qualityb, NULL, NULL, &mapping_quality, &mapping_flags, !halves_record -> is_reversed_list[i]);

		}
		else
		{
			if(ginp->file_type >= GENE_INPUT_SAM_SINGLE)
				geinput_readline_back(ginp, old_line);
			rl = geinput_next_read_sam(ginp, nameb, inb, qualityb, NULL, NULL, &mapping_quality, &mapping_flags, halves_record -> is_reversed_list[i]);

			/*
			if(halves_record -> is_reversed_list[i])
			{
				int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb[0]>='A' && inb[0]<='Z')?1:0;
				reverse_read(inb+ rev_offset, rl- rev_offset, ginp->space_type);
				if(qualityb[0])
					reverse_quality(qualityb, rl- rev_offset);
			}*/

		}



		if (rl<0){
			break;
		}

		if(!qualityb[0])
		{
			int xx;
			EXON_FASTQ_FORMAT = FASTQ_PHRED33;
			for (xx=0; xx<rl; xx++) qualityb[xx]='J';
			qualityb[xx]='\0';
		}

		if(EXON_FASTQ_FORMAT == FASTQ_PHRED64)
			fastq_64_to_33(qualityb);

		if (ginp->file_type == GENE_INPUT_PLAIN)
			sprintf (nameb, "Raw-%d", ginp2?(i/2):i);






		if((halves_record -> best_vote1_list[i]>=1 && halves_record -> best_vote2_list[i]>=1 ))
		{
			char /** chro_name, * chro_name2,*/ *chro_name_min;
			unsigned int chro_pos;//, chro_pos_small, chro_pos_large;
			char cigar_buf[100], *cigar_print;

			unsigned int pos = halves_record -> best_pos1_list[i];

			locate_gene_position(pos, &_global_offsets, &chro_name_min, &chro_pos);

			gene_quality_score_t mapping_quality;

			if (halves_record -> cigar_string_buffer[i * EXON_MAX_CIGAR_LEN])
			{
				cigar_print = halves_record -> cigar_string_buffer + i * EXON_MAX_CIGAR_LEN;
			}
			else
			{
				sprintf (cigar_buf, "%dM",  rl);
				cigar_print = cigar_buf;
			}
			mapping_quality = halves_record -> final_quality[i] ;


			fprintf(out_fp,"%s\t%d\t%s\t%u\t%d\t%s\t*\t0\t0\t%s\t%s",  nameb, halves_record -> is_reversed_list[i]?16:0 , chro_name_min, chro_pos+1, (int)(0.5+mapping_quality) , cigar_print , inb, qualityb);

			fprintf(out_fp,"\n");

		}
		else if(!EXON_JUNCTION_READS_ONLY)
		{
			#ifdef QUALITY_KILL_SUBREAD
			{
				int k=0, field=0, old_quality=0;
				//char old_cigar[100];
				unsigned int is_ID =0;
				char cc;
				while( (cc = old_line[k]) )
				{
					if(cc=='\t')
					{
						field++;
						k++;
						continue;
					}
					if(field == 4)
						old_quality = old_quality * 10 + (cc-'0');
					else if(field == 5)
					{
						if(cc=='I'||cc=='D') is_ID = 1;
					}
					k++;

				}
				
				if(old_quality<=QUALITY_KILL_SUBREAD + (is_ID?0:0))
					fprintf(out_fp,"%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",  nameb,inb, qualityb);
				else
					fprintf(out_fp, "%s\n",old_line);
			}
			#else
			fprintf(out_fp, "%s\n",old_line);
			#endif	

		}
		i++;

	}
	SUBREADprintf("\n");
}





void print_exon_res_paired(gene_value_index_t *array_index , halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads)
{
	int i=0, ic=0 , j;

	if(ftello(out_fp)<1)
	{
		unsigned int last_offset = 0;
		i=0;
		while(_global_offsets.read_offsets[i])
		{
			fprintf(out_fp, "@SQ\tSN:%s\tLN:%u\n", _global_offsets.read_names + i*MAX_READ_NAME_LEN, _global_offsets.read_offsets[i] - last_offset+16);
			last_offset = _global_offsets.read_offsets[i];
			i++;
		}
	}

	i=0;

	SUBREADprintf("%u fragments were processed. Saving the results for them:\n", processed_reads);

	while (1)
	{
		char nameb1[1201],nameb2[1201], inb1[1201], qualityb1[1201], inb2[1201], qualityb2[1201];
		char old_line1[3001],old_line2[3001];
		char old_chro1[120], old_chro2[120], old_cigar1[100], old_cigar2[100];
		unsigned int old_pos1=0, old_pos2=0;
		unsigned int old_flag1=0, old_flag2=0;
		unsigned int old_quality1=0, old_quality2=0;
		char old_inb1_reversed = 0;
		char old_inb2_reversed = 0;

		old_chro1 [0] = 0;
		old_chro2 [0] = 0;

		int rl1=0, rl2=0;
		if(i >= processed_reads*(1+(ginp2!=NULL)))break;


		if(i % 10000 ==0 && i>1)
			print_text_scrolling_bar("Saving results", i*1./(1+(ginp2!=NULL))/processed_reads, 80, &ic);

		old_line1[0]=0;
		old_line2[0]=0;

		if(ginp->file_type >= GENE_INPUT_SAM_SINGLE)
		{
			//long long int file_pos = ftello(ginp -> input_fp);
			geinput_readline(ginp, old_line1,0);
			geinput_readline(ginp, old_line2,0);
			//fseeko(ginp -> input_fp , file_pos , SEEK_SET);

			for(j=0; j<2; j++)
			{
				char * old_line = j?old_line2: old_line1, cc;
				char * old_chro = j?old_chro2:old_chro1;
				char * old_cigar = j?old_cigar2:old_cigar1;
				char * nameb = j?nameb2:nameb1;
				char * inb = j?inb2:inb1;
				char * qualityb = j?qualityb2:qualityb1;
				int k=0, field=0, ci=0,old_flag=0, old_quality=0;
				unsigned int old_pos = 0;
				while( (cc = old_line[k]) )
				{
					if(cc=='\t')
					{
						field++;
						k++;
						if(field == 1)nameb[ci]=0;
						else if(field == 3)old_chro[ci]=0;
						else if(field == 6)old_cigar[ci]=0;
						else if(field == 10)
						{
							inb[ci]=0;
							if(j) rl2=ci; else rl1 =ci;
						}
						else if(field == 11)qualityb[ci]=0;
						ci=0;
						continue;
					}
					if(field == 9)
						inb[ci++] = cc;
					else if(field == 10)
						qualityb[ci++] = cc;
					else if(field == 0)
						nameb[ci++] = cc;
					else if(field == 1)
						old_flag = old_flag*10 + (cc-'0');
					else if(field == 2)
						old_chro[ci++] = cc;
					else if(field == 3)
						old_pos = old_pos * 10 + (cc-'0');
					else if(field == 4)
						old_quality = old_quality * 10 + (cc-'0');
					else if(field == 5)
						old_cigar[ci++] = cc;

					k++;

				}
				if(ci>0)qualityb[ci]=0;

				if(j)
				{
					old_pos2 = old_pos-1;
					old_flag2 = old_flag;
					old_quality2 = old_quality;
				}
				else
				{
					old_pos1 = old_pos-1;
					old_flag1 = old_flag;
					old_quality1 = old_quality;
				}
			}
			if(old_chro1[0]=='*')old_chro1[0]=0;
			if(old_chro2[0]=='*')old_chro2[0]=0;

			old_inb1_reversed = (old_flag1 & 0x10)?1:0;
			old_inb2_reversed = (old_flag2 & 0x10)?1:0;

			#ifdef QUALITY_KILL_SUBREAD
				if(old_quality1 <= QUALITY_KILL_SUBREAD || (old_flag1 & SAM_FLAG_UNMAPPED))
				{
					/*
					if(old_inb1_reversed)
					{
						int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb1[0]>='A' && inb1[0]<='Z')?1:0;
						reverse_read(inb1+ rev_offset, rl1- rev_offset, ginp->space_type);
						if(qualityb1[0])
							reverse_quality(qualityb1, rl1- rev_offset);
					}
					*/

					old_chro1[0]=0;
					old_flag1=SAM_FLAG_UNMAPPED;
					old_quality1=0;
					old_cigar1[0]=0;
					old_pos1=0;
				}


				if(old_quality2 <= QUALITY_KILL_SUBREAD || (old_flag2 & SAM_FLAG_UNMAPPED))
				{
					/*

					if(old_inb2_reversed)
					{
						int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb2[0]>='A' && inb2[0]<='Z')?1:0;
						reverse_read(inb2+ rev_offset, rl2- rev_offset, ginp->space_type);
						if(qualityb2[0])
							reverse_quality(qualityb2, rl2- rev_offset);
					}

					*/


					old_chro2[0]=0;
					old_flag2=SAM_FLAG_UNMAPPED;
					old_quality2=0;
					old_cigar2[0]=0;
					old_pos2=0;
				}
			#endif
		}
		else
		{
			rl1 = geinput_next_read(ginp, nameb1, inb1, qualityb1);
			rl2 = geinput_next_read(ginp, nameb2, inb2, qualityb2);
		}


		if (rl1<0){
			break;
		}
		if (rl2<0){
			break;
		}



		if(!qualityb1[0])
		{
			int xx;
			EXON_FASTQ_FORMAT = FASTQ_PHRED33;
			for (xx=0; xx<rl1; xx++) qualityb1[xx]='J';
			qualityb1[xx]='\0';

			for (xx=0; xx<rl2; xx++) qualityb2[xx]='J';
			qualityb2[xx]='\0';
		}

		if(EXON_FASTQ_FORMAT == FASTQ_PHRED64)
		{
			fastq_64_to_33(qualityb1);
			fastq_64_to_33(qualityb2);
		}

		if (ginp->file_type == GENE_INPUT_PLAIN)
		{
			sprintf (nameb1, "Raw-%d", i/2);
			sprintf (nameb2, "Raw-%d", i/2);
		}
		else
		{
			remove_backslash(nameb1);
			remove_backslash(nameb2);
		}


		for(j = 0; j<2; j++)
		{
			int rl = j?rl2:rl1;
			int use_old_line = 1;
			if((halves_record -> best_vote1_list[i]>=1 && halves_record -> best_vote2_list[i]>=1 ))
			{
				char * chro_name;
				unsigned int chro_pos;
				char *cigar_print;

				unsigned int pos = halves_record -> best_pos1_list[i];

				locate_gene_position(pos, &_global_offsets, &chro_name, &chro_pos);

				pos = halves_record -> is_reversed_list[i]? (pos + rl):(pos);

				gene_quality_score_t mapping_quality;

				if (halves_record -> cigar_string_buffer[i * EXON_MAX_CIGAR_LEN])
				{
					cigar_print = halves_record -> cigar_string_buffer + i * EXON_MAX_CIGAR_LEN;
					mapping_quality = halves_record -> final_quality[i] ;



					if(j)
					{
						strcpy(old_chro2, chro_name);
						strcpy(old_cigar2, cigar_print); 
						old_pos2 = chro_pos;
						old_quality2 = mapping_quality;
						old_flag2 = halves_record -> is_reversed_list[i]?0:16;


						if((halves_record -> is_reversed_list[i]?1:0) == old_inb2_reversed)
						{
							int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb2[0]>='A' && inb2[0]<='Z')?1:0;
							reverse_read(inb2+ rev_offset, rl2- rev_offset, ginp->space_type);
							if(qualityb2[0])
								reverse_quality(qualityb2, rl2- rev_offset);
						}


					}else
					{
						strcpy(old_chro1, chro_name);
						strcpy(old_cigar1, cigar_print); 
						old_pos1 = chro_pos;
						old_quality1 = mapping_quality;
						old_flag1 = halves_record -> is_reversed_list[i]?16:0;


						if((halves_record -> is_reversed_list[i]?1:0) != old_inb1_reversed)
						{
							int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb1[0]>='A' && inb1[0]<='Z')?1:0;
							reverse_read(inb1+ rev_offset, rl1- rev_offset, ginp->space_type);
							if(qualityb1[0])
								reverse_quality(qualityb1, rl1- rev_offset);
						}
					}

					use_old_line = 0;
				}
			}
			if(use_old_line)
			{
				if(old_line1[0])
				{
					if(j)
					{
						if(old_inb2_reversed && !old_chro2[0])
						{
							int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb2[0]>='A' && inb2[0]<='Z')?1:0;
							reverse_read(inb2+ rev_offset, rl2- rev_offset, ginp->space_type);
							if(qualityb2[0])
								reverse_quality(qualityb2, rl2- rev_offset);
						}
					}
					else
					{
						if(old_inb1_reversed && !old_chro1[0])
						{
							int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb1[0]>='A' && inb1[0]<='Z')?1:0;
							reverse_read(inb1+ rev_offset, rl1- rev_offset, ginp->space_type);
							if(qualityb1[0])
								reverse_quality(qualityb1, rl1- rev_offset);
						}
					}
				}
				else
				{

					char * chro_name;
					unsigned int chro_pos;
					unsigned int pos = halves_record -> best_pos1_list[i];

					locate_gene_position(pos, &_global_offsets, &chro_name, &chro_pos);
					if (halves_record -> best_vote1_list[i]>=7)
					{
						if(j)
						{
							strcpy(old_chro2, chro_name);
							sprintf(old_cigar2, "%dM", rl2);
							old_pos2 = pos;
							old_flag1 = halves_record -> is_reversed_list[i]?16:0 ;


							if(halves_record -> is_reversed_list[i])
							{
								int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb2[0]>='A' && inb2[0]<='Z')?1:0;
								reverse_read(inb2+ rev_offset, rl2- rev_offset, ginp->space_type);
								if(qualityb2[0])
									reverse_quality(qualityb2, rl2- rev_offset);
							}

						}
						else
						{
							strcpy(old_chro2, chro_name);
							sprintf(old_cigar2, "%dM", rl2);
							old_pos2 = pos;
							old_flag1 = halves_record -> is_reversed_list[i]?16:0 ;

							if(!halves_record -> is_reversed_list[i])
							{
								int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb1[0]>='A' && inb1[0]<='Z')?1:0;
								reverse_read(inb1+ rev_offset, rl1- rev_offset, ginp->space_type);
								if(qualityb1[0])
									reverse_quality(qualityb1, rl1- rev_offset);
							}

						}
					}
				}
			}

			++i;
		}



		{
			char mate_for_1[120], mate_for_2[120];
			int tlen = 0;
			if(old_chro1[0] && (strcmp(old_chro1, old_chro2) == 0))
			{
				long long int pair_dist = old_pos1;
				pair_dist -= old_pos2;
				tlen = abs(pair_dist) + (tlen>0?rl1:rl2);
				if(tlen >= EXON_MIN_PAIRED_DISTANCE && tlen <= EXON_MAX_PAIRED_DISTANCE)
				{
					old_flag1 |= SAM_FLAG_MATCHED_IN_PAIR;
					old_flag2 |= SAM_FLAG_MATCHED_IN_PAIR;
				}
				else
				{
					old_flag1 &= ~SAM_FLAG_MATCHED_IN_PAIR;
					old_flag2 &= ~SAM_FLAG_MATCHED_IN_PAIR;
				}
				strcpy(mate_for_1,"=");
				strcpy(mate_for_2,"=");
			}
			else
			{
				if(old_chro2[0])
					strcpy(mate_for_1,old_chro2);
				else{
					strcpy(mate_for_1,"*");
					old_pos2=0;
				}
				
				if(old_chro1[0])
					strcpy(mate_for_2,old_chro1);
				else
				{
					strcpy(mate_for_2,"*");
					old_pos1=0;
				}

				old_flag1 &= ~SAM_FLAG_MATCHED_IN_PAIR;
				old_flag2 &= ~SAM_FLAG_MATCHED_IN_PAIR;
			}

			if(old_chro2[0])
				old_flag1 &= ~SAM_FLAG_MATE_UNMATCHED;
			else
				old_flag1 |= SAM_FLAG_MATE_UNMATCHED;

			if(old_chro1[0])
				old_flag2 &= ~SAM_FLAG_MATE_UNMATCHED;
			else
				old_flag2 |= SAM_FLAG_MATE_UNMATCHED;

			if(old_flag1 & SAM_FLAG_REVERSE_STRAND_MATCHED)
				old_flag2 |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
			else
				old_flag2 &= ~SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;

			if(old_flag2 & SAM_FLAG_REVERSE_STRAND_MATCHED)
				old_flag1 |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
			else
				old_flag1 &= ~SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;

			old_flag1 |= SAM_FLAG_FIRST_READ_IN_PAIR | SAM_FLAG_PAIRED_TASK;
			old_flag2 |= SAM_FLAG_SECOND_READ_IN_PAIR| SAM_FLAG_PAIRED_TASK;

			
			if(old_chro1[0]) old_pos1++;
			else
			{
				strcpy(old_chro1,"*");
				tlen = 0;
			}

			if(old_chro2[0]) old_pos2++;
			else
			{
				strcpy(old_chro2,"*");
				tlen = 0;
			}

			if(!old_cigar1[0])
				strcpy(old_cigar1,"*");
			if(!old_cigar2[0])
				strcpy(old_cigar2,"*");

			if((!EXON_JUNCTION_READS_ONLY) || strchr(old_cigar1, 'N'))
				fprintf(out_fp,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%d\t%s\t%s\n",  nameb1, old_flag1, old_chro1, old_pos1, old_quality1 , old_cigar1,mate_for_1, old_pos2, tlen,  inb1, qualityb1);
			if((!EXON_JUNCTION_READS_ONLY) || strchr(old_cigar2, 'N'))
				fprintf(out_fp,"%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%d\t%s\t%s\n",  nameb2, old_flag2, old_chro2, old_pos2, old_quality2 , old_cigar2,mate_for_2, old_pos1, -tlen,  inb2, qualityb2);

		}

	}

	SUBREADprintf("\n");

}




struct gene_thread_data
{
	gehash_t table;
	int * offsets;
	int number_t;
};


struct gene_thread_data_transport
{
	gehash_t * my_table;
	gene_value_index_t * my_value_array_index;
	int table_no;
	int all_tables;
	char * index_prefix;
	halves_record_t * halves_record;
	gene_input_t * ginp;
	gene_input_t * ginp2;
	unsigned long long int base_number;
	unsigned int * processed_reads;
	unsigned int section_length;
	int all_threads;
	int this_thread;
	HashTable * bed_table;
	HashTable * connection_table;
	HashTable * pos_table;

	int tolerable_scan;

	pthread_spinlock_t * input_data_lock;
	pthread_spinlock_t * init_lock;
};

void fragile_junction_voting(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table, gehash_t * my_table, gene_value_index_t * my_value_array_index , int table_no,  halves_record_t * halves_record, char * read, char * qual, unsigned int full_rl, int negative_strand, int color_space, unsigned int low_border, unsigned int high_border, gene_vote_t *vote_p1)
{
	int windows = full_rl / EXON_LARGE_WINDOW +1;
	float overlap = (1.0*windows * EXON_LARGE_WINDOW - full_rl) / (windows-1);

	int ww;
	int window_cursor = 0;
	for(ww=0; ww<windows;ww++)
	{
		window_cursor = (int)(ww * EXON_LARGE_WINDOW - ww * overlap);
		int read_len = EXON_LARGE_WINDOW;
		if(ww == windows-1)
			read_len = full_rl -window_cursor;

		float subread_step = 3.00001;

		int i;
		int subread_no;
		char * InBuff;
		InBuff = read + window_cursor;
		char tmp_char = InBuff[read_len];
		InBuff[read_len] = 0;
		
		init_gene_vote(vote_p1);
		for(subread_no=0; ; subread_no++)
		{
			int subread_offset1 = (int)(subread_step * (subread_no+1));
			subread_offset1 -= subread_offset1%GENE_SLIDING_STEP;
			subread_offset1 += GENE_SLIDING_STEP-1;

			for(i=0; i<GENE_SLIDING_STEP ; i++)
			{
				int subread_offset = (int)(subread_step * subread_no); 
				subread_offset -= subread_offset%GENE_SLIDING_STEP -i;

				char * subread_string = InBuff + subread_offset;

				gehash_key_t subread_integer = genekey2int(subread_string, color_space);

				//SUBREADprintf("TQ: POS=%d, TOL=%d, INT=%u, SR=%s\n", subread_offset  , tolerable_scan , subread_integer, subread_string);

				if(EXON_MAX_METHYLATION_C_NUMBER)
					gehash_go_q_CtoT(my_table, subread_integer , subread_offset, read_len,negative_strand, vote_p1, 1, 1, 21.9, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no, EXON_MAX_METHYLATION_C_NUMBER,  low_border, high_border - read_len);
				else
					gehash_go_q(my_table, subread_integer , subread_offset, read_len,negative_strand, vote_p1, 1, 1, 21.9, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no,  low_border, high_border - read_len);
			}
			if(subread_offset1 >= read_len -16)
				break;
		}


		if(1)
		{
			finalise_vote(vote_p1);
			//print_votes(&vote_p1, index_prefix);
			unsigned int best_pos1=0;
			unsigned int best_pos2=0;
			int best_vote1=0;
			int best_vote2=0;
			char is_abnormal=0;
			short half_marks=0;
			int is_reversed_halves=0, max_cover_start=0, max_cover_end=0;
			char indel_in_p1=0, indel_in_p2=0;
			short read_coverage_start =0, read_coverage_end=0;

			int splice_point = select_best_matching_halves(vote_p1, &best_pos1, &best_pos2, &best_vote1, &best_vote2, &is_abnormal ,&half_marks, &is_reversed_halves, accepted_support_rate, read_len, -1,  0, &read_coverage_start, &read_coverage_end, &indel_in_p1, &indel_in_p2, &max_cover_start, &max_cover_end, read_len, -1 , 0, NULL , 0xffffffff);
			if (splice_point>0 && best_vote1 >= 1 && best_vote2>=1)
			{
				int test_real_break_point = -1, test_donor_score=-1;
				int is_GTAG = 0;
				int is_accepted = test_donor(InBuff, read_len, min(best_pos1, best_pos2), max(best_pos1,best_pos2), splice_point, negative_strand, read_len/4, 0, EXON_INDEL_TOLERANCE, &test_real_break_point, my_value_array_index, 0, 0, negative_strand, color_space, halves_record ->best_vote2_list[i], &test_donor_score, &is_GTAG);

				if (is_accepted ){
					unsigned int pos_R1 = test_real_break_point+ best_pos1;
					unsigned int pos_R2 = test_real_break_point+ best_pos2;

					pos_R1 = get_grouped_position(pos_table, pos_R1);
					pos_R2 = get_grouped_position(pos_table, pos_R2);

					unsigned int pos_small = min(pos_R1, pos_R2);
					unsigned int pos_big   = max(pos_R1, pos_R2);

					paired_exon_key search_key; 
					search_key.small_key = pos_small;
					search_key.big_key = pos_big;
					put_connection_table(connection_table, pos_R1, pos_R2, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1)?1:0, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R2)?1:0);
					put_connection_table(connection_table, pos_R2, pos_R1, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R2)?1:0, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1)?1:0);

					exon_junction_t *search_res = (exon_junction_t*) HashTableGet(bed_table, &search_key);
					if(search_res)
						search_res -> feed_supporting_reads ++;
					else
					{
						//SUBREADprintf("FRAG: the %d-th window: %d : %d \t\tFOUND: %u-%u\n", ww, window_cursor , read_len,pos_R1, pos_R2 );

						search_res = (exon_junction_t *)malloc(sizeof(exon_junction_t));
						search_res->supporting_reads = 0;
						search_res -> feed_supporting_reads =0;
						search_res->strand = !is_GTAG;
						search_res -> left_extend = 0;
						search_res -> right_extend = 0;

						paired_exon_key * new_key = (paired_exon_key *) malloc(sizeof(paired_exon_key));
						new_key->small_key = pos_small;
						new_key->big_key = pos_big;
						HashTablePut(bed_table, new_key, search_res);
					}
				}
			}
		}
		InBuff[read_len] = tmp_char;
	}
}

int run_exon_search(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table , gehash_t * my_table, gene_value_index_t * my_value_array_index , int table_no,  halves_record_t * halves_record, gene_input_t * ginp,gene_input_t * ginp2, char * index_prefix, unsigned int * processed_reads, long long int base_number, int all_tables, pthread_spinlock_t * input_lock, int my_thread_no, unsigned int section_length, int tolerable_scan)
{
	//	FILE * fp;
	char BuffMemory [2500];
	char * InBuff = NULL , * InBuff2 = NULL;
	char BuffMemory2 [2500];
	char * QualityBuff = NULL, * QualityBuff2 = NULL;

	int all_reads = halves_record -> max_len;
	int queries = 0, i;
	double t0=miltime();
	int is_reversed;
	int good_match = 0;
	float subread_step; // = EXON_SUBREAD_GAP+.0001;
	struct stat read_fstat;
	stat (ginp->filename, &read_fstat);
	long long int read_fsize = read_fstat.st_size;

	double local_begin_ftime = miltime();
	unsigned int sam_pos1=0, sam_pos2=0;
	int read_len = 0, read2_len = 0, sam_qual1=0, sam_qual2=0, sam_flag1=0, sam_flag2=0;
	is_reversed = 0;

	if(my_thread_no==0)
		SUBREADprintf("Processing %s:\n", ginp2?"fragments":"reads");
	//SUBREADprintf ("I'm the %d-th thread RRRR\n", my_thread_no);return 0;

	if (ginp2)
		all_reads /=2;


	gene_vote_t * vote_p1, * vote_p2, *vote_p3;

	unsigned int low_border = my_value_array_index -> start_base_offset;
	unsigned int high_border = my_value_array_index -> start_base_offset + my_value_array_index -> length; 
	unsigned int index_valid_range =  my_value_array_index -> start_base_offset + my_value_array_index -> length;

	vote_p1 = (gene_vote_t *)malloc(sizeof(gene_vote_t));
	vote_p2	= (gene_vote_t *)malloc(sizeof(gene_vote_t));
	vote_p3	= (gene_vote_t *)malloc(sizeof(gene_vote_t));

	if(!vote_p1 || !vote_p2 || !vote_p3)
	{
		SUBREADprintf("Unable to allocate memory for voting\n");
		return -1;
	}

	// if this block of index is not the last one

	if(_global_offsets.read_offsets[ _global_offsets.total_offsets ] > high_border+10)
		for(i=0; i< _global_offsets.total_offsets; i++)
		{
			if(_global_offsets.read_offsets[i] > index_valid_range)
			{
				if(i>0 && _global_offsets.read_offsets[i-1] > index_valid_range-2000000)
					index_valid_range = _global_offsets.read_offsets[i-1] ;
				else
					index_valid_range = index_valid_range-2000000 ;
				break;
			}
		}

	while (1)
	{

		char namebuf[200];


		if (my_thread_no==0 && queries % (all_reads / 14) == 0 && !is_reversed)
		{
			if(table_no == 0)
			{
				long long int current_reads = base_number + queries;
				long long int fpos = ftello(ginp->input_fp);
				reads_density = fpos*1.0/current_reads; 
				//SUBREADprintf ("\nDENS=%.5f, POS=%llu, fsize=%llu\n", reads_density, fpos, read_fsize);
			}
			if(IS_DEBUG && queries % (all_reads/14)==0)
				SUBREADprintf("@LOG Done %d/%d, good %d, last time %f\n",queries, table_no, good_match, miltime() - t0);
			else
			{
				long long int all_steps = (read_fsize*1.0/reads_density) * all_tables;
				int remaining_load_libs = (all_tables - table_no);
				remaining_load_libs +=  all_tables * (int)(((read_fsize*1.0/reads_density) - base_number)/all_reads);
				long long int finished_steps =  ((section_length)*table_no + base_number*all_tables+queries);
				double finished_rate = finished_steps*1.0 / all_steps;
				double reads_per_second =  queries*1.0 / all_tables / (miltime()- local_begin_ftime);
				double expected_seconds = ((1.-finished_rate) * all_steps)/all_tables / reads_per_second + remaining_load_libs * 50 + (int)(((read_fsize*1.0/reads_density) - base_number)/all_reads)*100;
				if(queries>1)
					print_running_log(finished_rate, reads_per_second, expected_seconds, (unsigned long long int)all_steps / all_tables, ginp2 != NULL);
			}
			SUBREADfflush(stdout);
			t0 = miltime();
		}

		if (is_reversed)
		{

			reverse_read(InBuff, read_len, ginp->space_type);	
			reverse_quality(QualityBuff, read_len);
			if (ginp2)
			{
				reverse_read(InBuff2, read2_len, ginp->space_type);
				reverse_quality(QualityBuff2, read2_len);
			}

			if(!EXON_FUSION_DETECTION)
			{
				init_gene_vote(vote_p1);
				if(ginp2)
					init_gene_vote(vote_p2);
			}
		}
		else
		{
			if(input_lock!=NULL)
				pthread_spin_lock(input_lock);

			if(*processed_reads < all_reads)
			{
				InBuff = BuffMemory;
				QualityBuff = BuffMemory2;

				read_len = geinput_next_read_sam(ginp, namebuf, InBuff, QualityBuff, &_global_offsets, &sam_pos1, &sam_qual1, &sam_flag1, EXON_FIRST_READ_REVERSE);

				if(read_len>0){
					queries = *processed_reads;
					(*processed_reads ) ++;
				}
				if(ginp2)
				{
					InBuff2 = BuffMemory + 1250;
					QualityBuff2 = BuffMemory2 + 1250;

					read2_len = geinput_next_read_sam(ginp, namebuf, InBuff2, QualityBuff2, &_global_offsets, &sam_pos2, &sam_qual2, &sam_flag2, EXON_SECOND_READ_REVERSE);
				}

				if(input_lock!=NULL)
					pthread_spin_unlock(input_lock);
			}
			else
			{
				if(input_lock!=NULL)
					pthread_spin_unlock(input_lock);
				break;
			}

			if (read_len<0)
				break;

			if(sam_qual1 > MAX_QUALITY_TO_CALL_JUNCTION && ((!ginp2)|| sam_qual2 > MAX_QUALITY_TO_CALL_JUNCTION))
			{
				int qindex = queries;
				if(ginp2)qindex = qindex*2;
				halves_record -> best_vote1_list[qindex] = sam_qual1 ;
				halves_record -> best_vote2_list[qindex] = 0 ;
				halves_record -> best_pos1_list[qindex] = sam_pos1 ;
				halves_record -> is_reversed_list[qindex] = (sam_flag1 & SAM_FLAG_REVERSE_STRAND_MATCHED) ?1:0;
				halves_record -> half_marks_list[qindex] = 0;
				halves_record -> best1_read_end_pos[qindex] = read_len * 0.9;
				halves_record -> best1_read_start_pos[qindex] = read_len * 0.1;

				halves_record -> indel_in_piece1[qindex] = 0;
				if(ginp2)
				{
					qindex++;
					halves_record -> best_vote1_list[qindex] = sam_qual2 ;
					halves_record -> best_vote2_list[qindex] = 0 ;
					halves_record -> best_pos1_list[qindex] = sam_pos2 ;
					halves_record -> is_reversed_list[qindex] = (!(sam_flag2 & SAM_FLAG_REVERSE_STRAND_MATCHED)) ?1:0;
					halves_record -> half_marks_list[qindex] = 0;
					halves_record -> best1_read_end_pos[qindex] = read2_len * 0.9;
					halves_record -> best1_read_start_pos[qindex] = read2_len * 0.1;
					halves_record -> indel_in_piece1[qindex] = 0;
				}
				continue;
			}


			init_gene_vote(vote_p1);
			if(ginp2)
				init_gene_vote(vote_p2);

			/*

			if (ginp->space_type == GENE_SPACE_COLOR && InBuff[0]>='A' && InBuff[0]<='Z')
			{
				InBuff ++;
				if(QualityBuff[0])
					QualityBuff++;
				read_len --;
			}

			if(EXON_FIRST_READ_REVERSE)
			{
				reverse_read(InBuff, read_len, ginp->space_type);
				reverse_quality(QualityBuff,read_len);
			}

			if(ginp2)
			{
				if (ginp->space_type == GENE_SPACE_COLOR && InBuff2[0]>='A' && InBuff2[0]<='Z')
				{
					InBuff2 ++;
					if(QualityBuff2[0])
						QualityBuff2++;
					read2_len --;
				}

				if(EXON_SECOND_READ_REVERSE)
				{
					reverse_quality(QualityBuff2, read2_len);
					reverse_read(InBuff2, read2_len, ginp->space_type);
				}
			}*/

		}


		int repeated_record_pos_base = table_no>1?-1:( 12 * (table_no*2 + is_reversed));


		//if((sam_flag1 & SAM_FLAG_UNMAPPED) && ((!ginp2)|| (sam_flag2 & SAM_FLAG_UNMAPPED)))
		//	continue;

		if (tolerable_scan )
		{
			if (ginp2 && ( (halves_record->half_marks_list[queries*2] & IS_FINALISED_PROCESSING) && (halves_record->half_marks_list[queries*2+1] & IS_FINALISED_PROCESSING)))
			{
				is_reversed = 0;
				continue;
			}
			if (!ginp2 && ( (halves_record->half_marks_list[queries] & IS_FINALISED_PROCESSING)))
			{
				is_reversed = 0;
				continue;
			}
		}


		if (ginp2)
		{

			int i;
			int subread_no;
			int is_second_read;
			long long int hint_pos;
			for (is_second_read = 0; is_second_read <2; is_second_read ++)
			{
				gene_vote_t * current_vote = is_second_read?vote_p2: vote_p1;
				char * current_read =  is_second_read?InBuff2 : InBuff;
				int current_rlen = is_second_read?read2_len:read_len;
				if(tolerable_scan  && (!EXON_NO_TOLERABLE_SCAN)) subread_step = 3.00001;
				else
					subread_step =  (int)((read_len - 18)*1./(TOTAL_SUBREADS-1))+0.000001;
				for(subread_no=0; ; subread_no++)
				{
					int subread_offset1 = (int)(subread_step * (subread_no+1));
					subread_offset1 -= subread_offset1%GENE_SLIDING_STEP;
					subread_offset1 += GENE_SLIDING_STEP-1;

					for(i=0; i<GENE_SLIDING_STEP ; i++)
					{
						int subread_offset = (int)(subread_step * subread_no); 
						subread_offset -= subread_offset%GENE_SLIDING_STEP -i;

						char * subread_string = current_read + subread_offset;

						gehash_key_t subread_integer = genekey2int(subread_string, ginp->space_type);

						if(EXON_MAX_METHYLATION_C_NUMBER)
							gehash_go_q_CtoT(my_table, subread_integer , subread_offset, current_rlen, is_reversed, current_vote, 1, 1, 22, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no, EXON_MAX_METHYLATION_C_NUMBER,  low_border, high_border - current_rlen);
						else
							gehash_go_q(my_table, subread_integer , subread_offset, current_rlen, is_reversed, current_vote, 1, 1, 22, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no,  low_border, high_border - current_rlen);
					}
					if(subread_offset1 >= current_rlen -16)
						break;
				}
			}

			gene_vote_number_t numvote_read1, numvote_read2;
			gehash_data_t pos_read1, pos_read2;
			gene_quality_score_t sum_quality, qual_r1, qual_r2;
			char record_index1[48], record_index2[48];
			int is_breakeven = 0;


			if((!EXON_FUSION_DETECTION) || is_reversed)
			{
				finalise_vote(vote_p1);
				finalise_vote(vote_p2);
				int is_paired_match = select_positions_exons(vote_p1, vote_p2, &numvote_read1, &numvote_read2, &sum_quality, &qual_r1, &qual_r2, &pos_read1, &pos_read2, record_index1, record_index2, EXON_MAX_PAIRED_DISTANCE, EXON_MIN_PAIRED_DISTANCE, 3, ACCEPT_MINOR_SUBREADS, is_reversed, EXON_NUMBER_OF_ANCHORS_PAIRED,  EXON_INDEL_TOLERANCE, &is_breakeven, read_len, read2_len);
				for (is_second_read = 0; is_second_read <2; is_second_read ++)
				{
					gene_vote_t * current_vote = is_second_read?vote_p2: vote_p1;
					int current_rlen = is_second_read?read2_len:read_len;

					unsigned int best_pos1=0;
					unsigned int best_pos2=0;
					hint_pos = is_second_read?pos_read2:pos_read1;
					if (!is_paired_match ) hint_pos = -1;
					int best_vote1=0;
					int best_vote2=0;
					char is_abnormal=0;
					short half_marks=0;
					char indel_in_p1=0, indel_in_p2=0;
					int is_reversed_halves=0, max_cover_start=0, max_cover_end=0;
					short read_coverage_start = 0, read_coverage_end = 0;

					int splice_point = select_best_matching_halves(current_vote, &best_pos1, &best_pos2, &best_vote1, &best_vote2, &is_abnormal ,&half_marks, &is_reversed_halves, accepted_support_rate, current_rlen, hint_pos, tolerable_scan, &read_coverage_start, &read_coverage_end, &indel_in_p1, &indel_in_p2, &max_cover_start, &max_cover_end, current_rlen,  repeated_record_pos_base, is_reversed, halves_record -> cigar_string_buffer+(queries*2+is_second_read) * EXON_MAX_CIGAR_LEN, index_valid_range);

					if (splice_point>0)
					{
					//	if (is_second_read)printf ("SECOND!\n");
						add_best_matching_halves(halves_record, best_pos1, best_pos2, best_vote1, best_vote2,is_abnormal, is_reversed_halves, splice_point, half_marks, 2*queries+is_second_read, max_cover_start, max_cover_end, read_coverage_start, read_coverage_end, indel_in_p1, indel_in_p2);
					}
					else 
					{
						is_reversed_halves = (current_vote -> max_mask & IS_NEGATIVE_STRAND)?1:0;
						if (current_vote->max_vote > (halves_record -> best_vote1_list[queries*2+is_second_read] + halves_record -> best_vote2_list[queries*2+is_second_read]))
						{
							halves_record -> best_vote1_list[queries*2+is_second_read] = current_vote->max_vote ;
							halves_record -> best_vote2_list[queries*2+is_second_read] = 0 ;
							halves_record -> best_pos1_list[queries*2+is_second_read] = current_vote->max_position;
							halves_record -> is_reversed_list[queries*2+is_second_read] = is_reversed_halves;
							halves_record -> half_marks_list[queries*2+is_second_read] = (halves_record -> half_marks_list[queries*2+is_second_read]) & ~(IS_PAIRED_HINTED);

							halves_record -> best1_read_start_pos[queries*2+is_second_read] = current_vote -> max_coverage_start;
							halves_record -> best1_read_end_pos[queries*2+is_second_read] = current_vote -> max_coverage_end;
							halves_record -> read_coverage_start[queries*2+is_second_read] = current_vote -> max_coverage_start;
							halves_record -> read_coverage_end[queries*2+is_second_read] = current_vote -> max_coverage_end;
							if (is_paired_match)
								halves_record -> half_marks_list[queries*2+is_second_read] = (halves_record -> half_marks_list[queries*2+is_second_read]) | IS_PAIRED_HINTED;
	
						}
					}

					if( current_rlen >= EXON_LONG_READ_LENGTH)
						fragile_junction_voting(bed_table, pos_table, connection_table, my_table, my_value_array_index , table_no, halves_record, is_second_read?InBuff2:InBuff,  is_second_read?QualityBuff2:QualityBuff, current_rlen , is_reversed, ginp->space_type,  low_border,  high_border, vote_p3);
				}
			}
		}
		else
		{

			int i;
			int subread_no;
		
			if(tolerable_scan  && (!EXON_NO_TOLERABLE_SCAN) ) subread_step = 3.00001;
			else
				subread_step =  (int)((read_len - 18)*1./(TOTAL_SUBREADS-1)) + 0.00001;
			for(subread_no=0; ; subread_no++)
			{
				int subread_offset1 = (int)(subread_step * (subread_no+1));
				subread_offset1 -= subread_offset1%GENE_SLIDING_STEP;
				subread_offset1 += GENE_SLIDING_STEP-1;

				for(i=0; i<GENE_SLIDING_STEP ; i++)
				{
					int subread_offset = (int)(subread_step * subread_no); 
					subread_offset -= subread_offset%GENE_SLIDING_STEP -i;

					char * subread_string = InBuff + subread_offset;

					gehash_key_t subread_integer = genekey2int(subread_string, ginp->space_type);

					//SUBREADprintf("TQ: POS=%d, TOL=%d, INT=%u, SR=%s\n", subread_offset  , tolerable_scan , subread_integer, subread_string);
					if(EXON_MAX_METHYLATION_C_NUMBER)
						gehash_go_q_CtoT(my_table, subread_integer , subread_offset, read_len, is_reversed, vote_p1, 1, 1, 22, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no, EXON_MAX_METHYLATION_C_NUMBER,  low_border, high_border - read_len);
					else
						gehash_go_q(my_table, subread_integer , subread_offset, read_len, is_reversed, vote_p1, 1, 1, 22, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no,  low_border, high_border - read_len);
				}
				if(subread_offset1 >= read_len -16)
					break;
			}

			if((!EXON_FUSION_DETECTION) || is_reversed)
			{
				finalise_vote(vote_p1);
				//print_votes(&vote_p1, index_prefix);
				unsigned int best_pos1=0;
				unsigned int best_pos2=0;
				int best_vote1=0;
				int best_vote2=0;
				char is_abnormal=0;
				short half_marks=0;
				int is_reversed_halves=0, max_cover_start=0, max_cover_end=0;
				char indel_in_p1=0, indel_in_p2=0;
				short read_coverage_start =0, read_coverage_end=0;

				int splice_point = select_best_matching_halves(vote_p1, &best_pos1, &best_pos2, &best_vote1, &best_vote2, &is_abnormal ,&half_marks, &is_reversed_halves, accepted_support_rate, read_len, -1,  tolerable_scan, &read_coverage_start, &read_coverage_end, &indel_in_p1, &indel_in_p2, &max_cover_start, &max_cover_end, read_len, repeated_record_pos_base, is_reversed, halves_record -> cigar_string_buffer+queries * EXON_MAX_CIGAR_LEN, index_valid_range);
	
				if (splice_point>0 && (vote_p1->max_vote >= halves_record -> best_vote1_list[queries]))
				{
					#ifdef DEBUG
					SUBREADprintf("SPP=%d, v1=%d, v2=%d, POS=%u, %%=%s\nRead_Start=%d, Read_End=%d\n", splice_point, best_vote1, best_vote2, best_pos1, InBuff, read_coverage_start, read_coverage_end);
					#endif

					add_best_matching_halves(halves_record, best_pos1, best_pos2, best_vote1, best_vote2,is_abnormal, is_reversed_halves, splice_point, half_marks, queries, max_cover_start, max_cover_end, read_coverage_start, read_coverage_end, indel_in_p1, indel_in_p2);
				}
				else 
				{
					//SUBREADprintf("EDD=%d, v1=%d, v2=%d, R=%s, RL=%d, MM=%d, SPS=%.5f\n", splice_point, best_vote1, best_vote2, InBuff, read_len, vote_p1->max_vote, subread_step);
					is_reversed_halves = (vote_p1->max_mask & IS_NEGATIVE_STRAND)?1:0;
					if (vote_p1->max_vote > (halves_record -> best_vote1_list[queries] + 0* halves_record -> best_vote2_list[queries]))
					{
						halves_record -> best_vote1_list[queries] = vote_p1->max_vote ;
						halves_record -> best_vote2_list[queries] = 0 ;
						halves_record -> best_pos1_list[queries] = vote_p1->max_position;
						halves_record -> is_reversed_list[queries] = is_reversed_halves;
						halves_record -> half_marks_list[queries] = vote_p1->max_mask;
						halves_record -> best1_read_start_pos[queries] = vote_p1->max_coverage_start;
						halves_record -> best1_read_end_pos[queries] = vote_p1->max_coverage_end;

						halves_record -> read_coverage_start[queries] = vote_p1->max_coverage_start;
						halves_record -> read_coverage_end[queries] = vote_p1->max_coverage_end;
					}
				}
			}
			if(read_len >= EXON_LONG_READ_LENGTH)
				fragile_junction_voting(bed_table, pos_table, connection_table, my_table, my_value_array_index , table_no, halves_record, InBuff,  QualityBuff, read_len , is_reversed, ginp->space_type,  low_border,  high_border, vote_p3);

		}

	
		if (is_reversed)
		{
			if(queries >= all_reads)
				break;
		}

		is_reversed = !is_reversed;
	}

	free(vote_p1);
	free(vote_p2);
	free(vote_p3);
	if(my_thread_no==0)
		SUBREADprintf("\n");
	return 0;
}


void * run_exon_search_thread(void * parameters)
{
	struct gene_thread_data_transport * data_param = parameters;
	int thid = data_param->this_thread;

	pthread_spin_unlock(data_param -> init_lock);

	run_exon_search(data_param-> bed_table, data_param->pos_table, data_param->connection_table, data_param->my_table, data_param->my_value_array_index, data_param->table_no, data_param->halves_record, data_param->ginp , data_param->ginp2, data_param->index_prefix, data_param->processed_reads, data_param->base_number, data_param->all_tables, data_param->input_data_lock, thid, data_param-> section_length, data_param -> tolerable_scan);
	return NULL;
}







// This function search a segment of reads (length = read_number) 
// It returns the number of reads that were really processed;
int run_exon_search_index_tolerable(gene_input_t * ginp, gene_input_t * ginp2, char * index_prefix, halves_record_t * halves_record, FILE * out_fp, unsigned long long int base_number, int all_tables, unsigned long long int *succeed_reads, HashTable * overall_exon_bed, HashTable * pos_table, HashTable * connection_table, int tolerable_scan)
{
	unsigned int tabno=0;
	unsigned int processed_reads = 0, section_length = 0;
	unsigned long long int current_fp = ftello(ginp -> input_fp);
	unsigned long long int last_fp_pos =  0;
	struct stat filestat;
	int last_table = 0;
	int stat_ret;
	char table_fn [300];
	gehash_t my_raw_table;
	gehash_t * my_table = &my_raw_table;
	gene_value_index_t value_array_index;

	while (1)
	{
		sprintf(table_fn, "%s.%02d.%c.tab", index_prefix, tabno, ginp->space_type==GENE_SPACE_COLOR?'c':'b');

// Test and load the index partition
		stat_ret = stat(table_fn, &filestat);
		if (stat_ret !=0)
			break;

		if (IS_DEBUG)
			SUBREADprintf ("@LOG Loading table from %s\n", table_fn);
		else
			SUBREADprintf ("Loading the %02d-th index file ...					      \n", tabno+1);
		SUBREADfflush(stdout);

		if(gehash_load(my_table, table_fn)) return -1;
		if(EXON_USE_VALUE_ARRAY_INDEX)
		{
			
			sprintf(table_fn, "%s.%02d.%c.array", index_prefix, tabno, ginp->space_type==GENE_SPACE_COLOR?'c':'b');
			stat_ret = stat(table_fn, &filestat);
			if (stat_ret !=0)
			{
				SUBREADprintf("The index does not contain any raw base data which is required in detecting junctions. Please use the -b ooption while building the index.\n");
				return -1;
			}
			if (tabno>0)gvindex_destory(&value_array_index);

			if(gvindex_load(&value_array_index,table_fn)) return -1;
		}
		processed_reads = 0;
		last_table = tabno;

// Run the search algorithm on a part of the index
		if(EXON_ALL_THREADS <2)
			run_exon_search(overall_exon_bed, pos_table,connection_table,  my_table, &value_array_index, tabno, halves_record, ginp , ginp2, index_prefix, &processed_reads, base_number, all_tables, NULL /*the data lock is null*/, 0  /*I'm the 0-th thread*/, section_length, tolerable_scan);
		else
		{
			int i; 
			struct gene_thread_data_transport data_param;
			pthread_t runners [EXON_ALL_THREADS];
			pthread_spinlock_t  data_lock;
			pthread_spinlock_t  init_lock;

			data_param.my_table = my_table;
			data_param.my_value_array_index = &value_array_index;
			data_param.table_no = tabno;
			data_param.all_tables = all_tables;
			data_param.index_prefix = index_prefix;
			data_param.halves_record = halves_record;
			data_param.base_number = base_number;
			data_param.all_threads = EXON_ALL_THREADS;
			data_param.input_data_lock = &data_lock;
			data_param.init_lock = &init_lock;
			data_param.processed_reads = &processed_reads;
			data_param.section_length = section_length;
			data_param.ginp = ginp;
			data_param.tolerable_scan = tolerable_scan;
			data_param.ginp2 = ginp2;
			data_param.bed_table = overall_exon_bed;
			data_param.pos_table = pos_table;
			data_param.connection_table = connection_table;

			pthread_spin_init(&data_lock, PTHREAD_PROCESS_PRIVATE);
			pthread_spin_init(&init_lock, PTHREAD_PROCESS_PRIVATE);
			pthread_spin_lock(&init_lock);
			for (i=0; i< EXON_ALL_THREADS; i++)
			{
				data_param.this_thread = i;
				pthread_create(runners+i, NULL, run_exon_search_thread, &data_param);
				pthread_spin_lock(&init_lock);
			}

			for (i=0; i< EXON_ALL_THREADS; i++)
				pthread_join(*(runners+i), NULL);
			pthread_spin_destroy(&data_lock);
			pthread_spin_destroy(&init_lock);
		}
		tabno ++;

		if (section_length < 1)
			section_length = processed_reads;

		gehash_destory_fast(my_table);

		last_fp_pos = ftello(ginp -> input_fp);


		fseeko(ginp -> input_fp, current_fp, SEEK_SET);

	}

	tabno = 0;
	if(1)	// do two iterations
	{
		int i;
		gene_value_index_t * value_array_index_set[99];

		for(i=0; i<all_tables-1; i++)	
		{
			value_array_index_set[i]=(gene_value_index_t *)malloc(sizeof(gene_value_index_t));

			sprintf(table_fn, "%s.%02d.%c.array", index_prefix, i, ginp->space_type==GENE_SPACE_COLOR?'c':'b');
			if(gvindex_load(value_array_index_set[i],table_fn)) return -1;
		}
		value_array_index_set[i] = &value_array_index;

		if(EXON_IS_STEP1_RUN)
		{
			feed_exonbed_maybe_threads(overall_exon_bed, pos_table, connection_table, halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads,value_array_index_set, all_tables, tolerable_scan);

			remove_neighbours(overall_exon_bed ,  connection_table, pos_table);
		}

		
		if(EXON_IS_STEP2_RUN)
		{
			fseeko(ginp -> input_fp, current_fp, SEEK_SET);

			explorer_junc_exonbed_maybe_threads(overall_exon_bed,pos_table, connection_table, halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads, value_array_index_set, all_tables , tolerable_scan);
		}

		for(i=0; i<all_tables; i++)
		{
			gvindex_destory(value_array_index_set[i]);
			if(i< all_tables-1) free(value_array_index_set[i]);
		}

		
		fseeko(ginp -> input_fp, current_fp, SEEK_SET);
	}




	if (out_fp  && tolerable_scan && EXON_IS_STEP2_RUN)
	{
		if(ginp2)
			print_exon_res_paired(NULL , halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads);
		else
			print_exon_res_single(NULL , halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads);
	}


	fseeko(ginp -> input_fp, last_fp_pos, SEEK_SET);


	return processed_reads;
	
}


int run_exon_search_index(gene_input_t * ginp, gene_input_t * ginp2, char * index_prefix, halves_record_t * halves_record, FILE * out_fp, unsigned long long int base_number, int all_tables, unsigned long long int *succeed_reads, HashTable * overall_exon_bed, HashTable * pos_table, HashTable * connection_table)
{

	long long int tmp_pos;
	int ret = 0;

	
	if(!EXON_NO_TOLERABLE_SCAN)
	{
		tmp_pos = ftello(ginp -> input_fp);

		ret =run_exon_search_index_tolerable(ginp, ginp2, index_prefix, halves_record, out_fp, base_number, all_tables, succeed_reads,  overall_exon_bed, pos_table,  connection_table,0);
		clear_processed_marks(halves_record);

		fseeko(ginp -> input_fp, tmp_pos, SEEK_SET);
	}
	ret = run_exon_search_index_tolerable(ginp, ginp2, index_prefix, halves_record, out_fp, base_number, all_tables, succeed_reads,  overall_exon_bed, pos_table,  connection_table,1);
	return ret;
}


void exon_usage(char * execname)
{
	SUBREADprintf("Version %s\n\n", SUBREAD_VERSION);
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs(" ./subjunc [options] -i <index_name> -r <input> -o <output>");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("    -i --index     <index>  base name of the index.");
	SUBREADputs("");
	SUBREADputs("    -r --read      <input>  name of the input file(FASTQ/FASTA format). Both ");
	SUBREADputs("                            base-space and color-space read data are supported. ");
	SUBREADputs("                            For paired-end reads, this gives the first read file");
	SUBREADputs("                            and the other read file should be specified using");
	SUBREADputs("                            the -R option.");
	SUBREADputs("");
	SUBREADputs("    -o --output    <output> name of the output file(SAM format).");
	SUBREADputs("");
	SUBREADputs("Optional general arguments:");
	SUBREADputs("");
	SUBREADputs("    -n --subreads  <int>    number of selected subreads, 14 by default.");
	SUBREADputs("");
	SUBREADputs("       --singleSAM <input>  using the input file as a SAM file which includes");
	SUBREADputs("                            mapping results for single-end reads (e.g. 'subread-");
	SUBREADputs("                            align' output).");
	SUBREADputs("");
	SUBREADputs("       --pairedSAM <input>  using the input file as a SAM file which includes");
	SUBREADputs("                            mapping results for paired-end reads.");
	SUBREADputs("");
	SUBREADputs("    -T --threads   <int>    number of threads/CPUs used, 1 by default.");
	SUBREADputs("");
	SUBREADputs("    -I --indel     <int>    number of INDEL bases allowed, 5 by default.");
	SUBREADputs("");
	SUBREADputs("    -P --phred     <3:6>    the format of Phred scores used in input files, '3'");
	SUBREADputs("                            for phred+33 and '6' for phred+64. '3' by default.");
	SUBREADputs("");
	SUBREADputs("    -v                      displaying the version number.");
	SUBREADputs("");
	SUBREADputs("Optional arguments for paired-end reads:");
	SUBREADputs("");
	SUBREADputs("    -R --read2     <input>  name of the second input file from paired-end data. ");
	SUBREADputs("                            The program will then be switched to paired-end read");
	SUBREADputs("                            mapping mode.");
	SUBREADputs("");
	SUBREADputs("    -d --mindist   <int>    minimum fragment/template length, 50bp by default.");
	SUBREADputs("");
	SUBREADputs("    -D --maxdist   <int>    maximum fragment/template length, 600bp by default.");
	SUBREADputs("");
	SUBREADputs("    -S --order     <ff:fr:rf>  specifying if the first/second reads are forward");
	SUBREADputs("                            or reversed, 'fr' by default");
	SUBREADputs("");
	SUBREADputs("For more information about these arguments, please refer to the User Manual.");
	SUBREADputs("");
}

static struct option long_options[] =
{
	{"basewise", no_argument, &EXON_USE_VALUE_ARRAY_INDEX, 1},
	{"fusion", no_argument, &EXON_FUSION_DETECTION, 1},
	{"index", required_argument, 0, 'i'},
	{"read",  required_argument, 0, 'r'},
	{"read2", required_argument, 0, 'R'},
	{"bisulfite ", required_argument, 0, 'b'},
	{"indel", required_argument, 0, 'I'},
	{"nosam", required_argument, 0, 'A'},
	{"mindist", required_argument, 0, 'd'},
	{"maxdist", required_argument, 0, 'D'},
	{"minhalf", required_argument, 0, 'H'},
	{"subreads", required_argument, 0, 'n'},
	{"minmatch", required_argument, 0, 'm'},
	{"minmatch2", required_argument, 0, 'p'},
	{"quality", required_argument, 0, 'Q'},
	{"threads", required_argument, 0, 'T'},
	{"index-threshold", required_argument, 0, 'f'},
	{"output", required_argument, 0, 'o'},
	{"order",  required_argument,0, 'S'},
	{"halflen",  required_argument,0, 'L'},
	{"halfmatch",  required_argument,0, 'l'},
	{"singleSAM",  required_argument,0, '1'},
	{"pairedSAM",  required_argument,0, '2'},
	{"nofull",  no_argument, &EXON_NO_TOLERABLE_SCAN, 1},
	{"extending",  no_argument, &EXON_EXTENDING_SCAN, 1},
	{"direct",  no_argument, &EXON_DIRECT_READS, 1},
	{"jreadonly",  no_argument, &EXON_JUNCTION_READS_ONLY, 1},
	{0, 0, 0, 0}
};
void print_bed_table(HashTable * bed_table, char * out_fn, unsigned long long int * junction_number , unsigned long long int * support_number)
{
	int bucket;
	char fn2 [310];
	KeyValuePair * cursor;

	snprintf(fn2, 309, "%s.bed", out_fn);
	FILE * ofp = fopen(fn2, "w");

	for(bucket=0; bucket<bed_table -> numOfBuckets; bucket++)
	{
		cursor = bed_table -> bucketArray[bucket];
		while (1)
		{
			char * chro_name, * chro_name2;
			unsigned int chro_pos, chro_pos_big;
			if (!cursor) break;
			paired_exon_key * p = (paired_exon_key * ) cursor -> key;
			exon_junction_t *counter = (exon_junction_t*) cursor ->value;

			if(counter -> supporting_reads>0 || ! EXON_IS_STEP2_RUN)
			{
				locate_gene_position( p->small_key , &_global_offsets, &chro_name, &chro_pos);
				locate_gene_position( p->big_key , &_global_offsets, &chro_name2, &chro_pos_big);


				(*junction_number)++;
				(*support_number)+= (counter -> supporting_reads);
			//	fprintf(ofp,"%s\t%u\t%s\t%u\t%d\t%s,%s,\n", chro_name, chro_pos, chro_name2, chro_pos_big, counter -> supporting_reads, is_fusion, is_recovered);
				unsigned int feature_start = max(0,chro_pos-counter ->left_extend);
				unsigned int feature_end = chro_pos_big+counter ->right_extend;

				fprintf(ofp,"%s\t%u\t%u\tJUNC%08llu\t%d\t%c\t%u\t%u\t%d,0,%d\t2\t%d,%d\t0,%u\n", chro_name, feature_start,  feature_end,
													*junction_number, counter -> supporting_reads, counter -> strand?'-':'+',
													feature_start,  feature_end, counter -> strand?0:255, counter -> strand?255:0,
													counter ->left_extend, counter ->right_extend, feature_end-feature_start- counter ->right_extend);
			}

			free(counter);
			cursor = cursor->next;
		}
	}
	fclose(ofp);
}







void load_bed_table(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table, char * filename, gene_offset_t * offsets)
{
	FILE * fp = fopen(filename,"r");
	if(!fp)
	{
		SUBREADprintf("Warning: specified BED file '%s' is not found. Subjunc will run without preloaded BED table.\n", filename);
		return;
	}


	int step = 0;
	char chro [300];
	int chro_pnt = 0;
	unsigned int chro_pos1=0;
	unsigned int chro_pos2=0;
	int l_ext=0, r_ext = 0, lr_part = 0;

	int is_reverse_strand=0;
	int rows = 0;

	while(1)
	{
		if(feof(fp)) break;
		char nch = fgetc(fp);
		if(nch=='\r') continue;
		else if(nch == '\n')
		{

			if(step >0)
			{
				unsigned int linear_pos1;
				unsigned int linear_pos2;

				chro_pos1 += l_ext;
				chro_pos2 -= r_ext;


				linear_pos1 = linear_gene_position(offsets , chro, chro_pos1);
				linear_pos2 = linear_gene_position(offsets , chro, chro_pos2);

				if(linear_pos1==0xffffffff)
				{
					SUBREADprintf("Warning! Unknown chromosome: %s\n", chro);
				}
				else
				{

					linear_pos1 = get_grouped_position(pos_table, linear_pos1);
					linear_pos2 = get_grouped_position(pos_table, linear_pos2);

					put_connection_table(connection_table, linear_pos1, linear_pos2, is_reverse_strand, is_reverse_strand);
					put_connection_table(connection_table, linear_pos2, linear_pos1, is_reverse_strand, is_reverse_strand);

					exon_junction_t *search_res = (exon_junction_t *)malloc(sizeof(exon_junction_t));
					search_res->supporting_reads = 0;
					search_res->feed_supporting_reads = 0;
					search_res->strand = is_reverse_strand;
					search_res->left_extend=0;
					search_res->right_extend=0;

					paired_exon_key * new_key = (paired_exon_key *) malloc(sizeof(paired_exon_key));
					new_key->small_key = linear_pos1;
					new_key->big_key = linear_pos2;
					HashTablePut(bed_table, new_key, search_res);
					rows++;
				}
			}
			step = 0;
			chro_pnt = 0;
			chro_pos1=0;
			chro_pos2=0;
			l_ext = 0;
			r_ext = 0;
			lr_part =0;
			is_reverse_strand=0;
		}
		else if(nch=='\t')
		{
			step++;
		}
		else if(step == 0)
		{
			if((nch == '@' || nch == '#') &&(chro_pnt == 0))
			{
				chro_pnt=0;
				step = -999;
				continue;
			}
			else
			{
				chro[chro_pnt++]=nch;
				chro[chro_pnt]=0;
			}
		}
		else if(step == 1)
			chro_pos1 = chro_pos1*10 + (nch-'0');
		else if(step == 2)
			chro_pos2 = chro_pos2*10 + (nch-'0');
		else if(step == 5)
		{
			if(nch=='-') is_reverse_strand=1;
		}
		else if(step == 10)
		{
			if(nch==',')lr_part++;
			else
			{
				if(lr_part) r_ext = r_ext*10+(nch-'0');
				else l_ext = l_ext*10+(nch-'0');
			}
		}


	}
	fclose(fp);
	SUBREADprintf("%d rows have been loaded from the junction table.\n", rows);
	remove_neighbours(bed_table ,  connection_table, pos_table);

}



#ifdef MAKE_STANDALONE
int main(int argc,char ** argv)
#else
int main_junction(int argc,char ** argv)
#endif
{
	char read_file [300], read2_file [300];
	char output_file [300];
	char bed_input_file[300];

	char tmpfile[300];
	sprintf(tmpfile, "./subjunc-temp-sam-%06u-XXXXXX", getpid());

	char index_prefix [300];
	unsigned int all_reads, all_tables;
	halves_record_t halves_record;
	unsigned long long int processed_reads = 0, succeed_reads = 0;
	gene_input_t ginp, ginp2;
	HashTable * bed_index ;
	HashTable * pos_index ;
	HashTable * connection_index ;
	//gene_flat_t my_flat ;
	//create_flat_strip(&my_flat);

	int c;
	int option_index = 0;
	unsigned long long int junction_number=0, support_number=0; 

	bed_input_file[0]=0;

	EXON_IS_STEP1_RUN=1;
	EXON_IS_STEP2_RUN=1;
	EXON_MAJOR_HALF_VOTES_RATE = .2;

	IS_SAM_INPUT=0;
	TOTAL_SUBREADS = 14;
	ACCEPT_MINOR_SUBREADS = 1;
	EXON_MAX_METHYLATION_C_NUMBER = 0;
	INDEX_THRESHOLD = 24;
	read_file[0]=0;
	read2_file[0]=0;
	index_prefix[0]=0;
	output_file[0]=0;
	all_reads = 14*1024*1024;///128;
//	all_reads = 300000;
	int using_base_distance = 0;

	SUBREADprintf("\n");

	//	SUBREADprintf ("I'm the %d-th thread RRRKK\n", -9999);return 0;


	while ((c = getopt_long (argc, argv, "vxS:L:AH:d:D:n:m:p:f:P:R:r:i:l:o:T:Q:I:1:2:t:B:b:F?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'v':
				print_version_info();
				return 0;
				break;
			case 'A':
				REPORT_SAM_FILE = 0;
				break;
			case 'x':
				EXON_FUSION_DETECTION=1;
				break;
			case 'H':
				using_base_distance = 1;
				break;
			case 'S':
				EXON_FIRST_READ_REVERSE = optarg[0]=='r'?1:0;
				EXON_SECOND_READ_REVERSE = optarg[1]=='f'?0:1;
				break;
			case 'b':
				EXON_MAX_METHYLATION_C_NUMBER = atoi(optarg);
				break;
			case 'D':
				EXON_MAX_PAIRED_DISTANCE = atoi(optarg);
				break;
			case 'd':
				EXON_MIN_PAIRED_DISTANCE = atoi(optarg);
				break;
			case 'n':
				TOTAL_SUBREADS = atoi(optarg);
			//	SUBREADprintf(" === WARNING ===\n You cannot set subread numbers in exon detection. It is automatically set to the max number.");
				break;
			case 'f':
				INDEX_THRESHOLD  = atoi(optarg);
				break;
			case 'm':
				EXON_MAJOR_HALF_VOTES_RATE = atof(optarg);
				break;
				
			case 'T':
				EXON_ALL_THREADS = atoi(optarg);
				if(EXON_ALL_THREADS <1) EXON_ALL_THREADS=1;

				break;
			case 'r':
				if(IS_SAM_INPUT)
					SUBREADputs("\n!!WARNING!!\nYou should not specify both FASTQ/FASTA files and SAM files at same time!\nOnly the SAM file is used in the program.");
				else
					strncpy(read_file, optarg, 299);
				break;
			case 'R':
				strncpy(read2_file, optarg, 299);
				break;
			case 'i':
				strncpy(index_prefix, optarg, 299);
				break;
			case 'o':
				strncpy(output_file, optarg, 299);
				break;
			case 'I':
				EXON_INDEL_TOLERANCE = atoi(optarg);
				if( EXON_INDEL_TOLERANCE >5)EXON_INDEL_TOLERANCE=5;
				EXON_INDEL_TOLERANCE ++;
				break ;
			case 'Q':
				if(optarg[0]=='l')
					EXON_QUALITY_SCALE = QUALITY_SCALE_LINEAR;
				if(optarg[0]=='e')
					EXON_QUALITY_SCALE = QUALITY_SCALE_LOG;
				if(optarg[0]=='n')
					EXON_QUALITY_SCALE = QUALITY_SCALE_NONE;
				break;
			case 'P':
				if (optarg[0]=='3')
					EXON_FASTQ_FORMAT = FASTQ_PHRED33;
				else
					EXON_FASTQ_FORMAT = FASTQ_PHRED64;
				break;
			case 'p':
				ACCEPT_MINOR_SUBREADS = atoi(optarg);
				break;
			case 'L':
				EXON_MIN_HALF_LENGTH = atoi(optarg);
				break;
			//case 'l':
			//	EXON_HALF_MATCH_PERCENTAGE = atof(optarg);
			//	break;
			case '1':
				if(read_file[0] || IS_SAM_INPUT)
				{
					SUBREADputs("\n!!WARNING!!\nYou should not specify both FASTQ/FASTA files and SAM files at same time!\nOnly the SAM file is used in the program.");
				}
				IS_SAM_INPUT=1;

				strncpy(read_file, optarg, 299);
				EXON_FASTQ_FORMAT = FASTQ_PHRED33;
				break;
			case '2':
				if(read_file[0]|| IS_SAM_INPUT)
				{
					SUBREADputs("\n!!WARNING!!\nYou should not specify both FASTQ/FASTA files and SAM files at same time!\nOnly the SAM file is used in the program.");
				}
				IS_SAM_INPUT=2;
				strncpy(read_file, optarg, 299);
				EXON_FASTQ_FORMAT = FASTQ_PHRED33;
				strncpy(read2_file, optarg, 299);
				break;
			case 't':
				sprintf(tmpfile, "%s/subjunc-tem-sum-%06u-XXXXXX", optarg, getpid());
				break;
			case 'F':
				EXON_IS_STEP2_RUN = 0;
				REPORT_SAM_FILE = 0;
				break;
			case 'B':
				strcpy(bed_input_file, optarg);
				EXON_IS_STEP1_RUN = 0;
				break;
			case '?':
				return -1 ;
		}
	if (!read_file[0] || !index_prefix[0] || !output_file[0])
	{
		exon_usage(argv[0]);

		return -1 ;
	}

	if((!EXON_DIRECT_READS) && (!IS_SAM_INPUT))
	{
		int xx=0;
		char command[1000];
		char cwd[300];
		int is_successful = 0;

		is_successful = mkstemp(tmpfile);

		strcpy(cwd, argv[0]);
		for(xx=strlen(cwd);xx>=0; xx--)
		{
			if(cwd[xx]=='/')
			{
				cwd[xx]=0;
				break;
			}
		}

		if(xx>=0) strcat(cwd,"/");
		else cwd[0]=0;

		SUBREADputs("Call subread-align to map reads...");

		if(read2_file[0])
			sprintf(command, "%ssubread-align -J --allow-repeating  -T %d -i '%s' -r '%s' -R '%s' -o '%s' -P %d -d %d -D %d %s %s -I %d -b %d ",cwd, EXON_ALL_THREADS, index_prefix, read_file, read2_file, tmpfile, EXON_FASTQ_FORMAT == FASTQ_PHRED33?3:6, EXON_MIN_PAIRED_DISTANCE, EXON_MAX_PAIRED_DISTANCE, EXON_QUALITY_SCALE?"-Q":"", using_base_distance?"-H":"", EXON_INDEL_TOLERANCE-1, EXON_MAX_METHYLATION_C_NUMBER);
		else
			sprintf(command, "%ssubread-align  -J --allow-repeating  -T %d -i '%s' -r '%s' -o '%s' -P %d %s %s -I %d -b %d",cwd, EXON_ALL_THREADS, index_prefix, read_file, tmpfile, EXON_FASTQ_FORMAT == FASTQ_PHRED33?3:6 , EXON_QUALITY_SCALE?"-Q":"", using_base_distance?"-H":"", EXON_INDEL_TOLERANCE-1, EXON_MAX_METHYLATION_C_NUMBER);
		//SUBREADputs(command);
		is_successful = system(command);

		strcpy(read_file, tmpfile);
		IS_SAM_INPUT=1;

		if(read2_file[0])
		{
			strcpy(read2_file, tmpfile);
			IS_SAM_INPUT=2;
		}
		EXON_FASTQ_FORMAT = FASTQ_PHRED33;
		reads_density = guess_reads_density(read_file, IS_SAM_INPUT);
		if(reads_density<0)
		{
			SUBREADputs("Subjunc is terminated because subread could not map the reads.");
			return -1;
		}
	}
	else
	{
		tmpfile[0]=0;
		reads_density = guess_reads_density(read_file, IS_SAM_INPUT);
  
		if(reads_density<0)
			SUBREADprintf("Input file '%s' is not found or is in an incorrect format.\n", read_file);
	}

	SUBREADputs("Detect exon-exon junctions and map reads...");

	if(IS_SAM_INPUT==0)
	{
		if(geinput_open(read_file, &ginp))
			return -1;
	}
	else if(IS_SAM_INPUT==1)
	{
		if(geinput_open_sam(read_file, &ginp,0))
			return -1;
	}
	else if(IS_SAM_INPUT==2)
	{
		if(geinput_open_sam(read_file, &ginp,0))
			return -1;
		//if(geinput_open_sam(read_file, &ginp2,2))
		//	return -1;
	}

	if(ginp.space_type==GENE_SPACE_COLOR)
	{
		SUBREADprintf("Subjunc currently does not support color-space junction detection.\n");
		return -1;
	}

	for(all_tables=0; ; all_tables++)
	{
		char fn[300];
		struct stat fstatbuf;
		sprintf(fn, "%s.%02d.%c.tab", index_prefix, all_tables, ginp.space_type==GENE_SPACE_COLOR?'c':'b');
		if(stat(fn, &fstatbuf))
			break;
	}

	if(all_tables==0)
	{
		SUBREADprintf("Unable to open the index files in the %s space.\n",  ginp.space_type==GENE_SPACE_COLOR?"color":"base");
		return -1;
	}

	FILE * out_fp = NULL;
	if(REPORT_SAM_FILE)
		out_fp= fopen(output_file, "w");
	if (REPORT_SAM_FILE && !out_fp)
	{
		SUBREADprintf("Unable to open the output file at '%s'.\n", output_file);
		return -1;
	}


	//SUBREADprintf("Number of subreads selected for each read=%d\n", TOTAL_SUBREADS);
	SUBREADprintf("Threshold on number of subreads for a successful mapping=%f\n", EXON_MAJOR_HALF_VOTES_RATE);
	SUBREADprintf("Number of threads=%d\n", EXON_ALL_THREADS);
	if (EXON_INDEL_TOLERANCE)
		SUBREADprintf("Tolerance for Indel=%d\n", EXON_INDEL_TOLERANCE-1);
	if (EXON_QUALITY_SCALE==QUALITY_SCALE_LINEAR)
		SUBREADputs("Quality scale=linear\n\n");
	else if (EXON_QUALITY_SCALE==QUALITY_SCALE_LOG)
		SUBREADputs("Quality scale=exponential\n\n");
	else 	SUBREADputs("\n");


	if (read2_file[0] || IS_SAM_INPUT==2)
	{
		if (EXON_MAX_PAIRED_DISTANCE <= EXON_MIN_PAIRED_DISTANCE)
		{
			SUBREADprintf ("The value of the '-D' option must be greater than that of the '-d' option. \n");
			return -1;
		}

		SUBREADprintf ("Performing paired-end alignment:\n");
		SUBREADprintf ("Maximum fragment length=%d\n", EXON_MAX_PAIRED_DISTANCE);
		SUBREADprintf ("Minimum fragment length=%d\n", EXON_MIN_PAIRED_DISTANCE);
		SUBREADprintf ("Threshold on number of subreads for a successful mapping (the minor end in the pair)=%d\n", ACCEPT_MINOR_SUBREADS);
		SUBREADprintf ("The directions of the two input files are: %s, %s\n\n", EXON_FIRST_READ_REVERSE?"reversed":"forward", EXON_SECOND_READ_REVERSE?"reversed":"forward");
	}

#ifdef REPORT_ALL_THE_BEST
	SUBREADprintf("***** WARNING: the REPORT_ALL_THE_BEST switch is turned on. You need an extra 1 GBytes of RAM space for saving the temporary results. *****\n");
#endif

	if(init_halves_record(&halves_record, all_reads)) return -1;
	if(read2_file[0])
		all_reads/=2;
		
	SUBREADfflush(stdout);

	begin_ftime = miltime();

	load_offsets (&_global_offsets, index_prefix);

	bed_index = HashTableCreate(399997);
	pos_index = HashTableCreate(399997);
	connection_index = HashTableCreate(399997);

	HashTableSetKeyComparisonFunction(connection_index, pointercmp_forpos);
	HashTableSetHashFunction(connection_index, pointerHashFunction_forpos);

	HashTableSetKeyComparisonFunction(pos_index, pointercmp_forpos);
	HashTableSetHashFunction(pos_index, pointerHashFunction_forpos);
	
	HashTableSetKeyComparisonFunction(bed_index, pointercmp_forbed);
	HashTableSetHashFunction(bed_index, pointerHashFunction_forbed);

	if(bed_input_file[0])
	{
		load_bed_table(bed_index, pos_index, connection_index, bed_input_file, &_global_offsets);
	}

	while (1)
	{
		char inbuff[1201];
		int i;

		for(i=0; i<EXON_ALL_THREADS; i++)
		{
			thread_block_locations[i]=-1;
		}

		int new_processed_reads = run_exon_search_index(&ginp, read2_file[0] ? (&ginp2):NULL, index_prefix, &halves_record, out_fp, processed_reads, all_tables, &succeed_reads, bed_index, pos_index, connection_index);
		if(new_processed_reads<0)break;

		processed_reads += new_processed_reads;
		clear_halve_record(&halves_record);

		// test if there no anyreads remaining
		unsigned long long int current_fp = ftello(ginp.input_fp);
		int rl = geinput_next_read(&ginp, NULL, inbuff, NULL);
		//SUBREADprintf("RL=%d\nfp=%llu\npr=%u",rl, current_fp, processed_reads);
		if (rl<0)
			break;
		fseeko(ginp.input_fp, current_fp, SEEK_SET);
	}


	geinput_close(&ginp);
	print_bed_table(bed_index, output_file, &junction_number, &support_number);
	HashTableDestroy(bed_index);
	HashTableDestroy(pos_index);
	HashTableDestroy(connection_index);


	if(IS_DEBUG)
		SUBREADprintf("@LOG THE END. \n");
	else
		SUBREADprintf("\n\n %llu %s were processed in %.1f seconds.\n There are %llu junction pairs found, supported by %llu reads.\n\n", processed_reads, read2_file[0]?"fragments":"reads", miltime()-begin_ftime, junction_number, support_number );

	if(out_fp)
		fclose(out_fp);
	SUBREADprintf("\n\nCompleted successfully.\n");

	if(tmpfile[0])
		unlink(tmpfile);
	destory_halves_record(&halves_record);
	destroy_offsets (&_global_offsets);
	return 0;
}
