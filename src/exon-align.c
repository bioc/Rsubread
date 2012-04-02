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
int TOTAL_SUBREADS;
float EXON_MAJOR_HALF_VOTES_RATE = 0.2;
float EXON_MIN_HALF_VOTES_RATE = 0.15;


int MIN_VOTE2_TMP_VAR = 1;
int EXON_SUBREAD_GAP = 6;
int EXON_LARGE_WINDOW = 60;
int EXON_LONG_READ_LENGTH = 120;


int ACCEPT_MINOR_SUBREADS;
int INDEX_THRESHOLD;
int EXON_MAX_PAIRED_DISTANCE = 600;
int EXON_MIN_PAIRED_DISTANCE = 50;
int EXON_INDEL_TOLERANCE = 4;
int EXON_QUALITY_SCALE = QUALITY_SCALE_NONE;
int EXON_USE_VALUE_ARRAY_INDEX = 1;
int EXON_FIRST_READ_REVERSE = 0;
int EXON_SECOND_READ_REVERSE = 1;
int EXON_NUMBER_OF_ANCHORS_PAIRED = 50;
int EXON_DIRECT_READS = 0;
int EXON_NO_TOLERABLE_SCAN = 1;
int EXON_MAX_CIGAR_LEN = 48;

int EXON_FASTQ_FORMAT = FASTQ_PHRED33;

int IS_SAM_INPUT = 0;
int EXON_MIN_HALF_LENGTH = 0;
float EXON_HALF_MATCH_PERCENTAGE =.7f;
double reads_density;

int EXON_FUSION_DETECTION = 0;

#define EXON_DONOR_TEST_WINDOW 17
#define EXON_GROUPING_SIZE 1
#define MAX_EXON_CONNECTIONS 10
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

void init_halves_record(halves_record_t* halves_record, int items)
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

	memset(halves_record -> best_vote1_list, 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> best_vote2_list, 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> cigar_string_buffer , 0 , sizeof(char)*items*EXON_MAX_CIGAR_LEN);
	memset(halves_record -> half_marks_list, 0, sizeof(short)*halves_record -> max_len);
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
	memset(halves_record -> cigar_string_buffer , 0 , sizeof(char)*halves_record -> max_len*EXON_MAX_CIGAR_LEN);
}

void add_best_matching_halves(halves_record_t * halves_record, unsigned int best_pos1, unsigned int best_pos2, unsigned char best_vote1 , unsigned char best_vote2, char is_abnormal, char is_reversed, short splice_point, short half_marks, int read_number, int best_pos_start, int best_pos_end, short read_coverage_start, short read_coverage_end, char indel_in_p1, char indel_in_p2)
{
	//printf ("\nADDED_BEST_HALVES #%d V1=%d  V2=%d;  oldv = %d\n", read_number, best_vote1, best_vote2, halves_record -> best_vote1_list[read_number] );
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

		//printf("ADD_BEST=%d, v1=%d, v2=%d\n", splice_point, best_vote1, best_vote2);
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
			int long_overlap = 0;
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
			//printf("MAX: %d-%d   OTHER %d-%d    COV=%d   OVLP=%d\n", max_start, max_end, vote->coverage_start[i][j], vote->coverage_end[i][j], coverage_len, overlapped_len);



			if(overlapped_len >=14)
				continue;
			else if (overlapped_len >=7)
				long_overlap = IS_LONG_OVERLAP | IS_SHORT_OVERLAP;
			else if (overlapped_len >=3)
				long_overlap = IS_SHORT_OVERLAP;


			long long int dist = vote->pos[i][j];
			dist -= max_pos;

			//printf ("D=%lld\n", abs(dist));
			if (abs(dist)<6)
				continue;

			//int support_r1 = (int)(max(EXON_MIN_HALF_VOTES+0.01, (n1*1.-15)/3*accept_rate*1.0001));
			//int support_r2 = (int)(max(EXON_MIN_HALF_VOTES+0.01, (n2*1.-15)/3*accept_rate*1.0001));

			int support_r1 = (int) (/*(rl -15)/EXON_SUBREAD_GAP*/ TOTAL_SUBREADS * EXON_MAJOR_HALF_VOTES_RATE); 
			int support_r2 = 1;//(int) ((rl -15)/3 * EXON_MIN_HALF_VOTES_RATE); 


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
			if( is_reversed && ((max_pos > vote->pos[i][j]) + (vote->coverage_start[i][j] < max_start) != 1))is_fusion = 1;
			if((! is_reversed) && ((max_pos > vote->pos[i][j]) + (vote->coverage_start[i][j] > max_start) != 1)) is_fusion = 1;

			if(abs(dist) > 500000 || chro_name != best_chro_name)
				is_fusion = 1;

			if (is_fusion && !EXON_FUSION_DETECTION) continue;

			int test_vote_value ;
	//		test_vote_value = 8888888 +  vote->votes[i][j]* 4000 - abs(dist);
			test_vote_value = 8888888 +  vote->votes[i][j]* 1000000 - abs(dist);
			if (hint_pos>=0)
			{
				long long int hint_dist = hint_pos;
				hint_dist -= vote->pos[i][j];
				if (abs (hint_dist) < 1000000)
					test_vote_value += 100;
				if (abs (hint_dist) < 10000)
					test_vote_value += 100;
			}

			//printf ("\n%u and %u : V1 %d >= S1 %d\t\tV2 %d >= S2 %d  ;  TVV=%d <=SEL=%d; BRK=%d\n", vote->max_position, vote->pos[i][j] , vote->max_vote, support_r1,  vote->votes[i][j], support_r2, test_vote_value , selected_max_votes, ((vote->coverage_start[i][j] < vote->max_coverage_start)? (vote->coverage_end[i][j]):(vote->max_coverage_end)) + ((vote->coverage_start[i][j] < vote->max_coverage_start)? (vote->max_coverage_start):(vote->coverage_start[i][j]))/2);

			if (test_vote_value<selected_max_votes)continue;
			// Conditions of order of R3 and R5
			*half_marks &= ~IS_REVERSED_HALVES;
			if (vote->coverage_start[i][j] < max_start && (((max_pos < vote->pos[i][j]) && !is_reversed) || ((max_pos > vote->pos[i][j]) && is_reversed) ) )
				*half_marks |= IS_REVERSED_HALVES;
			if (vote->coverage_start[i][j] >= max_end  &&  (((max_pos > vote->pos[i][j]) && !is_reversed) || ((max_pos < vote->pos[i][j]) && is_reversed) ) )
				*half_marks |= IS_REVERSED_HALVES;

			if (vote->coverage_start[i][j] < max_start)
			{
				*half_marks = (*half_marks) & ~IS_R1_CLOSE_TO_5;
			}
			else
			{
				*half_marks |= IS_R1_CLOSE_TO_5;
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

			*half_marks = (*half_marks) & ~(IS_LONG_OVERLAP|IS_SHORT_OVERLAP);
			if (long_overlap)
				*half_marks |= long_overlap;

			if (is_fusion)
				*half_marks = (*half_marks)    | IS_FUSION;
			else
				*half_marks = (*half_marks) & ~( IS_FUSION);
	

			selected_max_votes = test_vote_value; 

		}
	*best_select_max_votes = selected_max_votes ;
	return best_splicing_point;
}



int select_best_matching_halves(gene_vote_t * vote, unsigned int * best_pos1, unsigned int * best_pos2, int * best_vote1, int * best_vote2, char * is_abnormal, short * half_marks, int * is_reversed_halves, float accept_rate, int read_len, long long int hint_pos, int tolerable_bases, short * read_coverage_start, short * read_coverage_end, char * indel_in_p1, char * indel_in_p2 , int * max_cover_start, int * max_cover_end, int rl)
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
	int encountered = 0;
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
			if(vote->votes[i][j] >=  vote->max_vote -1)encountered++;

	if (encountered>1)
	{
		return 0;
	}

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			test_select_votes = 0;
			if(vote->votes[i][j] < vote->max_vote) continue;
			//if(( vote-> coverage_start[i][j] <15 )&& (vote-> coverage_end[i][j] > rl -15 )) continue;
			//if(vote->pos[i][j] != vote->max_position) continue;

			ret = select_best_matching_halves_maxone(vote, &tmp_best_pos1, &tmp_best_pos2, &tmp_best_vote1, &tmp_best_vote2,  &tmp_is_abnormal,half_marks, &tmp_is_reversed_halves, accept_rate, read_len, hint_pos,  tolerable_bases, &tmp_read_coverage_start, &tmp_read_coverage_end, &tmp_indel_in_p1, &tmp_indel_in_p2, vote -> pos[i][j],  vote->votes[i][j], vote-> coverage_start[i][j], vote-> coverage_end[i][j],  vote-> masks[i][j], vote->indel_recorder[i][j], &test_select_votes, rl);
			test_select_votes += vote->votes[i][j]*1000000;
			//printf("TSV=%d\n",test_select_votes);

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

				* max_cover_start = vote-> coverage_start[i][j];
				* max_cover_end = vote-> coverage_end[i][j];
				best_ret = ret;
			}
		}
	return best_ret;
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
int test_donor(char *read, int read_len, unsigned int pos1, unsigned int pos2, int guess_break_point, char negative_strand, int test_range, char is_soft_condition, int EXON_INDEL_TOLERANCE, int* real_break_point, gene_value_index_t * my_value_array_index, int indel_offset1, int indel_offset2, int is_reversed, int space_type, int confidence, int * best_donor_score)
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
		//printf("C1=%s @%u, C2=%s @%u\n",h1_2ch, pos1 + bps_pos_x, h2_2ch,pos2 - 2 + indel_offset + bps_pos_x);
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
				//printf("DONOR_CONF_LEN=%d\n", donar_conf_len);

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
					printf("DONOR TEST STR=%s, %s ; pos=%d    %d %d ; M=%d %d ; X=%d %d\n", h1_2ch, h2_2ch, bps_pos_x, indel_offset1, indel_offset2, m1, m2, x1, x2);
				}
				#endif
	
				int threshold = 4;
				if (paired_score == 1)
					threshold = 4;

				if (m1 >= donar_conf_len - 1 && m2>=donar_conf_len - 1)
					if(x1<donar_conf_len - threshold && x2<donar_conf_len - threshold)
					{
						int score = 1000 * paired_score - (x1 + x2) + (m1+ m2) ;
						if (min_x < score)
						{
							min_x = score;
							best_break = bps_pos_x;
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
					printf("SELECRED!!!_BREAKPOINT=%d, RAW POS=%u,%u, R=%s\n",  best_break, pos1 , pos2, read);
				#endif
		//printf ("FINAL BREAK: %d   ; REV = %d\n ", best_break, is_reversed);
		*real_break_point = best_break;
		return 1;
	}
	else
	{
				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
					printf("KILLED!!!_BREAKPOINT=%d, R=%s\n",  best_break+ pos1, read);
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
//	else	printf("REUSE:%d\n", grouped_pos);
	return grouped_pos;
}

#define EXON_EXPLAIN_DEPTH 5
//#define DEBUG

void junction_tree_b_explorer(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table, char * read , int rl, int full_rl ,unsigned int read_tail_pos, int number_of_piece, unsigned int * cigar_recorder, int total_matched_bases, int * result_total_matched_bases,unsigned int * result_cigar_recorder, int *result_number_pieces, gene_value_index_t * my_value_array_index, gene_input_t * ginp,int subread_votes_cover_start, char * quality_str, int quality_format, float match_score)
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


	int matched_bases =  match_indel_chro_to_back(read, my_value_array_index, read_tail_pos, rl,&indels, &indel_point, EXON_INDEL_TOLERANCE);
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
			
			matched_bases = match_indel_chro_to_back(read+rl-test_piece_len, my_value_array_index, read_tail_pos - test_piece_len, test_piece_len,&indels, &indel_point, EXON_INDEL_TOLERANCE);

			#ifdef DEBUG
			printf("BSEARCH: INDEL=%d\n",indels);
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
		printf("KKK1: %d, %d, %u\n", max_indel_point + cigar_recorder [number_of_piece*4+1], max_indels, cigar_recorder [number_of_piece * 4+3] );
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

			#ifdef DEBUG
			printf ("DIG-IN B-SEARCH REMAIN_LEN=%d, POS_TAIL=%u\n", rl + max_indels- max_piece_len, connect_to_iii);
			#endif
			junction_tree_b_explorer(bed_table, pos_table, connection_table, read, rl + max_indels - max_piece_len, full_rl, connect_to_iii , 1+number_of_piece , cigar_recorder, total_matched_bases + max_matched_bases, result_total_matched_bases , result_cigar_recorder, result_number_pieces,  my_value_array_index, ginp, subread_votes_cover_start, quality_str , quality_format, match_score);
			digged = 1;
		}
	}
		// the "read_tail_pos" is the exon that matches through this read; it has to be well matched to accept the read.
	//if (!digged)
	{
		indels = 0;
		indel_point =0;


		int max_piece_quality_good = 1;

		int matched_bases = match_indel_chro_to_back(read, my_value_array_index, read_tail_pos - rl, rl, &indels, &indel_point, EXON_INDEL_TOLERANCE);
		if ( matched_bases > -1965 * rl )
		{
			//put this explaination to the result, if the matched bases is the max.
			if (matched_bases +total_matched_bases > (*result_total_matched_bases))
			{
				cigar_recorder [number_of_piece * 4 + 1] = 0;
				cigar_recorder [number_of_piece * 4 + 3] = ((indel_point + cigar_recorder [number_of_piece*4+1]) << 16) + (0xffff&(indels << 4)) + (max_piece_quality_good << 3);

				memcpy(result_cigar_recorder, cigar_recorder, 138);
				(*result_total_matched_bases) = matched_bases + total_matched_bases;
				(*result_number_pieces) = number_of_piece+1;
			}
		}
	}
}


void junction_tree_f_explorer(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table, char * read , int rl, int full_rl ,unsigned int read_head_pos, int number_of_piece, unsigned int * cigar_recorder, int total_matched_bases, int * result_total_matched_bases,unsigned int * result_cigar_recorder, int *result_number_pieces, gene_value_index_t * my_value_array_index, gene_input_t * ginp, int next_head_modify, char * quality_str, int quality_format, float match_score)
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



	int matched_bases =  match_indel_chro_to_front(read, my_value_array_index, read_head_pos, rl,&indels, &indel_point, EXON_INDEL_TOLERANCE);
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

		//printf("TESTPOS=%u > %u\n", grouped_pos, read_head_pos);

		if (grouped_pos && grouped_pos >= read_head_pos)
		{
			int test_piece_len = grouped_pos - read_head_pos;
			if(test_piece_len < 8) continue;
			int matched_bases;

			matched_bases = match_indel_chro_to_front(read, my_value_array_index, read_head_pos, test_piece_len, &indels, &indel_point, EXON_INDEL_TOLERANCE);


			#ifdef DEBUG
			printf ("TEST match = %d ; piece_len = %d; GRP_POS=%u; \n", matched_bases, test_piece_len, grouped_pos);
			printf("FSEARCH: INDEL=%d\n",indels);
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
		printf("KKK3: %d, %d, %u\n", max_indel_point + cigar_recorder [number_of_piece*4+1], max_indels, cigar_recorder [number_of_piece * 4+3] );
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
#ifdef DEBUG
			printf ("DIG-IN:match/piece_len = %d/%d, rl=%d head_pos=%u  READ=%s\n", max_matched_bases, max_piece_len, rl - max_piece_len, connect_to_iii, read+max_piece_len);
#endif
			//printf("DIGIN: TOTAL=%d NEWADD=%d\n", total_matched_bases, max_matched_bases);
			//printf("F-MAXINDEL=%d ; NEXT_HEAD_MODIFY=%d ; #PIECES=%d\n", max_indels, next_head_modify, number_of_piece );
			junction_tree_f_explorer(bed_table, pos_table, connection_table, read+max_piece_len-max_indels, rl - max_piece_len + max_indels, full_rl , connect_to_iii , 1+number_of_piece , cigar_recorder, total_matched_bases + max_matched_bases, result_total_matched_bases , result_cigar_recorder,result_number_pieces, my_value_array_index, ginp, 0 , quality_str + max_piece_len-max_indels , quality_format, match_score);
			digged = 1;
		}
	}
		// the "read_head_pos" is the exon that matches through this read; it has to be well matched to accept the read.
	//if(!digged)
	{
		indels = 0;
		indel_point =0;



		int max_piece_quality_good = 1;

		int matched_bases = match_indel_chro_to_front(read, my_value_array_index, read_head_pos, rl, &indels, &indel_point, EXON_INDEL_TOLERANCE);

		if ( matched_bases > -1965 * rl )
		{
			//put this explaination to the result, if the matched bases is the max.
			if (matched_bases + total_matched_bases > (*result_total_matched_bases))
			{
				cigar_recorder [number_of_piece * 4 + 2] = full_rl;
				cigar_recorder [number_of_piece * 4 + 3] = ((indel_point + cigar_recorder [number_of_piece*4+1]) << 16) +  ((0xfff&indels) << 4)  + (max_piece_quality_good << 3);

#ifdef DEBUG
				printf("KKK4: %d+%d, %d, %u\n", indel_point,cigar_recorder [number_of_piece*4+1], indels, cigar_recorder [number_of_piece * 4+3] );
#endif

				memcpy(result_cigar_recorder, cigar_recorder, 138);
				(*result_total_matched_bases) = matched_bases + total_matched_bases;
				(*result_number_pieces) = number_of_piece+1;
			}
		}
	}
}

void explorer_junc_exonbed(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table, halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads, gene_value_index_t * my_value_array_index, int first_index, int last_index, int tolerable_scan)
{
	int i=0, ic=0, j, result_head_modify=-1;
	//return;
	while (1)
	{
		char nameb[1201], inb[1201], qualityb[1201];
		int rl;
		int is_settle = 1;
		if(i >= processed_reads*(1+(ginp2!=NULL)))break;


		if(i % 10000 ==0 && i>1)
			print_text_scrolling_bar("Second Iteration", i*1./(1+(ginp2!=NULL))/processed_reads, 80, &ic);

		unsigned int pos = halves_record -> best_pos1_list[i];

		if (ginp2 && (i % 2))
			rl = geinput_next_read(ginp2, nameb, inb, qualityb);
		else
			rl = geinput_next_read(ginp, nameb, inb, qualityb);
		if (rl<0){
			break;
		}

		char negative_strand = halves_record -> is_reversed_list[i];

		if ((ginp2 && (i % 2)) + negative_strand == 1)
		{
			reverse_read(inb, rl, ginp->space_type);
			reverse_quality(qualityb, rl);
		}
		//printf("INPUT=%d , %s\n", ginp->file_type, nameb);

#ifdef TEST_TARGET
		if(memcmp(inb, TEST_TARGET, 15)==0)
			printf("P0=%u, P1=%d, P1_INDEL=%d, P2_INDEL=%d, BPOINT=%d\n", pos, halves_record -> best_pos2_list[i] , halves_record -> splice_point_offset_1[i], halves_record -> splice_point_offset_2[i],  halves_record -> splice_point_list[i]);
#endif


		if(halves_record -> best_vote1_list[i]< /*(rl - 15)/EXON_SUBREAD_GAP*/ TOTAL_SUBREADS * EXON_MAJOR_HALF_VOTES_RATE)
		{
			i++; continue;
		}

		if (!first_index)
			if (pos< my_value_array_index -> start_base_offset + 1000)
			{
				i++;
				continue;
			}
		if (!last_index)
			if (pos > my_value_array_index -> start_base_offset+ my_value_array_index ->length - 1000)
			{
				i++; continue;
			}
		if (IS_PROCESSED_READ_R2 & halves_record -> half_marks_list[i])
		{
			i++; continue;
		}

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

#ifdef DEBUG
		printf ("\n%s %s\nB-SEARCH RANGE [%d - %d] ; POS=%u ; P1_Indel=%d ; P2_Indel = %d\n",nameb, inb ,f_search_head , b_search_tail,pos + b_search_tail ,  halves_record ->indel_in_piece1[i],  halves_record ->indel_in_piece2[i]);
#endif

		junction_tree_b_explorer(bed_table, pos_table, connection_table, inb , b_search_tail, rl , pos + b_search_tail + halves_record ->indel_in_piece1[i] , 0, explain_buff, 0, &result_total_matched_bases , explain_result,  &result_number_pieces, my_value_array_index, ginp, f_search_head, qualityb, EXON_FASTQ_FORMAT, match_score);
	//	printf(" POS BBB RL=%d, MATCHED_TOTAL=%d LEN=%d\n", rl, result_total_matched_bases, b_search_tail + halves_record ->indel_in_piece1[i]);

#ifdef DEBUG

		int xx;
		printf("BSEARCH MATCH BASES=%d/%d; PIECES=%d\n", result_total_matched_bases,b_search_tail, result_number_pieces );
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

			printf ("BS-RES: [%d ~ %d  %c %s,%d] \n", piece_start, piece_end, ((ginp2 && (i % 2)) + negative_strand == 1)?'~':'@', chro_name , chro_pos);
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
				//			printf("indels_in_section=%d\n",indels_in_section);
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


#ifdef DEBUG
			printf ("\n%s %s\nF-SEARCH RANGE [%d - END] ; POS=%u ; P1_Indel=%d ; P2_Indel = %d\n",nameb, inb ,f_search_head ,explain_result[(total_pieces-1)*4],  halves_record ->indel_in_piece1[i],  halves_record ->indel_in_piece2[i]);
#endif
		}
		{
			result_number_pieces = 0,  result_total_matched_bases = 0;
			junction_tree_f_explorer(bed_table, pos_table, connection_table, inb + f_search_head , rl - f_search_head, rl , explain_result[(total_pieces-1)*4], 0, explain_buff, 0, &result_total_matched_bases , explain_result+ 4 * (total_pieces>0?(total_pieces-1):0),  &result_number_pieces, my_value_array_index, ginp, result_head_modify, qualityb+ f_search_head,  EXON_FASTQ_FORMAT, match_score);

			total_total_matched_bases = result_total_matched_bases;

			total_pieces --;
		}

		if (result_number_pieces >0)
			total_pieces += result_number_pieces ;
		else{	
			total_pieces=0;
			is_settle=0;
		}


	//	printf(" POS AAA RL=%d, MATCHED_TOTAL=%d TEST_LEB=%d\n", rl, total_total_matched_bases, rl-f_search_head);

#ifdef DEBUG



		printf ("\nFSEARCH FROM %d MATCH BASES=%d/%d; PIECES=%d\n" , f_search_head, result_total_matched_bases, rl - f_search_head, result_number_pieces );
		for (xx=0; xx<total_pieces;xx++)
		{
			char * chro_name;
			unsigned int chro_pos;

			unsigned int piece_start_pos = explain_result[xx*4];
			int piece_start = explain_result[xx*4+1];
			int piece_end = explain_result[xx*4+2];

			locate_gene_position(piece_start_pos, &_global_offsets, &chro_name, &chro_pos);
			chro_pos -= piece_start ;

			printf ("FS-RES(ALL): [%d ~ %d  %c %s,%d] \n", piece_start, piece_end, ((ginp2 && (i % 2)) + negative_strand == 1)?'~':'@', chro_name , chro_pos);
		}

#endif


		if (total_pieces>4)
			is_settle = 0;

		char cigar_buf[100];

		if(halves_record -> best_vote1_list[i]< /*(rl - 15)/EXON_SUBREAD_GAP*/TOTAL_SUBREADS * EXON_MAJOR_HALF_VOTES_RATE ){ total_pieces=0 ; is_settle = 0;}

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
							halves_record -> best_vote2_list[i]=0;
						}
					}
					else{
						if(explain_result[ 4*j+2]- explain_result[ 4*j+1] >0)
							sprintf (cigar_piece, "%dM", explain_result[ 4*j+2]- explain_result[ 4*j+1] );
						else
						{
							halves_record -> best_vote1_list[i]=0;
							halves_record -> best_vote2_list[i]=0;
						}
					}

					strcat(cigar_buf, cigar_piece);
					if (j<total_pieces-1)
					{
#ifdef TEST_TARGET
						if(memcmp(inb, TEST_TARGET, 15)==0)
							printf("CIGAR- %dN p1=%u p2=%u pos_def=%d, half1=%d, indel=%d\n", explain_result[4*j+4]- (explain_result[ 4*j] + explain_result[ 4*j+2] - explain_result[4*j+1] + indel_length), explain_result[4*j+4]  , explain_result[4*j]  , explain_result[4*j+4] - explain_result[4*j],  explain_result[ 4*j+2] - explain_result[4*j+1] + indel_length,  indel_length);
#endif
						sprintf (cigar_piece, "%dN",explain_result[4*j+4]- (explain_result[ 4*j] + explain_result[ 4*j+2] - explain_result[4*j+1] + indel_length));
						strcat(cigar_buf, cigar_piece);
						last_indel = indel_length;

					}

				//	printf("CIG=%s\n", cigar_buf);

				}
			}/*
			    else
			    {
			    unsigned int pos_small = min(halves_record -> best_pos1_list[i], halves_record -> best_pos2_list[i]);
			    int h1_len = (halves_record -> is_reversed_list[i]) ? (rl - halves_record ->splice_point_list[i]) : halves_record ->splice_point_list[i];
			    long long int dist =  pos_small;
			    dist -= pos_small;
			    sprintf (cigar_buf, "%dM%dN%dM",  h1_len , (int)abs(dist), rl - h1_len);
			    }*/
			strncpy(halves_record -> cigar_string_buffer + i * EXON_MAX_CIGAR_LEN, cigar_buf, EXON_MAX_CIGAR_LEN-1);


			if (tolerable_scan)
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

					int *search_res = (int*) HashTableGet(bed_table, &search_key);
					//assert(search_res);
					if(search_res)
						(*search_res)++;
#ifdef DEBUG
					printf("FINAL RECOVERED-IS-FINISH\n");
#endif

				}
				if(is_settle)
					halves_record -> best_pos1_list[i] = explain_result[0];
			}

			compress_cigar(cigar_buf, rl, inb);
			halves_record -> final_quality [i] = final_mapping_quality(my_value_array_index, explain_result[0], inb, qualityb[0]?qualityb:NULL, cigar_buf, EXON_FASTQ_FORMAT);
			if(halves_record -> final_quality [i] < 10)	
				printf("QUAL=%.4f POS=%u READ=%s", halves_record -> final_quality [i] ,  explain_result[0] , inb);
		}
		else if(is_settle)
		{
			halves_record -> best_vote1_list[i] = max(9, halves_record -> best_vote1_list[i]);
			halves_record -> best_vote2_list[i] = 0;
			halves_record -> best_pos1_list[i] = explain_result[0];

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
			halves_record -> best_vote2_list[i]=0;
		}
		halves_record -> half_marks_list[i] |= IS_PROCESSED_READ_R2;

		if (is_settle)
			halves_record -> half_marks_list[i] |= IS_FINALISED_PROCESSING;
		i++;

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
	if (i == MAX_EXON_CONNECTIONS)return;	// too many connections

	connect_to ->connect_to[i] = p2;
	connect_to ->is_opposite_reversed[i] = is_p2_reversed;

	if (i<MAX_EXON_CONNECTIONS-1) connect_to -> connect_to[1+i]=0;

}


#define SHORT_EXON_MIN_LENGTH 16
#define SHORT_EXON_WINDOW 6 
#define SHORT_EXON_EXTEND 5000



void search_short_exons(char * inb0, char * qualityb0, int rl, HashTable * connection_table, HashTable * bed_table, HashTable * pos_table, unsigned int P1_Pos, unsigned int P2_Pos, short read_coverage_start, short read_coverage_end,  gene_value_index_t *base_index, int space_type, int is_negative, int tolerable_scan)
{
	char inb[1201], qualityb[1201];
	if (!EXON_EXTENDING_SCAN) return;

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
		//	printf("FOUND IN STEP 4: %u %u\n", pos_R1, pos_R2);
		int *search_res = (int*) HashTableGet(bed_table, &search_key);
		if(!search_res)
		{
			search_res = (int *)malloc(sizeof(int));
			*search_res = tolerable_scan? 0x80000000 : 0;
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
		//	printf("FOUND IN STEP 4: %u %u\n", pos_R1, pos_R2);
		int *search_res = (int*) HashTableGet(bed_table, &search_key);
		if(!search_res)
		{
			search_res = (int *)malloc(sizeof(int));
			*search_res = tolerable_scan? 0x80000000 : 0;
			paired_exon_key * new_key = (paired_exon_key *) malloc(sizeof(paired_exon_key));
			new_key->small_key = pos_small_x;
			new_key->big_key = pos_big_x;
			HashTablePut(bed_table, new_key, search_res);
		}
	}
}

void search_short_exons_old(char * inb0, char * qualityb0, int rl, HashTable * connection_table, HashTable * bed_table, HashTable * pos_table, unsigned int P1_Pos, unsigned int P2_Pos, short read_coverage_start, short read_coverage_end,  gene_value_index_t *base_index, int space_type, int is_negative, int tolerable_scan)
{
	char inb[1201], qualityb[1201];
	if (!EXON_EXTENDING_SCAN) return;

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

	float max_score , test_score;
	unsigned int best_j1_edge=0 , best_j2_edge=0;
	if (read_coverage_start > SHORT_EXON_MIN_LENGTH)// && bad_quality_base_number(qualityb, read_coverage_start, EXON_FASTQ_FORMAT)< 2)
	{
		max_score = -1;
		int window_end = read_coverage_start + SHORT_EXON_WINDOW + 3;
		if(8 > match_chro(inb, base_index, pos_small, SHORT_EXON_MIN_LENGTH , 0, space_type))
		{
			for(; window_end>=SHORT_EXON_WINDOW+SHORT_EXON_MIN_LENGTH; window_end --)
			{
				int matched_in_exon , matched_out_exon;
				matched_in_exon = match_chro(inb + window_end - SHORT_EXON_WINDOW, base_index, pos_small + window_end - SHORT_EXON_WINDOW , SHORT_EXON_WINDOW , 0, space_type);
				matched_out_exon = match_chro(inb + window_end - SHORT_EXON_WINDOW*2, base_index, pos_small + window_end - 2*SHORT_EXON_WINDOW , SHORT_EXON_WINDOW , 0, space_type);
				//		printf("\n MAT_IN=%d  MAT_OUT=%d\n",matched_in_exon,matched_out_exon );
				test_score = 100000 + matched_in_exon - matched_out_exon;
				if(matched_in_exon >= SHORT_EXON_WINDOW && matched_out_exon <= SHORT_EXON_WINDOW -1 )
				{
					// test donor
					char cc[3];
					cc[2]=0;
					get_chro_2base(cc, base_index, pos_small + window_end - SHORT_EXON_WINDOW-2, 0);
					if(is_donar_chars_part(cc))
					{
						//				printf("%s FOUND!\n", cc);
						unsigned int new_pos = pos_small + window_end - SHORT_EXON_WINDOW-3;
						while(1)
						{
							new_pos = match_chro_range(inb,  base_index, new_pos, window_end - SHORT_EXON_WINDOW , SHORT_EXON_EXTEND - ( pos_small + window_end - SHORT_EXON_WINDOW-3 - new_pos), SEARCH_BACK);
							//printf("%s BBT=%u\n", inb , new_pos);
							if(new_pos==0xffffffff) break;
							char cc2[3];
							cc2[2]=0;
								//						printf("%u FOUND!\n", i);
							get_chro_2base(cc2, base_index, new_pos + window_end - SHORT_EXON_WINDOW, 0);
							if(is_donar_chars_part(cc2) && paired_chars_part(cc2 , cc, 0)) 
							{
								if(test_score > max_score + 0.0001)
								{
			//						printf("%s ACCBBT=%u ; S=%.5f\n", inb , new_pos , test_score);
									best_j1_edge = new_pos + window_end - SHORT_EXON_WINDOW;
									best_j2_edge = pos_small + window_end - SHORT_EXON_WINDOW;
									max_score = test_score;
								}
							}
							//new_pos-=4;
						}
					}
				}
			}
		}
		if(max_score>0)
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
			//	printf("FOUND IN STEP 4: %u %u\n", pos_R1, pos_R2);
			int *search_res = (int*) HashTableGet(bed_table, &search_key);
			if(!search_res)
			{
				search_res = (int *)malloc(sizeof(int));
				*search_res = tolerable_scan? 0x80000000 : 0;
				paired_exon_key * new_key = (paired_exon_key *) malloc(sizeof(paired_exon_key));
				new_key->small_key = pos_small_x;
				new_key->big_key = pos_big_x;
				HashTablePut(bed_table, new_key, search_res);
			}
		}
	}

	if (read_coverage_end < rl - SHORT_EXON_MIN_LENGTH)// && bad_quality_base_number(qualityb + read_coverage_end, rl - read_coverage_start, EXON_FASTQ_FORMAT )< 2)
	{
		max_score = -1;
		int window_start = read_coverage_end - SHORT_EXON_WINDOW - 3;
		if(8 > match_chro(inb + rl - SHORT_EXON_MIN_LENGTH, base_index, pos_big + rl - SHORT_EXON_MIN_LENGTH, SHORT_EXON_MIN_LENGTH , 0, space_type))
		{

			for(; window_start<= rl -SHORT_EXON_WINDOW - SHORT_EXON_MIN_LENGTH; window_start++)
			{
				int matched_in_exon , matched_out_exon;
				matched_in_exon = match_chro(inb + window_start, base_index, pos_big + window_start , SHORT_EXON_WINDOW , 0, space_type);
				matched_out_exon = match_chro(inb + window_start + SHORT_EXON_WINDOW, base_index, pos_big + window_start +SHORT_EXON_WINDOW , SHORT_EXON_WINDOW , 0, space_type);
				test_score = 100000 + matched_in_exon - matched_out_exon;
				//		printf("\n MAT_IN=%d  MAT_OUT=%d\n",matched_in_exon,matched_out_exon );
				if(matched_in_exon >= SHORT_EXON_WINDOW  && matched_out_exon <= SHORT_EXON_WINDOW - 1)
				{
					// test donor
					char cc[3];
					cc[2]=0;
					get_chro_2base(cc, base_index, pos_big + window_start + SHORT_EXON_WINDOW, 0);
					if(is_donar_chars_part(cc))
					{

						//				printf("%s FOUND!\n", cc);

						unsigned int new_pos = pos_big + window_start + SHORT_EXON_WINDOW ;
						while(1)
						{

							//printf("%s FFT=%u\n", inb , new_pos);
							new_pos = match_chro_range(inb +  window_start + SHORT_EXON_WINDOW,  base_index, new_pos, rl - window_start - SHORT_EXON_WINDOW , SHORT_EXON_EXTEND - (new_pos - pos_big + window_start), SEARCH_FRONT);
							if(new_pos==0xffffffff) break;

							//printf("%s FFT=%u\n", inb , new_pos);
							char cc2[3];
							cc2[2]=0;
							//						printf("%u FOUND!\n", pos_small- new_pos);
							get_chro_2base(cc2, base_index, new_pos -2, 0);
							if(is_donar_chars_part(cc2) && paired_chars_part(cc , cc2, 0))
							{

								if(test_score > max_score+0.0001)
								{
		//							printf("%s ACCFFT=%u ; S=%.5f\n", inb , new_pos , test_score);
									best_j1_edge = pos_big + window_start + SHORT_EXON_WINDOW;
									best_j2_edge = new_pos;
									max_score = test_score;
								}
							}
							//new_pos+=4;
						}
					}
				}
			}
		}

		if(max_score>0)
		{


			pos_small = get_grouped_position(pos_table, best_j1_edge);
			pos_big = get_grouped_position(pos_table, best_j2_edge);

			paired_exon_key search_key; 
			search_key.small_key = pos_small;
			search_key.big_key = pos_big;



			put_connection_table(connection_table, pos_small, pos_big, is_negative, is_negative);
			put_connection_table(connection_table, pos_big, pos_small, is_negative, is_negative);

			//if (tolerable_scan && !(halves_record ->half_marks_list[i] & IS_FINALISED_PROCESSING))
			//	printf("FOUND IN STEP 4: %u %u\n", pos_R1, pos_R2);

			int *search_res = (int*) HashTableGet(bed_table, &search_key);
			if(!search_res)
			{
				search_res = (int *)malloc(sizeof(int));
				*search_res = tolerable_scan? 0x80000000 : 0;
				paired_exon_key * new_key = (paired_exon_key *) malloc(sizeof(paired_exon_key));
				new_key->small_key = pos_small;
				new_key->big_key = pos_big;
				HashTablePut(bed_table, new_key, search_res);
			}
		}
	}

}


void feed_exonbed(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table,  halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads, gene_value_index_t * my_value_array_index, int first_index, int last_index, int tolerable_scan)
{
	int i=0, ic=0;
	while (1)
	{
		char nameb[1201], inb[1201], qualityb[1201];
		int rl;
		if(i >= processed_reads*(1+(ginp2!=NULL)))break;

		if(i % 10000 ==0 && i>1)
			print_text_scrolling_bar("First Iteration", i*1./(1+(ginp2!=NULL))/processed_reads, 80, &ic);

		unsigned int pos = ( IS_R1_CLOSE_TO_5 & halves_record -> half_marks_list[i] ) ?halves_record -> best_pos1_list[i]:halves_record -> best_pos2_list[i];
		unsigned int pos2 = ( IS_R1_CLOSE_TO_5 & halves_record -> half_marks_list[i] ) ?halves_record -> best_pos2_list[i]:halves_record -> best_pos1_list[i];

		if (ginp2 && (i % 2))
			rl = geinput_next_read(ginp2, nameb, inb, qualityb);
		else
			rl = geinput_next_read(ginp, nameb, inb, qualityb);

		if (ginp2 && (i % 2))
			reverse_read(inb, rl, ginp2->space_type);

		if (rl<0){
			break;
		}

#ifdef TEST_TARGET
		if(memcmp(TEST_TARGET, inb,15)==0)
		{
			printf("RAW i=%d PRE TEST DONOR = %u, %u\n", i, halves_record -> best_pos1_list[i] ,halves_record -> best_pos2_list[i]);
		}
#endif

		/*
		   if (tolerable_scan && (halves_record -> half_marks_list[i] & IS_FINALISED_PROCESSING))
		   {
		   i++;
		   continue;
		   }
		 */

		if (!first_index)
		{
			if ((halves_record -> best_vote2_list[i]>=1) && min(pos, pos2) < my_value_array_index -> start_base_offset + 1000)
			{
				i++; continue;
			}
			else if((halves_record -> best_vote2_list[i]==0) && halves_record -> best_pos1_list[i] < my_value_array_index -> start_base_offset + 1000)
			{
				i++; continue;
			}
		}

		if (!last_index)
		{
			if ((halves_record -> best_vote2_list[i]>=1) && max(pos, pos2)  > my_value_array_index -> start_base_offset + my_value_array_index ->length - 1000)
			{
				i++; continue;
			}
			else if((halves_record -> best_vote2_list[i]==0) && halves_record -> best_pos1_list[i] > my_value_array_index -> start_base_offset + my_value_array_index ->length - 1000)
			{
				i++; continue;
			}
		}

		if (IS_PROCESSED_READ & halves_record -> half_marks_list[i])
		{
			i++;
			continue;
		}

		//printf("\n P1=%u  P2=%u \n",  halves_record -> best_pos1_list[i],  halves_record -> best_pos2_list[i]);

		//printf("FEED-P2 v1=%d, v2=%d\n", halves_record -> best_vote1_list[i] ,halves_record -> best_vote2_list[i] );
		long long int dist = halves_record -> best_pos1_list[i] ;
		dist-= halves_record -> best_pos2_list[i];

		if(halves_record -> best_vote1_list[i]>= (int)( /*(rl -15)/EXON_SUBREAD_GAP*/ TOTAL_SUBREADS* EXON_MAJOR_HALF_VOTES_RATE) && halves_record -> best_vote2_list[i]>=MIN_VOTE2_TMP_VAR )
		{
			//printf("FEED-P3\n");
			//char * chro_name;
			//unsigned int chro_pos, chro_pos_small, chro_pos_large;
			//char is_strong = (min(halves_record ->best_vote1_list[i],halves_record ->best_vote2_list[i])>5);

			search_short_exons(inb, qualityb, rl, connection_table, bed_table, pos_table, halves_record -> best_pos1_list[i], halves_record -> best_pos2_list[i], halves_record -> read_coverage_start[i], halves_record -> read_coverage_end[i], my_value_array_index, ginp->space_type, halves_record -> is_reversed_list[i], tolerable_scan);
			//char short_overlap = !(halves_record -> half_marks_list[i] & IS_LONG_OVERLAP);
			char negative_strand = halves_record -> is_reversed_list[i];
			int real_break_point;
			char is_soft_condition = 0;
			int test_range = rl / 4;
			char is_accepted = 0;
			int indel_offset_pos1=0, indel_offset_pos2=0;
			int guess_break_point = (halves_record -> is_reversed_list[i]) ? (rl - halves_record ->splice_point_list[i]) : halves_record ->splice_point_list[i];
			int best_donor_score =-1;

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

						int indel_offset2 = ((indel_x2 %2)?1:-1) * ((indel_x2+1)>>1);
						int test_donor_score=-1;
						int test_real_break_point;

						is_accepted = test_donor(inb, rl, min(pos, pos2), max(pos,pos2), guess_break_point, negative_strand, test_range, is_soft_condition, EXON_INDEL_TOLERANCE, &test_real_break_point, my_value_array_index, indel_offset1, indel_offset2, negative_strand, ginp->space_type, halves_record ->best_vote2_list[i], &test_donor_score);

						test_donor_score -= abs(indel_offset2)*1 + abs(indel_offset1)*1;

						if (is_accepted  && (test_donor_score > best_donor_score)){
							//if(best_donor_score >0)
							//	printf("TEST SCORE=%d MAX=%d O1=%d O2=%d\n", test_donor_score , best_donor_score , indel_offset1, indel_offset2);
							best_indel_p1 = indel_offset1;
							best_indel_p2 = indel_offset2;
							best_donor_score = test_donor_score;
							real_break_point = test_real_break_point;
						}
					}
				}

				if(best_donor_score >0)
				{
					//printf("FINAL BEST SCORE=%d\n", best_donor_score);
					indel_offset_pos1 = best_indel_p1;
					indel_offset_pos2 = best_indel_p2;
					halves_record ->splice_point_offset_1[i] = best_indel_p1;
					halves_record ->splice_point_offset_2[i] = best_indel_p2;
					is_accepted =1;
				}
			}
			else
			{
				is_accepted = test_donor(inb, rl, min(pos, pos2), max(pos,pos2), guess_break_point, negative_strand, test_range, is_soft_condition, EXON_INDEL_TOLERANCE, &real_break_point, my_value_array_index,0,0, negative_strand, ginp->space_type, halves_record ->best_vote2_list[i], &best_donor_score);
				halves_record ->splice_point_offset_2[i]=0;
				halves_record ->splice_point_offset_1[i]=0;
			}

			//if ((i % 2) && is_accepted)printf ("22 SECOND!\n");

			//printf ("\nPRE-ACCEPT=%d\n", is_accepted);

			if((! is_accepted) && EXON_FUSION_DETECTION && (halves_record -> half_marks_list[i] & IS_FUSION))
			{
				int p1_offset = (halves_record ->half_marks_list[i] & IS_R1_CLOSE_TO_5)? 20:0; 
				int p2_offset = (halves_record ->half_marks_list[i] & IS_R1_CLOSE_TO_5)? 0:20;
				int c1_offset = (!(halves_record ->half_marks_list[i] & IS_R1_CLOSE_TO_5) ^ !(halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1))? 20:0; 
				int c2_offset = (!(halves_record ->half_marks_list[i] & IS_R1_CLOSE_TO_5) ^ !(halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R2))? 0:20; 
				int m1 = match_chro(inb+p1_offset, my_value_array_index, halves_record -> best_pos1_list[i]+c1_offset , rl -20, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1)?1:0, ginp->space_type);
				int m2 = match_chro(inb+p2_offset, my_value_array_index, halves_record -> best_pos2_list[i]+c2_offset , rl -20, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R2)?1:0, ginp->space_type);
				//printf ("M1=%d   M2=%d    H=%f\n", m1, m2, rl * .85);
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

				//printf("ACCP: %u, %u / indel = %d, %d?\n", pos_R1, pos_R2, indel_offset_pos1, indel_offset_pos2);

				halves_record -> splice_point_list[i] = real_break_point;

				pos_R1 = get_grouped_position(pos_table, pos_R1);
				pos_R2 = get_grouped_position(pos_table, pos_R2);

				unsigned int pos_small = min(pos_R1, pos_R2);
				unsigned int pos_big   = max(pos_R1, pos_R2);

#ifdef TEST_TARGET
				if(memcmp(TEST_TARGET, inb, 15)==0) 
					printf("ACCEPTED: %u, %u\nRAW: %u, %u\nBP:%d\tOFFSETS=%d %d\n", pos_small, pos_big, halves_record -> best_pos1_list[i] , halves_record -> best_pos2_list[i], real_break_point, indel_offset_pos1, indel_offset_pos2);
#endif

				paired_exon_key search_key; 
				search_key.small_key = pos_small;
				search_key.big_key = pos_big;
				put_connection_table(connection_table, pos_R1, pos_R2, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1)?1:0, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R2)?1:0);
				put_connection_table(connection_table, pos_R2, pos_R1, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R2)?1:0, (halves_record ->half_marks_list[i] & IS_NEGATIVE_STRAND_R1)?1:0);

				//if (tolerable_scan && !(halves_record ->half_marks_list[i] & IS_FINALISED_PROCESSING))
				//	printf("FOUND IN STEP 4: %u %u\n", pos_R1, pos_R2);

				int *search_res = (int*) HashTableGet(bed_table, &search_key);
				if(!search_res)
				{
					search_res = (int *)malloc(sizeof(int));
					*search_res = tolerable_scan? 0x80000000 : 0;
					paired_exon_key * new_key = (paired_exon_key *) malloc(sizeof(paired_exon_key));
					new_key->small_key = pos_small;
					new_key->big_key = pos_big;
					//printf("INSERT0: POS=%u , %u\n",  pos_small, pos_big);
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
		halves_record -> half_marks_list[i] |= IS_PROCESSED_READ;
		i++;

	}
}

void print_exon_res(gene_value_index_t *array_index , halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads)
{
	int i=0, ic=0 ;

	if(ftello(out_fp)<1)
	{
		unsigned int last_offset = 0;
		i=0;
		while(_global_offsets.read_offset[i])
		{
			fprintf(out_fp, "@SQ\tSN:%s\tLN:%u\n", _global_offsets.read_name[i], _global_offsets.read_offset[i] - last_offset+16);
			last_offset = _global_offsets.read_offset[i];
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


		if(i % 10000 ==0 && i>1)
			print_text_scrolling_bar("Saving results", i*1./(1+(ginp2!=NULL))/processed_reads, 80, &ic);

		old_line[0]=0;

		if (ginp2 && (i % 2))
		{
			if(ginp->file_type >= GENE_INPUT_SAM_SINGLE)
				geinput_readline_back(ginp2, old_line);
			rl = geinput_next_read(ginp2, nameb, inb, qualityb);

			if(!halves_record -> is_reversed_list[i])
			{
				int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb[0]>='A' && inb[0]<='Z')?1:0;
				reverse_read(inb+ rev_offset, rl- rev_offset, ginp->space_type);
				if(qualityb[0])
					reverse_quality(qualityb, rl- rev_offset);
			}

		}
		else
		{
			if(ginp->file_type >= GENE_INPUT_SAM_SINGLE)
				geinput_readline_back(ginp, old_line);
			rl = geinput_next_read(ginp, nameb, inb, qualityb);

			if(halves_record -> is_reversed_list[i])
			{
				int rev_offset = (ginp->space_type==GENE_SPACE_COLOR && inb[0]>='A' && inb[0]<='Z')?1:0;
				reverse_read(inb+ rev_offset, rl- rev_offset, ginp->space_type);
				if(qualityb[0])
					reverse_quality(qualityb, rl- rev_offset);
			}

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






		long long int dist = halves_record -> best_pos1_list[i] ;
		dist-= halves_record -> best_pos2_list[i];

		if((halves_record -> best_vote1_list[i]>=1 && halves_record -> best_vote2_list[i]>=1 ))
		{
			char * chro_name, * chro_name2, *chro_name_min;
			unsigned int chro_pos, chro_pos_small, chro_pos_large;
			char cigar_buf[100], *cigar_print;

			unsigned int pos = ( IS_R1_CLOSE_TO_5 & halves_record -> half_marks_list[i] ) ?halves_record -> best_pos1_list[i]:halves_record -> best_pos2_list[i];
			unsigned int pos2 = ( IS_R1_CLOSE_TO_5 & halves_record -> half_marks_list[i] ) ?halves_record -> best_pos2_list[i]:halves_record -> best_pos1_list[i];

			locate_gene_position(min(pos,pos2), &_global_offsets, &chro_name_min, &chro_pos);

			pos = halves_record -> is_reversed_list[i]? (pos + rl):(pos);
			pos2 = halves_record -> is_reversed_list[i]? (pos2):(pos2+rl);

			int h1_len = (halves_record -> is_reversed_list[i]) ? (rl - halves_record ->splice_point_list[i]) : halves_record ->splice_point_list[i];

			gene_quality_score_t mapping_quality;

			if (halves_record -> cigar_string_buffer[i * EXON_MAX_CIGAR_LEN])
			{
				cigar_print = halves_record -> cigar_string_buffer + i * EXON_MAX_CIGAR_LEN;
			}else
			{
				sprintf (cigar_buf, "%dM%dN%dM",  h1_len , (int)abs(dist), rl - h1_len);
				cigar_print = cigar_buf;
			}
			mapping_quality = halves_record -> final_quality[i] ;


			fprintf(out_fp,"%s\t%d\t%s\t%u\t%d\t%s\t*\t0\t0\t%s\t%s",  nameb, halves_record -> is_reversed_list[i]?16:0 , chro_name_min, chro_pos+1, (int)(0.5+mapping_quality) , cigar_print , inb, qualityb);

			if(0 && halves_record -> half_marks_list[i] & IS_FUSION)
			{
				unsigned int pos_R1 = halves_record -> best_pos1_list[i];
				unsigned int pos_R2 = halves_record -> best_pos2_list[i];
				locate_gene_position(pos_R1, &_global_offsets, &chro_name, &chro_pos_small) ;
				locate_gene_position(pos_R2, &_global_offsets, &chro_name2, &chro_pos_large) ;

				fprintf(out_fp,"\t%c%s,%d:%c%s,%d", (halves_record -> half_marks_list[i] & IS_NEGATIVE_STRAND_R1)?'~':'@' ,chro_name, chro_pos_small, (halves_record -> half_marks_list[i] & IS_NEGATIVE_STRAND_R2)?'~':'@' ,chro_name2,chro_pos_large);
			}

			fprintf(out_fp,"\n");

		}
		else if(!old_line[0])
		{
			if (halves_record -> best_vote1_list[i]>=8)
			{
				char * chro_name;
				unsigned int chro_pos;
				unsigned int pos = halves_record -> best_pos1_list[i];

				locate_gene_position(pos, &_global_offsets, &chro_name, &chro_pos);
				char * pair_info_refed = (IS_PAIRED_HINTED & halves_record -> half_marks_list[i]) ? "PA7":"SG7";
				if (IS_RECOVERED_JUNCTION_READ & halves_record -> half_marks_list[i] ) pair_info_refed = "RC7";
				if (IS_RECOVERED_JUNCTION_READ_STEP4 & halves_record -> half_marks_list[i] ) pair_info_refed = "RF7";

				fprintf(out_fp,"%s\t%d\t%s\t%u\t%d\t%dM\t*\t0\t0\t%s\t%s\n",  nameb, halves_record -> is_reversed_list[i]?16:0 , chro_name, chro_pos, halves_record -> best_vote1_list[i] , rl, inb, qualityb);
			}
			else
			{
				fprintf(out_fp,"%s\t4\t*\t*\t0\t*\t*\t0\t0\t%s\t%s\n",  nameb,inb, qualityb);
			}
		}
		else
			fprintf(out_fp, "%s\n",old_line);
		i++;

	}

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

void fragile_junction_voting(HashTable * bed_table, HashTable * pos_table, HashTable * connection_table, gehash_t * my_table, gene_value_index_t * my_value_array_index , int table_no,  halves_record_t * halves_record, char * read, char * qual, unsigned int full_rl, int negative_strand, int color_space, unsigned int low_border, unsigned int high_border)
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
		gene_vote_t vote_p1;
		
		init_gene_vote(&vote_p1);
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

				//printf("TQ: POS=%d, TOL=%d, INT=%u, SR=%s\n", subread_offset  , tolerable_scan , subread_integer, subread_string);
				gehash_go_q(my_table, subread_integer , subread_offset, read_len,negative_strand, &vote_p1, 1, 1, 21.9, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no,  low_border, high_border - read_len);
			}
			if(subread_offset1 >= read_len -16)
				break;
		}


		if(1)
		{
			finalise_vote(&vote_p1);
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

			int splice_point = select_best_matching_halves(&vote_p1, &best_pos1, &best_pos2, &best_vote1, &best_vote2, &is_abnormal ,&half_marks, &is_reversed_halves, accepted_support_rate, read_len, -1,  0, &read_coverage_start, &read_coverage_end, &indel_in_p1, &indel_in_p2, &max_cover_start, &max_cover_end, read_len);
			if (splice_point>0 && best_vote1 >= 1 && best_vote2>=1)
			{
				int test_real_break_point = -1, test_donor_score=-1;
				int is_accepted = test_donor(InBuff, read_len, min(best_pos1, best_pos2), max(best_pos1,best_pos2), splice_point, negative_strand, read_len/4, 0, EXON_INDEL_TOLERANCE, &test_real_break_point, my_value_array_index, 0, 0, negative_strand, color_space, halves_record ->best_vote2_list[i], &test_donor_score);

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

					int *search_res = (int*) HashTableGet(bed_table, &search_key);
					if(!search_res)
					{
						//printf("FRAG: the %d-th window: %d : %d \t\tFOUND: %u-%u\n", ww, window_cursor , read_len,pos_R1, pos_R2 );

						search_res = (int *)malloc(sizeof(int));
						*search_res = 0;
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
	int queries = 0;
	double t0=miltime();
	int is_reversed;
	int good_match = 0;
	float subread_step; // = EXON_SUBREAD_GAP+.0001;
	struct stat read_fstat;
	stat (ginp->filename, &read_fstat);
	long long int read_fsize = read_fstat.st_size;

	double local_begin_ftime = miltime();
	int read_len = 0, read2_len = 0;
	is_reversed = 0;


	if (ginp2)
		all_reads /=2;

	//	printf ("I'm the %d-th thread\n", my_thread_no);
	gene_vote_t vote_p1, vote_p2;

	unsigned int low_border = my_value_array_index -> start_base_offset;
	unsigned int high_border = my_value_array_index -> start_base_offset + my_value_array_index -> length; 


	while (1)
	{

		char namebuf[200];


		if (my_thread_no==0 && queries % 10000 == 0 && !is_reversed)
		{
			if(table_no == 0)
			{
				long long int current_reads = base_number + queries;
				long long int fpos = ftello(ginp->input_fp);
				reads_density = fpos*1.0/current_reads; 
				//printf ("\nDENS=%.5f, POS=%llu, fsize=%llu\n", reads_density, fpos, read_fsize);
			}
			if(IS_DEBUG && queries % 100000==0)
				printf("@LOG Done %d/%d, good %d, last time %f\n",queries, table_no, good_match, miltime() - t0);
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
			fflush(stdout);
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
				init_gene_vote(&vote_p1);
				if(ginp2)
					init_gene_vote(&vote_p2);
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

				read_len = geinput_next_read(ginp, namebuf, InBuff, QualityBuff);

				if(read_len>0){
					queries = *processed_reads;
					(*processed_reads ) ++;
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

				}
				if(ginp2)
				{
					InBuff2 = BuffMemory + 1250;
					QualityBuff2 = BuffMemory2 + 1250;

					read2_len = geinput_next_read(ginp2, namebuf, InBuff2, QualityBuff2);

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

			init_gene_vote(&vote_p1);
			if(ginp2)
				init_gene_vote(&vote_p2);
		}




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
				gene_vote_t * current_vote = is_second_read?&vote_p2: &vote_p1;
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

						if (tolerable_scan && (!EXON_NO_TOLERABLE_SCAN))
							gehash_go_q_tolerable(my_table, subread_integer , subread_offset, current_rlen, is_reversed, current_vote, 1, 1, 21.9, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no, tolerable_scan, low_border, high_border - current_rlen);
						else
							gehash_go_q(my_table, subread_integer , subread_offset, current_rlen, is_reversed, current_vote, 1, 1, 21.9, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no,  low_border, high_border - current_rlen);
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
				finalise_vote(&vote_p1);
				finalise_vote(&vote_p2);
				int is_paired_match = select_positions(&vote_p1, &vote_p2, &numvote_read1, &numvote_read2, &sum_quality, &qual_r1, &qual_r2, &pos_read1, &pos_read2, record_index1, record_index2, EXON_MAX_PAIRED_DISTANCE, EXON_MIN_PAIRED_DISTANCE, 3, ACCEPT_MINOR_SUBREADS, is_reversed, EXON_NUMBER_OF_ANCHORS_PAIRED,  EXON_INDEL_TOLERANCE, &is_breakeven);
				for (is_second_read = 0; is_second_read <2; is_second_read ++)
				{
					gene_vote_t * current_vote = is_second_read?&vote_p2: &vote_p1;
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

					int splice_point = select_best_matching_halves(current_vote, &best_pos1, &best_pos2, &best_vote1, &best_vote2, &is_abnormal ,&half_marks, &is_reversed_halves, accepted_support_rate, current_rlen, hint_pos, tolerable_scan, &read_coverage_start, &read_coverage_end, &indel_in_p1, &indel_in_p2, &max_cover_start, &max_cover_end, current_rlen);
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
						fragile_junction_voting(bed_table, pos_table, connection_table, my_table, my_value_array_index , table_no, halves_record, is_second_read?InBuff2:InBuff,  is_second_read?QualityBuff2:QualityBuff, current_rlen , is_reversed, ginp->space_type,  low_border,  high_border);
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

					//printf("TQ: POS=%d, TOL=%d, INT=%u, SR=%s\n", subread_offset  , tolerable_scan , subread_integer, subread_string);

					if (tolerable_scan && (!EXON_NO_TOLERABLE_SCAN))
						gehash_go_q_tolerable(my_table, subread_integer , subread_offset, read_len, is_reversed, &vote_p1, 1, 1, 21.9, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no, tolerable_scan,  low_border, high_border - read_len);
					else
						gehash_go_q(my_table, subread_integer , subread_offset, read_len, is_reversed, &vote_p1, 1, 1, 21.9, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no,  low_border, high_border - read_len);
				}
				if(subread_offset1 >= read_len -16)
					break;
			}

			if((!EXON_FUSION_DETECTION) || is_reversed)
			{
				finalise_vote(&vote_p1);
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

				int splice_point = select_best_matching_halves(&vote_p1, &best_pos1, &best_pos2, &best_vote1, &best_vote2, &is_abnormal ,&half_marks, &is_reversed_halves, accepted_support_rate, read_len, -1,  tolerable_scan, &read_coverage_start, &read_coverage_end, &indel_in_p1, &indel_in_p2, &max_cover_start, &max_cover_end, read_len);
				if (splice_point>0 && (vote_p1.max_vote >= halves_record -> best_vote1_list[queries]))
				{
					#ifdef DEBUG
					printf("SPP=%d, v1=%d, v2=%d, POS=%u, %%=%s\nRead_Start=%d, Read_End=%d\n", splice_point, best_vote1, best_vote2, best_pos1, InBuff, read_coverage_start, read_coverage_end);
					#endif

					add_best_matching_halves(halves_record, best_pos1, best_pos2, best_vote1, best_vote2,is_abnormal, is_reversed_halves, splice_point, half_marks, queries, max_cover_start, max_cover_end, read_coverage_start, read_coverage_end, indel_in_p1, indel_in_p2);
				}
				else 
				{
					//printf("EDD=%d, v1=%d, v2=%d, R=%s, RL=%d, MM=%d, SPS=%.5f\n", splice_point, best_vote1, best_vote2, InBuff, read_len, vote_p1.max_vote, subread_step);
					is_reversed_halves = (vote_p1.max_mask & IS_NEGATIVE_STRAND)?1:0;
					if (vote_p1.max_vote > (halves_record -> best_vote1_list[queries] + 0* halves_record -> best_vote2_list[queries]))
					{
						halves_record -> best_vote1_list[queries] = vote_p1.max_vote ;
						halves_record -> best_vote2_list[queries] = 0 ;
						halves_record -> best_pos1_list[queries] = vote_p1.max_position;
						halves_record -> is_reversed_list[queries] = is_reversed_halves;
						halves_record -> half_marks_list[queries] = vote_p1.max_mask;
						halves_record -> best1_read_start_pos[queries] = vote_p1.max_coverage_start;
						halves_record -> best1_read_end_pos[queries] = vote_p1.max_coverage_end;

						halves_record -> read_coverage_start[queries] = vote_p1.max_coverage_start;
						halves_record -> read_coverage_end[queries] = vote_p1.max_coverage_end;
					}
				}
			}
			if(read_len >= EXON_LONG_READ_LENGTH)
				fragile_junction_voting(bed_table, pos_table, connection_table, my_table, my_value_array_index , table_no, halves_record, InBuff,  QualityBuff, read_len , is_reversed, ginp->space_type,  low_border,  high_border);

		}

	
		if (is_reversed)
		{
			if(queries >= all_reads)
				break;
		}

		is_reversed = !is_reversed;
	}

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
	unsigned long long int current_fp2 = 0, last_fp_pos = 0, last_fp2_pos = 0;
	struct stat filestat;
	int last_table = 0;
	int stat_ret;
	char table_fn [300];
	gehash_t my_raw_table;
	gehash_t * my_table = &my_raw_table;
	gene_value_index_t value_array_index;

	if(ginp2)
		current_fp2 = ftello(ginp2 -> input_fp);

	while (1)
	{
		sprintf(table_fn, "%s.%02d.%c.tab", index_prefix, tabno, ginp->space_type==GENE_SPACE_COLOR?'c':'b');

// Test and load the index partition
		stat_ret = stat(table_fn, &filestat);
		if (stat_ret !=0)
			break;

		if (IS_DEBUG)
			printf ("@LOG Loading table from %s\n", table_fn);
		else
			printf ("Loading the %02d-th index file ...					      \r", tabno+1);
		fflush(stdout);

		gehash_load(my_table, table_fn);
		if(EXON_USE_VALUE_ARRAY_INDEX)
		{
			
			sprintf(table_fn, "%s.%02d.%c.array", index_prefix, tabno, ginp->space_type==GENE_SPACE_COLOR?'c':'b');
			stat_ret = stat(table_fn, &filestat);
			if (stat_ret !=0)
			{
				printf("The index does not contain any raw base data which is required in detecting junctions. Please use the -b ooption while building the index.\n");
				return -1;
			}
			if (tabno>0)gvindex_destory(&value_array_index);

			gvindex_load(&value_array_index,table_fn);
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
		if(ginp2)
			last_fp2_pos = ftello(ginp2-> input_fp);



		fseeko(ginp -> input_fp, current_fp, SEEK_SET);
		if (ginp2)
			fseeko(ginp2 -> input_fp, current_fp2, SEEK_SET);


	}

	if (tabno == 1)
	{
		feed_exonbed(overall_exon_bed,pos_table, connection_table, halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads, &value_array_index, 1,1, tolerable_scan);

		fseeko(ginp -> input_fp, current_fp, SEEK_SET);

		if (ginp2)
			fseeko(ginp2 -> input_fp, current_fp2, SEEK_SET);

		explorer_junc_exonbed(overall_exon_bed,pos_table, connection_table, halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads, &value_array_index, 1,1, tolerable_scan);
	}
	else
	{
		//gvindex_destory(&value_array_index);
		tabno = 0;
		while (1)
		{
			sprintf(table_fn, "%s.%02d.%c.array", index_prefix, tabno, ginp->space_type==GENE_SPACE_COLOR?'c':'b');

			stat_ret = stat(table_fn, &filestat);
			if (stat_ret !=0)
				break;

			gvindex_destory(&value_array_index);
			gvindex_load(&value_array_index,table_fn);
			feed_exonbed(overall_exon_bed, pos_table, connection_table, halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads, &value_array_index, tabno==0, tabno== last_table, tolerable_scan);

			fseeko(ginp -> input_fp, current_fp, SEEK_SET);
			if (ginp2)
				fseeko(ginp2 -> input_fp, current_fp2, SEEK_SET);

			explorer_junc_exonbed(overall_exon_bed,pos_table, connection_table, halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads, &value_array_index,tabno==0, tabno== last_table , tolerable_scan);
			tabno ++;

			fseeko(ginp -> input_fp, current_fp, SEEK_SET);
			if (ginp2)
				fseeko(ginp2 -> input_fp, current_fp2, SEEK_SET);
			
		}
	}

	fseeko(ginp -> input_fp, current_fp, SEEK_SET);
	if (ginp2)
		fseeko(ginp2 -> input_fp, current_fp2, SEEK_SET);



	if (out_fp  && tolerable_scan)
		print_exon_res(NULL , halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads);
	gvindex_destory(&value_array_index);

	fseeko(ginp -> input_fp, last_fp_pos, SEEK_SET);
	if (ginp2)
		fseeko(ginp2 -> input_fp, last_fp2_pos, SEEK_SET);



	return processed_reads;
	
}


int run_exon_search_index(gene_input_t * ginp, gene_input_t * ginp2, char * index_prefix, halves_record_t * halves_record, FILE * out_fp, unsigned long long int base_number, int all_tables, unsigned long long int *succeed_reads, HashTable * overall_exon_bed, HashTable * pos_table, HashTable * connection_table)
{

	long long int tmp_pos, tmp_pos2 = 0;
	int ret = 0;

	
	if(!EXON_NO_TOLERABLE_SCAN)
	{
		tmp_pos = ftello(ginp -> input_fp);
		if (ginp2)
			tmp_pos2 = ftello(ginp2 -> input_fp);

		ret =run_exon_search_index_tolerable(ginp, ginp2, index_prefix, halves_record, out_fp, base_number, all_tables, succeed_reads,  overall_exon_bed, pos_table,  connection_table,0);
		clear_processed_marks(halves_record);

		fseeko(ginp -> input_fp, tmp_pos, SEEK_SET);
		if (ginp2)
			fseeko(ginp2 -> input_fp, tmp_pos2, SEEK_SET);
	}
	ret = run_exon_search_index_tolerable(ginp, ginp2, index_prefix, halves_record, out_fp, base_number, all_tables, succeed_reads,  overall_exon_bed, pos_table,  connection_table,1);
	return ret;
}


void exon_usage(char * execname)
{
	puts("Usage:");
	puts(" ./subjunc [options] -i <index_name> -r <input> -o <output>");
	puts("");
	puts("Required arguments:");
	puts("    -i --index     <index>\t base name of the index.");
	puts("    -r --read      <input>\t name of the input file(FASTQ/FASTA format). Both base-space and color-space read data are supported. For paired-end reads, this gives the first read file and the other read file should be specified using the -R option.");
	puts("    -o --output    <output>\t name of the output file(SAM format).");
	puts("");
	puts("Optional arguments:");
        puts("    -n --subreads  <int>\t number of selected subreads, 14 by default.");
        puts("       --singleSAM <input>\t using as input file a SAM file which includes mapping results for single-end reads (e.g. 'subread-align' output).");
        puts("       --pairedSAM <input>\t using as input file a SAM file which includes mapping results for paired-end reads.");
	puts("    -x --fusion           \t enabling the detection of fusion events.");
	puts("    -A --nosam           \t disabling the SAM output for the reads. Only discovered exon junction locations will be reported (BED file).");
	//puts("    -H --hamming        \t passed to subread-align program.");
        //puts("    -Q --quality        \t passed to subread-align program.");
        //puts("    -u --unique        \t passed to subread-align program.");
	//puts("    -L --halflen   <int>\t minimal distance allowed between the junction location and the end base (first or last) of the read or the segment, 0 by default.");
	puts("    -T --threads   <int>\t number of threads/CPUs used, 1 by default.");
	puts("    -I --indel     <int>\t number of INDEL bases allowed, 5 by default.");
        puts("    -P --phred     <3:6>\t the format of Phred scores used in input files, '3' for phred+33 and '6' for phred+64. '3' by default.");
	puts("");
	puts("Arguments for paired-end reads:");
	puts("    -R --read2     <input>\t name of the second input file from paired-end data. The program will then be switched to paired-end read mapping mode.");
	puts("    -d --mindist   <int>\t the minimum distance between two reads in a pair, 50 by default");
	puts("    -D --maxdist   <int>\t the maximum distance between two reads in a pair, 600 by default");
	puts("    -S --order     <ff:fr:rf> \t specifying if the first/second reads are forward or reversed, 'fr' by default.");
	puts("");
	puts("Example:");
	puts(" ./subjunc -i my_index -r reads.fastq -o my_result.sam ");
	puts("");
	puts("");
	puts("For more information about these arguments, please refer to the User Manual.\n");

}

static struct option long_options[] =
{
	{"basewise", no_argument, &EXON_USE_VALUE_ARRAY_INDEX, 1},
	{"fusion", no_argument, &EXON_FUSION_DETECTION, 1},
	{"index", required_argument, 0, 'i'},
	{"read",  required_argument, 0, 'r'},
	{"read2", required_argument, 0, 'R'},
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
	{0, 0, 0, 0}
};


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
			int *counter = (int*) cursor ->value;

			locate_gene_position( p->small_key , &_global_offsets, &chro_name, &chro_pos);
			locate_gene_position( p->big_key , &_global_offsets, &chro_name2, &chro_pos_big);
			char * is_fusion =(chro_name2 == chro_name)?"SP0":"FS0";
			char * is_recovered = (*counter & 0x80000000)?"RC1":"NM1";
			if ((*counter& 0x7fffffff) >0) 
			{
				(*junction_number)++;
				(*support_number)+= (*counter& 0x7fffffff);
				fprintf(ofp,"%s\t%u\t%s\t%u\t%d\t%s,%s,\n", chro_name, chro_pos, chro_name2, chro_pos_big, (*counter& 0x7fffffff), is_fusion, is_recovered);
			}
			cursor = cursor->next;
		}
	}
	fclose(ofp);
}

int main_junction(int argc,char ** argv)
{
	char read_file [300], read2_file [300];
	char output_file [300];

	char tmpfile[100];
	strcpy(tmpfile, "./subjunc-sam-XXXXXX");

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

	TOTAL_SUBREADS = 14;
	EXON_MAJOR_HALF_VOTES_RATE = .2;
	ACCEPT_MINOR_SUBREADS = 1;
	INDEX_THRESHOLD = 24;
	read_file[0]=0;
	read2_file[0]=0;
	index_prefix[0]=0;
	output_file[0]=0;
	all_reads = 14*1024*1024;

	int using_base_distance = 0;

	printf("\n");



	while ((c = getopt_long (argc, argv, "xbS:L:AH:d:D:n:m:p:f:P:R:r:i:l:o:T:Q:I:1:2:t:?", long_options, &option_index)) != -1)
		switch(c)
		{
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
			//case 'b':
			//	EXON_USE_VALUE_ARRAY_INDEX = 1;
			//	break;
			case 'D':
				EXON_MAX_PAIRED_DISTANCE = atoi(optarg);
				break;
			case 'd':
				EXON_MIN_PAIRED_DISTANCE = atoi(optarg);
				break;
			case 'n':
				TOTAL_SUBREADS = atoi(optarg);
			//	printf(" === WARNING ===\n You cannot set subread numbers in exon detection. It is automatically set to the max number.");
				break;
			case 'f':
				INDEX_THRESHOLD  = atoi(optarg);
				break;
			case 'm':
				EXON_MAJOR_HALF_VOTES_RATE = atof(optarg);
				break;
				
			case 'T':
				EXON_ALL_THREADS = atoi(optarg);
				break;
			case 'r':
				if(IS_SAM_INPUT)
				{
					puts("You cannot specify the input files in FASTQ/FASTA formats when specifying a SAM file as input");
					return -1;
				}
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
			//	if(read_file[0] || IS_SAM_INPUT)
			//	{
			//		puts("You cannot specify the input files in FASTQ/FASTA formats when specifying a SAM file as input.");
			//		return -1;
			//	}
				strncpy(read_file, optarg, 299);
				IS_SAM_INPUT=1;
				EXON_FASTQ_FORMAT = FASTQ_PHRED33;
				break;
			case '2':
			//	if(read_file[0]|| IS_SAM_INPUT)
			//	{
			//		puts("You cannot specify the input files in FASTQ/FASTA formats when specifying a SAM file as input");
			//		return -1;
			//	}
				IS_SAM_INPUT=2;
				strncpy(read_file, optarg, 299);
				EXON_FASTQ_FORMAT = FASTQ_PHRED33;
				strncpy(read2_file, optarg, 299);
				break;
			case 't':
				sprintf(tmpfile, "%s/XXXXXX", optarg);
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


		if(read2_file[0])
			sprintf(command, "%ssubread-align -J -T %d -i '%s' -r '%s' -R '%s' -o '%s' -P %d -d %d -D %d %s %s -I %d ",cwd, EXON_ALL_THREADS, index_prefix, read_file, read2_file, tmpfile, EXON_FASTQ_FORMAT == FASTQ_PHRED33?3:6, EXON_MIN_PAIRED_DISTANCE, EXON_MAX_PAIRED_DISTANCE, EXON_QUALITY_SCALE?"-Q":"", using_base_distance?"-H":"", EXON_INDEL_TOLERANCE-1);
		else
			sprintf(command, "%ssubread-align -J -T %d -i '%s' -r '%s' -o '%s' -P %d %s %s -I %d",cwd, EXON_ALL_THREADS, index_prefix, read_file, tmpfile, EXON_FASTQ_FORMAT == FASTQ_PHRED33?3:6 , EXON_QUALITY_SCALE?"-Q":"", using_base_distance?"-H":"", EXON_INDEL_TOLERANCE-1);
		puts(command);
		is_successful = system(command);

		strcpy(read_file, tmpfile);
		IS_SAM_INPUT=1;

		if(read2_file[0])
		{
			strcpy(read2_file, tmpfile);
			IS_SAM_INPUT=2;
		}
		EXON_FASTQ_FORMAT = FASTQ_PHRED33;
	}
	else
		tmpfile[0]=0;



	if(!EXON_USE_VALUE_ARRAY_INDEX)
	{
		printf("Detecting junction reads must reference to the base-wise indel. Please enable the '-b' option while building the index.\n");
		return -1;
	}
  
	reads_density = guess_reads_density(read_file, IS_SAM_INPUT);
	if(reads_density<0)
		printf("Input file '%s' is not found or is in an incorrect format.\n", read_file);

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
		if(geinput_open_sam(read_file, &ginp,1))
			return -1;
		if(geinput_open_sam(read_file, &ginp2,2))
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
		printf("Unable to open the index files in the %s space.\n",  ginp.space_type==GENE_SPACE_COLOR?"color":"base");
		return -1;
	}

	FILE * out_fp = NULL;
	if(REPORT_SAM_FILE)
		out_fp= fopen(output_file, "w");
	if (REPORT_SAM_FILE && !out_fp)
	{
		printf("Unable to open the output file at '%s'.\n", output_file);
		return -1;
	}


	//printf("Number of subreads selected for each read=%d\n", TOTAL_SUBREADS);
	printf("Threshold on number of subreads for a successful mapping=%f\n", EXON_MAJOR_HALF_VOTES_RATE);
	printf("Number of threads=%d\n", EXON_ALL_THREADS);
	if (EXON_INDEL_TOLERANCE)
		printf("Tolerance for Indel=%d\n", EXON_INDEL_TOLERANCE-1);
	if (EXON_QUALITY_SCALE==QUALITY_SCALE_LINEAR)
		puts("Quality scale=linear\n\n");
	else if (EXON_QUALITY_SCALE==QUALITY_SCALE_LOG)
		puts("Quality scale=exponential\n\n");
	else 	puts("\n");

	if (read2_file[0] || IS_SAM_INPUT==2)
	{
		if (EXON_MAX_PAIRED_DISTANCE <= EXON_MIN_PAIRED_DISTANCE)
		{
			printf ("The value of the '-D' option must be greater than that of the '-d' option. \n");
			return -1;
		}

		printf ("Performing paired-end alignment:\n");
		printf ("Maximum distance between reads=%d\n", EXON_MAX_PAIRED_DISTANCE);
		printf ("Minimum distance between reads=%d\n", EXON_MIN_PAIRED_DISTANCE);
		printf ("Threshold on number of subreads for a successful mapping (the minor end in the pair)=%d\n", ACCEPT_MINOR_SUBREADS);
		printf ("The directions of the two input files are: %s, %s\n\n", EXON_FIRST_READ_REVERSE?"reversed":"forward", EXON_SECOND_READ_REVERSE?"reversed":"forward");
	}

#ifdef REPORT_ALL_THE_BEST
	printf("***** WARNING: the REPORT_ALL_THE_BEST switch is turned on. You need an extra 1 GBytes of RAM space for saving the temporary results. *****\n");
#endif

	init_halves_record(&halves_record, all_reads);
	if(read2_file[0])
		all_reads/=2;
		
	if((!IS_SAM_INPUT) && read2_file[0] && geinput_open(read2_file, &ginp2))
	{
		printf("Input file '%s' is not found or is in an incorrect format.\n", read2_file);
		return -1;
	}

	fflush(stdout);

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


	while (1)
	{
		char inbuff[1201];

		int new_processed_reads = run_exon_search_index(&ginp, read2_file[0] ? (&ginp2):NULL, index_prefix, &halves_record, out_fp, processed_reads, all_tables, &succeed_reads, bed_index, pos_index, connection_index);
		if(new_processed_reads<0)break;

		processed_reads += new_processed_reads;
		clear_halve_record(&halves_record);

		// test if there no anyreads remaining
		unsigned long long int current_fp = ftello(ginp.input_fp);
		int rl = geinput_next_read(&ginp, NULL, inbuff, NULL);
		//printf("RL=%d\nfp=%llu\npr=%u",rl, current_fp, processed_reads);
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
		printf("@LOG THE END. \n");
	else
		printf("\n\n %llu reads were processed in %.1f seconds.\n There are %llu junction pairs found, supported by %llu reads.\n\n", processed_reads, miltime()-begin_ftime, junction_number, support_number );

	if(out_fp)
		fclose(out_fp);
	printf("\n\nCompleted successfully.\n");

	if(tmpfile[0])
		unlink(tmpfile);
	return 0;
}
