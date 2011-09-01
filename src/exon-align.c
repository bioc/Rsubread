#include <stdio.h>
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
float accepted_support_rate = 0.300000;
int EXON_ALL_THREADS=1;
int TOTAL_SUBREADS;
int ACCEPT_SUBREADS;
int ACCEPT_MINOR_SUBREADS;
int INDEX_THRESHOLD;
int EXON_MAX_PAIRED_DISTANCE = 600;
int EXON_MIN_PAIRED_DISTANCE = 50;
int EXON_INDEL_TOLERANCE = 0;
int EXON_QUALITY_SCALE = 0;
int EXON_USE_VALUE_ARRAY_INDEX = 1;
int EXON_FIRST_READ_REVERSE = 0;
int EXON_SECOND_READ_REVERSE = 1;
int EXON_NUMBER_OF_ANCHORS_PAIRED = 50;
double reads_density;


#define IS_LONG_OVERLAP 4
#define IS_SHORT_OVERLAP 8
#define IS_PAIRED_HINTED 16
#define IS_R1_CLOSE_TO_5 1
#define IS_REVERSED_HALVES 2
#define	IS_PROCESSED_READ 32

#define MIN_HALF_VOTES 2
#define DONAR_CONFIRM_SIZE 13
#define EXON_GROUPING_SIZE 12

typedef struct {
	int max_len;

	unsigned int * best_pos1_list;
	unsigned int * best_pos2_list;
	unsigned char * best_vote1_list;
	unsigned char * best_vote2_list;
	char * is_abnormal_list;
	char * is_reversed_list;
	char * half_marks_list;
	short * splice_point_list;
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
	halves_record -> half_marks_list = (char * )malloc(sizeof(char)*items);
	halves_record -> splice_point_list = (short * )malloc(sizeof(short)*items);

	memset(halves_record -> best_vote1_list, 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> best_vote2_list, 0, sizeof(char)*halves_record -> max_len);
}

void clear_halve_record(halves_record_t* halves_record)
{
	memset(halves_record -> best_vote1_list, 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> best_vote2_list, 0, sizeof(char)*halves_record -> max_len);
	memset(halves_record -> half_marks_list, 0, sizeof(char)*halves_record -> max_len);
}

void add_best_matching_halves(halves_record_t * halves_record, unsigned int best_pos1, unsigned int best_pos2, unsigned char best_vote1 , unsigned char best_vote2, char is_abnormal, char is_reversed, short splice_point, char half_marks, int read_number)
{
	if (halves_record -> best_vote1_list[read_number] < best_vote1)
	{
		halves_record -> best_pos1_list[read_number] = best_pos1;
		halves_record -> best_pos2_list[read_number] = best_pos2;
		halves_record -> best_vote1_list[read_number] = best_vote1;
		halves_record -> best_vote2_list[read_number] = best_vote2;
		halves_record -> is_abnormal_list[read_number] = is_abnormal;
		halves_record -> is_reversed_list[read_number] = is_reversed;
		halves_record -> half_marks_list[read_number] = half_marks;
		halves_record -> splice_point_list[read_number] = splice_point;
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

int select_best_matching_halves(gene_vote_t * vote, unsigned int * best_pos1, unsigned int * best_pos2, int * best_vote1, int * best_vote2, char * is_abnormal, char * half_marks, char is_reversed, float accept_rate, int read_len, long long int hint_pos)
{
	int best_splicing_point = -1, i,j;
	char * best_chro_name;
	unsigned int best_chro_pos;
	int selected_max_votes = -1;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			int n1, n2;
			char * chro_name;
			int long_overlap = 0;
			unsigned int chro_pos;

			// All logical conditions
			int overlap_start = max(vote->max_coverage_start , vote->coverage_start[i][j]);
			int overlap_end   = min(vote->max_coverage_end   , vote->coverage_end[i][j]  );
			int overlapped_len = overlap_end - overlap_start;
		//	printf("\nTest Cov : %d ~ %d vs %d ~ %d = %d", vote->max_coverage_start, vote->max_coverage_end ,  vote->coverage_start[i][j],  vote->coverage_end[i][j] , overlapped_len);

			if(overlapped_len >=14)
				continue;
			else if (overlapped_len >=7)
				long_overlap = IS_LONG_OVERLAP | IS_SHORT_OVERLAP;
			else if (overlapped_len >=3)
				long_overlap = IS_SHORT_OVERLAP;


			long long int dist = vote->pos[i][j];
			dist -= vote->max_position;

			//printf ("D=%lld\n", abs(dist));
			if (abs(dist)<6)
				continue;

			// All quality conditions 
			if (vote->coverage_start[i][j] < vote->max_coverage_start)
				n1 = read_len - vote->max_coverage_start;
			else
				n1 = vote->max_coverage_end;
			n2 = read_len - n1;

			int support_r1 = (int)(max(MIN_HALF_VOTES+0.01, (n1*1.-15)/3*accept_rate*1.0001));
			int support_r2 = (int)(max(MIN_HALF_VOTES+0.01, (n2*1.-15)/3*accept_rate*1.0001));

			//printf ("%u and %u : V1 %d >= S1 %d\t\tV2 %d >= S2 %d\n", vote->max_position, vote->pos[i][j] , vote->max_vote, support_r1,  vote->votes[i][j], support_r2);

			if (vote->max_vote < support_r1 || vote->votes[i][j]<support_r2)
				continue;

			// Same chromosome
			if ((vote->coverage_start[i][j] < vote->max_coverage_start) + is_reversed == 1)
			{
				locate_gene_position(vote->max_position + read_len, &_global_offsets, &best_chro_name, &best_chro_pos);
				locate_gene_position(vote->pos[i][j] , &_global_offsets, &chro_name, &chro_pos);
			}else
			{
				locate_gene_position(vote->max_position , &_global_offsets, &best_chro_name, &best_chro_pos);
				locate_gene_position(vote->pos[i][j] +read_len, &_global_offsets, &chro_name, &chro_pos);
			}
			if (chro_name != best_chro_name)	// The pointers can be compared because they can be the same.
				continue;


			int test_vote_value = vote->votes[i][j];
			if (hint_pos>=0)
			{
				long long int hint_dist = hint_pos;
				hint_dist -= vote->pos[i][j];
				if (abs (hint_dist) < 1000000)
					test_vote_value += 100;
			}
			if (test_vote_value<=selected_max_votes)continue;
			// Conditions of order of R3 and R5
			*half_marks &= ~IS_REVERSED_HALVES;
			if (vote->coverage_start[i][j] < vote->max_coverage_start && (((vote->max_position < vote->pos[i][j]) && !is_reversed) || ((vote->max_position > vote->pos[i][j]) && is_reversed) ) )
				*half_marks |= IS_REVERSED_HALVES;
			if (vote->coverage_start[i][j] >= vote->max_coverage_end  &&  (((vote->max_position > vote->pos[i][j]) && !is_reversed) || ((vote->max_position < vote->pos[i][j]) && is_reversed) ) )
				*half_marks |= IS_REVERSED_HALVES;

			if (vote->coverage_start[i][j] < vote->max_coverage_start)
				*half_marks = (*half_marks) & ~IS_R1_CLOSE_TO_5;
			else
				*half_marks |= IS_R1_CLOSE_TO_5;
			
			best_splicing_point = ((vote->coverage_start[i][j] < vote->max_coverage_start)? (vote->coverage_end[i][j]):(vote->max_coverage_end)) + ((vote->coverage_start[i][j] < vote->max_coverage_start)? (vote->max_coverage_start):(vote->coverage_start[i][j]));
			best_splicing_point /=2;
			/*printf("B %u, %u\n", vote->max_coverage_start, vote->max_coverage_end);
			printf("S %u, %u\n", vote->coverage_start[i][j], vote->coverage_end[i][j]);
			printf("SP = %d\n\n", best_splicing_point);
*/
			* best_pos1 = vote->max_position ;
			* best_pos2 = vote->pos[i][j] ;
			* best_vote1 = vote->max_vote ;
			* best_vote2 = vote->votes[i][j] ;

			*half_marks = (*half_marks) & ~(IS_PAIRED_HINTED);
			if (test_vote_value >=100)
				*half_marks = (*half_marks) | IS_PAIRED_HINTED;
			*half_marks = (*half_marks) & ~(IS_LONG_OVERLAP|IS_SHORT_OVERLAP);
			if (long_overlap)
				*half_marks |= long_overlap;
			selected_max_votes = test_vote_value; 

		}

	return best_splicing_point;


}


#define ceq(c,t) ((c)[0]==(t)[0] && (c)[1]==(t)[1])
#define c2eq(ch1, ch2, tg1, tg2) ((ceq(ch1, tg1) && ceq(ch2, tg2)) || (ceq(ch1, tg2) && ceq(ch2, tg1)) )

int paired_chars(char * ch1, char * ch2, int is_reverse)
{
	if (c2eq(ch1, ch2, "GT", "AG") || c2eq(ch1, ch2, "CT", "AC")|| c2eq(ch1, ch2,"GC","AG") || c2eq(ch1, ch2,"GC","CT") || c2eq(ch1, ch2,"AT","AC") || c2eq(ch1, ch2,"GT","AT"))
	{
		if (is_reverse) if (ceq(ch1, "AG") || ceq(ch1, "AC") ||  ceq(ch1, "GC") || ceq(ch1, "AT")) return 1;
		if (!is_reverse) if (ceq(ch1, "CT") || ceq(ch1, "GT") || ceq(ch1, "GC") || ceq(ch1, "AT")) return 1;
	}
	return 0;
}

#define is_donar_chars(cc) (((cc)[0]=='G' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='G') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') || \
			    ((cc)[0]=='C' && (cc)[1]=='T') || \
			    ((cc)[0]=='G' && (cc)[1]=='C') || \
			    ((cc)[0]=='A' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') ) 

int match_chro(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand)
{
	int ret = 0;
	int i;
	if (is_negative_strand)
		for (i=test_len -1;i>=0;i--)
		{
			char tt = gvindex_get (index, pos+test_len-1-i);
			switch(tt)
			{
				case 'A': ret += read[i] == 'T'; break;
				case 'T': ret += read[i] == 'A'; break;
				case 'G': ret += read[i] == 'C'; break;
				case 'C': ret += read[i] == 'G'; break;
			}
		}
	else
		for (i=0;i<test_len;i++)
			ret +=read[i] == gvindex_get (index, pos +i);
	return ret;
}



void get_chro(char *buf, gene_value_index_t * index, unsigned int pos, int is_negative_strand)
{
	int i;
	if (is_negative_strand)
		for (i=1;i>=0;i--)
		{
			buf[i] = gvindex_get (index, pos+1-i);
			switch(buf[i])
			{
				case 'A': buf[i] = 'T'; break;
				case 'T': buf[i] = 'A'; break;
				case 'G': buf[i] = 'C'; break;
				case 'C': buf[i] = 'G'; break;
			}
		}
	else
		for (i=0;i<2;i++)
			buf[i] = gvindex_get (index, pos +i);
}

int test_donar(char *read, int read_len, unsigned int pos1, unsigned int pos2, int guess_break_point, char negative_strand, int test_range, char is_soft_condition, int EXON_INDEL_TOLERANCE, int* real_break_point, gene_value_index_t * my_value_array_index, int indel_offset, int is_reversed)
{
	int bps_pos_x;
	int search_start = guess_break_point - test_range ;
	int search_end   = guess_break_point + test_range ;
	char h1_2ch[3], h2_2ch[3];

	h1_2ch[2] = h2_2ch[2]=0;
	search_start=max(16, search_start);
	search_end = min(read_len-16, search_end);
	int best_break = -1;
	int min_x = -9099;

	for (bps_pos_x = search_start; bps_pos_x < search_end ; bps_pos_x ++)
	{
		get_chro(h1_2ch, my_value_array_index, pos1 + bps_pos_x, is_reversed);
		get_chro(h2_2ch, my_value_array_index, pos2 - 2 + indel_offset + bps_pos_x, is_reversed);

		//if(!is_reversed)
		//printf("C1=%s @%u, C2=%s @%u\n",h1_2ch, pos1 + bps_pos_x, h2_2ch,pos2 - 2 + indel_offset + bps_pos_x);
		if(h1_2ch[0]==h2_2ch[0] && h1_2ch[1]==h2_2ch[1]) continue;

		if(is_donar_chars(h1_2ch) && is_donar_chars(h2_2ch))
			if(paired_chars(h1_2ch, h2_2ch, is_reversed))
			{
				int m1, m2, x1, x2;
				int break_point_half = is_reversed?(read_len - bps_pos_x):bps_pos_x;
				int first_exon_end,second_half_start;

				if (is_reversed)
				{
					first_exon_end = pos2 + bps_pos_x + indel_offset;
					second_half_start = pos1 + bps_pos_x;

					m1 = match_chro(read + break_point_half - DONAR_CONFIRM_SIZE -3, my_value_array_index, first_exon_end, DONAR_CONFIRM_SIZE +3, is_reversed);
					m2 = match_chro(read + break_point_half , my_value_array_index, second_half_start-DONAR_CONFIRM_SIZE -3, DONAR_CONFIRM_SIZE +3, is_reversed);

					x1 = match_chro(read + break_point_half ,  my_value_array_index, first_exon_end - DONAR_CONFIRM_SIZE, DONAR_CONFIRM_SIZE , is_reversed);
					x2 = match_chro(read + break_point_half -DONAR_CONFIRM_SIZE ,  my_value_array_index, second_half_start , DONAR_CONFIRM_SIZE , is_reversed);
				}
				else
				{
					first_exon_end = pos1 + bps_pos_x + indel_offset;
					second_half_start = pos2 + bps_pos_x;

					m1 = match_chro(read + break_point_half - DONAR_CONFIRM_SIZE -3, my_value_array_index, first_exon_end-DONAR_CONFIRM_SIZE -3, DONAR_CONFIRM_SIZE +3, is_reversed);
					m2 = match_chro(read + break_point_half , my_value_array_index, second_half_start, DONAR_CONFIRM_SIZE +3 , is_reversed);

					x1 = match_chro(read + break_point_half ,  my_value_array_index, first_exon_end, DONAR_CONFIRM_SIZE , is_reversed);
					x2 = match_chro(read + break_point_half -DONAR_CONFIRM_SIZE ,  my_value_array_index, second_half_start - DONAR_CONFIRM_SIZE, DONAR_CONFIRM_SIZE , is_reversed);
					//x1 = match_chro(read + break_point_half ,  my_value_array_index, first_exon_end-DONAR_CONFIRM_SIZE, DONAR_CONFIRM_SIZE , is_reversed);
					//x2 = match_chro(read + break_point_half -DONAR_CONFIRM_SIZE ,  my_value_array_index, second_half_start , DONAR_CONFIRM_SIZE , is_reversed);
				}
		//		printf("M1=%d M2=%d X1=%d X2=%d\n", m1,m2,x1,x2);
				if (m1 >= DONAR_CONFIRM_SIZE - 1  && m2>=DONAR_CONFIRM_SIZE - 1)
					if(x1<=DONAR_CONFIRM_SIZE - 5 && x2<=DONAR_CONFIRM_SIZE - 5)
					{
						int score = 40 - max(x1, x2) + min(m1, m2);
						if (min_x < score)
						{
							min_x = score;
							best_break = bps_pos_x;
						}
					}

			}

	}
	if (best_break>0)
	{
		//printf ("FINAL BREAK: %d   ; REV = %d\n ", best_break, is_reversed);
		*real_break_point = best_break;
		return 1;
	}
	return 0;
}

unsigned int get_grouped_position(HashTable * pos_table, unsigned int pos)
{

	int delta_pos;
	int group_anchor = pos / EXON_GROUPING_SIZE - 1;
	unsigned int grouped_pos = 0;
	for (delta_pos = 0 ; delta_pos < 3; delta_pos ++)
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

typedef struct{
	unsigned int small_key;
	unsigned int big_key;
} paired_exon_key;

void feed_exonbed(HashTable * bed_table, HashTable * pos_table,  halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads, gene_value_index_t * my_value_array_index, int first_index, int last_index)
{
	int i=0;
	while (1)
	{
		char nameb[1201], inb[1201], qualityb[1201];
		int rl;
		if(i >= processed_reads*2)break;
		unsigned int pos = ( IS_R1_CLOSE_TO_5 & halves_record -> half_marks_list[i] ) ?halves_record -> best_pos1_list[i]:halves_record -> best_pos2_list[i];
		unsigned int pos2 = ( IS_R1_CLOSE_TO_5 & halves_record -> half_marks_list[i] ) ?halves_record -> best_pos2_list[i]:halves_record -> best_pos1_list[i];


		if (ginp2 && (i % 2))
			rl = geinput_next_read(ginp2, nameb, inb, qualityb);
		else
			rl = geinput_next_read(ginp, nameb, inb, qualityb);
		if (rl<0){
			break;
		}

		if (!first_index)
		if (min(pos, pos2) < my_value_array_index -> start_base_offset + 1000)
		{
			i++;
			continue;
		}

		if (!last_index)
		if (max(pos, pos2)  > my_value_array_index -> start_base_offset+ my_value_array_index ->length - 1000)
		{
			i++; continue;
		}

		if (IS_PROCESSED_READ & halves_record -> half_marks_list[i])
		{
			i++;
			continue;
		}

		long long int dist = halves_record -> best_pos1_list[i] ;
		dist-= halves_record -> best_pos2_list[i];

		if(halves_record -> best_vote1_list[i]>=2 && halves_record -> best_vote2_list[i]>=1 )
		{
			//char * chro_name;
			//unsigned int chro_pos, chro_pos_small, chro_pos_large;
			//char is_strong = (min(halves_record ->best_vote1_list[i],halves_record ->best_vote2_list[i])>5);
			char short_overlap = !(halves_record -> half_marks_list[i] & IS_LONG_OVERLAP);
			char negative_strand = halves_record -> is_reversed_list[i];
			int real_break_point;
			char is_soft_condition = 0;
			int test_range = rl / 4;
			char is_accepted = 0;
			int guess_break_point = (halves_record -> is_reversed_list[i]) ? (rl - halves_record ->splice_point_list[i]) : halves_record ->splice_point_list[i];

		//	pos = halves_record -> is_reversed_list[i]? (pos + rl):(pos);
		//	pos2 = halves_record -> is_reversed_list[i]? (pos2):(pos2+rl);

		//	printf ("%s\n", inb);

			if (min(halves_record ->best_vote1_list[i],halves_record ->best_vote2_list[i]) >=4 && short_overlap)
			{
				int indel_x;
				for(indel_x = 0; indel_x < 2*EXON_INDEL_TOLERANCE +1 ; indel_x ++)
				{
					int indel_offset = ((indel_x %2)?1:-1) * ((indel_x+1)>>1);
					is_accepted = test_donar(inb, rl, min(pos, pos2), max(pos,pos2), guess_break_point, negative_strand, test_range, is_soft_condition, EXON_INDEL_TOLERANCE, &real_break_point, my_value_array_index, indel_offset, negative_strand);
					if (is_accepted) break;
				}
			}
			else
				is_accepted = test_donar(inb, rl, min(pos, pos2), max(pos,pos2), guess_break_point, negative_strand, test_range, is_soft_condition, EXON_INDEL_TOLERANCE, &real_break_point, my_value_array_index,0, negative_strand);

			if (is_accepted)
			{
				// real_break_point is on chromosome: pos+real_break_point / pos2 + real_break_point are end/start of exons
				unsigned int pos_small = min(pos, pos2) + real_break_point;
				unsigned int pos_big   = max(pos, pos2) + real_break_point;
				pos_small = get_grouped_position(pos_table, pos_small);
				pos_big = get_grouped_position(pos_table, pos_big);

				paired_exon_key search_key; 
				search_key.small_key = pos_small;
				search_key.big_key = pos_big;

				int *search_res = (int*) HashTableGet(bed_table, &search_key);
				if(!search_res)
				{
					search_res = (int *)malloc(sizeof(int));
					*search_res = 0;
					paired_exon_key * new_key = (paired_exon_key *) malloc(sizeof(paired_exon_key));
					new_key->small_key = pos_small;
					new_key->big_key = pos_big;
					HashTablePut(bed_table, new_key, search_res);
				}
				
				(*search_res)++;
				
			}
		}
		halves_record -> half_marks_list[i] |= IS_PROCESSED_READ;
		i++;
	}
}

void print_exon_res(halves_record_t * halves_record,  gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads)
{
	int i=0;


	while (1)
	{
		char nameb[1201], inb[1201], qualityb[1201];
		int rl;
		if(i >= processed_reads*2)break;
		if (ginp2 && (i % 2))
			rl = geinput_next_read(ginp2, nameb, inb, qualityb);
		else
			rl = geinput_next_read(ginp, nameb, inb, qualityb);
		if (rl<0){
			break;
		}


		long long int dist = halves_record -> best_pos1_list[i] ;
		dist-= halves_record -> best_pos2_list[i];

		if(halves_record -> best_vote1_list[i]>=2 && halves_record -> best_vote2_list[i]>=1 )
		{
			char * chro_name;
			unsigned int chro_pos, chro_pos_small, chro_pos_large;
			char * diststr =  abs(dist)<1000000?(abs(dist)<200000?(abs(dist)<5000?"NN1":"NM1"):"NF1"):"FF1"; 
			char * reversed_halves = (IS_REVERSED_HALVES &  halves_record -> half_marks_list[i] )? "RE2":"NM2";
			char * strong_votes = (min(halves_record ->best_vote1_list[i],halves_record ->best_vote2_list[i])>5)?"ST3":"WK3";
			char * long_overlap = (halves_record -> half_marks_list[i] & IS_SHORT_OVERLAP) ? ((halves_record -> half_marks_list[i] & IS_LONG_OVERLAP)?"LO4":"SO4"):"NO4";
			char * negative_strand = halves_record -> is_reversed_list[i]? "NE5":"PO5";
			char * pair_info_refed = (IS_PAIRED_HINTED & halves_record -> half_marks_list[i]) ? ((i % 2)?"PB7" : "PA7"):"SG7";

			unsigned int pos = ( IS_R1_CLOSE_TO_5 & halves_record -> half_marks_list[i] ) ?halves_record -> best_pos1_list[i]:halves_record -> best_pos2_list[i];
			unsigned int pos2 = ( IS_R1_CLOSE_TO_5 & halves_record -> half_marks_list[i] ) ?halves_record -> best_pos2_list[i]:halves_record -> best_pos1_list[i];

			locate_gene_position(min(pos,pos2), &_global_offsets, &chro_name, &chro_pos);

			pos = halves_record -> is_reversed_list[i]? (pos + rl):(pos);
			pos2 = halves_record -> is_reversed_list[i]? (pos2):(pos2+rl);

			locate_gene_position(pos, &_global_offsets, &chro_name, &chro_pos_small) ;
			locate_gene_position(pos2, &_global_offsets, &chro_name, &chro_pos_large) ;

			int h1_len = (halves_record -> is_reversed_list[i]) ? (rl - halves_record ->splice_point_list[i]) : halves_record ->splice_point_list[i];


			fprintf(out_fp,"%s\t%d\t%s\t%u\t%d\t%dM%dN%dM\t*\t0\t0\t%s\t%s\t%s,%s,%s,%s,%s,JN6,%s\t%u\t%u\n",  nameb, halves_record -> is_reversed_list[i]?16:0 , chro_name, chro_pos, min(halves_record ->best_vote1_list[i],halves_record ->best_vote2_list[i]) , h1_len , (int)abs(dist), rl - h1_len , inb, qualityb, diststr,reversed_halves,strong_votes,long_overlap,negative_strand, pair_info_refed, chro_pos_small, chro_pos_large);

		}
		else if (halves_record -> best_vote1_list[i]>=8)
		{
			char * chro_name;
			unsigned int chro_pos;
			unsigned int pos = halves_record -> best_pos1_list[i];
			locate_gene_position(pos, &_global_offsets, &chro_name, &chro_pos);
			//fprintf(out_fp,"%s\t%s\t%u\t%s\tEXONIC\n", nameb, chro_name, chro_pos, halves_record -> is_reversed_list[i]?"NEG":"POS");
			char * pair_info_refed = (IS_PAIRED_HINTED & halves_record -> half_marks_list[i]) ? "PA7":"SG7";
			fprintf(out_fp,"%s\t%d\t%s\t%u\t%d\t%dM\t*\t0\t0\t%s\t%s\tEX6,%s\n",  nameb, halves_record -> is_reversed_list[i]?16:0 , chro_name, chro_pos, halves_record -> best_vote1_list[i] , rl, inb, qualityb, pair_info_refed);
		}
		else
		{
			//fprintf(out_fp,"%s\tUNMAPPED\n",nameb);

			fprintf(out_fp,"%s\t4\t*\t*\t0\t*\t*\t0\t0\t%s\t%s\tNM6\n",  nameb,inb, qualityb);
		}
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

	pthread_spinlock_t * input_data_lock;
	pthread_spinlock_t * init_lock;
};

int run_exon_search(gehash_t * my_table, gene_value_index_t * my_value_array_index , int table_no,  halves_record_t * halves_record, gene_input_t * ginp,gene_input_t * ginp2, char * index_prefix, unsigned int * processed_reads, long long int base_number, int all_tables, pthread_spinlock_t * input_lock, int my_thread_no, unsigned int section_length)
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
	float subread_step = 3.0001;
	struct stat read_fstat;
	stat (ginp->filename, &read_fstat);
	long long int read_fsize = read_fstat.st_size;
	
	double local_begin_ftime = miltime();
	int read_len = 0, read2_len = 0;
	is_reversed = 0;


	if (ginp2)
		all_reads /=2;

//	printf ("I'm the %d-th thread\n", my_thread_no);


	while (1)
	{

		char namebuf[200];

		if (is_reversed)
		{
			
			reverse_read(InBuff, read_len, ginp->space_type);	
			reverse_quality(QualityBuff, read_len);
			if (ginp2)
			{
				reverse_read(InBuff2, read2_len, ginp->space_type);
				reverse_quality(QualityBuff2, read2_len);
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
//			subread_step = max(3.00001, (read_len-16-GENE_SLIDING_STEP)*1.0/(TOTAL_SUBREADS-1)+0.00001); 
		}
		if (ginp2)
		{

			int i;
			int subread_no;
			int is_second_read;
			gene_vote_t vote_p1, vote_p2;
			long long int hint_pos;
			init_gene_vote(&vote_p1);
			init_gene_vote(&vote_p2);

			for (is_second_read = 0; is_second_read <2; is_second_read ++)
			{
				gene_vote_t * current_vote = is_second_read?&vote_p2: &vote_p1;
				char * current_read =  is_second_read?InBuff2 : InBuff;
				int current_rlen = is_second_read?read2_len:read_len;
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

						gehash_go_q(my_table, subread_integer , subread_offset, read_len, is_reversed, current_vote, 1, 1, 21.9, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no);
					}
					if(subread_offset1 >= current_rlen -16)
						break;
				}
			}

			gene_vote_number_t numvote_read1, numvote_read2;
			gehash_data_t pos_read1, pos_read2;
			gene_quality_score_t sum_quality, qual_r1, qual_r2;
			char record_index1[48], record_index2[48];

			int is_paired_match = select_positions(&vote_p1, &vote_p2, &numvote_read1, &numvote_read2, &sum_quality, &qual_r1, &qual_r2, &pos_read1, &pos_read2, record_index1, record_index2, EXON_MAX_PAIRED_DISTANCE, EXON_MIN_PAIRED_DISTANCE, ACCEPT_SUBREADS, ACCEPT_MINOR_SUBREADS, is_reversed, EXON_NUMBER_OF_ANCHORS_PAIRED,  EXON_INDEL_TOLERANCE);
			for (is_second_read = 0; is_second_read <2; is_second_read ++)
			{
				gene_vote_t * current_vote = is_second_read?&vote_p2: &vote_p1;
				unsigned int best_pos1=0;
				unsigned int best_pos2=0;
				hint_pos = is_second_read?pos_read2:pos_read1;
				if (!is_paired_match ) hint_pos = -1;
				int best_vote1=0;
				int best_vote2=0;
				char is_abnormal=0;
				char half_marks=0;
		
				int splice_point = select_best_matching_halves(current_vote, &best_pos1, &best_pos2, &best_vote1, &best_vote2, &is_abnormal ,&half_marks, is_reversed, accepted_support_rate, read_len, hint_pos);
				if (splice_point>0)
					add_best_matching_halves(halves_record, best_pos1, best_pos2, best_vote1, best_vote2,is_abnormal, is_reversed, splice_point, half_marks, 2*queries+is_second_read);
				else 
				{
					if (current_vote->max_vote > (halves_record -> best_vote1_list[queries*2+is_second_read] + halves_record -> best_vote2_list[queries*2+is_second_read]))
					{
						halves_record -> best_vote1_list[queries*2+is_second_read] = current_vote->max_vote ;
						halves_record -> best_pos1_list[queries*2+is_second_read] = current_vote->max_position;
						halves_record -> is_reversed_list[queries*2+is_second_read] = is_reversed;
						halves_record -> half_marks_list[queries*2+is_second_read] = (halves_record -> half_marks_list[queries*2+is_second_read]) & ~(IS_PAIRED_HINTED);
						if (is_paired_match)
							halves_record -> half_marks_list[queries*2+is_second_read] = (halves_record -> half_marks_list[queries*2+is_second_read]) | IS_PAIRED_HINTED;
	
					}
				}
			}
		
		}
		else
		{

			int i;
			int subread_no;
			gene_vote_t vote;
			init_gene_vote(&vote);
		
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

					gehash_go_q(my_table, subread_integer , subread_offset, read_len, is_reversed, &vote, 1, 1, 21.9, INDEX_THRESHOLD, EXON_INDEL_TOLERANCE, subread_no);
				}
				if(subread_offset1 >= read_len -16)
					break;
			}

			unsigned int best_pos1=0;
			unsigned int best_pos2=0;
			int best_vote1=0;
			int best_vote2=0;
			char is_abnormal=0;
			char half_marks=0;
	
			int splice_point = select_best_matching_halves(&vote, &best_pos1, &best_pos2, &best_vote1, &best_vote2, &is_abnormal ,&half_marks, is_reversed, accepted_support_rate, read_len, -1);
			if (splice_point>0)
			{
				add_best_matching_halves(halves_record, best_pos1, best_pos2, best_vote1, best_vote2,is_abnormal, is_reversed, splice_point, half_marks, queries);
			}
			else 
			{
				if (vote.max_vote > (halves_record -> best_vote1_list[queries] + halves_record -> best_vote2_list[queries]))
				{
					halves_record -> best_vote1_list[queries] = vote.max_vote ;
					halves_record -> best_pos1_list[queries] = vote.max_position;
					halves_record -> is_reversed_list[queries] = is_reversed;
				}
			}
		
		}
		if (my_thread_no==0 && queries % 10000 == 0 && !is_reversed)
		{
			if(table_no == 0)
			{
				long long int current_reads = base_number + queries;
				long long int fpos = ftello(ginp->input_fp);
				reads_density = fpos*1.0/current_reads; 
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

	run_exon_search(data_param->my_table, data_param->my_value_array_index, data_param->table_no, data_param->halves_record, data_param->ginp , data_param->ginp2, data_param->index_prefix, data_param->processed_reads, data_param->base_number, data_param->all_tables, data_param->input_data_lock, thid, data_param-> section_length);
	return NULL;
}



// This function search a segment of reads (length = read_number) 
// It returns the number of reads that were really processed;
int run_exon_search_index(gene_input_t * ginp, gene_input_t * ginp2, char * index_prefix, halves_record_t * halves_record, FILE * out_fp, unsigned long long int base_number, int all_tables, unsigned long long int *succeed_reads, HashTable * overall_exon_bed, HashTable * pos_table)
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
			run_exon_search(my_table, &value_array_index, tabno, halves_record, ginp , ginp2, index_prefix, &processed_reads, base_number, all_tables, NULL /*the data lock is null*/, 0  /*I'm the 0-th thread*/, section_length);
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
			data_param.ginp2 = ginp2;

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

	print_exon_res(halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads);
	fseeko(ginp -> input_fp, current_fp, SEEK_SET);
	if (ginp2)
		fseeko(ginp2 -> input_fp, current_fp2, SEEK_SET);

	if (tabno == 1)
		feed_exonbed(overall_exon_bed,pos_table, halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads, &value_array_index, 1,1);
	else
	{
		gvindex_destory(&value_array_index);
		tabno = 0;
		while (1)
		{
			sprintf(table_fn, "%s.%02d.%c.array", index_prefix, tabno, ginp->space_type==GENE_SPACE_COLOR?'c':'b');

			stat_ret = stat(table_fn, &filestat);
			if (stat_ret !=0)
				break;

			gvindex_load(&value_array_index,table_fn);
			feed_exonbed(overall_exon_bed, pos_table, halves_record, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads, &value_array_index, tabno==0, tabno== last_table);
			tabno ++;

			fseeko(ginp -> input_fp, current_fp, SEEK_SET);
			if (ginp2)
				fseeko(ginp2 -> input_fp, current_fp2, SEEK_SET);

		}
	}
	if (tabno>0)gvindex_destory(&value_array_index);


	fseeko(ginp -> input_fp, last_fp_pos, SEEK_SET);
	if (ginp2)
		fseeko(ginp2 -> input_fp, last_fp2_pos, SEEK_SET);

	return processed_reads;
	
}

void exon_usage(char * execname)
{
	puts("Usage:");
	puts(" ./subjunc [options] -i <index_name> -r <input> -o <output>");
	puts("");
	puts("Basic arguments:");
	puts("    -i --index     <index>\t name of the index, same as that of the index builder.");
	puts("    -r --read      <input>\t name of an input file(FASTQ/FASTA format), either in the base-space or in the color-space.");
	puts("    -o --output    <output>\t name of the output file(SAM format)");
	puts("");
	puts("Optional arguments:");
//	puts("    -n --subreads  <int>\t optional, number of subreads selected from each read for mapping, 10 by default");
//	puts("    -m --minmatch  <int>\t optional, minimal number of subreads which have the consensus mapping location, 3 by default");
	puts("    -T --threads   <int>\t optional, number of threads/CPUs used for mapping the reads, 1 by default");
	puts("    -I --indel     <int>\t optional, the maximum number of bases for insertion/deletion, 0 by default");
//	puts("    -Q --quality   <l:e:n>\t optional, quality scale: l for linear, e for exponential, n for none.");
//	puts("    -a --basewise       \t optional, using the base-wise quality index; the index must be built with a -a option.");
	puts("");
	puts("Paired-end alignment arguments:");
	puts("    -R --read2     <input>\t optional, the second input file; using this argument to activate paired-end alignment");
//	puts("    -p --minmatch2 <int>\t optional, the `-m' option for the read receiving less votes in a pair to be accepted, 1 by default");
	puts("    -d --mindist   <int>\t optional, the minimum distance between two reads in a pair, 50 by default");
	puts("    -D --maxdist   <int>\t optional, the maximum distance between two reads in a pair, 600 by default");
	puts("    -S --order     <ff:fr:rf> \t optional, specifying if the first/second reads are forward or reversed, 'fr' by default.");
	puts("");
	puts("Example:");
	puts(" ./subjunc -i my_index -r reads.fastq -o my_result.sam ");
	puts("");
	puts("");

}

static struct option long_options[] =
{
	{"basewise", no_argument, &EXON_USE_VALUE_ARRAY_INDEX, 1},
	{"index", required_argument, 0, 'i'},
	{"read",  required_argument, 0, 'r'},
	{"read2", required_argument, 0, 'R'},
	{"indel", required_argument, 0, 'I'},
	{"mindist", required_argument, 0, 'd'},
	{"maxdist", required_argument, 0, 'D'},
	{"subreads", required_argument, 0, 'n'},
	{"minmatch", required_argument, 0, 'm'},
	{"minmatch2", required_argument, 0, 'p'},
	{"quality", required_argument, 0, 'Q'},
	{"threads", required_argument, 0, 'T'},
	{"index-threshold", required_argument, 0, 'f'},
	{"output", required_argument, 0, 'o'},
	{"order",  required_argument,0, 'S'},
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

void print_bed_table(HashTable * bed_table, char * out_fn)
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
			char * chro_name;
			unsigned int chro_pos, chro_pos_big;
			if (!cursor) break;
			paired_exon_key * p = (paired_exon_key * ) cursor -> key;
			int *counter = (int*) cursor ->value;

			locate_gene_position( p->small_key , &_global_offsets, &chro_name, &chro_pos);
			locate_gene_position( p->big_key , &_global_offsets, &chro_name, &chro_pos_big);
			fprintf(ofp,"%s\t%u\t%u\t%d\n", chro_name, chro_pos, chro_pos_big, *counter);
			cursor = cursor->next;
		}
	}
	fclose(ofp);
}

int main_junction(int argc,char ** argv)
{
	char read_file [300], read2_file [300];
	char output_file [300];

	char index_prefix [300];
	unsigned int all_reads, all_tables;
	halves_record_t halves_record;
	unsigned long long int processed_reads = 0, succeed_reads = 0;
	gene_input_t ginp, ginp2;
	HashTable * bed_index ;
	HashTable * pos_index ;
	//gene_flat_t my_flat ;
	//create_flat_strip(&my_flat);

	int c;
	int option_index = 0;

	TOTAL_SUBREADS = 20;
	ACCEPT_SUBREADS = 3;
	ACCEPT_MINOR_SUBREADS = 1;
	INDEX_THRESHOLD = 12;
	read_file[0]=0;
	read2_file[0]=0;
	index_prefix[0]=0;
	output_file[0]=0;
	all_reads = 14*1024*1024;

	printf("\n");



	while ((c = getopt_long (argc, argv, "bSd:D:n:m:p:f:R:r:i:o:T:Q:I:?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'S':
				EXON_FIRST_READ_REVERSE = optarg[0]=='r'?1:0;
				EXON_SECOND_READ_REVERSE = optarg[1]=='f'?0:1;
				break;
			case 'b':
				EXON_USE_VALUE_ARRAY_INDEX = 1;
				break;
			case 'D':
				EXON_MAX_PAIRED_DISTANCE = atoi(optarg);
				break;
			case 'd':
				EXON_MIN_PAIRED_DISTANCE = atoi(optarg);
				break;
			case 'n':
			//	TOTAL_SUBREADS = atoi(optarg);
				printf(" === WARNING ===\n You cannot set subread numbers in exon detection. It is automatically set to the max number.");
				break;
			case 'f':
				INDEX_THRESHOLD  = atoi(optarg);
				break;
			case 'm':
				ACCEPT_SUBREADS = atoi(optarg);
				break;
			case 'T':
				EXON_ALL_THREADS = atoi(optarg);
				break;
			case 'r':
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
				break ;
			case 'Q':
				if(optarg[0]=='l')
					EXON_QUALITY_SCALE = QUALITY_SCALE_LINEAR;
				if(optarg[0]=='e')
					EXON_QUALITY_SCALE = QUALITY_SCALE_LOG;
				if(optarg[0]=='n')
					EXON_QUALITY_SCALE = QUALITY_SCALE_NONE;
				break;
			case 'p':
				ACCEPT_MINOR_SUBREADS = atoi(optarg);
				break;
			case '?':
				return -1 ;
		}

	if (!read_file[0] || !index_prefix[0] || !output_file[0])
	{
		exon_usage(argv[0]);

		return -1 ;
	}

	if(!EXON_USE_VALUE_ARRAY_INDEX)
	{
		printf("Detecting junction reads must reference to the base-wise indel. Please enable the '-b' option while building the index.\n");
		return -1;
	}
  
	reads_density = guess_reads_density(read_file);
	if(reads_density<0)
		printf("Input file '%s' is not found or is in an incorrect format.\n", read_file);

	if(geinput_open(read_file, &ginp))
		return -1;

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

	FILE * out_fp = fopen(output_file, "w");
	if (!out_fp)
	{
		printf("Unable to open the output file at '%s'.\n", output_file);
		return -1;
	}


	printf("Number of subreads selected for each read=%d\n", TOTAL_SUBREADS);
	printf("Threshold on number of subreads for a successful mapping=%d\n", ACCEPT_SUBREADS);
	printf("Number of threads=%d\n", EXON_ALL_THREADS);
	printf("Tolerance for Indel=%d\n", EXON_INDEL_TOLERANCE);
	if (EXON_QUALITY_SCALE==QUALITY_SCALE_LINEAR)
		puts("Quality scale=linear\n\n");
	else if (EXON_QUALITY_SCALE==QUALITY_SCALE_LOG)
		puts("Quality scale=exponential\n\n");
	else 	puts("\n");

	if (read2_file[0])
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
		
	if(read2_file[0] && geinput_open(read2_file, &ginp2))
	{
		printf("Input file '%s' is not found or is in an incorrect format.\n", read2_file);
		return -1;
	}

	fflush(stdout);

	begin_ftime = miltime();

	load_offsets (&_global_offsets, index_prefix);

	bed_index = HashTableCreate(399997);
	pos_index = HashTableCreate(399997);

	HashTableSetKeyComparisonFunction(pos_index, pointercmp_forpos);
	HashTableSetHashFunction(pos_index, pointerHashFunction_forpos);

	HashTableSetKeyComparisonFunction(bed_index, pointercmp_forbed);
	HashTableSetHashFunction(bed_index, pointerHashFunction_forbed);


	while (1)
	{
		char inbuff[1201];

		int new_processed_reads = run_exon_search_index(&ginp, read2_file[0] ? (&ginp2):NULL, index_prefix, &halves_record, out_fp, processed_reads, all_tables, &succeed_reads, bed_index, pos_index);
		if(new_processed_reads<0)break;

		processed_reads += new_processed_reads;
		clear_halve_record(&halves_record);

		// test if there no anyreads remaining
		unsigned long long int current_fp = ftello(ginp.input_fp);
		int rl = geinput_next_read(&ginp, NULL, inbuff, NULL);
		if (rl<0)
			break;
		fseeko(ginp.input_fp, current_fp, SEEK_SET);
	}

	geinput_close(&ginp);
	print_bed_table(bed_index, output_file);
	HashTableDestroy(bed_index);
	HashTableDestroy(pos_index);


	if(IS_DEBUG)
		printf("@LOG THE END. \n");
	else
		printf("\n\n %llu reads were processed in %.1f seconds.\nPercentage of successfully mapped reads is %0.2f%%.\n\n", processed_reads, miltime()-begin_ftime, succeed_reads*100.0/processed_reads/(read2_file[0]?2:1));

	printf("\n\nCompleted successfully.\n");

	return 0;
}
