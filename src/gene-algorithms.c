#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/timeb.h>
#include "input-files.h"
#include "gene-algorithms.h"
#include "sorted-hashtable.h"

void non_func(const char * fmt, ...)
{
}


int MAX_CIGAR_LEN = 30;


//#define indel_debug

#ifdef indel_debug
#define ddprintf printf
#else
#define ddprintf non_func 
#endif

#define get_base_quality_score(quality_chr, quality_scale)  ((quality_scale)==QUALITY_SCALE_NONE?-0.01: ((quality_scale) == QUALITY_SCALE_LINEAR?(	(quality_chr) - '@')*0.01-0.03: correct_rate_table[(quality_chr) - '@']))


double begin_ftime;

unsigned char get_next_char(FILE * fp)
{
	int find_br = 0;
	while (!feof(fp))
	{
		char nch;
		nch = fgetc(fp);
		if (find_br)
		{
			if (nch == '\n')
				find_br=0;
		}
		else if (nch == '>')
			find_br = 1;
		else if (nch > 32)
			return nch;
	}
	return 0;
}

int is_valid_subread(const char * read_str)
{
	int i;
	for (i=0; i<16; i++)
		if(read_str[i] == 'N' || read_str[i] == '.')
			return 0;
	return 1;
}

double correct_rate_table [] ={ -1.58147375341 , -0.99684304401 , -0.69552447133 , -0.50767587370 , -0.38013040807 , -0.28926818720 , -0.22255151597 , -0.17255657291 , -0.13455196029 , -0.10536051566 , -0.08276530267 , -0.06517417320 , -0.05141827416 , -0.04062484422 , -0.03213357402 , -0.02543972753 , -0.02015436476 , -0.01597586925 , -0.01266917021 , -0.01005033585 , -0.00797499828 , -0.00632956293 , -0.00502447389 , -0.00398901727 , -0.00316728823 , -0.00251504651 , -0.00199725550 , -0.00158615046 , -0.00125971852 , -0.00100050033 , -0.00079464388 , -0.00063115648 , -0.00050131287 , -0.00039818644 , -0.00031627778 , -0.00025122020 , -0.00019954614 , -0.00015850188 , -0.00012590047 , -0.00010000500 , -0.00007943598 , -0.00006309773 , -0.00005011998 , -0.00003981151 , -0.00003162328 , -0.00002511918 , -0.00001995282 , -0.00001584906 , -0.00001258933 , -0.00001000005 , -0.00000794331 , -0.00000630959 , -0.00000501188 , -0.00000398108 , -0.00000316228 , -0.00000251189 , -0.00000199526 , -0.00000158489 , -0.00000125893 , -0.00000100000 , -0.00000079433 , -0.00000063096 , -0.00000050119 , -0.00000039811 , -0.00000031623 , -0.00000025119 , -0.00000019953 , -0.00000015849 , -0.00000012589 , -0.00000010000 , -0.00000007943 , -0.00000006310 , -0.00000005012 , -0.00000003981 , -0.00000003162 , -0.00000002512 , -0.00000001995 , -0.00000001585 , -0.00000001259 , -0.00000001000 , 0., 0.,0.,0.,0.,0.,0.,0.};


#define get_base_error_prob64(a) ((a) < '@'-1?1:pow(10., -0.1*((a)-'@')))
#define get_base_error_prob33(a) ((a) < '!'-1?1:pow(10., -0.1*((a)-'!'))) 

int reported_version_error = 0;
gene_quality_score_t get_subread_quality(const char * quality_str, const char * read_str, int quality_scale, int phred_version)
{
	gene_quality_score_t ret =0;
	int i;
/*
	for (i=0; i<16; i++)
		if(read_str[i] == 'N' || read_str[i] == '.')
			return -22.;
*/

	//for(i=0;i<16;i++) ret += get_base_quality_score(quality_str[i] , quality_scale);
	if(FASTQ_PHRED64 == phred_version)
		for(i=0;i<16;i++) ret += (1. - get_base_error_prob64(quality_str[i])); 
	else
		for(i=0;i<16;i++) ret += (1. - get_base_error_prob33(quality_str[i])); 

	if (ret <0 && !reported_version_error)
	{
		printf("\nWARNING: negative Phred quality score! Please verify the version of the Phred scores.\n");
		reported_version_error=1;
	}
	
	return ret/16;
}

/*int get_base_phred(char quality_chr)
{
	return quality_chr - '@';
}
*/

void print_running_log(double finished_rate, double read_per_second, double expected_seconds, unsigned long long int total_reads, int is_pair)
{
        char outbuff[99]; int i;
        snprintf(outbuff, 98,"completed=%0.2f%%; time used=%.1fs; rate=%.1fk reads/s; time left=%.1fs; total=%lluk %s", finished_rate*100, miltime()-begin_ftime, read_per_second/1000 ,expected_seconds, total_reads/1000, is_pair?"pairs":"reads");
        fputs(outbuff, stdout);
        for(i=strlen(outbuff); i<105; i++)
                printf(" ");
        printf("\r");
}

void print_text_scrolling_bar(char * hint, float percentage, int width, int * internal_counter)
{
	char fan = '-';
	int bar_width = width - 7 - strlen(hint) , i;
	int dash_width = (int)(bar_width * percentage+0.5);
	dash_width = min(dash_width, bar_width - 1);
	int space_width = bar_width - 1 - dash_width;

	for (i=0; i<width; i++);
		putchar(' ');
	printf("\r");

 
	switch ((*internal_counter) % 4)
	{
		case 0:
			fan='-';
			break;
		case 1:
			fan='\\';
			break;
		case 2:
			fan='|';
			break;
		case 3:
			fan='/';
			break;
	}

	(*internal_counter) ++;

	printf (" %c %s [", fan, hint);
	for(i = 0; i< dash_width; i++)
		putchar('=');
	putchar('>');
	for(i = 0; i< space_width; i++)
		putchar(' ');
	printf("]\r");
	
	fflush(stdout);
}


int remove_repeated_reads(gehash_t * table, gehash_t * huge_table, int index_threshold)
{
	int i;
	int vals[200000];
	int val_len[200000];
	int scroll_count = 0;
	unsigned int all_removed = 0;


	for (i=0; i <table->buckets_number ; i++)
	{
		struct gehash_bucket * cb = table->buckets + i;
		int j;
		int val_c = 0;

		if(i % 300 == 0)
			print_text_scrolling_bar ("Removing non-informative subreads", 1.0*i/table->buckets_number, 80, &scroll_count);

		for(j=0; j<cb->current_items; j++)
		{
			int k, found = 0;
			for(k=0;k<val_c;k++)
			{
				if (vals[k]==cb->item_keys[j])
				{
					val_len[k]++;
					found = 1;
					break;
				}
			}

			if(!found)
			{
				if(val_c>199999)
				{
					printf("\nThere are too many items in a bucket; you may reduce the threshold of non-informative subreads to eliminate this problem.\n");
					exit(-1);
				}
				else
				{
					vals[val_c] = cb->item_keys[j];
					val_len[val_c] = 1;
					val_c++;
				}
			}
		}

//			printf ("F3\n");

		for(j=0; j<val_c; j++)
		{

			if (gehash_exist(huge_table, vals[j]))
				 gehash_remove(table, vals[j]);
			else if(val_len[j]> index_threshold)
			{
				gehash_remove(table, vals[j]);
				gehash_insert(huge_table, vals[j], 1);
				all_removed += val_len[j];
			}
		}
	}

	if(IS_DEBUG)
		printf ("@LOG There are %u subreads removed from the index.\n", all_removed);

	return all_removed;
}

unsigned int get_gene_linear(int chrono, int offset, const unsigned int offsets [])
{
	if (chrono>1)return offsets[chrono-1]+offset;
	return offset;
}

int locate_gene_position(unsigned int linear, const gene_offset_t* offsets , char ** chro_name, unsigned int * pos)
{
	int n;

	for (n=0; offsets->read_offset[n]; n++)
	{
		if (offsets->read_offset[n] > linear)
		{
			if (n==0)
				*pos = linear;
			else
				*pos = linear - offsets->read_offset[n-1];

			*chro_name = (char *)offsets->read_name[n];

			return 0;
		}
	}
	return 1;
}

#define _index_vote(key) (((unsigned int)key)%GENE_VOTE_TABLE_SIZE)

int vv=0;
inline void add_gene_vote(gene_vote_t* vote, int key , int add_new)
{
	add_gene_vote_weighted(vote, key, add_new, 1);
}

inline void add_gene_vote_weighted(gene_vote_t* vote, int key , int add_new, int w)
{
	int offset = _index_vote(key);
	int datalen = vote -> items[offset];
	unsigned int * dat = vote -> pos[offset];
	int i;

	for (i=0;i<datalen;i++)
	{
		if (dat[i] == key)
		{

			int test_max = (vote->votes[offset][i]);
			test_max += w;
			vote->votes[offset][i] = test_max;

			if(test_max > vote->max_vote){
				vote->max_vote = test_max;
				vote->max_position = key;
			}
			return;
		}
	}

	if(!add_new || datalen >=GENE_VOTE_SPACE)
		return;

//	w += 16;

	vote -> items[offset] ++;
	dat[i] = key;
	vote->votes[offset][i]=w;
}

int max_gene_vote(gene_vote_t* vote, int * position_result, int query)
{
	int n, i, max_index=0;
	int max_votes = -1;
	for (n=0; n<GENE_VOTE_TABLE_SIZE; n++)
	{
		int itms = vote->items[n];
		gene_vote_number_t * vots = vote->votes[n];
		for (i=0; i<itms; i++)
		{

			if (vots[i] > max_votes)
			{
				max_index = n<<16|i;
				max_votes = vots[i];
			}
		}
	}

	if(max_votes == -1)
	{
		*position_result = -1;
		return 0;
	}
	else
	{
		//printf("err:%X\n",max_index);
		*position_result = vote->pos[max_index>>16][max_index&0xffff];
		return max_votes;
	}
}


void clear_allvote(gene_allvote_t* allvote)
{
	bzero(allvote -> max_votes,  allvote -> max_len* sizeof(*allvote -> max_votes));
	bzero(allvote -> masks,  allvote -> max_len* sizeof(*allvote -> masks));
}

void init_allvote(gene_allvote_t* allvote, int expected_len, int allowed_indels)
{
	allvote -> max_len = expected_len; 
	allvote -> max_positions = (unsigned int *) malloc(sizeof(int)*expected_len);
	allvote -> max_votes = (gene_vote_number_t *) calloc(sizeof(gene_vote_number_t), expected_len);
	allvote -> max_quality = (gene_quality_score_t *) calloc(sizeof(gene_quality_score_t), expected_len);
	allvote -> masks = (unsigned char *) calloc(1, expected_len);
#ifdef REPORT_ALL_THE_BEST
	allvote -> best_records = (gene_best_record_t *) malloc(sizeof(gene_best_record_t)* expected_len);
#endif
	allvote -> is_counterpart = (unsigned char *) malloc(expected_len);

	allvote -> max_indel_tolerance = allowed_indels;
	allvote	-> indel_recorder_length = max(3*allowed_indels+1 , MAX_CIGAR_LEN+2);

	if(allowed_indels)
		allvote -> max_indel_recorder = (char *)malloc( allvote -> indel_recorder_length *expected_len);
	else	allvote -> max_indel_recorder =  NULL;
}

void add_allvote_q(gene_allvote_t* allvote,int qid , int pos, gene_vote_number_t votes, gene_quality_score_t quality, int is_counterpart, char mask, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads)
{


	if((votes > allvote -> max_votes[qid] + 0.1) || (votes > allvote -> max_votes[qid] - 0.1 && quality >  allvote -> max_quality[qid]))
	{
#ifdef REPORT_ALL_THE_BEST
		if(votes > allvote -> max_votes[qid])
			allvote -> best_records [qid].best_len = 0;

		if(allvote -> best_records [qid].best_len < BEXT_RESULT_LIMIT)
		{
			int bestlen = allvote -> best_records [qid].best_len;
			allvote -> best_records [qid].offsets[bestlen] = pos;
			allvote -> best_records [qid].is_reverse[bestlen] = is_counterpart;
			allvote -> best_records [qid].best_len ++;
		}
#endif

		allvote -> is_counterpart[qid] = is_counterpart;
		allvote -> max_positions[qid] = pos;
		allvote -> max_votes[qid] = votes;
		allvote -> max_quality[qid] =  quality;
		allvote -> masks[qid] = mask;

		if(allvote -> max_indel_recorder)
		{
			//PNT111
			if (array_index)
				*(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length) = 0xff ;
			else
				memcpy(allvote -> max_indel_recorder + qid *allvote -> indel_recorder_length, max_indel_recorder, allvote -> max_indel_tolerance *3* sizeof(char));

//			printf("%s\n", read_txt);

			if (array_index && max_indel_recorder[3]){
				explain_indel(allvote, qid , pos, max_indel_recorder,  array_index, read_txt, read_len,  max_indel, total_subreads);
			}


			if(max_indel_recorder[3] && array_index && 0)
			{
				
				char indel_operations[1500];
				//printf ("\n%c v=%d ; r=%s\n", is_counterpart?'~':'@' , votes, read_txt);
				int k;
				for (k=0; max_indel_recorder[k]; k+=3);
				int moves = dynamic_align(read_txt, read_len, array_index, pos, max_indel, indel_operations, max_indel_recorder[k-1], 0, read_len);
				int xx;

				int ixx = 0;
				for (xx=0; xx<moves; xx++)
				{
					if (indel_operations[xx] == 2)
						printf (".");
					else
						printf("%c", gvindex_get(array_index, (ixx++) + pos));
				}
				printf ("\n");
	
				for (xx=0; xx<moves; xx++)
				{
					if (indel_operations[xx]==0) printf("|");
					else if (indel_operations[xx]==1) printf ("<");
					else if (indel_operations[xx]==2) printf (">");
					else if (indel_operations[xx]==3) printf (".");
				}
				
				printf("\n");
				ixx = 0;
				for (xx = 0; xx<moves;xx++)
				{
					if(indel_operations[xx] == 1)
						printf (".");
					else if(read_txt[ixx])
						printf("%c",read_txt[ixx++]);
					else    printf (".");
				}
				printf("\n");

			}
		}

		//printf("%d\t%.11f\n", votes , good_votes);
	}
}

void explain_indel(gene_allvote_t* allvote, int qid , int pos, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads)
{

	char indel_operations[1500];
	int i,xx;
	char tmp_cigar [MAX_CIGAR_LEN+1];

	int current_pos = 0; 
	int explain_cursor = 0;
	char last_operation = 0;
	int del_number = 0;
	int dynamic_delta = 0;

	tmp_cigar[0]=0;

	for (i=3; max_indel_recorder[i]; i+=3)
	{
		int last_dist = max_indel_recorder[i-1];
		int black_subread_start = max_indel_recorder[i-2]-1;
		int black_subread_end = max_indel_recorder[i]-1;
		if (black_subread_end < black_subread_start+1) black_subread_end = black_subread_start+1;
		while(max_indel_recorder[i+3])
		{
			if (max_indel_recorder[i+3]- black_subread_end>2) break;
			black_subread_end = max_indel_recorder[i+3]-1;
			i+=3;
		}
		int black_base_start = find_subread_end(read_len, total_subreads,black_subread_start )- 3 + max(0,-last_dist);
		int black_base_end = find_subread_end(read_len, total_subreads, black_subread_end)-15 + 3 + max(0,-max_indel_recorder[i+2]);

		int blackref_base_start = black_base_start - max(0,-last_dist) + max(0,last_dist);
//		printf("baseref_offset=%d\n", max(0,last_dist));

		int gap_end_read =  black_base_end + max(0,max_indel_recorder[i+2]);

		blackref_base_start = max(0,blackref_base_start);
		black_base_start = max(0,black_base_start);
		black_base_end = min(read_len, black_base_end);

		int vpos = strlen(tmp_cigar);
		if (vpos > MAX_CIGAR_LEN -6) break;

		int exp_indel = last_dist - max_indel_recorder[i+2];

	//	int moves = window_indel_align (read_txt + black_base_start, black_base_end-black_base_start , array_index, pos + blackref_base_start, max_indel, indel_operations, exp_indel,  -10, black_base_end - black_base_start+5);
		int moves = dynamic_align (read_txt + black_base_start, black_base_end-black_base_start , array_index, pos + blackref_base_start, max_indel, indel_operations, exp_indel,  -10, black_base_end - black_base_start+5);

#ifdef indel_debug
		char tt = *(read_txt + black_base_end);
		read_txt [black_base_end] = 0;
		ddprintf ("%s\n", read_txt + black_base_start);
		read_txt [black_base_end] = tt;
		ddprintf("\n TESTING SUBR %d to %d BASE %d to %d CHRO %d to %d exp offset %d\n",black_subread_start, black_subread_end, black_base_start, black_base_end, blackref_base_start, -max(0,exp_indel)+blackref_base_start +(black_base_end-black_base_start), exp_indel);
		ddprintf(" moves %d\n", moves);
#endif

		last_operation = 0;
		if(moves)// < (black_base_end-black_base_start) + 6 + 1.5 * exp_indel)
		{
			int pre_read_len = 0;
			current_pos = black_base_start ;
			xx=0;
			while (pre_read_len < 0 && xx < moves)
			{
				if(indel_operations[xx]!=2)pre_read_len++;
				xx++;
			}
			for (; xx<moves; xx++)
			{
				int current_operation = indel_operations[xx];
				if (current_operation == 3) current_operation = 0;

				if(current_operation != last_operation && ( current_pos < gap_end_read-1 || current_operation == 0))
				{
					int vpos = strlen(tmp_cigar);
					if (vpos>MAX_CIGAR_LEN-6) break;
					sprintf(tmp_cigar + vpos, "%d%c", last_operation==1?del_number:(current_pos - explain_cursor), last_operation==0?'M':(last_operation==1?'D':'I'));
					explain_cursor = current_pos ;
					del_number = 0;
				}
				last_operation = current_operation; 

				if (indel_operations[xx]==0) current_pos ++;
				else if (indel_operations[xx]==1) {del_number++; dynamic_delta++;}//"D"
				else if (indel_operations[xx]==2) {current_pos ++; dynamic_delta--;} //"I"
				else if (indel_operations[xx]==3) current_pos ++;

				if (current_pos >= gap_end_read)break;

/*
				if (indel_operations[xx]==1 && exp_indel <0 || indel_operations[xx]==2 && exp_indel >0)
				{
					int vpos = strlen(tmp_cigar);
					if (vpos>MAX_CIGAR_LEN-13) break;
					sprintf(tmp_cigar + vpos, "%dM", current_pos-explain_cursor);
					vpos = strlen(tmp_cigar);
					if (vpos>MAX_CIGAR_LEN-10) break;

					if (indel_operations[xx]==1)
					{
						sprintf(tmp_cigar + vpos, "%dD", -exp_indel);
						explain_cursor = current_pos;
					}
					else
					{
						sprintf(tmp_cigar + vpos, "%dI", exp_indel);
						explain_cursor = current_pos+exp_indel	;
					}
					del_number = 0;
					break;
				}
				else if (indel_operations[xx]==0) current_pos ++;
				else if (indel_operations[xx]==1) {del_number++; dynamic_delta++;}//"D"
				else if (indel_operations[xx]==2) {current_pos ++; dynamic_delta--;} //"I"
				else if (indel_operations[xx]==3) current_pos ++;
*/
			} 
		}
		else
		{
			int movement = last_dist - max_indel_recorder[i+2];
			int vpos = strlen(tmp_cigar);
			if (vpos>MAX_CIGAR_LEN-6) break;

			current_pos = find_subread_end(read_len, total_subreads,max_indel_recorder[i]-1)-15 + (max_indel_recorder[i+2] >0?max_indel_recorder[i+2]:-max_indel_recorder[i+2]);
			current_pos -= (black_base_end - black_base_start)/2 -3;

			sprintf(tmp_cigar + vpos, "%dM%d%c", current_pos - explain_cursor , abs(movement), movement<0?'D':'I' );

			explain_cursor = current_pos + (movement >0?movement:0);
			current_pos = find_subread_end(read_len, total_subreads,max_indel_recorder[i+1]-1);
			last_operation = 0;
			del_number = 0;
			moves = 0;
		}
		last_operation = 0;
		current_pos = black_base_end;
//		if (!moves)
//			strcat(tmp_cigar,"X");
		//	printf("Wrong explaination! dist=%d   dyn=%d\n", max_indel_recorder[i+2], dynamic_delta);

	}

	int vpos = strlen(tmp_cigar);
	if (vpos > MAX_CIGAR_LEN-10)
	{
		*(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length) = 0xfd; 
		memcpy(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length+1, max_indel_recorder, 3 * max_indel);
	}
	else
	{
		if (explain_cursor<read_len)	
		{
			sprintf(tmp_cigar+vpos,"%dM", read_len - explain_cursor);
		}

	//printf("\n%s\n", tmp_cigar);
		*(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length) = 0xfe ;
		strncpy(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length+1, tmp_cigar, allvote -> indel_recorder_length - 2);
	}
}

int evaluate_piece(char * piece_str, int chron, int offset, int is_counterpart, int start_pos, int end_pos)
{
	char fname[300];
	int inner_pos = 0, i;
	FILE * fp;
	char next_char=0;
	int ret = 0;

	if (chron == 0)
		sprintf(fname, "/opt/Work2001/Gene-Search/src/GENE-LIB/%02da.fa", chron);
	else
		sprintf(fname, "/opt/Work2001/Gene-Search/src/GENE-LIB/%02d.fa", chron);

	inner_pos = offset + offset / 70;

	fp = fopen(fname,"r");

	while(next_char!='\n')
		next_char=fgetc(fp);
	fseek(fp, inner_pos, SEEK_CUR);

	for(i=0 ; i<end_pos; i++)
	{
		next_char = get_next_char(fp);
		if(i < start_pos)
			continue;
		if(next_char == 'N')
			ret ++;
		else{
			if(is_counterpart)
			{
				if (piece_str[99-i] == 'A' && next_char == 'T')
					ret ++;
				else if (piece_str[99-i] == 'G' && next_char == 'C')
                                        ret ++;
				else if (piece_str[99-i] == 'T' && next_char == 'A')
                                        ret ++;
				else if (piece_str[99-i] == 'C' && next_char == 'G')
                                        ret ++;
			}
			else if(piece_str[i] == next_char)
				ret ++;
		}
	}

	fclose(fp);

	return ret;

}

void mark_votes_array_index(char * read_str, int read_len, gene_vote_t * dest, gene_vote_t * src, gene_value_index_t * my_array_index, int color_space, int indel_tolerance, int min_minor, const char quality_str [], int quality_scale)
{
	int i, j;

	dest -> max_vote = 0;
	dest -> max_quality = 0;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
	{
		dest -> items[i] = src -> items[i];
		for(j=0; j< src->items[i]; j++)
		{
			unsigned int potential_position = src-> pos[i][j];
			float matchingness_count = 0;
			if (src -> votes [i][j]>= min_minor)
				matchingness_count = match_read(read_str, read_len, potential_position, my_array_index, color_space, indel_tolerance, quality_str, quality_scale);
	
			dest -> pos[i][j] = potential_position;
			dest -> quality[i][j] = matchingness_count;
			dest -> votes [i][j] =src -> votes [i][j];
			dest -> masks[i][j] = src -> masks[i][j];
			// find for the `best positions'; only the best positions are replaced by the base machingness scores.
			if((matchingness_count >  dest -> max_quality && src -> votes [i][j] == dest -> max_vote) || (src -> votes [i][j] > dest -> max_vote))
			{
				dest -> max_vote  = src -> votes [i][j];
				dest -> max_mask  = src -> masks[i][j];
				dest -> max_quality = matchingness_count;
				dest -> max_position = potential_position;
			}
		}
	}
}

int select_positions_array(char * read1_str, int read1_len, char * read2_str, int read2_len,  gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2,  gene_quality_score_t * sum_quality, gene_quality_score_t * qual_r1, gene_quality_score_t * qual_r2, gehash_data_t * pos_read1, gehash_data_t * pos_read2, unsigned int max_pair_dest, unsigned int min_pair_dest, int min_major, int min_minor,int is_negative_strand, gene_value_index_t * my_array_index, int color_space, int indel_tolerance, int number_of_anchors_quality, const char quality_str1 [], const char quality_str2 [], int quality_scale, int max_indel_len)
{
	gene_vote_t base_vote_1 , base_vote_2;
	char r1_recorder [MAX_INDEL_TOLERANCE*3], r2_recorder[MAX_INDEL_TOLERANCE*3];

	mark_votes_array_index(read1_str, read1_len, &base_vote_1, vote_read1, my_array_index, color_space, indel_tolerance, min_minor, quality_str1, quality_scale);
	mark_votes_array_index(read2_str, read2_len, &base_vote_2, vote_read2, my_array_index, color_space, indel_tolerance, min_minor, quality_str2, quality_scale);

	return select_positions(&base_vote_1, &base_vote_2, numvote_read1, numvote_read2, sum_quality, qual_r1, qual_r2,  pos_read1, pos_read2, r1_recorder , r2_recorder ,max_pair_dest,  min_pair_dest,  min_major,  min_minor, is_negative_strand, number_of_anchors_quality,max_indel_len);
}


int select_positions(gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2, gene_quality_score_t * sum_quality, gene_quality_score_t * qual_r1, gene_quality_score_t * qual_r2 , gehash_data_t * pos_read1, gehash_data_t * pos_read2, char * read1_indel_recorder, char * read2_indel_recorder, unsigned int max_pair_dest, unsigned int min_pair_dest, int min_major, int min_minor,int is_negative_strand, int number_of_anchors_quality, int max_indel_len)
{

	int k, i, j, anchors = 0;
	gehash_data_t anchors_position [ANCHORS_NUMBER];
	unsigned char anchor_read [ANCHORS_NUMBER];

	gene_vote_number_t anchor_votes = max(vote_read1->max_vote, vote_read2->max_vote);
	gene_vote_number_t anchors_votes [ANCHORS_NUMBER];
	gene_quality_score_t anchors_quality [ANCHORS_NUMBER];
	int anchors_buckets [ANCHORS_NUMBER];
	int anchors_index [ANCHORS_NUMBER];

	gene_vote_number_t minor_votes [ANCHORS_NUMBER];
	gene_quality_score_t minor_quality [ANCHORS_NUMBER];
	gehash_data_t minor_position [ANCHORS_NUMBER];
	int minor_buckets [ANCHORS_NUMBER];
	int minor_index [ANCHORS_NUMBER];

	if(anchor_votes < min_major)
		return 0;

	for (k=0; k<2; k++)
	{
		gene_vote_t * current_vote = k?vote_read2:vote_read1;
		if (current_vote->max_vote < anchor_votes)
			continue;

		for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
			for(j=0; j< current_vote->items[i]; j++)
			{
				if(current_vote->votes[i][j] == anchor_votes)
				{
					anchors_position[anchors] = current_vote->pos[i][j];
					anchors_votes   [anchors] = current_vote->votes[i][j];
					anchors_quality [anchors] = current_vote->quality[i][j];
					anchor_read	[anchors] = k;
					anchors_buckets [anchors] = i;
					anchors_index   [anchors] = j;

					anchors++;

					if(anchors >= ANCHORS_NUMBER)
					{
	//					printf("Abnormally too many anchors [%d, %d] %f\n", i, j, anchor_votes);
						anchors --;
					}
				}
			}

	}

	bzero(minor_votes, ANCHORS_NUMBER*sizeof(gene_vote_number_t));
	bzero(minor_quality, ANCHORS_NUMBER*sizeof(gene_quality_score_t));

	for (k=0; k<2; k++)
	{
		gene_vote_t * current_vote = k?vote_read2:vote_read1;
		gene_vote_t * current_vote2 = k?vote_read1:vote_read2;

		if (current_vote2->max_vote < anchor_votes)
			continue;

		for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
			for(j=0; j< current_vote->items[i]; j++)
			{
				if(current_vote->votes[i][j] >= min_minor)
				{
					int l;
					for (l=0; l<anchors; l++)
					{

						if(anchor_read[l]==k)continue;

						long long int abdist = current_vote->pos[i][j];
						abdist -= anchors_position[l];

						if((anchor_read[l] && !is_negative_strand) ||
						    ((!anchor_read[l]) && is_negative_strand))
							abdist = -abdist;
						//printf ("\n\nDIST=%d, ANCHOR=%d, NEGATIVE=%d\n", abdist, anchor_read[l], is_negative_strand);

						if (abdist < 0)
						//#warning WGSIM generates paired-end data in a wrong way; abdist is given its absolute value in such a case.
							abdist = -abdist;
						//	continue;

						if (	(minor_votes[l] < current_vote->votes[i][j] ||
							(minor_votes[l] == current_vote->votes[i][j] && current_vote->quality [i][j] > minor_quality[l]) ) && 
							(abdist <= max_pair_dest) &&
							(abdist >= min_pair_dest) 
						   )
						{
							minor_votes[l] = current_vote->votes[i][j];
							minor_position[l] = current_vote->pos[i][j];
							minor_quality[l] = current_vote->quality [i][j];
							minor_buckets[l] = i;
							minor_index[l] = j;
						}
					}
				}
			}

	}

	gene_vote_number_t selection_minor_vote = 0;
	gene_vote_number_t selection_major_vote = 0;
	gehash_data_t selection_minor_pos = 0;
	gehash_data_t selection_major_pos = 0;

	int selected_major_index=0, selected_major_bucket=0, selected_minor_bucket=0, selected_minor_index=0;

	*sum_quality = 0;

	unsigned char selection_major_read_no = 0;

	for (k=0; k<anchors; k++)
	{
		if(minor_votes[k]==0)
			continue;

		if ((minor_votes[k]+anchors_votes[k]) > (selection_minor_vote+selection_major_vote) ||
	   	       ((minor_votes[k]+anchors_votes[k]) == (selection_minor_vote+selection_major_vote) &&
		 	 minor_quality[k] + anchors_quality[k] > *sum_quality))
		{
			selection_minor_pos = minor_position[k];
			selection_major_pos = anchors_position[k];
			selection_minor_vote = minor_votes[k];
			selection_major_vote = anchors_votes[k];
			selection_major_read_no = anchor_read[k];
			selected_major_bucket = anchors_buckets[k];
			selected_minor_bucket = minor_buckets[k];
			selected_major_index = anchors_index[k];
			selected_minor_index = minor_index[k];
			* qual_r1 = anchor_read[k] ? minor_quality[k]:anchors_quality[k];
			* qual_r2 = anchor_read[k] ? anchors_quality[k]: minor_quality[k];
			* sum_quality = minor_quality[k] + anchors_quality[k];
		}
	}

	if(selection_minor_vote > 0)
	{
		if(selection_major_read_no)	// major is on read 2
		{
			*numvote_read1 = selection_minor_vote;
			*numvote_read2 = selection_major_vote;
			*pos_read1 = selection_minor_pos;
			*pos_read2 = selection_major_pos ;
			memcpy(read1_indel_recorder, vote_read1 -> indel_recorder[selected_minor_bucket][selected_minor_index], max_indel_len*3*sizeof(char));
			memcpy(read2_indel_recorder, vote_read2 -> indel_recorder[selected_major_bucket][selected_major_index], max_indel_len*3*sizeof(char));
		}
		else
		{
			*numvote_read1 = selection_major_vote;
			*numvote_read2 = selection_minor_vote;
			*pos_read1 = selection_major_pos;
			*pos_read2 = selection_minor_pos;
			memcpy(read1_indel_recorder, vote_read1 -> indel_recorder[selected_major_bucket][selected_major_index], max_indel_len*3*sizeof(char));
			memcpy(read2_indel_recorder, vote_read2 -> indel_recorder[selected_minor_bucket][selected_minor_index], max_indel_len*3*sizeof(char));
		}
		return 1;
	}

	return 0;
}

int load_offsets(gene_offset_t* offsets , const char index_prefix [])
{
	char fn[300];
	FILE * fp;
	int n=0;

	sprintf(fn, "%s.reads", index_prefix);

	fp = fopen(fn, "r");

	if(!fp)
	{
		printf("file not found :%s\n", fn);
		return 1;
	}

	while (!feof(fp))
	{
		int i=0, step = 0, j=0;

		read_line(fp, fn, 0);
		if (strlen(fn)<2)continue;
		while (fn[i])
		{
			if (fn[i] == '\t')
			{
				fn[i]=0;
				offsets->read_offset[n] = (unsigned int)atoll(fn);
				step = 1;
			}
			else if (step)
			{
				if(j<47)
				{
					offsets->read_name[n][j++] = fn[i];
					offsets->read_name[n][j]=0;
				}
			}
			i++;
		}
		n++;

		offsets->read_offset[n] = 0;
	}

	fclose(fp);
	return 0;
}

double miltime(){
        struct timeb trp;
        double ret;

        ftime(&trp);

        ret = trp.time*1.0+(trp.millitm*1.0/1000.0);
        return ret;
}

#define front2(str, bias)	(*((str)+(bias))+*((str)+1+(bias)))
#define front4(str, bias)	(front2(str, bias)+front2(str, bias+2))
#define front8(str, bias)	(front4(str, bias)+front4(str, bias+4))
#define front16(str, bias)	(front8(str, bias)+front8(str, bias+8))

float match_read(const char read_str[], int read_len, unsigned int potential_position,  gene_value_index_t * my_array_index, int space_type, int indel_tolerance, const char quality_str [], int quality_scale)
{
	float ret = 0;
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

/*			if(i % 3 == 2)
			{
				int b2_match = (read_matchingness[bias+indel_tolerance][i -2]  && read_matchingness[bias+indel_tolerance][i -1] && read_matchingness[bias+indel_tolerance][i]);
				if(b2_match)
				{
					read_matchingness[bias+indel_tolerance][i]  = 2;
					read_matchingness[bias+indel_tolerance][i-1]  = 2;
					read_matchingness[bias+indel_tolerance][i-2]  = 2;
				}
			}*/
		}
	}

	for(i=0; i<read_len-4; i+=4)
	{
	
		int j;
		int best_movement = 0; 
		float max_matchness = -1;
		for (j=-indel_tolerance; j<=indel_tolerance; j++)
		{
			int m = front4(read_matchingness[j],i);
			if(m > max_matchness)
			{
				max_matchness = m*1.;
				best_movement = j;
			}
		}

//		printf("best-match=%f\n", max_matchness);
//

		if (quality_str[0])
		{
			max_matchness = 0;
		//	float qs0 = max_matchness;
			for (j=0; j<4; j++)
			{
				if(read_matchingness[best_movement][j+i])
				{
					// penalty
					max_matchness += 1.03+get_base_quality_score(quality_str[j+i], quality_scale);
				}
				else
				{
					// relief
					//max_matchness += 0.2-get_base_phred(quality_str[j+i])*0.01;
				}
			}
		//	printf("Delta matchingness = %f\n", qs0 - max_matchness);
		}
		
		ret += max_matchness; //front8(read_matchingness[movement], i);  // (front4(read_matchingness[movement], i)==4?4:0);
	}

	//printf ("Read matchingness = %f\n", (ret));

	return  ret;
}


void final_matchingness_scoring(const char read_str[], const char quality_str[], int read_len, gene_vote_t * vote, gehash_data_t * max_position, gene_vote_number_t * max_vote, char *max_mask, gene_quality_score_t * max_quality, gene_value_index_t * my_array_index, int space_type, int indel_tolerance, int quality_scale)
{
	int i, j;
	gene_quality_score_t max_matching = -1.0;

	*max_vote = vote -> max_vote; //vote -> votes[i][j] + matchingness_count * 3;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			unsigned int potential_position = vote -> pos[i][j];
			gene_quality_score_t matchingness_count = match_read(read_str, read_len, potential_position, my_array_index, space_type, indel_tolerance, quality_str, quality_scale);
			if(vote->votes[i][j] == vote -> max_vote && matchingness_count > max_matching)
			{
				max_matching = matchingness_count+1E-5;
		//		printf ("MAX MATCHING = %f / %s\n", max_matching, read_str);
				*max_position = potential_position;
				*max_mask = vote -> masks[i][j];
			}
		}
}

#define MAX_READ_LENGTH 1208
int DPALIGN_CREATEGAP_PENALTY = 0 ;
int DPALIGN_EXTENDGAP_PENALTY = 0 ;
int DPALIGN_MISMATCH_PENALTY = 0;
int DPALIGN_MATCH_SCORE = 2;


int backup_dynamic_align(char * read, int read_len, gene_value_index_t * index, unsigned int begin_position, int max_indel, char * movement_buffer)
// read must be converted to the positive strand.
// movement buffer: 0:match, 1: read-insert, 2: gene-insert, 3:mismatch
// the size of the movement buffer must be equal to the length of the read plus max_indel * 3.
{
//	max_indel = 30; 
	short table[MAX_READ_LENGTH][MAX_READ_LENGTH];
	int i,j;
	for (i=0; i<read_len; i++)
		for(j=0; j<read_len + max_indel; j++)
		{

			if (j < i - max_indel || j > max_indel + i)
			{
				table[i][j]=-9999;
				continue;
			}

			short from_upper;

			if (i>0) from_upper = table[i-1][j] + DPALIGN_CREATEGAP_PENALTY;
			else     from_upper = DPALIGN_CREATEGAP_PENALTY;

			short from_left;

			if (j>0) from_left = table[i][j-1] + DPALIGN_CREATEGAP_PENALTY;
			else     from_left = DPALIGN_CREATEGAP_PENALTY;


			char is_matched_ij = gvindex_get(index, begin_position + j) == read[i]?DPALIGN_MATCH_SCORE :DPALIGN_MISMATCH_PENALTY;
			short from_upperleft;
			if (i>0 && j>0) from_upperleft = table[i-1][j-1] + is_matched_ij;
			else            from_upperleft = is_matched_ij; 

			table[i][j]=max(from_upper, max(from_left, from_upperleft));
		}

	short current_score = table[read_len-1] [read_len + max_indel-1];
	short path_i = read_len -1;
	int out_pos = 1499;
	char out_tmp [1500];
	j = read_len + max_indel-1;
	int last_move =-1;
	while (1)
	{
	//	printf("(%d, %d) = %d [%c]\n", path_i,j, current_score, read[path_i]);

		short upper_score;
		if (path_i>0) upper_score = table[path_i-1][j];
		else	      upper_score = 0;

		short left_score;
		if (j>0)      left_score = table[path_i][j-1];
		else	      left_score = 0;

		short upperleft_score;
		if (j>0 && path_i>0) upperleft_score = table[path_i-1][j-1];
		else	      upperleft_score = 0;

                char is_matched_ij = DPALIGN_MATCH_SCORE * (gvindex_get(index, begin_position + j) == read[path_i]?2:-1);

	//	printf("LS=%d, US=%d, LUS=%d, CURS=%d, MATCH=%d\n", left_score, upper_score, upperleft_score,  current_score , is_matched_ij);

		char finished=0;

		if (last_move ==0 || last_move == 3)
		{
	                if (is_matched_ij + (upperleft_score - current_score) == 0)
	                {
        	                path_i --;
                	        j --;
                        	out_tmp[out_pos] = is_matched_ij==2?0:3;
				current_score = upperleft_score;
				finished = 1;
                	}
		}
		else if (last_move == 1)
		{

			if (left_score + DPALIGN_CREATEGAP_PENALTY == current_score) 
			{
				out_tmp[out_pos] = 1;
				current_score = left_score;
				j --;
				finished = 1;
			}
		}
		else if (last_move == 2)
		{
			if (upper_score + DPALIGN_CREATEGAP_PENALTY == current_score)
			{
				out_tmp[out_pos] = 2;
				path_i --;
				current_score = upper_score;
				finished = 1;
			}
		}

		if (!finished){
			if (is_matched_ij + (upperleft_score - current_score) == 0)
			{
				path_i --;
				j --;
				out_tmp[out_pos] = is_matched_ij==2?0:3;
				current_score = upperleft_score;
				last_move = 0;
			}
			else if (left_score + DPALIGN_CREATEGAP_PENALTY == current_score) 
			{
				out_tmp[out_pos] = 1;
				current_score = left_score;
				j --;
				last_move = 1;
			}
			else if (upper_score + DPALIGN_CREATEGAP_PENALTY == current_score)
			{
				out_tmp[out_pos] = 2;
				path_i --;
				current_score = upper_score;
				last_move = 2;
			}
 			else  printf("ERROR!\n");
		}

		out_pos --;

		if(j < 0 || path_i < 0)
		{
			for (i=0; i<=j; i++) out_tmp[out_pos--] = 1;
			for (i=0; i<=path_i; i++) out_tmp[out_pos--] = 2;
			break;
		}
	}
	memcpy(movement_buffer, out_tmp  + out_pos +1, 1499 - out_pos);
	return 1499 - out_pos;
}


#define MAX_READ_LENGTH 1208 


int search_DP_branch(char * read, int read_len, gene_value_index_t * index, unsigned int begin_position, int path_i, int path_j, short table[MAX_READ_LENGTH][MAX_READ_LENGTH],char table_mask[MAX_READ_LENGTH][MAX_READ_LENGTH], int max_indel, char * movement_buffer, int expected_offset, int current_score, int out_pos, int current_offset, int init_read_offset, int shutdown_read_offset, int * all_steps)
{

	if (1499 - out_pos > (read_len << 2) || (*all_steps) > 3000 + (read_len << 5) || (*all_steps) > 57000){
		ddprintf("\nTOO MANY STEPS: rl = %d    len=%d    steps=%d\n", read_len, 1499 - out_pos, *all_steps);
//		if ((*all_steps) > 199999)
//			exit(0);
		return 0;
	}
	if(path_j < 0 || path_i < 0)
	{
		int i;
		ddprintf("\nFINAL TEST exp %d  real-offset %d\n", expected_offset, current_offset);
		if (expected_offset !=current_offset)return 0;
		for (i=0; i<=path_j; i++) movement_buffer[out_pos--] = 1;
		for (i=0; i<=path_i; i++) movement_buffer[out_pos--] = 2;
		return  out_pos ;
	}
	
	short upper_score;
	if (path_i>0) upper_score = table[path_i-1][path_j];
	else	      upper_score = 0;

	short left_score;
	if (path_j>0)      left_score = table[path_i][path_j-1];
	else	      left_score = 0;

	short upperleft_score;
	if (path_j>0 && path_i>0) upperleft_score = table[path_i-1][path_j-1];
	else	      upperleft_score = 0;

	char is_matched_ij =  gvindex_get(index, begin_position + path_j) == read[path_i]?DPALIGN_MATCH_SCORE :DPALIGN_MISMATCH_PENALTY;

	int found = 0;
	ddprintf ("CS = %d ; ", current_score);
	if (is_matched_ij + (upperleft_score - current_score) == 0)
	{
		movement_buffer[out_pos] = is_matched_ij==2?0:3;
		ddprintf ("Path_i = %d ; Path_j = %d ; Offset = %d => 0\n", path_i, path_j, current_offset);
		(*all_steps) ++;
		found = search_DP_branch (read, read_len, index, begin_position, path_i -1 , path_j -1, table , table_mask, max_indel, movement_buffer, expected_offset, upperleft_score, out_pos -1, current_offset, init_read_offset, shutdown_read_offset, all_steps); 
	}
	if ((!found ) && (left_score + (table_mask[path_i][path_j-1]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) == current_score || (current_score == 0 && left_score + (table_mask[path_i][path_j-1]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) < 0)))
	{
		movement_buffer[out_pos] = 1;
		(*all_steps) ++;
		ddprintf ("Path_i = %d ; Path_j = %d ; Offset = %d => 1\n", path_i, path_j, current_offset);
		found =search_DP_branch (read, read_len, index, begin_position, path_i , path_j -1, table  , table_mask, max_indel, movement_buffer, expected_offset, left_score, out_pos -1, current_offset - ((path_i >= init_read_offset && path_i <= shutdown_read_offset)?1:0),  init_read_offset, shutdown_read_offset, all_steps); 
	}
	if ((!found ) && (upper_score + (table_mask[path_i-1][path_j]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) == current_score || (current_score == 0 && upper_score + (table_mask[path_i-1][path_j]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) < 0 )))
	{
		movement_buffer[out_pos] = 2;
		(*all_steps) ++;
		ddprintf ("Path_i = %d ; Path_j = %d ; Offset = %d => 2\n", path_i, path_j, current_offset);
		found =search_DP_branch (read, read_len, index, begin_position, path_i -1 , path_j , table , table_mask , max_indel, movement_buffer, expected_offset, upper_score, out_pos -1, current_offset + ((path_i >= init_read_offset && path_i <= shutdown_read_offset)?1:0),  init_read_offset, shutdown_read_offset, all_steps); 
	}

	return found;
}

#define INDEL_WINDOW_WIDTH 4

int window_indel_align(char * read, int read_len, gene_value_index_t * index, unsigned int begin_position, int max_indel, char * movement_buffer, int expected_offset, int init_read_offset, int shutdown_read_offset)
{
	short indel_windows[32];
	int scores [32][MAX_READ_LENGTH];
	int i,j,windows_number;
	char chro_str[200];

	memset(indel_windows,0, 32*sizeof(short));

	windows_number = abs(expected_offset)+1;
	ddprintf ("\nWindow size = %d ; ExpOffset = %d\n", windows_number , expected_offset);
	for (i=0; i< read_len; i++)
		for (j=0; j< windows_number; j++)
		{
			int matchingness ;
			if (j==0) chro_str[i] = gvindex_get(index, begin_position+i);
			if (expected_offset<0)
				matchingness = *(read + i) == gvindex_get(index, begin_position + j + i) ;
			else
				matchingness = *(read + i) == gvindex_get(index, begin_position - windows_number + j + 1 + i) ;
			indel_windows [j] += matchingness;

			if (i >= INDEL_WINDOW_WIDTH)// The result for i - INDEL_WINDOW_WIDTH is OK
			{
				scores[j][i - INDEL_WINDOW_WIDTH] = indel_windows [j] ;
			}
			if (i >= INDEL_WINDOW_WIDTH)// The result for i - INDEL_WINDOW_WIDTH is OK
			{

				if (expected_offset<0)
					matchingness = *(read - INDEL_WINDOW_WIDTH + i) == gvindex_get(index, begin_position + j + i - INDEL_WINDOW_WIDTH) ;
				else
					matchingness = *(read - INDEL_WINDOW_WIDTH + i) == gvindex_get(index, begin_position -windows_number + j +1 + i - INDEL_WINDOW_WIDTH) ;
				indel_windows [j] -= matchingness;
			}
		}
	chro_str[i]=0;
	j = read[read_len];
	read[read_len]=0;

	ddprintf ("CHRO=%s\nREAD=%s\n", chro_str, read);

	/*
	float contingency_list [MAX_READ_LENGTH];
	for (i=0; i< read_len - windows_number-INDEL_WINDOW_WIDTH; i++)
	{
		float contingency;
		if (expected_offset >0)
			contingency  = ((scores[0][i]*1. + 1)-(scores[windows_number-1][i]*1.+1));
		else
			contingency  = ((scores[windows_number-1][i]*1. + 1)-(scores[0][i]*1.+1));
		contingency_list[i]=contingency;
	}*/

	int max_score=-1;
	int max_pos = -1;
	if (expected_offset>0)
	{
		//insertion
		for (i=read_len-INDEL_WINDOW_WIDTH-1; i>=0; i--)
		{
			if (scores[windows_number - 1][i - expected_offset]>=2 && scores[0][i] >= max_score)
			{
				max_pos = i - windows_number+1;
				max_score = scores[0][i];
			}
		}
	}
	else
	{
		//deletion
		for (i=read_len-INDEL_WINDOW_WIDTH-1; i>=0; i--)
		{
			if (scores[windows_number -1][i] >= max_score && scores[0][i + expected_offset] >=2)
			{
				max_score = scores[windows_number -1][i];
				max_pos = i ;
			}
		}
	
	}

	max_pos = min(read_len, max(0, max_pos));

	int move_number = 0;

	for (i=0; i< read_len -INDEL_WINDOW_WIDTH; i++)
	{
		if (i==max_pos)
			for (j = 0; j < windows_number-1; j++)
				movement_buffer[move_number++] = 1+(expected_offset>0);	// 1=del; 2=ins
		if (i!=max_pos || expected_offset <0)
			movement_buffer [move_number++] = 0; // matched

		#ifdef indel_debug
		for (j=0; j< windows_number; j++)
			printf("%d ",scores[j][i]);
		printf ("     %c\n", (i==max_pos)?'*':' ') ;
		#endif
	}
	for(; i< read_len; i++)
		movement_buffer [move_number++] = 0; // matched



	read[read_len] = j;
	return move_number;
}

int dynamic_align(char * read, int read_len, gene_value_index_t * index, unsigned int begin_position, int max_indel, char * movement_buffer, int expected_offset, int init_read_offset, int shutdown_read_offset)
// read must be converted to the positive strand.
// movement buffer: 0:match, 1: read-insert, 2: gene-insert, 3:mismatch
// the size of the movement buffer must be equal to the length of the read plus max_indel * 3.
{

	short table[MAX_READ_LENGTH][MAX_READ_LENGTH];
	char table_mask[MAX_READ_LENGTH][MAX_READ_LENGTH];
	int i,j;
	for (i=0; i<read_len; i++)
		for(j=0; j<read_len - expected_offset; j++)
		{
			table_mask[i][j]=0;
			if (j < i - max_indel || j > max_indel + i)
			{
				table[i][j]=-9999;
				continue;
			}

			short from_upper;

			if (i>0) from_upper = table[i-1][j] + (table_mask[i-1][j]?DPALIGN_EXTENDGAP_PENALTY:DPALIGN_CREATEGAP_PENALTY);
			else     from_upper = DPALIGN_CREATEGAP_PENALTY;

			short from_left;

			if (j>0) from_left = table[i][j-1] + (table_mask[i][j-1]?DPALIGN_EXTENDGAP_PENALTY:DPALIGN_CREATEGAP_PENALTY);
			else     from_left = DPALIGN_CREATEGAP_PENALTY;

			char is_matched_ij = gvindex_get(index, begin_position + j) == read[i]?DPALIGN_MATCH_SCORE:DPALIGN_MISMATCH_PENALTY;

			short from_upperleft;
			if (i>0 && j>0) from_upperleft = table[i-1][j-1] + is_matched_ij;
			else            from_upperleft = is_matched_ij; 

			if (from_upperleft <from_left || from_upperleft < from_upper)
				table_mask[i][j]=1;

			table[i][j]=max(0,max(from_upper, max(from_left, from_upperleft)));
		}

	short current_score = table[read_len-1] [read_len - expected_offset-1];
	short path_i = read_len -1;
	int out_pos = 1499;
	char out_tmp [1500];
	int all_steps = 0;
	j = read_len -expected_offset-1;

	out_pos = search_DP_branch(read, read_len, index, begin_position, path_i, j, table, table_mask, max_indel,  out_tmp , expected_offset,  current_score , out_pos, 0, init_read_offset,shutdown_read_offset, &all_steps); 

	if (out_pos)
	{
		memcpy(movement_buffer, out_tmp  + out_pos +1, 1499 - out_pos);
		return 1499 - out_pos;
	}
	else return 0;

}
