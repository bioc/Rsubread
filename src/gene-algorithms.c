#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/timeb.h>
#include "subread.h"
#include "input-files.h"
#include "gene-algorithms.h"
#include "sorted-hashtable.h"

void non_func(const char * fmt, ...)
{
}


int MAX_CIGAR_LEN = 26;


//#define indel_debug

#ifdef indel_debug
#define ddprintf printf
#define ddfflush fflush
#else
#define ddprintf non_func
#define ddfflush(a) 
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

float read_quality_score(char * qualityb, int rl , int format)
{
	int i;
	int qual = 0;
	int testlen = 0;

	char base;

	if(format==FASTQ_PHRED64)
		base = 'B';
	else	base = '#';

	for(i=0; i<rl; i++)
	{
		int testv = qualityb[i] - base;
		if(testv > 1)
		{
			qual += testv;
			testlen ++;
		}
	}
	return qual*1./testlen;
}

int min_matched_bases(char * qualityb, int rl , int format, float match_score)
{
	if(!(qualityb) || !qualityb[0]) return 0;
	int i;

	char base;
	int ret = 0;

	if(format==FASTQ_PHRED64)
		base = 'B';
	else	base = '#';

	for(i=0; i<rl; i++)
		ret += (qualityb[i] - base)<=5;

	return (rl - ret*3/4)*match_score;
}
int bad_quality_base_number(char * qualityb, int rl , int format)
{
	if(!(qualityb) || !qualityb[0]) return 0;
	int ret = 0, i;
	if(FASTQ_PHRED64==format)
	{
		for(i=0; i<rl; i++)
			if(qualityb[i] <='F') ret ++;
	}
	else
		for(i=0; i<rl; i++)
			if(qualityb[i] <='#'+4) ret ++;

		
	return ret;
}


double correct_rate_table [] ={ -1.58147375341 , -0.99684304401 , -0.69552447133 , -0.50767587370 , -0.38013040807 , -0.28926818720 , -0.22255151597 , -0.17255657291 , -0.13455196029 , -0.10536051566 , -0.08276530267 , -0.06517417320 , -0.05141827416 , -0.04062484422 , -0.03213357402 , -0.02543972753 , -0.02015436476 , -0.01597586925 , -0.01266917021 , -0.01005033585 , -0.00797499828 , -0.00632956293 , -0.00502447389 , -0.00398901727 , -0.00316728823 , -0.00251504651 , -0.00199725550 , -0.00158615046 , -0.00125971852 , -0.00100050033 , -0.00079464388 , -0.00063115648 , -0.00050131287 , -0.00039818644 , -0.00031627778 , -0.00025122020 , -0.00019954614 , -0.00015850188 , -0.00012590047 , -0.00010000500 , -0.00007943598 , -0.00006309773 , -0.00005011998 , -0.00003981151 , -0.00003162328 , -0.00002511918 , -0.00001995282 , -0.00001584906 , -0.00001258933 , -0.00001000005 , -0.00000794331 , -0.00000630959 , -0.00000501188 , -0.00000398108 , -0.00000316228 , -0.00000251189 , -0.00000199526 , -0.00000158489 , -0.00000125893 , -0.00000100000 , -0.00000079433 , -0.00000063096 , -0.00000050119 , -0.00000039811 , -0.00000031623 , -0.00000025119 , -0.00000019953 , -0.00000015849 , -0.00000012589 , -0.00000010000 , -0.00000007943 , -0.00000006310 , -0.00000005012 , -0.00000003981 , -0.00000003162 , -0.00000002512 , -0.00000001995 , -0.00000001585 , -0.00000001259 , -0.00000001000 , 0., 0.,0.,0.,0.,0.,0.,0.};

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


unsigned int linear_gene_position(const gene_offset_t* offsets , char *chro_name, unsigned int chro_pos)
{
	unsigned int ret = 0 ;
	int n;
	for (n=0; offsets->read_offset[n]; n++)
	{
		if (strcmp(offsets->read_name[n], chro_name) == 0)
			return ret + chro_pos;
		else
			ret = offsets->read_offset[n];
	}
	return 0xffffffff;
}

unsigned int get_gene_linear(int chrono, int offset, const unsigned int offsets [])
{
	if (chrono>1)return offsets[chrono-1]+offset;
	return offset;
}

int locate_gene_position(unsigned int linear, const gene_offset_t* offsets , char ** chro_name, unsigned int * pos)
{
	int n = 0;

	int total_offsets = offsets -> total_offsets;

	#define GENE_LOCATE_JUMP 15

	while(n+GENE_LOCATE_JUMP < total_offsets &&  offsets->read_offset[n+GENE_LOCATE_JUMP] <= linear)
		n+=GENE_LOCATE_JUMP;

	for (; offsets->read_offset[n]; n++)
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
	bzero(allvote -> span_coverage ,  allvote -> max_len* sizeof(char));
}

void destory_allvote(gene_allvote_t* allvote)
{
	free(allvote -> max_positions);
	free(allvote -> max_votes);
	free(allvote -> max_quality);
	free(allvote -> max_final_quality);
	free(allvote -> masks);
	#ifdef REPORT_ALL_THE_BEST
	free(allvote -> best_records);
	#endif
	free(allvote -> is_counterpart);
	free(allvote -> span_coverage);
	if(allvote -> max_indel_recorder)
		free(allvote -> max_indel_recorder);
}

int init_allvote(gene_allvote_t* allvote, int expected_len, int allowed_indels)
{
	int is_OK = 0;
	allvote -> max_len = expected_len; 
	allvote -> max_positions = (unsigned int *) malloc(sizeof(int)*expected_len);
	allvote -> max_votes = (gene_vote_number_t *) calloc(sizeof(gene_vote_number_t), expected_len);
	allvote -> max_quality = (gene_quality_score_t *) calloc(sizeof(gene_quality_score_t), expected_len);
	allvote -> max_final_quality = (gene_quality_score_t *) calloc(sizeof(gene_quality_score_t), expected_len);
	allvote -> masks = (short *) calloc(sizeof(short), expected_len);
#ifdef REPORT_ALL_THE_BEST
	allvote -> best_records = (gene_best_record_t *) malloc(sizeof(gene_best_record_t)* expected_len);
#endif
	allvote -> is_counterpart = (unsigned char *) malloc(expected_len);

	allvote -> max_indel_tolerance = allowed_indels;
	allvote	-> indel_recorder_length = max(3*(allowed_indels+1)+1, MAX_CIGAR_LEN+2);
	allvote -> span_coverage = (char *)calloc(sizeof(char), expected_len);

	if((allvote -> max_quality &&  allvote -> max_positions  && allvote -> max_votes  && allvote -> max_final_quality && allvote -> masks && allvote -> is_counterpart && allvote -> span_coverage))
		is_OK = 1;

	if(allowed_indels && is_OK)
	{
		allvote -> max_indel_recorder = (char *)malloc( allvote -> indel_recorder_length *expected_len);
		is_OK = (allvote -> max_indel_recorder!=NULL);
	}
	else	allvote -> max_indel_recorder =  NULL;

	if(is_OK)
		return 0;
	else
	{
		puts(MESSAGE_OUT_OF_MEMORY);
		return 1;
	}
}


void compress_cigar(char *cigar, int total_length, char * read)
{
	char tmp[200];
	char cigar_piece [10];

	int cigar_len = strnlen(cigar,MAX_CIGAR_LEN);
	cigar[cigar_len]=0;

	int i;
	int tmpv = 0;

	char last_operation = 'X';
	int last_tmpv = 0;
	tmp[0]=0;
	int cigar_length=0;
	
	for(i=0; i < cigar_len; i++)
	{
		char cc = cigar[i];
		if(isdigit(cc))
		{
			tmpv=tmpv*10+(cc-'0');
		}
		else if(cc=='-')
		{
			last_tmpv = 0;
			cigar_length = 0;
			break;
		}
		else
		{
			if((cc!=last_operation) && (last_operation!='X'))
			{
				sprintf(cigar_piece,"%d%c", last_tmpv, last_operation); 
				if(last_operation == 'M' || last_operation == 'S' || last_operation == 'I')
					cigar_length += last_tmpv;
				strcat(tmp, cigar_piece);
				last_tmpv = 0;
			}
			last_tmpv += tmpv;
		
			tmpv = 0;
			last_operation = cc;
		
		}
	
	}

	if(last_tmpv)
	{
		sprintf(cigar_piece,"%d%c", tmpv+last_tmpv, last_operation); 
		strcat(tmp, cigar_piece);
		if(last_operation == 'M' || last_operation == 'S' || last_operation == 'I')
			cigar_length += tmpv+last_tmpv;
	}
	//printf("\nCIG=%s\n", cigar);
	if(cigar_length == total_length)
		strcpy(cigar, tmp);
	else
	{
		//printf("CIG=%s; read=%s\n", cigar, read);
		sprintf(cigar, "%dM", total_length);
	}
}

void show_cigar(char * info, int len, int is_reversed_map, char * buf, int indel_tolerance, int total_subreads, char *read)
{
	int i, is_error = 0;
	int last_end=0, last_offset = 0, cursor = 0;

	//printf("L0=%d\n", info[0]);
	if(info[0]==-1){
		sprintf(buf, "%dM", len);
		return;
	}
	else if(info[0]==-2){
		if (strchr( info+1, '-'))
			sprintf(buf, "%dM", len);
		else
		{
			strncpy(buf, info+1, 98);
			compress_cigar(buf, len, read);
		}
		return;
	}
	else if(info[0]==-3){
		info ++;
		is_error = 1;
	}

	for(i=0; i<indel_tolerance*3; i+=3)
	{
		if (!info[i])break;
		int dist = info[i+2];
		//int subread_start = info[i]-1;
		int subread_end = info[i+1]-1;

		//int base_start = (last_end==0)?0:(find_subread_end(len, TOTAL_SUBREADS, subread_start)-16);
		int base_end = (i < indel_tolerance*3-3 && info[i+3])?find_subread_end(len, total_subreads, subread_end):len;

//			printf("BE=%d ; II+3=%d ; len=%d\n", base_end, info[i+3], len);

		int offset = last_offset-dist;
		if (base_end - cursor - (offset>0?offset:0) <0) 
		{
			buf[0]=0;

			cursor = 0;
			break;
		}
		if (i >0)
			sprintf(buf+strlen(buf), "%d%c%dM", abs(offset), offset>0?'I':'D', base_end - cursor - (offset>0?offset:0));
		else
			sprintf(buf+strlen(buf), "%dM", base_end);
		last_offset = dist;
		last_end = base_end;
		cursor = base_end;
		//if (dist >0)cursor += dist;
	}
	if(cursor < len)
		sprintf(buf+strlen(buf), "%dM", len-cursor);
//	if(is_error)
//		strcat(buf, "XX");
}

void add_allvote_q(gene_allvote_t* allvote,int qid , int pos, gene_vote_number_t votes, gene_quality_score_t quality, int is_counterpart, short mask, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads, int space_type, int report_junction, int is_head_high_quality, char * qual_txt, int phred_version, char span_coverage)
{

	int is_add_new = 0;
	if((votes > allvote -> max_votes[qid] ) || (votes  == allvote -> max_votes[qid]  && span_coverage > allvote -> span_coverage[qid]) || (votes == allvote -> max_votes[qid] && (span_coverage == allvote -> span_coverage[qid]) && quality > allvote -> max_quality[qid]))
	{
		is_add_new = 1;
	}
	else if (votes == allvote -> max_votes[qid] && quality == allvote -> max_quality[qid] && span_coverage == allvote -> span_coverage[qid])
	{
		//printf("ADD1\n");
		allvote -> masks[qid]|= IS_BREAKEVEN_READ;
		if(pos < allvote -> max_positions[qid] ) is_add_new=1;
	}

	if(is_add_new)
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
		//printf("ADD ALL: POS=%u > %u VOTES=%d > %d ; Q=%f > %f ; SPAN=%d > %d\n", pos, allvote -> max_positions[qid], votes,  allvote->max_votes[qid], quality, allvote -> max_quality[qid], span_coverage, allvote -> span_coverage[qid]);
		allvote -> is_counterpart[qid] = is_counterpart;
		allvote -> max_positions[qid] = pos;
		allvote -> max_votes[qid] = votes;
		allvote -> max_quality[qid] =  quality;
		allvote -> masks[qid] = mask;
		allvote -> span_coverage[qid] = span_coverage;

		/*if(mask & IS_BREAKEVEN_READ)		printf("ADD2\n");
		else printf("NOADD3\n");*/
	

		int mismatch=0;

		if(max_indel>0)
		{
			find_and_explain_indel(allvote, qid, pos, votes, quality, is_counterpart, mask, max_indel_recorder, array_index, read_txt, read_len, max_indel, total_subreads, space_type, report_junction, is_head_high_quality, qual_txt, phred_version);

			if(allvote -> max_indel_recorder)
			{
				char cigar_str [100];
				cigar_str[0]=0;
				show_cigar(allvote->max_indel_recorder + qid * allvote -> indel_recorder_length, read_len, is_counterpart, cigar_str, max_indel, total_subreads, read_txt);
				allvote->max_final_quality[qid] = final_mapping_quality(array_index, allvote -> max_positions[qid], read_txt, qual_txt, cigar_str, phred_version, &mismatch);
			}
		}
		else
		{
			char cigar_str [100];
			sprintf(cigar_str,"%dM", read_len);	
			allvote->max_final_quality[qid] = final_mapping_quality(array_index, allvote -> max_positions[qid], read_txt, qual_txt, cigar_str, phred_version, &mismatch);
		}


	}

}


float EXON_RECOVER_MATCHING_RATE = 0.9;
float EXON_INDEL_MATCHING_RATE_TAIL = 0.8;
float EXON_INDEL_MATCHING_RATE_HEAD = 0.8;

void find_and_explain_indel(gene_allvote_t* allvote,int qid , int pos, gene_vote_number_t votes, gene_quality_score_t quality, int is_counterpart, char mask, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads, int space_type, int report_junction,  int is_head_high_quality, char * qual_txt, int phred_version)
{
	//unsigned int pos0 = pos;
	if(allvote -> max_indel_recorder)
	{
		//PNT111
		short head_indel_pos=-1 , tail_indel_pos=-1;
		int head_indel_movement=0, tail_indel_movement=0;
		if (array_index)
		{
			int k;
			if (max_indel_recorder[0])
			{
				int cover_start = find_subread_end(read_len, total_subreads, max_indel_recorder[0]-1)- 15;

				for (k=0; max_indel_recorder[k]; k+=3);
				int cover_end = find_subread_end(read_len, total_subreads, max_indel_recorder[k-2]-1) + max(0,-max_indel_recorder[k-1]);

				float head_qual = is_head_high_quality?999:0, tail_qual = is_head_high_quality?0:999;


				if(qual_txt && qual_txt[0])
				{
					int j;
					head_qual = 0;
					tail_qual = 0;
					for (j = 0; j<cover_start; j++)
					{
						if(FASTQ_PHRED64 == phred_version)
						{
							head_qual += (1-get_base_error_prob64(qual_txt[j]));
						}
						else
						{
							head_qual += (1-get_base_error_prob33(qual_txt[j]));
						}
					}
					for (j = 0; j<read_len - cover_end; j++)
					{
						if(FASTQ_PHRED64 == phred_version)
						{
							tail_qual += (1-get_base_error_prob64(qual_txt[read_len - j -1]));
						}
						else
						{
							tail_qual += (1-get_base_error_prob33(qual_txt[read_len - j -1]));
						}
					}
				}

				int head_must_correct = 4;
				EXON_INDEL_MATCHING_RATE_HEAD = 0.92;
				if (head_qual / cover_start < 0.95 )
				{
					EXON_INDEL_MATCHING_RATE_HEAD = 0.85;
					head_must_correct =3;
				}
				if (head_qual / cover_start < 0.85)
				{
					EXON_INDEL_MATCHING_RATE_HEAD = 0.75;
					head_must_correct =2;
				}

				int tail_must_correct = 4;
				EXON_INDEL_MATCHING_RATE_TAIL = 0.92;
				if (tail_qual / (read_len-cover_end) < 0.95)
				{
					EXON_INDEL_MATCHING_RATE_TAIL = 0.85;
					tail_must_correct =3;
				}
				if (tail_qual / (read_len-cover_end) < 0.85)
				{
					EXON_INDEL_MATCHING_RATE_TAIL = 0.75;
					tail_must_correct =2;
				}

				//printf("HQ=%.4f in %d ; TQ=%.4f in %d\n%s\n", head_qual, cover_start, tail_qual, (read_len-cover_end), qual_txt);

				

				int is_full_covered = 0;
				is_full_covered = extend_covered_region(array_index, pos, read_txt, read_len, cover_start, cover_end, 4, head_must_correct, tail_must_correct,  max_indel, space_type, max_indel_recorder[k-1], &head_indel_pos, &head_indel_movement, &tail_indel_pos, &tail_indel_movement, is_head_high_quality, qual_txt, phred_version);
				if(head_indel_movement)
				{
					pos += head_indel_movement;
					allvote -> max_positions[qid] = pos; 
				}
				if(is_full_covered == 3)
				{
					allvote -> masks[qid] &= ~IS_RECOVERED_JUNCTION_READ;
				}
				else	allvote -> masks[qid] |= IS_RECOVERED_JUNCTION_READ;
			}
		}

		if (array_index)
			*(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length) = 0xff ;
		else
			memcpy(allvote -> max_indel_recorder + qid *allvote -> indel_recorder_length, max_indel_recorder, allvote -> max_indel_tolerance *3* sizeof(char));


		if (array_index && (max_indel_recorder[3] || head_indel_pos>=0 || tail_indel_pos>0)){
			int head_pos = 0;
			if((head_indel_pos>=0 && report_junction) || head_indel_movement)head_pos = head_indel_pos ;//- max(0,head_indel_movement);
			int tail_pos = read_len;
			if((tail_indel_pos>0 && report_junction) || tail_indel_movement)tail_pos = tail_indel_pos ;//- min(0, tail_indel_movement);
			//printf ("\nEE H=%d T=%d HM=%d TM=%d\n",head_pos, tail_pos, head_indel_movement, tail_indel_movement);
			explain_indel_in_middle(allvote, qid , pos, max_indel_recorder,  array_index, read_txt, read_len,  max_indel, total_subreads, head_pos, tail_pos , head_indel_movement, tail_indel_movement, report_junction);
			//if(*(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length) != 0xfd ) allvote -> max_positions[qid] = pos0;
		}

	}

}

void explain_indel_in_middle(gene_allvote_t* allvote, int qid , int pos, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads, int head_start_point, int tail_end_point, int head_indel_movement, int tail_indel_movement, int report_junction)
{

	char indel_operations[1500];
	int i,xx;
	char tmp_cigar [MAX_CIGAR_LEN+1];

	int current_pos = head_start_point; 
	int explain_cursor = head_start_point;
	char last_operation = 0;
	int del_number = 0;
	int dynamic_delta = 0;

	tmp_cigar[0]=0;

	if(current_pos >0)
	{
		if (head_indel_movement != 0)
			sprintf(tmp_cigar, "%dM%d%c", head_start_point  - max(head_indel_movement, 0), abs(head_indel_movement), head_indel_movement>0?'I':'D');
		else if (report_junction)
			sprintf(tmp_cigar, "%dS", head_start_point);
	}

	ddprintf ("R=%s; REC[8]=%d [9]=%d\n", read_txt, max_indel_recorder[8], max_indel_recorder[9]);

	for (i=3; max_indel_recorder[i]; i+=3)
	{
		if (i >= max_indel*3) break;

		int last_dist = max_indel_recorder[i-1];
		int black_subread_start = max_indel_recorder[i-2]-1;
		int black_subread_end = max_indel_recorder[i]-1;
		if (black_subread_end < black_subread_start+1) black_subread_end = black_subread_start+1;
		while(max_indel_recorder[i+3])
		{
			if (max_indel_recorder[i+3]- black_subread_end>2 || i >= max_indel*3-3) break;
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

		
		int moves = 0;
		if(black_base_end > black_base_start)
			moves = dynamic_align (read_txt + black_base_start, black_base_end-black_base_start , array_index, pos + blackref_base_start, max_indel, indel_operations, exp_indel,  -10, black_base_end - black_base_start+5);

#ifdef indel_debug
		char tt = *(read_txt + black_base_end);
		read_txt [black_base_end] = 0;
		ddprintf ("%s\n", read_txt + black_base_start);
		read_txt [black_base_end] = tt;
		ddprintf("\n TESTING SUBR %d to %d BASE %d to %d CHRO %d to %d exp offset %d\n",black_subread_start, black_subread_end, black_base_start, black_base_end, blackref_base_start, -max(0,exp_indel)+blackref_base_start +(black_base_end-black_base_start), exp_indel);
		ddprintf(" moves %d\n", moves);
		ddfflush(stdout) ;
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
				if(current_operation == 3) current_operation = 0;

				if(current_operation != last_operation && (current_pos < gap_end_read-1+1 || current_operation == 0))// && current_pos!=explain_cursor)
				{
					int vpos = strlen(tmp_cigar);
					if (vpos>MAX_CIGAR_LEN-6) break;
					sprintf(tmp_cigar + vpos, "%d%c", last_operation==1?del_number:(current_pos - explain_cursor), last_operation==0?'M':(last_operation==1?'D':'I'));
					explain_cursor = current_pos ;
					//printf("EXP=%d of %d;  XX=%d of %d;  MOV=%d->%d ; vpos=%d>%d cpos=%d > gapend=%d\nS=%s\n", explain_cursor, gap_end_read , xx , moves, last_operation,current_operation, vpos, MAX_CIGAR_LEN-6, current_pos, gap_end_read-1, tmp_cigar);
					del_number = 0;
				}
				if(current_pos >= gap_end_read-1 || xx == moves - 1){
					if(current_operation != 0 && current_pos!=explain_cursor)
					{
						int vpos = strlen(tmp_cigar);
						sprintf(tmp_cigar + vpos, "%d%c", current_operation==1?del_number:(current_pos - explain_cursor), current_operation==1?'D':'I');
						explain_cursor = current_pos;
						del_number=0;
					}
					break;
				}

				if (indel_operations[xx]==0) current_pos ++;
				else if (indel_operations[xx]==1) {del_number++; dynamic_delta++;}//"D"
				else if (indel_operations[xx]==2) {current_pos ++; dynamic_delta--;} //"I"
				else if (indel_operations[xx]==3) current_pos ++;


				last_operation = current_operation; 
			} 
		}
		else
		{
			int movement = last_dist - max_indel_recorder[i+2];
			int vpos = strlen(tmp_cigar);
			if (vpos>MAX_CIGAR_LEN-6) break;

			current_pos = find_subread_end(read_len, total_subreads,max_indel_recorder[i]-1)-15 + (max_indel_recorder[i+2] >0?max_indel_recorder[i+2]:-max_indel_recorder[i+2]);
			current_pos -= (black_base_end - black_base_start)/2 -3;

			if(movement)
				sprintf(tmp_cigar + vpos, "%dM%d%c", current_pos - explain_cursor , abs(movement), movement<0?'D':'I' );
			else if (current_pos - explain_cursor>0)
				sprintf(tmp_cigar + vpos, "%dM", current_pos - explain_cursor );

			explain_cursor = current_pos + (movement >0?movement:0);
			current_pos = find_subread_end(read_len, total_subreads,max_indel_recorder[i+1]-1);
			last_operation = 0;
			del_number = 0;
			moves = 0;
		}
		last_operation = 0;
		current_pos = black_base_end;
	//	if (!moves)
			//strcat(tmp_cigar,"X");
			//printf("Wrong explaination! dist=%d   dyn=%d R=%s\n", max_indel_recorder[i+2], dynamic_delta, read_txt);
			
		
	}

	int vpos = strlen(tmp_cigar);
	if (vpos > MAX_CIGAR_LEN-9 )
	{
		*(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length) = 0xfd; 
		memcpy(allvote -> max_indel_recorder + qid * allvote -> indel_recorder_length+1, max_indel_recorder, 3 * max_indel);
	}
	else
	{
		if (explain_cursor<read_len)	
		{
			if(tail_end_point > explain_cursor)
				sprintf(tmp_cigar+vpos,"%dM", tail_end_point - explain_cursor);
			vpos = strlen(tmp_cigar);
			if(tail_end_point< read_len)
			{
				if (tail_indel_movement!=0) 
				{
					int tail_len_m = read_len -(tail_end_point - min(tail_indel_movement, 0));
					if(tail_len_m)
						sprintf(tmp_cigar+vpos, "%d%c%dM", abs(tail_indel_movement), tail_indel_movement<0?'I':'D', tail_len_m);
					else
						sprintf(tmp_cigar+vpos, "%d%c", abs(tail_indel_movement), tail_indel_movement<0?'I':'D');
				}
				else if (report_junction)
					sprintf(tmp_cigar+vpos, "%dS",read_len - tail_end_point );
			}
		}
	//	if (head_start_point >0 || tail_end_point < read_len)
	//		strcat(tmp_cigar,"X");

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

int match_chro_indel(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type, int indel_size)
{

	int i;
	int ret = 0;
	for (i = -indel_size; i<=indel_size ; i++)
	{
		if (pos +i + test_len < index -> start_base_offset + index -> length && pos > -i)
			ret += match_chro(read, index, pos + i, test_len, is_negative_strand, space_type);
	}
	return ret;

}



void mark_votes_array_index(char * read_str, int read_len, gene_vote_t * dest, gene_vote_t * src, gene_value_index_t * my_array_index, int color_space, int indel_tolerance, int min_minor, const char quality_str [], int quality_scale, int is_negative_strand)
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
				//matchingness_count = match_read(read_str, read_len, potential_position, my_array_index, color_space, indel_tolerance, quality_str, quality_scale);
				matchingness_count = 1.*match_chro_indel(read_str, my_array_index, potential_position, read_len, 0, color_space, indel_tolerance);
	
			dest -> pos[i][j] = potential_position;
			dest -> quality[i][j] = matchingness_count;
			//printf ("Q=%.4f\n", matchingness_count);
			dest -> votes [i][j] =src -> votes [i][j];
			dest -> masks[i][j] = src -> masks[i][j];
			dest -> coverage_start[i][j] = src -> coverage_start[i][j];
			dest -> coverage_end[i][j] = src -> coverage_end[i][j];
			memcpy(dest -> indel_recorder[i][j], src -> indel_recorder[i][j], MAX_INDEL_TOLERANCE*3*sizeof(char));
			// find for the `best positions'; only the best positions are replaced by the base machingness scores.
			if((matchingness_count >  dest -> max_quality && src -> votes [i][j] == dest -> max_vote) || (src -> votes [i][j] > dest -> max_vote))
			{
				memcpy(dest -> max_indel_recorder, src -> indel_recorder[i][j], MAX_INDEL_TOLERANCE*3*sizeof(char));
				dest -> max_vote  = src -> votes [i][j];
				dest -> max_mask  = src -> masks[i][j];
				dest -> max_quality = matchingness_count;
				dest -> max_position = potential_position;
				dest -> max_coverage_start = src -> coverage_start[i][j];
				dest -> max_coverage_end = src -> coverage_end[i][j];
			}
		}
	}
}

int select_positions_array(char * read1_str, int read1_len, char * read2_str, int read2_len,  gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2,  gene_quality_score_t * sum_quality, gene_quality_score_t * qual_r1, gene_quality_score_t * qual_r2, char * r1_recorder, char * r2_recorder, gehash_data_t * pos_read1, gehash_data_t * pos_read2, unsigned int max_pair_dest, unsigned int min_pair_dest, int min_major, int min_minor,int is_negative_strand, gene_value_index_t * my_array_index, int color_space, int indel_tolerance, int number_of_anchors_quality, const char quality_str1 [], const char quality_str2 [], int quality_scale, int max_indel_len, int * is_break_even)
{
	gene_vote_t base_vote_1 , base_vote_2;

	mark_votes_array_index(read1_str, read1_len, &base_vote_1, vote_read1, my_array_index, color_space, indel_tolerance, min_minor, quality_str1, quality_scale, is_negative_strand);
	mark_votes_array_index(read2_str, read2_len, &base_vote_2, vote_read2, my_array_index, color_space, indel_tolerance, min_minor, quality_str2, quality_scale, is_negative_strand);

//	printf("Q=%.3f\n", *qual_r2);

	return select_positions(&base_vote_1, &base_vote_2, numvote_read1, numvote_read2, sum_quality, qual_r1, qual_r2,  pos_read1, pos_read2, r1_recorder , r2_recorder ,max_pair_dest,  min_pair_dest,  min_major,  min_minor, is_negative_strand, number_of_anchors_quality,max_indel_len, is_break_even);
}


int select_positions(gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2, gene_quality_score_t * sum_quality, gene_quality_score_t * qual_r1, gene_quality_score_t * qual_r2 , gehash_data_t * pos_read1, gehash_data_t * pos_read2, char * read1_indel_recorder, char * read2_indel_recorder, unsigned int max_pair_dest, unsigned int min_pair_dest, int min_major, int min_minor,int is_negative_strand, int number_of_anchors_quality, int max_indel_len, int * is_breakeven)
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
	char is_minor_breakeven [ANCHORS_NUMBER];

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
	bzero(is_minor_breakeven, ANCHORS_NUMBER);
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

						long long int abdist_old = minor_position[l];
						abdist_old -=  anchors_position[l];
						if(abdist_old <0) abdist_old=-abdist_old;
						if (	(minor_votes[l] < current_vote->votes[i][j] ||
							(minor_votes[l] == current_vote->votes[i][j] && current_vote->quality [i][j] > minor_quality[l]) ||
							(minor_votes[l] == current_vote->votes[i][j] && abs(current_vote->quality [i][j] - minor_quality[l])<0.0001 && abdist < abdist_old))&& 
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
	//int reason = 0;

	for (k=0; k<anchors; k++)
	{
		if(minor_votes[k]==0)
			continue;

		int need_replace = 0;

		if ((minor_votes[k]+anchors_votes[k]) > (selection_minor_vote+selection_major_vote) ||
	   	       ((minor_votes[k]+anchors_votes[k]) == (selection_minor_vote+selection_major_vote) &&
		 	 minor_quality[k] + anchors_quality[k] > *sum_quality))
		{
			need_replace = 1;
		}
		else if (minor_votes[k]+anchors_votes[k] == selection_minor_vote+selection_major_vote &&
		 	 minor_quality[k] + anchors_quality[k] == *sum_quality)
		{
			if (selection_minor_pos != minor_position[k] && minor_position[k] != selection_major_pos)
			{
				* is_breakeven = 1;
				if (selection_major_pos >  anchors_position[k]) need_replace = 1;
			}
		}

		if(need_replace){
			selection_minor_pos = minor_position[k];
			selection_major_pos = anchors_position[k];
			selection_minor_vote = minor_votes[k];
			selection_major_vote = anchors_votes[k];
			selection_major_read_no = anchor_read[k];
			selected_major_bucket = anchors_buckets[k];
			selected_minor_bucket = minor_buckets[k];
			selected_major_index = anchors_index[k];
			selected_minor_index = minor_index[k];
			* is_breakeven = is_minor_breakeven[k];
			* qual_r1 = anchor_read[k] ? minor_quality[k]:anchors_quality[k];
			* qual_r2 = anchor_read[k] ? anchors_quality[k]: minor_quality[k];
			* sum_quality = minor_quality[k] + anchors_quality[k];
	
		}
	}
//	if (*is_breakeven)printf("BKF3 R=%d\n", reason);

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

		read_line(299,fp, fn, 0);
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

		if( n >= OFFSET_TABLE_SIZE -1) break;
	}

	offsets->total_offsets=n;
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


void final_matchingness_scoring(const char read_str[], const char quality_str[], int read_len, gene_vote_t * vote, gehash_data_t * max_position, gene_vote_number_t * max_vote, short *max_mask, gene_quality_score_t * max_quality, gene_value_index_t * my_array_index, int space_type, int indel_tolerance, int quality_scale, char * max_indel_recorder, int * max_coverage_start, int * max_coverage_end)
{
	int i, j;
	gene_quality_score_t max_matching = -1.0;

	*max_vote = vote -> max_vote; 

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if( vote->votes[i][j]  < vote -> max_vote ) continue;

			unsigned int potential_position = vote -> pos[i][j];
			gene_quality_score_t matchingness_count = 1.*match_chro_indel((char *)read_str, my_array_index, potential_position, read_len, 0, space_type, indel_tolerance);
			//printf("TQ=%.4f VS MQ=%.4f\n", matchingness_count, max_matching);
			if(matchingness_count > max_matching)
			{
				max_matching = matchingness_count;
				*max_position = potential_position;
				*max_mask = vote -> masks[i][j];
				*max_coverage_start = vote-> coverage_start[i][j];
				*max_coverage_end = vote -> coverage_end[i][j];
				memcpy(max_indel_recorder, vote-> indel_recorder[i][j], MAX_INDEL_TOLERANCE*3);

				* max_quality = max_matching;
			}
			else if (matchingness_count == max_matching)
			{
//				printf("\nBREAK EVEN DETECTED AT SORTED TABLE: %u (kept) and %u\n", vote->max_position, vote -> pos[i][j]);
				(*max_mask) |= IS_BREAKEVEN_READ ;
			}
		}
}

int DPALIGN_CREATEGAP_PENALTY = -2 ;
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

 


int search_DP_branch(char * read, int read_len, gene_value_index_t * index, unsigned int begin_position, int path_i, int path_j, short table[MAX_READ_LENGTH][MAX_READ_LENGTH],char table_mask[MAX_READ_LENGTH][MAX_READ_LENGTH], int max_indel, char * movement_buffer, int expected_offset, int current_score, int out_pos, int current_offset, int init_read_offset, int shutdown_read_offset, int * all_steps)
{

	if (1499 - out_pos > (read_len << 2) || (*all_steps) > 3000 + (read_len << 5) || (*all_steps) > 57000){
		ddprintf("\nTOO MANY STEPS: rl = %d    len=%d    steps=%d\n", read_len, 1499 - out_pos, *all_steps);
		ddfflush(stdout) ;
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
	if ((!found ) && (left_score + (table_mask[path_i][path_j-1]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) == current_score || (current_score == 0 && left_score + (table_mask[path_i][path_j-1]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) < 0)))
	{
		movement_buffer[out_pos] = 1;
		(*all_steps) ++;
		ddprintf ("Path_i = %d ; Path_j = %d ; Offset = %d => 1\n", path_i, path_j, current_offset);
		ddfflush(stdout) ;
		found =search_DP_branch (read, read_len, index, begin_position, path_i , path_j -1, table  , table_mask, max_indel, movement_buffer, expected_offset, left_score, out_pos -1, current_offset - ((path_i >= init_read_offset && path_i <= shutdown_read_offset)?1:0),  init_read_offset, shutdown_read_offset, all_steps); 
	}
	if ((!found ) && (upper_score + (table_mask[path_i-1][path_j]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) == current_score || (current_score == 0 && upper_score + (table_mask[path_i-1][path_j]? DPALIGN_EXTENDGAP_PENALTY: DPALIGN_CREATEGAP_PENALTY) < 0 )))
	{
		movement_buffer[out_pos] = 2;
		(*all_steps) ++;
		ddprintf ("Path_i = %d ; Path_j = %d ; Offset = %d => 2\n", path_i, path_j, current_offset);
		ddfflush(stdout) ;
		found =search_DP_branch (read, read_len, index, begin_position, path_i -1 , path_j , table , table_mask , max_indel, movement_buffer, expected_offset, upper_score, out_pos -1, current_offset + ((path_i >= init_read_offset && path_i <= shutdown_read_offset)?1:0),  init_read_offset, shutdown_read_offset, all_steps); 
	}
	if ((!found ) && (is_matched_ij + (upperleft_score - current_score) == 0))
	{
		movement_buffer[out_pos] = is_matched_ij==2?0:3;
		ddprintf ("Path_i = %d ; Path_j = %d ; Offset = %d => 0\n", path_i, path_j, current_offset);
		(*all_steps) ++;
		found = search_DP_branch (read, read_len, index, begin_position, path_i -1 , path_j -1, table , table_mask, max_indel, movement_buffer, expected_offset, upperleft_score, out_pos -1, current_offset, init_read_offset, shutdown_read_offset, all_steps); 
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
	ddfflush(stdout) ;
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

	#ifdef indel_debug
	int ii, jj;
	for(ii = 0; ii<read_len - expected_offset; ii++)
		putchar(gvindex_get(index, begin_position + ii));

	printf ("\n%s\n", read);
	#endif

	for (i=0; i<read_len ; i++)
	{
		for(j=0; j<read_len - expected_offset; j++)
		{
			table_mask[i][j]=0;
			if (j < i - max_indel || j > max_indel + i)
			{
				table[i][j]=-9999;
			#ifdef indel_debug
				putchar('\t');
			#endif
				continue;
			}

			short from_upper;

			if (i>0) from_upper = table[i-1][j] + (table_mask[i-1][j]?DPALIGN_EXTENDGAP_PENALTY:DPALIGN_CREATEGAP_PENALTY);
			else     from_upper = DPALIGN_CREATEGAP_PENALTY;

			short from_left;

			if (j>0) from_left = table[i][j-1] + (table_mask[i][j-1]?DPALIGN_EXTENDGAP_PENALTY:DPALIGN_CREATEGAP_PENALTY);
			else     from_left = DPALIGN_CREATEGAP_PENALTY;

			char chromo_ch = gvindex_get(index, begin_position + j);
			char is_matched_ij = (chromo_ch == read[i])?DPALIGN_MATCH_SCORE:DPALIGN_MISMATCH_PENALTY;

			short from_upperleft;
			if (i>0 && j>0) from_upperleft = table[i-1][j-1] + is_matched_ij;
			else            from_upperleft = is_matched_ij; 
			if ((i ==0 || j ==0) && (i+j>0)) from_upperleft += DPALIGN_CREATEGAP_PENALTY;

			if (from_upperleft <=from_left || from_upperleft <= from_upper)
				table_mask[i][j]=1;

			table[i][j]=max(0,max(from_upper, max(from_left, from_upperleft)));
			#ifdef indel_debug
			printf("%c%c\t", chromo_ch, read[i]);
			//printf("%d,%d,%d\t", from_left,from_upperleft,from_upper);
			#endif
		}
		#ifdef indel_debug
		puts("");
		#endif

	}
	#ifdef indel_debug
	//puts("");
	//puts("");

	#endif
	//printf ("EXP_OFFSET=%d    RLEN=%d\n", read_len, expected_offset);

	short current_score = table[read_len-1] [read_len - expected_offset-1];
	short path_i = read_len-1 ;
	int out_pos = 1499;
	char out_tmp [1500];
	int all_steps = 0;
	j = read_len-expected_offset -1;

	#ifdef indel_debug

	for(ii=0;ii< 20; ii++)
	{
		for(jj=0; jj<20; jj++)
			printf("%d\t",table[ii][jj]);
		puts("");
	}
		puts("");
		puts("");

	for(ii=0;ii< 20; ii++)
	{
		for(jj=0; jj<20; jj++)
			printf("%d\t",table_mask[ii][jj]);
		puts("");
	}
	#endif

	out_pos = search_DP_branch(read, read_len, index, begin_position, path_i, j, table, table_mask, max_indel,  out_tmp , expected_offset,  current_score , out_pos, 0, init_read_offset,shutdown_read_offset, &all_steps); 

	if (out_pos)
	{
		memcpy(movement_buffer, out_tmp  + out_pos +1, 1499 - out_pos);
		return 1499 - out_pos;
	}
	else return 0;

}

int extend_covered_region(gene_value_index_t *array_index, unsigned int read_start_pos, char * read, int read_len, int cover_start, int cover_end, int window_size, int req_match_5end , int req_match_3end, int indel_tolerance, int space_type, int tail_indel, short * head_indel_pos, int * head_indel_movement, short * tail_indel_pos, int * tail_indel_movement, int is_head_high_quality, char * qual_txt, int qual_format)
{
	int ret = 0;
	*head_indel_pos = -1;
	*tail_indel_pos = -1;
	if (cover_start >= window_size && EXON_INDEL_MATCHING_RATE_HEAD < 1.0001) 
	{
		int head_test_len =  cover_start;
		int roughly_mapped = match_chro(read, array_index, read_start_pos, head_test_len , 0, space_type);
		if (roughly_mapped >= head_test_len  * EXON_RECOVER_MATCHING_RATE - 0.0001)
		{
			//printf ("HEAD: ROUGHLY MAPPED %d in %d!\n", roughly_mapped,head_test_len );
			ret |= 1;
		}
		else
		{
			int window_end_pos = cover_start + window_size-1;
			int not_too_bad = 1;
			int right_match_number=0;

			//printf ("HEAD: ROUGHLY UNMAPPED %d in %d!\n", roughly_mapped,head_test_len );
			while (1)
			{
				int matched_bases_in_window = match_chro_wronglen(read+ window_end_pos - window_size, array_index, read_start_pos + window_end_pos - window_size, window_size, space_type, NULL, &right_match_number);
				int best_indel_movement = -1000;
				int best_indel_pos = -1;

				if (matched_bases_in_window >= req_match_5end) // this window is not bad enough so that an indel is considered
					window_end_pos--;
				else
				{
					roughly_mapped = match_chro(read, array_index, read_start_pos, window_end_pos - right_match_number , 0, space_type);
					//printf("\nMATCH: %d in %d\n", roughly_mapped , (window_end_pos - window_size/2));
					if (roughly_mapped < (int)(0.5+ (window_end_pos - right_match_number)  * EXON_RECOVER_MATCHING_RATE))
					{
						// the quality of this window is too low (at most 1 base is matched). I must consider if it is an indel.
						int indel_movement_i ;

						int best_matched_bases_after_indel = -1;

						not_too_bad = 0;
						for (indel_movement_i = 0; indel_movement_i < 2* indel_tolerance-1 ; indel_movement_i ++)
						{
							int indel_movement = (indel_movement_i+1)/2 * (indel_movement_i %2?1:-1);
							int test_length = window_end_pos /*- 1*/ - max(0, indel_movement) -  right_match_number;

							if (test_length < window_size) continue;
							//if (test_length <= 1+ abs(indel_movement)/4) continue;
							//test_length = min(10, test_length);
							if (abs(indel_movement) > indel_tolerance) continue;

							int test_start =0;// window_end_pos - max(0, indel_movement) -  right_match_number - test_length;

							int matched_bases_after_indel = match_chro_support(read +test_start, array_index, read_start_pos + indel_movement +test_start, test_length,0, space_type, qual_txt, qual_format);
							//printf("MOV=%d ; MATCH=%d ; TLEN=%d\n", indel_movement , matched_bases_after_indel, test_length);
							float test_rate = EXON_INDEL_MATCHING_RATE_HEAD;
							
							if(test_length < 3) test_rate = 1;

							if(best_matched_bases_after_indel <matched_bases_after_indel  && matched_bases_after_indel >= (int)(0.5+test_length * test_rate))
							{
								not_too_bad = 1;
								//printf ("0WSP=%d, RHT_MAT=%d\n", window_end_pos, right_match_number);

								best_matched_bases_after_indel = matched_bases_after_indel;
								best_indel_movement = indel_movement;
								best_indel_pos = window_end_pos - right_match_number;

								//printf("HEAD: RECOVERED : %d/%d  ; BINDEL=%d , MOVE=%d\n",matched_bases_after_indel , test_length, best_indel_pos , best_indel_movement);

								
								*head_indel_pos = best_indel_pos;
								*head_indel_movement = best_indel_movement;
							}
						}
						if(best_indel_pos<0) *head_indel_pos =  window_end_pos - right_match_number;
						break;
					}else window_end_pos--;
				}
				if (window_end_pos - window_size <= 0) break;
			}
			if(not_too_bad)
				ret |=1;
			//else
				//*head_indel_pos = -1;
		}
	}
	else ret |= 1;



	if (cover_end <= read_len - window_size && EXON_INDEL_MATCHING_RATE_TAIL < 1.0001) 
	{
		int tail_test_len = read_len - cover_end;
		int roughly_mapped = match_chro(read + cover_end, array_index, read_start_pos + cover_end + tail_indel, tail_test_len , 0, space_type);
		if (roughly_mapped >= tail_test_len  * EXON_RECOVER_MATCHING_RATE - 0.0001)
		{
			//printf ("TAIL: ROUGHLY MAPPED %d in %d! T=%f\n", roughly_mapped,tail_test_len,tail_test_len  * EXON_RECOVER_MATCHING_RATE - 0.0001 );
			ret |= 2;
		}
		else
		{
			int window_start_pos = cover_end - window_size +1;
			int not_too_bad = 1;
			//printf ("TAIL: ROUGHLY UNMAPPED %d in %d!\n", roughly_mapped,tail_test_len );

			while (1)
			{
				int left_match_number = 0;
				int matched_bases_in_window = match_chro_wronglen(read+ window_start_pos, array_index, read_start_pos + window_start_pos + tail_indel, window_size, space_type, &left_match_number, NULL);
				int best_indel_movement = -1000;
				int best_indel_pos = -1;

				if (matched_bases_in_window >= req_match_3end) // this window is not bad enough so that an indel is considered
					window_start_pos++;
				else
				{
					roughly_mapped = match_chro(read+window_start_pos + left_match_number, array_index, read_start_pos + window_start_pos + tail_indel + left_match_number, read_len - window_start_pos - left_match_number , 0, space_type);
					//printf("\nMATCH2: %d in %d\n", roughly_mapped , (int)(0.5+(read_len - window_start_pos -left_match_number )  * EXON_RECOVER_MATCHING_RATE));

					if (roughly_mapped < (int)(0.5+(read_len - window_start_pos -left_match_number )  * EXON_RECOVER_MATCHING_RATE))
					{
						// the quality of this window is too low (at most 1 base is matched). I must consider if it is an indel.
						int indel_movement_i ;
						int best_matched_bases_after_indel = -1;
						not_too_bad = 0;
						for (indel_movement_i = 0; indel_movement_i < 2* indel_tolerance; indel_movement_i ++)
						{
							
							int indel_adjustment = (indel_movement_i+1)/2 * (indel_movement_i %2?1:-1);
							int indel_movement = tail_indel + indel_adjustment;
							int test_length = read_len - window_start_pos  - left_match_number + min(0, indel_adjustment);
							//test_length = min(10, test_length);

							//printf("TAIL: RECOVERING0 : indel_adjustment=%d;   indel_movement=%d;  %d\n",indel_adjustment ,indel_movement,test_length);

							if (test_length < window_size) continue;
							//if (test_length <= 1 + abs(indel_movement)/4) continue;
							if (abs(indel_movement) > indel_tolerance) continue;

							int matched_bases_after_indel = match_chro_support(read + window_start_pos - min(0, indel_adjustment) + left_match_number, array_index, read_start_pos + window_start_pos  + max(0,indel_movement) +left_match_number , test_length,0, space_type, qual_txt +  window_start_pos - min(0, indel_adjustment) + left_match_number , qual_format);

							float test_rate = EXON_INDEL_MATCHING_RATE_TAIL;
							
							if(test_length < 3) test_rate = 1;

							if(best_matched_bases_after_indel <matched_bases_after_indel  && matched_bases_after_indel >= (int)(0.5+test_length * test_rate))
							{
								not_too_bad = 1;
								best_matched_bases_after_indel = matched_bases_after_indel;
								best_indel_movement = indel_movement;
								//printf ("WSP=%d, LFT_MAT=%d\n", window_start_pos, left_match_number);
								best_indel_pos = window_start_pos + left_match_number ;//-1;
								*tail_indel_movement = best_indel_movement ;
							}
						}
					
						if(best_indel_pos<0)
							*tail_indel_pos =  window_start_pos + left_match_number ;
						else
							*tail_indel_pos = best_indel_pos;
						break;
					}else window_start_pos++;
				}
				if (window_start_pos + window_size >= read_len) break;
			}
			if(not_too_bad)
				ret |=2;
			//else
			//	*tail_indel_pos =-1;
		}
	}
	else ret |= 2;


	//printf("\nEXRET=%d\n", ret);
	return ret;
}


int extend_covered_region_backup(gene_value_index_t *array_index, unsigned int read_start_pos, char * read, int read_len, int cover_start, int cover_end, int window_size, short * output_tail, int * head_indel_movement, short * output_head, int * tail_indel_movement, int indel_tolerance, int space_type, int tail_indel)
{
	int ret = 0;
	if (cover_start >= window_size) 
	{
		int head_test_len =  cover_start;
		//printf ("COVERAGE:%d ~ %d\n", cover_start,cover_end);
		int roughly_mapped = match_chro(read, array_index, read_start_pos, head_test_len , 0, space_type);
		if (roughly_mapped >= head_test_len  * 0.75)
		{
			printf ("ROUGHLY MAPPED %d in %d!\n", roughly_mapped,head_test_len );
			ret = 1;
		}
		else
		{

			printf ("ROUGHLY UNMAPPED %d in %d!\n", roughly_mapped,head_test_len );
			
			int window_end_pos = cover_start + window_size-1;
			int current_indel_movements = 0;
			while (1)
			{
			//	printf ("WF=%d\n",window_end_pos);
				if (gvindex_get(array_index, read_start_pos + window_end_pos - window_size + current_indel_movements) == read [window_end_pos - window_size])
					window_end_pos --;
				else
				{
					int search_indel_movement_i;
					int best_window_end_pos = -1;
					int best_indel_movements = 0;
					int best_matching = -1;
					int search_window_end_pos = window_end_pos - window_size+1;
					for (; search_window_end_pos - window_size >=0; search_window_end_pos--)
					{
						for (search_indel_movement_i = 0; search_indel_movement_i < indel_tolerance *2 ; search_indel_movement_i++)
						{
							int search_indel_movement = (search_indel_movement_i+2) / 2*(search_indel_movement_i%2? 1:-1);
							int chromosome_test_pos = read_start_pos+search_window_end_pos - window_size + search_indel_movement + current_indel_movements;
							int read_test_pos = search_window_end_pos - window_size;
							int matched_bases_in_window = match_chro(read+read_test_pos , array_index, chromosome_test_pos, window_size, 0, space_type);
							//printf ("TM: @%d : %d  M=%d\n", read_test_pos, search_indel_movement, matched_bases_in_window);
							if (matched_bases_in_window > max(window_size -2, window_size * 0.8) && best_matching < matched_bases_in_window)
							{
								best_window_end_pos = search_window_end_pos;
								best_indel_movements =  current_indel_movements + search_indel_movement ; 
								best_matching = matched_bases_in_window;
							}
						}
					}
					if (best_window_end_pos >0)
					{
						int toli = 0;
						while(output_head[toli*3])toli++;
						if (toli < 6){
							output_head[toli*3] = best_window_end_pos - window_size ;
							output_head[toli*3+1] = best_indel_movements;
							output_head[toli*3+3] = 0;
						}
						window_end_pos = best_window_end_pos;
						current_indel_movements = best_indel_movements;
						//printf ("ROUGHLY UNMAPPED %d in %d!  SEARCH RANGE %d ~ %d\n", roughly_mapped,head_test_len , cover_start, cover_end);
						//printf ("AT %d INDEL:%d,  BEST MATCHING=%d\n", best_window_end_pos,best_indel_movements , best_matching);
						//printf ("READ=%s , POS=%u\n", read, read_start_pos);
					}
					else{
						//printf("????????????????? @ %d\n", window_end_pos);

						window_end_pos = window_size;
					}
				}

				if (window_end_pos - window_size <= 0) break;
			}
			int toli = 0;
			int remap_begin, remap_end;
			int last_movement = 0;
			int all_mapped = match_chro(read + output_head[0] + window_size, array_index, read_start_pos + output_head[0] +window_size, cover_start - output_head[0]-window_size, 0, space_type);
			int char_skipped = 0;

			while(output_head[toli*3])
			{
				//printf("REC: %d : %d\n", output_head[toli*3], output_head[toli*3+1]);
				if (output_head[toli*3+3])
					remap_begin = output_head[toli*3+3]+ window_size;
				else
					remap_begin = 0;
				int current_indel_movement = output_head[toli*3+1];

				char_skipped  += max(0, current_indel_movement - last_movement);
				remap_end = output_head[toli*3] - max(0, current_indel_movement- last_movement)+ window_size;
				int number_mapped = match_chro(read + remap_begin , array_index, read_start_pos + remap_begin + current_indel_movement , remap_end- remap_begin, 0, space_type);
				all_mapped  += number_mapped;
				//printf("REC: %d (%d) : %d (%d)     REMAP=%d/%d;   ALLMAP=%d/%d CS=%d   RL=%d\n", output_head[toli*3], remap_begin, output_head[toli*3+1], remap_end, number_mapped, remap_end- remap_begin, all_mapped, cover_start - char_skipped, cover_start, read_len);

				toli++;
				last_movement = current_indel_movement ;
			}

			float new_accept_rate = all_mapped*1./(cover_start - char_skipped);
			if(char_skipped < indel_tolerance && new_accept_rate > 0.75)
				ret = 1;
			//printf("NEW_ACCEPT_HEAD=%.4f\n", new_accept_rate);
		}
	}
	else ret = 1;



	if (cover_end <= read_len - window_size ) 
	{
		int tail_test_len =  read_len - cover_end;
		//printf ("COVERAGE:%d ~ %d\n", cover_start,cover_end);
		int current_indel_movements =  tail_indel;
		int roughly_mapped = match_chro(read + read_len - tail_test_len, array_index, read_start_pos+ read_len - tail_test_len + current_indel_movements , tail_test_len , 0, space_type);
		if (roughly_mapped >= tail_test_len  * 0.75)
		{
			//printf ("ROUGHLY MAPPED %d in %d!\n", roughly_mapped,tail_test_len );
			ret |=2;
		}
		else
		{

			//printf ("ROUGHLY UNMAPPED %d in %d!\n", roughly_mapped,tail_test_len );
			
			int window_start_pos = cover_end - window_size +1;
			while (1)
			{
			//	printf ("WF=%d\n",window_end_pos);
				if (gvindex_get(array_index, read_start_pos + window_start_pos + window_size + current_indel_movements) == read [window_start_pos + window_size])
					window_start_pos ++;
				else
				{
					int search_indel_movement_i;
					int best_window_start_pos = -1;
					int best_indel_movements = 0;
					int best_matching = -1;
					int search_window_start_pos = window_start_pos +window_size-1;
					for (; search_window_start_pos + window_size <= read_len; search_window_start_pos++)
					{
						for (search_indel_movement_i = 0; search_indel_movement_i < indel_tolerance *2 ; search_indel_movement_i++)
						{
							int search_indel_movement = (search_indel_movement_i+2) / 2*(search_indel_movement_i%2? 1:-1);
							int chromosome_test_pos = read_start_pos+search_window_start_pos +search_indel_movement + current_indel_movements;
							int read_test_pos = search_window_start_pos;
							int matched_bases_in_window = match_chro(read+read_test_pos , array_index, chromosome_test_pos, window_size, 0, space_type);
							//printf ("TM: @%d : %d  M=%d\n", read_test_pos, search_indel_movement, matched_bases_in_window);
							if (matched_bases_in_window > max(window_size -2, window_size * 0.8) && best_matching < matched_bases_in_window)
							{
								best_window_start_pos = search_window_start_pos;
								best_indel_movements =  current_indel_movements + search_indel_movement ; 
								best_matching = matched_bases_in_window;
							}
						}
					}
					if (best_window_start_pos >0)
					{
						int toli = 0;
						while(output_tail[toli*3])toli++;
						if (toli < 6){
							output_tail[toli*3] = best_window_start_pos ;
							output_tail[toli*3+1] = best_indel_movements;
							output_tail[toli*3+3] = 0;

						}
						window_start_pos = best_window_start_pos;
						current_indel_movements = best_indel_movements;
						//printf ("ROUGHLY UNMAPPED %d in %d!  SEARCH RANGE %d ~ %d\n", roughly_mapped,tail_test_len , cover_start, cover_end);
						//printf ("AT %d INDEL:%d,  BEST MATCHING=%d\n", best_window_start_pos,best_indel_movements , best_matching);
						//printf ("READ=%s , POS=%u\n", read, read_start_pos);
					}
					else{
						//printf("????????????????? @ %d\n", window_end_pos);

						window_start_pos = read_len -window_size;
					}
				}

				if (window_start_pos + window_size >= read_len) break;
			}
			int toli = 0;
			int remap_begin, remap_end;
			int last_movement = tail_indel;
			int all_mapped = 0;
			int char_skipped = 0;

			if(output_tail[0])
				all_mapped=match_chro(read + cover_end , array_index, read_start_pos + tail_indel + cover_end , output_tail[0]-cover_end-max(0,tail_indel-output_tail[1]), 0, space_type);

			while(output_tail[toli*3])
			{
				//printf("REC: %d : %d\n", output_head[toli*3], output_head[toli*3+1]);
				int current_indel_movement = output_tail[toli*3+1];
				if (output_tail[toli*3+3])
					remap_end = output_tail[toli*3+3] -max(0,last_movement-current_indel_movement);
				else
					remap_end = read_len;

				char_skipped  += max(0, last_movement-current_indel_movement);
				remap_begin = output_tail[toli*3];
				//printf("RMBEGIN=%d , RL=%d, TL=%d\n", remap_begin, read_len, remap_end- remap_begin);
				int number_mapped = match_chro(read + remap_begin , array_index, read_start_pos + remap_begin + current_indel_movement , remap_end- remap_begin, 0, space_type);
				all_mapped  += number_mapped;
				//printf("REC: %d (%d) : %d (%d)    REMAP=%d/%d;   ALLMAP=%d/%d; CE=%d ; RL=%d ; LASTM=%d\n", output_tail[toli*3], remap_begin, output_tail[toli*3+1], remap_end, number_mapped, remap_end- remap_begin, all_mapped, read_len - cover_end - char_skipped, cover_end, read_len, last_movement);

				toli++;
				last_movement = current_indel_movement ;
			}

			float new_accept_rate = all_mapped*1./(read_len - cover_end - char_skipped);
			if(char_skipped < indel_tolerance && new_accept_rate > 0.75)
				ret |= 2;
			//printf("NEW_ACCEPT_TAIL=%.4f\n", new_accept_rate);
		}
	}
	else ret |=2;
	return ret;

}

float match_base_quality(gene_value_index_t *array_index, char * read_txt,  unsigned int pos, char * qual_txt, int read_len, int phred_version, int * high_qual_unmatch)
{
	int i;
	float ret =0;
	if(pos < array_index -> start_base_offset || pos + read_len >= array_index -> start_base_offset + array_index -> length){
		printf("WARNING: BASE INDEX OUT OF LIMIT: %u < %u < %u\n%s\n", array_index -> start_base_offset , pos, array_index -> start_base_offset + array_index -> length, read_txt);
	//	exit(-1);
		return 100;
	}
	for(i=0; i<read_len; i++)
	{
		char true_chr = gvindex_get(array_index, pos + i);
		//printf("%c VS %c\n", true_chr, read_txt[i]);
		if (true_chr == read_txt[i])
		{
			if(!qual_txt)
				ret += 1;
			else if(FASTQ_PHRED64 == phred_version)
				ret += (1-get_base_error_prob64(qual_txt[i]));
			else
				ret += (1-get_base_error_prob33(qual_txt[i]));
		}
		else
		{
			if(!qual_txt)
			{
				ret -= 1;
				(*high_qual_unmatch)++;
			}
			else
			{
				float ql ;
				if(FASTQ_PHRED64 == phred_version)
					ql = get_base_error_prob64(qual_txt[i]);
				else
					ql = get_base_error_prob33(qual_txt[i]);
				#ifdef QUALITY_KILL
					#if QUALITY_KILL > 196
						#define QL_MIN 0.999
					#endif
				#endif
	
				#ifndef QL_MIN
					#define QL_MIN 0.2
				#endif

				if( ql < QL_MIN) (*high_qual_unmatch)++;

				ret += ql-1;
			}
		}
	}
	//printf ("SECTION QUAL = %.4f LEN = %d\n", ret, read_len);
	return ret;
}

float final_mapping_quality(gene_value_index_t *array_index, unsigned int pos, char * read_txt, char * qual_txt, char * cigar_txt, int phred_version, int * mismatch)
{
	int cigar_cursor = 0;
	int read_cursor = 0;
	int chromosome_cursor = pos;
	int cigar_length = strlen(cigar_txt);
	int i, x;
	float ret = 0.;

	x= 0;
	while(cigar_cursor < cigar_length)
	{
		if(cigar_txt[cigar_cursor] =='X')
		{
			cigar_cursor++;
			continue;
		}
		if(cigar_txt[cigar_cursor]>='0' && cigar_txt[cigar_cursor]<='9')
			x = x*10+ (cigar_txt[cigar_cursor]-'0');
		else
		{
			if(cigar_txt[cigar_cursor] == 'M' || cigar_txt[cigar_cursor] == 'S') 
			{
				float nret = match_base_quality(array_index, read_txt + read_cursor, chromosome_cursor , (qual_txt && qual_txt[0])?qual_txt + read_cursor:NULL, x, phred_version, mismatch);
				ret += nret;
				chromosome_cursor +=x;
				read_cursor +=x;
			}
			else if(cigar_txt[cigar_cursor] == 'I')
			{
				if(!qual_txt)
				{
					ret += x;
					read_cursor +=x;
				}
				else if(FASTQ_PHRED64 == phred_version)
				{
					for(i = 0; i<x; i++)
					{
						char qchar = qual_txt[read_cursor++];
						ret += (1-get_base_error_prob64(qchar));
					}
				}
				else
				{
					for(i = 0; i<x; i++)
					{
						char qchar = qual_txt[read_cursor++];
						ret += (1-get_base_error_prob33(qchar));
					}
				}
			}
			else if(cigar_txt[cigar_cursor] == 'D' || cigar_txt[cigar_cursor] == 'N')
			{
				chromosome_cursor +=x;
			}
			x= 0;
		}
		cigar_cursor++;
	}

	//printf("S=%.5f, LEN=%d\n", ret , read_cursor);
	return (ret*100) / read_cursor+100.;
}


void print_votes(gene_vote_t * vote, char *index_prefix)
{

	gene_offset_t offsets;
	int i,j;
	char * chrname = NULL;
	unsigned int chrpos = 0;
	

	load_offsets (&offsets, index_prefix);

	locate_gene_position(vote -> max_position, &offsets, &chrname, &chrpos);

	printf("Max votes = %d , Position is %s,%u\n", vote->max_vote, chrname, chrpos );
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
                for(j=0; j< vote->items[i]; j++)
                {
			locate_gene_position(vote -> pos[i][j], &offsets, &chrname, &chrpos);
			printf("\tVote = %d , Position is %s,%u (+%u)\n", vote->votes[i][j] , chrname, chrpos, vote -> pos[i][j]);
		}
	

}


