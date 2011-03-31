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

gene_vote_number_t get_subread_quality(const char * quality_str, const char * read_str, int quality_scale)
{
	gene_vote_number_t ret;
	int i;
/*
	for (i=0; i<16; i++)
		if(read_str[i] == 'N' || read_str[i] == '.')
			return -22.;
*/
	if(quality_scale == QUALITY_SCALE_LOG)
	{
		ret = 0.;
		for(i=0; i<16; i++)
		{
			int ix = quality_str[i] - '@';
			ret += correct_rate_table[ix-1];
		}
	}
	else if(quality_scale == QUALITY_SCALE_LINEAR)
	{
		ret = -.8;
		for(i=0; i<16; i++)
		{
			int ix = quality_str[i] - '@';
			ret += (ix*0.001);
		}
	}
	else ret = -0.1;
	
	return ret;
}

void print_running_log(double finished_rate, double read_per_second, double expected_seconds, unsigned long long int total_reads)
{
        char outbuff[99]; int i;
        snprintf(outbuff, 98,"completed=%0.2f%%; time used=%.1fs; rate=%.1fk reads/s; time left=%.1fs; total=%lluk reads", finished_rate*100, miltime()-begin_ftime, read_per_second/1000 ,expected_seconds, total_reads/1000);
        printf(outbuff);
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
}

void init_allvote(gene_allvote_t* allvote, int expected_len)
{
	allvote -> max_len = expected_len; 
	allvote -> max_positions = (unsigned int *) malloc(sizeof(int)*expected_len);
	allvote -> max_votes = (gene_vote_number_t *) calloc(sizeof(gene_vote_number_t), expected_len);
	allvote -> masks = (unsigned char *) calloc(1, expected_len);
#ifdef REPORT_ALL_THE_BEST
	allvote -> best_records = (gene_best_record_t *) malloc(sizeof(gene_best_record_t)* expected_len);
#endif
	allvote -> is_counterpart = (unsigned char *) malloc(expected_len);
}

void add_allvote(gene_allvote_t* allvote,int qid , int pos, gene_vote_number_t votes, int is_counterpart, char mask)
{
	if(votes >= allvote -> max_votes[qid])
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
		allvote -> max_votes[qid] = (unsigned char)votes;
		allvote -> masks[qid] = mask;
	}
}


void add_allvote_q(gene_allvote_t* allvote,int qid , int pos, gene_vote_number_t votes, int is_counterpart, char mask)
{
	if(votes > allvote -> max_votes[qid])
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
		allvote -> max_votes[qid] = (unsigned char)votes;
		allvote -> masks[qid] = mask;

		//printf("%d\t%.11f\n", votes , good_votes);
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

void mark_votes_array_index(char * read_str, int read_len, gene_vote_t * dest, gene_vote_t * src, gene_value_index_t * my_array_index, int color_space, int indel_tolerance)
{
	int i, j;
	dest -> max_vote = -1;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
	{
		dest -> items[i] = src -> items[i];
		for(j=0; j< src->items[i]; j++)
		{
			unsigned int potential_position = src-> pos[i][j];
			float matchingness_count = match_read(read_str, read_len, potential_position, my_array_index, color_space, indel_tolerance);
	
			dest -> pos[i][j] = potential_position;
			dest -> votes[i][j] = matchingness_count;
			dest -> masks[i][j] = src -> masks[i][j];
			if(matchingness_count >  dest -> max_vote )
			{
				dest -> max_vote  = matchingness_count + 1E-5;
				dest -> max_mask  = src -> masks[i][j];
				dest -> max_position   = potential_position;
			}
		}
	}
}

int select_positions_array(char * read1_str, int read1_len, char * read2_str, int read2_len,  gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2, gehash_data_t * pos_read1, gehash_data_t * pos_read2, unsigned int max_pair_dest, unsigned int min_pair_dest, int min_major, int min_minor,int is_negative_strand, gene_value_index_t * my_array_index, int color_space, int indel_tolerance)
{
	gene_vote_t base_vote_1 , base_vote_2;

	mark_votes_array_index(read1_str, read1_len, &base_vote_1, vote_read1, my_array_index, color_space, indel_tolerance);
	mark_votes_array_index(read2_str, read2_len, &base_vote_2, vote_read2, my_array_index, color_space, indel_tolerance);

	return select_positions(&base_vote_1, &base_vote_2, numvote_read1, numvote_read2, pos_read1, pos_read2, max_pair_dest,  min_pair_dest,  min_major,  min_minor, is_negative_strand);
}

#define ANCHORS_NUMBER 259

int select_positions(gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2, gehash_data_t * pos_read1, gehash_data_t * pos_read2, unsigned int max_pair_dest, unsigned int min_pair_dest, int min_major, int min_minor,int is_negative_strand)
{

	int k, i, j, anchors = 0;
	gehash_data_t anchors_position [ANCHORS_NUMBER];
	unsigned char anchor_read [ANCHORS_NUMBER];
	gene_vote_number_t anchor_votes = max(vote_read1->max_vote, vote_read2->max_vote);
	gene_vote_number_t anchors_votes [ANCHORS_NUMBER];
	gene_vote_number_t minor_votes [ANCHORS_NUMBER];
	gehash_data_t minor_position [ANCHORS_NUMBER];

	if(anchor_votes < min_major*22-21.99999)
		return 0;

	anchor_votes = -0.00001+22.*(int)((anchor_votes+0.000001)/22);

	for (k=0; k<2; k++)
	{
		gene_vote_t * current_vote = k?vote_read2:vote_read1;
		if (current_vote->max_vote < anchor_votes)
			continue;

		for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
			for(j=0; j< current_vote->items[i]; j++)
			{
				if(current_vote->votes[i][j] >= anchor_votes)
				{
					anchors_position[anchors] = current_vote->pos[i][j];
					anchors_votes   [anchors] = current_vote->votes[i][j];
					anchor_read	[anchors] = k;
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

	for (k=0; k<2; k++)
	{
		gene_vote_t * current_vote = k?vote_read2:vote_read1;
		gene_vote_t * current_vote2 = k?vote_read1:vote_read2;

		if (current_vote2->max_vote < anchor_votes)
			continue;

		for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
			for(j=0; j< current_vote->items[i]; j++)
			{
				if(current_vote->votes[i][j] >= min_minor*22-21.99999)
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
							continue;
						if (	(minor_votes[l] < current_vote->votes[i][j]) && 
							(abdist <= max_pair_dest) &&
							(abdist >= min_pair_dest) 
						   )
						{
							minor_votes[l] = current_vote->votes[i][j];
							minor_position[l] = current_vote->pos[i][j];
						}
					}
				}
			}

	}

	gene_vote_number_t selection_minor_vote = 0;
	gene_vote_number_t selection_major_vote = 0;
	gehash_data_t selection_minor_pos = 0;
	gehash_data_t selection_major_pos = 0;
	unsigned char selection_major_read_no = 0;

	for (k=0; k<anchors; k++)
	{
		if(minor_votes[k]<0.00001)
			continue;

		if ((minor_votes[k]+anchors_votes[k]) > (selection_minor_vote+selection_major_vote))
		{
			selection_minor_pos = minor_position[k];
			selection_major_pos = anchors_position[k];
			selection_minor_vote = minor_votes[k];
			selection_major_vote = anchors_votes[k];
			selection_major_read_no = anchor_read[k];
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
		}
		else
		{
			*numvote_read1 = selection_major_vote;
			*numvote_read2 = selection_minor_vote;
			*pos_read1 = selection_major_pos;
			*pos_read2 = selection_minor_pos;
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

float match_read(const char read_str[], int read_len, unsigned int potential_position,  gene_value_index_t * my_array_index, int space_type, int indel_tolerance)
{
	int ret = 0;
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
			if(i % 3 == 2)
			{
				int b2_match = (read_matchingness[bias+indel_tolerance][i -2]  && read_matchingness[bias+indel_tolerance][i -1] && read_matchingness[bias+indel_tolerance][i]);
				if(b2_match)
				{
					read_matchingness[bias+indel_tolerance][i]  = 2;
					read_matchingness[bias+indel_tolerance][i-1]  = 2;
					read_matchingness[bias+indel_tolerance][i-2]  = 2;
				}
			}
		}
	}

	for(i=0; i<read_len; i+=8)
	{
	
		int j;
		int movement = 0; 
		int max_matchness = -1;
		for (j=-indel_tolerance; j<=indel_tolerance; j++)
		{
			int m = front8(read_matchingness[j],i);
			if(m > max_matchness)
			{
				max_matchness = m;
				movement = j;
			}

		}
		ret += max_matchness; //front8(read_matchingness[movement], i);  // (front4(read_matchingness[movement], i)==4?4:0);
	}

	return 21.99 * (int)(ret / 14) ;
}


void final_matchingness_scoring(const char read_str[], const char quality_str[], int read_len, gene_vote_t * vote, gehash_data_t * max_position, gene_vote_number_t * max_vote, char *max_mask, gene_value_index_t * my_array_index, int space_type, int indel_tolerance)
{
	int i, j;
	float max_matching = 0;
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			unsigned int potential_position = vote -> pos[i][j];
			float matchingness_count = match_read(read_str, read_len, potential_position, my_array_index, space_type, indel_tolerance);
			if(matchingness_count > max_matching)
			{
				max_matching = matchingness_count+1E-5;
		//		printf ("MAX MATCHING = %f / %s\n", max_matching, read_str);
				*max_vote = matchingness_count; //vote -> votes[i][j] + matchingness_count * 3;
				*max_position = potential_position;
				*max_mask = vote -> masks[i][j];
			}
		}
}

