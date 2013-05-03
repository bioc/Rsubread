#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifndef MAKE_STANDALONE
  #include <R.h>
#endif

#include <zlib.h>
#include <pthread.h>
#include <getopt.h>
#include "subread.h"
#include "sambam-file.h"
#include "hashtable.h"
#include "HelperFunctions.h"

/********************************************************************/
/********************************************************************/
/********************************************************************/
//  NEW FUNCTION FOR MULTI-THREADING
/********************************************************************/
/********************************************************************/
/********************************************************************/
#define FEATURE_NAME_LENGTH 48
#define CHROMOSOME_NAME_LENGTH 48
#define MAX_LINE_LENGTH 100000
#define FILE_TYPE_RSUBREAD 10
#define FILE_TYPE_GTF 100

typedef struct
{
	unsigned char feature_name[FEATURE_NAME_LENGTH];
	unsigned char chro[CHROMOSOME_NAME_LENGTH];
	
	unsigned int start;
	unsigned int end;

	unsigned int sorted_order;

	char is_negative_strand;
} fc_feature_info_t;

typedef struct
{
	unsigned short thread_id;
	char * line_buffer1;
	char * line_buffer2;
	unsigned long long int nreads_mapped_to_exon;
	unsigned short current_read_length1;
	unsigned short current_read_length2;
	unsigned int * count_table;
	pthread_t thread_object;

	char * input_buffer;
	unsigned int input_buffer_remainder;
	unsigned int input_buffer_write_ptr;	
	pthread_spinlock_t input_buffer_lock;
} fc_thread_thread_context_t;


typedef struct
{
	int is_gene_level;
	int is_paired_end_data;
	int is_multi_overlap_allowed;
	int is_strand_checked;
	int min_paired_end_distance;
	int max_paired_end_distance;
	int read_length;
	int line_length;

	unsigned short thread_number;
	fc_thread_thread_context_t * thread_contexts;
	int is_all_finished;
	unsigned int input_buffer_max_size;

	HashTable * gene_name_table;	// gene_name -> gene_internal_number
	unsigned char ** gene_name_array;	// gene_name -> gene_internal_number
	int exontable_nchrs;
	int exontable_exons;
	int * exontable_geneid;
	char * exontable_strand;
	char ** exontable_chr;
	long * exontable_start;
	long * exontable_stop;

	char ** exontable_anno_chrs;
	long * exontable_anno_chr_heads;

	FILE * SAM_output_fp;
	
} fc_thread_global_context_t;

unsigned int tick_time = 1000;


// This function loads annotations from the file.
// It returns the number of featres loaded, or -1 if something is wrong. 
// Memory will be allowcated in this function. The pointer is saved in *loaded_features.
// The invoker must release the memory itself.

int load_feature_info(const char * annotation_file, int file_type, fc_feature_info_t ** loaded_features)
{
	unsigned int features = 0, xk1 = 0;
	char * file_line = malloc(MAX_LINE_LENGTH);
	FILE * fp = fopen(annotation_file,"r"); 
	if(!fp) return -1;
	if(file_type == FILE_TYPE_RSUBREAD)
		fgets(file_line, MAX_LINE_LENGTH, fp);

	while(1)
	{
		char * fgets_ret = fgets(file_line, MAX_LINE_LENGTH, fp);
		if(!fgets_ret) break;
		if(file_type == FILE_TYPE_GTF)
		{
			char * token_temp;
			strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp);
			char * feature_type = strtok_r(NULL,"\t", &token_temp);
			if(strcmp(feature_type, "exon")==0)
				features++;
		}
		else
			features++;
	}

	fseek(fp,0,SEEK_SET);
	if(file_type == FILE_TYPE_RSUBREAD)
		fgets(file_line, MAX_LINE_LENGTH, fp);

	fc_feature_info_t * ret_features = malloc(sizeof(fc_feature_info_t) * features);

	while(xk1 < features)
	{
		fgets(file_line, MAX_LINE_LENGTH, fp);
		char * token_temp;

		if(file_type == FILE_TYPE_RSUBREAD)
		{
			char * feature_name = strtok_r(file_line,"\t",&token_temp);
			if(!feature_name) break;
			strncpy((char *)ret_features[xk1].feature_name, (char *)feature_name, FEATURE_NAME_LENGTH);
			char * seq_name = strtok_r(NULL,"\t", &token_temp);
			if(!seq_name) break;
			strncpy((char *)ret_features[xk1].chro, (char *)seq_name, CHROMOSOME_NAME_LENGTH);
			ret_features[xk1].start = atoi(strtok_r(NULL,"\t", &token_temp));// start 
			ret_features[xk1].end = atoi(strtok_r(NULL,"\t", &token_temp));//end 
			char * strand_str = strtok_r(NULL,"\t", &token_temp); 
			if(strand_str == NULL)
				ret_features[xk1].is_negative_strand = 0;
			else
				ret_features[xk1].is_negative_strand = ('-' ==strand_str[0]);
			ret_features[xk1].sorted_order = xk1;
			xk1++;
		}
		else if(file_type == FILE_TYPE_GTF)
		{
			sprintf((char *)ret_features[xk1].feature_name, "UK%06u", xk1 + 1);
			char * seq_name = strtok_r(file_line,"\t",&token_temp);
			strncpy((char *)ret_features[xk1].chro, (char *)seq_name, CHROMOSOME_NAME_LENGTH);

			strtok_r(NULL,"\t", &token_temp);// source
			char * feature_type = strtok_r(NULL,"\t", &token_temp);// feature_type
			if(strcmp(feature_type, "exon")==0)
			{
				ret_features[xk1].start = atoi(strtok_r(NULL,"\t", &token_temp));// start 
				ret_features[xk1].end = atoi(strtok_r(NULL,"\t", &token_temp));//end 
				strtok_r(NULL,"\t", &token_temp);// score 
				ret_features[xk1].is_negative_strand = ('-' == (strtok_r(NULL,"\t", &token_temp)[0]));//strand 
				ret_features[xk1].sorted_order = xk1;
				strtok_r(NULL,"\t",&token_temp);	// "frame"
				char * extra_attrs = strtok_r(NULL,"\t",&token_temp);	// name_1 "val1"; name_2 "val2"; ... 
				if(extra_attrs && (strlen(extra_attrs)>6))
				{
					#define FC_GTF_GENE_ID_NAME	"gene_id"
					int attr_state = 0;	// 0: name; 1:value
					int name_start = 0;
					int val_start = 0;
					int is_gene_id = 0;
					int xk2, exch;
					for(xk2 = 0; exch = extra_attrs[xk2]; xk2++)
					{
						if(exch == ' ' && attr_state == 0)
							continue;
						else if(exch == '\"')
						{
							if(xk2<1) break;
							if(attr_state == 0){
								val_start = xk2 + 1;
								extra_attrs[xk2 - 1] = 0;
								if(strcmp(extra_attrs + name_start, FC_GTF_GENE_ID_NAME)==0)
									is_gene_id = 1;
							}
							else if(attr_state == 1)
							{
								if(is_gene_id)
								{
									extra_attrs[xk2]=0;
									strncpy((char *)ret_features[xk1].feature_name, extra_attrs + val_start, FEATURE_NAME_LENGTH);
									//printf("N=%s\n", ret_features[xk1].feature_name);
									is_gene_id = 2;
									break;
								}
								is_gene_id = 0 ;
								name_start = xk2 + 2;	// "; NAME2
							}
							
							attr_state = !attr_state;
							continue;
						}
						
					}
				}
				xk1++;
			}
		}
	}
	fclose(fp);

	(*loaded_features) = ret_features;
	SUBREADprintf("There are %d features loaded from the annotation file.\n", features);
	return features;
}

int find_or_insert_gene_name(fc_thread_global_context_t * global_context, unsigned char * feature_name)
{
	HashTable * genetable = global_context -> gene_name_table;

	long long int gene_number = HashTableGet(genetable, feature_name) - NULL;
	if(gene_number>0)
		return gene_number-1;
	else
	{
		gene_number = genetable -> numOfElements; 
		HashTablePut(genetable, feature_name, NULL+gene_number+1);
		global_context -> gene_name_array[gene_number] = feature_name;
			// real memory space of feature_name is in the "loaded_features" data structure.
			// now we only save its pointer.

		return gene_number;
	}
}

int fc_strcmp(const void * s1, const void * s2)
{
	return strcmp((char*)s1, (char*)s2);
}


void sort_feature_info(fc_thread_global_context_t * global_context, unsigned int features, fc_feature_info_t * loaded_features, char *** sorted_chr_names, int ** sorted_entrezid, long ** sorted_start, long ** sorted_end, unsigned char ** sorted_strand, int * chro_number, char *** anno_chrs, long ** anno_chr_head)
{
	unsigned int chro_pnt;
	unsigned char * is_chro_name_explored = malloc(features);
	int * ret_entrez = malloc(sizeof(int) * features);
	long * ret_start = malloc(sizeof(long) * features);
	long * ret_end = malloc(sizeof(long) * features);
	unsigned char * ret_strand = malloc(features);
	char ** ret_char_name = malloc(sizeof(void *) * features);
	fc_feature_info_t ** old_info_ptr = malloc(sizeof(void *) * features);
	(*anno_chrs) = malloc(sizeof(void *) * 500);
	(*anno_chr_head) = malloc(sizeof(long) * 500);

	global_context -> gene_name_array = malloc(sizeof(char *) * features);	// there should be much less identical names.
	global_context -> gene_name_table = HashTableCreate(5000);
	HashTableSetHashFunction(global_context -> gene_name_table, HashTableStringHashFunction);
	HashTableSetKeyComparisonFunction(global_context -> gene_name_table, fc_strcmp);


	unsigned int this_chro_start = 0;
	unsigned int this_chro_end = 0;
	unsigned int ret_array_pointer = 0;

	memset(is_chro_name_explored, 0 , features);

	(*sorted_chr_names) = ret_char_name;
	(*sorted_entrezid) = ret_entrez;
	(*sorted_start) = ret_start;
	(*sorted_end) = ret_end;
	(*sorted_strand) = ret_strand;
	(*chro_number) = 0;

	for(chro_pnt=0; chro_pnt < features; chro_pnt++)
	{
		if(!is_chro_name_explored[chro_pnt])
		{
			unsigned char * current_chro = loaded_features[chro_pnt].chro;
			unsigned int xk1;
			// aggregate all features belonging to this chromosome.
			this_chro_start = ret_array_pointer;
			(*anno_chrs) [(*chro_number)] = (char *)current_chro;
			(*anno_chr_head) [(*chro_number)] = this_chro_start;
			(*chro_number) ++;
			for(xk1 = chro_pnt; xk1<features; xk1++)
			{
				if(strcmp((char *)current_chro, (char *)loaded_features[xk1].chro)==0)
				{
					is_chro_name_explored[xk1]=1;
					//ret_entrez[ret_array_pointer] = atoi((char *)loaded_features[xk1].feature_name);
					ret_entrez[ret_array_pointer] = find_or_insert_gene_name(global_context, loaded_features[xk1].feature_name);
					ret_start[ret_array_pointer] = loaded_features[xk1].start;
					ret_end[ret_array_pointer] = loaded_features[xk1].end;
					ret_strand[ret_array_pointer] = loaded_features[xk1].is_negative_strand;
					ret_char_name[ret_array_pointer] = (char *)loaded_features[xk1].chro;
					old_info_ptr[ret_array_pointer] = &loaded_features[xk1];
					ret_array_pointer++;
				}
			}
			this_chro_end = ret_array_pointer;
			(*anno_chr_head) [(*chro_number)] = this_chro_end;

			// sorting the features.
			unsigned int sort_i, sort_j;
			for(sort_i = this_chro_start; sort_i< this_chro_end; sort_i++)
			{
				unsigned int min_j = sort_i;
				unsigned int min_start = ret_start[min_j];
				for(sort_j = sort_i+1; sort_j < this_chro_end; sort_j++)
				{
					// put the minimun start to *sort_i 
					if(ret_start[sort_j] < min_start)
					{
						min_j = sort_j;
						min_start = ret_start[sort_j];
					}
				}
				if(min_j != sort_i)
				{
					unsigned int tmp;
					void * tmp_ptr;
					long tmpl;

					tmpl = ret_start[sort_i];
					ret_start[sort_i] = ret_start[min_j];
					ret_start[min_j] = tmpl;

					tmpl = ret_end[sort_i];
					ret_end[sort_i] = ret_end[min_j];
					ret_end[min_j] = tmpl;

					tmp = ret_strand[sort_i];
					ret_strand[sort_i] = ret_strand[min_j];
					ret_strand[min_j] = tmp;

					tmp = ret_entrez[sort_i];
					ret_entrez[sort_i] = ret_entrez[min_j];
					ret_entrez[min_j] = tmp;

					tmp_ptr = old_info_ptr[sort_i];
					old_info_ptr[sort_i] = old_info_ptr[min_j];
					old_info_ptr[min_j] = tmp_ptr;
				}
				old_info_ptr[sort_i]->sorted_order = sort_i;
			}
		}
	}
	SUBREADprintf("The %u features are sorted.\n", ret_array_pointer);
	free(old_info_ptr);
}

#define MAX_HIT_NUMBER 200

void process_line_buffer(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context)
{

	char * mate_chr = NULL, * read_chr, * line = thread_context -> line_buffer1, *tmp_tok_ptr, *CIGAR_str ;
	long read_pos, mate_pos = 0, fragment_length = 0, search_start, search_end, pos_leftmost;
	int is_chro_found, i, j, nhits, flag_overlap, prev_gid, alignment_masks;
	long hits_indices[MAX_HIT_NUMBER];
	char is_gene_explored[MAX_HIT_NUMBER];

	//printf("L=%s (%d)\n", line, strlen(line));

	// process the current read or read pair
	char * read_name = strtok_r(line,"\t", &tmp_tok_ptr);	// read name

	alignment_masks = atoi(strtok_r(NULL,"\t", &tmp_tok_ptr));

	read_chr = strtok_r(NULL,"\t", &tmp_tok_ptr);

	//skip the read if unmapped (its mate will be skipped as well if paired-end)
	if(*read_chr == '*' || (alignment_masks & SAM_FLAG_UNMAPPED) == 4){
		if(global_context -> SAM_output_fp)
			fprintf(global_context -> SAM_output_fp,"%s\tUNMAPPED\n", read_name);
		return;	// do nothing if a read is unmapped, or the first read in a pair of reads is unmapped.
	}
	int is_first_read_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
	int is_second_read_negative_strand = (alignment_masks & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)?1:0; 


	//get mapping location of the read (it could be the first read in a pair)
	read_pos = atoi(strtok_r(NULL,"\t", &tmp_tok_ptr));

	strtok_r(NULL,"\t", &tmp_tok_ptr);	// mapping quality
	CIGAR_str = strtok_r(NULL,"\t", &tmp_tok_ptr);	// CIGAR string

	//remove reads which are not properly paired if paired-end reads are used (on different chromsomes or paired-end distance is too big or too small)
	if(global_context -> is_paired_end_data){
		mate_chr = strtok_r(NULL,"\t", &tmp_tok_ptr); //get chr which the mate read is mapped to
		mate_pos = atoi(strtok_r(NULL,"\t", &tmp_tok_ptr));
		char * frag_len_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
		fragment_length = abs(atoi(frag_len_str)); //get the fragment length

		//printf("MM=%s ; FL=%d; LL=%d ; CK=%d   %d~%d\n", mate_chr, fragment_length,  thread_context->current_read_length1, global_context->is_strand_checked, global_context -> min_paired_end_distance, global_context -> max_paired_end_distance);
		if(strcmp(mate_chr,"=") != 0 || fragment_length > (global_context -> max_paired_end_distance + thread_context->current_read_length1 -1) || fragment_length < (global_context -> min_paired_end_distance + thread_context->current_read_length1 - 1) || (global_context->is_strand_checked && (is_first_read_negative_strand==is_second_read_negative_strand))){
		//	printf("PMQ\n");
			//the two reads are not properly paired and are skipped
			if(global_context -> SAM_output_fp)
				fprintf(global_context -> SAM_output_fp,"%s\tPAIR_CONDITION\n", read_name);
			return;
		}
	} //end if(isPE==1)

	is_chro_found = 0;
	for(i=0;i<global_context -> exontable_nchrs;i++)
		if(strcmp(read_chr,global_context -> exontable_anno_chrs[i])==0){
			//get chr to which the current read or fragment is mapped and also the searching range
			is_chro_found = 1;
			search_start = global_context -> exontable_anno_chr_heads[i];
			search_end = global_context -> exontable_anno_chr_heads[i+1] - 1;
			break;
		}


	//printf("FX=%d\n", is_chro_found);
	if(is_chro_found)
	{
		nhits = 0;

		if(global_context -> is_paired_end_data)
		{
			if(read_pos < mate_pos)
				pos_leftmost = read_pos;
			else
				pos_leftmost = mate_pos;

			for(i=search_start;i<=search_end;i++){
				if (global_context -> exontable_start[i] > (pos_leftmost + fragment_length - 1)) break;
				if (global_context -> exontable_stop[i] >= pos_leftmost && ((!global_context->is_strand_checked)||is_first_read_negative_strand == global_context -> exontable_strand[i])){
					hits_indices[nhits] = i;
					is_gene_explored[nhits] = 0;
					nhits++;
					if(nhits>=MAX_HIT_NUMBER) break;
				} 
			}
		}
		else{
			int cigar_section_id, cigar_sections;
			unsigned int Staring_Points[6];
			unsigned short Section_Lengths[6];

			cigar_sections = RSubread_parse_CIGAR_string(CIGAR_str, Staring_Points, Section_Lengths);
			for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
			{
				long section_begin_pos = read_pos + Staring_Points[cigar_section_id];
				long section_length = Section_Lengths[cigar_section_id];
		//		printf("CIGAR_str=%s; cigar_sections=%d; base=%ld; pos[%d]=%u ; len[%d]=%d\n", CIGAR_str, cigar_sections, read_pos, cigar_section_id, Staring_Points[cigar_section_id],cigar_section_id, Section_Lengths[cigar_section_id]);
				for(i=search_start;i<=search_end;i++){
					if (global_context -> exontable_start[i] > (section_begin_pos + section_length -1)) break;
					if (global_context -> exontable_stop[i] >= section_begin_pos && ((!global_context->is_strand_checked)||is_first_read_negative_strand == global_context -> exontable_strand[i])){
						hits_indices[nhits] = i;
						is_gene_explored[nhits] = 0;
		//				printf("HIT: %s, %ld ~ %ld\n",global_context -> exontable_chr[i] , global_context -> exontable_start[i] ,global_context -> exontable_stop[i] );
						nhits++;
						if(nhits>=MAX_HIT_NUMBER) break;
					} 
				}
				if(nhits>=MAX_HIT_NUMBER) break;
			}
		}

		//printf("HHITS=%d\n", nhits);
		if (nhits > 0){
			if (nhits == 1){
				thread_context->nreads_mapped_to_exon++;
				thread_context->count_table[hits_indices[0]]++; 
				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tACCEPTED\t%d\n", read_name, nhits);
			}
			else { // nhits greater than 1	    		
				if (global_context -> is_multi_overlap_allowed){
					if (global_context -> is_gene_level){
						//gives only one vote to each gene.
						int xk1;

						for(xk1=0; xk1<nhits; xk1++)
						{
							if(is_gene_explored[xk1])continue;
							prev_gid = global_context -> exontable_geneid[hits_indices[xk1]];
							thread_context->count_table[hits_indices[xk1]]++;
							int xk2;
							for(xk2=xk1; xk2<nhits; xk2++)
								if(prev_gid == global_context -> exontable_geneid[hits_indices[xk2]]) is_gene_explored[xk2]=1;
						}
					}
					else
						for (j=0;j<nhits;j++){
							thread_context->count_table[hits_indices[j]]++;
						}
					if(global_context -> SAM_output_fp)
						fprintf(global_context -> SAM_output_fp,"%s\tACCEPTED\t%d\n", read_name, nhits);
					thread_context->nreads_mapped_to_exon++;
				}
				else { // multi-overlap is not allowed		
					if (global_context -> is_gene_level){
						prev_gid = global_context -> exontable_geneid[hits_indices[0]];
						flag_overlap = 0;
						for (j=1;j<nhits;j++){
							if (global_context -> exontable_geneid[hits_indices[j]] != prev_gid){
								flag_overlap = 1;
								break;
							} 
						}

						if (flag_overlap){
							if(global_context -> SAM_output_fp)
								fprintf(global_context -> SAM_output_fp,"%s\tOVERLAPPED_GENES\n", read_name);
						}else{ //overlap multiple exons from a single gene
							thread_context->nreads_mapped_to_exon++;
							// only one vote for a gene. the vote is given to the first exon.
							thread_context->count_table[hits_indices[0]]++;
							if(global_context -> SAM_output_fp)
								fprintf(global_context -> SAM_output_fp,"%s\tACCEPTED_SAME_GENE\t%d\n", read_name, nhits);
						}
					}
					else
						if(global_context -> SAM_output_fp)
							fprintf(global_context -> SAM_output_fp,"%s\tOVERLAPPED_EXONS\n", read_name);
				}
			}
		}
		else	if(global_context -> SAM_output_fp)
				fprintf(global_context -> SAM_output_fp,"%s\tNO_FEATURE\n", read_name);
	}
	else
		if(global_context -> SAM_output_fp)
			fprintf(global_context -> SAM_output_fp,"%s\tUNKNOWN_CHRO\t%s\n", read_name, read_chr);

	// Note that we actually make NO use of the second read in a pair.
	// All information we need is in the first line.
}

void * feature_count_worker(void * vargs)
{
	void ** args = (void **) vargs;

	fc_thread_global_context_t * global_context = args[0];
	fc_thread_thread_context_t * thread_context = args[1];

	free(vargs);

	while (1)
	{
		while(1)
		{
			int is_retrieved = 0;
			pthread_spin_lock(&thread_context->input_buffer_lock);
			if(thread_context->input_buffer_remainder)
			{
				int is_second_read;
				unsigned int buffer_read_bytes ;
				unsigned int buffer_read_ptr;
				if(thread_context->input_buffer_remainder <= thread_context->input_buffer_write_ptr)
					buffer_read_ptr = thread_context->input_buffer_write_ptr - thread_context->input_buffer_remainder; 
				else
					buffer_read_ptr = thread_context->input_buffer_write_ptr + global_context->input_buffer_max_size - thread_context->input_buffer_remainder;

				for(is_second_read = 0; is_second_read < (global_context->is_paired_end_data ? 2:1); is_second_read++)
				{
					char * curr_line_buff = is_second_read?thread_context -> line_buffer2:thread_context -> line_buffer1;
					//printf("R=%u; WPTR=%u ;RPTR=%u\n", thread_context->input_buffer_remainder, thread_context->input_buffer_write_ptr, buffer_read_ptr);
					//if(buffer_read_ptr % 7 == 0)
					//	fflush(stdout);
					
					for(buffer_read_bytes=0; ; buffer_read_bytes++)
					{
						char nch =  thread_context->input_buffer[buffer_read_ptr ++];
						curr_line_buff[buffer_read_bytes] = nch;
						if(buffer_read_ptr == global_context->input_buffer_max_size)
							buffer_read_ptr = 0; 
						if(nch=='\n'){
							curr_line_buff[buffer_read_bytes+1]=0;
							break;
						}
					}

					//if(buffer_read_bytes + 1 > thread_context->input_buffer_remainder)
					//	(*(int*)NULL) = 1;
					thread_context->input_buffer_remainder -= buffer_read_bytes + 1;
				}
				is_retrieved = 1;

			}

			pthread_spin_unlock(&thread_context->input_buffer_lock);
			if(global_context->is_all_finished && !is_retrieved) return NULL;

			if(is_retrieved) break;
			else
				usleep(tick_time);
		}


		thread_context -> current_read_length1 = global_context -> read_length;
		thread_context -> current_read_length2 = global_context -> read_length;

		process_line_buffer(global_context, thread_context);

	}
}

void fc_thread_merge_results(fc_thread_global_context_t * global_context, int * nreads , int *nreads_mapped_to_exon)
{
	int xk1, xk2;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		for(xk2=0; xk2<global_context -> exontable_exons; xk2++)
		{
			nreads[xk2]+=global_context -> thread_contexts[xk1].count_table[xk2];
		}
		SUBREADprintf("The %d-th thread processed %llu reads\n", xk1, global_context -> thread_contexts[xk1].nreads_mapped_to_exon);
		(*nreads_mapped_to_exon) += global_context -> thread_contexts[xk1].nreads_mapped_to_exon;
	}
}

void fc_thread_init_global_context(fc_thread_global_context_t * global_context, unsigned int buffer_size, unsigned short threads, int line_length , int is_PE_data, int min_pe_dist, int max_pe_dist, int is_gene_level, int is_overlap_allowed, int is_strand_checked, char * output_fname, int is_sam_out)
{

	global_context -> input_buffer_max_size = buffer_size;

	global_context -> is_multi_overlap_allowed = is_overlap_allowed;
	global_context -> is_paired_end_data = is_PE_data;
	global_context -> is_gene_level = is_gene_level;
	global_context -> is_strand_checked = is_strand_checked;


	if(is_sam_out)
	{
		char tmp_fname[350];
		sprintf(tmp_fname, "%s.reads",output_fname);
		global_context -> SAM_output_fp = fopen(tmp_fname, "w");
	}
	else
		global_context -> SAM_output_fp = NULL;

	global_context -> min_paired_end_distance = min_pe_dist;
	global_context -> max_paired_end_distance = max_pe_dist;
	global_context -> thread_number = threads;
	global_context -> line_length = line_length;

}
int fc_thread_start_threads(fc_thread_global_context_t * global_context, int et_nchrs, int et_exons, int * et_geneid, char ** et_chr, long * et_start, long * et_stop, unsigned char * et_strand, char ** et_anno_chrs, long * et_anno_chr_heads, int read_length)
{
	int xk1;

	global_context -> read_length = read_length;

	global_context -> exontable_nchrs = et_nchrs;
	global_context -> exontable_exons = et_exons;
	global_context -> exontable_geneid = et_geneid;
	global_context -> exontable_chr = et_chr;
	global_context -> exontable_start = et_start;
	global_context -> exontable_stop = et_stop;
	global_context -> exontable_strand = (char *)et_strand;
	global_context -> exontable_anno_chrs = et_anno_chrs;
	global_context -> exontable_anno_chr_heads = et_anno_chr_heads;
	global_context -> is_all_finished = 0;
	global_context -> thread_contexts = malloc(sizeof(fc_thread_thread_context_t) * global_context -> thread_number);
	for(xk1=0; xk1<global_context -> thread_number; xk1++)
	{
		pthread_spin_init(&global_context->thread_contexts[xk1].input_buffer_lock, PTHREAD_PROCESS_PRIVATE);
		global_context -> thread_contexts[xk1].input_buffer_remainder = 0;
		global_context -> thread_contexts[xk1].input_buffer_write_ptr = 0;
		global_context -> thread_contexts[xk1].input_buffer = malloc(global_context -> input_buffer_max_size);
		global_context -> thread_contexts[xk1].thread_id = xk1;
		global_context -> thread_contexts[xk1].count_table = calloc(sizeof(unsigned int), et_exons);
		global_context -> thread_contexts[xk1].nreads_mapped_to_exon = 0;
		global_context -> thread_contexts[xk1].line_buffer1 = malloc(global_context -> line_length + 2);
		global_context -> thread_contexts[xk1].line_buffer2 = malloc(global_context -> line_length + 2);
		if(!global_context ->  thread_contexts[xk1].count_table) return 1;
		void ** thread_args = malloc(sizeof(void *)*2);
		thread_args[0] = global_context;
		thread_args[1] = & global_context -> thread_contexts[xk1];
		pthread_create(&global_context -> thread_contexts[xk1].thread_object, NULL, feature_count_worker, thread_args);
	}

	return 0;
}

void fc_thread_destroy_global_context(fc_thread_global_context_t * global_context)
{
	int xk1;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		free(global_context -> thread_contexts[xk1].count_table);	
		free(global_context -> thread_contexts[xk1].line_buffer1);	
		free(global_context -> thread_contexts[xk1].line_buffer2);	
		free(global_context -> thread_contexts[xk1].input_buffer);
		pthread_spin_destroy(&global_context -> thread_contexts[xk1].input_buffer_lock);
	}
	HashTableDestroy(global_context -> gene_name_table);
	free(global_context -> thread_contexts);
	free(global_context -> gene_name_array);
	if(global_context -> SAM_output_fp) fclose(global_context -> SAM_output_fp);
}
void fc_thread_wait_threads(fc_thread_global_context_t * global_context)
{
	int xk1;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
		pthread_join(global_context -> thread_contexts[xk1].thread_object, NULL);
}


void fc_write_final_gene_results(fc_thread_global_context_t * global_context, int * et_geneid, long * et_start, long * et_stop, const char * out_file, int features, int * nreads, fc_feature_info_t * loaded_features)
{
	int xk1;
	int genes = global_context -> gene_name_table -> numOfElements;
	int * gene_nreads = calloc(sizeof(int) , genes);
	int * gene_lengths = calloc(sizeof(int) , genes);

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1];
		gene_nreads[gene_id] += nreads[xk1];
		gene_lengths[gene_id] += et_stop[xk1] - et_start[xk1] + 1;
	}
	FILE * fp_out = fopen(out_file,"w");
	if(!fp_out){
		SUBREADprintf("Failed to create file %s\n", out_file);
		return;
	}

	fprintf(fp_out,"geneid\tlength\tnreads\n");
	for(xk1 = 0 ; xk1 < genes; xk1++)
	{
		unsigned char * gene_symbol = global_context -> gene_name_array [xk1];
		//printf("%s\tXK1=%d\tGENEs=%d\n", gene_symbol, xk1, genes);
		fprintf(fp_out,"%s\t%d\t%d\n", gene_symbol , gene_lengths[xk1], gene_nreads[xk1]);
	}

	free(gene_nreads);
	free(gene_lengths);
	fclose(fp_out);
}
void fc_write_final_results(fc_thread_global_context_t * global_context, const char * out_file, int features, int * nreads, fc_feature_info_t * loaded_features)
{
	/* save the results */
	FILE * fp_out;
	int i;
	fp_out = fopen(out_file,"w");
	if(!fp_out){
		SUBREADprintf("Failed to create file %s\n", out_file);
			return;
		}
	fprintf(fp_out,"geneid\tchr\tstart\tend\tstrand\tnreads\n");
	for(i=0;i<features;i++)
	{
		fprintf(fp_out,"%s\t%s\t%u\t%u\t%c\t%d\n",loaded_features[i].feature_name,loaded_features[i].chro,loaded_features[i].start, loaded_features[i].end, loaded_features[i].is_negative_strand?'-':'+', nreads[loaded_features[i].sorted_order]);
	}

	fclose(fp_out);
}

static struct option long_options[] =
{
	{0, 0, 0, 0}
};

void print_usage()
{
	SUBREADprintf("\nUsage: readSummary -i <input_file> -o <output_file> -a <annotation_file> { -T <n_threads> } {-b} {-O} {-p} {-S} {-G} {-m} {-d <min_PE_dist>} {-D <max_PE_dist>} {-F <annot_format>} \n\n  Optional parameters:\n   -T n\tThe number of threads, 1 by default.\n   -b  \tRead input file as BAM, false by default.\n   -O \tAllow overlapped exons, false by default.\n   -p \tUse paired-end information, false by default.\n   -F f\tThe format of annotations, either 'GTF' or 'Customized', 'GTF' by default\n   -m  \tDo gene-level aggregation, false by default.\n   -S  \tWrite the read classification result into output_file.reads.\n   -s  \tMatch the strand of features with the strands of reads mapped to.\n\n");
}

int readSummary(int argc,char *argv[]){

	/*
	   This function counts the number of reads falling into each exon region.
	   The order of exons in the output is the same as that of exons included in the annotation.
	   The annotation, if provided as a file, should be sorted by chromosome name.

	   Parameters passed from the featureCounts R function:
	0: "readSummary"
	1: ann
	2: files[i]
	3: fout
	4: as.numeric(isPairedEnd)
	5: min.distance
	6: max.distance
	7: as.numeric(tolower(file.type)=="sam")
	8: as.numeric(allowMultiOverlap)
	9: as.numeric(isGeneLevel)
	10: as.numeric(nthreads)
	11: as.numeric(isGTFannotation)
	12: as.numeric(isStrandChecked)
	13: as.numeric(isReadSummaryReported)
	 */

	FILE *fp_in = NULL;

	int read_length = 0, isStrandChecked, isCVersion;
	char **chr;
	long *start, *stop;
	int *geneid, *nreads;

	char * line = NULL;
	long nexons;


	int nreads_mapped_to_exon = 0;

	long * anno_chr_head;
	char ** anno_chrs;
	long curchr, curpos;
	char * curchr_name;
	unsigned char * sorted_strand;
	int nchr;
	curchr = 0;
	curpos = 0;
	curchr_name = "";

	int isPE, minPEDistance, maxPEDistance, isReadSummaryReport;

	int isSAM, isGTF;
	char * ret = NULL;
	SamBam_FILE * fp_in_bam = NULL;

	int isMultiOverlapAllowed, isGeneLevel;
	double time_start = miltime();

	line = (char*)calloc(MAX_LINE_LENGTH, 1);

	isCVersion = ((argv[0][0]=='C')?1:0);

	isPE = atoi(argv[4]);
	minPEDistance = atoi(argv[5]);
	maxPEDistance = atoi(argv[6]);

	isSAM = atoi(argv[7]);
	isMultiOverlapAllowed = atoi(argv[8]);
	isGeneLevel = atoi(argv[9]);
	unsigned short thread_number;
	if(argc > 10)
		thread_number = atoi(argv[10]);
	else	thread_number = 4;
	if(argc > 11)
		isGTF = atoi(argv[11]);
	else	isGTF = 0;
	if(argc > 12)
		isStrandChecked = atoi(argv[12]);
	else	isStrandChecked = 0;
	if(argc > 13)
		isReadSummaryReport = atoi(argv[13]);
	else	isReadSummaryReport = 0;

	if(thread_number<1) thread_number=1;
	if(thread_number>16)thread_number=16;

	unsigned int buffer_size = 1024*1024*6;


	// Open the SAM/BAM file
	// Nothing is done if the file does not exist.
	if (isSAM == 1)
	{
		#ifdef MAKE_STANDALONE
		if(strcmp("STDIN",argv[2])==0)
			fp_in = stdin;
		else
			fp_in = fopen(argv[2],"r");
		#else
			fp_in = fopen(argv[2],"r");
		#endif
		if(!fp_in){
			SUBREADprintf("Failed to open file %s. Please check if the file name and specified file type are correct.\n", argv[2]); 
			return -1;
		}
	}
	else
	{
		fp_in_bam = SamBam_fopen(argv[2], SAMBAM_FILE_BAM);
		if(!fp_in_bam){
			SUBREADprintf("Failed to open file %s. Please check if the file name and specified file type are correct.\n", argv[2]); 
			return -1;
		}
	}


	// Loading the annotations.
	// Nothing is done if the annotation does not exist.
	fc_feature_info_t * loaded_features;
	nexons = load_feature_info(argv[1], isGTF?FILE_TYPE_GTF:FILE_TYPE_RSUBREAD, &loaded_features);
	if(nexons<1){
		SUBREADprintf("Failed to open the annotation file %s\n",argv[1]);
		return -1;
	}


	fc_thread_global_context_t global_context;
	fc_thread_init_global_context(& global_context, buffer_size, thread_number, MAX_LINE_LENGTH, isPE, minPEDistance, maxPEDistance,isGeneLevel, isMultiOverlapAllowed, isStrandChecked, (char *)argv[3] , isReadSummaryReport);

	sort_feature_info(&global_context, nexons, loaded_features, &chr, &geneid, &start, &stop, &sorted_strand, &nchr, &anno_chrs, &anno_chr_head);
	nreads = (int *) calloc(nexons,sizeof(int));

	SUBREADprintf("Number of chromosomes included in the annotation is \%d\n",nchr);

	int thread_ret = 0;
	thread_ret |= fc_thread_start_threads(& global_context, nchr, nexons, geneid, chr, start, stop, sorted_strand, anno_chrs, anno_chr_head, read_length);

	int buffer_pairs = 16;
	char * preload_line = malloc(sizeof(char) * (2+MAX_LINE_LENGTH)*(isPE?2:1)*buffer_pairs);
	int preload_line_ptr;
	int current_thread_id = 0;

	while (1){
		int pair_no;
		int is_second_read;
		int fresh_read_no = 0;
		preload_line[0] = 0;
		preload_line_ptr = 0;
		

		for(pair_no=0; pair_no < buffer_pairs; pair_no++)
		{
			for(is_second_read=0;is_second_read<(isPE?2:1);is_second_read++)
			{
				while(1)
				{
					if (isSAM == 1)
						ret = fgets(line, MAX_LINE_LENGTH, fp_in);
					else
						ret = SamBam_fgets(fp_in_bam, line, MAX_LINE_LENGTH);  
					if(!ret) break;
					if(line[0] != '@') break;
				}

				if(!ret) break;
				int curr_line_len = strlen(line);

				// Try to figure out the length of the first read.
				// It is assumed that all reads have the same length.
				if(read_length < 1)
				{
					int tab_no = 0;
					int read_len_tmp=0, read_cursor;
					for(read_cursor=0; read_cursor<curr_line_len; read_cursor++)
					{
						if(line[read_cursor] == '\t')
							tab_no++;
						else
						{
							if(tab_no == 9)	// SEQ
								read_len_tmp++;
						}
					}
					read_length = read_len_tmp;
					global_context.read_length = read_length;
				}

				strcpy(preload_line+preload_line_ptr, line);
				preload_line_ptr += curr_line_len;
				if(line[curr_line_len-1]!='\n')
				{
					strcpy(preload_line+preload_line_ptr, "\n");
					preload_line_ptr++;
				}
				fresh_read_no++;
			}
			if(!ret) break;
		}


		int line_length = preload_line_ptr;

		if(isPE && (fresh_read_no%2>0))
		{
			// Safegarding -- it should not happen if the SAM file has a correct format.
			SUBREADprintf("WARNING! THE LAST READ IS UNPAIRED!!\nIGNORED!\n");
			line_length = 0;
		}

		if(line_length > 0)
		{
			while(1)
			{
				int is_finished = 0;
				fc_thread_thread_context_t * thread_context = global_context.thread_contexts+current_thread_id;

				pthread_spin_lock(&thread_context->input_buffer_lock);
				unsigned int empty_bytes = global_context.input_buffer_max_size -  thread_context->input_buffer_remainder; 
				if(empty_bytes > line_length)
				{
					unsigned int tail_bytes = global_context.input_buffer_max_size -  thread_context->input_buffer_write_ptr; 
					unsigned int write_p1_len = (tail_bytes > line_length)?line_length:tail_bytes;
					unsigned int write_p2_len = (tail_bytes > line_length)?0:(line_length - tail_bytes);
					memcpy(thread_context->input_buffer + thread_context->input_buffer_write_ptr, preload_line, write_p1_len);
					if(write_p2_len)
					{
						memcpy(thread_context->input_buffer, preload_line + write_p1_len, write_p2_len);
						thread_context->input_buffer_write_ptr = write_p2_len;
					}
					else	thread_context->input_buffer_write_ptr += write_p1_len;
					if(thread_context->input_buffer_write_ptr == global_context.input_buffer_max_size) 
						thread_context->input_buffer_write_ptr=0;


					thread_context->input_buffer_remainder += line_length;
					is_finished = 1;
				}

				pthread_spin_unlock(&thread_context->input_buffer_lock);

				current_thread_id++;
				if(current_thread_id >= thread_number) current_thread_id = 0;

				if(is_finished) break;
				else usleep(tick_time);
			}
		}

		if(!ret) break;

	}

	free(preload_line);
	global_context.is_all_finished = 1;
	fc_thread_wait_threads(&global_context);

	fc_thread_merge_results(&global_context, nreads , &nreads_mapped_to_exon);

	double time_end = miltime();
	if(isPE == 1)
		SUBREADprintf("Number of fragments mapped to the features is: %d\nTime cost = %.1f seconds\n\n", nreads_mapped_to_exon, time_end - time_start);
	else
		SUBREADprintf("Number of reads mapped to the features is: %d\nTime cost = %.1f seconds\n\n", nreads_mapped_to_exon, time_end - time_start);


	if(isCVersion && isGeneLevel)
		fc_write_final_gene_results(&global_context, geneid, start, stop, argv[3], nexons, nreads, loaded_features);
	else
		fc_write_final_results(&global_context, argv[3], nexons, nreads, loaded_features);

	if (isSAM == 1)
	{
		#ifdef MAKE_STANDALONE
		if(strcmp("STDIN",argv[2])!=0)
		#endif
			fclose(fp_in);
	}
	else
		SamBam_fclose(fp_in_bam);


	fc_thread_destroy_global_context(&global_context);
	free(line);
	free(loaded_features);
	free(geneid);
	free(chr);
	free(start);
	free(sorted_strand);
	free(anno_chrs);
	free(anno_chr_head);
	free(stop);
	free(nreads);

	return 0;
}


#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int feature_count_main(int argc, char ** argv)
#endif
{
	char * Rargv[14];
	char annot_name[300];
	char sam_name[300];
	char out_name[300];
	int min_dist = 50;
	int max_dist = 600;
	char min_dist_str[11];
	char max_dist_str[11];
	int is_PE = 0;
	int is_SAM = 1;
	int is_GeneLevel = 0;
	int is_Overlap = 0;
	int is_Strand_Sensitive = 0;
	int is_ReadSummary_Report = 0;
	int threads = 1;
	int isGTF = 1;
	char nthread_str[4];
	int option_index = 0;
	int c;

	annot_name[0]=0;sam_name[0]=0;out_name[0]=0;

	while ((c = getopt_long (argc, argv, "T:i:o:a:d:D:pbF:msOS?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'T':
				threads = atoi(optarg);
				break;
			case 'd':
				min_dist = atoi(optarg);
				break;
			case 'D':
				max_dist = atoi(optarg);
				break;
			case 'p':
				is_PE = 1;
				break;
			case 'b':
				is_SAM = 0;
				break;
			case 'm':
				is_GeneLevel = 1;
				break;
			case 'F':
				isGTF = (optarg[0]=='G');
				break;
			case 'O':
				is_Overlap = 1;
				break;
			case 'S':
				is_ReadSummary_Report = 1;
				break;
			case 's':
				is_Strand_Sensitive = 1;
				break;
			case 'i':
				strncpy(sam_name, optarg,299);
				break;
			case 'o':
				strncpy(out_name, optarg,299);
				break;
			case 'a':
				strncpy(annot_name, optarg,299);
				break;
			case '?':
			default :
				print_usage();
				return -1;
				break;
		}

	if(sam_name[0]==0||out_name[0]==0 || annot_name[0]==0)
	{
		SUBREADprintf("The input file, out put file and annotation file must be specified.\n");
		print_usage();
		return -1;
	}

	sprintf(nthread_str,"%d", threads);
	sprintf(min_dist_str,"%d",min_dist);
	sprintf(max_dist_str,"%d",max_dist);
	Rargv[0] = "CreadSummary";
	Rargv[1] = annot_name;
	Rargv[2] = sam_name;
	Rargv[3] = out_name;
	Rargv[4] = is_PE?"1":"0";
	Rargv[5] = min_dist_str;
	Rargv[6] = max_dist_str;
	Rargv[7] = is_SAM?"1":"0";
	Rargv[8] = is_Overlap?"1":"0";
	Rargv[9] = is_GeneLevel?"1":"0";
	Rargv[10] = nthread_str;
	Rargv[11] = isGTF?"1":"0";
	Rargv[12] = is_Strand_Sensitive?"1":"0";
	Rargv[13] = is_ReadSummary_Report?"1":"0";
	readSummary(14, Rargv);
	return 0;
}


