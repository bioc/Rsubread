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
#include "gene-algorithms.h"
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
	int is_both_end_required;
	int is_chimertc_disallowed;
	int is_PE_distance_checked;
	int min_paired_end_distance;
	int max_paired_end_distance;
	int read_length;
	int line_length;

	unsigned short thread_number;
	fc_thread_thread_context_t * thread_contexts;
	int is_all_finished;
	unsigned int input_buffer_max_size;

	HashTable * gene_name_table;	// gene_name -> gene_internal_number
	unsigned char ** gene_name_array;	// gene_internal_number -> gene_name 
	int exontable_nchrs;
	int exontable_exons;
	int * exontable_geneid;
	char * exontable_strand;
	char ** exontable_chr;
	long * exontable_start;
	long * exontable_stop;
	char feature_name_column[66];
	char gene_id_column[66];

	long * exontable_block_end_index;
	long * exontable_block_max_end;
	long * exontable_block_min_start;

	char ** exontable_anno_chrs;
	char * exontable_anno_chr_2ch;
	long * exontable_anno_chr_heads;

	FILE * SAM_output_fp;
	
} fc_thread_global_context_t;

unsigned int tick_time = 1000;


int is_comment_line(const char * l, int file_type)
{
	int tabs = 0, xk1 = 0;
	if(l[0]=='#') return 1;
	while(l[xk1]) tabs += (l[xk1++] == '\t');

	return tabs < (file_type == FILE_TYPE_GTF?8:4);
}

// This function loads annotations from the file.
// It returns the number of featres loaded, or -1 if something is wrong. 
// Memory will be allowcated in this function. The pointer is saved in *loaded_features.
// The invoker must release the memory itself.

int load_feature_info(fc_thread_global_context_t *global_context, const char * annotation_file, int file_type, fc_feature_info_t ** loaded_features)
{
	unsigned int features = 0, xk1 = 0;
	char * file_line = malloc(MAX_LINE_LENGTH+1);
	FILE * fp = fopen(annotation_file,"r"); 
	int is_GFF_warned = 0;
	if(!fp) return -1;
	if(file_type == FILE_TYPE_RSUBREAD)
		fgets(file_line, MAX_LINE_LENGTH, fp);

	while(1)
	{
		char * fgets_ret = fgets(file_line, MAX_LINE_LENGTH, fp);
		if(!fgets_ret) break;
		if(is_comment_line(file_line, file_type))continue;
		if(file_type == FILE_TYPE_GTF)
		{
			char * token_temp;
			strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp);
			char * feature_type = strtok_r(NULL,"\t", &token_temp);
			if(strcmp(feature_type, global_context -> feature_name_column)==0)
				features++;
		}
		else
			features++;
	}

	fseek(fp,0,SEEK_SET);
	if(file_type == FILE_TYPE_RSUBREAD)
		fgets(file_line, MAX_LINE_LENGTH, fp);

	fc_feature_info_t * ret_features = malloc(sizeof(fc_feature_info_t) * features);

	int lineno = 0;
	while(xk1 < features)
	{
		int is_gene_id_found = 0;
		fgets(file_line, MAX_LINE_LENGTH, fp);
		lineno++;
		char * token_temp;
		if(is_comment_line(file_line, file_type))continue;

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
			is_gene_id_found = 1;
			xk1++;
		}
		else if(file_type == FILE_TYPE_GTF)
		{
			sprintf((char *)ret_features[xk1].feature_name, "LINE_%07u", xk1 + 1);
			char * seq_name = strtok_r(file_line,"\t",&token_temp);
			strncpy((char *)ret_features[xk1].chro, (char *)seq_name, CHROMOSOME_NAME_LENGTH);

			strtok_r(NULL,"\t", &token_temp);// source
			char * feature_type = strtok_r(NULL,"\t", &token_temp);// feature_type
			if(strcmp(feature_type, global_context -> feature_name_column)==0)
			{
				ret_features[xk1].start = atoi(strtok_r(NULL,"\t", &token_temp));// start 
				ret_features[xk1].end = atoi(strtok_r(NULL,"\t", &token_temp));//end 

				if(ret_features[xk1].start < 1 || ret_features[xk1].end<1 || ret_features[xk1].start >= ret_features[xk1].end)
					SUBREADprintf("WARNING: the feature on the %d-th line has zero coordinate or zero lengths\n", lineno);

				strtok_r(NULL,"\t", &token_temp);// score 
				ret_features[xk1].is_negative_strand = ('-' == (strtok_r(NULL,"\t", &token_temp)[0]));//strand 
				ret_features[xk1].sorted_order = xk1;
				strtok_r(NULL,"\t",&token_temp);	// "frame"
				char * extra_attrs = strtok_r(NULL,"\t",&token_temp);	// name_1 "val1"; name_2 "val2"; ... 
				if(extra_attrs && (strlen(extra_attrs)>6))
				{
					int attr_state = 0;	// 0: name; 1:value
					int name_start = 0;
					int val_start = 0;
					int is_gene_id = 0;
					int xk2, exch;
					for(xk2 = 0;  extra_attrs[xk2]; xk2++)
					{
						exch = extra_attrs[xk2];
						if(exch == ' ' && attr_state == 0)
							continue;
						else if(exch == '#' && attr_state == 0)
							break;
						else if(exch == '\"')
						{
							if(xk2<1) break;
							if(attr_state == 0){
								val_start = xk2 + 1;
								extra_attrs[xk2 - 1] = 0;
								if(strcmp(extra_attrs + name_start, global_context -> gene_id_column)==0)
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
									is_gene_id_found = 1;
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

		if(!is_gene_id_found)
		{
			if(!is_GFF_warned)
				SUBREADprintf("**********\n**********\n  WARNING\n**********\nNo meta-feature id is found on the %d-th line. If it is a GTF file, you may need to check the name of the gene_id field and specify a correct field name using a '-g' option.\n**********\n**********\n", lineno);
			is_GFF_warned++;
		}
	}
	fclose(fp);
	free(file_line);

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


#define FC_FAST_BLOCK_SIZE 70
#define FC_MAX_CHROMOSOME_NUMBER 500

void sort_feature_info(fc_thread_global_context_t * global_context, unsigned int features, fc_feature_info_t * loaded_features, char *** sorted_chr_names, int ** sorted_entrezid, long ** sorted_start, long ** sorted_end, unsigned char ** sorted_strand, int * chro_number, char ** anno_chr_2ch, char *** anno_chrs, long ** anno_chr_head, long ** block_end_index, long ** block_min_start_pos, long ** block_max_end_pos)
{
	unsigned int chro_pnt;
	unsigned char * is_chro_name_explored = malloc(features);
	int * ret_entrez = malloc(sizeof(int) * features);
	long * ret_start = malloc(sizeof(long) * features);
	long * ret_end = malloc(sizeof(long) * features);
	long * ret_block_end_index = malloc(sizeof(long) * features);
	long * ret_block_min_start = malloc(sizeof(long) * features);
	long * ret_block_max_end = malloc(sizeof(long) * features);
	unsigned char * ret_strand = malloc(features);
	char ** ret_char_name = malloc(sizeof(void *) * features);
	fc_feature_info_t ** old_info_ptr = malloc(sizeof(void *) * features);
	(*anno_chrs) = malloc(sizeof(void *) * FC_MAX_CHROMOSOME_NUMBER);
	(*anno_chr_head) = malloc(sizeof(long) * FC_MAX_CHROMOSOME_NUMBER);
	(*anno_chr_2ch) = malloc(sizeof(char) * FC_MAX_CHROMOSOME_NUMBER*2); 

	global_context -> gene_name_array = malloc(sizeof(char *) * features);	// there should be much less identical names.
	global_context -> gene_name_table = HashTableCreate(5000);
	HashTableSetHashFunction(global_context -> gene_name_table, HashTableStringHashFunction);
	HashTableSetKeyComparisonFunction(global_context -> gene_name_table, fc_strcmp);


	unsigned int this_chro_start = 0;
	unsigned int this_chro_end = 0;
	unsigned int ret_array_pointer = 0;
	int current_block_id = 0;

	memset(is_chro_name_explored, 0 , features);

	(*sorted_chr_names) = ret_char_name;
	(*sorted_entrezid) = ret_entrez;
	(*sorted_start) = ret_start;
	(*sorted_end) = ret_end;
	(*sorted_strand) = ret_strand;
	(*block_end_index) = ret_block_end_index;
	(*block_min_start_pos) = ret_block_min_start;
	(*block_max_end_pos) = ret_block_max_end;
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
			int chro_name_len = strlen((char *)current_chro);
			(*anno_chr_2ch)[(*chro_number)*2+1] = current_chro[chro_name_len-1];
			if(chro_name_len>1)
				(*anno_chr_2ch)[(*chro_number)*2] = current_chro[chro_name_len-2];
			else
				(*anno_chr_2ch)[(*chro_number)*2] = 0;

			(*anno_chr_head) [(*chro_number)] = current_block_id;
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

			// sorting the features.
			unsigned int sort_i, sort_j, this_block_items = 0;
			long this_block_min_start = 0x7fffffff, this_block_max_end = 0;
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
				this_block_max_end = max(this_block_max_end, ret_end[sort_i]);
				this_block_min_start = min(this_block_min_start, ret_start[sort_i]);
				this_block_items ++;
				if(this_block_items >= FC_FAST_BLOCK_SIZE)
				{
					ret_block_end_index[current_block_id] = sort_i+1;
					ret_block_min_start[current_block_id] = this_block_min_start;
					ret_block_max_end[current_block_id] = this_block_max_end;
					current_block_id++;
					this_block_max_end = 0;
					this_block_items = 0;
					this_block_min_start = 0x7fffffff;
				}
			}
			if(this_block_items)
			{
				ret_block_end_index[current_block_id] = this_chro_end;
				ret_block_min_start[current_block_id] = this_block_min_start;
				ret_block_max_end[current_block_id] = this_block_max_end;
				current_block_id++;
			}

			(*anno_chr_head) [(*chro_number)] = current_block_id; 
		}
	}
	SUBREADprintf("The %u features are sorted.\n", ret_array_pointer);
	free(old_info_ptr);
}

#define MAX_HIT_NUMBER 80

void process_line_buffer(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context)
{

	char * read_chr, *tmp_tok_ptr, *CIGAR_str , *read_name = NULL;
	long read_pos, fragment_length = 0, search_start, search_end;
	int is_chro_found, nhits1 = 0, nhits2 = 0, alignment_masks, search_block_id, search_item_id;
	long hits_indices1[MAX_HIT_NUMBER], hits_indices2[MAX_HIT_NUMBER];

	int is_second_read;

	for(is_second_read = 0 ; is_second_read < 2; is_second_read++)
	{
		if(is_second_read && !global_context -> is_paired_end_data) break;

		char * line = is_second_read? thread_context -> line_buffer2:thread_context -> line_buffer1;
	
		read_name = strtok_r(line,"\t", &tmp_tok_ptr);	// read name
		alignment_masks = atoi(strtok_r(NULL,"\t", &tmp_tok_ptr));

		if(is_second_read == 0)
		{
			//skip the read if unmapped (its mate will be skipped as well if paired-end)
			if( ((!global_context -> is_paired_end_data) &&  (alignment_masks & SAM_FLAG_UNMAPPED) ) ||
			    ((alignment_masks & SAM_FLAG_UNMAPPED)   &&  (alignment_masks & SAM_FLAG_MATE_UNMATCHED) && global_context -> is_paired_end_data) ||
			    (((alignment_masks & SAM_FLAG_UNMAPPED) || (alignment_masks & SAM_FLAG_MATE_UNMATCHED)) && global_context -> is_paired_end_data && global_context -> is_both_end_required)
			  ){
				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUNMAPPED\n", read_name);
				return;	// do nothing if a read is unmapped, or the first read in a pair of reads is unmapped.
			}
		}



		read_chr = strtok_r(NULL,"\t", &tmp_tok_ptr);
		read_pos = atoi(strtok_r(NULL,"\t", &tmp_tok_ptr));
		strtok_r(NULL,"\t", &tmp_tok_ptr);	// mapping quality
		CIGAR_str = strtok_r(NULL,"\t", &tmp_tok_ptr);	// CIGAR string

		if(is_second_read == 0 && global_context -> is_paired_end_data && 
	   	  (global_context -> is_PE_distance_checked || global_context -> is_chimertc_disallowed)
		  )
		{
			int is_half_mapped = (alignment_masks & SAM_FLAG_UNMAPPED) || (alignment_masks & SAM_FLAG_MATE_UNMATCHED);

			if(!is_half_mapped)
			{
				char * mate_chr = strtok_r(NULL,"\t", &tmp_tok_ptr); //get chr which the mate read is mapped to
				atoi(strtok_r(NULL,"\t", &tmp_tok_ptr));	// mate_pos
				char * frag_len_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
				fragment_length = abs(atoi(frag_len_str)); //get the fragment length


				int is_first_read_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
				int is_second_read_negative_strand = (alignment_masks & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)?1:0; 

				if(mate_chr[0]=='=' && is_first_read_negative_strand!=is_second_read_negative_strand)
				{
					if((global_context -> is_PE_distance_checked && fragment_length > global_context -> max_paired_end_distance + thread_context->current_read_length1 -1) || (fragment_length < global_context -> min_paired_end_distance + thread_context->current_read_length1 - 1))
					{
						if(global_context -> SAM_output_fp)
							fprintf(global_context -> SAM_output_fp,"%s\tPAIR_DISTANCE\t%ld\n", read_name, fragment_length);
						return;
					}
				}
				else
				{
					if(global_context -> is_chimertc_disallowed)
					{
						if(global_context -> SAM_output_fp)
							fprintf(global_context -> SAM_output_fp,"%s\tPAIR_CHIMERIC\n", read_name);
						return;
					}
				}
			}
		}

		if(SAM_FLAG_UNMAPPED & alignment_masks) continue;

		int is_this_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
		if(is_second_read) is_this_negative_strand = !is_this_negative_strand;


		is_chro_found = 0;
		char last_chr_chr, second_chr_chr;
		int chro_name_len = strlen(read_chr);
		last_chr_chr = read_chr[chro_name_len-1];
		if(chro_name_len>1)
			second_chr_chr = read_chr[chro_name_len-2];
		else	second_chr_chr = 0 ;

		for(search_item_id=0;search_item_id<global_context -> exontable_nchrs;search_item_id++)
		{
			if(last_chr_chr != global_context -> exontable_anno_chr_2ch[search_item_id*2+1]||second_chr_chr != global_context -> exontable_anno_chr_2ch[search_item_id*2])
				continue;

			if(memcmp(read_chr,global_context -> exontable_anno_chrs[search_item_id], chro_name_len+1)==0){
				//get chr to which the current read or fragment is mapped and also the searching range
				is_chro_found = 1;
				search_start = global_context -> exontable_anno_chr_heads[search_item_id];
				search_end = global_context -> exontable_anno_chr_heads[search_item_id+1] - 1;
				break;
			}
		}

		//printf("FX=%d\n", is_chro_found);
		if(is_chro_found)
		{
			int nhits = 0;

			int cigar_section_id, cigar_sections;
			unsigned int Staring_Points[6];
			unsigned short Section_Lengths[6];
			long * hits_indices = (is_second_read?hits_indices2:hits_indices1);

			cigar_sections = RSubread_parse_CIGAR_string(CIGAR_str, Staring_Points, Section_Lengths);
			for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
			{
				long section_begin_pos = read_pos + Staring_Points[cigar_section_id];
				long section_length = Section_Lengths[cigar_section_id];
		//		printf("CIGAR_str=%s; cigar_sections=%d; base=%ld; pos[%d]=%u ; len[%d]=%d\n", CIGAR_str, cigar_sections, read_pos, cigar_section_id, Staring_Points[cigar_section_id],cigar_section_id, Section_Lengths[cigar_section_id]);
				for(search_block_id=search_start;search_block_id<=search_end;search_block_id++){
					if (global_context -> exontable_block_min_start[search_block_id] > (section_begin_pos + section_length -1)) break;
					if (global_context -> exontable_block_max_end[search_block_id] < section_begin_pos) continue;

					int search_item_start = 0, search_item_end = global_context -> exontable_block_end_index[search_block_id];
					if(search_block_id>0)search_item_start = global_context -> exontable_block_end_index[search_block_id-1];
					for(search_item_id = search_item_start ; search_item_id < search_item_end; search_item_id++)
					{
						if (global_context -> exontable_start[search_item_id] > (section_begin_pos + section_length -1)) break;
						if (global_context -> exontable_stop[search_item_id] >= section_begin_pos && ((!global_context->is_strand_checked)||is_this_negative_strand == global_context -> exontable_strand[search_item_id])){
							hits_indices[nhits] = search_item_id;
							nhits++;
							if(nhits>=MAX_HIT_NUMBER) break;
						} 
					}
				}
				if(nhits>=MAX_HIT_NUMBER) break;
			}

			if(is_second_read) nhits2 = nhits;
			else	nhits1 = nhits;
		}
	}	// loop for is_second_read


	if(nhits2+nhits1==1)
	{
		long hit_exon_id = nhits2?hits_indices2[0]:hits_indices1[0];
		thread_context->count_table[hit_exon_id]++;
		thread_context->nreads_mapped_to_exon++;
		if(global_context -> SAM_output_fp)
		{
			int final_gene_number = global_context -> exontable_geneid[hit_exon_id];
			unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
			fprintf(global_context -> SAM_output_fp,"%s\tACCEPTED_%s\t%s\n", read_name, global_context -> is_gene_level?"GENE":"EXON", final_feture_name);
		}
	}
	else if(nhits2 == 1 && nhits1 == 1 && hits_indices2[0]==hits_indices1[0])
	{
		long hit_exon_id = hits_indices1[0];
		thread_context->count_table[hit_exon_id]++;
		thread_context->nreads_mapped_to_exon++;
		if(global_context -> SAM_output_fp)
		{
			int final_gene_number = global_context -> exontable_geneid[hit_exon_id];
			unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
			fprintf(global_context -> SAM_output_fp,"%s\tACCEPTED_%s\t%s\n", read_name, global_context -> is_gene_level?"GENE":"EXON", final_feture_name);
		}
	}
	else
	{

		long decision_table_ids[MAX_HIT_NUMBER];
		unsigned char decision_table_votes[MAX_HIT_NUMBER];
		long decision_table_exon_ids[MAX_HIT_NUMBER];
		int decision_table_items = 0, xk1, xk2;

		for(is_second_read = 0; is_second_read < 2; is_second_read++)
		{
			if(is_second_read && !global_context -> is_paired_end_data) break;
			long * hits_indices = is_second_read?hits_indices2:hits_indices1;
			int nhits = is_second_read?nhits2:nhits1;
			if (nhits<1) continue;
			if(global_context -> is_gene_level)
			{
				long uniq_gene_table[MAX_HIT_NUMBER];
				long uniq_gene_exonid_table[MAX_HIT_NUMBER];
				int uniq_genes = 0;
				for(xk1=0;xk1<nhits;xk1++)
				{
					int gene_id = global_context -> exontable_geneid[hits_indices[xk1]];
					int is_unique = 1;
					for(xk2=0; xk2<uniq_genes; xk2++)
					{
						if(gene_id == uniq_gene_table[xk2])
						{
							is_unique = 0;
							break;
						}
					}
					if(is_unique){
						uniq_gene_exonid_table[uniq_genes] = hits_indices[xk1];
						uniq_gene_table[uniq_genes++] = gene_id;
					}
				}

				for(xk1=0;xk1<uniq_genes; xk1++)
				{
					long gene_id = uniq_gene_table[xk1];
					int is_fresh = 1;
					if(decision_table_items >= MAX_HIT_NUMBER) break;
					for(xk2=0; xk2<decision_table_items; xk2++)
					{
						if(gene_id == decision_table_ids[xk2])
						{
							decision_table_votes[xk2]++;
							is_fresh = 0;
							break;
						}
						
					}
					if(is_fresh)
					{
						decision_table_votes[decision_table_items] = 1;
						decision_table_exon_ids[decision_table_items] = uniq_gene_exonid_table[xk1];
						decision_table_ids[decision_table_items++] = gene_id;
					}
				}
			}
			else
			{
				for(xk1=0;xk1<nhits;xk1++)
				{
					long exon_id = global_context -> exontable_geneid[hits_indices[xk1]];
					int is_fresh = 1;
					if(decision_table_items >= MAX_HIT_NUMBER) break;
					for(xk2=0; xk2<decision_table_items; xk2++)
					{
						if(exon_id == decision_table_ids[xk2])
						{
							decision_table_votes[xk2]++;
							is_fresh = 0;
							break;
						}
					}
					if(is_fresh)
					{
						decision_table_votes[decision_table_items] = 1;
						decision_table_ids[decision_table_items++] = exon_id;
					}

				}
			}

		}
		if(decision_table_items>0)
		{
			int max_votes = 0;
			int top_voters = 0;
			long top_voter_id = 0;

			for(xk1 = 0; xk1 < decision_table_items; xk1++)
			{
				if(decision_table_votes[xk1] > max_votes)
				{
					max_votes = decision_table_votes[xk1];
					top_voters = 1;
					top_voter_id = (global_context -> is_gene_level)?decision_table_exon_ids[xk1]:decision_table_ids[xk1];
				}
				else
					if(decision_table_votes[xk1] == max_votes) top_voters++;
			}

			if(top_voters == 1)
			{
				thread_context->count_table[top_voter_id]++;
				thread_context->nreads_mapped_to_exon++;
				if(global_context -> SAM_output_fp)
				{
					int final_gene_number = global_context -> exontable_geneid[top_voter_id];
					unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
					fprintf(global_context -> SAM_output_fp,"%s\tACCEPTED_%s\t%s\n", read_name, global_context -> is_gene_level?"GENE":"EXON", final_feture_name);
				}
			}
			else if(top_voters >1)
			{
				if(global_context -> is_multi_overlap_allowed)
				{
					char final_feture_names[1000];
					final_feture_names[0]=0;
					for(xk1 = 0; xk1 < decision_table_items; xk1++)
					{
						if(decision_table_votes[xk1] == max_votes)
						{
							long tmp_voter_id = (global_context -> is_gene_level)?decision_table_exon_ids[xk1]:decision_table_ids[xk1];
							thread_context->count_table[tmp_voter_id]++;

							if(global_context -> SAM_output_fp)
							{
								int final_gene_number = global_context -> exontable_geneid[tmp_voter_id];
								unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
								strncat(final_feture_names, (char *)final_feture_name, 999);
								strncat(final_feture_names, ",", 999);
							}
						}
					}
					thread_context->nreads_mapped_to_exon++;
					if(global_context -> SAM_output_fp)
					{
						int ffnn = strlen(final_feture_names);
						if(ffnn>0) final_feture_names[ffnn-1]=0;
						fprintf(global_context -> SAM_output_fp,"%s\tACCEPTED_MULTI_%sS\t%s\n", read_name, global_context -> is_gene_level?"GENE":"EXON", final_feture_names);
					}
				}
				else if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tOVERLAPPED_%sS\t%d\n", read_name, global_context -> is_gene_level?"GENE":"EXON", top_voters);
			}
		}
		else if(global_context -> SAM_output_fp)
			fprintf(global_context -> SAM_output_fp,"%s\tNOTFOUND_%s\n", read_name, global_context -> is_gene_level?"GENE":"EXON");
	}

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

void fc_thread_init_global_context(fc_thread_global_context_t * global_context, unsigned int buffer_size, unsigned short threads, int line_length , int is_PE_data, int min_pe_dist, int max_pe_dist, int is_gene_level, int is_overlap_allowed, int is_strand_checked, char * output_fname, int is_sam_out, int is_both_end_required, int is_chimertc_disallowed, int is_PE_distance_checked, char *feature_name_column, char * gene_id_column)
{

	global_context -> input_buffer_max_size = buffer_size;

	global_context -> is_multi_overlap_allowed = is_overlap_allowed;
	global_context -> is_paired_end_data = is_PE_data;
	global_context -> is_gene_level = is_gene_level;
	global_context -> is_strand_checked = is_strand_checked;
	global_context -> is_both_end_required = is_both_end_required;
	global_context -> is_chimertc_disallowed = is_chimertc_disallowed;
	global_context -> is_PE_distance_checked = is_PE_distance_checked;
	strcpy(global_context -> feature_name_column,feature_name_column);
	strcpy(global_context -> gene_id_column,gene_id_column);

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
int fc_thread_start_threads(fc_thread_global_context_t * global_context, int et_nchrs, int et_exons, int * et_geneid, char ** et_chr, long * et_start, long * et_stop, unsigned char * et_strand, char * et_anno_chr_2ch, char ** et_anno_chrs, long * et_anno_chr_heads, long * et_bk_end_index, long * et_bk_min_start, long * et_bk_max_end, int read_length)
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
	global_context -> exontable_anno_chr_2ch = et_anno_chr_2ch;
	global_context -> exontable_anno_chrs = et_anno_chrs;
	global_context -> exontable_anno_chr_heads = et_anno_chr_heads;
	global_context -> exontable_block_end_index = et_bk_end_index;
	global_context -> exontable_block_max_end = et_bk_max_end;
	global_context -> exontable_block_min_start = et_bk_min_start;


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

		if(global_context -> thread_number>1)
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
	SUBREADputs("\nUsage: featureCounts -a <annotation_file> -i <input_file> -o <output_file> {optional parameters} \n");
	SUBREADputs("    Required parameters:\n"); 
	SUBREADputs("    -a <input>\tGive the name of the annotation file. The program assumes"); 
	SUBREADputs("              \tthat the provided annotation file is in GTF format. Use -F"); 
	SUBREADputs("              \toption to specify other annotation formats."); 
	SUBREADputs("    "); 
	SUBREADputs("    -i <input>\tGive the name of input file including the read data. The "); 
	SUBREADputs("              \tfile can have a SAM or BAM format. -b option needs to be"); 
	SUBREADputs("              \tspecified if the file is in BAM format. "); 
	SUBREADputs("              "); 
	SUBREADputs("    -o <input>\tGive the name of the output file. The output file contains"); 
	SUBREADputs("              \tthe number of reads assigned to each meta-feature (or each"); 
	SUBREADputs("              \tfeature if -f is specified). A meta-feature is the aggregation");
	SUBREADputs("              \tof features, grouped by using gene identifiers. Please refer");
	SUBREADputs("              \tto the users guide for more details."); 
	SUBREADputs("    "); 
	SUBREADputs("    Optional parameters:"); 
	SUBREADputs("    "); 
	SUBREADputs("    -F        \tSpecify the format of the annotation file. Acceptable formats");
	SUBREADputs("              \tinclude `GTF' and `SAF'. `GTF' by default. Please refer to the");
	SUBREADputs("              \tusers guide for SAF annotation format."); 
	SUBREADputs("    "); 
	SUBREADputs("    -t <input>\tSpecify the feature type. Only rows which have the matched"); 
	SUBREADputs("              \tmatched feature type in the provided GTF annotation file"); 
	SUBREADputs("              \t will be included for read counting. `exon' by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -g <input>\tSpecify the attribute type used to group features (eg. exons)");
	SUBREADputs("              \tinto meta-features (eg. genes), when GTF annotation is provided.");
	SUBREADputs("              \t`gene_id' by default. This attribute type is usually the gene");
	SUBREADputs("              \tidentifier. This argument is useful for the meta-feature level");
	SUBREADputs("              \tsummarization.");
	SUBREADputs("    "); 
	SUBREADputs("    -b        \tIndicate that the input file is in BAM format."); 
	SUBREADputs("    "); 
	SUBREADputs("    -f        \tIf specified, read summarization will be performed at the "); 
	SUBREADputs("              \tfeature level. By default (-f is not specified), the read"); 
	SUBREADputs("              \tsummarization is performed at the meta-feature level."); 
	SUBREADputs("    "); 
	SUBREADputs("    -O        \tIf specified, reads (or fragments if -p is specified) will"); 
	SUBREADputs("              \tbe allowed to be assigned to more than one matched meta-"); 
	SUBREADputs("              \tfeature (or matched feature if -f is specified). "); 
	SUBREADputs("    "); 
	SUBREADputs("    -s        \tIf specified, strand-specific read assignment will be "); 
	SUBREADputs("              \tperformed."); 
	SUBREADputs("    "); 
	SUBREADputs("    -T <int>  \tNumber of the threads. 1 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -R        \tOutput the read assignment result for each read."); 
	SUBREADputs("    "); 
	SUBREADputs("    Optional paired-end parameters:"); 
	SUBREADputs("    "); 
	SUBREADputs("    -p        \tIf specified, fragments (or templates) will be counted "); 
	SUBREADputs("              \tinstead of reads. This option is only applicable for "); 
	SUBREADputs("              \tpaired-end reads. "); 
	SUBREADputs("    "); 
	SUBREADputs("    -P        \tIf specified, paired-end distance will be checked when "); 
	SUBREADputs("              \tassigning fragments to meta-features or features. This "); 
	SUBREADputs("              \toption is only applicable when -p is specified. The "); 
	SUBREADputs("              \tdistance thresholds should be specified using -d and -D "); 
	SUBREADputs("              \toptions."); 
	SUBREADputs("    "); 
	SUBREADputs("    -d <int>  \tMinimal allowed paired-end distance. 50 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -D <int>  \tMaximal allowed paired-end distance. 600 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -B        \tIf specified, only fragments that have both ends "); 
	SUBREADputs("              \tsuccessfully aligned will be considered for summarization."); 
	SUBREADputs("              \tThis option is only applicable for paired-end reads."); 
	SUBREADputs("    "); 
	SUBREADputs("    -C        \tIf specified, the chimeric fragments (those fragments that "); 
	SUBREADputs("              \thave their two ends aligned to different chromosomes) will"); 
	SUBREADputs("              \tNOT be included for summarization. This option is only "); 
	SUBREADputs("              \tapplicable for paired-end read data."); 
	SUBREADputs("    "); 


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
	14: as.numeric(isBothEndMapped)
	15: as.numeric(isChimericDisallowed)
	16: as.numeric(isPEDistChecked)
	17: nameFeatureTypeColumn 
	18: nameGeneIDColumn
	 */

	FILE *fp_in = NULL;

	int read_length = 0, isStrandChecked, isCVersion, isChimericDisallowed, isPEDistChecked;
	char **chr;
	long *start, *stop;
	int *geneid, *nreads;

	char * line = NULL, *nameFeatureTypeColumn, *nameGeneIDColumn;
	long nexons;


	int nreads_mapped_to_exon = 0;

	long * anno_chr_head, * block_min_start, *block_max_end, *block_end_index;
	char ** anno_chrs, * anno_chr_2ch;
	long curchr, curpos;
	char * curchr_name;
	unsigned char * sorted_strand;
	int nchr;
	curchr = 0;
	curpos = 0;
	curchr_name = "";

	int isPE, minPEDistance, maxPEDistance, isReadSummaryReport, isBothEndRequired;

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
	if(argc > 14)
		isBothEndRequired = atoi(argv[14]);
	else	isBothEndRequired = 0;
	if(argc > 15)
		isChimericDisallowed = atoi(argv[15]);
	else	isChimericDisallowed = 0;
	if(argc > 16)
		isPEDistChecked = atoi(argv[16]);
	else	isPEDistChecked = 0;
	if(argc > 17)
		nameFeatureTypeColumn = argv[17];
	else	nameFeatureTypeColumn = "exon";
	if(argc > 18)
		nameGeneIDColumn = argv[18];
	else	nameGeneIDColumn = "gene_id";

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
			SUBREADprintf("Failed to open file %s. Please check if the file name is correct and it is a SAM file.\n", argv[2]); 
			return -1;
		}
	}
	else
	{
		fp_in_bam = SamBam_fopen(argv[2], SAMBAM_FILE_BAM);
		if(!fp_in_bam){
			SUBREADprintf("Failed to open file %s. Please check if the file name is correct and it is a BAM file.\n", argv[2]); 
			return -1;
		}
	}


	fc_thread_global_context_t global_context;
	fc_thread_init_global_context(& global_context, buffer_size, thread_number, MAX_LINE_LENGTH, isPE, minPEDistance, maxPEDistance,isGeneLevel, isMultiOverlapAllowed, isStrandChecked, (char *)argv[3] , isReadSummaryReport, isBothEndRequired, isChimericDisallowed, isPEDistChecked, nameFeatureTypeColumn, nameGeneIDColumn);


	// Loading the annotations.
	// Nothing is done if the annotation does not exist.
	fc_feature_info_t * loaded_features;
	nexons = load_feature_info(&global_context, argv[1], isGTF?FILE_TYPE_GTF:FILE_TYPE_RSUBREAD, &loaded_features);
	if(nexons<1){
		SUBREADprintf("Failed to open the annotation file %s, or its format is incorrect, or it contains no '%s' features.\n",argv[1], nameFeatureTypeColumn);
		return -1;
	}

	sort_feature_info(&global_context, nexons, loaded_features, &chr, &geneid, &start, &stop, &sorted_strand, &nchr, &anno_chr_2ch, &anno_chrs, &anno_chr_head, & block_end_index, & block_min_start, & block_max_end);
	nreads = (int *) calloc(nexons,sizeof(int));

	SUBREADprintf("Number of chromosomes included in the annotation is \%d\n",nchr);

	int thread_ret = 0;
	thread_ret |= fc_thread_start_threads(& global_context, nchr, nexons, geneid, chr, start, stop, sorted_strand, anno_chr_2ch, anno_chrs, anno_chr_head, block_end_index, block_min_start , block_max_end, read_length);

	int buffer_pairs = thread_number>1?20:1;
	char * preload_line = malloc(sizeof(char) * (2+MAX_LINE_LENGTH)*(isPE?2:1)*buffer_pairs);
	int preload_line_ptr;
	int current_thread_id = 0;
	fc_thread_thread_context_t * one_thread_context = global_context.thread_contexts;

	while (1){
		int pair_no;
		int is_second_read;
		int fresh_read_no = 0;
		preload_line[0] = 0;
		preload_line_ptr = 0;


		
		if(thread_number==1)
		{
			int is_second_read;

			for(is_second_read=0;is_second_read<(isPE?2:1);is_second_read++)
			{
				char * lbuf = is_second_read?one_thread_context -> line_buffer2:one_thread_context -> line_buffer1;
				while(1)
				{
					if (isSAM == 1)
						ret = fgets(lbuf, MAX_LINE_LENGTH, fp_in);
					else
						ret = SamBam_fgets(fp_in_bam, lbuf, MAX_LINE_LENGTH);  
					if(!ret) break;
					if(lbuf[0] != '@') break;
				}

				if(!ret) break;
				if(read_length < 1)
				{
					int tab_no = 0;
					int read_len_tmp=0, read_cursor;
					int curr_line_len = strlen(lbuf);
					for(read_cursor=0; read_cursor<curr_line_len; read_cursor++)
					{
						if(lbuf[read_cursor] == '\t')
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
			}

			if(!ret) break;
			
			one_thread_context -> current_read_length1 = global_context.read_length;
			one_thread_context -> current_read_length2 = global_context.read_length;

			process_line_buffer(&global_context, one_thread_context);
		}
		else
		{
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
	}

	free(preload_line);
	global_context.is_all_finished = 1;

	if(thread_number > 1)
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
	free(anno_chr_2ch);
	free(anno_chrs);
	free(anno_chr_head);
	free(block_min_start);
	free(block_max_end);
	free(block_end_index);
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
	char * Rargv[19];
	char annot_name[300];
	char sam_name[300];
	char out_name[300];
	char nameFeatureTypeColumn[66];
	char nameGeneIDColumn[66];
	int min_dist = 50;
	int max_dist = 600;
	char min_dist_str[11];
	char max_dist_str[11];
	int is_PE = 0;
	int is_SAM = 1;
	int is_GeneLevel = 1;
	int is_Overlap = 0;
	int is_Both_End_Mapped = 0;
	int is_Strand_Sensitive = 0;
	int is_ReadSummary_Report = 0;
	int is_Chimeric_Disallowed = 0;
	int is_PE_Dist_Checked = 0;
	int threads = 1;
	int isGTF = 1;
	char nthread_str[4];
	int option_index = 0;
	int c;

	strcpy(nameFeatureTypeColumn,"exon");
	strcpy(nameGeneIDColumn,"gene_id");
	annot_name[0]=0;sam_name[0]=0;out_name[0]=0;

	while ((c = getopt_long (argc, argv, "g:t:T:i:o:a:d:D:pbF:fsCBPORv?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'v':
				print_version_info();
				return 0;
				break;
			case 't':
				strcpy(nameFeatureTypeColumn, optarg);
				break;
			case 'g':
				strcpy(nameGeneIDColumn, optarg);
				break;
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
			case 'C':
				is_Chimeric_Disallowed = 1;
				break;
			case 'P':
				is_PE_Dist_Checked = 1;
				break;
			case 'B':
				is_Both_End_Mapped = 0;
				break;
			case 'f':
				is_GeneLevel = 0;
				break;
			case 'F':
				isGTF = 1;
				if(strcmp("SAF", optarg)==0) isGTF=0;
				else if(strcmp("GTF", optarg)==0) isGTF=1;
				else SUBREADprintf("WARNING: Unknown annotation format: %s. GTF format is used.\n", optarg); 
				break;
			case 'O':
				is_Overlap = 1;
				break;
			case 'R':
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
	Rargv[14] = is_Both_End_Mapped?"1":"0";
	Rargv[15] = is_Chimeric_Disallowed?"1":"0";
	Rargv[16] = is_PE_Dist_Checked?"1":"0";
	Rargv[17] = nameFeatureTypeColumn;
	Rargv[18] = nameGeneIDColumn;
	readSummary(19, Rargv);
	return 0;
}


