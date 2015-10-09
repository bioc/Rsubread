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
  
  
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>


#ifndef MAKE_STANDALONE
  #include <R.h>
#endif

#include <zlib.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>
#include "subread.h"
#include "interval_merge.h"
#include "core.h"
#include "gene-algorithms.h"
#include "sambam-file.h"
#include "input-files.h"
#include "hashtable.h"
#include "HelperFunctions.h"

/********************************************************************/
/********************************************************************/
/********************************************************************/
//  NEW FUNCTION FOR MULTI-THREADING
/********************************************************************/
/********************************************************************/
/********************************************************************/
#define FEATURE_NAME_LENGTH  256 
#define CHROMOSOME_NAME_LENGTH 256 
#define MAX_LINE_LENGTH 3000
#define FILE_TYPE_RSUBREAD 10
#define FILE_TYPE_GTF 100

#define ALLOW_ALL_MULTI_MAPPING 1
#define ALLOW_PRIMARY_MAPPING 2

#define MAX_HIT_NUMBER 3000

typedef struct {
	char chromosome_name_left[CHROMOSOME_NAME_LENGTH + 1];
	char chromosome_name_right[CHROMOSOME_NAME_LENGTH + 1];
	unsigned int last_exon_base_left;
	unsigned int first_exon_base_right;
} fc_junction_info_t;


typedef struct{
	char gene_name[FEATURE_NAME_LENGTH];
	unsigned int pos_first_base;
	unsigned int pos_last_base;
} fc_junction_gene_t;

typedef struct {
	unsigned int feature_name_pos;
	unsigned int start;
	unsigned int end;
	unsigned int sorted_order;

	unsigned short chro_name_pos_delta;
	char is_negative_strand;
} fc_feature_info_t;

typedef struct {
	unsigned long long assigned_reads;
	unsigned long long unassigned_ambiguous;
	unsigned long long unassigned_multimapping;
	unsigned long long unassigned_nofeatures;
	unsigned long long unassigned_unmapped;
	unsigned long long unassigned_mappingquality;
	unsigned long long unassigned_fragmentlength;
	unsigned long long unassigned_chimericreads;
	unsigned long long unassigned_secondary;
	unsigned long long unassigned_nonjunction;
	unsigned long long unassigned_duplicate;
} fc_read_counters;

typedef unsigned long long read_count_type_t;

typedef struct {
	unsigned short thread_id;
	char * line_buffer1;
	char * line_buffer2;
	unsigned long long int nreads_mapped_to_exon;
	unsigned long long int all_reads;
	//unsigned short current_read_length1;
	//unsigned short current_read_length2;
	read_count_type_t * count_table;
	read_count_type_t unpaired_fragment_no;
	unsigned int chunk_read_ptr;
	pthread_t thread_object;

	char * input_buffer;
	unsigned int input_buffer_remainder;
	unsigned int input_buffer_write_ptr;	
	pthread_spinlock_t input_buffer_lock;


	unsigned short hits_read_start_base1[MAX_HIT_NUMBER];
	unsigned short hits_read_start_base2[MAX_HIT_NUMBER];

	short hits_read_len1[MAX_HIT_NUMBER];
	short hits_read_len2[MAX_HIT_NUMBER];

	long hits_indices1 [MAX_HIT_NUMBER];
	long hits_indices2 [MAX_HIT_NUMBER];
	long decision_table_ids [MAX_HIT_NUMBER];
	unsigned short decision_table_votes [MAX_HIT_NUMBER];
	long decision_table_exon_ids [MAX_HIT_NUMBER];
	char decision_table_read1_used[MAX_HIT_NUMBER];
	char decision_table_read2_used[MAX_HIT_NUMBER];
	long uniq_gene_exonid_table [MAX_HIT_NUMBER];
	long uniq_gene_table [MAX_HIT_NUMBER];
	char * read_coverage_bits;

	char * chro_name_buff;
	z_stream * strm_buffer;

	HashTable * junction_counting_table;   // key: string chro_name \t last_base_previous_exont \t first_base_next_exon
	HashTable * splicing_point_table;
	fc_read_counters read_counters;

	SamBam_Alignment aln_buffer;
} fc_thread_thread_context_t;

#define REVERSE_TABLE_BUCKET_LENGTH 131072
#define REDUCE_TO_5_PRIME_END 5
#define REDUCE_TO_3_PRIME_END 3

typedef struct {
	unsigned int chro_number;
	unsigned int chro_features;
	unsigned int chro_feature_table_start;
	unsigned int chro_block_table_start;
	unsigned int chro_block_table_end;
	unsigned int chro_possible_length;

	unsigned short chro_reverse_table_current_size;
	unsigned int * reverse_table_start_index;
	//unsigned int * reverse_table_end_index;
} fc_chromosome_index_info;

typedef struct {
	int is_gene_level;
	int is_paired_end_input_file;
	int is_paired_end_mode_assign;
	int is_multi_overlap_allowed;
	int is_strand_checked;
	int is_both_end_required;
	int is_chimertc_disallowed;
	int is_PE_distance_checked;
	int is_multi_mapping_allowed;
	int is_input_file_resort_needed;
	int is_SAM_file;
	int is_read_details_out;
	int is_SEPEmix_warning_shown;
	int is_unpaired_warning_shown;
	int is_stake_warning_shown;
	int is_split_alignments_only;
	int is_duplicate_ignored;
	int is_first_read_reversed;
	int is_second_read_straight;
	int do_not_sort;
	int reduce_5_3_ends_to_one;
	int isCVersion;
	int use_fraction_multi_mapping;
	int do_junction_counting;

	int min_mapping_quality_score;
	int min_paired_end_distance;
	int max_paired_end_distance;
	int feature_block_size;
	int read_length;
	int line_length;
	int longest_chro_name;
	int five_end_extension;
	int three_end_extension;
	int fragment_minimum_overlapping;
	int calculate_overlapping_lengths;
	int use_overlapping_break_tie;

	unsigned long long int all_reads;

	unsigned short thread_number;
	fc_thread_thread_context_t * thread_contexts;
	int is_all_finished;
	unsigned int input_buffer_max_size;
	SamBam_Reference_Info * sambam_chro_table;

	char * debug_command;
	char * unistr_buffer_space;
	unsigned int unistr_buffer_size;
	unsigned int unistr_buffer_used;
	HashTable * junction_features_table;
	fasta_contigs_t * fasta_contigs;
	HashTable * gene_name_table;	// gene_name -> gene_number
	HashTable * annot_chro_name_alias_table;	// name in annotation file -> alias name
	char alias_file_name[300];
	char input_file_name[300];
	char raw_input_file_name[300];
	char output_file_name[300];
	unsigned char ** gene_name_array;	// gene_internal_number -> gene_name 

	HashTable * exontable_chro_table;	// gene_name -> fc_chromosome_index_info structure (contains chro_number, feature_number, block_start, block_end, etc) 
	int exontable_nchrs;
	int exontable_exons;
	int * exontable_geneid;
	char * exontable_strand;
	char ** exontable_chr;
	long * exontable_start;
	long * exontable_stop;
	char feature_name_column[100];
	char gene_id_column[100];

	long * exontable_block_end_index;
	long * exontable_block_max_end;
	long * exontable_block_min_start;

	char ** exontable_anno_chrs;
	char * exontable_anno_chr_2ch;
	long * exontable_anno_chr_heads;

	FILE * SAM_output_fp;
	double start_time;

	char * cmd_rebuilt;

	char   redo;

	fc_read_counters read_counters;
	
} fc_thread_global_context_t;

unsigned int tick_time = 1000;


int fetch_boundaries(char * chroname,char * cigar, unsigned int pos, char strand, int *has_left, unsigned short *left_on_read, unsigned int *left_pos, int *has_right, unsigned short *right_on_read, unsigned int *right_pos, fc_junction_info_t *  result_junctions, int junction_space){

	int cigar_cursor = 0, nch, read_len = 0, ret = 0;
	unsigned int chro_cursor = pos, tmpi = 0;
	unsigned int right_boundary = 0;
	unsigned short left_clipped = 0;
	unsigned short right_clipped = 0;
	*has_right = 0;
	*has_left = 0;

	for(; (nch = cigar[cigar_cursor])!=0 ; cigar_cursor++){
		if(isdigit(nch)){
			tmpi = tmpi*10 + (nch - '0');
		} else {
			if (nch == 'S'){
				if(chro_cursor == pos) left_clipped = tmpi;else right_clipped=tmpi;
				read_len += tmpi;
			} else if(nch == 'M' || nch == 'D'){
				if(nch == 'M')read_len += tmpi;

				chro_cursor += tmpi;
				right_boundary = chro_cursor -1;
			} else if(nch == 'N'){
				unsigned int last_exon_last_base = chro_cursor - 1;
				unsigned int next_exon_first_base = chro_cursor + tmpi;
				chro_cursor += tmpi;

				if(ret < junction_space){
					result_junctions[ret].last_exon_base_left = last_exon_last_base;
					result_junctions[ret].first_exon_base_right = next_exon_first_base;
					strcpy(result_junctions[ret].chromosome_name_left, chroname);
					strcpy(result_junctions[ret].chromosome_name_right, chroname);

					ret ++;
				}


			} else if(nch == 'I') read_len += tmpi;
			tmpi = 0;
		}
	}
	if(left_clipped){
		*has_left = 1;
		*left_on_read = left_clipped;
		*left_pos = pos;
	}
	if(right_clipped){
		*has_right = 1;
		*right_on_read = read_len - right_clipped - 1;
		*right_pos = right_boundary;
	}
	return ret;
}

// This function parses the cigar string and returns the number of exon-exon junctions found in the cigar.
// It returns 0 if no junctions are found.
int calc_junctions_from_cigar(fc_thread_global_context_t * global_context, int flag, char * chroname, unsigned int pos, char * cigar , char * extra_tags, fc_junction_info_t * result_junctions){
	unsigned short boundaries_inclusive_base_on_read[FC_CIGAR_PARSER_ITEMS];
	unsigned int boundaries_inclusive_base_pos[FC_CIGAR_PARSER_ITEMS];
	char boundaries_chromosomes[FC_CIGAR_PARSER_ITEMS][MAX_CHROMOSOME_NAME_LEN];
	char boundaries_extend_to_left_on_read[FC_CIGAR_PARSER_ITEMS];
	int boundaries = 0;

	int cigar_cursor = 0, nch, ret = 0, read_len = 0, x1, x2;
	unsigned int chro_cursor = pos, tmpi = 0;
	unsigned int right_boundary = 0;
	unsigned short left_clipped = 0;
	unsigned short right_clipped = 0;

	for(; (nch = cigar[cigar_cursor])!=0 ; cigar_cursor++){
		if(isdigit(nch)){
			tmpi = tmpi*10 + (nch - '0');
		} else {
			if (nch == 'S'){
				if(chro_cursor == pos) left_clipped = tmpi;else right_clipped=tmpi;
				read_len += tmpi;
			} else if(nch == 'M' || nch == 'D'){
				if(nch == 'M')read_len += tmpi;

				chro_cursor += tmpi;
				right_boundary = chro_cursor -1;
			} else if(nch == 'N'){
				unsigned int last_exon_last_base = chro_cursor - 1;
				unsigned int next_exon_first_base = chro_cursor + tmpi;
				if(ret <= FC_CIGAR_PARSER_ITEMS - 1){
					result_junctions[ret].last_exon_base_left = last_exon_last_base;
					result_junctions[ret].first_exon_base_right = next_exon_first_base;
					strcpy(result_junctions[ret].chromosome_name_left, chroname);
					strcpy(result_junctions[ret].chromosome_name_right, chroname);

					ret ++;
				}
				chro_cursor += tmpi;
			} else if(nch == 'I') read_len += tmpi;
			tmpi = 0;
		}
	}
	if(left_clipped){
		strcpy(boundaries_chromosomes[boundaries] , chroname);
		boundaries_extend_to_left_on_read[boundaries] = 0;
		boundaries_inclusive_base_pos[boundaries] = pos;
		boundaries_inclusive_base_on_read[boundaries++] = left_clipped;
	}
	if(right_clipped){
		strcpy(boundaries_chromosomes[boundaries] , chroname);
		boundaries_extend_to_left_on_read[boundaries] = 1;
		boundaries_inclusive_base_pos[boundaries] = chro_cursor - 1;
		boundaries_inclusive_base_on_read[boundaries++] = read_len - right_clipped - 1;
	}

	int tag_cursor=0;

	//if(strstr(extra_tags, "CG:Z")) {
	//	SUBREADprintf("CIGAR=%s, EXTRA=%s\n", cigar, extra_tags);
	//}
	int status = PARSE_STATUS_TAGNAME;
	char tag_name[2], typechar=0;
	int tag_inner_cursor=0;

	char read_main_strand = (((flag & 0x10) == 0x10) == ((flag & 0x40)==0x40))?'-':'+';
	char current_fusion_char[MAX_CHROMOSOME_NAME_LEN];
	unsigned int current_fusion_pos = 0;
	char current_fusion_strand = 0;
	char current_fusion_cigar[FC_CIGAR_PARSER_ITEMS * 15];
	current_fusion_cigar [0] =0;
	current_fusion_char [0]=0;

	while(1){
		int nch = extra_tags[tag_cursor];
		if(status == PARSE_STATUS_TAGNAME){
			tag_name[tag_inner_cursor++] = nch;
			if(tag_inner_cursor == 2){
				status = PARSE_STATUS_TAGTYPE;
				tag_cursor += 1;
				assert(extra_tags[tag_cursor] == ':');
			}
		}else if(status == PARSE_STATUS_TAGTYPE){
			typechar = nch;
			tag_cursor +=1;
			assert(extra_tags[tag_cursor] == ':');
			tag_inner_cursor = 0;
			status = PARSE_STATUS_TAGVALUE;
		}else if(status == PARSE_STATUS_TAGVALUE){
			if(nch == '\t' || nch == 0){
				if(current_fusion_cigar[0] && current_fusion_char[0] && current_fusion_pos && current_fusion_strand){

					unsigned int left_pos = 0, right_pos = 0;
					unsigned short left_on_read = 0, right_on_read = 0;
					int has_left = 0, has_right = 0;

					unsigned int start_pos = current_fusion_pos;
					if(current_fusion_strand!=read_main_strand)
						start_pos = find_left_end_cigar(current_fusion_pos, current_fusion_cigar);

					ret += fetch_boundaries(current_fusion_char, current_fusion_cigar, start_pos, current_fusion_strand, &has_left, &left_on_read, &left_pos, &has_right, &right_on_read, &right_pos, result_junctions + ret, FC_CIGAR_PARSER_ITEMS - ret );

					if(has_left){
						strcpy(boundaries_chromosomes[boundaries] , current_fusion_char);
						boundaries_extend_to_left_on_read[boundaries] = 0;
						boundaries_inclusive_base_pos[boundaries] = left_pos;
						boundaries_inclusive_base_on_read[boundaries++] = left_on_read;
					}
					if(has_right){
						strcpy(boundaries_chromosomes[boundaries] , current_fusion_char);
						boundaries_extend_to_left_on_read[boundaries] = 1;
						boundaries_inclusive_base_pos[boundaries] = right_pos;
						boundaries_inclusive_base_on_read[boundaries++] = right_on_read;
					}
	

			//		SUBREADprintf("BOUND_EXT: %s:%u (at %u) (%c)  ~  %s:%u (at %u) (%c)\n", current_fusion_char, left_pos, left_on_read, has_left?'Y':'X' , current_fusion_char, right_pos, right_on_read,  has_right?'Y':'X');

					current_fusion_pos = 0;
					current_fusion_strand = 0;
					current_fusion_cigar [0] =0;
					current_fusion_char [0]=0;
				}

				tag_inner_cursor = 0;
				status = PARSE_STATUS_TAGNAME;
			}else{
				if(tag_name[0]=='C' && tag_name[1]=='C' && typechar == 'Z'){
					current_fusion_char[tag_inner_cursor++]=nch;
					current_fusion_char[tag_inner_cursor]=0;
				}else if(tag_name[0]=='C' && tag_name[1]=='G' && typechar == 'Z'){
					current_fusion_cigar[tag_inner_cursor++]=nch;
					current_fusion_cigar[tag_inner_cursor]=0;
				}else if(tag_name[0]=='C' && tag_name[1]=='P' && typechar == 'i'){
					current_fusion_pos = current_fusion_pos * 10 + (nch - '0');
				}else if(tag_name[0]=='C' && tag_name[1]=='T' && typechar == 'Z'){
					current_fusion_strand = nch;
				}
			}
		}

		if(nch == 0){
			assert(status == PARSE_STATUS_TAGNAME);
			break;
		}

		tag_cursor++;
	}


	//for(x1 = 0; x1 < boundaries; x1++)
	//	SUBREADprintf("HAS: LR:%d, READ:%d\n", boundaries_extend_to_left_on_read[x1], boundaries_inclusive_base_on_read[x1]);

	for(x1 = 0; x1 < boundaries; x1++)
		for(x2 = 0; x2 < boundaries; x2++){
			if(x1==x2) continue;
			if(boundaries_chromosomes[x1][0]==0 || boundaries_chromosomes[x2][0]==0) continue;
			if(boundaries_extend_to_left_on_read[x1] == 1 && boundaries_extend_to_left_on_read[x2] == 0){
				if( boundaries_inclusive_base_on_read[x1] == boundaries_inclusive_base_on_read[x2]-1 ){

					if(ret <= FC_CIGAR_PARSER_ITEMS - 1){
						result_junctions[ret].last_exon_base_left = boundaries_inclusive_base_pos[x1];
						result_junctions[ret].first_exon_base_right = boundaries_inclusive_base_pos[x2];
						strcpy(result_junctions[ret].chromosome_name_left, boundaries_chromosomes[x1]);
						strcpy(result_junctions[ret].chromosome_name_right, boundaries_chromosomes[x2]);
						ret++;
					}


	//				SUBREADprintf("MATCH: %d ~ %d\n", boundaries_inclusive_base_on_read[x1], boundaries_inclusive_base_on_read[x2]);
					boundaries_chromosomes[x1][0]=0;
					boundaries_chromosomes[x2][0]=0;
				}
			}
		}

	//for(x1 = 0; x1 < boundaries; x1++)
	//	if(boundaries_chromosomes[x1][0])
	//		SUBREADprintf("LEFT: LR:%d, READ:%d\n", boundaries_extend_to_left_on_read[x1], boundaries_inclusive_base_on_read[x1]);
	return ret;
}


unsigned int unistr_cpy(fc_thread_global_context_t * global_context, char * str, int strl)
{
	unsigned int ret;
	if(global_context->unistr_buffer_used + strl >= global_context->unistr_buffer_size-1)
	{
		if( global_context->unistr_buffer_size < 3435973835u) // 4G / 5 * 4 - 5
		{
			global_context -> unistr_buffer_size = global_context->unistr_buffer_size /4 *5;
			global_context -> unistr_buffer_space = realloc(global_context -> unistr_buffer_space, global_context->unistr_buffer_size);
		}
		else
		{
			SUBREADprintf("Error: exceed memory limit (4GB) for storing annotation data.\n");
			return 0xffffffffu;
		}
	}

	strcpy(global_context -> unistr_buffer_space + global_context->unistr_buffer_used, str);
	ret = global_context->unistr_buffer_used;

	global_context->unistr_buffer_used += strl +1;

	return ret;
}

void print_FC_configuration(fc_thread_global_context_t * global_context, char * annot, char * sam, char * out, int is_sam, int is_GTF, int *n_input_files, int isReadSummaryReport)
{
	char * tmp_ptr1 = NULL , * next_fn, *sam_used = malloc(strlen(sam)+1), sam_ntxt[30],bam_ntxt[30], next_ntxt[50];
	int nfiles=1, nBAMfiles = 0, nNonExistFiles = 0;

	strcpy(sam_used, sam);

	SUBREADputs("");
	print_subread_logo();
	SUBREADputs("");
	print_in_box(80,1,1,"featureCounts setting");
	print_in_box(80,0,0,"");
	
	nfiles = 0;

	while(1)
	{
		next_fn = strtok_r(nfiles==0?sam_used:NULL, ";", &tmp_ptr1);
		if(next_fn == NULL || strlen(next_fn)<1) break;
		nfiles++;

		int file_probe = is_certainly_bam_file(next_fn, NULL);
		if(file_probe==-1) nNonExistFiles++;
		if(file_probe == 1) nBAMfiles++;		
	}

	sam_ntxt[0]=0;
	bam_ntxt[0]=0;
	next_ntxt[0]=0;

	if(nNonExistFiles)
		sprintf(next_ntxt, "%d unknown file%s", nNonExistFiles, nNonExistFiles>1?"s":"");
	if(nBAMfiles)
		sprintf(bam_ntxt, "%d BAM file%s  ", nBAMfiles, nBAMfiles>1?"s":"");
	if(nfiles-nNonExistFiles-nBAMfiles)
		sprintf(sam_ntxt, "%d SAM file%s  ", nfiles-nNonExistFiles-nBAMfiles , (nfiles-nNonExistFiles-nBAMfiles)>1?"s":"");


	strcpy(sam_used, sam);

	print_in_box(80,0,0,"            Input files : %s%s%s", sam_ntxt, bam_ntxt, next_ntxt);
	nfiles=0;

	while(1)
	{
		next_fn = strtok_r(nfiles==0?sam_used:NULL, ";", &tmp_ptr1);
		if(next_fn == NULL || strlen(next_fn)<1) break;
		int is_first_read_PE = 0 , file_probe = is_certainly_bam_file(next_fn, &is_first_read_PE);

		char file_chr = 'S';
		if(file_probe == -1) file_chr = '?';
		else if(is_first_read_PE == 1) file_chr = 'P';
		//file_chr = 'o';

		print_in_box(94,0,0,"                          %c[32m%c%c[36m %s%c[0m",CHAR_ESC, file_chr,CHAR_ESC, next_fn,CHAR_ESC);
		nfiles++;
	}

	(*n_input_files) = nfiles;
	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"            Output file : %s", out);
	print_in_box(80,0,0,"            Annotations : %s (%s)", annot, is_GTF?"GTF":"SAF");
	if(isReadSummaryReport)
		print_in_box(80,0,0,"     Assignment details : <input_file>.featureCounts");

	if(global_context -> alias_file_name[0])
		print_in_box(80,0,0,"  Chromosome alias file : %s", global_context -> alias_file_name);

	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"                Threads : %d", global_context->thread_number);
	print_in_box(80,0,0,"                  Level : %s level", global_context->is_gene_level?"meta-feature":"feature");
	print_in_box(80,0,0,"             Paired-end : %s", global_context->is_paired_end_mode_assign?"yes":"no");
	if(global_context -> do_not_sort && global_context->is_paired_end_mode_assign)
	{
		print_in_box(80,0,0,"       Sorting PE Reads : never");
		print_in_box(80,0,0,"");
	}

	print_in_box(80,0,0,"        Strand specific : %s", global_context->is_strand_checked?(global_context->is_strand_checked==1?"yes":"inversed"):"no");
	char * multi_mapping_allow_mode = "not counted";
	if(global_context->is_multi_mapping_allowed == ALLOW_PRIMARY_MAPPING)
		multi_mapping_allow_mode = "primary only";
	else if(global_context->is_multi_mapping_allowed == ALLOW_ALL_MULTI_MAPPING)
		multi_mapping_allow_mode = global_context -> use_fraction_multi_mapping?"counted (as fractions)": "counted (as integer)";

	print_in_box(80,0,0,"     Multimapping reads : %s", multi_mapping_allow_mode);
	print_in_box(80,0,0,"Multi-overlapping reads : %s", global_context->is_multi_overlap_allowed?"counted":"not counted");
	if(global_context -> is_split_alignments_only)
		print_in_box(80,0,0,"       Split alignments : required");
	if(global_context -> fragment_minimum_overlapping !=1)
		print_in_box(80,0,0,"      Overlapping bases : %d", global_context -> fragment_minimum_overlapping);
	if(global_context -> five_end_extension || global_context -> three_end_extension)
		print_in_box(80,0,0,"        Read extensions : %d on 5' and %d on 3' ends", global_context -> five_end_extension , global_context -> three_end_extension);
	if(global_context -> reduce_5_3_ends_to_one)
		print_in_box(80,0,0,"      Read reduction to : %d' end" , global_context -> reduce_5_3_ends_to_one == REDUCE_TO_5_PRIME_END ?5:3);
	if(global_context -> is_duplicate_ignored)
		print_in_box(80,0,0,"       Duplicated Reads : ignored");
	print_in_box(80,0,0,"      Read orientations : %c%c", global_context->is_first_read_reversed?'r':'f', global_context->is_second_read_straight?'f':'r' );

	if(global_context->is_paired_end_mode_assign)
	{
		print_in_box(80,0,0,"");
		print_in_box(80,0,0,"         Chimeric reads : %s", global_context->is_chimertc_disallowed?"not counted":"counted");
		print_in_box(80,0,0,"       Both ends mapped : %s", global_context->is_both_end_required?"required":"not required");

		if(global_context->is_PE_distance_checked)
			print_in_box(80,0,0,"        Fragment length : %d - %d", global_context -> min_paired_end_distance, global_context -> max_paired_end_distance);
	}

	print_in_box(80,0,0,"");
	print_in_box(80,2,1,"http://subread.sourceforge.net/");
	SUBREADputs("");
	print_in_box(80,1,1,"Running");
	print_in_box(80,0,0,"");
	if(global_context->annot_chro_name_alias_table)
		print_in_box(80,0,0,"%ld chromosome name aliases are loaded.", global_context -> annot_chro_name_alias_table ->numOfElements);

	free(sam_used);
}

void print_FC_results(fc_thread_global_context_t * global_context)
{
	print_in_box(89,0,1,"%c[36mRead assignment finished.%c[0m", CHAR_ESC, CHAR_ESC);
	print_in_box(80,0,0,"");
	print_in_box(80,2,1,"http://subread.sourceforge.net/");
	SUBREADputs("");
	return;


	if(0){
		print_in_box(80,1,1,"Summary");
		print_in_box(80,0,0,"");
		if(global_context->is_paired_end_mode_assign)
			print_in_box(80,0,0,"        All fragments : %llu", global_context -> all_reads);
		else
			print_in_box(80,0,0,"            All reads : %llu", global_context -> all_reads);

		if(global_context->is_gene_level)
			print_in_box(80,0,0,"        Meta-features : %lu", global_context -> gene_name_table -> numOfElements);
		else
			print_in_box(80,0,0,"             Features : %u", global_context -> exontable_exons);

		if(global_context->is_paired_end_mode_assign)
			print_in_box(80,0,0,"   Assigned fragments : %llu", global_context -> read_counters.assigned_reads);
		else
			print_in_box(80,0,0,"       Assigned reads : %llu", global_context -> read_counters.assigned_reads);

		print_in_box(80,0,0,"            Time cost : %.3f minutes", (miltime() - global_context -> start_time)/60);
		print_in_box(80,0,0,"");
		print_in_box(80,2,1,"http://subread.sourceforge.net/");
	}
	SUBREADputs("");
}

int fc_strcmp(const void * s1, const void * s2)
{
	return strcmp((char*)s1, (char*)s2);
}


int is_comment_line(const char * l, int file_type, unsigned int lineno)
{
	int tabs = 0, xk1 = 0;
	if(l[0]=='#') return 1;

	if(isalpha(l[0]) && file_type == FILE_TYPE_RSUBREAD)
	{
		char target_chr[16];
		memcpy(target_chr, l, 16);
		for(xk1=0; xk1<16; xk1++)
			target_chr[xk1] = tolower(target_chr[xk1]);

		if(memcmp(target_chr, "geneid\tchr\tstart",16)==0) return 1;
	}

	xk1=0;
	while(l[xk1]) tabs += (l[xk1++] == '\t');

	return tabs < ((file_type == FILE_TYPE_GTF)?8:4);
}

void register_junc_feature(fc_thread_global_context_t *global_context, char * feature_name, char * chro, unsigned int start, unsigned int stop){
	HashTable * gene_table = HashTableGet(global_context -> junction_features_table, chro);
	//SUBREADprintf("REG %s : %p\n", chro, gene_table);
	if(NULL == gene_table){
		gene_table = HashTableCreate(48367);
		HashTableSetDeallocationFunctions(gene_table, NULL, free);
		HashTableSetKeyComparisonFunction(gene_table, fc_strcmp);
		HashTableSetHashFunction(gene_table, fc_chro_hash);

		char * new_name = malloc(strlen(chro)+1);
		strcpy(new_name, chro);
		HashTablePut(global_context -> junction_features_table, new_name, gene_table);
	}
	fc_junction_gene_t * gene_info = HashTableGet(gene_table, feature_name);
	if(NULL == gene_info){
		gene_info = malloc(sizeof(fc_junction_gene_t));
		strcpy(gene_info -> gene_name, feature_name);
		gene_info -> pos_first_base = start;
		gene_info -> pos_last_base = stop;

		HashTablePut(gene_table, gene_info -> gene_name, gene_info);
	}else{
		gene_info -> pos_first_base = min(start, gene_info -> pos_first_base);
		gene_info -> pos_last_base = max(stop, gene_info -> pos_last_base);
	}
}

int locate_junc_features(fc_thread_global_context_t *global_context, char * chro, unsigned int pos, fc_junction_gene_t ** ret_info, int max_ret_info_size){
	HashTable * gene_table = NULL;

	if(global_context -> annot_chro_name_alias_table) {
		char * anno_chro_name = HashTableGet( global_context -> annot_chro_name_alias_table , chro);
		if(anno_chro_name)
			gene_table = HashTableGet( global_context -> junction_features_table , anno_chro_name);
	}
	if(gene_table == NULL)
		gene_table = HashTableGet(global_context -> junction_features_table, chro);

	if(gene_table == NULL && strlen(chro)>3 && memcmp(chro, "chr", 3)==0){
		gene_table = HashTableGet(global_context -> junction_features_table, chro + 3);
	}

	if(gene_table == NULL){
		char new_name [FEATURE_NAME_LENGTH];

		strcpy(new_name, "chr");
		strcat(new_name, chro);
		gene_table = HashTableGet(global_context -> junction_features_table, new_name);
	}


	//SUBREADprintf("GET %s : %p\n", chro, gene_table);
	if(gene_table == NULL) return 0;

	int bucket, ret = 0;
	KeyValuePair * cursor;
	for(bucket = 0; bucket < gene_table -> numOfBuckets; bucket++){
		cursor = gene_table -> bucketArray[bucket];
		while(cursor){
			fc_junction_gene_t * gene_info = cursor -> value;
			//SUBREADprintf("COMP: %u ~ %u:%u\n", pos, gene_info -> pos_first_base , gene_info -> pos_last_base);
			if(gene_info -> pos_first_base <= pos && gene_info -> pos_last_base >= pos){
				if(ret < max_ret_info_size)
					ret_info [ret ++] = gene_info;
			}
			cursor = cursor -> next;
		}
	}
	return ret;
}

// This function loads annotations from the file.
// It returns the number of featres loaded, or -1 if something is wrong. 
// Memory will be allowcated in this function. The pointer is saved in *loaded_features.
// The invoker must release the memory itself.

int load_feature_info(fc_thread_global_context_t *global_context, const char * annotation_file, int file_type, fc_feature_info_t ** loaded_features)
{
	unsigned int features = 0, xk1 = 0, lineno=0;
	char * file_line = malloc(MAX_LINE_LENGTH+1);
	FILE * fp = f_subr_open(annotation_file,"r"); 
	int is_GFF_warned = 0;
	if(!fp) return -1;

	HashTable * chro_name_table = HashTableCreate(1603);
	HashTableSetHashFunction(chro_name_table, fc_chro_hash);
	HashTableSetKeyComparisonFunction(chro_name_table, fc_strcmp_chro);
	global_context -> longest_chro_name = 0;

	if(global_context -> do_junction_counting){
		global_context -> junction_features_table = HashTableCreate(1603);
		HashTableSetDeallocationFunctions(global_context -> junction_features_table, free, (void (*)(void *))HashTableDestroy);
		HashTableSetKeyComparisonFunction(global_context -> junction_features_table, fc_strcmp);
		HashTableSetHashFunction(global_context -> junction_features_table, fc_chro_hash);
	}


	// first scan: get the chromosome size, etc
	while(1)
	{
		char * fgets_ret = fgets(file_line, MAX_LINE_LENGTH, fp);
		char * token_temp, *chro_name;
		fc_chromosome_index_info * chro_stab;
		unsigned int feature_pos = 0;
		if(!fgets_ret) break;

		lineno++;
		if(is_comment_line(file_line, file_type, lineno-1))continue;
		if(file_type == FILE_TYPE_GTF)
		{
			chro_name = strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp); // lib_name (not needed)
			char * feature_type = strtok_r(NULL,"\t", &token_temp);
			if(strcmp(feature_type, global_context -> feature_name_column)==0)
			{
				strtok_r(NULL,"\t", &token_temp); // feature_start
				feature_pos = atoi(strtok_r(NULL,"\t", &token_temp));// feature_end
				features++;
			}
			else chro_name = NULL;
		}
		else
		{
			strtok_r(file_line,"\t", &token_temp);
			chro_name = strtok_r(NULL,"\t",&token_temp);
			strtok_r(NULL,"\t",&token_temp);	// feature_start
			feature_pos = atoi(strtok_r(NULL,"\t", &token_temp));// feature_end

			features++;
		}

		if(chro_name)
		{
			if(strlen(chro_name)>=CHROMOSOME_NAME_LENGTH) 
				chro_name[CHROMOSOME_NAME_LENGTH-1]=0;
			chro_stab = HashTableGet(chro_name_table, chro_name);

			if(chro_stab)
			{
				chro_stab -> chro_possible_length = max(chro_stab -> chro_possible_length , feature_pos+1);
			}else
			{
				char * tmp_chro_name = malloc(CHROMOSOME_NAME_LENGTH);
				term_strncpy(tmp_chro_name, chro_name, CHROMOSOME_NAME_LENGTH);
				chro_stab = calloc(sizeof(fc_chromosome_index_info),1);
				chro_stab -> chro_number = chro_name_table->numOfElements;
				chro_stab -> chro_possible_length = feature_pos+1;
				chro_stab -> reverse_table_start_index = NULL;
				HashTablePut(chro_name_table, tmp_chro_name, chro_stab);
			}

			chro_stab -> chro_features ++;
		}
	}

	fseek(fp,0,SEEK_SET);

	fc_feature_info_t * ret_features = malloc(sizeof(fc_feature_info_t) * features);

	lineno = 0;
	while(xk1 < features)
	{
		int is_gene_id_found = 0;
		fgets(file_line, MAX_LINE_LENGTH, fp);
		lineno++;
		char * token_temp;
		if(is_comment_line(file_line, file_type, lineno-1))continue;

		if(file_type == FILE_TYPE_RSUBREAD)
		{
			char * feature_name = strtok_r(file_line,"\t",&token_temp);
			int feature_name_len = strlen(feature_name);
			if(feature_name_len > FEATURE_NAME_LENGTH) feature_name[FEATURE_NAME_LENGTH -1 ] = 0;
			ret_features[xk1].feature_name_pos = unistr_cpy(global_context, (char *)feature_name, feature_name_len);

			char * seq_name = strtok_r(NULL,"\t", &token_temp);
			int chro_name_len = strlen(seq_name);
			if(chro_name_len > CHROMOSOME_NAME_LENGTH) seq_name[CHROMOSOME_NAME_LENGTH -1 ] = 0;
			unsigned int chro_name_pos = unistr_cpy(global_context, (char *)seq_name, chro_name_len);
			global_context -> longest_chro_name = max(chro_name_len, global_context -> longest_chro_name);

			ret_features[xk1].chro_name_pos_delta = chro_name_pos - ret_features[xk1].feature_name_pos;
			ret_features[xk1].start = atoi(strtok_r(NULL,"\t", &token_temp));// start 
			if(ret_features[xk1].start<0 || ret_features[xk1].start>0x7fffffff)
			{
				ret_features[xk1].start = 0;
				print_in_box(80,0,0,"WARNING the %d-th line has a negative chro coordinate.", lineno);
			}

			ret_features[xk1].end = atoi(strtok_r(NULL,"\t", &token_temp));//end 
			if(ret_features[xk1].end<0 || ret_features[xk1].end>0x7fffffff)
			{
				ret_features[xk1].end = 0;
				print_in_box(80,0,0,"WARNING the %d-th line has a negative chro coordinate.", lineno);
			}




			char * strand_str = strtok_r(NULL,"\t", &token_temp); 
			if(strand_str == NULL)
				ret_features[xk1].is_negative_strand = 0;
			else
				ret_features[xk1].is_negative_strand = ('-' ==strand_str[0]);
			ret_features[xk1].sorted_order = xk1;

			int bin_location = ret_features[xk1].start / REVERSE_TABLE_BUCKET_LENGTH;
			
			fc_chromosome_index_info * chro_stab = HashTableGet(chro_name_table, seq_name);
			if(!chro_stab -> reverse_table_start_index)
			{
				chro_stab -> reverse_table_start_index = malloc(sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
				memset(chro_stab -> reverse_table_start_index, 0 , sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
			}
			chro_stab -> reverse_table_start_index[bin_location]++;

			is_gene_id_found = 1;

			assert(feature_name);
			if(global_context -> do_junction_counting)
				register_junc_feature(global_context , feature_name, seq_name, ret_features[xk1].start, ret_features[xk1].end);

			xk1++;
		}
		else if(file_type == FILE_TYPE_GTF)
		{
			char feature_name_tmp[FEATURE_NAME_LENGTH];
			sprintf(feature_name_tmp, "LINE_%07u", xk1 + 1);
			char * seq_name = strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp);// source
			char * feature_type = strtok_r(NULL,"\t", &token_temp);// feature_type
			if(strcmp(feature_type, global_context -> feature_name_column)==0)
			{
				ret_features[xk1].start = atoi(strtok_r(NULL,"\t", &token_temp));// start 
				ret_features[xk1].end = atoi(strtok_r(NULL,"\t", &token_temp));//end 

				if(ret_features[xk1].start < 1 || ret_features[xk1].end<1 ||  ret_features[xk1].start > 0x7fffffff ||  ret_features[xk1].end > 0x7fffffff || ret_features[xk1].start > ret_features[xk1].end)
					SUBREADprintf("\nWarning: the feature on the %d-th line has zero coordinate or zero lengths\n\n", lineno);


				strtok_r(NULL,"\t", &token_temp);// score 
				ret_features[xk1].is_negative_strand = ('-' == (strtok_r(NULL,"\t", &token_temp)[0]));//strand 
				ret_features[xk1].sorted_order = xk1;
				strtok_r(NULL,"\t",&token_temp);	// "frame"
				char * extra_attrs = strtok_r(NULL,"\t",&token_temp);	// name_1 "val1"; name_2 "val2"; ... 
				if(extra_attrs && (strlen(extra_attrs)>2))
				{
					int attr_val_len = GTF_extra_column_value(extra_attrs , global_context -> gene_id_column , feature_name_tmp, FEATURE_NAME_LENGTH);
					if(attr_val_len>0) is_gene_id_found=1;
			//		printf("V=%s\tR=%d\n", extra_attrs , attr_val_len);
				}

				if(is_gene_id_found)
				{
				}
				else
				{
					if(!is_GFF_warned)
					{
						int ext_att_len = strlen(extra_attrs);
						if(extra_attrs[ext_att_len-1] == '\n') extra_attrs[ext_att_len-1] =0;
						SUBREADprintf("\nWarning: failed to find the gene identifier attribute in the 9th column of the provided GTF file.\nThe specified gene identifier attribute is '%s' \nThe attributes included in your GTF annotation are '%s' \n\n",  global_context -> gene_id_column, extra_attrs);
					}
					is_GFF_warned++;
				}

				int feature_name_len = strlen(feature_name_tmp);
				if(feature_name_len > FEATURE_NAME_LENGTH) feature_name_tmp[FEATURE_NAME_LENGTH -1 ] = 0;
				ret_features[xk1].feature_name_pos = unistr_cpy(global_context, (char *)feature_name_tmp, feature_name_len);

				int chro_name_len = strlen(seq_name);
				if(chro_name_len > CHROMOSOME_NAME_LENGTH) seq_name[CHROMOSOME_NAME_LENGTH -1 ] = 0;
				unsigned int chro_name_pos = unistr_cpy(global_context, (char *)seq_name, chro_name_len);
				global_context -> longest_chro_name = max(chro_name_len, global_context -> longest_chro_name);

				ret_features[xk1].chro_name_pos_delta = chro_name_pos - ret_features[xk1].feature_name_pos;

				int bin_location = ret_features[xk1].start / REVERSE_TABLE_BUCKET_LENGTH;
				fc_chromosome_index_info * chro_stab = HashTableGet(chro_name_table, seq_name);
				if(!chro_stab -> reverse_table_start_index)
				{
					chro_stab -> reverse_table_start_index = malloc(sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
					memset(chro_stab -> reverse_table_start_index, 0 , sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
				}
				chro_stab -> reverse_table_start_index[bin_location]++;

				if(global_context -> do_junction_counting)
					register_junc_feature(global_context , feature_name_tmp, seq_name, ret_features[xk1].start, ret_features[xk1].end);

				xk1++;
			}
		}
	}
	fclose(fp);
	free(file_line);

	(*loaded_features) = ret_features;
	global_context -> exontable_nchrs = (int)chro_name_table-> numOfElements;
	global_context -> exontable_chro_table = chro_name_table;

	print_in_box(80,0,0,"   Features : %d\n", features);
	if(features < 1)
	{
		print_in_box(80,0,0,"WARNING no features were loaded in format %s.", file_type == FILE_TYPE_GTF?"GTF":"SAF");
		print_in_box(80,0,0,"        annotation format can be specified using '-F'.");
	}
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

void register_reverse_table(int block_no, long this_block_min_start, long this_block_max_end, fc_chromosome_index_info * chro_inf)
{

	unsigned int reversed_bucket_start = this_block_min_start /  REVERSE_TABLE_BUCKET_LENGTH;
	unsigned int reversed_bucket_end = this_block_max_end / REVERSE_TABLE_BUCKET_LENGTH;
	assert(this_block_min_start <= this_block_max_end);
	assert(reversed_bucket_end < chro_inf -> chro_possible_length);
	int x1;
	for(x1 = reversed_bucket_start; x1 <= reversed_bucket_end; x1++)
	{
		chro_inf->reverse_table_start_index[x1] = min(chro_inf->reverse_table_start_index[x1], block_no);
		//chro_inf->reverse_table_end_index[x1] = max(chro_inf->reverse_table_end_index[x1], block_no+1);
	}

}

void feature_merge(void * arrv, int start, int items, int items2)
{

	void ** arr = (void **) arrv;

	long * ret_start = (long *) arr[0];
	long * ret_end = (long *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	int total_items = items+items2;
	long * tmp_start = malloc(sizeof(long) * total_items);
	long * tmp_end = malloc(sizeof(long) * total_items);
	unsigned char * tmp_strand = malloc(sizeof(char) * total_items);
	int * tmp_entyrez = malloc(sizeof(int) * total_items);
	fc_feature_info_t ** tmp_info_ptr = malloc(sizeof(fc_feature_info_t*) * total_items);

	int read_1_ptr = start;
	int read_2_ptr = start+items;
	int write_ptr;

	for(write_ptr=0; write_ptr<total_items; write_ptr++)
	{
		if((read_1_ptr >= start+items)||(read_2_ptr < start+total_items && ret_start[read_1_ptr] >= ret_start[read_2_ptr]))
		{
			tmp_start[write_ptr] = ret_start[read_2_ptr];
			tmp_end[write_ptr] = ret_end[read_2_ptr];
			tmp_strand[write_ptr] = ret_strand[read_2_ptr];
			tmp_entyrez[write_ptr] = ret_entyrez[read_2_ptr];
			tmp_info_ptr[write_ptr] = old_info_ptr[read_2_ptr];
			read_2_ptr++;
		}
		else
		{
			tmp_start[write_ptr] = ret_start[read_1_ptr];
			tmp_end[write_ptr] = ret_end[read_1_ptr];
			tmp_strand[write_ptr] = ret_strand[read_1_ptr];
			tmp_entyrez[write_ptr] = ret_entyrez[read_1_ptr];
			tmp_info_ptr[write_ptr] = old_info_ptr[read_1_ptr];
			read_1_ptr++;
		}
	}

	memcpy(ret_start+ start, tmp_start, sizeof(long) * total_items);
	memcpy(ret_end+ start, tmp_end, sizeof(long) * total_items);
	memcpy(ret_strand+ start, tmp_strand, sizeof(char) * total_items);
	memcpy(ret_entyrez+ start, tmp_entyrez, sizeof(int) * total_items);
	memcpy(old_info_ptr+ start, tmp_info_ptr, sizeof(fc_feature_info_t*) * total_items);

	free(tmp_start);
	free(tmp_end);
	free(tmp_strand);
	free(tmp_entyrez);
	free(tmp_info_ptr);
}


int feature_sort_compare(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	long * ret_start = (long *)arr[0];
	long ll = ret_start[l];
	long rl = ret_start[r];

	if(ll==rl) return 0;
	else if(ll>rl) return 1;
	else return -1;
}

void feature_sort_exchange(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	long tmp;
	fc_feature_info_t * tmpptr;

	long * ret_start = (long *) arr[0];
	long * ret_end = (long *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	
	tmp = ret_start[r];
	ret_start[r]=ret_start[l];
	ret_start[l]=tmp;

	tmp = ret_end[r];
	ret_end[r]=ret_end[l];
	ret_end[l]=tmp;

	tmp = ret_strand[r];
	ret_strand[r]=ret_strand[l];
	ret_strand[l]=tmp;

	tmp = ret_entyrez[r];
	ret_entyrez[r]=ret_entyrez[l];
	ret_entyrez[l]=tmp;

	tmpptr = old_info_ptr[r];
	old_info_ptr[r]=old_info_ptr[l];
	old_info_ptr[l]=tmpptr;

}



void sort_feature_info(fc_thread_global_context_t * global_context, unsigned int features, fc_feature_info_t * loaded_features, char *** sorted_chr_names, int ** sorted_entrezid, long ** sorted_start, long ** sorted_end, unsigned char ** sorted_strand, char ** anno_chr_2ch, char *** anno_chrs, long ** anno_chr_head, long ** block_end_index, long ** block_min_start_pos, long ** block_max_end_pos)
{
	unsigned int chro_pnt;
	unsigned int xk1,xk2;
	int * ret_entrez = malloc(sizeof(int) * features);
	long * ret_start = malloc(sizeof(long) * features);
	long * ret_end = malloc(sizeof(long) * features);
	int current_block_buffer_size = 2000;

	long * ret_block_end_index = malloc(sizeof(long) * current_block_buffer_size);
	long * ret_block_min_start = malloc(sizeof(long) * current_block_buffer_size);
	long * ret_block_max_end = malloc(sizeof(long) * current_block_buffer_size);
	unsigned char * ret_strand = malloc(features);
	char ** ret_char_name = malloc(sizeof(void *) * features);
	fc_feature_info_t ** old_info_ptr = malloc(sizeof(void *) * features);
	(*anno_chrs) = malloc(sizeof(void *) * global_context -> exontable_nchrs);
	(*anno_chr_head) = malloc(sizeof(long) * global_context -> exontable_nchrs);
	(*anno_chr_2ch) = malloc(sizeof(char) * global_context -> exontable_nchrs*2); 
	unsigned int * chro_feature_ptr = calloc(sizeof(int) * global_context -> exontable_nchrs,1);
	fc_chromosome_index_info ** tmp_chro_info_ptrs = malloc(global_context -> exontable_nchrs * sizeof(fc_chromosome_index_info *));

	global_context -> gene_name_array = malloc(sizeof(char *) * features);	// there should be much less identical names.
	global_context -> gene_name_table = HashTableCreate(5000);
	HashTableSetHashFunction(global_context -> gene_name_table, HashTableStringHashFunction);
	HashTableSetKeyComparisonFunction(global_context -> gene_name_table, fc_strcmp);

	// init start positions of each chromosome block.
	if(1)
	{
		KeyValuePair * cursor;
		int bucket;
		unsigned int sum_ptr = 0;
		for(bucket=0; bucket < global_context -> exontable_chro_table  -> numOfBuckets; bucket++)
		{
			cursor = global_context -> exontable_chro_table -> bucketArray[bucket];
			while (1)
			{
				if (!cursor) break;
				fc_chromosome_index_info * tmp_chro_inf = cursor -> value;
				cursor = cursor->next;
				//tmp_chro_inf -> reverse_table_end_index = calloc(sizeof(int), tmp_chro_inf->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2);
				chro_feature_ptr [tmp_chro_inf -> chro_number] = tmp_chro_inf -> chro_features;
				tmp_chro_info_ptrs[tmp_chro_inf -> chro_number] = tmp_chro_inf;
			}
		}

		for(xk1 = 0; xk1 < global_context -> exontable_nchrs; xk1++)
		{
			unsigned int tmpv = chro_feature_ptr[xk1];
			chro_feature_ptr[xk1] = sum_ptr;
			tmp_chro_info_ptrs[xk1] -> chro_feature_table_start = sum_ptr;
		//		printf("SII=%u  +  %u\n", sum_ptr, tmpv);
			sum_ptr += tmpv;
		}

	}
	int current_block_id = 0, sort_i = 0;

	(*sorted_chr_names) = ret_char_name;
	(*sorted_entrezid) = ret_entrez;
	(*sorted_start) = ret_start;
	(*sorted_end) = ret_end;
	(*sorted_strand) = ret_strand;
	int curr_chro_number = 0;

	for(chro_pnt=0; chro_pnt < features; chro_pnt++)
	{
		char * this_chro_name = global_context -> unistr_buffer_space + loaded_features[chro_pnt].feature_name_pos + loaded_features[chro_pnt].chro_name_pos_delta;
		fc_chromosome_index_info * this_chro_info = HashTableGet(global_context -> exontable_chro_table , this_chro_name);
		assert(this_chro_info);
		unsigned int this_chro_number = this_chro_info -> chro_number;
		unsigned int this_chro_table_ptr = chro_feature_ptr[this_chro_number];

		ret_char_name[this_chro_table_ptr] = this_chro_name;// (char *)loaded_features[chro_pnt].chro;
		ret_entrez[this_chro_table_ptr] = find_or_insert_gene_name(global_context, (unsigned char *)(global_context -> unistr_buffer_space + loaded_features[chro_pnt].feature_name_pos));
		ret_start[this_chro_table_ptr] = loaded_features[chro_pnt].start;
		ret_end[this_chro_table_ptr] = loaded_features[chro_pnt].end;
		ret_strand[this_chro_table_ptr] = loaded_features[chro_pnt].is_negative_strand;
		old_info_ptr[this_chro_table_ptr] = &loaded_features[chro_pnt];

		chro_feature_ptr[this_chro_number]++;
	}

	for(xk1 = 0; xk1 < global_context -> exontable_nchrs; xk1++)
	{
		fc_chromosome_index_info * tmp_chro_inf = tmp_chro_info_ptrs[xk1];
		int bins_in_chr = ( tmp_chro_inf->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2);
		short * features_per_block_bins = malloc(sizeof(short)*bins_in_chr);
		for(xk2=0; xk2<bins_in_chr; xk2++)
		{
			features_per_block_bins[xk2] = max(1,min(1000,(int)(0.9999999+sqrt(tmp_chro_inf -> reverse_table_start_index[xk2]))));
			//printf("CHR%d : SQR[%d]=%d (%d)\n",  tmp_chro_inf -> chro_number,xk2, features_per_block_bins[xk2], tmp_chro_inf -> reverse_table_start_index[xk2] );
		}

		memset(tmp_chro_inf -> reverse_table_start_index, 0xff, sizeof(int) *bins_in_chr);

		tmp_chro_inf -> chro_block_table_start = current_block_id; 
		unsigned int this_block_items = 0;
		long this_block_min_start = 0x7fffffff, this_block_max_end = 0;
		unsigned int this_chro_tab_end =  tmp_chro_inf -> chro_features + tmp_chro_inf -> chro_feature_table_start;

		void * in_array[5];
		in_array[0] = ret_start + tmp_chro_inf -> chro_feature_table_start; 
		in_array[1] = ret_end + tmp_chro_inf -> chro_feature_table_start; 
		in_array[2] = ret_strand + tmp_chro_inf -> chro_feature_table_start; 
		in_array[3] = ret_entrez + tmp_chro_inf -> chro_feature_table_start; 
		in_array[4] = old_info_ptr + tmp_chro_inf -> chro_feature_table_start; 

		merge_sort(in_array, this_chro_tab_end - tmp_chro_inf -> chro_feature_table_start, feature_sort_compare, feature_sort_exchange, feature_merge);

		for(sort_i = tmp_chro_inf -> chro_feature_table_start; sort_i< this_chro_tab_end ; sort_i++)
		{
			// NOW THE FEATURES (ret_start, ret_end, ret_strand, ret_entrez, old_info_ptr) ARE ALL SORTED!
			//printf("NT=%lu\tCHRO=%d\n", ret_start[sort_i], tmp_chro_inf->chro_number);
			old_info_ptr[sort_i]->sorted_order = sort_i;

			int feature_bin_location = ret_start[sort_i] / REVERSE_TABLE_BUCKET_LENGTH;
			int block_bin_location = this_block_min_start / REVERSE_TABLE_BUCKET_LENGTH;

			if(this_block_items && (this_block_items > features_per_block_bins[block_bin_location] || feature_bin_location != block_bin_location))//global_context -> feature_block_size)
			{

				if(current_block_id >= current_block_buffer_size - 1)
				{
					current_block_buffer_size *= 1.3;
					ret_block_min_start = realloc(ret_block_min_start, sizeof(long)*current_block_buffer_size);
					ret_block_max_end = realloc(ret_block_max_end, sizeof(long)*current_block_buffer_size);
					ret_block_end_index = realloc(ret_block_end_index, sizeof(long)*current_block_buffer_size);
				}


				ret_block_end_index[current_block_id] = sort_i;	// FIRST UNWANTED ID
				ret_block_min_start[current_block_id] = this_block_min_start;
				ret_block_max_end[current_block_id] = this_block_max_end;
				register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
				//printf("B=%d; ST=%ld, END=%ld, ITM=%d\n", current_block_id, this_block_min_start, this_block_max_end, this_block_items);
				current_block_id++;
				this_block_max_end = 0;
				this_block_items = 0;
				this_block_min_start = 0x7fffffff;
			}

			this_block_max_end = max(this_block_max_end, ret_end[sort_i]);
			this_block_min_start = min(this_block_min_start, ret_start[sort_i]);
			this_block_items ++;
		
		}
		if(this_block_items)
		{
			if(current_block_id >= current_block_buffer_size)
			{
				current_block_buffer_size *= 1.3;
				ret_block_min_start = realloc(ret_block_min_start, sizeof(long)*current_block_buffer_size);
				ret_block_max_end = realloc(ret_block_max_end, sizeof(long)*current_block_buffer_size);
				ret_block_end_index = realloc(ret_block_end_index, sizeof(long)*current_block_buffer_size);
			}

			ret_block_end_index[current_block_id] = this_chro_tab_end;	// FIRST UNWANTED ID
			ret_block_min_start[current_block_id] = this_block_min_start;
			ret_block_max_end[current_block_id] = this_block_max_end;
			register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
			current_block_id++;
		}

		(*anno_chr_head) [curr_chro_number] = current_block_id; 
		tmp_chro_inf -> chro_block_table_end = current_block_id; 
		free(features_per_block_bins);
	}

	(*block_end_index) = ret_block_end_index;
	(*block_min_start_pos) = ret_block_min_start;
	(*block_max_end_pos) = ret_block_max_end;

	//print_in_box(80, 0,0,"The %u features are sorted.\n", sort_i);
	free(old_info_ptr);
	free(tmp_chro_info_ptrs);
	free(chro_feature_ptr);
}

int strcmp_slash(char * s1, char * s2)
{
	char nch;
	while(0!=(nch = *(s1++))){
		if(nch == '/') break;
		if(nch != (*s2)) return 1;
		s2++;
	}
	return nch != *s2;
}

#define NH_FRACTION_INT 65536

unsigned int calc_fixed_fraction(int nh){
	if(nh==1) return NH_FRACTION_INT;
	else if(nh == 2) return NH_FRACTION_INT>>1;
	else return NH_FRACTION_INT / nh; 
}


int calc_float_fraction(read_count_type_t score, read_count_type_t * integer_count, double * float_count){
	if(score % NH_FRACTION_INT == 0){
		(*integer_count) = score / NH_FRACTION_INT;
		return 0;
	}else{
		(*float_count) = score * 1./NH_FRACTION_INT;
		return 1;
	}
}


void print_read_wrapping(char * rl, int is_second){
	int refill_spaces = 3;

	int read_length = 0, x1 = 0, spaces=0;

	for(x1 = 0; x1 < 3100; x1++){
		if(rl[x1]==0 && rl[x1+1]==0)break;
		if(rl[x1]=='0' || rl[x1]=='\t') spaces++;
		read_length ++;
	}

	char *out_buf1 = malloc(read_length + spaces * refill_spaces + 1), out_buf2[100];
	int ox=0;

	for(x1 = 0; x1 < 3000; x1++){
		if(rl[x1]=='\n' || (rl[x1]==0 && rl[x1+1]==0)){
			out_buf1[ox]=0;
			break;
		} else if((rl[x1]==0 && rl[x1+1]!=0) || rl[x1] == '\t'){
			int x2;
			for(x2 = 0; x2 < refill_spaces ; x2++){
				out_buf1[ox]=' ';
				ox++;
			}
		} else {
			out_buf1[ox]=rl[x1];
			ox++;
		}
	}
	out_buf1[ox] = 0;

	x1=0;

	while(1){
		int x2;
		for(x2 = 0; x2 < 67 ; x2 ++){
			char nch = out_buf1[x1];
			if(nch == 0) break;
			out_buf2[x2] = nch;
			x1++;
		}
		out_buf2[x2] = 0;

		print_in_box(80,0,PRINT_BOX_NOCOLOR_FOR_COLON,"      %s", out_buf2);
		if(out_buf1[x1] == 0)break;
	}

	free(out_buf1);

}

void report_unpair_warning(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, int * this_noproperly_paired_added){
	//printf("WARN:%d [%d]\n", global_context->is_unpaired_warning_shown, thread_context -> thread_id);
	if(!global_context->is_unpaired_warning_shown)
	{
		global_context->is_unpaired_warning_shown=1;
		print_in_box(80,0,0,"   Found reads that are not properly paired.");
		print_in_box(80,0,0,"   (missing mate or the mate is not the next read)");

		if(global_context -> do_not_sort){
			print_in_box(85,0,0,"   %c[31mHowever, the reads will not be re-ordered.", 27);
		}else{
			global_context->redo = 1;
		}
		print_in_box(80,0,0,"   Below are the two reads that are not properly paired:");
		print_read_wrapping(thread_context -> line_buffer1,0);
		print_read_wrapping(thread_context -> line_buffer2,1);

	}
	if(0==(*this_noproperly_paired_added))thread_context -> unpaired_fragment_no++;
	(*this_noproperly_paired_added) = 1;
}


void vote_and_add_count(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context,
			long * hits_indices1, unsigned short * hits_read_start_base1, short * hits_read_len1, int nhits1, unsigned short rl1,
			long * hits_indices2, unsigned short * hits_read_start_base2, short * hits_read_len2, int nhits2, unsigned short rl2,
			int fixed_fractional_count, char * read_name, fc_junction_info_t * supported_junctions1, int njunc1, fc_junction_info_t * supported_junctions2, int njunc2);


void process_line_buffer(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context)
{

	char * read_chr, *read_1_chr = NULL, *tmp_tok_ptr= NULL, *CIGAR_str , *read_name = NULL, *read_name1 = NULL;
	long read_pos, fragment_length = 0, read_1_pos = 0;
	unsigned int search_start = 0, search_end;
	int nhits1 = 0, nhits2 = 0, alignment_masks, search_block_id, search_item_id;
	long * hits_indices1 = thread_context -> hits_indices1, * hits_indices2 = thread_context -> hits_indices2;
	unsigned short * hits_read_start_base1 = thread_context -> hits_read_start_base1 ,  * hits_read_start_base2 = thread_context -> hits_read_start_base2;
	short * hits_read_len1 = thread_context -> hits_read_len1, * hits_read_len2 = thread_context -> hits_read_len2;
	unsigned short cigar_read_len1 = 0, cigar_read_len2 = 0;
	fc_junction_info_t supported_junctions1[FC_CIGAR_PARSER_ITEMS], supported_junctions2[FC_CIGAR_PARSER_ITEMS];
	int njunc1, njunc2;

	int is_second_read;
	int maximum_NH_value = 1;
	int skipped_for_exonic = 0;
	int first_read_quality_score = 0;
	int this_noproperly_paired_added = 0;

	thread_context->all_reads++;
	//if(thread_context->all_reads>1000000) printf("TA=%llu\n%s\n",thread_context->all_reads, thread_context -> line_buffer1);

	if(global_context -> do_junction_counting){
		njunc1 = njunc2 = 0;
	}

	for(is_second_read = 0 ; is_second_read < 2; is_second_read++)
	{
		if(is_second_read && !global_context -> is_paired_end_mode_assign) break;

		char * line = is_second_read? thread_context -> line_buffer2:thread_context -> line_buffer1;
//		if(strstr(line, "CG:Z:"))
	//	SUBREADprintf("LINE_BUF=%s\n",line);

		read_name = strtok_r(line,"\t", &tmp_tok_ptr);	// read name
		if(!read_name)return;

		if(is_second_read)
		{
			if(read_name)
			{
				int x1;
				for(x1=0; read_name[x1]; x1++)
				{
					if(read_name[x1]=='/')
					{
						read_name[x1]=0;
						read_name1[x1]=0;
						break;
					}
				}
				//printf("R1=%s; R2=%s\n",read_name,read_name1 );
				if(strcmp_slash(read_name,read_name1)!=0)
					report_unpair_warning(global_context, thread_context, &this_noproperly_paired_added);
			}
		}
		else
				read_name1 = read_name;

		char * mask_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
		if((!mask_str) || !isdigit(mask_str[0])) return;

		alignment_masks = atoi(mask_str);

		if(is_second_read == 0)
		{
			//skip the read if unmapped (its mate will be skipped as well if paired-end)
			if( ((!global_context -> is_paired_end_mode_assign) &&  (alignment_masks & SAM_FLAG_UNMAPPED) ) ||
			    ((alignment_masks & SAM_FLAG_UNMAPPED)   &&  (alignment_masks & SAM_FLAG_MATE_UNMATCHED) && global_context -> is_paired_end_mode_assign) ||
			    (((alignment_masks & SAM_FLAG_UNMAPPED) || (alignment_masks & SAM_FLAG_MATE_UNMATCHED)) && global_context -> is_paired_end_mode_assign && global_context -> is_both_end_required)
			  ){
				thread_context->read_counters.unassigned_unmapped ++;

				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Unmapped\t*\t*\n", read_name);

				if(global_context -> is_paired_end_mode_assign){
					char * read_name2 = strtok_r(thread_context -> line_buffer2,"\t", &tmp_tok_ptr);
					if(strcmp_slash(read_name,read_name2)!=0)
						report_unpair_warning(global_context, thread_context, &this_noproperly_paired_added);
				}
				return;	// do nothing if a read is unmapped, or the first read in a pair of reads is unmapped.
			}
		}

		if(global_context -> is_paired_end_mode_assign && (!global_context ->is_SEPEmix_warning_shown)){
			if(((!global_context -> is_paired_end_input_file)  && ( alignment_masks & SAM_FLAG_PAIRED_TASK )) || ((global_context -> is_paired_end_input_file)  && 0 == ( alignment_masks & SAM_FLAG_PAIRED_TASK ))){
				print_in_box(85,0,0,"   %c[31mBoth single-end and paired-end reads were found.", 27);
				global_context ->is_SEPEmix_warning_shown = 1;
			}
		}



		read_chr = strtok_r(NULL,"\t", &tmp_tok_ptr);
		if(!read_chr) return;
		char * read_pos_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
		if(!read_pos_str) return;

		read_pos = atoi(read_pos_str);
		if(read_pos < 1 && read_pos_str[0]!='0') return;

		char * mapping_qual_str = strtok_r(NULL,"\t", &tmp_tok_ptr);

		CIGAR_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
		if(!CIGAR_str)
			continue;

		if(global_context -> min_mapping_quality_score>0)
		{
			int mapping_qual =atoi(mapping_qual_str);

			//printf("SECOND=%d; FIRST=%d; THIS=%d; Q=%d\n", is_second_read, first_read_quality_score, mapping_qual, );
			if(( mapping_qual < global_context -> min_mapping_quality_score  && ! global_context -> is_paired_end_mode_assign)||( is_second_read  && max( first_read_quality_score, mapping_qual ) < global_context -> min_mapping_quality_score))
			{
				thread_context->read_counters.unassigned_mappingquality ++;

				if(global_context -> SAM_output_fp)
				{
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_MappingQuality\t*\tMapping_Quality=%d,%d\n", read_name, first_read_quality_score, mapping_qual);
				}
				return;
			}
			if(is_second_read==0 && global_context -> is_paired_end_mode_assign)
			{
				first_read_quality_score = mapping_qual;
			}
		}

		long mate_pos = 0;
		char * mate_chr = NULL;
		char * mate_chr_abs_pos = tmp_tok_ptr;

		if(is_second_read)
		{
			mate_chr = strtok_r(NULL,"\t", &tmp_tok_ptr);// mate_chr
			if(mate_chr[0]=='=') mate_chr = read_chr;
			char * mate_pos_str = strtok_r(NULL,"\t", &tmp_tok_ptr);	// mate_pos
			mate_pos = atol(mate_pos_str);

		}

		if(is_second_read == 0 && global_context -> is_paired_end_mode_assign && 
	   	  (global_context -> is_PE_distance_checked || global_context -> is_chimertc_disallowed)
		  )
		{
			int is_half_mapped = (alignment_masks & SAM_FLAG_UNMAPPED) || (alignment_masks & SAM_FLAG_MATE_UNMATCHED);

			if(!is_half_mapped)
			{
				char * mate_chrx = strtok_r(NULL,"\t", &tmp_tok_ptr); //get chr which the mate read is mapped to
				if(!mate_chrx) return;
				strtok_r(NULL,"\t", &tmp_tok_ptr);
				if(!tmp_tok_ptr) return;
				char * frag_len_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
				if(!tmp_tok_ptr) return;

				fragment_length = abs(atoi(frag_len_str)); //get the fragment length

				int is_first_read_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
				int is_second_read_negative_strand = (alignment_masks & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)?1:0; 

				if(mate_chrx[0]=='=' && is_first_read_negative_strand!=is_second_read_negative_strand)
				{
					if(global_context -> is_PE_distance_checked && ((fragment_length > global_context -> max_paired_end_distance) || (fragment_length < global_context -> min_paired_end_distance)))
					{

						if(global_context -> is_paired_end_mode_assign && 0==is_second_read){
							char * read_name2 = strtok_r(thread_context -> line_buffer2,"\t", &tmp_tok_ptr);
							if(strcmp_slash(read_name,read_name2)!=0)
								report_unpair_warning(global_context, thread_context, &this_noproperly_paired_added);
						}

						thread_context->read_counters.unassigned_fragmentlength ++;

						if(global_context -> SAM_output_fp)
							fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_FragmentLength\t*\tLength=%ld\n", read_name, fragment_length);
						return;
					}
				}
				else
				{
					if(global_context -> is_chimertc_disallowed)
					{

						if(global_context -> is_paired_end_mode_assign && 0==is_second_read){
							char * read_name2 = strtok_r(thread_context -> line_buffer2,"\t", &tmp_tok_ptr);
							if(strcmp_slash(read_name,read_name2)!=0)
								report_unpair_warning(global_context, thread_context, &this_noproperly_paired_added);
						}

						thread_context->read_counters.unassigned_chimericreads ++;

						if(global_context -> SAM_output_fp)
							fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Chimera\t*\t*\n", read_name);
						return;
					}
				}
			}
		}

		if(!tmp_tok_ptr) return;


		// This filter has to be put here because the 0x400 FLAG is not about mapping but about sequencing.
		// A unmapped read with 0x400 FLAG should be able to kill the mapped mate which may have no 0x400 FLAG. 
		if(global_context -> is_duplicate_ignored)
		{
			if(alignment_masks & SAM_FLAG_DUPLICATE)
			{
				if(global_context -> is_paired_end_mode_assign && 0==is_second_read){
					char * read_name2 = strtok_r(thread_context -> line_buffer2,"\t", &tmp_tok_ptr);
					if(strcmp_slash(read_name,read_name2)!=0)
						report_unpair_warning(global_context, thread_context, &this_noproperly_paired_added);
				}

				thread_context->read_counters.unassigned_duplicate ++;
				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Duplicate\t*\t*\n", read_name);

				return;
			}

		}

		if(SAM_FLAG_UNMAPPED & alignment_masks) continue;

		int NH_value = 1;
		char * NH_pos = strstr(tmp_tok_ptr,"\tNH:i:");
		if(NH_pos)
		{
			if(NH_pos[6]>'1' || isdigit(NH_pos[7]))
			{

				if(is_second_read && read_1_chr)
				{
					if((strcmp(read_1_chr, mate_chr)!=0 || mate_pos!=read_1_pos) && read_1_chr[0] != '*'  && mate_chr[0]!='*')
						report_unpair_warning(global_context, thread_context, &this_noproperly_paired_added);
				}
				else
				{
					read_1_chr = read_chr;
					read_1_pos = read_pos;
				}


				if(global_context -> is_multi_mapping_allowed == 0)
				{
					// now it is a NH>1 read!
					// not allow multimapping -> discard!
					thread_context->read_counters.unassigned_multimapping ++;

					if(global_context -> SAM_output_fp)
						fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_MultiMapping\t*\t*\n", read_name);

					if(global_context -> is_paired_end_mode_assign && is_second_read == 0){
						char * read_name2 = strtok_r(thread_context -> line_buffer2,"\t", &tmp_tok_ptr);
						if(strcmp_slash(read_name,read_name2)!=0)
							report_unpair_warning(global_context, thread_context, &this_noproperly_paired_added);
					}
					return;
				}
			}
			int nh_i, NHtmpi=0;
			for(nh_i = 6; nh_i < 15; nh_i++){
				char nch = NH_pos[nh_i];
				if(isdigit(nch)){
					NHtmpi = NHtmpi * 10 + (nch-'0');
				}else{break; }
			}
			NH_value = NHtmpi;
		}


		maximum_NH_value = max(maximum_NH_value, NH_value);

		// if a pair of reads have one secondary, the entire fragment is seen as secondary.
		if((alignment_masks & SAM_FLAG_SECONDARY_MAPPING) && (global_context -> is_multi_mapping_allowed == ALLOW_PRIMARY_MAPPING))
		{
			if(global_context -> is_paired_end_mode_assign && is_second_read == 0){
				char * read_name2 = strtok_r(thread_context -> line_buffer2,"\t", &tmp_tok_ptr);
				if(strcmp_slash(read_name,read_name2)!=0)
					report_unpair_warning(global_context, thread_context, &this_noproperly_paired_added);
			}

			thread_context->read_counters.unassigned_secondary ++;

			if(global_context -> SAM_output_fp)
				fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Secondary\t*\t*\n", read_name);
			return;
		}

		int is_this_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
		int is_fragment_negative_strand = is_this_negative_strand;

		if( global_context -> is_paired_end_mode_assign ){
			int is_second_read_in_pair = alignment_masks & SAM_FLAG_SECOND_READ_IN_PAIR;
			//is_fragment_negative_strand = is_second_read_in_pair?(!is_this_negative_strand):is_this_negative_strand;
			if(is_second_read_in_pair)
				is_fragment_negative_strand = global_context -> is_second_read_straight?is_this_negative_strand:(!is_this_negative_strand);
			else
				is_fragment_negative_strand = global_context -> is_first_read_reversed?(!is_this_negative_strand):is_this_negative_strand;
		}

		int nhits = 0;

		int cigar_section_id, cigar_sections, is_junction_read = 0;
		unsigned int Starting_Chro_Points[FC_CIGAR_PARSER_ITEMS];
		unsigned short Starting_Read_Points[FC_CIGAR_PARSER_ITEMS];
		unsigned short Section_Lengths[FC_CIGAR_PARSER_ITEMS];
		char * ChroNames[FC_CIGAR_PARSER_ITEMS];
		long * hits_indices = (is_second_read?hits_indices2:hits_indices1);
		unsigned short * hits_read_start_base = is_second_read?hits_read_start_base2:hits_read_start_base1;
		short * hits_read_len = is_second_read?hits_read_len2:hits_read_len1;
		unsigned short * cigar_read_len = is_second_read?&cigar_read_len2:&cigar_read_len1;

		if(global_context -> do_junction_counting)
		{
			int this_col = 6;
			while(this_col < 11){
				char nch = *(mate_chr_abs_pos++);
				if(nch == 0 || nch == '\t') this_col++;
			}
			cigar_sections = RSubread_parse_CIGAR_Extra_string(alignment_masks, read_chr, read_pos, CIGAR_str, mate_chr_abs_pos, ChroNames, Starting_Chro_Points, Starting_Read_Points, Section_Lengths, &is_junction_read);
		} else  cigar_sections = RSubread_parse_CIGAR_string(read_chr, read_pos, CIGAR_str, ChroNames, Starting_Chro_Points, Starting_Read_Points, Section_Lengths, &is_junction_read);

		(*cigar_read_len) = Starting_Read_Points[cigar_sections-1] + Section_Lengths[cigar_sections-1];

		if(is_junction_read || !global_context->is_split_alignments_only) 
		{
			if(global_context -> do_junction_counting){
				int * njunc_current = is_second_read?&njunc2:&njunc1;
				fc_junction_info_t * junctions_current = is_second_read?supported_junctions2:supported_junctions1;
				(*njunc_current) = calc_junctions_from_cigar(global_context, alignment_masks , read_chr, read_pos , CIGAR_str, mate_chr_abs_pos, junctions_current);
			}

		//#warning "=================== COMMENT THESE 2 LINES ================================"
		//for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
		//	SUBREADprintf("BCCC: %llu , sec[%d] %s: %u ~ %u ; secs=%d ; flags=%d ; second=%d\n", read_pos, cigar_section_id , ChroNames[cigar_section_id] , Starting_Chro_Points[cigar_section_id], Section_Lengths[cigar_section_id], cigar_sections, alignment_masks, is_second_read);

			if(global_context -> reduce_5_3_ends_to_one)
			{
				if((REDUCE_TO_5_PRIME_END == global_context -> reduce_5_3_ends_to_one) + is_this_negative_strand == 1) // reduce to 5' end (small coordinate if positive strand / large coordinate if negative strand)
				{
					Section_Lengths[0]=1;
				}
				else
				{
					Starting_Chro_Points[0] = Starting_Chro_Points[cigar_sections-1] + Section_Lengths[cigar_sections-1] - 1;
					Section_Lengths[0]=1;
				}

				cigar_sections = 1;
			}

			// Extending the reads to the 3' and 5' ends. (from the read point of view) 
			if(global_context -> five_end_extension)
			{
				if(is_this_negative_strand){
					Section_Lengths [cigar_sections - 1] += global_context -> five_end_extension;
				}else{
					//SUBREADprintf("5-end extension: %d [%d]\n", Starting_Chro_Points[0], Section_Lengths[0]);
					if( read_pos > global_context -> five_end_extension)
					{
						Section_Lengths [0] += global_context -> five_end_extension;
						Starting_Chro_Points [0] -= global_context -> five_end_extension;
					}
					else
					{
						Section_Lengths [0] += read_pos-1;
						Starting_Chro_Points [0] -= read_pos-1;
					}
				}
			}

			if(global_context -> three_end_extension)
			{

				if(is_this_negative_strand){
					if( read_pos > global_context -> three_end_extension)
					{
						Section_Lengths [0] += global_context -> three_end_extension;
						Starting_Chro_Points [0] -= global_context -> three_end_extension;
					}
					else
					{
						Section_Lengths [0] += read_pos - 1;
						Starting_Chro_Points [0] -= read_pos - 1;
					}
				}
				else	Section_Lengths [cigar_sections - 1] += global_context -> three_end_extension;

			}

			for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
			{
				long section_begin_pos = Starting_Chro_Points[cigar_section_id];
				long section_end_pos = Section_Lengths[cigar_section_id] + section_begin_pos - 1;

				
				int start_reverse_table_index = section_begin_pos / REVERSE_TABLE_BUCKET_LENGTH;
				int end_reverse_table_index = (1+section_end_pos) / REVERSE_TABLE_BUCKET_LENGTH;

				fc_chromosome_index_info * this_chro_info = HashTableGet(global_context -> exontable_chro_table, ChroNames[cigar_section_id]);
				if(this_chro_info == NULL)
				{
					if(global_context -> annot_chro_name_alias_table)
					{
						char * anno_chro_name = HashTableGet( global_context -> annot_chro_name_alias_table , ChroNames[cigar_section_id]);
						if(anno_chro_name)
							this_chro_info = HashTableGet(global_context -> exontable_chro_table, anno_chro_name);
					}
					if(this_chro_info == NULL && memcmp(ChroNames[cigar_section_id], "chr", 3)==0)
					{
						this_chro_info = HashTableGet(global_context -> exontable_chro_table, ChroNames[cigar_section_id]+3);
					//	SUBREADprintf("INQ: %p : '%s'\n", this_chro_info , ChroNames[cigar_section_id]+3);
					}
					if(this_chro_info == NULL && strlen(ChroNames[cigar_section_id])<=2)
					{
						strcpy(thread_context -> chro_name_buff, "chr");
						strcpy(thread_context -> chro_name_buff+3, ChroNames[cigar_section_id]);
						this_chro_info = HashTableGet(global_context -> exontable_chro_table, thread_context -> chro_name_buff);
					}
				}

				//SUBREADprintf("INF: %p : %s\n", this_chro_info , ChroNames[cigar_section_id]);

				if(this_chro_info)
				{
					start_reverse_table_index = min(start_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH);
					end_reverse_table_index = min(end_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH+ 1);

					while(start_reverse_table_index<=end_reverse_table_index)
					{
						search_start = this_chro_info -> reverse_table_start_index [start_reverse_table_index];
						if(search_start<0xffffff00)break;
						start_reverse_table_index++;
					}
					if(search_start>0xffffff00) continue;

					//search_start = this_chro_info -> chro_block_table_start;

					search_end = this_chro_info -> chro_block_table_end;//reverse_table_end_index [end_reverse_table_index];
		
					for(search_block_id=search_start;search_block_id<search_end;search_block_id++){
						if (global_context -> exontable_block_min_start[search_block_id] > section_end_pos) break;
						if (global_context -> exontable_block_max_end[search_block_id] < section_begin_pos) continue;

						int search_item_start = 0, search_item_end = global_context -> exontable_block_end_index[search_block_id];
						if(search_block_id>0)search_item_start = global_context -> exontable_block_end_index[search_block_id-1];

						// search_item_id is the inner number of the exons.
						// the exontables in global_index has search_item_id as the index.

						for(search_item_id = search_item_start ; search_item_id < search_item_end; search_item_id++)
						{
							if (global_context -> exontable_stop[search_item_id] >= section_begin_pos)
							{
								if (global_context -> exontable_start[search_item_id] > section_end_pos) break;
								// there is an overlap >=1 between read and feature.
								// the overlap length is min(end_r, end_F) - max(start_r, start_F) + 1
								
								int is_strand_ok =1;

								if(global_context->is_strand_checked){
									if(global_context->is_strand_checked == 1)
										is_strand_ok = (is_fragment_negative_strand == global_context -> exontable_strand[search_item_id]);
									else// if(global_context->is_strand_checked == 2)
										is_strand_ok = (is_fragment_negative_strand != global_context -> exontable_strand[search_item_id]);
									//SUBREADprintf("%d = %d == %d\n", is_strand_ok, is_fragment_negative_strand, global_context -> exontable_strand[search_item_id]);
								}

								if(is_strand_ok){
									if(nhits<=MAX_HIT_NUMBER - 1)
									{
										hits_indices[nhits] = search_item_id;

										if(global_context -> calculate_overlapping_lengths)
										{
											int section_overlapped = min(global_context -> exontable_stop[search_item_id] , section_end_pos) 
													       - max(global_context -> exontable_start[search_item_id] , section_begin_pos) + 1;
											hits_read_len[nhits] = (short)section_overlapped;
											hits_read_start_base[nhits] = (section_begin_pos > global_context -> exontable_start[search_item_id]?(Starting_Read_Points[cigar_section_id]):(Starting_Read_Points[cigar_section_id] + (global_context -> exontable_start[search_item_id] - section_begin_pos)));
										}

										nhits++;
									}
									else break;
								}
							} 
						}
					}
				}
			}
		}else if(global_context->is_split_alignments_only) // must be true.
		{
			skipped_for_exonic ++;
			if((is_second_read && skipped_for_exonic == 2) || (!global_context -> is_paired_end_mode_assign) || (alignment_masks & 0x8))
			{
				if(global_context -> is_paired_end_mode_assign && is_second_read == 0){
					char * read_name2 = strtok_r(thread_context -> line_buffer2,"\t", &tmp_tok_ptr);
					if(strcmp_slash(read_name,read_name2)!=0)
						report_unpair_warning(global_context, thread_context, &this_noproperly_paired_added);
				}

				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Nonjunction\t*\t*\n", read_name);

				thread_context->read_counters.unassigned_nonjunction ++;
				return;
			}
		}

		if(is_second_read) nhits2 = nhits;
		else	nhits1 = nhits;
	}	// loop for is_second_read


	int fixed_fractional_count = global_context -> use_fraction_multi_mapping ?calc_fixed_fraction(maximum_NH_value): NH_FRACTION_INT;

	// we have hits_indices1 and hits_indices2 and nhits1 and nhits2 here
	// we also have fixed_fractional_count which is the value to add

	vote_and_add_count(global_context, thread_context,
			   hits_indices1, hits_read_start_base1, hits_read_len1, nhits1, cigar_read_len1,
			   hits_indices2, hits_read_start_base2, hits_read_len2, nhits2, cigar_read_len2,
			   fixed_fractional_count, read_name, supported_junctions1, njunc1, supported_junctions2, njunc2);
	return;
}

void add_bitmap_overlapping(char * x1_bitmap, short start_base, short len){
	int x1;
	int rl16 = start_base+len-16;
	for(x1 = start_base; x1 < start_base+len; x1++){
		int bit = x1 % 8;
		int byte = x1 / 8;
		if(bit == 0 && x1 < rl16){
			x1_bitmap[byte]=-1;
			x1_bitmap[byte+1]=-1;
			x1+=15;
		}else{
			x1_bitmap[byte] |= (1<<bit);
		}
	}
}

int count_bitmap_overlapping(char * x1_bitmap, unsigned short rl){

	int x1;
	int ret = 0;
	for(x1 = 0; x1 < rl; x1++){
		int byte = x1 / 8;
		int bit = x1 % 8;

		if(bit == 0 && x1_bitmap[byte]==-1){
			x1 += 7;
			ret += 8;
		}else if(x1_bitmap[byte] &  (1<<bit)) ret ++;
	}
	return ret;
}

void add_fragment_supported_junction(	fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, fc_junction_info_t * supported_junctions1,
					int njunc1, fc_junction_info_t * supported_junctions2, int njunc2){
	int x1,x2, in_total_junctions = njunc2 + njunc1;
	for(x1 = 0; x1 < in_total_junctions - 1; x1 ++){
		fc_junction_info_t * j_one = (x1 >= njunc1)?supported_junctions2+(x1-njunc1):supported_junctions1+x1;
		if(j_one->chromosome_name_left[0]==0) continue;

		for(x2 = x1+1; x2 < in_total_junctions ; x2 ++){
			fc_junction_info_t * j_two = (x2 >= njunc1)?supported_junctions2+(x2-njunc1):supported_junctions1+x2;
			if(j_two->chromosome_name_left[0]==0) continue;
			if(
				j_one -> last_exon_base_left == j_two -> last_exon_base_left &&
				j_one -> first_exon_base_right == j_two -> first_exon_base_right &&
				strcmp(j_one -> chromosome_name_left, j_two -> chromosome_name_left) == 0 &&
				strcmp(j_one -> chromosome_name_right, j_two -> chromosome_name_right) == 0
			)
				j_two -> chromosome_name_left[0]=0;
		}
		if(j_one->chromosome_name_left[0]){
			char * this_key = malloc(strlen(j_one->chromosome_name_left) + strlen(j_one->chromosome_name_right)  + 27);
			sprintf(this_key, "%s\t%u\t%s\t%u", j_one->chromosome_name_left, j_one -> last_exon_base_left, j_one->chromosome_name_right, j_one -> first_exon_base_right);
			char * left_key = malloc(strlen(j_one->chromosome_name_left) + 13);
			char * right_key = malloc(strlen(j_one->chromosome_name_right) + 13);
			sprintf(left_key, "%s\t%u", j_one->chromosome_name_left, j_one -> last_exon_base_left);
			sprintf(right_key, "%s\t%u", j_one->chromosome_name_right, j_one -> first_exon_base_right);

			void * count_ptr = HashTableGet(thread_context -> junction_counting_table, this_key);
			unsigned long long count_junc = count_ptr - NULL;
			if(count_ptr != NULL){
				count_junc ++;
				HashTablePut(thread_context -> junction_counting_table, this_key, NULL+count_junc);
			}else	HashTablePut(thread_context -> junction_counting_table, this_key, NULL+1);

			for( x2 = 0 ; x2 < 2 ; x2++ ){
				char * lr_key = x2?right_key:left_key;
				count_ptr = HashTableGet(thread_context -> splicing_point_table, lr_key);
				count_junc = count_ptr - NULL;
				HashTablePut(thread_context -> splicing_point_table, lr_key, NULL + count_junc + 1);
			}
		}
	}
}

void vote_and_add_count(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context,
			long * hits_indices1, unsigned short * hits_read_start_base1, short * hits_read_len1, int nhits1, unsigned short rl1,
			long * hits_indices2, unsigned short * hits_read_start_base2, short * hits_read_len2, int nhits2, unsigned short rl2,
			int fixed_fractional_count, char * read_name, fc_junction_info_t * supported_junctions1, int njunc1, fc_junction_info_t * supported_junctions2, int njunc2)
{
	if(global_context -> calculate_overlapping_lengths == 0 && nhits2+nhits1==1)
	{
		long hit_exon_id = nhits2?hits_indices2[0]:hits_indices1[0];
		thread_context->count_table[hit_exon_id] += fixed_fractional_count;
		thread_context->nreads_mapped_to_exon++;
		if(global_context -> SAM_output_fp)
		{
			int final_gene_number = global_context -> exontable_geneid[hit_exon_id];
			unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
			fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\t*\n", read_name, final_feture_name);
		}
		thread_context->read_counters.assigned_reads ++;

		if(global_context -> do_junction_counting)
			add_fragment_supported_junction(global_context, thread_context, supported_junctions1, njunc1, supported_junctions2, njunc2);
	}
	else if(global_context -> calculate_overlapping_lengths == 0 && nhits2 == 1 && nhits1 == 1 && hits_indices2[0]==hits_indices1[0])
	{
		long hit_exon_id = hits_indices1[0];
		thread_context->count_table[hit_exon_id] += fixed_fractional_count;
		thread_context->nreads_mapped_to_exon++;
		if(global_context -> SAM_output_fp)
		{
			int final_gene_number = global_context -> exontable_geneid[hit_exon_id];
			unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
			fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\t*\n", read_name, final_feture_name);
		}
		thread_context->read_counters.assigned_reads ++;

		if(global_context -> do_junction_counting)
			add_fragment_supported_junction(global_context, thread_context, supported_junctions1, njunc1, supported_junctions2, njunc2);
	}
	else
	{
		int ends;
		// Build a voting table.
		// The voting table should be:
		//      total_length [nhit_final] = total_length_overlapping
		//      final_id [nhit_final] = final_exon_id

		// if is_gene_leven, then decision_table_exon_ids[nhit_final] is the exon id where the count is added.

		// After all, the count is added to all hits where total_length has the maximum value.
		// If there are more than one locations having the same total_length, then the fragment is ambiguous. 
		// Count is added when "-O" is specified.

		// merge feature : if a read overlaps with an EXON twice or more times (by >=2 segments in cigar),
		//                 then the total length of the overlapped bases is calculated.
		// 
		// two ends in a fragment is considered individually; the overlapping bases are not added up.
		//
		char * read_coverage_bits = thread_context -> read_coverage_bits;

		for(ends = 0; ends < global_context -> is_paired_end_mode_assign + 1; ends++)
		{
			long * hits_indices = ends?hits_indices2:hits_indices1;
			unsigned short * hits_read_start_base = ends?hits_read_start_base2:hits_read_start_base1;
			short * hits_read_len = ends?hits_read_len2:hits_read_len1;
			int  nhits = ends?nhits2:nhits1;
			int x1;

			if(read_coverage_bits)
			{
				int covergae_bitmap_size = max(rl1,rl2)/8+2; 
				for(x1 = 0; x1 < nhits; x1++)
					memset(read_coverage_bits +  x1 * (MAX_READ_LENGTH / 8 + 1), 0, covergae_bitmap_size);
			}

			// calculating the summed lengths of overlapping exons
			for(x1=0; x1<nhits;x1++)
			{
				long exon_no = hits_indices[x1];
				if(exon_no>=0x7fffffff) continue;	//already removed

				char * x1_bitmap = read_coverage_bits + (MAX_READ_LENGTH /8 +1) * x1;

				//if(FIXLENstrcmp("V0112_0155:7:1308:19321:196983", read_name)==0)
				//	SUBREADprintf("CREATE bitmap: for x1 = %d, on read %d len = %d \n", x1, hits_read_start_base[x1], hits_read_len[x1]);

				if(read_coverage_bits)
					add_bitmap_overlapping(x1_bitmap, hits_read_start_base[x1], hits_read_len[x1]);

				long merge_key = global_context -> is_gene_level? global_context -> exontable_geneid[exon_no] : exon_no;
				int x2;
				for(x2=x1+1; x2<nhits; x2++)
				{
					long tomerge_exon_no = hits_indices[x2];
					if(tomerge_exon_no >=0x7fffffff) continue;

					long tomerge_key = global_context -> is_gene_level? global_context -> exontable_geneid[tomerge_exon_no] : tomerge_exon_no;
					if(tomerge_key==merge_key)
					{
						if(read_coverage_bits)
							add_bitmap_overlapping(x1_bitmap, hits_read_start_base[x2], hits_read_len[x2]);

						//if(FIXLENstrcmp("V0112_0155:7:1308:19321:196983", read_name)==0)
						//	SUBREADprintf("APPEND bitmap: for x1 = %d, on read %d len = %d \n", x1, hits_read_start_base[x2], hits_read_len[x2]);

						//TODO: change this part to uniquely-overlapping length. Not simply adding.
						//hits_read_start_base[x1]+=hits_read_start_base[x2];
						hits_indices[x2]=0x7fffffff;
						
					}
				}
			}

			// remove the exons in the hits table when it is marked as removed (0x7fffffff)
			int new_hits=0;
			for(x1 = 0; x1< nhits; x1++)
			{
				if(hits_indices[x1]>=0x7fffffff) continue;

				if(new_hits != x1)
					hits_indices[new_hits]=hits_indices[x1];
				
				if(global_context -> calculate_overlapping_lengths)
				{
					char * x1_bitmap = read_coverage_bits + (MAX_READ_LENGTH /8 +1) * x1;
					hits_read_len[new_hits]=count_bitmap_overlapping(x1_bitmap, ends?rl2:rl1);//hits_read_len[x1];
				}

				new_hits++;
			}
			if(ends) nhits2 = new_hits;
			else nhits1 = new_hits;
		}

		unsigned short * decision_total_lengths = thread_context ->decision_table_votes;
		int decision_number = 0, decision_table_no;
		long * decision_table_exon_ids = thread_context -> decision_table_exon_ids;
		long * decision_table_keys = thread_context -> decision_table_ids;
		char * read1_used =  thread_context -> decision_table_read1_used;
		char * read2_used =  thread_context -> decision_table_read2_used;

		for(ends =0 ; ends < global_context -> is_paired_end_mode_assign + 1 ; ends++){
			int nhits = ends?nhits2:nhits1;
			short * hits_read_len = ends?hits_read_len2:hits_read_len1;
			long * hits_indices = ends?hits_indices2:hits_indices1;
			int hit_x1;
			for(hit_x1 = 0; hit_x1 < nhits ; hit_x1++){
				int decision_i;
				long decision_key; // gene_id if is_gene_level, or exon_id.
				long decision_exon_id;
				if (global_context -> is_gene_level ){
					decision_exon_id = hits_indices[hit_x1];
					decision_key = global_context -> exontable_geneid[decision_exon_id];
				}else
					decision_exon_id = decision_key = hits_indices[hit_x1];

				// get decision_table_no from decision_key
				// add into decision_tables if decision_key is unfound.
				decision_table_no = -1;
				for(decision_i = 0; decision_i < decision_number ; decision_i ++){
					if(decision_key == decision_table_keys[decision_i])
					{
						decision_table_no = decision_i;
						break;
					}
				}

				if(decision_table_no<0){
					decision_table_exon_ids[decision_number] = decision_exon_id;
					decision_table_keys[decision_number] = decision_key;
					decision_total_lengths[decision_number] = 0;
					read1_used[decision_number]=0;
					read2_used[decision_number]=0;

					decision_table_no = decision_number;
					decision_number ++;
				}

				if(ends)read2_used[decision_table_no] = 1;
				else read1_used[decision_table_no] = 1;

				assert(read2_used[decision_table_no]<2);
				assert(read1_used[decision_table_no]<2);

				if(global_context -> calculate_overlapping_lengths)
					decision_total_lengths[decision_table_no] += hits_read_len[hit_x1];

			}
		}

		int maximum_decision_score = 0;
		int maximum_total_count = 0;
		int maximum_decision_no = 0;
		for(decision_table_no = 0; decision_table_no < decision_number; decision_table_no++)
		{
			if(global_context -> fragment_minimum_overlapping == 1 || decision_total_lengths[decision_table_no] >= global_context -> fragment_minimum_overlapping)
			{
				int this_decision_score = global_context -> use_overlapping_break_tie? decision_total_lengths[decision_table_no] :(read1_used[decision_table_no] + read2_used[decision_table_no]);
				if(!global_context -> use_overlapping_break_tie)assert(this_decision_score<3);
				//SUBREADprintf("%s LEN[%d] = %d , score=%d\n", read_name, decision_table_no ,  decision_total_lengths[decision_table_no], this_decision_score );
				if(maximum_decision_score < this_decision_score){
					maximum_total_count = 1;
					maximum_decision_no = decision_table_no;
					maximum_decision_score = this_decision_score;
				}else if(maximum_decision_score == this_decision_score)
					maximum_total_count ++;
			}
		}

		if(maximum_total_count == 0){
			if(global_context -> SAM_output_fp)
				fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_NoFeatures\t*\t*\n", read_name);

			thread_context->read_counters.unassigned_nofeatures ++;
		}else{

			// final adding votes.
			if(1 == maximum_total_count && !global_context -> is_multi_overlap_allowed) {
				// simple add to the exon ( EXON_ID = decision_table_exon_ids[maximum_decision_no])
				long max_exon_id = decision_table_exon_ids[maximum_decision_no];
				thread_context->count_table[max_exon_id] += fixed_fractional_count;
				thread_context->nreads_mapped_to_exon++;
				if(global_context -> SAM_output_fp)
				{
					int final_gene_number = global_context -> exontable_geneid[max_exon_id];
					unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
					if(decision_number>1)
						fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\t%s/Targets=%d/%d\n", read_name, final_feture_name, global_context -> use_overlapping_break_tie? "MaximumOverlapping":"Votes", maximum_decision_score, decision_number);
					else
						fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\t*\n", read_name, final_feture_name);
				}
				thread_context->read_counters.assigned_reads ++;
				if(global_context -> do_junction_counting)
					add_fragment_supported_junction(global_context, thread_context, supported_junctions1, njunc1, supported_junctions2, njunc2);

			}else if(global_context -> is_multi_overlap_allowed) {
				char final_feture_names[1000];
				int assigned_no = 0, xk1;
				final_feture_names[0]=0;
				for(xk1 = 0; xk1 < decision_number; xk1++)
				{
					long tmp_voter_id = decision_table_exon_ids[xk1];
					thread_context->count_table[tmp_voter_id] += fixed_fractional_count;

					if(global_context -> SAM_output_fp)
					{
						if(strlen(final_feture_names)<700)
						{
							int final_gene_number = global_context -> exontable_geneid[tmp_voter_id];
							unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
							strncat(final_feture_names, (char *)final_feture_name, 999);
							strncat(final_feture_names, ",", 999);
							assigned_no++;
						}
					}
				}
				final_feture_names[999]=0;
				thread_context->nreads_mapped_to_exon++;
				if(global_context -> SAM_output_fp)
				{
					int ffnn = strlen(final_feture_names);
					if(ffnn>0) final_feture_names[ffnn-1]=0;
					// overlapped but still assigned 
					fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\tTotal=%d\n", read_name, final_feture_names, assigned_no);
				}
				thread_context->read_counters.assigned_reads ++;
				if(global_context -> do_junction_counting)
					add_fragment_supported_junction(global_context, thread_context, supported_junctions1, njunc1, supported_junctions2, njunc2);

			} else {
				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Ambiguit\t*\tNumber_Of_Overlapped_Genes=%d\n", read_name, maximum_total_count);

				thread_context->read_counters.unassigned_ambiguous ++;
			}
		}
	}
}


void * feature_count_worker(void * vargs)
{
	void ** args = (void **) vargs;

	fc_thread_global_context_t * global_context = args[0];
	fc_thread_thread_context_t * thread_context = args[1];

	free(vargs);


	//printf("QQQ0:T%d\n", thread_context->thread_id);
	//Rprintf("QQQ1:T%d\n", thread_context->thread_id);
	//printf("QQQ2:T%d\n", thread_context->thread_id);

	if(global_context -> is_SAM_file)
	{
		//thread_context -> current_read_length1 = global_context -> read_length;
		//thread_context -> current_read_length2 = global_context -> read_length;
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

					//if(buffer_read_ptr>= global_context->input_buffer_max_size)
					//	if(buffer_read_ptr>6*1024*1024) printf("REALLY BIG PTR:%u = %u + %u - %u\n", buffer_read_ptr, thread_context->input_buffer_write_ptr , global_context->input_buffer_max_size, thread_context->input_buffer_remainder);

					for(is_second_read = 0; is_second_read < (global_context->is_paired_end_mode_assign ? 2:1); is_second_read++)
					{
						char * curr_line_buff = is_second_read?thread_context -> line_buffer2:thread_context -> line_buffer1;
						//printf("R=%u; WPTR=%u ;RPTR=%u\n", thread_context->input_buffer_remainder, thread_context->input_buffer_write_ptr, buffer_read_ptr);
						//if(buffer_read_ptr % 7 == 0)
						//	fflush(stdout);
						
						for(buffer_read_bytes=0; ; buffer_read_bytes++)
						{
							//printf("%p + %d\n", thread_context->input_buffer, buffer_read_ptr);
							//if(buffer_read_ptr>6*1024*1024) printf("VERY BIG PTR:%u > %u\n", buffer_read_ptr , global_context->input_buffer_max_size);
							char nch =  thread_context->input_buffer[buffer_read_ptr ++];
							curr_line_buff[buffer_read_bytes] = nch;
							if(buffer_read_ptr >= global_context->input_buffer_max_size)
								buffer_read_ptr = 0; 
							if(nch=='\n' || buffer_read_bytes>2998){
								curr_line_buff[buffer_read_bytes+1]=0;
								curr_line_buff[buffer_read_bytes+2]=0;
								break;
							}
						}

						//printf("%s\n", curr_line_buff);

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



			process_line_buffer(global_context, thread_context);

		}
	}
	else
	{	// if is BAM: decompress the chunk and process reads.
		char * PDATA = malloc(2*70000);
		SamBam_Alignment * aln = &thread_context->aln_buffer;

		//thread_context -> current_read_length1 = global_context -> read_length;
		//thread_context -> current_read_length2 = global_context -> read_length;
		while(1)
		{
			int PDATA_len = 0;
			while(1)
			{
				int is_retrieved = 0;
				PDATA_len = 0;
				//retrieve the next chunk.

				pthread_spin_lock(&thread_context->input_buffer_lock);
				if(thread_context->input_buffer_remainder)
				{
					assert(thread_context->input_buffer_remainder>4);
					unsigned int tail_bytes = global_context->input_buffer_max_size - thread_context -> chunk_read_ptr ;
					if(tail_bytes<4)
					{
						thread_context -> chunk_read_ptr = 0;
						thread_context -> input_buffer_remainder -= tail_bytes;
						memcpy(&PDATA_len, thread_context->input_buffer + thread_context -> chunk_read_ptr , 4);
					}
					else
					{
						memcpy(&PDATA_len, thread_context->input_buffer + thread_context -> chunk_read_ptr , 4);
						if(PDATA_len==0)
						{
							thread_context -> chunk_read_ptr = 0;
							thread_context -> input_buffer_remainder -= tail_bytes;
							memcpy(&PDATA_len, thread_context->input_buffer , 4);
						}
					}
					thread_context -> chunk_read_ptr+=4;
					thread_context -> input_buffer_remainder -= 4;

					//fprintf(stderr,"chunk_read_ptr=%d , input_buffer_remainder = %d\n", thread_context -> chunk_read_ptr , thread_context -> input_buffer_remainder);
					if(PDATA_len<0 || PDATA_len > 140000)
					{
						SUBREADprintf("THREAD ABNORMALLY QUIT\n");
						return NULL;
					}

					memcpy(PDATA, thread_context -> input_buffer + thread_context -> chunk_read_ptr , PDATA_len);
					thread_context -> chunk_read_ptr += PDATA_len;
					thread_context -> input_buffer_remainder -= PDATA_len;

					if( PDATA_len > 0 )
						is_retrieved = 1;
				}


				pthread_spin_unlock(&thread_context->input_buffer_lock);
				if(global_context->is_all_finished && !is_retrieved){
					free(PDATA);
					return NULL;
				}

				if(is_retrieved) break;
				else
					usleep(tick_time);

			}
			
			// convert binary reads into sam lines and process;
			int processed_reads = 0, PDATA_ptr = 0;
			while(PDATA_ptr < PDATA_len)
			{
				int is_second_read;
				for(is_second_read = 0; is_second_read <= global_context -> is_paired_end_mode_assign; is_second_read++)
				{
					int binary_read_len, local_PDATA_ptr = PDATA_ptr;
					char * curr_line_buff = is_second_read?thread_context -> line_buffer2:thread_context -> line_buffer1;

					memcpy(&binary_read_len, PDATA + PDATA_ptr, 4);
					int ret = PBam_chunk_gets(PDATA, &local_PDATA_ptr, PDATA_len, global_context -> sambam_chro_table, curr_line_buff, 2999, aln,0);
					//printf("LL=%s\n", curr_line_buff);
					if(ret<0)
						SUBREADprintf("READ DECODING ERROR!\n");

					PDATA_ptr += 4+binary_read_len;
					processed_reads++;
				}

				process_line_buffer(global_context, thread_context);
				//printf("LE\n\n");
			}
		}
	}
}

void fc_thread_merge_results(fc_thread_global_context_t * global_context, read_count_type_t * nreads , unsigned long long int *nreads_mapped_to_exon, fc_read_counters * my_read_counter, HashTable * junction_global_table, HashTable * splicing_global_table)
{
	int xk1, xk2;

	long long int total_input_reads = 0 ;
	read_count_type_t unpaired_fragment_no = 0;

	(*nreads_mapped_to_exon)=0;

	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		for(xk2=0; xk2<global_context -> exontable_exons; xk2++)
		{
			nreads[xk2]+=global_context -> thread_contexts[xk1].count_table[xk2];
		}
		total_input_reads += global_context -> thread_contexts[xk1].all_reads;
		(*nreads_mapped_to_exon) += global_context -> thread_contexts[xk1].nreads_mapped_to_exon;
		unpaired_fragment_no += global_context -> thread_contexts[xk1].unpaired_fragment_no;

		global_context -> read_counters.unassigned_ambiguous += global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous;
		global_context -> read_counters.unassigned_nofeatures += global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures;
		global_context -> read_counters.unassigned_unmapped += global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped;
		global_context -> read_counters.unassigned_mappingquality += global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality;
		global_context -> read_counters.unassigned_fragmentlength += global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength;
		global_context -> read_counters.unassigned_chimericreads += global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads;
		global_context -> read_counters.unassigned_multimapping += global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping;
		global_context -> read_counters.unassigned_secondary += global_context -> thread_contexts[xk1].read_counters.unassigned_secondary;
		global_context -> read_counters.unassigned_nonjunction += global_context -> thread_contexts[xk1].read_counters.unassigned_nonjunction;
		global_context -> read_counters.unassigned_duplicate += global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate;
		global_context -> read_counters.assigned_reads += global_context -> thread_contexts[xk1].read_counters.assigned_reads;

		my_read_counter->unassigned_ambiguous += global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous;
		my_read_counter->unassigned_nofeatures += global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures;
		my_read_counter->unassigned_unmapped += global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped;
		my_read_counter->unassigned_mappingquality += global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality;
		my_read_counter->unassigned_fragmentlength += global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength;
		my_read_counter->unassigned_chimericreads += global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads;
		my_read_counter->unassigned_multimapping += global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping;
		my_read_counter->unassigned_secondary += global_context -> thread_contexts[xk1].read_counters.unassigned_secondary;
		my_read_counter->unassigned_nonjunction += global_context -> thread_contexts[xk1].read_counters.unassigned_nonjunction;
		my_read_counter->unassigned_duplicate += global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate;
		my_read_counter->assigned_reads += global_context -> thread_contexts[xk1].read_counters.assigned_reads;

		if(global_context -> do_junction_counting){
			int bucket_i;
			for(bucket_i = 0 ; bucket_i < global_context -> thread_contexts[xk1].junction_counting_table -> numOfBuckets; bucket_i++){
				KeyValuePair * cursor;
				cursor = global_context -> thread_contexts[xk1].junction_counting_table -> bucketArray[bucket_i];
				while(cursor){
					char * junckey = (char *) cursor -> key;

					void * globval = HashTableGet(junction_global_table, junckey);
					char * new_key = malloc(strlen(junckey)+1);
					strcpy(new_key, junckey);
					globval += (cursor -> value - NULL);
					HashTablePut(junction_global_table, new_key, globval);
						// new_key will be freed when it is replaced next time or when the global table is destroyed.

					cursor = cursor->next;
				}
			}

			for(bucket_i = 0 ; bucket_i < global_context -> thread_contexts[xk1].splicing_point_table -> numOfBuckets; bucket_i++){
				KeyValuePair * cursor;
				cursor = global_context -> thread_contexts[xk1].splicing_point_table -> bucketArray[bucket_i];
				while(cursor){
					char * junckey = (char *) cursor -> key;
					void * globval = HashTableGet(splicing_global_table, junckey);
					char * new_key = malloc(strlen(junckey)+1);
					strcpy(new_key, junckey);

					//if(xk1>0)
					//SUBREADprintf("MERGE THREAD-%d : %s    VAL=%u, ADD=%u\n", xk1, junckey, globval - NULL, cursor -> value - NULL);

					globval += (cursor -> value - NULL);
					HashTablePut(splicing_global_table, new_key, globval);
					cursor = cursor->next;
				}
			}
		}
	}

	char pct_str[10];
	if(total_input_reads>0)
		sprintf(pct_str,"(%.1f%%%%)", (*nreads_mapped_to_exon)*100./total_input_reads);
	else	pct_str[0]=0;

	if(unpaired_fragment_no){
		print_in_box(80,0,0,"   Not properly paired fragments : %llu", unpaired_fragment_no);
	}
	print_in_box(80,0,0,"   Total %s : %llu", global_context -> is_paired_end_mode_assign?"fragments":"reads", total_input_reads); 
	print_in_box(pct_str[0]?81:80,0,0,"   Successfully assigned %s : %llu %s", global_context -> is_paired_end_mode_assign?"fragments":"reads", *nreads_mapped_to_exon,pct_str); 
	print_in_box(80,0,0,"   Running time : %.2f minutes", (miltime() - global_context -> start_time)/60);
	print_in_box(80,0,0,"");
}

HashTable * load_alias_table(char * fname)
{
	FILE * fp = f_subr_open(fname, "r");
	if(!fp)
	{
		print_in_box(80,0,0,"WARNING unable to open alias file '%s'", fname);
		return NULL;
	}

	char * fl = malloc(2000);

	HashTable * ret = HashTableCreate(1013);
	HashTableSetDeallocationFunctions(ret, free, free);
	HashTableSetKeyComparisonFunction(ret, fc_strcmp);
	HashTableSetHashFunction(ret, fc_chro_hash);
	
	while (1)
	{
		char *ret_fl = fgets(fl, 1999, fp);
		if(!ret_fl) break;
		if(fl[0]=='#') continue;
		char * sam_chr = NULL;
		char * anno_chr = strtok_r(fl, ",", &sam_chr);
		if((!sam_chr)||(!anno_chr)) continue;

		sam_chr[strlen(sam_chr)-1]=0;
		char * anno_chr_buf = malloc(strlen(anno_chr)+1);
		strcpy(anno_chr_buf, anno_chr);
		char * sam_chr_buf = malloc(strlen(sam_chr)+1);
		strcpy(sam_chr_buf, sam_chr);
		
		//printf("ALIAS: %s -> %s\n", sam_chr, anno_chr);
		HashTablePut(ret, sam_chr_buf, anno_chr_buf);
	}


	free(fl);
	return ret;
}

void fc_thread_init_global_context(fc_thread_global_context_t * global_context, unsigned int buffer_size, unsigned short threads, int line_length , int is_PE_data, int min_pe_dist, int max_pe_dist, int is_gene_level, int is_overlap_allowed, int is_strand_checked, char * output_fname, int is_sam_out, int is_both_end_required, int is_chimertc_disallowed, int is_PE_distance_checked, char *feature_name_column, char * gene_id_column, int min_map_qual_score, int is_multi_mapping_allowed, int is_SAM, char * alias_file_name, char * cmd_rebuilt, int is_input_file_resort_needed, int feature_block_size, int isCVersion, int fiveEndExtension,  int threeEndExtension, int minFragmentOverlap, int is_split_alignments_only, int reduce_5_3_ends_to_one, char * debug_command, int is_duplicate_ignored, int is_not_sort, int use_fraction_multimapping, int useOverlappingBreakTie, char * pair_orientations, int do_junction_cnt)
{

	global_context -> input_buffer_max_size = buffer_size;
	global_context -> all_reads = 0;
	global_context -> redo = 0;
	global_context -> SAM_output_fp = NULL;


	global_context -> isCVersion = isCVersion;
	global_context -> is_read_details_out = is_sam_out;
	global_context -> is_multi_overlap_allowed = is_overlap_allowed;
	global_context -> is_paired_end_mode_assign = is_PE_data;
	global_context -> is_gene_level = is_gene_level;
	global_context -> is_strand_checked = is_strand_checked;
	global_context -> is_both_end_required = is_both_end_required;
	global_context -> is_chimertc_disallowed = is_chimertc_disallowed;
	global_context -> is_PE_distance_checked = is_PE_distance_checked;
	global_context -> is_multi_mapping_allowed = is_multi_mapping_allowed;
	global_context -> is_split_alignments_only = is_split_alignments_only;
	global_context -> is_duplicate_ignored = is_duplicate_ignored;
	global_context -> is_first_read_reversed = (pair_orientations[0]=='r');
	global_context -> is_second_read_straight = (pair_orientations[1]=='f');

	global_context -> reduce_5_3_ends_to_one = reduce_5_3_ends_to_one;
	global_context -> do_not_sort = is_not_sort;
	global_context -> is_SAM_file = is_SAM;
	global_context -> use_fraction_multi_mapping = use_fraction_multimapping;
	global_context -> do_junction_counting = do_junction_cnt;

	global_context -> thread_number = threads;
	global_context -> min_mapping_quality_score = min_map_qual_score;
	global_context -> unistr_buffer_size = 1024*1024*2;
	global_context -> unistr_buffer_used = 0;
	global_context -> unistr_buffer_space = malloc(global_context -> unistr_buffer_size);
	global_context -> annot_chro_name_alias_table = NULL;
	global_context -> cmd_rebuilt = cmd_rebuilt;
	global_context -> is_input_file_resort_needed = is_input_file_resort_needed;
	global_context -> feature_block_size = feature_block_size;
	global_context -> five_end_extension = fiveEndExtension;
	global_context -> three_end_extension = threeEndExtension;
	global_context -> fragment_minimum_overlapping = minFragmentOverlap;
	global_context -> use_overlapping_break_tie = useOverlappingBreakTie;
	global_context -> calculate_overlapping_lengths = (global_context -> fragment_minimum_overlapping > 1) || global_context -> use_overlapping_break_tie;
	global_context -> debug_command = debug_command;

	global_context -> read_counters.unassigned_ambiguous=0;
	global_context -> read_counters.unassigned_nofeatures=0;
	global_context -> read_counters.unassigned_unmapped=0;
	global_context -> read_counters.unassigned_mappingquality=0;
	global_context -> read_counters.unassigned_fragmentlength=0;
	global_context -> read_counters.unassigned_chimericreads=0;
	global_context -> read_counters.unassigned_multimapping=0;
	global_context -> read_counters.unassigned_secondary=0;
	global_context -> read_counters.unassigned_nonjunction=0;
	global_context -> read_counters.unassigned_duplicate=0;
	global_context -> read_counters.assigned_reads=0;
	
	if(alias_file_name && alias_file_name[0])
	{
		strcpy(global_context -> alias_file_name,alias_file_name);
		global_context -> annot_chro_name_alias_table = load_alias_table(alias_file_name);
	}
	else	global_context -> alias_file_name[0]=0;

	strcpy(global_context -> feature_name_column,feature_name_column);
	strcpy(global_context -> gene_id_column,gene_id_column);
	strcpy(global_context -> output_file_name, output_fname);

	global_context -> min_paired_end_distance = min_pe_dist;
	global_context -> max_paired_end_distance = max_pe_dist;
	global_context -> thread_number = threads;
	global_context -> line_length = line_length;


}
int fc_thread_start_threads(fc_thread_global_context_t * global_context, int et_exons, int * et_geneid, char ** et_chr, long * et_start, long * et_stop, unsigned char * et_strand, char * et_anno_chr_2ch, char ** et_anno_chrs, long * et_anno_chr_heads, long * et_bk_end_index, long * et_bk_min_start, long * et_bk_max_end, int read_length)
{
	int xk1;

	global_context -> read_length = read_length;
	global_context -> is_unpaired_warning_shown = 0;
	global_context -> is_stake_warning_shown = 0;


	if(global_context -> is_read_details_out)
	{
		char tmp_fname[350], *modified_fname;
		int i=0;
		sprintf(tmp_fname, "%s.featureCounts", global_context -> raw_input_file_name);
		modified_fname = tmp_fname;
		while(modified_fname[0]=='/' || modified_fname[0]=='.' || modified_fname[0]=='\\'){
			modified_fname ++;
		}
		while(modified_fname[i]){
			if(modified_fname[i]=='\\' || modified_fname[i]=='/'||modified_fname[i]==' ')modified_fname[i]='.';
			i++;
		}
		global_context -> SAM_output_fp = f_subr_open(modified_fname, "w");
		if(!global_context -> SAM_output_fp)
		{
			SUBREADprintf("Unable to create file '%s'; the read assignment details are not written.\n", tmp_fname);
		}
	}
	else
		global_context -> SAM_output_fp = NULL;

	global_context -> redo = 0;
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
	//	printf("CHRR_MALLOC\n");
		pthread_spin_init(&global_context->thread_contexts[xk1].input_buffer_lock, PTHREAD_PROCESS_PRIVATE);
		global_context -> thread_contexts[xk1].input_buffer_remainder = 0;
		global_context -> thread_contexts[xk1].input_buffer_write_ptr = 0;
		global_context -> thread_contexts[xk1].input_buffer = malloc(global_context -> input_buffer_max_size);
		global_context -> thread_contexts[xk1].thread_id = xk1;
		global_context -> thread_contexts[xk1].chunk_read_ptr = 0;
		global_context -> thread_contexts[xk1].count_table = calloc(sizeof(read_count_type_t), et_exons);
		global_context -> thread_contexts[xk1].nreads_mapped_to_exon = 0;
		global_context -> thread_contexts[xk1].all_reads = 0;
		global_context -> thread_contexts[xk1].line_buffer1 = malloc(global_context -> line_length + 2);
		global_context -> thread_contexts[xk1].line_buffer2 = malloc(global_context -> line_length + 2);
		global_context -> thread_contexts[xk1].chro_name_buff = malloc(CHROMOSOME_NAME_LENGTH);
		global_context -> thread_contexts[xk1].strm_buffer = malloc(sizeof(z_stream));
		if(global_context -> calculate_overlapping_lengths)
			global_context -> thread_contexts[xk1].read_coverage_bits = malloc((MAX_READ_LENGTH / 8 + 1)*MAX_HIT_NUMBER);
		else global_context -> thread_contexts[xk1].read_coverage_bits = NULL;

		global_context -> thread_contexts[xk1].unpaired_fragment_no = 0;
		global_context -> thread_contexts[xk1].read_counters.assigned_reads = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_secondary = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_nonjunction = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate = 0;

		if(global_context -> do_junction_counting)
		{
			global_context -> thread_contexts[xk1].junction_counting_table = HashTableCreate(131317);
			HashTableSetHashFunction(global_context -> thread_contexts[xk1].junction_counting_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(global_context -> thread_contexts[xk1].junction_counting_table, free, NULL);
			HashTableSetKeyComparisonFunction(global_context -> thread_contexts[xk1].junction_counting_table, fc_strcmp_chro);
			
			global_context -> thread_contexts[xk1].splicing_point_table = HashTableCreate(131317);
			HashTableSetHashFunction(global_context -> thread_contexts[xk1].splicing_point_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(global_context -> thread_contexts[xk1].splicing_point_table, free, NULL);
			HashTableSetKeyComparisonFunction(global_context -> thread_contexts[xk1].splicing_point_table, fc_strcmp_chro);
		}

		if(!global_context ->  thread_contexts[xk1].count_table) return 1;
		void ** thread_args = malloc(sizeof(void *)*2);
		thread_args[0] = global_context;
		thread_args[1] = & global_context -> thread_contexts[xk1];

		if(global_context -> thread_number>1 || ! global_context -> is_SAM_file)
			pthread_create(&global_context -> thread_contexts[xk1].thread_object, NULL, feature_count_worker, thread_args);
	}

	return 0;
}

void fc_thread_destroy_thread_context(fc_thread_global_context_t * global_context)
{
	int xk1;
	if(global_context -> is_read_details_out)
	{
		fclose(global_context -> SAM_output_fp);
		global_context -> SAM_output_fp = NULL;
	}

	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		//printf("CHRR_FREE\n");
		if(global_context -> thread_contexts[xk1].read_coverage_bits)
			free(global_context -> thread_contexts[xk1].read_coverage_bits);	
		free(global_context -> thread_contexts[xk1].count_table);	
		free(global_context -> thread_contexts[xk1].line_buffer1);	
		free(global_context -> thread_contexts[xk1].line_buffer2);	
		free(global_context -> thread_contexts[xk1].input_buffer);
		free(global_context -> thread_contexts[xk1].chro_name_buff);
		free(global_context -> thread_contexts[xk1].strm_buffer);
		pthread_spin_destroy(&global_context -> thread_contexts[xk1].input_buffer_lock);
		if(global_context -> do_junction_counting){
			HashTableDestroy(global_context -> thread_contexts[xk1].junction_counting_table);
			HashTableDestroy(global_context -> thread_contexts[xk1].splicing_point_table);
		}
	}
	free(global_context -> thread_contexts);
}
void fc_thread_wait_threads(fc_thread_global_context_t * global_context)
{
	int xk1;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
		pthread_join(global_context -> thread_contexts[xk1].thread_object, NULL);
}

int resort_input_file(fc_thread_global_context_t * global_context)
{
	char * temp_file_name = malloc(300), * fline = malloc(3000);

	if(!global_context->redo)
		print_in_box(80,0,0,"   Resort the input file ...");
	sprintf(temp_file_name, "./temp-core-%06u-%08X.sam", getpid(), rand());

	unsigned long long unpaired_in_sort = 0;
	if( global_context-> is_SAM_file ){
		SamBam_FILE * sambam_reader ;
		sambam_reader = SamBam_fopen(global_context-> input_file_name, global_context-> is_SAM_file?SAMBAM_FILE_SAM:SAMBAM_FILE_BAM);

		if(!sambam_reader){
			SUBREADprintf("Unable to open %s.\n", global_context-> input_file_name);
			return -1;
		}
		SAM_sort_writer writer;
		int ret = sort_SAM_create(&writer, temp_file_name, ".");
		if(ret)
		{
			SUBREADprintf("Unable to sort input file because temporary file '%s' cannot be created.\n", temp_file_name);
			return -1;
		}
		int is_read_len_warned = 0;

		while(1)
		{
			//#warning "Find out why I ndded the sequences, then replace '1 - 1' with 0"
			char * is_ret = SamBam_fgets(sambam_reader, fline, 2999, 1 - 1);
			if(!is_ret) break;
			int ret = sort_SAM_add_line(&writer, fline, strlen(fline));
			if(ret<0) 
			{
				if(!is_read_len_warned)
					print_in_box(80,0,0,"WARNING: reads with very long names were found.");
				is_read_len_warned = 1;
			//	break;
			}
		//printf("N1=%llu\n",  writer.unpaired_reads);
		}

		sort_SAM_finalise(&writer);
		global_context->is_SAM_file = 1;
		SamBam_fclose(sambam_reader);
	}else{
		char temp_temp_name[300];
		sprintf(temp_temp_name, "./fspr-FC-%06u-%08X",  getpid(), rand());
		SAM_pairer_context_t pairer;
		SAM_pairer_writer_main_t writer_main;
		SAM_pairer_writer_create(&writer_main, global_context -> thread_number, 1, 1, Z_NO_COMPRESSION, temp_file_name);
		SAM_pairer_create(&pairer, global_context -> thread_number, 64, 1, 1, 0, global_context-> input_file_name, SAM_pairer_writer_reset, SAM_pairer_multi_thread_header, SAM_pairer_multi_thread_output, temp_temp_name, &writer_main);
		SAM_pairer_run(&pairer);
		SAM_pairer_destroy(&pairer);
		SAM_pairer_writer_destroy(&writer_main);
		unpaired_in_sort = pairer.total_orphan_reads;
		global_context->is_SAM_file = 0;
	}
	print_in_box(80,0,0,"   %llu read%s ha%s missing mates.", unpaired_in_sort, unpaired_in_sort>1?"s":"", unpaired_in_sort>1?"ve":"s");
	print_in_box(80,0,0,"   Input was converted to a format accepted by featureCounts.");

	strcpy(global_context-> input_file_name, temp_file_name);
	free(temp_file_name);
	free(fline);
	return 0;
}


void BUFstrcat(char * targ, char * src, char ** buf){
	int srclen = strlen(src);
	if( (*buf) == NULL){
		(*buf) = targ;
	}
	memcpy((*buf), src, srclen);
	(*buf) += srclen;
	(**buf) = 0;
}

void fc_write_final_gene_results(fc_thread_global_context_t * global_context, int * et_geneid, char ** et_chr, long * et_start, long * et_stop, unsigned char * et_strand, const char * out_file, int features, read_count_type_t ** column_numbers, char * file_list, int n_input_files, fc_feature_info_t * loaded_features, int header_out)
{
	int xk1;
	int genes = global_context -> gene_name_table -> numOfElements;
	read_count_type_t *gene_columns;

	FILE * fp_out = f_subr_open(out_file,"w");
	if(!fp_out){
		SUBREADprintf("Failed to create file %s\n", out_file);
		return;
	}

	if(header_out)
	{
		fprintf(fp_out, "# Program:featureCounts v%s", SUBREAD_VERSION);
		if(global_context->cmd_rebuilt)
			fprintf(fp_out, "; Command:%s", global_context->cmd_rebuilt);
		fprintf(fp_out, "\n");
	}

	char * tmp_ptr = NULL, * next_fn;
	int non_empty_files = 0, i_files=0;
	fprintf(fp_out,"Geneid\tChr\tStart\tEnd\tStrand\tLength");
	next_fn = strtok_r(file_list, ";", &tmp_ptr);
	while(1){
		if(!next_fn||strlen(next_fn)<1) break;
		if(column_numbers[i_files])
		{
			fprintf(fp_out,"\t%s", next_fn);
			non_empty_files ++;
		}
		next_fn = strtok_r(NULL, ";", &tmp_ptr);
		i_files++;
	}
	fprintf(fp_out,"\n");

	gene_columns = calloc(sizeof(read_count_type_t) , genes * non_empty_files);
	unsigned int * gene_exons_number = calloc(sizeof(unsigned int) , genes);
	unsigned int * gene_exons_pointer = calloc(sizeof(unsigned int) , genes);
	unsigned int * gene_exons_start = malloc(sizeof(unsigned int) * features);
	unsigned int * gene_exons_end = malloc(sizeof(unsigned int) * features);
	char ** gene_exons_chr = malloc(sizeof(char *) * features);
	char * gene_exons_strand = malloc(features);

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1];
		gene_exons_number[gene_id]++;
	}

	unsigned int accumulative_no = 0;
	unsigned longest_gene_exons = 0;
	for(xk1 = 0 ; xk1 < genes; xk1++)
	{
		unsigned int tmpv = gene_exons_number[xk1];
		longest_gene_exons = max(longest_gene_exons, tmpv);
		accumulative_no += gene_exons_number[xk1];
		gene_exons_number[xk1] = accumulative_no - tmpv;
	}

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1];
		int gene_write_ptr = gene_exons_number[gene_id] + gene_exons_pointer[gene_id];

		gene_exons_chr[gene_write_ptr] = et_chr[xk1];
		gene_exons_start[gene_write_ptr] = et_start[xk1]; 
		gene_exons_end[gene_write_ptr] = et_stop[xk1]; 
		gene_exons_strand[gene_write_ptr] = et_strand[xk1]; 

		gene_exons_pointer[gene_id]++;
	}

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1], k_noempty = 0;
		for(i_files=0;i_files < n_input_files; i_files++)
		{
			if(column_numbers[i_files]==NULL) continue;
			gene_columns[gene_id * non_empty_files + k_noempty ] += column_numbers[i_files][xk1];
			k_noempty++;
		}
	}


	char *is_occupied = malloc(longest_gene_exons);
	unsigned int * input_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);
	unsigned int * output_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);

	char * out_chr_list = malloc(longest_gene_exons * (1+global_context -> longest_chro_name) + 1), * tmp_chr_list = NULL;
	char * out_start_list = malloc(11 * longest_gene_exons + 1), * tmp_start_list = NULL;
	char * out_end_list = malloc(11 * longest_gene_exons + 1), * tmp_end_list = NULL;
	char * out_strand_list = malloc(2 * longest_gene_exons + 1), * tmp_strand_list = NULL;

	for(xk1 = 0 ; xk1 < genes; xk1++)
	{
		int xk2;
		
		memset(is_occupied,0,gene_exons_pointer[xk1]);
		tmp_chr_list = NULL;
		tmp_start_list = NULL;
		tmp_end_list = NULL;
		tmp_strand_list = NULL;
		out_chr_list[0]=0;
		out_start_list[0]=0;
		out_end_list[0]=0;
		out_strand_list[0]=0;
		int gene_nonoverlap_len =0;

		for(xk2=0; xk2<gene_exons_pointer[xk1]; xk2++)
		{
			if(!is_occupied[xk2])
			{
				int xk3;
				char * matched_chr = gene_exons_chr[xk2 + gene_exons_number[xk1]];
				char matched_strand = gene_exons_strand[xk2 + gene_exons_number[xk1]];

				memset(input_start_stop_list, 0, gene_exons_pointer[xk1] * sizeof(int) * 2);
				int gap_merge_ptr = 1;
				input_start_stop_list[0] = gene_exons_start[xk2 + gene_exons_number[xk1]];
				input_start_stop_list[1] = gene_exons_end[xk2 + gene_exons_number[xk1]] + 1;

				for(xk3 = xk2+1; xk3 < gene_exons_pointer[xk1]; xk3++)
				{
					if((!is_occupied[xk3]) && strcmp(matched_chr, gene_exons_chr[xk3+gene_exons_number[xk1]])==0 && matched_strand == gene_exons_strand[xk3 + gene_exons_number[xk1]])
					{
						is_occupied[xk3]=1;
						input_start_stop_list[gap_merge_ptr*2] = gene_exons_start[xk3+gene_exons_number[xk1]]; 
						input_start_stop_list[gap_merge_ptr*2+1] = gene_exons_end[xk3+gene_exons_number[xk1]]+1;
						gap_merge_ptr++;
					}
				}

				{
						int merged_gaps = mergeIntervals(input_start_stop_list, output_start_stop_list, gap_merge_ptr);

						for(xk3=0; xk3<merged_gaps; xk3++)
						{
							char numbbuf[12];
							BUFstrcat(out_chr_list, matched_chr, &tmp_chr_list);
							BUFstrcat(out_chr_list, ";", &tmp_chr_list);

							sprintf(numbbuf,"%u;", output_start_stop_list[xk3 * 2]);
							BUFstrcat(out_start_list, numbbuf, &tmp_start_list);
							sprintf(numbbuf,"%u;", output_start_stop_list[xk3 * 2 + 1] - 1);
							BUFstrcat(out_end_list, numbbuf, &tmp_end_list);
							sprintf(numbbuf,"%c;", matched_strand?'-':'+');
							BUFstrcat(out_strand_list, numbbuf, &tmp_strand_list);

							gene_nonoverlap_len += output_start_stop_list[xk3 * 2 + 1] - output_start_stop_list[xk3 * 2];
						}
				}	
			}
		}

		unsigned char * gene_symbol = global_context -> gene_name_array [xk1];

		#define _cut_tail(x) (x)[strlen(x)-1]=0

		_cut_tail(out_chr_list);
		_cut_tail(out_start_list);
		_cut_tail(out_end_list);
		_cut_tail(out_strand_list);

		fprintf(fp_out, "%s\t%s\t%s\t%s\t%s\t%d"    , gene_symbol, out_chr_list, out_start_list, out_end_list, out_strand_list, gene_nonoverlap_len);

		// all exons: gene_exons_number[xk1] : gene_exons_pointer[xk1]
		for(i_files=0; i_files< non_empty_files; i_files++)
		{
			if(column_numbers[i_files])
			{
				read_count_type_t longlong_res = 0;
				double double_res = 0;
				int is_double_number = calc_float_fraction(gene_columns[i_files+non_empty_files*xk1], &longlong_res, &double_res);
				if(is_double_number){
					fprintf(fp_out,"\t%.2f", double_res);
				}else{
					fprintf(fp_out,"\t%llu", longlong_res);
				}
			}
		}
		fprintf(fp_out,"\n");

	}
	free(is_occupied);
	free(input_start_stop_list);
	free(output_start_stop_list);
	free(out_chr_list);
	free(out_strand_list);
	free(out_start_list);
	free(out_end_list);

	free(gene_exons_number);
	free(gene_exons_pointer);
	free(gene_columns);
	free(gene_exons_chr);
	free(gene_exons_start);
	free(gene_exons_end);
	free(gene_exons_strand);
	fclose(fp_out);
}

void fc_write_final_counts(fc_thread_global_context_t * global_context, const char * out_file, int nfiles, char * file_list, read_count_type_t ** column_numbers, fc_read_counters *read_counters, int isCVersion)
{
	char fname[300];
	int i_files, xk1;

	sprintf(fname, "%s.summary", out_file);
	FILE * fp_out = f_subr_open(fname,"w");

	if(!fp_out){
		SUBREADprintf("Unable to create summary file '%s'\n", fname);
		return;
	}

	fprintf(fp_out,"Status");
	char * next_fn = file_list;
	
	for(i_files=0; i_files<nfiles; i_files++)
	{
		if(!next_fn||strlen(next_fn)<1) break;
		if(column_numbers[i_files])
			fprintf(fp_out,"\t%s", next_fn);

		next_fn += strlen(next_fn)+1;
	}

	fprintf(fp_out,"\n");
	char * keys [] ={ "Assigned" , "Unassigned_Ambiguity", "Unassigned_MultiMapping" ,"Unassigned_NoFeatures", "Unassigned_Unmapped", "Unassigned_MappingQuality", "Unassigned_FragmentLength", "Unassigned_Chimera", "Unassigned_Secondary", "Unassigned_Nonjunction", "Unassigned_Duplicate"};

	for(xk1=0; xk1<11; xk1++)
	{
		fprintf(fp_out,"%s", keys[xk1]);
		for(i_files = 0; i_files < nfiles; i_files ++)
		{
			unsigned long long * array_0 = (unsigned long long *)&(read_counters[i_files]);
			unsigned long long * cntr = array_0 + xk1;
			if(column_numbers[i_files])
				fprintf(fp_out,"\t%llu", *cntr);
		}
		fprintf(fp_out,"\n");
	}


	fclose(fp_out);
}
void fc_write_final_results(fc_thread_global_context_t * global_context, const char * out_file, int features, read_count_type_t ** column_numbers, char * file_list, int n_input_files, fc_feature_info_t * loaded_features, int header_out)
{
	/* save the results */
	FILE * fp_out;
	int i, i_files = 0;
	fp_out = f_subr_open(out_file,"w");
	if(!fp_out){
		SUBREADprintf("Failed to create file %s\n", out_file);
			return;
		}

	if(header_out)
	{
		fprintf(fp_out, "# Program:featureCounts v%s", SUBREAD_VERSION);
		if(global_context->cmd_rebuilt)
			fprintf(fp_out, "; Command:%s", global_context->cmd_rebuilt);
		fprintf(fp_out, "\n");
	}



	char * tmp_ptr = NULL, * next_fn;
	fprintf(fp_out,"Geneid\tChr\tStart\tEnd\tStrand\tLength");
	next_fn = strtok_r(file_list, ";", &tmp_ptr);
	while(1){
		if(!next_fn||strlen(next_fn)<1) break;
		if(column_numbers[i_files])
			fprintf(fp_out,"\t%s", next_fn);
		next_fn = strtok_r(NULL, ";", &tmp_ptr);
		i_files++;
	}
	fprintf(fp_out,"\n");
	for(i=0;i<features;i++)
	{
		fprintf(fp_out,"%s\t%s\t%u\t%u\t%c\t%d", global_context -> unistr_buffer_space + loaded_features[i].feature_name_pos,
 							   global_context -> unistr_buffer_space + loaded_features[i].feature_name_pos + loaded_features[i].chro_name_pos_delta,
						   	   loaded_features[i].start, loaded_features[i].end, loaded_features[i].is_negative_strand?'-':'+',loaded_features[i].end-loaded_features[i].start+1);
		for(i_files=0; i_files<n_input_files; i_files++)
		{
			if(column_numbers[i_files])
			{
				int sorted_exon_no = loaded_features[i].sorted_order;
				unsigned long long count_frac_raw =  column_numbers[i_files][sorted_exon_no], longlong_res = 0;

				double double_res = 0;
				int is_double_number = calc_float_fraction(count_frac_raw, &longlong_res, &double_res);
				if(is_double_number){
					fprintf(fp_out,"\t%.2f", double_res);
				}else{
					fprintf(fp_out,"\t%llu", longlong_res);
				}


			}
		}
		fprintf(fp_out,"\n");
	}

	fclose(fp_out);
}

static struct option long_options[] =
{
	{"primary",no_argument, 0, 0},
	{"readExtension5", required_argument, 0, 0},
	{"readExtension3", required_argument, 0, 0},
	{"read2pos", required_argument, 0, 0},
	{"minOverlap", required_argument, 0, 0},
	{"countSplitAlignmentsOnly", no_argument, 0, 0},
	{"debugCommand", required_argument, 0, 0},
	{"ignoreDup", no_argument, 0, 0},
	{"donotsort", no_argument, 0, 0},
	{"fraction", no_argument, 0, 0},
	{"order", required_argument, 0, 'S'},
	{"fasta", required_argument, 0, 0},
	{"largestOverlap", no_argument, 0,0},
	{0, 0, 0, 0}
};

void print_usage()
{
	SUBREADprintf("\nVersion %s\n\n", SUBREAD_VERSION);

	SUBREADputs("\nUsage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... \n");
	SUBREADputs("    Required parameters:\n"); 
	SUBREADputs("    -a <input>\tGive the name of the annotation file. The program assumes"); 
	SUBREADputs("              \tthat the provided annotation file is in GTF format. Use -F"); 
	SUBREADputs("              \toption to specify other annotation formats."); 
	SUBREADputs("    "); 
	SUBREADputs("    -o <input>\tGive the name of the output file. The output file contains"); 
	SUBREADputs("              \tthe number of reads assigned to each meta-feature (or each"); 
	SUBREADputs("              \tfeature if -f is specified). A meta-feature is the aggregation");
	SUBREADputs("              \tof features, grouped by using gene identifiers. Please refer");
	SUBREADputs("              \tto the users guide for more details."); 
	SUBREADputs("    "); 
	SUBREADputs("   input_files\tGive the names of input read files that include the read");
	SUBREADputs("              \tmapping results. Format of input files is automatically");
	SUBREADputs("              \tdetermined (SAM or BAM). Paired-end reads will be");
	SUBREADputs("              \tautomatically re-ordered if it is found that reads from the");
	SUBREADputs("              \tsame pair are not adjacent to each other. Multiple files can");
	SUBREADputs("              \tbe provided at the same time."); 
	SUBREADputs("    "); 
	SUBREADputs("    Optional parameters:"); 
	SUBREADputs("    "); 
	SUBREADputs("    -A <input>\tSpecify the name of a file including aliases of chromosome");
	SUBREADputs("              \tnames. The file should be a comma delimited text file that");
	SUBREADputs("              \tincludes two columns. The first column gives the chromosome");
	SUBREADputs("              \tnames used in the annotation and the second column gives the");
	SUBREADputs("              \tchromosome names used by reads. This file should not contain");
	SUBREADputs("              \theader lines. Names included in this file are case sensitive.");
	SUBREADputs("    "); 
	SUBREADputs("    -F <input>\tSpecify the format of the annotation file. Acceptable formats");
	SUBREADputs("              \tinclude `GTF' and `SAF'. `GTF' by default. Please refer to the");
	SUBREADputs("              \tusers guide for SAF annotation format."); 
	SUBREADputs("    "); 
	SUBREADputs("    -t <input>\tSpecify the feature type. Only rows which have the matched"); 
	SUBREADputs("              \tmatched feature type in the provided GTF annotation file"); 
	SUBREADputs("              \twill be included for read counting. `exon' by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -g <input>\tSpecify the attribute type used to group features (eg. exons)");
	SUBREADputs("              \tinto meta-features (eg. genes), when GTF annotation is provided.");
	SUBREADputs("              \t`gene_id' by default. This attribute type is usually the gene");
	SUBREADputs("              \tidentifier. This argument is useful for the meta-feature level");
	SUBREADputs("              \tsummarization.");
	SUBREADputs("    "); 
	SUBREADputs("    -f        \tIf specified, read summarization will be performed at the "); 
	SUBREADputs("              \tfeature level (eg. exon level). Otherwise, it is performed at");
	SUBREADputs("              \tmeta-feature level (eg. gene level).");
	SUBREADputs("    "); 
	SUBREADputs("    -O        \tIf specified, reads (or fragments if -p is specified) will"); 
	SUBREADputs("              \tbe allowed to be assigned to more than one matched meta-"); 
	SUBREADputs("              \tfeature (or feature if -f is specified). "); 
	SUBREADputs("    "); 
	SUBREADputs("    -s <int>  \tIndicate if strand-specific read counting should be performed.");
	SUBREADputs("              \tIt has three possible values:  0 (unstranded), 1 (stranded) and");
	SUBREADputs("              \t2 (reversely stranded). 0 by default.");
	SUBREADputs("    "); 
	SUBREADputs("    -M        \tIf specified, multi-mapping reads/fragments will also be counted");
	SUBREADputs("              \t(ie. a multi-mapping read will be counted up to N times if it");
	SUBREADputs("              \thas N reported mapping locations). The program uses the `NH' tag");
	SUBREADputs("              \tto find multi-mapping reads.");
	SUBREADputs("    "); 
	SUBREADputs("    -Q <int>  \tThe minimum mapping quality score a read must satisfy in order");
	SUBREADputs("              \tto be counted. For paired-end reads, at least one end should");
	SUBREADputs("              \tsatisfy this criteria. 0 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -T <int>  \tNumber of the threads. 1 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -v        \tOutput version of the program.");
	SUBREADputs("    "); 
	SUBREADputs("    -R        \tOutput read counting result for each read/fragment. For each");
	SUBREADputs("              \tinput read file, read counting results for reads/fragments will");
	SUBREADputs("              \tbe saved to a tab-delimited file that contains four columns");
	SUBREADputs("              \tincluding read name, status(assigned or the reason if not");
	SUBREADputs("              \tassigned), name of target feature/meta-feature and number of");
	SUBREADputs("              \thits if the read/fragment is counted multiple times. Name of");
	SUBREADputs("              \tthe file is the same as name of the input read file except a");
	SUBREADputs("              \tsuffix `.featureCounts' is added.");
	SUBREADputs("    "); 
	SUBREADputs("    --largestOverlap            If specified, reads (or fragments) will be");
	SUBREADputs("              \tassigned to the target that has the largest number of");
	SUBREADputs("              \toverlapping bases.");
	SUBREADputs("    "); 
	SUBREADputs("    --minOverlap <int>          Specify the minimum required number of");
	SUBREADputs("              \toverlapping bases between a read (or a fragment) and a feature.");
	SUBREADputs("              \t1 by default. If a negative value is provided, the read will be");
	SUBREADputs("              \textended from both ends.");
	SUBREADputs("    "); 
	SUBREADputs("    --readExtension5 <int>      Reads are extended upstream by <int> bases from");
	SUBREADputs("              \ttheir 5' end."); 
	SUBREADputs("    "); 
	SUBREADputs("    --readExtension3 <int>      Reads are extended upstream by <int> bases from");
	SUBREADputs("              \ttheir 3' end."); 
	SUBREADputs("    "); 
	SUBREADputs("    --read2pos <5:3>            The read is reduced to its 5' most base or 3'");
	SUBREADputs("              \tmost base. Read summarization is then performed based on the");
	SUBREADputs("              \tsingle base which the read is reduced to."); 
	SUBREADputs("    "); 
	SUBREADputs("    --fraction\tIf specified, a fractional count 1/n will be generated for each");
	SUBREADputs("              \tmulti-mapping read, where n is the number of alignments (indica-");
	SUBREADputs("              \tted by 'NH' tag) reported for the read. This option must be used");
	SUBREADputs("              \ttogether with the '-M' option.");
	SUBREADputs("    "); 
	SUBREADputs("    --primary \tIf specified, only primary alignments will be counted. Primary");
	SUBREADputs("              \tand secondary alignments are identified using bit 0x100 in the");
	SUBREADputs("              \tFlag field of SAM/BAM files. All primary alignments in a dataset");
	SUBREADputs("              \twill be counted no matter they are from multi-mapping reads or");
	SUBREADputs("              \tnot ('-M' is ignored). ");
	SUBREADputs("    "); 
	SUBREADputs("    --ignoreDup                 If specified, reads that were marked as");
	SUBREADputs("              \tduplicates will be ignored. Bit Ox400 in FLAG field of SAM/BAM");
	SUBREADputs("              \tfile is used for identifying duplicate reads. In paired end");
	SUBREADputs("              \tdata, the entire read pair will be ignored if at least one end");
	SUBREADputs("              \tis found to be a duplicate read.");
	SUBREADputs("    "); 
	SUBREADputs("    --countSplitAlignmentsOnly  If specified, only split alignments (CIGAR");
	SUBREADputs("              \tstrings containing letter `N') will be counted. All the other");
	SUBREADputs("              \talignments will be ignored. An example of split alignments is");
	SUBREADputs("              \tthe exon-spanning reads in RNA-seq data.");
	SUBREADputs("    "); 
	SUBREADputs("    Optional paired-end parameters:"); 
	SUBREADputs("    "); 
	SUBREADputs("    -p        \tIf specified, fragments (or templates) will be counted instead");
	SUBREADputs("              \tof reads. This option is only applicable for paired-end reads.");
	SUBREADputs("              \tThe two reads from the same fragment must be adjacent to each");
	SUBREADputs("              \tother in the provided SAM/BAM file.");
	SUBREADputs("    "); 
	SUBREADputs("    -P        \tIf specified, paired-end distance will be checked when assigning");
	SUBREADputs("              \tfragments to meta-features or features. This option is only");
	SUBREADputs("              \tapplicable when -p is specified. The distance thresholds should");
	SUBREADputs("              \tbe specified using -d and -D options."); 
	SUBREADputs("    "); 
	SUBREADputs("    -d <int>  \tMinimum fragment/template length, 50 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -D <int>  \tMaximum fragment/template length, 600 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -B        \tIf specified, only fragments that have both ends "); 
	SUBREADputs("              \tsuccessfully aligned will be considered for summarization."); 
	SUBREADputs("              \tThis option is only applicable for paired-end reads."); 
	SUBREADputs("    "); 
	SUBREADputs("    -S <ff:fr:rf> Orientation of the two read from the same pair, 'fr' by");
	SUBREADputs("              \tby default.");
	SUBREADputs("    "); 
	SUBREADputs("    -C        \tIf specified, the chimeric fragments (those fragments that "); 
	SUBREADputs("              \thave their two ends aligned to different chromosomes) will"); 
	SUBREADputs("              \tNOT be included for summarization. This option is only "); 
	SUBREADputs("              \tapplicable for paired-end read data."); 
	SUBREADputs("    "); 
	SUBREADputs("    --donotsort   If specified, paired end reads will not be reordered even if");
        SUBREADputs("              \treads from the same pair were found not to be next to each other");
        SUBREADputs("              \tin the input. ");
	SUBREADputs("    "); 
}

int junckey_sort_compare(void * inptr, int i, int j){
	char ** inp = (char **) inptr;
	return strcmp(inp[i], inp[j]);
}

void junckey_sort_exchange(void * inptr, int i, int j){

	char ** inp = (char **) inptr;
	char * tmpp = inp[j];
	inp[j]=inp[i];
	inp[i]=tmpp;
}

void junckey_sort_merge(void * inptr, int start, int items1, int items2){
	char ** inp = (char **) inptr;
	char ** tmpp = malloc(sizeof(char *) * (items1+items2));
	int read_1_ptr = start, read_2_ptr = start+items1, outptr = 0;
	while(1){
		if(read_1_ptr == start+items1 && read_2_ptr == start+items1+items2) break;
		if((read_1_ptr == start+items1)||(read_2_ptr < start+items1+items2 &&  junckey_sort_compare(inptr, read_1_ptr, read_2_ptr) > 0 )) {
			// select 2
			tmpp[outptr++]=inp[read_2_ptr++];
		} else {
			// select 1
			tmpp[outptr++]=inp[read_1_ptr++];
		}
	}
	memcpy(inp + start, tmpp, sizeof(char *)*(items1+items2));
	free(tmpp);
}

int junccmp(fc_junction_gene_t * j1, fc_junction_gene_t * j2){
	if(strcmp( j1 -> gene_name, j2 -> gene_name ) == 0)
		return 0;
	return 1;
}


void fc_write_final_junctions(fc_thread_global_context_t * global_context,  char * output_file_name,  read_count_type_t ** table_columns, char * input_file_names, int n_input_files, HashTable ** junction_global_table_list, HashTable ** splicing_global_table_list){
	int infile_i;
	char * buffered_names = malloc(strlen(input_file_names)+1);
	strcpy(buffered_names, input_file_names);

	HashTable * merged_junction_table = HashTableCreate(156679);

	HashTableSetHashFunction(merged_junction_table,HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(merged_junction_table, NULL, NULL);
	HashTableSetKeyComparisonFunction(merged_junction_table, fc_strcmp_chro);

	HashTable * merged_splicing_table = HashTableCreate(156679);

	HashTableSetHashFunction(merged_splicing_table,HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(merged_splicing_table, NULL, NULL);
	HashTableSetKeyComparisonFunction(merged_splicing_table, fc_strcmp_chro);


	for(infile_i = 0 ; infile_i < n_input_files ; infile_i ++){
		if(!table_columns[infile_i]) continue;	// bad input file
		KeyValuePair * cursor;
		int bucket;
		for(bucket=0; bucket < splicing_global_table_list[infile_i]  -> numOfBuckets; bucket++)
		{
			cursor = splicing_global_table_list[infile_i] -> bucketArray[bucket];
			while (cursor)
			{
				char * ky = (char *)cursor -> key;
				unsigned int old_supp = HashTableGet(merged_splicing_table, ky) - NULL;
				old_supp += (cursor -> value - NULL);
				HashTablePut(merged_splicing_table, ky, NULL+old_supp);
				cursor = cursor -> next;
			}	
		}
	}

	for(infile_i = 0 ; infile_i < n_input_files ; infile_i ++){
		if(!table_columns[infile_i]) continue;	// bad input file
		KeyValuePair * cursor;
		int bucket;
		for(bucket=0; bucket < junction_global_table_list[infile_i]  -> numOfBuckets; bucket++)
		{
			cursor = junction_global_table_list[infile_i] -> bucketArray[bucket];
			while (cursor)
			{
				char * ky = (char *)cursor -> key;

				if(HashTableGet(merged_junction_table, ky)==NULL)
					HashTablePut(merged_junction_table, ky, NULL+1);
				cursor = cursor -> next;
			}	
		}
	}

	char ** key_list;
	key_list = malloc(sizeof(char *) * merged_junction_table -> numOfElements);

	KeyValuePair * cursor;
	int bucket, ky_i = 0;
	for(bucket=0; bucket < merged_junction_table -> numOfBuckets; bucket++){
		cursor = merged_junction_table -> bucketArray[bucket];
		while (cursor){
			char * ky = (char *)cursor -> key;

			key_list[ky_i ++] = ky;
			cursor = cursor -> next;
		}
	}

	merge_sort(key_list,  merged_junction_table -> numOfElements , junckey_sort_compare, junckey_sort_exchange, junckey_sort_merge);

	char outfname[300];
	sprintf(outfname, "%s.junctions", output_file_name);

	int max_junction_genes = 3000;
	char * gene_names = malloc(max_junction_genes * FEATURE_NAME_LENGTH), * gene_name_tail;
	fc_junction_gene_t ** ret_juncs_small = malloc(sizeof(fc_junction_gene_t *) * max_junction_genes);
	fc_junction_gene_t ** ret_juncs_large = malloc(sizeof(fc_junction_gene_t *) * max_junction_genes);
	fc_junction_gene_t ** junction_key_list = malloc(sizeof(fc_junction_gene_t *)* max_junction_genes * 2);
	unsigned int * junction_support_list = malloc(sizeof(int)* max_junction_genes * 2);
	unsigned char * junction_source_list = malloc(sizeof(char)* max_junction_genes * 2 );

	int ky_i1, ky_i2;
	FILE * ofp = fopen(outfname, "w");
	char * tmpp = NULL;

	//SUBREADprintf("%s\n", buffered_names);
	fprintf(ofp, "#PrimaryGene\tSecondaryGenes\tChro1\tSplicePoint1\tStrand1\tChro2\tSplicePoint2\tStrand2");
	for(infile_i = 0 ; infile_i < n_input_files ; infile_i ++){
		char *fname;
		if(0 == infile_i) fname = strtok_r(buffered_names, ";", &tmpp);
		else fname = strtok_r(NULL, ";", &tmpp);
		
		if(!table_columns[infile_i]) continue;

		fprintf(ofp, "\t%s", fname);
	}
	fprintf(ofp, "\n");

	for(ky_i = 0; ky_i < merged_junction_table -> numOfElements ; ky_i ++){

		//SUBREADprintf("KY=%s\n", key_list[ky_i]);

		int unique_junctions = 0;
		char * chro_small = strtok_r( key_list[ky_i] , "\t", &tmpp);
		char * pos_small_str = strtok_r( NULL, "\t", &tmpp);
		char * chro_large = strtok_r( NULL, "\t", &tmpp);
		char * pos_large_str = strtok_r( NULL, "\t", &tmpp);

		unsigned int pos_small = atoi(pos_small_str);
		unsigned int pos_large = atoi(pos_large_str);

		int found_features_small = locate_junc_features(global_context, chro_small, pos_small, ret_juncs_small , max_junction_genes); 
		int found_features_large = locate_junc_features(global_context, chro_large, pos_large, ret_juncs_large , max_junction_genes);

		char strand = '?';
		if(global_context -> fasta_contigs){
			char donor[3], receptor[3];
			donor[2]=receptor[2]=0;
			int has = !get_contig_fasta(global_context -> fasta_contigs, chro_small, pos_small, 2, donor);
			has = has && !get_contig_fasta(global_context -> fasta_contigs, chro_large, pos_large-3, 2, receptor);
			if(has){
				if(donor[0]=='G' && donor[1]=='T' && receptor[0]=='A' && receptor[1]=='G') strand = '+';
				else if(donor[0]=='C' && donor[1]=='T' && receptor[0]=='A' && receptor[1]=='C') strand = '-';
			}
		}

		//SUBREADprintf("FOUND=%d, %d\n", found_features_small, found_features_large);

		gene_name_tail = gene_names;
		gene_names[0]=0;

		// rules to choose the primary gene:
		// (1) if some genes have one support but the other have multiple supporting reads: remove the lowly supported genes
		// (2) if all genes have only one support but from different ends of the fragment, then remove the genes that are assigned to the end having lower supporting fragments
		// (3) choose the gene that have the smallest coordinate.

		int max_supp = 0;
		for(ky_i1 = 0; ky_i1 < found_features_small + found_features_large; ky_i1++){
			int is_duplicate = 0;
			fc_junction_gene_t * tested_key = (ky_i1 < found_features_small)?ret_juncs_small[ky_i1] :ret_juncs_large[ky_i1 - found_features_small];
			for(ky_i2 = 0; ky_i2 < unique_junctions; ky_i2 ++){
				if(junccmp( tested_key, junction_key_list[ky_i2]  )==0){
					junction_support_list[ ky_i2 ] ++;
					junction_source_list[ky_i2] |= ( (ky_i1 < found_features_small)? 1 : 2 );
					is_duplicate = 1;
					break;
				}
			}

			if(!is_duplicate){
				junction_key_list[unique_junctions] = tested_key;
				junction_support_list[unique_junctions] = 1;
				junction_source_list[unique_junctions] = ( (ky_i1 < found_features_small)? 1 : 2 );
				max_supp = max(junction_support_list[unique_junctions], max_supp);
				unique_junctions++;
			}
		}

		if(1 == max_supp){
			if(found_features_small > 0 && found_features_large > 0){
				char junc_key [FEATURE_NAME_LENGTH + 15]; 
				sprintf(junc_key, "%s\t%u", chro_small, pos_small);
				unsigned int supp_small = HashTableGet(merged_splicing_table, junc_key) - NULL;
				sprintf(junc_key, "%s\t%u", chro_large, pos_large);
				unsigned int supp_large = HashTableGet(merged_splicing_table, junc_key) - NULL;

				if(supp_small !=supp_large){
					for(ky_i2 = 0; ky_i2 < unique_junctions; ky_i2 ++){
						if(supp_small > supp_large && junction_source_list[ky_i2] == 1) junction_key_list[ky_i2] = NULL;
						else if(supp_small < supp_large && junction_source_list[ky_i2] == 2) junction_key_list[ky_i2] = NULL;
					}
				} 
			}
		}

		int smallest_coordinate_gene = 0x7fffffff;
		fc_junction_gene_t * primary_gene = NULL;
		
		for(ky_i2 = 0; ky_i2 < unique_junctions; ky_i2 ++){
			fc_junction_gene_t * tested_key = junction_key_list[ky_i2];
			if(tested_key != NULL && tested_key -> pos_first_base < smallest_coordinate_gene){
				primary_gene = tested_key;
				smallest_coordinate_gene = tested_key -> pos_first_base;
			}
		}

		if(primary_gene == NULL){
			strcpy(gene_names, "NONE");
		}else{
			strcpy(gene_names, primary_gene -> gene_name);
		}

		*(pos_small_str-1)='\t';
		*(pos_large_str-1)='\t';

		fprintf(ofp, "%s", gene_names);
	
		gene_name_tail = gene_names;
		gene_names[0]=0;
		for(ky_i2 = 0; ky_i2 < unique_junctions; ky_i2 ++){
			fc_junction_gene_t * tested_key = junction_key_list[ky_i2];
			if(tested_key && tested_key != primary_gene)
				gene_name_tail += sprintf(gene_name_tail, "%s,", tested_key -> gene_name);
		}
		if( gene_names[0] ) gene_name_tail[-1]=0;
		else strcpy(gene_names, "NONE");
		fprintf(ofp, "\t%s", gene_names);

		fprintf(ofp, "\t%s\t%c\t%s\t%c", chro_small, strand, chro_large, strand);

		chro_large[-1]='\t';

		for(infile_i = 0 ; infile_i < n_input_files ; infile_i ++){
			if(!table_columns[infile_i]) continue;
			unsigned long count = HashTableGet(junction_global_table_list[infile_i]  , key_list[ky_i]) - NULL;
			fprintf(ofp,"\t%lu", count);
		}
		fprintf(ofp, "\n");
	}
	fclose(ofp);
	free(junction_key_list);
	free(gene_names);
	free(ret_juncs_small);
	free(ret_juncs_large);
	free(buffered_names);
	free(junction_support_list);
	free(key_list);
	free(junction_source_list);
	HashTableDestroy(merged_junction_table);
	HashTableDestroy(merged_splicing_table);
}

int readSummary_single_file(fc_thread_global_context_t * global_context, read_count_type_t * column_numbers, int nexons,  int * geneid, char ** chr, long * start, long * stop, unsigned char * sorted_strand, char * anno_chr_2ch, char ** anno_chrs, long * anno_chr_head, long * block_end_index, long * block_min_start , long * block_max_end, fc_read_counters * my_read_counter, HashTable * junc_glob_tab, HashTable * splicing_glob_tab);

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
	19: min.MappingQualityScore
	20: as.numeric(isMultiMappingAllowed)
	21: Annotation Chromosome Alias Name File. If the file is not specified, set this value to NULL or a zero-length string.
	22: Command line for CfeatureCounts header output; RfeatureCounts should set this value to NULL or a zero-length string or a space (' ').
	23: as.numeric(isInputFileResortNeeded)
	24: NOT IN USE: as.numeric(feature_block_size) # This parameter is no longer used. Give "14" for safe. 
	25: as.numeric(Five_End_Extension_Length)  # 5' end extension
	26: as.numeric(Three_End_Extension_Length)  # 3' end extension
	27: as.numeric(Minimum_Overlap_Between_Read_And_Feature) # 1 by default
	28: as.numeric(is_Split_Alignment_Only) # 0 by default
	29: as.numeric(reduce_5_3_ends_to_one) # 0= no reduction; 1= reduce to 5' end; 2= reduce to 3' end
	30: debug_command # This is for debug only; RfeatureCounts should pass a space (" ") to this parameter, disabling the debug command.
	31: as.numeric(is_duplicate_ignored) # 0 = INCLUDE DUPLICATE READS; 1 = IGNORE DUPLICATE READS (0x400 FLAG IS SET) ; "0" by default.
	32: as.numeric(do_not_sort)   # 1 = NEVER SORT THE PE BAM/SAM FILES; 0 = SORT THE BAM/SAM FILE IF IT IS FOUND NOT SORTED.
	33: as.numeric(fractionMultiMapping) # 1 = calculate fraction numbers if a read overlaps with multiple features or meta-features. "-M" must be specified when fractions are caculated.
	34: as.numeric(useOverlappingBreakTie) # 1 = Select features or meta-features with a longer overlapping length; 0 = just use read-voting strategy: one overlapping read = 1 vote
	35: Pair_Orientations # FF, FR, RF or RR. This parameter matters only if "-s" option is 1 or 2.
	36: as.numeric(doJunctionCounting)  # 1 = count the number of junction reads spaining each exon-exon pairs;  0 = do not.
	37: file name of fasta (for determine the strandness of junctions by looking for GT/AG or CT/AC).
	 */

	int isStrandChecked, isCVersion, isChimericDisallowed, isPEDistChecked, minMappingQualityScore=0, isInputFileResortNeeded, feature_block_size = 20, reduce_5_3_ends_to_one;
	char **chr;
	long *start, *stop;
	int *geneid;

	char *nameFeatureTypeColumn, *nameGeneIDColumn,*debug_command, *pair_orientations;
	long nexons;


	long * anno_chr_head, * block_min_start, *block_max_end, *block_end_index;
	char ** anno_chrs, * anno_chr_2ch;
	long curchr, curpos;
	char * curchr_name, * fasta_contigs_fname;
	unsigned char * sorted_strand;
	curchr = 0;
	curpos = 0;
	curchr_name = "";

	int isPE, minPEDistance, maxPEDistance, isReadSummaryReport, isBothEndRequired, isMultiMappingAllowed, fiveEndExtension, threeEndExtension, minFragmentOverlap, isSplitAlignmentOnly, is_duplicate_ignored, doNotSort, fractionMultiMapping, useOverlappingBreakTie, doJuncCounting;

	int isSAM, isGTF, n_input_files=0;
	char *  alias_file_name = NULL, * cmd_rebuilt = NULL;

	int isMultiOverlapAllowed, isGeneLevel;

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
	if(argc > 19)
		minMappingQualityScore = atoi(argv[19]);
	else	minMappingQualityScore = 0;
	if(argc > 20)
		isMultiMappingAllowed = atoi(argv[20]);
	else	isMultiMappingAllowed = 1;
	if(argc > 21)
	{
		alias_file_name = argv[21];
		if(alias_file_name == NULL || alias_file_name[0]==' ' || alias_file_name[0]==0)
			alias_file_name = NULL;
	}
	else	alias_file_name = NULL;
	if(argc > 22)
	{
		cmd_rebuilt = argv[22];
		if(cmd_rebuilt == NULL || cmd_rebuilt[0]==' '||cmd_rebuilt[0]==0)
			cmd_rebuilt=NULL;
	}
	else	cmd_rebuilt = NULL;
	if(argc>23)
		isInputFileResortNeeded = atoi(argv[23]);
	else	isInputFileResortNeeded = 0;
	if(thread_number<1) thread_number=1;
	if(thread_number>16)thread_number=16;

	int Param_fiveEndExtension, Param_threeEndExtension;
	if(argc>25)
		Param_fiveEndExtension = atoi(argv[25]);
	else    Param_fiveEndExtension = 0;

	if(argc>26)
		Param_threeEndExtension = atoi(argv[26]);
	else    Param_threeEndExtension = 0;

	if(argc>27)
		minFragmentOverlap = atoi(argv[27]);
	else    minFragmentOverlap = 1;

	if(minFragmentOverlap <1){
		fiveEndExtension = 1 - minFragmentOverlap;
		threeEndExtension = 1 - minFragmentOverlap;
		minFragmentOverlap = 1;
	}else{
		fiveEndExtension = Param_fiveEndExtension;
		threeEndExtension = Param_threeEndExtension;
	}

	if(argc>28)
		isSplitAlignmentOnly = atoi(argv[28]);
	else	isSplitAlignmentOnly = 0;

	if(argc>29)
		reduce_5_3_ends_to_one = atoi(argv[29]);	// 0 : no reduce; 1: reduce to 5' end; 2: reduce to 3' end.
	else	reduce_5_3_ends_to_one = 0;


	if(argc>30 && strlen(argv[30])>0 && argv[30][0]!=' ')
		debug_command = argv[30];
	else
		debug_command = " ";

	if(argc>31)
		is_duplicate_ignored = atoi(argv[31]);
	else
		is_duplicate_ignored = 0;

	if(argc>32)
		doNotSort = atoi(argv[32]);
	else
		doNotSort = 0;

	if(argc>33)
		fractionMultiMapping = atoi(argv[33]);
	else
		fractionMultiMapping = 0;

	if(argc>34)
		useOverlappingBreakTie = atoi(argv[34]);
	else	useOverlappingBreakTie = 0;


	if(argc>35)
		pair_orientations = argv[35];
	else	pair_orientations = "FR";

	if(argc>36)
		doJuncCounting = atoi(argv[36]);
	else	doJuncCounting = 0;

	fasta_contigs_fname = NULL;
	if(argc>37)
		if(argv[37][0] != 0 && argv[37][0]!=' ')
			fasta_contigs_fname = argv[37];
	
	unsigned int buffer_size = 1024*1024*6;


	fc_thread_global_context_t global_context;
	fc_thread_init_global_context(& global_context, buffer_size, thread_number, MAX_LINE_LENGTH, isPE, minPEDistance, maxPEDistance,isGeneLevel, isMultiOverlapAllowed, isStrandChecked, (char *)argv[3] , isReadSummaryReport, isBothEndRequired, isChimericDisallowed, isPEDistChecked, nameFeatureTypeColumn, nameGeneIDColumn, minMappingQualityScore,isMultiMappingAllowed, isSAM, alias_file_name, cmd_rebuilt, isInputFileResortNeeded, feature_block_size, isCVersion, fiveEndExtension, threeEndExtension , minFragmentOverlap, isSplitAlignmentOnly, reduce_5_3_ends_to_one, debug_command, is_duplicate_ignored, doNotSort, fractionMultiMapping, useOverlappingBreakTie, pair_orientations, doJuncCounting);

	if( global_context.is_multi_mapping_allowed != ALLOW_ALL_MULTI_MAPPING && global_context.use_fraction_multi_mapping)
	{
		SUBREADprintf("ERROR: '--fraction' option should be used together with '-M'.\nThis option should not be used when multi-mapping reads are disallowed or the primary mapping is only counted.\n");
		return -1;
	}
	print_FC_configuration(&global_context, argv[1], argv[2], argv[3], global_context.is_SAM_file, isGTF, & n_input_files, isReadSummaryReport);


	//print_in_box(80,0,0,"IG=%d, IS=%d", isGeneLevel, isSplitAlignmentOnly);
	if(isSplitAlignmentOnly && ( isGeneLevel || !isMultiOverlapAllowed) )
	{
		print_in_box(80,0,0,"NOTICE --countSplitAlignmentsOnly is specified, but '-O' and '-f' are not");
		print_in_box(80,0,0,"       both specified. Please read the manual for details.");
		print_in_box(80,0,0,"");
	}

	// Loading the annotations.
	// Nothing is done if the annotation does not exist.
	fc_feature_info_t * loaded_features;
	print_in_box(84,0,0,"Load annotation file %s %c[0m...", argv[1], CHAR_ESC);
	nexons = load_feature_info(&global_context,argv[1], isGTF?FILE_TYPE_GTF:FILE_TYPE_RSUBREAD, &loaded_features);
	if(nexons<1){
		SUBREADprintf("Failed to open the annotation file %s, or its format is incorrect, or it contains no '%s' features.\n",argv[1], nameFeatureTypeColumn);
		return -1;
	}

	sort_feature_info(&global_context, nexons, loaded_features, &chr, &geneid, &start, &stop, &sorted_strand, &anno_chr_2ch, &anno_chrs, &anno_chr_head, & block_end_index, & block_min_start, & block_max_end);
	print_in_box(80,0,0,"   Meta-features : %d", global_context . gene_name_table -> numOfElements);
	print_in_box(80,0,0,"   Chromosomes/contigs : %d", global_context . exontable_nchrs);

	print_in_box(80,0,0,"");


	if(fasta_contigs_fname){
		print_in_box(80,0,0,"Loading FASTA contigs : %s", fasta_contigs_fname);
		global_context.fasta_contigs = malloc(sizeof(fasta_contigs_t));
		read_contig_fasta(global_context.fasta_contigs, fasta_contigs_fname);
		print_in_box(80,0,0,"");
	}else	global_context.fasta_contigs = NULL;
	


	global_context.exontable_exons = nexons;
	unsigned int * nreads = (unsigned int *) calloc(nexons,sizeof(int));









	char * tmp_pntr = NULL;
	char * file_list_used = malloc(strlen(argv[2])+1);
	strcpy(file_list_used, argv[2]);
	char * next_fn = strtok_r(file_list_used,";", &tmp_pntr);
	read_count_type_t ** table_columns = calloc( n_input_files , sizeof(read_count_type_t *)), i_files=0;
	fc_read_counters * read_counters = calloc(n_input_files , sizeof(fc_read_counters)); 
	HashTable ** junction_global_table_list = NULL;
	HashTable ** splicing_global_table_list = NULL;

	if(global_context.do_junction_counting){
		junction_global_table_list = calloc(n_input_files, sizeof(HashTable *));
		splicing_global_table_list = calloc(n_input_files, sizeof(HashTable *));
	}

	for(;;){
		int redoing, original_sorting = global_context.is_input_file_resort_needed, orininal_isPE = global_context.is_paired_end_mode_assign;
		if(next_fn==NULL || strlen(next_fn)<1) break;

		read_count_type_t * column_numbers = calloc(nexons, sizeof(read_count_type_t));
		HashTable * junction_global_table = NULL;
		HashTable * splicing_global_table = NULL;

		strcpy(global_context.input_file_name, next_fn);
		strcpy(global_context.raw_input_file_name, next_fn);
		global_context.redo=0;
		
		for(redoing = 0; redoing < 1 + !original_sorting; redoing++)
		{

			if(global_context.do_junction_counting){
				junction_global_table = HashTableCreate(156679);
				splicing_global_table = HashTableCreate(156679);

				HashTableSetHashFunction(junction_global_table,HashTableStringHashFunction);
				HashTableSetDeallocationFunctions(junction_global_table, free, NULL);
				HashTableSetKeyComparisonFunction(junction_global_table, fc_strcmp_chro);

				HashTableSetHashFunction(splicing_global_table,HashTableStringHashFunction);
				HashTableSetDeallocationFunctions(splicing_global_table, free, NULL);
				HashTableSetKeyComparisonFunction(splicing_global_table, fc_strcmp_chro);
			}

			fc_read_counters * my_read_counter = &(read_counters[i_files]);
			memset(my_read_counter, 0, sizeof(fc_read_counters));

			int ret_int = readSummary_single_file(& global_context, column_numbers, nexons, geneid, chr, start, stop, sorted_strand, anno_chr_2ch, anno_chrs, anno_chr_head, block_end_index, block_min_start, block_max_end, my_read_counter, junction_global_table, splicing_global_table);
			if(ret_int!=0 || (global_context.redo && redoing)){
				// give up this file.

				table_columns[i_files] = NULL;
				if(global_context.do_junction_counting){
					HashTableDestroy(junction_global_table);
					HashTableDestroy(splicing_global_table);
				}
				free(column_numbers);
				break;
			} else {
				// maybe finished or need redoing
				table_columns[i_files] = column_numbers;
				if(global_context.do_junction_counting){
					junction_global_table_list[ i_files ] = junction_global_table;
					splicing_global_table_list[ i_files ] = splicing_global_table;
				}
			}

			if(redoing || !global_context.redo) break;
			
			global_context.is_input_file_resort_needed = 1;
			memset(column_numbers, 0, nexons * sizeof(read_count_type_t));

			if(global_context.do_junction_counting){
				HashTableDestroy(junction_global_table);
				HashTableDestroy(splicing_global_table);
			}
		}
		global_context.is_SAM_file = isSAM;
		global_context.is_input_file_resort_needed = original_sorting;
		global_context.is_paired_end_mode_assign = orininal_isPE;

		i_files++;
		next_fn = strtok_r(NULL, ";", &tmp_pntr);
	}

	free(file_list_used);

	if(isGeneLevel)
		fc_write_final_gene_results(&global_context, geneid, chr, start, stop, sorted_strand, argv[3], nexons,  table_columns, argv[2], n_input_files , loaded_features, isCVersion);
	else
		fc_write_final_results(&global_context, argv[3], nexons, table_columns, argv[2], n_input_files ,loaded_features, isCVersion);

	if(global_context.do_junction_counting)
		fc_write_final_junctions(&global_context, argv[3], table_columns, argv[2], n_input_files , junction_global_table_list, splicing_global_table_list);

	fc_write_final_counts(&global_context, argv[3], n_input_files, argv[2], table_columns, read_counters, isCVersion);

	int total_written_coulmns = 0;
	for(i_files=0; i_files<n_input_files; i_files++)
		if(table_columns[i_files]){
			free(table_columns[i_files]);
			if(global_context.do_junction_counting){
				HashTableDestroy(junction_global_table_list[i_files]);
				HashTableDestroy(splicing_global_table_list[i_files]);
			}

			total_written_coulmns++;

		}
	free(table_columns);


	print_FC_results(&global_context);
	KeyValuePair * cursor;
	int bucket;
	for(bucket=0; bucket < global_context.exontable_chro_table  -> numOfBuckets; bucket++)
	{
		cursor = global_context.exontable_chro_table -> bucketArray[bucket];
		while (1)
		{
			if (!cursor) break;
			fc_chromosome_index_info * del_chro_info = cursor->value;
			free(del_chro_info->reverse_table_start_index);
			//free(del_chro_info->reverse_table_end_index);
			free((void *)cursor -> key);
			free(del_chro_info);
			cursor = cursor->next;
		}
	}

	if(global_context.SAM_output_fp) fclose(global_context. SAM_output_fp);
	HashTableDestroy(global_context.gene_name_table);
	free(global_context.gene_name_array);

	HashTableDestroy(global_context.exontable_chro_table);
	if(global_context.fasta_contigs){
		destroy_contig_fasta(global_context.fasta_contigs);
		free(global_context.fasta_contigs);
	}
	if(global_context.annot_chro_name_alias_table)
		HashTableDestroy(global_context.annot_chro_name_alias_table);
	if(global_context.do_junction_counting){
		HashTableDestroy(global_context.junction_features_table);
		free(junction_global_table_list);
		free(splicing_global_table_list);
	}

	free(global_context.unistr_buffer_space);
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


	return total_written_coulmns?0:-1;
}













int readSummary_single_file(fc_thread_global_context_t * global_context, read_count_type_t * column_numbers, int nexons,  int * geneid, char ** chr, long * start, long * stop, unsigned char * sorted_strand, char * anno_chr_2ch, char ** anno_chrs, long * anno_chr_head, long * block_end_index, long * block_min_start , long * block_max_end, fc_read_counters * my_read_counter, HashTable * junction_global_table, HashTable * splicing_global_table)
{
	FILE *fp_in = NULL;
	int read_length = 0;
	int is_first_read_PE=0;
	char * line = (char*)calloc(MAX_LINE_LENGTH, 1);
	char * file_str = "";

	if(strcmp( global_context->input_file_name,"STDIN")!=0)
	{
		int file_probe = is_certainly_bam_file(global_context->input_file_name, &is_first_read_PE);

		
		global_context -> is_paired_end_input_file = is_first_read_PE;
		// a Singel-end SAM/BAM file cannot be assigned as a PE SAM/BAM file;
		// but a PE SAM/BAM file may be assigned as a SE file if the user wishes to do so.
		if(is_first_read_PE==0)
				global_context -> is_paired_end_mode_assign = 0;

		if(file_probe == 1){
			global_context->is_SAM_file = 0;
		}
		else if(file_probe == 0) global_context->is_SAM_file = 1;

		global_context -> start_time = miltime();

		file_str = "SAM";
		if(file_probe == 1) file_str = "BAM" ;
		if(file_probe == -1) file_str = "Unknown";

		if(!global_context->redo)
		{
			print_in_box(80,0,0,"Process %s file %s...", file_str, global_context->input_file_name);
			if(is_first_read_PE)
				print_in_box(80,0,0,"   Paired-end reads are included.");
			else
				print_in_box(80,0,0,"   Single-end reads are included.");
		}
		
	}

	int isInputFileResortNeeded = global_context->is_input_file_resort_needed;

	if(strcmp( global_context->input_file_name,"STDIN")!=0)
	{
		FILE * exist_fp = f_subr_open( global_context->input_file_name,"r");
		if(!exist_fp)
		{
			print_in_box(80,0,0,"Failed to open file %s",  global_context->input_file_name);
			print_in_box(80,0,0,"No counts were generated for this file.");
			print_in_box(80,0,0,"");
			return -1;
		}
		fclose(exist_fp);
	}

	if(strcmp(global_context->input_file_name,"STDIN")!=0)
		if(warning_file_type(global_context->input_file_name, global_context->is_SAM_file?FILE_TYPE_SAM:FILE_TYPE_BAM))
			global_context->is_unpaired_warning_shown=1;
	if(strcmp(global_context->input_file_name,"STDIN")!=0 && isInputFileResortNeeded)
		if(resort_input_file( global_context)) return -1;
	int isSAM = global_context->is_SAM_file;
	// Open the SAM/BAM file
	// Nothing is done if the file does not exist.

	#ifdef MAKE_STANDALONE
	if(strcmp("STDIN",global_context->input_file_name)==0)
		fp_in = stdin;
	else
		fp_in = f_subr_open(global_context->input_file_name,"r");
	#else
		fp_in = f_subr_open(global_context->input_file_name,"r");
	#endif


	// begin to load-in the data.
	if(!global_context->redo)
	{
		if( global_context->is_paired_end_mode_assign)
		{
			print_in_box(80,0,0,"   Assign fragments (read pairs) to features...");
//				print_in_box(80,0,0,"   Each fragment is counted no more than once.");
		}
		else
			print_in_box(80,0,0,"   Assign reads to features...");
	}



	fc_thread_start_threads(global_context, nexons, geneid, chr, start, stop, sorted_strand, anno_chr_2ch, anno_chrs, anno_chr_head, block_end_index, block_min_start , block_max_end, read_length);

	int buffer_pairs = global_context -> thread_number>1?512:1;
	int isPE = global_context->is_paired_end_mode_assign;
	char * preload_line = malloc(sizeof(char) * (2+MAX_LINE_LENGTH)*(isPE?2:1)*buffer_pairs);
	int preload_line_ptr;
	int current_thread_id = 0;
	fc_thread_thread_context_t * one_thread_context = global_context->thread_contexts;

	SamBam_Reference_Info * sb_header_tab = NULL;
	
	unsigned long long int chunk_id = 0;
	int binary_remainder = 0, binary_read_ptr = 0;
	char * chunk_in_buff = malloc(70000);
	char * binary_in_buff = malloc(80000 * 2);

	if(!isSAM)
	{
		int remainder_read_data_len = 0;

		PBum_load_header(fp_in, &sb_header_tab, binary_in_buff,  & remainder_read_data_len);
		//printf("RMD=%d\n", remainder_read_data_len);

		if(remainder_read_data_len)
		{
			binary_remainder = remainder_read_data_len;
		}
		global_context->sambam_chro_table = sb_header_tab;
	}

	while (1){
		int pair_no;
		int is_second_read;
		int fresh_read_no = 0;
		preload_line[0] = 0;
		preload_line_ptr = 0;

		char * ret = NULL;
		
		// one-thread BAM is not supported.
		if( isSAM && global_context->thread_number==1)
		{
			int is_second_read;

			for(is_second_read=0;is_second_read<(isPE?2:1);is_second_read++)
			{
				char * lbuf = is_second_read?one_thread_context -> line_buffer2:one_thread_context -> line_buffer1;
				while(1)
				{
					ret = fgets(lbuf, MAX_LINE_LENGTH, fp_in);
					if(global_context -> redo) ret = NULL;
					if(!ret) break;
					if(lbuf[0] == '@') 
					{
						int retlen = strlen(ret);
						if(ret[retlen-1]!='\n')
						{
							while(1){
								int nch = getc(fp_in);
								if(nch == EOF || nch == '\n') break;
							}
						}
					}
					else break;
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
					global_context->read_length = read_length;
				}
				lbuf[strlen(lbuf)+1]=0;
			}

			//printf("RRR=%d\n",ret);
			
			//one_thread_context -> current_read_length1 = global_context->read_length;
			//one_thread_context -> current_read_length2 = global_context->read_length;

			if(is_second_read == 1 && isPE){
				print_in_box(85,0,0,"   %c[31mThere are odd number of reads in the paired-end data.", CHAR_ESC);
				print_in_box(80,0,0,"   Please make sure that the format is correct.");
			}

			if(ret)
			{
				global_context->all_reads ++;
				process_line_buffer(global_context, one_thread_context);
			}

			if(!ret)break;
		}
		else if(!isSAM)
		{
			int no_of_reads = 0;
			unsigned int real_len = 0;
			// most of the data must have been given out before this step.

			int cdata_size = 0;

			if(binary_remainder > 70000)
				SUBREADprintf("SOMETHING IS WRONG!\n");

			if(global_context -> redo)
				cdata_size = -1;
			else{
				if(binary_remainder<10000)
					cdata_size = PBam_get_next_zchunk(fp_in, chunk_in_buff, 65537, & real_len);
			}

			if(cdata_size>0 || binary_remainder>0)
			{
				int x1;


				if(binary_read_ptr>0)
				{
					for(x1=0; x1< binary_remainder; x1++)
						binary_in_buff[x1] = binary_in_buff [x1 + binary_read_ptr];
					binary_read_ptr = 0;
				}

				//fprintf(stderr,"NBN=%d, OBN=%d\n", cdata_size , binary_remainder);
				if(cdata_size>0)
				{
					int new_binary_bytes = SamBam_unzip(binary_in_buff + binary_remainder , chunk_in_buff , cdata_size);
					if(new_binary_bytes>=0)
						binary_remainder += new_binary_bytes;
					else	SUBREADprintf("ERROR: BAM GZIP FORMAT ERROR.\n");
				//	fprintf(stderr,"BBN=%d\n", new_binary_bytes);
				}

				while(binary_remainder>4)
				{
					unsigned int binary_read_len = 0;
					memcpy(& binary_read_len , binary_in_buff + binary_read_ptr , 4);
					//printf("RLEN=%d; PTR=%d; RMD=%d\n", binary_read_len , binary_read_ptr, binary_remainder);
					if(binary_read_len > 10000)
					{
						binary_remainder = -1;
						//SUBREADprintf("FATAL ERROR: BAM RECORD SIZE = %u ; PTR=%d ; REM=%d.\n", binary_read_len, binary_read_ptr , binary_remainder);
						print_in_box(80,0,0,"   A format error was detected in this BAM file.");
						print_in_box(80,0,0,"   The remaining part in the file is skipped.");
						print_in_box(80,0,0,"   Please check the file format using samtools.");
						print_in_box(80,0,0,"");
						break;
					}
					// if the program runs on PE mode, no_of_reads must be even.

					if(isPE)
					{
						if(4 + binary_read_len + 4 < binary_remainder)
						{
							int binary_read2_len=0;
							memcpy(&binary_read2_len , binary_in_buff + binary_read_ptr + 4 + binary_read_len, 4);
							if(4 + binary_read_len + 4 + binary_read2_len <= binary_remainder)
							{
								no_of_reads +=2;
								binary_read_ptr += 4 + binary_read_len + 4 + binary_read2_len;
								binary_remainder -= 4 + binary_read_len + 4 + binary_read2_len;
							}
							else break;
						}
						else break;
					}
					else
					{
						if(binary_read_len + 4<= binary_remainder)
						{
							no_of_reads ++;
							binary_read_ptr  += 4 + binary_read_len;
							binary_remainder -= 4 + binary_read_len;
						}
						else break;
					}
				}
			}

			if(binary_remainder <0)break;

			if(no_of_reads>0) 
			{
				while(1)
				{
					int is_finished = 0;

					fc_thread_thread_context_t * thread_context = global_context->thread_contexts+current_thread_id;

					pthread_spin_lock(&thread_context->input_buffer_lock);

					// the number of bytes can be utilised given the two_chunk_len.
					int empty_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_remainder;
					int tail_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_write_ptr;
					
					if(thread_context->input_buffer_remainder > global_context->input_buffer_max_size)
					{
						SUBREADprintf("RMD=%d\n", thread_context->input_buffer_remainder );
						assert(0);
					}

					if(tail_bytes < binary_read_ptr + 4)
						empty_bytes -= tail_bytes;

					// copy the new buffer to thread buffer.
					// format: read_number=n, read_chunk1, read_chunk2, ..., read_chunk_n
					if(empty_bytes >= binary_read_ptr + 4)
					{

						if(tail_bytes < binary_read_ptr + 4)
						{
							if(tail_bytes>=4)
								memset(thread_context->input_buffer + thread_context->input_buffer_write_ptr, 0, 4);
							thread_context->input_buffer_write_ptr = 0;
							thread_context->input_buffer_remainder += tail_bytes;
						}

						memcpy( thread_context->input_buffer + thread_context->input_buffer_write_ptr, & binary_read_ptr, 4);
						memcpy( thread_context->input_buffer + thread_context->input_buffer_write_ptr + 4, binary_in_buff , binary_read_ptr);
						thread_context->input_buffer_write_ptr += 4 + binary_read_ptr;
						thread_context->input_buffer_remainder += 4 + binary_read_ptr;
						is_finished = 1;
					}

					pthread_spin_unlock(&thread_context->input_buffer_lock);
					current_thread_id++;
					if(current_thread_id >= global_context->thread_number) current_thread_id = 0;

					if(is_finished) break;
					else usleep(tick_time);

				}

				chunk_id++;
			}
			else if(cdata_size<0)
				break;
	
		}
		else
		{
			for(pair_no=0; pair_no < buffer_pairs; pair_no++)
			{
				for(is_second_read=0;is_second_read<(isPE?2:1);is_second_read++)
				{
					while(1)
					{
						ret = fgets(preload_line+preload_line_ptr, MAX_LINE_LENGTH, fp_in);
						if(global_context -> redo ) ret = NULL;
						if(!ret) break;
						if(preload_line[preload_line_ptr] == '@'){
							int retlen = strlen(ret);
							if(ret[retlen-1]!='\n')
							{
								while(1){
									int nch = getc(fp_in);
									if(nch == EOF || nch == '\n') break;
								}
							}
						}else break;
					}

					if(!ret) break;

					int curr_line_len = strlen(preload_line+preload_line_ptr);
					if(curr_line_len >= MAX_LINE_LENGTH || preload_line[preload_line_ptr + curr_line_len-1]!='\n')
					{
						print_in_box(80,0,0,"ERROR: the lines are too long. Please check the input format!!\n");
						ret = NULL;
						preload_line_ptr = 0;
						break;
					}
					preload_line_ptr += curr_line_len;

					fresh_read_no++;
				}
				if(!ret) break;
				else
					global_context->all_reads ++;
			}

			int line_length = preload_line_ptr;
			if(line_length >= global_context->input_buffer_max_size-1)
			{
				SUBREADprintf("ERROR: the lines are too long. Please check the input format!!\n");
				break;
			}
			if(isPE && (fresh_read_no%2>0))
			{
				// Safegarding -- it should not happen if the SAM file has a correct format.
				//line_length = 0;
				if( (!global_context -> redo)){
					print_in_box(85,0,0,"   %c[31mThere are odd number of reads in the paired-end data.", CHAR_ESC);
					print_in_box(80,0,0,"   Please make sure that the format is correct.");
				}
				if(line_length > 0){
					int xx1, enters = 0;
					for(xx1 = line_length; xx1 >=0; xx1--){
						if( preload_line[xx1]=='\n' ) enters ++;
						if(2 == enters){
							line_length = xx1+1;
							break;
						}
					}
					if(xx1 <= 0) line_length = 0;
				}
			}

			//printf("FRR=%d\n%s\n", fresh_read_no, preload_line);

			if(line_length > 0)
			{
				while(1)
				{
					int is_finished = 0;
					fc_thread_thread_context_t * thread_context = global_context->thread_contexts+current_thread_id;
					//printf("WRT_THR_IBUF_REM [%d]=%d\n", current_thread_id , thread_context->input_buffer_remainder);

					pthread_spin_lock(&thread_context->input_buffer_lock);
					unsigned int empty_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_remainder; 
					if(empty_bytes > line_length)
					{
						unsigned int tail_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_write_ptr; 
						unsigned int write_p1_len = (tail_bytes > line_length)?line_length:tail_bytes;
						unsigned int write_p2_len = (tail_bytes > line_length)?0:(line_length - tail_bytes);
						memcpy(thread_context->input_buffer + thread_context->input_buffer_write_ptr, preload_line, write_p1_len);
						if(write_p2_len)
						{
							memcpy(thread_context->input_buffer, preload_line + write_p1_len, write_p2_len);
							thread_context->input_buffer_write_ptr = write_p2_len;
						}
						else	thread_context->input_buffer_write_ptr += write_p1_len;
						if(thread_context->input_buffer_write_ptr == global_context->input_buffer_max_size) 
							thread_context->input_buffer_write_ptr=0;


						thread_context->input_buffer_remainder += line_length;
						//printf("WRT_THR_IBUF_REM [%d] + %d =%d\n", current_thread_id, line_length , thread_context->input_buffer_remainder);
						is_finished = 1;
					}

					pthread_spin_unlock(&thread_context->input_buffer_lock);

					current_thread_id++;
					if(current_thread_id >= global_context->thread_number) current_thread_id = 0;

					if(is_finished) break;
					else usleep(tick_time);
				}
			}
			if(!ret) break;
		}
	}

	free(chunk_in_buff);
	free(binary_in_buff);
	free(preload_line);
	global_context->is_all_finished = 1;

	if(global_context->thread_number > 1 || !isSAM)
		fc_thread_wait_threads(global_context);

	unsigned long long int nreads_mapped_to_exon = 0;


	if(!global_context->redo)
		fc_thread_merge_results(global_context, column_numbers , &nreads_mapped_to_exon, my_read_counter, junction_global_table, splicing_global_table);

	fc_thread_destroy_thread_context(global_context);

	//global_context .read_counters.assigned_reads = nreads_mapped_to_exon;

	#ifdef MAKE_STANDALONE
	if(strcmp("STDIN",global_context->input_file_name)!=0)
	#endif
		fclose(fp_in);

	if(sb_header_tab) free(sb_header_tab);
	if(strcmp(global_context->input_file_name,"STDIN")!=0 && isInputFileResortNeeded)
		unlink(global_context->input_file_name);
	free(line);
	return 0;
}


#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int feature_count_main(int argc, char ** argv)
#endif
{
	char * Rargv[38];
	char annot_name[300];
	char * out_name = malloc(300);
	char * fasta_contigs_name = malloc(300);
	char * alias_file_name = malloc(300);
	int cmd_rebuilt_size = 200;
	char * cmd_rebuilt = malloc(cmd_rebuilt_size);
	char nameFeatureTypeColumn[66];
	char nameGeneIDColumn[66];
	int min_qual_score = 0;
	int min_dist = 50;
	int max_dist = 600;
	char debug_command[10];
	char min_dist_str[11];
	char max_dist_str[11];
	char min_qual_score_str[11];
	char feature_block_size_str[11];
	char Strand_Sensitive_Str[11];
	char Pair_Orientations[3];
	char * very_long_file_names;
	int is_Input_Need_Reorder = 0;
	int is_PE = 0;
	int is_SAM = 1;
	int is_GeneLevel = 1;
	int is_Overlap = 0;
	int is_Both_End_Mapped = 0;
	int feature_block_size = 14;
	int Strand_Sensitive_Mode = 0;
	int is_ReadSummary_Report = 0;
	int is_Chimeric_Disallowed = 0;
	int is_PE_Dist_Checked = 0;
	int is_Multi_Mapping_Allowed = 0;
	int is_Split_Alignment_Only = 0;
	int is_duplicate_ignored = 0;
	int do_not_sort = 0;
	int do_junction_cnt = 0;
	int reduce_5_3_ends_to_one = 0;
	int use_fraction_multimapping = 0;
	int threads = 1;
	int isGTF = 1;
	int use_overlapping_length_break_tie = 0;
	char nthread_str[4];
	int option_index = 0;
	int c;
	int very_long_file_names_size = 200;
	int fiveEndExtension = 0, threeEndExtension = 0, minFragmentOverlap = 1;
	char strFiveEndExtension[11], strThreeEndExtension[11], strMinFragmentOverlap[11];
	very_long_file_names = malloc(very_long_file_names_size);
	very_long_file_names [0] = 0;
	fasta_contigs_name[0]=0;

	alias_file_name[0]=0;
	debug_command[0] = 0;

	strcpy(nameFeatureTypeColumn,"exon");
	strcpy(nameGeneIDColumn,"gene_id");
	annot_name[0]=0;out_name[0]=0;
	

	cmd_rebuilt[0]=0;
	for(c = 0; c<argc;c++)
	{
		if(strlen(cmd_rebuilt) + 300 > cmd_rebuilt_size)
		{
			cmd_rebuilt_size*=2;
			cmd_rebuilt = realloc(cmd_rebuilt, cmd_rebuilt_size);
		}
		sprintf(cmd_rebuilt+strlen(cmd_rebuilt), "\"%s\" ", argv[c]);
	}

	optind=0;
	opterr=1;
	optopt=63;

	while ((c = getopt_long (argc, argv, "A:g:t:T:o:a:d:D:L:Q:pbF:fs:S:CBJPMORv?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'S':
				if(strlen(optarg)!=2 || (strcmp(optarg, "ff")!=0 && strcmp(optarg, "rf")!=0 && strcmp(optarg, "fr")!=0)){
					SUBREADprintf("The order parameter can only be ff, fr or rf.\n");
					print_usage();
					return -1;
				}
				Pair_Orientations[0]=(optarg[0]=='r'?'r':'f');
				Pair_Orientations[1]=(optarg[1]=='f'?'f':'r');
				Pair_Orientations[2]=0;

				break;
			case 'J':
				do_junction_cnt = 1;
				break;
			case 'A':
				strcpy(alias_file_name, optarg);
				break;
			case 'M':
				if(0 == is_Multi_Mapping_Allowed)
					is_Multi_Mapping_Allowed = ALLOW_ALL_MULTI_MAPPING;
				break;
			case 'v':
				core_version_number("featureCounts");
				return 0;
			case 'Q':
				min_qual_score = atoi(optarg);
				break;
			case 't':
				strcpy(nameFeatureTypeColumn, optarg);
				break;
			case 'g':
				while((*optarg) == ' ') optarg++;
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
				SUBREADprintf("The '-b' option has been deprecated.\n FeatureCounts will automatically examine the file format.\n");
				is_SAM = 0;
				break;
			case 'C':
				is_Chimeric_Disallowed = 1;
				break;
			case 'P':
				is_PE_Dist_Checked = 1;
				break;
			case 'B':
				is_Both_End_Mapped = 1;
				break;
			case 'f':
				is_GeneLevel = 0;
				break;
			case 'F':
				isGTF = 1;
				if(strcmp("SAF", optarg)==0) isGTF=0;
				else if(strcmp("GTF", optarg)==0) isGTF=1;
				else SUBREADprintf("\nWarning: Unknown annotation format: %s. GTF format is used.\n\n", optarg); 
				break;
			case 'O':
				is_Overlap = 1;
				break;
			case 'R':
				is_ReadSummary_Report = 1;
				break;
			case 's':
				Strand_Sensitive_Mode = atoi(optarg);
				break;
//			case 'i':
//				term_strncpy(sam_name, optarg,299);
//				break;
			case 'o':
				term_strncpy(out_name, optarg,299);
				break;
			case 'a':
				term_strncpy(annot_name, optarg,299);
				break;
			case 'L':
				feature_block_size = atoi(optarg);
				break;
			case 0 :	// long options

				if(strcmp("primary", long_options[option_index].name)==0)
				{
					is_Multi_Mapping_Allowed = ALLOW_PRIMARY_MAPPING;
				}

				if(strcmp("readExtension5", long_options[option_index].name)==0)
				{
					fiveEndExtension = atoi(optarg);
					fiveEndExtension = max(0, fiveEndExtension);
				}

				if(strcmp("readExtension3", long_options[option_index].name)==0)
				{
					threeEndExtension = atoi(optarg);
					threeEndExtension = max(0, threeEndExtension);
				}

				if(strcmp("minOverlap", long_options[option_index].name)==0)
				{
					minFragmentOverlap = atoi(optarg);
				}

				if(strcmp("debugCommand", long_options[option_index].name)==0)
				{
					strcpy(debug_command, optarg);
				}


				if(strcmp("ignoreDup", long_options[option_index].name)==0)
				{
					is_duplicate_ignored = 1 ;
				}

				if(strcmp("fraction", long_options[option_index].name)==0)
				{
					use_fraction_multimapping = 1;
				}

				if(strcmp("fasta", long_options[option_index].name)==0)
				{
					strcpy(fasta_contigs_name , optarg);
				}
				if(strcmp("read2pos", long_options[option_index].name)==0)
				{
					if(optarg[0]=='3')
						reduce_5_3_ends_to_one = REDUCE_TO_3_PRIME_END;
					else if(optarg[0]=='5')
						reduce_5_3_ends_to_one = REDUCE_TO_5_PRIME_END;
						
				}				

				if(strcmp("largestOverlap", long_options[option_index].name)==0)
				{
					use_overlapping_length_break_tie = 1;
				}

				if(strcmp("donotsort", long_options[option_index].name)==0)
				{
					do_not_sort = 1;
				}

				if(strcmp("countSplitAlignmentsOnly", long_options[option_index].name)==0)
				{
					is_Split_Alignment_Only = 1;
				}

				break;
			case '?':
			default :
				print_usage();
				return -1;
				break;
		}


	if(minFragmentOverlap<1)
	{
		fiveEndExtension = - minFragmentOverlap + 1;
		threeEndExtension =  - minFragmentOverlap + 1;
		minFragmentOverlap = 1;
	}

	if(out_name[0]==0 || annot_name[0]==0||argc == optind)
	{
		print_usage();
		return -1;
	}

	for(; optind < argc; optind++)
	{
		int curr_strlen = strlen(very_long_file_names);
		if( very_long_file_names_size - curr_strlen <300)
		{
			very_long_file_names_size *=2;
			//printf("CL=%d ; NS=%d\n", curr_strlen , very_long_file_names_size);
			very_long_file_names=realloc(very_long_file_names , very_long_file_names_size);
		}

		strcat(very_long_file_names, argv[optind]);
		strcat(very_long_file_names, ";");
	}

	very_long_file_names[strlen(very_long_file_names)-1]=0;

	sprintf(strFiveEndExtension, "%d", fiveEndExtension);
	sprintf(strThreeEndExtension, "%d", threeEndExtension);
	sprintf(strMinFragmentOverlap, "%d", minFragmentOverlap);
	sprintf(nthread_str,"%d", threads);
	sprintf(min_dist_str,"%d",min_dist);
	sprintf(max_dist_str,"%d",max_dist);
	sprintf(min_qual_score_str,"%d", min_qual_score);
	sprintf(feature_block_size_str,"%d", feature_block_size);
	sprintf(Strand_Sensitive_Str,"%d", Strand_Sensitive_Mode);
	Rargv[0] = "CreadSummary";
	Rargv[1] = annot_name;
	Rargv[2] = very_long_file_names;
	Rargv[3] = out_name;
	Rargv[4] = is_PE?"1":"0";
	Rargv[5] = min_dist_str;
	Rargv[6] = max_dist_str;
	Rargv[7] = is_SAM?"1":"0";
	Rargv[8] = is_Overlap?"1":"0";
	Rargv[9] = is_GeneLevel?"1":"0";
	Rargv[10] = nthread_str;
	Rargv[11] = isGTF?"1":"0";
	Rargv[12] = Strand_Sensitive_Str;
	Rargv[13] = is_ReadSummary_Report?"1":"0";
	Rargv[14] = is_Both_End_Mapped?"1":"0";
	Rargv[15] = is_Chimeric_Disallowed?"1":"0";
	Rargv[16] = is_PE_Dist_Checked?"1":"0";
	Rargv[17] = nameFeatureTypeColumn;
	Rargv[18] = nameGeneIDColumn;
	Rargv[19] = min_qual_score_str;
	Rargv[20] = is_Multi_Mapping_Allowed == ALLOW_PRIMARY_MAPPING?"2":(is_Multi_Mapping_Allowed == ALLOW_ALL_MULTI_MAPPING?"1":"0");
	Rargv[21] = alias_file_name;
	Rargv[22] = cmd_rebuilt;
	Rargv[23] = is_Input_Need_Reorder?"1":"0";
	Rargv[24] = feature_block_size_str;
	Rargv[25] = strFiveEndExtension;
	Rargv[26] = strThreeEndExtension;
	Rargv[27] = strMinFragmentOverlap;
	Rargv[28] = is_Split_Alignment_Only?"1":"0";
	Rargv[29] = (reduce_5_3_ends_to_one == 0?"0":(reduce_5_3_ends_to_one==REDUCE_TO_3_PRIME_END?"3":"5"));
	Rargv[30] = debug_command;
	Rargv[31] = is_duplicate_ignored?"1":"0";
	Rargv[32] = do_not_sort?"1":"0";
	Rargv[33] = use_fraction_multimapping?"1":"0";
	Rargv[34] = use_overlapping_length_break_tie?"1":"0";
	Rargv[35] = Pair_Orientations;
	Rargv[36] = do_junction_cnt?"1":"0";
	Rargv[37] = fasta_contigs_name;
	int retvalue = readSummary(38, Rargv);

	free(very_long_file_names);
	free(out_name);
	free(alias_file_name);
	free(fasta_contigs_name);
	free(cmd_rebuilt);

	return retvalue;

}


