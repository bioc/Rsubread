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
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>


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
#include "input-blc.h"
#include "hashtable.h"
#include "seek-zlib.h"
#include "HelperFunctions.h"

/********************************************************************/
/********************************************************************/
/********************************************************************/
//  NEW FUNCTION FOR MULTI-THREADING
/********************************************************************/
/********************************************************************/
/********************************************************************/
#define CHROMOSOME_NAME_LENGTH 256 
#define MAX_FC_READ_LENGTH 10001
#define MAX_HIT_NUMBER (1000*1000*1000)
#define MAX_EXTRA_COLS 15
#define FC_FLIST_SPLITOR "\026"

typedef struct{
	int is_negative_strand;
	ArrayList * exons_in_transcript; // the items in IVT tree obj is deallocated when this ArrayList is destroyed. 
	char * gene_name;
} fc_junction_transcript_t;

typedef struct{
	char * transcript_id; // This should be the memory address in the fc_junction_transcript_t table. Not a dedicated memory space for exons.
	int chro_start;
	int chro_stop;
} fc_junction_exon_in_transcript_t;

typedef struct{
	char chromosome_name [CHROMOSOME_NAME_LENGTH+1];
	char gene_name [FEATURE_NAME_LENGTH+1];
	ArrayList * exon_list; // values of fc_junction_exon_in_transcript_t; 
} fc_junction_genebody_t;


#define MAXIMUM_INSERTION_IN_SECTION 8

typedef struct {
	char * chro;
	unsigned int start_pos;
	unsigned int chromosomal_length;
	short insertions;
	unsigned int insertion_start_pos[ MAXIMUM_INSERTION_IN_SECTION ];
	unsigned short insertion_lengths[ MAXIMUM_INSERTION_IN_SECTION ];
} CIGAR_interval_t;


typedef struct {
	char chromosome_name_left[CHROMOSOME_NAME_LENGTH + 1];
	char chromosome_name_right[CHROMOSOME_NAME_LENGTH + 1];
	unsigned int last_exon_base_left;
	unsigned int first_exon_base_right;
} fc_junction_info_t;

typedef struct {
	srInt_64 feature_name_pos;
	unsigned int start;
	unsigned int end;
	unsigned int sorted_order;

	unsigned short chro_name_pos_delta;
	char is_negative_strand;
	char * extra_columns;
} fc_feature_info_t;

typedef struct {
	srInt_64 assigned_reads;

	srInt_64 unassigned_unmapped;
	srInt_64 unassigned_read_type;
	srInt_64 unassigned_singleton;
	srInt_64 unassigned_mappingquality;
	srInt_64 unassigned_chimericreads;
	srInt_64 unassigned_fragmentlength;
	srInt_64 unassigned_duplicate;
	srInt_64 unassigned_multimapping;
	srInt_64 unassigned_secondary;
	srInt_64 unassigned_junction_condition;
	srInt_64 unassigned_nofeatures;
	srInt_64 unassigned_overlapping_length;
	srInt_64 unassigned_ambiguous;
} fc_read_counters;

typedef srInt_64 read_count_type_t;

typedef struct {
	unsigned short thread_id;
	srInt_64 nreads_mapped_to_exon;
	srInt_64 all_reads;
	unsigned int count_table_size;
	read_count_type_t * count_table;
	unsigned int chunk_read_ptr;
	pthread_t thread_object;

	int hits_number_capacity;
	unsigned int * hits_start_pos1;
	unsigned int * hits_start_pos2;

	unsigned short * hits_length1;
	unsigned short * hits_length2;

	char ** hits_chro1;
	char ** hits_chro2;

	srInt_64 * hits_indices1;
	srInt_64 * hits_indices2;

	unsigned int proc_Starting_Chro_Points_1BASE[65536];
	unsigned short proc_Starting_Read_Points[65536];
	unsigned short proc_Section_Read_Lengths[65536];
	char * proc_ChroNames[65536];
	char proc_Event_After_Section[65536];
	CIGAR_interval_t proc_CIGAR_intervals_R1[65536], proc_CIGAR_intervals_R2[65536];

	char ** scoring_buff_gap_chros;
	unsigned int * scoring_buff_gap_starts;
	unsigned short * scoring_buff_gap_lengths;
	char * read_details_buff;
	char * bam_compressed_buff;
	int read_details_buff_used;

	fc_junction_info_t * memory_supported_junctions1;
	fc_junction_info_t * memory_supported_junctions2;

	unsigned int * scoring_buff_numbers;
	unsigned int * scoring_buff_flags;
	unsigned int * scoring_buff_overlappings;
	srInt_64 * scoring_buff_exon_ids;
	srInt_64 del4_added_reads;

	char * chro_name_buff;
	z_stream bam_file_output_stream;

	HashTable * junction_counting_table;   // key: string chro_name \t last_base_previous_exont \t first_base_next_exon
	HashTable * splicing_point_table;
	HashTable * RG_table;	// rg_name -> [ count_table, sum_fc_read_counters, junction_counting_table,  splicing_point_table]
				// NOTE: some reads have no RG tag. These reads are put into the tables in this object but not in the RG_table -> tables.
	fc_read_counters read_counters;
} fc_thread_thread_context_t;

#define READ_SHIFT_UPSTREAM 10
#define READ_SHIFT_DOWNSTREAM 20
#define READ_SHIFT_LEFT 30
#define READ_SHIFT_RIGHT 40
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
	unsigned int reverse_table_start_index_size;
	//unsigned int * reverse_table_end_index;
} fc_chromosome_index_info;

typedef struct {
	int is_gene_level;
	int is_paired_end_mode_assign;
	int is_paired_end_reads_expected;
	int is_multi_overlap_allowed;
	int restricted_no_multi_overlap;
	char * strand_check_mode;
	int is_strand_checked;
	int is_both_end_required;
	int is_chimertc_disallowed;
	int is_PE_distance_checked;
	int is_multi_mapping_allowed;
	int is_primary_alignment_only;
	int is_SAM_file;
	int is_read_details_out;
	int is_junction_no_chro_shown;
	int is_unpaired_warning_shown;
	int is_cigar_overflow_warning_shown;
	int is_stake_warning_shown;
	int is_read_too_long_to_SAM_BAM_shown;
	int is_split_or_exonic_only;
	int is_duplicate_ignored;
	int is_first_read_reversed;
	int is_second_read_straight;
	int is_verbose;
	int long_read_minimum_length;
	int assign_reads_to_RG;
	int use_stdin_file;
	int is_mixed_PE_SE;
	int disk_is_full;
	int do_not_sort;
	int reduce_5_3_ends_to_one;
	int isCVersion;
	int use_fraction_multi_mapping;
	int do_junction_counting;
	int do_detection_call;
	int this_input_number;

	int need_calculate_overlap_len;
	int need_calculate_fragment_len;

	int min_mapping_quality_score;
	int min_paired_end_distance;
	int max_paired_end_distance;
	int max_M;
	int feature_block_size;
	int read_length;
	int line_length;
	int longest_chro_name;
	int five_end_extension;
	int three_end_extension;
	int read_shift_type;
	int read_shift_size;
	int fragment_minimum_overlapping;
	float fractional_minimum_overlapping; 
	float fractional_minimum_feature_overlapping;
	int only_want_largest_overlap;
	int max_missing_bases_in_read, max_missing_bases_in_feature;

	srInt_64 all_reads;

	unsigned short thread_number;
	fc_thread_thread_context_t * thread_contexts;
	int sambam_chro_table_items;
	int is_input_bad_format;
	int any_reads_are_PE;
	SamBam_Reference_Info * sambam_chro_table;
	pthread_spinlock_t read_details_out_lock;

	SAM_pairer_context_t read_pairer;

	char * debug_command;
	char * unistr_buffer_space;
	srInt_64 max_BAM_header_size;
	srInt_64 unistr_buffer_size;
	srInt_64 unistr_buffer_used;
	HashTable * lineno_2_sortedno_tab;
	int known_cell_barcode_length;
//	HashTable * junction_features_table;
//	HashTable * junction_bucket_table;
	HashTable * junction_transcript_table;
	HashTable * junction_genebody_table;
	HashTable * junction_GenebodyTree_table;
	HashTable * junction_ExonTree_table;
	fasta_contigs_t * fasta_contigs;

	HashTable * gene_name_table;	// gene_name -> gene_number
	HashTable * BAM_chros_to_anno_table;	// name in annotation file -> alias name
	HashTable * GCcontent_table; // gene_name -> "qc_content_frac"

	char * RGnames_set;
	int RGnames_capacity;
	int RGnames_ptr;

	char alias_file_name[MAX_FILE_NAME_LENGTH];
	char input_file_name[MAX_FILE_NAME_LENGTH];
	char * input_file_short_name;
	char raw_input_file_name[MAX_FILE_NAME_LENGTH];
	char output_file_name[MAX_FILE_NAME_LENGTH];
	char output_file_path[MAX_FILE_NAME_LENGTH];
	char temp_file_dir[MAX_FILE_NAME_LENGTH];
	char read_details_path[MAX_FILE_NAME_LENGTH];
	char annotation_file_screen_output[MAX_FILE_NAME_LENGTH];
	unsigned char ** gene_name_array;	// gene_internal_number -> gene_name 
	int input_file_unique;

	char * reported_extra_columns;
	HashTable * exontable_chro_table;	// gene_name -> fc_chromosome_index_info structure (contains chro_number, feature_number, block_start, block_end, etc) 
	int exontable_nchrs;
	int exontable_exons;
	int * exontable_geneid;
	char * exontable_strand;
	char ** exontable_chr;
	srInt_64 * exontable_start;
	srInt_64 * exontable_stop;
	char feature_name_column[2000];
	char gene_id_column[100];
	char transcript_id_column[100];

	srInt_64 * exontable_block_end_index;
	srInt_64 * exontable_block_max_end;
	srInt_64 * exontable_block_min_start;

	char ** exontable_anno_chrs;
	char * exontable_anno_chr_2ch;
	srInt_64 * exontable_anno_chr_heads;

	FILE * read_details_out_FP;
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
	unsigned short boundaries_inclusive_base_on_read[global_context -> max_M];
	unsigned int boundaries_inclusive_base_pos[global_context -> max_M];
	char boundaries_chromosomes[global_context -> max_M][MAX_CHROMOSOME_NAME_LEN];
	char boundaries_extend_to_left_on_read[global_context -> max_M];
	int boundaries = 0;

	int cigar_cursor = 0, nch, ret = 0, read_len = 0, x1, x2;
	unsigned int chro_cursor = pos, tmpi = 0;
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
			} else if(nch == 'N'){
				unsigned int last_exon_last_base = chro_cursor - 1;
				unsigned int next_exon_first_base = chro_cursor + tmpi;
				if(ret <= global_context -> max_M - 1){
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
	char current_fusion_cigar[global_context -> max_M * 15];
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

					ret += fetch_boundaries(current_fusion_char, current_fusion_cigar, start_pos, current_fusion_strand, &has_left, &left_on_read, &left_pos, &has_right, &right_on_read, &right_pos, result_junctions + ret, global_context -> max_M - ret );

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

					if(ret <= global_context -> max_M - 1){
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


srInt_64 unistr_cpy(fc_thread_global_context_t * global_context, char * str, int strl)
{
	srInt_64 ret;
	if(global_context->unistr_buffer_used + strl >= global_context->unistr_buffer_size-1)
	{
		if( global_context->unistr_buffer_size < (1000llu*1000u*1000u*32)) // 32GB
		{
			global_context -> unistr_buffer_size = global_context->unistr_buffer_size /2 *3;
			global_context -> unistr_buffer_space = realloc(global_context -> unistr_buffer_space, global_context->unistr_buffer_size);
		}
		else
		{
			SUBREADprintf("Error: exceed memory limit (32GB) for storing feature names.\n");
			return 0xffffffffu;
		}
	}

	strcpy(global_context -> unistr_buffer_space + global_context->unistr_buffer_used, str);
	ret = global_context->unistr_buffer_used;

	global_context->unistr_buffer_used += strl +1;

	return ret;
}

int print_FC_configuration(fc_thread_global_context_t * global_context, char * annot, char * sam, char * out, int is_sam, int is_GTF, int *n_input_files, int isReadSummaryReport, char * PE_exp, char * PE_ass)
{
	char * tmp_ptr1 = NULL , * next_fn, *sam_used = malloc(strlen(sam)+MAX_FILE_NAME_LENGTH), sam_ntxt[30],bam_ntxt[30], next_ntxt[50];
	int nfiles=1, nBAMfiles = 0, nNonExistFiles = 0, x1;
	char MAC_or_random[13];
	mac_or_rand_str(MAC_or_random);

	SUBreadSprintf(sam_used, strlen(sam)+MAX_FILE_NAME_LENGTH, "%s/featureCounts_test_file_writable-%06d-%s.tmp", global_context -> temp_file_dir, getpid(), MAC_or_random);
	FILE * fp = fopen(sam_used,"w");
	if(fp){
		fclose(fp);
		unlink(sam_used);
	}else{
		SUBREADprintf("\nERROR: temporary directory is not writable: '%s'\n\n", global_context -> temp_file_dir);
		return 1;
	}

	strcpy(sam_used, sam);
	nfiles = 0;
	while(1)
	{
		next_fn = strtok_r(nfiles==0?sam_used:NULL, FC_FLIST_SPLITOR, &tmp_ptr1);
		if(next_fn == NULL || strlen(next_fn)<1) break;
		nfiles++;

		srInt_64 BAM_header_size = -1;
		int file_probe = is_certainly_bam_file(next_fn, NULL, &BAM_header_size);
		if(BAM_header_size>0) global_context -> max_BAM_header_size = max( global_context -> max_BAM_header_size , BAM_header_size + 180000);
		if(file_probe==-1){
			nNonExistFiles++;
			if(global_context -> use_stdin_file){
				SUBREADprintf("\nERROR: no valid SAM or BAM file is received from <STDIN>\n\n");
			}else{
				SUBREADprintf("\nERROR: invalid parameter: '%s'\n\n", next_fn);
			}
			return 1;
		}
		if(file_probe == 1) nBAMfiles++;		
	}

	SUBREADputs("");
	print_subread_logo();
	SUBREADputs("");
	print_in_box(80,1,1,"featureCounts setting");
	print_in_box(80,0,0,"");

	sam_ntxt[0]=0;
	bam_ntxt[0]=0;
	next_ntxt[0]=0;

	if(nNonExistFiles)
		SUBreadSprintf(next_ntxt,30, "%d unknown file%s", nNonExistFiles, nNonExistFiles>1?"s":"");
	if(nBAMfiles)
		SUBreadSprintf(bam_ntxt, 30, "%d BAM file%s  ", nBAMfiles, nBAMfiles>1?"s":"");
	if(nfiles-nNonExistFiles-nBAMfiles)
		SUBreadSprintf(sam_ntxt, 30, "%d SAM file%s  ", nfiles-nNonExistFiles-nBAMfiles , (nfiles-nNonExistFiles-nBAMfiles)>1?"s":"");


	strcpy(sam_used, sam);

	print_in_box(80,0,0,"            Input files : %s%s%s", sam_ntxt, bam_ntxt, next_ntxt);
	print_in_box(80,0,0,"");
	nfiles=0;

	while(1){
		next_fn = strtok_r(nfiles==0?sam_used:NULL, FC_FLIST_SPLITOR, &tmp_ptr1);
		if(next_fn == NULL || strlen(next_fn)<1) break;
		//int is_first_read_PE = 0 , file_probe = is_certainly_bam_file(next_fn, &is_first_read_PE, NULL);
		print_in_box(89,0,0,"                          %c[36m%s%c[0m",CHAR_ESC, global_context -> use_stdin_file?"<STDIN>":get_short_fname(next_fn),CHAR_ESC);
		nfiles++;
	}

	(*n_input_files) = nfiles;
	print_in_box(80,0,0,"");

	if(global_context -> annotation_file_screen_output[0]==0){
		print_in_box(80,0,0,"            Output file : %s", get_short_fname(out));
		print_in_box(80,0,0,"                Summary : %s.summary", get_short_fname(out));
	}

	char * PEassignStr = malloc(nfiles * 6);
	char * PEexpectStr = malloc(nfiles * 6);

	int exp_all_same = 1, ass_all_same = 1;

	SUBreadSprintf(PEexpectStr, nfiles * 6,"%s, ", (PE_exp[0]=='1')?"yes":"no");
	char * Ystr = nfiles>5?"Y":"yes";
	char * Nstr = nfiles>5?"N":"no";
	if(PE_exp[1]){
		for(x1=1; PE_exp[x1]; x1++)
			if(PE_exp[x1]!=PE_exp[0]) exp_all_same=0;

		if(!exp_all_same){
			PEexpectStr[0]=0;
			for(x1=0; PE_exp[x1]; x1++)
				SUBreadSprintf(PEexpectStr+strlen(PEexpectStr), nfiles * 6 - strlen(PEexpectStr), "%s, ", (PE_exp[x1]=='1')?Ystr:Nstr);
		}
	}
	PEexpectStr[strlen(PEexpectStr)-2]=0;

	SUBreadSprintf(PEassignStr, nfiles * 6,"%s, ", (PE_ass[0]=='1')?"yes":"no");
	if(PE_ass[1]){
		for(x1=1; PE_ass[x1]; x1++)
			if(PE_ass[x1]!=PE_ass[0]) ass_all_same=0;

		if(!ass_all_same){
			PEassignStr[0]=0;
			for(x1=0; PE_ass[x1]; x1++)
				SUBreadSprintf(PEassignStr+strlen(PEassignStr), nfiles * 6 - strlen(PEassignStr), "%s, ", (PE_ass[x1]=='1')?Ystr:Nstr);
		}
	}
	PEassignStr[strlen(PEassignStr)-2]=0;

	print_in_box(80,0,0,"             Paired-end : %s",PEexpectStr);
	print_in_box(80,0,0,"       Count read pairs : %s",PEassignStr);
	free(PEassignStr);
	free(PEexpectStr);

	if(global_context -> annotation_file_screen_output[0])
		print_in_box(80,0,0,"             Annotation : %s",global_context -> annotation_file_screen_output);
	else
		print_in_box(80,0,0,"             Annotation : %s (%s)", get_short_fname(annot), is_GTF?"GTF":"SAF");
	print_in_box(80,0,0,"     Dir for temp files : %s", global_context->temp_file_dir);

	if(isReadSummaryReport){
		print_in_box(80,0,0,"     Assignment details : <input_file>.featureCounts%s", isReadSummaryReport == FILE_TYPE_BAM?".bam":(isReadSummaryReport == FILE_TYPE_SAM?".sam":""));
		if(global_context -> read_details_path[0])
			print_in_box(80,0,0,"    Details output path : %s", global_context ->read_details_path);
		else
			print_in_box(80,0,0,"                     (Note that files are saved to the output directory)");
		print_in_box(80,0,0,"");
	}

	if(global_context -> do_junction_counting)
		print_in_box(80,0,0,"      Junction Counting : <output_file>.jcounts");
	#ifdef MAKE_STANDALONE	
	#endif

	if(global_context -> alias_file_name[0])
		print_in_box(80,0,0,"  Chromosome alias file : %s", get_short_fname(global_context -> alias_file_name));

	#ifdef MAKE_STANDALONE	
	print_in_box(80,0,0,"");
	#endif
	print_in_box(80,0,0,"                Threads : %d", global_context->thread_number);
	print_in_box(80,0,0,"                  Level : %s level", global_context->is_gene_level?"meta-feature":"feature");
//	print_in_box(80,0,0,"             Paired-end : %s", global_context->is_paired_end_mode_assign?"yes":"no");
	if(global_context -> do_not_sort && global_context->is_paired_end_mode_assign) {
		print_in_box(80,0,0,"       Sorting PE Reads : never");
		print_in_box(80,0,0,"");
	}

	char * multi_mapping_allow_mode = "not counted";
	if(global_context->is_multi_mapping_allowed)
		multi_mapping_allow_mode = global_context -> use_fraction_multi_mapping?"counted (fractional)": "counted";

	print_in_box(80,0,0,"     Multimapping reads : %s", multi_mapping_allow_mode);

	if(global_context-> is_primary_alignment_only)
		print_in_box(80,0,0,"    Multiple alignments : primary alignment only");

	print_in_box(80,0,0,"Multi-overlapping reads : %s", global_context->is_multi_overlap_allowed?"counted":"not counted");
	if(global_context -> is_split_or_exonic_only)
		print_in_box(80,0,0,"       Split alignments : %s", (1 == global_context -> is_split_or_exonic_only)?"only split alignments":"only exonic alignments");
	print_in_box(80,0,0,"  Min overlapping bases : %d", global_context -> fragment_minimum_overlapping);
	if(global_context -> max_missing_bases_in_read >= 0)
		print_in_box(80,0,0,"      Max missing bases : %d in reads", global_context -> max_missing_bases_in_read);
	if(global_context -> max_missing_bases_in_feature >= 0)
		print_in_box(80,0,0,"      Max missing bases : %d in features", global_context -> max_missing_bases_in_feature);
	if(global_context -> fractional_minimum_overlapping > 0.000001)
		print_in_box(81,0,0,"  Min overlapping frac. : %0.1f%%%% to reads", global_context -> fractional_minimum_overlapping*100);
	if(global_context -> fractional_minimum_feature_overlapping > 0.000001)
		print_in_box(81,0,0,"  Min overlapping frac. : %0.1f%%%% to features", global_context -> fractional_minimum_feature_overlapping*100);
	if(global_context -> read_shift_size >0)
		print_in_box(80,0,0,"             Read shift : %d to %s", global_context -> read_shift_size, global_context -> read_shift_type==READ_SHIFT_UPSTREAM?"upstream":( global_context -> read_shift_type==READ_SHIFT_DOWNSTREAM?"downstream":( global_context -> read_shift_type==READ_SHIFT_LEFT?"left":"right")));
	if(global_context -> five_end_extension || global_context -> three_end_extension)
		print_in_box(80,0,0,"         Read extension : %d on 5' and %d on 3' ends", global_context -> five_end_extension , global_context -> three_end_extension);
	if(global_context -> reduce_5_3_ends_to_one)
		print_in_box(80,0,0,"         Read reduction : to %d' end" , global_context -> reduce_5_3_ends_to_one == REDUCE_TO_5_PRIME_END ?5:3);
	if(global_context -> is_duplicate_ignored)
		print_in_box(80,0,0,"       Duplicated Reads : ignored");
	if(global_context -> long_read_minimum_length < 5000)
		print_in_box(80,0,0,"         Long read mode : yes");
	//print_in_box(80,0,0,"      Read orientations : %c%c", global_context->is_first_read_reversed?'r':'f', global_context->is_second_read_straight?'f':'r' );

	if(global_context->is_paired_end_mode_assign)
	{
		print_in_box(80,0,0,"");
		print_in_box(80,0,0,"         Chimeric reads : %s", global_context->is_chimertc_disallowed?"not counted":"counted");
		print_in_box(80,0,0,"       Both ends mapped : %s", global_context->is_both_end_required?"required":"not required");

		if(global_context->is_PE_distance_checked)
			print_in_box(80,0,0,"        Fragment length : %d - %d", global_context -> min_paired_end_distance, global_context -> max_paired_end_distance);
	}

	print_in_box(80,0,0,"");
	print_in_box(80,2,1,"");
	SUBREADputs("");
	print_in_box(80,1,1,"Running");
	print_in_box(80,0,0,"");
	if( global_context -> max_BAM_header_size > 32 * 1024 * 1024 ){
	}
	if(global_context->BAM_chros_to_anno_table)
		print_in_box(80,0,0,"%ld chromosome name aliases are loaded.", global_context -> BAM_chros_to_anno_table ->numOfElements);

	free(sam_used);
	return 0;
}

void print_FC_results(fc_thread_global_context_t * global_context, char * out)
{
	//print_in_box(89,0,1,"%c[36mAlignment assignment finished.%c[0m", CHAR_ESC, CHAR_ESC);
	print_in_box(80,0,0,"");
	#ifdef MAKE_STANDALONE
	print_in_box(80,0,PRINT_BOX_WRAPPED,"Summary of counting results can be found in file \"%s.summary\"", out);
	print_in_box(80,0,0,"");
	#endif
	print_in_box(80,2,1,"");
	SUBREADputs("");
	return;

	SUBREADputs("");
}

int fc_strcmp(const void * s1, const void * s2)
{
	return strcmp((char*)s1, (char*)s2);
}

void junc_transcript_free(fc_junction_transcript_t * txp){
	free(txp -> gene_name);
	ArrayListDestroy(txp -> exons_in_transcript);
	free(txp);
}


void sort_junc_feature_make_gaps(void *k, void *v, HashTable * tab){
	fc_junction_genebody_t * jg = v;
	unsigned int min_start = 0xffffffffu;
	unsigned int max_stop = 0x0u;
	int xk1;
	ArrayList * exl = jg -> exon_list; 
	for(xk1=0; xk1<exl->numOfElements; xk1++){
		fc_junction_exon_in_transcript_t * one_ex = ArrayListGet(exl,xk1);
		min_start = min(min_start, one_ex->chro_start);
		max_stop = max(max_stop, one_ex->chro_stop);
	}
	if(0)if(strcmp(k,"768215\tchrX")==0){
		fprintf(stderr,"POS %u %u of %s\n", min_start,max_stop,k);
	}
	fc_thread_global_context_t *global_context = tab -> appendix1;
	IVT_IntervalTreeNode * IVT_rootnode = HashTableGet(global_context -> junction_GenebodyTree_table, jg -> chromosome_name);
	IVT_rootnode = IVT_insert(IVT_rootnode, min_start, max_stop, jg); 
	HashTablePutReplaceEx(global_context -> junction_GenebodyTree_table, jg -> chromosome_name, IVT_rootnode,0,0,0); // use old key ptr in table; don't free key and value.
}

// this is called after all the features are loaded.
void sort_junc_feature(fc_thread_global_context_t *global_context){
	global_context -> junction_genebody_table -> appendix1 = global_context;
	HashTableIteration(global_context -> junction_genebody_table, sort_junc_feature_make_gaps);
//	fprintf(stderr,"MERGED GeneBodyTree = %ld\n", global_context -> junction_GenebodyTree_table -> numOfElements);
}

void register_junc_feature(fc_thread_global_context_t *global_context, char * feature_name, char * transcript_id, char * chro, unsigned int start, unsigned int stop){
	if(transcript_id==NULL)return;	// SAF format: don't do anything about junction assignment to genes.
	fc_junction_exon_in_transcript_t * new_item = calloc(sizeof(sizeof(fc_junction_exon_in_transcript_t)),1);
	new_item -> transcript_id = HashTableGetKey(global_context -> junction_transcript_table, transcript_id); 
	if(NULL == new_item -> transcript_id){
		new_item -> transcript_id = strdup(transcript_id);
		fc_junction_transcript_t * new_txp = calloc(sizeof(fc_junction_transcript_t),1);
		new_txp -> exons_in_transcript = ArrayListCreate(10); 
		ArrayListSetDeallocationFunction(new_txp -> exons_in_transcript, free);
		ArrayListPush(new_txp -> exons_in_transcript, new_item);
		new_txp -> gene_name = strdup(feature_name);
		HashTablePut(global_context -> junction_transcript_table, new_item->transcript_id, new_txp);
	}else{
		fc_junction_transcript_t * has_txp = HashTableGet(global_context -> junction_transcript_table, transcript_id);
		ArrayListPush(has_txp -> exons_in_transcript, new_item);
	} 

	new_item -> chro_start = start;
	new_item -> chro_stop = stop;

	IVT_IntervalTreeNode * IVT_rootnode = HashTableGet(global_context -> junction_ExonTree_table, chro);
	if(NULL == IVT_rootnode){
		IVT_rootnode = IVT_insert(NULL, start, stop, new_item); // JCount module uses 1-based coordinates.
		char * new_name = strdup(chro);
		HashTablePut(global_context -> junction_ExonTree_table, new_name, IVT_rootnode);
	}else{
		IVT_rootnode = IVT_insert(IVT_rootnode, start, stop, new_item); // JCount module uses 1-based coordinates.
		HashTablePutReplaceEx(global_context -> junction_ExonTree_table, chro, IVT_rootnode,0,0,0); // use old key ptr in table; don't free key and value.
	}

	char gene_body_key[FEATURE_NAME_LENGTH + CHROMOSOME_NAME_LENGTH+10];
	sprintf(gene_body_key,"%s\t%s", feature_name, chro);
	fc_junction_genebody_t * jgbody = HashTableGet(global_context -> junction_genebody_table, gene_body_key);
	if(!jgbody){
		jgbody = malloc(sizeof(fc_junction_genebody_t));
		memset(jgbody,0,sizeof(fc_junction_genebody_t));
		strcpy(jgbody -> gene_name, feature_name);
		strcpy(jgbody -> chromosome_name, chro);
		jgbody -> exon_list = ArrayListCreate(10); 
		// the items are destroyed in  global_context -> junction_transcript_table[txid] => exons_in_transcript
		HashTablePut(global_context -> junction_genebody_table, strdup(gene_body_key), jgbody);
	}

	ArrayListPush(jgbody -> exon_list, new_item); 
}

int match_feature_name_column(char * infile, char * needed){
	char * ptt = NULL;
	char lneeded[strlen(needed)+1];
	strcpy(lneeded, needed);
	char * t1 = strtok_r(lneeded, ",", &ptt);
	while(t1){
		if(strcmp(t1, infile)==0) return 1;
		t1 = strtok_r(NULL,",", &ptt);
	}
	return 0;
}

void junction_genebody_table_free(void * jcv){
	fc_junction_genebody_t * jc = jcv;
	ArrayListDestroy(jc -> exon_list);
 	free(jc);
}

// This function loads annotations from the file.
// It returns the number of featres loaded, or -1 if something is wrong. 
// Memory will be allowcated in this function. The pointer is saved in *loaded_features.
// The invoker must release the memory itself.

#define MAX_ANNOT_LINE_LENGTH 1500000
int load_feature_info(fc_thread_global_context_t *global_context, const char * annotation_file, int file_type, fc_feature_info_t ** loaded_features)
{
	unsigned int features = 0, xk1 = 0, lineno=0;
	char * file_line = malloc(MAX_ANNOT_LINE_LENGTH+1);
	autozip_fp anno_fp;
	int apret = autozip_open(annotation_file, &anno_fp);
	int is_GFF_warned = 0;
	if(apret < 0) return -1;

	HashTable * chro_name_table = HashTableCreate(1603);
	HashTableSetHashFunction(chro_name_table, fc_chro_hash);
	HashTableSetKeyComparisonFunction(chro_name_table, fc_strcmp_chro);
	global_context -> longest_chro_name = 0;

	if(global_context -> do_junction_counting){
		global_context -> junction_genebody_table = StringTableCreate(1603);
		HashTableSetDeallocationFunctions(global_context -> junction_genebody_table, free, (void (*)(void *))junction_genebody_table_free); 

		global_context -> junction_transcript_table = StringTableCreate(1603);
		HashTableSetDeallocationFunctions(global_context -> junction_transcript_table, free, (void (*)(void *))junc_transcript_free); 

		global_context -> junction_ExonTree_table = StringTableCreate(1603);
		HashTableSetDeallocationFunctions(global_context -> junction_ExonTree_table, free, (void (*)(void *))IVT_freeTree);

		global_context -> junction_GenebodyTree_table = StringTableCreate(1603);
		HashTableSetDeallocationFunctions(global_context -> junction_GenebodyTree_table, NULL, (void (*)(void *))IVT_freeTree);
	}

	// first scan: get the chromosome size (that have exons), total number of features
	// also create chro_name_table : chro_name => fc_chromosome_index_info 
	while(0)
	{
		int rchars = autozip_gets(&anno_fp, file_line, MAX_ANNOT_LINE_LENGTH);
		char * token_temp = NULL, *chro_name;
		fc_chromosome_index_info * chro_stab;
		unsigned int feature_pos = 0;
		if(rchars < 1) break;

		lineno++;
		if(is_comment_line(file_line, file_type, lineno-1))continue;
		if(file_type == FILE_TYPE_GTF)
		{
			chro_name = strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp); // lib_name (not needed)
			char * feature_type = strtok_r(NULL,"\t", &token_temp);
			if(match_feature_name_column(feature_type, global_context -> feature_name_column))
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
				chro_stab -> reverse_table_start_index_size = 5000000;
				chro_stab -> reverse_table_start_index = NULL;
				HashTablePut(chro_name_table, tmp_chro_name, chro_stab);
			}

			chro_stab -> chro_features ++;
		}
	}

	//autozip_rewind(&anno_fp);

	unsigned int ret_features_size = 400000;
	fc_feature_info_t * ret_features = malloc(sizeof(fc_feature_info_t) * ret_features_size);
	char * tmpnameex = malloc(50001);

	lineno = 0;
	while(1)
	{
		int is_gene_id_found = 0, is_txn_id_found = 0;
		int rchars = autozip_gets(&anno_fp, file_line, MAX_ANNOT_LINE_LENGTH);
		if(rchars < 1) break;
		if(rchars >= MAX_ANNOT_LINE_LENGTH - 1){
			SUBREADprintf("\nERROR: the %u-th line in your GTF file is extremely long (longer than %d bytes).\nThe program cannot parse this line.\n", lineno+1, MAX_ANNOT_LINE_LENGTH-1);
			return -2;
		}

		lineno++;
		char * token_temp = NULL;
		if(is_comment_line(file_line, file_type, lineno-1))continue;

		if(file_type == FILE_TYPE_RSUBREAD){
			if(xk1 >= ret_features_size) {
				ret_features_size *=2;
				ret_features = realloc(ret_features, sizeof(fc_feature_info_t) * ret_features_size);
			}
			char * feature_name = strtok_r(file_line,"\t",&token_temp);
			int feature_name_len = strlen(feature_name);
			if(feature_name_len > FEATURE_NAME_LENGTH-2){
				SUBREADprintf("WARNING: feature name on the %d-th line is longer than %d bytes. The name is truncated\n", lineno, FEATURE_NAME_LENGTH -2);
				feature_name[FEATURE_NAME_LENGTH -2 ] = 0;
			}

			srInt_64 genename_pos = unistr_cpy(global_context, (char *)feature_name, feature_name_len);
			
	//		SUBREADprintf("REALL: '%s'=%d  [%d] %p  POS=%d\t\tOLD_NAME_POS=%d\n" , feature_name, feature_name_len , xk1, ret_features+xk1, genename_pos, xk1>0?ret_features[xk1-1].feature_name_pos:-1);
			ret_features[xk1].feature_name_pos = genename_pos;

			char * seq_name = strtok_r(NULL,"\t", &token_temp);
			int chro_name_len = strlen(seq_name);
			if(chro_name_len > CHROMOSOME_NAME_LENGTH) seq_name[CHROMOSOME_NAME_LENGTH -1 ] = 0;
			srInt_64 chro_name_pos = unistr_cpy(global_context, (char *)seq_name, chro_name_len);
			global_context -> longest_chro_name = max(chro_name_len, global_context -> longest_chro_name);


			char * start_ptr = strtok_r(NULL,"\t", &token_temp);
			char * end_ptr = strtok_r(NULL,"\t", &token_temp);

			if(start_ptr == NULL || end_ptr == NULL){
				SUBREADprintf("\nWarning: the format on the %d-th line is wrong.\n", lineno);
			}
			srInt_64 tv1 = atoll(start_ptr);
			srInt_64 tv2 = atoll(end_ptr);

			if( isdigit(start_ptr[0]) && isdigit(end_ptr[0]) ){
				if(strlen(start_ptr) > 10 || strlen(end_ptr) > 10 || tv1 > 0x7fffffff || tv2> 0x7fffffff){
					SUBREADprintf("\nError: Line %d contains a coordinate greater than 2^31.\n", lineno);
					return -2;
				}

				if(tv1 >tv2){
					SUBREADprintf("\nError: Line %d contains a feature that do not have a positive length.\n", lineno);
					return -2;
				}
			}else{
				SUBREADprintf("\nError: Line %d contains a format error. The expected annotation format is SAF.\n", lineno);
				return -2;
			}

			ret_features[xk1].chro_name_pos_delta = chro_name_pos - ret_features[xk1].feature_name_pos;
			ret_features[xk1].start = atoi( start_ptr );// start 
			if(ret_features[xk1].start>0x7fffffff)
			{
				ret_features[xk1].start = 0;
				print_in_box(80,0,0,"WARNING the %d-th line has a negative chro coordinate.", lineno);
			}

			ret_features[xk1].end = atoi( end_ptr );//end 
			if(ret_features[xk1].end>0x7fffffff)
			{
				ret_features[xk1].end = 0;
				print_in_box(80,0,0,"WARNING the %d-th line has a negative chro coordinate.", lineno);
			}

			char * strand_str = strtok_r(NULL,"\t", &token_temp); 
			if(strand_str == NULL)
				ret_features[xk1].is_negative_strand = 0;
			else
				ret_features[xk1].is_negative_strand = ('+' ==strand_str[0])?0:(('-' ==strand_str[0])?1:-1);

			if(global_context -> do_detection_call){
				char * GCcontent = strtok_r(NULL,"\t", &token_temp);
				if(GCcontent){
					int gclen = strlen(GCcontent);
					if(gclen>0)GCcontent[gclen-1]=0;
					HashTablePut(global_context -> GCcontent_table, strdup(feature_name) , strdup(GCcontent));
				}
			}

			ret_features[xk1].sorted_order = xk1;

			int bin_location = ret_features[xk1].start / REVERSE_TABLE_BUCKET_LENGTH;
			
			fc_chromosome_index_info * chro_stab = HashTableGet(chro_name_table, seq_name);
			int feature_pos = ret_features[xk1].end;
			if(NULL == chro_stab){
				char * tmp_chro_name = malloc(CHROMOSOME_NAME_LENGTH);
				term_strncpy(tmp_chro_name, seq_name, CHROMOSOME_NAME_LENGTH);
				chro_stab = calloc(sizeof(fc_chromosome_index_info),1);
				chro_stab -> chro_number = chro_name_table->numOfElements;
				chro_stab -> chro_possible_length = feature_pos+1;
				chro_stab -> reverse_table_start_index_size = 5000000;
				chro_stab -> reverse_table_start_index = calloc( chro_stab -> reverse_table_start_index_size  / REVERSE_TABLE_BUCKET_LENGTH +2, sizeof(int));
				HashTablePut(chro_name_table, tmp_chro_name, chro_stab);
			}else chro_stab -> chro_possible_length = max(feature_pos+1, chro_stab -> chro_possible_length);
			chro_stab -> chro_features ++;

			if( chro_stab -> chro_possible_length >= chro_stab -> reverse_table_start_index_size ) {
				unsigned int old_end = sizeof(int) *( chro_stab -> reverse_table_start_index_size  / REVERSE_TABLE_BUCKET_LENGTH +2);
				chro_stab -> reverse_table_start_index_size = max(chro_stab -> reverse_table_start_index_size * 2, (int)(chro_stab -> chro_possible_length * 1.3));
				unsigned int new_size = sizeof(int) *( chro_stab -> reverse_table_start_index_size  / REVERSE_TABLE_BUCKET_LENGTH +2);
				chro_stab -> reverse_table_start_index = realloc( chro_stab -> reverse_table_start_index , new_size);
				memset(chro_stab -> reverse_table_start_index + old_end / sizeof(int), 0, new_size - old_end);
			}
			chro_stab -> reverse_table_start_index[bin_location]++;
			is_gene_id_found = 1;

			assert(feature_name);
			if(global_context -> do_junction_counting){
				register_junc_feature(global_context , feature_name, NULL, seq_name, ret_features[xk1].start, ret_features[xk1].end);
			}

			xk1++;
		} else if(file_type == FILE_TYPE_GTF) {
			char feature_name_tmp[FEATURE_NAME_LENGTH];
			char transcript_id_tmp[FEATURE_NAME_LENGTH];
			SUBreadSprintf(feature_name_tmp, FEATURE_NAME_LENGTH, "LINE_%07u", xk1 + 1);
			transcript_id_tmp[0]=0;
			char * seq_name = strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp);// source
			char * feature_type = strtok_r(NULL,"\t", &token_temp);// feature_type
			if(match_feature_name_column(feature_type, global_context -> feature_name_column)) {
				if(xk1 >= ret_features_size) {
					ret_features_size *=2;
					ret_features = realloc(ret_features, sizeof(fc_feature_info_t) * ret_features_size);
				}

				char * start_ptr = strtok_r(NULL,"\t", &token_temp);
				char * end_ptr = strtok_r(NULL,"\t", &token_temp);

				if(start_ptr == NULL || end_ptr == NULL){
					SUBREADprintf("\nWarning: the format on the %d-th line is wrong.\n", lineno);
				}
				srInt_64 tv1 = atoll(start_ptr);
				srInt_64 tv2 = atoll(end_ptr);

				if( isdigit(start_ptr[0]) && isdigit(end_ptr[0]) ){
					if(strlen(start_ptr) > 10 || strlen(end_ptr) > 10 || tv1 > 0x7fffffff || tv2> 0x7fffffff){
						SUBREADprintf("\nError: Line %d contains a coordinate greater than 2^31.\n", lineno);
						return -2;
					}
				}else{
					SUBREADprintf("\nError: Line %d contains a format error. The expected annotation format is GTF/GFF.\n", lineno);
					return -2;
				}
				ret_features[xk1].start = atoi(start_ptr);// start 
				ret_features[xk1].end = atoi(end_ptr);//end 

				if(ret_features[xk1].start < 1 || ret_features[xk1].end<1 ||  ret_features[xk1].start > 0x7fffffff ||  ret_features[xk1].end > 0x7fffffff || ret_features[xk1].start > ret_features[xk1].end){
					SUBREADprintf("\Error: the feature on the %d-th line has zero coordinate or zero lengths\n\n", lineno);
					return -2;
				}


				strtok_r(NULL,"\t", &token_temp);// score 
				char * strand_str = strtok_r(NULL,"\t", &token_temp);
				ret_features[xk1].is_negative_strand = ('-' == strand_str[0])?1:( ('+' == strand_str[0])?0:-1 );//strand 
				ret_features[xk1].sorted_order = xk1;
				strtok_r(NULL,"\t",&token_temp);	// "frame"
				char * extra_attrs = strtok_r(NULL,"\t",&token_temp);	// name_1 "val1"; name_2 "val2"; ... 
				ret_features[xk1].extra_columns = NULL;
				if(extra_attrs && (strlen(extra_attrs)>2))
				{
					int attr_val_len = GTF_extra_column_value(extra_attrs , global_context -> gene_id_column , feature_name_tmp, FEATURE_NAME_LENGTH);
					int txpid_val_len = GTF_extra_column_value(extra_attrs , global_context -> transcript_id_column , transcript_id_tmp, FEATURE_NAME_LENGTH);
					if(attr_val_len>0) is_gene_id_found=1;
					if(txpid_val_len>0) is_txn_id_found=1;
			//		printf("V=%s\tR=%d\n", extra_attrs , attr_val_len);

					if(global_context -> reported_extra_columns){
						char * extcols = malloc(30);
						int extcols_size = 30, extcols_len = 0;

						char * this_exname_ptr=global_context -> reported_extra_columns;
						while(1){
							int padd0, is_last=1;
							for(padd0=0; this_exname_ptr[padd0]; padd0++)
								if(this_exname_ptr[padd0]=='\t'){
									this_exname_ptr[padd0]=0;
									is_last=0;
									break;
								}

							attr_val_len = GTF_extra_column_value(extra_attrs , this_exname_ptr , tmpnameex, 50000);

							if(attr_val_len<0){
								attr_val_len=2;
								strcpy(tmpnameex,"NA");
							}
							if(attr_val_len + extcols_len + 2 > extcols_size){
								extcols_size = max(extcols_size*2, attr_val_len + extcols_len+2);
								extcols = realloc(extcols, extcols_size);
							}
							memcpy(extcols+extcols_len, tmpnameex, attr_val_len);
							extcols_len += attr_val_len;
							extcols[extcols_len]='\t';
							extcols_len += 1;
							
							if(is_last)break;
							this_exname_ptr[padd0]='\t';
							this_exname_ptr += padd0+1;
						}
						extcols[extcols_len-1]=0;
						ret_features[xk1].extra_columns = extcols;
					}
				}

				if(!is_gene_id_found) {
					if(!is_GFF_warned) {
						int ext_att_len = strlen(extra_attrs);
						if(extra_attrs[ext_att_len-1] == '\n') extra_attrs[ext_att_len-1] =0;
						SUBREADprintf("\nERROR: failed to find the gene identifier attribute in the 9th column of the provided GTF file.\nThe specified gene identifier attribute is '%s' \nAn example of attributes included in your GTF annotation is '%s'.\n\n",  global_context -> gene_id_column, extra_attrs);
					}
					is_GFF_warned++;
				}

				int feature_name_len = strlen(feature_name_tmp);
				if(feature_name_len > FEATURE_NAME_LENGTH-2){
					SUBREADprintf("WARNING: feature name on the %d-th line is longer than %d bytes. The name is truncated\n", lineno, FEATURE_NAME_LENGTH-2);
					feature_name_tmp[FEATURE_NAME_LENGTH -2 ] = 0;
				}
				ret_features[xk1].feature_name_pos = unistr_cpy(global_context, (char *)feature_name_tmp, feature_name_len);

				int chro_name_len = strlen(seq_name);
				if(chro_name_len > CHROMOSOME_NAME_LENGTH) seq_name[CHROMOSOME_NAME_LENGTH -1 ] = 0;
				srInt_64 chro_name_pos = unistr_cpy(global_context, (char *)seq_name, chro_name_len);
				global_context -> longest_chro_name = max(chro_name_len, global_context -> longest_chro_name);

				ret_features[xk1].chro_name_pos_delta = chro_name_pos - ret_features[xk1].feature_name_pos;

				int bin_location = ret_features[xk1].start / REVERSE_TABLE_BUCKET_LENGTH;
				fc_chromosome_index_info * chro_stab = HashTableGet(chro_name_table, seq_name);
				int feature_pos = ret_features[xk1].end;

				if(NULL == chro_stab){
					char * tmp_chro_name = malloc(CHROMOSOME_NAME_LENGTH);
					term_strncpy(tmp_chro_name, seq_name, CHROMOSOME_NAME_LENGTH);
					chro_stab = calloc(sizeof(fc_chromosome_index_info),1);
					chro_stab -> chro_number = chro_name_table->numOfElements;
					chro_stab -> chro_possible_length = feature_pos+1;
					chro_stab -> reverse_table_start_index_size = 5000000;
					chro_stab -> reverse_table_start_index = calloc( chro_stab -> reverse_table_start_index_size  / REVERSE_TABLE_BUCKET_LENGTH +2 , sizeof(int));
					HashTablePut(chro_name_table, tmp_chro_name, chro_stab);
				}else chro_stab -> chro_possible_length = max(feature_pos+1, chro_stab -> chro_possible_length);
				chro_stab -> chro_features ++;
	
				if( chro_stab -> chro_possible_length >= chro_stab -> reverse_table_start_index_size ) {
					int old_end = sizeof(int) *( chro_stab -> reverse_table_start_index_size  / REVERSE_TABLE_BUCKET_LENGTH +2);
					chro_stab -> reverse_table_start_index_size = max(chro_stab -> reverse_table_start_index_size * 2, (int)(chro_stab -> chro_possible_length * 1.3));
					int new_size = sizeof(int) *( chro_stab -> reverse_table_start_index_size  / REVERSE_TABLE_BUCKET_LENGTH +2);
					chro_stab -> reverse_table_start_index = realloc(chro_stab -> reverse_table_start_index, new_size);
					memset(chro_stab -> reverse_table_start_index + old_end / sizeof(int), 0, new_size - old_end);
				}
	
				chro_stab -> reverse_table_start_index[bin_location]++;

				if(global_context -> do_junction_counting){
					if(!is_txn_id_found){
						int ext_att_len = strlen(extra_attrs);
						if(extra_attrs[ext_att_len-1] == '\n') extra_attrs[ext_att_len-1] =0;
						SUBREADprintf("\nERROR: failed to find the transcript identifier attribute in the 9th column of the provided GTF file.\nThe specified gene identifier attribute is '%s' \nAn example of attributes included in your GTF annotation is '%s'.\n\n",  global_context -> transcript_id_column, extra_attrs);
					}
					register_junc_feature(global_context , feature_name_tmp, transcript_id_tmp, seq_name, ret_features[xk1].start, ret_features[xk1].end);
				}

				xk1++;
			}
		}
	}
	features = xk1;
	autozip_close(&anno_fp);
	free(file_line);
	free(tmpnameex);

	(*loaded_features) = ret_features;
	global_context -> exontable_nchrs = (int)chro_name_table-> numOfElements;
	global_context -> exontable_chro_table = chro_name_table;


	if(is_GFF_warned) return -2;
	if(features < 1){
		if(global_context -> annotation_file_screen_output[0]){
			SUBREADprintf("ERROR: no features were loaded in format %s. The annotation format can be specified by the 'isGTFAnnotationFile' option%s.\n", file_type == FILE_TYPE_GTF?"GTF":"SAF", file_type == FILE_TYPE_GTF?", and the required feature type can be specified by the 'GTF.featureType' option":"");
		}else{
			SUBREADprintf("ERROR: no features were loaded in format %s. The annotation format can be specified by the '-F' option%s.\n", file_type == FILE_TYPE_GTF?"GTF":"SAF", file_type == FILE_TYPE_GTF?", and the required feature type can be specified by the '-t' option.":"");
		}
		SUBREADprintf("\n\n");
		return -2;
	}

	print_in_box(80,0,0,"   Features : %d\n", features);
	return features;
}

int find_or_insert_gene_name(fc_thread_global_context_t * global_context, unsigned char * feature_name)
{
	HashTable * genetable = global_context -> gene_name_table;

	srInt_64 gene_number = HashTableGet(genetable, feature_name) - NULL;
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

void register_reverse_table(int block_no, srInt_64 this_block_min_start, srInt_64 this_block_max_end, fc_chromosome_index_info * chro_inf)
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

	srInt_64 * ret_start = (srInt_64 *) arr[0];
	srInt_64 * ret_end = (srInt_64 *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	int total_items = items+items2;
	srInt_64 * tmp_start = malloc(sizeof(srInt_64) * total_items);
	srInt_64 * tmp_end = malloc(sizeof(srInt_64) * total_items);
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

	memcpy(ret_start+ start, tmp_start, sizeof(srInt_64) * total_items);
	memcpy(ret_end+ start, tmp_end, sizeof(srInt_64) * total_items);
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
	srInt_64 * ret_start = (srInt_64 *)arr[0];
	srInt_64 ll = ret_start[l];
	srInt_64 rl = ret_start[r];

	if(ll==rl) return 0;
	else if(ll>rl) return 1;
	else return -1;
}

void feature_sort_exchange(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	srInt_64 tmp;
	fc_feature_info_t * tmpptr;

	srInt_64 * ret_start = (srInt_64 *) arr[0];
	srInt_64 * ret_end = (srInt_64 *) arr[1];
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



void sort_feature_info(fc_thread_global_context_t * global_context, unsigned int features, fc_feature_info_t * loaded_features, char *** sorted_chr_names, int ** sorted_entrezid, srInt_64 ** sorted_start, srInt_64 ** sorted_end, unsigned char ** sorted_strand, char ** anno_chr_2ch, char *** anno_chrs, srInt_64 ** anno_chr_head, srInt_64 ** block_end_index, srInt_64 ** block_min_start_pos, srInt_64 ** block_max_end_pos)
{
	unsigned int chro_pnt;
	unsigned int xk1,xk2;
	int * ret_entrez = malloc(sizeof(int) * features);
	srInt_64 * ret_start = malloc(sizeof(srInt_64) * features);
	srInt_64 * ret_end = malloc(sizeof(srInt_64) * features);
	int current_block_buffer_size = 2000;

	srInt_64 * ret_block_end_index = malloc(sizeof(srInt_64) * current_block_buffer_size);
	srInt_64 * ret_block_min_start = malloc(sizeof(srInt_64) * current_block_buffer_size);
	srInt_64 * ret_block_max_end = malloc(sizeof(srInt_64) * current_block_buffer_size);
	unsigned char * ret_strand = malloc(features);
	char ** ret_char_name = malloc(sizeof(void *) * features);
	fc_feature_info_t ** old_info_ptr = malloc(sizeof(void *) * features);
	(*anno_chrs) = malloc(sizeof(void *) * global_context -> exontable_nchrs);
	(*anno_chr_head) = malloc(sizeof(srInt_64) * global_context -> exontable_nchrs);
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
		srInt_64 this_block_min_start = 0x7fffffff, this_block_max_end = 0;
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
					ret_block_min_start = realloc(ret_block_min_start, sizeof(srInt_64)*current_block_buffer_size);
					ret_block_max_end = realloc(ret_block_max_end, sizeof(srInt_64)*current_block_buffer_size);
					ret_block_end_index = realloc(ret_block_end_index, sizeof(srInt_64)*current_block_buffer_size);
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
				ret_block_min_start = realloc(ret_block_min_start, sizeof(srInt_64)*current_block_buffer_size);
				ret_block_max_end = realloc(ret_block_max_end, sizeof(srInt_64)*current_block_buffer_size);
				ret_block_end_index = realloc(ret_block_end_index, sizeof(srInt_64)*current_block_buffer_size);
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

	if(global_context -> do_junction_counting) sort_junc_feature(global_context);
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

unsigned int calculate_multi_overlap_fraction(fc_thread_global_context_t * global_context, unsigned int fixed_fractional_count, int passed_filters_count){
	if(global_context -> use_fraction_multi_mapping) return fixed_fractional_count / passed_filters_count;
	else return fixed_fractional_count;
}

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


void disallocate_RG_tables(void * pt){
	void ** t4 = pt;
	free(t4[0]);
	free(t4[1]);
	if(t4[2]){
		HashTableDestroy(t4[2]);
		HashTableDestroy(t4[3]);
	}
	free(pt);
}

void fc_thread_start_open_report_fp(fc_thread_global_context_t * global_context);

void process_pairer_reset(void * pairer_vp){
//	SUBREADprintf("RESET COUNTERS\n");
	SAM_pairer_context_t * pairer = (SAM_pairer_context_t *) pairer_vp;
	fc_thread_global_context_t * global_context = (fc_thread_global_context_t * )pairer -> appendix1;
	if(global_context -> sambam_chro_table) free(global_context -> sambam_chro_table);
	global_context -> sambam_chro_table = NULL;
	global_context -> sambam_chro_table_items = 0;
	if(global_context -> assign_reads_to_RG) free(global_context -> RGnames_set);

	int xk1, xk2;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		for(xk2=0; xk2<global_context -> exontable_exons; xk2++)
		{
			global_context -> thread_contexts[xk1].count_table[xk2] = 0;
		}

		global_context -> thread_contexts[xk1].del4_added_reads = 0;

		global_context -> thread_contexts[xk1].all_reads = 0;
		global_context -> thread_contexts[xk1].nreads_mapped_to_exon = 0;

		global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_singleton = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_read_type = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_secondary = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_junction_condition = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_overlapping_length = 0;
		global_context -> thread_contexts[xk1].read_counters.assigned_reads = 0;
		global_context -> thread_contexts[xk1].read_details_buff_used = 0;

		if(global_context -> do_junction_counting)
		{
			HashTableDestroy(global_context -> thread_contexts[xk1].junction_counting_table);
			global_context -> thread_contexts[xk1].junction_counting_table = HashTableCreate(131317);
			HashTableSetHashFunction(global_context -> thread_contexts[xk1].junction_counting_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(global_context -> thread_contexts[xk1].junction_counting_table, free, NULL);
			HashTableSetKeyComparisonFunction(global_context -> thread_contexts[xk1].junction_counting_table, fc_strcmp_chro);
			
			HashTableDestroy(global_context -> thread_contexts[xk1].splicing_point_table);
			global_context -> thread_contexts[xk1].splicing_point_table = HashTableCreate(131317);
			HashTableSetHashFunction(global_context -> thread_contexts[xk1].splicing_point_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(global_context -> thread_contexts[xk1].splicing_point_table, free, NULL);
			HashTableSetKeyComparisonFunction(global_context -> thread_contexts[xk1].splicing_point_table, fc_strcmp_chro);
		}
		
		if(global_context -> assign_reads_to_RG){
			HashTableDestroy(global_context -> thread_contexts[xk1].RG_table);
			global_context -> thread_contexts[xk1].RG_table = HashTableCreate(97);
			HashTableSetHashFunction(global_context -> thread_contexts[xk1].RG_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(global_context -> thread_contexts[xk1].RG_table, free, disallocate_RG_tables);
			HashTableSetKeyComparisonFunction(global_context -> thread_contexts[xk1].RG_table, fc_strcmp_chro);
		}


	}

	if(global_context -> read_details_out_FP) fc_thread_start_open_report_fp(global_context);
}

int is_value_contig_name(char * n, int l){
	int x;
	for(x=0; x<l; x++){
		if(n[x]==0)continue;
		if(n[x]>'~' || n[x]<'!') return 0;
	}
	return 1;
}


#ifdef __MINGW32__
#define this_memmem windows_memmem
#else
#define this_memmem memmem
#endif


void ** get_RG_tables(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * rg_name);
int compress_read_detail_BAM(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, int write_start, int write_end, char * bam_buf);
int process_pairer_header (void * pairer_vp, int thread_no, int is_text, unsigned int items, char * bin, unsigned int bin_len){
	SAM_pairer_context_t * pairer = (SAM_pairer_context_t *) pairer_vp;
	fc_thread_global_context_t * global_context = (fc_thread_global_context_t * )pairer -> appendix1;
	fc_thread_thread_context_t * thread_context = global_context -> thread_contexts;

	//SUBREADprintf("ENTER PROCESS (THRD %d): IS_TXT=%d,  ITEMS = %d, CURRENT_ITEMS=%d\n", thread_no, is_text, items, global_context -> sambam_chro_table_items);

	if(global_context -> is_read_details_out == FILE_TYPE_BAM){
		int write_cursor;
		int first_block = 1;
		for(write_cursor = 0; write_cursor < bin_len; write_cursor += 55000){
			int wlen = min(55000, bin_len - write_cursor);
		
			if( first_block ){
				if(is_text)memcpy(thread_context -> read_details_buff, "BAM\1", 4);
				memcpy(thread_context -> read_details_buff + (is_text?4:0), is_text?(&bin_len):(&items), 4);
			}

			memcpy(thread_context -> read_details_buff + (first_block?4*(1+is_text):0), bin + write_cursor, wlen);
			int blen = compress_read_detail_BAM(global_context, thread_context, 0, wlen + (first_block?4*(1+is_text):0), thread_context -> bam_compressed_buff);
			fwrite( thread_context -> bam_compressed_buff, 1, blen, global_context -> read_details_out_FP);
			first_block = 0;
		} 
	}else if( global_context -> is_read_details_out == FILE_TYPE_SAM && is_text ){
		fwrite( bin, 1, bin_len, global_context -> read_details_out_FP);
	}
	if(is_text ){
		if( global_context -> assign_reads_to_RG ){
			global_context->RGnames_capacity = 10000;
			global_context->RGnames_ptr = 0;
			global_context->RGnames_set =  malloc( global_context->RGnames_capacity );

			int rcursor=0;
			for(;rcursor<bin_len; rcursor++){
				assert(bin[rcursor] == '@'&& bin[rcursor+3] == '\t');
				if(bin[rcursor+1]=='R' && bin[rcursor+2]=='G'){
					int id_start = -1, id_end = -1;
					for(; rcursor < bin_len; rcursor++){
						if(bin[rcursor]=='I' && bin[rcursor+1]=='D'){
							id_start = rcursor + 3;
							id_end = 0;
						}
						for(; rcursor < bin_len; rcursor++){
							if(bin[rcursor]=='\t' || bin[rcursor]=='\n'){
								if(id_end < 1)id_end = rcursor;
								break;
							}
						}
						if(bin[rcursor]=='\n') break;
					}
	
					if(id_start > 0){
						int id_len = id_end - id_start;
						if(global_context->RGnames_capacity < global_context->RGnames_ptr + id_len + 3){
							global_context->RGnames_capacity = global_context->RGnames_capacity * 17 / 10;
							global_context->RGnames_set = realloc( global_context->RGnames_set , global_context->RGnames_capacity );
						}
						memcpy(global_context->RGnames_set + global_context->RGnames_ptr, bin + id_start, id_len);
						global_context->RGnames_set[global_context->RGnames_ptr+id_len]='\t';
						global_context->RGnames_ptr += id_len+1;
					}
				}
				for( ;rcursor<bin_len; rcursor++ ) if(bin[rcursor] == '\n')break;
			}
			if(global_context->RGnames_ptr>0){
				global_context->RGnames_set[global_context->RGnames_ptr-1]=0;
				global_context->RGnames_ptr--;
			}
			//SUBREADprintf("RGList: %s\n", global_context->RGnames_set);

			int thread_no;
			for(thread_no = 0; thread_no < global_context -> thread_number; thread_no ++){
				fc_thread_thread_context_t * RGthread_context = global_context -> thread_contexts + thread_no;
				int RGcursor = 0;
				char *lastRGptr = global_context->RGnames_set;
				for(; RGcursor < global_context->RGnames_ptr+1; RGcursor++){
					if(global_context->RGnames_set[ RGcursor ] == '\t' || global_context->RGnames_set[ RGcursor ] == 0){
						global_context->RGnames_set[ RGcursor ] = 0;
						if(strlen(lastRGptr)>0){
					//		SUBREADprintf("PUT 4Tab:'%s'\n", lastRGptr);
							get_RG_tables(global_context, RGthread_context, lastRGptr);
							lastRGptr = global_context->RGnames_set + RGcursor +1;
							if(RGcursor < global_context->RGnames_ptr)
								global_context->RGnames_set[ RGcursor ] = '\t';
						}
					}
				}
			}
		}
	}else{
		if(global_context -> sambam_chro_table)
			global_context -> sambam_chro_table = delay_realloc(global_context -> sambam_chro_table, global_context -> sambam_chro_table_items * sizeof(SamBam_Reference_Info), (items + global_context -> sambam_chro_table_items) * sizeof(SamBam_Reference_Info));
		else global_context -> sambam_chro_table = malloc(items * sizeof(SamBam_Reference_Info));

		int x1, bin_ptr = 0;
		for(x1 =  global_context -> sambam_chro_table_items; x1 <  global_context -> sambam_chro_table_items+items; x1++){
			int l_name;
			memcpy(&l_name, bin + bin_ptr, 4);
			bin_ptr += 4;

			if( !is_value_contig_name(bin + bin_ptr, l_name)){
				SUBREADprintf("The chromosome name contains unexpected characters: \"%s\" (%d chars)\nfeatureCounts has to stop running\n", bin + bin_ptr, l_name);
				return -1;
			}
			if(l_name >= MAX_CHROMOSOME_NAME_LEN){
				SUBREADprintf("The chromosome name of \"%s\" contains %d characters, longer than the upper limit of %d\nfeatureCounts has to stop running\n",  bin + bin_ptr , l_name, MAX_CHROMOSOME_NAME_LEN - 1);
				return -1;
			}
			memcpy(global_context -> sambam_chro_table[x1].chro_name ,  bin + bin_ptr, l_name);
			//SUBREADprintf("The %d-th is '%s'\n", x1, global_context -> sambam_chro_table[x1].chro_name);
			bin_ptr += l_name;
			memcpy(&global_context -> sambam_chro_table[x1].chro_length ,  bin + bin_ptr, 4);
			bin_ptr += 4;
		}
		global_context -> sambam_chro_table_items += items;
	}
	return 0;
}

void process_line_buffer(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * bin1, char * bin2);

int reverse_flag(int mf){
	int ret = mf & 3;
	if(mf & 4) ret |= 8;
	if(mf & 8) ret |= 4;
	if((mf & 1)==0) ret |= 4;

	if(mf & 0x10) ret |= 0x20;
	if(mf & 0x20) ret |= 0x10;

	if(mf & 0x40) ret |= 0x80;
	if(mf & 0x80) ret |= 0x40;
	return ret;
}

int calc_total_frag_one_len(CIGAR_interval_t * intvs, int intvn, char * read_name){
	int ret = 0, x1;
	for(x1 = 0; x1 < intvn; x1++){
		int x2;
		//#warning "=========== DEBUG OUT =============="
		if(0 && FIXLENstrcmp("V0112_0155:7:1101:20072:12961#ATCAC", read_name)==0){
			SUBREADprintf("READ %s SINGLE: chro_len = %d, secs = %d\n" , read_name, intvs[x1].chromosomal_length, intvs[x1].insertions);
		}
		for(x2 = 0; x2 < intvs[x1].insertions; x2++) ret += intvs[x1].insertion_lengths[x2];
		ret += intvs[x1].chromosomal_length;
	}
	return ret;
}

int calc_total_has_overlap(unsigned int r1_start, unsigned int r1_end, unsigned int r2_start, unsigned int r2_end, unsigned int * overlap_start, unsigned int * overlap_end){
	if((r1_start <= r2_start && r1_end > r2_start) || (r2_start <= r1_start && r2_end > r1_start) ){
		(*overlap_start) = max( r1_start, r2_start );
		(*overlap_end) = min( r1_end, r2_end );
		return 1;
	}
	return 0;
}

int calc_total_frag_len( fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, CIGAR_interval_t * CIGAR_intervals_R1, int CIGAR_intervals_R1_sections, CIGAR_interval_t * CIGAR_intervals_R2, int CIGAR_intervals_R2_sections, char * read_name){
	if     ( CIGAR_intervals_R1_sections == 0 && CIGAR_intervals_R2_sections > 0) return calc_total_frag_one_len( CIGAR_intervals_R2,CIGAR_intervals_R2_sections , read_name);
	else if( CIGAR_intervals_R1_sections  > 0 && CIGAR_intervals_R2_sections== 0) return calc_total_frag_one_len( CIGAR_intervals_R1,CIGAR_intervals_R1_sections , read_name);
	else if( CIGAR_intervals_R1_sections == 0 && CIGAR_intervals_R2_sections== 0) return 0;

	if(CIGAR_intervals_R1_sections > 0 && CIGAR_intervals_R2_sections > 0 && strcmp(CIGAR_intervals_R1[0].chro, CIGAR_intervals_R2[0].chro )!=0 )
		// two reads are from different chromosomes
		return calc_total_frag_one_len( CIGAR_intervals_R2,CIGAR_intervals_R2_sections , read_name) + calc_total_frag_one_len( CIGAR_intervals_R1,CIGAR_intervals_R1_sections , read_name);

	if(0 && FIXLENstrcmp("V0112_0155:7:1101:20072:12961#ATCAC", read_name)==0){
		int xx;
		for(xx = 0; xx < CIGAR_intervals_R1_sections; xx++)
				SUBREADprintf("R1 SEC %d: %u + %d\n", xx, CIGAR_intervals_R1[xx].start_pos,  CIGAR_intervals_R1[xx].chromosomal_length );
		for(xx = 0; xx < CIGAR_intervals_R2_sections; xx++)
				SUBREADprintf("R2 SEC %d: %u + %d\n", xx, CIGAR_intervals_R2[xx].start_pos,  CIGAR_intervals_R2[xx].chromosomal_length );
	}
	
	unsigned int merged_section_count = 0;
	unsigned short merged_section_lengths[ MAXIMUM_INSERTION_IN_SECTION * 3 ];
	unsigned int merged_section_indel_counts[ MAXIMUM_INSERTION_IN_SECTION * 3 ];
	unsigned short merged_section_indel_lengths[ MAXIMUM_INSERTION_IN_SECTION * 3 ][ MAXIMUM_INSERTION_IN_SECTION ];

	int R1_i = 0 , R2_i = 0;
	while (1){
		//SUBREADprintf("FRAGDEBUG %s : %d < %d & %d < %d; MC=%d; INS1=%d INS2=%d\n", read_name, R1_i,CIGAR_intervals_R1_sections,R2_i,CIGAR_intervals_R2_sections, merged_section_count, CIGAR_intervals_R1[R1_i].insertions, CIGAR_intervals_R2[R2_i].insertions);
		if( R1_i >= CIGAR_intervals_R1_sections &&  R2_i >= CIGAR_intervals_R2_sections ) break;

		if( R1_i < CIGAR_intervals_R1_sections && R2_i < CIGAR_intervals_R2_sections){
			// see if R1 and R2 overlap
			// if not: add R2 to specific sction; R2_i ++
			// elif overlap: add the R1 first_half and/or R2 first_half or zero to specific section, and add overlapping part to overlapping section; DO NOT add the second specific half!
			// 	if R1_end > R2_end: R1_section_start = overlapping_end; R2_i ++
			// 	elif R2_end > R1_end: R2_section_start = overlapping_end; R1_i ++ 
			// 	elif R2_end == R1_end: R1_i++; R2_i++

			unsigned int overlapping_start= 0 ,  overlapping_end = 0;

			int is_r1r2_overlap = 0;

			is_r1r2_overlap = calc_total_has_overlap( CIGAR_intervals_R1[R1_i].start_pos, CIGAR_intervals_R1[R1_i].start_pos + CIGAR_intervals_R1[R1_i].chromosomal_length , CIGAR_intervals_R2[R2_i].start_pos, CIGAR_intervals_R2[R2_i].start_pos + CIGAR_intervals_R2[R2_i].chromosomal_length , & overlapping_start , & overlapping_end);

			if( is_r1r2_overlap ){
				if (CIGAR_intervals_R1[R1_i].start_pos > CIGAR_intervals_R2[R2_i].start_pos ){
					//first half_R2 add special
					merged_section_lengths[merged_section_count] = overlapping_start - CIGAR_intervals_R2[R2_i].start_pos;

					int indel_i;
					for(indel_i = 0; indel_i < min(MAXIMUM_INSERTION_IN_SECTION,CIGAR_intervals_R2[R2_i].insertions); indel_i++){
						if( CIGAR_intervals_R2[R2_i].insertion_start_pos[indel_i] >= overlapping_start ){
							if(indel_i>0){
								int insmov_i, ins_dist_i = 0;
								for(insmov_i = indel_i ; insmov_i < CIGAR_intervals_R2[R2_i].insertions; insmov_i++){
									assert(MAXIMUM_INSERTION_IN_SECTION > ins_dist_i);
									assert(MAXIMUM_INSERTION_IN_SECTION > insmov_i);
									CIGAR_intervals_R2[R2_i].insertion_start_pos[ins_dist_i] = CIGAR_intervals_R2[R2_i].insertion_start_pos[insmov_i];
									CIGAR_intervals_R2[R2_i].insertion_lengths[ins_dist_i] = CIGAR_intervals_R2[R2_i].insertion_lengths[insmov_i];
									ins_dist_i++;
								}
								CIGAR_intervals_R2[R2_i].insertions = ins_dist_i;
							}
							break;
						}
						merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R2[R2_i].insertion_lengths[indel_i];
					}
					merged_section_indel_counts[merged_section_count] = indel_i;

					merged_section_count ++;
	
				}else if( CIGAR_intervals_R1[R1_i].start_pos < CIGAR_intervals_R2[R2_i].start_pos ){
					//first half_R1 add special
					merged_section_lengths[merged_section_count] = overlapping_start - CIGAR_intervals_R1[R1_i].start_pos;

					int indel_i;
					for(indel_i = 0; indel_i < min(MAXIMUM_INSERTION_IN_SECTION,CIGAR_intervals_R1[R1_i].insertions); indel_i++){
						if( CIGAR_intervals_R1[R1_i].insertion_start_pos[indel_i] >= overlapping_start ){
							if(indel_i>0){
								int insmov_i, ins_dist_i = 0;
								for(insmov_i = indel_i ; insmov_i < CIGAR_intervals_R1[R1_i].insertions; insmov_i++){
									assert(MAXIMUM_INSERTION_IN_SECTION > insmov_i);
									CIGAR_intervals_R1[R1_i].insertion_start_pos[ins_dist_i] = CIGAR_intervals_R1[R1_i].insertion_start_pos[insmov_i];
									CIGAR_intervals_R1[R1_i].insertion_lengths[ins_dist_i] = CIGAR_intervals_R1[R1_i].insertion_lengths[insmov_i];
									ins_dist_i++;
								}
								CIGAR_intervals_R1[R1_i].insertions = ins_dist_i;
							}
							break;
						}

						merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R1[R1_i].insertion_lengths[indel_i];
					}
					merged_section_indel_counts[merged_section_count] = indel_i;

					merged_section_count ++;
				}

				merged_section_lengths[merged_section_count] = overlapping_end - overlapping_start;
				merged_section_indel_counts[merged_section_count] = 0;


				int indel_i_R1 = 0, indel_i_R2 = 0;
				while(1){
					//SUBREADprintf("FRAGDEBUG: CC[%d] = %d ; II1=%d < %d; II2=%d < %d\n", merged_section_count,  merged_section_indel_counts[merged_section_count], indel_i_R1, CIGAR_intervals_R1[R1_i].insertions , indel_i_R2, CIGAR_intervals_R2[R2_i].insertions);

					if( indel_i_R1 >= CIGAR_intervals_R1[R1_i].insertions ||  indel_i_R2 >= CIGAR_intervals_R2[R2_i].insertions ){
						if(indel_i_R1 > 0){
							int insmov_i, ins_dist_i = 0;
							for(insmov_i = indel_i_R1 ; insmov_i < min(MAXIMUM_INSERTION_IN_SECTION,CIGAR_intervals_R1[R1_i].insertions); insmov_i++){
								assert(MAXIMUM_INSERTION_IN_SECTION > insmov_i);
								CIGAR_intervals_R1[R1_i].insertion_start_pos[ins_dist_i] = CIGAR_intervals_R1[R1_i].insertion_start_pos[insmov_i];
								CIGAR_intervals_R1[R1_i].insertion_lengths[ins_dist_i] = CIGAR_intervals_R1[R1_i].insertion_lengths[insmov_i];
								ins_dist_i++;
							}
							CIGAR_intervals_R1[R1_i].insertions = ins_dist_i;
						}
						if(indel_i_R2 > 0){
							int insmov_i, ins_dist_i = 0;
							for(insmov_i = indel_i_R2 ; insmov_i < min(CIGAR_intervals_R2[R2_i].insertions,MAXIMUM_INSERTION_IN_SECTION); insmov_i++){
								assert(MAXIMUM_INSERTION_IN_SECTION > insmov_i);
								CIGAR_intervals_R2[R2_i].insertion_start_pos[ins_dist_i] = CIGAR_intervals_R2[R2_i].insertion_start_pos[insmov_i];
								CIGAR_intervals_R2[R2_i].insertion_lengths[ins_dist_i] = CIGAR_intervals_R2[R2_i].insertion_lengths[insmov_i];
								ins_dist_i++;
							}
							CIGAR_intervals_R2[R2_i].insertions = ins_dist_i;
						}
						break;
					}

					if( CIGAR_intervals_R1[R1_i].insertion_start_pos[indel_i_R1] > CIGAR_intervals_R2[R2_i].insertion_start_pos[indel_i_R2] ) indel_i_R2 ++;
					else if( CIGAR_intervals_R1[R1_i].insertion_start_pos[indel_i_R1] < CIGAR_intervals_R2[R2_i].insertion_start_pos[indel_i_R2]  ) indel_i_R1 ++;
					else{
						if( CIGAR_intervals_R1[R1_i].insertion_lengths[ indel_i_R1 ] == CIGAR_intervals_R2[R2_i].insertion_lengths[ indel_i_R2 ] ){
							merged_section_indel_lengths[merged_section_count][ merged_section_indel_counts[merged_section_count] ] = CIGAR_intervals_R1[R1_i].insertion_lengths[indel_i_R1];
							merged_section_indel_counts[merged_section_count] ++;
						}
						indel_i_R2++;
						indel_i_R1++;
					}
				}

				merged_section_count ++;

				// add common

				if(CIGAR_intervals_R1[R1_i].start_pos + CIGAR_intervals_R1[R1_i].chromosomal_length > CIGAR_intervals_R2[R2_i].start_pos + CIGAR_intervals_R2[R2_i].chromosomal_length){
					CIGAR_intervals_R1[R1_i].chromosomal_length -= ( overlapping_end - CIGAR_intervals_R1[R1_i].start_pos );
					CIGAR_intervals_R1[R1_i].start_pos = overlapping_end;
					R2_i ++;
				}else if(CIGAR_intervals_R1[R1_i].start_pos + CIGAR_intervals_R1[R1_i].chromosomal_length < CIGAR_intervals_R2[R2_i].start_pos + CIGAR_intervals_R2[R2_i].chromosomal_length){
					CIGAR_intervals_R2[R2_i].chromosomal_length -= ( overlapping_end - CIGAR_intervals_R2[R2_i].start_pos );
					CIGAR_intervals_R2[R2_i].start_pos = overlapping_end;
					R1_i ++;
				}else{
					R1_i ++;
					R2_i ++;
				}

			}else if(CIGAR_intervals_R1[R1_i].start_pos >  CIGAR_intervals_R2[R2_i].start_pos){
				merged_section_lengths[merged_section_count] = CIGAR_intervals_R2[R2_i].chromosomal_length;

				int indel_i;
				for(indel_i = 0; indel_i < CIGAR_intervals_R2[R2_i].insertions; indel_i++){
					merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R2[R2_i].insertion_lengths[indel_i];
				}
				merged_section_indel_counts[merged_section_count] = CIGAR_intervals_R2[R2_i].insertions;

				merged_section_count ++;
				R2_i ++;
			}else{
				merged_section_lengths[merged_section_count] = CIGAR_intervals_R1[R1_i].chromosomal_length;

				int indel_i;
				for(indel_i = 0; indel_i < CIGAR_intervals_R1[R1_i].insertions; indel_i++){
					merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R1[R1_i].insertion_lengths[indel_i];
				}
				merged_section_indel_counts[merged_section_count] = CIGAR_intervals_R1[R1_i].insertions;

				merged_section_count ++;
				R1_i ++;
			}
		}else if(R1_i < CIGAR_intervals_R1_sections){ 
			// add R1 section to specific section
			// R1_i ++
			merged_section_lengths[merged_section_count] = CIGAR_intervals_R1[R1_i].chromosomal_length;

			int indel_i;
			for(indel_i = 0; indel_i < CIGAR_intervals_R1[R1_i].insertions; indel_i++){
				merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R1[R1_i].insertion_lengths[indel_i];
			}
			merged_section_indel_counts[merged_section_count] = CIGAR_intervals_R1[R1_i].insertions;

			merged_section_count ++;
			R1_i ++;
		}else if(R2_i < CIGAR_intervals_R2_sections){
			merged_section_lengths[merged_section_count] = CIGAR_intervals_R2[R2_i].chromosomal_length;

			int indel_i;
			for(indel_i = 0; indel_i < CIGAR_intervals_R2[R2_i].insertions; indel_i++){
				merged_section_indel_lengths[merged_section_count][indel_i] = CIGAR_intervals_R2[R2_i].insertion_lengths[indel_i];
			}
			merged_section_indel_counts[merged_section_count] = CIGAR_intervals_R2[R2_i].insertions;

			merged_section_count ++;
			R2_i ++;
		}
	}

	int ret = 0, x1, x2;
	for(x1 = 0; x1 < merged_section_count ; x1++){
		ret += merged_section_lengths[x1];
		for(x2 = 0; x2 < merged_section_indel_counts[x1]; x2++)
			ret += merged_section_indel_lengths[x1][x2];
//		SUBREADprintf("FRAGDEBUG %s [%d] : len = %d , indels = %d\n" , read_name, x1, merged_section_lengths[x1] , merged_section_indel_counts[x1]);
	}

	return ret; 
}

void get_readname_from_bin(char * bin, char ** read_name){
	(*read_name) = bin + 36;
}

#define INC_intervals_i_AND_CLEAN {(*intervals_i) ++;memset(intervals_buffer+ (*intervals_i), 0, sizeof(CIGAR_interval_t));}

void parse_bin(SamBam_Reference_Info * sambam_chro_table, char * bin, char * bin2, char ** read_name, int * flag, char ** chro, srInt_64 * pos, int * mapq, char ** mate_chro, srInt_64 * mate_pos, srInt_64 * tlen, int * is_junction_read, int * cigar_sect, int * cigar_overflow, unsigned int * Starting_Chro_Points_1BASE, unsigned short * Starting_Read_Points, unsigned short * Section_Read_Lengths, char ** ChroNames, char * Event_After_Section, int * NH_value, int max_M, CIGAR_interval_t * intervals_buffer, int * intervals_i, int assign_reads_to_RG, char ** RG_ptr, int * ret_me_refID, int * ret_mate_refID){
	int x1, len_of_S1 = 0;
	*cigar_sect = 0;
	*NH_value = 1;
	*flag = 0;
	*is_junction_read = 0;
	assert(bin||bin2);

	if(bin){
		(*read_name) = bin + 36;
		memcpy(flag, bin + 16, 4);
		int cigar_opts = (*flag) & 0xffff;
		(*flag) = (*flag) >> 16;
		int refID, mate_refID;
		memcpy(&refID, bin + 4, 4);
		if(refID >= 0) (*chro) = sambam_chro_table[refID].chro_name;
		else (*chro) = NULL;

		(*pos) = 0;
		memcpy(pos, bin+8, 4);
		(*pos) ++;
		
		memcpy(mapq, bin+12, 4);
		int l_read_name = (*mapq)& 0xff;
		(*mapq) = ((*mapq)>>8)&0xff;

		int seq_len;
		memcpy(&seq_len, bin + 20,4);
		memcpy(&mate_refID, bin+24, 4);
		if(mate_refID>=0) (*mate_chro) = sambam_chro_table[mate_refID].chro_name;
		else	(*mate_chro) = NULL;

		*ret_mate_refID = mate_refID;
		*ret_me_refID = refID;

		(*mate_pos)=0;
		memcpy(mate_pos, bin+28, 4);
		(*mate_pos)++;

		int tlen_int;
		memcpy(&tlen_int, bin+32, 4);
		(*tlen) = tlen_int;

		int * cigar_opt_ints = (int *)(bin + 36 + l_read_name);
		unsigned int chro_cursor = (*pos), section_start_chro = (*pos);
		unsigned short read_cursor = 0, this_section_length = 0, section_start_read = 0;

		if(intervals_buffer){
			intervals_buffer[ *intervals_i ].start_pos = chro_cursor;
			intervals_buffer[ *intervals_i ].chro = *chro;
		}

		for(x1 = 0 ; x1 < cigar_opts; x1++){
			int optype = cigar_opt_ints[x1]&0xf;
			int optval = (cigar_opt_ints[x1]>>4)& 0xfffffff;
			if(optype == 0 || optype == 7 || optype == 8){ // 'M' , '=', 'X'
				chro_cursor += optval;
				read_cursor += optval;
				this_section_length += optval;
			}else if(optype == 1 || optype == 2 || optype == 3){ // 'I', 'D' or 'N'
				if(3 == optype)
					(*is_junction_read) = 1;
				char event_char=0;
				if(optype == 3) event_char = 'N';
				if(optype == 2) event_char = 'D';
				else if(optype == 1){
					if(intervals_buffer && intervals_buffer[ *intervals_i ].insertions < MAXIMUM_INSERTION_IN_SECTION){
						intervals_buffer[ *intervals_i ].insertion_start_pos[  intervals_buffer[ *intervals_i ].insertions  ] = chro_cursor;
						intervals_buffer[ *intervals_i ].insertion_lengths[ intervals_buffer[ *intervals_i ].insertions ] = optval;
						intervals_buffer[ *intervals_i ].insertions ++;
					}
					event_char = 'I';
				}

				if( (*cigar_sect) < max_M){
					Event_After_Section[*cigar_sect] = event_char;
					Starting_Chro_Points_1BASE[*cigar_sect] = section_start_chro; 
					Starting_Read_Points[*cigar_sect] = section_start_read;
					Section_Read_Lengths[*cigar_sect] = this_section_length;
					ChroNames[*cigar_sect] = (*chro);
					(*cigar_sect)++;

					if(intervals_buffer){
						intervals_buffer[ *intervals_i ].chromosomal_length = chro_cursor - intervals_buffer[ *intervals_i ].start_pos;
						INC_intervals_i_AND_CLEAN;
					}
				} else (*cigar_overflow)=1;

				if(optype == 2 || optype == 3)// N or D
					chro_cursor += optval;
				else
					read_cursor += optval;

				if(intervals_buffer && (*cigar_sect) < max_M){
					intervals_buffer[ *intervals_i ].start_pos = chro_cursor;
					intervals_buffer[ *intervals_i ].chro = *chro;
				}

				section_start_chro = chro_cursor;
				section_start_read = read_cursor;
				this_section_length = 0;
			}else if(optype == 4){ // 'S'
				if(read_cursor==0)
				{
					read_cursor += optval;
					section_start_read = read_cursor;

					if(intervals_buffer){
						if(intervals_buffer[ *intervals_i ].start_pos > optval) intervals_buffer[ *intervals_i ].start_pos -= optval;
						else intervals_buffer[ *intervals_i ].start_pos = 0;
					}
				}else	len_of_S1 = optval;
			}	// H and P do not have effect on cigar parsing.
		}
		if(this_section_length>0){
			// add new section
			if( (*cigar_sect) < max_M){
				if(intervals_buffer){
					intervals_buffer[ *intervals_i ].chromosomal_length = chro_cursor - intervals_buffer[ *intervals_i ].start_pos + len_of_S1;
					INC_intervals_i_AND_CLEAN;
					
				}
				Starting_Chro_Points_1BASE[*cigar_sect] = section_start_chro; 
				Starting_Read_Points[*cigar_sect] = section_start_read;
				Section_Read_Lengths[*cigar_sect] = this_section_length ;
				ChroNames[*cigar_sect] = (*chro);
				(*cigar_sect)++;
			} else (*cigar_overflow)=1;
		}
		
		int bin_ptr = 36 + l_read_name + seq_len + (seq_len+1)/2 + 4 * cigar_opts;
		int block_len;
		memcpy(&block_len, bin, 4);
		int found_NH = SAM_pairer_iterate_int_tags((unsigned char *)bin+bin_ptr, block_len + 4 - bin_ptr, "NH", NH_value);
		if(!found_NH) *(NH_value) = 1;

		if(assign_reads_to_RG){
			char RG_type = 0;
			SAM_pairer_iterate_tags((unsigned char *)bin+bin_ptr, block_len + 4 - bin_ptr, "RG", &RG_type, RG_ptr);
			if(RG_type != 'Z') (*RG_ptr) = NULL;
		}
		//SUBREADprintf("FOUND=%d, NH=%d, TAG=%.*s\n", found_NH, *(NH_value), 3 , bin+bin_ptr);
	}else{
		(*read_name) = bin2 + 36;
		int mate_flag;
		memcpy(&mate_flag, bin2 + 16, 4);
		mate_flag = mate_flag >> 16;
		(*flag) = reverse_flag(mate_flag);

		int refID, mate_refID;
		memcpy(&refID, bin2 + 24, 4);
		memcpy(&mate_refID, bin2 + 4, 4);
		if(refID < 0) *chro = NULL; 
		else (*chro) = sambam_chro_table[refID].chro_name;

		if(mate_refID < 0) *mate_chro = NULL;
		else (*mate_chro) = sambam_chro_table[mate_refID].chro_name; 
		*ret_mate_refID = mate_refID;
		*ret_me_refID = refID;

		*pos=0;
		memcpy(pos, bin2+28, 4);
		(*pos)++;

		*mate_pos=0;
		memcpy(mate_pos, bin2+8, 4);
		(*mate_pos)++;
	
		(*tlen) = 0;
		memcpy(tlen, bin2+32, 4);
		(*tlen) = -(*tlen);

		if(assign_reads_to_RG){
			char RG_type = 0;
			int block2_len = 0;
			memcpy(&block2_len, bin2, 4);
			int rname2len = 0, cigar2len = 0, seq2len = 0;
			memcpy(&rname2len, bin2+12, 1);
			memcpy(&cigar2len, bin2+16, 2);
			memcpy(&seq2len, bin2+20, 4);

			int bin2_ptr = 36 + rname2len + 4 * cigar2len + seq2len + (seq2len+1)/2;
			SAM_pairer_iterate_tags((unsigned char *)bin2+bin2_ptr, block2_len + 4 - bin2_ptr, "RG", &RG_type, RG_ptr);
			if(RG_type != 'Z') (*RG_ptr) = NULL;
		}

	}
}

/*
typedef struct {
	char chromosome_name_left[CHROMOSOME_NAME_LENGTH + 1];
	char chromosome_name_right[CHROMOSOME_NAME_LENGTH + 1];
	unsigned int last_exon_base_left;
	unsigned int first_exon_base_right;
} fc_junction_info_t;

*/
int calc_junctions_from_cigarInts(fc_thread_global_context_t * global_context, int alignment_masks , int cigar_sections, unsigned int * Starting_Chro_Points_1BASE, unsigned short * Starting_Read_Points, unsigned short * Section_Lengths, char ** ChroNames, char * Event_After_Section, fc_junction_info_t * junctions_current){
	int x1, ret = 0;
	unsigned int last_base_pos = Starting_Chro_Points_1BASE[0] + Section_Lengths[0] - 1;
	for(x1 = 1; x1 < cigar_sections; x1++){
		if(!ChroNames[x1]) continue; // NULL chro name for https://groups.google.com/forum/#!topic/subread/QDT6npjAZuE
		if(Event_After_Section[x1-1] == 'N'){
			unsigned int first_base_pos = Starting_Chro_Points_1BASE[x1];
			junctions_current[ret].last_exon_base_left = last_base_pos;
			junctions_current[ret].first_exon_base_right = first_base_pos;
			strcpy(junctions_current[ret].chromosome_name_left, ChroNames[x1]);
			strcpy(junctions_current[ret].chromosome_name_right, ChroNames[x1]);
			ret ++;
		}

		last_base_pos = Starting_Chro_Points_1BASE[x1] + Section_Lengths[x1] - 1;
	}
	return ret;
}

void add_fragment_supported_junction(	fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, fc_junction_info_t * supported_junctions1, int njunc1, fc_junction_info_t * supported_junctions2, int njunc2, char * RG_name);

void process_line_junctions(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * bin1, char * bin2) {
	fc_junction_info_t * supported_junctions1 = thread_context -> memory_supported_junctions1;
	fc_junction_info_t * supported_junctions2 = thread_context -> memory_supported_junctions2;
	int is_second_read, njunc1=0, njunc2=0, is_junction_read, cigar_sections, cigar_overflow=1;
	int alignment_masks, mapping_qual, NH_value;
	char *RG_ptr=NULL;

	unsigned int * Starting_Chro_Points_1BASE = malloc(sizeof(int)*global_context -> max_M);
	unsigned short * Starting_Read_Points = malloc(sizeof(short)*global_context -> max_M);
	unsigned short * Section_Read_Lengths = malloc(sizeof(short)*global_context -> max_M);
	char ** ChroNames = malloc(sizeof(char*)*global_context -> max_M);
 	char * Event_After_Section = malloc(sizeof(char)*global_context -> max_M);
	for(is_second_read = 0 ; is_second_read < 2; is_second_read++){
		char * read_chr, *read_name, *mate_chr;
		srInt_64 read_pos, fragment_length = 0, mate_pos;
		/*
		unsigned int Starting_Chro_Points_1BASE[global_context -> max_M];
		unsigned short Starting_Read_Points[global_context -> max_M];
		unsigned short Section_Read_Lengths[global_context -> max_M];
		char * ChroNames[global_context -> max_M];
		char Event_After_Section[global_context -> max_M];
		*/
		if(is_second_read && !global_context -> is_paired_end_mode_assign) break;
		char * RG_ptr_one = NULL;
		int me_refID, mate_refID;

		parse_bin(global_context -> sambam_chro_table, is_second_read?bin2:bin1, is_second_read?bin1:bin2 , &read_name,  &alignment_masks , &read_chr, &read_pos, &mapping_qual, &mate_chr, &mate_pos, &fragment_length, &is_junction_read, &cigar_sections, &cigar_overflow, Starting_Chro_Points_1BASE, Starting_Read_Points, Section_Read_Lengths, ChroNames, Event_After_Section, &NH_value, global_context -> max_M, NULL, NULL, global_context -> assign_reads_to_RG, &RG_ptr_one, &me_refID, &mate_refID);
		if(RG_ptr_one) RG_ptr = RG_ptr_one;

		int * njunc_current = is_second_read?&njunc2:&njunc1;
		fc_junction_info_t * junctions_current = is_second_read?supported_junctions2:supported_junctions1;
		(*njunc_current) = calc_junctions_from_cigarInts(global_context, alignment_masks , cigar_sections, Starting_Chro_Points_1BASE, Starting_Read_Points, Section_Read_Lengths, ChroNames, Event_After_Section, junctions_current);
	}
	if(njunc1 >0 || njunc2>0)
		add_fragment_supported_junction(global_context, thread_context, supported_junctions1, njunc1, supported_junctions2, njunc2, RG_ptr);
	free(Starting_Chro_Points_1BASE);
	free(Starting_Read_Points);
	free(Section_Read_Lengths);
	free(ChroNames);
	free(Event_After_Section);
}

void ** get_RG_tables(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * rg_name){
	void ** ret = HashTableGet(thread_context->RG_table, rg_name);
	if(ret) return ret;
	
	ret = malloc(sizeof(void *)*4);
	
	ret[0] = malloc(thread_context -> count_table_size * sizeof(read_count_type_t));
	ret[1] = malloc(sizeof(fc_read_counters));

	memset(ret[0], 0, thread_context -> count_table_size * sizeof(read_count_type_t));
	memset(ret[1], 0, sizeof(fc_read_counters));

	if(global_context -> do_junction_counting){
		HashTable * junction_counting_table = HashTableCreate(131317);
		HashTableSetHashFunction(junction_counting_table,HashTableStringHashFunction);
		HashTableSetDeallocationFunctions(junction_counting_table, free, NULL);
		HashTableSetKeyComparisonFunction(junction_counting_table, fc_strcmp_chro);
		
		HashTable * splicing_point_table = HashTableCreate(131317);
		HashTableSetHashFunction(splicing_point_table,HashTableStringHashFunction);
		HashTableSetDeallocationFunctions(splicing_point_table, free, NULL);
		HashTableSetKeyComparisonFunction(splicing_point_table, fc_strcmp_chro);

		ret [2] = junction_counting_table;
		ret [3] = splicing_point_table;
	}else ret[2] = NULL;
	
	char * rg_name_mem = strdup(rg_name);//malloc(strlen(rg_name)+1);
	if(NULL == rg_name_mem){
		SUBREADprintf("cannot allocate memory for %s (%d)\n", rg_name, strlen(rg_name));
		return NULL;
	}
	strcpy(rg_name_mem, rg_name);
	HashTablePut(thread_context->RG_table, rg_name_mem, ret);
	return ret;
}

int process_pairer_output(void * pairer_vp, int thread_no, char * bin1, char * bin2){
	SAM_pairer_context_t * pairer = (SAM_pairer_context_t *) pairer_vp;
	fc_thread_global_context_t * global_context = (fc_thread_global_context_t * )pairer -> appendix1;
	fc_thread_thread_context_t * thread_context = global_context -> thread_contexts + thread_no;

	if(pairer -> using_long_read_parser){
		if(0==global_context->is_read_too_long_to_SAM_BAM_shown &&(global_context -> is_read_details_out == FILE_TYPE_SAM || global_context -> is_read_details_out == FILE_TYPE_BAM)){
			global_context -> read_pairer.is_internal_error = 1;
			SUBREADprintf("ERROR: The read is too long to the SAM or BAM output.\nPlease use the 'CORE' mode for the assignment detail output. \n");
			global_context->is_read_too_long_to_SAM_BAM_shown = 1;
		}
	}

	process_line_buffer(global_context, thread_context, bin1, bin2);
	if(global_context -> do_junction_counting){// We decided to use all reads for junction counting on 08FEB2023. All reads with "N" in them will be counted against junctions, no matter whether they were used in read counting. 
		process_line_junctions(global_context, thread_context, bin1, bin2);
	}
	return 0;
}

void vote_and_add_count(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context,
			srInt_64 * hits_indices1, int nhits1, srInt_64 * hits_indices2, int nhits2, unsigned int total_frag_len,
			char ** hits_chro1, char ** hits_chro2, unsigned int * hits_start_pos1, unsigned int * hits_start_pos2, unsigned short * hits_length1, unsigned short * hits_length2, int fixed_fractional_count, char * read_name, char * RG_name, char * bin1, char * bin2);

void add_bin_new_tags(char * oldbin, char **newbin, char ** tags, char * types, void ** vals){
	int new_tags_length = 0;
	int tagi;

	add_bin_new_tags_reduce_longtag(oldbin);

	for(tagi = 0; tags[tagi]; tagi++){
		char type = types[tagi];
		if(type == 'i') new_tags_length += 7;
		else new_tags_length += 4 + strlen((char *)vals[tagi]);
	}

	int oldbin_len;
	memcpy(&oldbin_len, oldbin, 4);
	oldbin_len += 4;

	int newbin_len = oldbin_len + new_tags_length;
	(*newbin) = malloc(newbin_len);
	memcpy(*newbin, oldbin, oldbin_len);
	newbin_len -= 4;
	memcpy(*newbin, &newbin_len, 4);
	
	for(tagi = 0; tags[tagi]; tagi++){
		memcpy( (*newbin) + oldbin_len, tags[tagi] ,2);
		(*newbin)[oldbin_len+2] = types[tagi];
		if(types[tagi] == 'i'){
			int intv = vals[tagi] - NULL;
			memcpy((*newbin) + oldbin_len + 3, &intv, 4);
			oldbin_len += 7;
		}else{
			int vlen = strlen((char *)(vals[tagi]))+1;
			memcpy((*newbin) + oldbin_len + 3, vals[tagi], vlen);
			oldbin_len += 3 + vlen;
		}
	}
}


int compress_read_detail_BAM(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, int write_start, int write_end, char * bam_buf){
	if(global_context -> is_read_details_out == FILE_TYPE_SAM){
		// there MUST be only one read in the buffer.
		int write_ptr = write_start;
		int tmplen = 0 ;
		int sam_ptr = 0;
		while(1){
			if(write_ptr >= write_end) break;
			memcpy(&tmplen, thread_context -> read_details_buff + write_ptr, 4);
			tmplen +=4;
			int txtlen = convert_BAM_binary_to_SAM(global_context -> sambam_chro_table, thread_context -> read_details_buff + write_ptr, bam_buf + sam_ptr);
			bam_buf[sam_ptr + txtlen] = '\n';
			bam_buf[sam_ptr + txtlen + 1] = 0;
			sam_ptr += txtlen + 1;
			write_ptr += tmplen;
		}
		return sam_ptr;

	}else{
		// there may be multiple reads in the buffer.
			int bin_len = write_end - write_start;
			char * compressed_buff = bam_buf + 18;
					
			int compressed_size ; 
			unsigned int CRC32;
			thread_context -> bam_file_output_stream.avail_out = 66600;
			thread_context -> bam_file_output_stream.avail_in = bin_len;
			//SUBREADprintf("COMPRESS PTR=%p , LEN=%d\n", thread_context -> read_details_buff + write_start , bin_len);
			CRC32 = FC_CRC32(thread_context -> read_details_buff + write_start , bin_len);
					
			int Z_DEFAULT_MEM_LEVEL = 8;
			thread_context -> bam_file_output_stream.zalloc = Z_NULL;
			thread_context -> bam_file_output_stream.zfree = Z_NULL;
			thread_context -> bam_file_output_stream.opaque = Z_NULL;

			deflateInit2(&thread_context -> bam_file_output_stream, bin_len?Z_BEST_SPEED:Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
					
			thread_context -> bam_file_output_stream.next_in = (unsigned char*) thread_context -> read_details_buff + write_start;
			thread_context -> bam_file_output_stream.next_out = (unsigned char*) compressed_buff;

			deflate(&thread_context -> bam_file_output_stream, Z_FINISH);
			deflateEnd(&thread_context -> bam_file_output_stream);

			compressed_size = 66600 -thread_context -> bam_file_output_stream.avail_out;
					
			bam_buf[0]=31;
			bam_buf[1]=-117;
			bam_buf[2]=8;
			bam_buf[3]=4;
			memset(bam_buf+4, 0, 5);
			bam_buf[9] = 0xff;	// OS
				
			int tmpi = 6;
			memcpy(bam_buf+10, &tmpi, 2); //XLSN
			bam_buf[12]=66; // SI1
			bam_buf[13]=67; // SI2
			tmpi = 2;
			memcpy(bam_buf+14, &tmpi, 2); //BSIZE
			tmpi = compressed_size + 19 + 6;
			memcpy(bam_buf+16, &tmpi, 2); //BSIZE
				
			memcpy(bam_buf+18+compressed_size, &CRC32, 4);
			memcpy(bam_buf+18+compressed_size+4, &bin_len, 4);
			return 18+compressed_size+8;
	}
}

void write_read_detailed_remainder(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context){
	int write_bin_ptr = 0;
	int last_written_ptr = 0;
	int bam_compressed_buff_ptr = 0;

	if(thread_context -> read_details_buff_used <1)return;

	if(global_context -> is_read_details_out == FILE_TYPE_BAM && thread_context -> read_details_buff_used < 64000){
		bam_compressed_buff_ptr = compress_read_detail_BAM(global_context, thread_context, 0, thread_context -> read_details_buff_used, thread_context -> bam_compressed_buff);
	}else while(1){
		if(write_bin_ptr >= thread_context -> read_details_buff_used ) break;
		int tmplen = 0;
		memcpy(&tmplen, thread_context -> read_details_buff + write_bin_ptr, 4);
		if(tmplen < 9 || tmplen + 4 > FC_UNSAFE_BLOCK_SIZE){
			SUBREADprintf("ERROR: Format error : len = %d\n", tmplen);
			//oexit(-1);
			return ;
		}
		tmplen +=4;
		write_bin_ptr += tmplen;
		if(write_bin_ptr - last_written_ptr > 64000 || write_bin_ptr >= thread_context -> read_details_buff_used || global_context -> is_read_details_out == FILE_TYPE_SAM){
			bam_compressed_buff_ptr += compress_read_detail_BAM(global_context, thread_context, last_written_ptr, write_bin_ptr, thread_context -> bam_compressed_buff + bam_compressed_buff_ptr);
			last_written_ptr = write_bin_ptr;
		}
	}
	pthread_spin_lock(&global_context -> read_details_out_lock);
	fwrite(thread_context -> bam_compressed_buff, 1, bam_compressed_buff_ptr , global_context -> read_details_out_FP);
	pthread_spin_unlock(&global_context -> read_details_out_lock);
	thread_context -> read_details_buff_used =0;
}

int add_read_detail_bin_buff(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context,  char * bin, int do_write){
	int binlen=0;

	memcpy(&binlen, bin, 4);
	binlen += 4;
	if(binlen > FC_UNSAFE_BLOCK_SIZE){
		if(!global_context->is_read_too_long_to_SAM_BAM_shown){
			global_context -> read_pairer.is_internal_error = 1;
			SUBREADprintf("ERROR: The read is too long to the SAM or BAM output.\nPlease use the 'CORE' mode for the assignment detail output.\n");
			global_context->is_read_too_long_to_SAM_BAM_shown = 1;
		}
		return -1;
	}

	memcpy(thread_context -> read_details_buff + thread_context -> read_details_buff_used, bin, binlen);
	thread_context -> read_details_buff_used  += binlen;

	if(do_write){
		if(global_context -> is_read_details_out == FILE_TYPE_SAM || thread_context -> read_details_buff_used >= 65000 - FC_UNSAFE_BLOCK_SIZE) write_read_detailed_remainder(global_context, thread_context);
	}
	return 0;
}

int write_read_details_FP(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * status, int feature_count, char * features, char * bin1, char * bin2){
	int ret = 1;

	char * read_name;
	
	if(global_context -> is_read_details_out == FILE_TYPE_RSUBREAD){
		get_readname_from_bin(bin1?bin1:bin2, &read_name);
		fprintf(global_context -> read_details_out_FP, "%s\t%s\t%d\t%s\n", read_name, status, feature_count, features?features:"NA");
	}else{
		char * out_bin1 = NULL, *out_bin2 = NULL;
		char * tags[4];
		char types[4];
		void * vals[4];

		tags[0]="XS";
		tags[1]=feature_count >0?"XN":NULL;
		tags[2]=feature_count >0?"XT":NULL;
		tags[3]=NULL;
		types[0]='Z';
		types[1]='i';
		types[2]='Z';
		vals[0]=status;
		vals[1]=NULL+feature_count;
		vals[2]=features;
		
		if(bin1){
			add_bin_new_tags(bin1, &out_bin1, tags, types, vals);
			add_read_detail_bin_buff(global_context, thread_context, out_bin1, bin2 == NULL);
			free(out_bin1);
		}

		if(bin2){
			add_bin_new_tags(bin2, &out_bin2, tags, types, vals);
			add_read_detail_bin_buff(global_context, thread_context, out_bin2, 1);
			free(out_bin2);
		}
	}
	if(ret < 1) global_context -> disk_is_full = 1;
	return ret;
}

void warning_anno_BAM_chromosomes(fc_thread_global_context_t * global_context){
	int x1;
	HashTable * BAM_chro_tab = HashTableCreate(1117);
	HashTableSetHashFunction(BAM_chro_tab,HashTableStringHashFunction);
	HashTableSetKeyComparisonFunction(BAM_chro_tab,fc_strcmp_chro );

	for(x1 = 0; x1 < global_context -> sambam_chro_table_items; x1++){
		char * BAM_chro = global_context -> sambam_chro_table[x1].chro_name;
		if( global_context -> BAM_chros_to_anno_table){
			char * tmp_chro = HashTableGet(global_context -> BAM_chros_to_anno_table, global_context -> sambam_chro_table[x1].chro_name);
			if(tmp_chro) BAM_chro = tmp_chro;
		}
		HashTablePut(BAM_chro_tab, BAM_chro, NULL+1);
	}

	HashTable * ANNO_chro_tab = HashTableCreate(1117);
	HashTableSetHashFunction(ANNO_chro_tab,HashTableStringHashFunction);
	HashTableSetKeyComparisonFunction(ANNO_chro_tab,fc_strcmp_chro );

	for(x1 = 0 ; x1 < global_context -> exontable_exons ; x1++)
		HashTablePut(ANNO_chro_tab, global_context -> exontable_chr[x1], NULL+1);

	if(global_context -> is_verbose){
		warning_hash_hash(ANNO_chro_tab, BAM_chro_tab, "Chromosomes/contigs in annotation but not in input file");
		warning_hash_hash(BAM_chro_tab, ANNO_chro_tab, "Chromosomes/contigs in input file but not in annotation");
	}
	HashTableDestroy(BAM_chro_tab);
	HashTableDestroy(ANNO_chro_tab);
}

void process_line_buffer(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, char * bin1, char * bin2)
{
	if(global_context -> is_input_bad_format) return;
	char * read_chr, *read_name, *mate_chr;
	srInt_64 read_pos, fragment_length = 0, mate_pos;
	unsigned int search_start = 0, search_end;
	int nhits1 = 0, nhits2 = 0, alignment_masks, search_block_id, search_item_id, mapping_qual;
	unsigned int  total_frag_len =0;

	int cigar_sections, is_junction_read;
	unsigned int * Starting_Chro_Points_1BASE = thread_context -> proc_Starting_Chro_Points_1BASE;
	unsigned short * Starting_Read_Points = thread_context -> proc_Starting_Read_Points;
	unsigned short * Section_Read_Lengths = thread_context -> proc_Section_Read_Lengths;
	char ** ChroNames = thread_context -> proc_ChroNames;
	char * Event_After_Section = thread_context -> proc_Event_After_Section;

	CIGAR_interval_t * CIGAR_intervals_R1 = thread_context -> proc_CIGAR_intervals_R1;
	CIGAR_interval_t * CIGAR_intervals_R2 = thread_context -> proc_CIGAR_intervals_R2;

	int is_second_read;
	int maximum_NH_value = 1, NH_value;
	int skipped_for_exonic = 0;
	int first_read_quality_score = 0, CIGAR_intervals_R1_sections = 0, CIGAR_intervals_R2_sections = 0;

	if(thread_context -> thread_id == 0 && thread_context -> all_reads < 1){
		warning_anno_BAM_chromosomes(global_context);
	}

	if(global_context -> need_calculate_overlap_len ){
		// we only need to clean the 1st item. The following items are cleaned in parse_bin.
		memset( CIGAR_intervals_R1, 0, sizeof(CIGAR_interval_t));
		memset( CIGAR_intervals_R2, 0, sizeof(CIGAR_interval_t));
	}

	thread_context->all_reads++;

	char * RG_ptr;
	int me_refID =-1, mate_refID =-1, this_is_inconsistent_read_type = 0;
	for(is_second_read = 0 ; is_second_read < 2; is_second_read++)
	{
		if(is_second_read && !global_context -> is_paired_end_mode_assign) break;

		RG_ptr = NULL;
		int cigar_overflow=0;
		parse_bin(global_context -> sambam_chro_table, is_second_read?bin2:bin1, is_second_read?bin1:bin2 , &read_name,  &alignment_masks , &read_chr, &read_pos, &mapping_qual, &mate_chr, &mate_pos, &fragment_length, &is_junction_read, &cigar_sections, &cigar_overflow, Starting_Chro_Points_1BASE, Starting_Read_Points, Section_Read_Lengths, ChroNames, Event_After_Section, &NH_value, global_context -> max_M , global_context -> need_calculate_overlap_len?(is_second_read?CIGAR_intervals_R2:CIGAR_intervals_R1):NULL, is_second_read?&CIGAR_intervals_R2_sections:&CIGAR_intervals_R1_sections, global_context -> assign_reads_to_RG, &RG_ptr, &me_refID, &mate_refID);
		//fprintf(stderr,"READ NAME %s FLAG %d\n",read_name, alignment_masks);

		if(cigar_overflow && 0==global_context -> is_cigar_overflow_warning_shown) {
			print_in_box(80,0,0,"");
			print_in_box(80,0,0,"WARNING: more than %d 'M' segments were found in CIGAR string.", global_context -> max_M);
			print_in_box(80,0,0,"The exceeding sections are not useed. Please consider to increase maxMOp.");
			print_in_box(80,0,0,"");
			global_context -> is_cigar_overflow_warning_shown = 1;
		}

		// this will be done in the other function.
		if(global_context -> is_paired_end_mode_assign && (alignment_masks&1)==0) alignment_masks|=8;

		if(global_context -> assign_reads_to_RG && NULL == RG_ptr)return;

		if(  ( alignment_masks & SAM_FLAG_PAIRED_TASK ) && !global_context -> any_reads_are_PE ) global_context -> any_reads_are_PE=1;
		if(((!global_context -> is_paired_end_reads_expected)  && ( alignment_masks & SAM_FLAG_PAIRED_TASK )) || ((global_context -> is_paired_end_reads_expected)  && 0 == ( alignment_masks & SAM_FLAG_PAIRED_TASK ))){
			if(global_context -> is_mixed_PE_SE == 0) global_context -> is_mixed_PE_SE =1;
			if(!global_context -> is_paired_end_reads_expected){
				SUBREADprintf("ERROR: Paired-end reads were detected in single-end read library : %s\n", global_context -> input_file_name);
				global_context -> is_input_bad_format = 1;
				return;
			}
			this_is_inconsistent_read_type = 1;
		}

		if(is_second_read == 0)
		{
			//skip the read if unmapped (its mate will be skipped as well if paired-end)
			if( ((!global_context -> is_paired_end_mode_assign) &&  (alignment_masks & SAM_FLAG_UNMAPPED) ) ||
			    ((alignment_masks & SAM_FLAG_UNMAPPED)   &&  (alignment_masks & SAM_FLAG_MATE_UNMATCHED) && global_context -> is_paired_end_mode_assign)) { 
				if(RG_ptr){
					void ** tab4s = get_RG_tables(global_context, thread_context, RG_ptr);
					fc_read_counters * sumtab = tab4s[1];
					sumtab -> unassigned_unmapped++;
				}else
					thread_context->read_counters.unassigned_unmapped ++;

				if(global_context -> read_details_out_FP)
					write_read_details_FP(global_context , thread_context ,"Unassigned_Unmapped",0, NULL, bin1, bin2);

				return;	// do nothing if a read is unmapped, or the first read in a pair of reads is unmapped.
			}
		}

		if(((alignment_masks & SAM_FLAG_UNMAPPED) || (alignment_masks & SAM_FLAG_MATE_UNMATCHED)) && global_context -> is_paired_end_mode_assign && global_context -> is_both_end_required){
				if(RG_ptr){
					void ** tab4s = get_RG_tables(global_context, thread_context, RG_ptr);
					fc_read_counters * sumtab = tab4s[1];
					sumtab -> unassigned_singleton++;
				}else
					thread_context->read_counters.unassigned_singleton ++;

				if(global_context -> read_details_out_FP)
					write_read_details_FP(global_context , thread_context ,"Unassigned_Singleton",0, NULL, bin1, bin2);

				return;
		}

		if(this_is_inconsistent_read_type){
			if(global_context -> is_both_end_required && 0 == ( alignment_masks & SAM_FLAG_PAIRED_TASK )){
				if(global_context -> read_details_out_FP)
					write_read_details_FP(global_context , thread_context ,"Unassigned_Singleton",0, NULL, bin1, bin2);
				if(RG_ptr){
					void ** tab4s = get_RG_tables(global_context, thread_context, RG_ptr);
					fc_read_counters * sumtab = tab4s[1];
					sumtab -> unassigned_singleton++;
				}else thread_context->read_counters.unassigned_singleton ++;

				return; // when running on PE mode, SE reads are seen as "only one end mapped"
			}
		}

		if(global_context -> min_mapping_quality_score>0)
		{
			//printf("SECOND=%d; FIRST=%d; THIS=%d; Q=%d\n", is_second_read, first_read_quality_score, mapping_qual, );
			if(( mapping_qual < global_context -> min_mapping_quality_score  && ! global_context -> is_paired_end_mode_assign)||( is_second_read  && max( first_read_quality_score, mapping_qual ) < global_context -> min_mapping_quality_score))
			{
				thread_context->read_counters.unassigned_mappingquality ++;

				if(global_context -> read_details_out_FP)
				{
					write_read_details_FP(global_context, thread_context, "Unassigned_MappingQuality", 0, NULL, bin1, bin2);
				}

				return;
			}
			if(is_second_read==0 && global_context -> is_paired_end_mode_assign)
			{
				first_read_quality_score = mapping_qual;
			}
		}

		if(is_second_read == 0 && global_context -> is_paired_end_mode_assign && 
	   	  (global_context -> is_PE_distance_checked || global_context -> is_chimertc_disallowed)
		  )
		{
			int is_half_mapped = (alignment_masks & SAM_FLAG_UNMAPPED) || (alignment_masks & SAM_FLAG_MATE_UNMATCHED);

			if(!is_half_mapped)
			{
				fragment_length = abs( fragment_length ); //get the fragment length

				int is_first_read_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
				int is_second_read_negative_strand = (alignment_masks & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)?1:0; 

				if(mate_chr == read_chr && is_first_read_negative_strand!=is_second_read_negative_strand) {
				 //^^^^^^^^^^^^^^^^^^^^ They are directly compared because they are both pointers in the same contig name table.
				 //
					if(global_context -> is_PE_distance_checked && ((fragment_length > global_context -> max_paired_end_distance) || (fragment_length < global_context -> min_paired_end_distance))) {
						if(RG_ptr){
							void ** tab4s = get_RG_tables(global_context, thread_context, RG_ptr);
							fc_read_counters * sumtab = tab4s[1];
							sumtab -> unassigned_fragmentlength++;
						}else
							thread_context->read_counters.unassigned_fragmentlength ++;

						if(global_context -> read_details_out_FP)
							write_read_details_FP(global_context, thread_context, "Unassigned_FragmentLength", -1, NULL, bin1, bin2);
						return;
					}
				} else {
					if(global_context -> is_chimertc_disallowed) {
						if(RG_ptr){
							void ** tab4s = get_RG_tables(global_context, thread_context, RG_ptr);
							fc_read_counters * sumtab = tab4s[1];
							sumtab -> unassigned_chimericreads++;
						}else
							thread_context->read_counters.unassigned_chimericreads ++;

						if(global_context -> read_details_out_FP)
							write_read_details_FP(global_context, thread_context, "Unassigned_Chimera", -1, NULL, bin1, bin2);
						return;
					}
				}
			}
		}

		// This filter has to be put here because the 0x400 FLAG is not about mapping but about sequencing.
		// A unmapped read with 0x400 FLAG should be able to kill the mapped mate which may have no 0x400 FLAG. 
		if(global_context -> is_duplicate_ignored)
		{
			if(alignment_masks & SAM_FLAG_DUPLICATE)
			{
				if(RG_ptr){
					void ** tab4s = get_RG_tables(global_context, thread_context, RG_ptr);
					fc_read_counters * sumtab = tab4s[1];
					sumtab -> unassigned_duplicate++;
				}else thread_context->read_counters.unassigned_duplicate ++;
				if(global_context -> read_details_out_FP)
					write_read_details_FP(global_context, thread_context, "Unassigned_Duplicate", -1, NULL, bin1, bin2);

				return;
			}

		}

		if(SAM_FLAG_UNMAPPED & alignment_masks) continue;

		if( NH_value > 1 ) {
			if(global_context -> is_multi_mapping_allowed == 0) {
				// now it is a NH>1 read!
				// not allow multimapping -> discard!
				if(RG_ptr){
					void ** tab4s = get_RG_tables(global_context, thread_context, RG_ptr);
					fc_read_counters * sumtab = tab4s[1];
					sumtab -> unassigned_multimapping++;
				}else thread_context->read_counters.unassigned_multimapping ++;

				if(global_context -> read_details_out_FP)
					write_read_details_FP(global_context, thread_context, "Unassigned_MultiMapping", -1, NULL, bin1, bin2);

				return;
			}
		}

		maximum_NH_value = max(maximum_NH_value, NH_value);

		// if a pair of reads have one secondary, the entire fragment is seen as secondary.
		if((alignment_masks & SAM_FLAG_SECONDARY_MAPPING) && (global_context -> is_primary_alignment_only)) {
			if(RG_ptr){
				void ** tab4s = get_RG_tables(global_context, thread_context, RG_ptr);
				fc_read_counters * sumtab = tab4s[1];
				sumtab -> unassigned_secondary++;
			}else thread_context->read_counters.unassigned_secondary ++;

			if(global_context -> read_details_out_FP)
				write_read_details_FP(global_context, thread_context, "Unassigned_Secondary", -1, NULL, bin1, bin2);
			return;
		}

		int is_this_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
		int is_fragment_negative_strand = is_this_negative_strand;

		if(1 || global_context -> is_paired_end_mode_assign){ // On 20 JULY 2020: If strand-specific counting is on, isPairedEnd = TRUE, countReadPairs = FALSE and the BAM file contains mixed reads, then the single-end reads and the R1 reads in read-pairs will be directly compared with the strand of the gene, but the R2 reads in read-pairs will be compared with the opposite strand of the gene. A read-pair will be counted twice no matter if the strand-specific mode is on or off. If the argument to the strand-specific option is "1", then R1 must have the same strand of the gene and R2 must have the opposite strand of the gene to be counted. 
			int is_second_read_in_pair = alignment_masks & SAM_FLAG_SECOND_READ_IN_PAIR;
			//is_fragment_negative_strand = is_second_read_in_pair?(!is_this_negative_strand):is_this_negative_strand;
			if(is_second_read_in_pair)
				is_fragment_negative_strand = global_context -> is_second_read_straight?is_this_negative_strand:(!is_this_negative_strand);
			else
				is_fragment_negative_strand = global_context -> is_first_read_reversed?(!is_this_negative_strand):is_this_negative_strand;
		}

		int nhits = 0;

		int cigar_section_id;
		srInt_64 * hits_indices = is_second_read?thread_context -> hits_indices2:thread_context -> hits_indices1;
		unsigned int * hits_start_pos = is_second_read?thread_context -> hits_start_pos2:thread_context -> hits_start_pos1;
		unsigned short * hits_length = is_second_read?thread_context -> hits_length2:thread_context -> hits_length1;
		char ** hits_chro = is_second_read?thread_context -> hits_chro2:thread_context -> hits_chro1;

		if(global_context->is_split_or_exonic_only == 1 && !is_junction_read) {
			skipped_for_exonic ++;

			if(skipped_for_exonic == 1 + global_context -> is_paired_end_mode_assign){
				if(global_context -> read_details_out_FP)
					write_read_details_FP(global_context, thread_context, (global_context->is_split_or_exonic_only == 2)?"Unassigned_Split":"Unassigned_NonSplit", -1, NULL, bin1, bin2);

				if(RG_ptr){
					void ** tab4s = get_RG_tables(global_context, thread_context, RG_ptr);
					fc_read_counters * sumtab = tab4s[1];
					sumtab -> unassigned_junction_condition++;
				}else thread_context->read_counters.unassigned_junction_condition ++;

				return;
			}
		}


		if(global_context->is_split_or_exonic_only == 2 && is_junction_read) {
			if(global_context -> read_details_out_FP)
				write_read_details_FP(global_context, thread_context,(global_context->is_split_or_exonic_only == 2)?"Unassigned_Split":"Unassigned_NonSplit", -1, NULL, bin1, bin2);
			if(RG_ptr){
				void ** tab4s = get_RG_tables(global_context, thread_context, RG_ptr);
				fc_read_counters * sumtab = tab4s[1];
				sumtab -> unassigned_junction_condition++;
			}else thread_context->read_counters.unassigned_junction_condition ++;

			return;
		}

		if(1) {

			if(0)SUBREADprintf("MAPPED R_%d to %s : CHR_POS=%u + %u, CHR_LEN=%u\n", is_second_read+1, global_context -> sambam_chro_table[me_refID]. chro_name, Starting_Chro_Points_1BASE[0], Section_Read_Lengths[0], global_context -> sambam_chro_table[me_refID] .chro_length );

			if(global_context -> read_shift_size>0){
				int shifting_applied_length = 0;
				int shifting_i;

				if((global_context -> read_shift_type == READ_SHIFT_UPSTREAM   && (!is_this_negative_strand))||
					(global_context -> read_shift_type == READ_SHIFT_DOWNSTREAM &&   is_this_negative_strand ))
					shifting_applied_length = -global_context -> read_shift_size;

				if((global_context -> read_shift_type == READ_SHIFT_UPSTREAM   &&   is_this_negative_strand)||
					(global_context -> read_shift_type == READ_SHIFT_DOWNSTREAM && (!is_this_negative_strand)))
					shifting_applied_length = global_context -> read_shift_size;

				if(global_context -> read_shift_type == READ_SHIFT_LEFT) shifting_applied_length = -global_context -> read_shift_size;
				if(global_context -> read_shift_type == READ_SHIFT_RIGHT) shifting_applied_length = global_context -> read_shift_size;

				if(shifting_applied_length < 0 && Starting_Chro_Points_1BASE[0] <=-shifting_applied_length) shifting_applied_length = - Starting_Chro_Points_1BASE[0]+1;
				if(shifting_applied_length > 0 && Starting_Chro_Points_1BASE[cigar_sections-1] + Section_Read_Lengths[cigar_sections-1] + shifting_applied_length > global_context -> sambam_chro_table[me_refID].chro_length +1 )
					shifting_applied_length =  global_context -> sambam_chro_table[me_refID].chro_length +1 - (Starting_Chro_Points_1BASE[cigar_sections-1] + Section_Read_Lengths[cigar_sections-1]);

				for(shifting_i = 0; shifting_i < cigar_sections ; shifting_i++)
						Starting_Chro_Points_1BASE[shifting_i] += shifting_applied_length;
			}

			if(global_context -> five_end_extension)
			{
				if(is_this_negative_strand){
					int applied_ext = global_context -> five_end_extension;

					if( Starting_Chro_Points_1BASE[cigar_sections-1] + Section_Read_Lengths[cigar_sections-1] + applied_ext > global_context -> sambam_chro_table[me_refID].chro_length +1  )
						applied_ext =  global_context -> sambam_chro_table[me_refID].chro_length +1 - (Starting_Chro_Points_1BASE[cigar_sections-1] + Section_Read_Lengths[cigar_sections-1]);

					Section_Read_Lengths [cigar_sections - 1] += applied_ext;
				}else{
					//SUBREADprintf("5-end extension: %d [%d]\n", Starting_Chro_Points_1BASE[0], Section_Lengths[0]);
					if( read_pos > global_context -> five_end_extension)
					{
						Section_Read_Lengths [0] += global_context -> five_end_extension;
						Starting_Chro_Points_1BASE [0] -= global_context -> five_end_extension;
					}
					else
					{
						Section_Read_Lengths [0] += read_pos-1;
						Starting_Chro_Points_1BASE [0] -= read_pos-1;
					}
				}
			}

			if(global_context -> three_end_extension) {

				if(is_this_negative_strand){
					if( read_pos > global_context -> three_end_extension)
					{
						Section_Read_Lengths [0] += global_context -> three_end_extension;
						Starting_Chro_Points_1BASE [0] -= global_context -> three_end_extension;
					}
					else
					{
						Section_Read_Lengths [0] += read_pos - 1;
						Starting_Chro_Points_1BASE [0] -= read_pos - 1;
					}
				} else{
					int applied_ext = global_context -> three_end_extension;
					if( Starting_Chro_Points_1BASE[cigar_sections-1] + Section_Read_Lengths[cigar_sections-1] + applied_ext > global_context -> sambam_chro_table[me_refID].chro_length +1 )
						applied_ext = global_context -> sambam_chro_table[me_refID].chro_length +1 - (Starting_Chro_Points_1BASE[cigar_sections-1] + Section_Read_Lengths[cigar_sections-1]);
					Section_Read_Lengths [cigar_sections - 1] += applied_ext;
				}

			}

			if(global_context -> reduce_5_3_ends_to_one) {
				if((REDUCE_TO_5_PRIME_END == global_context -> reduce_5_3_ends_to_one) + is_this_negative_strand == 1) // reduce to 5' end (small coordinate if positive strand / large coordinate if negative strand)
				{
					Section_Read_Lengths[0]=1;
				}
				else
				{
					Starting_Chro_Points_1BASE[0] = Starting_Chro_Points_1BASE[cigar_sections-1] + Section_Read_Lengths[cigar_sections-1] - 1;
					Section_Read_Lengths[0]=1;
				}
				cigar_sections = 1;
			}

			for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
			{
				
				if(!ChroNames[ cigar_section_id ]) continue; // NULL chro name for https://groups.google.com/forum/#!topic/subread/QDT6npjAZuE
				srInt_64 section_begin_pos = Starting_Chro_Points_1BASE[cigar_section_id];
				srInt_64 section_end_pos = Section_Read_Lengths[cigar_section_id] + section_begin_pos - 1;

				
				int start_reverse_table_index = section_begin_pos / REVERSE_TABLE_BUCKET_LENGTH;
				int end_reverse_table_index = (1+section_end_pos) / REVERSE_TABLE_BUCKET_LENGTH;

				/*if(ChroNames[cigar_section_id] < (char *)NULL + 0xfffff){
					unsigned char * tbbin = is_second_read?bin2:bin1;
					int * refid = (int*)(tbbin);
					
					SUBREADprintf("DANGEROUS! RNAME=%s, REC_LEN=%d,  CNAME=[%d]%p,  LEN_P=%d,  SECID=%d\n", read_name, refid[0], refid[1], ChroNames[cigar_section_id], Section_Read_Lengths[cigar_section_id], cigar_section_id);
				}*/

				fc_chromosome_index_info * this_chro_info = HashTableGet(global_context -> exontable_chro_table, ChroNames[cigar_section_id]);
				if(this_chro_info == NULL)
				{
					if(global_context -> BAM_chros_to_anno_table)
					{
						char * anno_chro_name = HashTableGet( global_context -> BAM_chros_to_anno_table , ChroNames[cigar_section_id]);
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

									if(nhits >= thread_context -> hits_number_capacity - 1){
										//SUBREADprintf("RESIZE hits: %d\n", thread_context -> hits_number_capacity);
										thread_context -> hits_number_capacity = thread_context -> hits_number_capacity/2 * 3;
										thread_context -> hits_number_capacity = max(10, thread_context -> hits_number_capacity);
										thread_context -> hits_start_pos1 = realloc(thread_context -> hits_start_pos1 , sizeof(int) * thread_context -> hits_number_capacity);
										thread_context -> hits_start_pos2 = realloc(thread_context -> hits_start_pos2 , sizeof(int) * thread_context -> hits_number_capacity);

										thread_context -> hits_length1 = realloc(thread_context -> hits_length1, sizeof(short) * thread_context -> hits_number_capacity);
										thread_context -> hits_length2 = realloc(thread_context -> hits_length2, sizeof(short) * thread_context -> hits_number_capacity);

										thread_context -> hits_chro1 = realloc(thread_context -> hits_chro1, sizeof(char *) * thread_context -> hits_number_capacity);
										thread_context -> hits_chro2 = realloc(thread_context -> hits_chro2, sizeof(char *) * thread_context -> hits_number_capacity);

										thread_context -> hits_indices1 = realloc(thread_context -> hits_indices1, sizeof(srInt_64) * thread_context -> hits_number_capacity);
										thread_context -> hits_indices2 = realloc(thread_context -> hits_indices2, sizeof(srInt_64) * thread_context -> hits_number_capacity);

										thread_context -> scoring_buff_numbers = realloc(thread_context -> scoring_buff_numbers, sizeof(int)*2*thread_context -> hits_number_capacity);
										thread_context -> scoring_buff_flags = realloc(thread_context -> scoring_buff_flags, sizeof(int)*2*thread_context -> hits_number_capacity);
										thread_context -> scoring_buff_overlappings = realloc(thread_context -> scoring_buff_overlappings, sizeof(int)*2*thread_context -> hits_number_capacity);
										thread_context -> scoring_buff_exon_ids = realloc(thread_context -> scoring_buff_exon_ids, sizeof(srInt_64)*2*thread_context -> hits_number_capacity);

										if(global_context -> need_calculate_overlap_len){
											thread_context -> scoring_buff_gap_chros = realloc(thread_context -> scoring_buff_gap_chros, sizeof(char *) * 2 * global_context -> max_M *2 * thread_context -> hits_number_capacity);
											thread_context -> scoring_buff_gap_starts = realloc(thread_context -> scoring_buff_gap_starts, sizeof(int) * 2 * global_context -> max_M *2 * thread_context -> hits_number_capacity);
											thread_context -> scoring_buff_gap_lengths = realloc(thread_context -> scoring_buff_gap_lengths, sizeof(short) * 2 * global_context -> max_M *2 * thread_context -> hits_number_capacity);
										}

										hits_indices = is_second_read?thread_context -> hits_indices2:thread_context -> hits_indices1;
										hits_start_pos = is_second_read?thread_context -> hits_start_pos2:thread_context -> hits_start_pos1;
										hits_length = is_second_read?thread_context -> hits_length2:thread_context -> hits_length1;
										hits_chro = is_second_read?thread_context -> hits_chro2:thread_context -> hits_chro1;
										//SUBREADprintf("RESIZE hits2: %d\n", thread_context -> hits_number_capacity);
									}

									if(nhits <= MAX_HIT_NUMBER - 1) {
										hits_indices[nhits] = search_item_id;

										if(global_context -> need_calculate_overlap_len) {
											hits_start_pos[nhits] = max(Starting_Chro_Points_1BASE[cigar_section_id], global_context -> exontable_start[search_item_id]);
											hits_length[nhits] =  min(global_context -> exontable_stop[search_item_id] , section_end_pos)+1 - hits_start_pos[nhits] ;
											hits_chro[nhits] = ChroNames[cigar_section_id];
											if(0 && FIXLENstrcmp("V0112_0155:7:1101:10214:3701", read_name)==0)
												SUBREADprintf("QNAME: [%d] %s %d ~ %d\n", nhits, hits_chro[nhits],  hits_start_pos[nhits],  hits_start_pos[nhits]+hits_length[nhits]);
										}

										nhits++;
									} else {
										SUBREADprintf("ERROR: the read overlapped with more than %d features.\n", nhits);
										global_context -> is_input_bad_format = 1;
										return ;
									}
								}
							} 
						}
					}
				}
			}
		}


		if(is_second_read) nhits2 = nhits;
		else	nhits1 = nhits;
	}	// loop for is_second_read


	if(0&&global_context -> do_junction_counting)// junction reads that passed the basic filters will be considered with the junction counting. Filters: Unmapped, Singleton, MAPQ, TemplateLength, Chimeric, Duplicate, Multimapping, Secondary alignment, Junction-containing status. But this was reversed on 08FEB2023. Hence it isn't executed here. 
	        process_line_junctions(global_context, thread_context, bin1, bin2);

	if(global_context -> need_calculate_fragment_len )
		total_frag_len = calc_total_frag_len( global_context, thread_context, CIGAR_intervals_R1, CIGAR_intervals_R1_sections, CIGAR_intervals_R2, CIGAR_intervals_R2_sections , read_name);

	//SUBREADprintf("FRAGLEN: %s %d; CIGARS=%d,%d\n", read_name, total_frag_len, CIGAR_intervals_R1_sections,CIGAR_intervals_R2_sections);

	int fixed_fractional_count = ( global_context -> use_fraction_multi_mapping && ! global_context -> is_primary_alignment_only )?calc_fixed_fraction(maximum_NH_value): NH_FRACTION_INT;

	// we have hits_indices1 and hits_indices2 and nhits1 and nhits2 here
	// we also have fixed_fractional_count which is the value to add

	vote_and_add_count(global_context, thread_context,
			    thread_context -> hits_indices1,  nhits1, thread_context -> hits_indices2,  nhits2, total_frag_len,
			    thread_context -> hits_chro1, thread_context -> hits_chro2,
				thread_context -> hits_start_pos1, thread_context -> hits_start_pos2,
				thread_context -> hits_length1, thread_context ->hits_length2,
			    fixed_fractional_count, read_name, RG_ptr, bin1, bin2);
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

void add_fragment_supported_junction(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context, fc_junction_info_t * supported_junctions1, int njunc1, fc_junction_info_t * supported_junctions2, int njunc2, char * RG_name){
	assert(njunc1 >= 0 && njunc1 <= global_context -> max_M -1 );
	assert(njunc2 >= 0 && njunc2 <= global_context -> max_M -1 );
	int x1,x2, in_total_junctions = njunc2 + njunc1;
	
	HashTable * junction_counting_table, *splicing_point_table;
	
	if(RG_name){
		void ** tab4s = get_RG_tables(global_context, thread_context, RG_name);
		junction_counting_table = tab4s[2];
		splicing_point_table = tab4s[3];
	}else{
		junction_counting_table = thread_context -> junction_counting_table;
		splicing_point_table = thread_context -> splicing_point_table;
	}
	
	for(x1 = 0; x1 < in_total_junctions; x1 ++){
		fc_junction_info_t * j_one = (x1 >= njunc1)?supported_junctions2+(x1-njunc1):(supported_junctions1+x1);
		if(j_one->chromosome_name_left[0]==0) continue;

		for(x2 = x1+1; x2 < in_total_junctions ; x2 ++){
			fc_junction_info_t * j_two = (x2 >= njunc1)?supported_junctions2+(x2-njunc1):(supported_junctions1+x2);
			if(j_two->chromosome_name_left[0]==0) continue;
			if(
				j_one -> last_exon_base_left == j_two -> last_exon_base_left &&
				j_one -> first_exon_base_right == j_two -> first_exon_base_right &&
				strcmp(j_one -> chromosome_name_left, j_two -> chromosome_name_left) == 0 &&
				strcmp(j_one -> chromosome_name_right, j_two -> chromosome_name_right) == 0
			) j_two -> chromosome_name_left[0]=0;
		}

		size_t key_maxlen = strlen(j_one->chromosome_name_left) + strlen(j_one->chromosome_name_right)  + 36;
		char * this_key = malloc(key_maxlen);
		SUBreadSprintf(this_key, key_maxlen, "%s\t%u\t%s\t%u", j_one->chromosome_name_left, j_one -> last_exon_base_left, j_one->chromosome_name_right, j_one -> first_exon_base_right);
		void * count_ptr = HashTableGet(junction_counting_table, this_key);
		srInt_64 count_junc = count_ptr - NULL;
		HashTablePut(junction_counting_table, this_key, NULL+count_junc + 1);
		size_t left_keysize = strlen(j_one->chromosome_name_left) + 16;
		size_t right_keysize = strlen(j_one->chromosome_name_right) + 16;
		char * left_key = malloc(left_keysize);
		char * right_key = malloc(right_keysize);
		SUBreadSprintf(left_key, left_keysize , "%s\t%u", j_one->chromosome_name_left, j_one -> last_exon_base_left);
		SUBreadSprintf(right_key, right_keysize, "%s\t%u", j_one->chromosome_name_right, j_one -> first_exon_base_right);

		for( x2 = 0 ; x2 < 2 ; x2++ ){
			char * lr_key = x2?right_key:left_key;
			count_ptr = HashTableGet(splicing_point_table, lr_key);
			count_junc = count_ptr - NULL;
			HashTablePut(splicing_point_table, lr_key, NULL + count_junc + 1);
		}
	}
}

int overlap_compare(void * arr, int L, int R){
	unsigned int * pos = (unsigned int *)arr;
	return pos[ L*2 ] -  pos[R*2];
}

void overlap_exchange(void * arr, int L, int R){
	unsigned int * pos = (unsigned int *)arr, tt;
	tt=pos[L*2];
	pos[L*2] = pos[R*2];
	pos[R*2] = tt;

	tt=pos[L*2+1];
	pos[L*2+1] = pos[R*2+1];
	pos[R*2+1] = tt;
}

unsigned int calc_score_overlaps(fc_thread_global_context_t * global_context,  fc_thread_thread_context_t * thread_context, char ** chros, unsigned int * start_poses, unsigned short * lens, int sections, char * read_name){
	unsigned int in_intervals[ 2*sections ];
	unsigned int out_intervals[ 2*sections ], x1;
	char used_interval[ sections ];

	memset(used_interval, 0 , sections);
	unsigned int ret = 0;

	for(x1 = 0  ; x1 < sections ; x1++){
		if( used_interval [x1] )continue;

		in_intervals[0] = start_poses[x1];
		in_intervals[1] = start_poses[x1] + lens[x1];
		used_interval[x1]=1;
	
		int x2, this_sections = 1;
		for(x2 = x1 + 1; x2 < sections; x2++){
			if(strcmp( chros[x2], chros[x1] ) == 0){
				in_intervals[this_sections*2] = start_poses[x2];
				in_intervals[this_sections*2 + 1] = start_poses[x2] + lens[x2];
				used_interval[x2]=1;
				this_sections++;
			}
		}

		basic_sort( in_intervals, this_sections, overlap_compare, overlap_exchange );

		int merged_secs = mergeIntervals( in_intervals, out_intervals, this_sections );
		for(x2 = 0; x2 < merged_secs; x2++)
			ret += ( out_intervals[x2*2+1] - out_intervals[x2*2] );
	}
	return ret;
}


void vote_and_add_count(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context,
			srInt_64 * hits_indices1, int nhits1, srInt_64 * hits_indices2, int nhits2, unsigned int total_frag_len,
			char ** hits_chro1, char ** hits_chro2, unsigned int * hits_start_pos1, unsigned int * hits_start_pos2, unsigned short * hits_length1, unsigned short * hits_length2, int fixed_fractional_count, char * read_name, char * RG_name, char * bin1, char * bin2){
	if(global_context -> need_calculate_overlap_len == 0 && nhits2+nhits1==1) {
		srInt_64 hit_exon_id = nhits2?hits_indices2[0]:hits_indices1[0];

		//SUBREADprintf("V_AND_A: '%p'\n", RG_name);

		if(RG_name){
			void ** tab4s = get_RG_tables(global_context, thread_context, RG_name);
			fc_read_counters * sumtab = tab4s[1];
			sumtab -> assigned_reads++;
			
			read_count_type_t * count_table = tab4s[0];
			count_table[hit_exon_id] += fixed_fractional_count;				
		}else{
			thread_context->count_table[hit_exon_id] += fixed_fractional_count;
			thread_context->read_counters.assigned_reads ++;
		}
		thread_context->nreads_mapped_to_exon++;
		if(global_context -> read_details_out_FP){
			int final_gene_number = global_context -> exontable_geneid[hit_exon_id];
			char * final_feture_name = (char *)global_context -> gene_name_array[final_gene_number];
			write_read_details_FP(global_context, thread_context, "Assigned", 1, final_feture_name, bin1, bin2);
		}
	} else if(global_context -> need_calculate_overlap_len == 0 && nhits2 == 1 && nhits1 == 1 && hits_indices2[0]==hits_indices1[0]) {
		srInt_64 hit_exon_id = hits_indices1[0];
		
		if(RG_name){
			void ** tab4s = get_RG_tables(global_context, thread_context, RG_name);
			fc_read_counters * sumtab = tab4s[1];
			sumtab -> assigned_reads++;
			
			read_count_type_t * count_table = tab4s[0];
			count_table[hit_exon_id] += fixed_fractional_count;	
		}else{
			thread_context->count_table[hit_exon_id] += fixed_fractional_count;
			thread_context->read_counters.assigned_reads ++;
		}
		thread_context->nreads_mapped_to_exon++;
		if(global_context -> read_details_out_FP)
		{
			int final_gene_number = global_context -> exontable_geneid[hit_exon_id];
			char * final_feture_name = (char *)global_context -> gene_name_array[final_gene_number];
			write_read_details_FP(global_context, thread_context, "Assigned", 1, final_feture_name, bin1, bin2);
		}

	} else {
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

		
		unsigned int * scoring_numbers = thread_context -> scoring_buff_numbers;	// size is : MAX_HIT_NUMBER *2
		unsigned int * scoring_flags = thread_context -> scoring_buff_flags;		// size is : MAX_HIT_NUMBER *2
		unsigned int * scoring_overlappings = thread_context -> scoring_buff_overlappings;		// size is : MAX_HIT_NUMBER *2
		srInt_64 * scoring_exon_ids = thread_context -> scoring_buff_exon_ids;		// size is : MAX_HIT_NUMBER *2
		int scoring_count = 0,  score_x1, longest_overlap_score;

		if( global_context -> need_calculate_overlap_len ){
			int end1, end2, hit_x1, hit_x2;
			char ** scoring_gap_chros = thread_context -> scoring_buff_gap_chros;
			unsigned int * scoring_gap_starts = thread_context -> scoring_buff_gap_starts; // size is : MAX_HIT_NUMBER *2;
			unsigned short * scoring_gap_lengths = thread_context -> scoring_buff_gap_lengths; 	// size is : MAX_HIT_NUMBER *2*  global_context -> max_M*2
			longest_overlap_score = 0;

			char used_hit1 [nhits1];
			char used_hit2 [nhits2];

			if( global_context ->  fractional_minimum_feature_overlapping > 1E-10 || global_context -> max_missing_bases_in_feature >= 0){
				memset(used_hit1 , 0 , nhits1);
				memset(used_hit2 , 0 , nhits2);
				for(end1 = 0; end1 < global_context -> is_paired_end_mode_assign + 1 ; end1++){
					int allhits = end1?nhits2:nhits1;
					srInt_64 * hits_indices_X1 = end1?hits_indices2:hits_indices1;
					char * used_hit_X1 = end1?used_hit2:used_hit1;

					for(hit_x1 = 0; hit_x1 < allhits; hit_x1++){
						if(used_hit_X1[hit_x1])continue;

						srInt_64 tested_exon_id = hits_indices_X1[hit_x1];
						srInt_64 exon_span = global_context -> exontable_stop[tested_exon_id] +1;
						exon_span -= global_context -> exontable_start[tested_exon_id];

						srInt_64 applied_overlapping_threshold_frac = 0, applied_overlapping_threshold_missing = 0; 
						if(global_context -> max_missing_bases_in_feature >= 0){
							if(exon_span <= global_context -> max_missing_bases_in_feature) applied_overlapping_threshold_missing = 0;
							else applied_overlapping_threshold_missing = exon_span - global_context -> max_missing_bases_in_feature;
						}

						double applied_overlapping_threshold_frac_float = exon_span * global_context ->  fractional_minimum_feature_overlapping;
						applied_overlapping_threshold_frac = (srInt_64)(applied_overlapping_threshold_frac_float);
						// use the first 3 digits after decimal point to round-up: if the first 3 digits is 001, 002, ..., 999 then add 1.
						if(applied_overlapping_threshold_frac_float - applied_overlapping_threshold_frac >=0.001) applied_overlapping_threshold_frac++;

						srInt_64 applied_overlapping_threshold = max(applied_overlapping_threshold_frac , applied_overlapping_threshold_missing);

						scoring_gap_chros[0 ] = (end1?hits_chro2:hits_chro1)[hit_x1];
						scoring_gap_starts[0 ] = (end1?hits_start_pos2:hits_start_pos1)[hit_x1];
						scoring_gap_lengths[0 ] = (end1?hits_length2:hits_length1)[hit_x1];
						int gaps=1;

						for(end2 = 0; end2 < global_context -> is_paired_end_mode_assign + 1 ; end2++){
							int allhits2 = end2?nhits2:nhits1;
							char * used_hit_X2 = end2?used_hit2:used_hit1;
							srInt_64 * hits_indices_X2 = end2?hits_indices2:hits_indices1;


							for(hit_x2 = 0; hit_x2 < allhits2; hit_x2++){
								if(used_hit_X2[hit_x2]) continue;
								srInt_64 other_exon_id =  hits_indices_X2[hit_x2];
								if(other_exon_id == tested_exon_id){
									used_hit_X2[ hit_x2 ]=1;
									scoring_gap_chros[ gaps ] = (end2?hits_chro2:hits_chro1)[hit_x2];
									scoring_gap_starts[ gaps ] = (end2?hits_start_pos2:hits_start_pos1)[hit_x2];
									scoring_gap_lengths[ gaps ] = (end2?hits_length2:hits_length1)[hit_x2];
									gaps ++;
								}
							}
						}


						srInt_64 tested_exon_overlap_any_read = calc_score_overlaps(global_context, thread_context, scoring_gap_chros, scoring_gap_starts, scoring_gap_lengths, gaps, read_name);
						if(applied_overlapping_threshold > tested_exon_overlap_any_read){
							// remove this exon from lists
							for(end2 = 0; end2 < global_context -> is_paired_end_mode_assign + 1 ; end2++){
								int allhits2 = end2?nhits2:nhits1;
								srInt_64 * hits_indices_X2 = end2?hits_indices2:hits_indices1;

								for(hit_x2 = 0; hit_x2 < allhits2; hit_x2++){
									srInt_64 other_exon_id =  hits_indices_X2[hit_x2];
									if(other_exon_id == tested_exon_id){
										hits_indices_X2[hit_x2] = -1;
									}
								}
							}
						}
					}
				}
			}

			memset(used_hit1 , 0 , nhits1);
			memset(used_hit2 , 0 , nhits2);

			for(end1 = 0; end1 < global_context -> is_paired_end_mode_assign + 1 ; end1++){
				srInt_64 * hits_indices_X1 = end1?hits_indices2:hits_indices1;
				char * used_hit_X1 = end1?used_hit2:used_hit1;
				int nhit_X1 = end1?nhits2:nhits1;

				for( hit_x1 = 0 ; hit_x1 < nhit_X1; hit_x1 ++ ){
					if(used_hit_X1[hit_x1])continue;

					int gaps = 0;
					srInt_64 tmp_exon_id = hits_indices_X1[hit_x1];
					if(tmp_exon_id < 0) continue;
					srInt_64 score_merge_key;
					if (global_context -> is_gene_level )
						score_merge_key = global_context -> exontable_geneid[tmp_exon_id];
					else	score_merge_key = tmp_exon_id;


					scoring_gap_chros[0 ] = (end1?hits_chro2:hits_chro1)[hit_x1]; 
					scoring_gap_starts[0 ] = (end1?hits_start_pos2:hits_start_pos1)[hit_x1]; 
					scoring_gap_lengths[0 ] = (end1?hits_length2:hits_length1)[hit_x1]; 

					gaps=1;

					scoring_flags[scoring_count] = end1?2:1;
					scoring_numbers[scoring_count] =1;
					scoring_exon_ids[scoring_count] = tmp_exon_id;

					used_hit_X1[ hit_x1 ]=1;

					for(end2 = 0; end2 < global_context -> is_paired_end_mode_assign + 1 ; end2++){
						srInt_64 * hits_indices_X2 = end2?hits_indices2:hits_indices1;
						char * used_hit_X2 = end2?used_hit2:used_hit1;
						int nhit_X2 = end2?nhits2:nhits1;

						for( hit_x2 = 0 ; hit_x2 < nhit_X2; hit_x2 ++ ){
							if(used_hit_X2[hit_x2])continue;
							if(hits_indices_X2[hit_x2] < 0) continue;

							srInt_64 X2_merge_key;
							if (global_context -> is_gene_level )
								X2_merge_key = global_context -> exontable_geneid[ hits_indices_X2[hit_x2] ];
							else	X2_merge_key = hits_indices_X2[hit_x2];

							if( X2_merge_key == score_merge_key ){
								used_hit_X2[ hit_x2 ]=1;
								scoring_gap_chros[ gaps ] = (end2?hits_chro2:hits_chro1)[hit_x2]; 
								scoring_gap_starts[ gaps ] = (end2?hits_start_pos2:hits_start_pos1)[hit_x2]; 
								scoring_gap_lengths[ gaps ] = (end2?hits_length2:hits_length1)[hit_x2]; 

								if(!global_context -> is_multi_overlap_allowed)if((scoring_flags[scoring_count] & (end2?2:1))== 0 ){// when "-O" is specified: no preference toward two-end overlapping.
									scoring_flags[scoring_count] |= end2?2:1;
									scoring_numbers[scoring_count] ++;
								}
								gaps ++;
							}
						}
					}

					scoring_overlappings [scoring_count] = calc_score_overlaps(global_context, thread_context, scoring_gap_chros, scoring_gap_starts, scoring_gap_lengths, gaps, read_name);
					if( global_context -> only_want_largest_overlap ) {
						int this_overlap_score = scoring_overlappings [scoring_count] *2 + (scoring_flags[scoring_count] ==3?1:0); // after a discussion on 15FEB2023, we still prefer paired-end overlapping if largest-overlap is enabled.
						longest_overlap_score = max(longest_overlap_score,this_overlap_score);
						scoring_numbers[scoring_count] = this_overlap_score;
					}
					scoring_count++;
				}
			}
		}else{
			int ends;
			for(ends =0 ; ends < global_context -> is_paired_end_mode_assign + 1 ; ends++){
				int nhits = ends?nhits2:nhits1;
				srInt_64 * hits_indices = ends?hits_indices2:hits_indices1;

				int hit_x1;
				for(hit_x1 = 0; hit_x1 < nhits; hit_x1++){
					srInt_64 tmp_exon_id = hits_indices[hit_x1], score_merge_key;
					int found = 0;
					if (global_context -> is_gene_level )
						score_merge_key = global_context -> exontable_geneid[tmp_exon_id];
					else	score_merge_key = tmp_exon_id;

					for(score_x1 = 0; score_x1 < scoring_count; score_x1 ++){
						srInt_64 score_x1_key ; 
						if (global_context -> is_gene_level )
							score_x1_key = global_context -> exontable_geneid[ scoring_exon_ids[score_x1] ];
						else	score_x1_key = scoring_exon_ids[score_x1] ;

						if( score_x1_key == score_merge_key ){
							if(!global_context -> is_multi_overlap_allowed)if((scoring_flags[score_x1] & ( ends?2:1 )) == 0) {
								// when "-O" is specified: no preference toward two-end overlapping. SEE the Email on 15FEB2023
								scoring_flags[score_x1] |= (ends?2:1);
								scoring_numbers[score_x1] ++;
							}

							found = 1;
							break;
						}
					}

					if(0 == found){
						scoring_exon_ids[scoring_count] = tmp_exon_id;
						scoring_flags[scoring_count] = ends?2:1;
						scoring_numbers[scoring_count] = 1;

						scoring_count++;
					}
				}
			}
		}


		int passed_filters_score = 0;
		int passed_filters_count = 0;
		int passed_filters_first_index = 0;
		if(scoring_count == 0){
			if(global_context -> read_details_out_FP)
				write_read_details_FP(global_context, thread_context,"Unassigned_NoFeatures",-1, NULL, bin1, bin2);
			if(RG_name){
				void ** tab4s = get_RG_tables(global_context, thread_context, RG_name);
				fc_read_counters * sumtab = tab4s[1];
				sumtab -> unassigned_nofeatures++;
			}else thread_context->read_counters.unassigned_nofeatures ++;

		}else{
			if(global_context -> need_calculate_overlap_len) { 
				// Apply filters against overlap lengths. 
				srInt_64 applied_fragment_minimum_overlapping_overlap = 1, applied_fragment_minimum_overlapping_missing = 1;
				srInt_64 applied_fragment_minimum_overlapping = 1;

				if(global_context -> max_missing_bases_in_read >=0){
					if(total_frag_len <= global_context -> max_missing_bases_in_read) applied_fragment_minimum_overlapping_missing = 0;
					else applied_fragment_minimum_overlapping_missing = total_frag_len - global_context -> max_missing_bases_in_read;
				}

				double fractional_minimum_overlapping_float = global_context -> fractional_minimum_overlapping * total_frag_len;
				srInt_64 applied_fractional_minimum_overlapping = (srInt_64) fractional_minimum_overlapping_float;
				// use the first 3 digits after decimal point to round-up: if the first 3 digits is 001, 002, ..., 999 then add 1.
				if(fractional_minimum_overlapping_float - applied_fractional_minimum_overlapping >= 0.001) applied_fractional_minimum_overlapping++;

				applied_fragment_minimum_overlapping_overlap = max( global_context -> fragment_minimum_overlapping, applied_fractional_minimum_overlapping);
				applied_fragment_minimum_overlapping = max(applied_fragment_minimum_overlapping_overlap , applied_fragment_minimum_overlapping_missing);

				for(score_x1 = 0; score_x1 < scoring_count ; score_x1++){
					if(applied_fragment_minimum_overlapping > scoring_overlappings[score_x1]) scoring_numbers[score_x1] = 0;
					if(longest_overlap_score > scoring_numbers[score_x1]) scoring_numbers[score_x1] = 0;
						//"It is the largest overlap" is a filter step. Discussed on 15FEB2023.
						// longest_overlap_score is 0 if largestOverlap isn't enabled.
				}
			}

			for(score_x1 = 0; score_x1 < scoring_count ; score_x1++){
				if(scoring_numbers[score_x1]<1) continue;
				if(passed_filters_score < scoring_numbers[score_x1]){
					passed_filters_count = 1;
					passed_filters_score = scoring_numbers[score_x1];
					passed_filters_first_index = score_x1;
				}else if(passed_filters_score == scoring_numbers[score_x1])
					passed_filters_count++;
			}
		
			if(passed_filters_count == 0){
				if(global_context -> read_details_out_FP)
					write_read_details_FP(global_context, thread_context,"Unassigned_Overlapping_Length", -1, NULL, bin1, bin2);
				
				if(RG_name){
					void ** tab4s = get_RG_tables(global_context, thread_context, RG_name);
					fc_read_counters * sumtab = tab4s[1];
					sumtab -> unassigned_overlapping_length++;
				}else thread_context->read_counters.unassigned_overlapping_length ++;

			}else{
				if(1 == passed_filters_count && !global_context -> is_multi_overlap_allowed) {
					// some undesired scoring_numbers may be >0.
					// but we saved passed_filters_first_index so we just add score to it.
					srInt_64 max_exon_id = scoring_exon_ids[passed_filters_first_index];
					
					if(RG_name){
						void ** tab4s = get_RG_tables(global_context, thread_context, RG_name);
						fc_read_counters * sumtab = tab4s[1];
						sumtab -> assigned_reads++;
						
						read_count_type_t * count_table = tab4s[0];
						count_table[max_exon_id] += fixed_fractional_count;
					}else{
						thread_context->count_table[max_exon_id] += fixed_fractional_count;
						thread_context->read_counters.assigned_reads ++;
					}
					thread_context->nreads_mapped_to_exon++;
					if(global_context -> read_details_out_FP) {
						int final_gene_number = global_context -> exontable_geneid[max_exon_id];
						char * final_feture_name = (char *)global_context -> gene_name_array[final_gene_number];
						write_read_details_FP(global_context, thread_context,"Assigned", 1, final_feture_name, bin1, bin2);
					}

				}else if(global_context -> is_multi_overlap_allowed) {
					for(score_x1 = 0; score_x1 < scoring_count ; score_x1++)
						if(scoring_numbers[score_x1]<passed_filters_score) scoring_numbers[score_x1] = 0;
	
					#define GENE_NAME_LIST_BUFFER_SIZE (FEATURE_NAME_LENGTH * 50) 

					char final_feture_names[GENE_NAME_LIST_BUFFER_SIZE];
					int assigned_no = 0, xk1;
					int is_etc = 0;

					ArrayList * assigned_list = NULL;
					if(global_context -> read_details_out_FP) assigned_list = ArrayListCreate(nhits1+nhits2);
					for(xk1 = 0; xk1 < scoring_count; xk1++) {
						// This change was made on 31/MAR/2016
						if( scoring_numbers[xk1] < 1 ) continue ;

						srInt_64 tmp_voter_id = scoring_exon_ids[xk1];
						srInt_64 assignment_target_number = tmp_voter_id;
						if(global_context->is_gene_level) assignment_target_number = global_context -> exontable_geneid[tmp_voter_id];

						if(RG_name){
							void ** tab4s = get_RG_tables(global_context, thread_context, RG_name);
							read_count_type_t * count_table = tab4s[0];
							count_table[tmp_voter_id] += calculate_multi_overlap_fraction(global_context, fixed_fractional_count, passed_filters_count);
						}else thread_context->count_table[tmp_voter_id] += calculate_multi_overlap_fraction(global_context, fixed_fractional_count, passed_filters_count);

						if(global_context -> read_details_out_FP) {
							int final_gene_number = global_context -> exontable_geneid[tmp_voter_id];
							unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
							ArrayListPush(assigned_list, final_feture_name);
							assigned_no++;
						}
					}

					if(global_context -> read_details_out_FP){
						final_feture_names[0]=0;
						if(assigned_list -> numOfElements>1) ArrayListSort(assigned_list,ArrayListStringComparison); 
						long gene_cat_i;
						for(gene_cat_i=0; gene_cat_i<assigned_list -> numOfElements; gene_cat_i++){
							void * final_feture_name = ArrayListGet(assigned_list,gene_cat_i);
							if(strlen(final_feture_names)< (GENE_NAME_LIST_BUFFER_SIZE - 40 - FEATURE_NAME_LENGTH)) {
								strncat(final_feture_names, (char *)final_feture_name, GENE_NAME_LIST_BUFFER_SIZE-1);
								strncat(final_feture_names, ",", GENE_NAME_LIST_BUFFER_SIZE-1);
							}else{
								is_etc ++;
							}
						}
						ArrayListDestroy(assigned_list);

						if(is_etc) SUBreadSprintf(final_feture_names + strlen(final_feture_names), GENE_NAME_LIST_BUFFER_SIZE-strlen(final_feture_names), "... (%d names ommited),", is_etc);
						final_feture_names[GENE_NAME_LIST_BUFFER_SIZE-1]=0;
					}

					if(RG_name){
						void ** tab4s = get_RG_tables(global_context, thread_context, RG_name);
						fc_read_counters * sumtab = tab4s[1];
						sumtab -> assigned_reads++;
					}else{
						thread_context->read_counters.assigned_reads ++;
					}
					thread_context->nreads_mapped_to_exon++;
					
					if(global_context -> read_details_out_FP) {
						int ffnn = strlen(final_feture_names);
						if(ffnn>0) final_feture_names[ffnn-1]=0;
						// overlapped but still assigned 
						write_read_details_FP(global_context, thread_context, "Assigned", assigned_no, final_feture_names, bin1, bin2);
					}
				} else {
					if(global_context -> read_details_out_FP)
						write_read_details_FP(global_context, thread_context,"Unassigned_Ambiguity", -1, NULL, bin1, bin2);
					if(RG_name){
						fc_read_counters * sumtab = get_RG_tables(global_context, thread_context, RG_name)[1];
						sumtab -> unassigned_ambiguous++;
					}else{
						thread_context->read_counters.unassigned_ambiguous ++;
					}

				}
			}
		}
	}
}

// return the number of RG result sets
int fc_thread_merge_results(fc_thread_global_context_t * global_context, read_count_type_t * nreads , srInt_64 *nreads_mapped_to_exon, fc_read_counters * my_read_counter, HashTable * junction_global_table, HashTable * splicing_global_table, HashTable * RGmerged_table, fc_feature_info_t * loaded_features, srInt_64 nexons)
{
	int xk1, xk2, ret = 0, sample_i;

	srInt_64 total_input_reads = 0 ;
	(*nreads_mapped_to_exon)=0;
	SAM_pairer_destroy(&global_context -> read_pairer);

	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		if(global_context -> assign_reads_to_RG){
			HashTable * thread_rg_tab = global_context -> thread_contexts[xk1].RG_table;
			int buck_i;
			for(buck_i = 0; buck_i < thread_rg_tab -> numOfBuckets; buck_i++){
				KeyValuePair *cursor = thread_rg_tab -> bucketArray[buck_i];
				while(cursor){
					char * rg_name = (char *)cursor -> key;
					void ** rg_thread_tabs = cursor -> value;
					void ** rg_old_tabs = HashTableGet(RGmerged_table, rg_name);
					if(!rg_old_tabs){
						rg_old_tabs = malloc(sizeof(char *)*4); // all_counts, sum_counts , junc_table, split_table
						rg_old_tabs[0] = calloc(global_context -> thread_contexts[xk1].count_table_size, sizeof(srInt_64));
						rg_old_tabs[1] = calloc(1, sizeof(fc_read_counters));
						if(global_context -> do_junction_counting){
							HashTable * junction_counting_table = HashTableCreate(131317);
							HashTableSetHashFunction(junction_counting_table,HashTableStringHashFunction);
							HashTableSetDeallocationFunctions(junction_counting_table, free, NULL);
							HashTableSetKeyComparisonFunction(junction_counting_table, fc_strcmp_chro);
							
							HashTable * splicing_point_table = HashTableCreate(131317);
							HashTableSetHashFunction(splicing_point_table,HashTableStringHashFunction);
							HashTableSetDeallocationFunctions(splicing_point_table, free, NULL);
							HashTableSetKeyComparisonFunction(splicing_point_table, fc_strcmp_chro);

							rg_old_tabs[2] = junction_counting_table;
							rg_old_tabs[3] = splicing_point_table;
						}else rg_old_tabs[2] = NULL;
						
						HashTablePut(RGmerged_table, memstrcpy(rg_name), rg_old_tabs);
					}
					srInt_64 * rg_counts = rg_old_tabs[0];
					fc_read_counters * rg_sum_reads = rg_old_tabs[1];
					HashTable * rg_junc_tab = rg_old_tabs[2];
					HashTable * rg_split_tab = rg_old_tabs[3];
					
					srInt_64 * rg_thread_counts = rg_thread_tabs[0];
					fc_read_counters * rg_thread_sum_reads = rg_thread_tabs[1];
					HashTable * rg_thread_junc_table = rg_thread_tabs[2];
					HashTable * rg_thread_split_table = rg_thread_tabs[3];
					
					for(xk2=0; xk2<global_context -> exontable_exons; xk2++)
						rg_counts[xk2] += rg_thread_counts[xk2];
					
					rg_sum_reads->unassigned_ambiguous += rg_thread_sum_reads->unassigned_ambiguous;
					rg_sum_reads->unassigned_nofeatures += rg_thread_sum_reads->unassigned_nofeatures;
					rg_sum_reads->unassigned_overlapping_length += rg_thread_sum_reads->unassigned_overlapping_length;
					rg_sum_reads->unassigned_unmapped += rg_thread_sum_reads->unassigned_unmapped;
					rg_sum_reads->unassigned_singleton += rg_thread_sum_reads->unassigned_singleton;
					rg_sum_reads->unassigned_read_type += rg_thread_sum_reads->unassigned_read_type;
					rg_sum_reads->unassigned_mappingquality += rg_thread_sum_reads->unassigned_mappingquality;
					rg_sum_reads->unassigned_fragmentlength += rg_thread_sum_reads->unassigned_fragmentlength;
					rg_sum_reads->unassigned_chimericreads += rg_thread_sum_reads->unassigned_chimericreads;
					rg_sum_reads->unassigned_multimapping += rg_thread_sum_reads->unassigned_multimapping;
					rg_sum_reads->unassigned_secondary += rg_thread_sum_reads->unassigned_secondary;
					rg_sum_reads->unassigned_junction_condition += rg_thread_sum_reads->unassigned_junction_condition;
					rg_sum_reads->unassigned_duplicate += rg_thread_sum_reads->unassigned_duplicate;
					rg_sum_reads->assigned_reads += rg_thread_sum_reads->assigned_reads;

					if(global_context -> do_junction_counting){
						int bucket_i;
						for(bucket_i = 0 ; bucket_i < rg_thread_junc_table -> numOfBuckets; bucket_i++){
							KeyValuePair * cursor;
							cursor = rg_thread_junc_table -> bucketArray[bucket_i];
							while(cursor){
								char * junckey = (char *) cursor -> key;
								void * globval = HashTableGet(rg_junc_tab, junckey);
								char * new_key = memstrcpy(junckey);

								globval += (cursor -> value - NULL);
								HashTablePut(rg_junc_tab, new_key, globval);
									// new_key will be freed when it is replaced next time or when the global table is destroyed.

								cursor = cursor->next;
							}
						}

						for(bucket_i = 0 ; bucket_i < rg_thread_split_table -> numOfBuckets; bucket_i++){
							KeyValuePair * cursor;
							cursor = rg_thread_split_table -> bucketArray[bucket_i];
							while(cursor){
								char * junckey = (char *) cursor -> key;
								void * globval = HashTableGet(rg_split_tab, junckey);
								char * new_key = memstrcpy(junckey);

								//if(xk1>0)
								//SUBREADprintf("MERGE THREAD-%d : %s    VAL=%u, ADD=%u\n", xk1, junckey, globval - NULL, cursor -> value - NULL);
								globval += (cursor -> value - NULL);
								HashTablePut(rg_split_tab, new_key, globval);
								cursor = cursor->next;
							}
						}
					} // end : merge junc tables
					ret++;
					cursor = cursor -> next;
				}
			}
		}
		
		for(xk2=0; xk2<global_context -> exontable_exons; xk2++)
			nreads[xk2]+=global_context -> thread_contexts[xk1].count_table[xk2];

		total_input_reads += global_context -> thread_contexts[xk1].all_reads;
		(*nreads_mapped_to_exon) += global_context -> thread_contexts[xk1].nreads_mapped_to_exon;

		global_context -> read_counters.unassigned_ambiguous += global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous;
		global_context -> read_counters.unassigned_nofeatures += global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures;
		global_context -> read_counters.unassigned_overlapping_length += global_context -> thread_contexts[xk1].read_counters.unassigned_overlapping_length;
		global_context -> read_counters.unassigned_unmapped += global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped;
		global_context -> read_counters.unassigned_singleton += global_context -> thread_contexts[xk1].read_counters.unassigned_singleton;
		global_context -> read_counters.unassigned_read_type += global_context -> thread_contexts[xk1].read_counters.unassigned_read_type;
		global_context -> read_counters.unassigned_mappingquality += global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality;
		global_context -> read_counters.unassigned_fragmentlength += global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength;
		global_context -> read_counters.unassigned_chimericreads += global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads;
		global_context -> read_counters.unassigned_multimapping += global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping;
		global_context -> read_counters.unassigned_secondary += global_context -> thread_contexts[xk1].read_counters.unassigned_secondary;
		global_context -> read_counters.unassigned_junction_condition += global_context -> thread_contexts[xk1].read_counters.unassigned_junction_condition;
		global_context -> read_counters.unassigned_duplicate += global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate;
		global_context -> read_counters.assigned_reads += global_context -> thread_contexts[xk1].read_counters.assigned_reads;

		my_read_counter->unassigned_ambiguous += global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous;
		my_read_counter->unassigned_nofeatures += global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures;
		my_read_counter->unassigned_overlapping_length += global_context -> thread_contexts[xk1].read_counters.unassigned_overlapping_length;
		my_read_counter->unassigned_unmapped += global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped;
		my_read_counter->unassigned_singleton += global_context -> thread_contexts[xk1].read_counters.unassigned_singleton;
		my_read_counter->unassigned_read_type += global_context -> thread_contexts[xk1].read_counters.unassigned_read_type;
		my_read_counter->unassigned_mappingquality += global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality;
		my_read_counter->unassigned_fragmentlength += global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength;
		my_read_counter->unassigned_chimericreads += global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads;
		my_read_counter->unassigned_multimapping += global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping;
		my_read_counter->unassigned_secondary += global_context -> thread_contexts[xk1].read_counters.unassigned_secondary;
		my_read_counter->unassigned_junction_condition += global_context -> thread_contexts[xk1].read_counters.unassigned_junction_condition;
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



	if(0 == global_context -> is_input_bad_format){

		if(global_context -> is_paired_end_reads_expected){
			if(global_context -> is_mixed_PE_SE)
					print_in_box(80,0,0,"   WARNING: Single-end reads were found%s.", global_context -> is_strand_checked?" and excluded":"");
			else print_in_box(80,0,0,"   Paired-end reads are included.");
			if(!global_context -> is_paired_end_mode_assign)
				print_in_box(80,0,0, "   The reads are assigned on the single-end mode.");
		}else{
			// paired-end reads in a single-end lib will result in error.
			print_in_box(80,0,0,"   Single-end reads are included.");
		}

		char pct_str[10];
		if(total_input_reads>0)
			SUBreadSprintf(pct_str, 10,"(%.1f%%%%)", (*nreads_mapped_to_exon)*100./total_input_reads);
		else	pct_str[0]=0;
	
		int show_summary = 1;
		if(global_context -> assign_reads_to_RG){
			if(RGmerged_table -> numOfElements)
				print_in_box(80,0,0,"   Total read groups : %ld", RGmerged_table -> numOfElements);
			else{
				print_in_box(80,0,0,"   No read groups are found; no output is generated.");
				show_summary = 0;
			}
		}
		if(show_summary){
			print_in_box(80,0,0,"   Total alignments : %llu", total_input_reads); 
			print_in_box(pct_str[0]?81:80,0,0,"   Successfully assigned alignments : %llu %s", *nreads_mapped_to_exon,pct_str); 
		}
		print_in_box(80,0,0,"   Running time : %.2f minutes", (miltime() - global_context -> start_time)/60);
		print_in_box(80,0,0,"");
	}
	return ret;
}

void get_temp_dir_from_out(char * tmp, char * out){
	char * slash = strrchr(out,'/');
	if(NULL == slash){
		strcpy(tmp, "./");
	}else{
		memcpy(tmp, out, slash - out);
		tmp[slash - out]=0;
	}
}

void fc_thread_init_input_files(fc_thread_global_context_t * global_context, char * in_fnames, char ** out_ptr ){
	if(global_context -> use_stdin_file){
		#ifdef MAKE_STANDALONE

		char MAC_or_random[13];

		(*out_ptr) = malloc(MAX_FILE_NAME_LENGTH);
		mac_or_rand_str(MAC_or_random);
		SUBreadSprintf(*out_ptr, MAX_FILE_NAME_LENGTH, "%s/temp-core-%06u-%s.sam", global_context -> temp_file_dir, getpid(), MAC_or_random);

		SUBREADprintf("\nReading data from <STDIN> for featureCounts ...\n\n");
		
		FILE * ifp = fopen(*out_ptr,"w");
		while(1){
			char nchar[100];
			int rlen = fread(nchar, 1, 100, stdin);
			if(rlen > 0) fwrite(nchar, 1, rlen, ifp);
			else break;
			//if(rlen < 100)break;
		}
		fclose(ifp);

		#endif
	}else{
		(*out_ptr) = malloc(strlen(in_fnames)+1);
		strcpy((*out_ptr), in_fnames);
	}

}

void fc_NCfree(void * vv){
	char ** cc = vv;
	int i;
	for(i=0; cc[i]; i++) free(cc[i]);
	free(vv);
}


void fc_thread_init_global_context(fc_thread_global_context_t * global_context, unsigned int buffer_size, unsigned short threads, int line_length, int min_pe_dist, int max_pe_dist, int is_gene_level, int is_overlap_allowed, char * strand_check_mode, char * output_fname, int is_sam_out, int is_both_end_required, int is_chimertc_disallowed, int is_PE_distance_checked, char *feature_name_column, char * gene_id_column, char * transcript_id_column, int min_map_qual_score, int is_multi_mapping_allowed, int is_SAM, char * alias_file_name, char * cmd_rebuilt, int is_input_file_resort_needed, int feature_block_size, int isCVersion, int fiveEndExtension,  int threeEndExtension, int minFragmentOverlap, int is_split_or_exonic_only, int reduce_5_3_ends_to_one, char * debug_command, int is_duplicate_ignored, int is_not_sort, int use_fraction_multimapping, int useOverlappingBreakTie, char * pair_orientations, int do_junction_cnt, int max_M, int isRestrictlyNoOvelrapping, float fracOverlap, char * temp_dir, int use_stdin_file, int assign_reads_to_RG, int long_read_minimum_length, int is_verbose, float frac_feature_overlap, int do_detection_call, int max_missing_bases_in_read, int max_missing_bases_in_feature, int is_primary_alignment_only, char * Rpath, char * extra_column_names , char * annotation_file_screen_output, int read_shift_type, int read_shift_size, int is_dual_index ) {
	int x1;
	myrand_srand(time(NULL));

	memset(global_context, 0, sizeof(fc_thread_global_context_t));
	global_context -> max_BAM_header_size = buffer_size;
	global_context -> all_reads = 0;
	global_context -> redo = 0;
	global_context -> read_details_out_FP = NULL;

	global_context -> reported_extra_columns = extra_column_names;
	global_context -> isCVersion = isCVersion;
	global_context -> is_read_details_out = is_sam_out;
	global_context -> is_multi_overlap_allowed = is_overlap_allowed;
	global_context -> restricted_no_multi_overlap = isRestrictlyNoOvelrapping;
	global_context -> is_gene_level = is_gene_level;
	global_context -> strand_check_mode = strand_check_mode;
	global_context -> is_both_end_required = is_both_end_required;
	global_context -> is_chimertc_disallowed = is_chimertc_disallowed;
	global_context -> is_PE_distance_checked = is_PE_distance_checked;
	global_context -> is_multi_mapping_allowed = is_multi_mapping_allowed;
	global_context -> is_primary_alignment_only = is_primary_alignment_only;
	global_context -> is_split_or_exonic_only = is_split_or_exonic_only;
	global_context -> is_duplicate_ignored = is_duplicate_ignored;
	global_context -> use_stdin_file = use_stdin_file;
	global_context -> assign_reads_to_RG = assign_reads_to_RG;
	global_context -> long_read_minimum_length = long_read_minimum_length;
	global_context -> is_verbose = is_verbose;
	global_context -> do_detection_call = do_detection_call;
	//global_context -> is_first_read_reversed = (pair_orientations[0]=='r');
	//global_context -> is_second_read_straight = (pair_orientations[1]=='f');

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
	global_context -> BAM_chros_to_anno_table = NULL;
	global_context -> cmd_rebuilt = cmd_rebuilt;
	global_context -> feature_block_size = feature_block_size;
	global_context -> five_end_extension = fiveEndExtension;
	global_context -> three_end_extension = threeEndExtension;
	global_context -> read_shift_type = read_shift_type;
	global_context -> read_shift_size = read_shift_size;
	global_context -> fragment_minimum_overlapping = minFragmentOverlap;
	global_context -> fractional_minimum_overlapping = fracOverlap;
	global_context -> fractional_minimum_feature_overlapping = frac_feature_overlap;
	global_context -> max_missing_bases_in_read = max_missing_bases_in_read;
	global_context -> max_missing_bases_in_feature = max_missing_bases_in_feature;
	global_context -> only_want_largest_overlap = useOverlappingBreakTie;
	global_context -> need_calculate_fragment_len = ( global_context -> fractional_minimum_overlapping > 1E-10 ) || (global_context -> fractional_minimum_feature_overlapping > 1E-10) || ( global_context -> max_missing_bases_in_read >= 0 ) || ( global_context -> max_missing_bases_in_feature >= 0 );
	global_context -> need_calculate_overlap_len = (global_context -> fractional_minimum_overlapping > 1E-10) || (global_context -> fragment_minimum_overlapping > 1) || global_context -> only_want_largest_overlap || (global_context -> fractional_minimum_feature_overlapping > 1E-10) || ( global_context -> max_missing_bases_in_read >= 0 ) || ( global_context -> max_missing_bases_in_feature >= 0 );
	global_context -> debug_command = debug_command;
	global_context -> max_M = max_M;
	global_context -> max_BAM_header_size = buffer_size;

	global_context -> read_counters.unassigned_ambiguous=0;
	global_context -> read_counters.unassigned_nofeatures=0;
	global_context -> read_counters.unassigned_overlapping_length=0;
	global_context -> read_counters.unassigned_unmapped=0;
	global_context -> read_counters.unassigned_read_type=0;
	global_context -> read_counters.unassigned_singleton=0;
	global_context -> read_counters.unassigned_mappingquality=0;
	global_context -> read_counters.unassigned_fragmentlength=0;
	global_context -> read_counters.unassigned_chimericreads=0;
	global_context -> read_counters.unassigned_multimapping=0;
	global_context -> read_counters.unassigned_secondary=0;
	global_context -> read_counters.unassigned_junction_condition=0;
	global_context -> read_counters.unassigned_duplicate=0;
	global_context -> read_counters.assigned_reads=0;

	global_context -> GCcontent_table = HashTableCreate(20000);
	HashTableSetHashFunction(global_context -> GCcontent_table, HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(global_context -> GCcontent_table, free, free);
	HashTableSetKeyComparisonFunction(global_context -> GCcontent_table, fc_strcmp_chro);

	if(annotation_file_screen_output) strcpy(global_context -> annotation_file_screen_output, annotation_file_screen_output);
	else global_context ->annotation_file_screen_output[0]=0;

	if(alias_file_name && alias_file_name[0])
	{
		strcpy(global_context -> alias_file_name,alias_file_name);
		global_context -> BAM_chros_to_anno_table = load_alias_table(alias_file_name);
	}
	else	global_context -> alias_file_name[0]=0;

	global_context -> read_details_path[0]=0;
	if(Rpath)strcpy(global_context -> read_details_path, Rpath);

	strcpy(global_context -> feature_name_column,feature_name_column);
	strcpy(global_context -> gene_id_column,gene_id_column);
	strcpy(global_context -> transcript_id_column,transcript_id_column);
	strcpy(global_context -> output_file_name, output_fname);
	global_context -> output_file_path[0]=0;
	for( x1 = strlen(output_fname)-1; x1 >= 0; x1 --){
		if(output_fname[x1]=='/'){
			memcpy(global_context -> output_file_path, output_fname, x1);
			global_context -> output_file_path[x1]=0;
			break;
		}
	}
	if(0 == global_context -> output_file_path[0]){
		strcpy(global_context -> output_file_path, ".");
	}

	if(temp_dir == NULL)get_temp_dir_from_out(global_context -> temp_file_dir, output_fname);
	else strcpy(global_context -> temp_file_dir, temp_dir);
	//SUBREADprintf("OFPP:%s, OFNN:%s\n", global_context -> output_file_path, global_context -> output_file_name);

	global_context -> min_paired_end_distance = min_pe_dist;
	global_context -> max_paired_end_distance = max_pe_dist;
	global_context -> thread_number = threads;
	global_context -> line_length = line_length;
}


void fc_thread_start_open_report_fp(fc_thread_global_context_t * global_context){
	if(global_context -> read_details_out_FP) fclose( global_context -> read_details_out_FP );
	if(global_context -> is_read_details_out){
		char tmp_fname[MAX_FILE_NAME_LENGTH+20], *modified_fname;
		int i=0;
		char * applied_detail_path = global_context -> output_file_path;
		if(global_context -> read_details_path[0]) applied_detail_path = global_context -> read_details_path;

		if( global_context -> input_file_unique ){
			SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+20, "%s/%s.featureCounts%s", applied_detail_path, global_context -> input_file_short_name, global_context -> is_read_details_out == FILE_TYPE_BAM?".bam":(global_context -> is_read_details_out == FILE_TYPE_SAM?".sam":""));
			global_context -> read_details_out_FP = f_subr_open(tmp_fname, "w");
			//SUBREADprintf("FCSSF=%s\n", tmp_fname);
		} else {
			SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+20, "%s.featureCounts%s", global_context -> raw_input_file_name, global_context -> is_read_details_out == FILE_TYPE_BAM?".bam":(global_context -> is_read_details_out == FILE_TYPE_SAM?".sam":""));
			modified_fname = tmp_fname;
			while(modified_fname[0]=='/' || modified_fname[0]=='.' || modified_fname[0]=='\\'){
				modified_fname ++;
			}
			while(modified_fname[i]){
				if(modified_fname[i]=='\\' || modified_fname[i]=='/'||modified_fname[i]==' ')modified_fname[i]='.';
				i++;
			}
			char tmp_fname2[MAX_FILE_NAME_LENGTH*2+100];
			SUBreadSprintf(tmp_fname2, MAX_FILE_NAME_LENGTH*2+100, "%s/%s", applied_detail_path, modified_fname);
			global_context -> read_details_out_FP = f_subr_open(tmp_fname2, "w");
			//SUBREADprintf("FCSSF=%s\n", tmp_fname2);
		}
		if(global_context -> read_details_out_FP){
			pthread_spin_init(&global_context -> read_details_out_lock, 1);
		}else{
			SUBREADprintf("Unable to create file '%s'; the read assignment details are not written.\n", tmp_fname);
		}
	}else global_context -> read_details_out_FP = NULL;
}

int fc_thread_start_threads(fc_thread_global_context_t * global_context, int et_exons, int * et_geneid, char ** et_chr, srInt_64 * et_start, srInt_64 * et_stop, unsigned char * et_strand, char * et_anno_chr_2ch, char ** et_anno_chrs, srInt_64 * et_anno_chr_heads, srInt_64 * et_bk_end_index, srInt_64 * et_bk_min_start, srInt_64 * et_bk_max_end, int read_length)
{
	int xk1;

	global_context -> read_length = read_length;
	global_context -> is_unpaired_warning_shown = 0;
	global_context -> is_stake_warning_shown = 0;
	global_context -> is_read_too_long_to_SAM_BAM_shown = 0;

	fc_thread_start_open_report_fp(global_context);

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
	global_context -> sambam_chro_table_items = 0;
	global_context -> sambam_chro_table = NULL;


	global_context -> thread_contexts = malloc(sizeof(fc_thread_thread_context_t) * global_context -> thread_number);
	for(xk1=0; xk1<global_context -> thread_number; xk1++)
	{
	//	printf("CHRR_MALLOC\n");
		global_context -> thread_contexts[xk1].thread_id = xk1;
		global_context -> thread_contexts[xk1].chunk_read_ptr = 0;
		global_context -> thread_contexts[xk1].count_table = calloc(sizeof(read_count_type_t), et_exons);
		global_context -> thread_contexts[xk1].count_table_size = et_exons;
		global_context -> thread_contexts[xk1].nreads_mapped_to_exon = 0;
		global_context -> thread_contexts[xk1].all_reads = 0;
		global_context -> thread_contexts[xk1].chro_name_buff = malloc(CHROMOSOME_NAME_LENGTH);

		global_context -> thread_contexts[xk1].read_counters.assigned_reads = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_singleton = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_read_type = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_secondary = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_junction_condition = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_overlapping_length = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate = 0;
		global_context -> thread_contexts[xk1].read_details_buff_used = 0;
		global_context -> thread_contexts[xk1].hits_number_capacity = 300 ;

		global_context -> thread_contexts[xk1].hits_start_pos1 = malloc(sizeof(int)* global_context -> thread_contexts[xk1].hits_number_capacity);
		global_context -> thread_contexts[xk1].hits_start_pos2 = malloc(sizeof(int)* global_context -> thread_contexts[xk1].hits_number_capacity);
		global_context -> thread_contexts[xk1].hits_length1 = malloc(sizeof(short)* global_context -> thread_contexts[xk1].hits_number_capacity);
		global_context -> thread_contexts[xk1].hits_length2 = malloc(sizeof(short)* global_context -> thread_contexts[xk1].hits_number_capacity);
		global_context -> thread_contexts[xk1].hits_chro1 = malloc(sizeof(char*)* global_context -> thread_contexts[xk1].hits_number_capacity);
		global_context -> thread_contexts[xk1].hits_chro2 = malloc(sizeof(char*)* global_context -> thread_contexts[xk1].hits_number_capacity);
		global_context -> thread_contexts[xk1].hits_indices1 = malloc(sizeof(srInt_64)* global_context -> thread_contexts[xk1].hits_number_capacity);
		global_context -> thread_contexts[xk1].hits_indices2 = malloc(sizeof(srInt_64)* global_context -> thread_contexts[xk1].hits_number_capacity);

		global_context -> thread_contexts[xk1].scoring_buff_numbers = malloc(sizeof(int)* global_context -> thread_contexts[xk1].hits_number_capacity * 2);
		global_context -> thread_contexts[xk1].scoring_buff_flags = malloc(sizeof(int)* global_context -> thread_contexts[xk1].hits_number_capacity * 2);
		global_context -> thread_contexts[xk1].scoring_buff_overlappings = malloc(sizeof(int)* global_context -> thread_contexts[xk1].hits_number_capacity * 2);
		global_context -> thread_contexts[xk1].scoring_buff_exon_ids =malloc(sizeof(srInt_64)* global_context -> thread_contexts[xk1].hits_number_capacity * 2);

		if(global_context -> read_details_out_FP){
			global_context -> thread_contexts[xk1].read_details_buff = malloc(70000 + 2 * MAX_FC_READ_LENGTH * 3);
			global_context -> thread_contexts[xk1].bam_compressed_buff = malloc(70000 + 2 * MAX_FC_READ_LENGTH * 3);
		}

		if(global_context -> need_calculate_overlap_len){
			global_context -> thread_contexts[xk1].scoring_buff_gap_chros = malloc( sizeof(char *) * global_context -> thread_contexts[xk1].hits_number_capacity * 2 * global_context -> max_M *2);
			global_context -> thread_contexts[xk1].scoring_buff_gap_starts = malloc( sizeof(unsigned int ) * global_context -> thread_contexts[xk1].hits_number_capacity * 2 * global_context -> max_M *2);
			global_context -> thread_contexts[xk1].scoring_buff_gap_lengths = malloc( sizeof(unsigned short) * global_context -> thread_contexts[xk1].hits_number_capacity * 2 * global_context -> max_M *2);
		} else global_context -> thread_contexts[xk1].scoring_buff_gap_chros = NULL;

		if(global_context -> do_junction_counting){
			global_context -> thread_contexts[xk1].junction_counting_table = HashTableCreate(131317);
			HashTableSetHashFunction(global_context -> thread_contexts[xk1].junction_counting_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(global_context -> thread_contexts[xk1].junction_counting_table, free, NULL);
			HashTableSetKeyComparisonFunction(global_context -> thread_contexts[xk1].junction_counting_table, fc_strcmp_chro);
			
			global_context -> thread_contexts[xk1].splicing_point_table = HashTableCreate(131317);
			HashTableSetHashFunction(global_context -> thread_contexts[xk1].splicing_point_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(global_context -> thread_contexts[xk1].splicing_point_table, free, NULL);
			HashTableSetKeyComparisonFunction(global_context -> thread_contexts[xk1].splicing_point_table, fc_strcmp_chro);

			global_context ->  thread_contexts[xk1].memory_supported_junctions1 = malloc(sizeof(fc_junction_info_t) * 65536);
			global_context ->  thread_contexts[xk1].memory_supported_junctions2 = malloc(sizeof(fc_junction_info_t) * 65536);
		}
		
		if(global_context -> assign_reads_to_RG){
			global_context -> thread_contexts[xk1].RG_table = HashTableCreate(97);
			HashTableSetHashFunction(global_context -> thread_contexts[xk1].RG_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(global_context -> thread_contexts[xk1].RG_table, free, disallocate_RG_tables);
			HashTableSetKeyComparisonFunction(global_context -> thread_contexts[xk1].RG_table, fc_strcmp_chro);
		}

		if(!global_context ->  thread_contexts[xk1].count_table) return 1;
	}

	char new_fn[MAX_FILE_NAME_LENGTH+10];
	char MAC_or_random[13];
	mac_or_rand_str(MAC_or_random);
	char rand_prefix[MAX_FILE_NAME_LENGTH+100];
	SUBreadSprintf(rand_prefix, MAX_FILE_NAME_LENGTH+100, "%s/temp-core-%06u-%s.sam", global_context -> temp_file_dir, getpid(), MAC_or_random);
	if(global_context -> use_stdin_file) SUBreadSprintf(new_fn, MAX_FILE_NAME_LENGTH+10, "<%s",  global_context -> input_file_name );
	else SUBreadSprintf(new_fn, MAX_FILE_NAME_LENGTH+10, "%s",  global_context -> input_file_name );

	int pairer_error = SAM_pairer_create(&global_context -> read_pairer, global_context -> thread_number , global_context -> max_BAM_header_size/1024/1024+2, !global_context-> is_SAM_file, !( global_context -> is_read_details_out == FILE_TYPE_BAM ||global_context -> is_read_details_out == FILE_TYPE_SAM ) , !global_context -> is_paired_end_mode_assign, global_context ->is_paired_end_mode_assign && global_context -> do_not_sort, global_context -> assign_reads_to_RG ,0, new_fn, process_pairer_reset, process_pairer_header, process_pairer_output, rand_prefix, global_context,  global_context -> long_read_minimum_length);

	return pairer_error;
}

void fc_thread_destroy_thread_context(fc_thread_global_context_t * global_context)
{
	int xk1;

	if(global_context -> is_read_details_out)for(xk1=0; xk1<global_context-> thread_number; xk1++)
		write_read_detailed_remainder(global_context, global_context -> thread_contexts+xk1);

	if(global_context -> is_read_details_out) {
		if( global_context -> is_read_details_out == FILE_TYPE_BAM ){
			char bam_tail_block[1000];
			int tail_size = compress_read_detail_BAM( global_context, global_context -> thread_contexts, 0,0,bam_tail_block);
			assert(tail_size > 0);
			//SUBREADprintf("TAIL SIZE=%d\n", tail_size);
			fwrite(bam_tail_block, 1, tail_size, global_context -> read_details_out_FP);
		}
		fclose(global_context -> read_details_out_FP);
		global_context -> read_details_out_FP = NULL;
		pthread_spin_destroy(&global_context -> read_details_out_lock);
	}

	for(xk1=0; xk1<global_context-> thread_number; xk1++) {
		//printf("CHRR_FREE\n");
		free(global_context -> thread_contexts[xk1].count_table);	
		free(global_context -> thread_contexts[xk1].chro_name_buff);
		free(global_context -> thread_contexts[xk1].hits_start_pos1);
		free(global_context -> thread_contexts[xk1].hits_start_pos2);
		free(global_context -> thread_contexts[xk1].hits_length1);
		free(global_context -> thread_contexts[xk1].hits_length2);
		free(global_context -> thread_contexts[xk1].hits_chro1);
		free(global_context -> thread_contexts[xk1].hits_chro2);
		free(global_context -> thread_contexts[xk1].hits_indices1);
		free(global_context -> thread_contexts[xk1].hits_indices2);
		free(global_context -> thread_contexts[xk1].scoring_buff_numbers);
		free(global_context -> thread_contexts[xk1].scoring_buff_flags);
		free(global_context -> thread_contexts[xk1].scoring_buff_overlappings);
		free(global_context -> thread_contexts[xk1].scoring_buff_exon_ids);

		if(global_context -> thread_contexts[xk1].scoring_buff_gap_chros){
			free(global_context -> thread_contexts[xk1].scoring_buff_gap_chros);
			free(global_context -> thread_contexts[xk1].scoring_buff_gap_starts);
			free(global_context -> thread_contexts[xk1].scoring_buff_gap_lengths);
		}
		if(global_context -> do_junction_counting){
			HashTableDestroy(global_context -> thread_contexts[xk1].junction_counting_table);
			HashTableDestroy(global_context -> thread_contexts[xk1].splicing_point_table);
			free(global_context ->  thread_contexts[xk1].memory_supported_junctions1);
			free(global_context ->  thread_contexts[xk1].memory_supported_junctions2);
		}
		if(global_context -> assign_reads_to_RG)
			HashTableDestroy(global_context -> thread_contexts[xk1].RG_table);
		if(global_context -> is_read_details_out ){
			free(global_context -> thread_contexts[xk1].read_details_buff);
			free(global_context -> thread_contexts[xk1].bam_compressed_buff);
		}

	}

	free(global_context -> thread_contexts);
}
void fc_thread_wait_threads(fc_thread_global_context_t * global_context)
{
	int assign_ret = SAM_pairer_run(&global_context -> read_pairer);
	if(0 && assign_ret){
		print_in_box(80,0,0,"");
		print_in_box(80,0,0,"   format error found in this file.");
	}
	global_context -> is_input_bad_format |= assign_ret;
}

void merge_repeated_extra_columns(char * cols){
	if(cols[0]!=';')return;

	int is_diff = 0;
	int seglen = -1, laststart = 0;
	int xx;
	for(xx=0; ; xx++){
		if(cols[xx]==';' || cols[xx]==0){
			if(seglen <0)seglen = xx -1;
			else{
				is_diff = (xx-laststart != seglen )|| memcmp(cols+laststart, cols+1, seglen);
				if(is_diff)break;
			}
			laststart = xx+1;
		}
		if(cols[xx]==0)break;
	}

	if(seglen>0 && !is_diff) cols[seglen+1]=0;
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

void fc_write_final_gene_results(fc_thread_global_context_t * global_context, int * et_geneid, char ** et_chr, srInt_64 * et_start, srInt_64 * et_stop, unsigned char * et_strand, char ** et_extra_columns, const char * out_file, int features, ArrayList * column_numbers, ArrayList * column_names, fc_feature_info_t * loaded_features, int header_out)
{
	int xk1,xk4;
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

	int i_files;
	fprintf(fp_out,"Geneid\t%sChr\tStart\tEnd\tStrand\tLength%s%s", global_context->do_detection_call?"GCfraction\t":"", global_context -> reported_extra_columns?"\t":"", global_context -> reported_extra_columns?global_context -> reported_extra_columns:"");
	for(i_files=0; i_files<column_names->numOfElements; i_files++)
	{
		char * next_fn = ArrayListGet(column_names, i_files);
		fprintf(fp_out,"\t%s", global_context -> use_stdin_file?"STDIN":next_fn);
	}
	
	fprintf(fp_out,"\n");

	gene_columns = calloc(sizeof(read_count_type_t) , genes * column_names->numOfElements);
	unsigned int * gene_exons_number = calloc(sizeof(unsigned int) , genes);
	unsigned int * gene_exons_pointer = calloc(sizeof(unsigned int) , genes);
	unsigned int * gene_exons_start = malloc(sizeof(unsigned int) * features);
	unsigned int * gene_exons_end = malloc(sizeof(unsigned int) * features);
	char ** gene_exons_chr = malloc(sizeof(char *) * features);
	char ** gene_exons_extra_columns = malloc(sizeof(char *) * features);
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
		unsigned int this_gene_exons = gene_exons_number[xk1];
		longest_gene_exons = max(longest_gene_exons, this_gene_exons);
		gene_exons_number[xk1] = accumulative_no;
		accumulative_no += this_gene_exons;
	}

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1];
		int gene_write_ptr = gene_exons_number[gene_id] + gene_exons_pointer[gene_id];

		gene_exons_chr[gene_write_ptr] = et_chr[xk1];
		gene_exons_start[gene_write_ptr] = et_start[xk1]; 
		gene_exons_end[gene_write_ptr] = et_stop[xk1]; 
		gene_exons_strand[gene_write_ptr] = et_strand[xk1]; 
		if(global_context -> reported_extra_columns!=NULL)gene_exons_extra_columns[gene_write_ptr] = et_extra_columns[xk1];

		gene_exons_pointer[gene_id]++;
	}

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1], k_noempty = 0;
		for(i_files=0;i_files < column_names->numOfElements; i_files++)
		{
			srInt_64 * this_col = ArrayListGet(column_numbers, i_files);
			gene_columns[gene_id * column_names->numOfElements + k_noempty ] += this_col[xk1];
			k_noempty++;
		}
	}


	char *is_occupied = malloc(longest_gene_exons);
	unsigned int * input_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);
	unsigned int * output_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);
	int disk_is_full = 0;

	char * out_chr_list = malloc(longest_gene_exons * (1+global_context -> longest_chro_name) + 1), * tmp_chr_list = NULL;
	char * out_start_list = malloc(11 * longest_gene_exons + 1), * tmp_start_list = NULL;
	char * out_end_list = malloc(11 * longest_gene_exons + 1), * tmp_end_list = NULL;
	char * out_strand_list = malloc(2 * longest_gene_exons + 1), * tmp_strand_list = NULL;

	char * out_extra_columns[MAX_EXTRA_COLS];
	int out_extra_column_size[MAX_EXTRA_COLS];
	int total_extra_cols = 0;
	if(global_context -> reported_extra_columns){
		char * tnamep = global_context -> reported_extra_columns;
		total_extra_cols =1;
		while(*(tnamep++))
			total_extra_cols += '\t' ==(*tnamep);
		for(xk1=0; xk1<total_extra_cols; xk1++){
			out_extra_columns[xk1] = malloc(220);
			out_extra_column_size[xk1] = 220;
		}
	}


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
		for(xk4=0; xk4<total_extra_cols; xk4++)
			out_extra_columns[xk4][0]=0;
		int gene_nonoverlap_len =0;

		unsigned char * gene_symbol = global_context -> gene_name_array [xk1];
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

				for(xk3 = xk2; xk3 < gene_exons_pointer[xk1]; xk3++)
				{
					if( global_context -> reported_extra_columns &&  (xk3==xk2 || (0 == is_occupied[xk3] && strcmp(matched_chr, gene_exons_chr[xk3+gene_exons_number[xk1]])==0 && matched_strand == gene_exons_strand[xk3 + gene_exons_number[xk1]] ))){
						char * this_col_ptr = NULL;
						char * this_col = strtok_r(gene_exons_extra_columns[xk3+gene_exons_number[xk1]], "\t", &this_col_ptr);
						for(xk4 = 0; xk4 < total_extra_cols; xk4++){
							int exlen = strlen( this_col), ollen = strlen(out_extra_columns[xk4]);
							if(ollen + exlen +2 > out_extra_column_size[xk4]){
								out_extra_column_size[xk4] = max(ollen + exlen +2, out_extra_column_size[xk4]);
								out_extra_columns[xk4] = realloc(out_extra_columns[xk4], out_extra_column_size[xk4]);
							}
							SUBreadSprintf(out_extra_columns[xk4]+ollen, out_extra_column_size[xk4]-ollen,";%s", this_col);
							this_col = strtok_r(NULL, "\t", &this_col_ptr);
						}
					}

					if(xk3==xk2)continue;

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

						for(xk3=0; xk3<gap_merge_ptr; xk3++)
						{
							char numbbuf[12];
							BUFstrcat(out_chr_list, matched_chr, &tmp_chr_list);
							BUFstrcat(out_chr_list, ";", &tmp_chr_list);

							SUBreadSprintf(numbbuf,12,"%u;", input_start_stop_list[xk3 * 2]);
							BUFstrcat(out_start_list, numbbuf, &tmp_start_list);
							SUBreadSprintf(numbbuf,12,"%u;", input_start_stop_list[xk3 * 2 + 1] - 1);
							BUFstrcat(out_end_list, numbbuf, &tmp_end_list);
							SUBreadSprintf(numbbuf,12,"%c;", (matched_strand==1)?'-':( ( matched_strand==0 )? '+':'.'));
							BUFstrcat(out_strand_list, numbbuf, &tmp_strand_list);

						}
						for(xk3=0; xk3<merged_gaps; xk3++)
							gene_nonoverlap_len += output_start_stop_list[xk3 * 2 + 1] - output_start_stop_list[xk3 * 2];
				}
			}
		}
		#define _cut_tail(x) (x)[strlen(x)-1]=0

		_cut_tail(out_chr_list);
		_cut_tail(out_start_list);
		_cut_tail(out_end_list);
		_cut_tail(out_strand_list);

		char * QCcontent = "";
		char * QCtab = "";
		if(global_context -> do_detection_call){
			QCcontent = HashTableGet(global_context -> GCcontent_table, gene_symbol);
			QCtab = "\t";
			if(!QCcontent)QCcontent="nan";
		}

		int wlen = fprintf(fp_out, "%s\t%s%s%s\t%s\t%s\t%s\t%d", gene_symbol, QCcontent, QCtab, out_chr_list, out_start_list, out_end_list, out_strand_list, gene_nonoverlap_len);
		for(xk4 = 0; xk4<total_extra_cols; xk4++){
			merge_repeated_extra_columns(out_extra_columns[xk4]);
			fprintf(fp_out, "\t%s", out_extra_columns[xk4]+1);
		}

		for(i_files=0; i_files< column_names->numOfElements; i_files++)
		{
			read_count_type_t longlong_res = 0;
			double double_res = 0;
			int is_double_number = calc_float_fraction(gene_columns[i_files + column_names->numOfElements*xk1], &longlong_res, &double_res);
			if(is_double_number){
				fprintf(fp_out,"\t%.2f", double_res);
			}else{
				#ifdef __MINGW32__
				fprintf(fp_out,"\t%" PRIu64, (srInt_64)longlong_res);
				#else
				fprintf(fp_out,"\t%lld", (srInt_64)longlong_res);
				#endif
			}
		}
		fprintf(fp_out,"\n");
		if(wlen < 6)disk_is_full = 1;
	}

	for(xk1=0; xk1<total_extra_cols; xk1++) free(out_extra_columns[xk1]);
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
	free(gene_exons_extra_columns);
	free(gene_exons_start);
	free(gene_exons_end);
	free(gene_exons_strand);
	fclose(fp_out);

	if(disk_is_full){
		SUBREADprintf("ERROR: disk is full; the count file cannot be generated.\n");
		unlink(out_file);
	}
}

void fc_write_final_counts(fc_thread_global_context_t * global_context, const char * out_file, ArrayList * column_names, ArrayList * read_counters, int isCVersion)
{
	char fname[MAX_FILE_NAME_LENGTH];
	int i_files, xk1, disk_is_full = 0;

	SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH, "%s.summary", out_file);
	FILE * fp_out = f_subr_open(fname,"w");

	if(!fp_out){
		SUBREADprintf("Unable to create summary file '%s'\n", fname);
		return;
	}

	fprintf(fp_out,"Status");
	
	for(i_files=0; i_files<column_names->numOfElements; i_files++)
	{
		char * next_fn = ArrayListGet(column_names, i_files);
		fprintf(fp_out,"\t%s", global_context -> use_stdin_file?"STDIN":next_fn);
	}

	fprintf(fp_out,"\n");
	char * keys [] ={ "Assigned" ,  "Unassigned_Unmapped", "Unassigned_Read_Type", "Unassigned_Singleton", "Unassigned_MappingQuality", "Unassigned_Chimera", "Unassigned_FragmentLength", "Unassigned_Duplicate", "Unassigned_MultiMapping" , "Unassigned_Secondary",  (global_context->is_split_or_exonic_only == 2)?"Unassigned_Split":"Unassigned_NonSplit", "Unassigned_NoFeatures", "Unassigned_Overlapping_Length", "Unassigned_Ambiguity"};

	for(xk1=0; xk1<14; xk1++)
	{
		fprintf(fp_out,"%s", keys[xk1]);
		for(i_files = 0; i_files < column_names->numOfElements; i_files ++)
		{
			srInt_64 * array_0 = ArrayListGet(read_counters,i_files);
			srInt_64 * cntr = array_0 + xk1;
			#ifdef __MINGW32__
			fprintf(fp_out,"\t%" PRIu64, (srInt_64)*cntr);
			#else
			fprintf(fp_out,"\t%lld", (srInt_64)*cntr);
			#endif
		}
		int wlen = fprintf(fp_out,"\n");
		if(wlen < 1)disk_is_full = 1;
	}


	fclose(fp_out);

	if(disk_is_full){
		SUBREADprintf("ERROR: disk is full; the count file cannot be generated.\n");
		unlink(out_file);
	}

}
void fc_write_final_results(fc_thread_global_context_t * global_context, const char * out_file, int features, ArrayList* column_numbers, ArrayList * column_names,fc_feature_info_t * loaded_features, int header_out)
{
	/* save the results */
	FILE * fp_out;
	int i, i_files = 0, disk_is_full =0;
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



	char * next_fn;
	fprintf(fp_out,"Geneid\tChr\tStart\tEnd\tStrand\tLength");
	if(global_context -> reported_extra_columns)fprintf(fp_out,"\t%s", global_context -> reported_extra_columns);

	for(i_files = 0; i_files < column_names -> numOfElements; i_files++){
		next_fn = ArrayListGet(column_names, i_files);
		fprintf(fp_out,"\t%s", global_context -> use_stdin_file?"STDIN":next_fn);
	}
	fprintf(fp_out,"\n");
	for(i=0;i<features;i++)
	{
		fprintf(fp_out,"%s\t%s\t%u\t%u\t%c\t%d%s%s", global_context -> unistr_buffer_space + loaded_features[i].feature_name_pos,
 							   global_context -> unistr_buffer_space + loaded_features[i].feature_name_pos + loaded_features[i].chro_name_pos_delta,
						   	   loaded_features[i].start, loaded_features[i].end, loaded_features[i].is_negative_strand == 1?'-':(  loaded_features[i].is_negative_strand ==  0? '+':'.'),
							loaded_features[i].end-loaded_features[i].start+1, global_context -> reported_extra_columns ?"\t":"", global_context -> reported_extra_columns ?loaded_features[i].extra_columns:"");
		for(i_files=0; i_files < column_names -> numOfElements; i_files++)
		{
			srInt_64 * this_list = ArrayListGet(column_numbers, i_files);			
			int sorted_exon_no = loaded_features[i].sorted_order;
			srInt_64 count_frac_raw = this_list[sorted_exon_no], longlong_res = 0;

			double double_res = 0;
			int is_double_number = calc_float_fraction(count_frac_raw, &longlong_res, &double_res);
			if(is_double_number){
				fprintf(fp_out,"\t%.2f", double_res);
			}else{
				#ifdef __MINGW32__
				fprintf(fp_out,"\t%" PRId64, (srInt_64)longlong_res);
				#else
				fprintf(fp_out,"\t%lld", (srInt_64)longlong_res);
				#endif
			}
		}
		int wlen = fprintf(fp_out,"\n");
		if(wlen < 1)disk_is_full = 1;
	}

	fclose(fp_out);
	if(disk_is_full){
		SUBREADprintf("ERROR: disk is full; unable to write into the output file.\n");
		unlink(out_file);
	}
}

static struct option long_options[] =
{
	{"primary",no_argument, 0, 0},
	{"readShiftSize", required_argument, 0, 0},
	{"readShiftType", required_argument, 0, 0},
	{"readExtension5", required_argument, 0, 0},
	{"readExtension5", required_argument, 0, 0},
	{"readExtension3", required_argument, 0, 0},
	{"read2pos", required_argument, 0, 0},
	{"minOverlap", required_argument, 0, 0},
	{"fracOverlap", required_argument, 0, 0},
	{"nonOverlap", required_argument, 0, 0},
	{"nonOverlapFeature", required_argument, 0, 0},
	{"fracOverlapFeature", required_argument, 0, 0},
	{"splitOnly", no_argument, 0, 0},
	{"nonSplitOnly", no_argument, 0, 0},
	{"debugCommand", required_argument, 0, 0},
	{"ignoreDup", no_argument, 0, 0},
	{"donotsort", no_argument, 0, 0},
	{"restrictedlyNoOverlap", no_argument, 0, 0},
	{"fraction", no_argument, 0, 0},
	{"order", required_argument, 0, 'S'},
	{"genome", required_argument, 0, 'G'},
	{"maxMOp", required_argument, 0, 0},
	{"tmpDir", required_argument, 0, 0},
	{"extraAttributes", required_argument, 0, 0},
	{"largestOverlap", no_argument, 0,0},
	{"countReadPairs", no_argument, 0, 0},
	{"byReadGroup", no_argument, 0,0},
	{"verbose", no_argument, 0,0},
	{"detectionCall", no_argument, 0,0},
	{"Rpath", required_argument, 0, 0},
	{"scSampleSheet", required_argument, 0, 0},
	{"scInputMode", required_argument, 0, 0},
	{"scCellBarcodeFile", required_argument, 0, 0},
	{0, 0, 0, 0}
};

void print_usage()
{
	SUBREADprintf("\nVersion %s\n\n", SUBREAD_VERSION);

	SUBREADputs("Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... \n");
	SUBREADputs("## Mandatory arguments:");
	SUBREADputs("");
	SUBREADputs("  -a <string>         Name of an annotation file. GTF/GFF format by default. See");
	SUBREADputs("                      -F option for more format information. Inbuilt annotations");
	SUBREADputs("                      (SAF format) is available in 'annotation' directory of the");
	SUBREADputs("                      package. Gzipped file is also accepted.");
	SUBREADputs("");
	SUBREADputs("  -o <string>         Name of output file including read counts. A separate file");
	SUBREADputs("                      including summary statistics of counting results is also");
	SUBREADputs("                      included in the output ('<string>.summary'). Both files");
	SUBREADputs("                      are in tab delimited format."); 
	SUBREADputs("");
	SUBREADputs("  input_file1 [input_file2] ...   A list of SAM or BAM format files. They can be");
	SUBREADputs("                      either name or location sorted. If no files provided,");
	SUBREADputs("                      <stdin> input is expected. Location-sorted paired-end reads");
	SUBREADputs("                      are automatically sorted by read names.");
	SUBREADputs("");
	
	SUBREADputs("## Optional arguments:");
	SUBREADputs("# Annotation");
	SUBREADputs("");
	SUBREADputs("  -F <string>         Specify format of the provided annotation file. Acceptable");
	SUBREADputs("                      formats include 'GTF' (or compatible GFF format) and");
	SUBREADputs("                      'SAF'. 'GTF' by default.  For SAF format, please refer to");
	SUBREADputs("                      Users Guide.");
	SUBREADputs("");
	SUBREADputs("  -t <string>         Specify feature type(s) in a GTF annotation. If multiple");
	SUBREADputs("                      types are provided, they should be separated by ',' with");
	SUBREADputs("                      no space in between. 'exon' by default. Rows in the");
	SUBREADputs("                      annotation with a matched feature will be extracted and");
	SUBREADputs("                      used for read mapping. ");
	SUBREADputs("");
	SUBREADputs("  -g <string>         Specify attribute type in GTF annotation. 'gene_id' by ");
	SUBREADputs("                      default. Meta-features used for read counting will be ");
	SUBREADputs("                      extracted from annotation using the provided value.");
	SUBREADputs("");
	SUBREADputs("  --extraAttributes   Extract extra attribute types from the provided GTF");
	SUBREADputs("                      annotation and include them in the counting output. These");
	SUBREADputs("                      attribute types will not be used to group features. If");
	SUBREADputs("                      more than one attribute type is provided they should be");
	SUBREADputs("                      separated by comma.");
	SUBREADputs("");
	SUBREADputs("  -A <string>         Provide a chromosome name alias file to match chr names in");
	SUBREADputs("                      annotation with those in the reads. This should be a two-");
	SUBREADputs("                      column comma-delimited text file. Its first column should");
	SUBREADputs("                      include chr names in the annotation and its second column");
	SUBREADputs("                      should include chr names in the reads. Chr names are case");
	SUBREADputs("                      sensitive. No column header should be included in the");
	SUBREADputs("                      file.");
	SUBREADputs("");

	SUBREADputs("# Level of summarization");
	SUBREADputs("");
	SUBREADputs("  -f                  Perform read counting at feature level (eg. counting ");
	SUBREADputs("                      reads for exons rather than genes).");
	SUBREADputs("");

	SUBREADputs("# Overlap between reads and features");
	SUBREADputs("");
	SUBREADputs("  -O                  Assign reads to all their overlapping meta-features (or ");
	SUBREADputs("                      features if -f is specified).");
	SUBREADputs("");
	SUBREADputs("  --minOverlap <int>  Minimum number of overlapping bases in a read that is");
	SUBREADputs("                      required for read assignment. 1 by default. Number of");
	SUBREADputs("                      overlapping bases is counted from both reads if paired");
	SUBREADputs("                      end. If a negative value is provided, then a gap of up");
	SUBREADputs("                      to specified size will be allowed between read and the");
	SUBREADputs("                      feature that the read is assigned to.");
	SUBREADputs("");
	SUBREADputs("  --fracOverlap <float> Minimum fraction of overlapping bases in a read that is");
	SUBREADputs("                      required for read assignment. Value should be within range");
	SUBREADputs("                      [0,1]. 0 by default. Number of overlapping bases is");
	SUBREADputs("                      counted from both reads if paired end. Both this option");
	SUBREADputs("                      and '--minOverlap' option need to be satisfied for read");
	SUBREADputs("                      assignment.");
	SUBREADputs("");
	SUBREADputs("  --fracOverlapFeature <float> Minimum fraction of overlapping bases in a");
	SUBREADputs("                      feature that is required for read assignment. Value");
	SUBREADputs("                      should be within range [0,1]. 0 by default.");
	SUBREADputs("");
	SUBREADputs("  --largestOverlap    Assign reads to a meta-feature/feature that has the ");
	SUBREADputs("                      largest number of overlapping bases.");
	SUBREADputs("");
	SUBREADputs("  --nonOverlap <int>  Maximum number of non-overlapping bases in a read (or a");
	SUBREADputs("                      read pair) that is allowed when being assigned to a");
	SUBREADputs("                      feature. No limit is set by default.");
	SUBREADputs("");
	SUBREADputs("  --nonOverlapFeature <int> Maximum number of non-overlapping bases in a feature");
	SUBREADputs("                      that is allowed in read assignment. No limit is set by");
	SUBREADputs("                      default.");
	SUBREADputs("");
	SUBREADputs("  --readExtension5 <int> Reads are extended upstream by <int> bases from their");
	SUBREADputs("                      5' end.");
	SUBREADputs("");
	SUBREADputs("  --readExtension3 <int> Reads are extended upstream by <int> bases from their");
	SUBREADputs("                      3' end.");
	SUBREADputs("");
	SUBREADputs("  --read2pos <5:3>    Reduce reads to their 5' most base or 3' most base. Read");
	SUBREADputs("                      counting is then performed based on the single base the ");
	SUBREADputs("                      read is reduced to.");
	SUBREADputs("");

	SUBREADputs("# Multi-mapping reads");
	SUBREADputs("");
	SUBREADputs("  -M                  Multi-mapping reads will also be counted. For a multi-");
	SUBREADputs("                      mapping read, all its reported alignments will be ");
	SUBREADputs("                      counted. The 'NH' tag in BAM/SAM input is used to detect ");
	SUBREADputs("                      multi-mapping reads.");
	SUBREADputs("");
	SUBREADputs("# Fractional counting");
	SUBREADputs("");
	SUBREADputs("  --fraction          Assign fractional counts to features. This option must");
	SUBREADputs("                      be used together with '-M' or '-O' or both. When '-M' is");
	SUBREADputs("                      specified, each reported alignment from a multi-mapping");
	SUBREADputs("                      read (identified via 'NH' tag) will carry a fractional");
	SUBREADputs("                      count of 1/x, instead of 1 (one), where x is the total");
	SUBREADputs("                      number of alignments reported for the same read. When '-O'");
	SUBREADputs("                      is specified, each overlapping feature will receive a");
	SUBREADputs("                      fractional count of 1/y, where y is the total number of");
	SUBREADputs("                      features overlapping with the read. When both '-M' and");
	SUBREADputs("                      '-O' are specified, each alignment will carry a fractional");
	SUBREADputs("                      count of 1/(x*y).");
	SUBREADputs("");
	

	SUBREADputs("# Read filtering");
	SUBREADputs("");
	SUBREADputs("  -Q <int>            The minimum mapping quality score a read must satisfy in");
	SUBREADputs("                      order to be counted. For paired-end reads, at least one");
	SUBREADputs("                      end should satisfy this criteria. 0 by default.");
	SUBREADputs("");
	SUBREADputs("  --splitOnly         Count split alignments only (ie. alignments with CIGAR");
	SUBREADputs("                      string containing 'N'). An example of split alignments is");
	SUBREADputs("                      exon-spanning reads in RNA-seq data.");
	SUBREADputs("");
	SUBREADputs("  --nonSplitOnly      If specified, only non-split alignments (CIGAR strings do");
	SUBREADputs("                      not contain letter 'N') will be counted. All the other");
	SUBREADputs("                      alignments will be ignored.");
	SUBREADputs("");
	SUBREADputs("  --primary           Count primary alignments only. Primary alignments are ");
	SUBREADputs("                      identified using bit 0x100 in SAM/BAM FLAG field.");
	SUBREADputs("");
	SUBREADputs("  --ignoreDup         Ignore duplicate reads in read counting. Duplicate reads ");
	SUBREADputs("                      are identified using bit Ox400 in BAM/SAM FLAG field. The ");
	SUBREADputs("                      whole read pair is ignored if one of the reads is a ");
	SUBREADputs("                      duplicate read for paired end data.");
	SUBREADputs("");

	SUBREADputs("# Strandness");
	SUBREADputs("");
	SUBREADputs("  -s <int or string>  Perform strand-specific read counting. A single integer");
	SUBREADputs("                      value (applied to all input files) or a string of comma-");
	SUBREADputs("                      separated values (applied to each corresponding input");
	SUBREADputs("                      file) should be provided. Possible values include:");
	SUBREADputs("                      0 (unstranded), 1 (stranded) and 2 (reversely stranded).");
	SUBREADputs("                      Default value is 0 (ie. unstranded read counting carried");
	SUBREADputs("                      out for all input files).");
	SUBREADputs("");

	SUBREADputs("# Exon-exon junctions");
	SUBREADputs("");
	SUBREADputs("  -J                  Count number of reads supporting each exon-exon junction.");
	SUBREADputs("                      Junctions were identified from all the exon-spanning reads");
	SUBREADputs("                      in the input (containing 'N' in CIGAR string). Counting");
	SUBREADputs("                      results are saved to a file named '<output_file>.jcounts'");
	SUBREADputs("");
	SUBREADputs("  -G <string>         Provide the name of a FASTA-format file that contains the");
	SUBREADputs("                      reference sequences used in read mapping that produced the");
	SUBREADputs("                      provided SAM/BAM files. This optional argument can be used");
	SUBREADputs("                      with '-J' option to improve read counting for junctions.");
	SUBREADputs("");

	SUBREADputs("# Parameters specific to paired end reads");
	SUBREADputs("");
	SUBREADputs("  -p                  If specified, libraries are assumed to contain paired-end");
	SUBREADputs("                      reads. For any library that contains paired-end reads, the");
	SUBREADputs("                      'countReadPairs' parameter controls if read pairs or reads");
	SUBREADputs("                      should be counted.");
	SUBREADputs("");
	SUBREADputs("  --countReadPairs    If specified, fragments (or templates) will be counted");
	SUBREADputs("                      instead of reads. This option is only applicable for");
	SUBREADputs("                      paired-end reads. For single-end data, it is ignored.");
	SUBREADputs("");
	SUBREADputs("  -B                  Only count read pairs that have both ends aligned.");
	SUBREADputs("");
	SUBREADputs("  -P                  Check validity of paired-end distance when counting read ");
	SUBREADputs("                      pairs. Use -d and -D to set thresholds.");
	SUBREADputs("");
	SUBREADputs("  -d <int>            Minimum fragment/template length, 50 by default.");
	SUBREADputs("");
	SUBREADputs("  -D <int>            Maximum fragment/template length, 600 by default.");
	SUBREADputs("");
	SUBREADputs("  -C                  Do not count read pairs that have their two ends mapping ");
	SUBREADputs("                      to different chromosomes or mapping to same chromosome ");
	SUBREADputs("                      but on different strands.");
	SUBREADputs("");
	SUBREADputs("  --donotsort         Do not sort reads in BAM/SAM input. Note that reads from ");
	SUBREADputs("                      the same pair are required to be located next to each ");
	SUBREADputs("                      other in the input.");
	SUBREADputs("");

	SUBREADputs("# Number of CPU threads");
	SUBREADputs("");
	SUBREADputs("  -T <int>            Number of the threads. 1 by default.");
	SUBREADputs("");

	SUBREADputs("# Read groups");
	SUBREADputs("");
	SUBREADputs("  --byReadGroup       Assign reads by read group. \"RG\" tag is required to be");
	SUBREADputs("                      present in the input BAM/SAM files.");
	SUBREADputs("                      ");
	SUBREADputs("");

	SUBREADputs("# Long reads");
	SUBREADputs("");
	SUBREADputs("  -L                  Count long reads such as Nanopore and PacBio reads. Long");
	SUBREADputs("                      read counting can only run in one thread and only reads");
	SUBREADputs("                      (not read-pairs) can be counted. There is no limitation on");
	SUBREADputs("                      the number of 'M' operations allowed in a CIGAR string in");
	SUBREADputs("                      long read counting.");
	SUBREADputs("");

	SUBREADputs("# Assignment results for each read");
	SUBREADputs("");
	SUBREADputs("  -R <format>         Output detailed assignment results for each read or read-");
	SUBREADputs("                      pair. Results are saved to a file that is in one of the");
	SUBREADputs("                      following formats: CORE, SAM and BAM. See Users Guide for");
	SUBREADputs("                      more info about these formats.");
	SUBREADputs("");
	SUBREADputs("  --Rpath <string>    Specify a directory to save the detailed assignment");
	SUBREADputs("                      results. If unspecified, the directory where counting");
	SUBREADputs("                      results are saved is used.");
	SUBREADputs("");

	SUBREADputs("# Miscellaneous");
	SUBREADputs("");
	SUBREADputs("  --tmpDir <string>   Directory under which intermediate files are saved (later");
	SUBREADputs("                      removed). By default, intermediate files will be saved to");
	SUBREADputs("                      the directory specified in '-o' argument.");
	SUBREADputs("");
	SUBREADputs("  --verbose           Output verbose information for debugging, such as un-");
	SUBREADputs("                      matched chromosome/contig names.");
	SUBREADputs("");
	SUBREADputs("  -v                  Output version of the program.");
	SUBREADputs("");

}

int junckey_sort_compare(void * inptr, int i, int j){
	char ** inp = (char **) inptr;
	int x1;

	int chrI=-1, chrJ=-1;

	if(atoi(inp[i])>0) chrI = atoi(inp[i]);
	if(atoi(inp[j])>0) chrJ = atoi(inp[j]);

	if(inp[i][0]=='X' && !isdigit(inp[i][1])&& !isalpha(inp[i][1])) chrI = 90;
	if(inp[i][0]=='Y' && !isdigit(inp[i][1])&& !isalpha(inp[i][1])) chrI = 91;
	if(inp[i][0]=='M' && !isdigit(inp[i][1])&& !isalpha(inp[i][1])) chrI = 99;
	if(inp[j][0]=='X' && !isdigit(inp[j][1])&& !isalpha(inp[j][1])) chrJ = 90;
	if(inp[j][0]=='Y' && !isdigit(inp[j][1])&& !isalpha(inp[j][1])) chrJ = 91;
	if(inp[j][0]=='M' && !isdigit(inp[j][1])&& !isalpha(inp[j][1])) chrJ = 99;



	if(memcmp(inp[i], "chr", 3)==0){
		chrI=atoi(inp[i]+3);
		if(0 == chrI && inp[i][3] == 'X') chrI = 90;
		if(0 == chrI && inp[i][3] == 'Y') chrI = 91;
		if(0 == chrI && inp[i][3] == 'M') chrI = 99;
	}
	if(memcmp(inp[j], "chr", 3)==0){
		chrJ=atoi(inp[j]+3);
		if(0 == chrJ && inp[j][3] == 'X') chrJ = 90;
		if(0 == chrJ && inp[j][3] == 'Y') chrJ = 91;
		if(0 == chrJ && inp[j][3] == 'M') chrJ = 99;
	}

	int len_I_long = 9;
	for(x1 = 0 ; x1 < FEATURE_NAME_LENGTH + 15 ; x1++){
		int c1 = inp[i][x1];
		int c2 = inp[j][x1];
		if(c1 == '\t' && c2 != '\t')
			len_I_long = -1;
		else if(c1 != '\t' && c2 == '\t')
			len_I_long = 1;
		else if(c1 == '\t' && c2 == '\t')
			len_I_long = 0;

		if(len_I_long != 9) break;
	}

	if(chrI != chrJ || len_I_long != 0){
		return (chrI * 100 + len_I_long) - (chrJ * 100);
	}

	for(x1 = 0 ; x1 < FEATURE_NAME_LENGTH + 15 ; x1++){
		int c1 = inp[i][x1];
		int c2 = inp[j][x1];
		if(c1 != c2){
			return c1 - c2;
		}else if(c1 == '\t' && c1 == c2){
			int pos1 = atoi(inp[i]+x1+1);
			int pos2 = atoi(inp[j]+x1+1);
			if( pos1 == pos2)
				return strcmp(inp[i], inp[j]);
			else
				return pos1 - pos2;
		}

		if(c1 == 0 || c2 == 0)return c1 - c2;
	}
	return 0;
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


#define JC_STATUS_KNOWN 1
#define JC_STATUS_NOVEL 2
#define JC_STATUS_NOVEL_FUSION 3
#define MAX_OVERLAP_EDGE_NUMBER 1000
#define JC_OUT_GENE_COLUMNS_LENGTH (MAX_OVERLAP_EDGE_NUMBER * (4+CHROMOSOME_NAME_LENGTH)) 
int determine_jcount_gene_transcript_report(fc_thread_global_context_t * global_context, int side_small, int side_large, int junc_olay_genebody_left_result_no,int junc_olay_genebody_right_result_no, IVT_Interval ** junc_genebody_olayleft, IVT_Interval ** junc_genebody_olayright, char * gene_ids_str, char * transcript_ids_str, int * num_exons, char * dist_to_nearest_splice_side_str){
	// test "exact hit edges" situation
	int xk1,xk2;
	HashTable * match1_txn_table = StringTableCreate(100);
	HashTable * edge1P1_table = StringTableCreate(100);
	HashTable * edge2P1_table = StringTableCreate(100);
	for(xk1=0; xk1<junc_olay_genebody_left_result_no; xk1++){
		srInt_64 edge1_dist = 0xffffffffu;
		fc_junction_genebody_t * jg_ptr = junc_genebody_olayleft[xk1]->attr;
		ArrayList * exon_list = jg_ptr -> exon_list;
		for(xk2=0; xk2< exon_list->numOfElements; xk2++){
			fc_junction_exon_in_transcript_t * jte_ptr = ArrayListGet(exon_list,xk2);
			int exon_end_known = jte_ptr -> chro_stop;

			// NULL+1 : txn has an exon that matched one edge.
			// NULL+2 : txn has no exon that can match any edge.
			if(side_small == exon_end_known) HashTablePut(match1_txn_table, jte_ptr -> transcript_id, NULL+1);
			else if( HashTableGet(match1_txn_table, jte_ptr -> transcript_id)!= NULL+1 ) HashTablePut(match1_txn_table, jte_ptr -> transcript_id, NULL+2);
			if(side_small - exon_end_known < edge1_dist && side_small >= exon_end_known){
				edge1_dist = side_small - exon_end_known;
			}
		}
		HashTablePut(edge1P1_table, jg_ptr -> gene_name, NULL+1+edge1_dist);
	}
	
	HashTable * two_matched_txn_ids_tab = StringTableCreate(10);
	for(xk1=0; xk1<junc_olay_genebody_right_result_no; xk1++){
		srInt_64 edge2_dist = 0xffffffffu;
		fc_junction_genebody_t * jg_ptr = junc_genebody_olayright[xk1]->attr;
		ArrayList * exon_list = jg_ptr -> exon_list;
		for(xk2=0; xk2< exon_list->numOfElements; xk2++){
			fc_junction_exon_in_transcript_t * jte_ptr = ArrayListGet(exon_list,xk2);
			char * txn_id = jte_ptr -> transcript_id;
			int exon_start_known = jte_ptr -> chro_start; 
			if(exon_start_known == side_large && NULL+1==HashTableGet(match1_txn_table, txn_id)) HashTablePut(two_matched_txn_ids_tab, txn_id, NULL+1);
			if(exon_start_known - side_large < edge2_dist && side_large <= exon_start_known){
				edge2_dist = exon_start_known - side_large;
			}
		}
		HashTablePut(edge2P1_table, jg_ptr -> gene_name, NULL+1+edge2_dist);
	}
//	fprintf(stderr,"COMM1  %lld  HASHS %lld  HITS %d %d\n", two_matched_txn_ids_tab -> numOfElements , match1_txn_table -> numOfElements, junc_olay_genebody_left_result_no, junc_olay_genebody_right_result_no);
	ArrayList * two_matched_txn_ids_list = HashTableKeys(two_matched_txn_ids_tab); 
	HashTableDestroy(two_matched_txn_ids_tab);

	int retv = 0;

	// mode 1 : the two sides of the junction matches exactly two exons in the same transcript.
	if(two_matched_txn_ids_list -> numOfElements >0){
		HashTable * gene_tab = StringTableCreate(100);
		ArrayList * txn_id_list = ArrayListCreate(10);
		for(xk1=0; xk1< two_matched_txn_ids_list -> numOfElements; xk1++){
			char * txn_id = ArrayListGet(two_matched_txn_ids_list, xk1);
			fc_junction_transcript_t * txn_1 = HashTableGet(global_context -> junction_transcript_table, txn_id); 
			HashTablePut(gene_tab, txn_1 -> gene_name, NULL+1);
			ArrayListPush(txn_id_list, txn_id);
		}
		ArrayListSort(txn_id_list, ArrayListStringComparison);
		ArrayList * gene_name_list = HashTableKeys(gene_tab);
		ArrayListSort(gene_name_list, ArrayListStringComparison);

		ArrayListStringJoin(gene_name_list, gene_ids_str, JC_OUT_GENE_COLUMNS_LENGTH);
		ArrayListStringJoin(txn_id_list,transcript_ids_str, JC_OUT_GENE_COLUMNS_LENGTH);
		
		strcpy(dist_to_nearest_splice_side_str,"0/0");
		ArrayListDestroy(gene_name_list);
		ArrayListDestroy(txn_id_list);
		*num_exons = 2;
		retv = JC_STATUS_KNOWN;
	}
	ArrayListDestroy(two_matched_txn_ids_list);

	if(!retv){
		// mode 2: two sides matched the same gene body.
		HashTable * two_mismatched_gene_tab = StringTableCreate(10);
		HashTable * two_fusion_gene_tab = StringTableCreate(10);
		ArrayList * mis2_genes = HashTableKeys(edge2P1_table);
		ArrayList * mis1_genes = HashTableKeys(edge1P1_table);
		for(xk1 = 0; xk1< mis1_genes -> numOfElements; xk1++){
			char * mis1_gene_name = ArrayListGet(mis1_genes,xk1); 
			if(HashTableGet(edge2P1_table, mis1_gene_name))
				HashTablePut(two_mismatched_gene_tab, mis1_gene_name, NULL+1);
		}

		if(two_mismatched_gene_tab -> numOfElements>0){
			ArrayList * gene_name_list = HashTableKeys(two_mismatched_gene_tab),
				* gene_name_list_best = ArrayListCreate(10),
				* gene_dist_best = ArrayListCreate(10);
			ArrayListSort(gene_name_list, ArrayListStringComparison);
			ArrayListSetDeallocationFunction(gene_dist_best,free);
			transcript_ids_str[0]=0;
			int best_nexons = 0;
			for(xk1 = 0; xk1 < gene_name_list-> numOfElements; xk1++){
				char * gene_name = ArrayListGet(gene_name_list,xk1);
				srInt_64 dist1 = HashTableGet(edge1P1_table, gene_name) - NULL - 1 ;
				srInt_64 dist2 = HashTableGet(edge2P1_table, gene_name) - NULL - 1 ;
				int this_exons=0;
				if(dist1 == 0 && dist2 == 0) this_exons=2;
				else if(dist1 == 0 || dist2 == 0) this_exons=1;
				best_nexons = max(best_nexons, this_exons);
			}

			for(xk1 = 0; xk1 < gene_name_list-> numOfElements; xk1++){
				char * gene_name = ArrayListGet(gene_name_list,xk1);
				int dist1 = HashTableGet(edge1P1_table, gene_name) - NULL - 1 ;
				int dist2 = HashTableGet(edge2P1_table, gene_name) - NULL - 1 ;
//if(side_small==539704 && side_large==546211) fprintf(stderr,"TEST_ONEDEGE %d %d\n", dist1, dist2);
				int this_exons=0;
				if(dist1 == 0 && dist2 == 0) this_exons=2;
				else if(dist1 == 0 || dist2 == 0) this_exons=1;
				if(this_exons==best_nexons){
					char* gene_dist_str = malloc(30);
					if(dist1<0 && dist2>=0) sprintf(gene_dist_str,"NA/%d", dist2);
					else if(dist1>=0 && dist2<0) sprintf(gene_dist_str,"%d/NA", dist1);
					else if(dist1<0 && dist2<0) strcpy(gene_dist_str,"NA/NA");
					else sprintf(gene_dist_str,"%d/%d", dist1,dist2);
					ArrayListPush(gene_name_list_best, gene_name);
					ArrayListPush(gene_dist_best, gene_dist_str);
				}
			}
	
			*num_exons = best_nexons;
			ArrayListStringJoin(gene_name_list_best, gene_ids_str, JC_OUT_GENE_COLUMNS_LENGTH);
			ArrayListStringJoin(gene_dist_best, dist_to_nearest_splice_side_str, JC_OUT_GENE_COLUMNS_LENGTH);
			ArrayListDestroy(gene_name_list_best);
			ArrayListDestroy(gene_dist_best);
			ArrayListDestroy(gene_name_list);
			retv = JC_STATUS_NOVEL;
		}
		HashTableDestroy(two_mismatched_gene_tab);

		ArrayListDestroy(mis1_genes);
		ArrayListDestroy(mis2_genes);
//
		HashTableDestroy(two_fusion_gene_tab);
	}
	if( edge1P1_table->numOfElements && edge2P1_table -> numOfElements && !retv){
		// mode 3: two sides matched different gene bodies.
		ArrayList * gene_name_list = HashTableKeys(edge1P1_table),
				* edge2_gene_list = HashTableKeys(edge2P1_table),
				* gene_name_list_best = ArrayListCreate(10),
				* gene_dist_best = ArrayListCreate(10);
		ArrayListExtend(gene_name_list,edge2_gene_list);

		int best_nexons1 = 0, best_nexons2 = 0;
		for(xk1 = 0; xk1 < gene_name_list-> numOfElements; xk1++){
			char * gene_name = ArrayListGet(gene_name_list,xk1);
			srInt_64 dist1 = HashTableGet(edge1P1_table, gene_name) - NULL - 1 ;
			srInt_64 dist2 = HashTableGet(edge2P1_table, gene_name) - NULL - 1 ;
			if(dist1 == 0 && best_nexons1 ==0) best_nexons1=1;
			if(dist2 == 0 && best_nexons2 ==0) best_nexons2=1;
		}

		ArrayListSort(gene_name_list, ArrayListStringComparison);
		for(xk1 = 0; xk1 < gene_name_list-> numOfElements; xk1++){
			char * gene_name = ArrayListGet(gene_name_list,xk1);
			void * dist1ptr = HashTableGet(edge1P1_table, gene_name) ;
			void * dist2ptr = HashTableGet(edge2P1_table, gene_name) ;
			int dist1 = dist1ptr - NULL - 1 ;
			int dist2 = dist2ptr - NULL - 1 ;
			if(((best_nexons1 == 0 || dist1 == 0) && dist1ptr) || ((best_nexons2 == 0 || dist2 == 0)&& dist2ptr)){
				char * gene_dist_str = malloc(30);
				if(dist1<0 && dist2>=0) sprintf(gene_dist_str,"NA/%d", dist2);
				else if(dist1>=0 && dist2<0) sprintf(gene_dist_str,"%d/NA", dist1);
				else if(dist1<0 && dist2<0) strcpy(gene_dist_str,"NA/NA");
				else sprintf(gene_dist_str,"%d/%d", dist1,dist2);
				ArrayListPush(gene_name_list_best, gene_name);
				ArrayListPush(gene_dist_best, gene_dist_str);
			}
		}
		ArrayListStringJoin(gene_name_list_best, gene_ids_str, JC_OUT_GENE_COLUMNS_LENGTH);
		ArrayListStringJoin(gene_dist_best, dist_to_nearest_splice_side_str, JC_OUT_GENE_COLUMNS_LENGTH);
	
//fprintf(stderr,"G_FUSE_LISTS\t%d  %d\t%s\t%s\t%s\n", side_small,side_large , gene_ids_str, dist_to_nearest_splice_side_str, transcript_ids_str);
		ArrayListDestroy(gene_dist_best);
		ArrayListDestroy(gene_name_list_best);
		ArrayListDestroy(gene_name_list);
		ArrayListDestroy(edge2_gene_list);

		*num_exons = best_nexons1 + best_nexons2;
		retv = JC_STATUS_NOVEL_FUSION;
	}

	HashTableDestroy(match1_txn_table);
	HashTableDestroy(edge1P1_table);
	HashTableDestroy(edge2P1_table);

	return retv;
}

void fc_write_final_junctions(fc_thread_global_context_t * global_context,  char * output_file_name, ArrayList * column_names, ArrayList * junction_global_table_list, ArrayList * splicing_global_table_list){
	int infile_i, disk_is_full = 0;

	HashTable * merged_junction_table = HashTableCreate(156679);

	HashTableSetHashFunction(merged_junction_table,HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(merged_junction_table, NULL, NULL);
	HashTableSetKeyComparisonFunction(merged_junction_table, fc_strcmp_chro);

	HashTable * merged_splicing_table = HashTableCreate(156679);

	HashTableSetHashFunction(merged_splicing_table,HashTableStringHashFunction);
	HashTableSetDeallocationFunctions(merged_splicing_table, NULL, NULL);
	HashTableSetKeyComparisonFunction(merged_splicing_table, fc_strcmp_chro);


	for(infile_i = 0 ; infile_i < column_names -> numOfElements ; infile_i ++){
		KeyValuePair * cursor;
		int bucket;
		HashTable * spl_table = ArrayListGet(splicing_global_table_list, infile_i);
		for(bucket=0; bucket < spl_table -> numOfBuckets; bucket++)
		{
			cursor = spl_table -> bucketArray[bucket];
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

	for(infile_i = 0 ; infile_i < column_names -> numOfElements ; infile_i ++){
		KeyValuePair * cursor;
		int bucket;
		HashTable *  junc_table = ArrayListGet(junction_global_table_list, infile_i);
		for(bucket=0; bucket < junc_table -> numOfBuckets; bucket++)
		{
			cursor = junc_table -> bucketArray[bucket];
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

	char outfname[MAX_FILE_NAME_LENGTH];
	SUBreadSprintf(outfname, MAX_FILE_NAME_LENGTH, "%s.jcounts", output_file_name);

	int max_junction_genes = 3000;
	char * gene_names = malloc(max_junction_genes * FEATURE_NAME_LENGTH), * gene_name_tail;

	int ky_i1, ky_i2;
	FILE * ofp = fopen(outfname, "w");
	char * tmpp = NULL;

	fprintf(ofp, "GeneName\tTranscriptID\tStatus\tnExons\tDist.to.nearest.splicing.site\tSite1_chr\tSite1_location\tSite1_strand\tSite2_chr\tSite2_location\tSite2_strand");

	for(infile_i=0; infile_i < column_names -> numOfElements; infile_i++)
	{
		char * next_fn = ArrayListGet(column_names, infile_i);		
		fprintf(ofp,"\t%s", global_context -> use_stdin_file?"STDIN":next_fn);
	}
	fprintf(ofp, "\n");

	IVT_Interval ** junc_genebody_olayleft = malloc(sizeof(void*) * MAX_OVERLAP_EDGE_NUMBER);
	IVT_Interval ** junc_genebody_olayright = malloc(sizeof(void*) * MAX_OVERLAP_EDGE_NUMBER);
	for(ky_i = 0; ky_i < merged_junction_table -> numOfElements ; ky_i ++){
		int unique_junctions = 0;
		char * chro_small = strtok_r( key_list[ky_i] , "\t", &tmpp);
		char * pos_small_str = strtok_r( NULL, "\t", &tmpp);
		char * chro_large = strtok_r( NULL, "\t", &tmpp);
		char * pos_large_str = strtok_r( NULL, "\t", &tmpp);

		unsigned int pos_small = atoi(pos_small_str);
		unsigned int pos_large = atoi(pos_large_str);

		char * strand = "NA";
		if(global_context -> fasta_contigs){
			char donor[3], receptor[3];
			donor[2]=receptor[2]=0;
			int has = !get_contig_fasta(global_context -> fasta_contigs, chro_small, pos_small, 2, donor);
			has = has && !get_contig_fasta(global_context -> fasta_contigs, chro_large, pos_large-3, 2, receptor);
			if(has){
				if(donor[0]=='G' && donor[1]=='T' && receptor[0]=='A' && receptor[1]=='G') strand = "+";
				else if(donor[0]=='C' && donor[1]=='T' && receptor[0]=='A' && receptor[1]=='C') strand = "-";
			}else if(!global_context ->is_junction_no_chro_shown){
				global_context ->is_junction_no_chro_shown = 1;
				print_in_box(80,0,0, "   WARNING contig '%s' is not found in the", chro_small);
				print_in_box(80,0,0, "   provided genome file or its length is too short in it.");
				print_in_box(80,0,0,"");

			}
		}
		assert(0==strcmp(chro_small, chro_large));
		IVT_IntervalTreeNode * IVT_gbody_root = HashTableGet(global_context -> junction_GenebodyTree_table, chro_small);
		char gene_ids_str[JC_OUT_GENE_COLUMNS_LENGTH], transcript_ids_str[JC_OUT_GENE_COLUMNS_LENGTH], dist_to_nearest_splice_side_str[JC_OUT_GENE_COLUMNS_LENGTH];
		strcpy(gene_ids_str,"NA");
		strcpy(transcript_ids_str,"NA");
		strcpy(dist_to_nearest_splice_side_str,"NA");
		int num_exons=0;

//fprintf(stderr,"IVT_PTR %p  %s\n", IVT_gbody_root, chro_small);
		char * jc_retv_str = "NA";
		int jc_gene_status = -1;
		if(IVT_gbody_root){
			int junc_olay_genebody_left_result_no, junc_olay_genebody_right_result_no;
			junc_olay_genebody_left_result_no = IVT_query(IVT_gbody_root, pos_small, junc_genebody_olayleft, MAX_OVERLAP_EDGE_NUMBER);
			junc_olay_genebody_right_result_no = IVT_query(IVT_gbody_root, pos_large, junc_genebody_olayright, MAX_OVERLAP_EDGE_NUMBER);

			if(junc_olay_genebody_left_result_no == MAX_OVERLAP_EDGE_NUMBER|| junc_olay_genebody_right_result_no==MAX_OVERLAP_EDGE_NUMBER){
				SUBREADprintf("WARNING: Your annotation file contains very many exons that start/end at the same location. Consider to increase MAX_OVERLAP_EDGE_NUMBER in readSummary.c to accomodate these exons for junction counting.\n");
			}

			jc_gene_status = determine_jcount_gene_transcript_report(global_context , pos_small, pos_large, junc_olay_genebody_left_result_no, junc_olay_genebody_right_result_no, junc_genebody_olayleft, junc_genebody_olayright, gene_ids_str, transcript_ids_str, &num_exons, dist_to_nearest_splice_side_str);
		}
		if(jc_gene_status == JC_STATUS_KNOWN) jc_retv_str = "KNOWN";
		if(jc_gene_status == JC_STATUS_NOVEL) jc_retv_str = "NOVEL";
		if(jc_gene_status == JC_STATUS_NOVEL_FUSION) jc_retv_str = "NOVEL_FUSION";

		fprintf(ofp, "%s\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%s", gene_ids_str, transcript_ids_str, jc_retv_str, num_exons, dist_to_nearest_splice_side_str, chro_small, pos_small, strand, chro_large,pos_large, strand);

		*(pos_small_str-1)='\t';
		*(pos_large_str-1)='\t';
		chro_large[-1]='\t';

		for(infile_i = 0 ; infile_i < column_names -> numOfElements ; infile_i ++){
			HashTable * junc_table = ArrayListGet(junction_global_table_list, infile_i);
			srInt_64 count = HashTableGet(junc_table, key_list[ky_i]) - NULL;
			#ifdef __MINGW32__
			fprintf(ofp,"\t%" PRId64, count);
			#else
			fprintf(ofp,"\t%lld", count);
			#endif
		}
		int wlen = fprintf(ofp, "\n");
		if(wlen < 1) disk_is_full = 1;
	}
	fclose(ofp);
	free(gene_names);
	free(key_list);
	free(junc_genebody_olayleft);
	free(junc_genebody_olayright);

	//print_in_box(80,0,PRINT_BOX_CENTER,"Found %llu junctions in all the input files.", merged_junction_table -> numOfElements);
	//print_in_box(80,0,0,"");

	HashTableDestroy(merged_junction_table);
	HashTableDestroy(merged_splicing_table);
	if(disk_is_full){
		unlink(outfname);
		SUBREADprintf("ERROR: disk is full; no junction counting table is generated.\n");
	}
}

int readSummary_single_file(fc_thread_global_context_t * global_context, read_count_type_t * column_numbers, srInt_64 nexons,  int * geneid, char ** chr, srInt_64 * start, srInt_64 * stop, unsigned char * sorted_strand, char * anno_chr_2ch, char ** anno_chrs, srInt_64 * anno_chr_head, srInt_64 * block_end_index, srInt_64 * block_min_start , srInt_64 * block_max_end, fc_read_counters * my_read_counter, HashTable * junc_glob_tab, HashTable * splicing_glob_tab, HashTable * merged_RG_table, fc_feature_info_t * loaded_features);

int Input_Files_And_Strand_Mode_Pair(char * fnames, char * smodes){
	int ret = 0, ch, bad_fmt = 0, numbs = 0;
	//SUBREADputs(fnames);
	//SUBREADputs(smodes);
	if(strstr(smodes, ".")==NULL){
		bad_fmt = smodes[0]<'0' || smodes[0]>'2';
	}else{
		while('\0'!=(ch=*(fnames++)))if(ch == FC_FLIST_SPLITOR[0])ret++;
		while('\0'!=(ch=*(smodes++))){
			if(ch == '.'){
				if(numbs != 1) bad_fmt = 1;
				numbs = 0;
				ret--;
			}else if(ch >= '0' && ch <= '2') numbs++;
		}
		if(numbs != 1) bad_fmt = 1;
	}
	if(bad_fmt) SUBREADputs("Error: The strand mode list has a wrong format.");
	if(ret) SUBREADputs("Error: The length of strand mode list differs from the length of input file list");
	ret |= bad_fmt;
	return ret;
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
	12: isStrandChecked
	13: as.numeric(isReadSummaryReported)
	14: as.numeric(isBothEndMapped)
	15: as.numeric(isChimericDisallowed)
	16: as.numeric(isPEDistChecked)
	17: nameFeatureTypeColumn 
	18: nameGeneIDColumn
	19: min.MappingQualityScore
	20: as.numeric(isMultiMappingAllowed) # "1" : NH tag > 1 is allowed  ; "0" : not allowd (by default)
	21: Annotation Chromosome Alias Name File. If the file is not specified, set this value to NULL or a zero-length string.
	22: Command line for CfeatureCounts header output; RfeatureCounts should set this value to NULL or a zero-length string or a space (' ').
	23: as.numeric(isInputFileResortNeeded)
	24: NOT IN USE: as.numeric(feature_block_size) # This parameter is no longer used. Give "14" for safe. 
	25: as.numeric(Five_End_Extension_Length)  # 5' end extension
	26: as.numeric(Three_End_Extension_Length)  # 3' end extension
	27: as.numeric(Minimum_Overlap_Between_Read_And_Feature) # 1 by default
	28: as.numeric(is_Split_or_Exonic_Only) # 0 by default; 0: all reads are counted ; 1: only split (Cigar has "N") reads are counted ; 2: only exonic (no "N" in Cigar) are counted.
	29: as.numeric(reduce_5_3_ends_to_one) # 0= no reduction; 1= reduce to 5' end; 2= reduce to 3' end
	30: debug_command # This is for debug only; RfeatureCounts should pass a space (" ") to this parameter, disabling the debug command.
	31: as.numeric(is_duplicate_ignored) # 0 = INCLUDE DUPLICATE READS; 1 = IGNORE DUPLICATE READS (0x400 FLAG IS SET) ; "0" by default.
	32: as.numeric(do_not_sort)   # 1 = NEVER SORT THE PE BAM/SAM FILES; 0 = SORT THE BAM/SAM FILE IF IT IS FOUND NOT SORTED.
	33: as.numeric(fractionMultiMapping) # 1 = calculate fraction numbers if a read overlaps with multiple features or meta-features. "-M" must be specified when fractions are caculated.
	34: as.numeric(useOverlappingBreakTie) # 1 = Select features or meta-features with a longer overlapping length; 0 = just use read-voting strategy: one overlapping read = 1 vote
	35: Pair_Orientations # FF, FR, RF or RR. This parameter matters only if "-s" option is 1 or 2.
	36: as.numeric(doJunctionCounting)  # 1 = count the number of junction reads spaining each exon-exon pairs;  0 = do not.
	37: file name of genome fasta (for determine the strandness of junctions by looking for GT/AG or CT/AC).
	38: as.numeric(max_M_Ops) # maximum "M" sections allowed in the CIGAR string. This parameter is needed in parse_BIN()
	39: as.numeric(is_Restrictly_No_Overlapping) # when "1", disable the voting-based tie breaking (e.g., when the reads are paired-end and one gene receives two votes but the other gene only has one.). "0" by default.
	40: as.numeric(min_Fractional_Overlap) # A fractioal number.  0.00 : at least 1 bp overlapping
	41: temp_directory # the directory to put temp files. "<use output directory>" by default, namely find it from the output file dir.
	42: as.numeric(use_stdin_stdout) # only for CfeatureCounts. When use_stdin_stdout & 0x01 > 0, the input file is from stdin (stored in a temporary file); when use_stdin_stdout & 0x02 > 0, the output should be written to STDOUT instead of a file.
	43: as.numeric(assign_reads_to_RG) # 1: reads with "RG" tags will be assigned to read groups' 0: default setting
	44: as.numeric(long_read_minimum_length) # "1": treat the input BAM or SAM files as containing long reads. No multi-threading. "0": classic behaviour.
	45: as.numeric(is_verbose) # 1: show the mismatched chromosome names on screet; 0: don't do so
	46: as.numeric(frac_feature_overlap) # fraction of the feature to be overlapped with a read
	47: as.numeric(do_detection_call) # do detectionCalls : put the GC fraction into the 2nd column.
	48: as.numeric(max_missing_bases_in_read) # maximum # of bases in a read or fragment not overlapping with an exon ; efault value: "-1" means no limit
	49: as.numeric(max_missing_bases_in_feature) # maximum # of bases in an exon not overlapping with a read or fragment ; default value: "-1" means no limit
	50: as.numeric(is_Primary_Alignment_only) # "1" : only count the primary alignment (FLAG doesn't have 0x100 bit); "0" : count alignments no metter the 0x100 bit (by default)
	51: Rpath : the path where the assignment details per read are stored.
	52: AdditionalColumnList: the names of additional column names written after "Length". Comma deliminated.
	53: annotation_file_screen_output : just for displaying the annotation file name or inbuilt (mm10/hg39/...) or R data.frame.
	54: read_shift_type : how to shift reads? "upstream" : to the 5' end; "downstream" : to the 3' end; "left" : to the smaller coordinates in chromosome ; "right" : to the larger coordinates in chromosome.
	55: as.numeric(read_shift_size) : how many bases to shift. Mush be a positive number or zero.
	 */

	int isCVersion, isChimericDisallowed, isPEDistChecked, minMappingQualityScore=0, isInputFileResortNeeded, feature_block_size = 20, reduce_5_3_ends_to_one, useStdinFile, assignReadsToRG, long_read_minimum_length, is_verbose, do_detectionCall, max_missing_bases_in_feature, max_missing_bases_in_read, is_Primary_Alignment_only, read_shift_size, read_shift_type;
	float fracOverlap, fracOverlapFeature;
	char **chr;
	srInt_64 *start, *stop;
	int *geneid;

	char *nameFeatureTypeColumn, *nameTranscriptIDColumn, *nameGeneIDColumn,*debug_command, *pair_orientations="fr", *temp_dir, *file_name_ptr =NULL, *strand_check_mode = NULL, *extra_column_names = NULL;
	srInt_64 nexons;


	srInt_64 * anno_chr_head, * block_min_start, *block_max_end, *block_end_index;
	char ** anno_chrs, * anno_chr_2ch;
	char * fasta_contigs_fname, *annotation_file_screen_output;
	unsigned char * sorted_strand;

	int minPEDistance, maxPEDistance, isReadSummaryReport, isBothEndRequired, isMultiMappingAllowed, fiveEndExtension, threeEndExtension, minFragmentOverlap, isSplitOrExonicOnly, is_duplicate_ignored, doNotSort, fractionMultiMapping, useOverlappingBreakTie, doJuncCounting, max_M, isRestrictlyNoOvelrapping;
	char * isPEassign,  *is_paired_end_reads_expected;

	int  isGTF, n_input_files=0,is_dual_index, scrna_total_BAM_no;
	char * alias_file_name = NULL, * cmd_rebuilt = NULL, * Rpath = NULL;

	int isMultiOverlapAllowed, isGeneLevel;

	isCVersion = ((argv[0][0]=='C')?1:0);

	isPEassign = argv[4];
	minPEDistance = atoi(argv[5]);
	maxPEDistance = atoi(argv[6]);

	//  isSAM = atoi(argv[7]);
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
		strand_check_mode = argv[12];
	else	strand_check_mode = NULL;
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


	if(isPEDistChecked && 0==isBothEndRequired){
		#ifdef MAKE_STANDALONE
		SUBREADprintf("ERROR: when the '-P' option is specified for checking fragment lengths, the '-B' option must also be specified to require both ends mapped.\n");
		#else
		SUBREADprintf("ERROR: when parameter checkFragLength is set to TRUE, parameter requireBothEndMapped also needs to be set to TRUE.\n");
		#endif
		return -1;
	}

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
	else	isMultiMappingAllowed = 0;
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
	if(thread_number>FC_MAX_THREADS)thread_number=FC_MAX_THREADS;

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
		isSplitOrExonicOnly = atoi(argv[28]);
	else	isSplitOrExonicOnly = 0;

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


	/*if(argc>35) "-S" is depreciated.
		pair_orientations = argv[35];
	else	pair_orientations = "FR";
	*/

	if(argc>36)
		doJuncCounting = atoi(argv[36]);
	else	doJuncCounting = 0;

	fasta_contigs_fname = NULL;
	if(argc>37)
		if(argv[37][0] != 0 && argv[37][0]!=' ')
			fasta_contigs_fname = argv[37];

	if(argc>38)
		max_M = min(65536, max(1,atoi(argv[38])));
	else	max_M = 65536;

	if(argc>39)
		isRestrictlyNoOvelrapping = atoi(argv[39]);
	else	isRestrictlyNoOvelrapping = 0;

	if(argc>40)
		fracOverlap = atof(argv[40]);
	else	fracOverlap= 0.0;
	
	if(argc>41){
		if(strcmp("<use output directory>", argv[41])!=0)temp_dir = argv[41];
		else temp_dir = NULL;
	}
	else	temp_dir = NULL;//	get_temp_dir_from_out(temp_dir, (char *)argv[3]);

	if(argc>42){
		useStdinFile = (atoi(argv[42]) & 1)!=0;
	}else	useStdinFile = 0;
	
	if(argc>43)
		assignReadsToRG = (argv[43][0]=='1');
	else  assignReadsToRG = 0;

	if(argc>44)
		long_read_minimum_length = atoi(argv[44])?1:1999999999;
	else  long_read_minimum_length = 1999999999;
	if(long_read_minimum_length == 1) max_M = 65536; 

	if(long_read_minimum_length < 2 && isPEassign[0]=='1'){
		// this is because long read pairser cannot pair PE reads by names.
		SUBREADputs("ERROR: long read assignment can only be done on single-end mode");
		return -1;
	}

	if(argc>45)
		is_verbose = (argv[45][0]=='1'); 
	else  is_verbose = 0;

	if(argc>46)
		fracOverlapFeature = atof(argv[46]);
	else	fracOverlapFeature = 0.0;

	if(argc>47)
		do_detectionCall = (argv[47][0]=='1'); 
	else  do_detectionCall = 0;

	if(argc>48) max_missing_bases_in_read = atoi(argv[48]);
	else  max_missing_bases_in_read = -1;

	if(argc>49) max_missing_bases_in_feature = atoi(argv[49]);
	else  max_missing_bases_in_feature = -1;
		
	if(argc>50) is_Primary_Alignment_only = atoi(argv[50]);
	else is_Primary_Alignment_only = 0;

	if(argc>51 && argv[51]!=NULL && argv[51][0]!=0 && argv[51][0]!=' ') Rpath = argv[51];
	else Rpath = NULL;

	if(argc>52 && argv[52]!=NULL && argv[52][0]!=0 && argv[52][0]!=' ') extra_column_names = argv[52];
	else extra_column_names = NULL;

	annotation_file_screen_output = NULL;
#ifndef MAKE_STANDALONE
	if(argc>53) annotation_file_screen_output = argv[53];
#endif

	if(argc>54){
		read_shift_type = -1;
		if(strcmp(argv[54], "upstream")==0)read_shift_type = READ_SHIFT_UPSTREAM;
		if(strcmp(argv[54], "downstream")==0) read_shift_type = READ_SHIFT_DOWNSTREAM;
		if(strcmp(argv[54], "left")==0) read_shift_type = READ_SHIFT_LEFT;
		if(strcmp(argv[54], "right")==0) read_shift_type = READ_SHIFT_RIGHT;
	} else read_shift_type = READ_SHIFT_UPSTREAM;

	if(argc>55) read_shift_size = atoi(argv[55]);
	else read_shift_size = 0;

	if(argc>58 && strlen(argv[58])>0 && argv[58][0]!=' ') is_paired_end_reads_expected = argv[58];
	else is_paired_end_reads_expected = "0";

	if(argc>63) is_dual_index = (strcmp(argv[63],"Dual-Index")==0);
	else is_dual_index =0; 

	if(argc>64) scrna_total_BAM_no = atoi(argv[64]);
	else scrna_total_BAM_no =0; 

	if(argc > 65)
		nameTranscriptIDColumn = argv[65];
	else	nameTranscriptIDColumn = "transcript_id";

	if(read_shift_size<0){
		SUBREADprintf("ERROR: why the value for read_shift_size is negative?\n");
		return -1;
	}

	if(read_shift_type<0){
		SUBREADprintf("ERROR: why the value for read_shift_type is %s?\n", argv[54]);
		return -1;
	}

	if(SAM_pairer_warning_file_open_limit()) return -1;
	if(strand_check_mode != NULL && Input_Files_And_Strand_Mode_Pair(argv[2],strand_check_mode)) return -1;
	if(extra_column_names){
		if(!isGTF){
			SUBREADputs("ERROR: only GTF files contain additional attributes");
			return -1;
		}
		int xk1, total_cols =1;
		for(xk1=0; extra_column_names[xk1]; xk1++)
			if(extra_column_names[xk1] == ';' || extra_column_names[xk1]==',' || extra_column_names[xk1]=='\t'){
				extra_column_names[xk1]='\t';
				total_cols ++;
			}
		if(total_cols>MAX_EXTRA_COLS){
			SUBREADprintf("ERROR: there are more than %d additional attributes required\n", MAX_EXTRA_COLS);
			return -1;
		}
	}

	fc_thread_global_context_t global_context;

	fc_thread_init_global_context(& global_context, FEATURECOUNTS_BUFFER_SIZE, thread_number, MAX_LINE_LENGTH, minPEDistance, maxPEDistance,isGeneLevel, isMultiOverlapAllowed, strand_check_mode, (char *)argv[3] , isReadSummaryReport, isBothEndRequired, isChimericDisallowed, isPEDistChecked, nameFeatureTypeColumn, nameGeneIDColumn, nameTranscriptIDColumn, minMappingQualityScore,isMultiMappingAllowed, 0, alias_file_name, cmd_rebuilt, isInputFileResortNeeded, feature_block_size, isCVersion, fiveEndExtension, threeEndExtension , minFragmentOverlap, isSplitOrExonicOnly, reduce_5_3_ends_to_one, debug_command, is_duplicate_ignored, doNotSort, fractionMultiMapping, useOverlappingBreakTie, pair_orientations, doJuncCounting, max_M, isRestrictlyNoOvelrapping, fracOverlap, temp_dir, useStdinFile, assignReadsToRG, long_read_minimum_length, is_verbose, fracOverlapFeature, do_detectionCall, max_missing_bases_in_read, max_missing_bases_in_feature, is_Primary_Alignment_only, Rpath, extra_column_names, annotation_file_screen_output, read_shift_type, read_shift_size, is_dual_index );

	fc_thread_init_input_files( & global_context, argv[2], &file_name_ptr );

	if( print_FC_configuration(&global_context, argv[1], file_name_ptr, argv[3], global_context.is_SAM_file, isGTF, & n_input_files, isReadSummaryReport, is_paired_end_reads_expected, isPEassign) )
		return -1;
	// Loading the annotations.
	// Nothing is done if the annotation does not exist.
	fc_feature_info_t * loaded_features;
	print_in_box(84,0,0,"Load annotation file %s %c[0m...", get_short_fname(argv[1]), CHAR_ESC);
	nexons = load_feature_info(&global_context,argv[1], isGTF?FILE_TYPE_GTF:FILE_TYPE_RSUBREAD, &loaded_features);
	if(nexons<1){
		if(nexons >= -1) SUBREADprintf("Failed to open the annotation file %s, or its format is incorrect, or it contains no '%s' features.\n",argv[1], nameFeatureTypeColumn);
		return -1;
	}

	sort_feature_info(&global_context, nexons, loaded_features, &chr, &geneid, &start, &stop, &sorted_strand, &anno_chr_2ch, &anno_chrs, &anno_chr_head, & block_end_index, & block_min_start, & block_max_end);
	global_context.lineno_2_sortedno_tab = NULL;

	print_in_box(80,0,0,"   Meta-features : %d", global_context . gene_name_table -> numOfElements);
	print_in_box(80,0,0,"   Chromosomes/contigs : %d", global_context . exontable_nchrs);

	print_in_box(80,0,0,"");

	if(fasta_contigs_fname){
		print_in_box(80,0,0,"Load FASTA contigs from %s...", get_short_fname(fasta_contigs_fname));
		global_context.fasta_contigs = malloc(sizeof(fasta_contigs_t));
		int ret_fq = read_contig_fasta(global_context.fasta_contigs, fasta_contigs_fname);
		if(ret_fq){
			print_in_box(80,0,0,"   WARNING unable to open the FASTA file.");
			print_in_box(80,0,0,"");
			free(global_context.fasta_contigs);
			global_context.fasta_contigs = NULL;
		}else{
			print_in_box(80,0,0,"   %lu contigs were loaded", global_context.fasta_contigs -> contig_table -> numOfElements);
			print_in_box(80,0,0,"");
		}
	}else	global_context.fasta_contigs = NULL;
	

	global_context.exontable_exons = nexons;
	unsigned int x1, total_written_coulmns=0;




	char * tmp_pntr = NULL, *tmp_smode_ptr = NULL;
	char * strand_mode_list = strdup(global_context.strand_check_mode);
	char * file_list_used = malloc(strlen(file_name_ptr)+1);
	char * file_list_used2 = malloc(strlen(file_name_ptr)+1);
	char * is_unique = malloc(strlen(file_name_ptr)+1);
	strcpy(file_list_used, file_name_ptr);
	for(x1 = 0;;x1++){
		char * test_fn = strtok_r(x1?NULL:file_list_used, FC_FLIST_SPLITOR, &tmp_pntr);
		if(NULL == test_fn) break; 
		char * short_fname = get_short_fname(test_fn);
		strcpy(file_list_used2, file_name_ptr);

		is_unique[x1]=1;
		char * loop_ptr = NULL;
		int x2;
		for(x2 = 0;;x2++){
			char * test_loopfn = strtok_r(x2?NULL:file_list_used2, FC_FLIST_SPLITOR, &loop_ptr);
			if(NULL == test_loopfn) break;
			if(x1==x2)continue;

			char * short_loop_fname = get_short_fname(test_loopfn);

			if(strcmp(short_loop_fname, short_fname)==0) {
				is_unique[x1] = 0;
				break;
			}
		}
	}
	free(file_list_used2);

	tmp_pntr = NULL;
	strcpy(file_list_used, file_name_ptr);
	char * next_fn = strtok_r(file_list_used, FC_FLIST_SPLITOR, &tmp_pntr);
	char * next_strand_mode = strtok_r(strand_mode_list, ".", &tmp_smode_ptr);
	int one_single_strand_mode = -1;
	if(NULL == strstr( global_context.strand_check_mode, "." )){
		one_single_strand_mode = next_strand_mode[0] - '0';
		assert(one_single_strand_mode >= 0 && one_single_strand_mode < 3);
	}

	ArrayList * table_columns = ArrayListCreate(n_input_files+1);
	ArrayList * table_column_names = ArrayListCreate(n_input_files+1);
	ArrayList * read_counters = ArrayListCreate(n_input_files+1);
	ArrayListSetDeallocationFunction(table_columns, free);
	ArrayListSetDeallocationFunction(table_column_names, free);
	ArrayListSetDeallocationFunction(read_counters, free);
	
	ArrayList * junction_global_table_list = NULL;
	ArrayList * splicing_global_table_list = NULL;

	if(global_context.do_junction_counting){
		junction_global_table_list = ArrayListCreate(n_input_files+1);
		splicing_global_table_list = ArrayListCreate(n_input_files+1);
		ArrayListSetDeallocationFunction(junction_global_table_list, (void (*)(void *))HashTableDestroy);
		ArrayListSetDeallocationFunction(splicing_global_table_list, (void (*)(void *))HashTableDestroy);
	}

	int ret_int = 0;

#ifdef MAKE_STANDALONE
	#define NO_SORT_OPTION_NAME "donotsort"
#else
	#define NO_SORT_OPTION_NAME "autosort"
#endif

	for(x1 = 0;;x1++){
		int orininal_isPE = global_context.is_paired_end_mode_assign;
		if(next_fn==NULL || strlen(next_fn)<1 || global_context.disk_is_full || ret_int) break;
		int this_file_isPEassign = isPEassign[1]?isPEassign[x1] == '1' :(isPEassign[0]=='1');
		int this_file_isPEexpected = is_paired_end_reads_expected[1]?is_paired_end_reads_expected[x1]=='1' :(is_paired_end_reads_expected[0]=='1');
		global_context.is_paired_end_reads_expected = this_file_isPEexpected;
		global_context.is_paired_end_mode_assign = this_file_isPEassign;
		if(global_context.do_not_sort && 0==this_file_isPEassign){
			print_in_box(80,0,0,"      WARNING the %s option is ignored when single-end reads", NO_SORT_OPTION_NAME);
			print_in_box(80,0,0,"              are being counted.");
		}

		read_count_type_t * column_numbers = calloc(nexons, sizeof(read_count_type_t));
		HashTable * junction_global_table = NULL;
		HashTable * splicing_global_table = NULL;

		strcpy(global_context.input_file_name, next_fn);
		strcpy(global_context.raw_input_file_name, next_fn);
		global_context.this_input_number = x1;
		global_context.input_file_unique = is_unique[x1];
		global_context.input_file_short_name = get_short_fname(next_fn);
		if(strstr( global_context.strand_check_mode, "." )){
			global_context.is_strand_checked = next_strand_mode[0]-'0';
			assert(global_context.is_strand_checked >=0 && global_context.is_strand_checked <=2);
		}else global_context.is_strand_checked = one_single_strand_mode;
		global_context.redo=0;

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

		HashTable * merged_RG_table = NULL;
		if(global_context.assign_reads_to_RG){
			merged_RG_table = HashTableCreate(97);
			HashTableSetHashFunction(merged_RG_table,HashTableStringHashFunction);
			HashTableSetDeallocationFunctions(merged_RG_table, NULL, free); // the names are put into the column_names table, but the 4-pointer arrays are not used anymore.
			HashTableSetKeyComparisonFunction(merged_RG_table, fc_strcmp_chro);
		}
		
		fc_read_counters * my_read_counter = calloc(1, sizeof(fc_read_counters));
		global_context.is_read_details_out = isReadSummaryReport;
		global_context.max_M = max_M;

		ret_int = ret_int || readSummary_single_file(& global_context, column_numbers, nexons, geneid, chr, start, stop, sorted_strand, anno_chr_2ch, anno_chrs, anno_chr_head, block_end_index, block_min_start, block_max_end, my_read_counter, junction_global_table, splicing_global_table, merged_RG_table, loaded_features);
		if(global_context.disk_is_full){
			SUBREADprintf("ERROR: disk is full. Please check the free space in the output directory.\n");
		}
		if(ret_int!=0){
			// give up this file.
			if(global_context.do_junction_counting){
				HashTableDestroy(junction_global_table);
				HashTableDestroy(splicing_global_table);
			}
			free(column_numbers);
		} else {
			// finished

			char * mem_file_name = memstrcpy(next_fn);
			if(!global_context.assign_reads_to_RG){
				ArrayListPush(table_columns, column_numbers);
				ArrayListPush(table_column_names, mem_file_name);
				ArrayListPush(read_counters, my_read_counter);
				if(global_context.do_junction_counting){
					ArrayListPush(junction_global_table_list,junction_global_table);
					ArrayListPush(splicing_global_table_list,splicing_global_table);
				}
			}
			
			if(global_context.assign_reads_to_RG){
				int rgcur;
				char * rg_name = global_context.RGnames_set;
				for(rgcur = 0; rgcur < global_context.RGnames_ptr+1;  rgcur ++){
					if(global_context.RGnames_set[rgcur] == '\t'||global_context.RGnames_set[rgcur] == '\0'){
						global_context.RGnames_set[rgcur] = 0;
						int rg_name_len = strlen(rg_name);
						if(rg_name_len > 0){
	//						SUBREADprintf("GET 4Tab:'%s'\n", rg_name);
							void ** tab4 = HashTableGet(merged_RG_table, rg_name);
							int file_len = strlen(mem_file_name);
		
							char * rg_file_name = malloc(rg_name_len + 3 + file_len);
							SUBreadSprintf(rg_file_name, rg_name_len + 3 + file_len, "%s:%s", mem_file_name, rg_name);
							
							ArrayListPush(table_column_names, rg_file_name);
							ArrayListPush(table_columns, tab4[0]);
							ArrayListPush(read_counters, tab4[1]);
							if(global_context.do_junction_counting){
								ArrayListPush(junction_global_table_list,tab4[2]);
								ArrayListPush(splicing_global_table_list,tab4[3]);
							}
							rg_name = global_context.RGnames_set + rgcur + 1;
						}
					}
				}
				free(mem_file_name);
			}
			total_written_coulmns ++;
		}
		global_context.is_paired_end_mode_assign = orininal_isPE;
		next_fn = strtok_r(NULL, FC_FLIST_SPLITOR, &tmp_pntr);

		if(strstr( global_context.strand_check_mode, "." )) next_strand_mode = strtok_r(NULL, ".", &tmp_smode_ptr);
		if(global_context.assign_reads_to_RG) free(global_context.RGnames_set);
		if(merged_RG_table) HashTableDestroy(merged_RG_table);
	}

	free(file_list_used);
	free(is_unique);

	if(global_context.is_input_bad_format){
	//	SUBREADprintf("\nEEROR: The program has to terminate and no counting file is generated.\n\n");
	}else if(!(global_context.disk_is_full||ret_int)){
		print_in_box(80,0,0,"Write the final count table.");
		if(isGeneLevel){
			char ** sorted_extra_columns = NULL;
			if(global_context.reported_extra_columns != NULL){
				sorted_extra_columns = malloc(sizeof(char**) * nexons);
				int ii;
				for(ii = 0; ii < nexons; ii++){
					sorted_extra_columns[loaded_features[ii].sorted_order] = loaded_features[ii].extra_columns;
					//SUBREADprintf("SSMQ: %d = %s\n", loaded_features[ii].sorted_order, loaded_features[ii].extra_columns);
				}
			}

			fc_write_final_gene_results(&global_context, geneid, chr, start, stop, sorted_strand, sorted_extra_columns, argv[3], nexons,  table_columns, table_column_names, loaded_features, isCVersion);

			if(sorted_extra_columns) free(sorted_extra_columns);
		} else
			fc_write_final_results(&global_context, argv[3], nexons, table_columns, table_column_names, loaded_features, isCVersion);
	}
	if(global_context.do_junction_counting && global_context.is_input_bad_format == 0 && !(ret_int||global_context.disk_is_full)){
		print_in_box(80,0,0,"Write the junction count table.");
		fc_write_final_junctions(&global_context, argv[3], table_column_names, junction_global_table_list, splicing_global_table_list);
	}

	if(global_context.is_input_bad_format == 0 && !(ret_int||global_context.disk_is_full)){
		print_in_box(80,0,0,"Write the read assignment summary.");
		fc_write_final_counts(&global_context, argv[3], table_column_names,  read_counters, isCVersion);
	}

	ArrayListDestroy(table_columns);
	ArrayListDestroy(table_column_names);
	ArrayListDestroy(read_counters);
	if(global_context.do_junction_counting){
		ArrayListDestroy(junction_global_table_list);
		ArrayListDestroy(splicing_global_table_list);
	}
	free(file_name_ptr);

	if(global_context.is_input_bad_format == 0) print_FC_results(&global_context, (char *)argv[3]/*out file name*/);
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

	if(global_context.read_details_out_FP) fclose(global_context. read_details_out_FP);
	HashTableDestroy(global_context.gene_name_table);
	HashTableDestroy(global_context.GCcontent_table);
	free(global_context.gene_name_array);

	HashTableDestroy(global_context.exontable_chro_table);
	if(global_context.fasta_contigs){
		destroy_contig_fasta(global_context.fasta_contigs);
		free(global_context.fasta_contigs);
	}
	if(global_context.BAM_chros_to_anno_table)
		HashTableDestroy(global_context.BAM_chros_to_anno_table);
	if(global_context.do_junction_counting){
		HashTableDestroy(global_context.junction_genebody_table);
		HashTableDestroy(global_context.junction_transcript_table);
		HashTableDestroy(global_context.junction_GenebodyTree_table);
		HashTableDestroy(global_context.junction_ExonTree_table);
	}


	free(global_context.unistr_buffer_space);

	if(global_context.reported_extra_columns){
		for(bucket = 0; bucket < nexons; bucket++)
			free(loaded_features[bucket].extra_columns);
	}
	if(global_context.lineno_2_sortedno_tab)HashTableDestroy(global_context.lineno_2_sortedno_tab);

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
	free(strand_mode_list);

	return total_written_coulmns?0:-1;
}

int readSummary_single_file(fc_thread_global_context_t * global_context, read_count_type_t * column_numbers, srInt_64 nexons,  int * geneid, char ** chr, srInt_64 * start, srInt_64 * stop, unsigned char * sorted_strand, char * anno_chr_2ch, char ** anno_chrs, srInt_64 * anno_chr_head, srInt_64 * block_end_index, srInt_64 * block_min_start , srInt_64 * block_max_end, fc_read_counters * my_read_counter, HashTable * junction_global_table, HashTable * splicing_global_table, HashTable * merged_RG_table, fc_feature_info_t * loaded_features)
{
	int read_length = 0;
	int is_first_read_PE=0;
	char * line = (char*)calloc(MAX_LINE_LENGTH, 1);
	char * file_str = "";

	int file_probe = is_certainly_bam_file(global_context->input_file_name, &is_first_read_PE, NULL);
		
		// a Singel-end SAM/BAM file cannot be assigned as a PE SAM/BAM file;
		// but a PE SAM/BAM file may be assigned as a SE file if the user wishes to do so.

	global_context->is_SAM_file = 1;
	if(file_probe == 1) global_context->is_SAM_file = 0;
	global_context->is_mixed_PE_SE = 0;
	global_context->any_reads_are_PE = 0;
	global_context -> start_time = miltime();

	file_str = "SAM";
	if(file_probe == 1) file_str = "BAM" ;
	if(file_probe == -1) file_str = "Unknown";

	if(!global_context->redo)
	{
		print_in_box(80,0,0,"Process %s file %s...", file_str, global_context -> use_stdin_file? "<STDIN>":get_short_fname(global_context->input_file_name));
		if(global_context->is_strand_checked)
			print_in_box(80,0,0,"   Strand specific : %s", global_context->is_strand_checked==1?"stranded":"reversely stranded");
	}
		
	// Open the SAM/BAM file
	// Nothing is done if the file does not exist.

	int thread_error = fc_thread_start_threads(global_context, nexons, geneid, chr, start, stop, sorted_strand, anno_chr_2ch, anno_chrs, anno_chr_head, block_end_index, block_min_start , block_max_end, read_length);
	if(thread_error) return -1;
	fc_thread_wait_threads(global_context);
	if(global_context -> is_paired_end_reads_expected && !global_context -> any_reads_are_PE){
		SUBREADprintf("ERROR: No paired-end reads were detected in paired-end read library : %s\n", global_context -> input_file_name);
		global_context -> is_input_bad_format=1;
		return -1;
	}

	srInt_64 nreads_mapped_to_exon = 0;
	fc_thread_merge_results(global_context, column_numbers , &nreads_mapped_to_exon, my_read_counter, junction_global_table, splicing_global_table, merged_RG_table, loaded_features, nexons);
	fc_thread_destroy_thread_context(global_context);

	if(global_context -> sambam_chro_table) free(global_context -> sambam_chro_table);
	global_context -> sambam_chro_table = NULL;

	free(line);
	if(global_context -> is_input_bad_format) return -1;
	return 0;
}


#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int feature_count_main(int argc, char ** argv)
#endif
{
	char * Rargv[65];
	char annot_name[MAX_FILE_NAME_LENGTH];
	char temp_dir[MAX_FILE_NAME_LENGTH];
	char * out_name = malloc(MAX_FILE_NAME_LENGTH);
	char * fasta_contigs_name = malloc(MAX_FILE_NAME_LENGTH);
	char * alias_file_name = malloc(MAX_FILE_NAME_LENGTH);
	char * Rpath = malloc(MAX_FILE_NAME_LENGTH);

	int cmd_rebuilt_size = 2000;
	char * cmd_rebuilt = malloc(cmd_rebuilt_size);
	char max_M_str[8];
	char nameFeatureTypeColumn[2000];
	char nameGeneIDColumn[66];
	char nameTranscriptIDColumn[66];
	int min_qual_score = 0;
	int min_dist = 50;
	int max_dist = 600;
	int read_shift_size = 0;
	char debug_command[15];
	char max_missing_bases_in_read_str[15];
	char max_missing_bases_in_feature_str[15];
	char min_dist_str[15];
	char max_dist_str[15];
	char read_shift_size_str[15];
	char read_shift_type[15];
	char min_qual_score_str[15];
	char feature_block_size_str[15];
	char * Strand_Sensitive_Str = "0";
	char * old_zero_smode = Strand_Sensitive_Str;
	char strFeatureFracOverlap[20];
	char Pair_Orientations[3];
	char * extra_column_names = NULL;
	char * very_long_file_names;
	char is_paired_end_reads_expected[2];
	int is_Input_Need_Reorder = 0;
	int is_PE = 0;
	int is_SAM = 1;
	int is_primary_alignment_only = 0;
	int is_GeneLevel = 1;
	int is_Overlap = 0;
	int is_Both_End_Mapped = 0;
	int is_Restrictedly_No_Overlap = 0;
	int feature_block_size = 14;
	int is_ReadSummary_Report = 0;
	int is_Chimeric_Disallowed = 0;
	int is_PE_Dist_Checked = 0;
	int is_Multi_Mapping_Allowed = 0;
	int is_Split_or_Exonic_Only = 0;
	int is_duplicate_ignored = 0;
	int assign_reads_to_RG = 0;
	int do_not_sort = 0;
	int do_junction_cnt = 0;
	int do_detection_call = 0;
	int reduce_5_3_ends_to_one = 0;
	int use_fraction_multimapping = 0;
	int threads = 1;
	int isGTF = 1;
	int only_want_largest_overlap = 0;
	char nthread_str[4];
	int option_index = 0;
	int max_missing_bases_in_feature = -1;
	int max_missing_bases_in_read = -1;
	int c;
	int very_long_file_names_size = 200;
	int fiveEndExtension = 0, threeEndExtension = 0, minFragmentOverlap = 1;
	float fracOverlap = 0.0, fracOverlapFeature = 0.0;
	int std_input_output_mode = 0, long_read_mode = 0, is_verbose = 0;
	char strFiveEndExtension[11], strThreeEndExtension[11], strMinFragmentOverlap[11], fracOverlapStr[20], std_input_output_mode_str[16], long_read_mode_str[16];
	very_long_file_names = malloc(very_long_file_names_size);
	very_long_file_names [0] = 0;
	fasta_contigs_name[0]=0;
	is_paired_end_reads_expected[0]='0';
	is_paired_end_reads_expected[1]='\0';

	alias_file_name[0]=0;
	debug_command[0] = 0;

	strcpy(read_shift_type,"upstream");
	strcpy(nameFeatureTypeColumn,"exon");
	strcpy(nameGeneIDColumn,"gene_id");
	strcpy(nameTranscriptIDColumn,"transcript_id");
	strcpy(temp_dir, "<use output directory>");
	annot_name[0]=0;out_name[0]=0;Rpath[0]=0;
	

	cmd_rebuilt[0]=0;
	for(c = 0; c<argc;c++)
	{
		if(strlen(cmd_rebuilt) + 1000 > cmd_rebuilt_size)
		{
			cmd_rebuilt_size*=2;
			cmd_rebuilt = realloc(cmd_rebuilt, cmd_rebuilt_size);
		}
		SUBreadSprintf(cmd_rebuilt+strlen(cmd_rebuilt), cmd_rebuilt_size-strlen(cmd_rebuilt), "\"%s\" ", argv[c]);
	}

	optind=0;
	opterr=1;
	optopt=63;
	strcpy(max_M_str, "65536");
	strcpy(Pair_Orientations,"fr");

	while ((c = getopt_long (argc, argv, "G:A:g:t:T:o:a:d:D:LQ:pbF:fs:S:CBJPMOR:v?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'S':
				/*
				if(strlen(optarg)!=2 || (strcmp(optarg, "ff")!=0 && strcmp(optarg, "rf")!=0 && strcmp(optarg, "fr")!=0)){
					SUBREADprintf("The order parameter can only be ff, fr or rf.\n");
					print_usage();
					return -1;
				}
				Pair_Orientations[0]=(optarg[0]=='r'?'r':'f');
				Pair_Orientations[1]=(optarg[1]=='f'?'f':'r');
				Pair_Orientations[2]=0;
				*/
				SUBREADprintf("The \"-S\" option has been depreciated.\n");

				break;
			case 'G':
				strcpy(fasta_contigs_name , optarg);
				break;
			case 'J':
				do_junction_cnt = 1;
				break;
			case 'A':
				strcpy(alias_file_name, optarg);
				break;
			case 'M':
				is_Multi_Mapping_Allowed = 1;
				break;
			case 'v':
				core_version_number("featureCounts");
				return 0;
			case 'Q':
				if(!is_valid_digit_range(optarg, "Q", 0 , 255))
					STANDALONE_exit(-1);

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
				if(!is_valid_digit_range(optarg, "T", 1, FC_MAX_THREADS))
					STANDALONE_exit(-1);

				threads = atoi(optarg);
				break;
			case 'd':
				if(!is_valid_digit(optarg, "d"))
					STANDALONE_exit(-1);

				min_dist = atoi(optarg);
				break;
			case 'D':
				if(!is_valid_digit(optarg, "D"))
					STANDALONE_exit(-1);

				max_dist = atoi(optarg);
				break;
			case 'p':
				is_paired_end_reads_expected[0]='1';
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
				if(strcmp(optarg, "SAM")==0) is_ReadSummary_Report = FILE_TYPE_SAM;
				else if(strcmp(optarg, "BAM")==0) is_ReadSummary_Report = FILE_TYPE_BAM;
				else if(strcmp(optarg, "CORE")==0) is_ReadSummary_Report = FILE_TYPE_RSUBREAD;
				else{
					SUBREADprintf("\nERROR: unknown output format: '%s'\n\n", optarg);
					STANDALONE_exit(-1);
				}
				break;
			case 's':
				Strand_Sensitive_Str = strdup(optarg);
				int xx;
				for(xx =0; Strand_Sensitive_Str[xx]!='\0'; xx++) if(Strand_Sensitive_Str[xx]==',') Strand_Sensitive_Str[xx]='.';
				break;
//			case 'i':
//				term_strncpy(sam_name, optarg,299);
//				break;
			case 'o':
				term_strncpy(out_name, optarg,MAX_FILE_NAME_LENGTH-1);
				break;
			case 'a':
				term_strncpy(annot_name, optarg,MAX_FILE_NAME_LENGTH-1);
				break;
			case 'L':
				long_read_mode = 1;
				break;
			case 0 :	// long options

				if(strcmp("countReadPairs", long_options[option_index].name)==0){
					is_PE=1;
				}

				if(strcmp("primary", long_options[option_index].name)==0)
				{
					is_primary_alignment_only = 1;
				}

				if(strcmp("readExtension5", long_options[option_index].name)==0)
				{
					if(!is_valid_digit_range(optarg, "readExtension5", 0, 0x7fffffff))
						STANDALONE_exit(-1);
					fiveEndExtension = atoi(optarg);
					fiveEndExtension = max(0, fiveEndExtension);
				}

				if(strcmp("readExtension3", long_options[option_index].name)==0)
				{
					if(!is_valid_digit_range(optarg, "readExtension3", 0, 0x7fffffff))
						STANDALONE_exit(-1);
					threeEndExtension = atoi(optarg);
					threeEndExtension = max(0, threeEndExtension);
				}

				if(strcmp("fracOverlap", long_options[option_index].name)==0)
				{
					if(!is_valid_float(optarg, "fracOverlap"))
						STANDALONE_exit(-1);
					fracOverlap = atof(optarg);
				}


				if(strcmp("fracOverlapFeature", long_options[option_index].name)==0)
				{
					if(!is_valid_float(optarg, "fracOverlapFeature"))
						STANDALONE_exit(-1);
					fracOverlapFeature = atof(optarg);
				}

				if(strcmp("nonOverlapFeature", long_options[option_index].name)==0){
					if(!is_valid_digit_range(optarg, "nonOverlapFeature", 0, 0x7fffffff))
						STANDALONE_exit(-1);
					max_missing_bases_in_feature = atoi(optarg);
				}

				if(strcmp("nonOverlap", long_options[option_index].name)==0){
					if(!is_valid_digit_range(optarg, "nonOverlap", 0, 0x7fffffff))
						STANDALONE_exit(-1);
					max_missing_bases_in_read = atoi(optarg);
				}

				if(strcmp("extraAttributes", long_options[option_index].name)==0)
				{
					extra_column_names = strdup(optarg);
				}

				if(strcmp("Rpath", long_options[option_index].name)==0)
				{
					strcpy(Rpath, optarg);
				}

				if(strcmp("minOverlap", long_options[option_index].name)==0)
				{
					if(!is_valid_digit(optarg, "minOverlap"))
						STANDALONE_exit(-1);
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
				if(strcmp("tmpDir", long_options[option_index].name)==0){
					strcpy(temp_dir, optarg);
				}
				if(strcmp("maxMOp", long_options[option_index].name)==0){
					if(!is_valid_digit_range(optarg, "maxMOp", 1 , 65555))
						STANDALONE_exit(-1);
					strcpy(max_M_str, optarg);
				}
				if(strcmp("read2pos", long_options[option_index].name)==0)
				{
					if(optarg[0]=='3')
						reduce_5_3_ends_to_one = REDUCE_TO_3_PRIME_END;
					else if(optarg[0]=='5')
						reduce_5_3_ends_to_one = REDUCE_TO_5_PRIME_END;
					else{
						SUBREADprintf("Invalide parameter to the --read2pos option: %s\n", optarg);
						STANDALONE_exit(-1);
					}		
				}				

				if(strcmp("largestOverlap", long_options[option_index].name)==0)
				{
					only_want_largest_overlap = 1;
				}

				if(strcmp("detectionCall", long_options[option_index].name)==0)
				{
					do_detection_call = 1;
				}

				if(strcmp("donotsort", long_options[option_index].name)==0)
				{
					do_not_sort = 1;
				}

				if(strcmp("readShiftSize", long_options[option_index].name)==0)
				{
					if(!is_valid_digit_range(optarg, "readShiftSize", 1 , 0x7fffffff))
						STANDALONE_exit(-1);
					read_shift_size = atoi(optarg);
				}

				if(strcmp("readShiftType", long_options[option_index].name)==0)
				{
					if(strcmp(optarg,"upstream")!=0 && strcmp(optarg,"downstream")!=0 && strcmp(optarg,"left")!=0 && strcmp(optarg,"right")!=0){
						SUBREADprintf("Error: the readShiftType parameter can only be 'upstream', 'downstream', 'left' or 'right'\n");
						STANDALONE_exit(-1);
					}
					strcpy(read_shift_type, optarg);
				}

				if(strcmp("splitOnly", long_options[option_index].name)==0)
				{
					if(is_Split_or_Exonic_Only == 2) {
						SUBREADprintf("Error: You can not specify both splitOnly and nonSplitOnly\n");
						return -1;
					}
					is_Split_or_Exonic_Only = 1;
				}

				if(strcmp("restrictedlyNoOverlap", long_options[option_index].name)==0)
				{
					is_Restrictedly_No_Overlap = 1;
				}
				if(strcmp("nonSplitOnly", long_options[option_index].name)==0)
				{
					if(is_Split_or_Exonic_Only == 1) {
						SUBREADprintf("Error: You can not specify both splitOnly and nonSplitOnly\n");
						return -1;
					}
					is_Split_or_Exonic_Only = 2;
				}
				
				if(strcmp("verbose", long_options[option_index].name)==0){
					is_verbose = 1;
				}

				if(strcmp("byReadGroup", long_options[option_index].name)==0){
					assign_reads_to_RG = 1;
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

	if(out_name[0]==0 || annot_name[0]==0)
	{
		print_usage();
		return -1;
	}

	for(; optind < argc; optind++)
	{
		int curr_strlen = strlen(very_long_file_names);
		if( very_long_file_names_size - curr_strlen < MAX_FILE_NAME_LENGTH+1)
		{
			very_long_file_names_size *=2;
			//printf("CL=%d ; NS=%d\n", curr_strlen , very_long_file_names_size);
			very_long_file_names=realloc(very_long_file_names , very_long_file_names_size);
		}

		strcat(very_long_file_names, argv[optind]);
		strcat(very_long_file_names, FC_FLIST_SPLITOR);
	}

	very_long_file_names[strlen(very_long_file_names)-1]=0;
	std_input_output_mode = (strcmp(very_long_file_names, "") == 0?1:0);

	SUBreadSprintf(strFiveEndExtension, 11, "%d", fiveEndExtension);
	SUBreadSprintf(strThreeEndExtension,11, "%d", threeEndExtension);
	SUBreadSprintf(strMinFragmentOverlap,11, "%d", minFragmentOverlap);
	SUBreadSprintf(nthread_str,11,"%d", threads);
	SUBreadSprintf(min_dist_str,11,"%d",min_dist);
	SUBreadSprintf(max_dist_str,11,"%d",max_dist);
	SUBreadSprintf(min_qual_score_str,11,"%d", min_qual_score);
	SUBreadSprintf(feature_block_size_str,11,"%d", feature_block_size);
	SUBreadSprintf(fracOverlapStr,20, "%g", fracOverlap);
	SUBreadSprintf(std_input_output_mode_str,11,"%d",std_input_output_mode);
	SUBreadSprintf(long_read_mode_str,11, "%d", long_read_mode);
	SUBreadSprintf(strFeatureFracOverlap,20, "%g", fracOverlapFeature);
	SUBreadSprintf(max_missing_bases_in_feature_str,11, "%d", max_missing_bases_in_feature);
	SUBreadSprintf(max_missing_bases_in_read_str,11, "%d", max_missing_bases_in_read);
	SUBreadSprintf(read_shift_size_str,11, "%d", read_shift_size);

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
	Rargv[13] = is_ReadSummary_Report == 0 ? "0":(is_ReadSummary_Report == FILE_TYPE_RSUBREAD?"10":(is_ReadSummary_Report == FILE_TYPE_BAM?"500":"50"));
	Rargv[14] = is_Both_End_Mapped?"1":"0";
	Rargv[15] = is_Chimeric_Disallowed?"1":"0";
	Rargv[16] = is_PE_Dist_Checked?"1":"0";
	Rargv[17] = nameFeatureTypeColumn;
	Rargv[18] = nameGeneIDColumn;
	Rargv[19] = min_qual_score_str;
	Rargv[20] = is_Multi_Mapping_Allowed?"1":"0";
	Rargv[21] = alias_file_name;
	Rargv[22] = cmd_rebuilt;
	Rargv[23] = is_Input_Need_Reorder?"1":"0";
	Rargv[24] = feature_block_size_str;
	Rargv[25] = strFiveEndExtension;
	Rargv[26] = strThreeEndExtension;
	Rargv[27] = strMinFragmentOverlap;
	Rargv[28] = is_Split_or_Exonic_Only == 1?"1":(is_Split_or_Exonic_Only ==  2 ? "2":"0");
	Rargv[29] = (reduce_5_3_ends_to_one == 0?"0":(reduce_5_3_ends_to_one==REDUCE_TO_3_PRIME_END?"3":"5"));
	Rargv[30] = debug_command;
	Rargv[31] = is_duplicate_ignored?"1":"0";
	Rargv[32] = do_not_sort?"1":"0";
	Rargv[33] = use_fraction_multimapping?"1":"0";
	Rargv[34] = only_want_largest_overlap?"1":"0";
	Rargv[35] = Pair_Orientations;
	Rargv[36] = do_junction_cnt?"1":"0";
	Rargv[37] = fasta_contigs_name;
	Rargv[38] = max_M_str;
	Rargv[39] = is_Restrictedly_No_Overlap?"1":"0"; 
	Rargv[40] = fracOverlapStr;
	Rargv[41] = temp_dir;
	Rargv[42] = std_input_output_mode_str;
	Rargv[43] = assign_reads_to_RG?"1":"0"; 
	Rargv[44] = long_read_mode_str;
	Rargv[45] = is_verbose?"1":"0";
	Rargv[46] = strFeatureFracOverlap;
	Rargv[47] = do_detection_call?"1":"0";
	Rargv[48] = max_missing_bases_in_read_str;
	Rargv[49] = max_missing_bases_in_feature_str;
	Rargv[50] = is_primary_alignment_only?"1":"0";
	Rargv[51] = Rpath;
	Rargv[52] = extra_column_names;
	Rargv[54] = "NA"; // C featureCounts dosn't need the display_annotation_name.
	Rargv[54] = read_shift_type;
	Rargv[55] = read_shift_size_str;
	Rargv[56] = "NA";
	Rargv[57] = "NA";
	Rargv[58] = is_paired_end_reads_expected;
	Rargv[59] = "0";
	Rargv[60] = "3";
	Rargv[61] = "0";
	Rargv[62] = "-1";
	Rargv[63] = "Single-Index";
	Rargv[64] = nameTranscriptIDColumn;

	int retvalue = -1;
	if(is_ReadSummary_Report && (std_input_output_mode & 1)==1) SUBREADprintf("ERROR: no detailed assignment results can be written when the input is from STDIN. Please remove the '-R' option.\n");
	else retvalue = readSummary(65, Rargv);

	free(very_long_file_names);
	free(out_name);
	free(alias_file_name);
	free(fasta_contigs_name);
	if(old_zero_smode != Strand_Sensitive_Str)free(Strand_Sensitive_Str);
	free(cmd_rebuilt);
	free(Rpath);
	if(extra_column_names)free(extra_column_names);

	return retvalue;

}


