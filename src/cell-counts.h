#include "subread.h"
#include "core.h"
#include "seek-zlib.h"
#ifndef __CELL_COUNTS_H
#define __CELL_COUNTS_H

#define SCRNA_VBUFF_SIZE (32*1024*1024)
#define MAX_FC_READ_LENGTH 10001
#define READ_BIN_BUF_SIZE 1000 // sufficient for a <=150bp read.
#define CELLBC_BATCH_NUMBER 149
#define MAX_UMI_LEN 14 // cannot be higher than 16: must be able to encode into a 32-bit integer.
#define MAX_CELLBC_LEN 24
#define MAX_SCRNA_SAMPLE_NUMBER 64 
#define MAX_SUBREADS_PER_READ 32
#define SCRNA_SUBREADS_HARD_LIMIT 20 
#define REVERSED_READ_BIN_OFFSET ( MAX_SCRNA_READ_LENGTH /4+1 )
#define MIN_LEN_VISIUM_HD_CELLBC 14


#define JUNCTION_REALIGNMENT_MAX_TRIES 2000
#define JUNCTION_REALIGNMENT_MAX_DEPTH 5
#define JUNCTION_MAX_COLOCATION 10000
#define JUNCTION_MAX_MISMATCHING_BASES_IN_REALIGNMENT 1
#define JUNCTION_MAX_CHRO_DISTANCE 500000
#define JUNCTION_WIDDEN_GAP_LEN 4
#define JUNCTION_MAX_MISMA_MEET max(1,(gaplen * .1) )
#define JUNCTION_MINIMUM_MAIN_HALF ( (int)(all_subreads / 6.0 - 0.0001) )

// the configurations used in our cellCounts paper
#define GENE_SCRNA_VOTE_SPACE 3

#define cellCounts_lock_occupy pthread_mutex_lock
#define cellCounts_lock_release pthread_mutex_unlock

#define GENE_SCRNA_VOTE_TABLE_SIZE 17 

#define NOT__DEBUG_NO_LOOK

#define chroEvent_t_TYPE_INDEL 1
#define chroEvent_t_TYPE_JUNCTION 2
#define chroEvent_t_TYPE_EXON 3
typedef struct{
	int event_type;	// 1==indel; 2==junction
	unsigned int left_edge;	// linear chromosome location for the left edge (included in the read).
	int n_events;	// two edges can define multiple insertions.
	int * length;	// Positive: deletion or junction; negative: insertion.
			// For saving space: if no_events == 1: length = length - NULL; 
			// Othervise, treat length as array of lengths: [0 ~ no_events];
			// The length is the number in Cigar, e.g., 10M2000N10M has length=2000.
	char * inserted_bases;
	int step1_supported_reads;
	int * step2_supported_reads;
	int * step2_non_supported_reads;
	int from_truth;
} chroEvent_t;


// These 3 structure types are copied from featureCounts junction reporting.
typedef struct{
	ArrayList * exons_in_transcript; // the items in IVT tree obj is deallocated when this ArrayList is destroyed.
	char * transcript_id; // this piece of memory is owned by this struct.
	char * gene_name;  // this piece of memory is owned by this struct.
} cct_junction_transcript_t;

typedef struct{
	char * transcript_id; // This should be the memory address in the fc_junction_transcript_t table. Not a dedicated memory space for exons.
	int chro_start;
	int chro_stop;
	int is_negative;
} cct_junction_exon_in_transcript_t;

typedef struct{
	char chromosome_name [MAX_CHROMOSOME_NAME_LEN+1];
	char gene_name [FEATURE_NAME_LENGTH+1];
	int is_negative;
	ArrayList * transcript_list; // values of cct_junction_transcript_t
} cct_junction_genebody_t;

typedef struct{
	gene_vote_number_t max_vote;
	int max_vote_IJ;
	gehash_data_t max_position;
	gene_quality_score_t max_quality;
	gene_vote_number_t max_indel_recorder[MAX_INDEL_TOLERANCE*3];
	gene_vote_number_t * max_tmp_indel_recorder;
	int max_mask;
	gene_vote_number_t noninformative_subreads;

	unsigned short items[GENE_SCRNA_VOTE_TABLE_SIZE];
	unsigned int pos [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	int masks [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	int marked_shift_indel[GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	gene_vote_number_t votes [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	gene_quality_score_t quality [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	gene_vote_number_t last_subread_cluster [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	gene_vote_number_t indel_recorder [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE][MAX_INDEL_TOLERANCE*3];
	char current_indel_cursor[GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	char toli[GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	int topK_votes[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	int topK_IJ[SCRNA_HIGHEST_REPORTED_ALIGNMENTS ];

	short coverage_start [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	short coverage_end [GENE_SCRNA_VOTE_TABLE_SIZE][GENE_SCRNA_VOTE_SPACE];
	short max_coverage_start;
	short max_coverage_end;
} gene_sc_vote_t;

typedef struct {
	unsigned int selected_position;
	short result_flags;
	short read_length;
	// 4 bytes
	gene_vote_number_t selected_votes;
	gene_vote_number_t used_subreads_in_vote;
	char indels_in_confident_coverage;
	char is_fully_covered;
	gene_vote_number_t selected_indel_record [MAX_INDEL_SECTIONS*3 + 1];
	unsigned short confident_coverage_start;
	unsigned short confident_coverage_end;
} voting_location_t;

typedef struct {
	int votes[MAX_SUBREADS_PER_READ *2];
	int offsets[MAX_SUBREADS_PER_READ *2];
	unsigned int * start_location_in_index [MAX_SUBREADS_PER_READ *2]; // * 2 because positive and negative strand votes are evaluated together.
} temp_votes_per_read_t;

typedef struct {
	unsigned int linear_l, linear_r;
	short insertion_idx;
	short matching_bases_in_alignment, mismatching_bases_in_alignment;
	short read_covered_first_base, read_covered_last_base;
	chroEvent_t *event_details;
} realignment_event_stack_item_t;

typedef struct{
	FILE * fp;
	unsigned char rle_buffer[31];
	unsigned char rle_buffer_used;
	unsigned char rle_run_byte;
	unsigned char rle_run_repeats;
	unsigned char rle_run_active;
} cellcounts_temp_file_point_t;

typedef struct{
	int thread_no;
	pthread_t thread;

	int event_space_capacity;
	int total_events;
	HashTable * event_entry_table;

	topK_buffer_t topKbuff;
	short * final_reads_mismatches_array;
	short * final_counted_reads_array;

	int hits_number_capacity;
	int * hits_start_pos;
	int * hits_length;
	char ** hits_chro;
	srInt_64 * hits_indices;

	srInt_64 mapped_reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER];
	srInt_64 assigned_reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER];
	srInt_64 reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER+1];
	srInt_64 bcl_input_local_start_no;
	int bcl_input_local_filled, bcl_input_local_cached;
	char bcl_input_local_readbin[BCL_READBIN_ITEMS_LOCAL][BCL_READBIN_SIZE];
	int bcl_input_local_readlane[BCL_READBIN_ITEMS_LOCAL];

	srInt_64 hiconf_map;
	srInt_64 loconf_map;
	int populating_voteIJ_buf_index;
	int total_voteIJs_to_write, writing_voteID_buf_index;
	srInt_64 reporting_scores[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	srInt_64 reporting_flags[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	unsigned int reporting_positions[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	int reporting_mapq[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	char reporting_cigars[SCRNA_HIGHEST_REPORTED_ALIGNMENTS][MAX_SCRNA_READ_LENGTH+20];
	int reporting_editing_distance[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	int reporting_vote_for_aln[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	srInt_64 reporting_ma_misma_ins_Sclip[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];

	int temp_realign_record_capacity;
	int temp_realign_work_capacity;
	unsigned char tempbin_v_buffer[SCRNA_VBUFF_SIZE];
	unsigned char * temp_realign_record_buf;
	unsigned char * temp_realign_work_buf;

	char ** dynamic_align_buffers[4];
	int dynamic_align_penalties[4];

	// realignment stack related data
	cellcounts_temp_file_point_t realign_temp_fp;
	int realignment_event_stack_runcount;
	int realignment_event_stack_current_depth;
	int realignment_event_stack_best_score;
	realignment_event_stack_item_t realignment_event_current_stack[JUNCTION_REALIGNMENT_MAX_DEPTH];
	char * realignment_event_read_name;
	ArrayList * best_3end_stack_list;
	ArrayList * best_5end_stack_list;
	HashTable * alignment_repating_table;
	HashTable * junction_to_cell_umi_table[ MAX_SCRNA_SAMPLE_NUMBER+1 ];

	// realignment candidature related data
	int realignment_alignment_candidature_score;
} cellcounts_align_thread_t;

typedef struct{
	HashTable       * chroEvent_entry_table;
	HashTable	* chroEvent_detail_table;
} junction_index_t;

typedef struct{
	int total_threads;
	cellcounts_align_thread_t * all_thread_contexts;
	int reads_per_chunk;
	int allow_multi_overlapping_reads;
	int max_candidate_voteIJ_per_read;
	int max_reported_alignments_per_read;
	int max_indel_length;
	int max_distinct_top_vote_numbers;
	int max_differential_from_top_vote_number;
	int max_mismatching_bases_in_reads;
	int min_mapped_length_for_mapped_read;
	int min_votes_per_mapped_read;
	int enable_soft_clipping;
	int total_subreads_per_read;
	int report_multi_mapping_reads;
	int is_BAM_and_FQ_out_generated;
	int current_dataset_no;

	int processed_reads_in_chunk;
	int running_processed_reads_in_chunk;

	srInt_64 mapped_reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER];
	srInt_64 assigned_reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER];
	srInt_64 reads_per_sample[MAX_SCRNA_SAMPLE_NUMBER];

	srInt_64 hiconf_map;
	srInt_64 loconf_map;
	srInt_64 all_processed_reads_before_chunk;
	double program_start_time;
	int is_final_voting_run;
	int longest_chro_name;
	int output_binfiles_are_full;
	gene_inputfile_position_t current_circle_start_position, current_circle_end_position;
	int last_written_fragment_number;

	char index_prefix[MAX_FILE_NAME_LENGTH];
	char output_prefix[MAX_FILE_NAME_LENGTH];
	char temp_file_dir[MAX_FILE_NAME_LENGTH];
        char read_assignment_detail_file[MAX_FILE_NAME_LENGTH];
	char input_dataset_name[MAX_FILE_NAME_LENGTH * MAX_SCRNA_FASTQ_FILES * 3];
	int input_mode;

	int total_index_blocks;
	int current_index_block_number;
	gene_value_index_t * value_index;
	gene_input_t input_dataset;
	gehash_t * current_index;
	cellCounts_lock_t input_dataset_lock;

	char cell_barcode_list_file[MAX_FILE_NAME_LENGTH];
	char bcl_sample_sheet_file[MAX_FILE_NAME_LENGTH];
	int visium_hd_barcodes;
	int known_cell_barcode_length;
	int is_dual_index;
	HashTable * cell_barcode_head_tail_table;
	ArrayList * cell_barcodes_array;
	HashTable * sample_sheet_table;
	ArrayList * sample_barcode_list;
	ArrayList * sample_id_to_name;
	HashTable * lineno1B_to_sampleno1B_tab;
	FILE * batch_files[CELLBC_BATCH_NUMBER+2];
	cellCounts_lock_t batch_file_locks[CELLBC_BATCH_NUMBER+2];
	HashTable * sample_BAM_writers;

	parallel_gzip_writer_t fastq_unassigned_writer[4];
	cellCounts_lock_t fastq_unassigned_lock;
	pthread_t thread_delete_files;
	cellCounts_lock_t read_assignment_detail_lock;
	FILE * read_assignment_detail_fp;

	int UMI_length;
	int barcode_batched_max_genes;
	int barcode_batched_max_Rbin_len;
	float umi_cutoff;
	int applied_umi_cut[MAX_SCRNA_SAMPLE_NUMBER];
	int do_one_batch_runner_current;
	int report_excluded_barcodes;
	int has_error;
	int need_check_strand;
	
	char features_annotation_file[MAX_FILE_NAME_LENGTH];
	char features_annotation_alias_file[MAX_FILE_NAME_LENGTH];
	int  features_annotation_file_type;
	char features_annotation_gene_id_column[MAX_READ_NAME_LEN];
	char features_annotation_feature_type[MAX_READ_NAME_LEN];
	srInt_64 * block_min_start, *block_max_end, *block_end_index;
	gene_offset_t chromosome_table;
	ArrayList * all_features_array;

	HashTable * chromosome_exons_table;
	unsigned char ** gene_name_array;
	HashTable * gene_name_table; 
	char * unistr_buffer_space;
	srInt_64 unistr_buffer_size;
	srInt_64 unistr_buffer_used;
	char * cmd_rebuilt;


	unsigned char 	* features_sorted_strand;
	srInt_64 	* features_sorted_start, * features_sorted_stop;
	int 		* features_sorted_geneid;
	char 		** features_sorted_chr;
	HashTable 	* sam_chro_to_anno_chr_alias;

	int do_cell_level_junction_detection;		// switch to enable step1 (pre-alignment)
	cellCounts_lock_t * read_assignment_counter_locks; 
	int 		    chroEvent_lock_number;

	cellCounts_lock_t chroEvent_entry_table_lock;
	HashTable	* junction_to_cell_umi_table[ MAX_SCRNA_SAMPLE_NUMBER+1 ];
	HashTable       * chroEvent_entry_table[MAX_SCRNA_SAMPLE_NUMBER+1];
	HashTable	* chroEvent_detail_table[MAX_SCRNA_SAMPLE_NUMBER+1];
	HashTable       * transcript_exon_table;
	HashTable       * transcript_to_gene_name_table;
	char		* exonic_region_bitmap;

	HashTable	* cluster_spec_junction_table; // sample_no << 56 | cluster_no => table : (junction linear left <<32|right) => 1 // i.e. this junction is for this cluster spec
	HashTable	* cluster_cell_map_table; // sample_no << 56 | cell_bc_no => cluster_no ; cell_bc_no and cluster_no are 0-based; sample_no is 1-based 
	char		cluster_junctions_file[MAX_FILE_NAME_LENGTH];
	char		cluster_map_file[MAX_FILE_NAME_LENGTH];

	// The following 4 variables were copied from featureCounts junction reporting.
	HashTable 	* junction_ExonEdgeTree_table[3];
	HashTable	* junction_GenebodyTree_table;
	HashTable	* junction_transcript_table;
	HashTable	* junction_genebody_table;
	int 		ignore_transcript_junction_assignment;
} cellcounts_global_t;

typedef struct{
	int thread_no;
} cellcounts_final_thread_t;

#define CHROMOSOME_NAME_LENGTH 256 

struct TempForRealign{
	int saved_alignments;
	unsigned short num_of_votes[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	unsigned int   voted_position[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	unsigned short coverage_start[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	unsigned short coverage_end[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];
	unsigned char  flags[SCRNA_HIGHEST_REPORTED_ALIGNMENTS];

	unsigned int sample_number;
	unsigned int cell_number;
	unsigned int read_number;
	unsigned int raw_umi_sequence;
	unsigned short sample_seq_length;
	unsigned short sample_qual_length;
	int read_length;
	unsigned char qual_values[MAX_SCRNA_READ_LENGTH];
	unsigned char read_bases[MAX_SCRNA_READ_LENGTH/2];
	unsigned char cellbc_umi_qual[MAX_UMI_LEN+MAX_CELLBC_LEN];
	unsigned char cellbc_umi_bases[(MAX_UMI_LEN+MAX_CELLBC_LEN)/2+1];
};

int cellCounts_select_and_write_temps(cellcounts_global_t * cct_context, int thread_no, int sample_i, gene_sc_vote_t * votetab, char * read_name, char * read_text, char * read_bin, char * read_qual, int read_len, gene_vote_number_t all_subreads);

#endif
