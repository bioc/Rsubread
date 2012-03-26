#ifndef _SUBREAD_H_
#define _SUBREAD_H_

#include <stdlib.h>
#include <stdio.h>

#define SAM_FLAG_PAIRED_TASK	0x01
#define SAM_FLAG_FIRST_READ_IN_PAIR 0x40
#define SAM_FLAG_SECOND_READ_IN_PAIR 0x80
#define SAM_FLAG_MATE_UNMATCHED 0x08
#define SAM_FLAG_MATCHED_IN_PAIR 0x02
#define SAM_FLAG_REVERSE_STRAND_MATCHED 0x10
#define SAM_FLAG_MATE_REVERSE_STRAND_MATCHED 0x20
#define SAM_FLAG_UNMAPPED 0x04

#define MAX_PIECE_JUNCTION_READ 7

#define IS_LONG_OVERLAP 4
#define IS_SHORT_OVERLAP 8
#define IS_PAIRED_HINTED 16
#define IS_R1_CLOSE_TO_5 1
#define IS_REVERSED_HALVES 2
#define	IS_PROCESSED_READ 32
#define	IS_PROCESSED_READ_R2 64
#define IS_PAIRED_MATCH 128
#define IS_NEGATIVE_STRAND_R1 256 
#define IS_NEGATIVE_STRAND_R2 512 
#define IS_FUSION 1024 
#define IS_NEGATIVE_STRAND 2048
#define IS_RECOVERED_JUNCTION_READ 4096
#define IS_FINALISED_PROCESSING 8192
#define IS_RECOVERED_JUNCTION_READ_STEP4 (8192*2)
#define	IS_BREAKEVEN_READ (8192*4)

//#define TEST_TARGET "GCAGGCCGAAGCCGACAAGAA"


typedef unsigned int gehash_key_t;
typedef unsigned int gehash_data_t;
typedef float gene_quality_score_t;
typedef char gene_vote_number_t;

#define GENE_TABLE_SIZE 500000000

#define ANCHORS_NUMBER 259

#define GENE_SLIDING_STEP 3
#define BEXT_RESULT_LIMIT 16

#define SEARCH_BACK 0
#define SEARCH_FRONT 1

#define GENE_VOTE_SPACE 64
#define GENE_VOTE_TABLE_SIZE 91

#define MAX_INDEL_TOLERANCE 16


#define base2int(c) ((c)=='G'?1:((c)=='A'?0:((c)=='C'?2:3)))
//#define int2base(c) ((c)==1?'G':((c)==0?'A':((c)==2?'C':'T')))
#define int2base(c) ("AGCT"[(c)]) 
#define color2int(c) ((c) - '0')
#define int2color(c) ("0123"[(c)])

#define get_base_error_prob64(a) ((a) < '@'-1?1:pow(10., -0.1*((a)-'@')))
#define get_base_error_prob33(a) ((a) < '!'-1?1:pow(10., -0.1*((a)-'!'))) 



#define FASTQ_PHRED33 1
#define FASTQ_PHRED64 0

#define IS_DEBUG 0

typedef struct{
	unsigned int start_base_offset;
	unsigned int start_point;
	unsigned int length;
	unsigned char * values;
	unsigned int values_bytes;
} gene_value_index_t;



struct gehash_bucket {
	int current_items;
	int space_size;
	gehash_key_t * item_keys;
	gehash_data_t * item_values;
};

typedef struct {
	unsigned long long int current_items;
	int buckets_number;
	char is_small_table;
	struct gehash_bucket * buckets;
} gehash_t;


typedef struct {
	gene_vote_number_t max_vote;
	gehash_data_t max_position;
	gene_quality_score_t max_quality;
	char max_indel_recorder[MAX_INDEL_TOLERANCE*3];
	char * max_tmp_indel_recorder;
	short max_mask;

        unsigned short items[GENE_VOTE_TABLE_SIZE];
        unsigned int pos [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
        gene_vote_number_t votes [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
        gene_quality_score_t quality [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short masks [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short last_offset [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	char indel_recorder [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE][MAX_INDEL_TOLERANCE*3];
	char current_indel_cursor[GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];

	#ifdef MAKE_FOR_EXON
	short coverage_start [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short coverage_end [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short max_coverage_start;
	short max_coverage_end;
	//#warning Switch "MAKE_FOR_EXON" is turned on. It may cost more time. Do not turn it on unless you want to detect junction reads.
	#endif
} gene_vote_t ;


typedef struct{
	unsigned char best_len;
	unsigned int offsets [BEXT_RESULT_LIMIT];
	unsigned char is_reverse [BEXT_RESULT_LIMIT];
} gene_best_record_t;



typedef struct{
	int max_len;
	unsigned int * max_positions;
	unsigned char * is_counterpart;
	gene_vote_number_t * max_votes;
	gene_quality_score_t * max_quality;
	gene_quality_score_t * max_final_quality;
	short * masks;
	char * max_indel_recorder;
	char * span_coverage;
#ifdef REPORT_ALL_THE_BEST
	gene_best_record_t * best_records;
#endif
	char max_indel_tolerance;
	short indel_recorder_length;

} gene_allvote_t;


typedef struct{
        char read_name[1000][48];
        unsigned int read_offset[1000];
} gene_offset_t;



typedef struct {
	char filename [300];
	int space_type ;
	int file_type ;
	FILE * input_fp;
} gene_input_t;


typedef struct{
	unsigned int small_key;
	unsigned int big_key;
} paired_exon_key;


double miltime();


#define abs(a) 	  ((a)>=0?(a):-(a))
#define max(a,b)  ((a)<(b)?(b):(a))
#define min(a,b)  ((a)>(b)?(b):(a))


#endif
