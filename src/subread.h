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

typedef unsigned int gehash_key_t;
typedef unsigned int gehash_data_t;
typedef double gene_vote_number_t;

#define GENE_TABLE_SIZE 500000000

#define GENE_SLIDING_STEP 3
#define BEXT_RESULT_LIMIT 16

#define IS_DELETION 1
#define IS_INSERTION 2
#define IS_PAIRED_MATCH 128

#define GENE_VOTE_SPACE 64
#define GENE_VOTE_TABLE_SIZE 91

#define base2int(c) ((c)=='G'?1:((c)=='A'?0:((c)=='C'?2:3)))
#define int2base(c) ((c)==1?'G':((c)==0?'A':((c)==2?'C':'T')))
#define color2int(c) ((c) - '0')

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
	char max_mask;

        unsigned char items[GENE_VOTE_TABLE_SIZE];
        unsigned int pos [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
        gene_vote_number_t votes [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	char masks [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short last_offset [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
//	unsigned char last_offset [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
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
	unsigned char * masks;
#ifdef REPORT_ALL_THE_BEST
	gene_best_record_t * best_records;
#endif

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


double miltime();


#define abs(a) 	  ((a)>=0?(a):-(a))
#define max(a,b)  ((a)<(b)?(b):(a))
#define min(a,b)  ((a)>(b)?(b):(a))


#endif
