#ifndef _GENE_ALGORITHMS_H_
#define _GENE_ALGORITHMS_H_

#include "subread.h"
#include "sorted-hashtable.h"

// Load the ".reads" file into the memory
// Return 0 if succeed or -1 if errors.
int load_offsets(gene_offset_t* offsets , const char index_prefix []);

// Locate the position of a linear address
// Return 0 if the linear position is in a reasonable range or -1 if it is out of range.
// The pointer to the name of the chromosome is put into chro_name, and the position in this chromosome is in pos.
int locate_gene_position(unsigned int linear, const gene_offset_t* offsets , char ** chro_name, unsigned int * pos);


int remove_repeated_reads(gehash_t * table, gehash_t * huge_table,int index_threshold);

#ifdef __MAX_ACCURACY_
 	#define init_gene_vote(a) bzero((a)->items, GENE_VOTE_TABLE_SIZE*2); (a)->max_vote = 0; (a) -> max_indel_recorder[0]=0; (a)->max_tmp_indel_recorder = NULL; (a)->max_mask = 0;
#else
	#define init_gene_vote(a) bzero((a)->items, GENE_VOTE_TABLE_SIZE);  (a)->max_vote = 0; (a) -> max_indel_recorder[0]=0;  (a)->max_tmp_indel_recorder = NULL;  (a)->max_mask = 0;
#endif
// return current votes for a given position
// if create_new_pos == 0 then do not take this position if it does not exist in the vote array
inline void add_gene_vote(gene_vote_t* vote, int position, int create_new_pos);
inline void add_gene_vote_weighted(gene_vote_t* vote, int position, int create_new_pos, int w);

// return the votes; the position is put into position_result 
int max_gene_vote(gene_vote_t* vote, int * position_result, int query_id);

// evaluate the number of matched characters in the piece_str
int evaluate_piece(char * piece_str, int chromosome, int offset, int is_counterpart, int start_p, int end_p);

// These two functions maintain the global vote for pieces
void init_allvote(gene_allvote_t* allvote, int expected_len, int allowed_indels);
void clear_allvote(gene_allvote_t* allvote);


void add_allvote_q(gene_allvote_t* allvote,int qid , int pos, gene_vote_number_t votes, gene_quality_score_t quality, int is_counterpart, short mask, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads, int space_type, int report_junction, int is_head_high_quality, char * qual_txt, int phred_version, char span_coverage);

unsigned char get_next_char(FILE * fp);

unsigned char * replica_index;

extern double begin_ftime;

void print_text_scrolling_bar(char * hint, float percentage, int width, int * internal_counter);

void print_running_log(double finished_rate, double read_per_second, double expected_seconds, unsigned long long int total_reads, int is_pair) ;


int select_positions(gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2, gene_quality_score_t * sum_quality, gene_quality_score_t * qual_r1,  gene_quality_score_t * qual_r2, gehash_data_t * pos_read1, gehash_data_t * pos_read2, char * read1_indel_recorder, char* read2_indel_recorder, unsigned int max_pair_dist, unsigned int min_pair_dist, int min_major, int min_minor, int is_negatve, int number_of_anchors_quality, int max_indel_len, int * is_breakeven);


int select_positions_array(char * read1_str, int read1_len, char * read2_str, int read2_len,  gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2, gene_quality_score_t * sum_quality, gene_quality_score_t * q_r1, gene_quality_score_t *q_r2, char * indel_rec_r1, char * indel_rec_r2, gehash_data_t * pos_read1, gehash_data_t * pos_read2, unsigned int max_pair_dest, unsigned int min_pair_dest, int min_major, int min_minor,int is_negative_strand, gene_value_index_t * my_array_index, int color_space, int indel_tolerance, int number_of_anchors_quality, const char quality_str1 [], const char quality_str2 [], int quality_scale, int max_indel_len, int * is_breakeven);


#define QUALITY_SCALE_NONE 0
#define QUALITY_SCALE_LINEAR 1
#define QUALITY_SCALE_LOG 2


// This function returns the lowest score of phred in the given 16-byte quality string
// Assuming that the score = char - '@'
gene_quality_score_t get_subread_quality(const char * quality_str, const char * read_str, int quality_scale, int phred_version);

gene_vote_number_t calculate_penalty_score(const char * quality_str, const char * read_str, int quality_scale, int base_matchness);

int is_valid_subread(const char * read_str);

void final_matchingness_scoring(const char read_str[], const char quality_str[], int read_len, gene_vote_t * vote, gehash_data_t * max_position, gene_vote_number_t * max_vote, short *max_mask, gene_quality_score_t * max_quality, gene_value_index_t * my_array_index, int space_type, int indel_tolerance, int quality_scale, char * max_indel_recorder, int * max_coverage_start, int * max_coverage_end);

float match_read(const char read_str[], int read_len, unsigned int potential_position,  gene_value_index_t * my_array_index, int space_type, int indel_tolerance, const char quality_str [], int quality_scale) ;

void show_cigar(char * info, int len, int is_reversed_map, char * buf, int max_indel, int total_subreads, char * read);

int dynamic_align(char * read, int read_len, gene_value_index_t * index, unsigned int begin_position, int max_indel, char * movement_buffer, int expected_offset,int begin_read_offset, int end_read_offset);

int window_indel_align(char * read, int read_len, gene_value_index_t * index, unsigned int begin_position, int max_indel, char * movement_buffer, int expected_offset,int begin_read_offset, int end_read_offset);

void explain_indel_in_middle(gene_allvote_t* allvote, int qid , int pos, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads, int head_start_point, int tail_end_point, int head_indel_movement, int tail_indel_movement, int find_junctions);

void find_and_explain_indel(gene_allvote_t* allvote,int qid , int pos, gene_vote_number_t votes, gene_quality_score_t quality, int is_counterpart, char mask, char * max_indel_recorder, gene_value_index_t * array_index, char * read_txt, int read_len, int max_indel, int total_subreads, int space_type,int find_junctions, int is_head_high_quality, char * qual_txt, int phred_version);

int extend_covered_region(gene_value_index_t *array_index, unsigned int read_start_pos, char * read, int read_len, int cover_start, int cover_end, int window_size, int match_req_5end, int match_req_3end, int indel_tolerance, int space_type, int tail_indel, short * head_indel_pos, int * head_indel_movement, short * tail_indel_pos, int * tail_indel_movement, int is_head_high_quality, char * qual_str, int qual_format);

float final_mapping_quality(gene_value_index_t *array_index, unsigned int pos, char * read_txt, char * qual_txt, char * cigar_txt, int phred_version);

int bad_quality_base_number(char * qualityb, int rl , int format);

int min_matched_bases(char * qualityb, int rl , int format, float match_score);

float read_quality_score(char * qualityb, int rl , int format);

void compress_cigar(char * cigar,  int read_len, char *read);

void print_votes(gene_vote_t * vote, char *index_prefix);
#endif

