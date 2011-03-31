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
 	#define init_gene_vote(a) bzero((a)->items, GENE_VOTE_TABLE_SIZE*2); (a)->max_vote = 0;
#else
	#define init_gene_vote(a) bzero((a)->items, GENE_VOTE_TABLE_SIZE);  (a)->max_vote = 0;
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
void init_allvote(gene_allvote_t* allvote, int expected_len);
void clear_allvote(gene_allvote_t* allvote);
void add_allvote(gene_allvote_t* allvote,int qid, int pos, gene_vote_number_t votes, int is_counterpart, char mask);
void add_allvote_q(gene_allvote_t* allvote,int qid, int pos, gene_vote_number_t votes, int is_counterpart, char mask);

unsigned char get_next_char(FILE * fp);

unsigned char * replica_index;

extern double begin_ftime;

void print_text_scrolling_bar(char * hint, float percentage, int width, int * internal_counter);

void print_running_log(double finished_rate, double read_per_second, double expected_seconds, unsigned long long int total_reads) ;


int select_positions(gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2, gehash_data_t * pos_read1, gehash_data_t * pos_read2, unsigned int max_pair_dist, unsigned int min_pair_dist, int min_major, int min_minor, int is_negatve);


int select_positions_array(char * read1_str, int read1_len, char * read2_str, int read2_len,  gene_vote_t * vote_read1, gene_vote_t * vote_read2, gene_vote_number_t * numvote_read1, gene_vote_number_t * numvote_read2, gehash_data_t * pos_read1, gehash_data_t * pos_read2, unsigned int max_pair_dest, unsigned int min_pair_dest, int min_major, int min_minor,int is_negative_strand, gene_value_index_t * my_array_index, int color_space, int indel_tolerance);


#define QUALITY_SCALE_NONE 0
#define QUALITY_SCALE_LINEAR 1
#define QUALITY_SCALE_LOG 2


// This function returns the lowest score of phred in the given 16-byte quality string
// Assuming that the score = char - '@'
gene_vote_number_t get_subread_quality(const char * quality_str, const char * read_str, int quality_scale);

gene_vote_number_t calculate_penalty_score(const char * quality_str, const char * read_str, int quality_scale, int base_matchness);

int is_valid_subread(const char * read_str);

void final_matchingness_scoring(const char read_str[], const char quality_str[], int read_len,  gene_vote_t * vote, gehash_data_t * max_position, gene_vote_number_t * max_vote, char *max_mask , gene_value_index_t * my_array_index, int color_space, int indel_tolerance);

float match_read(const char read_str[], int read_len, unsigned int potential_position,  gene_value_index_t * my_array_index, int space_type, int indel_tolerance) ;

#endif
