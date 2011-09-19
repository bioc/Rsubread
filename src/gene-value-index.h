#ifndef __GENE_VALUE_INDEX_
#define __GENE_VALUE_INDEX_

#include "subread.h"
#include "sorted-hashtable.h"

void gvindex_init(gene_value_index_t * index, unsigned int start_point, unsigned int base_number);

void gvindex_set (gene_value_index_t * index, gehash_data_t offset, gehash_key_t base_value);

void gvindex_dump(gene_value_index_t * index, const char filename []);

void gvindex_load(gene_value_index_t * index, const char filename []);

void gvindex_destory(gene_value_index_t * index);

void gvindex_baseno2offset(unsigned int base_number, gene_value_index_t * index, int * offset_byte, int * offset_bit);

// returns a 16-bit bitmap showing if each base is matched.
int gvindex_match(gene_value_index_t * index, gehash_data_t offset, gehash_key_t base_value);

int gvindex_match_base(gene_value_index_t * index, gehash_data_t offset, const char base_int_value);

int gvindex_get(gene_value_index_t * index, gehash_data_t offset);

void gvindex_get_string(char *buf, gene_value_index_t * index, unsigned int pos, int len, int is_negative_strand);
#endif
