#include <stdio.h>
#include<assert.h>
#include <stdlib.h>
#include "gene-value-index.h"


void gvindex_init(gene_value_index_t * index, unsigned int start_point, unsigned int base_number)
{
	index->start_point = start_point;
	index->length = base_number;
	index->values = malloc(base_number / 4 + 1);
	index -> start_base_offset = index -> start_point - index -> start_point%4;
}


void gvindex_baseno2offset(unsigned int base_number, gene_value_index_t * index, int * offset_byte, int * offset_bit)
{
	// the base number corrsponding to the 0-th bit in the whole value array;

	unsigned int offset = (base_number - index -> start_base_offset);

	* offset_byte = offset / 4;
	* offset_bit = base_number % 4 * 2;
}

// return 'A', 'G', 'T' and 'C'
int gvindex_get(gene_value_index_t * index, gehash_data_t offset)
{
	int offset_byte, offset_bit;
	gvindex_baseno2offset(offset, index , &offset_byte, &offset_bit);

	unsigned char mask = 0x3 << (offset_bit);
	unsigned int one_base_value = (index->values [offset_byte] & mask) >> (offset_bit);

	return int2base(one_base_value);
}

int gvindex_match(gene_value_index_t * index, gehash_data_t offset, gehash_key_t base_values)
{
	int offset_byte, offset_bit;

	gvindex_baseno2offset(offset, index , &offset_byte, &offset_bit);
	int i, ret = 0;

	for (i=0; i<16; i++)
	{
		unsigned char mask = 0x3 << (offset_bit);
		unsigned char one_base_value = (index->values [offset_byte] & mask) >> (8-offset_bit);
		if ( ((base_values >> (30 - i*2)) & 0x3) == one_base_value)
			ret |= 1 << i;

		offset_bit +=2;
		if(offset_bit >=8)
		{
			offset_bit = 0;
			offset_byte ++;
		}
	}

	return ret;

}

void gvindex_set (gene_value_index_t * index, gehash_data_t offset, gehash_key_t base_values)
{
	int offset_byte, offset_bit;
	gvindex_baseno2offset(offset, index , &offset_byte, &offset_bit);
	int i;

	for (i=0; i<16; i++)
	{
		// 11110011
		//     ^^ base
		unsigned char mask = 0xff << (offset_bit+2) | 0xff >> (8-offset_bit);
		index->values [offset_byte] &= mask;
		index->values [offset_byte] |= ((base_values >> (30 - i*2))&0x03) << (offset_bit);

		offset_bit +=2;
		if(offset_bit >=8)
		{
			offset_bit = 0;
			offset_byte ++;
		}
	}

	index -> length = offset + 16 - index -> start_point ;
}

void gvindex_dump(gene_value_index_t * index, const char filename [])
{
	FILE * fp = fopen(filename, "wb");

	fwrite(&index->start_point,4,1, fp);
	fwrite(&index->length, 4, 1, fp);

	int useful_bytes, useful_bits;
	gvindex_baseno2offset (index -> length+ index -> start_point, index,&useful_bytes,&useful_bits);

	fwrite(index->values, 1, useful_bytes, fp);

	fclose(fp);
}


void gvindex_load(gene_value_index_t * index, const char filename [])
{
	FILE * fp = fopen(filename, "rb");
	assert(0<fread(&index->start_point,4,1, fp));
	assert(0<fread(&index->length,4,1, fp));

	//printf ("\nBINDEX %s : %u ~ +%u\n",filename, index->start_point, index->length );

	int useful_bytes, useful_bits;
	index -> start_base_offset = index -> start_point - index -> start_point%4;
	gvindex_baseno2offset (index -> length+ index -> start_point, index ,&useful_bytes,&useful_bits);
	index -> values = malloc(useful_bytes);
	index -> values_bytes = useful_bytes;

	assert(0<fread(index->values, 1, useful_bytes, fp));

	fclose(fp);

}

int gvindex_match_base(gene_value_index_t * index, gehash_data_t offset, const char base_int_value)
{
	int offset_byte, offset_bit;

	gvindex_baseno2offset(offset, index, &offset_byte, &offset_bit);

	unsigned char mask = 0x3 << (offset_bit);

	if(offset_byte >= index->values_bytes)
		return 0;
//		printf("\nERROR: %u > %u\n", offset_byte, index->values_bytes);

	char reference_base = ((index->values [offset_byte] & mask) >> offset_bit);

	return  (reference_base == base_int_value)?1:0 ;
}

void gvindex_destory(gene_value_index_t * index)
{
	free(index -> values);
}


