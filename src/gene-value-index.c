#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "gene-value-index.h"
#include "input-files.h"


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
	int read_length;
	read_length = fread(&index->start_point,4,1, fp);
	assert(read_length>0);
	read_length = fread(&index->length,4,1, fp);
	assert(read_length>0);

	//printf ("\nBINDEX %s : %u ~ +%u\n",filename, index->start_point, index->length );

	int useful_bytes, useful_bits;
	index -> start_base_offset = index -> start_point - index -> start_point%4;
	gvindex_baseno2offset (index -> length+ index -> start_point, index ,&useful_bytes,&useful_bits);
	index -> values = malloc(useful_bytes);
	index -> values_bytes = useful_bytes;

	read_length =fread(index->values, 1, useful_bytes, fp);
	assert(read_length>0);

	fclose(fp);

}


int match_chro_wronglen(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int space_type, int * left_match_bases, int * right_match_bases)
{
	int ret = 0;
	int i;
	char last_char='A';
	int left_correct_end = 0;
	if(left_match_bases) *left_match_bases=0;
	if(right_match_bases) *right_match_bases=0;

	if (space_type == GENE_SPACE_COLOR)
		last_char = (pos <= index -> start_point)?'A': gvindex_get(index,pos-1);

	for (i=0;i<test_len;i++)
	{
		char tt = gvindex_get (index, pos +i);
		int newv;
		if(space_type == GENE_SPACE_COLOR)
		{

			newv = read[i] == '0'+chars2color(last_char, tt); 
			last_char = tt;
		}
		else
			newv =read[i] == tt; 

		//if(left_wrong_bases)
		//	printf("I=%d, *LWB=%d, LWE=%d\n", i, *left_wrong_bases, left_wrong_end);

		if(left_match_bases && (newv) && (!left_correct_end ))
			(*left_match_bases)++;
		else if (!newv)left_correct_end=1;

		if(right_match_bases && (newv))
			(*right_match_bases) ++;
		else if (right_match_bases)
			(*right_match_bases) =0;

		ret += newv;
	}

	return ret;
}

#define INDEL_TEST_WINDOW 3

int match_indel_chro_to_front(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int * indels, int * indel_point, int max_indel_number)
{
	int offset = 0;
	int i;
	int ret = 0;

	for(i=0; i < test_len+min(0,offset); i++)
	{
		char tt = gvindex_get (index, pos + i + max(0, offset));

		if(read[i-min(0,offset)]==tt) ret++;
		else if(i + offset < test_len - INDEL_TEST_WINDOW - 3 && i >0)
		{
			// if there is a base unmatched, it is potentially an indel from here.
			float bast_match_score_remailing=-1;
			int window_match = match_chro(read+i-min(0,offset), index, pos+i+ max(0, offset), INDEL_TEST_WINDOW ,0,GENE_SPACE_BASE);

			if(window_match < INDEL_TEST_WINDOW -1)
			{
				// if the window is badly matched, it is very likely to be an indel from this base.
				int indel_test_i;
				for(indel_test_i =0; indel_test_i < 7; indel_test_i++)
				{
					int indel_test = (indel_test_i+1)/2*(indel_test_i%2?1:-1);
					if(abs(indel_test)>max_indel_number) continue;

					if(indel_test > 0)	// deletion
					{
						int matched_tail = match_chro(read+i, index, pos+i+indel_test, test_len - i,0,GENE_SPACE_BASE);
						float matched_score = matched_tail * 1. / ( test_len - i);
						//printf("INDEL_DEL_TEST i=%d: Indel=%d, Score=%f\n",i, indel_test,matched_score  );
						if(matched_score >  bast_match_score_remailing &&  matched_score > 0.8)
						{
							offset = indel_test;
							bast_match_score_remailing = matched_score;
						}
					}else	// insertion
					{
						int matched_tail = match_chro(read+i - indel_test, index, pos+i, test_len - i + offset ,0,GENE_SPACE_BASE);
						float matched_score = matched_tail * 1. / (test_len - i + offset);
						//printf("INDEL_INS_TEST i=%d: Indel=%d, Score=%f\n",i, indel_test,matched_score  );
						if(matched_score >  bast_match_score_remailing &&  matched_score > 0.8)
						{
							offset = indel_test;
							bast_match_score_remailing = matched_score;
						}
					}
				}
			}
			if(bast_match_score_remailing>0)
			{
				if(offset > 0)//deletion
				{
					tt = gvindex_get (index, pos + i + offset);
					ret += read[i] == tt;
				}
				else
					ret += read[i - offset] == tt;
				*indel_point  = i;
				
			}
		}
	}
	*indels = offset;
	return ret;

}


// "pos" here is the expected position of the head of the read , but it is not the final position. If indels are found in the read, the head position must be offset.
// Only certain is pos+test_len is the EXACT position of the TAIL of the read.
int match_indel_chro_to_back(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int * indels, int * indel_point, int max_indel_number)
{
	//return  match_chro(read, index, pos, test_len, 0, 1);
	//printf("TEST_INDEL_CHRO %s VS %u LEN=%d\n", read, pos, test_len);
	int offset = 0;
	int i;
	int ret = 0;

	for(i=test_len-1 ; i >=max(0,offset); i--)
	{
		char tt = gvindex_get (index, pos + i - max(0, offset));
		#ifdef TEST_TARGET
		if(memcmp(read, TEST_TARGET, 15)==0)
			printf("%c=?=%c OFF=%d\n",read[i+min(0,offset)], tt, offset);
		#endif

		if(read[i+min(0,offset)]==tt) ret++;
		else if(i + offset >INDEL_TEST_WINDOW +3 && i < test_len-1)
		{
			// if there is a base unmatched, it is potentially an indel from here.
			float bast_match_score_remailing=-1;
			int window_match = match_chro(read+i-min(0,offset)-INDEL_TEST_WINDOW, index, pos+i+ max(0, offset)-INDEL_TEST_WINDOW, INDEL_TEST_WINDOW ,0,GENE_SPACE_BASE);

			if(window_match < INDEL_TEST_WINDOW-1)
			{
				// if the window is badly matched, it is very likely to be an indel from this base.
				int indel_test_i;
				for(indel_test_i =0; indel_test_i < 7; indel_test_i++)
				{
					int indel_test = (indel_test_i+1)/2*(indel_test_i%2?1:-1);
					if(abs(indel_test)>max_indel_number) continue;
					if(indel_test > 0)	//deletion 
					{
						int matched_tail = match_chro(read , index, pos-indel_test, i ,0,GENE_SPACE_BASE);
						float matched_score = matched_tail * 1. / ( i );

						#ifdef TEST_TARGET
						if(memcmp(read, TEST_TARGET, 15)==0)
							printf("INDEL_DEL_TEST i=%d: Indel=%d, Score=%f, HEADPOS=%u\n",i, indel_test,matched_score , pos-indel_test );
						#endif

						if(matched_score >  bast_match_score_remailing &&  matched_score > 0.8)
						{
							offset = indel_test;
							bast_match_score_remailing = matched_score;
						}
					}else	//insertion 
					{
						int matched_tail = match_chro(read, index, pos - indel_test, i ,0,GENE_SPACE_BASE);
						float matched_score = matched_tail * 1. / (i + indel_test);
						#ifdef TEST_TARGET
						if(memcmp(read, TEST_TARGET, 15)==0)
							printf("INDEL_INS_TEST i=%d: Indel=%d, Score=%f, HEADPOS=%u\n",i, indel_test,matched_score  , pos + indel_test);
						#endif
						if(matched_score >  bast_match_score_remailing &&  matched_score > 0.8)
						{
							offset = indel_test;
							bast_match_score_remailing = matched_score;
						}
					}
				}
			}
			if(bast_match_score_remailing>0)
			{
				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
					printf("PIECE_LEN=%d ; INDEL AT %d ; INDEL=%d\n", test_len , i + min(0, offset) , offset);
				#endif
				if(offset > 0)//insertion
					ret += read[i-offset] == tt;
				else	//deletion
				{
					tt = gvindex_get (index, pos + i - offset);
					ret += read[i] == tt;
				}
				*indel_point  = i + min(0, offset);
			}
		}
	}
	//if(memcmp(read, "AACCCCTTGCAGAAAA", 15)==0)
	//	printf("\nINDEL_SPOT %d\n", offset);
	*indels = offset;
	return ret;
}

int match_chro(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type)
{
	int ret = 0;
	int i;
	char last_char='A';
	

	if (is_negative_strand)
	{

		if (space_type == GENE_SPACE_COLOR)
		{
			pos++;
			last_char = (pos+test_len>= index -> length + index -> start_point)?'A': gvindex_get(index,pos+test_len);
		}
		for (i=test_len -1;i>=0;i--)
		{
			char tt = gvindex_get (index, pos+test_len-1-i);
			if(space_type == GENE_SPACE_COLOR)
			{
				ret += read[i] == '0'+chars2color(tt, last_char); 
				last_char = tt;
			}
			else
				switch(tt)
				{
					case 'A': ret += read[i] == 'T'; break;
					case 'T': ret += read[i] == 'A'; break;
					case 'G': ret += read[i] == 'C'; break;
					case 'C': ret += read[i] == 'G'; break;
				}
		}
	}
	else
	{
		if (space_type == GENE_SPACE_COLOR)
			last_char = (pos <= index -> start_point)?'A': gvindex_get(index,pos-1);
		for (i=0;i<test_len;i++)
		{
			char tt = gvindex_get (index, pos +i);
			if(space_type == GENE_SPACE_COLOR)
			{
				ret += read[i] == '0'+chars2color(last_char, tt);
				last_char = tt;
			}
			else
				ret +=read[i] == tt; 
		}
	}
	return ret;
}


int match_chro_maxerror(char * read, gene_value_index_t * index, unsigned int pos, int test_len, int is_negative_strand, int space_type, int max_error)
{
	int ret = 0;
	int i;
	char last_char='A';
	

	if (is_negative_strand)
	{

		if (space_type == GENE_SPACE_COLOR)
		{
			pos++;
			last_char = (pos+test_len>= index -> length + index -> start_point)?'A': gvindex_get(index,pos+test_len);
		}
		for (i=test_len -1;i>=0;i--)
		{
			char tt = gvindex_get (index, pos+test_len-1-i);
			if(space_type == GENE_SPACE_COLOR)
			{
				ret += read[i] != '0'+chars2color(tt, last_char); 
				last_char = tt;
			}
			else
				switch(tt)
				{
					case 'A': ret += read[i] != 'T'; break;
					case 'T': ret += read[i] != 'A'; break;
					case 'G': ret += read[i] != 'C'; break;
					case 'C': ret += read[i] != 'G'; break;
				}
			if(ret>max_error)return 0;
		}
	}
	else
	{
		if (space_type == GENE_SPACE_COLOR)
			last_char = (pos <= index -> start_point)?'A': gvindex_get(index,pos-1);
		for (i=0;i<test_len;i++)
		{
			char tt = gvindex_get (index, pos +i);
			if(space_type == GENE_SPACE_COLOR)
			{
				ret += read[i] != '0'+chars2color(last_char, tt);
				last_char = tt;
			}
			else
				ret +=read[i] != tt; 
			//printf("RET=%d\n",ret);
			if(ret>max_error)return 0;
		}
	}
	return test_len - ret;
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


void gvindex_get_string(char *buf, gene_value_index_t * index, unsigned int pos, int len, int is_negative_strand)
{
	int i;
	if (is_negative_strand)
		for (i=len-1;i>=0;i--)
		{
			buf[i] = gvindex_get (index, pos + len - 1 - i);
			switch(buf[i])
			{
				case 'A': buf[i] = 'T'; break;
				case 'T': buf[i] = 'A'; break;
				case 'G': buf[i] = 'C'; break;
				case 'C': buf[i] = 'G'; break;
			}
		}
	else
		for (i=0;i<len;i++)
			buf[i] = gvindex_get (index, pos +i);
}


