#include "sorted-hashtable.h"
#include "gene-algorithms.h"
#include<assert.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include<math.h>

#define _gehash_hash(k) ((unsigned int)(k))

int gehash_create(gehash_t * the_table, size_t expected_size, char is_small_table)
{
	int expected_bucket_number;
	int i;

	if(expected_size ==0)
		expected_size = GEHASH_DEFAULT_SIZE;

	// calculate the number of buckets for creating the data structure
	expected_bucket_number = expected_size / GEHASH_BUCKET_LENGTH; 

	for (;;expected_bucket_number++)
	{
		int j, valid_v;
		
		valid_v = 1;
		for(j=2; j<=13;j++)
		{
			if (expected_bucket_number % j == 0)
				valid_v = 0;
		}

		if (valid_v)
			break;
	}

	the_table -> current_items = 0;
	the_table -> is_small_table = is_small_table;
	the_table -> buckets_number = expected_bucket_number;
	the_table -> buckets = (struct gehash_bucket *)
			malloc(
			  expected_bucket_number *
  			  sizeof(struct gehash_bucket )
			);

	for(i=0; i<expected_bucket_number; i++){
		the_table -> buckets [i].item_keys = NULL;
		the_table -> buckets [i].current_items = 0;
		the_table -> buckets [i].space_size = 0;
	}

	return 0;
}

void _gehash_resize_bucket(gehash_t * the_table , int bucket_no, char is_small_table)
{
	int new_bucket_length;
	struct gehash_bucket * current_bucket = &(the_table -> buckets [bucket_no]);
	gehash_key_t * new_item_keys;
	gehash_data_t * new_item_values;

	if (current_bucket->space_size<1)
	{
		if(is_small_table)
			new_bucket_length = (int)max(15,current_bucket->space_size*1.5+1);
		else
			new_bucket_length = (int)(GEHASH_BUCKET_LENGTH*(1.+3.*pow((rand()*1./RAND_MAX),30)));


		new_item_keys = (gehash_key_t *) malloc(sizeof(gehash_key_t) * new_bucket_length);
		new_item_values = (gehash_data_t *) malloc(sizeof(gehash_data_t) * new_bucket_length);
		bzero(new_item_keys, sizeof(gehash_key_t) * new_bucket_length);
		bzero(new_item_values, sizeof(gehash_data_t) * new_bucket_length);

	//printf("REALLOCATE %d\n", new_bucket_length);

		if(!new_item_values || !new_item_keys)
		{
			printf("\nThe system cannot allocate virtual memory for the index. It seems that you specified a too large memory limit on a 32-bit computer. You may reduce the limit and try again.\n");
			exit(1);
		}

		current_bucket->item_keys = new_item_keys;
		current_bucket->item_values = new_item_values;
		current_bucket->space_size = new_bucket_length;

	}
	else
	{
		int test_start = rand() % the_table-> buckets_number;
		int i;
		int need_allowc = 1;
		for (i =test_start; i<test_start+10000; i++)
		{
			if (i==bucket_no)continue;
			if (i>=the_table-> buckets_number)
			{
				test_start = rand() % the_table-> buckets_number;
				i = test_start;
				continue;
			} 

			if(the_table-> buckets[i].space_size > current_bucket->current_items+1 && current_bucket->space_size > the_table-> buckets[i].current_items+1 && current_bucket->current_items > the_table-> buckets[i].current_items)
			{
				int j;
				gehash_key_t t_key;
				gehash_data_t t_d;

				for (j=0; j<current_bucket->current_items; j++)
				{
					if (j< the_table -> buckets[i].current_items)
					{
						t_key = current_bucket->item_keys[j];
						current_bucket->item_keys[j] = the_table-> buckets[i].item_keys[j];
						the_table-> buckets[i].item_keys[j] = t_key;


						t_d = current_bucket->item_values[j];
						current_bucket->item_values[j] = the_table-> buckets[i].item_values[j];
						the_table-> buckets[i].item_values[j] = t_d;
					}
					else
					{
						the_table-> buckets[i].item_keys[j] = current_bucket->item_keys[j];
						the_table-> buckets[i].item_values[j] = current_bucket->item_values[j];
					}
				}

				gehash_key_t * p_key;
				gehash_data_t * p_d;
				int i_size;

				p_key = the_table-> buckets[i].item_keys ;
				the_table-> buckets[i].item_keys = current_bucket->item_keys;
				current_bucket->item_keys = p_key;

				p_d = the_table-> buckets[i].item_values ;
				the_table-> buckets[i].item_values = current_bucket->item_values;
				current_bucket->item_values = p_d;

				i_size = the_table-> buckets[i].space_size;
				the_table-> buckets[i].space_size =current_bucket->space_size;
				current_bucket->space_size=i_size;
				
				need_allowc = 0;
				break;
			}
		} 
		if(need_allowc)
		{
			if(is_small_table)
				new_bucket_length = (int)max(15,current_bucket->space_size*1.5+1);
			else
				new_bucket_length = (int)(current_bucket->space_size*1.5);
			new_item_keys = (gehash_key_t *) malloc(sizeof(gehash_key_t) * new_bucket_length);
			new_item_values = (gehash_data_t *) malloc(sizeof(gehash_data_t) * new_bucket_length);
			bzero(new_item_keys, sizeof(gehash_key_t) * new_bucket_length);
			bzero(new_item_values, sizeof(gehash_data_t) * new_bucket_length);

			if(!new_item_values || !new_item_keys)
			{
				printf("\nThe system cannot allocate virtual memory for the index. It seems that you specified a too large memory limit on a 32-bit computer. You may reduce the limit and try again.\n");
				exit(1);
			}

			memcpy(new_item_keys, current_bucket->item_keys, current_bucket->current_items*sizeof(gehash_key_t));
			memcpy(new_item_values, current_bucket->item_values, current_bucket->current_items*sizeof(gehash_data_t));

			free(current_bucket->item_keys);
			free(current_bucket->item_values);

			current_bucket->item_keys = new_item_keys;
			current_bucket->item_values = new_item_values;
			current_bucket->space_size = new_bucket_length;


		}
	}
	//printf ("New length = %d\n", new_bucket_length);

}

struct gehash_bucket * _gehash_get_bucket(gehash_t * the_table, gehash_key_t key)
{
	int bucket_number;
	bucket_number = _gehash_hash(key) % the_table -> buckets_number;
	return  &(the_table -> buckets [bucket_number]);
}

void gehash_insert(gehash_t * the_table, gehash_key_t key, gehash_data_t data)
{
	struct gehash_bucket * current_bucket;

	current_bucket = _gehash_get_bucket (the_table, key);
	//printf("P=%lld\n",(long long)(the_table->mem_keys - current_bucket->item_keys));
	if (current_bucket->current_items >= current_bucket->space_size)
	{

	        int bucket_number;
        	bucket_number = _gehash_hash(key) % the_table -> buckets_number;

		_gehash_resize_bucket(the_table, bucket_number, the_table->is_small_table);
	}
	current_bucket->item_keys[current_bucket->current_items] = key;
	current_bucket->item_values[current_bucket->current_items] = data;
	current_bucket->current_items ++;
	the_table ->current_items ++;
}

#define _index_vote(key) (((unsigned int)(key))%GENE_VOTE_TABLE_SIZE)
#define _index_vote_tol(key) (((unsigned int)(key)/16)%GENE_VOTE_TABLE_SIZE)

#define JUMP_GAP 3

#define is_quality_subread(scr)	((scr)>15?1:0)

size_t gehash_go_q_tolerable(gehash_t * the_table, gehash_key_t key, int offset, int read_len, int is_reversed, gene_vote_t * vote, int is_add, gene_vote_number_t weight, gene_quality_score_t quality ,int max_match_number, int indel_tolerance, int subread_number, int max_error_bases, unsigned int low_border, unsigned int high_border)
{
	char error_pos_stack[10];	// max error bases = 10;

	gehash_key_t mutation_key;
	int ret = 0;

	ret+=gehash_go_q(the_table, key, offset, read_len, is_reversed, vote, is_add,  weight,  quality , max_match_number, indel_tolerance, subread_number, low_border, high_border);
	int error_bases ;
	for (error_bases=1; error_bases <= max_error_bases; error_bases++)
	{
		int i, j;
		for(i=0; i<error_bases; i++)
			error_pos_stack [i] = i;
		while (1)
		{

			char mutation_stack[10];
			memset(mutation_stack, 0 , error_bases);
			while(1)
			{
				int is_end2 = 1;
				for (i = error_bases-1; i>=0; i--)
					if (mutation_stack[i] <= 3)
					{
						int base_to_change=-1;
						mutation_key = key;
		
						for (j = 0; j < error_bases; j++)
						{
							base_to_change = error_pos_stack[j];
							int new_value = mutation_stack[j];
							mutation_key = mutation_key & ~(0x3 << (2*base_to_change));
							mutation_key = mutation_key | (new_value << (2*base_to_change));
						}
						for (j = i+1; j<error_bases; j++)
							mutation_stack[j] = 0;
						is_end2 = 0;
					//printf ("K=%u VS M=%u @J=%d , BC=%d\n", key, mutation_key, j, base_to_change );
						if(key != mutation_key)
							ret+=gehash_go_q(the_table, mutation_key, offset, read_len, is_reversed, vote, is_add,  weight,  quality , max_match_number, indel_tolerance, subread_number, low_border, high_border);
						mutation_stack[i] ++;
						break;
					}
				if (is_end2)break;
			}


			int is_end = 1;
			for (i = error_bases-1; i>=0; i--)
				if (error_pos_stack[i] < 16 - (error_bases - i))
				{
					error_pos_stack[i] ++;
					for (j = i+1; j<error_bases; j++)
						error_pos_stack[j] = error_pos_stack[i] + (j-i);
					is_end = 0;
					break;
				}

			if(is_end) break;
		}
	}
	return ret;
}


size_t gehash_go_q(gehash_t * the_table, gehash_key_t key, int offset, int read_len, int is_reversed, gene_vote_t * vote, int is_add, gene_vote_number_t weight, gene_quality_score_t quality ,int max_match_number, int indel_tolerance, int subread_number, unsigned int low_border, unsigned int high_border)
{
	struct gehash_bucket * current_bucket;
	int i, items;
	gehash_key_t * keyp, *endkp, *endp12;
	int match_start, match_end;

	current_bucket = _gehash_get_bucket (the_table, key);
	items = current_bucket -> current_items;


	if(1)
	{
		int citems = items/4;
		keyp = current_bucket -> item_keys + items/2;
		endkp = keyp;
		for(i=0; i<5; i++)
		{
			if(*(keyp+citems) < key)
				keyp += citems;
			else if(*(keyp) >= key)
				keyp -= citems;

			citems = citems/2;
		}
	//	if(*(keyp) >= key || keyp < current_bucket -> item_keys || keyp >= current_bucket -> item_keys + items)
		if(*(keyp) >= key)
			keyp = current_bucket -> item_keys;
	}
	else
		keyp = current_bucket -> item_keys;

	endkp = current_bucket -> item_keys  + items;
	endp12 = endkp - JUMP_GAP;

	while(keyp < endp12 && *(keyp+JUMP_GAP) < key)
		keyp += JUMP_GAP;

	while(*keyp < key && keyp < endkp)
		keyp ++;

	if(keyp == endkp)return 0;

	match_start = keyp - current_bucket -> item_keys;

	gehash_key_t * tk = keyp;

	while(*keyp == key && keyp < endkp)
	{
		keyp ++;
		if(keyp - tk > max_match_number)
			return 0;
	}

	match_end = keyp - current_bucket -> item_keys;

	endkp = current_bucket -> item_values+match_end;

	if (indel_tolerance <1)
		for (keyp = current_bucket -> item_values+match_start; keyp < endkp ; keyp++)
		{
		//add_gene_vote_weighted(vote, (*keyp) - offset, is_add, weight);
			unsigned int kv = (*keyp) - offset;
			int offsetX = _index_vote(kv);
			int datalen = vote -> items[offsetX];
			unsigned int * dat = vote -> pos[offsetX];
			#ifdef MAKE_FOR_EXON
			short offset_from_5 = is_reversed?(read_len - offset - 16):offset ; 
			#endif
		

			if ((*keyp ) < offset)
				continue;

			if (kv < low_border || kv > high_border)
				continue;

			for (i=0;i<datalen;i++)
			{
				if (dat[i] == kv)
				{
					gene_vote_number_t test_max = (vote->votes[offsetX][i]);
					test_max += weight;
					vote->votes[offsetX][i] = test_max;
					vote->quality[offsetX][i] += quality;
					#ifdef MAKE_FOR_EXON
					if (offset_from_5 <  vote->coverage_start [offsetX][i])
						vote->coverage_start [offsetX][i] = offset_from_5;
					if (offset_from_5 +16 > vote->coverage_end [offsetX][i])
						vote->coverage_end [offsetX][i] = offset_from_5+16;
					#endif


					if((test_max > vote->max_vote) || (test_max == vote->max_vote && vote->quality[offsetX][i] > vote->max_quality))
					{
						vote->max_mask = vote->masks[offsetX][i];
						vote->max_vote = test_max;
						vote->max_position = kv;
						vote->max_quality = vote->quality[offsetX][i];
						#ifdef MAKE_FOR_EXON
						vote->max_coverage_start = vote->coverage_start [offsetX][i];
						vote->max_coverage_end = vote->coverage_end [offsetX][i];
						#endif

					}
					else if (test_max == vote->max_vote && vote->quality[offsetX][i] == vote->max_quality)
					{
						//printf("\nBREAK EVEN DETECTED AT SORTED TABLE: %u and %u\n", vote->max_position , kv);
						vote->max_mask |= IS_BREAKEVEN_READ ;
					}
					i = 9999999;
				}
			}

			if (i < 9999999 && is_add && datalen<GENE_VOTE_SPACE)
			{
				vote -> items[offsetX] ++;
				dat[i] = kv;
				vote->votes[offsetX][i]=weight;
				vote->quality[offsetX][i]=quality;
				vote->masks[offsetX][i]= (is_reversed?IS_NEGATIVE_STRAND:0);
				#ifdef MAKE_FOR_EXON
				vote->coverage_start [offsetX][i] = offset_from_5;
				vote->coverage_end [offsetX][i] = offset_from_5+16;
				#endif

				if(vote->max_vote==0)
				{
					#ifdef MAKE_FOR_EXON
					vote->max_coverage_start = offset_from_5;
					vote->max_coverage_end = offset_from_5+16;
					#endif
					vote->max_mask =  (is_reversed?IS_NEGATIVE_STRAND:0);
					vote->max_vote = weight;
					vote->max_position = kv;
					vote->max_quality = quality;
				}
			}
		}
	else
		// We duplicated all codes for indel_tolerance >= 1 for the minimal impact to performance.
		for (keyp = current_bucket -> item_values+match_start; keyp < endkp ; keyp++)
		{
			unsigned int kv = (*keyp) - offset;
			int iix;
			#ifdef MAKE_FOR_EXON
			short offset_from_5 = is_reversed?(read_len - offset - 16):offset ; 
			#endif

			if ((*keyp ) < offset)
				continue;

			if (kv < low_border || kv > high_border)
				continue;
			//if (kv == 2720480913)
			//printf ("kv=%u ; OFF=%d\n", kv, offset);

			for(iix = 0; iix<=3; iix++)
			{
				int ii;
				if (iix) ii = iix*2 - 3;
				else ii = 0;

				int offsetX = _index_vote_tol(kv+16*ii);
				//printf("X=%d\tii=%d\n", offsetX, ii);
				int datalen = vote -> items[offsetX];
				unsigned int * dat = vote -> pos[offsetX];

				for (i=0;i<datalen;i++)
				{
					long long int dist0 = kv;
					dist0 -= dat[i];
					if (abs(dist0) <= indel_tolerance)
					{
						gene_vote_number_t test_max = (vote->votes[offsetX][i]);
						test_max += weight;
						vote -> votes[offsetX][i] = test_max;
						vote -> quality[offsetX][i] += quality;
						vote -> last_offset[offsetX][i]=offset;

						//printf("\nFOUND-OLD! %u, V=%d\n", kv, test_max);
						#ifdef MAKE_FOR_EXON
						if (offset_from_5 <  vote->coverage_start [offsetX][i])
							vote->coverage_start [offsetX][i] = offset_from_5;
						if (offset_from_5 +16 > vote->coverage_end [offsetX][i])
							vote->coverage_end [offsetX][i] = offset_from_5+16;
						#endif

						int toli;
						for(toli=0; toli<indel_tolerance*3-3; toli+=3)
						{
							if (! vote -> indel_recorder[offsetX][i][toli+3])break;
						}
	
						if (dist0!=0 || vote->current_indel_cursor[offsetX][i] !=0)
						{
							if (dist0 !=  vote->current_indel_cursor[offsetX][i])
							{
								toli +=3;
								//printf("NEW TOLI=%d  DIST=%d  SUBR=%d\n", toli, dist0, subread_number);
								if (toli < indel_tolerance*3)
								{
									vote -> indel_recorder[offsetX][i][toli] = subread_number+1; 
									vote -> indel_recorder[offsetX][i][toli+2] = dist0; 
									
									if (toli < indel_tolerance*3-3) vote -> indel_recorder[offsetX][i][toli+3]=0;
								}
								vote->current_indel_cursor [offsetX][i] = (char)dist0;
							}
						}
						//printf("OLD TOLI=%d  DIST=%d  SUBR=%d\n", toli, dist0, subread_number);
						vote -> indel_recorder[offsetX][i][toli+1] = subread_number+1;

						if(test_max > vote->max_vote || (test_max == vote->max_vote && vote -> quality[offsetX][i] > vote->max_quality))
						{
							#ifdef MAKE_FOR_EXON
							vote->max_coverage_start = vote->coverage_start [offsetX][i];
							vote->max_coverage_end = vote->coverage_end [offsetX][i];
							#endif

							memcpy(vote->max_indel_recorder, vote->indel_recorder[offsetX][i], sizeof(char)*3 * MAX_INDEL_TOLERANCE);
							vote->max_vote = test_max;
							vote->max_quality = vote -> quality[offsetX][i];

							//if (vote->max_mask & IS_BREAKEVEN_READ) printf("\nREPLACED A BREAK EVEN! %u to %u\n",vote->max_position , dat[i]);
							vote->max_position = dat[i];
							vote->max_mask = vote->masks[offsetX][i];
						//	printf ("SETMAX %d  %d~%d @%u\n", test_max, vote->max_coverage_start, vote->max_coverage_end, vote->max_position);
						}
						else if (test_max == vote->max_vote && vote->quality[offsetX][i] == vote->max_quality)
						{
							//printf("\nBREAK EVEN DETECTED AT SORTED TABLE: %u (kept) and %u\n", vote->max_position , kv);
							vote->max_mask |= IS_BREAKEVEN_READ ;
						}
						i = 9999999;
					}
				}
				if (i==9999999){
					break;
				}
			}

			if (i < 9999999)
			{
				//printf("\nFOUND-NEW! %u\n", kv);
				int offsetX2 = _index_vote_tol(kv);
				int datalen2 = vote -> items[offsetX2];
				unsigned int * dat2 = vote -> pos[offsetX2];

				if (is_add && datalen2<GENE_VOTE_SPACE)
				{
					vote -> items[offsetX2] ++;
					dat2[datalen2] = kv;
					vote -> masks[offsetX2][datalen2]=(is_reversed?IS_NEGATIVE_STRAND:0);
					vote -> quality[offsetX2][datalen2]=quality;
					vote -> votes[offsetX2][datalen2]=weight;
					vote -> last_offset[offsetX2][datalen2]=offset;

					// data structure of recorder:
					// {unsigned char subread_start; unsigned char subread_end, char indel_offset_from_start}
					// All subread numbers are added with 1 for not being 0.

					vote -> indel_recorder[offsetX2][datalen2][0] = subread_number+1;
					vote -> indel_recorder[offsetX2][datalen2][1] = subread_number+1;
					vote -> indel_recorder[offsetX2][datalen2][2] = 0;
					vote -> indel_recorder[offsetX2][datalen2][3] = 0;
					vote->current_indel_cursor [offsetX2][datalen2] = 0;
					#ifdef MAKE_FOR_EXON
					vote->coverage_start [offsetX2][datalen2] = offset_from_5;
					vote->coverage_end [offsetX2][datalen2] = offset_from_5+16;
					#endif

					if (vote->max_vote==0)
					{
						#ifdef MAKE_FOR_EXON
						vote->max_coverage_start = offset_from_5;
						vote->max_coverage_end = offset_from_5+16;
						#endif
	
						vote->max_vote = weight;
						vote->max_position = kv;
						vote->max_mask = (is_reversed?IS_NEGATIVE_STRAND:0);
						vote->max_quality = quality;
						memcpy(vote->max_indel_recorder, vote->indel_recorder[offsetX2][datalen2], sizeof(char)*4);
					}
				}
			}
		}
	
		

	return match_end-match_start;
}

int gehash_exist(gehash_t * the_table, gehash_key_t key)
{
	struct gehash_bucket * current_bucket;
	int  items;
	gehash_key_t * keyp, *endkp;

	current_bucket = _gehash_get_bucket (the_table, key);
	items = current_bucket -> current_items;

	if(items <1)return 0;

	keyp = current_bucket -> item_keys;
	endkp = keyp + items;

	while(1)
	{
		if(*keyp == key)
			return 1;
		if(++keyp >= endkp) break;
	}
	return 0;
}
size_t gehash_get(gehash_t * the_table, gehash_key_t key, gehash_data_t * data_result, size_t max_result_space)
{
	struct gehash_bucket * current_bucket;
	size_t matched;
	int items;
	gehash_key_t * keyp, *endkp;

	if(max_result_space<1)
		return 0;

	current_bucket = _gehash_get_bucket (the_table, key);

	matched = 0;
	items = current_bucket -> current_items;
	keyp = current_bucket -> item_keys;
	endkp = keyp + items;

	while(1)
	{
		if(*keyp == key)
		{
			data_result [matched] = current_bucket -> item_values[keyp-current_bucket -> item_keys];
			matched +=1;
			if(matched >= max_result_space)
				break;
		}
		keyp +=1;
		if(keyp >= endkp)
			break;
	}
	return matched;
}


size_t gehash_remove(gehash_t * the_table, gehash_key_t key)
{
	struct gehash_bucket * current_bucket;
	int i;
	size_t removed;

	current_bucket = _gehash_get_bucket (the_table, key);	

	if(current_bucket -> current_items < 1)
		return 0;

	removed = 0;
	for(i=0; ; i++)
	{
		while(current_bucket -> item_keys [i+removed] == key && i+removed < current_bucket -> current_items)
			removed += 1;

		if(i+removed >= current_bucket -> current_items)
			break;

		if(removed)
		{
			current_bucket -> item_keys [i] = 
				current_bucket -> item_keys [i + removed];

			current_bucket -> item_values [i] = 
				current_bucket -> item_values [i + removed];
		}

	}

	current_bucket -> current_items -= removed;
	the_table-> current_items -= removed;
/*
	if (current_bucket -> space_size - current_bucket -> current_items > 50)
	{
		gehash_key_t * new_item_keys;
		gehash_data_t * new_item_values;

		new_item_keys = (gehash_key_t*)malloc(sizeof(gehash_key_t) * current_bucket -> current_items);
		memcpy(new_item_keys, current_bucket -> item_keys, sizeof(gehash_key_t) * current_bucket -> current_items);
		free(current_bucket -> item_keys);
		current_bucket -> item_keys = new_item_keys;

		new_item_values = (gehash_data_t*)malloc(sizeof(gehash_data_t) * current_bucket -> current_items);
		memcpy(new_item_values, current_bucket -> item_values, sizeof(gehash_data_t) * current_bucket -> current_items);
		free(current_bucket -> item_values);
		current_bucket -> item_values = new_item_values;

		current_bucket -> space_size  = current_bucket -> current_items;
	}
*/

	return removed;
}

size_t gehash_get_hpc(gehash_t * the_table, gehash_key_t key, gehash_data_t * data_result, size_t max_result_space)
{
	return -1;
}


void gehash_insert_sorted(gehash_t * the_table, gehash_key_t key, gehash_data_t data)
{
/*
        struct gehash_bucket * current_bucket;
        int search_start = 0, search_end;

        current_bucket = _gehash_get_bucket (the_table, key);
        if (current_bucket->current_items >= current_bucket->space_size)
                _gehash_resize_bucket(current_bucket, , the_table->is_small_table);

        for (;search_start<current_bucket->current_items;search_start++)
        {
                if(current_bucket->item_keys[search_start] >= key)
                        break;
        }



        for (search_end = current_bucket->current_items;search_end>search_start;search_start--)
                current_bucket->item_keys[search_end] = current_bucket->item_keys[search_end-1];

        current_bucket->item_keys[search_start] = key;
        current_bucket->item_values[search_start] = data;

        current_bucket->current_items ++;
*/
	printf("UNIMPLEMENTED! gehash_insert_sorted \n");
}



// Data Struct of dumpping:
// {
//      size_t current_items;
//      size_t buckets_number;
//      struct 
//      {
//              size_t current_items;
//              size_t space_size;
//              gehash_key_t item_keys [current_items];
//              gehash_data_t item_values [current_items]
//      } [buckets_number];
// }
//

inline unsigned int load_int32(FILE * fp)
{
	int ret;
	int read_length;
	read_length = fread(&ret, sizeof(int), 1, fp);
	assert(read_length>0);
	return ret;
}

inline long long int load_int64(FILE * fp)
{
	long long int ret;
	int read_length;
	read_length = fread(&ret, sizeof(long long int), 1, fp);
	assert(read_length>0);
	return ret;
}




int gehash_load(gehash_t * the_table, const char fname [])
{
	int i, read_length;

	FILE * fp = fopen(fname, "rb");
	if (!fp)
	{
		printf ("Table file `%s' is not found.\n", fname);
		return -1;
	}

	the_table -> current_items = load_int64(fp);
	the_table -> buckets_number = load_int32(fp);
	the_table -> buckets = (struct gehash_bucket * )malloc(sizeof(struct gehash_bucket) * the_table -> buckets_number);

	for (i=0; i<the_table -> buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(the_table -> buckets[i]);
		current_bucket -> current_items = load_int32(fp);
		current_bucket -> space_size = load_int32(fp);
		current_bucket -> space_size = current_bucket -> current_items;
		current_bucket -> item_keys = (gehash_key_t *) malloc ( sizeof(gehash_key_t) * current_bucket -> space_size);
		current_bucket -> item_values = (gehash_data_t *) malloc ( sizeof(gehash_data_t) * current_bucket -> space_size);

		//printf("IKV=%d\n", (int)current_bucket -> space_size);

		if(current_bucket -> current_items > 0)
		{
			read_length = fread(current_bucket -> item_keys, sizeof(gehash_key_t), current_bucket -> current_items, fp);
			assert(read_length>0);
			read_length = fread(current_bucket -> item_values, sizeof(gehash_data_t), current_bucket -> current_items, fp);
			assert(read_length>0);
		}

/*
		int j;
		for (j=0; j<current_bucket -> current_items; j++)
		{
			if(current_bucket -> item_keys[j]==4131185942)
			{
				printf("FOUND %u @ %d",4131185942, current_bucket -> item_values[j]);
			}
		}
*/	}

	read_length = fread(&(the_table -> is_small_table), sizeof(char), 1, fp);
	assert(read_length>0);
	fclose(fp);
	return 0;
}

int gehash_dump(gehash_t * the_table, const char fname [])
{
	int i, scroll_counter = 0;
	FILE * fp = fopen(fname, "wb");
	if (!fp)
	{
		printf ("Table file `%s' is not able to open.\n", fname);
		return -1;
	}

	fwrite(& (the_table -> current_items ), sizeof(long long int), 1, fp);
	fwrite(& (the_table -> buckets_number), sizeof(int), 1, fp);


	for (i=0; i<the_table -> buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(the_table -> buckets[i]);
		int ii, jj;
		gehash_key_t tmp_key;
		gehash_data_t tmp_val;


                if(i % 200 == 0)
			print_text_scrolling_bar("Saving index", 1.0*i/the_table -> buckets_number, 80, &scroll_counter);


		if(current_bucket -> current_items>=1)
		{
			for(ii=0;ii<current_bucket -> current_items -1; ii++)
			{
				for (jj = ii+1; jj < current_bucket -> current_items; jj++)
				{
					if (current_bucket -> item_keys[ii] > current_bucket -> item_keys[jj])
					{
						tmp_key = current_bucket -> item_keys[ii];
						current_bucket -> item_keys[ii] = current_bucket -> item_keys[jj];
						current_bucket -> item_keys[jj] = tmp_key;

						tmp_val = current_bucket -> item_values[ii];
						current_bucket -> item_values[ii] = current_bucket -> item_values[jj];
						current_bucket -> item_values[jj] = tmp_val;
					}
				}
			}
		}

//		if (i % 1000 ==0)
//			printf("Sorted %d buckets in %d\n",i, the_table -> buckets_number);

		fwrite(& (current_bucket -> current_items), sizeof(int), 1, fp);
		fwrite(& (current_bucket -> space_size), sizeof(int), 1, fp);
		fwrite(current_bucket -> item_keys, sizeof(gehash_key_t), current_bucket -> current_items, fp);
		fwrite(current_bucket -> item_values, sizeof(gehash_data_t), current_bucket -> current_items, fp);
/*
		int j;
		for (j=0; j<current_bucket -> current_items; j++)
		{
			if(current_bucket -> item_keys[j]==4131185942)
			{
				printf("FOUND %u @ %d\n",4131185942, current_bucket -> item_values[j]);
			}
		}
*/	}

	fwrite(&(the_table -> is_small_table), sizeof(char), 1, fp);
	fclose(fp);
	return 0;
}


void gehash_destory(gehash_t * the_table)
{
	int i;

	for (i=0; i<the_table -> buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(the_table -> buckets[i]);
		if (current_bucket -> space_size > 0)
		{
			free (current_bucket -> item_keys);
			free (current_bucket -> item_values);
		}
	}

	free (the_table -> buckets);

	the_table -> current_items = 0;
	the_table -> buckets_number = 0;
}
