#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "gene-algorithms.h"
#include "gene-value-index.h"
#include "input-files.h"
#include "sorted-hashtable.h"

int ALL_THREADS=1;
int TOTAL_SUBREADS;
int ACCEPT_SUBREADS;
int ACCEPT_MINOR_SUBREADS;
int INDEX_THRESHOLD;
int MAX_PAIRED_DISTANCE = 300;
int MIN_PAIRED_DISTANCE = 200;
int INDEL_TOLERANCE = 0;
int IS_DEBUG = 0;
int QUALITY_SCALE = 0;
int USE_VALUE_ARRAY_INDEX = 0;
int FIRST_READ_REVERSE = 0;
int SECOND_READ_REVERSE = 1;
double reads_density;


struct gene_thread_data
{
	gehash_t table;
	int * offsets;
	int number_t;
};




void print_res(gene_allvote_t *av, gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads)
{
	int i, j, ic=0;
	char inb [1201], nameb[1201], qualityb[1201];
	gene_offset_t offsets;
	unsigned int paired_match = 0;

	load_offsets (&offsets, index_prefix);

	if(all_processed_reads ==0)
	{
		unsigned int last_offset = 0;
		i=0;
		while(offsets.read_offset[i])
		{
			fprintf(out_fp, "@SQ\tSN:%s\tLN:%u\n", offsets.read_name[i], offsets.read_offset[i] - last_offset);
			last_offset = offsets.read_offset[i];
			i++;
		}
	}
	
	i=0;

	gene_input_t * Curr_ginp=ginp2;
	while(1)
	{
		int isOK = 0;
		nameb[0]=0;

		if(ginp2) 
			Curr_ginp = (Curr_ginp==ginp2)?ginp:ginp2;
		else
			Curr_ginp = ginp;

		int rl = geinput_next_read(Curr_ginp, nameb, inb, qualityb);
		if (rl<0){
			break;
		}

		if (nameb[0]==0)
			sprintf(nameb, "read_%llu", i+1+all_processed_reads);

		int Curr_position;
		if(ginp2)
		{
			Curr_position = i*2+(Curr_ginp==ginp2);
		}
		else
			Curr_position = i;

		int flags = 0;
		if (ginp2) flags |= SAM_FLAG_PAIRED_TASK;

		unsigned char votes = av->max_votes[Curr_position];
		if ((votes >= ACCEPT_SUBREADS*22-21.9 ) || (av->masks[Curr_position] & IS_PAIRED_MATCH))
		{
#ifdef REPORT_ALL_THE_BEST
		int k;
		for(k=0; k<av->best_records[i].best_len; k++)
		{
			unsigned int lpos = av->best_records[i].offsets[k];
			char is_counterpart = av->best_records[i].is_reverse[k];
#else

			if(av->masks[Curr_position] & IS_PAIRED_MATCH)
				paired_match++;
			int lpos = av->max_positions[Curr_position];
			char is_counterpart = av->is_counterpart[Curr_position];
#endif


			char * read_name;
			unsigned int read_pos;

			if(is_counterpart + (Curr_ginp==ginp2) == 1)
			{
				int rev_offset = (Curr_ginp->space_type==GENE_SPACE_COLOR && inb[0]>='A' && inb[0]<='Z')?1:0;
				reverse_read(inb+ rev_offset, rl- rev_offset, ginp->space_type);
			}

			if(locate_gene_position(lpos, &offsets, &read_name, &read_pos))
				printf ("ERROR: position out of range:%u\n", lpos); 
			else{
				char cigar_str [2];

				cigar_str[0]=0;
				if (av->masks[Curr_position] & IS_INSERTION)strcat(cigar_str,"I");
				if (av->masks[Curr_position] & IS_DELETION)strcat(cigar_str,"D");
				if (av->masks[Curr_position] & IS_PAIRED_MATCH)
					flags |= SAM_FLAG_MATCHED_IN_PAIR;

				if (is_counterpart + (Curr_ginp==ginp2) == 1)
					flags |= SAM_FLAG_REVERSE_STRAND_MATCHED; 
				else
					flags |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;

				if (ginp2)
				{
					if (Curr_ginp==ginp2)
						flags |= SAM_FLAG_SECOND_READ_IN_PAIR;
					else
						flags |= SAM_FLAG_FIRST_READ_IN_PAIR;

					if (!(av->masks[Curr_position] & IS_PAIRED_MATCH))
					{
						int mate_index = i*2+(Curr_ginp!=ginp2);
						if (av->max_votes[mate_index] < ACCEPT_SUBREADS)
							flags |= SAM_FLAG_MATE_UNMATCHED;
					}
					
				}

				fprintf (out_fp, "%s\t%d\t%s\t%d\t%d\t%dM%s\t*\t0\t0\t%s\t",  nameb, flags, read_name, read_pos+1, votes, rl, cigar_str, inb);
				if(qualityb[0])
					fputs(qualityb, out_fp);
				else
					for(j=0; j<rl; j++)
						fputc('h', out_fp);
				fprintf(out_fp, "\n");
				isOK = 1;
				(*succeed_reads) ++;
			}

#ifdef REPORT_ALL_THE_BEST
		}
#endif

		}

		if(!isOK)
		{
			flags |= SAM_FLAG_UNMAPPED;

			fprintf (out_fp, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t", nameb, flags, inb);
			if(qualityb[0])
				fputs(qualityb, out_fp);
			else
				for(j=0; j<rl; j++)
					fputc('h', out_fp);
			fprintf(out_fp, "\n");
		}
//			printf ("@UNK Q#%d UNKNOW %s\n", i, inb);

		if((!ginp2) || ginp2 == Curr_ginp)
			i++;
		if(i >= processed_reads)break;
		if(i % 10000 ==0 && i>1)
			print_text_scrolling_bar("Saving results", i*1./processed_reads, 80, &ic);
	}
	if (paired_match)
		printf("\nReads referencing paired-end info: %u\n", paired_match);
}



struct gene_thread_data_transport
{
	gehash_t * my_table;
	gene_value_index_t * my_value_array_index;
	int table_no;
	int all_tables;
	char * index_prefix;
	gene_allvote_t * all_vote;
	gene_input_t * ginp;
	gene_input_t * ginp2;
	unsigned long long int base_number;
	unsigned int * processed_reads;
	unsigned int section_length;
	int all_threads;
	int this_thread;

	pthread_spinlock_t * input_data_lock;
	pthread_spinlock_t * init_lock;
};

int run_search(gehash_t * my_table, gene_value_index_t * my_value_array_index , int table_no,  gene_allvote_t * allvote, gene_input_t * ginp,gene_input_t * ginp2, char * index_prefix, unsigned int * processed_reads, long long int base_number, int all_tables, pthread_spinlock_t * input_lock, int my_thread_no, unsigned int section_length)
{
//	FILE * fp;
	char BuffMemory [2500];
	char * InBuff = NULL , * InBuff2 = NULL;
	char BuffMemory2 [2500];
	char * QualityBuff = NULL, * QualityBuff2 = NULL;

	int all_reads = allvote -> max_len;
	int queries = 0;
	double t0=miltime();
	int is_counterpart;
	int good_match = 0;
	float subread_step = 3;
	struct stat read_fstat;
	stat (ginp->filename, &read_fstat);
	long long int read_fsize = read_fstat.st_size;
	
	double local_begin_ftime = miltime();
	int read_len = 0, read2_len = 0;
	is_counterpart = 0;
	if (ginp2)
		all_reads /=2;

//	printf ("I'm the %d-th thread\n", my_thread_no);


	while (1)
	{

		if (is_counterpart)
		{
			
			reverse_read(InBuff, read_len, ginp->space_type);	
			reverse_quality(QualityBuff, read_len);
			if (ginp2)
			{
				reverse_read(InBuff2, read2_len, ginp->space_type);
				reverse_quality(QualityBuff2, read2_len);
			}
		}
		else
		{
			if(input_lock!=NULL)
				pthread_spin_lock(input_lock);
			

			if(*processed_reads < all_reads)
			{
				char namebuf[200];
				InBuff = BuffMemory;
				QualityBuff = BuffMemory2;


				read_len = geinput_next_read(ginp, namebuf, InBuff, QualityBuff);
				//printf ("\nMAJOR: %s\n", namebuf);
				if(read_len>0){
					queries = *processed_reads;
					(*processed_reads ) ++;
					if (ginp->space_type == GENE_SPACE_COLOR && InBuff[0]>='A' && InBuff[0]<='Z')
					{
						InBuff ++;
						if(QualityBuff[0])
							QualityBuff++;
						read_len --;
					}

					if(FIRST_READ_REVERSE)
					{
						reverse_read(InBuff, read_len, ginp->space_type);
						reverse_quality(QualityBuff,read_len);
					}
//						reverse_read(InBuff, read_len, ginp->space_type);
				}
				if(ginp2)
				{
					InBuff2 = BuffMemory + 1250;
					QualityBuff2 = BuffMemory2 + 1250;

					read2_len = geinput_next_read(ginp2, namebuf, InBuff2, QualityBuff2);
				//printf ("MINOR: %s\n", namebuf);

					if (ginp->space_type == GENE_SPACE_COLOR && InBuff2[0]>='A' && InBuff2[0]<='Z')
					{
						InBuff2 ++;
						if(QualityBuff2[0])
							QualityBuff2++;
						read2_len --;
					}
	
					if(SECOND_READ_REVERSE)
					{
						reverse_quality(QualityBuff2, read2_len);
						reverse_read(InBuff2, read2_len, ginp->space_type);
					}
				
				}

				if(input_lock!=NULL)
					pthread_spin_unlock(input_lock);
			}
			else
			{
				if(input_lock!=NULL)
					pthread_spin_unlock(input_lock);
				break;
			}

			if (read_len<0)
				break;
			subread_step = max(3.00001, (read_len-16-GENE_SLIDING_STEP)*1.0/(TOTAL_SUBREADS-1)+0.00001); 
		}

		int i, matched = 0;
		if (ginp2)
		{
			gene_vote_t vote_read1;
			gene_vote_t vote_read2;

			init_gene_vote(&vote_read1);
			init_gene_vote(&vote_read2);
			for(i=0; i<GENE_SLIDING_STEP ; i++)
			{
				int is_paired ;
				int subread_no;

				for (is_paired = 0; is_paired <2; is_paired ++)
				{
					char * CurrInBuff = is_paired?InBuff2:InBuff;
					int Curr_read_len = is_paired?read2_len:read_len;
					gene_vote_t * vote = is_paired?&vote_read2:&vote_read1;
					char * CurrQualityBuff = is_paired?QualityBuff2:QualityBuff;

					for(subread_no=0; ; subread_no++)
					{
						int subread_offset = (int)(subread_step * subread_no); 
						subread_offset -= subread_offset%GENE_SLIDING_STEP -i;
						int subread_offset1 = (int)(subread_step * (subread_no+1));
						subread_offset1 -= subread_offset1%GENE_SLIDING_STEP;
						subread_offset1 += GENE_SLIDING_STEP-1;

						//printf("Q:%d I:%d 2nd:%d SR:%d INV:%d :: Pos:%d Len:%d\n", queries, i, is_paired, subread_no, is_counterpart, subread_offset, Curr_read_len);

						char * subread_string = CurrInBuff + subread_offset;
						gehash_key_t subread_integer = genekey2int(subread_string, ginp->space_type);

						gene_vote_number_t subread_quality = 22.;

						if(CurrQualityBuff[0])
						{
							char * quality_string = CurrQualityBuff + subread_offset;
							subread_quality += get_subread_quality(quality_string, subread_string, QUALITY_SCALE);
						}
						else subread_quality -= .1;
						if(is_valid_subread(subread_string))
							gehash_go_q(my_table, subread_integer , subread_offset , vote , 1 , subread_quality, INDEX_THRESHOLD, INDEL_TOLERANCE);

						if(subread_offset1 >= Curr_read_len -16)
							break;
					}
				}
			}

			gene_vote_number_t numvote_read1, numvote_read2;
			gehash_data_t pos_read1, pos_read2;

			int is_paired_match = 0;

			if(USE_VALUE_ARRAY_INDEX)
				is_paired_match = select_positions_array(InBuff, read_len, InBuff2, read2_len,&vote_read1, &vote_read2, &numvote_read1, &numvote_read2, &pos_read1, &pos_read2, MAX_PAIRED_DISTANCE, MIN_PAIRED_DISTANCE, ACCEPT_SUBREADS, ACCEPT_MINOR_SUBREADS, is_counterpart, my_value_array_index, ginp -> space_type, INDEL_TOLERANCE);
			else
				is_paired_match = select_positions(&vote_read1, &vote_read2, &numvote_read1, &numvote_read2, &pos_read1, &pos_read2, MAX_PAIRED_DISTANCE, MIN_PAIRED_DISTANCE, ACCEPT_SUBREADS, ACCEPT_MINOR_SUBREADS, is_counterpart);
			int pair_selected = 0;

			if (is_paired_match)
			{
				if (	(!(allvote->masks[queries*2] & IS_PAIRED_MATCH)) ||
					((allvote->max_votes[queries*2]+allvote->max_votes[queries*2+1]) < (numvote_read1+numvote_read2)
  					)
				   )
				{
					allvote->max_votes[queries*2] = numvote_read1+1E-4;
					allvote->max_votes[queries*2+1] = numvote_read2+1E-4;
					allvote->is_counterpart[queries*2] = is_counterpart;
					allvote->is_counterpart[queries*2+1] = is_counterpart;
					allvote->max_positions[queries*2] = pos_read1;
					allvote->max_positions[queries*2+1] = pos_read2;
					allvote->masks[queries*2] = vote_read1.max_mask |IS_PAIRED_MATCH;
					allvote->masks[queries*2+1] = vote_read2.max_mask |IS_PAIRED_MATCH;
					pair_selected = 1;
				}
			}

			if ((!(allvote->masks[queries*2] & IS_PAIRED_MATCH)) && !pair_selected)
			{
				pos_read1 = vote_read1.max_position;
				pos_read2 = vote_read2.max_position;
				numvote_read1 = vote_read1.max_vote;
				numvote_read2 = vote_read2.max_vote;

				add_allvote_q(allvote, queries*2  ,  pos_read1 , numvote_read1 , is_counterpart, vote_read1.max_mask);
				add_allvote_q(allvote, queries*2+1,  pos_read2 , numvote_read2 , is_counterpart, vote_read2.max_mask);
			}
		}
		else if (INDEL_TOLERANCE < 1)
		{
			gene_vote_t vote;
			for(i=0; i<GENE_SLIDING_STEP ; i++)
			{
				int subread_no;
				int all_break = 0;
				int total_results=0;

				init_gene_vote(&vote);
	
				for(subread_no=0; ; subread_no++)
				{
					int subread_offset = (int)(subread_step * subread_no); 
					subread_offset -= subread_offset%GENE_SLIDING_STEP -i;
					int subread_offset1 = (int)(subread_step * (subread_no+1));
					subread_offset1 -= subread_offset1%GENE_SLIDING_STEP;
					subread_offset1 += GENE_SLIDING_STEP-1;

					char * subread_string = InBuff + subread_offset;
					gehash_key_t subread_integer = genekey2int(subread_string, ginp->space_type);

					gene_vote_number_t subread_quality = 22;

					if(QualityBuff[0])
					{
						char * quality_string = QualityBuff + subread_offset;
						subread_quality += get_subread_quality(quality_string, subread_string, QUALITY_SCALE);
					}
					else subread_quality -= .1;

					if(is_valid_subread(subread_string))
						total_results = gehash_go_q(my_table, subread_integer , subread_offset , &vote , (subread_no < (TOTAL_SUBREADS - ACCEPT_SUBREADS)) , subread_quality, INDEX_THRESHOLD, INDEL_TOLERANCE);

					if (vote.max_vote <= (subread_no - (TOTAL_SUBREADS - ACCEPT_SUBREADS))*22-20)
						break;

					if(subread_offset1 >= read_len -16)
					{
						if(USE_VALUE_ARRAY_INDEX)
						{
							gehash_data_t max_position;
							gene_vote_number_t max_vote;
							char max_mask;

							final_matchingness_scoring(InBuff, QualityBuff, read_len, &vote, &max_position, &max_vote, &max_mask, my_value_array_index, ginp -> space_type, INDEL_TOLERANCE);
							add_allvote_q(allvote, queries,  max_position , max_vote , is_counterpart, max_mask);
						}
						else
							add_allvote_q(allvote, queries,  vote.max_position , vote.max_vote , is_counterpart, vote.max_mask);
						matched =1;
						break;
					}
				}

				good_match += matched;

				if(all_break) break;
	
			}
		}
		else
		{
			gene_vote_t vote;
			int subread_no;
			init_gene_vote(&vote);
			for(subread_no=0; ; subread_no++)
			{
				int total_results=0;

				int subread_offset1 = (int)(subread_step * (subread_no+1));
				subread_offset1 -= subread_offset1%GENE_SLIDING_STEP;
				subread_offset1 += GENE_SLIDING_STEP-1;

				for(i=0; i<GENE_SLIDING_STEP ; i++)
				{
					int subread_offset = (int)(subread_step * subread_no); 
					subread_offset -= subread_offset%GENE_SLIDING_STEP -i;

					char * subread_string = InBuff + subread_offset;

					gehash_key_t subread_integer = genekey2int(subread_string, ginp->space_type);

					gene_vote_number_t subread_quality = 22;

					if(QualityBuff[0])
					{
						char * quality_string = QualityBuff + subread_offset;
						subread_quality += get_subread_quality(quality_string, subread_string, QUALITY_SCALE);
					}
					else subread_quality -= .1;

					if(is_valid_subread(subread_string))
						total_results = gehash_go_q(my_table, subread_integer , subread_offset , &vote , 1 , subread_quality, INDEX_THRESHOLD, INDEL_TOLERANCE);

				}
				if(subread_offset1 >= read_len -16)
					break;

			}

			if(vote.max_vote > (subread_no - (TOTAL_SUBREADS - ACCEPT_SUBREADS))*22.)
			{
				if(USE_VALUE_ARRAY_INDEX)
				{
					gehash_data_t max_position;
					gene_vote_number_t max_vote;
					char max_mask;

					final_matchingness_scoring(InBuff, QualityBuff, read_len, &vote, &max_position, &max_vote, &max_mask, my_value_array_index, ginp -> space_type, INDEL_TOLERANCE);
					add_allvote_q(allvote, queries,  max_position , max_vote , is_counterpart, max_mask);
				}
				else
					add_allvote_q(allvote, queries,  vote.max_position , vote.max_vote , is_counterpart, vote.max_mask);

				matched =1;

				good_match += matched;
			}
	
		}


		if (my_thread_no==0 && queries % 10000 == 0 && !is_counterpart)
		{
			if(table_no == 0)
			{
				long long int current_reads = base_number + queries;
				long long int fpos = ftello(ginp->input_fp);
				reads_density = fpos*1.0/current_reads; 
			}
			if(IS_DEBUG && queries % 100000==0)
				printf("@LOG Done %d/%d, good %d, last time %f\n",queries, table_no, good_match, miltime() - t0);
			else
			{
				long long int all_steps = (read_fsize*1.0/reads_density) * all_tables;
				int remaining_load_libs = (all_tables - table_no);
				remaining_load_libs +=  all_tables * (int)(((read_fsize*1.0/reads_density) - base_number)/all_reads);
				long long int finished_steps =  ((section_length)*table_no + base_number*all_tables+queries);
				double finished_rate = finished_steps*1.0 / all_steps;
				double reads_per_second =  queries*1.0 / all_tables / (miltime()- local_begin_ftime);
				double expected_seconds = ((1.-finished_rate) * all_steps)/all_tables / reads_per_second + remaining_load_libs * 50 + (int)(((read_fsize*1.0/reads_density) - base_number)/all_reads)*100;
				if(queries>1)
				print_running_log(finished_rate, reads_per_second, expected_seconds, (unsigned long long int)all_steps / all_tables);
			}
			fflush(stdout);
			t0 = miltime();
		}
	
		if (is_counterpart)
		{
			if(queries >= all_reads)
				break;
		}

		is_counterpart = !is_counterpart;
	}

	return 0;
}


void * run_search_thread(void * parameters)
{
	struct gene_thread_data_transport * data_param = parameters;
	int thid = data_param->this_thread;

	pthread_spin_unlock(data_param -> init_lock);

	run_search(data_param->my_table, data_param->my_value_array_index, data_param->table_no, data_param->all_vote, data_param->ginp , data_param->ginp2, data_param->index_prefix, data_param->processed_reads, data_param->base_number, data_param->all_tables, data_param->input_data_lock, thid, data_param-> section_length);
	return NULL;
}



// This function search a segment of reads (length = read_number) 
// It returns the number of reads that were really processed;
int run_search_index(gene_input_t * ginp, gene_input_t * ginp2, char * index_prefix, gene_allvote_t * all_vote, FILE * out_fp, unsigned long long int base_number, int all_tables, unsigned long long int *succeed_reads)
{
	unsigned int tabno=0;
	unsigned int processed_reads = 0, section_length = 0;
	unsigned long long int current_fp = ftello(ginp -> input_fp);
	unsigned long long int current_fp2 = 0;
	struct stat filestat;
	int stat_ret;
	char table_fn [300];
	gehash_t my_raw_table;
	gehash_t * my_table = &my_raw_table;
	gene_value_index_t value_array_index;

	if(ginp2)
		current_fp2 = ftello(ginp2 -> input_fp);

	while (1)
	{
		sprintf(table_fn, "%s.%02d.%c.tab", index_prefix, tabno, ginp->space_type==GENE_SPACE_COLOR?'c':'b');

// Test and load the index partition
		stat_ret = stat(table_fn, &filestat);
		if (stat_ret !=0)
			break;

		if (IS_DEBUG)
			printf ("@LOG Loading table from %s\n", table_fn);
		else
			printf ("Loading the %02d-th index file ...                                              \r", tabno+1);
		fflush(stdout);

		gehash_load(my_table, table_fn);
		if(USE_VALUE_ARRAY_INDEX)
		{
			
			sprintf(table_fn, "%s.%02d.%c.array", index_prefix, tabno, ginp->space_type==GENE_SPACE_COLOR?'c':'b');
			stat_ret = stat(table_fn, &filestat);
			if (stat_ret !=0)
			{
				printf("The specified index does not contain any value array data.\n");
				return -1;
			}

			gvindex_load(&value_array_index,table_fn);
		}
		processed_reads = 0;

// Run the search algorithm on a part of the index
		if(ALL_THREADS <2)
			run_search(my_table, &value_array_index, tabno, all_vote, ginp , ginp2, index_prefix, &processed_reads, base_number, all_tables, NULL /*the data lock is null*/, 0  /*I'm the 0-th thread*/, section_length);
		else
		{
			int i; 
			struct gene_thread_data_transport data_param;
			pthread_t runners [ALL_THREADS];
			pthread_spinlock_t  data_lock;
			pthread_spinlock_t  init_lock;

			data_param.my_table = my_table;
			data_param.my_value_array_index = &value_array_index;
			data_param.table_no = tabno;
			data_param.all_tables = all_tables;
			data_param.index_prefix = index_prefix;
			data_param.all_vote = all_vote;
			data_param.base_number = base_number;
			data_param.all_threads = ALL_THREADS;
			data_param.input_data_lock = &data_lock;
			data_param.init_lock = &init_lock;
			data_param.processed_reads = &processed_reads;
			data_param.section_length = section_length;
			data_param.ginp = ginp;
			data_param.ginp2 = ginp2;

			pthread_spin_init(&data_lock, PTHREAD_PROCESS_PRIVATE);
			pthread_spin_init(&init_lock, PTHREAD_PROCESS_PRIVATE);
			pthread_spin_lock(&init_lock);
			for (i=0; i< ALL_THREADS; i++)
			{
				data_param.this_thread = i;
				pthread_create(runners+i, NULL, run_search_thread, &data_param);
				pthread_spin_lock(&init_lock);
			}

			for (i=0; i< ALL_THREADS; i++)
				pthread_join(*(runners+i), NULL);
			pthread_spin_destroy(&data_lock);
			pthread_spin_destroy(&init_lock);
		}
		tabno ++;

		if (section_length < 1)
			section_length = processed_reads;

		gehash_destory_fast(my_table);
		fseeko(ginp -> input_fp, current_fp, SEEK_SET);
		if (ginp2)
			fseeko(ginp2 -> input_fp, current_fp2, SEEK_SET);

	}

	print_res(all_vote, ginp, ginp2, out_fp, index_prefix, processed_reads, base_number, succeed_reads);

	return processed_reads;
	
}

void usage(char * execname)
{
	puts("Usage:");
	puts(" ./subread-align [options] -i <index_name> -r <input> -o <output>");
	puts("");
	puts("Basic arguments:");
	puts("    -i --index     <index>\t name of the index, same as that for the index builder.");
	puts("    -r --read      <input>\t name of an input file(FASTQ/FASTA format), either in the base-space or in the color-space.");
	puts("    -o --output    <output>\t name of the output file(SAM format)");
	puts("");
	puts("Optional arguments:");
	puts("    -n --subreads  <int>\t optional, number of subreads selected from each read for mapping, 10 by default");
	puts("    -m --minmatch  <int>\t optional, minimal number of subreads which have the consensus mapping location, 3 by default");
	puts("    -T --threads   <int>\t optional, number of threads/CPUs used for mapping the reads, 1 by default");
	puts("    -I --indel     <int>\t optional, the maximum number of bases for insertion/deletion, 0 by default");
	puts("    -Q --quality   <l:e:n>\t optional, quality scale: l for linear, e for exponential, n for none.");
	puts("    -a --basewise       \t optional, using the base-wise quality index; the index must be built with a -a option.");
	puts("");
	puts("Paired-end alignment arguments:");
	puts("    -R --read2     <input>\t optional, the second input file; using this argument to activate paired-end alignment");
	puts("    -p --minmatch2 <int>\t optional, the `-m' option for the read receiving less votes in a pair to be accepted, 1 by default");
	puts("    -d --mindist   <int>\t optional, the minimum distance between two reads in a pair, 200 by default");
	puts("    -D --maxdist   <int>\t optional, the maximum distance between two reads in a pair, 300 by default");
	puts("    -S --order     <ff:fr:rf> \t optional, specifying if the first/second reads are forward or reversed, 'fr' by default.");
	puts("");
	puts("Example:");
	puts(" ./subread-align -i my_index -r reads.fastq -o my_result.sam ");
	puts("");
	puts("Description:");
	puts("  To map a read to the reference seqeunce, a number of equally spaced subreads of 16 bases long will be selected from the read. These reads will be mapped to the reference sequence to determine whether the read can be successfully mapped and the mapping location.");
	puts("  The number of subreads seletecd from the read is specified by the -n argument. When mapping each subread to the reference sequence, no mismatches are allowed.");
	puts("  The mapping locations of each subread are recorded. A consensus mapping location is determined from all candidate locations, to which the majority of selected subreads were mapped. To successfully map a read, its consensus mapping location has to be mapped by at least a threshold number of subreads, which is specified by the -m argument.");
	puts("  The running time for aligning 10 million 100-bp reads to human genome is estimated to be about half an hour on a desktop. When running with 4 threads, the running time will be reduced to around 10 minutes.");
	puts("  Paired-end alignment is activated by specifying both -r and -R arguments. The two files should contain the same number of reads; two reads at the same line in the two files are considered as a paired observation.");
	puts("");

}

static struct option long_options[] =
{
	{"basewise", no_argument, &USE_VALUE_ARRAY_INDEX, 1},
	{"index", required_argument, 0, 'i'},
	{"read",  required_argument, 0, 'r'},
	{"read2", required_argument, 0, 'R'},
	{"indel", required_argument, 0, 'I'},
	{"mindist", required_argument, 0, 'd'},
	{"maxdist", required_argument, 0, 'D'},
	{"subreads", required_argument, 0, 'n'},
	{"minmatch", required_argument, 0, 'm'},
	{"minmatch2", required_argument, 0, 'p'},
	{"quality", required_argument, 0, 'Q'},
	{"threads", required_argument, 0, 'T'},
	{"index-threshold", required_argument, 0, 'f'},
	{"output", required_argument, 0, 'o'},
	{"order",  required_argument,0, 'S'},
	{0, 0, 0, 0}
};


int main_align(int argc,char ** argv)
{
	char read_file [300], read2_file [300];
	char output_file [300];

	char index_prefix [300];
	unsigned int all_reads, all_tables;
	gene_allvote_t allvote;
	unsigned long long int processed_reads = 0, succeed_reads = 0;
	gene_input_t ginp, ginp2;
	//gene_flat_t my_flat ;
	//create_flat_strip(&my_flat);

	int c;
	int option_index = 0;

	TOTAL_SUBREADS = 10;
	ACCEPT_SUBREADS = 3;
	ACCEPT_MINOR_SUBREADS = 1;
	INDEX_THRESHOLD = 12;
	read_file[0]=0;
	read2_file[0]=0;
	index_prefix[0]=0;
	output_file[0]=0;
	all_reads = 14*1024*1024;

	printf("\n");



	while ((c = getopt_long (argc, argv, "aSd:D:n:m:p:f:R:r:i:o:T:Q:I:?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'S':
				FIRST_READ_REVERSE = optarg[0]=='r'?1:0;
				SECOND_READ_REVERSE = optarg[1]=='f'?0:1;
				break;
			case 'a':
				USE_VALUE_ARRAY_INDEX = 1;
				break;
			case 'D':
				MAX_PAIRED_DISTANCE = atoi(optarg);
				break;
			case 'd':
				MIN_PAIRED_DISTANCE = atoi(optarg);
				break;
			case 'n':
				TOTAL_SUBREADS = atoi(optarg);
				break;
			case 'f':
				INDEX_THRESHOLD  = atoi(optarg);
				break;
			case 'm':
				ACCEPT_SUBREADS = atoi(optarg);
				break;
			case 'T':
				ALL_THREADS = atoi(optarg);
				break;
			case 'r':
				strncpy(read_file, optarg, 299);
				break;
			case 'R':
				strncpy(read2_file, optarg, 299);
				break;
			case 'i':
				strncpy(index_prefix, optarg, 299);
				break;
			case 'o':
				strncpy(output_file, optarg, 299);
				break;
			case 'I':
				INDEL_TOLERANCE = atoi(optarg);
				if( INDEL_TOLERANCE >5)INDEL_TOLERANCE=5;
				break ;
			case 'Q':
				if(optarg[0]=='l')
					QUALITY_SCALE = QUALITY_SCALE_LINEAR;
				if(optarg[0]=='e')
					QUALITY_SCALE = QUALITY_SCALE_LOG;
				if(optarg[0]=='n')
					QUALITY_SCALE = QUALITY_SCALE_NONE;
				break;
			case 'p':
				ACCEPT_MINOR_SUBREADS = atoi(optarg);
				break;
			case '?':
				return -1 ;
		}


	if (!read_file[0] || !index_prefix[0] || !output_file[0])
	{
		usage(argv[0]);

		return -1 ;
	}
  
	reads_density = guess_reads_density(read_file);
	if(reads_density<0)
		printf("Input file '%s' is not found or is in an incorrect format.\n", read_file);

	if(geinput_open(read_file, &ginp))
		return -1;

	for(all_tables=0; ; all_tables++)
	{
		char fn[300];
		struct stat fstatbuf;
		sprintf(fn, "%s.%02d.%c.tab", index_prefix, all_tables, ginp.space_type==GENE_SPACE_COLOR?'c':'b');
		if(stat(fn, &fstatbuf))
			break;
	}

	if(all_tables==0)
	{
		printf("Unable to open the index files in the %s space.\n",  ginp.space_type==GENE_SPACE_COLOR?"color":"base");
		return -1;
	}

	FILE * out_fp = fopen(output_file, "w");
	if (!out_fp)
	{
		printf("Unable to open the output file at '%s'.\n", output_file);
		return -1;
	}


	printf("Number of subreads selected for each read=%d\n", TOTAL_SUBREADS);
	printf("Threshold on number of subreads for a successful mapping=%d\n", ACCEPT_SUBREADS);
	printf("Number of threads=%d\n", ALL_THREADS);
	printf("Tolerance for Indel=%d\n", INDEL_TOLERANCE);
	if (QUALITY_SCALE==QUALITY_SCALE_LINEAR)
		puts("Quality scale=linear\n\n");
	else if (QUALITY_SCALE==QUALITY_SCALE_LOG)
		puts("Quality scale=exponential\n\n");
	else 	puts("\n");

	if (read2_file[0])
	{
		if (MAX_PAIRED_DISTANCE <= MIN_PAIRED_DISTANCE)
		{
			printf ("The value of the '-D' option must be greater than that of the '-d' option. \n");
			return -1;
		}

		printf ("Performing paired-end alignment:\n");
		printf ("Maximum distance between reads=%d\n", MAX_PAIRED_DISTANCE);
		printf ("Minimum distance between reads=%d\n", MIN_PAIRED_DISTANCE);
		printf ("Threshold on number of subreads for a successful mapping (the minor end in the pair)=%d\n", ACCEPT_MINOR_SUBREADS);
		printf ("The directions of the two input files are: %s, %s\n\n", FIRST_READ_REVERSE?"reversed":"forward", SECOND_READ_REVERSE?"reversed":"forward");
	}

#ifdef REPORT_ALL_THE_BEST
	printf("***** WARNING: the REPORT_ALL_THE_BEST switch is turned on. You need an extra 1 GBytes of RAM space for saving the temporary results. *****\n");
#endif
	
	init_allvote(&allvote, all_reads);
	if(read2_file[0])
		all_reads/=2;
		
	if(read2_file[0] && geinput_open(read2_file, &ginp2))
	{
		printf("Input file '%s' is not found or is in an incorrect format.\n", read2_file);
		return -1;
	}

	fflush(stdout);

	begin_ftime = miltime();

	while (1)
	{
		char inbuff[1201];

		processed_reads += run_search_index(&ginp, read2_file[0] ? (&ginp2):NULL, index_prefix, &allvote, out_fp, processed_reads, all_tables, &succeed_reads);
		clear_allvote(&allvote);

		// test if there no anyreads remaining
		unsigned long long int current_fp = ftello(ginp.input_fp);
		int rl = geinput_next_read(&ginp, NULL, inbuff, NULL);
		if (rl<0)
			break;
		fseeko(ginp.input_fp, current_fp, SEEK_SET);
	}
	fclose(out_fp);
	geinput_close(&ginp);
	if(IS_DEBUG)
		printf("@LOG THE END. \n");
	else
		printf("\n\n %llu reads were processed in %.1f seconds.\nPercentage of successfully mapped reads is %0.2f%%.\n\n", processed_reads, miltime()-begin_ftime, succeed_reads*100.0/processed_reads/(read2_file[0]?2:1));

	printf("\n\nCompleted successfully.\n");

	return 0;
}


