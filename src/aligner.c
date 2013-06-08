
/***************************************************************

   The Subread and Rsubread software packages are free
   software packages:
 
   you can redistribute it and/or modify it under the terms
   of the GNU General Public License as published by the 
   Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Subread is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
   See the GNU General Public License for more details.

   Authors: Drs Yang Liao and Wei Shi

  ***************************************************************/

#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include "gene-algorithms.h"
#include "gene-value-index.h"
#include "input-files.h"
#include "sorted-hashtable.h"

int ALL_THREADS=1;
int TOTAL_SUBREADS;
int ACCEPT_SUBREADS;
int ACCEPT_MINOR_SUBREADS;
int INDEX_THRESHOLD;
int MAX_PAIRED_DISTANCE = 600;
int MIN_PAIRED_DISTANCE = 50;
int REPORT_ONLY_UNIQUE = 0;
int NUMBER_OF_ANCHORS_PAIRED = -1;
int INDEL_TOLERANCE = 6;
int QUALITY_SCALE = QUALITY_SCALE_NONE;
int USE_VALUE_ARRAY_INDEX = 1;
int USE_BASEINDEX_BREAK_TIE = 0;
int MAX_METHYLATION_C_NUMBER = 0;
int REPORT_POTENTIAL_JUNCTION_READS=0;
int FIRST_READ_REVERSE = 0;
int SECOND_READ_REVERSE = 1;
int APPLY_REPEATING_PENALTY = 1;

extern int DPALIGN_CREATEGAP_PENALTY;
extern int DPALIGN_EXTENDGAP_PENALTY;
extern int DPALIGN_MISMATCH_PENALTY;
extern int DPALIGN_MATCH_SCORE;
extern int MAX_CIGAR_LEN;

int FASTQ_FORMAT = FASTQ_PHRED33;
double reads_density;


struct gene_thread_data
{
	gehash_t table;
	int * offsets;
	int number_t;
};



void print_res(gene_value_index_t *array_index , gene_allvote_t *av, gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads)
{
	int i, j, ic=0;
	char inb [1201], nameb[1201], qualityb[1201];

	gene_offset_t offsets;
	unsigned int paired_match = 0;

	load_offsets (&offsets, index_prefix);
	//printf("NS=%d\n", offsets.total_offsets);

	if(all_processed_reads ==0)
	{
		unsigned int last_offset = 0;
		i=0;
		while(offsets.read_offsets[i])
		{
			fprintf(out_fp, "@SQ\tSN:%s\tLN:%u\n", offsets.read_names+i*MAX_READ_NAME_LEN, offsets.read_offsets[i] - last_offset+16);
			last_offset = offsets.read_offsets[i];
			i++;
		}
	}
	
	i=0;

	int rl1=0, rl2=0, rl=0;
	while(1)
	{
		int best_read_id, total_best_read_id;
		nameb[0]=0;
		long long int fpos1 = 0 , fpos2 = 0;
		gene_input_t * Curr_ginp=ginp2;

		for(total_best_read_id = 0 ; total_best_read_id <  av->multi_best_reads; total_best_read_id ++)
		{
			int Curr_position = (ginp2?(i*2):i);
			voting_result_t * Curr_result = &(av -> results[Curr_position * av->multi_best_reads + total_best_read_id]);

			if(total_best_read_id && (REPORT_ONLY_UNIQUE || (ginp2 && !(Curr_result -> masks & IS_PAIRED_MATCH)))) break;
			if(Curr_result -> vote_number <1) break;
		}

		for(best_read_id = 0 ; best_read_id <  av->multi_best_reads; best_read_id += (ginp2?(Curr_ginp==ginp2):1))
		{
			int Curr_position, isOK=0;
			if(ginp2)
				Curr_position = i*2+(Curr_ginp==ginp);
			else
				Curr_position = i;


			voting_result_t * Curr_result = &(av -> results[Curr_position * av->multi_best_reads + best_read_id]);
			if(best_read_id && (Curr_ginp==ginp2 || !ginp2) && (REPORT_ONLY_UNIQUE || Curr_result -> vote_number <1 || (ginp2 && !(Curr_result -> masks & IS_PAIRED_MATCH)))) break;
	
			voting_result_t * Curr_result0 = &(av -> results[Curr_position * av->multi_best_reads]);

			if(ginp2) 
			{
				Curr_ginp = (Curr_ginp==ginp2)?ginp:ginp2;
				if((!(Curr_result0 -> masks & IS_R1R2_EQUAL_LEN)) && (Curr_ginp == ginp) )
				{
					if((!feof(ginp->input_fp)) && !feof(ginp2->input_fp))
					{
						long long int tfpos = ftello(ginp->input_fp);
						rl1 = geinput_next_read(ginp, NULL, inb, NULL);
						fseeko(ginp->input_fp, tfpos, SEEK_SET);

						tfpos = ftello(ginp2->input_fp);
						rl2 = geinput_next_read(ginp2, NULL, inb, NULL);
						fseeko(ginp2->input_fp, tfpos, SEEK_SET);
					}
				}
			}
			else
				Curr_ginp = ginp;

			if(best_read_id) 
			{
				if(!feof(Curr_ginp->input_fp))
					fseeko(Curr_ginp->input_fp, Curr_ginp==ginp2?fpos2:fpos1, SEEK_SET);
			}
			else{
				if(Curr_ginp==ginp2)
					fpos2 = ftello(Curr_ginp->input_fp);
				else
					fpos1 = ftello(Curr_ginp->input_fp);
			}


			//printf("R %d V0=%d\n", Curr_position, Curr_result->vote_number);
			rl = geinput_next_read(Curr_ginp, nameb, inb, qualityb);
			if(Curr_result->masks & IS_R1R2_EQUAL_LEN)
			{
				rl1 = rl2 = rl;
			}

			if (rl<0){
				break;
			}

			if (nameb[0]==0)
				sprintf(nameb, "read_%llu", i+1+all_processed_reads);

			int flags = 0;
			unsigned char votes = Curr_result -> vote_number;

			int mate_index = 0, mate_rl = 0; 
			int is_mate_ok = 0;
			unsigned int mate_pos=0;

			int mate_offset_pos = 0;

			char * rnext_pnt, mate_cigar[100], rnext[1201];	// the name of the chromosome where the other read is on; "=" means the same
			unsigned int pnext; 		// the offset of the other read on the above chromosome

			voting_result_t * Mate_result = NULL; 
			if(ginp2)
			{
				flags |= SAM_FLAG_PAIRED_TASK;
				mate_index = i*2+(Curr_ginp!=ginp2);
				Mate_result = av->results + mate_index*av->multi_best_reads + best_read_id;
				is_mate_ok = ((Mate_result -> vote_number>0 && Curr_result -> vote_number>0 && (Curr_result->masks & IS_PAIRED_MATCH)) ||Mate_result->vote_number  >= ACCEPT_SUBREADS) &&  ((!REPORT_ONLY_UNIQUE) ||  !(Curr_result->masks & IS_BREAKEVEN_READ));


				if (Curr_ginp==ginp2)
					flags |= SAM_FLAG_SECOND_READ_IN_PAIR;
				else
					flags |= SAM_FLAG_FIRST_READ_IN_PAIR;

				if(is_mate_ok)
				{
					int mate_rl_adj=0;
					int is_safeguarded= 0;

					mate_rl = (Curr_ginp == ginp2)?rl2:rl1;

					long long int mate_rec_offset =  mate_index*av->multi_best_reads + best_read_id;
					mate_rec_offset *= av -> indel_recorder_length;

					show_cigar(av-> all_indel_recorder+ mate_rec_offset , mate_rl, 0, mate_cigar, INDEL_TOLERANCE, TOTAL_SUBREADS, NULL,&mate_offset_pos, &mate_rl_adj, &is_safeguarded);

					mate_pos = Mate_result -> read_pos;
					mate_pos += mate_offset_pos;

					if(locate_gene_position_max(mate_pos, &offsets, &rnext_pnt, &pnext, mate_rl+mate_rl_adj - mate_offset_pos)) is_mate_ok = 0;
				}

				if(!is_mate_ok)
					flags|=SAM_FLAG_MATE_UNMATCHED;

			}

			if (((votes >= ACCEPT_SUBREADS) || ((Curr_result->masks & IS_PAIRED_MATCH) && votes>0 && Mate_result -> vote_number>0) ) && ((!REPORT_ONLY_UNIQUE) ||  !(Curr_result -> masks & IS_BREAKEVEN_READ)))
			{

				if(Curr_result->masks & IS_PAIRED_MATCH)
					paired_match++;
				int lpos = Curr_result->read_pos;
				char is_counterpart = Curr_result->is_negative_strand; 

				char cigar_str [100];
				int offset_pos = 0;

				char * read_name;
				unsigned int read_pos;
				int rl_adj = 0;
				int is_safeguarded=0;


				cigar_str[0]=0;

				if(is_counterpart + (Curr_ginp==ginp2) == 1)
				{
					int rev_offset = (Curr_ginp->space_type==GENE_SPACE_COLOR && inb[0]>='A' && inb[0]<='Z')?1:0;
					reverse_read(inb+ rev_offset, rl- rev_offset, ginp->space_type);
					if(qualityb[0])
						reverse_quality(qualityb, rl- rev_offset);
				}

				long long int curr_rec_offset =  Curr_position*av->multi_best_reads + best_read_id;
				curr_rec_offset *= av -> indel_recorder_length;
				show_cigar(av-> all_indel_recorder+ curr_rec_offset, rl, is_counterpart, cigar_str, INDEL_TOLERANCE, TOTAL_SUBREADS, inb,&offset_pos, &rl_adj, &is_safeguarded);

				lpos += offset_pos;

				if(!locate_gene_position_max(lpos, &offsets, &read_name, &read_pos, rl + rl_adj - offset_pos))
				{

					if (best_read_id) flags |= SAM_FLAG_SECONDARY_ALIGNMENT;

					if (Curr_result -> masks & IS_PAIRED_MATCH)
						flags |= SAM_FLAG_MATCHED_IN_PAIR;

					if (is_counterpart + (Curr_ginp==ginp2) == 1)
						flags |= SAM_FLAG_REVERSE_STRAND_MATCHED; 

					if(ginp2)
					{
						if(is_mate_ok && Mate_result -> is_negative_strand + (Curr_ginp!=ginp2) == 1)
							flags |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
						remove_backslash(nameb);
					}

					if (ginp2 && is_mate_ok)
					{



						long long int tlen; 		// the total length of the whole segment, which is split into two readso.

						if (Curr_ginp==ginp2)
							tlen = ( (mate_pos > lpos) ? (mate_pos - mate_offset_pos - lpos + offset_pos + rl1) :(lpos - offset_pos - mate_pos + mate_offset_pos + rl2));
						else
							tlen = ( (mate_pos > lpos) ? (mate_pos - mate_offset_pos - lpos + offset_pos + rl2) :(lpos - offset_pos - mate_pos + mate_offset_pos + rl1));
						if (mate_pos < lpos) tlen=-tlen;

						if(strncmp(rnext_pnt, read_name, 1200) ==0)
							strcpy(rnext,"=");
						else
						{
							strcpy(rnext, rnext_pnt);
							flags &= (0xffffffff ^ SAM_FLAG_MATCHED_IN_PAIR);
							tlen = 0;
						}

						gene_quality_score_t mapping_quality = Curr_result -> final_quality;

						fprintf (out_fp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%lld\t%s\t",  nameb, flags, read_name,read_pos+1 , (int)mapping_quality , cigar_str, rnext, pnext+1, tlen, inb);

						if(qualityb[0])
						{
							if (FASTQ_FORMAT == FASTQ_PHRED64)
								fastq_64_to_33(qualityb);
							fputs(qualityb,out_fp);
						}
						else
							for(j=0; j<rl; j++)
								fputc('J', out_fp);
						fprintf(out_fp, "\tAS:i:%d\tNM:i:%d\tNH:i:%d", Curr_result -> vote_number, Curr_result->edit_distance, total_best_read_id);

					}
					else
					{
						gene_quality_score_t mapping_quality = Curr_result -> final_quality;

						fprintf (out_fp, "%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t",  nameb, flags, read_name, read_pos+1, (int)mapping_quality, cigar_str, inb);

						if(qualityb[0])
						{
							if (FASTQ_FORMAT == FASTQ_PHRED64)
								fastq_64_to_33(qualityb);
							fputs(qualityb, out_fp);
						}
						else
							for(j=0; j<rl; j++)
								fputc('J', out_fp);
						fprintf(out_fp, "\tAS:i:%d\tNM:i:%d\tNH:i:%d", Curr_result -> vote_number, Curr_result->edit_distance,total_best_read_id);
					}
					fputc('\n',out_fp);

					isOK = 1;
					(*succeed_reads) ++;
				}
			}

			if(!isOK)
			{

				char * rnext_pnt;	// the name of the chromosome where the other read is on; "=" means the same
				unsigned int pnext=0; 		// the offset of the other read on the above chromosome

				if(is_mate_ok)
					locate_gene_position(mate_pos, &offsets, &rnext_pnt, &pnext);
				else	rnext_pnt = "*";
					
				flags |= SAM_FLAG_UNMAPPED;

				if(ginp2)
					remove_backslash(nameb);

				fprintf (out_fp, "%s\t%d\t*\t0\t0\t*\t%s\t%u\t0\t%s\t", nameb, flags,rnext_pnt  , pnext+(is_mate_ok?1:0), inb);
				if(qualityb[0])
				{

					if (FASTQ_FORMAT == FASTQ_PHRED64)
						fastq_64_to_33(qualityb);
					fputs(qualityb, out_fp);
				}
				else
					for(j=0; j<rl; j++)
						fputc('J', out_fp);
				fputs("\tNH:i:0", out_fp);
				fprintf(out_fp, "\n");
			}
	//			SUBREADprintf ("@UNK Q#%d UNKNOW %s\n", i, inb);

		}

		if (rl<0){
			break;
		}

		//if((!ginp2) || ginp2 == Curr_ginp)
		i++;
		if(i >= processed_reads)break;


		if(i % 10000 ==0 && i>1)
			print_text_scrolling_bar("Saving results", i*1./processed_reads, 80, &ic);
	}
	if (paired_match)
		SUBREADprintf("\nReads referencing paired-end info: %u\n", paired_match);
	//else
	//	SUBREADprintf("\nNo reads are mapped in pairs.\b");

	destroy_offsets (&offsets);
}



struct gene_thread_data_transport
{
	gehash_t * my_table;
	gene_value_index_t * my_value_array_index;
	gene_value_index_t * array_index_set;
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
	int run_method;

	pthread_spinlock_t * input_data_lock;
	pthread_spinlock_t * init_lock;
};

void run_final_stage(gene_value_index_t * array_index_set ,  gene_allvote_t * allvote,   gene_input_t * ginp,gene_input_t * ginp2 , unsigned int * processed_reads, int all_tables, pthread_spinlock_t * input_lock, int my_thread_no, int segment_length)
{
	char read1_txt_buf[1201], qual1_txt_buf[1201];
	char read2_txt_buf[1201], qual2_txt_buf[1201], namebuf[200];
	int queries = 0, rl1 = 0, rl2 = 0, i, ic=0;
	int all_reads = allvote -> max_len;
	int array_index_set_len = 0;

	short ** dynamic_align_short;
	char ** dynamic_align_char ;

	dynamic_align_short = (short **) malloc(sizeof(short*) * MAX_READ_LENGTH);
	dynamic_align_char  = (char  **) malloc(sizeof(char *) * MAX_READ_LENGTH);

	for(i=0; i<MAX_READ_LENGTH; i++)
	{
		dynamic_align_short[i] = (short *)malloc(sizeof(short) * MAX_READ_LENGTH);
		dynamic_align_char [i] = (char * )malloc(sizeof(char ) * MAX_READ_LENGTH);
	}


	if (ginp2)
		all_reads /=2;

	qual1_txt_buf[0]=0;
	qual2_txt_buf[0]=0;
	
	gene_value_index_t * tmp_index;
	for(tmp_index = array_index_set; tmp_index -> length; tmp_index ++)
		array_index_set_len++;


	while(1)
	{
		char * read1_txt = read1_txt_buf;
		char * read2_txt = read2_txt_buf;
		char * qual1_txt = qual1_txt_buf;
		char * qual2_txt = qual2_txt_buf;
		
		if(input_lock!=NULL)
			pthread_spin_lock(input_lock);

		if(*processed_reads >= all_reads)
		{
			if(input_lock!=NULL)
				pthread_spin_unlock(input_lock);
			break;
		}



		rl1 = geinput_next_read(ginp, namebuf, read1_txt, qual1_txt);
		if(ginp2)
			rl2 = geinput_next_read(ginp2, namebuf, read2_txt, qual2_txt);

		queries = *processed_reads;
		if(rl1>0)
			(*processed_reads ) ++;

		if(input_lock!=NULL)
			pthread_spin_unlock(input_lock);

		if (rl1<0)
			break;

		if (ginp->space_type == GENE_SPACE_COLOR && read1_txt[0]>='A' && read1_txt[0]<='Z')
		{
			read1_txt ++; 
			if(qual1_txt[0])
				qual1_txt++;
			rl1 --;

			if(ginp2){
				read2_txt ++;
				if(qual2_txt[0])
					qual2_txt++;
				rl2 --;
			}
	
		}

		if(my_thread_no<1 && queries % 10000 ==0 && queries>1)
			print_text_scrolling_bar("Finalising", queries*1./segment_length, 80, &ic);

		int best_read_id;
		int r1_reversed = 0;
		int r2_reversed = 0;
		for(best_read_id = 0; best_read_id < allvote->multi_best_reads; best_read_id++)
		{
			if(best_read_id>0)
			{
				if(ginp2 && allvote->results[queries*2*allvote->multi_best_reads+best_read_id].vote_number==0) break;
				if(!ginp2 && allvote->results[queries*allvote->multi_best_reads+best_read_id].vote_number==0) break;
			}
			if(ginp2)
			{
				int r1_need_reverse = (FIRST_READ_REVERSE  + allvote->results[queries*2*allvote->multi_best_reads+best_read_id].is_negative_strand)==1;
				int r2_need_reverse = (SECOND_READ_REVERSE + allvote->results[(queries*2+1)*allvote->multi_best_reads+best_read_id].is_negative_strand)==1;
				int numvote_read1 = allvote->results[queries*2*allvote->multi_best_reads+best_read_id].vote_number;
				int numvote_read2 = allvote->results[(queries*2+1)*allvote->multi_best_reads+best_read_id].vote_number;
				if((r1_need_reverse ^ r1_reversed ) && numvote_read1)
				{
					reverse_read(read1_txt, rl1, ginp->space_type);
					reverse_quality(qual1_txt, rl1);
					r1_reversed = !r1_reversed;
				}

				if((r2_need_reverse ^ r2_reversed)&& numvote_read2)
				{
					reverse_read(read2_txt, rl2, ginp->space_type);
					reverse_quality(qual2_txt, rl2);
					r2_reversed = !r2_reversed;
				}
		
				int read_second = 0;

				for(read_second = 0; read_second <2; read_second ++)
				{
					gene_value_index_t * my_value_array_index = NULL;
					int array_index_no;
					int numvote_read = read_second?numvote_read2:numvote_read1;

					if(numvote_read <1) continue;

					unsigned int pos_read = allvote->results[(queries*2 + read_second)  * allvote -> multi_best_reads + best_read_id].read_pos;
					char * qual_txt = read_second?qual2_txt:qual1_txt;
					char * read_txt = read_second?read2_txt:read1_txt;
					int r_reverse = read_second?r2_need_reverse:r1_need_reverse;
					int read_mask = allvote->results[(queries*2 + read_second)  * allvote -> multi_best_reads + best_read_id].masks;
					int rl = read_second?rl2:rl1;

					char indel_recorder[MAX_INDEL_TOLERANCE*3];


					for(array_index_no = 0; array_index_no < array_index_set_len ; array_index_no ++) 
					{
						my_value_array_index = array_index_set + array_index_no;
						if((pos_read > my_value_array_index -> start_point + (array_index_no?10000:0))
							&& (pos_read < my_value_array_index -> start_point + my_value_array_index -> length - ((array_index_no == array_index_set_len-1)?0:10000)))
								break;
					}

					int mismatch =0;
					int is_safeguarded=0;

					long long int qur_rec_offset =(queries*2 + read_second)  * allvote -> multi_best_reads + best_read_id;
					qur_rec_offset*= allvote -> indel_recorder_length;
					indel_recorder_copy(indel_recorder, allvote->all_indel_recorder +  qur_rec_offset);

					find_and_explain_indel(allvote, (queries*2 + read_second) * allvote -> multi_best_reads + best_read_id  ,  pos_read , numvote_read , 0, r_reverse , read_mask|IS_PAIRED_MATCH, indel_recorder, my_value_array_index, read_txt , rl , INDEL_TOLERANCE, TOTAL_SUBREADS, ginp->space_type, REPORT_POTENTIAL_JUNCTION_READS, ! r_reverse, qual_txt , FASTQ_FORMAT, dynamic_align_short, dynamic_align_char);

					char cigar_str [100];
					cigar_str[0]=0;

					show_cigar( allvote->all_indel_recorder + qur_rec_offset, rl, r_reverse , cigar_str, INDEL_TOLERANCE, TOTAL_SUBREADS, read_txt, NULL, NULL, &is_safeguarded);
					if(is_safeguarded)
						allvote->results[(2*queries+read_second) * allvote -> multi_best_reads + best_read_id].vote_number = 0;
					else
					{
						unsigned int tmp_mapped_pos = allvote->results[(2*queries+read_second) * allvote -> multi_best_reads + best_read_id].read_pos;

						allvote->results[(2*queries+read_second) * allvote -> multi_best_reads + best_read_id].final_quality = final_mapping_quality(my_value_array_index, tmp_mapped_pos, read_txt, qual_txt, cigar_str, FASTQ_FORMAT, &mismatch, rl, APPLY_REPEATING_PENALTY);
						//allvote->results[(2*queries+read_second) * allvote -> multi_best_reads + best_read_id].final_quality = final_mapping_quality_edit(my_value_array_index, tmp_mapped_pos, read_txt, qual_txt, cigar_str, FASTQ_FORMAT, &mismatch, rl);
						//allvote->results[(2*queries+read_second) * allvote -> multi_best_reads + best_read_id].edit_distance = mismatch;
					}
				}
			}
			else
			{
				char tmp_indel_recorder[MAX_INDEL_TOLERANCE*3];
				int numvote_read = allvote->results[queries* allvote -> multi_best_reads + best_read_id].vote_number;
				int r_need_reverse = allvote->results[queries* allvote -> multi_best_reads + best_read_id].is_negative_strand;
				int array_index_no;
				gene_value_index_t  * my_value_array_index = NULL;


				if(numvote_read>0)
				{
					if(r_need_reverse ^ r1_reversed)
					{
						reverse_read(read1_txt, rl1, ginp->space_type);
						reverse_quality(qual1_txt, rl1);
						r1_reversed=!r1_reversed;
					}


					unsigned int pos_read = allvote->results[queries * allvote -> multi_best_reads + best_read_id].read_pos;


					for(array_index_no = 0; array_index_no < array_index_set_len ; array_index_no ++) 
					{
						my_value_array_index = array_index_set + array_index_no;
						if((pos_read > my_value_array_index -> start_point + (array_index_no?10000:0))
							&& (pos_read < my_value_array_index -> start_point + my_value_array_index -> length - ((array_index_no == array_index_set_len-1)?0:10000)))
								break;
					}


					int mismatch =0;
					short read_mask = allvote->results[queries * allvote -> multi_best_reads + best_read_id].masks;
					int is_safeguarded = 0;



					long long int qur_rec_offset =queries * allvote -> multi_best_reads + best_read_id;
					qur_rec_offset*= allvote -> indel_recorder_length;
					indel_recorder_copy(tmp_indel_recorder, allvote->all_indel_recorder + qur_rec_offset);

					find_and_explain_indel(allvote, queries* allvote -> multi_best_reads + best_read_id , pos_read , numvote_read , 0, r_need_reverse , read_mask|IS_PAIRED_MATCH, tmp_indel_recorder, my_value_array_index, read1_txt , rl1 , INDEL_TOLERANCE, TOTAL_SUBREADS, ginp->space_type, REPORT_POTENTIAL_JUNCTION_READS, ! r_need_reverse, qual1_txt , FASTQ_FORMAT, dynamic_align_short, dynamic_align_char);

					char cigar_str [100];
					cigar_str[0]=0;

					show_cigar(allvote->all_indel_recorder + qur_rec_offset, rl1, r_need_reverse , cigar_str, INDEL_TOLERANCE, TOTAL_SUBREADS, read1_txt, NULL, NULL, &is_safeguarded);
					unsigned int tmp_mapped_pos = allvote->results[queries * allvote -> multi_best_reads + best_read_id].read_pos;
					
					if(is_safeguarded)
					{
						//printf("SG! %s\n", cigar_str);
						allvote ->results[queries * allvote -> multi_best_reads + best_read_id].vote_number = 0;
					}
					else
					{
						int repeated_number = 0;

						if(APPLY_REPEATING_PENALTY)
							repeated_number = test_big_margin(allvote, queries);

						allvote->results[queries * allvote -> multi_best_reads + best_read_id].final_quality = final_mapping_quality(my_value_array_index, tmp_mapped_pos, read1_txt, qual1_txt, cigar_str, FASTQ_FORMAT, &mismatch, rl1, APPLY_REPEATING_PENALTY) / (1+repeated_number);
						//allvote->results[queries * allvote -> multi_best_reads + best_read_id].final_quality = final_mapping_quality_edit(my_value_array_index, tmp_mapped_pos, read1_txt, qual1_txt, cigar_str, FASTQ_FORMAT, &mismatch, rl1) / (1+repeated_number);
						allvote->results[queries * allvote -> multi_best_reads + best_read_id].edit_distance = mismatch;
					}
				}

			}

		}
	}
	for(i=0; i< MAX_READ_LENGTH; i++)
	{
		free(dynamic_align_short[i]);
		free(dynamic_align_char [i]);
	}
	free(dynamic_align_short);
	free(dynamic_align_char);


}

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
	char namebuf[200];

	gene_vote_t * vote_memblock_1, *vote_memblock_2;
	//short ** dynamic_align_short =NULL;
	//char ** dynamic_align_char =NULL;
	
	double local_begin_ftime = miltime();
	int read_len = 0, read2_len = 0;

	is_counterpart = 0;
	if (ginp2)
		all_reads /=2;

//	SUBREADprintf ("I'm the %d-th thread\n", my_thread_no);

	unsigned int low_border = my_value_array_index -> start_base_offset;
	unsigned int high_border = my_value_array_index -> start_base_offset + my_value_array_index -> length; 

	vote_memblock_1 = (gene_vote_t*) malloc(sizeof(gene_vote_t));
	vote_memblock_2 = (gene_vote_t*) malloc(sizeof(gene_vote_t));

	if(!vote_memblock_1 || !vote_memblock_2)
	{
		SUBREADprintf("Cannot allocate memory for voting\n");
		return -1;
	}	


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
				InBuff = BuffMemory;
				QualityBuff = BuffMemory2;
				read_len = geinput_next_read(ginp, namebuf, InBuff, QualityBuff);
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
			gene_vote_t *vote_read1 = vote_memblock_1;
			gene_vote_t *vote_read2 = vote_memblock_2;

			init_gene_vote(vote_read1);
			init_gene_vote(vote_read2);

			int is_paired ;
			short current_masks = 0;

			if(read2_len == read_len){
				current_masks = IS_R1R2_EQUAL_LEN;
			}

			for (is_paired = 0; is_paired <2; is_paired ++)
			{
				int subread_no;

				char * CurrInBuff = is_paired?InBuff2:InBuff;
				int Curr_read_len = is_paired?read2_len:read_len;
				gene_vote_t * vote = is_paired?vote_read2:vote_read1;
				char * CurrQualityBuff = is_paired?QualityBuff2:QualityBuff;

				for(subread_no=0; ; subread_no++)
				{
					int subread_offset1 = (int)(subread_step * (subread_no+1));
					subread_offset1 -= subread_offset1%GENE_SLIDING_STEP;
					subread_offset1 += GENE_SLIDING_STEP-1;

					for(i=0; i<GENE_SLIDING_STEP ; i++)
					{
						int subread_offset = (int)(subread_step * subread_no); 
						subread_offset -= subread_offset%GENE_SLIDING_STEP -i;


						char * subread_string = CurrInBuff + subread_offset;
						gehash_key_t subread_integer = genekey2int(subread_string, ginp->space_type);

						gene_quality_score_t subread_quality = 1;

						if(QUALITY_SCALE != QUALITY_SCALE_NONE)
						{
							subread_quality = 1.;
							char * quality_string = CurrQualityBuff + subread_offset;
							subread_quality = get_subread_quality(quality_string, subread_string, QUALITY_SCALE, FASTQ_FORMAT);
						}


						if(is_valid_subread(subread_string))
						{
							if(MAX_METHYLATION_C_NUMBER)
								gehash_go_q_CtoT(my_table, subread_integer , subread_offset, read_len, is_counterpart , vote , 1 , 1, subread_quality, INDEX_THRESHOLD, INDEL_TOLERANCE, subread_no, MAX_METHYLATION_C_NUMBER, low_border, high_border - read_len);
							else
								gehash_go_q(my_table, subread_integer , subread_offset, read_len, is_counterpart , vote , 1 , 1, subread_quality, INDEX_THRESHOLD, INDEL_TOLERANCE, subread_no, low_border, high_border - read_len);
						}

					}
					if(subread_offset1 >= Curr_read_len -16)
					{
						break;
					}
				}
			}

			int is_paired_match = 0;
			finalise_vote(vote_read1);
			finalise_vote(vote_read2);



		/*	
			SUBREADprintf("\n\nV1:======\n");
			print_votes(vote_read1, index_prefix);

			SUBREADprintf("\n\nV2:======\n");
			print_votes(vote_read2, index_prefix);
		*/

			reg_big_margin_votes(vote_read1, allvote, queries*2);
			reg_big_margin_votes(vote_read2, allvote, queries*2+1);

			is_paired_match = select_positions_array(USE_BASEINDEX_BREAK_TIE, queries, allvote, current_masks, InBuff, read_len, InBuff2, read2_len,vote_read1, vote_read2, MAX_PAIRED_DISTANCE, MIN_PAIRED_DISTANCE, ACCEPT_SUBREADS, ACCEPT_MINOR_SUBREADS, is_counterpart, my_value_array_index, ginp -> space_type, INDEL_TOLERANCE, QualityBuff, QualityBuff2, QUALITY_SCALE, TOTAL_SUBREADS, my_thread_no);

			if (!is_paired_match && !(allvote -> results[queries*2*allvote->multi_best_reads + 0].masks & IS_PAIRED_MATCH))
			{
				if(vote_read1->max_vote > 0)
					final_matchingness_scoring(USE_BASEINDEX_BREAK_TIE,InBuff, (QUALITY_SCALE == QUALITY_SCALE_NONE )?NULL:QualityBuff,  is_counterpart, read_len, vote_read1, 0, allvote, queries*2, my_value_array_index, ginp -> space_type, INDEL_TOLERANCE, QUALITY_SCALE, TOTAL_SUBREADS);
				if(vote_read2->max_vote > 0)
					final_matchingness_scoring(USE_BASEINDEX_BREAK_TIE,InBuff2, (QUALITY_SCALE == QUALITY_SCALE_NONE )?NULL:QualityBuff2, is_counterpart, read2_len, vote_read2, 0, allvote, queries*2+1, my_value_array_index, ginp -> space_type, INDEL_TOLERANCE, QUALITY_SCALE, TOTAL_SUBREADS);
			}

		}
		else
		{
			gene_vote_t * vote = vote_memblock_1;
			int subread_no;
			init_gene_vote(vote);
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

					gene_quality_score_t subread_quality = 1;

					if(QUALITY_SCALE != QUALITY_SCALE_NONE )
					{
						char * quality_string = QualityBuff + subread_offset;
						subread_quality = get_subread_quality(quality_string, subread_string, QUALITY_SCALE, FASTQ_FORMAT);
					}

					if(is_valid_subread(subread_string))
					{
						if(MAX_METHYLATION_C_NUMBER)
							total_results = gehash_go_q_CtoT (my_table, subread_integer , subread_offset , read_len, is_counterpart, vote , 1, 1 , subread_quality, INDEX_THRESHOLD, INDEL_TOLERANCE, subread_no, MAX_METHYLATION_C_NUMBER, low_border, high_border - read_len);
						else
							total_results = gehash_go_q(my_table, subread_integer , subread_offset , read_len, is_counterpart, vote , 1, 1 , subread_quality, INDEX_THRESHOLD, INDEL_TOLERANCE, subread_no, low_border, high_border - read_len);
					}

				}
				if(subread_offset1 >= read_len -16)
					break;

			}



			//if(vote.max_vote > max(0,subread_no - (TOTAL_SUBREADS - ACCEPT_SUBREADS)))
			//	print_votes(vote, index_prefix);

			if(vote->max_vote >= ACCEPT_SUBREADS)
			{
				finalise_vote(vote);

				reg_big_margin_votes(vote, allvote, queries);
				final_matchingness_scoring(USE_BASEINDEX_BREAK_TIE,InBuff, (QUALITY_SCALE == QUALITY_SCALE_NONE )?NULL:QualityBuff, is_counterpart, read_len, vote, 0, allvote, queries, my_value_array_index, ginp -> space_type, INDEL_TOLERANCE, QUALITY_SCALE, TOTAL_SUBREADS);
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
				SUBREADprintf("@LOG Done %d/%d, good %d, last time %f\n",queries, table_no, good_match, miltime() - t0);
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
				print_running_log(finished_rate, reads_per_second, expected_seconds, (unsigned long long int)all_steps / all_tables, ginp2!=NULL);
			}
			SUBREADfflush(stdout);
			t0 = miltime();
		}
	
		if (is_counterpart)
		{
			if(queries >= all_reads)
				break;
		}

		is_counterpart = !is_counterpart;
	}
	free(vote_memblock_1);
	free(vote_memblock_2);
	return 0;
}


void * run_search_thread(void * parameters)
{
	struct gene_thread_data_transport * data_param = parameters;
	int thid = data_param->this_thread;

	pthread_spin_unlock(data_param -> init_lock);

	if(data_param -> run_method == RUN_ALIGN)
	{
		run_search(data_param->my_table, data_param->my_value_array_index, data_param->table_no, data_param->all_vote, data_param->ginp , data_param->ginp2, data_param->index_prefix, data_param->processed_reads, data_param->base_number, data_param->all_tables, data_param->input_data_lock, thid, data_param-> section_length);
	}
	else
	{

		run_final_stage(data_param -> array_index_set , data_param->all_vote ,  data_param-> ginp,  data_param-> ginp2 , data_param-> processed_reads, data_param-> all_tables,  data_param-> input_data_lock, thid,  data_param-> section_length);
	}
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
	int i; 

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
			SUBREADprintf ("@LOG Loading table from %s\n", table_fn);
		else
			SUBREADprintf ("Loading the %02d-th index file ...                                              \r", tabno+1);
		SUBREADfflush(stdout);

		if(gehash_load(my_table, table_fn)) return -1;
		if(USE_VALUE_ARRAY_INDEX)
		{
			
			sprintf(table_fn, "%s.%02d.%c.array", index_prefix, tabno, ginp->space_type==GENE_SPACE_COLOR?'c':'b');
			stat_ret = stat(table_fn, &filestat);
			if (stat_ret !=0)
			{
				SUBREADprintf("\nFatal Error: The specified index does not contain any value array data. \nPlease do not specify the '-b' option or rebuild the index with a '-b' option.\n");
				return -1;
			}

			if(tabno>0)gvindex_destory(&value_array_index);
			if(gvindex_load(&value_array_index,table_fn)) return -1;
		}

		processed_reads = 0;

// Run the search algorithm on a part of the index
		if(ALL_THREADS <2)
		{
			run_search(my_table, USE_VALUE_ARRAY_INDEX? &value_array_index : NULL, tabno, all_vote, ginp , ginp2, index_prefix, &processed_reads, base_number, all_tables, NULL /*the data lock is null*/, 0  /*I'm the 0-th thread*/, section_length);
		}
		else
		{
			struct gene_thread_data_transport data_param;
			pthread_t runners [ALL_THREADS];
			pthread_spinlock_t  data_lock;
			pthread_spinlock_t  init_lock;

			data_param.my_table = my_table;
			data_param.my_value_array_index = USE_VALUE_ARRAY_INDEX? &value_array_index:NULL;
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
			data_param.run_method = RUN_ALIGN;
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

		if (section_length < 1)
			section_length = processed_reads;

		gehash_destory_fast(my_table);

		// FINAL STAGE
		if(tabno >= all_tables - 1)
		{
			processed_reads = 0;
			gene_value_index_t value_index_set[100];

			memset(value_index_set,0, sizeof(gene_value_index_t)*100);

			memcpy(value_index_set+all_tables - 1, &value_array_index, sizeof(gene_value_index_t));
			for(i=0; i<all_tables - 1; i++)
			{
				sprintf(table_fn, "%s.%02d.%c.array", index_prefix, i, ginp->space_type==GENE_SPACE_COLOR?'c':'b');
				if(gvindex_load(value_index_set + i,table_fn)) return -1;
			}
			

			fseeko(ginp -> input_fp, current_fp, SEEK_SET);
			if (ginp2)
				fseeko(ginp2 -> input_fp, current_fp2, SEEK_SET);


			if(ALL_THREADS <2)
			{
				run_final_stage(value_index_set , all_vote ,  ginp,  ginp2 , & processed_reads,  all_tables, NULL, 0, section_length);
			}
			else
			{
				struct gene_thread_data_transport data_param;
				pthread_t runners [ALL_THREADS];
				pthread_spinlock_t  data_lock;
				pthread_spinlock_t  init_lock;

				data_param.array_index_set = value_index_set;
				data_param.all_tables = all_tables;
				data_param.all_vote = all_vote;
				data_param.all_threads = ALL_THREADS;
				data_param.input_data_lock = &data_lock;
				data_param.init_lock = &init_lock;
				data_param.processed_reads = &processed_reads;
				data_param.section_length = section_length;
				data_param.ginp = ginp;
				data_param.ginp2 = ginp2;
				data_param.run_method = RUN_FINAL;

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
	

			for(i=0; i<all_tables - 1; i++)
				gvindex_destory(value_index_set + i);

		}
		fseeko(ginp -> input_fp, current_fp, SEEK_SET);
		if (ginp2)
			fseeko(ginp2 -> input_fp, current_fp2, SEEK_SET);

		tabno ++;
	}

	print_res(NULL , all_vote, ginp, ginp2, out_fp, index_prefix, processed_reads, base_number, succeed_reads);

	
	gvindex_destory(&value_array_index);
	return processed_reads;
	
}

void usage(char * execname)
{
	SUBREADputs("Usage:");
	SUBREADputs(" ./subread-align [options] -i <index_name> -r <input> -o <output>");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("    -i --index     <index>\t base name of the index.");
	SUBREADputs("    -r --read      <input>\t name of the input file(FASTQ/FASTA format). For paired-end read data, this will be the first read file and the other read file should be provided via the -R option.");
	SUBREADputs("    -o --output    <output>\t name of the output file(SAM format).");
	SUBREADputs("");
	SUBREADputs("Optional arguments:");
	SUBREADputs("    -n --subreads  <int>\t number of selected subreads, 10 by default.");
	SUBREADputs("    -m --minmatch  <int>\t consensus threshold (minimal number of consensus subreads required) for reporting a hit. If paired-end read data are provided, this gives the consensus threshold for the read which receives more votes than the other read from the same pair. 3 by default.");
	SUBREADputs("    -T --threads   <int>\t number of threads, 1 by default.");
	SUBREADputs("    -I --indel     <int>\t number of indels allowed, 5 by default. Up to 16 indels are allowed.");
	SUBREADputs("    -B --multi     <int>\t Specify the maximal number of equally-best mapping locations allowed to be reported for any read. The value has to be within the range of 1 to 16. 1 by default. ");
	SUBREADputs("    -P --phred     <3:6>\t the format of Phred scores in input files, '3' for phred+33 and '6' for phred+64. '3' by default.");
	SUBREADputs("    -u --unique         \t reporting uniquely mapped reads only.");
	SUBREADputs("    -Q --quality        \t using mapping quality scores to break ties when more than one best mapping locations are found.");
	SUBREADputs("    -H --hamming        \t using Hamming distance to break ties when more than one best mapping locations are found.");
	SUBREADputs("    -J --junction       \t mark those bases which can not be aligned together with other bases from the same read using the `S' operation in the CIGAR string (soft-clipping). This option is useful for marking exon-spanning reads and fusion reads.");
	SUBREADputs("    -v                  \t displaying the version number.");
	SUBREADputs("");
	SUBREADputs("Arguments for paired-end reads:");
	SUBREADputs("    -R --read2     <input>\t name of the second input file. The program will then be switched to the paired-end read mapping mode.");
	SUBREADputs("    -p --minmatch2 <int>\t consensus threshold for the read which receives less votes than the other read from the same pair, 1 by default");
	SUBREADputs("    -d --mindist   <int>\t minimum fragment/template length, 50bp by default.");
	SUBREADputs("    -D --maxdist   <int>\t maximum fragment/template length, 600bp by default.");
	SUBREADputs("    -S --order     <ff:fr:rf> \t orientation of the two read from the same pair, 'fr' by default");
	SUBREADputs("");
	SUBREADputs("Arguments for INDEL detection:");
	//SUBREADputs("    --cigar_len       <int>\t  maximum length of a cigar string");
	SUBREADputs("    -G --cre_gap_penalty <int>\t gap opening penalty, default: -2");
	SUBREADputs("    -E --ext_gap_penalty <int>\t gap extending penalty, default: 0");
	SUBREADputs("    -X --mis_penalty     <int>\t mismatch penalty, default: 0");
	SUBREADputs("    -Y --match_score     <int>\t match score, default: 2");
	SUBREADputs("");
        //SUBREADputs("Arguments for mapping bisulfite reads:");
        //SUBREADputs("    -b --bisulfite <int>\t Number of bases 'c' allowed to be converted to 't' in each extracted subread (16bp mer), default: 1.");
        //SUBREADputs("");
	SUBREADputs("Example:");
	SUBREADputs(" ./subread-align -i my_index -r reads.fastq -o my_result.sam ");
	SUBREADputs("");
	SUBREADputs("For more information about these arguments, please refer to the User Manual.\n");

}

static struct option long_options[] =
{
	{"junction",  no_argument, &REPORT_POTENTIAL_JUNCTION_READS , 1},
	{"unique", no_argument, &REPORT_ONLY_UNIQUE, 1},
	{"hamming", no_argument, &USE_BASEINDEX_BREAK_TIE, 1},
	{"allow-repeating", no_argument, &APPLY_REPEATING_PENALTY, 0},
	{"bisulfite ", required_argument, 0, 'b'},
	{"phred", required_argument, 0, 'P'},
	{"index", required_argument, 0, 'i'},
	{"read",  required_argument, 0, 'r'},
	{"read2", required_argument, 0, 'R'},
	{"indel", required_argument, 0, 'I'},
	{"mindist", required_argument, 0, 'd'},
	{"maxdist", required_argument, 0, 'D'},
	{"subreads", required_argument, 0, 'n'},
	{"minmatch", required_argument, 0, 'm'},
	{"minmatch2", required_argument, 0, 'p'},
	{"quality", no_argument, &QUALITY_SCALE, QUALITY_SCALE_LINEAR},
	{"anchors", required_argument, 0, 'A'},
	{"threads", required_argument, 0, 'T'},
	{"index-threshold", required_argument, 0, 'f'},
	{"output", required_argument, 0, 'o'},
	{"order",  required_argument,0, 'S'},
	{"cigar_len", required_argument, 0, 'C'},
	{"crtgap_penalty", required_argument, 0, 'G'},
	{"extgap_penalty", required_argument, 0, 'E'},
	{"mis_penalty", required_argument, 0, 'X'},
	{"match_score", required_argument, 0, 'Y'},
	{"multi", required_argument, 0, 'B'},
	{0, 0, 0, 0}
};



#ifdef MAKE_STANDALONE
int main(int argc,char ** argv)
#else
int main_align(int argc,char ** argv)
#endif
{
	char read_file [300], read2_file [300];
	char output_file [300];

	char index_prefix [300];
	unsigned int all_reads, all_tables;
	gene_allvote_t allvote;
	unsigned long long int processed_reads = 0, succeed_reads = 0;
	gene_input_t ginp, ginp2;
	unsigned int multi_best_reads;
	//gene_flat_t my_flat ;
	//create_flat_strip(&my_flat);

	int c;
	int option_index = 0;

	multi_best_reads = 1;
	TOTAL_SUBREADS = 10;
	ACCEPT_SUBREADS = 3;
	ACCEPT_MINOR_SUBREADS = 1;
	INDEX_THRESHOLD = 1024;
	INDEL_TOLERANCE = 6;
	MAX_METHYLATION_C_NUMBER = 0;
	REPORT_ONLY_UNIQUE = 0;
	APPLY_REPEATING_PENALTY = 1;
	read_file[0]=0;
	read2_file[0]=0;
	index_prefix[0]=0;
	output_file[0]=0;
//	all_reads = 1*1024*1024/10; //*14
	all_reads = 14*1024*1024;


	SUBREADprintf("\n");



	while ((c = getopt_long (argc, argv, "vQJuHB:S:d:D:n:m:p:f:R:r:i:o:T:E:I:A:P:G:X:b:Y:?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'v':
				print_version_info();
				return 0;
				break;
			case 'B':
				multi_best_reads = max(1, min(16, atoi(optarg)));
				break;
			case 'b':
				MAX_METHYLATION_C_NUMBER = atoi(optarg);
				break;
			case 'H':
				USE_BASEINDEX_BREAK_TIE = 1;
				break;
			case 'u':
				REPORT_ONLY_UNIQUE = 1;
				break;
			case 'J':
				REPORT_POTENTIAL_JUNCTION_READS = 1;
				break;
			case 'S':
				FIRST_READ_REVERSE = optarg[0]=='r'?1:0;
				SECOND_READ_REVERSE = optarg[1]=='f'?0:1;
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
				if(ALL_THREADS>64)ALL_THREADS=64;
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
				if(INDEL_TOLERANCE >= MAX_INDEL_TOLERANCE )INDEL_TOLERANCE=MAX_INDEL_TOLERANCE;
				if(INDEL_TOLERANCE>0)
					INDEL_TOLERANCE ++;
				break ;
			case 'Q':
				QUALITY_SCALE = QUALITY_SCALE_LINEAR;
				break;
			case 'P':
				if (optarg[0]=='3')
					FASTQ_FORMAT = FASTQ_PHRED33;
				else
					FASTQ_FORMAT = FASTQ_PHRED64;
				break;
			case 'p':
				ACCEPT_MINOR_SUBREADS = atoi(optarg);
				break;
			case 'G':
				DPALIGN_CREATEGAP_PENALTY = atoi(optarg);
				break;
			case 'E':
				DPALIGN_EXTENDGAP_PENALTY = atoi(optarg);
				break;
			case 'X':
				DPALIGN_MISMATCH_PENALTY = atoi(optarg);
				break;
			case 'Y':
				DPALIGN_MATCH_SCORE = atoi(optarg);
				break;
			case 'A':
				NUMBER_OF_ANCHORS_PAIRED = atoi(optarg);
				break;
			case '?':
				return -1 ;
		}

	if (QUALITY_SCALE == QUALITY_SCALE_NONE && USE_VALUE_ARRAY_INDEX ==0 )
	{
		if (NUMBER_OF_ANCHORS_PAIRED!=-1)
		{
			SUBREADprintf("Please do not specify '--anchors' or '-A' parameters without specifying '-a' and '-Q' options\n");
			return -1;
		}
	}
	else if (NUMBER_OF_ANCHORS_PAIRED==-1)
		NUMBER_OF_ANCHORS_PAIRED = 10;
	else if (NUMBER_OF_ANCHORS_PAIRED >= ANCHORS_NUMBER)
	{
		SUBREADprintf("The number of anchors must not be greater than %d.\n", (ANCHORS_NUMBER-1));
		return -1;
	}

	if (!read_file[0] || !index_prefix[0] || !output_file[0])
	{
		usage(argv[0]);

		return -1 ;
	}
  
	reads_density = guess_reads_density(read_file,0);
	if(reads_density<0)
		SUBREADprintf("Input file '%s' is not found or is in an incorrect format.\n", read_file);

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
		SUBREADprintf("Unable to open the index files in the %s space.\n",  ginp.space_type==GENE_SPACE_COLOR?"color":"base");
		return -1;
	}

	FILE * out_fp = fopen(output_file, "w");
	if (!out_fp)
	{
		SUBREADprintf("Unable to open the output file at '%s'.\n", output_file);
		return -1;
	}


	SUBREADprintf("Number of selected subreads = %d\n", TOTAL_SUBREADS);
	SUBREADprintf("Consensus threshold = %d\n", ACCEPT_SUBREADS);
	SUBREADprintf("Number of threads=%d\n", ALL_THREADS);
	if(multi_best_reads)
		SUBREADprintf ("Number of positions reported per multi-mapping read=%d\n", multi_best_reads);
	if(INDEL_TOLERANCE)
		SUBREADprintf("Number of indels allowed=%d\n", INDEL_TOLERANCE-1);
	if (QUALITY_SCALE==QUALITY_SCALE_LINEAR)
		SUBREADputs("Quality scale=linear\n\n");
	else if (QUALITY_SCALE==QUALITY_SCALE_LOG)
		SUBREADputs("Quality scale=exponential\n\n");
	else 	SUBREADputs("\n");

	if (read2_file[0])
	{
		if (MAX_PAIRED_DISTANCE <= MIN_PAIRED_DISTANCE)
		{
			SUBREADprintf ("The value of the '-D' option must be greater than that of the '-d' option. \n");
			return -1;
		}

		SUBREADprintf ("Performing paired-end alignment:\n");
		SUBREADprintf ("Maximum distance between reads=%d\n", MAX_PAIRED_DISTANCE);
		SUBREADprintf ("Minimum distance between reads=%d\n", MIN_PAIRED_DISTANCE);
		SUBREADprintf ("Threshold on number of subreads for a successful mapping (the minor end in the pair)=%d\n", ACCEPT_MINOR_SUBREADS);
		if (NUMBER_OF_ANCHORS_PAIRED > -1)
			SUBREADprintf ("Number of anchors=%d\n", NUMBER_OF_ANCHORS_PAIRED);
		SUBREADprintf ("The directions of the two input files are: %s, %s\n\n", FIRST_READ_REVERSE?"reversed":"forward", SECOND_READ_REVERSE?"reversed":"forward");
	}

#ifdef REPORT_ALL_THE_BEST
	SUBREADprintf("***** WARNING: the REPORT_ALL_THE_BEST switch is turned on. You need an extra 1 GBytes of RAM space for saving the temporary results. *****\n");
#endif
	
	if(init_allvote(&allvote, all_reads, multi_best_reads, INDEL_TOLERANCE )) return -1;
	if(read2_file[0])
		all_reads/=2;
		
	if(read2_file[0] && geinput_open(read2_file, &ginp2))
	{
		destory_allvote(&allvote );

		SUBREADprintf("Input file '%s' is not found or is in an incorrect format.\n", read2_file);
		return -1;
	}

	SUBREADfflush(stdout);

	begin_ftime = miltime();

	while (1)
	{
		char inbuff[1201];
		int processed_reads_block = 0;

		clear_allvote(&allvote);
		processed_reads_block = run_search_index(&ginp, read2_file[0] ? (&ginp2):NULL, index_prefix, &allvote, out_fp, processed_reads, all_tables, &succeed_reads);
		if (processed_reads_block<0)
			return -1;
		processed_reads += processed_reads_block ;

		// test if there no anyreads remaining
		unsigned long long int current_fp = ftello(ginp.input_fp);
		int rl = geinput_next_read(&ginp, NULL, inbuff, NULL);
		if (rl<0)
			break;
		fseeko(ginp.input_fp, current_fp, SEEK_SET);
	}

	geinput_close(&ginp);
	if(IS_DEBUG)
		SUBREADprintf("@LOG THE END. \n");
	else
		SUBREADprintf("\n\n %llu %s were processed in %.1f seconds.\nPercentage of mapped reads is %0.2f%%.\n\n", processed_reads, read2_file[0] ?"read pairs":"reads", miltime()-begin_ftime, succeed_reads*100.0/processed_reads/(read2_file[0]?2:1));

	SUBREADprintf("Done.\n");

	fclose(out_fp);
	destory_allvote(&allvote );
	return 0;
}
