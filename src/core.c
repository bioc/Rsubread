/***************************************************************

   The Subread software package is free software package: 
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

/***************************************************************

   The ASCII Art used in this file was generated using FIGlet and
   the big font, contributed by Glenn Chappell to FIGlet.
  
  ***************************************************************/
  
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/stat.h>
#include <ctype.h>



#include "subread.h"
#include "sublog.h"
#include "core.h"
#include "input-files.h"
#include "sorted-hashtable.h"

#include "core-indel.h"
#include "core-junction.h"

static struct option long_options[] =
{
	{"memory-optimisation",  required_argument, 0, 0},
	{0, 0, 0, 0}
};

int (*progress_report_callback)(int, int, int);

int is_result_in_PE(alignment_result_t *al)
{
	if(al->Score_H & 0x8000000000000000llu)return 1;
	return 0;
}

void core_version_number(char * program)
{
	SUBREADprintf("\n%s v%s\n\n" , program, SUBREAD_VERSION);
}

void warning_file_limit()
{
	struct rlimit limit_st;
	getrlimit(RLIMIT_NOFILE, & limit_st);

	{
		if(min(limit_st.rlim_cur , limit_st.rlim_max) < 400)
		{
			print_in_box(80,0,0,"WARNING This operation needs to open many files at same time,");
			print_in_box(80,0,0,"        but your OS only allows to open %d files.", min(limit_st.rlim_cur , limit_st.rlim_max));
			print_in_box(80,0,0,"        You can use command 'ulimit -n 500' to raise this limit");
			print_in_box(80,0,0,"        to 500, or the program may crash or terminate unexpectedly.");
			print_in_box(80,0,0,"");
		}
	}
}

void print_in_box(int line_width, int is_boundary, int is_center, char * pattern,...)
{
	va_list args;
	va_start(args , pattern);
	char is_R_linebreak=0, * content;

	content= malloc(1000);
	vsprintf(content, pattern, args);
	int is_R_code,x1,content_len = strlen(content), state, txt_len, is_cut = 0, real_lenwidth;

	is_R_code = 1;
	#ifdef MAKE_STANDALONE
		is_R_code = 0;
	#endif

	if(content_len>0&&content[content_len-1]=='\r'){
		content_len--;
		content[content_len] = 0;
		is_R_linebreak = 1;
	}

	if(content_len>0&&content[content_len-1]=='\n'){
		content_len--;
		content[content_len] = 0;
	}

	state = 0;
	txt_len = 0;
	real_lenwidth = line_width;
	for(x1 = 0; content [x1]; x1++)
	{
		char nch = content [x1];
		if(nch == CHAR_ESC)
			state = 1;
		if(state){
			real_lenwidth --;
		}else{
			txt_len++;
			
			if(txt_len == 80 - 6)
			{
				is_cut = 1;
			} 
		}

		if(nch == 'm' && state)
			state = 0;
	}

	if(is_cut)
	{
		state = 0;
		txt_len = 0;
		for(x1 = 0; content [x1]; x1++)
		{
			char nch = content [x1];
			if(nch == CHAR_ESC)
				state = 1;
			if(!state){
				txt_len++;
				if(txt_len == 80 - 9)
				{
					strcpy(content+x1, "\x1b[0m ...");
					content_len = line_width - 4;
					content_len = 80 - 4;
					line_width = 80;
					break;
				} 
			}
			if(nch == 'm' && state)
				state = 0;
		}
	}

	if(content_len==0 && is_boundary)
	{
		sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,is_boundary==1?"//":"\\\\");
		for(x1=0;x1<line_width-4;x1++)
			sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"=");
		sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,is_boundary==1?"\\\\":"//");
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"");

		free(content);
		return;
	}
	else if(is_boundary)
	{
		int left_stars = (line_width - content_len)/2 - 1;
		int right_stars = line_width - content_len - 2 - left_stars;
		sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,is_boundary==1?"//":"\\\\");
		for(x1=0;x1<left_stars-2;x1++) sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"=");
		sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"%c[36m", CHAR_ESC);
		sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," %s ", content);
		sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"%c[0m", CHAR_ESC);
		for(x1=0;x1<right_stars-2;x1++) sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"=");
		sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,is_boundary==1?"\\\\":"//");
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"");

		free(content);
		return;
	}

	int right_spaces, left_spaces;	
	if(is_center)
		left_spaces = (line_width - content_len)/2-2;
	else
		left_spaces = 1;

	right_spaces = line_width - 4 - content_len- left_spaces; 

	char spaces[81];
	memset(spaces , ' ', 80);
	spaces[0]='|';
	spaces[1]='|';
	spaces[80]=0;

	//sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"||");
	
	//for(x1=0;x1<left_spaces;x1++) sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," ");

	spaces[left_spaces+2] = 0;
	sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,spaces);

	if(is_R_code)
	{
		sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,content);
	}
	else
	{
		int col1w=-1;
		for(x1=0; content[x1]; x1++)
		{
			if(content[x1]==':')
			{
				col1w=x1;
				break;
			}
		}
		if(col1w>0 && col1w < content_len-1)
		{
			content[col1w+1]=0;
			sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,content);
			sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," ");
			sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"%c[36m", CHAR_ESC);
			sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,content+col1w+2);
			sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"%c[0m", CHAR_ESC);
		}
		else
			sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,content);
	}
//	for(x1=0;x1<right_spaces - 1;x1++) sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," ");
	
	memset(spaces , ' ', 80);
	spaces[79]='|';
	spaces[78]='|';
	
	right_spaces = max(1,right_spaces);
	if(is_R_linebreak)
		sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," %c[0m%s%c", CHAR_ESC, spaces + (78 - right_spaces + 1) ,CORE_SOFT_BR_CHAR);
	else
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," %c[0m%s", CHAR_ESC , spaces + (78 - right_spaces + 1));
	free(content);
}



int show_summary(global_context_t * global_context)
{

	if(progress_report_callback)
	{
		long long int all_reads_K = global_context -> all_processed_reads / 1000;
		float mapped_reads_percentage = global_context -> all_mapped_reads * 1./global_context -> all_processed_reads;
		if(global_context->input_reads.is_paired_end_reads) mapped_reads_percentage/=2;
		progress_report_callback(10, 900000, (int) (miltime()-global_context->start_time));
		progress_report_callback(10, 900010, (int) all_reads_K);
		progress_report_callback(10, 900011, (int) (10000.*mapped_reads_percentage));
	}

	print_in_box(80,0,1,"");
	print_in_box(89,0,1,"%c[36mCompleted successfully.%c[0m", CHAR_ESC, CHAR_ESC);
	print_in_box(80,0,1,"");
	print_in_box(80,2,1,"");
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "");
	print_in_box(80, 1,1,"Summary");
	print_in_box(80, 0,1,"");
	print_in_box(80, 0,0,"         Processed : %llu %s" , global_context -> all_processed_reads, global_context->input_reads.is_paired_end_reads?"fragments":"reads");
	print_in_box(81, 0,0,"            Mapped : %llu %s (%.1f%%%%)", global_context -> all_mapped_reads, global_context->input_reads.is_paired_end_reads?"fragments":"reads" ,  global_context -> all_mapped_reads*100.0 / global_context -> all_processed_reads);
	if(global_context->input_reads.is_paired_end_reads)
		print_in_box(80, 0,0,"  Correctly paired : %llu fragments", global_context -> all_correct_PE_reads);

	if(global_context->config.output_prefix[0])
	{
		if(global_context->config.entry_program_name == CORE_PROGRAM_SUBJUNC)
			print_in_box(80, 0,0,"         Junctions : %u", global_context -> all_junctions);
		if(global_context->config.do_fusion_detection)
			print_in_box(80, 0,0,"           Fusions : %u", global_context -> all_fusions);
		print_in_box(80, 0,0,"            Indels : %u", global_context -> all_indels);
	}
	
	print_in_box(80, 0,1,"");
	print_in_box(80, 0,0,"      Running time : %.1f minutes", (miltime()-global_context->start_time)*1./60);
	print_in_box(80, 0,1,"");
	print_in_box(80, 2,1,"http://subread.sourceforge.net/");
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "");

	return 0;
}


void show_progress(global_context_t * global_context, thread_context_t * thread_context, unsigned int current_read_no, int task)
{
	if(thread_context&&thread_context->thread_id)
	{
		SUBREADputs("show_progress can only be called by thread#0\n");
		return;
	}

	gene_input_t * ginp1 = thread_context?(thread_context->ginp1):(&global_context->input_reads.first_read_file);
	unsigned long long ginp1_file_pos = ftello(ginp1->input_fp);

	if(task == STEP_VOTING)
	{
		unsigned long long real_read_number = global_context -> all_processed_reads + current_read_no;// * global_context -> config.all_threads;
		if(real_read_number>1000)
			global_context -> input_reads . avg_read_length = (ginp1_file_pos - ginp1->read_chunk_start) * 1./real_read_number ;
	}

	unsigned long long total_file_size = global_context -> input_reads.first_read_file_size;
	unsigned long long guessed_all_reads = total_file_size / global_context -> input_reads . avg_read_length;
	//printf("FS=%llu; AVG=%f; GAR=%llu; CURRENT_NO=%u\n", total_file_size, global_context -> input_reads . avg_read_length , guessed_all_reads, current_read_no);

	unsigned long long current_block_start_file_offset = global_context -> current_circle_start_position_file1;

	unsigned long long guessed_this_chunk_all_reads = (total_file_size - current_block_start_file_offset) / global_context -> input_reads . avg_read_length ;
	if(guessed_this_chunk_all_reads > global_context ->config.reads_per_chunk) guessed_this_chunk_all_reads = global_context ->config.reads_per_chunk;

	unsigned long long guessed_all_reads_before_this_chunk = current_block_start_file_offset / global_context -> input_reads . avg_read_length ;
	unsigned long long reads_finished_in_this_chunk = (ginp1_file_pos - current_block_start_file_offset) / global_context -> input_reads . avg_read_length;//* global_context -> config.all_threads;
	if(task != STEP_VOTING)
	   reads_finished_in_this_chunk = (ginp1_file_pos - current_block_start_file_offset) / global_context -> input_reads . avg_read_length;
	

	unsigned long long finished_steps = guessed_all_reads_before_this_chunk * (global_context -> index_block_number * 6 + 4);
	
	if(task == STEP_VOTING)
		finished_steps += guessed_this_chunk_all_reads * global_context -> current_index_block_number * 6 ;//* global_context -> config.all_threads; 
	if(task >= STEP_ITERATION_ONE)
		finished_steps += guessed_this_chunk_all_reads * global_context -> index_block_number * 6 ;//*global_context -> config.all_threads;
	if(task > STEP_ITERATION_ONE)
		finished_steps += guessed_this_chunk_all_reads ;//* global_context -> config.all_threads;
	if(task > STEP_ITERATION_TWO)
		finished_steps += guessed_this_chunk_all_reads;

	if(task == STEP_VOTING)	finished_steps += reads_finished_in_this_chunk*5*global_context -> config.all_threads;
	if(task > STEP_ITERATION_TWO)
		finished_steps += reads_finished_in_this_chunk * 2;
	else
		finished_steps += reads_finished_in_this_chunk*global_context -> config.all_threads;

	unsigned long long guessed_all_steps = guessed_all_reads * (global_context -> index_block_number * 6 + 4);

	float finished_rate = finished_steps*1./guessed_all_steps;
	float reads_per_second = 0;

	if(task == STEP_VOTING)
		reads_per_second = finished_steps / (miltime() - global_context -> align_start_time) / (global_context -> index_block_number*6 + 4);
	else
		reads_per_second = finished_steps / (miltime() - global_context -> start_time) / (global_context -> index_block_number*6 + 4);
	//float exp_mins = (miltime() - global_context -> start_time) / finished_rate / 60; 

	//fprintf(stderr, "FINISHED=%llu, FINISHED_READS=%llu, ALL=%llu, ALLREADS=%llu, ALLCHUNKREADS=%llu; BEFORE_CHUK=%llu; CUR-BLK=%d; IND-BLK=%d\n", finished_steps, reads_finished_in_this_chunk, guessed_all_steps, guessed_all_reads,guessed_this_chunk_all_reads, guessed_all_reads_before_this_chunk,  global_context -> current_index_block_number , global_context -> index_block_number );

	if(current_read_no>1000 && !progress_report_callback)
		print_in_box(81,0,0, "%4d%%%% completed, %3d mins elapsed, total=%dk %s, rate=%2.1fk/s\r", (int)(finished_rate*100), (int)((miltime() - global_context -> start_time)/60),(int)(guessed_all_reads*1./1000), global_context -> input_reads.is_paired_end_reads?"frags":"reads", reads_per_second/1000, reads_finished_in_this_chunk);

	if(progress_report_callback)
	{
		progress_report_callback(10 ,task,(int)(10000*finished_rate));
		progress_report_callback(20 ,task,(int)(guessed_all_reads/1000));
	}
}


/*
int Xmain(int argc , char ** argv);

int main(int argc, char ** argv)
{
	int xk1=0;
	for(xk1=0; xk1<4; xk1++)
		Xmain(argc ,  argv);
	return 0;
}
*/


int parse_opts_core(int argc , char ** argv, global_context_t * global_context)
{
	int c;
	int option_index = 0;	

	optind = 1;
	opterr = 1;
	optopt = 63;

	while ((c = getopt_long (argc, argv, "ExsS:L:AHd:D:n:m:p:P:R:r:i:l:o:T:Q:I:t:B:b:Q:FcuUfM?", long_options, &option_index)) != -1)
	{
		switch(c)
		{
			case 'Q':
				global_context->config.multi_best_reads = atoi(optarg); 
				break;
			case 'H':
				global_context->config.use_hamming_distance_break_ties = 1;
				break;
			case 's':
				global_context->config.downscale_mapping_quality = 1;
				break;
			case 'M':
				global_context->config.do_big_margin_reporting = 1;
				global_context->config.do_big_margin_filtering_for_reads = 1;
				global_context->config.report_multi_mapping_reads = 0;
				break;

			case 'A':
				global_context->config.report_sam_file = 0;
				break;
			case 'E':
				global_context->config.max_mismatch_exonic_reads = 200;
				global_context->config.max_mismatch_junction_reads = 200;

				break;
			case 'f':
				global_context->config.max_mismatch_exonic_reads = 200;
				global_context->config.max_mismatch_junction_reads = 200;
				global_context->config.do_fusion_detection = 1;
				global_context->config.minimum_subread_for_first_read = 1;
				global_context->config.minimum_subread_for_second_read = 1;
				global_context->config.total_subreads = 28;
				global_context->config.report_no_unpaired_reads = 0;
				global_context->config.limited_tree_scan = 0;
				global_context->config.use_hamming_distance_in_exon = 1;
				break;
			case 'x':
				global_context->config.max_mismatch_exonic_reads = 10;
				global_context->config.max_mismatch_junction_reads = 1;
				global_context->config.ambiguous_mapping_tolerance = 39;
				global_context->config.extending_search_indels = 0;

				global_context->config.is_rna_seq_reads = 1;
				global_context->config.total_subreads = 14;
				global_context->config.minimum_subread_for_first_read = 3;
				global_context->config.minimum_subread_for_second_read = 1;
				global_context->config.high_quality_base_threshold = 990000;
				global_context->config.do_big_margin_filtering_for_junctions = 1;
				global_context->config.report_no_unpaired_reads = 0;
				global_context->config.limited_tree_scan = 1;
				global_context->config.use_hamming_distance_in_exon = 0;
				break;
			case 'S':
				global_context->config.is_first_read_reversed = optarg[0]=='r'?1:0;
				global_context->config.is_second_read_reversed = optarg[0]=='f'?0:1;
				break;
			case 'U':
				global_context->config.report_no_unpaired_reads = 1;
				break;
			case 'u':
				global_context->config.report_multi_mapping_reads = 0;
				break;
			case 'b':
				global_context->config.is_methylation_reads = 1;
				break;
			case 'D':
				global_context->config.maximum_pair_distance = atoi(optarg);
				break;
			case 'd':
				global_context->config.minimum_pair_distance = atoi(optarg);
				break;
			case 'n':
				global_context->config.total_subreads = atoi(optarg);
				break;
			case 'm':
				global_context->config.minimum_subread_for_first_read = atoi(optarg);
				break;
			case 'T':
				global_context->config.all_threads = atoi(optarg);
				if(global_context->config.all_threads <1) global_context->config.all_threads = 1;
				if(global_context->config.all_threads >32) global_context->config.all_threads = 32;

				break;
			case 'r':
				strncpy(global_context->config.first_read_file, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'R':
				global_context->input_reads.is_paired_end_reads = 1;
				strncpy(global_context->config.second_read_file, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'i':
				strncpy(global_context->config.index_prefix, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'o':
				strncpy(global_context->config.output_prefix, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'I':
				global_context->config.max_indel_length = atoi(optarg);

				if(global_context->config.max_indel_length <0)global_context->config.max_indel_length =0;
				if(global_context->config.max_indel_length > MAX_INSERTION_LENGTH)global_context->config.max_indel_length = MAX_INSERTION_LENGTH;
				if(global_context->config.max_indel_length > 16)
				{
					global_context->config.reassembly_subread_length = 12;
					global_context->config.reassembly_window_multiplex = 3;
					global_context->config.reassembly_start_read_number = 5;
					global_context->config.reassembly_tolerable_voting = 0;
					global_context->config.reassembly_window_alleles = 2;
					global_context->config.reassembly_key_length = 28;

					global_context->config.is_third_iteration_running = 1;
					global_context->config.max_mismatch_exonic_reads = 2;
					global_context->config.max_mismatch_junction_reads = 2;
					global_context->config.total_subreads = 28;
					global_context->config.do_big_margin_filtering_for_reads = 1;

					global_context->config.do_superlong_indel_detection = 0;
				}
				break;
			case 'P':
				if (optarg[0]=='3')
					global_context->config.phred_score_format = FASTQ_PHRED33;
				else
					global_context->config.phred_score_format = FASTQ_PHRED64;
				break;
			case 'p':
				global_context->config.minimum_subread_for_second_read = atoi(optarg);
				break;
			case 't':
				sprintf(global_context->config.temp_file_prefix, "%s/core-temp-sum-%06u-%05u", optarg, getpid(), rand());
				break;
			case 'F':
				global_context->config.is_second_iteration_running = 0;
				global_context->config.report_sam_file = 0;
				break;
			case 'B':
				global_context->config.is_first_iteration_running = 0;
				strcpy(global_context->config.medium_result_prefix, optarg);
				break;
			case 'c':
				global_context->config.space_type = GENE_SPACE_COLOR; 
				break;
				
			case 0:
				//if(strcmp("memory-optimisation",long_options[option_index].name)==0)
				//	memory_optimisation = atoi(optarg);
				break;
			case '?':
			default:
				return -1 ;
		}
	}


	return 0;
}

#ifdef MAKE_STANDALONE
int main_back(int argc , char ** argv)
{
	progress_report_callback = NULL;
#else
int subread_core_main(int argc , char ** argv)
{
	if(progress_report_callback) progress_report_callback(0,0,1);
#endif

	return core_main(argc, argv, parse_opts_core);
}


int check_configuration(global_context_t * global_context)
{
	int expected_type = FILE_TYPE_FAST_;
	if(global_context -> config.is_SAM_file_input && global_context -> config.is_BAM_input)
		expected_type = FILE_TYPE_BAM;
	else if(global_context -> config.is_SAM_file_input && !global_context -> config.is_BAM_input)
		expected_type = FILE_TYPE_SAM;
	else if(global_context -> config.is_gzip_fastq)
		expected_type = FILE_TYPE_GZIP_FAST_;
	
	if(global_context -> config.max_indel_length > 16)
		warning_file_limit();

	warning_file_type(global_context -> config.first_read_file, expected_type);
	if(global_context -> config.second_read_file[0])
	{
		if(expected_type==FILE_TYPE_FAST_ || expected_type==FILE_TYPE_GZIP_FAST_)
			warning_file_type(global_context -> config.second_read_file, expected_type);
		else
			print_in_box(80,0,0,"You should specify only one input SAM or BAM file.");
	}

	return 0;

}

int core_main(int argc , char ** argv, int (parse_opts (int , char **, global_context_t * )))
{
	//int memory_optimisation = 0;
	global_context_t * global_context;
	global_context = (global_context_t*)malloc(sizeof(global_context_t));
	init_global_context(global_context);


	int ret = parse_opts(argc , argv, global_context);
	if(ret) return ret;
	//global_context->config.reads_per_chunk = 200*1024;

	if(global_context->config.max_indel_length > 20 && !global_context->input_reads.is_paired_end_reads)
	{
		global_context->config.total_subreads = 28;
		global_context->config.reassembly_start_read_number = 3;
		global_context->config.do_superlong_indel_detection = 1;
	}


	ret = print_configuration(global_context);

	ret = ret || check_configuration(global_context);
	ret = ret || load_global_context(global_context);
	ret = ret || init_modules(global_context);


	ret = ret || read_chunk_circles(global_context);
	ret = ret || write_final_results(global_context);
	ret = ret || destroy_modules(global_context);
	ret = ret || destroy_global_context(global_context);
	ret = ret || show_summary(global_context);

	free(global_context);

	return ret;
}

// the new file name is written into fname then.

int convert_BAM_to_SAM(global_context_t * global_context, char * fname, int is_bam)
{
	char temp_file_name[200], *fline=malloc(3000), tmp_readname[MAX_READ_NAME_LEN+1];
	short tmp_flags;
	SamBam_FILE * sambam_reader;


	int is_file_sorted = 1;
	unsigned long long int read_no = 0;
	sambam_reader = SamBam_fopen(fname, is_bam?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);
	if(!sambam_reader)
	{
		free(fline);
		return -1;
	}
	while(1)
	{
		char * is_ret = SamBam_fgets(sambam_reader, fline, 2999, 1);
		if(!is_ret) break;
		if(fline[0]=='@')continue;
		if(is_SAM_unsorted(fline, tmp_readname, &tmp_flags , read_no)){
			if(tmp_flags & 1) global_context->input_reads.is_paired_end_reads = 1;
			is_file_sorted = 0;
			break;
		}
		if(tmp_flags & 1) global_context->input_reads.is_paired_end_reads = 1;
		else global_context->input_reads.is_paired_end_reads = 0;

		if(! global_context->input_reads.is_paired_end_reads) break;
		read_no++;
	}
	SamBam_fclose(sambam_reader);
	print_in_box(80,0,0,"The input %s file contains %s-end reads.", is_bam?"BAM":"SAM",  global_context->input_reads.is_paired_end_reads?"paired":"single");
	if(!is_file_sorted)
		print_in_box(80,0,0,"The input %s file is unsorted. Reorder it...", is_bam?"BAM":"SAM");
	else if(is_bam)
		print_in_box(80,0,0,"Convert the input BAM file...");

	if(is_bam || (global_context->input_reads.is_paired_end_reads && !is_file_sorted))
	{
		sprintf(temp_file_name, "%s.sam", global_context->config.temp_file_prefix);
		sambam_reader = SamBam_fopen(fname, is_bam?SAMBAM_FILE_BAM: SAMBAM_FILE_SAM);
		if(!sambam_reader){
			SUBREADprintf("Unable to open %s.\n", fname);
			return -1;
		}
		FILE * sam_fp = NULL;
		SAM_sort_writer writer;

		if(is_file_sorted) sam_fp = f_subr_open(temp_file_name,"w");
		else sort_SAM_create(&writer, temp_file_name, ".");

		while(1)
		{
			char * is_ret = SamBam_fgets(sambam_reader, fline, 2999, 1);
			if(!is_ret) break;
			if(is_file_sorted)
				fputs(fline, sam_fp);
			else{
				int ret = sort_SAM_add_line(&writer, fline, strlen(fline));
				if(ret<0) {
					print_in_box(80,0,0,"ERROR: read name is too long; check input format.");
					break;
				}
			}
		}

		if(is_file_sorted) 
			fclose(sam_fp);
		else{
			sort_SAM_finalise(&writer);
			if(writer.unpaired_reads)
				print_in_box(80,0,0,"%llu single-end mapped reads in reordering.", writer.unpaired_reads);
		}

		SamBam_fclose(sambam_reader);
		strcpy(fname, temp_file_name);
		global_context -> will_remove_input_file = 1;
	}



	free(fline);
	return 0;
}

int convert_GZ_to_FQ(global_context_t * global_context, char * fname, int half_n)
{
	int is_OK = 0;
	char temp_file_name[200];
	char * linebuff=malloc(3001);
	gzFile * rawfp = gzopen(fname, "r");
	
	if(rawfp)
	{
		print_in_box(80,0,0,"Decompress %s...", fname);
		sprintf(temp_file_name, "%s-%d.fq", global_context->config.temp_file_prefix, half_n);
		FILE * outfp = fopen(temp_file_name, "w");
		if(outfp)
		{
			while(1)
			{
				char * bufr =gzgets(rawfp, linebuff, 3000);
				if(!bufr) break;
				fputs(bufr, outfp);
			}
			is_OK = 1;
			fclose(outfp);
		}

		gzclose(rawfp);
	}

	strcpy(fname, temp_file_name);
	global_context -> will_remove_input_file |= (1<< (half_n-1));

	return !is_OK;
}

int core_geinput_open(global_context_t * global_context, gene_input_t * fp, int half_number, int is_init)
{
	char *fname;
	if(global_context->config.is_SAM_file_input)
	{

		fname = is_init?global_context ->config.first_read_file:global_context -> input_reads.first_read_file.filename;
		if(is_init && half_number == 1)
			if(convert_BAM_to_SAM(global_context, global_context ->config.first_read_file, global_context ->config.is_BAM_input)) return -1;
		if(!global_context->input_reads.is_paired_end_reads) half_number=0;
		return geinput_open_sam(fname, fp, half_number);

	}
	else
	{
		if(is_init)
		{
			if(global_context -> config.is_gzip_fastq)
				if(convert_GZ_to_FQ(global_context, (half_number==2)? global_context ->config.second_read_file : global_context ->config.first_read_file, half_number)) return -1;
			fname = (half_number == 2)?global_context -> config.second_read_file:global_context -> config.first_read_file;
		}
		else
			fname = (half_number == 2)?global_context -> input_reads.second_read_file.filename:global_context -> input_reads.first_read_file.filename;
		return geinput_open(fname, fp);
	}
}

void relocate_geinputs(global_context_t * global_context, thread_context_t * thread_context)
{
	if(thread_context)
	{
		thread_context -> reads_to_be_done = global_context -> input_reads.reads_in_blocks[thread_context -> thread_id];
		thread_context -> read_block_start = global_context -> input_reads.start_read_number_blocks[thread_context -> thread_id];

		thread_context -> ginp1 = (gene_input_t *)malloc(sizeof(gene_input_t));
		core_geinput_open(global_context, thread_context -> ginp1,1, 0);
		fseeko(thread_context -> ginp1 -> input_fp, global_context -> input_reads.first_file_blocks[thread_context -> thread_id], SEEK_SET);

		if(global_context -> input_reads.is_paired_end_reads)
		{
			thread_context -> ginp2 = (gene_input_t *)malloc(sizeof(gene_input_t));
			core_geinput_open(global_context, thread_context -> ginp2, 2, 0);
			fseeko(thread_context -> ginp2-> input_fp, global_context -> input_reads.second_file_blocks[thread_context -> thread_id], SEEK_SET);
		}
	}
}


int fetch_next_read_pair(global_context_t * global_context, thread_context_t * thread_context, gene_input_t* ginp1, gene_input_t* ginp2, int *read_len_1, int *read_len_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, char * qual_text_1, char *qual_text_2, int remove_color_head)
{
	int rl1, rl2=0;

	rl1 = geinput_next_read_trim(ginp1, read_name_1, read_text_1 , qual_text_1, global_context->config.read_trim_5, global_context->config.read_trim_3);
	if(global_context->config.space_type == GENE_SPACE_COLOR && remove_color_head)
	{
		if(isalpha(read_text_1[0]))
		{
			int xk1;
			for(xk1=2; read_text_1[xk1]; xk1++)
				read_text_1[xk1-2]=read_text_1[xk1];
			read_text_1[xk1-2]=0;
		}
	}

	if(ginp2)
	{
		rl2 = geinput_next_read_trim(ginp2, read_name_2, read_text_2 , qual_text_2, global_context->config.read_trim_5, global_context->config.read_trim_3);
		if(global_context->config.space_type == GENE_SPACE_COLOR && remove_color_head)
		{
			if(isalpha(read_text_2[0]))
			{
				int xk1;
				for(xk1=2; read_text_2[xk1]; xk1++)
					read_text_2[xk1-2]=read_text_2[xk1];
				read_text_2[xk1-2]=0;
			}
		}
	}

	if( global_context->config.space_type == GENE_SPACE_COLOR)
	{
		rl1-=1;rl2-=1;
	}


	if(rl1>0 && (rl2>0 || !ginp2))
	{
		if(global_context->config.is_first_read_reversed)
		{
			reverse_read(read_text_1, rl1, global_context->config.space_type);
			if(qual_text_1)
				reverse_quality(qual_text_1, rl1);
		}

		if(ginp2 && global_context->config.is_second_read_reversed)
		{
			reverse_read(read_text_2, rl2, global_context->config.space_type);
			if(qual_text_2)
				reverse_quality(qual_text_2, rl2);
		}

		*read_len_1 = rl1;
		if(ginp2)
			*read_len_2 = rl2;
		return 0;
	}
	else	return 1;
}

int write_final_results(global_context_t * context)
{
	if(context -> config.output_prefix[0])
	{
		write_indel_final_results(context);

		if(context -> config.entry_program_name == CORE_PROGRAM_SUBJUNC)
			write_junction_final_results(context);

		if(context -> config.do_fusion_detection)
			write_fusion_final_results(context);
	}
	
	return 0;
}

int get_soft_clipping_length(char* CIGAR)
{
	int nch;
	int cigar_cursor, tmp_int = 0;

	for(cigar_cursor = 0; (nch = CIGAR[cigar_cursor]) > 0;cigar_cursor++)
	{
		if(isdigit(nch)) tmp_int = tmp_int*10+(nch-'0');
		else
		{
			if(nch=='S') return tmp_int;
			return 0;
		}
	}
	return 0;
}

#define CIGAR_PERFECT_SECTIONS 12

int write_chunk_results(global_context_t * global_context)
{
	unsigned int read_number, sqr_read_number = 0, sqr_interval;
	gene_input_t * ginp1, * ginp2=NULL;
	char * additional_information = malloc(1800);
	short current_display_offset = 0, current_display_tailgate = 0;

	unsigned int out_poses[CIGAR_PERFECT_SECTIONS+1], xk1;
	char * out_cigars[CIGAR_PERFECT_SECTIONS+1], *out_mate_cigars[CIGAR_PERFECT_SECTIONS+1];
	char out_strands[CIGAR_PERFECT_SECTIONS+1];
	short out_read_lens[CIGAR_PERFECT_SECTIONS+1];

	for(xk1 = 0; xk1 < CIGAR_PERFECT_SECTIONS+1; xk1++) out_cigars[xk1] = malloc(100);
	for(xk1 = 0; xk1 < CIGAR_PERFECT_SECTIONS+1; xk1++) out_mate_cigars[xk1] = malloc(100);

	//if(global_context -> config.space_type == GENE_SPACE_COLOR && !global_context -> config.convert_color_to_base)
	//	current_display_offset = 1;


	ginp1 = &global_context->input_reads.first_read_file;
	if(global_context->input_reads.is_paired_end_reads)
		ginp2 = &global_context->input_reads.second_read_file;
	
	sqr_interval = global_context -> processed_reads_in_chunk/10; 

	for(read_number = 0; read_number < global_context -> processed_reads_in_chunk ; read_number++)
	{
		char read_text_1[MAX_READ_LENGTH+1], read_text_2[MAX_READ_LENGTH+1];
		char qual_text_1[MAX_READ_LENGTH+1], qual_text_2[MAX_READ_LENGTH+1];
		char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
		int read_len_1, read_len_2=0;
		int best_read_id = 0;
		int total_best_reads = 0;
		int read1_has_reversed = 0;
		int read2_has_reversed = 0;
	
		int is_second_read;
		int applied_reverse_space;

		if(sqr_read_number > sqr_interval)
		{
			show_progress(global_context, NULL, read_number, STEP_WRITE_CHUNK_RESULTS);
			sqr_read_number = 0;
		}

		sqr_read_number ++;
		fetch_next_read_pair(global_context, NULL, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, 0);
		//printf("ORAW1=%s\nORAW2=%s\n\n", read_text_1,read_text_2);

		applied_reverse_space = global_context->config.space_type;
		if(global_context -> config.convert_color_to_base)
		{
			colorread2base(read_text_1, read_len_1+1);
			if(global_context->input_reads.is_paired_end_reads)
				colorread2base(read_text_2, read_len_2+1);
			applied_reverse_space = GENE_SPACE_BASE;

		}
		//printf("BAS1=%s\nBAS2=%s\n\n", read_text_1,read_text_2);

		alignment_result_t * prime_result = _global_retrieve_alignment_ptr(global_context  , read_number, 0, 0);
		for(total_best_reads=0; total_best_reads<global_context -> config.multi_best_reads; total_best_reads++)
		{
			alignment_result_t * current_result = _global_retrieve_alignment_ptr(global_context  , read_number, 0, total_best_reads);
			if(total_best_reads > 0 && current_result->selected_votes < 1) break;
			if(total_best_reads > 0 && global_context -> input_reads.is_paired_end_reads && !is_result_in_PE(prime_result)) break;
		}

		for(best_read_id=0; best_read_id<global_context -> config.multi_best_reads; best_read_id++)
		{
			char mate_cigar_decompress[100];
			char current_cigar_decompress[100];
			for(is_second_read = 0; is_second_read < 1+ global_context -> input_reads.is_paired_end_reads; is_second_read++)
			{
				alignment_result_t *current_result, *mate_result;
				current_result = _global_retrieve_alignment_ptr(global_context  , read_number, is_second_read, best_read_id);
				mate_result    = _global_retrieve_alignment_ptr(global_context  , read_number,!is_second_read, best_read_id);
				if(best_read_id>0 && current_result->selected_votes < 1 && (is_second_read == 0 || !global_context -> input_reads.is_paired_end_reads)) break;
				if(best_read_id>0 && ( global_context -> input_reads.is_paired_end_reads&& !is_result_in_PE(prime_result))) break;

				char * current_name      = is_second_read ? read_name_2 : read_name_1;
				char * current_read_text = is_second_read ? read_text_2 : read_text_1;
				char * current_qual_text = is_second_read ? qual_text_2 : qual_text_1;
				int * current_has_reversed = is_second_read ? (&read2_has_reversed):(&read1_has_reversed);
				int current_read_len     = is_second_read ? read_len_2  : read_len_1;
				int mate_read_len	= is_second_read ? read_len_1  : read_len_2;
				unsigned int mate_linear_pos=0, current_linear_pos=0;
				char mate_strand = 0, current_strand = 0;

				char * current_chro_name=NULL, * mate_chro_name = NULL;
				unsigned int current_chro_offset=0, mate_chro_offset=0;
				char * current_CIGAR;
				int second_char = -1;

				additional_information[0]=0;

				int mask = 0;
				int is_mate_ok = 0;
				int is_current_ok = (current_result -> result_flags & CORE_IS_FULLY_EXPLAINED)?1:0;
				int current_repeated_times = 0;
				float current_final_quality = current_result -> final_quality;
				int current_soft_clipping_movement  =0, mate_soft_clipping_movement = 0;

				current_linear_pos = current_result -> selected_position;
				current_strand = (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

				if(global_context -> input_reads.is_paired_end_reads)
				{
					remove_backslash(current_name);

					mask |= SAM_FLAG_PAIRED_TASK;
					is_mate_ok = (mate_result -> result_flags & CORE_IS_FULLY_EXPLAINED)?1:0;

					if((mate_result->result_flags & CORE_IS_BREAKEVEN) && !global_context -> config.report_multi_mapping_reads)
						is_mate_ok = 0;

					int mate_repeated_times = 0;

					mate_strand = (mate_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
					mate_linear_pos = mate_result -> selected_position;

					if(global_context -> config.do_big_margin_filtering_for_reads)
						mate_repeated_times = is_ambiguous_voting(global_context, read_number, !is_second_read, mate_result->selected_votes, mate_result->confident_coverage_start, mate_result->confident_coverage_end, mate_read_len, mate_strand);

					if( global_context -> config.do_big_margin_filtering_for_reads && mate_repeated_times > 1)
						is_mate_ok = 0;


					if(global_context->config.report_no_unpaired_reads && !is_result_in_PE(current_result))
					{
						is_mate_ok = mate_result->selected_votes > current_result->selected_votes?is_mate_ok:0;
						is_current_ok = mate_result->selected_votes < current_result->selected_votes?is_current_ok:0;
					}


					if(is_second_read)
						mask |= SAM_FLAG_SECOND_READ_IN_PAIR;
					else
						mask |= SAM_FLAG_FIRST_READ_IN_PAIR;



					if(is_mate_ok)
					{

						int is_jumped = 0;
						char * mate_CIGAR;
						if( mate_result -> cigar_string[0] == -1)
						{
							is_jumped = 1;
							bincigar2cigar(mate_cigar_decompress, 100, mate_result -> cigar_string + 1, CORE_MAX_CIGAR_LEN - 1, mate_read_len);

							mate_CIGAR = mate_cigar_decompress;
						}
						else
						{
							bincigar2cigar(mate_cigar_decompress, 100, mate_result -> cigar_string, CORE_MAX_CIGAR_LEN, mate_read_len);
							mate_CIGAR = mate_cigar_decompress;

						}

						if(global_context -> config.do_fusion_detection)
						{
							chimeric_cigar_parts(global_context, mate_linear_pos, is_jumped ^ mate_strand, is_jumped , mate_cigar_decompress , out_poses, out_mate_cigars, out_strands, mate_read_len, out_read_lens);

							mate_linear_pos = out_poses[0];
							mate_strand = out_strands[0]=='-';
							mate_CIGAR = out_mate_cigars[0];
						}

						if(locate_gene_position_max(mate_linear_pos, &global_context -> chromosome_table, & mate_chro_name, & mate_chro_offset, mate_read_len))
							is_mate_ok = 0;

						mate_soft_clipping_movement = get_soft_clipping_length(mate_CIGAR);
						if(global_context -> config.space_type == GENE_SPACE_COLOR)
						{
							if( (!is_second_read)  +  mate_strand == 1 )
							{
								//mate_chro_offset ++;
							}

						}
						mate_chro_offset += mate_soft_clipping_movement;

					}
					if(is_mate_ok)
					{
						if(mate_strand + (!is_second_read) == 1) mask |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
						mate_chro_offset++;
		
					}

				}
				if(!is_mate_ok)
				{
					mate_chro_name = "*";
					mate_chro_offset = 0;
				}
		

				if(global_context -> config.do_big_margin_reporting || global_context -> config.do_big_margin_filtering_for_reads)
					current_repeated_times = is_ambiguous_voting(global_context, read_number, is_second_read, current_result->selected_votes, current_result->confident_coverage_start, current_result->confident_coverage_end, current_read_len, current_strand);
		
				if(is_current_ok)
					if((current_result->result_flags & CORE_IS_BREAKEVEN) && !global_context -> config.report_multi_mapping_reads)
						is_current_ok = 0;

				if(is_current_ok)
					if(global_context -> config.do_big_margin_filtering_for_reads && current_repeated_times>1)
						is_current_ok = 0;

				if(is_current_ok)
				{
					int is_first_section_jumped = 0;
					if( current_result -> cigar_string[0] == -1)
					{
						bincigar2cigar(current_cigar_decompress, 100, current_result -> cigar_string + 1, CORE_MAX_CIGAR_LEN, current_read_len);
						//current_linear_pos = reverse_cigar(current_linear_pos , current_cigar_decompress, current_cigar_decompress_new);
						current_CIGAR  = current_cigar_decompress;
						is_first_section_jumped = 1;
					}
					else
					{
						bincigar2cigar(current_cigar_decompress, 100, current_result -> cigar_string, CORE_MAX_CIGAR_LEN, current_read_len);
						current_CIGAR  = current_cigar_decompress;
					}


					if(global_context -> config.do_fusion_detection)
					{

						int chimeric_sections = chimeric_cigar_parts(global_context, current_linear_pos, is_first_section_jumped ^ current_strand, is_first_section_jumped , current_CIGAR , out_poses, out_cigars, out_strands, current_read_len, out_read_lens);

						//sprintf(additional_information + strlen(additional_information), "\tXX:Z:%s", current_cigar_decompress);

						for(xk1=1; xk1<chimeric_sections; xk1++)
						{
							unsigned int chimeric_pos;
							char * chimaric_chr;

							if(0==locate_gene_position_max(out_poses[xk1],& global_context -> chromosome_table, & chimaric_chr, & chimeric_pos, 0+out_read_lens[xk1]))
							{
								int soft_clipping_movement = 0;
								soft_clipping_movement = get_soft_clipping_length( out_cigars[xk1]);
								char strand_xor = (out_strands[xk1] == '-')^ is_second_read;
								sprintf(additional_information + strlen(additional_information), "\tCG:Z:%s\tCP:i:%u\tCT:Z:%c\tCC:Z:%s", out_cigars[xk1] , chimeric_pos + soft_clipping_movement + 1, strand_xor?'-':'+' , chimaric_chr );
							}
						}
							

						current_linear_pos = out_poses[0];
						current_strand = out_strands[0]=='-';
						current_CIGAR = out_cigars[0];
					}

					//printf("CURCIGAR=%s ; OK=%d\n", current_cigar_decompress, is_current_ok);
					if(locate_gene_position_max(current_linear_pos,& global_context -> chromosome_table, & current_chro_name, & current_chro_offset, current_read_len))
						is_current_ok = 0;
					//printf("CURCIGAR=%s ; OK=%d\n", current_cigar_decompress, is_current_ok);
				}

				if(is_current_ok)
				{
					if(current_strand + is_second_read == 1)
						mask |= SAM_FLAG_REVERSE_STRAND_MATCHED;
					if(best_read_id>0) mask |= SAM_FLAG_SECONDARY_MAPPING;
					//if(1639 == read_number)
					//	printf("R0=%d ; NEG=%d; SEC=%d\n",(*current_has_reversed), (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0, is_second_read);
					int current_need_reverse = current_strand;

					//if(current_result -> cigar_string[0]==-1)
					//	current_need_reverse = ! current_need_reverse;

					if(current_need_reverse + (*current_has_reversed) == 1)
					{
						reverse_read(current_read_text, current_read_len + global_context->config.convert_color_to_base, applied_reverse_space);
						reverse_quality(current_qual_text , current_read_len);
						(*current_has_reversed)=!(*current_has_reversed);
					//	if(1639 == read_number)
					//		printf("RR=%d\n\n",(*current_has_reversed));
					}
					current_chro_offset++;

					current_soft_clipping_movement = get_soft_clipping_length(current_CIGAR);
					current_chro_offset += current_soft_clipping_movement;
				}
				else
				{
					mask |= SAM_FLAG_UNMAPPED;
					int this_should_nagetive = is_second_read;

					if(global_context -> input_reads.is_paired_end_reads && global_context -> config.report_unmapped_using_mate_pos&& is_mate_ok)
					{

						// DO NOT SHOW CORRDINATE IF IT IS UNMAPPED.
						//
						//current_chro_name = mate_chro_name;
						//current_chro_offset = mate_chro_offset;
						current_chro_name = "*";
						current_chro_offset = 0;

						/////////////////////////////////////////

						this_should_nagetive = mate_strand;
						current_strand = mate_strand;
						if(this_should_nagetive + is_second_read ==1)
							mask |= SAM_FLAG_REVERSE_STRAND_MATCHED;
						else
							mask &= ~SAM_FLAG_REVERSE_STRAND_MATCHED;
					}
					else
					{
						current_chro_name = "*";
						current_chro_offset = 0;
					}

					

					if(this_should_nagetive + (*current_has_reversed) == 1)
					{
						reverse_read(current_read_text, current_read_len + global_context->config.convert_color_to_base, applied_reverse_space);
						reverse_quality(current_qual_text , current_read_len);
						(*current_has_reversed)=!(*current_has_reversed);
					}

					current_CIGAR = "*";
					current_final_quality=0;
				}

				if(global_context -> config.space_type == GENE_SPACE_COLOR)
				{
					if( is_second_read  +  current_strand == 1 )
					{
						//if(is_current_ok)
						//	current_chro_offset ++;
						current_display_offset = 0;
						current_display_tailgate = 1;
					}
					else
					{
						if(!global_context -> config.convert_color_to_base)
						{
							// the first base was a fake prime base; the second base is the first meaningful base.
							second_char = current_read_text[1];
							current_read_text[1] = color2char(second_char, current_read_text[0]);
						}
						current_display_offset = 1;
						current_display_tailgate = 0;
					}
				}

				if(global_context -> input_reads.is_paired_end_reads && !is_mate_ok)
				{
					mask |= SAM_FLAG_MATE_UNMATCHED;
					mate_strand = current_strand;

					if(is_current_ok){

						int mate_should_nagetive = mate_strand;
						if(mate_should_nagetive + (!is_second_read) ==1)
							mask |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
						else
							mask &= ~SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
					}


				}
				long long int mate_distance = 0;

				// DO NOT SHOW CORRDINATE IF IT IS UNMAPPED.
				//
					if( 0 && is_current_ok && global_context -> config.report_unmapped_using_mate_pos && global_context -> input_reads.is_paired_end_reads &&!is_mate_ok)
					{
						mate_distance = 0;

						
						mate_chro_name = current_chro_name;
						mate_chro_offset = current_chro_offset;
					}
				//////////////////////////////////////////////

		//		printf("CURCIGAR=%s ; OK=%d\n", current_cigar_decompress, is_current_ok);

				if(is_current_ok && global_context -> input_reads.is_paired_end_reads && is_mate_ok)
				{
					mate_distance = mate_chro_offset - mate_soft_clipping_movement;
					mate_distance -= current_chro_offset - current_soft_clipping_movement;
					mate_distance = abs(mate_distance);
					if(current_chro_offset >mate_chro_offset)
						mate_distance += current_read_len; 
					else
						mate_distance += mate_read_len; 

					if(current_chro_offset - current_soft_clipping_movement > mate_chro_offset - mate_soft_clipping_movement) mate_distance = -mate_distance;

					if(mate_distance>0)
					{
						mate_distance = max(mate_distance, current_read_len);
						mate_distance = max(mate_distance, mate_read_len);
					}
					else
					{
						mate_distance = min(mate_distance, -current_read_len);
						mate_distance = min(mate_distance, -mate_read_len);
					}
				}

				if((best_read_id == 0 ) && current_qual_text[0] && global_context -> config.phred_score_format == FASTQ_PHRED64)
					fastq_64_to_33(current_qual_text);
				if(!current_qual_text[0])
				{
					int xi2;
					for(xi2=current_display_offset;current_read_text[xi2 + current_display_tailgate];xi2++) current_qual_text[xi2 - current_display_offset] = 'I';
					current_qual_text[xi2 - current_display_offset]=0;
				}

				if(current_repeated_times)
					current_final_quality /= current_repeated_times;
				if(global_context->config.downscale_mapping_quality)
					current_final_quality=current_final_quality/5;

				if(mate_chro_name[0]!='*' && mate_chro_name == current_chro_name && global_context -> input_reads.is_paired_end_reads)
				{
					mate_chro_name="=";
					if(is_mate_ok && is_current_ok && abs(mate_distance) >= global_context->config. minimum_pair_distance && abs(mate_distance) <= global_context-> config.maximum_pair_distance  && current_strand == mate_strand)
					{
						mask |= SAM_FLAG_MATCHED_IN_PAIR;
						if(is_second_read)
							global_context->all_correct_PE_reads  ++;
					}


				}
				if(mate_chro_name[0]!='=') 
					mate_distance = 0;

				int tailgate_0 = -1;
				if(current_display_tailgate)
				{
					tailgate_0 = current_read_text[strlen(current_read_text) -1];
					current_read_text[strlen(current_read_text) - 1] = 0;
				}

				//if(read_number==2461)
				//printf("CURCIGAR=%s ; FINAL_CIGAR=%s ; OK=%d\n", current_cigar_decompress, current_CIGAR, is_current_ok);

				if(is_current_ok || best_read_id==0)
				{
					int seq_len = strlen(additional_information);
					if(global_context->config.do_big_margin_reporting)
						seq_len += sprintf(additional_information+seq_len, "\tSA:i:%d", current_repeated_times);
					seq_len += sprintf(additional_information+seq_len, "\tSH:i:%d\tNH:i:%d", (int)((current_result -> Score_H >> 17) & 0xfff), is_current_ok?total_best_reads:0);
				//	seq_len += sprintf(additional_information+seq_len, "\tSG:Z:%s\tSB:i:%d\tSC:i:%d\tSD:i:%d\tSN:i:%u\tNH:i:%d",current_cigar_decompress, current_result -> used_subreads_in_vote, current_result -> selected_votes, current_result -> noninformative_subreads_in_vote, read_number, is_current_ok?total_best_reads:0); 
					if(global_context->config.read_group_id[0])
						seq_len += sprintf(additional_information+seq_len, "\tRG:Z:%s", global_context->config.read_group_id);

					if(global_context -> config.is_BAM_output)
						SamBam_writer_add_read(global_context -> output_bam_writer, current_name, mask, current_chro_name, current_chro_offset, (int)current_final_quality, current_CIGAR, mate_chro_name, mate_chro_offset, mate_distance, current_read_len, current_read_text + current_display_offset, current_qual_text, additional_information+1);
					else
						sambamout_fprintf(global_context -> output_sam_fp , "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%s\t%s%s\n", current_name, mask, current_chro_name, current_chro_offset, (int)current_final_quality, current_CIGAR, mate_chro_name, mate_chro_offset, mate_distance, current_read_text + current_display_offset, current_qual_text, additional_information);
				}

				if(current_display_tailgate)
				{
					current_read_text[strlen(current_read_text)] = tailgate_0;
				}

				if(second_char > 0)
					current_read_text[1] = second_char;

				if(global_context->input_reads.is_paired_end_reads)
				{
					if(is_second_read)
						if(is_current_ok || is_mate_ok)
							global_context -> all_mapped_reads++;
				}
				else
				{
					if(is_current_ok)
						global_context -> all_mapped_reads++;
				}
			} 
		}
	}

	free(additional_information);
	for(xk1 = 0; xk1 < CIGAR_PERFECT_SECTIONS; xk1++) free(out_cigars[xk1]);
	for(xk1 = 0; xk1 < CIGAR_PERFECT_SECTIONS; xk1++) free(out_mate_cigars[xk1]);
	return 0;
}
void init_chunk_scanning_parameters(global_context_t * global_context, thread_context_t * thread_context, gene_input_t ** ginp1, gene_input_t ** ginp2, unsigned int * read_block_start, unsigned int * reads_to_be_done)
{
	*ginp2 = NULL;
	*ginp1 = thread_context?thread_context->ginp1: & global_context->input_reads.first_read_file;
	if(global_context->input_reads.is_paired_end_reads)
		*ginp2 = thread_context?thread_context->ginp2:& global_context->input_reads.second_read_file;

	*read_block_start = thread_context?thread_context->read_block_start:0;
	*reads_to_be_done = thread_context?thread_context->reads_to_be_done:global_context -> config.reads_per_chunk;
}

gene_value_index_t * find_current_value_index(global_context_t * global_context, unsigned int pos, int len)
{
	int block_no;
	if(global_context->index_block_number<2)
	{
		unsigned index_begin = global_context -> all_value_indexes [0] . start_base_offset; 
		unsigned index_end = global_context -> all_value_indexes [0] . start_base_offset + global_context -> all_value_indexes [0] . length;

		if(pos>=index_begin && pos + len<index_end)
			return & global_context->all_value_indexes [0];
		else return NULL;
	}
	else
		for(block_no=0; block_no<global_context->index_block_number; block_no++)
		{
			unsigned index_begin = global_context -> all_value_indexes [block_no] . start_base_offset; 
			unsigned index_end = global_context -> all_value_indexes [block_no] . start_base_offset + global_context -> all_value_indexes [block_no] . length;
			if((block_no == 0 && pos >=  index_begin && pos < index_end - 1000000) ||
			   (block_no>0 && block_no<global_context->index_block_number-1 && pos >= index_begin+ 1000000 && pos < index_end - 1000000) ||
			   (block_no == global_context->index_block_number-1 && pos >= index_begin + 1000000 && pos < index_end ))
			{
				return & global_context -> all_value_indexes [block_no];
			}
		}
	return NULL;
}
//this function selects the correct all_value_indexes from global_context and put it to global_context or thread_context if thread_context is not NULL.
int locate_current_value_index(global_context_t * global_context, thread_context_t * thread_context, alignment_result_t * result, int rlen)
{
	int block_no;

	if(global_context->index_block_number<2)
	{
		unsigned index_begin = global_context -> all_value_indexes [0] . start_base_offset; 
		unsigned index_end = global_context -> all_value_indexes [0] . start_base_offset + global_context -> all_value_indexes [0] . length;

		if(result->selected_position>=index_begin && result->selected_position + rlen<index_end)
		{
			if(thread_context)thread_context->current_value_index = & global_context->all_value_indexes [0];
			else global_context->current_value_index =& global_context->all_value_indexes [0];
			return 0;
		}
		else return 1;
	}
	for(block_no=0; block_no<global_context->index_block_number; block_no++)
	{
		unsigned index_begin = global_context -> all_value_indexes [block_no] . start_base_offset; 
		unsigned index_end = global_context -> all_value_indexes [block_no] . start_base_offset + global_context -> all_value_indexes [block_no] . length;
		if((block_no == 0 && result->selected_position >=  index_begin && result->selected_position < index_end - 1000000) ||
		   (block_no>0 && block_no<global_context->index_block_number-1 && result->selected_position >= index_begin+ 1000000 && result->selected_position < index_end - 1000000) ||
		   (block_no == global_context->index_block_number-1 && result->selected_position >= index_begin + 1000000 && result->selected_position < index_end ))
		{
			if(thread_context)thread_context->current_value_index =& global_context -> all_value_indexes [block_no];
			else global_context->current_value_index =& global_context -> all_value_indexes [block_no];
			return 0;
		}
	}
	return 1;
}




int do_iteration_one(global_context_t * global_context, thread_context_t * thread_context)
{
	unsigned int reads_to_be_done = 0, read_block_start = 0;
	int ret;
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	unsigned int current_read_number;
	char read_text_1[MAX_READ_LENGTH+1], read_text_2[MAX_READ_LENGTH+1];
	char qual_text_1[MAX_READ_LENGTH+1], qual_text_2[MAX_READ_LENGTH+1];
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	int read_len_1, read_len_2=0;
	int need_junction_step = global_context -> config.is_rna_seq_reads || global_context -> config.do_fusion_detection;
	int sqr_interval, sqr_read_number = 0;

	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2, & read_block_start, & reads_to_be_done);
	sqr_interval = global_context -> processed_reads_in_chunk/10/ global_context -> config.all_threads;

	for(current_read_number = read_block_start; current_read_number < reads_to_be_done + read_block_start ; current_read_number++)
	{
		int is_second_read;

		sqr_read_number++;
		ret = fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, 1);
		// if no more reads
		if(ret)
			break;

		for (is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
		{
			int best_read_id, is_reversed_already = 0;
			for(best_read_id = 0; best_read_id < global_context -> config.multi_best_reads; best_read_id++)
			{
				alignment_result_t *current_result = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read, best_read_id); 

	//	if(strcmp("a4", read_name_1) == 0)
	//		printf("IDR=%u   VOT=%d  PAIR#=%u\n", current_result->selected_position, current_result->selected_votes, current_read_number);



				if(current_result -> selected_votes<1) break;
				if(!global_context->config.report_multi_mapping_reads)if(current_result -> result_flags & CORE_IS_BREAKEVEN) continue;

				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				char * current_read_name = is_second_read?read_name_2:read_name_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;

				if(current_result->selected_votes < global_context->config.minimum_subread_for_first_read)
					continue;
				int is_negative_strand = (current_result  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

				if(is_negative_strand + is_reversed_already == 1)
				{
					reverse_read(current_read, current_rlen ,  global_context->config.space_type);
					reverse_quality(current_qual, current_rlen);
					is_reversed_already = !is_reversed_already;
				}

				if(locate_current_value_index(global_context, thread_context, current_result, current_rlen))
				{
				//	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR, "Read position excesses index boundary : %u (%s : %s). V=%d", current_result -> selected_position, current_read_name, is_second_read?"SECOND":"FIRST", current_result -> selected_votes);
					continue;
				}

				find_new_indels(global_context, thread_context, current_read_number, current_read_name, current_read, current_qual, current_rlen, is_second_read, best_read_id);
				if(need_junction_step)
					find_new_junctions(global_context, thread_context, current_read_number, current_read, current_qual, current_rlen, is_second_read, best_read_id);
			}
		}
		
		if(!thread_context || thread_context->thread_id == 0)
		{
			if(sqr_read_number > sqr_interval)	
			{
				show_progress(global_context, thread_context, current_read_number, STEP_ITERATION_ONE);
				sqr_read_number = 0;
			}
		}

	}



	return 0;
}



int finish_iteration_three(global_context_t * global_context, thread_context_t * thread_context)
{
	return 0;
}
int do_iteration_three(global_context_t * global_context, thread_context_t * thread_context)
{
	unsigned int reads_to_be_done = 0, read_block_start = 0;
	int ret;
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	unsigned int current_read_number;
	char read_text_1[MAX_READ_LENGTH+1], read_text_2[MAX_READ_LENGTH+1];
	char qual_text_1[MAX_READ_LENGTH+1], qual_text_2[MAX_READ_LENGTH+1];
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	int read_len_1, read_len_2=0, sqr_interval, sqr_read_number = 0;

	//unsigned int low_index_border = global_context -> current_value_index -> start_base_offset;
	//unsigned int high_index_border = global_context -> current_value_index -> start_base_offset + global_context -> current_value_index -> length; 

	print_in_box(80,0,0,"Prepare for long indel deleteion...");
	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2, & read_block_start, & reads_to_be_done);
	sqr_interval = global_context -> processed_reads_in_chunk/10/ global_context -> config.all_threads;

	for(current_read_number = read_block_start; current_read_number < reads_to_be_done + read_block_start ; current_read_number++)
	{
		int is_second_read;

		sqr_read_number++;
		ret = fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, 1);
		if(ret)
			break;

		int best_read_id = 0;
		for (is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
		{
			int is_reversed_already = 0;
			for(best_read_id = 0; best_read_id < global_context -> config.multi_best_reads; best_read_id++)
			{
				alignment_result_t * current_result = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read, best_read_id);
				alignment_result_t * mate_result   = _global_retrieve_alignment_ptr(global_context, current_read_number,!is_second_read, best_read_id); 

				if(best_read_id && (current_result->selected_votes <1)) break;

				char * current_read_name =  is_second_read?read_name_2 : read_name_1;
				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;
				int mate_rlen = is_second_read?read_len_1:read_len_2;

				//if(!is_ambiguous_voting(global_context , current_read_number, is_second_read, current_result -> selected_votes, current_result -> confident_coverage_start, current_result -> confident_coverage_end, current_rlen, (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0))
				{
					// do local reassambly
					// a potential long-indel read has to have minimum supporting subreads, but not as many as total_subread - 1

					if(current_result->selected_votes >= global_context->config.minimum_subread_for_first_read)
					{
						int is_negative_strand = ((current_result  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0);

						if(is_negative_strand + is_reversed_already ==1)
						{
							reverse_read(current_read, current_rlen,  global_context->config.space_type);
							reverse_quality(current_qual, current_rlen);
							is_reversed_already=!is_reversed_already;
						}

						build_local_reassembly(global_context , thread_context , current_read_number, current_read_name , current_read , current_qual , current_rlen, 0 , is_second_read, best_read_id, 0);

					}
					else if(global_context -> input_reads.is_paired_end_reads && is_result_in_PE(current_result) && current_result -> selected_votes >= global_context->config.minimum_subread_for_second_read)
					{
						int is_negative_strand = (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
						if(is_negative_strand + is_reversed_already==1)
						{
							reverse_read(current_read, current_rlen,  global_context->config.space_type);
							reverse_quality(current_qual, current_rlen);
							is_reversed_already=!is_reversed_already;
						}

						build_local_reassembly(global_context , thread_context , current_read_number, current_read_name , current_read , current_qual , current_rlen , 0, is_second_read, best_read_id, 0);
					}
					else if(global_context -> input_reads.is_paired_end_reads && mate_result -> selected_votes >= global_context->config.minimum_subread_for_first_read)
					{
						int is_negative_strand = ((mate_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0);
						if(is_negative_strand+is_reversed_already==1)
						{
							reverse_read(current_read, current_rlen,  global_context->config.space_type);
							reverse_quality(current_qual, current_rlen);
							is_reversed_already=!is_reversed_already;
						}

						build_local_reassembly(global_context , thread_context , current_read_number , current_read_name , current_read , current_qual , current_rlen, mate_rlen , is_second_read, best_read_id, 1);
					}
				}	
			}
		}

		if(!thread_context || thread_context->thread_id == 0)
		{
			if(sqr_read_number>sqr_interval)
			{
				show_progress(global_context, thread_context, current_read_number, STEP_ITERATION_THREE);
				sqr_read_number = 0;
			}
		}
	}

	return 0;

}


int do_iteration_two(global_context_t * global_context, thread_context_t * thread_context)
{
	unsigned int reads_to_be_done = 0, read_block_start = 0;
	int ret;
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	unsigned int current_read_number;
	char read_text_1[MAX_READ_LENGTH+1], read_text_2[MAX_READ_LENGTH+1];
	char qual_text_1[MAX_READ_LENGTH+1], qual_text_2[MAX_READ_LENGTH+1];
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	int read_len_1, read_len_2=0;
	int sqr_interval, sqr_read_number=0;

	//unsigned int low_index_border = global_context -> current_value_index -> start_base_offset;
	//unsigned int high_index_border = global_context -> current_value_index -> start_base_offset + global_context -> current_value_index -> length; 

	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2, & read_block_start, & reads_to_be_done);
	sqr_interval = global_context -> processed_reads_in_chunk/10/ global_context -> config.all_threads;

	for(current_read_number = read_block_start; current_read_number < reads_to_be_done + read_block_start ; current_read_number++)
	{
		int is_second_read;

		sqr_read_number++;
		ret = fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2,1);
		// if no more reads
		if(ret)
			break;

		int best_read_id=0;
		for (is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
		{
			int is_reversed_already = 0;
			for(best_read_id = 0; best_read_id < global_context -> config.multi_best_reads; best_read_id++)
			{
				alignment_result_t *current_result = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read, best_read_id);
				if(best_read_id > 0 && current_result -> selected_votes==0)break;
				if(!global_context->config.report_multi_mapping_reads)if(current_result -> result_flags & CORE_IS_BREAKEVEN) continue;

				char * current_read_name =  is_second_read?read_name_2 : read_name_1;
				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;


				if(current_result->selected_votes < global_context->config.minimum_subread_for_second_read)
				{
					//printf("RESET0 %d %d\n", current_result->selected_votes , global_context->config.total_subreads * global_context->config.minimum_support_for_second_read);
					current_result->selected_votes = 0;
					continue;
				}
				if(locate_current_value_index(global_context, thread_context, current_result, current_rlen))
				{
					//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR, "Read position excesses index boundary.");
					continue;
				}

				int is_negative_strand = (current_result  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

				if(is_negative_strand + is_reversed_already == 1)
				{
					reverse_read(current_read, current_rlen,  global_context->config.space_type);
					reverse_quality(current_qual, current_rlen);
					is_reversed_already = !is_reversed_already;
				}

				explain_read(global_context, thread_context, current_read_number, current_rlen, current_read_name, current_read, current_qual, is_second_read, best_read_id, is_negative_strand);
			}
		}
		
		if(!thread_context || thread_context->thread_id == 0)
		{
			if(sqr_read_number>sqr_interval)
			{
				show_progress(global_context, thread_context, current_read_number, STEP_ITERATION_TWO);
				sqr_read_number=0;
			}
		}
	}

	return 0;
}

int core_get_subread_quality(global_context_t * global_context, thread_context_t * thread_context, char * qual, int qual_len)
{
	int x1, ret=1;

	if(!qual)return 1;
	if(!qual[0])return 1;

	int offset = (global_context->config.phred_score_format == FASTQ_PHRED33)?33:64; 

	for(x1=0; (x1 < qual_len) && qual[x1]; x1++)
		ret +=  (qual[x1] - offset);

	return  ret;
}

int do_voting(global_context_t * global_context, thread_context_t * thread_context)
{
	unsigned int reads_to_be_done = 0, read_block_start = 0;
	int ret, xk1;
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	unsigned int current_read_number;
	char * read_text_1, * read_text_2;
	char * qual_text_1, * qual_text_2;
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	int read_len_1, read_len_2=0;
	unsigned int processed_reads=0;
	int min_first_read_votes = global_context -> config.minimum_subread_for_first_read; 
	int min_second_read_votes = global_context -> config.minimum_subread_for_second_read; 
	int voting_max_indel_length = min(16, global_context->config.max_indel_length);
	int sqr_interval=10000, sqr_read_number = 0;

	read_text_1 = malloc(MAX_READ_LENGTH+1);
	read_text_2 = malloc(MAX_READ_LENGTH+1);
	qual_text_1 = malloc(MAX_READ_LENGTH+1);
	qual_text_2 = malloc(MAX_READ_LENGTH+1);

	gene_vote_t * vote_1, * vote_2, * vote_fg;

	vote_1 = (gene_vote_t *) malloc(sizeof(gene_vote_t));
	vote_2 = (gene_vote_t *) malloc(sizeof(gene_vote_t));
	vote_fg = (gene_vote_t *) malloc(sizeof(gene_vote_t));

	if(!(vote_1&&vote_2))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_FATAL,"Cannot allocate voting memory.");
		return 1;
	}

	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2, & read_block_start, & reads_to_be_done);

	unsigned int low_index_border = global_context -> current_value_index -> start_base_offset;
	unsigned int high_index_border = global_context -> current_value_index -> start_base_offset + global_context -> current_value_index -> length; 
	int has_second_read = 1 + global_context -> input_reads.is_paired_end_reads;
	//int need_junction_step = global_context -> config.is_rna_seq_reads || global_context -> config.do_fusion_detection;

	if(thread_context)
		thread_context -> current_value_index = global_context -> current_value_index;

	for(current_read_number = read_block_start; current_read_number < reads_to_be_done + read_block_start ; current_read_number++)
	{
		int is_second_read;
		int subread_no;
		int is_reversed;

		ret = fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2,1);
		if(ret)
			break;

		//printf("%s\t%d\n%s\t%d\n", read_name_1, thread_context -> thread_id, read_name_2,  thread_context -> thread_id);

		for(is_reversed = 0; is_reversed<2; is_reversed++)
		{
			//printf("MAP_REA = %s / %s\n", read_text_1, read_text_2);
			if(is_reversed==0 || !global_context->config.do_fusion_detection)
			{
				init_gene_vote(vote_1);
				if(global_context -> input_reads.is_paired_end_reads) init_gene_vote(vote_2);
			}


			for (is_second_read = 0; is_second_read < has_second_read; is_second_read ++)
			{
				gene_vote_t * current_vote = is_second_read?vote_2: vote_1;
				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;
				int subread_step;
				if(current_rlen<= EXON_LONG_READ_LENGTH)
				{
					subread_step = (((current_rlen - 18)<<16))/(global_context -> config.total_subreads -1);
					if(subread_step<(3<<16))subread_step = 3<<16;
				}else{
					subread_step = 6<<16;
					if(((current_rlen - 18)<<16) / subread_step > 62)
						subread_step = ((current_rlen - 18)<<16)/62;
				}

				#ifdef NEED_SUBREAD_STATISTIC
				int allsubreads_for_each_gap [3], noninformative_subreads_for_each_gap[3];
				#endif

				int applied_subreads = 1 + ((current_rlen - 18)<<16) / subread_step;
				unsigned int current_high_border = high_index_border -  current_rlen;

				if(global_context->config.is_rna_seq_reads && current_rlen > EXON_LONG_READ_LENGTH && global_context->config.all_threads<2)
					core_fragile_junction_voting(global_context, thread_context, current_read, current_qual, current_rlen, is_reversed, global_context->config.space_type, low_index_border, current_high_border, vote_fg);

				#ifdef NEED_SUBREAD_STATISTIC
				for(xk1=0;xk1<3;xk1++)
				{
					allsubreads_for_each_gap[xk1]=0;
					noninformative_subreads_for_each_gap[xk1]=0;
				}
				#endif

				for(subread_no=0; subread_no < applied_subreads ; subread_no++)
				{
					for(xk1=0; xk1<GENE_SLIDING_STEP ; xk1++)
					{

				#ifdef NEED_SUBREAD_STATISTIC
						current_vote -> noninformative_subreads = noninformative_subreads_for_each_gap[xk1];
						current_vote -> all_used_subreads = allsubreads_for_each_gap[xk1];
				#endif

						int subread_offset = ((subread_step * subread_no) >> 16);
						subread_offset -= subread_offset%3 - xk1; 

						int subread_quality = 1;
						char * subread_string = current_read + subread_offset;

						gehash_key_t subread_integer = genekey2int(subread_string, global_context->config.space_type);


						if(global_context -> config.use_quality_score_break_ties)
						{
							char * quality_string_subr = current_qual + subread_offset;
							subread_quality = core_get_subread_quality(global_context, thread_context, quality_string_subr, 16);
						}



						if(global_context->config.is_methylation_reads)
							gehash_go_q_CtoT(global_context->current_index, subread_integer , subread_offset, current_rlen, is_reversed, current_vote, 1, subread_quality, 0xffffff, voting_max_indel_length, subread_no, 1,  low_index_border, high_index_border - current_rlen);
						else
							gehash_go_q(global_context->current_index, subread_integer , subread_offset, current_rlen, is_reversed, current_vote, subread_quality, voting_max_indel_length, subread_no,  low_index_border, current_high_border);



				#ifdef NEED_SUBREAD_STATISTIC
						noninformative_subreads_for_each_gap[xk1] = current_vote -> noninformative_subreads;
						allsubreads_for_each_gap[xk1] = current_vote -> all_used_subreads;
				#endif


					}
				}

				//puts("");

				#ifdef NEED_SUBREAD_STATISTIC
				short max_noninformative_subreads = -1;
				unsigned char max_all_subreads = 0;

				for(xk1=0;xk1<3;xk1++)
					if(noninformative_subreads_for_each_gap[xk1] > max_noninformative_subreads)
					{
						max_noninformative_subreads = noninformative_subreads_for_each_gap[xk1];
						max_all_subreads = allsubreads_for_each_gap[xk1];
					}

				current_vote -> noninformative_subreads = max_noninformative_subreads;
				current_vote -> all_used_subreads = max_all_subreads;
				#endif
			}



			if(is_reversed==1 || !global_context->config.do_fusion_detection)
			{
				//if(strcmp(read_name_1,"b1")==0)

				//if(current_read_number == 17296){
				//print_votes(vote_1, global_context -> config.index_prefix);
				//print_votes(vote_2, global_context -> config.index_prefix);
				//}

				//finalise_vote(vote_1);
 
				/*
				if(407229 == current_read_number){
					fprintf(stderr, "TABLE_ITEMS=%llu\n", global_context->current_index->current_items);
					print_votes(vote_1, global_context -> config.index_prefix);
					//print_votes(vote_2, global_context -> config.index_prefix);
				}*/
				//if(global_context -> input_reads.is_paired_end_reads) finalise_vote(vote_2);

				if(global_context -> input_reads.is_paired_end_reads)
				{
					if((vote_1->max_vote >= min_first_read_votes || vote_2->max_vote >= min_first_read_votes)&& 
						(vote_1->max_vote >= min_second_read_votes && vote_2->max_vote >= min_second_read_votes))
							process_voting_junction(global_context, thread_context, current_read_number, vote_1, vote_2, read_name_1, read_name_2, read_text_1, read_text_2, read_len_1, read_len_2, is_reversed);
					else
						process_voting_junction(global_context, thread_context, current_read_number, vote_1, vote_2, read_name_1, read_name_2, read_text_1, read_text_2, read_len_1, read_len_2, is_reversed);
				}
				else{
					if(vote_1->max_vote >= min_first_read_votes)
						process_voting_junction(global_context, thread_context, current_read_number, vote_1, vote_2, read_name_1, NULL ,  read_text_1, NULL, read_len_1, 0, is_reversed);
					else if(_global_retrieve_alignment(global_context, current_read_number, 0,0).selected_votes < 1)
					{
						_global_retrieve_alignment_ptr(global_context, current_read_number, 0,0)->used_subreads_in_vote = max(_global_retrieve_alignment(global_context, current_read_number, 0,0).used_subreads_in_vote, vote_1 -> all_used_subreads);
						_global_retrieve_alignment_ptr(global_context, current_read_number, 0,0)->noninformative_subreads_in_vote = max(_global_retrieve_alignment(global_context, current_read_number, 0,0).noninformative_subreads_in_vote, vote_1 -> noninformative_subreads);
					}
				}


				//if(current_read_number == 0) print_votes(vote_1, global_context -> config.index_prefix);
				//if(current_read_number == 0) printf("V0=%d; %d\n", _global_retrieve_alignment_ptr(global_context, current_read_number, 0,0)->selected_votes, _global_retrieve_alignment_ptr(global_context, current_read_number, 1,0)->selected_votes);
			}
		
			if(is_reversed == 0)
			{
				reverse_read(read_text_1, read_len_1,  global_context->config.space_type);
				reverse_quality(qual_text_1, read_len_1);

				if(global_context -> input_reads.is_paired_end_reads)
				{
					reverse_read(read_text_2, read_len_2,  global_context->config.space_type);
					reverse_quality(qual_text_2, read_len_2);
				}
			}

		}
		
		if(!thread_context || thread_context->thread_id == 0)
		{
			if(sqr_read_number > sqr_interval)
			{
				show_progress(global_context, thread_context, processed_reads, STEP_VOTING);
				sqr_read_number = 0;
				unsigned long long total_file_size = global_context -> input_reads.first_read_file_size;
				unsigned long long guessed_all_reads = total_file_size / global_context -> input_reads . avg_read_length;// / (1+global_context -> config.is_SAM_file_input);
				sqr_interval = guessed_all_reads / global_context -> config.all_threads/10;
			}
		}




		sqr_read_number++;
		processed_reads++;

	}

	if(thread_context)
		thread_context -> processed_reads_in_chunk = processed_reads;
	else
		global_context -> processed_reads_in_chunk = processed_reads;

	free(vote_1);
	free(vote_2);
	free(vote_fg);
	free(read_text_1);
	free(read_text_2);
	free(qual_text_1);
	free(qual_text_2);

	return 0;
}

void * run_in_thread(void * pthread_param)
{
	void ** parameters = (void **)pthread_param;
	global_context_t * global_context = (global_context_t * ) parameters[0];
	thread_context_t * thread_context = (thread_context_t *) parameters[1];
	int task = *((int *)(parameters[2]));
	subread_lock_t * thread_init_lock = (subread_lock_t * ) parameters[3];
	int * ret_value_pointer = (int *)parameters[4];

	if(thread_init_lock)
		subread_lock_release(thread_init_lock);

	switch(task)
	{
		case STEP_VOTING:
			*ret_value_pointer = do_voting(global_context, thread_context);
		break;
		case STEP_ITERATION_ONE:
			*ret_value_pointer = do_iteration_one(global_context, thread_context);
		break;
		case STEP_ITERATION_TWO:
			*ret_value_pointer = do_iteration_two(global_context, thread_context);
		break;
	
	}

	//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_DETAILS, "finished running %d", task);

	return NULL;
}

int run_maybe_threads(global_context_t *global_context, int task)
{
	void * thr_parameters [5];
	int ret_value =0;

	if(task==STEP_VOTING)
		print_in_box(80,0,0, "Map %s...", global_context->input_reads.is_paired_end_reads?"fragments":"reads");
	else if(task == STEP_ITERATION_ONE)
		print_in_box(80,0,0, "Detect indels%s...", global_context->config.is_rna_seq_reads?" and junctions":"");
	else if(task == STEP_ITERATION_TWO)
		print_in_box(80,0,0, "Realign %s...", global_context->input_reads.is_paired_end_reads?"fragments":"reads");

	if(global_context->config.all_threads<2)
	{
		thr_parameters[0] = global_context;
		thr_parameters[1] = NULL;
		thr_parameters[2] = &task;
		thr_parameters[3] = NULL;
		thr_parameters[4] = &ret_value;

		run_in_thread(thr_parameters);
	}
	else
	{
		int current_thread_no ;
		thread_context_t thread_contexts[64];
		int ret_values[64];

		for(current_thread_no = 0 ; current_thread_no < global_context->config.all_threads ; current_thread_no ++)
		{
			thread_contexts[current_thread_no].thread_id = current_thread_no;
			init_indel_thread_contexts(global_context, thread_contexts+current_thread_no, task);
			if(global_context->config.is_rna_seq_reads || global_context->config.do_fusion_detection)
				init_junction_thread_contexts(global_context, thread_contexts+current_thread_no, task);

			relocate_geinputs(global_context, thread_contexts+current_thread_no);

			subread_lock_occupy(&global_context -> thread_initial_lock);
			thr_parameters[0] = global_context;
			thr_parameters[1] = thread_contexts+current_thread_no;
			thr_parameters[2] = &task;
			thr_parameters[3] = (void *)&global_context -> thread_initial_lock;
			thr_parameters[4] = ret_values + current_thread_no;

			pthread_create(&thread_contexts[current_thread_no].thread, NULL,  run_in_thread, &thr_parameters);
		}

		if(task == STEP_VOTING)
		{
			global_context -> processed_reads_in_chunk=0;
		}
		for(current_thread_no = 0 ; current_thread_no < global_context->config.all_threads ; current_thread_no ++)
		{
			pthread_join(thread_contexts[current_thread_no].thread, NULL);
			
			geinput_close(thread_contexts[current_thread_no].ginp1);
			free(thread_contexts[current_thread_no].ginp1);

			if(global_context->input_reads.is_paired_end_reads)
			{
				geinput_close(thread_contexts[current_thread_no].ginp2);
				free(thread_contexts[current_thread_no].ginp2);
			}

			ret_value += *(ret_values + current_thread_no);
			if(ret_value)break;

			if(task == STEP_VOTING)
			{
				//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_DEBUG, "The %d-th thread processed %u reads.", current_thread_no , thread_contexts[current_thread_no].processed_reads_in_chunk);
				global_context -> processed_reads_in_chunk += thread_contexts[current_thread_no].processed_reads_in_chunk;
			}

			finalise_indel_thread(global_context, thread_contexts+current_thread_no, task);
			finalise_junction_thread(global_context, thread_contexts+current_thread_no, task);
		}
	}

	if(CORE_SOFT_BR_CHAR == '\r')
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "");
	return ret_value;
}

void clean_context_after_chunk(global_context_t * context)
{
	memset(context -> chunk_alignment_records , 0 , sizeof(alignment_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);
	memset(context -> big_margin_record  , 0 , context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.big_margin_record_size);
	if(context ->chunk_subjunc_records)
		memset(context ->chunk_subjunc_records , 0 , sizeof(subjunc_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);
}

unsigned int split_read_files(global_context_t * global_context)
{
	unsigned int chunk_reads = global_context->config.reads_per_chunk;
	unsigned int processed_reads = 0;
	unsigned long long * read_position_1;
	unsigned long long * read_position_2 = NULL;
	char * read_line_buf = malloc(3002);
	read_position_1 = (unsigned long long*)malloc(global_context->config.reads_per_chunk * sizeof(long long));
	if(global_context->input_reads.is_paired_end_reads)
		read_position_2 = (unsigned long long*)malloc(global_context->config.reads_per_chunk * sizeof(long long));

	print_in_box(80,0,0, "Scan read files for multi-threaded alignment...");

	if(global_context->config.is_SAM_file_input)
	{
		unsigned long long fhead_pos1;
		unsigned long long fhead_pos2=0;

		while(1)
		{
			char * tok_tmp, * flag;
			if(processed_reads >= chunk_reads || feof(global_context->input_reads.first_read_file.input_fp))
				break;

			fhead_pos1 = ftello(global_context->input_reads.first_read_file.input_fp); 
			if(global_context->input_reads.is_paired_end_reads)
				fhead_pos2 = ftello(global_context->input_reads.second_read_file.input_fp); 

			/*char * is_ret = */fgets(read_line_buf, 3000, global_context->input_reads.first_read_file.input_fp);
			if(global_context->input_reads.is_paired_end_reads)
				fgets(read_line_buf, 3000, global_context->input_reads.second_read_file.input_fp);

			//if(!is_ret) break;

			flag = strtok_r(read_line_buf,"\t",&tok_tmp);
			if(!flag) break;

			flag = strtok_r(NULL,"\t",&tok_tmp);
			if(!flag) break;

			if((atoi(flag) & 0x100) == 0)
			{
				read_position_1[processed_reads] = fhead_pos1;
				if(global_context->input_reads.is_paired_end_reads)
					read_position_2[processed_reads] = fhead_pos2;
				processed_reads++;
			}
			if(global_context->input_reads.is_paired_end_reads)
			{
				fgets(read_line_buf, 3000, global_context->input_reads.first_read_file.input_fp);
				fgets(read_line_buf, 3000, global_context->input_reads.second_read_file.input_fp);
			}
		}
		//printf("PPPP=%llu\n", processed_reads);
		
	}
	else{
		while(1)
		{
			if(processed_reads >= chunk_reads || feof(global_context->input_reads.first_read_file.input_fp))
				break;

			read_position_1[processed_reads] = ftello(global_context->input_reads.first_read_file.input_fp);
			if(global_context->input_reads.is_paired_end_reads)
				read_position_2[processed_reads] = ftello(global_context->input_reads.second_read_file.input_fp);

			processed_reads++;

			geinput_jump_read(&global_context->input_reads.first_read_file);
			if(global_context->input_reads.is_paired_end_reads)
				geinput_jump_read(&global_context->input_reads.second_read_file);
		}
	}

	free(read_line_buf);

	int thread_no;
	for(thread_no = 0; thread_no < global_context->config.all_threads; thread_no++)
	{
		unsigned int my_start_read_no = processed_reads / global_context->config.all_threads * thread_no;
		unsigned int my_reads = (thread_no == global_context->config.all_threads-1)?(processed_reads - my_start_read_no):(processed_reads / global_context->config.all_threads);
		unsigned long long my_first_file_start = read_position_1[my_start_read_no]; 
		unsigned long long my_second_file_start = 0;
		if(global_context->input_reads.is_paired_end_reads)
			my_second_file_start = read_position_2[my_start_read_no];

		global_context -> input_reads.first_file_blocks[thread_no] = my_first_file_start;
		global_context -> input_reads.reads_in_blocks[thread_no] = my_reads;
		global_context -> input_reads.start_read_number_blocks[thread_no] = my_start_read_no;

		if(global_context->input_reads.is_paired_end_reads)
			global_context -> input_reads.second_file_blocks[thread_no] = my_second_file_start;
	}

	free(read_position_1);
	if(read_position_2)
		free(read_position_2);
	return processed_reads;
}

void locate_read_files(global_context_t * global_context, int type)
{
	if(type==SEEK_SET)
	{
		global_context -> current_circle_start_position_file1 = ftello(global_context -> input_reads.first_read_file.input_fp);
		if(global_context ->input_reads.is_paired_end_reads)
			global_context -> current_circle_start_position_file2 = ftello(global_context -> input_reads.second_read_file.input_fp);
	}
	else
	{
		global_context -> current_circle_end_position_file1 = ftello(global_context -> input_reads.first_read_file.input_fp);
		if(global_context ->input_reads.is_paired_end_reads)
			global_context -> current_circle_end_position_file2 = ftello(global_context -> input_reads.second_read_file.input_fp);
	
	}
}
void reward_read_files(global_context_t * global_context, int type)
{
	if(type==SEEK_SET)
	{
		fseeko(global_context -> input_reads.first_read_file.input_fp, global_context -> current_circle_start_position_file1, SEEK_SET);
		if(global_context ->input_reads.is_paired_end_reads)
			fseeko(global_context -> input_reads.second_read_file.input_fp, global_context -> current_circle_start_position_file2, SEEK_SET);
	}
	else
	{
		fseeko(global_context -> input_reads.first_read_file.input_fp, global_context -> current_circle_end_position_file1, SEEK_SET);
		if(global_context ->input_reads.is_paired_end_reads)
			fseeko(global_context -> input_reads.second_read_file.input_fp, global_context -> current_circle_end_position_file2, SEEK_SET);
	
	}
}


int read_chunk_circles(global_context_t *global_context)
{
	int block_no;

//	printf("GINP1 AT %llu\n", ftello(global_context -> input_reads.first_read_file.input_fp));
	
	while(1)
	{
		int ret;

		locate_read_files(global_context, SEEK_SET);
		if(global_context -> config.all_threads>1)
		{
			split_read_files(global_context);
			locate_read_files(global_context, SEEK_END);
			reward_read_files(global_context, SEEK_SET);
		}

		global_context -> current_index = (gehash_t*) malloc(sizeof(gehash_t));
		global_context -> current_value_index = (gene_value_index_t*) malloc(sizeof(gene_value_index_t));
		for(global_context->current_index_block_number = 0; global_context->current_index_block_number < global_context->index_block_number; global_context->current_index_block_number++)
		{
			char tmp_fname[MAX_FILE_NAME_LENGTH];
			sprintf(tmp_fname, "%s.%02d.%c.tab", global_context->config.index_prefix, global_context->current_index_block_number,  global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
			print_in_box(80,0,0, "Load the %d-th index block...",1+ global_context->current_index_block_number);


			if(gehash_load(global_context -> current_index, tmp_fname)) return -1;

			sprintf(tmp_fname, "%s.%02d.%c.array", global_context->config.index_prefix, global_context->current_index_block_number, global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
			if(gvindex_load(global_context -> current_value_index, tmp_fname)) return -1;


			if(global_context->current_index_block_number ==0 && global_context -> all_processed_reads==0)
				global_context->align_start_time = miltime();

			ret = run_maybe_threads(global_context, STEP_VOTING);

			if(global_context -> config.all_threads<2 && global_context->current_index_block_number ==0)
				locate_read_files(global_context, SEEK_END);
			if(global_context->current_index_block_number < global_context->index_block_number -1)
				reward_read_files(global_context, SEEK_SET);

			gehash_destory_fast(global_context -> current_index);
			gvindex_destory(global_context -> current_value_index);
			if(ret) break;
			if(!global_context -> processed_reads_in_chunk) break;
		}

		free(global_context -> current_index);
		free(global_context -> current_value_index);


		//sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "%d reads have been processed in this chunk.", global_context -> processed_reads_in_chunk);


		// after the voting step, all subread index blocks are released and all base index blocks are loaded at once.

		for(block_no = 0; block_no< global_context->index_block_number; block_no++)
		{
			char tmp_fname[MAX_FILE_NAME_LENGTH];
			sprintf(tmp_fname, "%s.%02d.%c.array", global_context->config.index_prefix, block_no,  global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
			if(gvindex_load(&global_context -> all_value_indexes[block_no], tmp_fname)) return -1;
		}

		if(!global_context -> processed_reads_in_chunk)
			// base value indexes loaded in the last circle are not destroyed and are used in writting the indel VCF.
			// the indexes will be destroyed in destroy_global_context
			break;
	
		reward_read_files(global_context, SEEK_SET);
		ret = run_maybe_threads(global_context, STEP_ITERATION_ONE);

		//HashTable * event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID])->event_entry_table;
		//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "There are %ld elements in the indel table before filtering.", event_table ->numOfElements);

		remove_neighbour(global_context);
		//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "There are only %ld elements in the indel table after filtering.", event_table ->numOfElements);

		reward_read_files(global_context, SEEK_SET);
		ret = ret || run_maybe_threads(global_context, STEP_ITERATION_TWO);

	//	printf("IBytes=%d+%d = %d\n", global_context -> all_value_indexes[0].start_base_offset, global_context -> all_value_indexes[0].values_bytes, global_context -> all_value_indexes[0].values [global_context -> all_value_indexes[0].values_bytes-1]);

		//gene_value_index_t * value_index = &global_context->all_value_indexes[0] ;
		
		//printf("=== I=%016llX   B=%016llX\n", (long long)value_index , (long long)value_index -> values);
		if(global_context -> config.is_third_iteration_running)
		{
			reward_read_files(global_context, SEEK_SET);
			ret = ret || do_iteration_three(global_context, NULL);
		}

		if(global_context -> config.report_sam_file)
		{
			reward_read_files(global_context, SEEK_SET);
			print_in_box(80, 0, 0, "%u %s were processed. Save the mapping results for them...", global_context ->processed_reads_in_chunk, global_context -> input_reads.is_paired_end_reads?"fragments":"reads");
			ret = ret || write_chunk_results(global_context);
			if('\r' == CORE_SOFT_BR_CHAR)
				sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"");
			
		}

		reward_read_files(global_context, SEEK_END);

		global_context -> all_processed_reads+= global_context ->processed_reads_in_chunk;

		if(ret) return ret;

		if(global_context -> processed_reads_in_chunk < global_context->config.reads_per_chunk)
			break; 
		else
			// base value indexes loaded in the last circle are not destroyed and are used in writting the indel VCF.
			// the indexes will be destroyed in destroy_global_context
			for(block_no = 0; block_no< global_context->index_block_number; block_no++)
				gvindex_destory(&global_context -> all_value_indexes[block_no]);

		
		clean_context_after_chunk(global_context);
	}

	// load all array index blocks at once.
	if(global_context -> config.is_third_iteration_running)
	{
		/*
		for(block_no = 0; block_no< global_context->index_block_number; block_no++)
		{
			char tmp_fname[MAX_FILE_NAME_LENGTH];
			sprintf(tmp_fname, "%s.%02d.%c.array", global_context->config.index_prefix, block_no,  global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
			if(gvindex_load(&global_context -> all_value_indexes[block_no], tmp_fname)) return -1;
		}
		*/

		finalise_long_insertions(global_context);

		/*
		for(block_no = 0; block_no< global_context->index_block_number; block_no++)
			gvindex_destory(&global_context -> all_value_indexes[block_no]);
		*/
	}
	return 0;
}

void char_strftime(char * tbuf){
	time_t rawtime;
	struct tm * timeinfo;

	time (&rawtime);
	timeinfo = localtime (&rawtime);
	strftime (tbuf,80,"%d-%b-%Y %H:%M:%S",timeinfo);

}

void print_subread_logo()
{
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m ========== %c[0m%c[36m    _____ _    _ ____  _____  ______          _____  ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m =====      %c[0m%c[36m   / ____| |  | |  _ \\|  __ \\|  ____|   /\\   |  __ \\ ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m   =====    %c[0m%c[36m  | (___ | |  | | |_) | |__) | |__     /  \\  | |  | |", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m     ====   %c[0m%c[36m   \\___ \\| |  | |  _ <|  _  /|  __|   / /\\ \\ | |  | |", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m       ==== %c[0m%c[36m   ____) | |__| | |_) | | \\ \\| |____ / ____ \\| |__| |", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m ========== %c[0m%c[36m  |_____/ \\____/|____/|_|  \\_\\______/_/    \\_\\_____/%c[0m", CHAR_ESC, CHAR_ESC, CHAR_ESC, CHAR_ESC);
	#ifdef MAKE_STANDALONE
	char * spaces = "";
	if(strlen(SUBREAD_VERSION) == 8) spaces = "";
	else if(strlen(SUBREAD_VERSION) == 5) spaces = "  ";
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"        %sv%s",spaces,SUBREAD_VERSION);
	#else
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %s",SUBREAD_VERSION);
	#endif
}
int print_configuration(global_context_t * context)
{
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"");
	print_subread_logo();
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"");
	print_in_box(80, 1, 1, context->config.entry_program_name == CORE_PROGRAM_SUBJUNC?"subjunc setting":"subread-align setting");
	print_in_box(80, 0, 1, "");

	if(context->config.is_rna_seq_reads)
	{
		if(context->config.do_fusion_detection)
		{
			print_in_box(80, 0, 0,         "          Function : Read alignment + Junction/Fusion detection%s", context->config.prefer_donor_receptor_junctions?" (RNA-Seq)":" (DNA-Seq)");
		}
		else
			print_in_box(80, 0, 0,         "          Function : Read alignment + Junction detection (RNA-Seq)");
	}
	else
		print_in_box(80, 0, 0,         "          Function : Read alignment");
	print_in_box(80, 0, 0,         "           Threads : %d", context->config.all_threads);
	if( context->config.second_read_file[0])
	{
		print_in_box(80, 0, 0, "      Input file 1 : %s", context->config.first_read_file);
		print_in_box(80, 0, 0, "      Input file 2 : %s", context->config.second_read_file);
	}
	else
		print_in_box(80, 0, 0, "        Input file : %s%s", context->config.first_read_file, context->config.is_SAM_file_input?(context->config.is_BAM_input?" (BAM)":" (SAM)"):"");

	if(context->config.output_prefix [0])
		print_in_box(80, 0, 0, "       Output file : %s (%s)", context->config.output_prefix, context->config.is_BAM_output?"BAM":"SAM");
	else
		print_in_box(80, 0, 0, "     Output method : STDOUT (%s)" , context->config.is_BAM_output?"BAM":"SAM");

	print_in_box(80, 0, 0,         "        Index name : %s", context->config.index_prefix);
	print_in_box(80, 0, 0,         "      Phred offset : %d", (context->config.phred_score_format == FASTQ_PHRED33)?33:64);
	//print_in_box(80, 0, 0,         "        Space type : %s", (context->config.space_type == GENE_SPACE_COLOR)?"color-space":"base-space");
	print_in_box(80, 0, 1, "");
	if( context->config.second_read_file[0])
	{
		print_in_box(80, 0, 0, "   Min read1 votes : %d", context->config.minimum_subread_for_first_read);
		print_in_box(80, 0, 0, "   Min read2 votes : %d", context->config.minimum_subread_for_second_read);
		print_in_box(80, 0, 0, " Max fragment size : %d", context->config.maximum_pair_distance);
		print_in_box(80, 0, 0, " Min fragment size : %d", context->config.minimum_pair_distance);
		print_in_box(80, 0, 1, "");
	}
	else
		print_in_box(80, 0, 0, "         Min votes : %d", context->config.minimum_subread_for_first_read);

	print_in_box(80, 0, 0,         "        Max indels : %d", context->config.max_indel_length);
	print_in_box(80, 0, 0,         " # of Best mapping : %d", context->config.multi_best_reads);
	print_in_box(80, 0, 0,         "    Unique mapping : %s", context->config.report_multi_mapping_reads?"no":"yes");
	print_in_box(80, 0, 0,         "  Hamming distance : %s", context->config.use_hamming_distance_break_ties?"yes":"no");
	print_in_box(80, 0, 0,         "    Quality scores : %s", context->config.use_quality_score_break_ties?"yes":"no");

	if(context->config.max_insertion_at_junctions)
		print_in_box(80, 0, 0,         "Insertions at junc : %d", context->config.max_insertion_at_junctions);

	if(context->config.read_group_id[0])
		print_in_box(80, 0, 0, "   Read group name : %s", context->config.read_group_id);

	print_in_box(80, 0, 1, "");
	print_in_box(80, 2, 1, "http://subread.sourceforge.net/");
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"");

	if(!context->config.first_read_file[0])
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"You have to specify at least one input file in the FASTQ/FASTA/PLAIN format using the '-r' option.\n");
		return -1;
	}

	if(0 && !context->config.output_prefix[0])
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"You have to specify the path of output using the '-o' option.\n");
		return -1;
	}

	if(!context->config.index_prefix[0])
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"You have to specify the prefix of the index files using the '-i' option.\n");
		return -1;
	}
	char tbuf[90];
	char_strftime(tbuf);

	print_in_box(80,1,1,"Running (%s)", tbuf);
	print_in_box(80,0,1,"");


	return 0;
}



int init_paired_votes(global_context_t *context)
{

	if(context -> config.is_rna_seq_reads)
		context -> chunk_subjunc_records = malloc(sizeof(subjunc_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);
	else	context -> chunk_subjunc_records = NULL;
	context -> chunk_alignment_records = malloc(sizeof(alignment_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);


	if(!context -> chunk_alignment_records)
	{	
		return 1;
	}

	context -> big_margin_record = malloc( (context->input_reads.is_paired_end_reads?2:1) * context -> config.big_margin_record_size * context ->config.reads_per_chunk);

	memset(context ->big_margin_record  , 0 , context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context -> config.big_margin_record_size);
	memset(context ->chunk_alignment_records , 0 , sizeof(alignment_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);

		//fprintf(stderr, "MALLOC=%llu = %d * %d * %d \n", sizeof(alignment_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads, sizeof(alignment_result_t), context ->config.reads_per_chunk, context->config.multi_best_reads);
		//sleep(10000);

	if(context -> chunk_subjunc_records)
		memset(context ->chunk_subjunc_records , 0 , sizeof(subjunc_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);


	return 0;
}

void write_sam_headers(global_context_t * context)
{
	if(context -> config.is_BAM_output)
	{
		SamBam_writer_add_header(context -> output_bam_writer,"@HD\tVN:1.0\tSO:unsorted", 0);
		int xk1;
		unsigned int last_offset = 0;
		char obuf[300];
		for(xk1=0; xk1< context->chromosome_table.total_offsets; xk1++)
		{
			SamBam_writer_add_chromosome(context -> output_bam_writer, context->chromosome_table.read_names+ xk1 * MAX_CHROMOSOME_NAME_LEN, context->chromosome_table.read_offsets[xk1] - last_offset+16, 1);
			last_offset = context->chromosome_table.read_offsets[xk1];
		}


		if(context->config.read_group_id[0])
		{
			snprintf(obuf,299, "@RG\tID:%s%s",context->config.read_group_id, context->config.read_group_txt);
			SamBam_writer_add_header(context -> output_bam_writer,obuf, 0);
		}
		snprintf(obuf,299, "@PG\tID:subread\tPN:subread\tVN:%s", SUBREAD_VERSION);
		SamBam_writer_add_header(context -> output_bam_writer,obuf, 0);
	}
	else
	{
		sambamout_fprintf(context -> output_sam_fp, "@HD\tVN:1.0\tSO:unsorted\n");
		int xk1;
		unsigned int last_offset = 0;
		for(xk1=0; xk1< context->chromosome_table.total_offsets; xk1++)
		{
			sambamout_fprintf(context -> output_sam_fp, "@SQ\tSN:%s\tLN:%u\n", context->chromosome_table.read_names+ xk1 * MAX_CHROMOSOME_NAME_LEN, context->chromosome_table.read_offsets[xk1] - last_offset+16);
			last_offset = context->chromosome_table.read_offsets[xk1];
		}

		if(context->config.read_group_id[0])
			sambamout_fprintf(context -> output_sam_fp, "@RG\tID:%s%s\n",context->config.read_group_id, context->config.read_group_txt);
		sambamout_fprintf(context -> output_sam_fp, "@PG\tID:subread\tPN:subread\tVN:%s\n", SUBREAD_VERSION);
		
	}
}

int load_global_context(global_context_t * context)
{
	char tmp_fname [MAX_FILE_NAME_LENGTH];
	int guess_phred_format = -1 ;
	

	if(core_geinput_open(context, &context->input_reads.first_read_file, 1,1))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable to open '%s' as input. Please check if it exists, you have the permission to read it, and it is in the correct format.\n", context->config.first_read_file);
		return -1;
	}

	context->config.space_type = context->input_reads.first_read_file.space_type;
	print_in_box(89,0,0,"The input file contains %c[36m%s%c[0m space reads.", CHAR_ESC, context->config.space_type == GENE_SPACE_COLOR?"color":"base", CHAR_ESC);
	if(context->config.space_type == GENE_SPACE_COLOR && context->config.is_BAM_output && !context->config.convert_color_to_base)
	{
		print_in_box(80,0,0,"The color-space bases will be converted into base space in the BAM output.");
		context->config.convert_color_to_base=1;
	}
	else if(context->config.space_type == GENE_SPACE_BASE && context->config.convert_color_to_base)
	{
		print_in_box(80,0,0,"The reads will not be converted into base space.");
		context->config.convert_color_to_base=0;
	}

	if(context->input_reads.is_paired_end_reads)
	{
		if(core_geinput_open(context, &context->input_reads.second_read_file, 2,1))
		{
			sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable to open '%s' as input. Please check if it exists, you have the permission to read it, and it is in the correct format.\n", context->config.second_read_file);
			return -1;
		}
	}


	if(context->input_reads.is_paired_end_reads)
		context->config.reads_per_chunk = 7*1024*1024*min(40,max(1,context->config.memory_use_multiplex));
	else
		context->config.reads_per_chunk = 14*1024*1024*min(40,max(1,context->config.memory_use_multiplex));

	struct stat ginp1_stat;
	stat(context->config.first_read_file , &ginp1_stat);
	context->input_reads.first_read_file_size = ginp1_stat.st_size;


	context -> input_reads.avg_read_length = guess_reads_density_format(context->config.first_read_file , context->config.is_SAM_file_input + context->input_reads.is_paired_end_reads, &guess_phred_format);
	if(context -> input_reads.avg_read_length<0 )context -> input_reads.avg_read_length = 250;


	if(context->config.report_sam_file && context -> config.output_prefix[0])
	{
		// ====== open output files ======
		// Only the sam file is opened here; other files like bed, indel and etc are opened in init_modules()
		sprintf(tmp_fname,"%s", context->config.output_prefix);

		if(context -> config.is_BAM_output)
		{
			context -> output_bam_writer = malloc(sizeof(SamBam_Writer));
			SamBam_writer_create(context -> output_bam_writer , tmp_fname);
			context -> output_sam_fp = NULL;
		}
		else
		{
			context -> output_sam_fp = f_subr_open(tmp_fname,"wb");
			context -> output_bam_writer = NULL;
		}
		if((!context -> output_bam_writer) && (!context->output_sam_fp))
		{
			sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable to open '%s' to write. Please check if the path exists and you have the permission to create/write this file", tmp_fname);
			return -1;
		}
	}
	else
	{
		if(context -> config.is_BAM_output)
		{
			context -> output_bam_writer = malloc(sizeof(SamBam_Writer));
			SamBam_writer_create(context -> output_bam_writer ,NULL);
		}
		context->output_sam_fp = NULL;
	}
	
	// ====== check index files, count blocks and load chro table ======
	sprintf(tmp_fname, "%s.reads", context->config.index_prefix);
	if(!does_file_exist(tmp_fname))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable top open index '%s'. Please make sure that the correct prefix is specified and you have the permission to read these files. For example, if there are files '/opt/my_index.reads', '/opt/my_index.files' and etc, the index prefix should be specified as '/opt/my_index' without any suffix. \n", context->config.index_prefix);
		return -1;
	}

	if(context->config.space_type == GENE_SPACE_COLOR)
		sprintf(tmp_fname, "%s.00.c.tab", context->config.index_prefix);
	else
		sprintf(tmp_fname, "%s.00.b.tab", context->config.index_prefix);
	if(!does_file_exist(tmp_fname))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Your reads are in the %s space but the index was not built in the same space. Unable to precess the reads.\n", context->config.space_type == GENE_SPACE_COLOR?"color":"base");
		return -1;
	}


	context->index_block_number = 0; 
	while(1)
	{
		sprintf(tmp_fname, "%s.%02d.%c.tab", context->config.index_prefix, context->index_block_number, context->config.space_type == GENE_SPACE_COLOR?'c':'b');
		if(!does_file_exist(tmp_fname))break;
		context->index_block_number ++;
		if(context->index_block_number>=2 && context->config.max_indel_length > 16)
		{
			print_in_box(80,0,0,"ERROR You cannot use multi-block index for very-long indel detection!");
			print_in_box(80,0,0,"Please set the maximum indel length <= 16.");
			return -1;
		}
	}

	context->current_index_block_number = 0;
	load_offsets(&context->chromosome_table, context->config.index_prefix);

	if(context->config.report_sam_file)
		write_sam_headers(context);

	// ====== init other variables ======
	if(context -> config.all_threads>1)
		subread_init_lock(&context -> thread_initial_lock);

	if(init_paired_votes(context))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Cannot initialise the voting space. You need at least 2GB of empty physical memory to run this program.\n");
		return 1;
	}
	context->all_processed_reads = 0;
	context->all_mapped_reads = 0;
	context->all_correct_PE_reads = 0;
	context->all_junctions = 0;
	context->all_fusions = 0;
	context->all_indels = 0;
	sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "load_global_context: finished");

	memset( context->all_value_indexes , 0 , 100 * sizeof(gene_value_index_t));

	return 0;
}

int init_modules(global_context_t * context)
{
	sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "init_modules: begin");
	int ret = init_indel_tables(context);
	if(context->config.is_rna_seq_reads || context->config.do_fusion_detection)
		ret = ret || init_junction_tables(context);

	sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "init_modules: finished: %d",ret);
	return ret;
}

int destroy_modules(global_context_t * context)
{
	destroy_indel_module(context);
	if(context->config.is_rna_seq_reads)
		destroy_junction_tables(context);
	return 0;
}

int destroy_global_context(global_context_t * context)
{
	int xk1, block_no;

	for(block_no = 0; block_no< context->index_block_number; block_no++)
		gvindex_destory(&context -> all_value_indexes[block_no]);

	if(context->output_sam_fp)
		fclose(context->output_sam_fp);
	if(context->output_bam_writer)
	{
		SamBam_writer_close(context->output_bam_writer);
		free(context->output_bam_writer);
		context->output_bam_writer=NULL;
	}
	free(context->chunk_alignment_records);
	free(context->big_margin_record);
	if(context->chunk_subjunc_records)
		free(context->chunk_subjunc_records);
	
	for(xk1=0; xk1<5; xk1++)
		if(context->module_contexts[xk1])free(context->module_contexts[xk1]);
	geinput_close(&context -> input_reads.first_read_file);
	if(context->input_reads.is_paired_end_reads) geinput_close(&context -> input_reads.second_read_file);
	destroy_offsets(&context->chromosome_table);
	if((context -> will_remove_input_file & 1) && (memcmp(context ->config.first_read_file, "./core-temp", 11) == 0)) unlink(context ->config.first_read_file);
	if((context -> will_remove_input_file & 2) && (memcmp(context ->config.second_read_file, "./core-temp", 11) == 0)) unlink(context ->config.second_read_file);

	return 0;
}


int write_bincigar_part(char * bincigar, int chropt, unsigned int optlen, int bincigar_len)
{
	int binopt, binbytes, x1;

	if(optlen<1) return -1;

	if(optlen < 4) binbytes=1;
	else if(optlen < 1024) binbytes=2;
	else if(optlen < 262144) binbytes=3;
	else if(optlen < 67108864) binbytes=4;
	else binbytes=5;

	if(bincigar_len<binbytes) return -1; 

	switch(chropt)
	{
		case 'S':
			binopt = CORE_CIGAR_OPT_S;
			break;
		case 'M':
			binopt = CORE_CIGAR_OPT_M;
			break;
		case 'I':
			binopt = CORE_CIGAR_OPT_I;
			break;
		case 'D':
			binopt = CORE_CIGAR_OPT_D;
			break;
		case 'N':
			binopt = CORE_CIGAR_OPT_N;
			break;
		case 'n':
			binopt = CORE_CIGAR_OPT_NN;
			break;
		case 'B':
			binopt = CORE_CIGAR_OPT_B;
			break;
		case 'b':
			binopt = CORE_CIGAR_OPT_BB;
			break;
		default:
			return -1;
	}

	bincigar[0]=binopt | (binbytes << 3) | ((optlen & 3)<< 6);
	optlen >>= 2;
	for(x1=1;x1<binbytes; x1++)
	{
		bincigar[x1] = optlen&0xff;
		optlen>>=8;
	}

	return binbytes;
}

// function returns the actual length of bincigar, or -1 if anything is wrong, e.g., bincigar_len is too short or unrecognized operations.
int cigar2bincigar(char *cigar, char *bincigar, int bincigar_len)
{
	int xk1=0;
	unsigned int tmpv=0, bincigar_cursor=0;
	while(1)
	{
		int nch = cigar[xk1];
		if(!nch) break;
		xk1++;

		if(isdigit(nch)) tmpv=tmpv*10+(nch-'0');
		else
		{
			int bincigar_sec_len = write_bincigar_part(bincigar+bincigar_cursor, nch, tmpv, bincigar_len-bincigar_cursor);
			if(bincigar_sec_len<0){
				bincigar[0]=0;
				return -1;
			}
			bincigar_cursor += bincigar_sec_len;
			tmpv=0;
		}
	}

	if(bincigar_cursor<bincigar_len) bincigar[bincigar_cursor] = 0;

	//printf("%s : BL=%d\n", cigar, bincigar_cursor);

	return bincigar_cursor;
}


int write_cigar_part(char *bincigar, char *cigar, int cigar_len , int * bincigar_move)
{
	int binbytes, x1, binopt, charopt;
	unsigned int tmpv = 0;
	char sec_buf[13];

	binbytes = 7& (bincigar[0] >> 3);
	binopt = 7 & bincigar[0];

	switch(binopt)
	{
		case CORE_CIGAR_OPT_D:
			charopt='D'; 
			break;
		case CORE_CIGAR_OPT_I:
			charopt='I'; 
			break;
		case CORE_CIGAR_OPT_M:
			charopt='M'; 
			break;
		case CORE_CIGAR_OPT_S:
			charopt='S'; 
			break;
		case CORE_CIGAR_OPT_B:
			charopt='B'; 
			break;
		case CORE_CIGAR_OPT_BB:
			charopt='b'; 
			break;
		case CORE_CIGAR_OPT_N:
			charopt='N'; 
			break;
		default:
			charopt='n'; 
			break;
	}

	tmpv = (bincigar[0]>>6) & 3;
	for(x1 = 1; x1 < binbytes; x1++)
	{
		unsigned int dtmpv = 0xff & bincigar[x1];
		dtmpv <<= (x1*8 - 6);
		tmpv += dtmpv; 
	}

	int added_len = sprintf(sec_buf, "%u%c", tmpv, charopt);
	if(added_len > cigar_len)
		return -1;
	memcpy(cigar, sec_buf, added_len);
	(*bincigar_move) = binbytes;

	return added_len;
}

int bincigar2cigar(char * cigar, int cigar_len, char * bincigar, int bincigar_max_len, int read_len)
{
	int cigar_cursor = 0, bincigar_cursor = 0;
	while(1)
	{
		int bincigar_move = 0;
		int cigar_sec_len = write_cigar_part(bincigar + bincigar_cursor, cigar+cigar_cursor, cigar_len-cigar_cursor-1, &bincigar_move);
		if(cigar_sec_len<0){
			sprintf(cigar,"%dM", read_len);
			return -1;
		}
		//printf("NPC=%s\n", cigar);
		cigar_cursor += cigar_sec_len;
		bincigar_cursor += bincigar_move;
		if(bincigar_cursor>=bincigar_max_len) break;
		if(bincigar[bincigar_cursor] == 0) break;
	}
	cigar[cigar_cursor] = 0;

	return cigar_cursor;
}

int term_strncpy(char * dst, char * src, int max_dst_mem)
{
	int i;

	for(i=0; i<max_dst_mem; i++)
	{
		if(!src[i]) break;
		dst[i]=src[i];
		if(i == max_dst_mem-1)
			SUBREADprintf("String out of memory limit: '%s'\n", src);
	}
	if(i >= max_dst_mem) i = max_dst_mem-1;
	dst[i] = 0;

	return 0;
}


// This assumes the first part of Cigar has differet strandness to the main part of the cigar.
// Pos is the LAST WANTED BASE location before the first strand jump (split by 'b' or 'n').
// The first base in the read actually has a larger coordinate than Pos. 
// new_cigar has to be at least 100 bytes.
unsigned int reverse_cigar(unsigned int pos, char * cigar, char * new_cigar)
{
	int cigar_cursor = 0;
	new_cigar[0]=0;
	unsigned int tmpi=0;
	int last_piece_end = 0;
	int last_sec_start = 0;
	unsigned int chro_pos = pos, this_section_start = pos, ret = pos;
	int is_positive_dir = 0;
	int read_cursor = 0;
	int section_no = 0;

	for(cigar_cursor = 0 ;  ; cigar_cursor++)
	{
		if( cigar [cigar_cursor] == 'n' ||  cigar [cigar_cursor] == 'b' ||  cigar [cigar_cursor] == 0)
		{
			int xk1, jmlen=0, nclen=strlen(new_cigar);
			char jump_mode [13];

			if(cigar [cigar_cursor] !=0)
			{
				sprintf(jump_mode, "%u%c", tmpi,  cigar [cigar_cursor] == 'b'?'n':'b');
				jmlen = strlen(jump_mode);
			}

			for(xk1=nclen-1;xk1>=0; xk1--)
				new_cigar[ xk1 +  last_piece_end + jmlen - last_sec_start ] = new_cigar[ xk1 ];
			new_cigar [nclen + jmlen + last_piece_end - last_sec_start ] = 0;

			memcpy(new_cigar , jump_mode, jmlen);
			memcpy(new_cigar + jmlen , cigar + last_sec_start, last_piece_end - last_sec_start);

			last_sec_start = cigar_cursor+1;

			if(is_positive_dir && cigar [cigar_cursor] !=0)
			{
				if(cigar [cigar_cursor] == 'b') chro_pos -= tmpi - read_cursor - 1;
				else	chro_pos += tmpi - read_cursor - 1;
			}
			if((!is_positive_dir) && cigar [cigar_cursor] !=0)
			{
				if(cigar [cigar_cursor] == 'b') chro_pos = this_section_start - tmpi - read_cursor - 1;
				else	chro_pos = this_section_start + tmpi - read_cursor - 1;
			}

			this_section_start = chro_pos;

			if(section_no == 0)
				ret = chro_pos;

			is_positive_dir = ! is_positive_dir;
			section_no++;
			tmpi=0;
		}
		else if(isalpha(cigar [cigar_cursor]))
		{
			if(cigar [cigar_cursor]=='M' || cigar [cigar_cursor] == 'S')
				read_cursor += tmpi;
			tmpi=0;
			last_piece_end = cigar_cursor+1;
		}
		else tmpi = tmpi*10 + (cigar [cigar_cursor] - '0');

		if(cigar [cigar_cursor] == 0)break;
	}

	//printf("REV CIGAR: %s  =>  %s\n", cigar, new_cigar);
	return ret;
}

int chimeric_cigar_parts(global_context_t * global_context, unsigned int sel_pos, int is_first_section_negative_strand, int is_first_section_reversed, char * in_cigar, unsigned int * out_poses, char ** out_cigars, char * out_strands, int read_len, short * perfect_lens)
{
	unsigned int current_perfect_map_start = sel_pos;
	int current_perfect_section_no = 0;
	int current_perfect_cursor = sel_pos;
	int is_reversed = is_first_section_reversed;
	int is_negative = is_first_section_negative_strand;
	int read_cursor = 0;
	int out_cigar_writer_ptr = 0;
	unsigned int tmpi = 0;

	short perfect_len = 0;

	int cigar_cursor;

	out_poses[0] = current_perfect_map_start;
	out_strands[0] = is_negative?'-':'+';
	char main_piece_strand = (is_first_section_negative_strand == is_first_section_reversed)?'+':'-';

	for(cigar_cursor=0;;cigar_cursor++)
	{
		char ncch = in_cigar[cigar_cursor];
		int is_chimeric_section_end = 0;

		if(!ncch){
			perfect_lens [current_perfect_section_no] = perfect_len ;
			current_perfect_section_no++;
			break;
		}

		if(toupper(ncch)=='N'||toupper(ncch)=='B')
		{

			unsigned int jummped_location;
			int is_chro_jump = 0, is_long_jump = 0;

			if(is_reversed)
			{
				if(toupper(ncch)=='N')
					jummped_location = current_perfect_map_start - 1 + tmpi;
				else
					jummped_location = current_perfect_map_start - 1 - tmpi;
			}
			else
			{
				if(toupper(ncch)=='N')
					jummped_location = current_perfect_cursor + tmpi;
				else
					jummped_location = current_perfect_cursor - tmpi;

			}

			if(ncch == 'N')
			{
				char * curr_chr, * new_chr;
				unsigned int curr_offset, new_offset;
				locate_gene_position_max(current_perfect_cursor, &global_context -> chromosome_table, & curr_chr, & curr_offset, 1);
				locate_gene_position_max(jummped_location      , &global_context -> chromosome_table, &  new_chr, &  new_offset, 1);
				is_chro_jump = (curr_chr != new_chr);

				long long int dist = current_perfect_cursor;
				dist -= jummped_location;
				if(abs(dist) >= 134217728)
					is_long_jump = 1;
				// A long jump is the jump longer than 2^27.
				// Picard does not like it!!
			}

			if(is_chro_jump || islower(ncch) || ncch == 'B' || is_long_jump)
			{
				current_perfect_cursor = jummped_location;

				if(islower(ncch)){
					is_reversed = !is_reversed;
					is_negative = !is_negative;
				}

				current_perfect_map_start = current_perfect_cursor;
				tmpi = 0;
				if(read_cursor<read_len)
					sprintf(out_cigars[current_perfect_section_no] + out_cigar_writer_ptr,"%dS", read_len - read_cursor);

				perfect_lens [current_perfect_section_no] = perfect_len ;
				perfect_len = 0;

				current_perfect_section_no++;
				if(current_perfect_section_no>CIGAR_PERFECT_SECTIONS)break;

				out_poses[current_perfect_section_no] = current_perfect_map_start - read_cursor;
				out_strands[current_perfect_section_no] = is_negative?'-':'+';
				out_cigar_writer_ptr = sprintf(out_cigars[current_perfect_section_no],"%dS", read_cursor);
				is_chimeric_section_end  = 1;
			}
		}

		if(!is_chimeric_section_end)
		{
			if(isalpha(ncch))
			{
				out_cigar_writer_ptr+=sprintf(out_cigars[current_perfect_section_no]+out_cigar_writer_ptr, "%u%c", tmpi, ncch);
			}
			if(ncch == 'M'|| ncch == 'S')
			{
				read_cursor += tmpi;
				if(ncch == 'M')
					perfect_len += tmpi;
				if(!is_reversed)
					current_perfect_cursor += tmpi;
				tmpi = 0;
			}
			else if(ncch == 'D' || ncch == 'N')
			{
				if(!is_reversed)
					current_perfect_cursor += tmpi;
				tmpi = 0;
			}
			else if(ncch == 'I')
			{
				read_cursor += tmpi;
				tmpi = 0;
			}
			else if(isdigit(ncch))
				tmpi = tmpi*10+(ncch-'0');
		}
	}

	int xk1 = 0, best_match = -9999;

	for(xk1=0; xk1<current_perfect_section_no;xk1++)
	{
		if(main_piece_strand == out_strands[xk1] && perfect_lens[xk1]>perfect_lens[best_match])
			best_match = xk1;
	}

	if(best_match>0)
	{
		unsigned int tmpv;
		char cigar_sec[100];
		tmpv = out_poses[0];
		out_poses[0]=out_poses[best_match];
		out_poses[best_match] = tmpv;

		tmpv = out_strands[0];
		out_strands[0] = out_strands[best_match];
		out_strands[best_match] = tmpv;

		strcpy(cigar_sec, out_cigars[0]);
		strcpy(out_cigars[0], out_cigars[best_match]);
		strcpy(out_cigars[best_match] , cigar_sec);

		tmpv = perfect_lens[0];
		perfect_lens[0] = perfect_lens[best_match];
		perfect_lens[best_match] = tmpv;
	}

	return current_perfect_section_no;
}

void quick_sort_run(void * arr, int spot_low,int spot_high, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r));

void quick_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r))
{
	quick_sort_run(arr, 0, arr_size-1, compare, exchange);
}
 
 
void quick_sort_run(void * arr, int spot_low,int spot_high, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r))
{
	int pivot,j,i;

	if(spot_high-spot_low<1) return;
	pivot = spot_low;
	i = spot_low;
	j = spot_high;

	while(i<=j)
	{
		if(compare(arr, i, pivot) <0)
		{
			i++;
			continue;
		}

		if(compare(arr, j, pivot)>0)
		{
			j--;
			continue;
		}

		if(i!=j)
			exchange(arr,i,j);
		i++;
		j--;
	}

	quick_sort_run(arr, spot_low, j, compare, exchange);
	quick_sort_run(arr, i, spot_high, compare, exchange);
	
}



void merge_sort_run(void * arr, int start, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	if(items > 11)
	{
		int half_point = items/2;
		merge_sort_run(arr, start, half_point, compare, exchange, merge);
		merge_sort_run(arr, start + half_point, items - half_point, compare, exchange, merge);
		merge(arr, start, half_point, items - half_point);
	}
	else
	{
		int i, j;
		for(i=start; i< start + items - 1; i++)
		{
			int min_j = i;
			for(j=i + 1; j< start + items; j++)
			{
				if(compare(arr, min_j, j) > 0)
					min_j = j;
			}
			if(i!=min_j)
				exchange(arr, i, min_j);
		}
	}
}
void merge_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	merge_sort_run(arr, 0, arr_size, compare, exchange, merge);
}
