#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>

#include "subread.h"
#include "input-files.h"
#include "core.h"


static struct option long_options[] =
{
	{"index",  required_argument, 0, 'i'},
	{"read",  required_argument, 0, 'r'},
	{"read2",  required_argument, 0, 'R'},
	{"output",  required_argument, 0, 'o'},
	{"subreads",  required_argument, 0, 'n'},
	{"singleSAM",  required_argument, 0, '1'},
	{"pairedSAM",  required_argument, 0, '2'},
	{"threads",  required_argument, 0, 't'},
	{"indel",  required_argument, 0, 'I'},
	{"phred",  required_argument, 0, 'P'},
	{"mindist",  required_argument, 0, 'd'},
	{"maxdist",  required_argument, 0, 'D'},
	{"order",  required_argument, 0, 'S'},
	{"trim5", required_argument, 0, '5'},
	{"trim3", required_argument, 0, '3'},
	{"color-convert",  no_argument, 0, 'b'},
	{"junctionIns", required_argument, 0, 0},
	{"rg",  required_argument, 0, 0},
	{"rg-id",  required_argument, 0, 0},
	{"BAMoutput", no_argument, 0, 0},
	{"BAMinput", no_argument, 0, 0},
	{"SAMinput", no_argument, 0, 0},
	{0, 0, 0, 0}
};

void print_usage_core_subjunc()
{

	SUBREADprintf("\nVersion %s\n\n", SUBREAD_VERSION);
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs(" ./subjunc [options] -i <index_name> -r <input> -o <output>");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("");
	SUBREADputs("    -i --index     <index>  base name of the index.");
	SUBREADputs("");
	SUBREADputs("    -r --read      <input>  name of the input file(FASTQ/FASTA format). Both ");
	SUBREADputs("                            base-space and color-space read data are supported. ");
	SUBREADputs("                            For paired-end reads, this gives the first read file");
	SUBREADputs("                            and the other read file should be specified using");
	SUBREADputs("                            the -R option.");
	SUBREADputs("");
	SUBREADputs("Optional general arguments:");
	SUBREADputs("");
	SUBREADputs("    -o --output    <output> name of the output file(SAM format by default). If");
	SUBREADputs("                            not provided, mapping results will be output to the");
	SUBREADputs("                            standard output (stdout).");
	SUBREADputs("");
	SUBREADputs("    -n --subreads  <int>    number of selected subreads, 14 by default.");
	SUBREADputs("");
	SUBREADputs("    -T --threads   <int>    number of threads/CPUs used, 1 by default.");
	SUBREADputs("");
	SUBREADputs("    -I --indel     <int>    number of INDEL bases allowed, 5 by default.");
	SUBREADputs("");
	SUBREADputs("    -P --phred     <3:6>    the format of Phred scores used in input files, '3'");
	SUBREADputs("                            for phred+33 and '6' for phred+64. '3' by default.");
	SUBREADputs("");
	SUBREADputs("    -b --color-convert      convert color-space read bases to base-space read");
	SUBREADputs("                            bases in the mapping output. Note that the mapping");
	SUBREADputs("                            itself will still be performed at color-space.");
	SUBREADputs("   ");
	SUBREADputs("    -v                      displaying the version number.");
	SUBREADputs("                                 ");
	SUBREADputs("       --trim5     <int>    trim off <int> number of bases from 5' end of each");
	SUBREADputs("                            read. 0 by default.");
	SUBREADputs("                                 ");
	SUBREADputs("       --trim3     <int>    trim off <int> number of bases from 3' end of each");
	SUBREADputs("                            read. 0 by default.");
	SUBREADputs("                                 ");
	SUBREADputs("       --rg-id     <string> specify the read group ID. It will be added to the");
	SUBREADputs("                            read group header field and also be added to each");
	SUBREADputs("                            read. ");
	SUBREADputs("                                 ");
	SUBREADputs("       --rg        <string> add a <tag:value> to  the read group (RG) header. ");
	SUBREADputs("                                 ");
	SUBREADputs("       --SAMinput           the input read data are in SAM format.");
	SUBREADputs("                                 ");
	SUBREADputs("       --BAMinput           the input read data are in BAM format.");
	SUBREADputs("                                 ");
	SUBREADputs("       --BAMoutput          mapping results will be saved into a BAM format file");
	SUBREADputs("                            instead of a SAM format file.");
	SUBREADputs("");
	SUBREADputs("");
	SUBREADputs("Optional arguments for paired-end reads:");
	SUBREADputs("");
	SUBREADputs("    -R --read2     <input>  name of the second input file from paired-end data. ");
	SUBREADputs("                            The program will then be switched to paired-end read");
	SUBREADputs("                            mapping mode.");
	SUBREADputs("");
	SUBREADputs("    -d --mindist   <int>    minimum fragment/template length, 50bp by default.");
	SUBREADputs("");
	SUBREADputs("    -D --maxdist   <int>    maximum fragment/template length, 600bp by default.");
	SUBREADputs("");
	SUBREADputs("    -S --order     <ff:fr:rf>  specifying if the first/second reads are forward");
	SUBREADputs("                            or reversed, 'fr' by default");
	SUBREADputs("");
	SUBREADputs("For more information about these arguments, please refer to the User Manual.");
	SUBREADputs("");


}

int parse_opts_subjunc(int argc , char ** argv, global_context_t * global_context)
{
	int c;
	int option_index = 0;	

	optind = 1;
	opterr = 1;
	optopt = 63;

	global_context->config.max_mismatch_exonic_reads = 10;
	global_context->config.max_mismatch_junction_reads = 1;
	global_context->config.ambiguous_mapping_tolerance = 39;
	global_context->config.extending_search_indels = 0;
	global_context->config.do_fusion_detection =0;
	global_context->config.use_dynamic_programming_indel = 0;

	global_context->config.is_rna_seq_reads = 1;
	global_context->config.total_subreads = 14;
	global_context->config.minimum_subread_for_first_read =1;
	global_context->config.minimum_subread_for_second_read = 1;
	global_context->config.high_quality_base_threshold = 990000;
	global_context->config.do_big_margin_filtering_for_junctions = 1;
	global_context->config.report_no_unpaired_reads = 0;
	global_context->config.limited_tree_scan = 1;
	global_context->config.use_hamming_distance_in_exon = 0;
	global_context->config.big_margin_record_size = 24;

	if(argc<2)
	{
		print_usage_core_subjunc();
		return -1;
	}
	while ((c = getopt_long (argc, argv, "ExsJ1:2:S:L:AHd:D:n:m:p:P:R:r:i:l:o:G:T:Q:I:t:B:bQ:FcuUfM3:5:9:?", long_options, &option_index)) != -1)
	{
		switch(c)
		{


			case '3':
				global_context->config.read_trim_3 = atoi(optarg); 
				break;
			case '5':
				global_context->config.read_trim_5 = atoi(optarg); 
				break;

			case '1':
			case '2':
			case 'r':
				strncpy(global_context->config.first_read_file, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'J':
				global_context->config.show_soft_cliping = 1;
				break;
			case 'Q':
				global_context->config.multi_best_reads = atoi(optarg); 
				if(global_context->config.multi_best_reads <1)
					global_context->config.multi_best_reads=1;
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
				global_context->config.convert_color_to_base = 1;
				break;
			case 'D':
				global_context->config.maximum_pair_distance = atoi(optarg);
				break;
			case 'd':
				global_context->config.minimum_pair_distance = atoi(optarg);
				break;
			case 'n':
				global_context->config.total_subreads = atoi(optarg);
				global_context->config.total_subreads = min(31,global_context->config.total_subreads);
				break;
			case 'm':
				global_context->config.minimum_subread_for_first_read = atoi(optarg);
				break;
			case 'T':
				global_context->config.all_threads = atoi(optarg);
				if(global_context->config.all_threads <1) global_context->config.all_threads = 1;
				if(global_context->config.all_threads >32) global_context->config.all_threads = 32;

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
					global_context->config.flanking_subread_indel_mismatch = 0;

					global_context->config.is_third_iteration_running = 1;
					global_context->config.max_mismatch_exonic_reads = 2;
					global_context->config.max_mismatch_junction_reads = 2;
					global_context->config.total_subreads = 28;
					global_context->config.minimum_subread_for_first_read = 3;
					global_context->config.minimum_subread_for_second_read = 1;
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
				if(strcmp("rg-id", long_options[option_index].name)==0) 
				{
					strcpy(global_context->config.read_group_id, optarg);
				}
				else if(strcmp("rg", long_options[option_index].name)==0) 
				{
					strcat(global_context->config.read_group_txt, "\t");
					strcat(global_context->config.read_group_txt, optarg);
				}
				else if(strcmp("BAMoutput", long_options[option_index].name)==0) 
				{
					global_context->config.is_BAM_output = 1;
				}
				else if(strcmp("BAMinput", long_options[option_index].name)==0) 
				{
					global_context->config.is_BAM_input = 1;
					global_context->config.is_SAM_file_input = 1;
				}
				else if(strcmp("SAMinput", long_options[option_index].name)==0) 
				{
					global_context->config.is_BAM_input = 0;
					global_context->config.is_SAM_file_input = 1;
				}
				else if(strcmp("junctionIns", long_options[option_index].name)==0)
				{
					global_context->config.check_donor_at_junctions=0;
					global_context->config.limited_tree_scan = 0;
					global_context->config.max_insertion_at_junctions = atoi(optarg);
				}
				break;


			case '?':
			default:
				print_usage_core_subjunc();
				return -1 ;
		}
	}

	if(global_context->config.is_SAM_file_input) global_context->config.phred_score_format = FASTQ_PHRED33;


	return 0;
}





#if defined MAKE_STANDALONE
int main(int argc , char ** argv)
{
#elif defined RUNNING_ENV_JAVA
int subread_subjunc_main(int argc , char ** argv)
{
#else
int main_junction(int argc , char ** argv)
{
#endif
	return core_main(argc, argv, parse_opts_subjunc);
}

