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
	{"minmatch",  required_argument, 0, 'm'},
	{"minmatch2",  required_argument, 0, 'p'},
	{"singleSAM",  required_argument, 0, '1'},
	{"pairedSAM",  required_argument, 0, '2'},
	{"threads",  required_argument, 0, 'T'},
	{"indel",  required_argument, 0, 'I'},
	{"phred",  required_argument, 0, 'P'},
	{"mindist",  required_argument, 0, 'd'},
	{"maxdist",  required_argument, 0, 'D'},
	{"order",  required_argument, 0, 'S'},
	{"trim5", required_argument, 0, '5'},
	{"trim3", required_argument, 0, '3'},
	{"color-convert",  no_argument, 0, 'b'},
	{"junctionIns", required_argument, 0, 0},
	{"multi",  required_argument, 0, 'B'},
	{"rg",  required_argument, 0, 0},
	{"rg-id",  required_argument, 0, 0},
	{"gzFASTQinput",  no_argument, 0, 0},
	{"unique",  no_argument, 0, 'u'},
	{"SAMoutput", no_argument, 0, 0},
	{"BAMinput", no_argument, 0, 0},
	{"SAMinput", no_argument, 0, 0},
	{"hamming",  no_argument, 0, 'H'},
	{"quality",  no_argument, 0, 'Q'},
	{"fast",  no_argument, 0, 0},
	{"DPMismatch",  required_argument, 0, 'X'},
	{"DPMatch",  required_argument, 0, 'Y'},
	{"DPGapOpen",  required_argument, 0, 'G'},
	{"DPGapExt",  required_argument, 0, 'E'},
	{"extendIndelDetection", no_argument, 0, 0},
	{"allJunctions",  no_argument, 0, 0},
	{"memoryMultiplex",  required_argument, 0, 0},
	{"ignoreUnmapped",  no_argument, 0, 0},
	{"extraColumns",  no_argument, 0, 0},
	{"disableBigMargin",  no_argument, 0, 0},
	{"relaxMismatchedBases",  no_argument, 0, 0},
	{"reportPairedMultiBest",  no_argument, 0, 0},
	{"maxMismatches",  required_argument, 0, 'M'},
	{"exonicSubreadFrac",  required_argument, 0, 0},
	{"SVdetection", no_argument, 0, 0},
	{"maxVoteSimples",  required_argument, 0, 0},
	{"maxRealignLocations",  required_argument, 0, 0},
	{"minVoteCutoff",  required_argument, 0, 0},
	{"minMappedFraction",  required_argument, 0, 0},
	{"disableBigMargin",  no_argument, 0, 0},
	{"complexIndels", no_argument, 0, 0},
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
	SUBREADputs("    -r --read      <input>  name of the input file. Several input formats can be");
	SUBREADputs("                            automatically detected, including gzipped fastq,");
	SUBREADputs("                            fastq, and fasta. If paired-end, this should give");
	SUBREADputs("                            the name of file including first reads.");
	SUBREADputs("");
	SUBREADputs("Optional general arguments:");
	SUBREADputs("");
	SUBREADputs("    -o --output    <output> name of the output file. By default, the output file");
	SUBREADputs("                            includes BAM-format mapping result.");
	SUBREADputs("");
	SUBREADputs("    -n --subreads  <int>    number of selected subreads, 14 by default.");
	SUBREADputs("");
	SUBREADputs("    -m --minmatch  <int>    consensus threshold (minimal number of consensus");
	SUBREADputs("                            subreads required) for reporting a hit. If paired-");
	SUBREADputs("                            end read data are provided, this gives the consensus");
	SUBREADputs("                            threshold for the read which receives more votes");
	SUBREADputs("                            than the other read from the same pair. 1 by default");
	SUBREADputs("");
	SUBREADputs("    -T --threads   <int>    number of CPU threads, 1 by default.");
	SUBREADputs("");
	SUBREADputs("    -I --indel     <int>    maximum length (in bp) of indels that can be");
	SUBREADputs("                            detected. 5 by default. The program can detect");
	SUBREADputs("                            indels of up to 200bp long.");
	SUBREADputs("");
	SUBREADputs("    -B --multi     <int>    maximal number of equally-best mapping locations to");
	SUBREADputs("                            be reported. 1 by default. Note that -u option takes");
	SUBREADputs("                            precedence over -B.  ");
	SUBREADputs("");
	SUBREADputs("    -P --phred     <3:6>    the format of Phred scores used in input files, '3'");
	SUBREADputs("                            for phred+33 and '6' for phred+64. '3' by default.");
	SUBREADputs("");
	SUBREADputs("    -u --unique             report uniquely mapped reads only. Number of matched");
	SUBREADputs("                            bases (for RNA-seq) or mis-matched bases(for genomic");
	SUBREADputs("                            DNA-seq) is used to break the tie.");
	SUBREADputs("");
	SUBREADputs("    -b --color-convert      convert color-space read bases to base-space read");
	SUBREADputs("                            bases in the mapping output. Note that the mapping");
	SUBREADputs("                            itself will still be performed at color-space.");
	SUBREADputs("   ");
	SUBREADputs("    -M --maxMismatches <int> specify the maximum number of mis-matched bases");
	SUBREADputs("                            allowed in the alignment. 3 by default. Mis-matches");
	SUBREADputs("                            found in soft-clipped bases are not counted.");
	SUBREADputs("   ");
	SUBREADputs("       --SAMinput           specify that the input read data is in SAM format.");
	SUBREADputs("");
	SUBREADputs("       --BAMinput           specify that the input read data is in BAM format.");
	SUBREADputs("");
	SUBREADputs("       --SAMoutput          specify that mapping results are saved in the SAM");
	SUBREADputs("                            format.");
	SUBREADputs("");
	SUBREADputs("       --trim5     <int>    trim off <int> number of bases from 5' end of each");
	SUBREADputs("                            read. 0 by default.");
	SUBREADputs("");
	SUBREADputs("       --trim3     <int>    trim off <int> number of bases from 3' end of each");
	SUBREADputs("                            read. 0 by default.");
	SUBREADputs("");
	SUBREADputs("       --rg-id     <string> specify the read group ID. If specified,the read");
	SUBREADputs("                            group ID will be added to the read group header");
	SUBREADputs("                            field and also to each read in the mapping output.");
	SUBREADputs("");
	SUBREADputs("       --rg        <string> add a <tag:value> to the read group (RG) header in");
	SUBREADputs("                            in the mapping output.");
	SUBREADputs("");
	SUBREADputs("       --DPGapOpen  <int>   penalty for gap opening in short indel detection. -1");
	SUBREADputs("                            by default.");
	SUBREADputs("");
	SUBREADputs("       --DPGapExt   <int>   penalty for gap extension in short indel detection.");
	SUBREADputs("                            0 by default.");
	SUBREADputs("");
	SUBREADputs("       --DPMismatch <int>   penalty for mis-matched bases in short indel");
	SUBREADputs("                            detection. 0 by default.");
	SUBREADputs("");
	SUBREADputs("       --DPMatch    <int>   score for matched bases in short indel detection. 2");
	SUBREADputs("                            by default.");
	SUBREADputs("   ");
	SUBREADputs("       --allJunctions       detect exon-exon junctions (both canonical and non-");
	SUBREADputs("                            canonical junctions) and fusion breakpoints in RNA-");
	SUBREADputs("                            seq data. Refer to Users Guide for reporting of");
	SUBREADputs("                            junctions and fusions.");
	SUBREADputs("   ");
	SUBREADputs("       --complexIndels      detect multiple short indels that occur concurrently");
	SUBREADputs("                            in a small genomic region (these indels could be as");
	SUBREADputs("                            close as 1bp apart).");
	SUBREADputs("");
	SUBREADputs("    -v                      output version of the program.");
	SUBREADputs("");
	SUBREADputs("");
	SUBREADputs("Optional arguments for paired-end reads:");
	SUBREADputs("");
	SUBREADputs("    -R --read2     <input>  name of the file including second reads.");
	SUBREADputs("");
	SUBREADputs("    -p --minmatch2 <int>    consensus threshold for the non-anchor read ");
	SUBREADputs("                            (receiving less votes than the anchor read from the");
	SUBREADputs("                            same pair). 1 by default.");
	SUBREADputs("");
	SUBREADputs("    -d --mindist   <int>    minimum fragment/template length, 50bp by default.");
	SUBREADputs("");
	SUBREADputs("    -D --maxdist   <int>    maximum fragment/template length, 600bp by default.");
	SUBREADputs("");
	SUBREADputs("    -S --order     <ff:fr:rf> orientation of first and second reads, 'fr' by");
	SUBREADputs("                            default (forward/reverse).");
	SUBREADputs("");
	SUBREADputs("");
	SUBREADputs("Advanced arguments:");
	SUBREADputs("");
	SUBREADputs("For more information about these arguments, please refer to the User Manual.");
	SUBREADputs("");


}

int parse_opts_subjunc(int argc , char ** argv, global_context_t * global_context)
{
	int c;
	int option_index = 0;	
	int is_64_bit_computer = sizeof(char *)>4;

	optind = 0;
	opterr = 1;
	optopt = 63;

	global_context->config.entry_program_name = CORE_PROGRAM_SUBJUNC;
	global_context->config.max_mismatch_exonic_reads = 3;
	global_context->config.max_mismatch_junction_reads = 3;
	global_context->config.ambiguous_mapping_tolerance = 39 - 20 ;
	global_context->config.extending_search_indels = 0;
	global_context->config.do_fusion_detection =0;
	global_context->config.use_dynamic_programming_indel = 1;

	global_context->config.do_breakpoint_detection = 1;
	global_context->config.total_subreads = 14;
	global_context->config.minimum_subread_for_first_read =1;
	global_context->config.minimum_subread_for_second_read = 1;
	global_context->config.minimum_exonic_subread_fraction = 0.3;
	global_context->config.high_quality_base_threshold = 990000;
	global_context->config.do_big_margin_filtering_for_junctions = 1;
	global_context->config.report_no_unpaired_reads = 0;
	global_context->config.experiment_type = CORE_EXPERIMENT_RNASEQ;

	//#warning " ========================= REMOVE ' + 1 ' FROM THE NEXT LINE !! =========================="
	global_context->config.limited_tree_scan = 0;
	global_context->config.use_hamming_distance_in_exon = 0;
	//global_context->config.big_margin_record_size = 24;

	if(argc<2)
	{
		print_usage_core_subjunc();
		return -1;
	}
	while ((c = getopt_long (argc, argv, "vxsJ1:2:S:L:AHd:D:n:m:p:P:R:r:i:l:o:G:Y:E:X:T:I:B:bQFcuUfM:3:5:9:?", long_options, &option_index)) != -1)
	{
		switch(c)
		{
			case 'v':
				core_version_number("Subjunc");
				return -1;
			case 'G':
				global_context->config.DP_penalty_create_gap = atoi(optarg);
				break;
			case 'Y':
				global_context->config.DP_match_score = atoi(optarg);
				break;
			case 'E':
				global_context->config.DP_penalty_extend_gap = atoi(optarg);
				break;
			case 'X':
				global_context->config.DP_mismatch_penalty = atoi(optarg);
				break;
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
				global_context->config.show_soft_cliping = 0;
				break;
			case 'Q':
				global_context->config.use_quality_score_break_ties = 1;
				break;
			case 'H':
				global_context->config.use_hamming_distance_break_ties = 1;
				break;
			case 's':
				global_context->config.downscale_mapping_quality = 1;
				break;
			case 'A':
				global_context->config.report_sam_file = 0;
				break;
			case 'S':
				global_context->config.is_first_read_reversed = optarg[0]=='r'?1:0;
				global_context->config.is_second_read_reversed = optarg[1]=='f'?0:1;
				break;
			case 'U':
				global_context->config.report_no_unpaired_reads = 1;
				break;
			case 'u':
				global_context->config.report_multi_mapping_reads = 0;
				global_context->config.use_hamming_distance_break_ties = 1;
				break;
			case 'b':
				global_context->config.convert_color_to_base = 1;
				break;
			case 'D':
				global_context->config.maximum_pair_distance = atoi(optarg);
				break;
			case 'd':
				global_context->config.minimum_pair_distance = atoi(optarg);
				if(global_context->config.minimum_pair_distance <0)
					global_context->config.restrected_read_order = 0;
				break;
			case 'n':
				global_context->config.total_subreads = atoi(optarg);
				//global_context->config.total_subreads = min(31,global_context->config.total_subreads);
				break;
			case 'm':
				global_context->config.minimum_subread_for_first_read = atoi(optarg);
				break;
			case 'T':
				global_context->config.all_threads = atoi(optarg);
				if(global_context->config.all_threads <1) global_context->config.all_threads = 1;
				if(global_context->config.all_threads > MAX_THREADS) global_context->config.all_threads = MAX_THREADS;

				break;
			case 'M':
				global_context->config.max_mismatch_exonic_reads = atoi(optarg);
				global_context->config.max_mismatch_junction_reads = atoi(optarg);
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

				if(!is_64_bit_computer) global_context->config.max_indel_length = min(global_context->config.max_indel_length , 16); 
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
					//global_context->config.max_mismatch_exonic_reads = 2;
					//global_context->config.max_mismatch_junction_reads = 2;
					//global_context->config.total_subreads = 28;
					//global_context->config.minimum_subread_for_first_read = 3;
					//global_context->config.minimum_subread_for_second_read = 1;
					//global_context->config.do_big_margin_filtering_for_reads = 1;

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
			case 'F':
				global_context->config.is_second_iteration_running = 0;
				global_context->config.report_sam_file = 0;
				break;
			case 'B':
				global_context->config.multi_best_reads = atoi(optarg); 

				if(global_context->config.multi_best_reads<1)
					global_context->config.multi_best_reads=1;

				global_context->config.reported_multi_best_reads = global_context->config.multi_best_reads;

				global_context->config.max_vote_combinations = max(global_context->config.max_vote_combinations, global_context->config.reported_multi_best_reads + 1);
				global_context->config.max_vote_simples = max(global_context->config.max_vote_simples, global_context->config.reported_multi_best_reads + 1);

				break;
			case 'c':
				global_context->config.space_type = GENE_SPACE_COLOR; 
				break;
				
			case 0:
				if(strcmp("memoryMultiplex", long_options[option_index].name)==0) 
				{
					global_context->config.memory_use_multiplex = atof(optarg);
				}
				else if(strcmp("ignoreUnmapped", long_options[option_index].name)==0) 
				{
					global_context->config.ignore_unmapped_reads = 1;
				}
				else if(strcmp("rg-id", long_options[option_index].name)==0) 
				{
					strcpy(global_context->config.read_group_id, optarg);
				}
				else if(strcmp("rg", long_options[option_index].name)==0) 
				{
					strcat(global_context->config.read_group_txt, "\t");
					strcat(global_context->config.read_group_txt, optarg);
				}
				else if(strcmp("SAMoutput", long_options[option_index].name)==0) 
				{
					global_context->config.is_BAM_output = 0;
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
				else if(strcmp("extraColumns", long_options[option_index].name)==0) 
				{
					global_context->config.SAM_extra_columns=1;
				}
				else if(strcmp("minMappedFraction", long_options[option_index].name)==0) 
				{
					global_context->config.min_mapped_fraction = atoi(optarg);
				}
				else if(strcmp("relaxMismatchedBases", long_options[option_index].name)==0) 
				{
					global_context->config.max_mismatch_junction_reads = 999;
					global_context->config.max_mismatch_exonic_reads = 999;
					global_context->config.min_mapped_fraction = 61;
				}
				else if(strcmp("reportPairedMultiBest", long_options[option_index].name)==0) 
				{
					global_context->config.report_multiple_best_in_pairs = 1;
				}
				else if(strcmp("exonicSubreadFrac", long_options[option_index].name)==0) 
				{
					if(atof(optarg)>0)
						global_context->config.minimum_exonic_subread_fraction = atof(optarg);
					else	SUBREADprintf("WARNING: unknown parameter: --exonicSubreadFrac '%s'\n", optarg);
				}
				else if(strcmp("fast", long_options[option_index].name)==0) 
				{
					global_context -> config.fast_run = 1;
				}
				else if(strcmp("SVdetection", long_options[option_index].name)==0) 
				{
					global_context -> config.do_structural_variance_detection = 1;
				}
				else if(strcmp("minDistanceBetweenVariants", long_options[option_index].name)==0)
				{
					int newdist = atoi(optarg);
					newdist = max(newdist, 1);
					newdist = min(newdist, MAX_READ_LENGTH);
					global_context->config.realignment_minimum_variant_distance = newdist;
				}
				else if(strcmp("disableBigMargin", long_options[option_index].name)==0) 
				{
					global_context->config.do_big_margin_filtering_for_junctions = 0;
					global_context->config.limited_tree_scan = 0;
				}
				else if(strcmp("maxVoteSimples", long_options[option_index].name)==0)
				{
					global_context->config.max_vote_simples = atoi(optarg);
				}
				else if(strcmp("maxRealignLocations", long_options[option_index].name)==0)
				{
					global_context->config.max_vote_combinations = atoi(optarg);
					global_context->config.multi_best_reads = atoi(optarg);
				}
				else if(strcmp("complexIndels", long_options[option_index].name)==0)
				{
					global_context->config.maximise_sensitivity_indel = 1;
					global_context->config.realignment_minimum_variant_distance = 1;
					global_context->config.max_indel_length = 16;
				}
				else if(strcmp("disableBigMargin", long_options[option_index].name)==0)
				{
					global_context->config.big_margin_record_size = 0;
				}
				else if(strcmp("extendIndelDetection", long_options[option_index].name)==0)
				{
					global_context->config.extending_search_indels = 1;
				}
				else if(strcmp("minVoteCutoff", long_options[option_index].name)==0)
				{
					global_context->config.max_vote_number_cutoff  = atoi(optarg);
				}
				else if(strcmp("allJunctions", long_options[option_index].name)==0)
				{
					global_context->config.do_fusion_detection = 1;
				}

				break;


			case '?':
			default:
				print_usage_core_subjunc();
				return -1 ;
		}
	}

	if(global_context->config.is_SAM_file_input) global_context->config.phred_score_format = FASTQ_PHRED33;
	global_context->config.more_accurate_fusions = global_context->config.more_accurate_fusions && global_context->config.do_fusion_detection;


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

