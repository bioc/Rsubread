#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "exon-algorithms.h"
#include "gene-algorithms.h"
#include "gene-value-index.h"
#include "input-files.h"
#include "sorted-hashtable.h"

int ALL_THREADS_junction=1;
int TOTAL_SUBREADS;
int ACCEPT_SUBREADS;
int ACCEPT_MINOR_SUBREADS;
int INDEX_THRESHOLD;
int MAX_PAIRED_DISTANCE_junction = 300;
int MIN_PAIRED_DISTANCE_junction = 200;
int INDEL_TOLERANCE_junction = 0;
int IS_DEBUG_junction = 0;
int QUALITY_SCALE_junction = 0;
int USE_VALUE_ARRAY_INDEX_junction = 1;
int FIRST_READ_REVERSE_junction = 0;
int SECOND_READ_REVERSE_junction = 1;
double reads_density;

void print_res_junction(gene_exon_allrecords_t *av, gene_input_t* ginp,gene_input_t * ginp2, FILE * out_fp, char * index_prefix, unsigned int processed_reads, unsigned long long int all_processed_reads,  unsigned long long int *succeed_reads, gene_value_index_t * my_value_array_index)
{
	int i=0;
	gene_offset_t offsets;

	load_offsets (&offsets, index_prefix);

	while (1)
	{
		char nameb[1201], inb[1201], qualityb[1201];
		if(i > processed_reads)break;
		int rl = geinput_next_read(ginp, nameb, inb, qualityb);
		if (rl<0){
			break;
		}

		unsigned int ex1p=0, ex2p=0;
		short ex1len=0, ex2len=0;
		char ex1_reversed=0, ex2_reversed=0, ex1gd=0, ex2gd=0;

		if(find_best_edges(rl, inb, &(av->arenas[i]), &ex1p,&ex2p,&ex1len,&ex2len,&ex1_reversed,&ex2_reversed, &ex1gd, &ex2gd))
		{
			char * e1name, * e2name;
			unsigned int e1chrpos, e2chrpos;
			if(locate_gene_position(ex1p, &offsets, &e1name, &e1chrpos) || locate_gene_position(ex2p, &offsets, &e2name, &e2chrpos))
				printf ("ERROR: position out of range.\n" );
			else{
				int value = 1000000;
				if (ex1p>ex2p)
					value -= ex1p-ex2p;
				else
					value -= ex2p-ex1p;

				fprintf(out_fp, "> %s\n", nameb);
				//fprintf(out_fp, "@%c %d + %d\n", (ex1len+ex2len > rl)?'!':'?', ex1len,ex2len);
				fprintf(out_fp, "@%c %c %d + %c %d\n", (value > 10) && (value < 1000000)?'!':'?',ex1gd?'G':'b', ex1len, ex2gd?'G':'b',ex2len);
				fprintf(out_fp, "Ex1: %s,%u %c\tEx2: %s,%u %c\t\n" ,e1name, e1chrpos, ex1_reversed?'~':'@', e2name, e2chrpos, ex2_reversed?'~':'@');
				fputs(inb, out_fp);
				fputs("\n", out_fp);
				reverse_read(inb, rl, ginp->space_type);	
				fputs(inb, out_fp);
				fputs("\n\n" , out_fp);
			}
		}
		else
		{
			fprintf(out_fp, "> %s\n", nameb);
			fprintf(out_fp, "@ 0 + 0\n");
		}

		i++;
	}


}

struct gene_thread_data
{
	gehash_t table;
	int * offsets;
	int number_t;
};


struct gene_thread_data_transport
{
	gehash_t * my_table;
	gene_value_index_t * my_value_array_index;
	int table_no;
	int all_tables;
	char * index_prefix;
	gene_exon_allrecords_t * all_vote;
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

int run_search_junction(gehash_t * my_table, gene_value_index_t * my_value_array_index , int table_no,  gene_exon_allrecords_t * allvote, gene_input_t * ginp,gene_input_t * ginp2, char * index_prefix, unsigned int * processed_reads, long long int base_number, int all_tables, pthread_spinlock_t * input_lock, int my_thread_no, unsigned int section_length)
{
//	FILE * fp;
	char BuffMemory [2500];
	char * InBuff = NULL , * InBuff2 = NULL;
	char BuffMemory2 [2500];
	char * QualityBuff = NULL, * QualityBuff2 = NULL;

	int all_reads = allvote -> max_len;
	int queries = 0;
	double t0=miltime();
	int is_reversed;
	int good_match = 0;
	float subread_step = 3;
	struct stat read_fstat;
	stat (ginp->filename, &read_fstat);
	long long int read_fsize = read_fstat.st_size;
	
	double local_begin_ftime = miltime();
	int read_len = 0, read2_len = 0;
	is_reversed = 0;
	if (ginp2)
		all_reads /=2;

//	printf ("I'm the %d-th thread\n", my_thread_no);


	while (1)
	{

		if (is_reversed)
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

					if(FIRST_READ_REVERSE_junction)
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
	
					if(SECOND_READ_REVERSE_junction)
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

//		printf("\n === NEW READ q=%d rev=%d === \n", queries, is_reversed);

		int i;
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

						//printf("Q:%d I:%d 2nd:%d SR:%d INV:%d :: Pos:%d Len:%d\n", queries, i, is_paired, subread_no, is_reversed, subread_offset, Curr_read_len);

						char * subread_string = CurrInBuff + subread_offset;
						gehash_key_t subread_integer = genekey2int(subread_string, ginp->space_type);

						gene_vote_number_t subread_quality = 22.;

						if(CurrQualityBuff[0])
						{
							char * quality_string = CurrQualityBuff + subread_offset;
							subread_quality += get_subread_quality(quality_string, subread_string, QUALITY_SCALE_junction);
						}
						else subread_quality -= .1;
						if(is_valid_subread(subread_string))
							gehash_go_q(my_table, subread_integer , subread_offset , vote , 1 , subread_quality, INDEX_THRESHOLD, INDEL_TOLERANCE_junction);

						if(subread_offset1 >= Curr_read_len -16)
							break;
					}
				}
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
						subread_quality += get_subread_quality(quality_string, subread_string, QUALITY_SCALE_junction);
					}
					else subread_quality -= .1;

					if(is_valid_subread(subread_string))
						total_results = gehash_go_q(my_table, subread_integer , subread_offset , &vote , 1 , subread_quality, INDEX_THRESHOLD, INDEL_TOLERANCE_junction);

				}
				if(subread_offset1 >= read_len -16)
					break;

			}

			if(vote.max_vote > (subread_no - (TOTAL_SUBREADS - ACCEPT_SUBREADS))*22.)
			{
				test_read_conjunction(InBuff, read_len, &vote, my_value_array_index, INDEL_TOLERANCE_junction, allvote, queries, is_reversed, my_value_array_index);
			}
	
		}


		if (my_thread_no==0 && queries % 10000 == 0 && !is_reversed)
		{
			if(table_no == 0)
			{
				long long int current_reads = base_number + queries;
				long long int fpos = ftello(ginp->input_fp);
				reads_density = fpos*1.0/current_reads; 
			}
			if(IS_DEBUG_junction && queries % 100000==0)
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
	
		if (is_reversed)
		{
			if(queries >= all_reads)
				break;
		}

		is_reversed = !is_reversed;
	}

	return 0;
}


void * run_search_thread_junction(void * parameters)
{
	struct gene_thread_data_transport * data_param = parameters;
	int thid = data_param->this_thread;

	pthread_spin_unlock(data_param -> init_lock);

	run_search_junction(data_param->my_table, data_param->my_value_array_index, data_param->table_no, data_param->all_vote, data_param->ginp , data_param->ginp2, data_param->index_prefix, data_param->processed_reads, data_param->base_number, data_param->all_tables, data_param->input_data_lock, thid, data_param-> section_length);
	return NULL;
}



// This function search a segment of reads (length = read_number) 
// It returns the number of reads that were really processed;
int run_search_index_junction(gene_input_t * ginp, gene_input_t * ginp2, char * index_prefix, gene_exon_allrecords_t * all_vote, FILE * out_fp, unsigned long long int base_number, int all_tables, unsigned long long int *succeed_reads)
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

		if (IS_DEBUG_junction)
			printf ("@LOG Loading table from %s\n", table_fn);
		else
			printf ("Loading the %02d-th index file ...					      \r", tabno+1);
		fflush(stdout);

		gehash_load(my_table, table_fn);
		if(USE_VALUE_ARRAY_INDEX_junction)
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
		if(ALL_THREADS_junction <2)
			run_search_junction(my_table, &value_array_index, tabno, all_vote, ginp , ginp2, index_prefix, &processed_reads, base_number, all_tables, NULL /*the data lock is null*/, 0  /*I'm the 0-th thread*/, section_length);
		else
		{
			int i; 
			struct gene_thread_data_transport data_param;
			pthread_t runners [ALL_THREADS_junction];
			pthread_spinlock_t  data_lock;
			pthread_spinlock_t  init_lock;

			data_param.my_table = my_table;
			data_param.my_value_array_index = &value_array_index;
			data_param.table_no = tabno;
			data_param.all_tables = all_tables;
			data_param.index_prefix = index_prefix;
			data_param.all_vote = all_vote;
			data_param.base_number = base_number;
			data_param.all_threads = ALL_THREADS_junction;
			data_param.input_data_lock = &data_lock;
			data_param.init_lock = &init_lock;
			data_param.processed_reads = &processed_reads;
			data_param.section_length = section_length;
			data_param.ginp = ginp;
			data_param.ginp2 = ginp2;

			pthread_spin_init(&data_lock, PTHREAD_PROCESS_PRIVATE);
			pthread_spin_init(&init_lock, PTHREAD_PROCESS_PRIVATE);
			pthread_spin_lock(&init_lock);
			for (i=0; i< ALL_THREADS_junction; i++)
			{
				data_param.this_thread = i;
				pthread_create(runners+i, NULL, run_search_thread_junction, &data_param);
				pthread_spin_lock(&init_lock);
			}

			for (i=0; i< ALL_THREADS_junction; i++)
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

	print_res_junction(all_vote, ginp, ginp2, out_fp, index_prefix, processed_reads, processed_reads,  succeed_reads, &value_array_index);
	return processed_reads;
	
}

void usage_junction(char * execname)
{
	puts("Usage:");
	puts(" ./subread-exons [options] -i <index_name> -r <input> -o <output>");
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
	{"basewise", no_argument, &USE_VALUE_ARRAY_INDEX_junction, 1},
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



int main_junction(int argc,char ** argv)
{
	char read_file [300], read2_file [300];
	char output_file [300];

	char index_prefix [300];
	unsigned int all_reads, all_tables;
	gene_exon_allrecords_t allvote;
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
				FIRST_READ_REVERSE_junction = optarg[0]=='r'?1:0;
				SECOND_READ_REVERSE_junction = optarg[1]=='f'?0:1;
				break;
			case 'a':
				USE_VALUE_ARRAY_INDEX_junction = 1;
				break;
			case 'D':
				MAX_PAIRED_DISTANCE_junction = atoi(optarg);
				break;
			case 'd':
				MIN_PAIRED_DISTANCE_junction = atoi(optarg);
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
				ALL_THREADS_junction = atoi(optarg);
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
				INDEL_TOLERANCE_junction = atoi(optarg);
				if( INDEL_TOLERANCE_junction >5)INDEL_TOLERANCE_junction=5;
				break ;
			case 'Q':
				if(optarg[0]=='l')
					QUALITY_SCALE_junction = QUALITY_SCALE_LINEAR;
				if(optarg[0]=='e')
					QUALITY_SCALE_junction = QUALITY_SCALE_LOG;
				if(optarg[0]=='n')
					QUALITY_SCALE_junction = QUALITY_SCALE_NONE;
				break;
			case 'p':
				ACCEPT_MINOR_SUBREADS = atoi(optarg);
				break;
			case '?':
				return -1 ;
		}

	if (!read_file[0] || !index_prefix[0] || !output_file[0])
	{
		usage_junction(argv[0]);

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
	printf("Number of threads=%d\n", ALL_THREADS_junction);
	printf("Tolerance for Indel=%d\n", INDEL_TOLERANCE_junction);
	if (QUALITY_SCALE_junction==QUALITY_SCALE_LINEAR)
		puts("Quality scale=linear\n\n");
	else if (QUALITY_SCALE_junction==QUALITY_SCALE_LOG)
		puts("Quality scale=exponential\n\n");
	else 	puts("\n");

	if (read2_file[0])
	{
		if (MAX_PAIRED_DISTANCE_junction <= MIN_PAIRED_DISTANCE_junction)
		{
			printf ("The value of the '-D' option must be greater than that of the '-d' option. \n");
			return -1;
		}

		printf ("Performing paired-end alignment:\n");
		printf ("Maximum distance between reads=%d\n", MAX_PAIRED_DISTANCE_junction);
		printf ("Minimum distance between reads=%d\n", MIN_PAIRED_DISTANCE_junction);
		printf ("Threshold on number of subreads for a successful mapping (the minor end in the pair)=%d\n", ACCEPT_MINOR_SUBREADS);
		printf ("The directions of the two input files are: %s, %s\n\n", FIRST_READ_REVERSE_junction?"reversed":"forward", SECOND_READ_REVERSE_junction?"reversed":"forward");
	}

#ifdef REPORT_ALL_THE_BEST
	printf("***** WARNING: the REPORT_ALL_THE_BEST switch is turned on. You need an extra 1 GBytes of RAM space for saving the temporary results. *****\n");
#endif

	init_exon_arena(&allvote, all_reads);
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

		processed_reads += run_search_index_junction(&ginp, read2_file[0] ? (&ginp2):NULL, index_prefix, &allvote, out_fp, processed_reads, all_tables, &succeed_reads);
		clear_exon_arena(&allvote);

		// test if there no anyreads remaining
		unsigned long long int current_fp = ftello(ginp.input_fp);
		int rl = geinput_next_read(&ginp, NULL, inbuff, NULL);
		if (rl<0)
			break;
		fseeko(ginp.input_fp, current_fp, SEEK_SET);
	}

	geinput_close(&ginp);
	if(IS_DEBUG_junction)
		printf("@LOG THE END. \n");
	else
		printf("\n\n %llu reads were processed in %.1f seconds.\nPercentage of successfully mapped reads is %0.2f%%.\n\n", processed_reads, miltime()-begin_ftime, succeed_reads*100.0/processed_reads/(read2_file[0]?2:1));

	printf("\n\nCompleted successfully.\n");

	return 0;
}
