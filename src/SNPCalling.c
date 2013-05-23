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
  
  
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "subread.h"
#include "gene-algorithms.h"
#include "SNPCalling.h"
#include "input-files.h"

struct SNP_Calling_Parameters{

	float supporting_read_rate;
	int max_supporting_read_number;
	int min_supporting_read_number;

	int neighbour_filter_testlen;
	float neighbour_filter_rate;
	int bases_ignored_head_tail;

	float fisher_exact_p_threshold;
	int fisher_exact_testlen;

	int min_phred_score;
	int test_two_strands;

};

#define PRECALCULATE_FACTORIAL 50000

double * precalculated_factorial;// [PRECALCULATE_FACTORIAL];

double factorial_float_real(int a)
{

	double ret = 0;
	while(a)
		ret += log(a--);
	return ret;
}


double factorial_float(int a)
{

	if(a<PRECALCULATE_FACTORIAL && (precalculated_factorial[a]!=0.))
		return precalculated_factorial[a]; 
	else
	{
		double ret = factorial_float_real(a);
		if(a<PRECALCULATE_FACTORIAL) precalculated_factorial[a]=ret;
		return ret;
	}
}

double fisherSub(int a, int b, int c, int d)
{
	double ret = factorial_float(a+b) + factorial_float(c+d) + factorial_float(a+c) + factorial_float(b+d) ;
	ret -= factorial_float(a) + factorial_float(b) + factorial_float(c) + factorial_float(d) + factorial_float(a+b+c+d);
	return pow(2.71828183, ret);
}




/**
 * See HELP string or run with no arguments for usage.
 * <p>
 * The code used to calculate a Fisher p-value comes originally from a
 * <a href="http://infofarm.affrc.go.jp/~kadasowa/fishertest.htm">JavaScript program</a>
 * by T. Kadosawa (kadosawa@niaes.affrc.go.jp).
 * Retrieved from http://www.users.zetnet.co.uk/hopwood/tools/StatTests.java on 3/Jul/2012
 *
 * @author David Hopwood
 * @date   2000/04/23
 */

float fisher_exact_test(int a, int b, int c, int d)
{

	if(a*1./c < b*1./d)	return 1.1;
		// the abnormal level at this base should be at least as the noise level.

	if (a * d > b * c) {
            a = a + b; b = a - b; a = a - b; 
            c = c + d; d = c - d; c = c - d;
        }

        if (a > d) { a = a + d; d = a - d; a = a - d; }
        if (b > c) { b = b + c; c = b - c; b = b - c; }

        double p_sum = 0.0;
        double p = fisherSub(a, b, c, d);


        while (a >= 0) {
            p_sum += p;
            if (a == 0) break;
            --a; ++b; ++c; --d;
            p = fisherSub(a, b, c, d);
        }

	return p_sum;
}

unsigned int fisher_test_size;

int process_snp_votes(FILE *out_fp, unsigned int offset , unsigned int reference_len, char * referenced_genome, char * chro_name , char * temp_prefix, struct SNP_Calling_Parameters * parameters)
{
	int block_no = (offset -1) / BASE_BLOCK_LENGTH, i;
	char temp_file_name[300];
	FILE *tmp_fp;
	unsigned int * snp_voting_table_Pos, *snp_voting_table_Neg;	// offset * 4 + "A/C/G/T"[0,1,2,3]
	float * snp_fisher_raws;

	sprintf(temp_file_name , "%s%s-%04u.bin", temp_prefix, chro_name, block_no);
	tmp_fp = fopen(temp_file_name, "rb");

	// if no temp file is here, do nothing 
	if(!tmp_fp)return 0;

	snp_voting_table_Pos = (unsigned int *)SUBREAD_malloc(sizeof(unsigned int) * reference_len*4); 
	snp_voting_table_Neg = (unsigned int *)SUBREAD_malloc(sizeof(unsigned int) * reference_len*4); 
	snp_fisher_raws = (float *)SUBREAD_malloc(sizeof(float)  * reference_len);
	if((!snp_voting_table_Neg) || (!snp_fisher_raws))
	{
		fatal_memory_size();
		return -1;
	}

	memset(snp_voting_table_Pos,0 ,sizeof(unsigned int) * reference_len*4);
	memset(snp_voting_table_Neg,0 ,sizeof(unsigned int) * reference_len*4);

	for(i=0; i<reference_len; i++)
		snp_fisher_raws[i]=-1.;

	while(!feof(tmp_fp))
	{
		base_block_temp_read_t read_rec;
		unsigned short read_len;
		unsigned int first_base_pos;
		char read[MAX_READ_LENGTH];
		char qual[MAX_READ_LENGTH];
		char base_neighbour_test[MAX_READ_LENGTH];
		char base_used[MAX_READ_LENGTH];

		fread(&read_rec, sizeof(read_rec), 1, tmp_fp);
		fread(&read_len, sizeof(short), 1, tmp_fp);
		fread(read, sizeof(char), read_len, tmp_fp);
		fread(qual, sizeof(char), read_len, tmp_fp);

		first_base_pos = read_rec.pos - block_no * BASE_BLOCK_LENGTH;


		if(first_base_pos + read_len -1> reference_len || first_base_pos<1)
		{
			SUBREADprintf("WARNING: read length %u+%d out of boundary: %u at the %d-th block.\n", first_base_pos, read_len, reference_len, block_no);
			continue;
		}

		if(parameters->neighbour_filter_testlen)
		{
			int j;
			memset(base_used,1,MAX_READ_LENGTH);
			for(i=0;i<read_len;i++)
			{
				char true_value = referenced_genome[i + first_base_pos -1];
				base_neighbour_test[i]=0;
				if(qual[i] >='!'+ parameters->min_phred_score)
					if(true_value == read[i])
						base_neighbour_test[i]=1;
			}

			for(i=0;i<read_len;i++)
			{
				int correct=0, wrong=0;
				for(j=max(0, i-parameters->neighbour_filter_testlen); j<min(read_len, i+parameters->neighbour_filter_testlen); j++)
				{
					if(base_neighbour_test[j])correct++;
					else wrong++;
				}
				if(correct*1./(correct+wrong)<parameters->neighbour_filter_rate)
					base_used[i]=0;
			}
		}
		for(i=0;i<read_len;i++)
		{
			char base_int = 3;

			if(qual[i] < '!'+parameters->min_phred_score)continue;
			if(parameters->neighbour_filter_testlen &&  !base_used[i])continue;

			if(i+first_base_pos > reference_len  || i+first_base_pos<1)
			{
				SUBREADprintf("Warning: read out of boundary: %u >= %u\n", i+first_base_pos, reference_len);
				break;
			}

			switch(read[i])
			{
				case 'A':
					base_int = 0;
					break;
				case 'C':
					base_int = 1;
					break;
				case 'G':
					base_int = 2;
			}
			int j;
			unsigned int supporting_bases = 0;
			for(j=0; j<4; j++)
			{
				supporting_bases += snp_voting_table_Neg[(first_base_pos+i-1)*4+j];
				supporting_bases += snp_voting_table_Pos[(first_base_pos+i-1)*4+j];
			}
			if(supporting_bases < parameters->max_supporting_read_number)
			{
				if(read_rec.strand)
					snp_voting_table_Neg[(first_base_pos+i-1)*4+base_int]++;
				else
					snp_voting_table_Pos[(first_base_pos+i-1)*4+base_int]++;
			}
		}
	}
	fclose(tmp_fp);

	unlink(temp_file_name);
	// Finding SNPs from the finished voting table

	char * base_is_reliable = NULL;


	if(parameters -> fisher_exact_testlen || parameters -> test_two_strands)
	{
		base_is_reliable = (char *)malloc(sizeof(char) * reference_len);
		memset(base_is_reliable,1 , sizeof(char) * reference_len);
	}

	if(parameters -> test_two_strands)
	{
		for(i= 0;i<reference_len; i++)
		{
			int j;
			int pos_match=0, pos_unmatch=0, neg_match=0, neg_unmatch=0;

			char true_value = referenced_genome[i + parameters -> fisher_exact_testlen];
			int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));

			for(j=0;j<4;j++)
			{
				if(j==true_value_int)
				{
					pos_match += snp_voting_table_Pos[i*4+j];
					neg_match += snp_voting_table_Neg[i*4+j];
				}
				else
				{
					pos_unmatch += snp_voting_table_Pos[i*4+j];
					neg_unmatch += snp_voting_table_Neg[i*4+j];
				}

			}

			float p_diff = fisher_exact_test(pos_match, neg_match, pos_unmatch, neg_unmatch);
			if(p_diff<0.1)
				base_is_reliable[i]=0;
		}
	}

	if(parameters -> fisher_exact_testlen)
	{
		int j;
		int a=0, b=0, c=0, d=0;
		long long int reference_len_long = reference_len;
		/**    | This | All_Window 
		 * ----+------+-------
		 * #mm |  a   |  b
		 * #pm |  c   |  d
		 **/


		for(i= - parameters -> fisher_exact_testlen;i<reference_len_long; i++)
		{
			a=0; c=0;
			for(j=0;j<4;j++)
			{
				if(i+parameters -> fisher_exact_testlen < reference_len_long)
				{
					char true_value = referenced_genome[i + parameters -> fisher_exact_testlen];
					int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));
					if(j == true_value_int)
						d += snp_voting_table_Pos[ (i+parameters -> fisher_exact_testlen) *4 + j ] + snp_voting_table_Neg[ (i+parameters -> fisher_exact_testlen) *4 + j ];
					else
						b += snp_voting_table_Pos[ (i+parameters -> fisher_exact_testlen) *4 + j ] + snp_voting_table_Neg[ (i+parameters -> fisher_exact_testlen) *4 + j ];
				}
				
				if(i>=0)
				{
					char true_value = referenced_genome[i];
					int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));

					if(j == true_value_int)
						c  =  snp_voting_table_Pos[ i * 4 + j ] + snp_voting_table_Neg[ i * 4 + j ];
					else
						a +=  snp_voting_table_Pos[ i * 4 + j ] + snp_voting_table_Neg[ i * 4 + j ];
				}

			}
		
			// test the middle base
			if(i>=0 && a > 0){
				float p_middle = fisher_exact_test(a, b-a, c, d-c);
				if(p_middle > parameters->fisher_exact_p_threshold)
					snp_fisher_raws [i] = - 2.0;
				else
					snp_fisher_raws [i] = p_middle;
				fisher_test_size ++;
			
			}

			for(j=0;j<4;j++)
				if(i >= parameters -> fisher_exact_testlen)
				{
					char true_value = referenced_genome[i - parameters -> fisher_exact_testlen];
					int  true_value_int =  (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3)));
					if(j == true_value_int)
						d -= snp_voting_table_Pos[ (i-parameters -> fisher_exact_testlen) *4 + j ] + snp_voting_table_Neg[ (i-parameters -> fisher_exact_testlen) *4 + j ];
					else
						b -= snp_voting_table_Pos[ (i-parameters -> fisher_exact_testlen) *4 + j ] + snp_voting_table_Neg[ (i-parameters -> fisher_exact_testlen) *4 + j ];
				}
	
		}
	}


	

	for(i=0;i<reference_len; i++)
	{
		char true_value = referenced_genome[i];
		int tested_int;
		int all_reads = 0;
		char base_list[10], supporting_list[55], snps=0;

		for(tested_int=0; tested_int <4; tested_int ++)
			all_reads += snp_voting_table_Pos[i*4+tested_int] + snp_voting_table_Neg[i*4+tested_int];

		if(all_reads<parameters->min_supporting_read_number) continue;
		if(base_is_reliable && !(base_is_reliable[i])) continue;
		base_list[0]=0;
		supporting_list[0]=0;

		for(tested_int=0; tested_int <4; tested_int ++)
		{
			if(tested_int != (true_value=='A'?0:(true_value=='C'?1:(true_value=='G'?2:3))))
				if((snp_voting_table_Pos[i*4+tested_int]+snp_voting_table_Neg[i*4+tested_int]) *1.0 / all_reads >= parameters->supporting_read_rate)
				{
					char int_buf[9];
					sprintf(int_buf, "%u", snp_voting_table_Pos[i*4+tested_int] + snp_voting_table_Neg[i*4+tested_int]);
					if(snps>0)
					{
						base_list[snps*2-1] = '/';
						base_list[snps*2] = tested_int==0?'A':(tested_int == 1?'C':(tested_int == 2?'G':'T'));
						base_list[snps*2+1] = 0;
						strcat(supporting_list,"/");
						strcat(supporting_list,int_buf);
					}
					else
					{
						base_list[0] = tested_int==0?'A':(tested_int == 1?'C':(tested_int == 2?'G':'T'));
						base_list[1] = 0;
						strcpy(supporting_list, int_buf);
					}
					snps++;
					
				}
		}
		if(snps)
			fprintf(out_fp, "%s\t%u\t%c\t%s\t%s\t%d\t%.9f\n", chro_name, BASE_BLOCK_LENGTH*block_no +1 + i, true_value,base_list, supporting_list , all_reads, snp_fisher_raws[i]);

	}

	if(base_is_reliable) free(base_is_reliable);
	free(snp_voting_table_Pos);
	free(snp_voting_table_Neg);
	free(snp_fisher_raws);
	return 0;
}


int run_chromosome_search(FILE *in_fp, FILE * out_fp, char * chro_name , char * temp_prefix, chromosome_t * chromosomes , struct SNP_Calling_Parameters* parameters)
{
	int chro_no;
	unsigned int offset=0, all_offset = 0;
	
	unsigned int chro_len=0;
	char * referenced_base;

	referenced_base = (char *) SUBREAD_malloc(sizeof(char)* BASE_BLOCK_LENGTH);
	if(!referenced_base)
	{
		fatal_memory_size();
		return -1;
	}


	for(chro_no=0;chro_no < OFFSET_TABLE_SIZE; chro_no++)
	{
		if(!chromosomes[chro_no].chromosome_name[0]) break;
		if(strcmp(chro_name , chromosomes[chro_no].chromosome_name)==0)
		{
			chro_len = chromosomes[chro_no].known_length;
			break;
		}
	}

	SUBREADprintf("Processing chromosome %s in FASTA file; expected length is %u.\n", chro_name, chro_len);
	fflush(stdout);
	if(!chro_len)
	{
		SUBREADprintf("Unknown chromosome name in FASTA file: %s\n", chro_name);
		free(referenced_base);
		return 1;
	}

	
	while(all_offset <= chro_len)
	{
		char nc = fgetc(in_fp);

		if(nc == ' ' || nc == '\r' || nc == '\n') continue;
		if(nc == '>')
			fseek(in_fp, -1, SEEK_CUR);

		if(nc != EOF && nc != '>')
		{
			nc = toupper(nc);
			referenced_base[offset++] = (nc=='A' || nc == 'G' || nc == 'C')?nc:'T';
			all_offset ++;
		}

		if((nc == '>'||nc == EOF) && all_offset < chro_len)
			SUBREADprintf("WARNING: Chromosome is shorter than expected: %s\n", chro_name);

		if(offset == BASE_BLOCK_LENGTH || nc == EOF || nc == '>')
		{
			process_snp_votes(out_fp, all_offset, offset, referenced_base, chro_name , temp_prefix, parameters);
			offset = 0;
		}

		if(nc == EOF || nc == '>') break;
	}

	free(referenced_base);
	

	return 0;
}



// This scan is driven by reading the FASTA file.
// While reading the FASTA file, if you see a new sequence, then find the base blocks for that sequence.
int parse_read_lists(char * in_FASTA_file, FILE * out_fp, char * temp_prefix, chromosome_t * chromosomes, struct SNP_Calling_Parameters * parameters, int all_threads, int thread_no)
{
	char line_buffer [3000];


	FILE *fp = fopen(in_FASTA_file,"r");

	if(!fp)
	{
		SUBREADprintf("Referenced Genome not found or inaccessible: '%s'.\n", in_FASTA_file);
		return -1;
	}
	
	while(!feof(fp))
	{
		int linelen = read_line(2999, fp, line_buffer, 0);
		if(line_buffer[0] == '>')
		{
			char chro_name [MAX_CHROMOSOME_NAME_LEN];
			int i;
			for(i=0; i< linelen -1; i++)
			{
				if(line_buffer[i+1] == ' ' || line_buffer[i+1] == '|' || line_buffer[i+1] == '\t')
					break;
				else chro_name[i] = line_buffer[i+1];
			}
			chro_name[i]=0;
			if(run_chromosome_search(fp, out_fp, chro_name , temp_prefix, chromosomes, parameters)) return -1;
		}
	}
	fclose(fp);

	return 0;
}

int parse_read_lists_maybe_threads(char * in_FASTA_file, char * out_BED_file, char * temp_prefix, chromosome_t * chromosomes, struct SNP_Calling_Parameters* parameters, int all_threads)
{
	FILE * out_fp = fopen(out_BED_file,"w");
	int ret;
	if(!out_fp){
		SUBREADprintf("Cannot open the output file: '%s'\n", out_BED_file);
	}
	fputs("chr\tpos\tref\talt\tfreq\tdepth\n", out_fp);
	ret =  parse_read_lists(in_FASTA_file, out_fp, temp_prefix, chromosomes, parameters , all_threads, 0);
	fprintf(out_fp, "## Fisher_Test_Size=%u\n",fisher_test_size);
	fclose(out_fp);
	return ret;
}

int SNP_calling(char * in_SAM_file, char * out_BED_file, char * in_FASTA_file, char * temp_location, unsigned int known_read_count, int threads, struct SNP_Calling_Parameters* parameters)
{
	char temp_file_prefix[300];
	chromosome_t * known_chromosomes;
	unsigned int real_read_count=0;

	double start_time = miltime();
	unsigned short rand48_seed[3];

	int i, fpos;

	fisher_test_size = 0;

	precalculated_factorial = (double*)malloc(sizeof(double)*PRECALCULATE_FACTORIAL);
	for(i=0; i<PRECALCULATE_FACTORIAL; i++)
		precalculated_factorial[i] = 0.; 
		

	known_chromosomes = (chromosome_t *) SUBREAD_malloc(sizeof(chromosome_t) * OFFSET_TABLE_SIZE);
	if(!known_chromosomes)
	{
		fatal_memory_size();
		return -1;
	}

	known_chromosomes[0].chromosome_name[0]=0;

	SUBREADprintf("Spliting SAM file\n");
	fflush(stdout);
	//Step 1:The SAM file is scanned to create a number of temp files "temp-snps-chrX-21-00000000-XXXXXX" and the related read positions/sequences are written into them (in 2-bit base coding). A read can contribute to two blocks if it crosses the border of the blocks. If there are indels in a read, each continuously mapped section is individually written into the temporary file. The quality scores are written companying the bases. The hierarchy of the data is: block -> read -> sections {bases, phred scores, CIGAR string, reported mapping quality}
	memcpy(rand48_seed, &start_time, 6);
	seed48(rand48_seed);
	sprintf(temp_file_prefix, "%s/temp-snps-%06u-%08lX-", temp_location==NULL?".":temp_location, getpid(), lrand48());


	fpos=0;
	while(1)
	{
		int fpos0 = fpos;
		char one_fn [300];
		while(in_SAM_file[fpos]!=',' && in_SAM_file[fpos]!=0)
			fpos++;
		strncpy(one_fn, in_SAM_file+fpos0, fpos-fpos0);
		one_fn[fpos-fpos0]=0;

		if(break_SAM_file(one_fn, temp_file_prefix, &real_read_count, known_chromosomes, 1 /* This 0 means that the sequence/quality/cigar fields are needed in the temp files*/, parameters -> bases_ignored_head_tail)) return -1;

		if(!in_SAM_file[fpos]) break;
		fpos++;
	}

	//Step 2:Each base blocks are load from the FASTA file into memory, then each temp file is scanned to create the SNP voting table. The sections in the temp files are scanned to create the SNP voting table. The voting table is compared with the sequences in the FASTA file to determine if each base is a SNP. The result is written into the bed file immediatly.
	if(parse_read_lists_maybe_threads(in_FASTA_file, out_BED_file, temp_file_prefix, known_chromosomes, parameters, threads)) return -1;

	free(known_chromosomes);

	free(precalculated_factorial);
	SUBREADprintf("Finished.\n");

	return 0;
}


void print_usage_snp(char * myname)
{
	SUBREADprintf("Usage: %s -i <input_SAM_file> -g <single_genome_FASTA_file> -o <output_BED_file> {-r threshold (float)} {-t temp_path} {-c max_read_number} {-n read_threshold_min} {-m read_threshold_max} {-q min_phred_score} {-N neighbour_test_len,matching_rate} {-S} {-p Len,Pvalue for Fisher's Exact Test} {-I number of first/last bases ignored from each read}\n\t-q\t ignoring low-quality bases in reads\n\t-N\t using neighbour filtering to remove wrongly mapped read segments\n", myname);
}



#ifdef MAKE_STANDALONE
int main(int argc,char ** argv)
#else
int main_snp_calling_test(int argc,char ** argv)
#endif
{
	char c;
	char in_SAM_file[5000];
	char out_BED_file[300];
	char temp_path[300];
	char in_FASTA_file[300];
	int threads;
	int t=0, k;
	struct SNP_Calling_Parameters parameters;
	unsigned int read_count;

	in_SAM_file [0] = 0;
	out_BED_file[0] = 0;
	temp_path[0] = 0;

	read_count = 0;
	threads = 0;
	parameters.supporting_read_rate = 0.5;
	parameters.min_supporting_read_number = 5;
	parameters.max_supporting_read_number = 10000;
	parameters.neighbour_filter_testlen = 0;
	parameters.neighbour_filter_rate = 0.5;
	parameters.min_phred_score = 0;
	parameters.fisher_exact_p_threshold = -1;
	parameters.fisher_exact_testlen = 0;
	parameters.test_two_strands = 0;

	parameters.bases_ignored_head_tail = 0;

	if(argc<2)
	{
		print_usage_snp(argv[0]);
		return 0;
	}
	while ((c = getopt (argc, argv, "N:I:p:q:i:o:r:t:g:n:c:m:?S")) != -1)
	{
		switch (c)
		{
			case 'I':
				parameters.bases_ignored_head_tail = atoi(optarg);
				break;
			case 'S':
				parameters.test_two_strands = 1;
				break;

			case 'p':
				if(parameters.neighbour_filter_testlen > 0)
				{
					SUBREADprintf("You cannot use both neighbour filtering and Fisher's exact test.\n");
					return -1;
				}

				k=strlen(optarg);
				for(t=0;t<k;t++)
					if(optarg[t]==',')
					{
						optarg[t]=0;
						break;
					}

				if(t==k)
					SUBREADprintf("Warning: the Fisher's exact test parameter is unparseable. It should be like \"-p 5,0.5\".\n");
				else
				{
					parameters.fisher_exact_testlen = atoi(optarg);
					parameters.fisher_exact_p_threshold = atof(optarg+t+1);
				}
				break;

			case 'q':
				parameters.min_phred_score = atoi(optarg);
				break;

			case 'N':

				if(parameters.fisher_exact_p_threshold >= -0.00001)
				{
					SUBREADprintf("You cannot use both neighbour filtering and Fisher's exact test.\n");
					return -1;
				}
	

				k=strlen(optarg);
				for(t=0;t<k;t++)
					if(optarg[t]==',')
					{
						optarg[t]=0;
						break;
					}
				if(t==k)
					SUBREADprintf("Warning: the neighbour filtering parameter is unparseable. It should be like \"-N 5,0.5\".\n");
				else
				{
					parameters.neighbour_filter_testlen = atoi(optarg);
					parameters.neighbour_filter_rate = atof(optarg+t+1);
				}
				break;

			case 'g':
				strncpy(in_FASTA_file, optarg,299);
				break;

			case 'i':
				strncpy(in_SAM_file, optarg,299);
				break;

			case 'o':
				strncpy(out_BED_file, optarg,299);
				break;

			case 't':
				strncpy(temp_path,  optarg,299);
				break;

			case 'T':
				threads = atoi(optarg);
				break;

			case 'r':
				parameters.supporting_read_rate = atof(optarg);
				break;

			case 'm':
				parameters.max_supporting_read_number = atof(optarg);
				break;

			case 'n':
				parameters.min_supporting_read_number = atof(optarg);
				break;

			case 'c':
				read_count = atoi(optarg);
				break;

			case '?':
			default:
				print_usage_snp(argv[0]);
		}
	}

	SUBREADprintf("Parameters: \nneighbour_filter_testlen = %d\nneighbour_filter_rate = %.4f\nmin_supporting_read_number = %d\nsupporting_read_rate = %.4f\nmin_phred_score = %d\nFisher's exact test len = %d\nFisher's exact test maximum p-valie = %.5f\n\n", parameters.neighbour_filter_testlen , parameters.neighbour_filter_rate, parameters.min_supporting_read_number, parameters.supporting_read_rate, parameters.min_phred_score, parameters.fisher_exact_testlen , parameters.fisher_exact_p_threshold);
	
	return SNP_calling(in_SAM_file, out_BED_file, in_FASTA_file, temp_path[0]?temp_path:NULL, read_count, threads, &parameters);
	
}
