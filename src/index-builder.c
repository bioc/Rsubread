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
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include "hashtable.h"
#include "gene-value-index.h"
#include "gene-algorithms.h"
#include "sorted-hashtable.h"
#include "input-files.h"

#define NO_GENE_DEBUG_
#define _GENE_DEBUG_SIZE_ 40000
#define MIN_READ_SPLICING 2000000


#define MAX_KEY_MATCH GENE_VOTE_SPACE 

int IS_COLOR_SPACE = 0;
int VALUE_ARRAY_INDEX = 1;
int QUICK_BUILD = 0;

#define NEXT_READ 1
#define NEXT_FILE 2
#define FULL_PARTITION 4


void print_build_log(double finished_rate, double read_per_second, double expected_seconds, unsigned long long int total_reads)
{
        char outbuff[81]; int i;
        snprintf(outbuff, 80,"completed=%0.2f%%; time used=%.1fs; rate=%.1fk bps/s; total=%llum bps", finished_rate*100, miltime()-begin_ftime, read_per_second/1000 ,total_reads/1000000);
        SUBREADprintf("%s",outbuff);
        for(i=strlen(outbuff); i<105; i++)
                SUBREADprintf(" ");
        SUBREADprintf("\r");
}

#define MAX_BASES_IN_INDEX 4294900000.0

int build_gene_index(const char index_prefix [], char ** chro_files, int chro_file_number, unsigned int memory_megabytes, int threshold, gehash_t * huge_table)
{
	int file_number, table_no;
	int status = NEXT_FILE;
	unsigned int offset, read_no;
	unsigned int segment_size = (int)(memory_megabytes * 1024.0 / 9.15) * 1024 ;
	long long int all_bases = guess_gene_bases(chro_files,chro_file_number);
	char fn[300];
	double local_begin_ftime = 0.;

	int chro_table_maxsize=100;
	unsigned int * read_offsets = malloc(sizeof(unsigned int) * chro_table_maxsize);
	char * read_names = malloc(MAX_READ_NAME_LEN * chro_table_maxsize);
	gehash_t table;
//	gehash_t huge_table;
	gene_value_index_t value_array_index;

	gene_input_t ginp;

	SUBREADprintf ("Index items per partition = %u\n\n", segment_size);



	if (chro_file_number > 199)
	{
		SUBREADprintf("There too many chromosome files. You may merge them into less than 199 files.\n");
		return -1;
	}

	if (strlen (index_prefix) > 290)
	{
		SUBREADprintf("The path is too long. It should not be longer than 290 chars.\n");
		return -1;
	}

	if(all_bases<1)
	{
		SUBREADprintf("File '%s' is inaccessible.\n", chro_files[-all_bases-1]);
		return -1;
	}

	if(gehash_create(& table, segment_size, 0)) return 1;

	unsigned int size_of_array_index = (unsigned int)(min(MAX_BASES_IN_INDEX, segment_size*4.35 + MIN_READ_SPLICING));

	if(VALUE_ARRAY_INDEX)
		if(gvindex_init(&value_array_index, 0, size_of_array_index)) return 1;

	file_number = 0;
	offset = 0;
	table_no = 0;
	read_no = 0;

	bzero(read_offsets, chro_table_maxsize*sizeof(int));
	sprintf(fn, "%s.files", index_prefix);
	unlink(fn);

	status = NEXT_FILE;

	SUBREADprintf("\nBuilding the index...\n");

	{
		char window [16], last_color_base=-1, last_last_color_base=-1;
		int i, read_len = 0;
		unsigned int int_key = 0, array_int_key = 0;
		int skips=0, all_skips = 0;

		//Pre-fill
		while(1)
		{
			//Subread Cycle
			char next_char;

			if (status == NEXT_FILE)
			{
				if(file_number == chro_file_number)
				{
					FILE * fp;

					geinput_close(&ginp);

					SUBREADprintf ("\n Processing chromosome files ...\n");

					sprintf (fn, "%s.%02d.%c.tab", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
					SUBREADfflush(stdout);

					gehash_dump(&table, fn);

					if(VALUE_ARRAY_INDEX)
					{
						sprintf (fn, "%s.%02d.%c.array", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
						gvindex_dump(&value_array_index, fn);
						gvindex_destory(&value_array_index) ;
					}



					gehash_destory(&table);

					read_offsets[read_no-1] = offset;

					for(i=table_no+1; i<100; i++)
					{
						sprintf(fn, "%s.%02d.%c.tab", index_prefix, i, IS_COLOR_SPACE?'c':'b');
						unlink(fn);
						sprintf(fn, "%s.%02d.%c.array", index_prefix, i, IS_COLOR_SPACE?'c':'b');
						unlink(fn);
					}

					sprintf (fn, "%s.reads", index_prefix);
					fp = fopen(fn, "w");
					for (i=0; i<read_no; i++)
						fprintf(fp, "%u\t%s\n", read_offsets[i], read_names+i*MAX_READ_NAME_LEN);

					fclose(fp);

					break;
				}
				else
				{
					if (file_number)
						geinput_close(&ginp);
					geinput_open(chro_files[file_number++], &ginp);
					status = NEXT_READ;
				}
			}
			if (status == NEXT_READ)
			{

				geinput_readline(&ginp, fn, 0);

				if(read_no>0)
					read_offsets[read_no-1] = offset;

				for(i=0;(fn[i+1] != ' ' && fn[i+1] != '\0' && fn[i+1] != '\t' && i<47); i++)
					*(read_names + MAX_READ_NAME_LEN*read_no + i) = fn[i+1];

				*(read_names + MAX_READ_NAME_LEN*read_no + i) = 0;

				sprintf(fn, "%s.files", index_prefix);
				FILE * fname_fp = fopen(fn, "a");
				fprintf(fname_fp, "%s\t%s\t%ld\n", read_names+read_no*MAX_READ_NAME_LEN, ginp.filename, ftell(ginp.input_fp));
				fclose(fname_fp);
				

				for (i=0; i<16; i++)
				{
					char nch = geinput_next_char(&ginp);
					if (nch == 'N') skips = 16;
					else if (skips>0) skips--;
					window[i] = nch;
				}
				read_len = 16;
				read_no ++;

				if(read_no >= chro_table_maxsize)
				{
					read_offsets = realloc(read_offsets, 2* chro_table_maxsize * sizeof(unsigned int));
					read_names = realloc(read_names, 2* chro_table_maxsize * MAX_READ_NAME_LEN);
					chro_table_maxsize *= 2;
				}

				if(IS_COLOR_SPACE)
				{
					int_key = genekey2color('A',window);
					if(VALUE_ARRAY_INDEX)
						array_int_key = genekey2int(window,GENE_SPACE_BASE);
					last_last_color_base = -1;
					last_color_base = window[15];
				}
				else
					array_int_key = int_key = genekey2int(window,GENE_SPACE_BASE);
	
				if(offset==0)
					local_begin_ftime = miltime();
			}
			if (status == FULL_PARTITION) 
			{
				int seek_back_reads ;

				if (read_len < MIN_READ_SPLICING)
				{
					seek_back_reads = read_len;
					offset -= seek_back_reads;
					offset += 16;
				}
				else
				{
					seek_back_reads = MIN_READ_SPLICING - 10;
					seek_back_reads -= seek_back_reads%3;
					seek_back_reads += 1;
					offset -= seek_back_reads;
					offset += 16;
				}

				sprintf(fn, "%s.%02d.%c.tab", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
				SUBREADfflush(stdout);

				gehash_dump(&table, fn);
				if(VALUE_ARRAY_INDEX)
				{
					sprintf(fn, "%s.%02d.%c.array", index_prefix, table_no, IS_COLOR_SPACE?'c':'b');
					gvindex_dump(&value_array_index, fn);
				}

				table_no ++;

				gehash_destory(&table);
				if(VALUE_ARRAY_INDEX)
					gvindex_destory(&value_array_index);

				read_len -= seek_back_reads;
				read_len += 16;

				while(seek_back_reads)
				{
					fseek(ginp.input_fp, -2, SEEK_CUR);
					char bnch = fgetc(ginp.input_fp);
					if ((bnch >='A' && bnch <= 'Z' ) || (bnch >='a' && bnch <= 'z' ) || bnch=='.')seek_back_reads--;
				}

				for (i=0; i<16; i++)
				{
					char nch = geinput_next_char(&ginp);
					if (nch == 'N' ) skips = 16;
					else if (skips>0) skips--;
					window[i] = nch;
				}
			
				if(IS_COLOR_SPACE)
				{
					int_key = genekey2color('A',window);
					if(VALUE_ARRAY_INDEX)
						array_int_key = genekey2int(window,GENE_SPACE_BASE);
					last_last_color_base = -1;
					last_color_base = window[15];
				}
				else
					array_int_key = int_key = genekey2int(window, GENE_SPACE_BASE);
				
				if(gehash_create(&table, segment_size, 0)) return 1;
				if(VALUE_ARRAY_INDEX)
					if(gvindex_init(&value_array_index, offset - (IS_COLOR_SPACE?0:0),(unsigned int)(min(MAX_BASES_IN_INDEX-offset + 2, size_of_array_index )))) return 1;
			}
	
			status = 0;

			if(skips || (IS_COLOR_SPACE && last_last_color_base<0))
				all_skips ++;
			else
			{
				int is_no_info = gehash_exist(huge_table, int_key);
				if(!is_no_info)
				{
					if(gehash_insert(&table, int_key, offset - (IS_COLOR_SPACE?1:0))) return 1;
				}
				if(VALUE_ARRAY_INDEX)
				{
					gvindex_set(&value_array_index, offset - (IS_COLOR_SPACE?0:0), array_int_key);
				}
			}

			if(table.current_items >= segment_size && (read_len > MIN_READ_SPLICING || read_len < 32))
			{
				status = FULL_PARTITION;
				continue;
			}

			for (i=0; i<GENE_SLIDING_STEP; i++)
			{
				next_char = geinput_next_char(&ginp);
				if(next_char < 0)
				{
					if (next_char == -1) status = NEXT_READ;
					if (next_char == -2) status = NEXT_FILE;
					if (next_char == -3) return 0;
					break;
				}

				if (next_char == 'N' )skips = 16;
				else if (skips>0){
					skips--;
					last_color_base = -1;
				}

				int_key = int_key << 2;


				if (IS_COLOR_SPACE)
				{
					if(last_color_base>0)
						int_key += chars2color(last_color_base, next_char);
					if(VALUE_ARRAY_INDEX)
					{
						array_int_key = array_int_key << 2;
						array_int_key += base2int (next_char);
					}

					last_last_color_base = last_color_base;
					last_color_base = next_char;
				}
				else
				{
					int_key += base2int (next_char); 
					array_int_key = int_key;
				}

				if (offset % 5000000 == 0)
				{
					if (offset>1)
					{
						double finished_rate = offset*1.0/all_bases;
						double base_per_second = offset / (miltime()-local_begin_ftime);
						double ETA = (all_bases-offset) / base_per_second;
						print_build_log(finished_rate,base_per_second,ETA, all_bases);
					}
					SUBREADfflush(stdout) ;
				}

				offset ++;
				read_len ++;

				if(offset > 0xFFFFFFFDU)	
				{
					SUBREADprintf("The chromosome data contains too many bases. The size of chromosome file should be less than 4G Bytes\n") ;
					return -1;
				}

			}

		}
	}

	free(read_names);
	free(read_offsets);
	return 0;
}

int add_repeated_subread(gehash_t * tab , unsigned int subr, unsigned char ** huge_index)
{
	unsigned int times;

	int huge_byte = (subr>>2) &0x3fffffff;
	int huge_offset = (subr % 4) * 2;
	unsigned int byte_value = huge_index[ (huge_byte >> 20) & 1023 ][huge_byte&0xfffff] ;

	int huge_value = (byte_value>> huge_offset) & 0x3;
	if(huge_value <3)
	{
		huge_value ++;
		huge_index[ (huge_byte >> 20) & 1023 ][huge_byte&0xfffff] = (byte_value & (~(0x3 << huge_offset))) | (huge_value << huge_offset);
		return 0;
	}

	int matched = gehash_get(tab, subr, &times, 1);
	if(matched)
	{
		gehash_update(tab, subr, times+1);
	}
	else
		if(gehash_insert(tab, subr,4)) return 1;
	return 0;
}


int scan_gene_index(const char index_prefix [], char ** chro_files, int chro_file_number, int threshold, gehash_t *huge_table)
{
	int file_number, table_no, i ,j;
	int status = NEXT_FILE;
	unsigned int offset, read_no;
	char fn[300];
	double local_begin_ftime = miltime();
	long long int all_bases = guess_gene_bases(chro_files,chro_file_number);

	gehash_t occurance_table;
	unsigned char * huge_index[1024];

	for(i=0;i<1024;i++)
	{
		huge_index[i] = (unsigned char *)malloc(1024*1024); 
		if(!huge_index[i])
		{
			for(j=0;j<i;j++) free(huge_index[j]);
			SUBREADprintf("You need at least one point five gigabytes of memory for building the index.\n");
			return -1;
		}
		memset(huge_index[i], 0 , 1024*1024);
	}


	if(gehash_create(&occurance_table , 100000000, 0)) return 1;


	gene_input_t ginp;

	SUBREADprintf("Scanning non-informative reads in the chromosomes...\n");

	if (chro_file_number > 199)
	{
		SUBREADprintf("There too many chromosome files. You may merge them into less than 199 files.\n");
		return -1;
	}

	if (strlen (index_prefix) > 290)
	{
		SUBREADprintf("The path is too long. It should not be longer than 290 chars.\n");
		return -1;
	}
	if(all_bases<1)
	{
		SUBREADprintf("File '%s' is inaccessible.\n", chro_files[-all_bases-1]);
		return -1;
	}


	file_number = 0;
	offset = 0;
	table_no = 0;
	read_no = 0;

	sprintf(fn, "%s.files", index_prefix);
	unlink(fn);

	status = NEXT_FILE;


	{
		char window [16], last_color_base=-1, last_last_color_base=-1;
		int i, read_len = 0;
		unsigned int int_key = 0, array_int_key = 0;
		int skips=0, all_skips = 0;

		//Pre-fill
		while(1)
		{
			//Subread Cycle
			char next_char;

			if (status == NEXT_FILE)
			{
				if(file_number == chro_file_number)
				{
					geinput_close(&ginp);

					break;
				}
				else
				{
					if (file_number)
						geinput_close(&ginp);
					geinput_open(chro_files[file_number++], &ginp);
					status = NEXT_READ;
				}
			}
			if (status == NEXT_READ)
			{

				geinput_readline(&ginp, fn, 0);

				for (i=0; i<16; i++)
				{
					char nch = geinput_next_char(&ginp);
					if (nch == 'N') skips = 16;
					else if (skips>0) skips--;
					window[i] = nch;
				}
				read_len = 16;
				read_no ++;

				if(IS_COLOR_SPACE)
				{
					int_key = genekey2color('A',window);
					if(VALUE_ARRAY_INDEX)
						array_int_key = genekey2int(window,GENE_SPACE_BASE);
					last_last_color_base = -1;
					last_color_base = window[15];
				}
				else
					array_int_key = int_key = genekey2int(window,GENE_SPACE_BASE);
	
			}
	
			status = 0;

			if(skips || (IS_COLOR_SPACE && last_last_color_base<0))
				all_skips ++;
			else
			{
				add_repeated_subread(&occurance_table, int_key, huge_index);
			}


			for (i=0; i<GENE_SLIDING_STEP; i++)
			{
				next_char = geinput_next_char(&ginp);
				if(next_char < 0)
				{
					if (next_char == -1) status = NEXT_READ;
					if (next_char == -2) status = NEXT_FILE;
					if (next_char == -3) return 0;
					break;
				}

				if (next_char == 'N' )skips = 16;
				else if (skips>0){
					skips--;
					last_color_base = -1;
				}

				int_key = int_key << 2;


				if (IS_COLOR_SPACE)
				{
					if(last_color_base>0)
						int_key += chars2color(last_color_base, next_char);
					if(VALUE_ARRAY_INDEX)
					{
						array_int_key = array_int_key << 2;
						array_int_key += base2int (next_char);
					}

					last_last_color_base = last_color_base;
					last_color_base = next_char;
				}
				else
				{
					int_key += base2int (next_char); 
					array_int_key = int_key;
				}

				offset ++;
				read_len ++;


				if (offset % 5000000 == 0)
				{
					if (offset>1)
					{
						double finished_rate = offset*1.0/all_bases;
						double base_per_second = offset / (miltime()-local_begin_ftime);
						double ETA = (all_bases-offset) / base_per_second;
						print_build_log(finished_rate,base_per_second,ETA, all_bases);
					}
					SUBREADfflush(stdout) ;
				}

				if(offset > 0xFFFFFFFDU)	
				{
					SUBREADprintf("The chromosome data contains too many bases. The size of chromosome file should be less than 4G Bytes\n") ;
					return -1;
				}

			}

		}
	}





	for (i=0; i<occurance_table.buckets_number; i++)
	{
		struct gehash_bucket * current_bucket = &(occurance_table.buckets[i]);

		if(current_bucket -> current_items>=1)
		{
			for (j=0; j<current_bucket -> current_items; j++)
			{
				if(current_bucket -> item_values [j] > threshold)
				{
					if(gehash_insert(huge_table, current_bucket -> item_keys[j], 1)) return 1;
				}
			}
		}
	}

	for(i=0;i<1024;i++)
		free(huge_index[i]);
	gehash_destory(&occurance_table);


	SUBREADprintf("There are %llu non-informative subreads found in the chromosomes.\n" , huge_table -> current_items);
	return 0;
}




#ifdef MAKE_STANDALONE
int main(int argc,char ** argv)
#else
int main_buildindex(int argc,char ** argv)
#endif
{
	int threshold = 24;
	int memory_limit = 8000;	// 8000 MBytes
	char output_file[300], c;
	output_file[0] = 0;

	SUBREADprintf("\n");
	while ((c = getopt (argc, argv, "vcqM:o:f:D?")) != -1)
		switch(c)
		{
			case 'v':
				print_version_info();
				return 0;
				break;
			case 'q':
				QUICK_BUILD = 1;
				break;
			case 'c':
				IS_COLOR_SPACE = 1;
				break;
			case 'M':
				memory_limit = atoi(optarg);
				break;
			case 'f':
				threshold = atoi(optarg);
				break;
			case 'o':
				strncpy(output_file, optarg, 299);
				break;
			case '?':
				return -1 ;
		}

	if (argc == optind || !output_file[0])
	{

                SUBREADprintf ("Usage:\n %s -o <basename> -M <int> {FASTA file1} [FASTA file2] ...\n\nArguments:\n    -o <basename>\t base name of the index to be created\n    -M <int>\t\t size of requested memory(RAM) in megabytes, 8000 by default\n    -f <int>\t\t specify the threshold for removing uninformative subreads (highly repetitive 16mers in the reference). 24 by default.\n    -c      \t\t build a color-space index\n    -v      \t\t display the version number.\n\nExample:\n %s -o my_index chr1.fa chr2.fa ...\n\n", argv[0], argv[0]);

		return -1 ;
	}


	if(IS_COLOR_SPACE)
		SUBREADprintf("Building a color-space index.\n");
	else
		SUBREADprintf("Building a base-space index.\n");

	SUBREADprintf("Size of memory used=%d MB\n",memory_limit);
	SUBREADprintf("Base name of the built index = %s\n",output_file);

	SUBREADfflush(stdout);
	begin_ftime = miltime();


	gehash_t huge_table;
	gehash_create(& huge_table, 50000000, 0);
	if(scan_gene_index(output_file, argv+optind , argc - optind, threshold, &huge_table)) return 1;
	if(build_gene_index(output_file, argv+optind , argc - optind,  memory_limit, threshold, &huge_table)) return 1;
	SUBREADprintf("\nIndex %s was successfully built.\n", output_file);

	return 0;
}

