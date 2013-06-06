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
  
  
#include <stdlib.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <getopt.h>
#include "subread.h"
#include "gene-value-index.h"
#include "input-files.h"

float MIN_REPORTING_RATIO = 0.8;
gene_offset_t _global_offsets;
unsigned int  SCAN_TOTAL_BASES=0;
int ic = 0;


void fullscan_usage()
{
	SUBREADputs("This program is for scanning a particular read in the whole chromosome space by completely scanning the index file.");

	SUBREADputs("Usage: subread-fullscan -i <index_base> -m <report_threshold> <read_string>");
}

void report_pos(unsigned int pos)
{
	char * chro_name;
	unsigned int chro_pos;
	locate_gene_position(pos, &_global_offsets, &chro_name, &chro_pos);
	SUBREADprintf ("%s,%u\n", chro_name,chro_pos);
}


int str_match_count(char * c1, char * c2, int rl)
{
	int i,ret =0;
	for (i=0; i<rl;i++)
		ret += (c1[i]==c2[i]);
	return ret;
}

void scan_test_match(char * read, char * read_rev, char * chro, int rl, unsigned int pos)
{
	int m = str_match_count(read, chro, rl);
	int mr = str_match_count(read_rev, chro, rl);
	int threshold = (int)(MIN_REPORTING_RATIO * rl - 0.001);
	if (m>=threshold)
	{
		SUBREADprintf("\nFound on positive strand (%0.2f%%): ", m*100./rl);
		report_pos(pos);
	}
	if (mr>=threshold)
	{
		SUBREADprintf("\nFound on negative strand (%0.2f%%): ", mr*100./rl);
		report_pos(pos);
	}

	//if (pos > 19999999)
	//SUBREADprintf ("m=%d  mr=%d  T=%d\n", m, mr, threshold);
}

void full_scan_read(char * index_name, char * read_str)
{
	int read_len = strlen(read_str);
	char read_rev_str[1208];
	char chro_str[1208];
	gene_value_index_t index;
	int tabno = 0;
	unsigned int current_pos =0xffffffffu;

	strcpy(read_rev_str, read_str);
	reverse_read(read_rev_str, read_len, GENE_SPACE_BASE);

	while (1){
		char table_fn[1250];
		struct stat filestat;
		sprintf(table_fn, "%s.%02d.b.array", index_name, tabno);

		int stat_ret = stat(table_fn, &filestat);
		if (stat_ret !=0 )
		{
			if(!tabno)
				SUBREADprintf("The index does not contain any raw base data which is required in scanning. Please use the -b option while building the index.\n");
			return ;
		}
		if (tabno>0)gvindex_destory(&index);

		gvindex_load(&index,table_fn);
		if (tabno==0)
		{
			gvindex_get_string (chro_str, &index, 0, read_len, 0);
			current_pos = 0;
		}

		for(; current_pos + read_len < index.start_point + index.length ; current_pos++)
		{
			scan_test_match(read_str, read_rev_str, chro_str, read_len, current_pos);
			char nch = gvindex_get(&index, current_pos + read_len);
			int i;
			for (i=0; i<read_len-1; i++)
				chro_str[i]= chro_str[i+1];
			chro_str[read_len-1] = nch;
			if (current_pos % 1000000 == 0)
				print_text_scrolling_bar("Scan:", current_pos*1./SCAN_TOTAL_BASES, 80, &ic);
		}
		tabno +=1;
	}
	if (tabno>0)gvindex_destory(&index);


}

int main (int argc , char ** argv)
{
	char index_name [1200];
	char read_str [1208];
	char c;
	int i;

	index_name[0]=0;

	while ((c = getopt (argc, argv, "i:m:?")) != -1)
		switch(c)
		{
			case 'i':
				strncpy(index_name,  optarg, 1199);
				break;
			case 'm':
				MIN_REPORTING_RATIO = atof(optarg);
				break;
			case '?':
				return -1 ;
		}

	if (!index_name[0] || optind == argc)
	{
		fullscan_usage();
		return -1;
	}
	i=0;
	while(argv[optind][i])
	{
		argv[optind][i] = toupper(argv[optind][i]);
		i++;
	}
	strncpy(read_str, argv[optind], 1199);

	load_offsets (&_global_offsets, index_name);
	SUBREADprintf ("Reporting threshold=%0.2f%%\n", MIN_REPORTING_RATIO*100);
	for(i=0;i<_global_offsets.total_offsets;i++)
	{
		if (!_global_offsets.read_offsets[i])break;
		SCAN_TOTAL_BASES = _global_offsets.read_offsets[i];
	}
	SUBREADprintf ("All bases =%u\n", SCAN_TOTAL_BASES);
	SUBREADprintf ("Scanning the full index for %s...\n\n", read_str);

	full_scan_read(index_name, read_str);

	SUBREADprintf ("\nFinished.\n");
	destroy_offsets(&_global_offsets);

	return 0;
}
