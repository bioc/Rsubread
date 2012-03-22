/*
 * remove_duplicate_reads.c
 *
 *  Created on: Jan 13, 2012
 *      Author: dai
 *  Modified by Wei Shi, March 21, 2012.
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>


/*	Constants for testing purpose but are variables for this function */

#define MAXCHRNUM 150
#define NUMREGCHR 24
#define STR 100	//maximum string length
#define MAXCHRLENGTH 250000000
#define LEN 1000



/* human chromosomes */
char *reg_chrs[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
				"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
				"chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};


/* Global Variable */
int debug = 0;	// whether output program progress to user	
int detail = 0;	// whether output duplicate removal result on chromosome basis
char *sam_file;
char *sam_header;
char *output_sam_file;
char *remove_logfile;
int count_seq = 0;
int32_t duplicate_threshold=0;
int32_t num_chr_reads = 0;
int32_t num_input_reads = 0;
int32_t num_output_reads = 0;
int number_of_site = 0;

char *chrs[MAXCHRNUM];
char *chr_reads_filenames[MAXCHRNUM];
FILE *chr_reads_files[MAXCHRNUM];
int32_t chr_length[MAXCHRNUM];
int site_removed[MAXCHRNUM];
int num_read[MAXCHRNUM];
int num_read_removed[MAXCHRNUM];

int depth[MAXCHRLENGTH];
int mapq[MAXCHRLENGTH];

/*
 *  Given the chromosome name(string), return chrosome order(int)
 *  For unidentified human chromosome, it returns -1
 * */
int find_chr_index(char *read_chr){
	int i = 1;
	while ((i <= count_seq) && (strcmp(read_chr, chrs[i]) != 0)){
		i++;
	}
	if (i <= count_seq){
		return i;
	} else {
		return -1;
	}
}



/*
 * For the input read mapping file(.sam file)
 * Group read according to their mapped chromosome
 * Each chromosome has its own read file
 * */
void group_reads(void){
	if (debug == 1) {
		printf("In grouping read ...\n");
		printf("Number of sequence: %d\n", count_seq);
	}
	
	
	int i;	//control variable

	for (i=1; i<=count_seq; i++){
		chr_reads_files[i] = fopen(chr_reads_filenames[i], "w");
	}
	

	/* Group reads into subgroups */
	FILE *fallreads;
	FILE *fout;
	char * line = NULL;
	char * line_copy = NULL;
	char * read_chr;
	int read_chr_index;
	int32_t pos;
	char *cigar;
	char *read_seq;

	char *readline;

	fallreads = fopen(sam_file,"r");
	fout = fopen(output_sam_file, "w");

	line = (char *)calloc(LEN+1,sizeof(char));
	line_copy = (char *)calloc(LEN+1,sizeof(char));
	
	while ((readline = fgets(line, LEN, fallreads)) != NULL){
		strcpy(line_copy, line);
		if(line[0] == '@'){
			fprintf(fout, "%s", line);
			continue;
		}
		num_input_reads++;
		strtok(line,"\t");
		strtok(NULL,"\t");
		read_chr = strtok(NULL,"\t");
		pos = atoi(strtok(NULL,"\t"));
		strtok(NULL,"\t");		
		cigar = strtok(NULL,"\t");
		strtok(NULL,"\t");
		strtok(NULL,"\t");
		strtok(NULL,"\t");
		read_seq = strtok(NULL, "\t");
		if( ((read_chr_index = find_chr_index(read_chr)) == -1) ){
			if ((read_chr[0] != '*') && (debug == 1)){
				printf("Attention unrecognized reference %s", read_chr);
			}
			fprintf(fout, "%s", line_copy);
			continue;
		} else {
			num_chr_reads++;
		}
		fprintf(chr_reads_files[read_chr_index], "%s", line_copy);
	}
	/* Close input file and read files for each chromosome */
  	fclose(fallreads);
  	fclose(fout);
	for (i=1; i<= count_seq; i++){
		fclose(chr_reads_files[i]);
	}
	
	if (debug == 1){
		printf("Finished grouping read ...\n");
	}
	
	/* If threshold not present, calculate as two times exonic base coverage*/
	if (duplicate_threshold == 0){
		double threshold = 2.0 * num_chr_reads * 49 / 0.02 / 3000000000;
		duplicate_threshold = (int)(threshold+0.5);
	}
	printf("Number of reads in input file: %d\n", num_input_reads);
	printf("Number of reads mapped to given reference: %d\n", num_chr_reads);
	printf("Mapping to reference rate is : %f\n\n", 1.0 * num_chr_reads/num_input_reads);
	printf("Maximal number of duplicates allowed is: %d\n", duplicate_threshold);
}



void filter_one_chr(int chr_index){
	if (debug == 1){
		printf("in processing chr:%s\n", chrs[chr_index]);
	}

	int i;
	for(i=0; i<= chr_length[chr_index]; i++){
		depth[i] = 0;
		mapq[i] = 0;
	}

	FILE *fin;
	fin = fopen(chr_reads_filenames[chr_index],"r");
	char * line = NULL;
	char * read_chr;
	int32_t pos;
	int mapqual;
	char *readline;

	line = (char *)calloc(LEN+1,sizeof(char));

	while ((readline = fgets(line, LEN, fin)) != NULL){
		if(line[0] == '@')
			continue;
		num_read[chr_index]++;
		strtok(line,"\t");
		strtok(NULL,"\t");
		read_chr = strtok(NULL,"\t");
		pos = atoi(strtok(NULL,"\t"));
		mapqual = atoi(strtok(NULL, "\t"));
		if ((mapqual > 255) || (mapqual <=0)){
			//printf("Is there something wrong with mapq? %d\n", mapqual);
		}
		depth[pos]++;
		if (mapqual > mapq[pos]){
			mapq[pos] = mapqual;
		}
	}
	fclose(fin);
	
	FILE *flog;
	flog = fopen(remove_logfile, "a");
	for(i=0; i<= chr_length[chr_index]; i++){
		if (depth[i] >= duplicate_threshold){
			number_of_site++;
			site_removed[chr_index]++;
			num_read_removed[chr_index] += depth[i];
			fprintf(flog,"%s\t%d\t%d\n", chrs[chr_index], i, depth[i]);
		}
	}
	fclose(flog);
	
	fin = fopen(chr_reads_filenames[chr_index],"r");
	FILE *fout;
	fout = fopen(output_sam_file, "a");

	char * line_copy = NULL;
	line_copy = (char *)calloc(LEN+1,sizeof(char));

	while ((readline = fgets(line, LEN, fin)) != NULL){
		strcpy(line_copy, line);
		if(line[0] == '@')
			continue;
		strtok(line,"\t");
		strtok(NULL,"\t");
		read_chr = strtok(NULL,"\t");
		pos = atoi(strtok(NULL,"\t"));
		mapqual = atoi(strtok(NULL, "\t"));
		if (duplicate_threshold != 1){
			if (depth[pos] < duplicate_threshold){
				fprintf(fout, "%s", line_copy);
				num_output_reads++;
			} 
		}
		if (duplicate_threshold == 1){
			if (mapqual == mapq[pos]){
				fprintf(fout, "%s", line_copy);
				/* modify best mapping quality for this site,  in case two reads exist with same best mapping quality. */
				mapq[pos] = 0;
				num_output_reads++;
			}
		}
	}

	fclose(fin);
	fclose(fout);
}


void output_clean_up(void){
	printf("Number of reads in output file %s: %d\n", output_sam_file, num_output_reads);
	printf("Number of site with reads removed is: %d\n", number_of_site);
	printf("Percentage of reads removed is : %f Percent\n\n", (1.0-(1.0*num_output_reads/num_chr_reads))*100);
	
	int p;
	
	if (detail == 1){
		printf("chr\tnreads\tnreads_removed\tnum_sites\n");
		
		for (p=1; p <= count_seq; p++){
			printf("%s\t%d\t%d\t%d\n", chrs[p], num_read[p], num_read_removed[p], site_removed[p]);
		}
	}

	
	for (p=1; p <= count_seq; p++){
		if (chr_reads_filenames[p] != NULL){
			remove(chr_reads_filenames[p]);
		}
	}
}


/* Check if reference sequence are present in SAM file 
 *	If not, regular chromosomes will be used. 
 */
void prepare_groups(void){
	
	FILE *fallreads;
	char * line = NULL;

	char *readline;
	char *filename;
	char *seq;
	char *len;
	int i;
	
	fallreads = fopen(sam_file,"r");

	line = (char *)calloc(LEN+1,sizeof(char));

	while ((readline = fgets(line, LEN, fallreads)) != NULL){
		if(line[0] == '@'){
			if ((line[1] == 'S') && (line[2]=='Q')){
				count_seq++;
				filename = (char *)calloc(STR,sizeof(char));
				strtok(line,"\t");
				seq = strtok(NULL,"\t");
				chrs[count_seq] = (char *)calloc(STR,sizeof(char));
				strcpy(chrs[count_seq], seq+3);
				len = strtok(NULL, "\t\n");
				chr_length[count_seq] = atoi(len+3);
				strcpy(filename, seq+3);
				strcat(filename, "_read_");
				//strcat(filename, sam_header);
				strcat(filename, ".txt");
				chr_reads_filenames[count_seq] = (char *)calloc(STR,sizeof(char));
				strcpy(chr_reads_filenames[count_seq], filename);
				site_removed[count_seq] = 0;
			}
		} else {
			break;
		}
	}

	fclose(fallreads);
	
	if (count_seq > 0){
		printf("Reference Sequence Presented: %d Sequences. \n", count_seq);
	} else {
	/* reference sequenct not provided, consider only regular chromosome*/
		printf("Reference Sequence Not Present, Use Regular Chromosomes.\n");
		for (i=0; i<NUMREGCHR; i++){
			strcpy(filename, chrs[i]);
			strcat(filename, "_read_");
			strcat(filename, sam_header);
			strcat(filename, ".txt");
			chr_reads_filenames[i] = (char *)calloc(STR,sizeof(char));
			strcpy(chr_reads_filenames[i], filename);
			site_removed[i] = 0;
		}
		count_seq = NUMREGCHR;
	}
	
	for(i=1; i<= count_seq; i++){
		num_read[i] = 0;
		num_read_removed[i] = 0;
	}

}



void duplicate_reads_remover(void){

	/* Step 1 - Check if Sequences are Present*/
	prepare_groups();
	
	/* Step 2 - group reads into chromosome groups */
	group_reads();

	int i_chr;
	/* Step 3 - for each chromosome, calculate read depth at each site and output less redundant reads */
	for (i_chr = 1; i_chr <= count_seq; i_chr++){
		filter_one_chr(i_chr);
	}
	/* Step 4 - Data Tidy Up*/
	output_clean_up();
}

int main_removeDuplicatedReads(int argc, char *argv[]){

	if (argc < 2){
		printf("Usage: ./remove_dup sample.sam (threshold) (progress) (detail)\n ");
		printf("\tIf threshold is not present, it will be calculated as 2*exonic base coverage.\n");
		printf("\tIf threshold >1, program with remove all reads at position whose depth is greater or equal to duplicate_threshold.\n");
		printf("\tIf threshold == 1: regardless of read depth, only one read per position is retained with best MAPQ.\n\n");
		printf("\tIf option progress is present, it is in debug mode and output program progess.\n");
		printf("\tIf option detail is present, program will output duplicate removal result on chromosome basis.\n\n");
		exit(1);
	}

	FILE *ftest_existence;
	if ( (ftest_existence = fopen ( argv[1], "r" ) ) == NULL ){
		printf("SAM file not found\n");
		exit(1);
	} else {
		sam_file = (char *)calloc(STR,sizeof(char));
		sam_header = (char *)calloc(STR,sizeof(char));
		output_sam_file = (char *)calloc(STR,sizeof(char));
		strcpy(sam_file, argv[1]);
		strcpy(sam_header, argv[1]);
		strcpy(output_sam_file, sam_header);
		strcat(output_sam_file, ".NoneDupReads");
		remove_logfile = (char *)calloc(STR,sizeof(char));
		strcpy(remove_logfile, sam_header);
		strcat(remove_logfile, ".DupReads.txt");

	}
	
	if (argc >= 3){
		/*Duplicate removal threshold is present */
		if (atoi(argv[2]) >= 0){
			duplicate_threshold = atoi(argv[2]);
		} else {
			printf("Invalid duplicate read removal threshold, it must be a positive integer. \n");
			exit(1);
		}
	}
	
	if (argc >= 4){
		if (strcmp(argv[3], "progress") == 0){
			debug = 1;
		}
		if (strcmp(argv[3], "detail") == 0){
			detail = 1;
		}
	}
	
	if (argc > 5){
		if (strcmp(argv[4], "progress") == 0){
			debug = 1;
		}
		if (strcmp(argv[5], "detail") == 0){
			detail = 1;
		}
	}
	duplicate_reads_remover();

}
