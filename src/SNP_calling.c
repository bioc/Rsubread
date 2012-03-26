/*
 * calling.c
 *
 *  Created on: Jan 3, 2012
 *      Author: dai
 *  Modified by Wei Shi, March 16 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

/*	Constants for testing purpose but are variables for this function */

#define MAXCHRNUM 150
#define NUMREGCHR 24
#define STR 100	//maximum string length
#define SEG_SIZE 100000000
#define LEN 1000
#define MAX_CIGAR_COMP 100
#define MAX_INDEL 50
#define MAX_INDEL_LENGTH 100
#define DEBUG 0
#define MAX_CALLED_DEPTH 10000000

int SNP_MIN_CALLED_DEPTH = 5;
double SNP_min_fraction = 0.5; 

char *numbers[] = {"0", "1", "2","3", "4", "5", "6", "7", "8", "9", "10"};


/* regular chromosomes */
char *SNP_reg_chrs[] = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
				"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
				"chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};
				

				
/* Global Variables - funcion arguments */
char *sam_file;
char *sam_header;
char *snp_indel_filename;
char *called_snp_filename;
char *fa_file;
int read_length = 0;
int num_segs = 0;
char *fully_match_cigar;

char *chrs[MAXCHRNUM];
char *chr_reads_filenames[MAXCHRNUM];
FILE *chr_reads_files[MAXCHRNUM];
int chr_length[MAXCHRNUM];

int count_seq; //number of reference sequence

/* Global Variables - Cigar components*/
char * NUCLEOTIDES = "ATGCN";
char * CIGAR_OPERATION = "MIDNSHP=X";
int cigar_comp_length[MAX_CIGAR_COMP];
int cigar_comp_type[MAX_CIGAR_COMP];
int cigar_comp = 0;


/* Global Variables - chromosome working unit */
typedef struct workingunit{
	int chr_index;
	char *segment_id;
	int32_t offset;
	char *seq_filename;
}unit;

unit *segments[100];
int num_segs;


/* Global Variable - sequence and base coount */
short ref_seq[SEG_SIZE+1];
short count[5][SEG_SIZE+1];	//count[5] represents number of read that has a deletion at that position


/* Global Variable - Indel node */
typedef struct insertiondeletion{
	int32_t start;
	int count;
	int freq[MAX_INDEL];
	char * seq[MAX_INDEL];
	int indel_length[MAX_INDEL];
	char flag[MAX_INDEL];

}node;

node * indel[SEG_SIZE+1];


/* Global Variable - Coverage Tesing */
int32_t positions[120000];
int pos_count[120000];
int num_pos = 0;
char *exonic_coverage_filename = "Exonic_sites_coverage.txt";



/*
 *  Given the chromosome name(string), return chrosome order(int)
 *  For unidentified human chromosome, it returns -1
 * */
int SNP_find_chr_index(char *read_chr){
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



/* Process Cigar string to retrieve each component
 * */
int process_cigar_info(char *cigar){
	/* translate CIGAR string into information */
	int length = 0;
	int cigar_base;
	cigar_comp = 0;
	for (cigar_base=0; cigar_base<strlen(cigar); cigar_base++){
		if (isdigit(cigar[cigar_base])){
			length = length * 10 + cigar[cigar_base] - 48;
		} else if (isalpha(cigar[cigar_base])){
			cigar_comp_length[cigar_comp] = length;
			char *indexof = strchr(CIGAR_OPERATION, cigar[cigar_base]);
			if (indexof != NULL){
				cigar_comp_type[cigar_comp] = indexof - CIGAR_OPERATION;
			} else {
				printf("Unseen cigar operation\n");
			}
			cigar_comp++;
			length = 0;
		}
	}
	/* Calculate read length represented by this CIGAR info */
	int i;
	int sum = 0;
	for (i=0; i<cigar_comp; i++){
		if ((cigar_comp_type[i] == 0) || (cigar_comp_type[i] == 1)){
			sum += cigar_comp_length[i];
		}
	}
	return sum;
}



/*
 * For the input read mapping file(.sam file)
 * Group read according to their mapped chromosome
 * Each chromosome has its own read file
 * */
void SNP_group_reads(void){
	if (DEBUG == 1) {
		printf("In grouping read ...\n");
	}
	/* control variables declared */
	int i;

	/* Create filenames and file pointer for each chromosome_reads_file */
	char *filename;
	filename = (char *)calloc(STR,sizeof(char));

	for (i=1; i<= count_seq; i++){
		strcpy(filename, chrs[i]);
		strcat(filename, "_read_");
		strcat(filename, sam_header);
		strcat(filename, ".txt");
		chr_reads_filenames[i] = (char *)calloc(STR,sizeof(char));
		strcpy(chr_reads_filenames[i], filename);
		chr_reads_files[i] = fopen(filename, "w");
	}
	
	if (DEBUG == 1){
		printf("finished generating filenames\n");
	}

	/* Group reads into subgroups */
	FILE *fallreads;
	char * line = NULL;
	char * line_copy = NULL;
	char * read_chr;
	int read_chr_index;
	int32_t pos;
	char *cigar;
	char *read_seq;

	char *readline;

	fallreads = fopen(sam_file,"r");

	line = (char *)calloc(LEN+1,sizeof(char));
	line_copy = (char *)calloc(LEN+1,sizeof(char));

	while ((readline = fgets(line, LEN, fallreads)) != NULL){
		strcpy(line_copy, line);
		if(line[0] == '@')
			continue;
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
		//printf("%s\t %d\t %d\n", read_chr, pos, SNP_find_chr_index(read_chr));
		if( ((read_chr_index = SNP_find_chr_index(read_chr)) == -1) ){
			continue;
		}
		fprintf(chr_reads_files[read_chr_index], "%s", line_copy);
	}

	/* Close input file and read files for each chromosome */
  	fclose(fallreads);
	for (i=1; i<=count_seq; i++){
		fclose(chr_reads_files[i]);
	}
	if (DEBUG == 1){
		printf("Finished grouping read ...\n");
	}
}


/*
 * Divide each chromosome into given SEG_SIZE size segments
 * prepare the index information and divide the reference genome
 * for future SNP calling purpose
 * */
void prepare_working_units(void){
	if (DEBUG == 1){
		printf("In prepare working units ...\n");
	}
	FILE *fin;
	fin = fopen(fa_file, "r");
	FILE *fseg;

	num_segs = -1;

	char line[200];
	char *chr;
	int chr_index = 0;
	int32_t base_length;
	int this_chr_seg;
	int i;
	int line_length;
	int non_read_flag = 0;
	
	//chr = (char *)calloc(STR,sizeof(char));
	
	while ( fgets ( line, sizeof(line), fin) != NULL ){
		if (line[0] == '>'){
			/*chr = (char *)calloc(STR,sizeof(char));
			strncpy(chr, line+1, (strlen(line) - 2));
			//printf("Notice\n");
			//printf("%s", line);
			//printf("%d\n", strlen(line) - 2);
			//printf("%s\n", chr);
			//printf("Notice\n");
			chr_index = SNP_find_chr_index(chr);*/
			chr_index ++;
			chr = chrs[chr_index];
			//printf("%s\t%s", chr, line);
			if (chr_index == -1){
				// This reference sequence does not have read mapped to it
				////////////////////////////////////////////////////////
				////////////////////////////////// TBC
				non_read_flag = 1;
				printf("reference sequence %s not in SAM header. \n", chr);
			} else {
				non_read_flag = 0;
			}
			if (non_read_flag == 1){
				continue;
			}
			if (DEBUG == 1){
				printf("%d\t%s\n", chr_index, chr);
			}
			num_segs++;
			segments[num_segs] = (unit *)calloc(1,sizeof(unit));
			segments[num_segs]->chr_index = chr_index;
			segments[num_segs]->offset = 0;
			segments[num_segs]->segment_id = (char *)calloc(STR,sizeof(char));
			strcpy(segments[num_segs]->segment_id, chr);
			strcat(segments[num_segs]->segment_id, "_seg_1");
			segments[num_segs]->seq_filename = (char *)calloc(STR,sizeof(char));
			strcpy(segments[num_segs]->seq_filename, chr);
			strcat(segments[num_segs]->seq_filename, "_seg_");
			strcat(segments[num_segs]->seq_filename, numbers[1]);
			strcat(segments[num_segs]->seq_filename, "_sequence.fa");
			fseg = fopen(segments[num_segs]->seq_filename, "w");
			//printf("%s\n",segments[num_segs]->seq_filename);
			this_chr_seg = 1;
			base_length = 0;
			continue;
		}
		line_length = strlen(line) - 1;
		base_length += line_length;
		if (base_length <= SEG_SIZE){
			fprintf(fseg, "%s", line);
		} else {
			// create a new segment
			int diff_old = SEG_SIZE - base_length + line_length;
			int diff_new = base_length - SEG_SIZE;
			for (i=0; i<diff_old; i++){
				fprintf(fseg, "%c", line[i]);
			}
			fclose(fseg);
			num_segs++;
			this_chr_seg++;
			segments[num_segs] = (unit *)calloc(1,sizeof(unit));
			segments[num_segs]->chr_index = chr_index;
			segments[num_segs]->offset = SEG_SIZE * (this_chr_seg - 1);
			segments[num_segs]->segment_id = (char *)calloc(STR,sizeof(char));
			strcpy(segments[num_segs]->segment_id, chr);
			strcat(segments[num_segs]->segment_id, "_seg_");
			strcat(segments[num_segs]->segment_id, numbers[this_chr_seg]);
			segments[num_segs]->seq_filename = (char *)calloc(STR,sizeof(char));
			strcpy(segments[num_segs]->seq_filename, chr);
			strcat(segments[num_segs]->seq_filename, "_seg_");
			strcat(segments[num_segs]->seq_filename, numbers[this_chr_seg]);
			strcat(segments[num_segs]->seq_filename, "_sequence.fa");
			fseg = fopen(segments[num_segs]->seq_filename, "w");
			for (i=diff_old; i<line_length; i++){
				fprintf(fseg, "%c", line[i]);
			}
			base_length = diff_new;
		}

	} // end-while
	num_segs++;
	fclose(fseg);
	fclose(fin);


}



void read_in_ref_seq(int seg_index){
	if (DEBUG == 1){
		printf("reading in reference sequence %s\n", segments[seg_index]->segment_id);
	}

	/* readin the reference sequence into array */
	FILE *fseq = fopen(segments[seg_index]->seq_filename, "r");
	int32_t index = 0;

	int i;
	char read_char;
	while ((read_char = fgetc(fseq)) != EOF){
		if (read_char == '\n'){
			continue;
		}
		index++;
		switch( read_char )
		{
		    case 'A':
		        ref_seq[index] = 0;
		        break;
		    case 'T':
		    	ref_seq[index] = 1;
		    	break;
		    case 'G':
		    	ref_seq[index] = 2;
		    	break;
		    case 'C':
		    	ref_seq[index] = 3;
		    	break;
		    case 'N':
		    	ref_seq[index] = 4;
		    	break;
		    default:
		    	ref_seq[index] = 4;
		    	if (DEBUG == 1){
					printf("at %s position %d: we encounter %c, it is treated as 'N'. \n", segments[seg_index]->segment_id, index, read_char);
		    	}
				break;
		    	//index--;
		}
		for (i=0; i<5; i++){
			count[i][index] = 0;
		}
	}

	fclose(fseq);
	if (DEBUG == 1){
		printf("Done reading in ref seq\n");
	}

}

void *
SNP_make_empty_node(void){
	node *new;
	new = (node *) calloc(1,sizeof(node));
	new->start = 0;
	new->count = 0;
	int i;
	for (i=0; i<MAX_INDEL; i++){
		new->freq[i] = 0;
	}
	return new;
}



void process_chr_reads(int seg_index){

	char *chr = chrs[segments[seg_index]->chr_index];

	FILE * fin;
	char *read_file = chr_reads_filenames[segments[seg_index]->chr_index];
	fin = fopen(read_file, "r");

	/* Group reads into subgroups */
	char * line = NULL;
	char * read_chr;
	int32_t pos;
	char *readline;
	char *read_seq;
	char *cigar;

	int num_indel_pos = 0;

	line = (char *)calloc(LEN+1,sizeof(char));

	int i;
	int32_t offset = segments[seg_index]->offset;
//printf("break1\n");

	while ((readline = fgets(line, LEN, fin)) != NULL){
//if(seg_index == 22) printf("%s",line);

		strtok(line,"\t");
		strtok(NULL,"\t");
		read_chr = strtok(NULL,"\t");
		if( (strcmp(read_chr, chr) != 0) )
			continue;
		pos = atoi(strtok(NULL,"\t"));
		/* NOTE HERE: pos is deducted by offset value */
		pos = pos - offset;
		if ((pos < 0) || (pos >= SEG_SIZE))
			continue;
		strtok(NULL,"\t");
		cigar = strtok(NULL,"\t");
		strtok(NULL,"\t");
		strtok(NULL,"\t");
		strtok(NULL,"\t");
		read_seq = strtok(NULL, "\t");

		if (strcmp(cigar, fully_match_cigar) == 0){
			/* This read is fully matched/mismatched */
			for (i=0; i<strlen(read_seq); i++){
				switch( read_seq[i] )
				{
					case 'A':
						count[0][pos+i] += 1;
						break;
					case 'T':
						count[1][pos+i] += 1;
						break;
					case 'G':
						count[2][pos+i] += 1;
						break;
					case 'C':
						count[3][pos+i] += 1;
						break;
				}
			}
		} else {
			//printf(" %d\t%s\n", pos, cigar);
			/* This read is NOT fully matched/mismatched */
			/* Process CIGAR information */
			process_cigar_info(cigar);
			int i_comp, rpos;	// iterator for read component and each component's read position
			int32_t refer_cursor; 	// cumulated position from up to now read components
			refer_cursor = pos;		// refer cursor starts from the read mapping location
			int32_t read_cursor = 0;	// cursor on read position

			for(i_comp=0; i_comp < cigar_comp; i_comp++){
				if (refer_cursor > SEG_SIZE){
						break;
				}
				if (DEBUG == 1){
					printf("in processing %c with length %d\n", CIGAR_OPERATION[cigar_comp_type[i_comp]],cigar_comp_length[i_comp]);
				}
				switch (cigar_comp_type[i_comp])
				{
					case 0:
						/* This component of read is a 'M': Match or Mismatch*/
						/* Modify both reference genome cursor and read cursor */
						for (rpos=0; rpos<cigar_comp_length[i_comp]; rpos++){
							switch( read_seq[read_cursor] )
							{
								case 'A':
									count[0][refer_cursor] += 1;
									break;
								case 'T':
									count[1][refer_cursor] += 1;
									break;
								case 'G':
									count[2][refer_cursor] += 1;
									break;
								case 'C':
									count[3][refer_cursor] += 1;
									break;
							}
							refer_cursor++;
							read_cursor++;
						}
						break;

					case 1:
						if (1==0){	// wrapper to single statement block
						/* This component of read is a 'I': Insertion*/
						/* modify the read cursor only */
						/* retrieve inserted sequence */
						char * inserted_seq;
						inserted_seq = (char *)calloc(MAX_INDEL_LENGTH,sizeof(char));
						strncpy(inserted_seq, (char *)(read_seq+read_cursor), cigar_comp_length[i_comp]);
						//printf("DEBUG: %s,%d\n",inserted_seq,cigar_comp_length[i_comp]);
						if (DEBUG == 1){
							//printf("read sequence is : %s\n", read_seq);
							//printf("inserted sequence is : %s\n", inserted_seq);
						}
						/* check if this is the first insertion at this position */
						if (indel[refer_cursor] == NULL){
							if (DEBUG == 1){
									printf("This position %s: %d first have indel(Insertion)\n", segments[seg_index]->segment_id, pos);
							}
							num_indel_pos++;	// another position on this chromosome has insertion
							indel[refer_cursor] = (node *)SNP_make_empty_node();
							node *insert_block = indel[refer_cursor];
							insert_block->start = pos;
							insert_block->count = 1;
							insert_block->freq[0] = 1;
							insert_block->flag[0] = 'I';
							insert_block->seq[0] = (char *)calloc(MAX_INDEL_LENGTH,sizeof(char));
							strcpy(insert_block->seq[0], inserted_seq); // store the insertion sequence
						} else {
							if (DEBUG == 1){
									printf("This position %s: %d had indel before", segments[seg_index]->segment_id, pos);
							}
							node *insert_block = indel[refer_cursor];
							/* check first if this inserted sequence has appear before */
							int occurance_flag = -1;
							int sym;
							for (sym=0; sym < insert_block->count; sym++){
								if (insert_block->flag[sym] == 'I'){	// it is an insertion node
									if (strcmp(inserted_seq, insert_block->seq[sym]) == 0){
										occurance_flag = sym;
									}
								}
							}	// end for
							if (occurance_flag != -1){
								/* This inserted sequence has appear before, increments its frequency */
								insert_block->freq[occurance_flag]++;
							} else {
							//	printf("creating new insertion at this position\n");
								/* create new insertion sample */
								insert_block->count++;
								insert_block->freq[insert_block->count-1] = 1;
								insert_block->flag[insert_block->count-1] = 'I';
								insert_block->seq[insert_block->count-1] = (char *)calloc(MAX_INDEL_LENGTH,sizeof(char));
								strcpy(insert_block->seq[insert_block->count-1], inserted_seq);
								//printf("new insertion created\n");
							}	// end if  (occurance_flag != -1)
						}	//end else (indel[refer_cursor-1] == NULL)
						read_cursor = read_cursor + cigar_comp_length[i_comp];
						free(inserted_seq);
						}	// closing bracket for wrapping if statement
						read_cursor = read_cursor + cigar_comp_length[i_comp];
						break;

					case 2:
						if (1 == 0){	// wrapper to single statement block/
							/* This component of read is a 'D': Deletion*/
							/* For deletion, only need to specify from current position, how many bases are deleted
							 * the actual deleted sequence can only be known when reference genome is inputed */

							/* check if this is the first indel at this position */
							if (indel[refer_cursor] == NULL){
								if (DEBUG == 1){
									printf("This position %s: %d first have indel(Deletion)\n", segments[seg_index]->segment_id, pos);
								}
								num_indel_pos++;	// another position on this chromosome has insertion
								indel[refer_cursor] = (node *)SNP_make_empty_node();
								node *del_block = indel[refer_cursor];
								del_block->start = pos;
								del_block->count = 1;
								del_block->freq[0] = 1;
								del_block->flag[0] = 'D';
								del_block->indel_length[0] = cigar_comp_length[i_comp];
							} else {
								if (DEBUG == 1){
									printf("This position %s: %d had indel before", segments[seg_index]->segment_id, pos);
								}
								/* check first if this inserted sequence has appear before */
								node *del_block = indel[refer_cursor];
								int occurance_flag = -1;
								int sym;
								for (sym=0; sym < del_block->count; sym++){
									if (del_block->flag[sym] == 'D'){	// it is an insertion node
										if (cigar_comp_length[i_comp] == del_block->indel_length[sym]){
											occurance_flag = sym;
										}
									}
								}	// end for
								if (occurance_flag != -1){
									/* This inserted sequence has appear before, increments its frequency */
									del_block->freq[occurance_flag]++;
								} else {
									/* create new deletion sample */
									del_block->count++;
									del_block->freq[del_block->count-1] = 1;
									del_block->flag[del_block->count-1] = 'D';
									del_block->indel_length[del_block->count-1] = cigar_comp_length[i_comp];
								}	// end if  (occurance_flag != -1)
							}	//end else (indel[refer_cursor-1] == NULL)
							/* modify the refer cursor only */
							refer_cursor = refer_cursor + cigar_comp_length[i_comp];
						}	// closing bracket for wrapping if statement
						refer_cursor = refer_cursor + cigar_comp_length[i_comp];
						break;
					case 3:
						/* This component is a 'N':skip */
						/* Modify the reference genome cursor only */
						refer_cursor = refer_cursor + cigar_comp_length[i_comp];
						break;

					case 4:
						/* This component is a 'S':soft clipping */
						/* Modify the reference genome cursor only */
						refer_cursor = refer_cursor + cigar_comp_length[i_comp];
						break;
				} // end-switch
			} // end-while cigar_component
		} // end-else not fully matched


	} // end-while
//printf("break2\n");

	fclose(fin);
	if (DEBUG == 1){
		printf("Finshed processing reads\n");
	}
}

void output_called_INDELs(int seg_index){

	FILE *fallsnv = fopen(snp_indel_filename,"a");
	/* Print header line indicating each column */
	//fprintf(fallsnv,"chr\tpos\tref\talt\tfreq\tdepth\n");

	char *chr = chrs[segments[seg_index]->chr_index];
	int i, j;
	int depth = 0;
	int count_other_nuc;
	char *alt;	//	string recording possible alternatives allele
	char *freq;	// string recording frequency of each alternative allele
	int offset = segments[seg_index]->offset;
	int INDEL_freq_largest;

	for (i=1; i<=SEG_SIZE; i++){
		/* Skip this position if reference sequences says N*/
		if (ref_seq[i] == 4){
			continue;
		}

		/* check whether this position generate an output line */
		count_other_nuc = 0;
		depth = 0;
		for (j=0; j<=4; j++){
			if (ref_seq[i] != j){
				count_other_nuc += count[j][i];
			}
			depth += count[j][i];
		}

		alt = (char *)calloc(STR,sizeof(char));
		strcpy(alt,"DUMMY");
		freq = (char *)calloc(STR,sizeof(char));
		strcpy(freq,"DUMMY");

		/* check whether this position has indel */
		if ((indel[i] != NULL) && (depth >= SNP_MIN_CALLED_DEPTH) && (depth <= MAX_CALLED_DEPTH)){
			INDEL_freq_largest = 0;
			//print_indel_at_pos(i+2);
			for (j=0; j<indel[i]->count; j++){
			    	char this_freq [ STR ];
			    	sprintf(this_freq, "%d", indel[i]->freq[j]);
				strcat(alt, "/");
				strcat(freq, "/");
				if  (indel[i]->flag[j] == 'I'){
					char *this_allele = indel[i]->seq[j];
					if(strlen(this_allele)+strlen(alt) > 50) 
						break;
					strcat(alt,"+");
					strcat(alt, this_allele);
                                        if(indel[i]->freq[j] > INDEL_freq_largest)
                                                INDEL_freq_largest = indel[i]->freq[j];
				} else if (indel[i]->flag[j] == 'D'){
					char this_allele[20];
					int p;
					for (p=0; p<indel[i]->indel_length[j]; p++){
						sprintf(this_allele, "%c", NUCLEOTIDES[ref_seq[i+p]]);
					}
					if(strlen(this_allele)+strlen(alt) > 50)
						break;
					strcat(alt,"-");
					strcat(alt, this_allele);
                                        if(indel[i]->freq[j] > INDEL_freq_largest)
                                                INDEL_freq_largest = indel[i]->freq[j];
				}
				if(strlen(this_freq)+strlen(freq) > 50)
					break;
				strcat(freq, this_freq);
			}
			if(INDEL_freq_largest >= SNP_min_fraction*depth)
				fprintf(fallsnv, "%s\t%d\t%c\t%s\t%s\t%d\n", chr, i+offset, NUCLEOTIDES[ref_seq[i]], alt+6, freq+6, depth);
		}

		if(alt) 
			free(alt);

		if(freq) 
			free(freq);

		/* In the final output line, both alt and freq should be refered with pointer +1
		 * in this case, regardless of how many possible alternatives are there, each is started with
		 * a slash, and the final output neglect the first slash by +1 */
	}
	fclose(fallsnv);
}



void output_called_SNPs(int seg_index){


	FILE *fsnp = fopen(called_snp_filename,"a");
	//fprintf(fsnp,"chr\tpos\tref\talt\tfreq\tdepth\n");

	char *chr = chrs[segments[seg_index]->chr_index];
	int i, j, isSNP;
	int depth = 0;
	int count_other_nuc;
	int offset = segments[seg_index]->offset;

	for (i=1; i<=SEG_SIZE; i++){
		/* Skip this position if reference sequences says N*/
		if (ref_seq[i] == 4){
			continue;
		}

		/* check whether this position generate an output line */
		count_other_nuc = 0;
		depth = 0;
		for (j=0; j<=4; j++){
			if (ref_seq[i] != j){
				count_other_nuc += count[j][i];
			}
			depth += count[j][i];
		}

		if ((count_other_nuc < SNP_min_fraction*depth) || (depth < SNP_MIN_CALLED_DEPTH) || (depth > MAX_CALLED_DEPTH)){
			continue;
		}

		
		int min_alter_count = MAX_CALLED_DEPTH;
		
		for (j=0; j<=3; j++){
			if ((ref_seq[i] != j) && (count[j][i] < min_alter_count)){
				min_alter_count = count[j][i];
			}
		}
		
		char * alt = (char *)calloc(STR,sizeof(char));
		strcpy(alt,"DUMMY");
		char * freq = (char *)calloc(STR,sizeof(char));
		strcpy(freq,"DUMMY");

		/* check whether position i has SNP */
		isSNP = 0;
		for (j=0; j<=3; j++){
			if ((ref_seq[i] != j) && (count[j][i] > min_alter_count)){
				isSNP = 1;
				/* Add in this allele's information to alt and freq */
				/* seperate each possible allele by a slash */
				char this_allele[3];
				sprintf(this_allele, "%c", NUCLEOTIDES[j]);
			    	char this_freq [ STR ];
				sprintf(this_freq, "%d", count[j][i]);
				strcat(alt, "/");
				strcat(freq, "/");
				strcat(alt, this_allele);
				strcat(freq, this_freq);
			}
		}

		if(isSNP == 1)
			fprintf(fsnp, "%s\t%d\t%c\t%s\t%s\t%d\n", chr, i+offset, NUCLEOTIDES[ref_seq[i]], alt+6, freq+6, depth);

		if(alt) free(alt);
		if(freq) free(freq);

		/* In the final output line, both alt and freq should be refered with pointer +1
		 * in this case, regardless of how many possible alternatives are there, each is started with
		 * a slash, and the final output neglect the first slash by +1 */
	}
	
	fclose(fsnp);


}


void clean_up_files(void){
	if (DEBUG == 1){
			printf("In cleaning up files...\n");
	}
	
	int p;
	for (p=0; p<MAXCHRNUM; p++){
		if (chr_reads_filenames[p] != NULL){
			remove(chr_reads_filenames[p]);
		}
	}
	for (p=0; p<num_segs; p++){
		if  ((segments[p] != NULL) && (segments[p]->seq_filename)){
			remove(segments[p]->seq_filename);
		}
	}
	if (DEBUG == 1){
			printf("Finished cleaning up files...\n");
	}
}


/*void coverage_count(char *chr){

	int chr_index;
	chr_index = SNP_find_chr_index(chr);


	FILE *freads;
	freads = fopen(chr_reads_filenames[chr_index], "r");

	char *readline;
	char *line;
	char *read_chr;
	int32_t pos;
	line = (char *)calloc(LEN+1,sizeof(char));

	int32_t start, end;
	int32_t p,q, mid;

	while (((readline = fgets(line, LEN, freads)) != NULL)){
		strtok(line,"\t");
		strtok(NULL,"\t");
		read_chr = strtok(NULL,"\t");
		pos = atoi(strtok(NULL,"\t"));

		/* Search the start and end covering sites */
		/*start = 0;
		end = num_pos;

		if ((pos < positions[0]) || (pos > positions[num_pos-1])){
			continue;
		}
		p = 0;
		q = num_pos-1;
		while (p < q){
			mid = (p+q)/2;
			if (positions[mid] < pos){
				p = mid+1;
			} else {
				q = mid-1;
			}
		}
		start = p;
		//printf("Located at position:%d,  %d\t%d\n", p, positions[start], num_pos);

		/*while (positions[start] < pos){
			start++;
		}*/

		/* Increment read depth */
		//if (p < 0){
		//	continue;
		//}
	/*	while ((p<num_pos) && (positions[p] <= (pos + read_length - 1))){
			if (positions[p] >= pos){
				pos_count[p]++;
				//printf("read position: %d is incrementing exonic position: %d\n", pos, positions[p]);
			}
			p++;
		}
	}

	
	FILE *fout;
	fout = fopen(exonic_coverage_filename, "a");
	for (p=0; p<num_pos; p++){
		fprintf(fout, "%s\t%d\t%d\n", read_chr, positions[p], pos_count[p]);
	}

	fclose(fout);
	fclose(freads);
}*/

/* Check if reference sequence are present in SAM file 
 *	If not, regular chromosomes will be used. 
 */
void SNP_prepare_groups(void){
	if (DEBUG == 1){
		printf("In preparing groups...\n");
	}
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
				strtok(line,"\t");
				seq = strtok(NULL,"\t");
				chrs[count_seq] = (char *)calloc(STR,sizeof(char));
				strcpy(chrs[count_seq], seq+3);
				len = strtok(NULL, "\t\n");
				chr_length[count_seq] = atoi(len+3);
				if (DEBUG == 1){
					printf("Reference Sequence: %s\n", chrs[count_seq]);
				}
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
			chrs[i+1] = (char *)calloc(STR,sizeof(char));
			strcpy(chrs[i+1], SNP_reg_chrs[i]);
		}
		count_seq = NUMREGCHR;
	}
	
}



/* This is the main working function that implements SNP calling.
 * It works with a segment of a chromosome at one time and output
 * SNP statistics to file.
 * */
void SNP_calling(void){

	int i_main;
	
	 // Step 1:  reads into chromosome groups, each chromosome has one reads file
	SNP_prepare_groups();
	SNP_group_reads();

	 // Step 2: Decide and divide working segments, it's part of a chromosome
	prepare_working_units();
	
	 // Step 3: For each working segment, calculate SNP statistics
	for (i_main = 0; i_main < num_segs; i_main++){
	//printf("DEBUG: Processing segment: %d, %s\n", i_main, segments[i_main]->segment_id);

	//printf("DEBUG: start of 3.1\n");

		/* Step 3.1 Read in reference genome sequence */
		read_in_ref_seq(i_main);
//printf("DEBUG: start of 3.2\n");

		/* Step 3.2 Process reads */
		process_chr_reads(i_main);
//printf("DEBUG: start of 3.3\n");

		/* Step 3.3 Output called INDELs */
		//output_called_INDELs(i_main);

//printf("DEBUG: start of 3.4\n");

		/* Step 3.4 Output called SNPs */
		output_called_SNPs(i_main);
//printf("DEBUG: end of 3.4\n");

	}

//printf("DEBUG: end of i_main for loop\n");

	clean_up_files();

}



/* Main Function checks whether program arguments are complete and legal
 * It calls SNP_calling function to do calculation
 * */
int main_SNPcalling(int argc, char *argv[]){

	if (argc == 1){
		printf("Usage: ./calling sample.sam read_length sample_reference.fa output_file minReadCoverage minAlleleFraction\n");
		exit(1);
	}

	SNP_MIN_CALLED_DEPTH = atoi(argv[5]); 
	SNP_min_fraction = atof(argv[6]);

	printf("Minimal read coverage required is: %d\n",SNP_MIN_CALLED_DEPTH);
	printf("Minimal fraction of supporting reads required is %f\n\n",SNP_min_fraction);

	FILE *ftest_existence;
	if ( (ftest_existence = fopen ( argv[1], "r" ) ) == NULL ){
		printf("SAM file not found\n");
		exit(1);
	} else {
		sam_file = (char *)calloc(STR,sizeof(char));
		sam_header = (char *)calloc(STR,sizeof(char));
		snp_indel_filename = (char *)calloc(STR,sizeof(char));
		called_snp_filename = (char *)calloc(STR,sizeof(char));
		strcpy(sam_file, argv[1]);
		strcpy(sam_header, argv[4]);
		strcpy(snp_indel_filename, sam_header);
		strcpy(called_snp_filename, sam_header);
		strcat(snp_indel_filename, ".INDEL");
		//strcat(called_snp_filename, ".SNP");
		//if ( (ftest_existence = fopen ( snp_indel_filename, "r" ) ) != NULL ){
		//	remove(snp_indel_filename);
		//}
		//if ( (ftest_existence = fopen ( called_snp_filename, "r" ) ) != NULL ){
		//	remove(called_snp_filename);
		//}
		FILE * ftmp = fopen(called_snp_filename,"w");
		fprintf(ftmp,"chr\tpos\tref\talt\tfreq\tdepth\n");		
		fclose(ftmp);
	}

	if (atoi(argv[2]) == 0){
		printf("Invalid read length\n");
		exit(1);
	} else {
		read_length = atoi(argv[2]);
		fully_match_cigar = (char *)calloc(STR,sizeof(char));
		strcpy(fully_match_cigar, argv[2]);
		strcat(fully_match_cigar, "M");
	}

	if ( (ftest_existence = fopen ( argv[3], "r" ) ) == NULL ){
		printf("Reference genome file not found\n");
		exit(1);
	} else {
		fa_file = (char *)calloc(STR,sizeof(char));
		strcpy(fa_file, argv[3]);
	}

	clock_t start = clock();
	SNP_calling();
	printf("Time elapsed for is: %f seconds\n",  ((double)clock() - start) / CLOCKS_PER_SEC);

}
