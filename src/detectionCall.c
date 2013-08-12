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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <R.h>
//typedef int int;

/*	INPUT: 	SAM FILE
 * 			GENE(EXON) ANNOTATION FILE BINSIZE=2000 WITH GC CONTENT
 * 			BG ANNOTATION FILE BINSIZE=2000 WITH GC CONTENT
 * 	OUTPUT: MAPPING RESULT FOR TWO ANNOTATION FILES
 *
 * 	PURPOSE: THE OUTPUT OF THIS FUNCTION IS THEN USED BY R FUNCTION
 * 			 TO CALCULATE THE P-VALUE FOR EACH GENE
 * */


/*	1. build up data structure for exon and integenic region
 *  2. simplify SAM, save just chr_id and read_position
 *  3. read_mapping procedure, map this SAM file reads to exon/ir
 * */


/* Constants*/
#define BINSIZE 2000
#define STR 300
#define MAXCHRNUM 24
#define MAXREADSPERCHR 2000000

/* chromosome names are in alphabic order that */
char *chrs_map[] = {"chr1", "chr10", "chr11","chr12", "chr13",
		"chr14","chr15","chr16", "chr17","chr18",
		"chr19", "chr2","chr20","chr21","chr22","chr3",
		"chr4", "chr5","chr6", "chr7", "chr8","chr9", "chrX","chrY"};
		
		
/* Variables*/
char *SAM_file;
char *simplified_SAM_file;
char *sorted_simplified_SAM_file;

char *annotation_exon_file;
char *annotation_ir_file;
char *mapping_result_exon;
char *mapping_result_ir;

/* Time Variables*/
time_t start_time;
time_t timer;
double timediff;
/*	Data Structure */

typedef struct an_node{
	int start, end, gene;
	int read, nnum, gcnum, atnum;
	struct an_node *next;
}node;


typedef struct a_chr{
	char *id;
	node *list;
}chr;

chr exon[MAXCHRNUM];
chr ir[MAXCHRNUM];

int chr_num;
char current_chr_id[STR];

void
generate_filenames(char *dataset, char *exon_file, char *ir_file, char *temp_header){

	SAM_file = (char *)malloc(STR);
	simplified_SAM_file = (char *)malloc(STR);
	sorted_simplified_SAM_file = (char *)malloc(STR);
	annotation_exon_file = (char *)malloc(STR);
	annotation_ir_file = (char *)malloc(STR);
	mapping_result_exon = (char *)malloc(STR);
	mapping_result_ir = (char *)malloc(STR);

	strcpy(SAM_file, dataset);

	strcpy(simplified_SAM_file, temp_header);
	strcat(simplified_SAM_file,"_simplified.txt");

	strcpy(sorted_simplified_SAM_file, temp_header);
	strcat(sorted_simplified_SAM_file,"_simplified_sorted.txt");

	strcpy(annotation_exon_file, exon_file);
	
	strcpy(annotation_ir_file, ir_file);

	strcpy(mapping_result_exon, temp_header);
	strcat(mapping_result_exon,"_mapping_exon_GC.txt");
	
	strcpy(mapping_result_ir, temp_header);
	strcat(mapping_result_ir,"_mapping_binned_integenic_region_GC.txt");

}

void *
make_empty_node_map(void){
	node *new;
	new = (node *) malloc(sizeof(node));
	new->next = NULL;
	new->start = 0;
	new->end = 0;
	new->gene = 0;
	new->read = 0;
	new->nnum = 0;
	new->gcnum = 0;
	new->atnum = 0;
	return new;
}


void
build_exon_data_structure_map(void){
	FILE *fin;
	int read_start, read_end, read_entrezid, read_n, read_gc, read_at;
	char chr_id[STR];
	node *cnode=NULL;
	char current_chr_id[STR];

	chr_num = 0;
	fin = fopen(annotation_exon_file,"r");
	while  (fscanf(fin, "%d %s %d %d %d %d %d",&read_entrezid, chr_id, &read_start, &read_end, &read_n, &read_gc, &read_at) != -1){
		if (strcmp(chr_id, current_chr_id) != 0){
			strcpy(current_chr_id, chr_id);
			chr_num++;
			exon[chr_num-1].id = (char *)malloc(STR);
			strcpy(exon[chr_num-1].id, chr_id);
			exon[chr_num-1].list = (node *)make_empty_node_map();
			cnode = exon[chr_num-1].list;
		}

		node *new_node;
		new_node = (node *)make_empty_node_map();
		new_node->start = read_start;
		new_node->end = read_end;
		new_node->gene = read_entrezid;
		new_node->nnum = read_n;
		new_node->gcnum = read_gc;
		new_node->atnum = read_at;
		cnode->next = new_node;
		cnode = cnode->next;
	}
	fclose(fin);
}


void
build_ir_data_structure_map(void){
	FILE *fin;
	int read_start, read_end, read_entrezid, read_n, read_gc, read_at;
	char chr_id[STR];
	node *cnode=NULL;
	char current_chr_id[STR];

	chr_num = 0;
	fin = fopen(annotation_ir_file,"r");
	while  (fscanf(fin, "%s %d %d %d %d %d",chr_id, &read_start, &read_end, &read_n, &read_gc, &read_at) != -1){
		if (strcmp(chr_id, current_chr_id) != 0){
			strcpy(current_chr_id, chr_id);
			chr_num++;
			ir[chr_num-1].id = (char *)malloc(STR);
			strcpy(ir[chr_num-1].id, chr_id);
			ir[chr_num-1].list = (node *)make_empty_node_map();
			cnode = ir[chr_num-1].list;
		}

		node *new_node;
		new_node = (node *)make_empty_node_map();
		new_node->start = read_start;
		new_node->end = read_end;
		new_node->nnum = read_n;
		new_node->gcnum = read_gc;
		new_node->atnum = read_at;
		cnode->next = new_node;
		cnode = cnode->next;
	}
	fclose(fin);
}


void
simplify_SAM_file(void){
	FILE *fin;
	FILE *fout;
	char * line = NULL;
	char * read_chr;
	size_t len = 1000;
	size_t z;
	int read_pos;
	char *readline;

	fin = fopen(SAM_file,"r");
	fout = fopen(simplified_SAM_file,"w");

	line = (char *)malloc(len+1);

	while ((readline = fgets(line, len, fin)) != NULL){
		if(line[0] == '@')
   		 	continue;
 		strtok(line,"\t");
  		strtok(NULL,"\t");
  		read_chr = strtok(NULL,"\t");
  		if(*read_chr == '*')
  			continue;
  		read_pos = atoi(strtok(NULL,"\t"));
  		fprintf(fout, "%s\t%d\n", read_chr, read_pos);
  	}
  	fclose(fin);
  	fclose(fout);
}


void
read_mapping(void){
	/* local variables */
	FILE *fread;
	char read_chr[STR];
	int read_pos, count;
	int chr_pos;
	node *cnode;


	/* read_mapping_exon*/
	chr_pos = 0;
	// count is the counter, every 2,000,000 for a chromosome,
	// we reset the current node to the head of list
	count = 0;

	cnode = exon[chr_pos].list->next;

	fread = fopen(sorted_simplified_SAM_file,"r");

	while ((fscanf(fread, "%s %d", read_chr, &read_pos)) != -1){
		// firstly find the chromosome
		if (strcmp(read_chr, exon[chr_pos].id) != 0){
			chr_pos++;
			cnode = exon[chr_pos].list->next;
			count = 1;
		} else {
			count++;
		}

		if (chr_pos == chr_num){
			// This read has an unknown chromosome id
			// do nothing
		} else {
			// go through the list, find the exon position.
			// increment read by 1
			if ((read_pos >= cnode->start) && (read_pos <= cnode->end)){
				cnode->read++;
			} else {
				// search the rest of list
				while ((cnode->next != NULL) && (cnode->next->start < read_pos)){
					cnode = cnode->next;
				}
				if ((cnode->start <= read_pos) && (cnode->end >= read_pos)){
					cnode->read++;
				}
			}
		}
		if (count == MAXREADSPERCHR){
			// reset the current node
			cnode = exon[chr_pos].list->next;
		}
	}
	fclose(fread);


	/* read_mapping_integenic_region*/
	chr_pos = 0;
	count = 0;
	cnode = ir[chr_pos].list->next;

	fread = fopen(sorted_simplified_SAM_file,"r");

	while ((fscanf(fread, "%s %d", read_chr, &read_pos)) != -1){
		// firstly find the chromosome
		if (strcmp(read_chr, ir[chr_pos].id) != 0){
			chr_pos++;
			cnode = ir[chr_pos].list->next;
			count = 1;
		} else {
			count++;
		}

		if (chr_pos == chr_num){
			// This read has an unknown chromosome id
			// do nothing
		} else {
			// go through the list, find the exon position.
			// increment read by 1
			if ((read_pos >= cnode->start) && (read_pos <= cnode->end)){
				cnode->read++;
			} else {
				// search the rest of list
				while ((cnode->next != NULL) && (cnode->next->start < read_pos)){
					cnode = cnode->next;
				}
				if ((cnode->start <= read_pos) && (cnode->end >= read_pos)){
					cnode->read++;
				}
			}
		}
		if (count == MAXREADSPERCHR){
			// reset the current node
			cnode = exon[chr_pos].list->next;
		}
	}
	fclose(fread);

}


void
output_mapping_result(void){

	/* Local Variables */
	int i;
	node *dh;

	/* Output exon result to file*/

	FILE *fout_exon = fopen(mapping_result_exon, "w");

	fprintf(fout_exon,"entrezid\tchr\tchr_start\tchr_stop\tnnum\tgcnum\tatnum\tnreads\n");
	for(i=0; i<chr_num; i++){
		dh = exon[i].list;
		while (dh->next != NULL){
			dh = dh->next;
			fprintf(fout_exon, "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", dh->gene, exon[i].id, dh->start, dh->end, dh->nnum, dh->gcnum, dh->atnum, dh->read);
		}
	}
	fclose(fout_exon);


	/* Output integenic result to file*/

	FILE *fout_ir = fopen(mapping_result_ir, "w");
	fprintf(fout_ir,"chr\tchr_start\tchr_stop\tnnum\tgcnum\tatnum\tnreads\n");
	for(i=0; i<chr_num; i++){
		dh = ir[i].list;
		while (dh->next != NULL){
			dh = dh->next;
			fprintf(fout_ir, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", ir[i].id, dh->start, dh->end, dh->nnum, dh->gcnum, dh->atnum, dh->read);
		}
	}
	fclose(fout_ir);

}

void q_sort(int numbers[], int left, int right)
{
  int pivot, l_hold, r_hold;

  l_hold = left;
  r_hold = right;
  pivot = numbers[left];
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      right--;
    }
  }
  numbers[left] = pivot;
  pivot = left;
  left = l_hold;
  right = r_hold;
  if (left < pivot)
    q_sort(numbers, left, pivot-1);
  if (right > pivot)
    q_sort(numbers, pivot+1, right);
}




/* From the simplified SAM file generated above
 * for each chromosome id, find all reads on that chromosome
 * sort the read positions and then output to a sorted_SAM_simplified file*/
void
sort_reads(void){
	FILE *freadin;
	FILE *freadout;
	int pos[MAXREADSPERCHR];
	int i;
	char read_chr[STR];
	int read_pos;
	int count;
	int j;
	freadout = fopen(sorted_simplified_SAM_file,"w");
	for (i=0; i<24; i++){
		count = 0;
		// extracting reads with chromosome id == chrs_map[i]
		freadin = fopen(simplified_SAM_file,"r");
		while ((fscanf(freadin, "%s %d", read_chr, &read_pos)) != -1){
			if (strcmp(read_chr, chrs_map[i]) == 0){
				// this read is on the target chromosome
				pos[count] = read_pos;
				count++;
				if (count == MAXREADSPERCHR){
					// sort these reads from the particular chromosome
					// then write them to sorted file
					q_sort(pos, 0, count-1);
					// export sorted list to file.
					for (j=0; j<count; j++){
						fprintf(freadout, "%s %d\n", chrs_map[i], pos[j]);
					}
					count = 0;
				}
			}
		}
		q_sort(pos, 0, count-1);
		// export sorted list to file.
		for (j=0; j<count; j++){
			fprintf(freadout, "%s %d\n", chrs_map[i], pos[j]);
		}
		fclose(freadin);
	}
	fclose(freadout);
}

void
detectionCall(char **dataset, char **exon_file, char **ir_file, char **temp_header){

	generate_filenames(*dataset, *exon_file, *ir_file, *temp_header);
	
	build_exon_data_structure_map();
	build_ir_data_structure_map();
	
	simplify_SAM_file();

	sort_reads();
	
	read_mapping();
	
	output_mapping_result();
	
	remove(simplified_SAM_file);
	remove(sorted_simplified_SAM_file);
}

