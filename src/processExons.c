#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <R.h>
#define MAX_GENE 1000000
#define MAX_CHR 200
#define STR 100

typedef struct an_exon{
	int32_t start, end;
	char orientation;
	struct an_exon *next;
}exon;

typedef struct a_chr{
	char *id;
	exon *exon_list;
}chr;

typedef struct a_gene{
	int32_t id;
	int chr_index;
	int chr_num;
	char *current_chr;
	chr chr_array[MAX_CHR];
}gene;


// global variables
int exon_num;
int gene_index;
int gene_num;
gene gene_array[MAX_GENE];
int32_t current_gene_id;



void
initialise_process_exons(){
	gene_num = 0;
	gene_index = 0;
	current_gene_id = 0;
	exon_num = 0;
}

// returns a pointer to a dummy head node to start off the linked list
// usage: 	exon *dh = (exon *)make_empty();
void *
make_empty(void)
{
	exon *head;
	head = (exon *) malloc(sizeof(*head));
	head->next = NULL;	// next = NULL implies end of list
	return head;
}

void print_list(FILE *fout, int x, int y){

	char *chromosome = gene_array[x].chr_array[y].id;
	int32_t gene = gene_array[x].id;
	exon *list = gene_array[x].chr_array[y].exon_list;;
	list = list->next;
	while (list !=NULL){
		fprintf(fout, "%d\t%s\t%d\t%d\n",gene, chromosome, list->start, list->end);
		exon_num++;
		list = list->next;
	}
}


void
insert_new_exon(void *pt, int32_t start, int32_t end){
	exon *node = (exon *)pt;
	exon *new_node;
	new_node = (exon *)(malloc(sizeof(exon)));
	new_node->start = start;
	new_node->end = end;
	new_node->next = node->next;
	node->next = new_node;
}


void
insert_exon(void *pt, int32_t start, int32_t end){
	exon *node = (exon *)pt;

	while (node ->next != NULL){
		//node = node->next;

		if ((start >= node-> next->start) && (start <= node->next->end)){
			if (end > node->next->end){
				node->next->end = end;
				return;
			} else {
				return;
			}
		}
		if ((end >= node->next->start) && (end <= node->next->end)){
			if (start < node->next->start){
				node->next->start = start;
				return;
			} else {
				return;
			}
		}
		if (start == ((node->next->end)+1)){
			node->next->end = end;
			return;
		}
		if (end == ((node->next->start)-1)){
			node->next->start = start;
			return;
		}
		if((node->end < start) && (node->next->start > end)){
			insert_new_exon(node, start, end);
			return;
		}
		node = node->next;
	}

	if (start == ((node->end)+1)){
		node->end = end;
		return;
	}
	if (end == ((node->start)-1)){
		node->start = start;
		return;
	}
	insert_new_exon(node,start,end);
}


void output_exons_to_file(){
	FILE *fout;
	fout = fopen("exon.txt","w");
	fprintf(fout, "entrezid	chromosome	chr_start	chr_stop\n");
	int j,i;
	int total_chr_num;
	for(j=0; j<gene_num; j++){
		total_chr_num = gene_array[j].chr_num;
		for(i=0; i<total_chr_num; i++){
			print_list(fout, j, i);
		}
		/*// inform uses if this gene appears on different chromosomes
		if (total_chr_num>1){
			Rprintf("GeneID: %d appears on chromosomes: ", gene_array[j].id);
			for(i=0; i<total_chr_num; i++){
				Rprintf("%s ",gene_array[j].chr_array[i].id);
			}
			Rprintf("\n");
		}*/
	}
	fclose(fout);
}



int
find_gene(int32_t target_gene_id){
	// already know that target_gene_id != current_gene_id
	gene_index = gene_num-1;

	while ((gene_index>=0) && (gene_array[gene_index].id != target_gene_id)){
		gene_index--;
	}
	// new gene_id, build a new gene record
	if ((gene_index < 0) && (gene_num < MAX_GENE)) {
		gene_num++;
		gene_index = gene_num - 1;
		gene_array[gene_index].id = target_gene_id;
		gene_array[gene_index].chr_index = 0;
		gene_array[gene_index].chr_num = 0;
		gene_array[gene_index].current_chr = "";
	}
	current_gene_id = target_gene_id;
	return gene_index;
}


int
find_chr(int gene_pos, char *chr_id){
	int current_chr_index;
	if (strcmp(chr_id, gene_array[gene_pos].current_chr) != 0){
		// current chromosome id doesn't match what we are processing.
		current_chr_index = gene_array[gene_pos].chr_num - 1;
		while ((current_chr_index >=0) && (strcmp(chr_id, gene_array[gene_pos].chr_array[current_chr_index].id) != 0)){
			current_chr_index--;
		}
		if ((current_chr_index < 0) && (gene_array[gene_pos].chr_num < MAX_CHR)){
			//build a new chromosome
			gene_array[gene_pos].chr_num++;
			gene_array[gene_pos].chr_index = gene_array[gene_pos].chr_num - 1;
			current_chr_index = gene_array[gene_pos].chr_index;
			gene_array[gene_pos].chr_array[current_chr_index].id = (char *)malloc(STR);
			char *pt = strcpy(gene_array[gene_pos].chr_array[current_chr_index].id,chr_id);
			gene_array[gene_pos].current_chr = gene_array[gene_pos].chr_array[current_chr_index].id;
			gene_array[gene_pos].chr_array[current_chr_index].exon_list =(exon *)make_empty();

		}
	} else {
		current_chr_index = gene_array[gene_pos].chr_index;
	}

	return(current_chr_index);
}




void *
find_list(int32_t gene_id, char *chr_id){
	int p_gene=0;
	int p_chr=0;
	// locate gene array index
	if (current_gene_id == gene_id){
		p_gene = gene_index;
	} else {
		p_gene = find_gene(gene_id);
	}
	if (p_gene < 0){
		printf("exceeding max number of genes that can be processed.\n");
		return(NULL);
	}
	// locate chromosome array index
	p_chr = find_chr(p_gene, chr_id);
	if (p_chr < 0) {
		printf("exceed the maximum number of chromosomes for parallel gene with id:%d\n",gene_id);
		return(NULL);
	}
	// now, we locate the linked list head.
	exon *dh = gene_array[p_gene].chr_array[p_chr].exon_list;
	return dh;
}

void
processExons(char **input){
	FILE *fin;
	int32_t start, end, gene_id;
	char chr_id[20];
	char ori;
	exon *dh;

	initialise_process_exons();
	fin = fopen(*input,"r");
	while  (fscanf(fin, "%s %d %d %c %d ",chr_id, &start, &end, &ori, &gene_id) != -1){
		dh = (exon *)find_list(gene_id, chr_id);
		if (dh != NULL){
			insert_exon(dh, start, end);
		}
	}
	fclose(fin);
	output_exons_to_file();
}
