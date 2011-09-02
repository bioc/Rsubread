#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <R.h>

/*	Constants */
#define MAXCHRNUM 24
#define STR 200


/*	Data Structure */

typedef struct an_node{
	int32_t start, end, gene;
	struct an_node *next;
}node;


typedef struct a_chr{
	char *id;
	node *list;
}chr;

/* Variables */
chr gene[MAXCHRNUM];

int chr_num;
char current_chr_id[STR];

void *
make_empty_node_generate_ir(void){
	node *new;
	new = (node *) malloc(sizeof(node));
	new->next = NULL;
	new->start = 0;
	new->end = 0;
	new->gene = 0;
	return new;
}


/* print out the region between meta-genes : integenic region*/
void
output_ir_file(char *ir_file){
	FILE *fir;
	fir = fopen(ir_file,"w");
	//fprintf(fir, "chromosome\tchr_start\tchr_stop\n");
	int i;
	node *dh;
	for(i=0; i<chr_num; i++){
		dh = gene[i].list->next;
		while (dh->next != NULL){
			fprintf(fir, "%s\t%d\t%d\n",gene[i].id, (dh->end+1), (dh->next->start-1));
			dh = dh->next;
		}
	}
	fclose(fir);
}


/*cheating function, build meta-gene data structure from gene file*/
void
generate_ir_file(char **gene_file, char **ir_file){
	FILE *fin;
	int32_t read_start, read_end, read_entrezid;
	char chr_id[STR];
	node *cnode=NULL;
	char current_chr_id[STR];

	chr_num = 0;
	fin = fopen(*gene_file,"r");
	while  (fscanf(fin, "%d %s %d %d",&read_entrezid, chr_id, &read_start, &read_end) != -1){
		if (strcmp(chr_id, current_chr_id) != 0){
			strcpy(current_chr_id, chr_id);
			chr_num++;
			gene[chr_num-1].id = (char *)malloc(STR);
			strcpy(gene[chr_num-1].id, chr_id);
			gene[chr_num-1].list = (node *)make_empty_node_generate_ir();
			cnode = gene[chr_num-1].list;
		}

		node *new_node;
		if (cnode->end >= read_start){
			if (read_end > cnode->end){
				cnode->end = read_end;
			}
		} else {
			new_node = (node *)make_empty_node_generate_ir();
			new_node->start = read_start;
			new_node->end = read_end;
			new_node->gene = read_entrezid;
			cnode->next = new_node;
			cnode = cnode->next;
		}
	}
	fclose(fin);
	output_ir_file(*ir_file);
}



