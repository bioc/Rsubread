#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <R.h>

/*	Constants */

#define MAXCHRNUM 24
#define STR 100	//maximum string length
#define NPERCENTAGE 0.3

/* chromosome names are in alphabic order that */
char *chrs[] = {"chr1", "chr10", "chr11","chr12", "chr13",
		"chr14","chr15","chr16", "chr17","chr18",
		"chr19", "chr2","chr20","chr21","chr22","chr3",
		"chr4", "chr5","chr6", "chr7", "chr8","chr9", "chrX","chrY"};


/*	Data Structure */

typedef struct an_node{
	int32_t start, end, gene;
	int read, nnum, gcnum, atnum;
	struct an_node *next;
}node;


typedef struct a_chr{
	char *id;
	node *list;
}chr;

/* Variables */
int binsize;
chr exon[MAXCHRNUM];
chr ir[MAXCHRNUM];

int chr_num;
char current_chr_id[STR];

char *utr_cds_file;
char *exon_file;
char *ir_file;
char *binned_ir_file;


void *
make_empty_node(void){
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

/*cheating function, build exon data structure from exon file*/
void
build_exon_data_structure(void){
	FILE *fin;
	int32_t read_start, read_end, read_entrezid;
	char chr_id[STR];
	node *cnode=NULL;
	char current_chr_id[STR];

	chr_num = 0;
	fin = fopen(exon_file,"r");
	while  (fscanf(fin, "%d %s %d %d",&read_entrezid, chr_id, &read_start, &read_end) != -1){
		if (strcmp(chr_id, current_chr_id) != 0){
			strcpy(current_chr_id, chr_id);
			chr_num++;
			exon[chr_num-1].id = (char *)malloc(STR);
			strcpy(exon[chr_num-1].id, chr_id);
			exon[chr_num-1].list = (node *)make_empty_node();
			cnode = exon[chr_num-1].list;
		}

		node *new_node;
		new_node = (node *)make_empty_node();
		new_node->start = read_start;
		new_node->end = read_end;
		new_node->gene = read_entrezid;
		cnode->next = new_node;
		cnode = cnode->next;
	}
	fclose(fin);
}


/* Build from binned integenic region */
void
build_ir_data_structure(void){
	FILE *fin;
	int32_t read_start, read_end, read_entrezid, read_length;
	char chr_id[STR];
	node *cnode=NULL;
	char current_chr_id[STR];

	chr_num = 0;
	fin = fopen(binned_ir_file,"r");
	while  (fscanf(fin, "%d %s %d %d",&read_entrezid, chr_id, &read_start, &read_end) != -1){
		if (strcmp(chr_id, current_chr_id) != 0){
			strcpy(current_chr_id, chr_id);
			chr_num++;
			ir[chr_num-1].id = (char *)malloc(STR);
			strcpy(ir[chr_num-1].id, chr_id);
			ir[chr_num-1].list = (node *)make_empty_node();
			cnode = ir[chr_num-1].list;
		}

		node *new_node;
		new_node = (node *)make_empty_node();
		new_node->start = read_start;
		new_node->end = read_end;
		new_node->gene = read_entrezid;
		cnode->next = new_node;
		cnode = cnode->next;
	}
	fclose(fin);
}


/* BREAK INTEGENIC REGION INTO BINSIZE BINS
 * IN THIS APPROACH, WE MAKE SURE ALL BINS ARE IN BINSIZE,
 * WE IGNORE BINS THAT ARE SMALLER THAN BINSIZE FOR INTEGENIC REGION*/
void
breakIntegenicRegion(){
	FILE *fin;
	FILE *fout;
	int32_t read_entrezid, read_start, read_end, read_length;
	char chr_id[STR];
	int32_t	binnum = 0;
	fin = fopen(ir_file, "r");
	fout = fopen(binned_ir_file, "w");

	while (fscanf(fin, "%d %s %d %d %d",&read_entrezid, chr_id, &read_start, &read_end, &read_length) != -1){

		if ((read_end - read_start +1) >= binsize){
			// need to break this integenic region
			// break it down into bins with bins being BINSIZE length and last being [binsize, 2*binsize)
			while ((read_end - read_start + 1) >= binsize){
				fprintf(fout, "%d\t%s\t%d\t%d\n",binnum, chr_id, read_start, (read_start+binsize-1));
				read_start = read_start + binsize;
				binnum++;
			}
		}
	}
	fclose(fin);
	fclose(fout);
}



void
calculateExonGCContent(void){
	int i,j;
	FILE *fseq;
	int32_t pos;
	char line[200];
	int len;
	node *bin;
	int linenum;
	char *filename;

	for (i=0; i<chr_num; i++){
		filename = (char *)malloc(STR*sizeof(char));
		strcpy(filename, "human_sequence_data/hs_ref_GRCh37_");
		strcat(filename,chrs[i]);
		strcat(filename,".fa");


		fseq = fopen(filename, "r");
		fgets(line, sizeof(line), fseq);
		bin = exon[i].list->next;
		pos = 0;
		linenum = 0;
		while ( fgets ( line, sizeof(line), fseq) != NULL ){
			linenum++;
			len = strlen(line)-1;
			for (j=0;j<len;j++){
				pos++;
				if (bin != NULL){
					if ((pos >= bin->start) && (pos <= bin->end)){
						switch (toupper(line[j])) {
							case 'N': 	bin->nnum++;
										break;
							case 'G': 	bin->gcnum++;
										break;
							case 'C': 	bin->gcnum++;
										break;
							case 'A': 	bin->atnum++;
										break;
							case 'T': 	bin->atnum++;
										break;
							default:	break;
						}
					} else {
						if (pos > bin->end){
							bin = bin->next;
						}
					} // end-else
				} // end-if
			} // end-for
		} // end-while
	} // end-for

	/*remove heavyN exons*/
	node *to_remove;
	for (i=0; i<chr_num; i++){
		bin = exon[i].list;
		while (bin->next != NULL){
			if ((1.0*bin->next->nnum/(bin->next->end - bin->next->start + 1)) > NPERCENTAGE){
				// this bin contains too much N, remove it
				to_remove = bin->next;
				bin->next = bin->next->next;
				free(to_remove);
			} else {
				bin = bin->next;
			}
		}
	}


	/* Output result to file*/
	char *exon_GC_file = "exon_GC.txt";
	FILE *f_annotation_exon = fopen(exon_GC_file, "w");
	node *dh;
	//fprintf(f_annotation_exon,"entrezid\tchr\tstart\tend\tnnum\tgcnum\tatnum\n");
	for(i=0; i<chr_num; i++){
		dh = exon[i].list;
		while (dh->next != NULL){
			dh = dh->next;
			fprintf(f_annotation_exon, "%d\t%s\t%d\t%d\t%d\t%d\t%d\n", dh->gene, exon[i].id, dh->start, dh->end, dh->nnum, dh->gcnum, dh->atnum);
		}
	}
	fclose(f_annotation_exon);
}

void
calculateIRGCContent(void){

	int i,j;
	FILE *fseq;
	int32_t pos;
	char line[200];
	int len;
	node *bin;
	int linenum;
	char *filename;

	for (i=0; i<chr_num; i++){
		filename = (char *)malloc(STR*sizeof(char));
		strcpy(filename, "human_sequence_data/hs_ref_GRCh37_");
		strcat(filename,chrs[i]);
		strcat(filename,".fa");

		fseq = fopen(filename, "r");
		fgets(line, sizeof(line), fseq);
		bin = ir[i].list->next;
		pos = 0;
		linenum = 0;
		while ( fgets ( line, sizeof(line), fseq) != NULL ){
			linenum++;
			len = strlen(line)-1;
			for (j=0;j<len;j++){
				pos++;
				if (bin != NULL){
					if ((pos >= bin->start) && (pos <= bin->end)){
						switch (toupper(line[j])) {
							case 'N': 	bin->nnum++;
										break;
							case 'G': 	bin->gcnum++;
										break;
							case 'C': 	bin->gcnum++;
										break;
							case 'A': 	bin->atnum++;
										break;
							case 'T': 	bin->atnum++;
										break;
							default:	break;
						}
					} else {
						if (pos > bin->end){
							bin = bin->next;
						}
					} // end-else
				} // end-if
			} // end-for
		} // end-while
	} // end-for

	/*remove heavyN integenic region bins*/
	node *to_remove;
	for (i=0; i<chr_num; i++){
		bin = ir[i].list;
		while (bin->next != NULL){
			if ((1.0*bin->next->nnum/(bin->next->end - bin->next->start + 1)) > NPERCENTAGE){
				// this bin contains too much N, remove it
				to_remove = bin->next;
				bin->next = bin->next->next;
				free(to_remove);
			} else {
				bin = bin->next;
			}
		}
	}


	/* Output result to file*/
	char *ir_GC_file = "binned_integenic_region_GC.txt";
	FILE *f_annotation_ir = fopen(ir_GC_file, "w");
	node *dh;
	//fprintf(f_annotation_ir,"chr\tstart\tend\tnnum\tgcnum\tatnum\n");
	for(i=0; i<chr_num; i++){
		dh = ir[i].list;
		while (dh->next != NULL){
			dh = dh->next;
			fprintf(f_annotation_ir, "%s\t%d\t%d\t%d\t%d\t%d\n", ir[i].id, dh->start, dh->end, dh->nnum, dh->gcnum, dh->atnum);
		}
	}
	fclose(f_annotation_ir);

}


void
DetectionCallAnnotationBody(void){

	breakIntegenicRegion();
	build_ir_data_structure();
	calculateIRGCContent();

	build_exon_data_structure();
	calculateExonGCContent();

}


void
detectionCallAnnotation(char **exonfile, char **irfile, int *user_binsize){

	// initialisation
	
	/* exon and ir files are pre-processed without unwanted chromosome data */
	exon_file = (char *)malloc(STR);
	strcpy(exon_file, *exonfile);
	ir_file = (char *)malloc(STR);
	strcpy(ir_file, *irfile);

	binned_ir_file = (char *)malloc(STR);
	strcpy(binned_ir_file, "binned_integenic_region.txt");
	binsize = *user_binsize;

	DetectionCallAnnotationBody();

	/*  printf("\nUSAGE: DetectionCallAnnotation(exonFilename, integenicRegionFilename, binsizse)\n\n");
		printf("BY DEFAULT:     exonFilename = exon.txt\n\t\tintegenicRegionFilename = integenic_region.txt\n\t\tbinsize = 2000\n\n");	*/
}
