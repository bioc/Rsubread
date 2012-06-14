#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int readSummary(int argc,char *argv[]){
if(argc == 1){
    printf("This function counts the number of reads falling into each exon region and saves the result to file output_file.\nThe order of exons in the output_file is the same as the exon order in the annotation file.\nThe input file should be in SAM format.\nThe annotation file should have already been sorted by chromosome name.\n\n"); 
    printf("Usage: readSummary annotation_file_sorted input_file output_file\n\n");
    exit(0);
}

FILE *fp_ann, *fp_in, *fp_out;

int read_length;
char **chr;
long *start, *stop;
int *geneid, *nreads;

char * line = NULL;
char * read_chr;
size_t len = 0;
ssize_t z;
long i,nexons,read_pos;

int nreads_mapped_to_exon = 0;

long anno_chr_head[500];
char * anno_chrs[500];
long curchr, curpos, search_start, search_end;
char * curchr_name;
int nchr, flag;
curchr = 0;
curpos = 0;
curchr_name = "";
	
/* read in annotation data */
fp_ann = fopen(argv[1],"r");
getline(&line, &len, fp_ann);
nexons = 0;
while ((z = getline(&line, &len, fp_ann)) != -1)
  nexons++;

geneid = (int *) calloc(nexons,sizeof(int));
chr = (char **) calloc(nexons,sizeof(char *));
start = (long *) calloc(nexons,sizeof(long));
stop = (long *) calloc(nexons,sizeof(long));
nreads = (int *) calloc(nexons,sizeof(int));
for(i=0;i<nexons;i++) nreads[i] = 0;

rewind(fp_ann);
getline(&line, &len, fp_ann);

for(i=0;i<nexons;i++){
 getline(&line, &len, fp_ann);
 geneid[i] = atoi(strtok(line,"\t"));
 chr[i] = malloc(41);
 strcpy(chr[i],strtok(NULL,"\t"));
 start[i] = atoi(strtok(NULL,"\t"));
 stop[i] = atoi(strtok(NULL,"\t"));
 
 if(strcmp(curchr_name,chr[i]) != 0){
  curchr_name = chr[i];
  anno_chrs[curchr] = chr[i];
  anno_chr_head[curchr] = curpos;
  curchr++;
 }
 curpos++;
}
nchr = curchr;
anno_chr_head[curchr] = nexons;
fclose(fp_ann);

printf("Number of chromosomes included in the annotation is %\d\n",nchr);

fp_in = fopen(argv[2],"r");
while ((z = getline(&line, &len, fp_in)) != -1){
  if(line[0] == '@')
    continue;
  strtok(line,"\t");
  for(i=0;i<8;i++) strtok(NULL,"\t");
  read_length = strlen(strtok(NULL,"\t"));
  break;
}
fclose(fp_in);

/* summarize reading mapping data */
fp_in = fopen(argv[2],"r");
while ((z = getline(&line, &len, fp_in)) != -1){
  if(line[0] == '@')
    continue;
  strtok(line,"\t");
  strtok(NULL,"\t");
  read_chr = strtok(NULL,"\t");
  if(*read_chr == '*')
    continue;
  read_pos = atoi(strtok(NULL,"\t"));
  
  flag = 0;
  for(i=0;i<nchr;i++)
    if(strcmp(read_chr,anno_chrs[i])==0){
	  flag = 1;
	  search_start = anno_chr_head[i];
	  search_end = anno_chr_head[i+1] - 1;
	  break;
	}
	
  if(flag == 1)
    for(i=search_start;i<=search_end;i++)
      if(read_pos >= (start[i]-read_length+1) && read_pos <= stop[i]){
		nreads_mapped_to_exon++;
        nreads[i]++;
        break;
	  }

}

printf("Number of reads mapped to features is: %d\n", nreads_mapped_to_exon);

/* save the results */
fp_out = fopen(argv[3],"w");
fprintf(fp_out,"entrezid\tchromosome\tchr_start\tchr_stop\tnreads\n");
for(i=0;i<nexons;i++)
  fprintf(fp_out,"%d\t%s\t%d\t%d\t%d\n",geneid[i],chr[i],start[i],stop[i],nreads[i]);

fclose(fp_in);
fclose(fp_out);

free(line);
for(i=0;i<nexons;i++) free(chr[i]);
free(geneid);
free(chr);
free(start);
free(stop);
free(nreads);
}

