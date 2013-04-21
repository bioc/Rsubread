#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <zlib.h>
#include "sambam-file.h"


int readSummary(int argc,char *argv[]){

/*
This function counts the number of reads falling into each exon region.
The order of exons in the output is the same as that of exons included in the annotation.
The annotation, if provided as a file, should be sorted by chromosome name.

Parameters passed from the featureCounts R function:
0: "readSummary"
1: ann
2: files[i]
3: fout
4: as.numeric(isPairedEnd)
5: min.distance
6: max.distance
7: as.numeric(tolower(file.type)=="sam")
8: as.numeric(allowMultiOverlap)
*/

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

int MAX_LINE_LENGTH = 100000;
  
int nreads_mapped_to_exon = 0;

long anno_chr_head[500];
char * anno_chrs[500];
long curchr, curpos, search_start, search_end;
char * curchr_name;
int nchr, flag;
curchr = 0;
curpos = 0;
curchr_name = "";

int isPE, minPEDistance, maxPEDistance;
char * mate_chr;
int mate_pos, pos_leftmost, fragment_length;

int isSAM, isMultiOverlapAllowed;
char * ret;
SamBam_FILE * fp_in_bam;

line = (char*)calloc(MAX_LINE_LENGTH, 1);

isPE = atoi(argv[4]);
minPEDistance = atoi(argv[5]);
maxPEDistance = atoi(argv[6]);

isSAM = atoi(argv[7]);
isMultiOverlapAllowed = atoi(argv[8]);
 	
/* read in annotation data */
fp_ann = fopen(argv[1],"r");
if(!fp_ann){
  Rprintf("Failed to open the annotation file %s\n",argv[1]);
  return -1;
}
fgets(line, MAX_LINE_LENGTH, fp_ann);
nexons = 0;
while (fgets(line, MAX_LINE_LENGTH, fp_ann))
  nexons++;

geneid = (int *) calloc(nexons,sizeof(int));
chr = (char **) calloc(nexons,sizeof(char *));
start = (long *) calloc(nexons,sizeof(long));
stop = (long *) calloc(nexons,sizeof(long));
nreads = (int *) calloc(nexons,sizeof(int));
for(i=0;i<nexons;i++) nreads[i] = 0;

rewind(fp_ann);
fgets(line, MAX_LINE_LENGTH, fp_ann);

for(i=0;i<nexons;i++){
 fgets(line, MAX_LINE_LENGTH, fp_ann);
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

Rprintf("Number of chromosomes included in the annotation is \%d\n",nchr);

/* get read length */
if (isSAM == 1){
  fp_in = fopen(argv[2],"r");
  if(!fp_in){
    Rprintf("Failed to open file %s. Please check if the file name and specified file type are correct.\n", argv[2]); 
	return -1;
  }
}
else{
  fp_in_bam = SamBam_fopen(argv[2], SAMBAM_FILE_BAM);
  if(!fp_in_bam){
    Rprintf("Failed to open file %s. Please check if the file name and specified file type are correct.\n", argv[2]); 
	return -1;
  }
}

while (1){
  if (isSAM == 1)
	ret = fgets(line, MAX_LINE_LENGTH, fp_in);
  else
	ret = SamBam_fgets(fp_in_bam, line, MAX_LINE_LENGTH);  

  if(!ret) break;
	
  if(line[0] != '@'){
	strtok(line,"\t");
	for(i=0;i<8;i++) strtok(NULL,"\t");
	read_length = strlen(strtok(NULL,"\t"));
	break;
  }
} //end while

if (isSAM == 1)
  fclose(fp_in);
else
  SamBam_fclose(fp_in_bam);

/* SUMMARIZE READ MAPPING DATA */
if (isSAM == 1)
  fp_in = fopen(argv[2],"r");
else
  fp_in_bam = SamBam_fopen(argv[2], SAMBAM_FILE_BAM);

while (1){
  if (isSAM == 1)
	ret = fgets(line, MAX_LINE_LENGTH, fp_in);
  else
	ret = SamBam_fgets(fp_in_bam, line, MAX_LINE_LENGTH);  

  if(!ret) break;
  if(line[0] == '@') continue;

  // process the current read or read pair
  strtok(line,"\t");
  strtok(NULL,"\t");
  read_chr = strtok(NULL,"\t");

  //skip the read if unmapped (its mate will be skipped as well if paired-end)
  if(*read_chr == '*'){
	if(isPE == 1){
	  if (isSAM == 1)
		fgets(line, MAX_LINE_LENGTH, fp_in);
	  else
		SamBam_fgets(fp_in_bam, line, MAX_LINE_LENGTH);
    }
	continue;
  }
  
  //get mapping location of the read (it could be the first read in a pair)
  read_pos = atoi(strtok(NULL,"\t"));
  
  //remove reads which are not properly paired if paired-end reads are used (on different chromsomes or paired-end distance is too big or too small)
  if(isPE == 1){
    strtok(NULL,"\t");
    strtok(NULL,"\t");
    mate_chr = strtok(NULL,"\t"); //get chr which the mate read is mapped to
    mate_pos = atoi(strtok(NULL,"\t"));
    fragment_length = abs(atoi(strtok(NULL,"\t"))); //get the fragment length
    if(strcmp(mate_chr,"=") != 0 || fragment_length > (maxPEDistance + read_length -1) || fragment_length < (minPEDistance + read_length - 1)){
      //the two reads are not properly paired and are skipped
	  if (isSAM == 1)
		fgets(line, MAX_LINE_LENGTH, fp_in);
	  else
		SamBam_fgets(fp_in_bam, line, MAX_LINE_LENGTH);
	continue;
    }
  } //end if(isPE==1)

  //assign reads or fragments to features 
  flag = 0;
  for(i=0;i<nchr;i++)
    if(strcmp(read_chr,anno_chrs[i])==0){
	  //get chr to which the current read or fragment is mapped and also the searching range
	  flag = 1;
	  search_start = anno_chr_head[i];
	  search_end = anno_chr_head[i+1] - 1;
	  break;
	}
	
  if(flag == 1)
    for(i=search_start;i<=search_end;i++){
      if(isPE == 1){
	  //get the mapping position of leftmost base of the fragment
        if(read_pos < mate_pos)
		  pos_leftmost = read_pos;
		else
	      pos_leftmost = mate_pos;

        if(pos_leftmost >= (start[i]-fragment_length+1) && pos_leftmost <= stop[i]){
          nreads_mapped_to_exon++;
          nreads[i]++;
          break;
        }   
      }
      else
        if(read_pos >= (start[i]-read_length+1) && read_pos <= stop[i]){
		  nreads_mapped_to_exon++;
          nreads[i]++;
          break;
        }
    } //end for

  //if paired end data are used, the current read pair is found properly paired and there is no need to process the second read in the pair 
  if(isPE == 1){
	if (isSAM == 1)
	  fgets(line, MAX_LINE_LENGTH, fp_in);
	else
	  SamBam_fgets(fp_in_bam, line, MAX_LINE_LENGTH);
  }

} //end while 


if(isPE == 1)
  Rprintf("Number of fragments mapped to the features is: %d\n\n", nreads_mapped_to_exon);
else
  Rprintf("Number of reads mapped to the features is: %d\n\n", nreads_mapped_to_exon);

/* save the results */
fp_out = fopen(argv[3],"w");
if(!fp_out){
  Rprintf("Failed to create file %s\n", argv[3]);
  return -1;
}
fprintf(fp_out,"geneid\tchr\tstart\tend\tnreads\n");
for(i=0;i<nexons;i++)
  fprintf(fp_out,"%d\t%s\t%ld\t%ld\t%d\n",geneid[i],chr[i],start[i],stop[i],nreads[i]);

if (isSAM == 1)
  fclose(fp_in);
else
  SamBam_fclose(fp_in_bam);

fclose(fp_out);

free(line);
for(i=0;i<nexons;i++) free(chr[i]);
free(geneid);
free(chr);
free(start);
free(stop);
free(nreads);
}

