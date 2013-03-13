#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int sam2bed(int argc,char *argv[]){

/*
  if(argc == 1){
    printf("Usage: sam2bed -n read_length sam_filename bed_filename\n");
    exit(0);
  }
*/

  FILE *fp, *fp_out;

  fp=fopen(argv[3],"r");
  fp_out=fopen(argv[4],"w");

  char * line = NULL;
  char * tok;
  size_t len = 0;
  ssize_t z;
  char * chr;
  int i, readlen, chr_start, chr_end;

  int MAX_LINE_LENGTH = 100000;
  
  readlen = atoi(argv[2]);
  
  line = (char*)calloc(MAX_LINE_LENGTH, 1);

  while (fgets(line, MAX_LINE_LENGTH, fp)) {
    if(line[0] == '@')
      continue;

    tok = strtok(line,"\t");
    for(i=0;i<2;i++)
       chr = strtok(NULL,"\t");
    if(chr[0] != '*'){
       chr_start = atoi(strtok(NULL,"\t")) - 1; 
       chr_end = chr_start + readlen;
       fprintf(fp_out,"%s\t%d\t%d\n", chr, chr_start, chr_end);
    }
  }

  if (line)
    free(line);
 
  fclose(fp);
  fclose(fp_out);
}
