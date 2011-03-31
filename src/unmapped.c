#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int unmapped(int argc,char *argv[]){

  if(argc == 1){
    printf("Usage: unmapped sam_file\n");
    exit(0);
  }

  FILE *fp;
  fp=fopen(argv[1],"r");
  char * line = NULL;
  size_t len = 0;
  ssize_t z;
  int unmapped = 0;
  int totalreads = 0;

  while ((z = getline(&line, &len, fp)) != -1) {
    if(line[0] == '@')
      continue;
    else
      totalreads++;

    if(*(index((index(line,'\t')+1),'\t')+1) == '*')
      unmapped++;
  }

  if (line)
    free(line);
 
  fclose(fp);

  printf("Total number of reads is: %d\n", totalreads);
  printf("Proportion of mapped reads is %f\n\n",1-(float)unmapped/totalreads);
}
