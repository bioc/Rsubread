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
  char strand;
  char * chr;
  int i, readlen, chr_start, chr_end, flag, mqs;

  int MAX_LINE_LENGTH = 100000;
  
  readlen = atoi(argv[2]);
  
  line = (char*)calloc(MAX_LINE_LENGTH, 1);

  while (fgets(line, MAX_LINE_LENGTH, fp)) {
    if(line[0] == '@')
      continue;

    tok = strtok(line,"\t");
	flag = atoi(strtok(NULL,"\t"));
		
	chr = strtok(NULL,"\t");
    if(chr[0] != '*'){
       chr_start = atoi(strtok(NULL,"\t")) - 1; 
       chr_end = chr_start + readlen;
	   mqs = atoi(strtok(NULL,"\t"));
	   
	   if ((flag & 0x10) == 0){ 
		 strand = '+';
	   }
	   else{
		 strand = '-';
	   }
       fprintf(fp_out,"%s\t%d\t%d\t%s\t%d\t%c\n", chr, chr_start, chr_end, ".",mqs,strand);
    }
  }

  if (line)
    free(line);
 
  fclose(fp);
  fclose(fp_out);
  
  return 0;
}
