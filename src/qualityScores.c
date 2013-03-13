#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <ctype.h>

void retrieve_scores(char ** input, int *offset_pt, int *size, char **output_temp, char ** output_sc){
  
  char * line = NULL;
  size_t len = 0;
  ssize_t z;
  int i;
  int offset;
  int score;
  int line_num = 0;
  int data_ready=1;
  offset = *offset_pt;
  int n;

  int MAX_LINE_LENGTH = 100000;
  
  char * input_name;
  char * output_temp_name;
  char * output_sc_name; 
  
  FILE *fin, *fout;

  
  input_name = (char *) calloc(1000,sizeof(char));
  output_temp_name = (char *) calloc(1000,sizeof(char));
  output_sc_name = (char *) calloc(1000,sizeof(char));

  strcpy(input_name,*input);
  strcpy(output_temp_name, *output_temp);
  strcpy(output_sc_name, *output_sc);

  fin = fopen(input_name, "r");
  if (!fin)
  {
     Rprintf("Unable to open file: '%s'\n", input_name);
     return;
  }

  fout = fopen(output_temp_name, "w");
  if (!fout)
  {
     Rprintf("Unable to create a temporary file at the current directory: '%s'\n", output_temp_name);
     return;
  }

  line = (char*)calloc(MAX_LINE_LENGTH, 1);
   
  while (fgets(line, MAX_LINE_LENGTH, fin)){
		line_num++;
		if(line[0] == '+'){
			fgets(line, MAX_LINE_LENGTH, fin);
			line_num++;
			score = toascii(*line)-offset;
			fprintf(fout, "%d", score);
			i=1;
			score = toascii(*(line+i))-offset;
			while ((score != (10-offset)) && (score != (32-offset))){
				fprintf(fout, ",%d", score);     
				i++;
				score = toascii(*(line+i))-offset;
			}
			fprintf(fout,"\n");
		}
  }
  

  n = *size;
  // only retrieve n lines out of m lines of input. 
  fclose(fin);
  fclose(fout);
  
  line_num = line_num / 4;
  int printed = 0;
  
  if (n <line_num) {
	int gap = line_num/n;
	int line_number=0;
	FILE *ffin, *ffout;
	ffin = fopen(output_temp_name, "r");
       	if (!ffin)
	{
      	  Rprintf("Unable to open file: '%s'\n", output_temp_name);
          return;
  	}

	ffout = fopen(output_sc_name, "w");
  	if (!ffout)
  	{
     	  Rprintf("Unable to create file: '%s'\n", output_sc_name);
     	  return;
  	}


	while (fgets(line, MAX_LINE_LENGTH, ffin)){
		line_number ++;
		if ((((line_number-1) % gap)==0) && (printed < n)){
		  fprintf(ffout, "%s", line);
		  printed++;
		}
	}

	free(line);
    free(input_name); 
    free(output_temp_name); 
    free(output_sc_name);  

    fclose(ffin);
    fclose(ffout);  

  }

}
