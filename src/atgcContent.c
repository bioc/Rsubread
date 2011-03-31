#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <ctype.h>

long int total_app[5];
long int app[5][1000]; 

int get_index(size_t pt){
  int index;
  char* a = (char *)pt;
  switch (*a) {
    case 'A':
	case 'a': index = 0;
              break;
    case 'T':
	case 't': index = 1;
              break;
    case 'G':
	case 'g': index = 2;
              break;
    case 'C':
	case 'c': index = 3;
              break;
    default:  index = 4;
              break;
  }
  return index; 
}


void initialise(){
  int p,q;
  for(p=0; p<5; p++){
    total_app[p] = 0;
    for(q=0;q<1000;q++){
      app[p][q] = 0;
    }
  }
}


void retrieve_sequence(char ** input, char ** output_seq){
  
  char * line = NULL;
  size_t len = 0;
  ssize_t z;
  int i;
  int readlen;
  int offset;
  int line_num = 0;
  int data_ready=1;
  
  FILE *fin,*fout;
  fin = fopen(*input,"r");
  fout = fopen(*output_seq,"w");


	while ((z = getline(&line, &len, fin)) != -1){
		line_num++;
		if(line[0] == '@'){
			z = getline(&line, &len, fin);
			line_num++;
			i=0;
			while ((toascii(*(line+i))!=10) && (toascii(*(line+i))!=32)){
				fprintf(fout, "%c",*(line+i)); 
				i++;
			}
			fprintf(fout,"\n");
		}
	}

	if (line)
		free(line);

  fclose(fin);
  fclose(fout);
}

void atgcContent(char ** input, char ** output, int *basewise){

  FILE *fin, *fout;
  
  char * line = NULL;
  size_t len = 0;
  ssize_t z;
  int i,j;
  int readlen;
  int line_num=0;
  int index;
  double freq[5];
  double total_freq[5];
  
  initialise();
 
  //corresponding index:
  // "A"-0 "T"-1 "G"-2 "C"-3 others-4
  fin = fopen(*input,"r");
  fout = fopen(*output, "w");
  
  fprintf(fout,"A,T,G,C,N\n");
  while ((z = getline(&line, &len, fin)) != -1){
    line_num++;
	i=0;
	while  ((toascii(*(line+i)) != 10) && (toascii(*(line+i)) != 32)){
	   index = get_index((size_t)(line+i));
       app[index][i]++;
       total_app[index]++;
	   i++;
	}
	readlen = i;
  }
  
  
  if (line)
    free(line);
  
  for(i=0; i<5; i++){
    total_freq[i] =(((double)total_app[i]/readlen)/line_num);
  }
  fprintf(fout, "%2.4f,%2.4f,%2.4f,%2.4f,%2.4f\n",total_freq[0], total_freq[1], total_freq[2], total_freq[3], total_freq[4]);
  
  if (*basewise==1){
	for(i=0; i<readlen; i++){
	  for(j=0; j<5; j++){
             freq[j] = ((double)app[j][i])/(double)line_num;
	  }	
	  fprintf(fout,"%2.4f,%2.4f,%2.4f,%2.4f,%2.4f\n",freq[0],freq[1],freq[2],freq[3],freq[4]);
	}
  }
  	
  fclose(fin);
  fclose(fout);
}
