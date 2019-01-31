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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <Rversion.h>
#if (R_VERSION >= R_Version(2,3,0))
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>
#endif 
#include "HelperFunctions.h"

// ========== This part of code is to run everything in thread threads; the main thread is only for screen output ============
struct R_child_thread_run_opt{
  int (*func)(int , char * []);
  char ** args;
  int n;
};


void * R_child_thread_child(void * aa){
  struct R_child_thread_run_opt *opts = aa;
  opts->func(opts->n, opts->args);
  free(opts);
  msgqu_notifyFinish();
  return NULL;
}

void R_child_thread_run(int (*func)(int , char *[]), int n, char **args){
  msgqu_init();
  struct R_child_thread_run_opt *opts = malloc(sizeof(struct R_child_thread_run_opt));
  opts -> func = func;
  opts -> n = n;
  opts -> args = args;
  pthread_t thread;
  pthread_create(&thread, NULL, R_child_thread_child , opts);
  msgqu_main_loop();
  pthread_join(thread, NULL);
  msgqu_destroy();
}

// ========== END: main thread screen output ============

int main_junction(int argc,char ** argv);
int main_align(int argc,char ** argv);
int main_buildindex(int argc,char ** argv);
int sam2bed(int argc,char *argv[]);
int propmapped(int argc,char *argv[]);
int readSummary(int argc,char *argv[]);
int main_snp_calling_test(int argc,char *argv[]);
int main_repeated_test(int argc, char *argv[]);
int main_read_repair(int argc,char *argv[]);
int main_qualityScores(int argc, char *argv[]);
int findCommonVariants(int argc, char *argv[]);
int longread_mapping_R(int argc, char *argv[]);
int TxUniqueMain(int argc, char *argv[]);
int R_flattenAnnotations(int argc, char *argv[]);
int gen_rnaseq_reads_main(int argc, char *argv[]);

void R_txUnique_wrapper(int * nargs, char ** argv){
	char * r_argv, ** c_argv;
	int i,n;

	n = *nargs;
	r_argv = (char *)calloc(15000, sizeof(char));
	strcpy(r_argv,*argv);

	c_argv = (char **) calloc(n+1,sizeof(char *));
	for(i=0;i<1+n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
	strcpy(c_argv[0],"R_txUnique");
	strcpy(c_argv[1],strtok(r_argv,"\t"));
	for(i=2;i<n+1;i++) strcpy(c_argv[i],strtok(NULL,"\t"));
	TxUniqueMain(n+1, c_argv);
	free(r_argv);
	for(i=0;i<n+1;i++) free(c_argv[i]);
	free(c_argv);
}

void R_mergeVCF(int * nargs, char ** argv)
{
	char * r_argv, ** c_argv;
	int i,n;

	n = *nargs;
	r_argv = (char *)calloc(15000, sizeof(char));
	strcpy(r_argv,*argv);

//	printf("N=%d; V=%s\n", n, r_argv);
	c_argv = (char **) calloc(n+1,sizeof(char *));
	for(i=0;i<1+n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
	strcpy(c_argv[0],"R_mergeVCF");
	strcpy(c_argv[1],strtok(r_argv,";"));
	for(i=2;i<n+1;i++) strcpy(c_argv[i],strtok(NULL,";"));

	findCommonVariants(n+1, c_argv);

	free(r_argv);
	for(i=0;i<n+1;i++) free(c_argv[i]);
	free(c_argv);
}

//#define SAVE_R_CS_STACK

void R_sublong_wrapper(int * nargs, char ** argv){
	char * r_argv, ** c_argv;
	int i,n;

	#ifdef SAVE_R_CS_STACK
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;
	#endif

	r_argv = (char *)calloc(1000, sizeof(char));
	strcpy(r_argv,*argv);

	n = *nargs;
	c_argv = (char **) calloc(n,sizeof(char *));
	for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
	strcpy(c_argv[0],strtok(r_argv,","));
	for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

	longread_mapping_R(n,c_argv);

	for(i=0;i<n;i++) free(c_argv[i]);
	free(c_argv);
	free(r_argv);

	#ifdef SAVE_R_CS_STACK
	R_CStackLimit = old_cstack_limit;
	#endif
}



void R_repair_wrapper(int * nargs, char ** argv){
	char * r_argv, ** c_argv;
	int i,n;

	r_argv = (char *)calloc(1000, sizeof(char));
	strcpy(r_argv,*argv);

	n = *nargs;
	c_argv = (char **) calloc(n,sizeof(char *));
	for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
	strcpy(c_argv[0],strtok(r_argv,","));
	for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

	main_read_repair(n,c_argv);

	for(i=0;i<n;i++) free(c_argv[i]);
	free(c_argv);
	free(r_argv);
}

void R_buildindex_wrapper(int * nargs, char ** argv)
{

	char * r_argv, ** c_argv;
	int i,n;
	
	r_argv = (char *)calloc(1000, sizeof(char));
	strcpy(r_argv,*argv);
	
	n = *nargs;
	c_argv = (char **) calloc(n,sizeof(char *));
	for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
	strcpy(c_argv[0],strtok(r_argv,","));
	for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

	R_child_thread_run(main_buildindex,n,c_argv);

	for(i=0;i<n;i++) free(c_argv[i]);
	free(c_argv);
	free(r_argv);

}

void R_align_wrapper(int * nargs, char ** argv)
{
	#ifdef SAVE_R_CS_STACK
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;
	#endif
	
        char * r_argv, ** c_argv;
        int i,n;
    
        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);
    	
        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

		R_child_thread_run(main_align,n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

	#ifdef SAVE_R_CS_STACK
	R_CStackLimit = old_cstack_limit;
	#endif
}

void R_junction_wrapper(int * nargs, char ** argv)
{
	#ifdef SAVE_R_CS_STACK
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;
	#endif

        char * r_argv, ** c_argv;
        int i,n;
    
        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);
    
        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

		R_child_thread_run(main_junction,n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

	#ifdef SAVE_R_CS_STACK
	R_CStackLimit = old_cstack_limit;
	#endif
}


void R_sam2bed_wrapper(int * nargs, char ** argv)
{

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        sam2bed(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}


void R_propmapped_wrapper(int * nargs, char ** argv)
{

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        propmapped(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}


void R_readSummary_wrapper(int * nargs, char ** argv)
{
  #ifdef SAVE_R_CS_STACK
  uintptr_t old_cstack_limit = R_CStackLimit;
  R_CStackLimit =(uintptr_t)-1;
  #endif

  //printf("RCL=%ld\n", R_CStackLimit);

  char * r_argv, ** c_argv;
  int i,n, arg_len;

  arg_len = strlen(*argv);
  r_argv = (char *)calloc(1+arg_len, sizeof(char));
  strcpy(r_argv,*argv);

  n = *nargs;
  c_argv = (char **) calloc(n,sizeof(char *));

  if(strstr(r_argv, ",,")==NULL )
  {
    for(i=0; i<n; i++)
    {
      char * current_arg = strtok(i?NULL:r_argv,",");
      if(current_arg == NULL)
        break;
      arg_len = strlen(current_arg);
      c_argv[i] = (char *)calloc(1+arg_len,sizeof(char));
      strcpy(c_argv[i], current_arg);
    }

    n=i;

    R_child_thread_run(readSummary,n,c_argv);

    for(i=0;i<n;i++) free(c_argv[i]);
  }
  else Rprintf("No input files are provided. \n");
    free(c_argv);
    free(r_argv);

  #ifdef SAVE_R_CS_STACK
  R_CStackLimit = old_cstack_limit;
  #endif

}

void R_SNPcalling_wrapper(int * nargs, char ** argv)
{
	#ifdef SAVE_R_CS_STACK
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;
	#endif

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

    	R_child_thread_run(main_snp_calling_test,n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

	#ifdef SAVE_R_CS_STACK
	R_CStackLimit = old_cstack_limit;
	#endif
}


void R_removeDupReads_wrapper(int * nargs, char ** argv)
{

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        main_repeated_test(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}


void R_qualityScores_wrapper(int * nargs, char ** argv)
{

	char * r_argv, ** c_argv;
	int i,n;
	
	r_argv = (char *)calloc(1000, sizeof(char));
	strcpy(r_argv,*argv);
	
	n = *nargs;
	c_argv = (char **) calloc(n,sizeof(char *));
	for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
	strcpy(c_argv[0],strtok(r_argv,","));
	for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

	main_qualityScores(n,c_argv);

	for(i=0;i<n;i++) free(c_argv[i]);
	free(c_argv);
	free(r_argv);

}

void R_generate_random_RNAseq_reads(int * nargs, char ** argv){
	char * r_argv, ** c_argv;
	int i,n;
	
	r_argv = (char *)calloc(1500, sizeof(char));
	strcpy(r_argv,*argv);
	//fprintf(stderr, "ARGSSS=%s\n", r_argv);
	
	n = *nargs;
	c_argv = (char **) calloc(n,sizeof(char *));
	for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
	strcpy(c_argv[0],strtok(r_argv,","));
	for(i=1;i<n;i++){
		strcpy(c_argv[i],strtok(NULL,","));
		//fprintf(stderr, "ARG_%d=%s\n",i,c_argv[i]);
	}

	gen_rnaseq_reads_main(n,c_argv);

	for(i=0;i<n;i++) free(c_argv[i]);
	free(c_argv);
	free(r_argv);
}
void R_flattenGTF_wrapper(int * nargs, char ** argv){
	char * r_argv, ** c_argv;
	int i,n;
	
	r_argv = (char *)calloc(1000, sizeof(char));
	strcpy(r_argv,*argv);
	
	n = *nargs;
	c_argv = (char **) calloc(n,sizeof(char *));
	for(i=0;i<n;i++) c_argv[i] = (char *)calloc(300,sizeof(char));
	strcpy(c_argv[0],strtok(r_argv,","));
	for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

	R_flattenAnnotations(n,c_argv);

	for(i=0;i<n;i++) free(c_argv[i]);
	free(c_argv);
	free(r_argv);
}
