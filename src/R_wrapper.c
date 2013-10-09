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

int main_junction(int argc,char ** argv);
int main_align(int argc,char ** argv);
int main_buildindex(int argc,char ** argv);
int sam2bed(int argc,char *argv[]);
int propmapped(int argc,char *argv[]);
int readSummary(int argc,char *argv[]);
int main_snp_calling_test(int argc,char *argv[]);
int main_repeated_test(int argc, char *argv[]);

void R_buildindex_wrapper(int * nargs, char ** argv)
{
	optind = 1;

	char * r_argv, ** c_argv;
	int i,n;
	
	r_argv = (char *)calloc(1000, sizeof(char));
	strcpy(r_argv,*argv);
	
	n = *nargs;
	c_argv = (char **) calloc(n,sizeof(char *));
	for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
	strcpy(c_argv[0],strtok(r_argv,","));
	for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

	main_buildindex(n,c_argv);

	for(i=0;i<n;i++) free(c_argv[i]);
	free(c_argv);
	free(r_argv);

}

void R_align_wrapper(int * nargs, char ** argv)
{
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;
    
        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);
    	
        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        main_align(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

	R_CStackLimit = old_cstack_limit;
}

void R_junction_wrapper(int * nargs, char ** argv)
{
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;
    
        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);
    
        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        main_junction(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

	R_CStackLimit = old_cstack_limit;
}


void R_sam2bed_wrapper(int * nargs, char ** argv)
{
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        sam2bed(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}


void R_propmapped_wrapper(int * nargs, char ** argv)
{
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        propmapped(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}

void R_readSummary_wrapper(int * nargs, char ** argv)
{
	//uintptr_t old_cstack_limit = R_CStackLimit;
	//R_CStackLimit =(uintptr_t)-1;
	//printf("RCL=%ld\n", R_CStackLimit);
        optind = 1;
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        readSummary(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

	//R_CStackLimit = old_cstack_limit;
}


void R_SNPcalling_wrapper(int * nargs, char ** argv)
{
	uintptr_t old_cstack_limit = R_CStackLimit;
	R_CStackLimit =(uintptr_t)-1;
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        main_snp_calling_test(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

	R_CStackLimit = old_cstack_limit;
}


void R_removeDupReads_wrapper(int * nargs, char ** argv)
{
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = (char *)calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = (char *)calloc(200,sizeof(char));
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        main_repeated_test(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}




