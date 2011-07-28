#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main_align(int argc,char ** argv);
int main_buildindex(int argc,char ** argv);
int sam2bed(int argc,char *argv[]);
int unmapped(int argc,char *argv[]);
int readSummary(int argc,char *argv[]);

void R_buildindex_wrapper(int * nargs, char ** argv)
{
	optind = 1;

	char * r_argv, ** c_argv;
	int i,n;
	
	r_argv = calloc(1000, sizeof(char));
	strcpy(r_argv,*argv);
	
	n = *nargs;
	c_argv = (char **) calloc(n,sizeof(char *));
	for(i=0;i<n;i++) c_argv[i] = malloc(200);
	strcpy(c_argv[0],strtok(r_argv,","));
	for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

	main_buildindex(n,c_argv);

	for(i=0;i<n;i++) free(c_argv[i]);
	free(c_argv);
	free(r_argv);

}

void R_align_wrapper(int * nargs, char ** argv)
{
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;
    
        r_argv = calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);
    	
        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = malloc(200);
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        main_align(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}

void R_sam2bed_wrapper(int * nargs, char ** argv)
{
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = malloc(200);
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        sam2bed(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}


void R_unmapped_wrapper(int * nargs, char ** argv)
{
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = malloc(200);
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        unmapped(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}

void R_readSummary_wrapper(int * nargs, char ** argv)
{
        optind = 1;

        char * r_argv, ** c_argv;
        int i,n;

        r_argv = calloc(1000, sizeof(char));
        strcpy(r_argv,*argv);

        n = *nargs;
        c_argv = (char **) calloc(n,sizeof(char *));
        for(i=0;i<n;i++) c_argv[i] = malloc(200);
        strcpy(c_argv[0],strtok(r_argv,","));
        for(i=1;i<n;i++) strcpy(c_argv[i],strtok(NULL,","));

        readSummary(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

}

