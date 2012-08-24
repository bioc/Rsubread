#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

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

}

void R_junction_wrapper(int * nargs, char ** argv)
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

        main_junction(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

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

}


void R_SNPcalling_wrapper(int * nargs, char ** argv)
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

        main_snp_calling_test(n,c_argv);

        for(i=0;i<n;i++) free(c_argv[i]);
        free(c_argv);
        free(r_argv);

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




