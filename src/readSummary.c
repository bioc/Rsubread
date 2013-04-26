#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <zlib.h>
#include <pthread.h>
#include "subread.h"
#include "sambam-file.h"

#define readSummary_multi_threads readSummary
//#define readSummary_single_thread readSummary

int readSummary_single_thread(int argc,char *argv[]){

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
9: as.numeric(isGeneLevel)
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

int isSAM;
char * ret;
SamBam_FILE * fp_in_bam;

int nhits, j, prev_gid, isMultiOverlapAllowed, isGeneLevel, flag_overlap;
long hits_indices[1000];

line = (char*)calloc(MAX_LINE_LENGTH, 1);

isPE = atoi(argv[4]);
minPEDistance = atoi(argv[5]);
maxPEDistance = atoi(argv[6]);

isSAM = atoi(argv[7]);
isMultiOverlapAllowed = atoi(argv[8]);
isGeneLevel = atoi(argv[9]);
 	
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
  if(*read_chr == '*'){ //better to use flag field to decide if it is mapped or not
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
	
  if(flag == 1){
    nhits = 0;
    for(i=search_start;i<=search_end;i++){
      if(isPE == 1){
        if(read_pos < mate_pos)
		  pos_leftmost = read_pos;
		else
	      pos_leftmost = mate_pos;

		if (start[i] > (pos_leftmost + fragment_length - 1)) break;
		if (stop[i] >= pos_leftmost){
			hits_indices[nhits] = i;
			nhits++;
		} 
		
        //if(pos_leftmost >= (start[i]-fragment_length+1) && pos_leftmost <= stop[i]){
        //  nreads_mapped_to_exon++;
        //  nreads[i]++;
        //  break;
        //}
		
      }
      else{
		if (start[i] > (read_pos + read_length -1)) break;
		if (stop[i] >= read_pos){
			hits_indices[nhits] = i;
			nhits++;
		} 
		
        //if(read_pos >= (start[i]-read_length+1) && read_pos <= stop[i]){
		//  nreads_mapped_to_exon++;
        //  nreads[i]++;
        //  break;
        //}
		
	  } //end else
    } //end for i from search start to search end
	
	if (nhits > 0){
	  if (nhits == 1){
	    nreads_mapped_to_exon++;
	    nreads[hits_indices[0]]++; 
	  }
	  else { // nhits greater than 1	    		
		if (isMultiOverlapAllowed == 1){
		  for (j=0;j<nhits;j++){
		    nreads[hits_indices[j]]++;
		  }
		  nreads_mapped_to_exon++;
		}
		else { // multi-overlap is not allowed		
		  if (isGeneLevel == 1){
		    prev_gid = geneid[hits_indices[0]];
			flag_overlap = 0;
			for (j=1;j<nhits;j++){
			  if (geneid[hits_indices[j]] != prev_gid){
			    flag_overlap = 1;
				break;
			  } 
			}
			
			if (flag_overlap == 0){ //overlap multiple exons from a single gene
			  nreads_mapped_to_exon++;
			  nreads[hits_indices[0]]++;
			}
		  } //end if isGeneLevel equal to 1
		} //end else multi-overlap not allowed
	  } //end else nhits greater than 1
	} // end if nhits greater than 0	
  } //end if flag equal to 1


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

















/********************************************************************/
/********************************************************************/
/********************************************************************/
//  NEW FUNCTION FOR MULTI-THREADING
/********************************************************************/
/********************************************************************/
/********************************************************************/

typedef struct
{
	unsigned short thread_id;
	char * line_buffer1;
	char * line_buffer2;
	unsigned long long int nreads_mapped_to_exon;
	unsigned short current_read_length1;
	unsigned short current_read_length2;
	unsigned int * count_table;
	pthread_t thread_object;

	char * input_buffer;
	unsigned int input_buffer_remainder;
	unsigned int input_buffer_write_ptr;	
	pthread_spinlock_t input_buffer_lock;
} fc_thread_thread_context_t;


typedef struct
{
	int is_gene_level;
	int is_paired_end_data;
	int is_multi_overlap_allowed;
	int min_paired_end_distance;
	int max_paired_end_distance;
	int read_length;
	int line_length;

	unsigned short thread_number;
	fc_thread_thread_context_t * thread_contexts;
	int is_all_finished;
	unsigned int input_buffer_max_size;

	int exontable_nchrs;
	int exontable_exons;
	int * exontable_geneid;
	char ** exontable_chr;
	long * exontable_start;
	long * exontable_stop;

	char ** exontable_anno_chrs;
	long * exontable_anno_chr_heads;
	
} fc_thread_global_context_t;

unsigned int tick_time = 1000;

void process_line_buffer(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context)
{

	char * mate_chr = NULL, * read_chr, * line = thread_context -> line_buffer1, *tmp_tok_ptr;
	long read_pos, mate_pos = 0, fragment_length = 0, search_start, search_end, pos_leftmost;
	int flag, i, j, nhits, flag_overlap, prev_gid;
	long hits_indices[1000];

	//printf("L=%s (%d)\n", line, strlen(line));

	// process the current read or read pair
	strtok_r(line,"\t", &tmp_tok_ptr);
	strtok_r(NULL,"\t", &tmp_tok_ptr);
	read_chr = strtok_r(NULL,"\t", &tmp_tok_ptr);

	//skip the read if unmapped (its mate will be skipped as well if paired-end)
	if(*read_chr == '*'){ //better to use flag field to decide if it is mapped or not
		return;	// do nothing if a read is unmapped, or the first read in a pair of reads is unmapped.
	}

	//get mapping location of the read (it could be the first read in a pair)
	read_pos = atoi(strtok_r(NULL,"\t", &tmp_tok_ptr));

	//remove reads which are not properly paired if paired-end reads are used (on different chromsomes or paired-end distance is too big or too small)
	if(global_context -> is_paired_end_data){
		strtok_r(NULL,"\t", &tmp_tok_ptr);
		strtok_r(NULL,"\t", &tmp_tok_ptr);
		mate_chr = strtok_r(NULL,"\t", &tmp_tok_ptr); //get chr which the mate read is mapped to
		mate_pos = atoi(strtok_r(NULL,"\t", &tmp_tok_ptr));
		char * frag_len_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
		fragment_length = abs(atoi(frag_len_str)); //get the fragment length
		//printf("PE_R : %s , '%s':%d > %d || %d < %d\n", mate_chr, frag_len_str, fragment_length, (global_context -> max_paired_end_distance + thread_context->current_read_length1 -1) , fragment_length, (global_context -> min_paired_end_distance + thread_context->current_read_length1 - 1));
		if(strcmp(mate_chr,"=") != 0 || fragment_length > (global_context -> max_paired_end_distance + thread_context->current_read_length1 -1) || fragment_length < (global_context -> min_paired_end_distance + thread_context->current_read_length1 - 1)){
			//the two reads are not properly paired and are skipped
			return;
		}
	} //end if(isPE==1)

	//assign reads or fragments to features 
	flag = 0;
	for(i=0;i<global_context -> exontable_nchrs;i++)
		if(strcmp(read_chr,global_context -> exontable_anno_chrs[i])==0){
			//get chr to which the current read or fragment is mapped and also the searching range
			flag = 1;
			search_start = global_context -> exontable_anno_chr_heads[i];
			search_end = global_context -> exontable_anno_chr_heads[i+1] - 1;
			break;
		}

	if(flag == 1){
		nhits = 0;
		for(i=search_start;i<=search_end;i++){
			if(global_context -> is_paired_end_data == 1){
				if(read_pos < mate_pos)
					pos_leftmost = read_pos;
				else
					pos_leftmost = mate_pos;

				if (global_context -> exontable_start[i] > (pos_leftmost + fragment_length - 1)) break;
				if (global_context -> exontable_stop[i] >= pos_leftmost){
					hits_indices[nhits] = i;
					nhits++;
				} 

				//if(pos_leftmost >= (start[i]-fragment_length+1) && pos_leftmost <= stop[i]){
				//  nreads_mapped_to_exon++;
				//  nreads[i]++;
				//  break;
				//}

			}
			else{
				if (global_context -> exontable_start[i] > (read_pos + thread_context->current_read_length1 -1)) break;
				if (global_context -> exontable_stop[i] >= read_pos){
					hits_indices[nhits] = i;
					nhits++;
				} 

				//if(read_pos >= (start[i]-read_length+1) && read_pos <= stop[i]){
				//  nreads_mapped_to_exon++;
				//  nreads[i]++;
				//  break;
				//}

			} //end else
		} //end for i from search start to search end

		if (nhits > 0){
			if (nhits == 1){
				thread_context->nreads_mapped_to_exon++;
				thread_context->count_table[hits_indices[0]]++; 
			}
			else { // nhits greater than 1	    		
				if (global_context -> is_multi_overlap_allowed == 1){
					for (j=0;j<nhits;j++){
						thread_context->count_table[hits_indices[j]]++;
					}
					thread_context->nreads_mapped_to_exon++;
				}
				else { // multi-overlap is not allowed		
					if (global_context -> is_gene_level == 1){
						prev_gid = global_context -> exontable_geneid[hits_indices[0]];
						flag_overlap = 0;
						for (j=1;j<nhits;j++){
							if (global_context -> exontable_geneid[hits_indices[j]] != prev_gid){
								flag_overlap = 1;
								break;
							} 
						}

						if (flag_overlap == 0){ //overlap multiple exons from a single gene
							thread_context->nreads_mapped_to_exon++;
							thread_context->count_table[hits_indices[0]]++;
						}
					} //end if isGeneLevel equal to 1
				} //end else multi-overlap not allowed
			} //end else nhits greater than 1
		} // end if nhits greater than 0	
	} //end if flag equal to 1

	// Note that we actually make NO use of the second read in a pair.
	// All information we need is in the first line.
}

void * feature_count_worker(void * vargs)
{
	void ** args = (void **) vargs;

	fc_thread_global_context_t * global_context = args[0];
	fc_thread_thread_context_t * thread_context = args[1];

	free(vargs);

	while (1)
	{
		while(1)
		{
			int is_retrieved = 0;
			pthread_spin_lock(&thread_context->input_buffer_lock);
			if(thread_context->input_buffer_remainder)
			{
				int is_second_read;
				unsigned int buffer_read_bytes ;
				unsigned int buffer_read_ptr;
				if(thread_context->input_buffer_remainder <= thread_context->input_buffer_write_ptr)
					buffer_read_ptr = thread_context->input_buffer_write_ptr - thread_context->input_buffer_remainder; 
				else
					buffer_read_ptr = thread_context->input_buffer_write_ptr + global_context->input_buffer_max_size - thread_context->input_buffer_remainder;

				for(is_second_read = 0; is_second_read < (global_context->is_paired_end_data ? 2:1); is_second_read++)
				{
					char * curr_line_buff = is_second_read?thread_context -> line_buffer2:thread_context -> line_buffer1;
					//printf("R=%llu + %u\n", global_context->input_buffer, buffer_read_ptr);
					for(buffer_read_bytes=0; ; buffer_read_bytes++)
					{
						char nch =  thread_context->input_buffer[buffer_read_ptr ++];
						curr_line_buff[buffer_read_bytes] = nch;
						if(buffer_read_ptr == global_context->input_buffer_max_size)
							buffer_read_ptr = 0; 
						if(nch=='\n'){
							curr_line_buff[buffer_read_bytes+1]=0;
							break;
						}
					}
					thread_context->input_buffer_remainder -= buffer_read_bytes + 1;
				}
				is_retrieved = 1;

			}

			pthread_spin_unlock(&thread_context->input_buffer_lock);
			if(global_context->is_all_finished && !is_retrieved) return NULL;

			if(is_retrieved) break;
			else
				usleep(tick_time);
		}


		thread_context -> current_read_length1 = global_context -> read_length;
		thread_context -> current_read_length2 = global_context -> read_length;

		process_line_buffer(global_context, thread_context);

	}
}

void fc_thread_merge_results(fc_thread_global_context_t * global_context, int * nreads , int *nreads_mapped_to_exon)
{
	int xk1, xk2;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		for(xk2=0; xk2<global_context -> exontable_exons; xk2++)
		{
			nreads[xk2]+=global_context -> thread_contexts[xk1].count_table[xk2];
		}
		printf("The %d-th thread processed %u reads\n", xk1, global_context -> thread_contexts[xk1].nreads_mapped_to_exon);
		(*nreads_mapped_to_exon) += global_context -> thread_contexts[xk1].nreads_mapped_to_exon;
	}
}

int fc_thread_init_global_context(fc_thread_global_context_t * global_context, unsigned int buffer_size, unsigned short threads, int et_nchrs, int et_exons, int * et_geneid, char ** et_chr, long * et_start, long * et_stop, char ** et_anno_chrs, long * et_anno_chr_heads, int line_length , int is_PE_data, int min_pe_dist, int max_pe_dist, int read_length, int is_gene_level, int is_overlap_allowed)
{
	int xk1;

	global_context -> input_buffer_max_size = buffer_size;

	global_context -> is_multi_overlap_allowed = is_overlap_allowed;
	global_context -> is_paired_end_data = is_PE_data;
	global_context -> is_gene_level = is_gene_level;
	global_context -> read_length = read_length;
	global_context -> line_length = read_length;

	global_context -> exontable_nchrs = et_nchrs;
	global_context -> exontable_exons = et_exons;
	global_context -> exontable_geneid = et_geneid;
	global_context -> exontable_chr = et_chr;
	global_context -> exontable_start = et_start;
	global_context -> exontable_stop = et_stop;
	global_context -> exontable_anno_chrs = et_anno_chrs;
	global_context -> exontable_anno_chr_heads = et_anno_chr_heads;
	global_context -> min_paired_end_distance = min_pe_dist;
	global_context -> max_paired_end_distance = max_pe_dist;

	global_context -> is_all_finished = 0;
	global_context -> thread_number = threads;
	global_context -> thread_contexts = malloc(sizeof(fc_thread_thread_context_t) * threads);
	for(xk1=0; xk1<threads; xk1++)
	{
		pthread_spin_init(&global_context->thread_contexts[xk1].input_buffer_lock, PTHREAD_PROCESS_PRIVATE);
		global_context -> thread_contexts[xk1].input_buffer_remainder = 0;
		global_context -> thread_contexts[xk1].input_buffer_write_ptr = 0;
		global_context -> thread_contexts[xk1].input_buffer = malloc(buffer_size);
		global_context -> thread_contexts[xk1].thread_id = xk1;
		global_context -> thread_contexts[xk1].count_table = calloc(sizeof(unsigned int), et_exons);
		global_context -> thread_contexts[xk1].nreads_mapped_to_exon = 0;
		global_context -> thread_contexts[xk1].line_buffer1 = malloc(line_length + 2);
		global_context -> thread_contexts[xk1].line_buffer2 = malloc(line_length);
		if(!global_context ->  thread_contexts[xk1].count_table) return 1;
		void ** thread_args = malloc(sizeof(void *)*2);
		thread_args[0] = global_context;
		thread_args[1] = & global_context -> thread_contexts[xk1];
		pthread_create(&global_context -> thread_contexts[xk1].thread_object, NULL, feature_count_worker, thread_args);
	}

	return 0;
}

void fc_thread_destroy_global_context(fc_thread_global_context_t * global_context)
{
	int xk1;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		free(global_context -> thread_contexts[xk1].count_table);	
		free(global_context -> thread_contexts[xk1].line_buffer1);	
		free(global_context -> thread_contexts[xk1].line_buffer2);	
		free(global_context -> thread_contexts[xk1].input_buffer);
		pthread_spin_destroy(&global_context -> thread_contexts[xk1].input_buffer_lock);
	}
	free(global_context -> thread_contexts);
}
void fc_thread_wait_threads(fc_thread_global_context_t * global_context)
{
	int xk1;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
		pthread_join(global_context -> thread_contexts[xk1].thread_object, NULL);
}

int readSummary_multi_threads(int argc,char *argv[]){

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
	9: as.numeric(isGeneLevel)
	10: as.numeric(nthreads)
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

	int isSAM;
	char * ret;
	SamBam_FILE * fp_in_bam;

	int nhits, j, prev_gid, isMultiOverlapAllowed, isGeneLevel, flag_overlap;
	double time_start = miltime();

	line = (char*)calloc(MAX_LINE_LENGTH, 1);

	isPE = atoi(argv[4]);
	minPEDistance = atoi(argv[5]);
	maxPEDistance = atoi(argv[6]);

	isSAM = atoi(argv[7]);
	isMultiOverlapAllowed = atoi(argv[8]);
	isGeneLevel = atoi(argv[9]);
	unsigned short thread_number;
	if(argc > 10)
		thread_number = atoi(argv[10]);
	else	thread_number = 4;

	if(thread_number<1) thread_number=1;
	if(thread_number>16)thread_number=16;

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


	
	/**********************************************************/
	/**********************************************************/
	// SO FAR THE CODES ARE UNCHANGED.
	// NOW THE NEW CODES: CREATE GLOBAL CONTEXT AND THREAD CONTEXTS.
	// THEN DO featureCount IN THREADS.
	/**********************************************************/
	/**********************************************************/

	unsigned int buffer_size = 1024*1024*6;
	int thread_ret = 0;

	fc_thread_global_context_t global_context;
	thread_ret |= fc_thread_init_global_context(& global_context, buffer_size, thread_number, nchr, nexons, geneid, chr, start, stop, anno_chrs, anno_chr_head, MAX_LINE_LENGTH, isPE, minPEDistance, maxPEDistance, read_length, isGeneLevel, isMultiOverlapAllowed);



	if (isSAM == 1)
		fp_in = fopen(argv[2],"r");
	else
		fp_in_bam = SamBam_fopen(argv[2], SAMBAM_FILE_BAM);

	int buffer_pairs = 16;
	char * preload_line = malloc(sizeof(char) * (2+MAX_LINE_LENGTH)*(isPE?2:1)*buffer_pairs);
	int preload_line_ptr;
	int current_thread_id = 0;

	while (1){
		int pair_no;
		int is_second_read;
		preload_line[0] = 0;
		preload_line_ptr = 0;

		for(pair_no=0; pair_no < buffer_pairs; pair_no++)
		{
			for(is_second_read=0;is_second_read<(isPE?2:1);is_second_read++)
			{
				while(1)
				{
					if (isSAM == 1)
						ret = fgets(line, MAX_LINE_LENGTH, fp_in);
					else
						ret = SamBam_fgets(fp_in_bam, line, MAX_LINE_LENGTH);  
					if(!ret) break;
					if(line[0] != '@') break;
				}

				if(!ret) break;
				int curr_line_len = strlen(line);

				//printf("L %d =%s\n", preload_line_ptr , line);
				strcpy(preload_line+preload_line_ptr, line);
				preload_line_ptr += curr_line_len;
				if(line[curr_line_len-1]!='\n')
				{
					strcpy(preload_line+preload_line_ptr, "\n");
					preload_line_ptr++;
				}
			}
			if(!ret) break;
		}


		int line_length = preload_line_ptr;
		//printf("DL=%s\n" , preload_line);

		if(line_length > 0)
		{
			while(1)
			{
				int is_finished = 0;
				fc_thread_thread_context_t * thread_context = global_context.thread_contexts+current_thread_id;

				pthread_spin_lock(&thread_context->input_buffer_lock);
				unsigned int empty_bytes = global_context.input_buffer_max_size -  thread_context->input_buffer_remainder; 
				if(empty_bytes > line_length)
				{
					unsigned int tail_bytes = global_context.input_buffer_max_size -  thread_context->input_buffer_write_ptr; 
					unsigned int write_p1_len = (tail_bytes > line_length)?line_length:tail_bytes;
					unsigned int write_p2_len = (tail_bytes > line_length)?0:(line_length - tail_bytes);
					memcpy(thread_context->input_buffer + thread_context->input_buffer_write_ptr, preload_line, write_p1_len);
					if(write_p2_len)
					{
						memcpy(thread_context->input_buffer, preload_line + write_p1_len, write_p2_len);
						thread_context->input_buffer_write_ptr = write_p2_len;
					}
					else	thread_context->input_buffer_write_ptr += write_p1_len;
					if(thread_context->input_buffer_write_ptr == global_context.input_buffer_max_size) 
						thread_context->input_buffer_write_ptr=0;


					thread_context->input_buffer_remainder += line_length;
					is_finished = 1;
				}

				pthread_spin_unlock(&thread_context->input_buffer_lock);

				current_thread_id++;
				if(current_thread_id >= thread_number) current_thread_id = 0;

				if(is_finished) break;
				else usleep(tick_time);
			}
		}

		if(!ret) break;

	}

	free(preload_line);
	global_context.is_all_finished = 1;

	fc_thread_wait_threads(&global_context);

	fc_thread_merge_results(&global_context, nreads , &nreads_mapped_to_exon);
	fc_thread_destroy_global_context(&global_context);



	/**********************************************************/
	/**********************************************************/
	// END OF THE NEW CODES
	/**********************************************************/
	/**********************************************************/

	double time_end = miltime();
	if(isPE == 1)
		Rprintf("Number of fragments mapped to the features is: %d\nTime cost = %.1f seconds\n\n", nreads_mapped_to_exon, time_end - time_start);
	else
		Rprintf("Number of reads mapped to the features is: %d\nTime cost = %.1f seconds\n\n", nreads_mapped_to_exon, time_end - time_start);

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
