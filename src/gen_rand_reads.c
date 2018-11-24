#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/stat.h>
#include <locale.h>
#include <ctype.h>
#include <math.h>
#include <zlib.h>

#ifndef MAKE_STANDALONE
  #include <R.h>
#endif



#include "subread.h"
#include "core.h"
#include "HelperFunctions.h"
#include "seek-zlib.h"
#include "input-files.h"
#include "gene-algorithms.h"

#define MAX_SIMULATION_READ_LEN 250

static struct option long_options[] =
{
	{"summarizeFasta",  no_argument, 0, 'M'},
	{"transcriptFasta",  required_argument, 0, 't'},
	{"totalReads",  required_argument, 0, 'r'},
	{"pairedEnd",  no_argument, 0, 'p'},
	{"expressionLevels",  required_argument, 0, 'e'},
	{"qualityRefFile",  required_argument, 0, 'q'},
	{"outputPrefix",  required_argument, 0, 'o'},
	{"randSeed",  required_argument, 0, 'S'},
	{"readLen",  required_argument, 0, 'L'},
	{"fragmentLenMean",  required_argument, 0, 'F'},
	{"fragmentLenMin",  required_argument, 0, 'N'},
	{"fragmentLenMax",  required_argument, 0, 'X'},
	{"fragmentLenSigma",  required_argument, 0, 'V'},
	{0, 0, 0, 0}
};

typedef struct {
	char random_seeds[16];
	char transcript_fasta_file[MAX_FILE_NAME_LENGTH];
	char output_prefix[MAX_FILE_NAME_LENGTH];
	char expression_level_file[MAX_FILE_NAME_LENGTH];
	char quality_string_file[MAX_FILE_NAME_LENGTH];

	unsigned long long total_fragments;
	int is_paired_end;
	float fragment_length_mean;
	int fragment_length_max;
	int fragment_length_min;
	float fragment_length_sigma;
	int read_length;

	ArrayList * quality_strings;
	ArrayList * transcript_hitting_space;
	ArrayList * transcript_names;
	HashTable * transcript_sequences;
	HashTable * transcript_lengths;
	HashTable * expression_levels;

	char fake_quality_string[MAX_SIMULATION_READ_LEN];
	char * cmd_line;
	gzFile out_fps[2];
	FILE * counts_out_fp;
} genRand_context_t;

void grc_incrand(genRand_context_t * grc){
	unsigned long long round_rand =0;
	memcpy(&round_rand, grc->random_seeds+8, sizeof(round_rand));
	round_rand++;
	memcpy(grc->random_seeds+8, &round_rand, sizeof(round_rand));
}

void grc_sequencing_error_read(char * seq, int qlen, char * qua){
	int b;
	for(b=0; b<qlen; b++){
		if(seq[b]=='N') continue;

		float randv = rand()*1./RAND_MAX;
		float errorp = pow(10,3.3 - qua[b]*.1); // Phred 33
		if(randv < errorp * 1.3333333333){// ATGC random can be the original
			seq[b]="ACGT"[rand()%4];
		}
	}
}

void gen_one_read_here(genRand_context_t * grc, char * seq, int is_PE_second, int trans_negative, unsigned long long rno){
	char read_seq [grc -> read_length+1];
	memcpy(read_seq, seq, grc -> read_length);
	read_seq[grc -> read_length]=0;
	if(trans_negative) reverse_read(read_seq , grc -> read_length , GENE_SPACE_BASE);

	char * qual_str = NULL;
	if(grc->quality_strings->numOfElements>0){
		qual_str = ArrayListRandom(grc->quality_strings);
		//SUBREADprintf("TESTQUAL\t%p\n", qual_str);
		grc_sequencing_error_read(read_seq, grc -> read_length, qual_str);
	}else{
		if(!grc->fake_quality_string[0]){
			int xx;
			for(xx = 0; xx < grc -> read_length; xx++) grc->fake_quality_string[xx]='X';
			grc->fake_quality_string[xx]=0;
		}
		qual_str = grc->fake_quality_string;
	}

	gzFile thisfp = (is_PE_second==1) ? grc -> out_fps[1] : grc -> out_fps[0];
	gzprintf(thisfp, "@R%09llu\n%s\n+\n%s\n", rno, read_seq, qual_str);
}

void gen_a_read_from_one_transcript(genRand_context_t * grc, long this_transcript_no, unsigned  long long rno){
	char * trans_name = ArrayListGet(grc->transcript_names, this_transcript_no);
	char * trans_seq = HashTableGet(grc->transcript_sequences, trans_name);
	int actual_transcript_len = HashTableGet(grc->transcript_lengths, trans_name) - NULL;
	int applied_fragment_maxlen = min(grc -> fragment_length_max, actual_transcript_len);
	double rand_01 = plain_txt_to_long_rand(grc->random_seeds, 16)*1./0xffffffffffffffffllu;
	int rand_01_int = (int)(rand_01*901267351);
	srand(rand_01_int); // for generating sequencing errors.
	grc_incrand(grc);

	if(grc -> is_paired_end){
		float fragment_len = inverse_sample_normal(rand_01) * grc-> fragment_length_sigma + grc -> fragment_length_mean;
		int fraglen = (int)(min(max(fragment_len, grc -> fragment_length_min), applied_fragment_maxlen));
		rand_01 = plain_txt_to_long_rand(grc->random_seeds, 16)*1./0xffffffffffffffffllu;
		grc_incrand(grc);
		int start_pos = (actual_transcript_len - fraglen) * rand_01;
		int is_first_end_negative = rand_01_int % 2;
		if(is_first_end_negative){
			gen_one_read_here(grc, trans_seq + start_pos + fraglen - grc -> read_length, 0, 1, rno);
			gen_one_read_here(grc, trans_seq + start_pos, 1, 0, rno);
		}else{
			gen_one_read_here(grc, trans_seq + start_pos, 0, 0, rno);
			gen_one_read_here(grc, trans_seq + start_pos + fraglen - grc -> read_length, 1, 1, rno);
		}
	}else{
		int start_pos = (actual_transcript_len - grc -> read_length)*rand_01;
		int is_negative = rand_01_int % 2;
		gen_one_read_here(grc, trans_seq + start_pos, -1, is_negative, rno);
	}
}

int grc_check_parameters(genRand_context_t * grc){
	int ret = 0;

	if(grc->read_length > MAX_SIMULATION_READ_LEN){
		SUBREADprintf("ERROR: the read length cannot be higher than  %d.\n", MAX_SIMULATION_READ_LEN);
		ret=1;
	}
	if(grc->is_paired_end){
		if(grc->fragment_length_min>grc->fragment_length_max){
			SUBREADprintf("ERROR: the minimum fragment length must be equal or higher than the maximum fragment length!\n");
			ret=1;
		}
	
		if(grc->fragment_length_min<grc->read_length){
			SUBREADprintf("ERROR: the minimum fragment length must be equal or higher than read length!\n");
			ret=1;
		}
	
		if(grc->fragment_length_max<1){
			SUBREADprintf("ERROR: the maximum fragment length must be a positive number!\n");
			ret=1;
		}
	}

	if(grc->read_length<1){
		SUBREADprintf("ERROR: the read length must be a positive number!\n");
		ret=1;
	}

	if(!grc->transcript_fasta_file[0]){
		SUBREADprintf("ERROR: a transcript file must be provide!\n");
		ret=1;
	}

	if(!grc->output_prefix[0]){
		SUBREADprintf("ERROR: the output prefix must be provide!\n");
		ret=1;
	}else{
		char outname[MAX_FILE_NAME_LENGTH+30];
		sprintf(outname, "%s.for_test.log",grc->output_prefix);
		FILE * test_out = fopen(outname, "w");
		if(test_out){
			fclose(test_out);
			unlink(outname);
		}else{
			SUBREADprintf("ERROR: cannot create the output file!\n");
			ret=1;
		}
	}

	if(!grc->expression_level_file[0]){
		SUBREADprintf("ERROR: the wanted expression levels must be provided!\n");
		ret=1;
	}

	if(grc->total_fragments < 1){
		SUBREADprintf("WARNING: no read number is specified. Generating one million read%s.\n", grc->is_paired_end?"-pairs":"s");
		grc->total_fragments = 1000000;
	}

	return ret;
}

// The 823,532,653,200th prime is 24,537,224,085,139.
//     -- https://primes.utm.edu/nthprime/index.php
#define A_LARGE_PRIME_FOR_MOD 24537224085139llu


int grc_gen( genRand_context_t *grc ){
	int ret = 0;
	unsigned long long read_i = 0;

	unsigned long long space_end = ArrayListGet(grc->transcript_hitting_space, grc->transcript_hitting_space->numOfElements -1)-NULL;
	ArrayList * num_of_frags_per_transcript = ArrayListCreate(100000);
	unsigned long long lastv = 0, current_total =0;
	ArrayList * rescure_hitting_space = ArrayListCreate(100000);
	unsigned long long to_rescure_read_top=0;
	int min_seq_len = grc->is_paired_end?grc->fragment_length_min:grc->read_length;

	for(read_i = 0; read_i < grc->transcript_hitting_space->numOfElements ; read_i++){
		char *seq_name = ArrayListGet(grc->transcript_names, read_i);
		int seq_len = HashTableGet(grc-> transcript_lengths, seq_name)-NULL;
		unsigned long long thisv = ArrayListGet(grc->transcript_hitting_space, read_i) - NULL;
		unsigned long long this_space_span = thisv - lastv;
		unsigned long long expected_reads =(unsigned long long )((this_space_span *1.0/space_end) * grc->total_fragments*0.99999999);
		unsigned long long to_rescure_reads = (unsigned long long)((this_space_span *1.0/space_end * grc->total_fragments- 1.*expected_reads)*100000.);
		if( seq_len < min_seq_len ){
			to_rescure_reads = 0;
			expected_reads = 0;
		}
		to_rescure_read_top+= to_rescure_reads;
		assert(to_rescure_read_top < 0x5fffffffffffffffllu);
		ArrayListPush(rescure_hitting_space, NULL+to_rescure_read_top);
		ArrayListPush(num_of_frags_per_transcript, NULL+expected_reads);
		current_total += expected_reads;

		lastv = thisv;
	}
	assert(current_total<=grc->total_fragments);
	//SUBREADprintf("TESTRESCURE_NO\t%llu\n", grc->total_fragments - current_total);

	for(read_i = current_total; read_i < grc->total_fragments; read_i++){
		unsigned long long longrand = plain_txt_to_long_rand(grc->random_seeds, 16);
		grc_incrand(grc);

		longrand = longrand % to_rescure_read_top;
		long this_transcript_no = ArrayListFindNextDent(rescure_hitting_space, longrand);
		unsigned long long expected_reads = ArrayListGet(num_of_frags_per_transcript, this_transcript_no)-NULL;
		expected_reads++;
		num_of_frags_per_transcript->elementList[this_transcript_no] = NULL+expected_reads;
	}

	ArrayList * per_transcript_reads_hitting_space = ArrayListCreate(100000);
	unsigned long long total_read_top =0;
	for(read_i =0; read_i < num_of_frags_per_transcript -> numOfElements; read_i++) {
		char *seq_name = ArrayListGet(grc->transcript_names, read_i);
		int seq_len = HashTableGet(grc-> transcript_lengths, seq_name)-NULL;
		unsigned long long expected_reads = ArrayListGet(num_of_frags_per_transcript, read_i)-NULL;
		if(seq_len >= min_seq_len)
			fprintf(grc->counts_out_fp, "%s\t%d\t%llu\n", seq_name, seq_len, expected_reads);
		else
			fprintf(grc->counts_out_fp, "%s\t%d\tNA\n", seq_name, seq_len);
		total_read_top+=expected_reads;
		ArrayListPush(per_transcript_reads_hitting_space, NULL+total_read_top);
	}
	assert(total_read_top == grc->total_fragments);

	if(0)
		for(read_i =0; read_i < num_of_frags_per_transcript -> numOfElements; read_i++) {
			char * trans_name = ArrayListGet(grc->transcript_names, read_i);
			unsigned long long expected_reads = ArrayListGet(num_of_frags_per_transcript, read_i)-NULL;
			long long int xx;
			for(xx =0; xx<expected_reads; xx++) SUBREADprintf("TESTGEN\t%s\n", trans_name);
		}
	else{
		unsigned long long mod_class = total_read_top/2; // an arbitratry starting point.
		for(read_i =0; read_i < grc->total_fragments; read_i++) {
			mod_class += A_LARGE_PRIME_FOR_MOD;
			mod_class = mod_class % grc->total_fragments;
			long this_transcript_no = ArrayListFindNextDent(per_transcript_reads_hitting_space, mod_class);
			//char * trans_name = ArrayListGet(grc->transcript_names, this_transcript_no);
			//SUBREADprintf("TESTGEN\t%s\n", trans_name);
			gen_a_read_from_one_transcript(grc, this_transcript_no, read_i);
		}
	}

	ArrayListDestroy(num_of_frags_per_transcript);
	ArrayListDestroy(rescure_hitting_space);
	ArrayListDestroy(per_transcript_reads_hitting_space);
	return ret;
}

int grc_finalize(genRand_context_t *grc){
	HashTableDestroy(grc->expression_levels);
	HashTableDestroy(grc->transcript_sequences);
	HashTableDestroy(grc->transcript_lengths);
	ArrayListDestroy(grc->quality_strings);
	ArrayListDestroy(grc->transcript_hitting_space);
	ArrayListDestroy(grc->transcript_names);
	gzclose(grc->out_fps[0]);
	if(grc->out_fps[1]) gzclose(grc->out_fps[1]);
	fclose(grc->counts_out_fp);
	free(grc->cmd_line);
	return 0;
}

#define TRANSCRIPT_FASTA_LINE_INIT 800
#define TRANSCRIPT_MAX_EXPRESSION_LEVEL 1000000.000001
#define TRANSCRIPT_FASTA_LINE_WIDTH 300



int grc_summary_fasta(genRand_context_t * grc){
	char outname[MAX_FILE_NAME_LENGTH+30];
	autozip_fp auto_FP;

	if(!grc->output_prefix[0]){
		SUBREADprintf("ERROR: the output prefix must be provide!\n");
		return -1;
	}

	sprintf(outname,"%s.faSummary", grc->output_prefix);
	int ret = autozip_open(grc->transcript_fasta_file, &auto_FP);
	if(ret<0){
		SUBREADprintf("ERROR: cannot open the fasta file as input\n");
		return -1;
	}

	FILE * sumfp = fopen(outname, "w");
	if(sumfp == NULL){
		SUBREADprintf("ERROR: cannot open the putput file\n");
		return -1;
	}
	fprintf(sumfp, "TranscriptID\tLength\n");

	char * seq_name = NULL;
	int seq_len = 0;
	while(1){
		char clinebuf[TRANSCRIPT_FASTA_LINE_WIDTH];
		int rlength = autozip_gets(&auto_FP, clinebuf, TRANSCRIPT_FASTA_LINE_WIDTH -1);
		if(rlength < 1)break;
		if(rlength >= TRANSCRIPT_FASTA_LINE_WIDTH -1 || clinebuf[rlength]!='\0' || clinebuf[rlength-1]!='\n'){
			SUBREADprintf("ERROR: The line width of the fasta file excessed %d bytes!\n", TRANSCRIPT_FASTA_LINE_WIDTH);
			ret = 1;
			break;
		}
		if(clinebuf[0]=='>'){
			if(seq_name){
				fprintf(sumfp, "%s\t%d\n", seq_name, seq_len);
				free(seq_name);
				seq_len = 0;
			}
			seq_name=malloc(rlength);
			clinebuf[rlength-1]=0;
			strcpy(seq_name, clinebuf+1);
		}else{
			seq_len += rlength-1; // no \n
		}
	}

	if(seq_name){
		fprintf(sumfp, "%s\t%d\n", seq_name, seq_len);
		free(seq_name);
	}

	autozip_close(&auto_FP);
	fclose(sumfp);
	return ret;
}

void grc_put_new_trans(genRand_context_t *grc, char * seq_name, char * seq_str, unsigned int seq_len, unsigned long long * linear_space_top){
	if(seq_len<1){
		SUBREADprintf("WARNING: a transcript, '%s', has a zero length. No read is generated from it!\n", seq_name);
	}
	HashTablePut(grc-> transcript_sequences,seq_name, seq_str);
	HashTablePut(grc-> transcript_lengths, seq_name, NULL+ seq_len);
	unsigned long long this_seq_exp_10000 = HashTableGet(grc->expression_levels, seq_name)-NULL;
	if(this_seq_exp_10000<1){
		SUBREADprintf("WARNING: a transcript, '%s', has no wanted expression level. No read is generated from it!\n", seq_name);
		this_seq_exp_10000=0;
	}else this_seq_exp_10000-=1;
	//SUBREADprintf("TESTLEN\t%s\t%d\n", seq_name, seq_len);
	(*linear_space_top) += this_seq_exp_10000 * seq_len;
	ArrayListPush(grc->transcript_names, seq_name);
	ArrayListPush(grc->transcript_hitting_space, NULL+*linear_space_top);
}

int grc_load_env(genRand_context_t *grc){

	grc->expression_levels = HashTableCreate(100000);
	HashTableSetDeallocationFunctions(grc->expression_levels, free, NULL);
	HashTableSetKeyComparisonFunction(grc->expression_levels, fc_strcmp_chro);
	HashTableSetHashFunction(grc->expression_levels, fc_chro_hash);

	grc->transcript_sequences = HashTableCreate(100000);
	HashTableSetDeallocationFunctions(grc->transcript_sequences, free, free);
	HashTableSetKeyComparisonFunction(grc->transcript_sequences, fc_strcmp_chro);
	HashTableSetHashFunction(grc->transcript_sequences, fc_chro_hash);

	grc->transcript_lengths = HashTableCreate(100000);
	HashTableSetKeyComparisonFunction(grc->transcript_lengths, fc_strcmp_chro);
	HashTableSetHashFunction(grc->transcript_lengths, fc_chro_hash);

	grc -> quality_strings = ArrayListCreate(100000);
	ArrayListSetDeallocationFunction(grc -> quality_strings, free);
	grc -> transcript_hitting_space = ArrayListCreate(100000);
	grc -> transcript_names = ArrayListCreate(100000); // the names are destroyed by destroying grc->transcript_sequences

	autozip_fp auto_FP;
	int xk1;
	int ret = autozip_open(grc->expression_level_file, &auto_FP);
	if(ret<0){
		ret = 1;
		SUBREADprintf("ERROR: unable to open the expression level file!\n");
	}else ret = 0;
	if(ret) return ret;

	unsigned long long total_tpm = 0;
	while(1){
		char linebuf[400], * tokbuf=NULL;
		int rline = autozip_gets(&auto_FP, linebuf, 399);
		if(rline<1) break;
		if(strstr(linebuf, "GeneID\tTPM"))continue;
		char * seqname = strtok_r(linebuf, "\t", &tokbuf);
		char * seqexp_str = tokbuf;
		if(NULL == seqexp_str){
			SUBREADprintf("ERROR: expression level file format error!\n");
			ret = 1;
		}
		double seqexp = atof(seqexp_str);
		if(seqexp > TRANSCRIPT_MAX_EXPRESSION_LEVEL){
			SUBREADprintf("ERROR: The transcript expression level shouldn't excess %.0f\n", TRANSCRIPT_MAX_EXPRESSION_LEVEL);
		}
		
		unsigned long long seqexp_int = (unsigned long long )(seqexp*10000.);
		total_tpm += seqexp_int;
		char * seqname_buf = malloc(strlen(seqname)+1);
		strcpy(seqname_buf, seqname);
		HashTablePut(grc->expression_levels, seqname_buf, NULL+seqexp_int+1);
	}
	autozip_close(&auto_FP);
	if(ret) return ret;

	#define ROUNDUP_TOLERANCE ( 1000llu * 10000llu )
	if(total_tpm > 1000000llu * 10000llu + ROUNDUP_TOLERANCE || total_tpm < 1000000llu * 10000llu - ROUNDUP_TOLERANCE) {
		SUBREADprintf("ERROR: total TPM is not 1,000,000\n");
		return 1;
	}




	if(grc->quality_string_file[0]){
		ret = autozip_open(grc->quality_string_file, &auto_FP);
		if(ret<0){
			ret = 1;
			SUBREADprintf("ERROR: unable to open the quality string file!\n");
		}else ret = 0;
		if(ret) return ret;
		while(1){
			char linebuf[400];
			int rline = autozip_gets(&auto_FP, linebuf, 399);
			if(rline<1) break;
			if(rline!=grc->read_length +1) {
				SUBREADprintf("ERROR: all your quality strings must be %d-byte long.\n", grc->read_length);
				ret = 1;
				break;
			}
			char * qstr = malloc(rline);
			memcpy(qstr, linebuf, rline-1);
			qstr[rline-1]=0;
			ArrayListPush(grc -> quality_strings, qstr);
		}
		autozip_close(&auto_FP);
	}

	if(ret) return ret;

	ret = autozip_open(grc->transcript_fasta_file, &auto_FP);
	if(ret<0){
		ret = 1;
		SUBREADprintf("ERROR: unable to open the transcript file!\n");
	} else ret = 0;
	if(ret) return ret;
	
	unsigned long long linear_space_top = 0;
	char * lbuf = NULL, * seq_name = NULL;
	unsigned int lbuf_cap = 0, lbuf_used = 0, this_seq_len = 0;
	while(1){
		
		char clinebuf[TRANSCRIPT_FASTA_LINE_WIDTH];
		int rlength = autozip_gets(&auto_FP, clinebuf, TRANSCRIPT_FASTA_LINE_WIDTH -1);
		if(rlength < 1)break;
		if(rlength >= TRANSCRIPT_FASTA_LINE_WIDTH -1 || clinebuf[rlength]!='\0' || clinebuf[rlength-1]!='\n'){
			SUBREADprintf("ERROR: The line width of the fasta file excessed %d bytes!\n", TRANSCRIPT_FASTA_LINE_WIDTH);
			ret = 1;
			break;
		}
		if(clinebuf[0]=='>'){
			if(NULL != seq_name)
				grc_put_new_trans(grc, seq_name, lbuf, this_seq_len, &linear_space_top);

			for(xk1 = 1; clinebuf[xk1]; xk1++){
				if(clinebuf[xk1]=='\r'||clinebuf[xk1]=='\n'){
					clinebuf[xk1]=0;
					break;
				}
			}

			seq_name = malloc(strlen(clinebuf));
			if( clinebuf[1]==0 ){
				SUBREADprintf("ERROR: Every transcript needs a name!\n");
				ret = 1;
				break;
			}
			strcpy(seq_name, clinebuf+1);
			lbuf_used = 0;
			lbuf = malloc(TRANSCRIPT_FASTA_LINE_INIT);
			lbuf_cap = TRANSCRIPT_FASTA_LINE_INIT;
		}else{
			if(NULL == seq_name){
				SUBREADprintf("ERROR: The fasta file did not start correctly! \n");
				ret = 1;
				break;
			}
			if(lbuf_cap - lbuf_used < rlength + 1){
				lbuf_cap = max(lbuf_cap *8/5, lbuf_cap + rlength);
				lbuf = realloc(lbuf, lbuf_cap);
			}
			//SUBREADprintf("STCP1 : %d used, %d len, %d cap\n", lbuf_used, strlen(clinebuf), lbuf_cap);
			strcpy(lbuf + lbuf_used, clinebuf );
			*(lbuf+lbuf_used+rlength-1)=0; // '\n' => 0

			lbuf_used += rlength -1;
			this_seq_len = lbuf_used;
		}
	}
	if(lbuf_used<1){
		SUBREADprintf("ERROR: The fasta file did not end correctly! \n");
		ret = 1;
	}
	if(NULL != seq_name && lbuf_used >0) grc_put_new_trans(grc, seq_name, lbuf, this_seq_len, &linear_space_top);
	
	autozip_close(&auto_FP);

	char outname[MAX_FILE_NAME_LENGTH+30];

	sprintf(outname,"%s.truthCounts", grc->output_prefix);
	grc->counts_out_fp = fopen(outname,"w");
	fprintf(grc->counts_out_fp, "## CMD :%s\nTranscriptID\tLength\tCount\n", grc->cmd_line);

	sprintf(outname,"%s_R1.fastq.gz", grc->output_prefix);
	grc->out_fps[0] = gzopen(outname, "wb");

	if(grc->is_paired_end){
		sprintf(outname,"%s_R2.fastq.gz", grc->output_prefix);
		grc->out_fps[1] = gzopen(outname, "wb");
	}else grc->out_fps[1]=NULL;

	return ret;
}

#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int gen_rnaseq_reads_main(int argc, char ** argv)
#endif
{
	int do_fasta_summary = 0;
	int c;
	int option_index = 0;
	genRand_context_t grc;
	memset(&grc,0,sizeof(grc));

	optind = 1;
	opterr = 1;
	optopt = 63;

	rebuild_command_line(&grc.cmd_line, argc, argv);
	// default settings of the read/fragment length: 100bp, a general illumina feeling 
	grc.fragment_length_sigma = 30.;
	grc.fragment_length_min = 110;
	grc.fragment_length_max = 400;
	grc.fragment_length_mean = 160.;
	grc.read_length = 100;

	long long seed = -1;

	while ((c = getopt_long (argc, argv, "S:V:N:X:F:L:q:r:t:e:o:pM?", long_options, &option_index)) != -1) {
		switch(c){
			case 'M':
				do_fasta_summary = 1;
				break;
			case 'V':
				grc.fragment_length_sigma = atof(optarg);
				break;
			case 'N':
				grc.fragment_length_min = atoi(optarg);
				break;
			case 'X':
				grc.fragment_length_max = atoi(optarg);
				break;
			case 'F':
				grc.fragment_length_mean = atof(optarg);
				break;
			case 'L':
				grc.read_length = atoi(optarg);
				break;
			case 'p':
				grc.is_paired_end = 1;
				break;
			case 'o':
				strcpy(grc.output_prefix, optarg);
				break;
			case 'e':
				strcpy(grc.expression_level_file, optarg);
				break;
			case 't':
				strcpy(grc.transcript_fasta_file, optarg);
				break;
			case 'r':
				grc.total_fragments = atoll(optarg);
				break;
			case 'S':
				seed = atoll(optarg);
				assert(seed>=0);
				break;
			case 'q':
				strcpy(grc.quality_string_file, optarg);
				break;
			default:
			case '?':
				break;
		}
	} 

	if(seed<0){
		double timemil = miltime();
		memcpy(&seed, &timemil, sizeof(seed));
	}
	memcpy(grc.random_seeds, &seed, sizeof(seed));

	if(do_fasta_summary){
		return grc_summary_fasta(&grc);
	}else{
		int ret = grc_check_parameters(&grc);
		ret =  ret || grc_load_env(&grc);
		ret =  ret || grc_gen(&grc);
		return ret || grc_finalize(&grc);
	}
}
