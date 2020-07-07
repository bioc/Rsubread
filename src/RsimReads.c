#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <sys/types.h>
#ifndef __MINGW32__
#include <sys/resource.h>
#endif
#include <unistd.h>
#include <sys/stat.h>
#include <locale.h>
#include <ctype.h>
#include <math.h>
#include <zlib.h>

#include <R.h>
#include "subread.h"
#include "HelperFunctions.h"
#include "seek-zlib.h"
#include "input-files.h"

#define TRANSCRIPT_FASTA_LINE_WIDTH 1000
#define MAX_SIMULATION_READ_LEN 250
#define TRANSCRIPT_FASTA_LINE_INIT 3000

#define warn_if_untrue(exp) if(!(exp)){\
  Rprintf("ERROR: unsatisfied assertion in %s at %d\n",__FILE__, __LINE__);\
}

typedef struct {
  int read_length;
  ArrayList * quality_strings;
  ArrayList * transcript_names;
  HashTable * transcript_lengths;
  HashTable * transcript_sequences;
  HashTable * expression_levels;
  char fake_quality_string[MAX_SIMULATION_READ_LEN+3];

  int is_paired_end;
  int simple_transcript_names;
  int truth_in_read_names;
  gzFile out_fps[2];
} RsimReads_context_t;


void Rgrc_sequencing_error_read(char * seq, int qlen, char * qua){
   int b;
   for(b=0; b<qlen; b++){
      if(seq[b]=='N') continue;

      int qub = qua[b];
      float randv = myrand_rand()*1./RAND_MAX;
      float errorp = pow(10,3.3 - qub*.1); // Phred 33
      if(randv < errorp * 1.3333333333){// ATGC random can be the original
         seq[b]="ACGT"[myrand_rand()%4];
      }
   }
}


void Rgen_one_read_here(RsimReads_context_t * grc, int seq_no, int seq_start_pos, int is_PE_second, int trans_negative, unsigned long long rno, int mate_pos){
  warn_if_untrue(seq_no>0);
  warn_if_untrue(seq_no<=grc -> transcript_names -> numOfElements);

  char * seq, *seq_name = ArrayListGet(grc -> transcript_names, seq_no-1);
  char read_seq [grc -> read_length+1];
  int trans_len = HashTableGet(grc->transcript_lengths, seq_name)-NULL-1;

  warn_if_untrue(trans_len>0);
  warn_if_untrue(seq_start_pos + grc -> read_length <= trans_len);
  if(seq_start_pos + grc -> read_length>trans_len){
    Rprintf("ERROR: seq %s has %d bases; wanted %d\n", seq_name, trans_len, seq_start_pos + grc -> read_length );
  }

  seq = HashTableGet( grc -> transcript_sequences, seq_name);
  warn_if_untrue(seq);
 // Rprintf("Extract: chro=%s; pos=%d; seq=%.8s; rseq=%.8s\n" , seq_name , seq_start_pos , seq, seq+seq_start_pos);
  memcpy(read_seq, seq + seq_start_pos, grc -> read_length);
  read_seq[grc -> read_length]=0;
  if(trans_negative) reverse_read(read_seq , grc -> read_length , GENE_SPACE_BASE);

  char * qual_str = NULL;
  if(grc->quality_strings->numOfElements>0){
    qual_str = ArrayListRandom(grc->quality_strings);
    //SUBREADprintf("TESTQUAL\t%p\n", qual_str);
    Rgrc_sequencing_error_read(read_seq, grc -> read_length, qual_str);
  }else{
    if(!grc->fake_quality_string[0]){
      int xx;
      for(xx = 0; xx < grc -> read_length; xx++) grc->fake_quality_string[xx]='H';
      grc->fake_quality_string[xx]=0;
    }
    qual_str = grc->fake_quality_string;
  }

  gzFile thisfp = (is_PE_second==1) ? grc -> out_fps[1] : grc -> out_fps[0];
  if(grc->truth_in_read_names){
    int R1_pos = is_PE_second?mate_pos:seq_start_pos;
    int R2_pos = is_PE_second?seq_start_pos:mate_pos;
    if(is_PE_second<0) gzprintf(thisfp, "@R%09llu:%s:%d\n%s\n+\n%s\n", rno, seq_name, 1+seq_start_pos, read_seq, qual_str);
    else gzprintf(thisfp, "@R%09llu:%s:%d:%d\n%s\n+\n%s\n", rno, seq_name, 1+R1_pos, 1+R2_pos, read_seq, qual_str);
  } else gzprintf(thisfp, "@R%09llu\n%s\n+\n%s\n", rno, read_seq, qual_str);
 // Rprintf("Extract - 2: chro=%s; pos=%d; rseq=%.8s; qual=%.8s\n" , seq_name , seq_start_pos , read_seq, qual_str);
}

void Rgrc_put_new_trans(RsimReads_context_t *grc, char * seq_name, char * lbuf, int lbuf_len){
 // Rprintf("Putting: chro=%s ; seq=%.8s\n", seq_name, lbuf);
  warn_if_untrue(lbuf_len>0);
  HashTablePut( grc->transcript_sequences , seq_name, lbuf );
  HashTablePut( grc->transcript_lengths, strdup(seq_name) , NULL+lbuf_len +1);
}

int init_grc_by_file(RsimReads_context_t *grc, char *fasta_name, char *output_name, char *qualstr_name, char ** input_order_transcript_names, int * input_order_transcript_per_read, int input_transcripts, int read_length, int all_reads, int simplify_names, int truth_in_rnames, int do_PE_reads){
  memset(grc, 0, sizeof(RsimReads_context_t));
  grc->transcript_sequences = StringTableCreate(100000);
  HashTableSetDeallocationFunctions(grc->transcript_sequences, free, free);

  grc -> expression_levels = StringTableCreate(100000);
  HashTableSetDeallocationFunctions(grc->expression_levels, free, NULL);

  grc -> transcript_lengths = StringTableCreate(100000);
  HashTableSetDeallocationFunctions(grc -> transcript_lengths ,free, NULL);
  grc -> transcript_names = ArrayListCreate(100000);
  ArrayListSetDeallocationFunction(grc -> transcript_names,free);

  grc -> quality_strings = ArrayListCreate(50000);
  ArrayListSetDeallocationFunction(grc -> quality_strings,free);
 
  grc -> simple_transcript_names = simplify_names;
  grc -> read_length = read_length;
  grc -> truth_in_read_names = truth_in_rnames;
  grc -> is_paired_end = do_PE_reads;
  autozip_fp auto_FP;
  int xk1, lbuf_cap=0, lbuf_used=0, this_seq_len=-1, total_dup=0;
  char * lbuf=NULL, * seq_name = NULL;

  for(xk1 = 0; xk1 < input_transcripts ; xk1++){
  //  Rprintf("INSERT : %d = %s\n", xk1, input_order_transcript_names[xk1]);
    ArrayListPush(grc -> transcript_names , strdup(input_order_transcript_names[xk1]));
    HashTablePut(grc -> expression_levels, strdup(input_order_transcript_names[xk1]), NULL +1);
    //if(xk1<100) Rprintf("TNAME [%d] = %s\n", xk1, input_order_transcript_names[xk1]);
  }

  for(xk1 = 0; xk1 < all_reads; xk1++){
    warn_if_untrue(input_order_transcript_per_read[xk1]>0);
    warn_if_untrue(input_order_transcript_per_read[xk1] <= input_transcripts);

    int to_display = 0;
    //if(input_order_transcript_per_read[xk1] > 206650) to_display = 1;
    if(to_display)  Rprintf("Inp %d . Req %d . StrP %p =%s ; xk1=%d/%d\n", input_transcripts , input_order_transcript_per_read[xk1], input_order_transcript_names[input_order_transcript_per_read[xk1]-1], input_order_transcript_names[input_order_transcript_per_read[xk1]-1], xk1, all_reads);
    int ov = HashTableGet( grc -> expression_levels, input_order_transcript_names[input_order_transcript_per_read[xk1]-1] ) -NULL;
    warn_if_untrue(ov>0);
    ov++;
    HashTablePutReplace(grc -> expression_levels, input_order_transcript_names[input_order_transcript_per_read[xk1]-1], NULL+ov, 0);
  }


  HelpFuncMD5_CTX md5ctx;
  HelpFuncMD5_Init(&md5ctx);
  HashTable * seq_duplicate_tab = StringTableCreate(100000);
  HashTableSetDeallocationFunctions(seq_duplicate_tab, free, NULL);

  int ret = autozip_open(fasta_name, &auto_FP);
  while(1){
    char clinebuf[TRANSCRIPT_FASTA_LINE_WIDTH];
    int rlength = autozip_gets(&auto_FP, clinebuf, TRANSCRIPT_FASTA_LINE_WIDTH -1);
    if(rlength < 1)break;
    if(rlength >= TRANSCRIPT_FASTA_LINE_WIDTH -1 || clinebuf[rlength]!='\0' || clinebuf[rlength-1]!='\n'){
      SUBREADprintf("Error: The line width of the fasta file excessed %d bytes.\n", TRANSCRIPT_FASTA_LINE_WIDTH);
      ret = 1;
      break;
    }
    if(clinebuf[0]=='>'){
      if(NULL != seq_name){
        char * md5mem = malloc(33); unsigned char md5res[16];
        HelpFuncMD5_Final(md5res, &md5ctx);
        int md5i;for(md5i=0;md5i<16;md5i++)sprintf(md5mem+2*md5i, "%02X", 0xff&(int)md5res[md5i]);
        char * had_tab = HashTableGet(seq_duplicate_tab, md5mem);
        long seq_exp = HashTableGet(grc->expression_levels, seq_name)-NULL;
        if(had_tab && seq_exp>1) total_dup++;//SUBREADprintf("Warning: duplicate sequence was found in '%s' and '%s'.\n", seq_name, had_tab);
        if(seq_exp>1)HashTablePut(seq_duplicate_tab, md5mem, 1+NULL);
        else free(md5mem);

        had_tab = HashTableGet(grc-> transcript_sequences, seq_name);
        if(had_tab){
          SUBREADprintf("Error: duplicate sequence names were found in the input: '%s'.n", seq_name);
          return -1;
        }

        Rgrc_put_new_trans(grc, seq_name, lbuf, this_seq_len);
      }

      clinebuf[rlength-1]=0;
      if(grc->simple_transcript_names)
        for(xk1=1; xk1<rlength-1; xk1++) if(clinebuf[xk1]=='|' || clinebuf[xk1]==' ') clinebuf[xk1]=0;

      seq_name = malloc(strlen(clinebuf));
      if( clinebuf[1]==0 ){
        SUBREADprintf("Error: Every transcript needs a name.\n");
        ret = 1;
        break;
      }
      strcpy(seq_name, clinebuf+1);
      lbuf_used = 0;
      lbuf = malloc(TRANSCRIPT_FASTA_LINE_INIT);
      lbuf_cap = TRANSCRIPT_FASTA_LINE_INIT;
      HelpFuncMD5_Init(&md5ctx);
    }else{
      if(NULL == seq_name){
        SUBREADprintf("Error: The fasta file did not start correctly! \n");
        ret = 1;
        break;
      }
      if(lbuf_cap - lbuf_used < rlength + 1){
        lbuf_cap = max(lbuf_cap *8/5, lbuf_cap + rlength);
        lbuf = realloc(lbuf, lbuf_cap);
      }

      int xx; for(xx=0; xx<rlength-1; xx++) clinebuf[xx] = toupper(clinebuf[xx]);
      //SUBREADprintf("STCP1 : %d used, %d len, %d cap\n", lbuf_used, strlen(clinebuf), lbuf_cap);
      HelpFuncMD5_Update(&md5ctx, clinebuf, rlength-1);
      strcpy(lbuf + lbuf_used, clinebuf );
      *(lbuf+lbuf_used+rlength-1)=0; // '\n' => 0

      lbuf_used += rlength -1;
      this_seq_len = lbuf_used;
    }
  }
  if(lbuf_used<1){
    SUBREADprintf("Error: The fasta file did not end correctly! \n");
    ret = 1;
  }
  if(NULL != seq_name && lbuf_used >0){
    char * md5mem = malloc(33); unsigned char md5res[16];
    HelpFuncMD5_Final(md5res, &md5ctx);
    int md5i;for(md5i=0;md5i<16;md5i++)sprintf(md5mem+2*md5i, "%02X", 0xff&(int)md5res[md5i]);
    long seq_exp = HashTableGet(grc->expression_levels, seq_name)-NULL;
    char * had_tab = HashTableGet(seq_duplicate_tab, md5mem);
    if(had_tab && seq_exp>1)total_dup++;// SUBREADprintf("Warning: duplicate sequence was found in '%s' and '%s'.\n", seq_name, had_tab);
    free(md5mem);

    had_tab = HashTableGet(grc-> transcript_sequences, seq_name);
    if(had_tab){
      SUBREADprintf("Error: duplicate sequence names were found in the input: '%s'.\n", seq_name);
      return -1;
    }

    Rgrc_put_new_trans(grc, seq_name, lbuf, this_seq_len);
  }

  if(total_dup)SUBREADprintf("Warning: there are %d transcripts that have replicate sequences and the wanted expression levels are non-zero. You may use scanFasta() to find their names.\n", total_dup);
  autozip_close(&auto_FP);
  HashTableDestroy(seq_duplicate_tab);

  if(qualstr_name && qualstr_name[0]){
    ret = autozip_open(qualstr_name, &auto_FP);
    if(ret){
      SUBREADprintf("Error: cannot open reference quality file '%s'\n", qualstr_name);
      return -1;
    }
    while(1){
      char linebuf[MAX_SIMULATION_READ_LEN+3];
      int rline = autozip_gets(&auto_FP, linebuf, MAX_SIMULATION_READ_LEN);
      if(rline<1) break;
      if(rline==grc->read_length +1 || ( rline==grc->read_length && linebuf[rline-1]!='\n' ) ) {
      }else{
         SUBREADprintf("Error: all your quality strings must be %d-byte long.\n",read_length);
         return -1;
      }
      char * qstr = malloc(rline);
      memcpy(qstr, linebuf, rline);
      if(qstr[rline-1]=='\n') qstr[rline-1]=0;

      ArrayListPush(grc -> quality_strings, qstr);
    }
    autozip_close(&auto_FP);
  }

  char outname[MAX_FILE_NAME_LENGTH+30];
  sprintf(outname,"%s_R1.fastq.gz", output_name);
  grc->out_fps[0] = gzopen(outname, "wb");

  if(grc->is_paired_end){
    sprintf(outname,"%s_R2.fastq.gz", output_name);
    grc->out_fps[1] = gzopen(outname, "wb");
  }else grc->out_fps[1]=NULL;

  if(grc->out_fps[0] == NULL || ( grc->out_fps[1]== NULL && grc->is_paired_end )){
    SUBREADprintf("Error: cannot create the output file.\n");
    return -1;
  }

  return 0;
}

int destroy_Rsim_context(RsimReads_context_t *grc){
  gzclose(grc->out_fps[0]);
  if(grc->out_fps[1]) gzclose(grc->out_fps[1]);

  HashTableDestroy(grc->transcript_sequences);
  HashTableDestroy(grc->expression_levels);
  HashTableDestroy(grc->transcript_lengths);
  ArrayListDestroy(grc->transcript_names);
  ArrayListDestroy(grc->quality_strings);
  return 0;
}

#define A_LARGE_PRIME_FOR_MOD 24537224085139llu

int simRead_at_main(char *fasta_name, char *output_name, char *qualstr_name, int all_transcripts, char ** trans_names_unique, int *trans_ids, int *start_poses, int *fra_lens, int read_length, int total_reads, int simplify_names, int truth_in_rnames,int do_paired_reads ){
  warn_if_untrue(read_length<=MAX_SIMULATION_READ_LEN);
  warn_if_untrue(total_reads>0);
  warn_if_untrue(all_transcripts>0);

  RsimReads_context_t grc;

  myrand_srand(0); // in R, the seed doesn't matter.
  unsigned long long read_pick_i = myrand_rand()&0xffff;
  read_pick_i = (0x10000llu * read_pick_i) +(myrand_rand()&0xffff);
  read_pick_i = (0x10000llu * read_pick_i) +(myrand_rand()&0xffff);
  read_pick_i = (0x10000llu * read_pick_i) +(myrand_rand()&0xffff);

  int retv = init_grc_by_file(&grc, fasta_name, output_name, qualstr_name, trans_names_unique, trans_ids, all_transcripts, read_length, total_reads, simplify_names, truth_in_rnames, do_paired_reads);
  if(retv==0){
    long long int read_i;
    for(read_i = 0 ; read_i < total_reads; read_i ++){
      read_pick_i %= 1llu * total_reads;

      warn_if_untrue(start_poses[read_pick_i] >0);
      int start_offset = start_poses[read_pick_i] -1; // it is 1-based from R!
      int end_offset = start_offset + fra_lens[read_pick_i];
      int is_R1_at_3End = myrand_rand() % 2;

      int pos_small = start_offset;
      int pos_large = end_offset - read_length;
      Rgen_one_read_here(&grc, trans_ids[read_pick_i], is_R1_at_3End?pos_large:pos_small, grc.is_paired_end?0:-1, is_R1_at_3End, read_i, is_R1_at_3End?pos_small:pos_large);
      if(grc.is_paired_end)Rgen_one_read_here(&grc, trans_ids[read_pick_i], is_R1_at_3End?pos_small:pos_large, 1,!is_R1_at_3End, read_i, is_R1_at_3End?pos_large:pos_small);

      read_pick_i += A_LARGE_PRIME_FOR_MOD;
    }
  }

  retv = retv || destroy_Rsim_context(&grc);
  return 0;
}


