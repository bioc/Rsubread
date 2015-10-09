/***************************************************************

   The Subread software package is free software package: 
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
#include <signal.h>
#include <dirent.h> 
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include <stdio.h>
#include <assert.h>
#include "input-files.h"
#include "seek-zlib.h"
#include "gene-algorithms.h"
#include "sublog.h"

unsigned int BASE_BLOCK_LENGTH = 15000000;
int is_R_warnned = 0;


FILE * f_subr_open(const char * fname, const char * mode)
{

#if defined(__LP64__) || defined(_LP64) || defined(MACOS)
		return fopen(fname, mode);
#else
		return fopen64(fname, mode);
#endif

}
void fastq_64_to_33(char * qs)
{
	int i=0;
	while(qs[i])
		qs[i++] -= 31;
}

double guess_reads_density(char * fname, int is_sam)
{
	return guess_reads_density_format(fname, is_sam, NULL, NULL, NULL);
}

unsigned long long geinput_file_offset( gene_input_t * input){
	if(input -> file_type == GENE_INPUT_GZIP_FASTQ){

		return ((seekable_zfile_t*)input -> input_fp) -> block_start_in_file_offset + ((seekable_zfile_t*)input -> input_fp) ->in_block_offset * 5/16; // compressed text ~= plain text * 28%
	}else{
		return ftello((FILE*)input ->input_fp);
	}
}

double guess_reads_density_format(char * fname, int is_sam, int * min_phred_score, int * max_phred_score, int * tested_reads)
{
	gene_input_t *ginp = malloc(sizeof(gene_input_t));
	long long int fpos =0, fpos2 = 0;
	int i;
	int max_qual_chr = -1, min_qual_chr = 127;
	char buff[MAX_READ_LENGTH] , qbuf[MAX_READ_LENGTH];

	float retv = 0;

	if(is_sam == 0)
	{
		if(geinput_open(fname, ginp))retv= -1.0;
	}else if(is_sam == 1)
	{
		if(geinput_open_sam(fname, ginp,0))retv= -1.0;
	}else if(is_sam == 2)
	{
		if(geinput_open_sam(fname, ginp,1))retv= -1.0;
	}

	if(retv > -0.1){
		geinput_next_read(ginp, NULL, buff, NULL);

		fpos = geinput_file_offset(ginp);
		for(i=0; i<3000; i++)
		{
			if(geinput_next_read(ginp, NULL, buff, qbuf)<0) break;
			if(qbuf[0])
			{
				int xk=0;
				while(qbuf[xk])
				{
					min_qual_chr = min(min_qual_chr,qbuf[xk]);
					max_qual_chr = max(max_qual_chr,qbuf[xk++]);
				}
			}
			if(tested_reads)
				(*tested_reads) ++;
				
		}

		if(min_phred_score)
		{
			(*min_phred_score) = min_qual_chr;
			(*max_phred_score) = max_qual_chr;

		}	
		fpos2 = geinput_file_offset(ginp) - fpos;
		geinput_close(ginp);

		retv= fpos2*1.0/i;
	}

	free(ginp);
	return retv;
}

int is_gene_char(char c)
{
	//if(c== 'M' || c == 'm' || c == 'U' || c == 'u' || c == 'A' || c=='a' || c=='G' || c=='g' || c=='C' || c=='c' || c=='T' || c=='t' || c=='N' || c=='n')
	if(c=='-' || c == '.' || c == 'N')
		return GENE_SPACE_BASE;
	if((c>='A' && c<'Z') || (c>='a' && c<='z'))
		return GENE_SPACE_BASE;
	if(c>='0' && c<'9')
		return GENE_SPACE_COLOR;
	return 0;
}

long long int guess_gene_bases(char ** files, int file_number)
{
	int i;
	long long int ret = 0;

	for(i=0; i<file_number; i++)
	{
		char * fname = files[i];
		struct stat statbuf;

		if (stat(fname , &statbuf))
		{
			//printf("NOT FOUND!!%s\n", fname);
			return -i-1;
		}

		ret += statbuf.st_size;
		ret -= 150;
		if(ret<0)ret=0;
	}
	return ret * 70 / 71;
}

int geinput_getc(gene_input_t * input){
	if(input -> file_type == GENE_INPUT_GZIP_FASTQ){
		return seekgz_next_char((seekable_zfile_t*)input -> input_fp);
	}else{
		return fgetc((FILE*)input -> input_fp);
	}
}


int read_line_noempty(int max_read_len, gene_input_t * input, char * buff, int must_upper)
{
	int ret =0;
	if(must_upper)
	{
		while(1)
		{
			char ch = geinput_getc(input);
			#ifdef WINDOWS
			if(ch == '\r') continue;
			#endif
			if(ch == EOF) break;
			if(ch == '\n'){
					if(ret)
						break;
			}
			else if(ret < max_read_len-1)
				buff[ret++] = toupper(ch);
		}
	}
	else
	{
		while(1)
		{
			char ch = geinput_getc(input);
			#ifdef WINDOWS
			if(ch == '\r') continue;
			#endif
			if (ch == EOF) break;
			if(ch == '\n'){
					if(ret)
						break;
			}
			else if(ret < max_read_len-1) buff[ret++] = ch;
		}
	
	}
	buff[ret]=0;
	return ret;
}



int read_line(int max_read_len, FILE * fp, char * buff, int must_upper)
{
	int ret =0;
	if(must_upper)
	{
		while(1)
		{
			char ch = fgetc(fp);
			#ifdef WINDOWS
			if(ch == '\r') continue;
			#endif
			if(ch == '\n' || ch == EOF) break;
			if(ret < max_read_len-1)
				buff[ret++] = toupper(ch);
		}
	}
	else
	{
		while(1)
		{
			char ch = fgetc(fp);
			#ifdef WINDOWS
			if(ch == '\r') continue;
			#endif
			if (ch == '\n' || ch == EOF) break;
			buff[ret++] = ch;
		}
	
	}
	buff[ret]=0;
	return ret;
}



int read_line_back(int max_read_len, FILE * fp, char * buff, int must_upper)
{
	int ret =0;
	int started = 0;
	if(must_upper)
	{
		while(1)
		{
			char ch = fgetc(fp);
			if (ch == '\n')
			{
				if (started)break;
				else continue;
			}
			else if(ch == EOF) break;
			else
				started = 1;
			if(ret <max_read_len && ch != '\r')
				if ((ch!=' ' && ch != '\t'))
					buff[ret++] = toupper(ch);
		}
	}
	else
	{
		while(1)
		{
			char ch = fgetc(fp);
			if (ch == '\n')
			{
				if (started)break;
				else continue;
			}
			else if(ch == EOF) break;
			else
				started = 1;
			
			if(ret <max_read_len && ch != '\r')
				buff[ret++] = ch;
		}
	
	}
	buff[ret]=0;
	return ret;
}

int geinput_readline(gene_input_t * input, char * buff, int conv_to_upper)
{
	return read_line(MAX_READ_LENGTH, input -> input_fp, buff, conv_to_upper);
}

int is_read(char * in_buff)
{
	int p=0;
	char c;
	int space_type = GENE_SPACE_BASE;
	while((c=in_buff[p++])!='\0')
	{
		if(c!='\r' && c!='\n'){
			int x = is_gene_char(c);
			if (x == GENE_SPACE_COLOR)
				space_type = GENE_SPACE_COLOR;
			else if(!x) 
				return 0;
		}
	}
	return space_type;
}

int geinput_open_sam(const char * filename, gene_input_t * input, int half_number)
{
	input->input_fp = f_subr_open(filename, "rb");

	strcpy(input->filename, filename);

	if(input->input_fp == NULL)	
		return 1;
	input -> file_type = half_number + GENE_INPUT_SAM_SINGLE;
	while(1){
		char in_buff[3001];
		long long int current_pos = ftello(input -> input_fp);
		int rlen = read_line(3000, input->input_fp, in_buff, 0);
		if(rlen < 1) return 1;

		if(in_buff[0] != '@')
		{
			int x, tab_no = 0;
			char *read_buf=NULL;
			for(x=0; x<rlen; x++)
			{
				if(in_buff[x]=='\t')
				{
					tab_no ++;
					if(tab_no ==9) read_buf = in_buff+x+1;
					if(tab_no ==10) in_buff[x]=0;
					continue;
				}
			}
			if (tab_no<10)return 1;
			input->space_type = is_read(read_buf);
			if (GENE_INPUT_SAM_PAIR_2 != input -> file_type) fseeko(input -> input_fp , current_pos, SEEK_SET);
			input -> read_chunk_start = ftell(input -> input_fp);
			break;
		}
	}	

	return 0;	
}

int geinput_open(const char * filename, gene_input_t * input)
{
	char in_buff[MAX_READ_LENGTH];
	int line_no = 0, ret = 0;
	if(strlen(filename)>298)
		return 1;

	strcpy(input->filename, filename);
	FILE * TMP_FP = f_subr_open(filename, "rb");

	if(TMP_FP == NULL)	
		return 1;

	int id1, id2;
	id1 = fgetc(TMP_FP);
	id2 = fgetc(TMP_FP);

	if(id1 == 31 && id2 == 139) {
		fclose(TMP_FP);
		input->input_fp = malloc(sizeof(seekable_zfile_t));
		input->file_type = GENE_INPUT_GZIP_FASTQ;
		ret = seekgz_open(filename, input->input_fp);
		if(ret == 0){
			int fq_stat = 0;
			for(line_no = 0; line_no < 1000; line_no++){
				int fl = seekgz_gets(input->input_fp, in_buff, 1000);
				if(fl < 1)break;	// EOF
				else if(fl == 1)continue;	// empty line 
				else{		// text line
					if(fq_stat%4 == 1) // read text
					{
						input->space_type = is_read(in_buff);
						break;
					}
					fq_stat ++;
				}
			}
			seekgz_close(input->input_fp);
			seekgz_open(filename, input->input_fp);
		}	
	}else{
		input->input_fp = TMP_FP;
		fseek(input->input_fp, 0, SEEK_SET);
		while (1){
			long long int last_pos = ftello(input->input_fp);
			int rlen = read_line_noempty(MAX_READ_LENGTH, input, in_buff, 0);
			if (rlen<=0){
				ret = 1;
				break;
			}else{
				if(line_no==0 && is_read(in_buff))
				{
					input->file_type = GENE_INPUT_PLAIN;
					input->space_type = is_read(in_buff);
					fseek(input->input_fp,last_pos,SEEK_SET);
					break;
				}
				if(in_buff[0]=='>')
				{
					input->file_type = GENE_INPUT_FASTA;
				//	printf("FILE %s OPENED AS FATSA.\n", filename);
					rlen += read_line(MAX_READ_LENGTH, input->input_fp, in_buff, 0);
					input->space_type = is_read(in_buff);
					
					fseek(input->input_fp,last_pos,SEEK_SET);
					break;
				}
				if(in_buff[0]=='@')
				{
					input->file_type = GENE_INPUT_FASTQ;
					rlen += read_line_noempty(MAX_READ_LENGTH, input, in_buff, 0);
					input->space_type = is_read(in_buff);
					fseek(input->input_fp, last_pos,SEEK_SET);
					break;
				}		
				line_no++;
			}
		}
	}
	input -> read_chunk_start = geinput_file_offset(input);

	if(0 == input->space_type)input->space_type = GENE_SPACE_BASE;
	return ret;
}


int geinput_next_char(gene_input_t * input)
{
	if(input->file_type == GENE_INPUT_FASTA)
	{
		int last_br = 0;
		while (1)
		{
			char nch = fgetc(input->input_fp);
			if (nch <0 && feof(input->input_fp))
				return -2;
			else if (nch < 0 || nch > 126)SUBREADprintf("\nUnrecognised char = #%d\n", nch);

			if (nch == '\r')
			{
				is_R_warnned = 1; 
				SUBREADprintf("The input FASTA file contains \\r characters. This should not result in any problem but we suggest to use UNIX-style line breaks.\n");
				continue;
			}
			if (nch == '\n')
			{
				last_br = 1;
				continue;
			}
			if (nch == ' ' || nch == '\t')
				continue;

			if (nch == '>' && last_br)
			{
				// if this is a new segment

				fseek(input->input_fp, -1 , SEEK_CUR);
				return -1;
			}

			if (is_gene_char(nch))
				return toupper(nch);
			else
			{
				long long int fpos = ftello(input->input_fp);
				int back_search_len =2;
				int is_empty_seq = 0;
				char *out_buf = malloc(2000);

				while( fpos >= back_search_len )
				{
					fseeko(input->input_fp, fpos - back_search_len, SEEK_SET);
					int bc_nch = fgetc(input->input_fp);
					if(bc_nch=='\n')
					{
						if(nch == '>' && back_search_len==2) is_empty_seq=1;
						break;
					}
					back_search_len++;
				}

				fgets(out_buf, 1999,input->input_fp);

				if(is_empty_seq)
				{
					if(strlen(out_buf)>0)
						out_buf[strlen(out_buf)-1]=0;
					SUBREADprintf ("\nEmpty chromosome sequence before '%s'. The file offset is %llu\n",out_buf, fpos);
					free(out_buf);
					return -1;
				}
				else
				{
					SUBREADprintf ("\nUnknown character in the chromosome data: '%c' (ASCII:%02X), ignored. The file offset is %llu\n", nch, nch, fpos);
					SUBREADprintf("%s", out_buf);
					for(; back_search_len>2; back_search_len--)
						SUBREADprintf(" ");
					SUBREADprintf("^\n");

					fseeko(input->input_fp, fpos, SEEK_SET);
					free(out_buf);
					return 'N';
				}
			}		
			last_br = 0;
		}
	}
	else
	{
		SUBREADprintf("Only the FASTA format is accepted for input chromosome data.\n");
		return -3;
	}

}


int geinput_readline_back(gene_input_t * input, char * linebuffer_3000) 
{
	long long int last_pos = ftello(input -> input_fp);
	int ret = read_line(3000, input->input_fp, linebuffer_3000, 0);
	if(ret<1) return -1;
	fseeko(input -> input_fp, last_pos, SEEK_SET);
	return ret;
}

#define SKIP_LINE { nch=' '; while(nch != EOF && nch != '\n') nch = geinput_getc(input); }
#define SKIP_LINE_NOEMPTY {int content_line_l = 0; nch=' '; while(nch != EOF && (nch != '\n' ||! content_line_l)){nch = geinput_getc(input); content_line_l += (nch != '\n');} }

void geinput_jump_read(gene_input_t * input)
{
	char nch=' ';
	if(input->file_type == GENE_INPUT_PLAIN)
		SKIP_LINE
	else if(input->file_type >= GENE_INPUT_SAM_SINGLE)
	{
		while(1)
		{
			nch = fgetc(input->input_fp); 
			if(nch=='@')
				SKIP_LINE
			else break;
		}
		
		SKIP_LINE
		if(input->file_type != GENE_INPUT_SAM_SINGLE)
			SKIP_LINE
	}
	else if(input->file_type == GENE_INPUT_FASTA)
	{
		SKIP_LINE
		while(1)
		{
			SKIP_LINE
			nch = fgetc(input->input_fp); 
			if(nch == EOF)
				break;
			if(nch=='>')
			{
				fseek(input->input_fp, -1, SEEK_CUR);
				break;
			}
		}
	}
	else if(input->file_type == GENE_INPUT_FASTQ)
	{
		SKIP_LINE_NOEMPTY
		SKIP_LINE_NOEMPTY
		SKIP_LINE_NOEMPTY
		SKIP_LINE_NOEMPTY
	}
}

unsigned int read_numbers(gene_input_t * input)
{
	unsigned int ret = 0;
	char nch;
	long long int fpos = ftello(input->input_fp);
	if(input->file_type >= GENE_INPUT_SAM_SINGLE)
	{
		while(1)
		{
			nch = fgetc(input->input_fp);
			if(nch=='@')
				SKIP_LINE
			else break;
		}
	}

	while(1)
	{
		SKIP_LINE
		if(nch==EOF) break;
		ret ++;
	}
	fseeko(input->input_fp, fpos, SEEK_SET);
	if (input->file_type == GENE_INPUT_FASTQ) return ret/4;
	if (input->file_type == GENE_INPUT_FASTA) return ret/2;
	return ret;
}

void geinput_tell(gene_input_t * input, gene_inputfile_position_t * pos){
	if(input -> file_type == GENE_INPUT_GZIP_FASTQ){
		seekgz_tell(( seekable_zfile_t *)input -> input_fp, &pos -> seekable_gzip_position);
	}else{
		pos -> simple_file_position = ftello((FILE *)input -> input_fp);
	}
}

void geinput_seek(gene_input_t * input, gene_inputfile_position_t * pos){
	if(input -> file_type == GENE_INPUT_GZIP_FASTQ){
		seekgz_seek(( seekable_zfile_t *)input -> input_fp, &pos -> seekable_gzip_position);
	}else{
		fseeko((FILE *)input -> input_fp, pos -> simple_file_position, SEEK_SET);
	}
}

int trim_read_inner(char * read_text, char * qual_text, int rlen, short t_5, short t_3)
{

	if(rlen > t_5)
	{
		int xk1;
		for(xk1 = 0; xk1 < rlen - t_5 ; xk1++)
			read_text[xk1] = read_text[xk1+t_5];

		if(qual_text)
			for(xk1 = 0; xk1 < rlen - t_5 ; xk1++)
				qual_text[xk1] = qual_text[xk1+t_5];
	}
	else{
		read_text[0]=0;
		if(qual_text)qual_text[0]=0;
		return 0;
	}

	if(rlen - t_5 > t_3)
	{
		read_text[rlen - t_5 - t_3]=0;
		if(qual_text)qual_text[rlen - t_5 - t_3]=0;
	}
	else{
		read_text[0]=0;
		if(qual_text)qual_text[0]=0;
		return 0;
	}



	return max(0, rlen - t_5 - t_3);
}

long long int tell_current_line_no(gene_input_t * input){
	long long int fpos = ftello(input->input_fp);
	fseeko(input->input_fp,0,SEEK_SET);
	long long ret = 0, fscanpos = 0;
	while(1)
	{
		char nch = fgetc(input->input_fp);
		if(nch == EOF) return -1;
		if(nch == '\n') ret ++;
		fscanpos ++;
		if(fscanpos >= fpos){
			fseeko(input->input_fp, fpos, SEEK_SET);
			return ret;
		}
	}
}

int geinput_next_read(gene_input_t * input, char * read_name, char * read_string, char * quality_string)
{
	return geinput_next_read_trim( input, read_name, read_string,  quality_string, 0, 0, NULL);
}
int geinput_next_read_trim(gene_input_t * input, char * read_name, char * read_string, char * quality_string, short trim_5, short trim_3, int * is_secondary)
{
	if(input -> file_type == GENE_INPUT_PLAIN) {
		int ret = read_line(MAX_READ_LENGTH, input->input_fp, read_string, 0);
		if(quality_string) *quality_string=0;

		if(ret <3)return -1;

		if(trim_5 || trim_3) ret = trim_read_inner(read_string, NULL, ret, trim_5, trim_3);
		return ret;
	} else if(input->file_type >= GENE_INPUT_SAM_SINGLE) {
		char in_buff [3001];
		int tabs;
		int current_str_pos;
		int i;
		int ret = -1;
		int need_reverse;
		char mask_buf[5];



		while(1)
		{
			//	int is_second_map = 0;
				int linelen = read_line(3000, input->input_fp, in_buff, 0);
				if(linelen <1)return -1;
				if(read_name)
					*read_name = 0;
				if(quality_string)
					*quality_string = 0;
				*read_string = 0;
				need_reverse = 0;
				current_str_pos = 0;
				ret = -1;
				tabs=0;

				for(i=0; i<linelen+1; i++)
				{
					if(in_buff[i]=='\t'|| i ==linelen)
					{
						if(tabs == 0 && read_name)read_name[current_str_pos] = 0;
						if(tabs == 1)
						{
							mask_buf[current_str_pos] = 0;
							int flags = atoi(mask_buf) ;
							if(is_secondary && (flags & SAM_FLAG_SECONDARY_MAPPING))
							{
								(*is_secondary) = 1;
							}
							need_reverse = ( flags & SAM_FLAG_REVERSE_STRAND_MATCHED )?1:0;
							
						}
						if(tabs == 9){
							read_string[current_str_pos] = 0;
							ret = current_str_pos;
						}
						if(tabs == 10 && quality_string){
							quality_string[current_str_pos] = 0;
							break;
						}

						current_str_pos = 0 ;
						tabs +=1;
					}
					else
					{
						if(tabs == 9)// read
							read_string[current_str_pos++] = in_buff[i];
						else if(tabs == 10 && quality_string)// quality string
							quality_string[current_str_pos++] = in_buff[i];
						else if(tabs == 0 && read_name)// name
							read_name[current_str_pos++] = in_buff[i];
						else if(tabs == 1)
							mask_buf[current_str_pos++] = in_buff[i];
					}
				}
				if(input->file_type > GENE_INPUT_SAM_SINGLE)
					// skip a line if not single-end
					read_line(1, input->input_fp, in_buff, 0);

				break;
				//printf("Repeated read skipped : %s\n", read_name);
		}

		if(need_reverse)
		{
			if(quality_string)
				reverse_quality(quality_string, ret);
			reverse_read(read_string, ret, input->space_type);
		}
		if(trim_5 || trim_3) ret = trim_read_inner(read_string, quality_string, ret, trim_5, trim_3);
		return ret;
	} else if(input->file_type == GENE_INPUT_FASTA) {
		int ret;
		if(quality_string) (*quality_string)=0;

		while(1) // fetch read name
		{
			ret = read_line(MAX_READ_LENGTH, input->input_fp, read_string, 0);
			if(ret <1)
			{
				sublog_printf(SUBLOG_STAGE_RELEASED,SUBLOG_LEVEL_DEBUG, "The input file normally exhausted.");
				return -1;
			}

			int cursor = 0;
			while(read_string[cursor])
			{
				if(cursor >=2 &&(read_string[cursor] == ' ' || read_string[cursor] == '\t'))
				{
					read_string [cursor] = 0;
					break;	
				}
				cursor++;
			}

			if(read_string[0]=='>'){
				if (read_name != NULL)
					strncpy(read_name, read_string+1, MAX_READ_NAME_LEN);
				break;
			}
			else
				sublog_printf(SUBLOG_STAGE_RELEASED,SUBLOG_LEVEL_FATAL,"The input file may be broken.");
		}
		ret = 0;
		while(1) // fetch read text
		{
			char nch;
			ret += read_line(MAX_READ_LENGTH-ret, input->input_fp, read_string+ret, 1);

			nch = fgetc(input->input_fp);

			if(nch!=EOF)
				fseek(input->input_fp, -1, SEEK_CUR);
		
			if(nch == '>'||nch<1 || nch == EOF)
				break;
		}
//		printf("LOAD R=|%s|\nRETV=%d\n", read_string, ret);
		if(ret <1)return -1;
		if(trim_5 || trim_3) ret = trim_read_inner(read_string, quality_string, ret, trim_5, trim_3);
		return ret;
		
	} else if(input->file_type == GENE_INPUT_FASTQ || input->file_type == GENE_INPUT_GZIP_FASTQ) {
		char nch;
		int ret;

		//READ NAME
		if (read_name == NULL)
		{
			SKIP_LINE_NOEMPTY;
			if(nch == EOF) return -1;
		}
		else
		{
			do{
				nch = geinput_getc(input);
			} while (nch == '\n');
			if(nch==EOF) return -1;
			
			if(nch != '@') {
				long long int lineno = tell_current_line_no(input);
				SUBREADprintf("ERROR: a format issue %c is found on the %lld-th line in input file '%s'!\nProgram aborted!\n", nch, lineno, input -> filename); 
				return -1;
			}

			read_line_noempty(MAX_READ_NAME_LEN, input, read_name, 0);

			int cursor = 1;
			while(read_name[cursor])
			{
				if(read_name[cursor] == ' ' || read_name[cursor] == '\t')
				{
					read_name [cursor] = 0;
					break;	
				}
				cursor++;
			}
		}
		// READ LINE 
		ret = read_line_noempty(MAX_READ_LENGTH, input, read_string, 1);

		// SKIP "+"
		do{
			nch = geinput_getc(input);
		} while( nch == '\n' );
		if(nch != '+'){
			long long int lineno = tell_current_line_no(input);
			SUBREADprintf("ERROR: a format issue %c is found on the %lld-th line in input file '%s'!\nProgram aborted!\n", nch, lineno, input -> filename); 
			return -1;
		}
		SKIP_LINE;

		// QUAL LINE 
		if (quality_string)
			read_line_noempty(MAX_READ_LENGTH, input, quality_string, 0);
		else
			SKIP_LINE_NOEMPTY;



		#ifdef MODIFIED_READ_LEN
		{
			int modified_start = 0;
			if(modified_start)
			{
				int i;
				for(i=0;i<MODIFIED_READ_LEN; i++)
				{
					read_string[i] = read_string[i+modified_start];
					if(quality_string)
						quality_string[i] = quality_string[i+modified_start];
				}
			}
			read_string[MODIFIED_READ_LEN]=0;
			if(quality_string)
				quality_string[MODIFIED_READ_LEN]=0;
			ret = MODIFIED_READ_LEN;
		}
		#endif

//		printf("LOAD R=|%s|\nRETV=%d\n", read_string, ret);
		
		if(trim_5 || trim_3) ret = trim_read_inner(read_string, quality_string, ret, trim_5, trim_3);
		return ret;
		
	}else return -1;
}
void geinput_close(gene_input_t * input)
{
	if(input -> file_type == GENE_INPUT_GZIP_FASTQ)
		seekgz_close((seekable_zfile_t * ) input->input_fp);
	else
		fclose((FILE*)input->input_fp);
}

char * __converting_char_table = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN  ";

void reverse_read(char * InBuff, int read_len, int space_type)
{
	int i;

	if(space_type == GENE_SPACE_COLOR)
	{
		int start_pos = 0;
		char last_base = InBuff[0];

		//printf("CLRLEN0=%d\nS0=%s\n", read_len, InBuff);
		if(isalpha(last_base))
		{
			read_len ++;

			for (i=1; i<read_len; i++)
			{
				int new_int = InBuff[i];
				int new_base = 0;
				if(new_int == '0')
					new_base=last_base;
				else if(new_int == '1')
				{
					if(last_base == 'A')new_base = 'C';
					else if(last_base == 'G')new_base = 'T';
					else if(last_base == 'T')new_base = 'G';
					else new_base = 'A';
				}
				else if(new_int == '2')
				{
					if(last_base == 'A')new_base = 'G';
					else if(last_base == 'G')new_base = 'A';
					else if(last_base == 'T')new_base = 'C';
					else new_base = 'T';
				}
				else
				{
					if(last_base == 'A')new_base = 'T';
					else if(last_base == 'G')new_base = 'C';
					else if(last_base == 'T')new_base = 'A';
					else new_base = 'G';
				}
				last_base = new_base;
			//	putchar(last_base);
			}	
			//puts("");
			InBuff[0] = *(__converting_char_table+last_base);
			start_pos = 1;
		}
		else read_len--;

		for (i=0; i<(read_len - start_pos)/2; i++)
		{
			int rll1 = read_len - 1 - i;
			char tmp = InBuff[rll1];
			InBuff[rll1] = InBuff[i + start_pos];
			InBuff[i + start_pos] = tmp;
		}
	}
	else
	{
		for (i=0; i<read_len/2; i++)
		{
			int rll1 = read_len - 1 - i;
			unsigned char tmp = InBuff[rll1];

			InBuff[rll1] = *(__converting_char_table+InBuff[i]);
			InBuff[i] = *(__converting_char_table+tmp);

		}
		if(i*2 == read_len-1)
		{
			InBuff[i] = *(__converting_char_table+InBuff[i]);
		}
	}

}



void reverse_quality(char * InBuff, int read_len)
{
	int i;
	if(!InBuff) return;
	if(!InBuff[0]) return;
	for (i=0; i<read_len/2; i++)
	{
		char tmp;
		tmp = InBuff[i];
		InBuff[i] = InBuff[read_len -1-i];
		InBuff[read_len -1-i] = tmp;
	}
}


int genekey2intX(char * key,int space_type)
{
	int i;
	int ret;

	ret = 0;
	if(space_type == GENE_SPACE_BASE)
		for (i=30; i>=0; i-=2)
		{
			char kv = *(key++);
			ret |= base2int(kv)<<i;
		}
	else
		for (i=0; i<16; i++)
		{
			ret = ret << 2;
			ret |= color2int (key[i]);
		}

//	printf("RET=%u\n",ret);
	
	return ret;
}


int genekey2int(char key [],int space_type)
{
	int i;
	int ret;

	ret = 0;
	if(space_type == GENE_SPACE_BASE)
		for (i=0; i<16; i++)
		{
			//SUBREADprintf("KV=%c\n", key[i]);
			//ret = ret << 2;
			ret |= (base2int(key[i]))<<(2*(15-i));
		}
	else
		for (i=0; i<16; i++)
		{
			ret = ret << 2;
			ret |= color2int (key[i]);
		}
	
	//SUBREADprintf("RET=%u\n",ret);
	return ret;
}

int genekey2color(char last_base, char key [])
{
	int i, ret = 0;
	char last_char = last_base;

	for (i=0; i<16; i++)
	{
		char next_char = key[i];

		ret = ret << 2;
		ret += chars2color(last_char, next_char);

		last_char = next_char;
	}

	return ret;
}

void colorread2base(char * read_buffer, int read_len)
{
	int i;
	char last_base = read_buffer[0];
	//printf("C2B:%s\n",read_buffer);
	for (i=1; i<read_len; i++)
	{
		int new_int = read_buffer[i];
		int new_base = 0;
		if(new_int == '0')
			new_base=last_base;
		else if(new_int == '1')
		{
			if(last_base == 'A')new_base = 'C';
			else if(last_base == 'G')new_base = 'T';
			else if(last_base == 'T')new_base = 'G';
			else new_base = 'A';
		}
		else if(new_int == '2')
		{
			if(last_base == 'A')new_base = 'G';
			else if(last_base == 'G')new_base = 'A';
			else if(last_base == 'T')new_base = 'C';
			else new_base = 'T';
		}
		else
		{
			if(last_base == 'A')new_base = 'T';
			else if(last_base == 'G')new_base = 'C';
			else if(last_base == 'T')new_base = 'A';
			else new_base = 'G';
		}
		read_buffer[i] = new_base;
		last_base = new_base;
	}
	//printf("CBX:%s\n",read_buffer);
}

char color2char(char clr, char c1)
{
	if(clr == '0')return c1;
	else if(clr == '1')
	{
		if(c1 == 'A') return 'C';
		else if(c1 == 'T') return 'G';
		else if(c1 == 'G') return 'T';
		else return 'A';
	}
	else if(clr == '2') 
	{
		if(c1 == 'A') return 'G';
		else if(c1 == 'T') return 'C';
		else if(c1 == 'G') return 'A';
		else return 'T';
	}
	else if(clr == '3') 
	{
		if(c1 == 'A') return 'T';
		else if(c1 == 'T') return 'A';
		else if(c1 == 'G') return 'C';
		else return 'G';
	}

	return 'N';	
}

int chars2color(char c1, char c2)
{
	if(c1 == 'A')
	{
		if (c2=='A') return 0;
		if (c2=='C') return 1;
		if (c2=='G') return 2;
		else return 3;
	}
	if (c1 == 'C')
	{
		if (c2=='A') return 1;
		if (c2=='C') return 0;
		if (c2=='G') return 3;
		else return 2;
	}
	if (c1 == 'G')
	{
		if (c2=='A') return 2;
		if (c2=='C') return 3;
		if (c2=='G') return 0;
		else return 1;
	}

	// if c1 == 'T', 'U'
	if (c2=='A') return 3;
	if (c2=='C') return 2;
	if (c2=='G') return 1;
	else return 0;



}

int find_subread_end(int len, int TOTAL_SUBREADS, int subread)
{
	if(len<= EXON_LONG_READ_LENGTH)
	{
		int subread_step =  ((len<<16) - (19<<16))/(TOTAL_SUBREADS -1);
		return ((subread_step*(subread))>>16)+15;
	}
	else
	{
		int subread_step;
		
		subread_step = 6<<16;
		if(((len - 18)<<16) / subread_step > 62)
			subread_step = ((len - 18)<<16)/62;
		return ((subread_step*(subread))>>16)+15;
	}
}


//This function returns 0 if the line is a mapped read; -1 if the line is in a wrong format and 1 if the read is unmapped.
int parse_SAM_line(char * sam_line, char * read_name, int * flags, char * chro, unsigned int * pos, char * cigar, int * mapping_quality, unsigned int * pair_dist, char * sequence , char * quality_string, int * rl, int * repeated)
{
	char cc;
	int ci = 0, k=0, field=0, ret_quality = 0, ret_flag = 0, ret_pairdist=0;
	unsigned int ret_pos = 0;
	int is_rep = 0;
	
	while( (cc = sam_line[k]) )
	{
		if(cc=='\t')
		{
			field++;
			k++;
			if(field == 1)read_name[ci]=0;
			else if(field == 3)chro[ci]=0;
			else if(field == 6)cigar[ci]=0;
			else if(field == 10)
			{
				sequence[ci]=0;
				(*rl) = ci;
			}
			else if(field == 11)quality_string[ci]=0;
			ci=0;
			is_rep = 0;
			continue;
		}
		if(field == 9)
			sequence[ci++] = cc;
		else if(field == 10)
			quality_string[ci++] = cc;
		else if(field == 0)
			read_name[ci++] = cc;
		else if(field == 1)
			ret_flag = ret_flag*10 + (cc-'0');
		else if(field == 8)
		{
			if(cc!='-')
				ret_pairdist = ret_pairdist*10 + (cc-'0');
		}
		else if(field == 2)
		{
			//if(ci == 0 && cc == '*') return 1;
			chro[ci++] = cc;
		}
		else if(field == 3)
			ret_pos = ret_pos * 10 + (cc-'0');
		else if(field == 4)
			ret_quality = ret_quality * 10 + (cc-'0');
		else if(field == 5)
			cigar[ci++] = cc;
		else if(field > 10)
		{
			if(cc == 'I' && ci==0) is_rep = 1;
			if(cc != 'H' && ci==1 ) is_rep = 0;
			if(is_rep && ci == 4) *repeated = 0;
			if(is_rep && ci>4)
				(*repeated)=(*repeated)*10+(cc-'0'); 
			ci++;
			
		}
		k++;

	}

	//printf("REP=%d\n", *repeated);

	if(field == 10 && ci>0)quality_string[ci]=0;
	else if(field < 10) return -1;
	
	if(ret_flag & 4)
		(*mapping_quality) = 0;
	else
		(*mapping_quality) = ret_quality;
	(*pos) = ret_pos;
	(*flags) = ret_flag;
	(*pair_dist) = ret_pairdist;
	//printf("FLAG=%d\n", (*flags));
	if(((*flags) & 4) == 4) return 1;
	return 0;
	
}


// This function returns 0 if the block is determined.
// The block is undeterminable if the chromosome name is not in known_chromosomes, or the position is larger than the known length.
// Pos is in terms of [1, ... , max_length]
int get_read_block(char *chro, unsigned int pos, char *temp_file_suffix, chromosome_t *known_chromosomes, unsigned int * max_base_position)
{
	int chro_no;
	unsigned int max_known_chromosome=0;

	for(chro_no=0;known_chromosomes[chro_no].chromosome_name[0]; chro_no++)
	{
		if(strcmp(chro , known_chromosomes[chro_no].chromosome_name) == 0)
		{
			max_known_chromosome = known_chromosomes[chro_no].known_length;
			break;
		}
		//if(chro_no > 1)
		//	printf("TOO MANY CHROS:%d\n", chro_no);
	}
	if(!known_chromosomes[chro_no].chromosome_name[0]) return 1;
	if(pos >= known_chromosomes[chro_no].known_length) return 1;

	int block_no = (pos-1) / BASE_BLOCK_LENGTH;
	sprintf(temp_file_suffix , "%s-%04u.bin", chro, block_no);
	if(max_base_position)*max_base_position=min((block_no+1)*BASE_BLOCK_LENGTH, max_known_chromosome);

	return 0;
}

FILE * get_temp_file_pointer(char *temp_file_name, HashTable* fp_table)
{
	FILE * temp_file_pointer = (FILE *) HashTableGet(fp_table, temp_file_name);
	if(!temp_file_pointer)
	{
		char *key_name;
		key_name = (char *)SUBREAD_malloc(300);
		if(!key_name)
			return NULL;
		strcpy(key_name, temp_file_name);
		temp_file_pointer = f_subr_open(key_name,"wb");

		if(0&&!temp_file_pointer)
		{
			struct rlimit limit_st;
			getrlimit(RLIMIT_NOFILE, & limit_st);
			if(limit_st.rlim_max>0 && limit_st.rlim_max <= 3000)
				limit_st.rlim_cur = min(limit_st.rlim_max, fp_table->numOfElements + 10);
			else
				limit_st.rlim_cur = max(limit_st.rlim_cur, fp_table->numOfElements + 10);
			setrlimit(RLIMIT_NOFILE, & limit_st);
			//if(rl==-1)
			//	SUBREADprintf("Cannot set limit: %d!\n", limit_st.rlim_cur);
			temp_file_pointer = f_subr_open(key_name,"wb");
		}


		if(!temp_file_pointer){
			SUBREADprintf("File cannot be opened: '%s' !!\nPlease increase the maximum open files by command 'ulimit -n'.\nThis number should be set to at least 500 for human genome, and more chromosomes require more opened files.\n\n", key_name);
			return NULL;
		}

		HashTablePut(fp_table, key_name ,temp_file_pointer);
	}

	return temp_file_pointer;
}

void my_fclose(void * fp)
{
	fclose((FILE *)fp);
}

int my_strcmp(const void * s1, const void * s2)
{
	return strcmp((char*)s1, (char*)s2);
}

void write_read_block_file(FILE *temp_fp , unsigned int read_number, char *read_name, int flags, char * chro, unsigned int pos, char *cigar, int mapping_quality, char *sequence , char *quality_string, int rl , int is_sequence_needed, char strand, unsigned short read_pos, unsigned short read_len)
{
	base_block_temp_read_t datum;
	datum.record_type = 100;
	datum.read_number = read_number;
	datum.pos = pos;
	datum.flags = flags;
	datum.strand = strand;
	datum.read_pos = read_pos;
	datum.read_len = read_len;
	datum.mapping_quality = mapping_quality;

	if(rl < 1|| rl > MAX_READ_LENGTH)
	{
		
		SUBREADprintf("READ IS TOO LONG:%d\n", rl);
		return;
	}

	fwrite(&datum, sizeof(datum), 1, temp_fp);
	if(is_sequence_needed)
	{
		unsigned short srl = rl&0xffff;
		fwrite(&srl, sizeof(short),1, temp_fp);
		fwrite(sequence , 1, rl,temp_fp );
		fwrite(quality_string , 1, rl,temp_fp );
	}
}


int get_known_chromosomes(char * in_SAM_file, chromosome_t * known_chromosomes)
{
	int i, is_first_read_PE;
	int is_BAM = is_certainly_bam_file(in_SAM_file,  &is_first_read_PE);
	SamBam_FILE * fp = SamBam_fopen(in_SAM_file,is_BAM?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);

	while(1)
	{
		char line_buffer [3000];
		char * is_ret = SamBam_fgets(fp, line_buffer, 2999, 0); 
		if(!is_ret) break;
		int linelen = strlen(line_buffer);

		if(line_buffer[0]=='@')
		{
			int chro_numb=0, field = 0, ci=0, ciw = 0;
			if(line_buffer[1]!='S' || line_buffer[2]!='Q' || line_buffer[3]!='\t' ) continue;

			while(known_chromosomes[chro_numb].chromosome_name[0]!=0) chro_numb++;
			if(chro_numb > XOFFSET_TABLE_SIZE-1)
			{
				SUBREADprintf("FATAL ERROR: the number of chromosomes excessed %d\n. Program terminated.\n", XOFFSET_TABLE_SIZE);
				return -1;
			}
			known_chromosomes[chro_numb].known_length = 0;
			for(i=0; i< linelen; i++)
			{
				char cc = line_buffer[i];

				if(cc == '\r' || cc=='\n') continue;

				if(cc == '\t')
				{
					if(field == 1)
						known_chromosomes[chro_numb].chromosome_name[ciw]=0;
					ci = 0;
					ciw = 0;
					field ++;
				}
				else if(field == 1)
				{
					if(ci >2)
						known_chromosomes[chro_numb].chromosome_name[ciw++]=cc;
					ci++;
				}
				else if(field == 2)
				{
					if(ci >2)
						known_chromosomes[chro_numb].known_length = known_chromosomes[chro_numb].known_length * 10 + (cc - '0');
					ci++;
				}
			}
		}
		else
			break;
	}
	SamBam_fclose(fp);
	return 0;
}

void add_cigar_indel_event(HashTable * event_table_ptr, char * chro, unsigned int chro_pos, int indels , char * ins_seq)
{
	if(abs(indels)>100) return;

	char event_token[100];
	snprintf(event_token, 99,"%s\t%u", chro, chro_pos);
	int x1;
	unsigned int indel_event_id = 0xffffffff, token_len;

	int exist_indel_count = HashTableGet(event_table_ptr, event_token) - NULL;
	unsigned short * app2_ptr =  event_table_ptr->appendix2;

	if(exist_indel_count)
		for(x1 = 0; x1< exist_indel_count; x1++)
		{
			snprintf(event_token, 99,"%s\t%u\t%d", chro, chro_pos, x1);
			long long int t64v =  (HashTableGet(event_table_ptr, event_token)-NULL);
			long long int indel_len = (t64v&0xff) - 0x80;
			if(indel_len == indels){
				indel_event_id = 0xffffff&(t64v >> 8) ;
				if(app2_ptr[indel_event_id]<65000)
					app2_ptr[indel_event_id] +=1;
				return;
			}
		}


	if(event_table_ptr->counter2<0xffff00)
	{
		unsigned int event_space_max_size = event_table_ptr-> counter1;
		indel_event_id = event_table_ptr->counter2 ++;

		if(indel_event_id >= event_space_max_size)
		{
			event_table_ptr->appendix1 = realloc(event_table_ptr->appendix1 , sizeof(char *) * event_space_max_size*2);
			event_table_ptr->appendix2 = realloc(event_table_ptr->appendix2 , sizeof(short) * event_space_max_size*2);
			memset(event_table_ptr->appendix2 + event_space_max_size * sizeof(short), 0, sizeof(short) * event_space_max_size);
			event_table_ptr-> counter1 = event_space_max_size*2;
			app2_ptr =  event_table_ptr->appendix2;
		}

		token_len=snprintf(event_token, 99,"%s\t%u", chro, chro_pos);
		if(exist_indel_count<1)
		{
			char * token_1 = malloc(token_len+1);
			strcpy(token_1, event_token);
			HashTablePut(event_table_ptr, token_1, NULL+1);
		}
		else
		{
			HashTablePutReplace(event_table_ptr, event_token, NULL+exist_indel_count+1, 0);
		}

		token_len=snprintf(event_token, 99,"%s\t%u\t%d", chro, chro_pos, exist_indel_count);
		char * token_2 = malloc(token_len+1);
		strcpy(token_2, event_token);
		long long int indel_event_id_long = indel_event_id;
		app2_ptr[indel_event_id] +=1;

		HashTablePut(event_table_ptr, token_2, NULL + ((0xff & (0x80 + indels)) | ((indel_event_id_long&0xffffff) << 8)));
		if(indels<0)
		{
			char * ins_seq_2 = malloc(-indels), ** app1_ptrptr = event_table_ptr->appendix1;
			memcpy(ins_seq_2, ins_seq, -indels);
			app1_ptrptr[indel_event_id] = ins_seq_2;
		}
	}
}

void destroy_cigar_event_table(HashTable * event_table)
{
	int bucket;
	KeyValuePair * cursor;
	char ** seq_tab = event_table->appendix1;
	for(bucket=0; bucket<event_table -> numOfBuckets; bucket++)
	{
		cursor = event_table -> bucketArray[bucket];
		while (1)
		{
			int xk1, tabs;
			if (!cursor) break;
			
			char * token = (char *)cursor -> key;
			tabs = 0;
			for(xk1=0; token[xk1]; xk1++) 
				if(token[xk1]=='\t') tabs++;
			long long int tmpv = cursor -> value - NULL;
			//printf("%s\t%lld\n", token, tmpv);

			if(tabs==3)
			{
				unsigned int event_id = (tmpv>>8)&0xffffff;
				free(seq_tab[event_id]);
			}
			free(token);
			cursor = cursor->next;
		}
	}

	free(event_table->appendix1);
	free(event_table->appendix2);
	HashTableDestroy(event_table);
}

void break_VCF_file(char * vcf_file, HashTable * fp_table, char * temp_file_prefix, chromosome_t* known_chromosomes)
{
	FILE * vfp=fopen(vcf_file, "r");
	char temp_file_suffix[MAX_CHROMOSOME_NAME_LEN+20];

	if(!vfp)
	{
		SUBREADprintf("The specified VCF does not exist.\n");
		return;
	}

	char * linebuf = malloc(3000);
	char * tmpfname = malloc(400);

	while(1)
	{
		char * tok_tmp;
		char * retstr = fgets(linebuf, 2999, vfp);
		if(!retstr) break;
		if(linebuf[0]=='#') continue;

		char * chro = strtok_r(linebuf, "\t", &tok_tmp);
		if(!tok_tmp) continue;
		char * pos_str = strtok_r(NULL, "\t", &tok_tmp);
		if(!tok_tmp) continue;

		strtok_r(NULL, "\t", &tok_tmp);// name
		if(!tok_tmp) continue;

		char * ref_seq = strtok_r(NULL, "\t", &tok_tmp);
		if(!tok_tmp) continue;
		char * alt_seq = strtok_r(NULL, "\t", &tok_tmp);
		if(!tok_tmp) continue;

		int is_snp = 0;
		if(strstr(alt_seq,","))
		{
			char * com_tmp = NULL;
			char * com_sec = strtok_r(alt_seq, ",", &com_tmp);
			while(com_sec)
			{
				if(strlen(com_sec)==strlen(ref_seq))
				{
					is_snp=1;
					break;
				}

				com_sec = strtok_r(NULL,  ",", &com_tmp);
			}

		}else if(strlen(ref_seq) == strlen(alt_seq)) is_snp=1;

		if(!is_snp)continue;
		unsigned int max_section_pos;

		if(get_read_block(chro, atoi(pos_str) , temp_file_suffix, known_chromosomes, &max_section_pos))continue;
		sprintf(tmpfname, "%s%s", temp_file_prefix , temp_file_suffix);
		FILE * temp_fp = get_temp_file_pointer(tmpfname, fp_table);
		if(temp_fp)
		{
			VCF_temp_read_t datum;
			datum.record_type = 200;
			datum.pos = atoi(pos_str);
			datum.type = CHRO_EVENT_TYPE_SNP;
			fwrite(&datum, sizeof(VCF_temp_read_t), 1, temp_fp);
		}
	}

	free(linebuf);
	free(tmpfname);
	fclose(vfp);
}

int break_SAM_file(char * in_SAM_file, int is_BAM_file, char * temp_file_prefix, unsigned int * real_read_count, int * block_count, chromosome_t * known_chromosomes, int is_sequence_needed, int base_ignored_head_tail, gene_value_index_t *array_index, gene_offset_t * offsets, unsigned long long int * all_mapped_bases, HashTable * event_table, char * VCF_file)
{
	int i, is_first_read=1;
	HashTable * fp_table;
	unsigned int read_number = 0;
	char line_buffer [3000];
	SamBam_FILE  * sambam_reader;

	sambam_reader = SamBam_fopen(in_SAM_file, is_BAM_file?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);
	 
	if(!sambam_reader){
		SUBREADprintf("SAM file does not exist or is not accessible: '%s'\n", in_SAM_file);
		return 1;
	}

	fp_table = HashTableCreate( 11011 );
	HashTableSetDeallocationFunctions(fp_table, free, my_fclose);
	HashTableSetKeyComparisonFunction(fp_table, my_strcmp);
	HashTableSetHashFunction(fp_table,HashTableStringHashFunction);

	if(event_table!=NULL && event_table->appendix1==NULL)
	{
		event_table->appendix1 = malloc(sizeof(char *) * 100);
		event_table->appendix2 = malloc(sizeof(unsigned short) * 100);
		memset(event_table->appendix2, 0, sizeof(unsigned short) * 100);
		event_table->counter1 = 100;
		event_table->counter2 = 0;
	}

	while(1)
	{
		//unsigned long long int file_position = ftello(fp);
		//int linelen = read_line(2999, fp, line_buffer, 0);
		char * is_ret = SamBam_fgets(sambam_reader, line_buffer, 2999, 1);

		if(!is_ret) break;

		if(line_buffer[0]=='@')
		{
			int chro_numb=0, field = 0, ci=0, ciw = 0;
			if(line_buffer[1]!='S' || line_buffer[2]!='Q' || line_buffer[3]!='\t' ) continue;

			while(known_chromosomes[chro_numb].chromosome_name[0]!=0) chro_numb++;

			if(chro_numb > XOFFSET_TABLE_SIZE-1)
			{
				SUBREADprintf("FATAL ERROR: the number of chromosomes excessed %d\n. Program terminated.\n", XOFFSET_TABLE_SIZE);
				return -1;
			}

			known_chromosomes[chro_numb].known_length = 0;
			for(i=0; ; i++)
			{
				char cc = line_buffer[i];
				if(!cc) break;

				if(cc == '\r' || cc=='\n') continue;

				if(cc == '\t')
				{
					if(field == 1)
						known_chromosomes[chro_numb].chromosome_name[ciw]=0;
					ci = 0;
					ciw = 0;
					field ++;
				}
				else if(field == 1)
				{
					if(ci >2)
						known_chromosomes[chro_numb].chromosome_name[ciw++]=cc;
					ci++;
				}
				else if(field == 2)
				{
					if(ci >2)
						known_chromosomes[chro_numb].known_length = known_chromosomes[chro_numb].known_length * 10 + (cc - '0');
					ci++;
				}
			}
			if(chro_numb < XOFFSET_TABLE_SIZE-1) known_chromosomes[chro_numb+1].chromosome_name[0]=0;
		}
		else
		{
			char read_name[MAX_READ_NAME_LEN], chro[MAX_CHROMOSOME_NAME_LEN], cigar[EXON_MAX_CIGAR_LEN], sequence[MAX_READ_LENGTH+1], quality_string[MAX_READ_LENGTH+1];
			int flags = 0, mapping_quality = 0, rl=0;
			char is_negative_strand = 0;
			unsigned int pos = 0, pairdist = 0;
			char temp_file_suffix[MAX_CHROMOSOME_NAME_LEN+20];
			char temp_file_name[MAX_CHROMOSOME_NAME_LEN+20+300];
			FILE * temp_fp;
			int repeated = -1;

			if(is_first_read)
			{
				is_first_read=0;

				if(VCF_file && VCF_file[0])
					break_VCF_file(VCF_file, fp_table, temp_file_prefix, known_chromosomes);
			}


			//SUBREADprintf("ARRI_0=%p ; OFFS=%p ; EVT=%p\n%s\n",array_index, offsets, event_table, line_buffer);
			int line_parse_result = parse_SAM_line(line_buffer, read_name, &flags, chro, &pos, cigar, & mapping_quality, &pairdist, sequence , quality_string, &rl, &repeated);
			//SUBREADprintf("ARRI_2=%p ; OFFS=%p ; EVT=%p\n",array_index, offsets, event_table);

			if(strlen(quality_string)<2)
			{
				int xk1;
				for(xk1=0; xk1<rl; xk1++)
				{
					quality_string[xk1]='I';
				}
				quality_string[xk1]=0;
			}

			if(line_parse_result || (flags & SAM_FLAG_UNMAPPED)|| (((flags & SAM_FLAG_PAIRED_TASK) && (pairdist ==0 || pairdist > 500000))))
			{

				read_number ++;
				continue;
			}

			if(array_index)
			{
				int mismatch = 0;

				unsigned int linear_pos = linear_gene_position(offsets , chro, pos)-1;
				float match_rate = final_mapping_quality(array_index, linear_pos, sequence, quality_string, cigar, FASTQ_PHRED33,  & mismatch,  rl, NULL, NULL);
				if(mismatch>8 || match_rate < 160)
				{
					read_number ++;
					continue;
				}
			}

			is_negative_strand = (flags & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0;

			if(is_sequence_needed == 2)
			{
				
			}
			else if(is_sequence_needed == 1)
			{
				int read_cursor = 0;
				int is_first_S = 1;
				unsigned int chromosome_cursor = pos;
				int j, tmpv=0;
				char cc;

				for(j=0; cigar[j]; j++)
				{
					cc = cigar[j];
					if(cc>='0' && cc<='9') tmpv= tmpv*10+(cc-'0');
					else if(cc == 'S'||cc == 'M')
					{
						if(cc == 'M') is_first_S = 0;

						if(cc == 'M')
						{
							unsigned int insertion_cursor = chromosome_cursor;
							// DO INSERTION
							while(insertion_cursor < (chromosome_cursor + tmpv) && read_cursor < (rl - base_ignored_head_tail))
							{
								unsigned int max_section_pos, insert_length;
								int need_write = 1;

								if(get_read_block(chro, insertion_cursor , temp_file_suffix, known_chromosomes, &max_section_pos))break;
								insert_length = min(max_section_pos + 1, chromosome_cursor + tmpv) - insertion_cursor;
								if(insert_length<1) break;

								if(base_ignored_head_tail)
								{
									if(read_cursor+insert_length < base_ignored_head_tail)
										need_write = 0;
									else if(read_cursor < base_ignored_head_tail)
									{
										int ignored_length = base_ignored_head_tail - read_cursor;
										insert_length = read_cursor + insert_length - base_ignored_head_tail;
										
										read_cursor = base_ignored_head_tail;
										insertion_cursor += ignored_length;
									}

									if(read_cursor >= (rl - base_ignored_head_tail))
										need_write = 0;
									else if(read_cursor +insert_length >= (rl - base_ignored_head_tail))
										insert_length = (rl - base_ignored_head_tail) - read_cursor;
								}
//								printf("INST: RL=%d; INSL=%d; READ_CUR=%d; IGNORE=%d\n", rl, insert_length, read_cursor , base_ignored_head_tail);

								if(need_write && insert_length >= 5)
								{
									sprintf(temp_file_name, "%s%s", temp_file_prefix , temp_file_suffix);
									temp_fp = get_temp_file_pointer(temp_file_name, fp_table);
									if(!temp_fp) return -1;
									if(all_mapped_bases)
										(*all_mapped_bases) += insert_length;
									write_read_block_file(temp_fp , read_number, read_name, flags, chro, insertion_cursor, cigar, mapping_quality, sequence + read_cursor , quality_string + read_cursor, insert_length , 1, is_negative_strand, read_cursor, rl);
								}
								insertion_cursor += insert_length;
								read_cursor += insert_length;
							}
						}
						else 
							read_cursor += tmpv;

						if(!is_first_S)
							chromosome_cursor += tmpv;

						tmpv=0;
					}
					else if(cc == 'D' || cc == 'N')
					{
						// the left edge ( last WANTED base ) is chromosome_cursor-1
						// the indel length is tmpv;
						// now we add this into the event table.
						if(event_table && cc=='D')
							add_cigar_indel_event(event_table, chro, chromosome_cursor-1, tmpv, NULL);
						chromosome_cursor += tmpv;
						tmpv = 0;
					}
					else if(cc == 'I' )
					{
						// the left edge ( last WANTED base ) is chromosome_cursor-1
						// the indel length is -tmpv;
						// now we add this into the event table.
						if(event_table)
							add_cigar_indel_event(event_table, chro, chromosome_cursor-1, -tmpv, sequence + read_cursor);
						read_cursor += tmpv;
						tmpv = 0;
					}
					else	tmpv = 0;

				}
				
			}else
			{
				if(get_read_block(chro, pos, temp_file_suffix, known_chromosomes, NULL))
				{
					read_number ++;
					continue;
				}
				sprintf(temp_file_name, "%s%s", temp_file_prefix , temp_file_suffix);
	
				temp_fp = get_temp_file_pointer(temp_file_name, fp_table);
				assert(temp_fp);

				write_read_block_file(temp_fp , read_number, read_name, flags, chro, pos, cigar, mapping_quality, sequence , quality_string, rl , is_sequence_needed, is_negative_strand, 0,rl);
			}
			read_number ++;
		}
	}

	if(block_count)
		(*block_count) = fp_table->numOfElements;
	HashTableDestroy(fp_table);
	SamBam_fclose(sambam_reader);
	if(real_read_count)
		(*real_read_count) = read_number;
	return 0;
}

int is_in_exon_annotations(gene_t *output_genes, unsigned int offset, int is_start)
{
	int i,j;

	for(i=0; i< MAX_ANNOTATION_EXONS; i++)
	{
		if(!output_genes[i].end_offset) break;
		if(output_genes[i].end_offset >= offset && output_genes[i].start_offset <= offset)
		{
			for(j=0; j< MAX_EXONS_PER_GENE; j++)
			{
				if(output_genes[i].exon_ends[j] >= offset && output_genes[i].exon_starts[j] <= offset)
				{
					if(output_genes[i].exon_starts[j] == offset && is_start) return 2;	// 2==exactly matched
					if(output_genes[i].exon_ends[j] == offset && !is_start)	return 2;
					return 1;	// 1==enclosed
				}
			}
		}
	}
	return 0;	//0==exon not found
}

int load_exon_annotation(char * annotation_file_name, gene_t ** output_genes, gene_offset_t* offsets)
{
	int line_len, gene_number = 0, exons = 0;
	char old_gene_name[MAX_GENE_NAME_LEN];
	FILE * fp = f_subr_open(annotation_file_name, "rb");

	if(!fp)
	{
		SUBREADprintf("Cannot open the exon annotation file: %s\n", annotation_file_name);
		return -1;
	}
	(*output_genes) = malloc(sizeof(gene_t)*MAX_ANNOTATION_EXONS);
	if(!*output_genes)
	{
		SUBREADprintf("Cannot allocate memory for the exon table. \n");
		return -1;
	}

	
	old_gene_name[0]=0;
	(*output_genes)[0].end_offset = 0;
	(*output_genes)[0].start_offset = 0xffffffff;
	while(gene_number < MAX_ANNOTATION_EXONS)
	{
		char buff[200], this_gene_name[MAX_GENE_NAME_LEN], chromosome_name[MAX_CHROMOSOME_NAME_LEN];
		int i = 0, j=0;
		unsigned int exon_location;

		line_len = read_line(200, fp, buff, 0);	

		if(line_len>0)	//Not EOF
		{
			if(!isdigit(buff[0]))	// it is a title line or something else
				continue;
		
			for(i=0; buff[i] != '\t' &&  buff[i] != '\n' && i < 200; i++)
				this_gene_name[i] = buff[i];
			this_gene_name[i] = 0;
		}
		
		if(line_len<=0 || (exons && old_gene_name[0] && strcmp(this_gene_name , old_gene_name)))	// it is a new gene
		{
			strncpy((*output_genes)[gene_number].gene_name , old_gene_name, MAX_GENE_NAME_LEN);
			(*output_genes)[gene_number].exon_ends[exons] = 0;
			gene_number++;
			exons = 0;
			(*output_genes)[gene_number].end_offset = 0;
			(*output_genes)[gene_number].start_offset = 0xffffffff;
		}

		if(line_len<=0) break;

	
		// copy chromosome name
		for(i++; buff[i] != '\t' &&  buff[i] != '\n' && i < 200; i++)
			chromosome_name[j++] = buff[i];
		chromosome_name[j] = 0;

		// start location
		exon_location = 0;
		for(i++; buff[i] != '\t' &&  buff[i] != '\n' && i < 200; i++)
			if(isdigit(buff[i]))
				exon_location = exon_location*10 + buff[i] - '0';

		(*output_genes)[gene_number].exon_starts[exons] = linear_gene_position(offsets, chromosome_name , exon_location-1); 
		if( (*output_genes)[gene_number].exon_starts[exons] == 0xffffffff)
			continue;

		if((*output_genes)[gene_number].start_offset > (*output_genes)[gene_number].exon_starts[exons])
			(*output_genes)[gene_number].start_offset = (*output_genes)[gene_number].exon_starts[exons];

		// end location
		exon_location = 0;
		for(i++; buff[i] != '\t' &&  buff[i] != '\n' && buff[i] && i < 200; i++)
			if(isdigit(buff[i]))
				exon_location = exon_location*10 + buff[i] - '0';

		(*output_genes)[gene_number].exon_ends[exons] = linear_gene_position(offsets, chromosome_name , exon_location); 

		if((*output_genes)[gene_number].end_offset <  (*output_genes)[gene_number].exon_ends[exons])
			(*output_genes)[gene_number].end_offset =  (*output_genes)[gene_number].exon_ends[exons];

		exons ++;
		if(exons >= MAX_EXONS_PER_GENE)
		{
			SUBREADprintf("The number of exons excesses the limit. Please increase the value of MAX_EXONS_PER_GENE in subread.h.\n");
			return -1;
		}

		strncpy(old_gene_name, this_gene_name , MAX_GENE_NAME_LEN);
	}
	fclose(fp);
	return 0;
}

int does_file_exist(char * path)
{
	int ret ;
	FILE * fp = f_subr_open(path, "rb");
	ret = fp!=NULL;
	if(fp)fclose(fp);

	return ret;
}

unsigned long long int sort_SAM_hash(char * str)
{
	unsigned long long int hash = 5381;
	int c, xk1=0;

	while (1)
	{
		c = str[xk1++];
		if(!c)break;
		hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
	}
	return hash;
}


void do_SIGINT_remove(char * prefix, int param);
char * _SAMSORT_SNP_delete_temp_prefix = NULL;
char * _REPAIRER_delete_temp_prefix = NULL;
void SAM_SORT_SIGINT_hook(int param) {
	do_SIGINT_remove(_SAMSORT_SNP_delete_temp_prefix,  param);
}
void REPAIR_SIGINT_hook(int param) {
	do_SIGINT_remove(_REPAIRER_delete_temp_prefix,  param);
}

void delete_with_prefix(char * prefix){
	if(prefix != NULL)
	{
		int xk1, last_slash = -1;
		char del2[300], del_suffix[200], del_name[400];
		for(xk1=0; prefix[xk1]; xk1++)
		{
			if(prefix[xk1]=='/') last_slash = xk1;
			else if(prefix[xk1]=='\\')
			{
				SUBREADprintf("The file name is unknown.\n");
				return;
			}
		}
		if(last_slash>=0)
		{
			memcpy(del2, prefix, last_slash);
			del2[last_slash] = 0;
			strcpy(del_suffix , prefix + last_slash + 1);
		}
		else
		{
			strcpy(del2,".");
			strcpy(del_suffix , prefix);
		}
	
		if(strlen(del_suffix)>8)
		{
			DIR           *d;
			struct dirent *dir;

			d = opendir(del2);
			if (d)
			{
				while ((dir = readdir(d)) != NULL)
				{
					if(strstr(dir->d_name, del_suffix))
					{
						//SUBREADprintf("DEL:%s\n", dir->d_name);
						strcpy(del_name, del2);
						strcat(del_name, "/");
						strcat(del_name, dir->d_name);
						unlink(del_name);
						//test fix
					}
				}
			}
		}
			
	}

}

void do_SIGINT_remove(char * prefix, int param) {
	#ifdef MAKE_STANDALONE
	delete_with_prefix(prefix);
	SUBREADprintf("\n\nReceived a terminal signal. The temporary files were removed.\n");
	exit(param);
	#endif
}


void * old_sig_TERM = NULL, * old_sig_INT = NULL;

#define PAIRER_GZIP_WINDOW_BITS -15
#define PAIRER_DEFAULT_MEM_LEVEL 8

int SAM_pairer_writer_create( SAM_pairer_writer_main_t * bam_main , int all_threads , int has_dummy, int BAM_input, int c_level, char * out_file){
	int x1;

	memset(bam_main, 0, sizeof(SAM_pairer_writer_main_t));
	bam_main -> bam_fp = f_subr_open(out_file, "wb");
	if(NULL == bam_main -> bam_fp) return 1;
	strcpy(bam_main -> bam_name, out_file);
	bam_main -> threads = malloc(all_threads * sizeof(SAM_pairer_writer_thread_t));
	bam_main -> all_threads = all_threads;
	bam_main -> has_dummy = has_dummy;
	bam_main -> compression_level = c_level;
	subread_init_lock(&bam_main -> output_fp_lock);

	for(x1 = 0; x1 < all_threads ; x1 ++){
		bam_main -> threads[x1].BIN_buffer_ptr = 0;
		bam_main -> threads[x1].strm.zalloc = Z_NULL;
		bam_main -> threads[x1].strm.zfree = Z_NULL;
		bam_main -> threads[x1].strm.opaque = Z_NULL;
		bam_main -> threads[x1].strm.avail_in = 0;
		bam_main -> threads[x1].strm.next_in = Z_NULL;

		deflateInit2(&bam_main -> threads[x1].strm, bam_main -> compression_level, Z_DEFLATED,
                PAIRER_GZIP_WINDOW_BITS, PAIRER_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
	}
	return 0;
}

void SAM_pairer_write_BAM_header(FILE * writer, int compressed_size)
{

	// the four magic characters
	fputc(31,  writer);
	fputc(139,  writer);
	fputc(8,  writer);
	fputc(4,  writer);

	time_t time_now = 0;
	fwrite(&time_now,4,1, writer);

	int tmp_i;
	// Extra flags and OS
	fputc(0,  writer);
	fputc(0xff,  writer); 

	// Extra length
	tmp_i = 6;
	fwrite(&tmp_i,2,1, writer);


	// SI1 and SI2 magic numbers, and SLEN
	fputc(66,  writer);
	fputc(67,  writer);
	tmp_i = 2;
	fwrite(&tmp_i,2,1, writer);
	tmp_i = compressed_size + 19 + 6;
	fwrite(&tmp_i,2,1, writer);
}



int SAM_pairer_multi_thread_compress(SAM_pairer_writer_main_t * bam_main ,  SAM_pairer_writer_thread_t * bam_thread)
{
	#define BAM_compressed_space 65536
	char * BAM_compressed = malloc(BAM_compressed_space);
	int ret, have;
	if(bam_thread -> BIN_buffer_ptr>0){
		deflateReset(&bam_thread -> strm);
		bam_thread -> strm.avail_in = bam_thread -> BIN_buffer_ptr;
		bam_thread -> strm.next_in = bam_thread -> BIN_buffer;
		bam_thread -> strm.avail_out = BAM_compressed_space;
		bam_thread -> strm.next_out = (unsigned char *)BAM_compressed;
		ret = deflate( &bam_thread -> strm , Z_FINISH);

		have = BAM_compressed_space - bam_thread -> strm.avail_out;
		assert(bam_thread -> strm.avail_in == 0);
	}else{
		z_stream nstrm;
		nstrm.zalloc = Z_NULL;
		nstrm.zfree = Z_NULL;
		nstrm.opaque = Z_NULL;
		nstrm.avail_in = 0;
		nstrm.next_in = Z_NULL;
	
		deflateInit2(&nstrm, SAMBAM_COMPRESS_LEVEL, Z_DEFLATED,
			PAIRER_GZIP_WINDOW_BITS, PAIRER_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);

		nstrm.avail_in = 0;
		nstrm.next_in = bam_thread -> BIN_buffer;
		nstrm.avail_out = BAM_compressed_space;
		nstrm.next_out = (unsigned char *)BAM_compressed;
		ret = deflate(&nstrm, Z_FINISH);
		deflateEnd(&nstrm);
		have = BAM_compressed_space - nstrm.avail_out;
	}
	if(ret == Z_OK || 1){

		//SUBREADprintf("Compress: %d -> %d  %p\n", bam_thread -> BIN_buffer_ptr, have, bam_main -> bam_fp);
		//if(bam_thread -> BIN_buffer_ptr == 0) have = 0;
		unsigned int crc0 = crc32(0, NULL, 0);
		unsigned int CRC32 = crc32(crc0, (unsigned char *)  bam_thread -> BIN_buffer ,bam_thread -> BIN_buffer_ptr);


		subread_lock_occupy( &bam_main -> output_fp_lock );
		SAM_pairer_write_BAM_header( bam_main -> bam_fp , have);
		fwrite(BAM_compressed,1, have, bam_main -> bam_fp );
		fwrite(&CRC32 , 4, 1, bam_main -> bam_fp);
		fwrite( &bam_thread -> BIN_buffer_ptr , 4, 1, bam_main -> bam_fp);

		subread_lock_release( &bam_main -> output_fp_lock );

		bam_thread -> BIN_buffer_ptr = 0;
	} else {
		SUBREADprintf("ERROR: Cannot compress a BAM block : %d\n", ret);
		return 1;
	}
	free(BAM_compressed);
	return 0;
}



void SAM_pairer_writer_destroy( SAM_pairer_writer_main_t * bam_main ) {
	int x1;
	for(x1 = 0; x1 < bam_main -> all_threads ; x1 ++){
		if(bam_main -> threads[x1].BIN_buffer_ptr>0){
			SAM_pairer_multi_thread_compress(bam_main, bam_main->threads+x1);
		}

		if(x1 == bam_main -> all_threads - 1){
			assert(0 == bam_main -> threads[x1].BIN_buffer_ptr);
			SAM_pairer_multi_thread_compress(bam_main, bam_main->threads+x1);
		}
		deflateEnd(&bam_main -> threads[x1].strm);
	}
	subread_destroy_lock(&bam_main -> output_fp_lock);
	fclose(bam_main -> bam_fp);
	free(bam_main -> threads);
}

// Tiny_Mode only write the following information:
// Name   Flag   Chro   Pos   Mapq   Cigar   MateChro   MatePos   Tlen  N  I  NH:i:xx  HI:i:xx
// Tiny_Mode does not work when output and input are both in BAM format
// in_format can be either 
// bin_buff_size_per_thread is in Mega-Bytes.
// It returns 0 if no error
int SAM_pairer_create(SAM_pairer_context_t * pairer, int all_threads, int bin_buff_size_per_thread, int BAM_input, int is_Tiny_Mode, int display_progress, char * in_file, void (* reset_output_function) (void * pairer), int (* output_header_function) (void * pairer, int thread_no, int is_text, unsigned int items, char * bin, unsigned int bin_len), int (* output_function) (void * pairer, int thread_no, char * readname, char * bin1, char * bin2), char * tmp_path, void * appendix1) {

	memset(pairer, 0, sizeof(SAM_pairer_context_t));
	pairer -> input_fp = f_subr_open(in_file, "rb");
	if(NULL == pairer -> input_fp) return 1;

	pairer -> input_is_BAM = BAM_input;
	pairer -> tiny_mode = is_Tiny_Mode;
	pairer -> reset_output_function = reset_output_function;
	pairer -> output_function = output_function;
	pairer -> output_header = output_header_function;
	pairer -> display_progress = display_progress;
	subread_init_lock(&pairer -> input_fp_lock);

	pairer -> total_threads = all_threads;
	pairer -> input_buff_SBAM_size = bin_buff_size_per_thread * 1024 * 1024;
	pairer -> input_buff_BIN_size = 1024*1024;

	pairer -> appendix1 = appendix1;

	_REPAIRER_delete_temp_prefix = tmp_path;
	old_sig_TERM = signal (SIGTERM, REPAIR_SIGINT_hook);
	old_sig_INT = signal (SIGINT, REPAIR_SIGINT_hook);
	
	strcpy(pairer -> tmp_file_prefix, tmp_path);
	pairer -> threads = malloc(all_threads * sizeof(SAM_pairer_thread_t));
	memset(pairer -> threads, 0, all_threads * sizeof(SAM_pairer_thread_t));

	int x1;

	for(x1 = 0; x1 < all_threads ; x1++){
		pairer -> threads[x1].thread_id = x1;
		pairer -> threads[x1].input_buff_SBAM = malloc(pairer -> input_buff_SBAM_size);
		pairer -> threads[x1].input_buff_BIN = malloc(pairer -> input_buff_BIN_size);
		pairer -> threads[x1].input_buff_BIN_used = 0;
		pairer -> threads[x1].orphant_table = HashTableCreate(pairer -> input_buff_SBAM_size / 100);
		HashTableSetHashFunction(pairer -> threads[x1].orphant_table, fc_chro_hash);
		HashTableSetKeyComparisonFunction(pairer -> threads[x1].orphant_table, fc_strcmp_chro);
		HashTableSetDeallocationFunctions(pairer -> threads[x1].orphant_table, free, free);
		pairer -> threads[x1].strm.zalloc = Z_NULL;
		pairer -> threads[x1].strm.zfree = Z_NULL;
		pairer -> threads[x1].strm.opaque = Z_NULL;
		pairer -> threads[x1].strm.avail_in = 0;
		pairer -> threads[x1].strm.next_in = Z_NULL;
		
		inflateInit2(&pairer -> threads[x1].strm, PAIRER_GZIP_WINDOW_BITS);
	}
	return 0;
}

void SAM_pairer_destroy(SAM_pairer_context_t * pairer){

	int x1;
	unsigned long long all_orphants = 0;
	for(x1 = 0; x1 < pairer -> total_threads ; x1++){
		inflateEnd(&pairer -> threads[x1].strm);
		free(pairer -> threads[x1].input_buff_SBAM);
		free(pairer -> threads[x1].input_buff_BIN);
		
		all_orphants += pairer -> threads[x1].orphant_table->numOfElements;
		HashTableDestroy(pairer -> threads[x1].orphant_table);
	}

	subread_destroy_lock(&pairer -> input_fp_lock);
	delete_with_prefix(pairer -> tmp_file_prefix);
	fclose(pairer -> input_fp);
	free(pairer -> threads);
	signal (SIGTERM, old_sig_TERM);
	signal (SIGINT, old_sig_INT);
	//SUBREADprintf("All orphans=%llu frags\n", all_orphants);
}

// always assume that fp is at the start of a BAM GZ block.
int SAM_pairer_read_BAM_block(FILE * fp, int max_read_len, char * inbuff) {
	unsigned char gz_header_12 [12];
	int read_len = fread(gz_header_12, 1, 12, fp );
	if(read_len < 12) return -1;
	if(gz_header_12[0]!=31 || gz_header_12[1]!=139){
		SUBREADprintf("BAD GZ BAM: %u, %u\n", gz_header_12[0], gz_header_12[1]);
		return -1;
	}
	unsigned short xlen = 0, bsize = 0;
	memcpy(&xlen, gz_header_12 + 10, 2);
	int xlen_read = 0;

	while( xlen_read < xlen ){
		unsigned char x_header_4[4];
		unsigned short slen = 0;
		read_len = fread(x_header_4, 1, 4, fp);
		if(read_len < 4){
			SUBREADprintf("BAD GZ BAM 6LEN\n");
			return -1;
		}
		memcpy(&slen, x_header_4+2 , 2);
		xlen_read += 4;
		if(x_header_4[0]==66 && x_header_4[1]==67 && slen == 2){
			read_len = fread(&bsize, 2, 1, fp);
			if(read_len < 1){
				SUBREADprintf("BAD GZ BAM XLEN\n");
				return -1;
			}
		}else{
			fseek(fp, slen, SEEK_CUR);
		}
		xlen_read += slen;
	}
	if(bsize < 1 || bsize < xlen + 19){
		SUBREADprintf("BAD GZ BAM BSIZE\n");
		return -1;
	}
	read_len = fread(inbuff, 1, bsize - xlen - 19, fp);
	//SUBREADprintf("GOOD GZ\n");

	// seek over CRC and ISIZE
	fseek(fp, 8, SEEK_CUR);
	if(read_len < bsize - xlen - 19) return -1;
	return read_len;
}

#define MIN_BAM_BLOCK_SIZE 66000

int SAM_pairer_read_SAM_line( FILE * fp, int max_read_len, char * inbuff ){
	return read_line(max_read_len, fp, inbuff, 0);
}

void SAM_pairer_fill_BIN_buff(SAM_pairer_context_t * pairer ,  SAM_pairer_thread_t * thread_contaxt , int * is_finished){
	// load continuous K blocks and decompress.
	int current_buffer_used = 0;
	int current_blocks = 0;
	while(1){
		if(feof(pairer -> input_fp)){
			*is_finished = 1;
			break;
		}
		if(pairer -> input_buff_SBAM_size - current_buffer_used < MIN_BAM_BLOCK_SIZE) {
			break;
		}
		int this_size = 0;

		if(pairer -> input_is_BAM){
			this_size = SAM_pairer_read_BAM_block( pairer -> input_fp , pairer -> input_buff_SBAM_size - current_buffer_used , thread_contaxt -> input_buff_SBAM + current_buffer_used);
		}else{
			this_size = SAM_pairer_read_SAM_line( pairer -> input_fp , pairer -> input_buff_SBAM_size - current_buffer_used , thread_contaxt -> input_buff_SBAM + current_buffer_used);
		}

		current_blocks ++;
		if(this_size >= 0) {
			current_buffer_used += this_size;
		} else {
			*is_finished = 1;
			break;
		}
	}

	//SUBREADprintf("PAPA:READ=%d by %d blocks  %p\n", current_buffer_used, current_blocks, thread_contaxt);
	thread_contaxt -> input_buff_SBAM_used = current_buffer_used;
	thread_contaxt -> input_buff_SBAM_ptr = 0;
}

#define BAM_next_nch { \
	while(thread_contaxt -> input_buff_BIN_ptr >= thread_contaxt -> input_buff_BIN_used){int retXX = SAM_pairer_fetch_BAM_block(pairer, thread_contaxt);  if(retXX) break;}\
	nch = thread_contaxt -> input_buff_BIN[thread_contaxt -> input_buff_BIN_ptr++];}

int SAM_pairer_fetch_BAM_block(SAM_pairer_context_t * pairer , SAM_pairer_thread_t * thread_contaxt){
	if(thread_contaxt -> input_buff_SBAM_used <=  thread_contaxt -> input_buff_SBAM_ptr){
		return 1;
	}

	if(thread_contaxt -> input_buff_BIN_ptr < thread_contaxt -> input_buff_BIN_used) {
		int x1;
		for(x1 = 0 ; x1 < thread_contaxt -> input_buff_BIN_used - thread_contaxt -> input_buff_BIN_ptr; x1++)
			thread_contaxt -> input_buff_BIN[x1] = thread_contaxt -> input_buff_BIN[x1+thread_contaxt -> input_buff_BIN_ptr];
		thread_contaxt -> input_buff_BIN_used -= thread_contaxt -> input_buff_BIN_ptr;
	} else thread_contaxt -> input_buff_BIN_used = 0;

	thread_contaxt -> input_buff_BIN_ptr = 0;

	thread_contaxt -> strm.zalloc = Z_NULL;
	thread_contaxt -> strm.zfree = Z_NULL;
	thread_contaxt -> strm.opaque = Z_NULL;
	thread_contaxt -> strm.avail_in = 0;
	thread_contaxt -> strm.next_in = Z_NULL;

	inflateReset(&thread_contaxt -> strm);

	thread_contaxt -> strm.avail_in = (unsigned int)(thread_contaxt -> input_buff_SBAM_used - thread_contaxt -> input_buff_SBAM_ptr);
	thread_contaxt -> strm.next_in = (unsigned char *)thread_contaxt -> input_buff_SBAM + thread_contaxt -> input_buff_SBAM_ptr;
	thread_contaxt -> strm.avail_out = pairer -> input_buff_BIN_size - thread_contaxt -> input_buff_BIN_used; 
	thread_contaxt -> strm.next_out = (unsigned char *)thread_contaxt -> input_buff_BIN + thread_contaxt -> input_buff_BIN_used;

	int ret = inflate(&thread_contaxt ->strm, Z_FINISH);
	if(ret == Z_OK || ret == Z_STREAM_END)
	{
		int have = pairer -> input_buff_BIN_size - thread_contaxt ->strm.avail_out - thread_contaxt -> input_buff_BIN_used;
		int used_BAM = (unsigned int)(thread_contaxt -> input_buff_SBAM_used - thread_contaxt -> input_buff_SBAM_ptr) - thread_contaxt -> strm.avail_in;
		
		thread_contaxt -> input_buff_BIN_used += have;
		thread_contaxt -> input_buff_SBAM_ptr += used_BAM;
		//SUBREADprintf("BAM USED=%d, BIN GENERATED=%d\n", used_BAM, have);
	} else {
		SUBREADprintf("GZIP ERROR:%d\n", ret);
		return 1;
	}
	return 0;
}

#define BAM_next_u32(v) {\
 (v) = 0; unsigned int poww = 1 ;  \
  BAM_next_nch; (v) += nch*poww; poww *= 256;\
  BAM_next_nch; (v) += nch*poww; poww *= 256;\
  BAM_next_nch; (v) += nch*poww; poww *= 256;\
  BAM_next_nch; (v) += nch*poww;\
}

void SAM_pairer_reduce_BAM_bin(SAM_pairer_context_t * pairer, SAM_pairer_thread_t * thread_contaxt,  unsigned char * bin_where, int * bin_len){
	unsigned int seq_len, name_len, cigar_ops;
	memcpy(&seq_len, bin_where + 20, 4);
	if(seq_len<=1) return;
	memcpy(&name_len, bin_where + 12, 4);
	name_len = name_len & 0xff;
	memcpy(&cigar_ops, bin_where + 16, 4);
	cigar_ops = cigar_ops & 0xffff;

	int targ_pos = 36+name_len+4*cigar_ops + 2;
	int src_pos = 36+name_len+4*cigar_ops + (1+seq_len) / 2 + seq_len;

	bin_where[targ_pos-2]=0xff;
	bin_where[targ_pos-1]=0xff;

	seq_len = 1;
	memcpy(bin_where + 20, &seq_len, 4);
	while(src_pos < (*bin_len)){
		bin_where[targ_pos++]=bin_where[src_pos++];
	}
	(* bin_len) = targ_pos - 4;
	memcpy(bin_where, bin_len, 4);
	(* bin_len) += 4;
	
}

int SAM_pairer_get_next_read_BIN( SAM_pairer_context_t * pairer , SAM_pairer_thread_t * thread_contaxt , unsigned char ** bin_where, int * bin_len ) {
	int read_is_load = 0;
	if( pairer -> input_is_BAM ){
		int nch = 0;
		while(!read_is_load){
			if(!pairer -> BAM_header_parsed){
				int x1;
				unsigned int bam_signature;
				BAM_next_u32(bam_signature);
				BAM_next_u32(pairer -> BAM_l_text);
				char * header_txt = malloc(pairer->BAM_l_text);

			//	SUBREADprintf("HEAD_TXT=%d\n", pairer -> BAM_l_text);
				for(x1 = 0 ; x1 < pairer -> BAM_l_text; x1++){
					BAM_next_nch;
					header_txt [x1] = nch;
				}
				pairer -> output_header(pairer, thread_contaxt -> thread_id, 1, pairer -> BAM_l_text , header_txt , pairer -> BAM_l_text );

				BAM_next_u32(pairer -> BAM_n_ref);
			//	SUBREADprintf("HEAD_CHRO=%d\n", pairer -> BAM_n_ref);
				unsigned int ref_bin_len = 0;
				for(x1 = 0; x1 < pairer -> BAM_n_ref; x1++) {
					unsigned int l_name, l_ref, x2;
					char ref_name[MAX_CHROMOSOME_NAME_LEN];
					BAM_next_u32(l_name);
					memcpy(header_txt + ref_bin_len, &l_name, 4);
					ref_bin_len += 4;
					for(x2 = 0; x2 < l_name; x2++){
						BAM_next_nch;
						header_txt[ref_bin_len++] = nch;
						ref_name[x2]=nch;
					}
					BAM_next_u32(l_ref);
					memcpy(header_txt + ref_bin_len, &l_ref, 4);
					ref_bin_len += 4;

					assert(ref_bin_len < pairer -> BAM_l_text);
					//SUBREADprintf("%d-th ref : %s [len=%u]\n", x1, ref_name, l_ref);
				}

				//exit(0);
				pairer -> output_header(pairer, thread_contaxt -> thread_id, 0, pairer -> BAM_n_ref , header_txt , ref_bin_len );

				free(header_txt);

				pairer -> BAM_header_parsed = 1;
				if(pairer -> display_progress)
					SUBREADprintf("\nThe header was parsed: %d reference sequences were found.\nScanning the input file.\n", pairer -> BAM_n_ref);
				SAM_pairer_fetch_BAM_block(pairer, thread_contaxt);
			}
			while( thread_contaxt -> input_buff_BIN_used == thread_contaxt -> input_buff_BIN_ptr){
				//SUBREADprintf("LOAD BLOCK: BINPTR = %d\n", thread_contaxt -> input_buff_BIN_ptr);
				int ret_fetch = SAM_pairer_fetch_BAM_block(pairer, thread_contaxt);
				//SUBREADprintf("LOAD BLOCK: BINLEN = %d\n", thread_contaxt -> input_buff_BIN_used);
				if(ret_fetch) 
					return 0;
			}		

			unsigned int record_len=0;

			if(thread_contaxt -> input_buff_BIN_ptr + 4 > thread_contaxt -> input_buff_BIN_used )
			{
				if(!pairer -> is_bad_format){
					pairer -> is_bad_format = 1;
					//SUBREADprintf("FORMAT ERROR 1 ; ptr=%u\n", thread_contaxt -> input_buff_BIN_ptr);
					return 0;
				}
			}
			memcpy(&record_len, thread_contaxt -> input_buff_BIN + thread_contaxt -> input_buff_BIN_ptr, 4);
			thread_contaxt -> input_buff_BIN_ptr += 4;

			if(thread_contaxt -> input_buff_BIN_ptr + record_len > thread_contaxt -> input_buff_BIN_used ){
				if(!pairer -> is_bad_format){
					pairer -> is_bad_format = 1;
					//SUBREADprintf("FORMAT ERROR 2 ; len=%u, ptr=%u\n", record_len, thread_contaxt -> input_buff_BIN_ptr);
					return 0;
				}
			}

			//SUBREADprintf("BIN-LEN=%u\n" , record_len);
			(* bin_where) = thread_contaxt -> input_buff_BIN + thread_contaxt -> input_buff_BIN_ptr - 4;
			(* bin_len) = record_len + 4;
			thread_contaxt -> input_buff_BIN_ptr += record_len;

			if( pairer -> tiny_mode )SAM_pairer_reduce_BAM_bin(pairer, thread_contaxt, *bin_where, bin_len);
			//SUBREADprintf("FORMAT OK 2 ; len=%u, ptr=%u\n", record_len, thread_contaxt -> input_buff_BIN_ptr);
			return 1;
		}
	} else {
		return 0;
	}
	return 0;
}

int SAM_pairer_iterate_int_tags(unsigned char * bin, int bin_len, char * tag_name, int * saved_value){
	int found = 0;
	int bin_cursor = 0;
	while(bin_cursor < bin_len){
		if(0){
			char outc[3];
			outc[0] = bin[bin_cursor];
			outc[1] = bin[bin_cursor+1];

			outc[2]=0;
			SUBREADprintf("TAG=%s, TYP=%c %d %c\n", outc, bin[bin_cursor+2],  bin[bin_cursor+3],  bin[bin_cursor+4]);
		}

                if(bin[bin_cursor] == tag_name[0] && bin[bin_cursor+1] == tag_name[1]){
                        int tag_int_val = 0;
                        if(bin[bin_cursor+2]=='i' || bin[bin_cursor+2]=='I'){
                                memcpy(&tag_int_val, bin+bin_cursor+3, 4);
                                found = 1;
                        } else if(bin[bin_cursor+2]=='s' || bin[bin_cursor+2]=='S'){
                                memcpy(&tag_int_val, bin+bin_cursor+3, 2);
                                found = 1;
                        } else if(bin[bin_cursor+2]=='c' || bin[bin_cursor+2]=='C'){
                                memcpy(&tag_int_val, bin+bin_cursor+3, 1);
                                found = 1;
                        }
                        if(found){
                                (* saved_value) = tag_int_val;
                                break;
                        }
                }
                int skip_content = 0;
		//SUBREADprintf("NextTag=%c; ", bin[bin_cursor+2]);
                if(bin[bin_cursor+2]=='i' || bin[bin_cursor+2]=='I' || bin[bin_cursor+2]=='f')
                        skip_content = 4;
                else if(bin[bin_cursor+2]=='s' || bin[bin_cursor+2]=='S')
                        skip_content = 2;
                else if(bin[bin_cursor+2]=='c' || bin[bin_cursor+2]=='C' ||  bin[bin_cursor+2]=='A')
                        skip_content = 1;
		else if(bin[bin_cursor+2]=='Z' || bin[bin_cursor+2]=='H'){
                        while(bin[bin_cursor+skip_content + 3]){
				//SUBREADprintf("ACHAR=%c\n", (bin[skip_content + 3]));
				skip_content++;
			}
			skip_content ++;
                } else if(bin[bin_cursor+2]=='B'){
                        char cell_type = tolower(bin[bin_cursor+3]);
			
                        memcpy(&skip_content, bin + bin_cursor + 4, 4);
		//	SUBREADprintf("Array Type=%c, cells=%d\n", cell_type, skip_content);
                        if(cell_type == 's')skip_content *=2;
                        else if(cell_type == 'i' || cell_type == 'f')skip_content *= 4;
			skip_content += 4 + 1;
                }else{
			assert(0);
		}
		//SUBREADprintf("SKIP=%d\n", skip_content);
                bin_cursor += skip_content + 3;
        }
        return found;
}

int SAM_pairer_get_read_full_name( SAM_pairer_context_t * pairer , SAM_pairer_thread_t * thread_contaxt , unsigned char * bin, int bin_len , char * full_name, int * this_flag){
        full_name[0]=0;
        int rlen = 0;
        if( pairer -> input_is_BAM ){
                unsigned int l_read_name = 0;
                unsigned int refID = 0;
                unsigned int next_refID = 0;
                unsigned int pos = 0, l_seq = 0, cigar_opts;
                unsigned int next_pos = 0, tmpi = 0;
                int FLAG;

                int HItag = -1;


                memcpy(&refID, bin + 4, 4);
                memcpy(&pos, bin + 8, 4);
                memcpy(&tmpi, bin + 12, 4);
		l_read_name = tmpi & 0xff;
		memcpy(&tmpi, bin + 16, 4);
		FLAG = (tmpi >> 16)&0xffff;
		(*this_flag) = FLAG;
		cigar_opts = tmpi & 0xffff;
		memcpy(&next_refID, bin + 24, 4);
		memcpy(&next_pos, bin + 28, 4);
		memcpy(full_name, bin+36, l_read_name);
		unsigned int r1_refID, r1_pos, r2_refID, r2_pos;

		if(FLAG & 4){
			refID = -1;
			pos = 0;
		}

		if(FLAG & 8){
			next_refID = -1;
			next_pos = 0;
		}

		if((FLAG & 0x40) == 0x40){
			r1_refID = refID;
			r1_pos = pos;
			r2_refID = next_refID;
			r2_pos = next_pos;
		} else {
			r2_refID = refID;
			r2_pos = pos;
			r1_refID = next_refID;
			r1_pos = next_pos;
		}


		memcpy(&l_seq, bin + 20, 4);
		//SUBREADprintf("LQ=%d, RL=%d, CIGAR_OPT=%d\n", l_seq, (l_seq+1)/2, cigar_opts);

		unsigned int tags_start = 36+l_read_name+4*cigar_opts+(l_seq+1)/2+l_seq;
		unsigned int tags_len = bin_len - tags_start;

		if(tags_len > 2){
			SAM_pairer_iterate_int_tags(bin + tags_start, tags_len, "HI", &HItag);
		}

		rlen = l_read_name-1 + sprintf(full_name+l_read_name-1, "\027%d\027%u\027%d\027%u\027%d", r1_refID, r1_pos, r2_refID, r2_pos, HItag);
	}
	return rlen;
}

int SAM_pairer_multi_thread_header (void * pairer_vp, int thread_no, int is_text, unsigned int items, char * bin, unsigned int bin_len){

	SAM_pairer_context_t * pairer = (SAM_pairer_context_t *) pairer_vp;
	SAM_pairer_writer_main_t * bam_main = (SAM_pairer_writer_main_t * )pairer -> appendix1; 
	SAM_pairer_writer_thread_t * bam_thread = bam_main -> threads + thread_no;
	unsigned int BIN_block_cursor = 0, bin_cursor = 0;
	if(is_text){
		memcpy( bam_thread -> BIN_buffer, "BAM\1", 4 );
		memcpy( bam_thread -> BIN_buffer + 4 , & items , 4 );
		BIN_block_cursor = 8;
	}else{
		memcpy( bam_thread -> BIN_buffer , & items , 4 );
		BIN_block_cursor = 4;
	}
	while( bin_cursor  < bin_len ){
		int write_text_len = min(SAM_PAIRER_WRITE_BUFFER - BIN_block_cursor, bin_len - bin_cursor);
		memcpy(bam_thread -> BIN_buffer + BIN_block_cursor , bin + bin_cursor, write_text_len);
		bam_thread -> BIN_buffer_ptr = write_text_len + BIN_block_cursor;

		SAM_pairer_multi_thread_compress(bam_main, bam_thread);
		bin_cursor += write_text_len;
		BIN_block_cursor = 0;
	}

	bam_thread -> BIN_buffer_ptr = 0;
	return 0;
}

void SAM_pairer_make_dummy(char * rname, char * bin1, char * out_bin2){
	char * tmptr = NULL;

	//SUBREADprintf("S=%s  ", rname);
	char * realname = strtok_r(rname, "\027", &tmptr);
	int len_name = strlen(realname);
	int r1_chro = atoi(strtok_r(NULL, "\027", &tmptr));
	int r1_pos = atoi(strtok_r(NULL, "\027", &tmptr));
	int r2_chro = atoi(strtok_r(NULL, "\027", &tmptr));
	int r2_pos = atoi(strtok_r(NULL, "\027", &tmptr));
	int HItag = atoi(strtok_r(NULL, "\027", &tmptr));
	int mate_FLAG = 0;
	memcpy(&mate_FLAG, bin1 + 16, 4);
	mate_FLAG = 0xffff&(mate_FLAG >>16);
	int mate_tlen = 0;
	memcpy(&mate_tlen, bin1 + 32, 4);

	if(r1_chro<0) r1_pos=-1;
	if(r2_chro<0) r2_pos=-1;

	int my_chro = (mate_FLAG&0x40)? r2_chro : r1_chro;
	int my_pos = (mate_FLAG&0x40)? r2_pos : r1_pos;
	int mate_chro = (mate_FLAG&0x40)? r1_chro : r2_chro;
	int mate_pos = (mate_FLAG&0x40)? r1_pos : r2_pos;

	int bin_mq_nl = (len_name+1);
	int my_flag = (mate_FLAG&0x40)? 0x80:0x40;
	my_flag |= 1;
	if(mate_FLAG & 8)my_flag |=4;
	if(mate_FLAG & 4)my_flag |=8;
	if(mate_FLAG & 0x10) my_flag |= 0x20;
	if(mate_FLAG & 0x20) my_flag |= 0x10;
	my_flag = my_flag << 16;

	memcpy(out_bin2+4, &my_chro,4);
	memcpy(out_bin2+8, &my_pos,4);
	memcpy(out_bin2+12, &bin_mq_nl, 4);
	memcpy(out_bin2+16, &my_flag, 4);

	my_flag = 1;
	memcpy(out_bin2+20, &my_flag, 4);
	memcpy(out_bin2+24, &mate_chro, 4); 
	memcpy(out_bin2+28, &mate_pos, 4); 

	mate_tlen = -mate_tlen;
	memcpy(out_bin2+32, &mate_tlen, 4);
	memcpy(out_bin2+36, realname, len_name+1);
	out_bin2[36 + len_name+1] = 0xff;
	out_bin2[36 + len_name+2] = 0x20;

	int all_len = 36 + len_name + 3 - 4;
	//SUBREADprintf("HI=%d\n", HItag);
	if(HItag>=0){
		out_bin2[36 + len_name+3]='H';
		out_bin2[36 + len_name+4]='I';
		if(HItag<128){
			out_bin2[36 + len_name+5]='C';
			memcpy(out_bin2 + 36 + len_name+6, &HItag, 1);
			all_len += 4;
		}else if(HItag<32767){
			out_bin2[36 + len_name+5]='S';
			memcpy(out_bin2 + 36 + len_name+6, &HItag, 2);
			all_len += 5;
		}else {
			out_bin2[36 + len_name+5]='I';
			memcpy(out_bin2 + 36 + len_name+6, &HItag, 4);
			all_len += 7;
		}
	}
	memcpy(out_bin2,&all_len,4);
}

void SAM_pairer_reset( SAM_pairer_context_t * pairer ) {
	int x1;
	pairer -> BAM_header_parsed = 0;
	pairer -> total_input_reads = 0;
	for(x1 = 0; x1 < pairer -> total_threads ; x1 ++){
		pairer -> threads[x1].input_buff_BIN_used = 0;
		pairer -> threads[x1].input_buff_BIN_ptr = 0;
		pairer -> threads[x1].input_buff_SBAM_used = 0;
		pairer -> threads[x1].input_buff_SBAM_ptr = 0;
		pairer -> threads[x1].orphant_block_no = 0;
		HashTableDestroy(pairer -> threads[x1].orphant_table);
		pairer -> threads[x1].orphant_table = HashTableCreate(pairer -> input_buff_SBAM_size / 100);
		HashTableSetHashFunction(pairer -> threads[x1].orphant_table, fc_chro_hash);
		HashTableSetKeyComparisonFunction(pairer -> threads[x1].orphant_table, fc_strcmp_chro);
		HashTableSetDeallocationFunctions(pairer -> threads[x1].orphant_table, free, free);
		inflateReset(&pairer -> threads[x1].strm);
	}

}
void SAM_pairer_writer_reset( void * pairer_vp ) {
	SAM_pairer_context_t * pairer = (SAM_pairer_context_t *) pairer_vp;
	SAM_pairer_writer_main_t * bam_main = (SAM_pairer_writer_main_t * )pairer -> appendix1;
	ftruncate(fileno(bam_main -> bam_fp), 0);
	fclose(bam_main -> bam_fp);
	bam_main -> bam_fp = f_subr_open(bam_main -> bam_name, "wb");
	int x1;
	for(x1 = 0; x1 < pairer -> total_threads ; x1 ++){
		bam_main -> threads[x1].BIN_buffer_ptr = 0;
		deflateReset(&bam_main -> threads[x1].strm);
	}


}

int SAM_pairer_multi_thread_output(void * pairer_vp, int thread_no, char * rname, char * bin1, char * bin2 ){
	SAM_pairer_context_t * pairer = (SAM_pairer_context_t *) pairer_vp;
	SAM_pairer_writer_main_t * bam_main = (SAM_pairer_writer_main_t * )pairer -> appendix1; 
	SAM_pairer_writer_thread_t * bam_thread = bam_main -> threads + thread_no;

	char dummy_bin2 [MAX_READ_NAME_LEN*2 + 180 ];
	if(bin2==NULL && bam_main -> has_dummy){
		SAM_pairer_make_dummy( rname, bin1, dummy_bin2 );
		bin2 = dummy_bin2;
	}

	int bin_len1, bin_len2 = 0;
	memcpy(&bin_len1, bin1, 4);
	bin_len1 +=4;

	if(bin2) {
		memcpy(&bin_len2, bin2, 4);
		bin_len2 +=4;
	}

	if( bin_len1 + bin_len2 >= SAM_PAIRER_WRITE_BUFFER){
		SUBREADprintf("ERROR: BAM Record larger than a BAM block!\n");
		return 1;
	}

	if(bin_len1 + bin_len2 + bam_thread -> BIN_buffer_ptr >= SAM_PAIRER_WRITE_BUFFER){
		int ret = SAM_pairer_multi_thread_compress(bam_main, bam_thread);
		if(ret)return 1;
	} 
	memcpy( bam_thread -> BIN_buffer + bam_thread -> BIN_buffer_ptr, bin1, bin_len1 );
	if(bin2)
		memcpy( bam_thread -> BIN_buffer + bam_thread -> BIN_buffer_ptr + bin_len1, bin2, bin_len2 );
	bam_thread -> BIN_buffer_ptr += bin_len1 + bin_len2;
	return 0;
}

int SAM_pairer_do_next_read( SAM_pairer_context_t * pairer , SAM_pairer_thread_t * thread_contaxt ){
	char read_full_name[ MAX_READ_NAME_LEN*2 +80 ];	// rname:chr_r1:pos_r1:chr_r2:pos_r2:HI_tag
	unsigned char * bin = NULL;
	int bin_len = 0, this_flags = 0;

	int has_next_read = SAM_pairer_get_next_read_BIN(pairer, thread_contaxt, &bin, &bin_len);
	if(has_next_read){
		int name_len = SAM_pairer_get_read_full_name(pairer, thread_contaxt, bin, bin_len, read_full_name, & this_flags);
		//SUBREADprintf("FNNM:%s\n", read_full_name);
		unsigned char * mate_bin = HashTableGet(thread_contaxt -> orphant_table, read_full_name);
		if(mate_bin){
			if(pairer -> output_function)
				pairer -> output_function(pairer, thread_contaxt -> thread_id, read_full_name, (char*) bin, (char*)mate_bin);
			HashTableRemove(thread_contaxt -> orphant_table, read_full_name);
			if(thread_contaxt -> orphant_space > bin_len)
				thread_contaxt -> orphant_space -= bin_len;
			else	thread_contaxt -> orphant_space = 0;
			//SUBREADprintf("Mate_found: %s\n", read_full_name);
		}else{
			char * mem_name = malloc(name_len + 1);
			memcpy(mem_name, read_full_name, name_len);
			mem_name[name_len] = 0;

			char * mem_bin = malloc(bin_len);
			memcpy(mem_bin, bin , bin_len);

			HashTablePut(thread_contaxt -> orphant_table, mem_name, mem_bin);
			thread_contaxt -> orphant_space += bin_len;
			//SUBREADprintf("Orphant_created [%d]: %s\n", thread_contaxt -> thread_id, read_full_name);
		}
		return 0;
	}
	return 1;
}


// all orphants are written into files, each has a size of buffer size.
// when the orphants are longer than buffer_size, then sort and save to disk.

void SAM_pairer_sort_exchange(void * arr, int l, int r){
	unsigned char *** sort_data = (unsigned char ***) arr;
	unsigned char * tmpc;

	tmpc = sort_data[0][r];
	sort_data[0][r] = sort_data[0][l];
	sort_data[0][l] = tmpc;

	tmpc = sort_data[1][r];
	sort_data[1][r] = sort_data[1][l];
	sort_data[1][l] = tmpc;
}

int SAM_pairer_sort_compare(void * arr, int l, int r){
	char *** sort_data = (char ***) arr;
	return strcmp(sort_data[0][l], sort_data[0][r]);
}

void SAM_pairer_sort_merge( void * arr, int start, int items, int items2 ){
	unsigned char *** sort_data = (unsigned char ***) arr;

	unsigned char ** tmp_name_list = malloc(sizeof(char *) * (items+items2));
	unsigned char ** tmp_bin_list = malloc(sizeof(char *) * (items+items2));

	int i1_cursor = start, i2_cursor = items + start;
	int tmp_cursor = 0;

	while(1){
		if(i1_cursor == items + start && i2_cursor == items + items2 + start )break;
		int select_items_1 = (i2_cursor == start + items + items2) || (i1_cursor < items + start && SAM_pairer_sort_compare(arr, i1_cursor, i2_cursor) <= 0);
		if(select_items_1){
			tmp_name_list[tmp_cursor] = sort_data[0][i1_cursor];
			tmp_bin_list[tmp_cursor ++] = sort_data[1][i1_cursor++];
		}else{
			tmp_name_list[tmp_cursor] = sort_data[0][i2_cursor];
			tmp_bin_list[tmp_cursor ++] = sort_data[1][i2_cursor++];
		}
	}
	assert(tmp_cursor == items + items2);

	memcpy( sort_data[0] + start, tmp_name_list, sizeof(char *) * (items+items2) );
	memcpy( sort_data[1] + start, tmp_bin_list, sizeof(char *) * (items+items2) );
	free(tmp_name_list);
	free(tmp_bin_list);
	
}

unsigned int SAM_pairer_osr_hash(char * st){
	int x1 = 0, nch;
	unsigned int ret = 0, ret2=0;
	while((nch = st[x1++])!=0){
		ret = (ret << 2) ^ nch;
		ret2 = (ret << 3) ^ nch;
	}
	return (ret^ret2) % 39846617;
}

int SAM_pairer_osr_next_name(FILE * fp , char * name, int thread_no, int all_threads){
	while(1){
		if(feof(fp)) return 0;
		int rlen =0;
		fread(&rlen, 1, 2, fp);
		if(rlen<1) return 0;
		assert(rlen < 1024);

		int rlen2 = fread(name, 1, rlen, fp);
		if(rlen2 != rlen) return 0;
		name[rlen]=0;
		if( SAM_pairer_osr_hash(name)% all_threads == thread_no  )
		{
			fseek(fp, -2-rlen, SEEK_CUR);
			return 1;
		}
		fread(&rlen, 1, 2, fp);
		assert(rlen < 65535);
		rlen +=4;
		fseek(fp, rlen, SEEK_CUR);
	}
	return 0;
}

void SAM_pairer_osr_next_bin(FILE * fp, char * bin){
	int rlen =0;
	fread(&rlen, 1, 2, fp);
	assert(rlen < 1024);
	fseek(fp, rlen, SEEK_CUR);
	rlen =0;
	fread(&rlen, 1, 2, fp);
	assert(rlen < 65535);
	rlen +=4;
	fread(bin, 1, rlen, fp);
}

void * SAM_pairer_rescure_orphants(void * params){

	void ** param_ptr = (void **) params;
	SAM_pairer_context_t * pairer = param_ptr[0];
	int thread_no = (int)(param_ptr[1]-NULL);
	free(params);

	int orphant_fp_size = 50, orphant_fp_no=0;
	FILE ** orphant_fps = malloc(sizeof(FILE *) * orphant_fp_size);
	int thno, bkno, x1;
	char * bin_tmp1 , * bin_tmp2;

	if(0 == thread_no && pairer -> display_progress)
		SUBREADprintf("Finished scanning the input file. Processing unpaired reads.\n");

	bin_tmp1 = malloc(66000);
	bin_tmp2 = malloc(66000);

	for( thno = 0 ; thno < pairer -> total_threads ; thno ++ ){
		for( bkno = 0 ; ; bkno++){
			char tmp_fname[MAX_FILE_NAME_LENGTH];
			sprintf(tmp_fname, "%s-TH%02d-BK%06d.tmp", pairer->tmp_file_prefix,  thno, bkno);

			FILE * in_fp = fopen(tmp_fname, "rb");
			if(NULL == in_fp) break;
			if(orphant_fp_no >= orphant_fp_size){
				orphant_fp_size *= 1.5;
				orphant_fps = realloc(orphant_fps, orphant_fp_size * sizeof(FILE *));
			}
			orphant_fps[orphant_fp_no++]=in_fp;
		}
	}

	int max_name_len = MAX_READ_NAME_LEN*2 +80;
	char * names = malloc( orphant_fp_no * max_name_len );
	memset(names, 0, orphant_fp_no * max_name_len );
	
	for(x1 = 0 ; x1 < orphant_fp_no; x1++)
	{
		int has = SAM_pairer_osr_next_name( orphant_fps[x1] , names + max_name_len*x1 , thread_no , pairer-> total_threads);
		if(!has) *(names + max_name_len*x1)=0;
	}

	unsigned long long rescured=0, died=0;

	while(1){
		int min_name_fileno = -1;
		int min2_name_fileno = -1;

		for(x1 = 0 ; x1 < orphant_fp_no; x1++){
			int has = *(names + max_name_len*x1);
			if(has){
				int strcv_12 = 1;
				if(min_name_fileno >=0) strcv_12 = strcmp(names+(min_name_fileno * max_name_len), names+(x1 * max_name_len));
				if(strcv_12 > 0){
					min_name_fileno = x1;
					min2_name_fileno = -1;
				}else if( strcv_12 == 0){
					min2_name_fileno = x1;
				}
			}

		}

		if(min_name_fileno >= 0){
			SAM_pairer_osr_next_bin( orphant_fps[ min_name_fileno ] , bin_tmp1);

			if( min2_name_fileno >=0){
				SAM_pairer_osr_next_bin( orphant_fps[ min2_name_fileno ] , bin_tmp2);
				pairer -> output_function(pairer, thread_no,  names + max_name_len*min_name_fileno , (char*) bin_tmp1, (char*)bin_tmp2);

				int read_has = SAM_pairer_osr_next_name( orphant_fps[min2_name_fileno],  names + max_name_len*min2_name_fileno, thread_no,  pairer-> total_threads);
				if(!read_has) *(names + max_name_len*min2_name_fileno)=0;
				rescured++;
			}else{
				//SUBREADprintf("FINAL_ORPHAN:%s\n" , names + max_name_len*min_name_fileno);
				pairer -> output_function(pairer, thread_no,  names + max_name_len*min_name_fileno, (char*) bin_tmp1, NULL);
				died++;
			}

			int read_has = SAM_pairer_osr_next_name( orphant_fps[min_name_fileno],  names + max_name_len*min_name_fileno, thread_no, pairer-> total_threads);
			if(!read_has) *(names + max_name_len*min_name_fileno)=0;
		} else break;
	}
	free(names);

	for(x1 = 0 ; x1 < orphant_fp_no; x1++)
	{
		fclose ( orphant_fps[x1] );
	}
	free( bin_tmp1 );
	free( bin_tmp2 );
	pairer -> total_orphan_reads += died;
	//SUBREADprintf("RESCURE THREAD %d   Rescured %llu, Died %llu\n", thread_no, rescured, died);
	return NULL;
}

void SAM_pairer_update_orphant_table(SAM_pairer_context_t * pairer , SAM_pairer_thread_t * thread_context){
	unsigned int x2 = 0;
	unsigned char ** name_list, ** bin_list;
	name_list = malloc(sizeof(char*) * thread_context->orphant_table->numOfElements);
	bin_list  = malloc(sizeof(char*) * thread_context->orphant_table->numOfElements);

	int x1;
	for(x1 = 0; x1 < thread_context->orphant_table->numOfBuckets; x1 ++){
		KeyValuePair *pair = thread_context->orphant_table->bucketArray[x1];
		while (pair != NULL) {
			KeyValuePair *nextPair = pair->next;
			name_list [x2] = (unsigned char *)pair -> key;
			bin_list [x2] = pair -> value;
			x2++;
			pair = nextPair;
		}
	}

	assert(x2 == thread_context->orphant_table->numOfElements);
	unsigned char ** sort_data[2];
	sort_data[0]=name_list;
	sort_data[1]=bin_list;
	merge_sort(sort_data, thread_context->orphant_table->numOfElements, SAM_pairer_sort_compare, SAM_pairer_sort_exchange, SAM_pairer_sort_merge);

	char tmp_fname[MAX_FILE_NAME_LENGTH];
	sprintf(tmp_fname, "%s-TH%02d-BK%06d.tmp", pairer->tmp_file_prefix, thread_context -> thread_id, thread_context -> orphant_block_no++);
	FILE * tmp_fp = fopen(tmp_fname, "wb");

	for(x1 = 0; x1 < x2;  x1 ++){
		unsigned int bin_len;

		memcpy(&bin_len, bin_list[x1] , 4);
		int namelen = strlen((char *)name_list[x1]);

		fwrite(&namelen,1,2,tmp_fp);
		fwrite(name_list[x1], 1, namelen, tmp_fp);
		fwrite(&bin_len,1,2,tmp_fp);
		fwrite(bin_list[x1],  1, bin_len + 4, tmp_fp);

		HashTableRemove(thread_context->orphant_table , name_list[x1]);
	}
	assert(thread_context -> orphant_table-> numOfElements == 0);
	fclose(tmp_fp);
	free(name_list);	
	free(bin_list);	
	thread_context -> orphant_space = 0;
}


#define PAIRER_WAIT_TICK_TIME 10000

void * SAM_pairer_thread_run( void * params ){
	void ** param_ptr = (void **) params;
	SAM_pairer_context_t * pairer = param_ptr[0];
	int thread_no = (int)(param_ptr[1]-NULL);
	free(params);

	SAM_pairer_thread_t * thread_contaxt = pairer -> threads + thread_no;
	int is_finished = 0;
	while(1){
		subread_lock_occupy(&pairer -> input_fp_lock);
		if(pairer -> BAM_header_parsed || thread_no == 0){
			SAM_pairer_fill_BIN_buff(pairer, thread_contaxt, &is_finished);
		}
		subread_lock_release(&pairer -> input_fp_lock);
		if(!pairer -> BAM_header_parsed && thread_no > 0) {
			usleep(PAIRER_WAIT_TICK_TIME);
		} else if(thread_contaxt -> input_buff_SBAM_used>0) {
			unsigned int processed_reads = 0;
			while(1){
				int has_no_more = SAM_pairer_do_next_read(pairer, thread_contaxt);
				if(has_no_more)break;
				processed_reads++;
			}

			pairer -> total_input_reads += processed_reads;
		}
		if(pairer -> is_bad_format) break;
		if(thread_contaxt -> orphant_space > pairer -> input_buff_SBAM_size)
			SAM_pairer_update_orphant_table(pairer, thread_contaxt);

		if(is_finished) break;
	}

	if(thread_contaxt -> orphant_table -> numOfElements > 0)
		SAM_pairer_update_orphant_table(pairer, thread_contaxt);

	return NULL;
}


// not only run, but also finalise.
// It returns 0 if no error.
int SAM_pairer_run_once( SAM_pairer_context_t * pairer){
	int x1;
	for(x1 = 0; x1 < pairer -> total_threads ; x1++){
		// this 16-byte memory block is freed in the thread worker.
		void ** init_params = malloc(sizeof(void *) * 2);

		init_params[0] = pairer;
		init_params[1] = (void *)(NULL+x1);
		pthread_create(&(pairer -> threads[x1].thread_stab), NULL, SAM_pairer_thread_run, init_params);
	}

	for(x1 = 0; x1 < pairer -> total_threads ; x1++){
		pthread_join(pairer -> threads[x1].thread_stab, NULL);
	}

	if(0 == pairer -> is_bad_format){
		for(x1 = 0; x1 < pairer -> total_threads ; x1++){
			// this 16-byte memory block is freed in the thread worker.

			void ** init_params = malloc(sizeof(void *) * 2);

			init_params[0] = pairer;
			init_params[1] = (void *)(NULL+x1);
			pthread_create(&(pairer -> threads[x1].thread_stab), NULL, SAM_pairer_rescure_orphants, init_params);
		}

		for(x1 = 0; x1 < pairer -> total_threads ; x1++){
			pthread_join(pairer -> threads[x1].thread_stab, NULL);
		}
	}

	return 0;
}

int fix_load_next_block(FILE * in, char * binbuf, z_stream * strm){
	char * bam_buf = malloc(70000);
	int x1, ret = 0;
	x1 = fgetc(in);
	if(x1 != 31) ret = -1;
	x1 = fgetc(in);
	if(x1 != 139) ret = -1;
	x1 = fgetc(in);
	if(x1 != 8) ret = -1;
	x1 = fgetc(in);
	if(x1 != 4) ret = -1;
	if(ret == 0){
		x1 = fgetc(in);
		x1 = fgetc(in);
		x1 = fgetc(in);
		x1 = fgetc(in);

		x1 = fgetc(in);//XFL

		x1 = fgetc(in);//OS
		int xlen;
		xlen = fgetc(in);
		xlen += fgetc(in) * 256;
		int bsize = -1, xlen_ptr = 0;

		while(xlen_ptr < xlen){
			int si1 = fgetc(in);
			int si2 = fgetc(in);
			int slen = fgetc(in);
			slen += fgetc(in) * 256;
			if(si1 == 66 && si2==67){
				bsize = fgetc(in);
				bsize += 256*fgetc(in);
			}else{
				fseek(in , slen, SEEK_CUR);
			}
			xlen_ptr += 4 + slen;
		}
		if(bsize > 0){
			fread(bam_buf, 1, bsize - xlen - 19, in);
		}
		fseek(in, 8, SEEK_CUR);

		strm -> avail_in = bsize - xlen - 19;
		strm -> next_in = (unsigned char*)bam_buf;
		strm -> avail_out = 70000;
		strm -> next_out = (unsigned char*)binbuf;
		int ret_inf = inflate(strm, Z_FINISH);
		if(ret_inf == Z_STREAM_END){
			ret = 70000 - strm -> avail_out;
		//	SUBREADprintf("FIX_DECOM: %d -> %d\n", bsize - xlen - 19, ret);
		}else{
			SUBREADprintf("FIX_DECOM_ERR:%d\n" , ret_inf);
			ret = -1;
		}
		inflateReset(strm);
	}
	free(bam_buf);
	return ret;
}

void fix_write_block(FILE * out, char * bin, int binlen, z_stream * strm){
	char * bam_buf = malloc(70000);
	int x1, bam_len = 0, retbam;

	if(binlen > 0){
		strm -> avail_in = binlen;
		strm -> next_in = (unsigned char*)bin;
		strm -> avail_out = 70000;
		strm -> next_out = (unsigned char*)bam_buf;
		retbam = deflate(strm , Z_FINISH);
		bam_len = 70000 - strm -> avail_out;
		deflateReset(strm);
	}else{
		z_stream nstrm;
		nstrm.zalloc = Z_NULL;
		nstrm.zfree = Z_NULL;
		nstrm.opaque = Z_NULL;
		nstrm.avail_in = 0;
		nstrm.next_in = Z_NULL;
	
		deflateInit2(&nstrm, SAMBAM_COMPRESS_LEVEL, Z_DEFLATED,
			PAIRER_GZIP_WINDOW_BITS, PAIRER_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);

		nstrm.avail_in = 0;
		nstrm.next_in = (unsigned char*)bin;
		nstrm.avail_out = 70000;
		nstrm.next_out = (unsigned char*)bam_buf;
		retbam = deflate(&nstrm, Z_FINISH);
		bam_len = 70000 - nstrm.avail_out;
		deflateEnd(&nstrm);
	}

	//SUBREADprintf("FIX_COMPR: %d -> %d  RET=%d\n", binlen , bam_len, retbam);

	unsigned int crc0 = crc32(0, NULL, 0);
	unsigned int crc = crc32(crc0, (unsigned char *) bin , binlen);

	fputc(31, out);
	fputc(139, out);
	fputc(8, out);
	fputc(4, out);
	fputc(0, out);
	fputc(0, out);
	fputc(0, out);
	fputc(0, out);

	fputc(0, out);//XFL
	fputc(0xff, out);//OS

	x1 = 6;
	fwrite( &x1, 2, 1 , out );
	fputc( 66, out );
	fputc( 67, out );
	x1 = 2;
	fwrite( &x1, 2, 1 , out );
	x1 = bam_len + 19 + 6;
	fwrite( &x1, 2, 1 , out );
	fwrite( bam_buf , 1,bam_len, out );
	
	fwrite( &crc, 4, 1, out );
	fwrite( &binlen, 4, 1, out );

	free(bam_buf);
}

#define FIX_GET_NEXT_NCH { while(in_bin_ptr == in_bin_size){ \
  in_bin_ptr = 0; in_bin_size = 0;\
  int newsize = fix_load_next_block(old_fp, in_bin, &in_strm);\
  if(newsize < 0){ break;}else{in_bin_size = newsize;}\
} if(in_bin_size>0){nch = in_bin[in_bin_ptr++];  if(nch < 0)nch += 256; } else nch = -1; } 

#define FIX_FLASH_OUT { if(out_bin_ptr > 0) fix_write_block(new_fp, out_bin, out_bin_ptr, &out_strm); out_bin_ptr = 0; }

#define FIX_APPEND_OUT(p, c) { if(out_bin_ptr > 60000){FIX_FLASH_OUT} ;  memcpy(out_bin + out_bin_ptr, p, c); out_bin_ptr +=c ; }
#define FIX_APPEND_READ(p, c){ memcpy(out_bin + out_bin_ptr, p, c); out_bin_ptr +=c ;  }

void SAM_pairer_fix_format(SAM_pairer_context_t * pairer){
	FILE * old_fp = pairer -> input_fp;
	fseek(old_fp, 0, SEEK_SET);
	char tmpfname [300];

	sprintf(tmpfname, "%s.fixbam", pairer -> tmp_file_prefix);

	FILE * new_fp = f_subr_open(tmpfname, "wb");
	char * in_bin = malloc(140000);
	char * out_bin = malloc(70000);

	z_stream in_strm;
	z_stream out_strm;
	in_strm.zalloc = Z_NULL;
	in_strm.zfree = Z_NULL;
	in_strm.opaque = Z_NULL;
	in_strm.avail_in = 0;
	in_strm.next_in = Z_NULL;
	
	inflateInit2(&in_strm, PAIRER_GZIP_WINDOW_BITS);

	out_strm.zalloc = Z_NULL;
	out_strm.zfree = Z_NULL;
	out_strm.opaque = Z_NULL;
	out_strm.avail_in = 0;
	out_strm.next_in = Z_NULL;
	
	deflateInit2(&out_strm, Z_NO_COMPRESSION, Z_DEFLATED,
                PAIRER_GZIP_WINDOW_BITS, PAIRER_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);

	int in_bin_ptr = 0;
	int out_bin_ptr = 0;
	int in_bin_size = 0;
	int content_count = 0;
	int content_size = 0;
	int x1, nch = 0;

	for(x1 = 0; x1 < 4; x1++){
		FIX_GET_NEXT_NCH; // BAM1
		FIX_APPEND_OUT(&nch, 1);
	}


	// ====== The header texts
	content_size = 0;
	for(x1 = 0; x1 < 4; x1++){
		FIX_GET_NEXT_NCH;
	//	SUBREADprintf("FIX: TLEN: %d\n", nch);
		content_size += (nch << (8 * x1));
	}
	FIX_APPEND_OUT(&content_size, 4);
	//SUBREADprintf("FIX: TXTLEN=%d\n", content_size);
	for(content_count = 0; content_count < content_size; content_count++){
		FIX_GET_NEXT_NCH;
		FIX_APPEND_OUT(&nch, 1);
	//	fputc(nch, stderr);
	}
	FIX_FLASH_OUT;

	// ====== The chromosome table
	content_size = 0;
	for(x1 = 0; x1 < 4; x1++){
		FIX_GET_NEXT_NCH;
		content_size += (nch << (8 * x1));
	}
	FIX_APPEND_OUT(&content_size, 4);
	//SUBREADprintf("FIX: CHROLEN=%d\n", content_size);
	for(content_count = 0; content_count < content_size; content_count++){
		int namelen = 0;
		for(x1 = 0; x1 < 4; x1++){
			FIX_GET_NEXT_NCH;
			namelen+= (nch << (8 * x1));
		}
		FIX_APPEND_READ(&namelen, 4);
		for(x1 = 0; x1 <  namelen + 4; x1++){ // inc. length
			FIX_GET_NEXT_NCH;
			FIX_APPEND_READ(&nch, 1);
		}

		if(out_bin_ptr > 60000){
			FIX_FLASH_OUT;
		}
	}
	FIX_FLASH_OUT;

	// ===== The reads
	unsigned long long reads =0;
	while(1){
		int block_size = 0;
		FIX_GET_NEXT_NCH;
		if(nch<0) break;
		FIX_APPEND_READ(&nch, 1);
		block_size = nch;
		for(x1 = 1; x1 < 4; x1++){
			FIX_GET_NEXT_NCH;
			block_size += (nch << (8 * x1));
			FIX_APPEND_READ(&nch, 1);
		}
		for(x1 = 0; x1 < block_size; x1++){
			FIX_GET_NEXT_NCH;
			FIX_APPEND_READ(&nch, 1);
		}
		reads ++;
		if(out_bin_ptr > 60000){
			FIX_FLASH_OUT;
		}
	}
	FIX_FLASH_OUT;
	//SUBREADprintf("FIX READS=%llu\n", reads);
	fix_write_block(new_fp, out_bin, 0, &out_strm);
	deflateEnd(&out_strm);
	inflateEnd(&in_strm);

	fclose(old_fp);
	fclose(new_fp);

	pairer -> input_fp = f_subr_open(tmpfname, "rb");
	free(in_bin);
	free(out_bin);
}

int SAM_pairer_run( SAM_pairer_context_t * pairer){
	int corrected_run;
	for(corrected_run = 0; corrected_run < 2 ; corrected_run ++){
		SAM_pairer_run_once(pairer);
		if(pairer -> is_bad_format){
			assert(0 == corrected_run);
			if(pairer -> display_progress)
				SUBREADprintf("Retrying with the corrected format...\n");
			delete_with_prefix(pairer -> tmp_file_prefix);
			SAM_pairer_fix_format(pairer);
			SAM_pairer_reset(pairer);
			pairer -> reset_output_function(pairer);
			pairer -> is_bad_format = 0;
		}else break;
	}

	return 0;
}

int sort_SAM_create(SAM_sort_writer * writer, char * output_file, char * tmp_path)
{
	char tmp_fname[MAX_FILE_NAME_LENGTH+40];
	memset(writer, 0, sizeof(SAM_sort_writer));

	old_sig_TERM = signal (SIGTERM, SAM_SORT_SIGINT_hook);
	old_sig_INT = signal (SIGINT, SAM_SORT_SIGINT_hook);

	sprintf(writer -> tmp_path, "%s/temp-sort-%06u-%08X-", tmp_path, getpid(), rand());
	_SAMSORT_SNP_delete_temp_prefix = writer -> tmp_path;

	sprintf(tmp_fname, "%s%s", writer -> tmp_path, "headers.txt");
	writer -> all_chunks_header_fp = f_subr_open(tmp_fname,"w");
	if(!writer -> all_chunks_header_fp) return -1;
	fclose(writer -> all_chunks_header_fp);
	unlink(tmp_fname);

	writer -> out_fp = f_subr_open(output_file,"w");
	if(!writer -> out_fp) return -1;

	return 0;
}

void find_tag_out(char * read_line_buf, char * tag, char * hi_tag_out)
{
	int hi_tag = -1;
	char tag_str[10];
	sprintf(tag_str , "\t%s:i:", tag);
	char * hi_tag_str = strstr(read_line_buf, tag_str);
	if(hi_tag_str)
	{


		hi_tag = 0;
		int line_cursor;
		for(line_cursor=6; ; line_cursor++)
		{
			char nch = hi_tag_str[line_cursor];
//								printf("HI:i=%s; nch [%d] ='%c'\n", hi_tag_str, line_cursor, nch);
			if(!isdigit(nch)) break;
			hi_tag = hi_tag*10 + (nch-'0');
		}
	}

	if(hi_tag >=0)
	{
		sprintf(hi_tag_out,"\t%s:i:%d", tag, hi_tag);
	}else hi_tag_out[0] = 0;


}

void sort_SAM_finalise(SAM_sort_writer * writer)
{
	int x1_chunk, x1_block;
	int xk1;
	for(xk1=0;xk1<SAM_SORT_BLOCKS;xk1++)
	{
		if(writer -> current_block_fp_array[xk1])
			fclose(writer -> current_block_fp_array[xk1]);
	}
	memset(writer -> current_block_fp_array, 0, sizeof(FILE *)*SAM_SORT_BLOCKS);
	writer -> current_chunk_size = 0;
	writer -> current_chunk++;

	for(x1_block = 0; x1_block <SAM_SORT_BLOCKS; x1_block++){  
		HashTable * first_read_name_table;
		first_read_name_table = HashTableCreate(SAM_SORT_BLOCK_SIZE / 100 );
		HashTableSetKeyComparisonFunction(first_read_name_table , fc_strcmp_chro);
		HashTableSetDeallocationFunctions(first_read_name_table , free, free);
		HashTableSetHashFunction(first_read_name_table, HashTableStringHashFunction);

		for(x1_chunk = 0; x1_chunk < writer -> current_chunk; x1_chunk++)
		{
			char tmpfname[MAX_FILE_NAME_LENGTH+40];
			sprintf(tmpfname, "%sCHK%08d-BLK%03d.bin", writer -> tmp_path, x1_chunk , x1_block);

			FILE * bbfp = f_subr_open(tmpfname,"rb");
			if(!bbfp) continue;

			while(!feof(bbfp))
			{
				char * read_name = NULL;
				short flags;
				short read_name_len;
				short read_len;
				int ret = fread(&flags, 2,1 , bbfp);
				if(ret<1) break;
				fread(&read_name_len, 2,1 , bbfp);
				if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
					fseek(bbfp, read_name_len, SEEK_CUR); 
				else
				{
					read_name = malloc(read_name_len+1);
					fread(read_name, 1, read_name_len, bbfp);
					read_name[read_name_len] = 0;
				}
				fread(&read_len,2,1,bbfp);
				if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
					fseek(bbfp, read_len, SEEK_CUR); 
				else
				{
					char * new_line_mem = malloc(read_len+1);
					fread(new_line_mem, 1, read_len, bbfp);
					new_line_mem[read_len] = 0;

					if(read_len<2)
					{
						SUBREADprintf("Cannot determain read length from the tmp file!\n");
						assert(0);
					}


					if( new_line_mem[0]==0 || new_line_mem[1]==0)
					{
						SUBREADprintf("Cannot load read part from the tmp file!\n");
						assert(0);
					}


					char * old_line_mem = HashTableGet(first_read_name_table, read_name);
					if(old_line_mem)
						old_line_mem[0]=0xff;
					else
						HashTablePut(first_read_name_table, read_name, new_line_mem);
					//if( first_read_name_table -> numOfElements<4)printf("RV=%s\n", read_name);
				}
			}

			fclose(bbfp);
		}

		//printf("BLK=%d; CKS=%d; READS=%llu\n", x1_block, x1_chunk, first_read_name_table -> numOfElements);
		unsigned long long int finished_second_reads = 0;

		for(x1_chunk = 0; x1_chunk < writer -> current_chunk; x1_chunk++)
		{
			char tmpfname[MAX_FILE_NAME_LENGTH+40];
			sprintf(tmpfname, "%sCHK%08d-BLK%03d.bin", writer -> tmp_path, x1_chunk , x1_block);

	//		printf("START_BLOCK: %s\n", tmpfname);

			FILE * bbfp = f_subr_open(tmpfname,"rb");
			if(!bbfp) continue;

			char * read_line_buf = malloc(3000);
			char * read_name_buf = malloc(MAX_READ_NAME_LEN + MAX_CHROMOSOME_NAME_LEN * 2 + 26);

			while(!feof(bbfp))
			{
				short flags;
				short read_name_len;
				short read_len;
				int ret = fread(&flags, 2,1 , bbfp);
				if(ret<1) break;

				fread(&read_name_len, 2,1 , bbfp);

				if(read_name_len>=MAX_READ_NAME_LEN + MAX_CHROMOSOME_NAME_LEN * 2 + 26)
					SUBREADprintf("VERY_LONG_NAME(%d)\n", read_name_len);
				if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
				{
					fread(read_name_buf, 1, read_name_len, bbfp);
					read_name_buf[read_name_len] = 0;
				}
				else	fseek(bbfp, read_name_len, SEEK_CUR);
				fread(&read_len, 2,1 , bbfp);
				if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
				{
					fread(read_line_buf, 1, read_len, bbfp);
					read_line_buf[read_len] = 0;
				}
				else	fseek(bbfp, read_len, SEEK_CUR);


				if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
				{
//					printf("RRNAME:%s\n", read_name_buf);

					char * first_read_text = HashTableGet(first_read_name_table, read_name_buf);
					strtok(read_name_buf,"\t");
					if(first_read_text && first_read_text[0]!=(char)0xff)
					{
						fputs(read_name_buf, writer->out_fp);
						putc('\t',  writer->out_fp);
						fputs(first_read_text, writer->out_fp);

						fputs(read_name_buf, writer->out_fp);
						putc('\t',  writer->out_fp);
						fputs(read_line_buf, writer->out_fp);

						read_name_buf[strlen(read_name_buf)]='\t';
						HashTableRemove(first_read_name_table, read_name_buf);
						finished_second_reads ++;
					}
					else{

						int dummy_flags = 4 | 1, mate_flags = 0;
						char * dummy_mate_chr = NULL;
						char dummy_mate_chr_buf[120];
						unsigned int dummy_mate_pos = 0, tmpi=0,dummy_char_strpos = 0;
						int tabs = 0;
						int read_cursor = 0;

						for(read_cursor = 0;; read_cursor++)
						{
							char nch = read_line_buf[read_cursor];
							if(!nch) break;
							if(nch == '\t')
							{
								if(tabs == 0){
									mate_flags = tmpi; 
									dummy_mate_chr = read_line_buf+read_cursor+1;
								}
								else if(tabs == 1)
									dummy_char_strpos = read_cursor;
								else if(tabs == 2)
								{
									dummy_mate_pos = tmpi;
									break;
								}
								tmpi=0;
								tabs++;
							}else{
								if(tabs==0 || tabs == 2) tmpi = tmpi * 10 + (nch - '0');
							}
						}


						dummy_flags |= SAM_FLAG_FIRST_READ_IN_PAIR;
						if(mate_flags & SAM_FLAG_UNMAPPED)  dummy_flags |= SAM_FLAG_MATE_UNMATCHED;
						if(mate_flags & SAM_FLAG_REVERSE_STRAND_MATCHED)  dummy_flags |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
						if(mate_flags & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)  dummy_flags |= SAM_FLAG_REVERSE_STRAND_MATCHED;

						memcpy(dummy_mate_chr_buf, dummy_mate_chr, read_line_buf +dummy_char_strpos - dummy_mate_chr);
						dummy_mate_chr_buf[read_line_buf +dummy_char_strpos - dummy_mate_chr]=0;

						char hi_tag_out[18];
						char nh_tag_out[18];

						find_tag_out(read_line_buf, "HI", hi_tag_out);
						find_tag_out(read_line_buf, "NH", nh_tag_out);

						// build a fake FIRST read for the mapped SECOND read.
						// note that the TLEN, MATE_POS and MATE_CHAR are incorrect for general use.
						fprintf(writer->out_fp, "%s\t%d\t*\t0\t0\t*\t%s\t%d\t0\tN\tI%s%s\n", read_name_buf, dummy_flags, dummy_mate_chr_buf, dummy_mate_pos, nh_tag_out, hi_tag_out);
						fputs(read_name_buf, writer->out_fp);
						putc('\t',  writer->out_fp);
						fputs(read_line_buf, writer->out_fp);
						writer -> unpaired_reads +=1;
					}

					//else SUBREADprintf("WARNING: Unpaired read found in file:%s\n", read_name_buf);
				}
			}

			fclose(bbfp);
			unlink(tmpfname);
			free(read_name_buf);
			free(read_line_buf);
		}



		if(1)
		{
			writer -> unpaired_reads += first_read_name_table -> numOfElements;

			KeyValuePair * cursor;
			int bucket;

			// go through the hash table and write correct FIRST lines and dummy SECOND lines.
			for(bucket=0; bucket< first_read_name_table -> numOfBuckets; bucket++)
			{
				cursor = first_read_name_table -> bucketArray[bucket];
				while(1)
				{
					if (!cursor) break;
					char * first_read_text = (char *)cursor -> value;
					char * first_read_name = (char *)cursor -> key;

					if(first_read_text[0]!=(char)0xff)
					{
						int dummy_flags = 4 | 1, mate_flags = 0;
						char * dummy_mate_chr = NULL;
						unsigned int dummy_mate_pos = 0, tmpi=0, dummy_char_strpos = 0;
						int tabs = 0;
						int read_cursor = 0;

						for(read_cursor = 0;; read_cursor++)
						{
							char nch = first_read_text[read_cursor];
							if(!nch) break;
							if(nch == '\t')
							{
								if(tabs == 0){
									mate_flags = tmpi; 
									dummy_mate_chr = first_read_text+read_cursor+1;
								}
								else if(tabs == 1)
									dummy_char_strpos = read_cursor;
								else if(tabs == 2)
								{
									dummy_mate_pos = tmpi;
									break;
								}
								tmpi=0;
								tabs++;
							}else{
								if(tabs==0 || tabs == 2) tmpi = tmpi * 10 + (nch - '0');
							}
						}

						dummy_flags |= SAM_FLAG_SECOND_READ_IN_PAIR;
						if(mate_flags & SAM_FLAG_UNMAPPED)  dummy_flags |= SAM_FLAG_MATE_UNMATCHED;
						if(mate_flags & SAM_FLAG_REVERSE_STRAND_MATCHED)  dummy_flags |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
						if(mate_flags & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)  dummy_flags |= SAM_FLAG_REVERSE_STRAND_MATCHED;

						if((!first_read_text[0])||(!first_read_text[1]))
						{
							SUBREADprintf("unable to recover the first read! : '%s' , flags = %d\n", first_read_name, mate_flags);
							assert(0);
						}

						char nh_tag_out[18];
						char hi_tag_out[18];
						find_tag_out(first_read_text, "NH", nh_tag_out);
						find_tag_out(first_read_text, "HI", hi_tag_out);

						strtok(first_read_name, "\t");
						fputs(first_read_name, writer->out_fp);
						putc('\t',  writer->out_fp);
						fputs(first_read_text, writer->out_fp);
						first_read_text[dummy_char_strpos] = 0;
						fprintf(writer->out_fp, "%s\t%d\t*\t0\t0\t*\t%s\t%d\t0\tN\tI%s%s\n", first_read_name, dummy_flags, dummy_mate_chr, dummy_mate_pos, nh_tag_out,hi_tag_out);
					}
					cursor = cursor->next;
				}
			}


		}

		HashTableDestroy(first_read_name_table);
	}
	fclose(writer -> out_fp);
	signal (SIGTERM, old_sig_TERM);
	signal (SIGINT, old_sig_INT);
}

void sort_SAM_check_chunk(SAM_sort_writer * writer)
{
	if(writer -> current_chunk_size > SAM_SORT_BLOCK_SIZE * SAM_SORT_BLOCKS)
	{
		int xk1;
		for(xk1=0;xk1<SAM_SORT_BLOCKS;xk1++)
		{
			if(writer -> current_block_fp_array[xk1])
				fclose(writer -> current_block_fp_array[xk1]);
		}
		memset(writer -> current_block_fp_array, 0, sizeof(FILE *)*SAM_SORT_BLOCKS);
		writer -> current_chunk_size = 0;
		writer -> current_chunk++;
	}
}

// the SAM_line includes "\n" at the tail!
// line_len = strlen(SAM_line)
int sort_SAM_add_line(SAM_sort_writer * writer, char * SAM_line, int line_len)
{
	assert(writer -> all_chunks_header_fp);
	if(line_len<3) return 0;
	if(SAM_line[0]=='@')
		fputs(SAM_line, writer -> out_fp);
	else
	{
		char read_name[MAX_READ_NAME_LEN + MAX_CHROMOSOME_NAME_LEN * 2 + 26];
		char chromosome_1_name[MAX_CHROMOSOME_NAME_LEN];
		char chromosome_2_name[MAX_CHROMOSOME_NAME_LEN];
		unsigned int pos_1, pos_2;
		int hi_tag,flags = 0, line_cursor = 0, field_cursor = 0, tabs=0;
		char * second_col_pos = NULL;

		chromosome_1_name[0]=0;
		chromosome_2_name[0]=0;
		pos_1 = 0;
		pos_2 = 0;
		hi_tag = -1;

		while(line_cursor < line_len)
		{
			char nch = SAM_line[line_cursor++];
			if(!nch)break;

			if(nch == '\t')
			{
				field_cursor = 0;
				tabs++;
				if(tabs == 1) second_col_pos = SAM_line + line_cursor;
				if(tabs>7) break;
			}
			else if(tabs == 0)
			{
				read_name[field_cursor++] = nch;
				if(MAX_READ_NAME_LEN<=field_cursor){
					return -1;
				}
				read_name[field_cursor] = 0;
			}
			else if(tabs == 1)
				flags = flags*10+(nch-'0');
			else if(tabs == 2)
			{
				chromosome_1_name[field_cursor++] = nch;
				chromosome_1_name[field_cursor]=0;
				if(MAX_CHROMOSOME_NAME_LEN - 1 <= field_cursor) return -1;
			}
			else if(tabs == 3)
				pos_1 = pos_1 * 10 + (nch-'0');
			else if(tabs == 6)
			{
				chromosome_2_name[field_cursor++] = nch;
				chromosome_2_name[field_cursor] = 0;
				if(MAX_CHROMOSOME_NAME_LEN - 1 <= field_cursor) return -1;
			}
			else if(tabs == 7)
				pos_2 = pos_2 * 10 + (nch-'0');

		}
		if(tabs <= 7) return -1;

		//if(memcmp("V0112_0155:7:1101:4561:132881", read_name, 27)==0)

		char * hi_tag_str = strstr(SAM_line,"\tHI:i:");
		if(hi_tag_str)
		{
			hi_tag = 0;
			for(line_cursor=6; ; line_cursor++)
			{
				char nch = hi_tag_str[line_cursor];
				if(!isdigit(nch)) break;
				hi_tag = hi_tag*10 + (nch-'0');
			}
		}

		line_len = strlen(second_col_pos);
		sort_SAM_check_chunk(writer);

		for(field_cursor = 0; read_name[field_cursor] ; field_cursor++)
			if(read_name[field_cursor] == '/') read_name[field_cursor] = 0;

		if(chromosome_2_name[0]=='=')
			strcpy(chromosome_2_name, chromosome_1_name);


		// new read name format: OLD_READ_NAME\tCHR_R1:POS_R1:CHR_R2:POS_R2


		if(flags & SAM_FLAG_MATE_UNMATCHED)
		{
			if(chromosome_2_name[0] != '*')
				strcpy(chromosome_2_name , "*");
			pos_2 = 0;
		}


		if(flags & SAM_FLAG_UNMAPPED)
		{
			if(chromosome_1_name[0] != '*')
				strcpy(chromosome_1_name , "*");
			pos_1 = 0;
		}

		char hi_key [13];
		if(hi_tag >=0)// && pos_1 && pos_2)
			sprintf(hi_key, ":%d", hi_tag);
		else
			hi_key[0]=0;

		if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)
			sprintf(read_name+strlen(read_name), "\t%s:%u:%s:%u%s",chromosome_2_name, pos_2, chromosome_1_name, pos_1, hi_key);
		else
			sprintf(read_name+strlen(read_name), "\t%s:%u:%s:%u%s",chromosome_1_name, pos_1, chromosome_2_name, pos_2, hi_key);

		//if(memcmp("V0112_0155:7:1101:4561:132881", read_name, 27)==0)
		//	printf("RRN=%s\n", read_name);
		
		int read_name_len = strlen(read_name);
		unsigned long long int read_line_hash = sort_SAM_hash(read_name);

		int block_id = read_line_hash % SAM_SORT_BLOCKS;
		if(!writer -> current_block_fp_array[block_id])
		{
			char tmpfname[MAX_FILE_NAME_LENGTH+40];
			sprintf(tmpfname,"%sCHK%08d-BLK%03d.bin", writer -> tmp_path , writer -> current_chunk , block_id);
			writer -> current_block_fp_array[block_id] = f_subr_open(tmpfname, "wb");
		}

		if(line_len < 2)
		{
			SUBREADprintf("unable to put the first read!\n");
			assert(0);
		}

		if(second_col_pos[0]==0 || second_col_pos[1]==0)
		{
			SUBREADprintf("unable to put the first read TEXT!\n");
			assert(0);
		}


//		printf("WRNAME:%s\n", read_name);

		fwrite(&flags, 2, 1, writer -> current_block_fp_array[block_id]);
		fwrite(&read_name_len, 2, 1, writer -> current_block_fp_array[block_id]);
		fwrite(read_name, 1, read_name_len, writer -> current_block_fp_array[block_id]);
		fwrite(&line_len, 2, 1, writer -> current_block_fp_array[block_id]);
		fwrite(second_col_pos, 1, line_len, writer -> current_block_fp_array[block_id]);

		writer -> output_file_size += line_len;
		writer -> current_chunk_size += line_len;
		writer -> written_reads ++;
	}

	return 0;
}

int is_SAM_unsorted(char * SAM_line, char * tmp_read_name, short * tmp_flag, unsigned long long int read_no)
{
	char read_name[MAX_READ_NAME_LEN];
	int flags = 0, line_cursor = 0, field_cursor = 0, tabs=0;

	while(1)
	{
		char nch = SAM_line[line_cursor++];
		if(!nch)break;
		if(nch == '\t')
		{
			field_cursor = 0;
			tabs++;
			if(tabs>1) break;
		}
		else if(tabs == 0)
		{
			read_name[field_cursor++] = nch;
			assert(MAX_READ_NAME_LEN>field_cursor);
			read_name[field_cursor] = 0;
		}
		else if(tabs == 1)
			flags = flags*10+(nch-'0');
	}

		//int is_second_read = (flags & 0x80) ? 1:0;
	for(field_cursor = 0; read_name[field_cursor] ; field_cursor++)
		if(read_name[field_cursor] == '/') read_name[field_cursor] = 0;


	(*tmp_flag) = flags;
	if(!(flags &1)) return 0;
	if(read_no % 2 == 0)
	{
		if(flags & SAM_FLAG_SECOND_READ_IN_PAIR)return 1;
		strcpy(tmp_read_name , read_name);
	}
	else
	{
		if(flags & SAM_FLAG_FIRST_READ_IN_PAIR) return 1;
		if(strcmp(tmp_read_name, read_name))return 1;
	}

	return 0;
}

int is_certainly_bam_file(char * fname, int * is_first_read_PE)
{

	int read_type = probe_file_type(fname, is_first_read_PE);
	if(read_type == FILE_TYPE_NONEXIST || read_type == FILE_TYPE_EMPTY || read_type == FILE_TYPE_UNKNOWN)
		return -1;
	if(read_type == FILE_TYPE_BAM)
		return 1;
	return 0;
}


int is_pipe_file(char * fname)
{
	FILE * fp = fopen(fname,"r");
	if(!fp) return 0;

	int seeked = fseek(fp, 0, SEEK_SET);
	fclose(fp);

	return (seeked != 0);
}

int warning_file_type(char * fname, int expected_type)
{
	int ret_pipe_file = is_pipe_file(fname);
	if(ret_pipe_file)
	{
		print_in_box(80,0,0,"WARNING file '%s' is not a regular file.", fname);
		print_in_box(80,0,0,"        No alignment can be done on a pipe file.");
		print_in_box(80,0,0,"        If the FASTQ file is gzipped, please use gzFASTQinput option.");
		print_in_box(80,0,0,"");
		return 1;
	}

	int read_type = probe_file_type(fname, NULL);

	if(read_type == FILE_TYPE_NONEXIST)
	{
		SUBREADprintf("ERROR: unable to open file '%s'. File name might be incorrect, or you do not have the permission to read the file.\n", fname);
		return -1;
	}
	else if(read_type == FILE_TYPE_EMPTY)
	{
		print_in_box(80,0,0,"WARNING file '%s' is empty.", fname);
		return 1;
	}

	else if((expected_type == FILE_TYPE_FAST_ && (read_type!= FILE_TYPE_FASTQ && read_type!= FILE_TYPE_FASTA && read_type!= FILE_TYPE_GZIP_FASTQ))||
		(expected_type == FILE_TYPE_GZIP_FAST_ && read_type!= FILE_TYPE_GZIP_FASTA) ||
		((  expected_type != FILE_TYPE_GZIP_FAST_ && expected_type != FILE_TYPE_FAST_) && expected_type != read_type))
	{
		char * req_fmt = "SAM";
		if(expected_type==FILE_TYPE_BAM) req_fmt = "BAM";
		else if(expected_type==FILE_TYPE_FAST_) req_fmt = "FASTQ or FASTA";
		else if(expected_type==FILE_TYPE_GZIP_FAST_) req_fmt = "gzip FASTQ or FASTA";

		char * real_fmt = "SAM";
		if(read_type==FILE_TYPE_BAM) real_fmt = "BAM";
		else if(read_type==FILE_TYPE_FASTA) real_fmt = "FASTA";
		else if(read_type==FILE_TYPE_FASTQ) real_fmt = "FASTQ";
		else if(read_type==FILE_TYPE_GZIP_FASTQ) real_fmt = "gzip FASTQ";
		else if(read_type==FILE_TYPE_GZIP_FASTA) real_fmt = "gzip FASTA";

		print_in_box(80,0,0,"WARNING format issue in file '%s':", fname);
		print_in_box(80,0,0,"        The required format is : %s", req_fmt); 
		if(read_type == FILE_TYPE_UNKNOWN)
			print_in_box(80,0,0,"        The file format is unknown.");
		else
			print_in_box(80,0,0,"        The real format seems to be : %s", real_fmt);
		print_in_box(80,0,0,"A wrong format may result in wrong results or crash the program.");
		print_in_box(80,0,0,"Please refer to the manual for file format options.");
		print_in_box(80,0,0,"If the file is in the correct format, please ignore this message.");
		print_in_box(80,0,0,"");

		return 1;
	}
	return 0;
}

char * gzgets_noempty(void * fp, char * buf, int maxlen)
{
	char * ret;
	while(1)
	{
		ret = gzgets(fp,buf, maxlen);
		if(!ret)return NULL;
		if(ret[0]!='\n') return ret;
	}
}


char * fgets_noempty(char * buf, int maxlen, FILE * fp)
{
	char * ret;
	while(1)
	{
		ret = fgets(buf, maxlen, fp);
		if(!ret)return NULL;
		if(ret[0]!='\n') return ret;
	}
}


int probe_file_type_fast(char * fname){
	FILE * fp = f_subr_open(fname, "rb");
	if(!fp) return FILE_TYPE_NONEXIST;

	int ret = FILE_TYPE_UNKNOWN; 
	int nch;
	char *test_buf=malloc(5000);

	nch = fgetc(fp);

	if(feof(fp))
		ret = FILE_TYPE_EMPTY;
	else
	{
		if(nch == '@')	// FASTQ OR SAM
		{
			char * rptr = fgets_noempty(test_buf, 4999, fp);
			int second_line_len = 0;
			if(rptr)
			{
				rptr = fgets_noempty(test_buf, 4999, fp);
				if(rptr)
				{
					second_line_len = strlen(test_buf);
					int tabs = 0, x1;
					for(x1=0;x1<4999;x1++)
					{
						if(test_buf[x1]=='\n' || !test_buf[x1]) break;
						if(test_buf[x1]=='\t'){
							tabs++;
							continue;
						}

						if(tabs == 1)
							if(!isdigit(test_buf[x1]))break;
					}
					if(rptr[0]=='@' || tabs>7)
						ret = FILE_TYPE_SAM;
				}
			}
			if(ret == FILE_TYPE_UNKNOWN)
			{
				rptr = fgets_noempty(test_buf, 4999, fp);
				if(rptr[0] == '+')
				{
					rptr = fgets_noempty(test_buf, 4999, fp);
					if(rptr && second_line_len == strlen(test_buf))
						ret = FILE_TYPE_FASTQ;
				}
			}
		}
		else if(nch == '>') // FASTA
		{
			ret = FILE_TYPE_FASTA;
		}
		else if(nch == 31) // BAM OR GZ_FASTQ
		{
			nch = fgetc(fp);
			if(nch == 139)
			{
				fclose(fp);
				fp=NULL;
				gzFile zfp = gzopen(fname, "rb");
				if(zfp)
				{
					int rlen = gzread(zfp, test_buf,4);
					if(rlen == 4 && memcmp(test_buf,"BAM\1",4)==0)
						ret = FILE_TYPE_BAM;
					if(rlen == 4 && test_buf[0]=='@')
						ret = FILE_TYPE_GZIP_FASTQ;
					if(rlen == 4 && test_buf[0]=='>')
						ret = FILE_TYPE_GZIP_FASTA;
					gzclose(zfp);
				}
			}
		}
		else if(nch >= 0x20 && nch <= 0x7f) // SAM without headers
		{
			int tabs = 0, x1;
			char * rptr = fgets(test_buf, 4999, fp);
			if(rptr)
				for(x1=0;x1<4999;x1++)
				{
					if(test_buf[x1]=='\n' || !test_buf[x1]) break;
					if(test_buf[x1]=='\t'){
						tabs++;
						continue;
					}
					if(tabs == 1)
						if(!isdigit(test_buf[x1]))break;
				}
			if(tabs>7)
				ret = FILE_TYPE_SAM;

		}
	}

	if(fp)fclose(fp);

	free(test_buf);
	return ret;

}
int probe_file_type(char * fname, int * is_first_read_PE)
{
	FILE * fp = f_subr_open(fname, "rb");
	if(!fp) return FILE_TYPE_NONEXIST;

	int ret = FILE_TYPE_UNKNOWN; 
	int nch;
	char *test_buf=malloc(5000);

	nch = fgetc(fp);

	if(feof(fp))
		ret = FILE_TYPE_EMPTY;
	
	else
	{
		if(nch == '@')	// FASTQ OR SAM
		{
			char * rptr = fgets_noempty(test_buf, 4999, fp);
			int second_line_len = 0;
			if(rptr)
			{
				rptr = fgets_noempty(test_buf, 4999, fp);
				if(rptr)
				{
					second_line_len = strlen(test_buf);
					int tabs = 0, x1;
					for(x1=0;x1<4999;x1++)
					{
						if(test_buf[x1]=='\n' || !test_buf[x1]) break;
						if(test_buf[x1]=='\t'){
							tabs++;
							continue;
						}

						if(tabs == 1)
							if(!isdigit(test_buf[x1]))break;
					}
					if(rptr[0]=='@' || tabs>7)
						ret = FILE_TYPE_SAM;
				}
			}
			if(ret == FILE_TYPE_UNKNOWN)
			{
				rptr = fgets_noempty(test_buf, 4999, fp);
				if(rptr[0] == '+')
				{
					rptr = fgets_noempty(test_buf, 4999, fp);
					if(rptr && second_line_len == strlen(test_buf))
						ret = FILE_TYPE_FASTQ;
				}
			}
		}
		else if(nch == '>') // FASTA
		{
			char * rptr = fgets(test_buf, 4999, fp);
			int x1;
			if(rptr)
			{
				ret = FILE_TYPE_FASTA;
				for(x1=0;x1<4999;x1++)
				{
					if(test_buf[x1]=='\n' || !test_buf[x1]) break;
					nch = toupper(test_buf[x1]);
					if(nch < ' ' || nch>127)
					{
						ret = FILE_TYPE_UNKNOWN;
						break;
					}
				}
				rptr = fgets(test_buf, 4999, fp);
				if(rptr && ret == FILE_TYPE_FASTA)
				{
					for(x1=0;x1<4999;x1++)
					{
						if(test_buf[x1]=='\n' || !test_buf[x1]) break;
						nch = toupper(test_buf[x1]);
						if(nch == 'A' || nch == 'T' || nch == 'G' || nch == 'C' || nch == 'N' || nch == '.' || (nch >='0' && nch <= '3'))
							;
						else
						{
							ret = FILE_TYPE_UNKNOWN;
							break;
						}
					}

					if(x1==0) ret = FILE_TYPE_UNKNOWN;
				}
			}
		}
		else if(nch == 31) // BAM OR GZ_FASTQ
		{
			nch = fgetc(fp);
			if(nch == 139)
			{
				fclose(fp);
				fp=NULL;
				gzFile zfp = gzopen(fname, "rb");
				if(zfp)
				{
					int rlen = gzread(zfp, test_buf,4);
					if(rlen == 4 && memcmp(test_buf,"BAM\1",4)==0)
						ret = FILE_TYPE_BAM;
					if(rlen == 4 && test_buf[0]=='@')
						ret = FILE_TYPE_GZIP_FASTQ;
					if(rlen == 4 && test_buf[0]=='>')
						ret = FILE_TYPE_GZIP_FASTA;
					gzclose(zfp);
				}
			}
		}
		else if(nch >= 0x20 && nch <= 0x7f) // SAM without headers
		{
			int tabs = 0, x1;
			char * rptr = fgets(test_buf, 4999, fp);
			if(rptr)
				for(x1=0;x1<4999;x1++)
				{
					if(test_buf[x1]=='\n' || !test_buf[x1]) break;
					if(test_buf[x1]=='\t'){
						tabs++;
						continue;
					}
					if(tabs == 1)
						if(!isdigit(test_buf[x1]))break;
				}
			if(tabs>7)
				ret = FILE_TYPE_SAM;

		}
	}

	if(fp)fclose(fp);

	if(FILE_TYPE_BAM == ret || FILE_TYPE_SAM == ret)
		if(is_first_read_PE)
		{
			SamBam_FILE * tpfp = SamBam_fopen(fname, (FILE_TYPE_BAM  == ret)?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);
			while(1)
			{
				char * tbr = SamBam_fgets(tpfp, test_buf, 4999, 0);
				if(!tbr){
					ret = FILE_TYPE_EMPTY;
					break;
				}
				if(tbr[0]=='\n' || tbr[0]=='\r'){
					ret = FILE_TYPE_UNKNOWN;
					break;
				}
				if(tbr[0]=='@') continue;
				char * rname_str, *tmpstr;
				rname_str = strtok_r(tbr, "\t", &tmpstr);
				if(!rname_str)
				{
					ret = FILE_TYPE_UNKNOWN;
					break;
				}
				rname_str = strtok_r(NULL, "\t", &tmpstr);
				if((!rname_str)|| (!isdigit(rname_str[0]))) 
				{
					ret = FILE_TYPE_UNKNOWN;
					break;
				}

				int flags = atoi(rname_str);
				(*is_first_read_PE) = flags &1;
				break;
			}
			SamBam_fclose(tpfp);
		}

	free(test_buf);
	return ret;
}

#ifdef MAKE_INPUTTEST
int main(int argc, char ** argv)
{
	FILE * ifp;
	unsigned long long int rno=0;
	short tmp_flags, is_sorted = 1;
	char buff[3000], tmp_rname[100];

	ifp = f_subr_open(argv[1],"r");
	while(1)
	{
		char * rr = fgets(buff,2999, ifp);
		if(!rr) break;
		if(buff[0]=='@')continue;
		if(is_SAM_unsorted(buff, tmp_rname, &tmp_flags, rno))
		{
			printf("The input file is unsorted.\n");
			is_sorted = 0;
			break;
		}
		rno++;
	}
	
	fclose(ifp);

	//if(is_sorted) return 0;

	ifp = f_subr_open(argv[1],"r");
	SAM_sort_writer writer;
	if(sort_SAM_create(&writer, argv[2], ".")){
		printf("ERROR: unable to create the writer!\n");
		return -1;
	}

	while(1)
	{
		char * rr = fgets(buff,2999, ifp);
		if(!rr) break;
		int line_len = strlen(buff);
		sort_SAM_add_line(&writer, buff, line_len);
	}
	fclose(ifp);
	sort_SAM_finalise(&writer);
	printf("WRITTEN=%llu\nUNPAIR=%llu\n", writer.written_reads, writer.unpaired_reads);
}
#endif
#ifdef MAKE_TYPETEST


int main(int argc, char ** argv)
{
	char * fn = argv[1];
	int type = probe_file_type(fn, NULL);
	switch(type)
	{
		case FILE_TYPE_FASTQ: printf("Type: FASTQ\n"); break;
		case FILE_TYPE_FASTA: printf("Type: FASTA\n"); break;
		case FILE_TYPE_SAM  : printf("Type: SAM\n"); break;
		case FILE_TYPE_BAM  : printf("Type: BAM\n"); break;
		default: printf("Unknown type!\n");
	}
}

#endif
