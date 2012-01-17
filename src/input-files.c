#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include "input-files.h"


void fastq_64_to_33(char * qs)
{
	int i=0;
	while(qs[i])
		qs[i++] -= 31;
}

double guess_reads_density(char * fname, int is_sam)
{
	gene_input_t ginp;
	long long int fpos =0, fpos2 = 0;
	int i;
	char buff[1200] ;

	if(is_sam == 0)
	{
		if(geinput_open(fname, &ginp))return -1.0;
	}else if(is_sam == 1)
	{
		if(geinput_open_sam(fname, &ginp,0))return -1.0;
	}else if(is_sam == 2)
	{
		if(geinput_open_sam(fname, &ginp,1))return -1.0;
	}

	geinput_next_read(&ginp, NULL, buff, NULL);

	fpos = ftello(ginp.input_fp);

	for(i=0; i<1000; i++)
	{
		if(geinput_next_read(&ginp, NULL, buff, NULL)<0) break;
	}
	fpos2 = ftello(ginp.input_fp) - fpos;
	geinput_close(&ginp);

	//printf("T0=%llu T1=%llu DENS=%.5f\n", fpos, fpos2+ fpos, fpos2*1.0/i);
	return fpos2*1.0/i;
}

int is_gene_char(char c)
{
	//if(c== 'M' || c == 'm' || c == 'U' || c == 'u' || c == 'A' || c=='a' || c=='G' || c=='g' || c=='C' || c=='c' || c=='T' || c=='t' || c=='N' || c=='n')
	if((c>='A' && c<'Z') || (c>='a' && c<='z'))
		return GENE_SPACE_BASE;
	if(c>='0' && c<'9')
		return GENE_SPACE_COLOR;
	if(c=='N' || c == '.')
		return GENE_SPACE_BASE;
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
			return -i-1;

		ret += statbuf.st_size;
		ret -= 150;
		if(ret<0)ret=0;
	}
	return ret * 70 / 71;
}

int read_line(int max_read_len, FILE * fp, char * buff, int must_upper)
{
	int ret =0;
	int started = 0;
	while(!feof(fp))
	{
		char ch = fgetc(fp);

/*		if(ret == 0 && ch == '#')
			while(!feof(fp))
				if(fgetc(fp) == '\n')break;
*/
		if (ch == '\n')
		{
			if (started)break;
			else continue;
		}

		started = 1;
		if(ret <max_read_len && ch != '\r')
			if ((ch!=' ' && ch != '\t') || !must_upper)
				buff[ret++] = must_upper?toupper(ch):ch;
	}
	buff[ret]=0;
	return ret;
}

int geinput_readline(gene_input_t * input, char * buff, int conv_to_upper)
{
	return read_line(1200, input -> input_fp, buff, conv_to_upper);
}

int is_read(char * in_buff)
{
	int p=0;
	char c;
	int space_type = GENE_SPACE_BASE;
	while((c=in_buff[p++])!='\0')
	{
		int x = is_gene_char(c);
		if (x == GENE_SPACE_COLOR)
			space_type = GENE_SPACE_COLOR;
		else if(!x) 
			return 0;
	}
	return space_type;
}

int geinput_open_sam(const char * filename, gene_input_t * input, int half_number)
{
	input->input_fp = fopen(filename, "r");

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
			break;
		}
	}	

	return 0;	
}

int geinput_open(const char * filename, gene_input_t * input)
{
	char in_buff[1201];
	int line_no = 0;
	if(strlen(filename)>298)
		return 1;

	strcpy(input->filename, filename);
	input->input_fp = fopen(filename, "r");

	if(input->input_fp == NULL)	
		return 1;

	while (1){
		int rlen = read_line(1200, input->input_fp, in_buff, 0);
		if (rlen<=0)
			return 1;

		if(line_no==0 && is_read(in_buff))
		{
			input->file_type = GENE_INPUT_PLAIN;
			input->space_type = is_read(in_buff);
			fseek(input->input_fp,0,SEEK_SET);
			break;
		}
		if(in_buff[0]=='>')
		{
			input->file_type = GENE_INPUT_FASTA;
			rlen += read_line(1200, input->input_fp, in_buff, 0);
			input->space_type = is_read(in_buff);

			fseek(input->input_fp,-rlen-2,SEEK_CUR);
			break;
		}
		if(in_buff[0]=='@')
		{
			input->file_type = GENE_INPUT_FASTQ;

			rlen += read_line(1200, input->input_fp, in_buff, 0);
			input->space_type = is_read(in_buff);

			fseek(input->input_fp,-rlen-2,SEEK_CUR);
			break;
		}		
		line_no++;
	}

	return 0;
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
			else if (nch < 0 || nch > 126)printf("\nUnrecognised char = #%d\n", nch);

			if (nch == '\r' || nch == '\n')
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
				printf ("\nUnknown character in the chromosome data: %d, ignored!\n", nch);
				return 'N';
			}		
			last_br = 0;
		}
	}
	else
	{
		printf("Only the FASTA format is accepted for input chromosome data.\n");
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

int geinput_next_read(gene_input_t * input, char * read_name, char * read_string, char * quality_string)
{
	if(input->file_type == GENE_INPUT_PLAIN)
	{
		int ret = read_line(1200, input->input_fp, read_string, 0);
		if(quality_string) *quality_string=0;

		if(ret <3)return -1;
		return ret;
	}
	else if(input->file_type >= GENE_INPUT_SAM_SINGLE)
	{
		char in_buff [3001];
		int tabs = 0;
		int current_str_pos = 0;
		int i;
		int ret = -1;
		int linelen = read_line(3000, input->input_fp, in_buff, 0);
		int need_reverse = 0;
		char mask_buf[5];
		if(linelen <1)return -1;
		if(read_name)
			*read_name = 0;
		if(quality_string)
			*quality_string = 0;
		*read_string = 0;

		for(i=0; i<linelen+1; i++)
		{
			if(in_buff[i]=='\t'|| i ==linelen)
			{
				if(tabs == 0 && read_name)read_name[current_str_pos] = 0;
				if(tabs == 1)
				{
					mask_buf[current_str_pos] = 0;
					need_reverse = (atoi(mask_buf) & 16 )?1:0;
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
		if(need_reverse)
		{
			if(quality_string)
				reverse_quality(quality_string, ret);
			reverse_read(read_string, ret, input->space_type);
		}
		if(input->file_type != GENE_INPUT_SAM_SINGLE)
			// skip a line if not single-end
			read_line(1, input->input_fp, in_buff, 0);
		return ret;
	}
	else if(input->file_type == GENE_INPUT_FASTA)
	{
		int ret;
		while(1)
		{
			ret = read_line(1200, input->input_fp, read_string, 0);
			if(ret <1)return -1;
			if(read_string[0]=='>'){
				if (read_name != NULL)
					strncpy(read_name, read_string+1, 100);
				break;
			}
		}
		ret = 0;
		while(1)
		{
			char nch;
			ret += read_line(1200, input->input_fp, read_string+ret, 1);

			while(1){
				nch = fgetc(input->input_fp);
				if (nch!='\n')break;
			}
			fseek(input->input_fp, -1, SEEK_CUR);
			if (nch<1)
				break;
			if(nch =='>') 
				break;
		}

		if(quality_string) (*quality_string)=0;

		if(ret <1)return -1;
		return ret;
		
	}
	else if(input->file_type == GENE_INPUT_FASTQ)
	{
		int ret;
		while(1)
		{
			ret = read_line(1200, input->input_fp, read_string, 0);
			if(ret <1)return -1;
			if(read_string[0]=='@'){
				if (read_name != NULL)
					strncpy(read_name, read_string+1,100);
				break;
			}
		}
		ret = 0;
		while(1)
		{
			char nch;
			ret += read_line(1200, input->input_fp, read_string+ret, 1);


			// test if the next line is the continued line of this read
			// read one more line, see if it is '+'
			nch = fgetc(input->input_fp);
			if (nch<1)
				break;
			fseek(input->input_fp, -1, SEEK_CUR);
			if(!is_gene_char(nch)) 
				break;
		}
		if(1)
		{
			int qret = 0;
			char nch;
			char fake_q_string [1200];
			if (!quality_string) quality_string = fake_q_string;
			// skip the line starting with '+'
			read_line(1200, input->input_fp, quality_string+qret, 0);

			while(1)
			{
				qret += read_line(1200, input->input_fp, quality_string+qret, 0);

				// test if the next line is the continued line of this read
				nch = fgetc(input->input_fp);
				if (nch<1)
					 break;
				fseek(input->input_fp, -1, SEEK_CUR);
				if(nch =='@')
					break;
			}
		}

		if(ret <1)return -1;
		return ret;
		
	}else return -1;
}
void geinput_close(gene_input_t * input)
{
	fclose(input->input_fp);
}


void reverse_read(char * InBuff, int read_len, int space_type)
{
	int i;

	if(space_type == GENE_SPACE_COLOR)
	{
		for (i=0; i<read_len/2; i++)
		{
			int rll1 = read_len - 1 - i;
			char tmp = InBuff[rll1];
			InBuff[rll1] = InBuff[i];
			InBuff[i] = tmp;
		}
	}
	else
	{
		for (i=0; i<read_len/2; i++)
		{
			int rll1 = read_len - 1 - i;
			char tmp = InBuff[rll1];

			if(InBuff[i]=='A')InBuff[rll1]='T';
			else if(InBuff[i]=='G')InBuff[rll1]='C';
			else if(InBuff[i]=='C')InBuff[rll1]='G';
			else if(InBuff[i]=='T' || InBuff[i]=='U')InBuff[rll1]='A';

			if(tmp=='A')InBuff[i]='T';
			else if(tmp=='G')InBuff[i]='C';
			else if(tmp=='C')InBuff[i]='G';
			else if(tmp=='T' || tmp=='U')InBuff[i]='A';
		}
		if(i*2 == read_len-1)
		{
			if(InBuff[i]=='A')InBuff[i]='T';
			else if(InBuff[i]=='G')InBuff[i]='C';
			else if(InBuff[i]=='C')InBuff[i]='G';
			else if(InBuff[i]=='T' || InBuff[i]=='U')InBuff[i]='A';
		}
	}

}



void reverse_quality(char * InBuff, int read_len)
{
	int i;
	if(!InBuff[0]) return;
	for (i=0; i<read_len/2; i++)
	{
		char tmp;
		tmp = InBuff[i];
		InBuff[i] = InBuff[read_len -1-i];
		InBuff[read_len -1-i] = tmp;		
	}
}

int genekey2int(char key [],int space_type)
{
	int i;
	int ret;

	ret = 0;
	if(space_type == GENE_SPACE_BASE)
		for (i=0; i<16; i++)
		{
			ret = ret << 2;
			ret += base2int (key[i]);
		}
	else
		for (i=0; i<16; i++)
		{
			ret = ret << 2;
			ret += color2int (key[i]);
		}
	
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
	float step = max(3.00001, (len-16-GENE_SLIDING_STEP)*1.0/(TOTAL_SUBREADS-1)+0.00001);
	return (int) (step * subread) + 15;
	//return (int)((1.*len-16.)/TOTAL_SUBREADS * subread+15);
}

