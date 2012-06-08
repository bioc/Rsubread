#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>
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
	if(must_upper)
	{
		while(!feof(fp))
		{
			char ch = fgetc(fp);
			if (ch == '\n')
			{
				if (started)break;
				else continue;
			}
			else
				started = 1;
			if(ret <max_read_len && ch != '\r')
				if ((ch!=' ' && ch != '\t'))
					buff[ret++] = toupper(ch);
		}
	}
	else
	{
		while(!feof(fp))
		{
			char ch = fgetc(fp);
			if (ch == '\n')
			{
				if (started)break;
				else continue;
			}
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

#define SKIP_LINE { nch=' '; while(nch != EOF && nch != '\n') nch = fgetc(input->input_fp); }

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
		//if(input->file_type != GENE_INPUT_SAM_SINGLE)
		//	SKIP_LINE
	}
	else if(input->file_type == GENE_INPUT_FASTA)
	{
		SKIP_LINE
		SKIP_LINE
	}
	else if(input->file_type == GENE_INPUT_FASTQ)
	{
		SKIP_LINE
		SKIP_LINE
		SKIP_LINE
		SKIP_LINE
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
		//if(input->file_type != GENE_INPUT_SAM_SINGLE)
			// skip a line if not single-end
		//	read_line(1, input->input_fp, in_buff, 0);
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


//This function returns 0 if the line is a mapped read; -1 if the line is in a wrong format and 1 if the read is unmapped.
int parse_SAM_line(char * sam_line, char * read_name, int * flags, char * chro, unsigned int * pos, char * cigar, int * mapping_quality, char * sequence , char * quality_string, int * rl)
{
	char cc;
	int ci = 0, k=0, field=0, ret_quality = 0, ret_flag = 0;
	unsigned int ret_pos = 0;
	
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
		else if(field == 2)
		{
			if(ci == 0 && cc == '*') return 1;
			chro[ci++] = cc;
		}
		else if(field == 3)
			ret_pos = ret_pos * 10 + (cc-'0');
		else if(field == 4)
			ret_quality = ret_quality * 10 + (cc-'0');
		else if(field == 5)
			cigar[ci++] = cc;
		k++;

	}

	if(field == 10 && ci>0)quality_string[ci]=0;
	else if(field < 10) return -1;
	
	(*mapping_quality) = ret_quality;
	(*pos) = ret_pos;
	(*flags) = ret_flag;
	return 0;
	
}


// This function returns 0 if the block is determined.
// The block is undeterminable if the chromosome name is not in known_chromosomes, or the position is larger than the known length.
int get_read_block(char *chro, unsigned int pos, char *temp_file_suffix, chromosome_t *known_chromosomes)
{
	int chro_no;

	for(chro_no=0;known_chromosomes[chro_no].chromosome_name[0]; chro_no++)
	{
		if(strcmp(chro , known_chromosomes[chro_no].chromosome_name) == 0)
			break;
	}
	if(!known_chromosomes[chro_no].chromosome_name[0]) return 1;
	if(pos >= known_chromosomes[chro_no].known_length) return 1;

	int block_no = pos / BASE_BLOCK_LENGTH;
	sprintf(temp_file_suffix , "%s-%04u.bin", chro, block_no);
	return 0;
}

FILE * get_temp_file_pointer(char *temp_file_name, HashTable* fp_table)
{
	FILE * temp_file_pointer = (FILE *) HashTableGet(fp_table, temp_file_name);
	if(!temp_file_pointer)
	{
		temp_file_pointer = fopen(temp_file_name,"w");
		assert(temp_file_pointer);

		HashTablePut(fp_table, temp_file_name ,temp_file_pointer);
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

void write_read_block_file(FILE *temp_fp , unsigned int read_number, char *read_name, int flags, char * chro, unsigned int pos, char *cigar, int mapping_quality, char *sequence , char *quality_string, int rl , int is_sequence_needed)
{
	if(is_sequence_needed)
	{
		assert(0);
	}
	else
	{
		base_block_temp_read_t datum;
		datum.read_number = read_number;
		datum.pos = pos;

		fwrite(&datum, sizeof(datum), 1, temp_fp);
	}
}

int break_SAM_file(char * in_SAM_file, char * temp_file_prefix, unsigned int * real_read_count, chromosome_t * known_chromosomes, int is_sequence_needed)
{
	FILE * fp = fopen(in_SAM_file,"r");
	int i;
	HashTable * fp_table;
	unsigned int read_number = 0;

	if(!fp){
		printf("SAM file does not exist or is not accessible: '%s'\n", in_SAM_file);
		return 1;
	}


	fp_table = HashTableCreate( OFFSET_TABLE_SIZE / 16 );
	HashTableSetDeallocationFunctions(fp_table, free, my_fclose);
	HashTableSetKeyComparisonFunction(fp_table, my_strcmp);

	while(!feof(fp))
	{
		char line_buffer [3000];
		int linelen = read_line(2999, fp, line_buffer, 0);

		if(line_buffer[0]=='@')
		{
			int chro_numb=0, field = 0, ci=0, ciw = 0;
			while(known_chromosomes[chro_numb].chromosome_name[0]!=0) chro_numb++;
			known_chromosomes[chro_numb].known_length = 0;
			for(i=0; i< linelen; i++)
			{
				char cc = line_buffer[i];

				if(cc == '\r' || cc=='\n') continue;

				if(cc == '\t')
				{
					if(field == 1)
						known_chromosomes[chro_numb].chromosome_name[ci]=0;
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

		else if((line_buffer[0] >='A' && line_buffer[0]<='Z') || (line_buffer[0] >='a' && line_buffer[0] <='z') || (line_buffer[0] >='0' && line_buffer[0] <='9') || line_buffer[0] =='_' || line_buffer[0] =='.')
		{
			char read_name[MAX_READ_NAME_LEN], chro[MAX_CHROMOSOME_NAME_LEN], cigar[EXON_MAX_CIGAR_LEN], sequence[MAX_READ_LENGTH+1], quality_string[MAX_READ_LENGTH+1];
			int flags = 0, mapping_quality = 0, rl=0;
			unsigned int pos = 0;
			char temp_file_suffix[MAX_CHROMOSOME_NAME_LEN+20];
			char temp_file_name[MAX_CHROMOSOME_NAME_LEN+20+300];
			FILE * temp_fp;

			int line_parse_result = parse_SAM_line(line_buffer, read_name, &flags, chro, &pos, cigar, & mapping_quality, sequence , quality_string, &rl);
			printf("R#%d P=%d L=%s\n", read_number, line_parse_result, line_buffer);
			if(line_parse_result)
			{
				read_number ++;
				continue;
			}

			// if the read block is not determinable, do nothing.
			if(get_read_block(chro, pos, temp_file_suffix, known_chromosomes))
			{
				read_number ++;
				continue;
			}

			sprintf(temp_file_name, "%s%s", temp_file_prefix , temp_file_suffix);

			temp_fp = get_temp_file_pointer(temp_file_name, fp_table);
			assert(temp_fp);

			write_read_block_file(temp_fp , read_number, read_name, flags, chro, pos, cigar, mapping_quality, sequence , quality_string, rl , is_sequence_needed);
			read_number ++;
		}
	}

	HashTableDestroy(fp_table);
	fclose(fp);
	(*real_read_count) = read_number;
	return 0;
}
