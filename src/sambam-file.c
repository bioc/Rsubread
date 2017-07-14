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

/***************************************************************

	The SamBam_reg2bin function was derived from the BAM
    specification. (The SAM Format Specication Working
    Group, September 7, 2011)

  ***************************************************************/
  
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "subread.h"
#include "core.h"
#include "gene-algorithms.h"
#include "sambam-file.h"

int SamBam_fetch_next_chunk(SamBam_FILE *fp)
{
	int x1, room =  SAMBAM_INPUT_STREAM_SIZE - ( fp -> input_binary_stream_write_ptr - fp -> input_binary_stream_read_ptr); 

	if(room < 65536)
		return -1;


	for(x1=0; x1 < fp->input_binary_stream_write_ptr - fp -> input_binary_stream_read_ptr; x1 ++)
	{
		fp -> input_binary_stream_buffer [x1] = fp -> input_binary_stream_buffer [x1 + fp->input_binary_stream_read_ptr - fp -> input_binary_stream_buffer_start_ptr];
	}
	fp -> input_binary_stream_buffer_start_ptr = fp->input_binary_stream_read_ptr;

	char * in_buff = malloc(65537);
	unsigned int real_len = 0;
	int ret, have = 0;
	
	while (1){
			int nchunk = 0;
			ret = PBam_get_next_zchunk(fp -> os_file, in_buff, 65536, & real_len);
			if(ret > 0)
				nchunk = SamBam_unzip(fp -> input_binary_stream_buffer + fp->input_binary_stream_write_ptr - fp -> input_binary_stream_read_ptr + have , in_buff , ret);
			else if(ret == -2){
				SUBREADputs("ERROR: BAM format is broken!");
				return -2;
			}

			//printf("RET=%d; CHK=%d\n", ret, nchunk);

			if(nchunk>0)
				have += nchunk; 
			if(have > 3000) break;
			if(feof(fp->os_file)){
				fp->is_eof=1;
				break;
			}
	}
	free(in_buff);

	fp -> input_binary_stream_write_ptr += have;

	return have;

}

SamBam_FILE * SamBam_fopen(char * fname , int file_type)
{
	SamBam_FILE * ret = (SamBam_FILE *)malloc(sizeof(SamBam_FILE));
	memset(ret, 0, sizeof(SamBam_FILE));
	ret -> file_type = file_type;

	if(file_type ==SAMBAM_FILE_SAM) 
	{
		ret -> os_file = f_subr_open(fname, "rb");
		if(!ret -> os_file)
		{
			free(ret);
			return NULL;
		}
		fseek(ret -> os_file,0,SEEK_SET);
	}
	else
	{
		ret -> os_file = f_subr_open(fname, "rb");
		if(ret -> os_file == NULL)
		{
			free(ret);
			return NULL;
		}
		unsigned char first_ch = fgetc(ret->os_file);
		unsigned char second_ch = fgetc(ret->os_file);

		if(first_ch!=31 || second_ch!=139)
		{
			free(ret);
			return NULL;
		}

		fseek(ret->os_file, 0, SEEK_SET);

		ret -> input_binary_stream_buffer = (char *)malloc(SAMBAM_INPUT_STREAM_SIZE);
		ret -> input_binary_stream_read_ptr = 0;
		ret -> input_binary_stream_write_ptr = 0;
		ret -> input_binary_stream_buffer_start_ptr = 0;

		ret -> bam_file_stage = BAM_FILE_STAGE_HEADER;
		ret -> is_eof = 0;
		
		SB_FETCH(ret);

		if(SB_EOF(ret))
		{
			free(ret->input_binary_stream_buffer);
			free(ret);
			SUBREADprintf("FEOF 0!\n");
			return NULL;

		}

		int magic_4 = 0;
		memcpy(&magic_4 , SB_READ(ret), 4);
		SB_RINC(ret, 4);

		if(magic_4 != 21840194) // this number is the four bytes of "BAM1"
		{
			free(ret->input_binary_stream_buffer);
			free(ret);
			SUBREADprintf("FEOF 4 == %d!\n", magic_4);
			return NULL;
		}


		int l_text = 0;
		memcpy(&l_text, SB_READ(ret), 4);
		SB_RINC(ret, 4);

		ret -> bam_file_next_section_start = ret -> input_binary_stream_read_ptr + l_text;
		ret -> header_length =  ret -> bam_file_next_section_start;
	}
	return ret;
}

char cigar_op_char(int ch)
{
	if(ch<9)
		return "MIDNSHP=X"[ch];
	else
	{
		SUBREADprintf("Unknwon opt was found in the CIGAR string: '%c'.\n", ch);
		return 'M';
	}
}

char read_int_char(int ch)
{
	assert(ch<16);
	return "=ACMGRSVTWYHKDBN"[ch];
}

void SamBam_fclose(SamBam_FILE * fp)
{
	if(fp->file_type==SAMBAM_FILE_SAM)
	{
		fclose(fp->os_file);
		free(fp);
	}
	else
	{
		fclose(fp->os_file);
		free(fp -> input_binary_stream_buffer);
		free(fp -> bam_chro_table);
		free(fp);
	}
}

int SamBam_feof(SamBam_FILE * fp)
{
	if(fp->file_type ==SAMBAM_FILE_SAM) 
		return feof(fp->os_file);
	else return SB_EOF(fp); 
}

void SamBam_read_ref_info(SamBam_FILE * ret)
{
	unsigned int ref_info_size = 0;
	ret ->bam_chro_table_size = 0;
	//printf("CKK0\n");

	SB_FETCH(ret);
	if(SB_EOF(ret))
		return;

	//printf("CKK1\n");

	memcpy(&ref_info_size, SB_READ(ret),4);
	SB_RINC(ret, 4);

	int xk1;
	ret -> bam_chro_table = malloc(sizeof(SamBam_Reference_Info) * ref_info_size);
	for(xk1=0;xk1<ref_info_size;xk1++)
	{
		SB_FETCH(ret);
	
		if(SB_EOF(ret))
			break;

		int ref_name_len = 0;
		memcpy(&ref_name_len, SB_READ(ret),4);
		SB_RINC(ret, 4);

		int ref_readin_len = min(ref_name_len, BAM_MAX_CHROMOSOME_NAME_LEN-1);
		int ref_skip_len = ref_name_len - ref_readin_len;

		memcpy(ret -> bam_chro_table[xk1].chro_name, SB_READ(ret), ref_readin_len);
		ret -> bam_chro_table[xk1].chro_name[ref_readin_len] = 0;
		SB_RINC(ret, ref_readin_len + ref_skip_len);

		memcpy(&(ret -> bam_chro_table[xk1].chro_length), SB_READ(ret),4);
		SB_RINC(ret, 4);

		//SUBREADprintf("CHRO[%d] : %s [%d]\n", xk1+1, ret -> bam_chro_table[xk1].chro_name , ret -> bam_chro_table[xk1].chro_length);
	}
	ret ->bam_chro_table_size = ref_info_size;
}

char * SamBam_fgets(SamBam_FILE * fp, char * buff , int buff_len, int seq_needed)
{
	if(fp->file_type==SAMBAM_FILE_SAM){
		char * ret = fgets(buff, buff_len, fp->os_file);
		int strlenbuff = strlen(buff);
		if(strlenbuff < 1 || ret == NULL) return NULL;
		else{
			if(ret[strlenbuff-1]!='\n')
			{
				while(1)
				{
					int ch = getc(fp->os_file);
					if(ch == EOF || ch == '\n')break;
				}
				ret[strlenbuff-1] = '\n';
			}
			if(fp -> is_paired_end < 10){
				if(buff[0]!='@'){
					int tabs = 0,x1=0, tmpi=0;
					for(x1 = 0; x1 < strlenbuff; x1++){
						if(buff[x1] == '\t'){
							if(tabs == 1){
		//						SUBREADprintf("TMPI_SAM = %d\n", tmpi);
								fp -> is_paired_end = 10 + (tmpi & 1);
								break;
							} else tabs ++;
						}else{
							if(tabs == 1)tmpi = tmpi * 10 + buff[x1]-'0';
						}
					}
				}
			}

			if(buff[0] == '@'){
				fp -> header_length = ftello(fp->os_file) + strlenbuff+1;
			}
			return ret;
		}
	}
	else
	{
		int xk1;
		// decrypt the BAM mess.
		if(fp-> bam_file_stage == BAM_FILE_STAGE_HEADER)
		{
			char nch;
			xk1=0;
			SB_FETCH(fp);
			if(SB_EOF(fp))
				return NULL;

			while(1)
			{
				if(fp -> input_binary_stream_read_ptr >= fp -> bam_file_next_section_start)
					break;

				SB_FETCH(fp);
				nch = *(SB_READ(fp));
				SB_RINC(fp,1);

				//printf("%c", nch);
				if(nch == '\r')continue;
				if(SB_EOF(fp)||nch == '\n' || nch <0) break;
				if(xk1 < buff_len-2)
				{
					buff[xk1]=nch;
					xk1++;
				}
			}

			buff[xk1]='\n';
			buff[xk1+1]=0;

			//SUBREADprintf("BUFF=%s\n========================================================================\n\n", buff);
			//printf("RL=%d , PTR %d >? RECORD_START %d\n\n\n", xk1, fp -> input_binary_stream_read_ptr , fp -> bam_file_next_section_start);

			if(fp -> input_binary_stream_read_ptr >= fp -> bam_file_next_section_start)
			{
				SamBam_read_ref_info(fp);
				fp -> bam_file_stage = BAM_FILE_STAGE_ALIGNMENT;
				fp -> header_length = fp-> input_binary_stream_read_ptr;
			}
			return buff;
		}
		else
		{
			SamBam_Alignment *aln = &fp->aln_buff;
			int chunk_ptr = 0;
			SB_FETCH(fp);
			if(SB_EOF(fp)) return NULL;

			fp -> is_paired_end = 10 + ((*(fp -> input_binary_stream_buffer + fp -> input_binary_stream_read_ptr - fp -> input_binary_stream_buffer_start_ptr + 18)) & 1);
			//SUBREADprintf("FLAG=%d\n",  *(fp -> input_binary_stream_buffer + fp -> input_binary_stream_read_ptr - fp -> input_binary_stream_buffer_start_ptr + 18) );
			int text_len = PBam_chunk_gets(SB_READ(fp) , &chunk_ptr, fp -> input_binary_stream_write_ptr - fp -> input_binary_stream_read_ptr , fp -> bam_chro_table, buff , buff_len, aln, seq_needed);
			SB_RINC(fp, chunk_ptr);

			if(text_len>0) return buff;
			return NULL;
		}
	}
}



int PBam_get_next_zchunk(FILE * bam_fp, char * buffer, int buffer_length, unsigned int * real_len)
{
	unsigned char ID1, ID2, CM, FLG;
	unsigned short XLEN;
	int BSIZE=-1, rlen, is_file_broken = 0;

	if(feof(bam_fp)) return -1;

	fread(&ID1, 1, 1, bam_fp);
	fread(&ID2, 1, 1, bam_fp);
	fread(&CM, 1, 1, bam_fp);
	fread(&FLG, 1, 1, bam_fp);
	if(feof(bam_fp)) return -1;

	if(ID1!=31 || ID2!=139 || CM!=8 || FLG!=4)
	{
		//SUBREADprintf("4CHR = %d, %d, %d, %d\n", ID1, ID2, CM, FLG);
		return -1;
	}
	fseeko(bam_fp, 6, SEEK_CUR);
	fread(&XLEN,1, 2, bam_fp );

	int XLEN_READ = 0;
	while(1)
	{
		unsigned char SI1, SI2;
		unsigned short SLEN, BSIZE_MID;
		
		fread(&SI1, 1, 1, bam_fp);
		fread(&SI2, 1, 1, bam_fp);
		rlen = fread(&SLEN, 2, 1, bam_fp);
		if(rlen < 1) is_file_broken = 1;

		if(SI1==66 && SI2== 67 && SLEN == 2)
		{
			fread(&BSIZE_MID, 1,2 , bam_fp);
			BSIZE = BSIZE_MID;
		}
		else	fseeko(bam_fp, SLEN, SEEK_CUR);
		XLEN_READ += SLEN + 4;
		if(XLEN_READ>=XLEN) break;
	}

	if(BSIZE>19)
	{
		int CDATA_LEN = BSIZE - XLEN - 19;
		int CDATA_READING = min(CDATA_LEN, buffer_length);
		fread(buffer, 1, CDATA_READING, bam_fp);
		if(CDATA_READING<CDATA_LEN)
			fseeko(bam_fp, CDATA_LEN-CDATA_READING, SEEK_CUR);
		fseeko(bam_fp, 4, SEEK_CUR);
		rlen = fread(&real_len, 4, 1, bam_fp);
		if(rlen < 1) is_file_broken = 1;

	//	SUBREADprintf("read_data=%u\n", CDATA_LEN);
		if(is_file_broken){
			SUBREADputs("ERROR: the input BAM file is broken.");
		}
		return is_file_broken?-2:CDATA_READING;
	}
	else
		return -1;
}


// returns 0 if the header finished.
// returns 1 if the header is going on.
// returns -1 if error.
int PBam_chunk_headers(char * chunk, int *chunk_ptr, int chunk_len, SamBam_Reference_Info ** bam_chro_table, int * table_size, int * table_items, int * state, int * header_txt_remainder, int * reminder_byte_len)
{

	if((*state)  == 0)
	{
		unsigned int header_txt_len ;
		if(0!=memcmp("BAM\x1",chunk + (*chunk_ptr),4))
			return -1;
		(*chunk_ptr)+=4;	// MAGIC
		(*state) = 1;

		memcpy(&header_txt_len, chunk + (*chunk_ptr),4);
		(*chunk_ptr)+=4;	
		if(header_txt_len + 8 < chunk_len)
		{
			(* state) = 2;
			(*chunk_ptr) += header_txt_len;
		}
		else
		{
			(* state) = 1;
			(* header_txt_remainder) = header_txt_len - (chunk_len - 8); 
			return 1;
		} 
	}

	if((*state) == 1)
	{
		if((*header_txt_remainder)<chunk_len)
		{
			(*state) = 2;
			(*chunk_ptr) += (*header_txt_remainder);
		}
		else if((*header_txt_remainder)==chunk_len)
		{
			(*state) = 2;
			return 1;
		}
		else	
		{
			(* header_txt_remainder) -= (chunk_len);
			return 1;
		}
	}

	if((*state) == 2 || (*state == 3))
	{
		int chrs, remainder_chrs;
		if((*state)==2)
		{
			memcpy(&chrs, chunk + (*chunk_ptr),4); 
			(*chunk_ptr)+=4;

			remainder_chrs = chrs;
		}
		else	remainder_chrs = (*header_txt_remainder);

		while((*chunk_ptr) < chunk_len && remainder_chrs>0)
		{
			int chro_name_len;
			unsigned int chro_len;
			(*reminder_byte_len) = chunk_len - (*chunk_ptr);

			if( (*chunk_ptr) < chunk_len-4)
			{
				memcpy(&chro_name_len, chunk + (*chunk_ptr),4);
				(*chunk_ptr)+=4;
				if( (*chunk_ptr) <= chunk_len-chro_name_len-4)
				{
					char * chro_name = chunk + (*chunk_ptr);
					(*chunk_ptr)+=chro_name_len;
					memcpy(&chro_len, chunk + (*chunk_ptr),4);
					(*chunk_ptr)+=4;

					(*reminder_byte_len) =0;

					//todo: insert item
					if(0==(* table_items))
					{
						(*table_size) = 50;
						(*bam_chro_table) = malloc(sizeof(SamBam_Reference_Info)*50);
					}
					else if((*table_size) <= (* table_items))
					{
						(*table_size) *= 2;
						(*bam_chro_table) = realloc((*bam_chro_table),sizeof(SamBam_Reference_Info)*(*table_size));
					}

					SamBam_Reference_Info * new_event = (*bam_chro_table) + (* table_items);
					strncpy(new_event->chro_name, chro_name, BAM_MAX_CHROMOSOME_NAME_LEN);
					new_event -> chro_length = chro_len;

					(* table_items)++;
					//SUBREADprintf("CHRO %d/%d added\n", (* table_items),(remainder_chrs));
					remainder_chrs --;
				}
				else break;
			}
			else break;

		}

		if(remainder_chrs)
		{
			(*state) = 3;
			(*header_txt_remainder) = remainder_chrs;
			return 1;
		}
		else{
			(*state) = 4;
			return 0;
		}
	}
	return -1;
}

int convert_BAM_binary_to_SAM( SamBam_Reference_Info * chro_table, char * bam_bin, char * sam_txt ){
	int bin_len = 0;
	memcpy(&bin_len, bam_bin, 4);
	bin_len += 4;

	int sam_ptr = 0, tmpint = 0;
	sam_ptr += sprintf(sam_txt + sam_ptr, "%s\t", bam_bin+36);

	memcpy(&tmpint, bam_bin + 16 ,4);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%d\t", (tmpint >> 16) & 0xffff);
	int cigar_opts = tmpint & 0xffff;

	memcpy(&tmpint, bam_bin + 4  ,4);
	int r1chro = tmpint;
	sam_ptr += sprintf(sam_txt + sam_ptr, "%s\t", tmpint<0?"*":chro_table[tmpint].chro_name);
	memcpy(&tmpint, bam_bin + 8  ,4);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%d\t", tmpint+1);
	memcpy(&tmpint, bam_bin + 12 ,4);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%d\t", (tmpint >> 8) & 0xff);
	int name_len = tmpint & 0xff;
	int cigar_i;
	for(cigar_i = 0; cigar_i < cigar_opts; cigar_i ++){
		unsigned int cigarint=0;
		memcpy(&cigarint, bam_bin + name_len + 36 + cigar_i * 4,4);
		sam_ptr += sprintf(sam_txt + sam_ptr, "%u%c", cigarint >> 4, "MIDNSHP=X"[cigarint&0xf]);
	}
	sam_ptr += sprintf(sam_txt + sam_ptr, "%s\t", cigar_i<1?"*":"");
	
	memcpy(&tmpint, bam_bin + 24, 4);
	//SUBREADprintf("CHRO_IDX=%d\n", tmpint);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%s\t", tmpint<0?"*":((tmpint == r1chro)?"=":chro_table[tmpint].chro_name));
	
	memcpy(&tmpint, bam_bin + 28, 4);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%d\t", tmpint+1);

	memcpy(&tmpint, bam_bin + 32, 4);
	sam_ptr += sprintf(sam_txt + sam_ptr, "%d\t", tmpint);

	int seq_len;
	memcpy(&seq_len, bam_bin + 20,4);
	int seqi, flex_ptr=name_len + 36 +  cigar_opts * 4;
	for(seqi=0; seqi<seq_len; seqi++){
		sam_txt[sam_ptr++]="=ACMGRSVTWYHKDBN"[ (bam_bin[flex_ptr] >> ( 4*!(seqi %2) )) & 15 ];
		if(seqi %2) flex_ptr++;
	}
	sam_txt[sam_ptr++]='\t';
	if(seqi %2) flex_ptr++;
	for(seqi=0; seqi<seq_len; seqi++){
		unsigned char nch = (unsigned char) bam_bin[flex_ptr++];
		if(nch!=0xff||seqi == 0)
				sam_txt[sam_ptr++]=nch==0xff?'*':(nch+33);
	}
	sam_txt[sam_ptr++]='\t';

	while(flex_ptr < bin_len){
		sam_txt[sam_ptr++]=bam_bin[flex_ptr++];
		sam_txt[sam_ptr++]=bam_bin[flex_ptr++];
		sam_txt[sam_ptr++]=':';
		char tagtype = bam_bin[flex_ptr++];

		if(tagtype == 'B'){
			char elemtype = bam_bin[flex_ptr++];
			int elem_no=0, is_int_type = 0, type_bytes = 0, is_signed = 0;
			memcpy(&elem_no, bam_bin + flex_ptr, 4);
			flex_ptr += 4;
			sam_txt[sam_ptr++]='B';
			sam_txt[sam_ptr++]=':';
			sam_txt[sam_ptr++]=elemtype;
			sam_txt[sam_ptr++]=',';

			if(elemtype == 'i' || elemtype == 'I'){
				is_int_type = 1;
				type_bytes = 4;
				is_signed = elemtype == 'i' ;
			}else if(elemtype == 's' || elemtype == 'S'){
				is_int_type = 1;
				type_bytes = 2;
				is_signed = elemtype == 's' ;
			}else if(elemtype == 'c' || elemtype == 'C'){
				is_int_type = 1;
				type_bytes = 1;
				is_signed = elemtype == 's' ;
			}else if(elemtype == 'f'){
				type_bytes = 4;
			}

			int elemi;
			for(elemi =0; elemi < elem_no; elemi++){
				if(is_int_type){
					int tagval = 0;
					memcpy(&tagval,  bam_bin + flex_ptr, type_bytes);
					long long printv = is_signed?tagval:( (unsigned int) tagval );
					sam_ptr += sprintf(sam_txt + sam_ptr, "%lld,", printv);
				}else{
					float tagval = 0;
					memcpy(&tagval,  bam_bin + flex_ptr, type_bytes);
					sam_ptr += sprintf(sam_txt + sam_ptr, "%f,", tagval);
				}
				flex_ptr += type_bytes;
			}

			sam_txt[sam_ptr-1] = '\t';
			sam_txt[sam_ptr] = 0;
		}else{
				int is_int_type = 0, is_float_type = 0, type_bytes = 0, is_string_type = 0, is_char_type = 0, is_signed = 0;
				if(tagtype == 'i' || tagtype == 'I'){
					is_int_type = 1;
					type_bytes = 4;
					is_signed = tagtype == 'i' ;
				}else if(tagtype == 's' || tagtype == 'S'){
					is_int_type = 1;
					type_bytes = 2;
					is_signed = tagtype == 's' ;
				}else if(tagtype == 'c' || tagtype == 'C'){
					is_int_type = 1;
					type_bytes = 1;
					is_signed = tagtype == 's' ;
				}else if(tagtype == 'f'){
					is_float_type = 1;
					type_bytes = 4;
				}else if(tagtype == 'Z' || tagtype == 'H'){
					is_string_type = 1;
					while(bam_bin[flex_ptr+(type_bytes ++)]);
				}else if(tagtype == 'A'){
					is_char_type = 1;
					type_bytes = 1;
				}


				sam_txt[sam_ptr++]=is_int_type?'i':tagtype;
				sam_txt[sam_ptr++]=':';

				if(is_int_type){
					int tagval = 0;
					memcpy(&tagval,  bam_bin + flex_ptr, type_bytes);
					long long printv = is_signed?tagval:( (unsigned int) tagval );
					sam_ptr += sprintf(sam_txt + sam_ptr, "%lld\t", printv);
				}else if(is_string_type){
					// type_bytes includes \0
					memcpy(sam_txt + sam_ptr, bam_bin + flex_ptr, type_bytes -1);
					sam_txt[ sam_ptr + type_bytes -1 ] = '\t';

					//sam_txt[ sam_ptr + type_bytes +1]=0;
					//SUBREADprintf("STR_LEN=%d\tSTR=%s\n", type_bytes-1, sam_txt + sam_ptr);
					sam_ptr += type_bytes;
				}else if(is_float_type){
					float tagval = 0;
					memcpy(&tagval,  bam_bin + flex_ptr, type_bytes);
					sam_ptr += sprintf(sam_txt + sam_ptr, "%f\t", tagval);
				}else if(is_char_type){
					sam_txt[ sam_ptr++ ] = bam_bin[flex_ptr];
					sam_txt[ sam_ptr++ ] = '\t';
				}
				flex_ptr += type_bytes;
		}
	}

	sam_txt[sam_ptr-1]=0; //last '\t' 
	return sam_ptr-1;
}

int PBam_chunk_gets(char * chunk, int *chunk_ptr, int chunk_limit, SamBam_Reference_Info * bam_chro_table, char * buff , int buff_len, SamBam_Alignment*aln, int seq_needed)
{
	int xk1;
	// decrypt the BAM mess.
	unsigned int block_size;
	if((*chunk_ptr) +4> chunk_limit) return -1;

	memcpy(&block_size, chunk+(*chunk_ptr), 4);
	//SUBREADprintf("PBSIZE=%u\n", block_size);
	(*chunk_ptr)+=4;
	unsigned int next_start = block_size+(*chunk_ptr);

	int ref_id;
	memcpy(&ref_id, chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	if(ref_id == -1) aln -> chro_name = NULL;
	else aln -> chro_name = bam_chro_table[ref_id].chro_name; 

	memcpy(&(aln -> chro_offset), chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	unsigned int comb1;
	memcpy(&comb1, chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	aln -> mapping_quality = 0xff & (comb1 >> 8);

	unsigned int comb2;
	memcpy(&comb2, chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	aln -> flags = 0xffff&(comb2 >> 16);

	unsigned int read_len;
	memcpy(&read_len, chunk+(*chunk_ptr), 4);

	(*chunk_ptr)+=4;

	unsigned int mate_ref_id;
	memcpy(&mate_ref_id, chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	if(mate_ref_id == -1) aln -> mate_chro_name = NULL;
	else aln -> mate_chro_name = bam_chro_table[mate_ref_id].chro_name; 

	memcpy(&(aln -> mate_chro_offset), chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	memcpy(&(aln -> templete_length), chunk+(*chunk_ptr), 4);
	(*chunk_ptr)+=4;

	int read_name_len = comb1 & 0xff;
	assert(read_name_len < BAM_MAX_READ_NAME_LEN);

	memcpy(aln -> read_name, chunk+(*chunk_ptr), read_name_len);
	aln -> read_name[read_name_len] = 0;
	(*chunk_ptr)+=read_name_len;

	int cigar_ops = comb2 & 0xffff;
	aln -> cigar[0]=0; 
	for(xk1=0; xk1<cigar_ops;xk1++)
	{
		char cigar_piece_buf[BAM_MAX_CIGAR_LEN];
		unsigned int cigar_piece;

		if((*chunk_ptr) +4 > chunk_limit) return -1;
		memcpy(&cigar_piece,  chunk+(*chunk_ptr),4);
		(*chunk_ptr)+=4;

		sprintf(cigar_piece_buf, "%u%c", cigar_piece>>4, cigar_op_char(cigar_piece&0xf));
		if(strlen(cigar_piece_buf)+strlen(aln->cigar)<BAM_MAX_CIGAR_LEN-1)
			strcat(aln->cigar, cigar_piece_buf);
		else
		{
			SUBREADprintf("WARNING: cigar string is too long to the buffer.\n");
			SUBREADprintf("Please only use the compressed BAM format.\n");
			assert(0);
			return -1;
		}
	}

	char read_2_seq = 0;
	int seq_qual_bytes = read_len + (read_len /2)+(read_len%2);

	if(seq_needed)
		memcpy( aln-> buff_for_seq, chunk+(*chunk_ptr), seq_qual_bytes);
	(*chunk_ptr) += seq_qual_bytes;

	char extra_tags [CORE_ADDITIONAL_INFO_LENGTH];
	extra_tags[0]=0;
	int extra_len = 0;
	while( (*chunk_ptr) < next_start)
	{
		char extag[2];
		char extype;
		int delta, need_tag = 1;
		memcpy(extag,  chunk+(*chunk_ptr), 2);
		extype = chunk[2+(*chunk_ptr)];
		(*chunk_ptr)+=3;
		//fprintf(stderr, "COL_EXTYPE: %c\n", extype);
		if(extype == 'Z' || extype == 'H')
		{
			delta = 0;
			// 'Z' columns are NULL-terminated.
			while(chunk[delta + (*chunk_ptr)]) delta++;
			delta += 1;
		}
		else if(extype == 'A' || extype == 'c' || extype=='C') delta=1;
		else if(extype == 'i' || extype=='I' || extype == 'f') delta=4;
		else if(extype == 's' || extype=='S') delta=2;
		else if(extype == 'B') 
		{
			extype = chunk[(*chunk_ptr)];
		//	fprintf(stderr, "B_EXTYPE: %c\n", extype);

			(*chunk_ptr)++;
			if(extype == 'A' || extype=='Z') delta=1;
			else if(extype == 'c' || extype=='C') delta=1;
			else if(extype == 'i' || extype=='I' || extype == 'f') delta=4;
			else if(extype == 's' || extype=='S') delta=2;
			else break;

			int array_len;
			need_tag = 0;
			memcpy(&array_len, chunk+(*chunk_ptr), 4);
			(*chunk_ptr)+=4;
			delta *= array_len;
		}
		else{
		//	fprintf(stderr, "NO_EXTYPE: %c\n", extype);
			break;
		}

		if(need_tag){
			if(extype == 'c' || extype=='C' || extype == 'i' || extype=='I' || extype == 's' || extype=='S'){
				int tmpi = 0;
				memcpy(&tmpi, chunk+(*chunk_ptr),delta);
				if(tmpi >= 0 && extra_len < CORE_ADDITIONAL_INFO_LENGTH - 18){
					int sret = sprintf(extra_tags + strlen(extra_tags), "\t%c%c:i:%d", extag[0], extag[1], tmpi);
					extra_len += sret;
				}
			}else if(extype == 'Z'){
				if(extra_len < CORE_ADDITIONAL_INFO_LENGTH - 7 - delta){
					sprintf(extra_tags + strlen(extra_tags), "\t%c%c:Z:", extag[0], extag[1]);
					extra_len += 6;
					*(extra_tags + strlen(extra_tags)+delta-1) = 0;
					memcpy(extra_tags + strlen(extra_tags), chunk + (*chunk_ptr), delta - 1);
					extra_len += delta - 1;
				}
			}else if(extype == 'A'){
				if(extra_len < CORE_ADDITIONAL_INFO_LENGTH - 8){
					int sret = sprintf(extra_tags + strlen(extra_tags), "\t%c%c:A:%c", extag[0], extag[1], *(chunk + *chunk_ptr) );
					extra_len += sret;
				}
			}
		}

		if((*chunk_ptr) + delta > chunk_limit) return -1;
		(*chunk_ptr)+=delta;
		
	}

	if(next_start > chunk_limit) return -1;
	(*chunk_ptr) = next_start;

	if(seq_needed)
	{
		for(xk1=0;xk1<read_len;xk1++)
		{
			if(xk1 %2 == 0){
				read_2_seq = aln-> buff_for_seq[xk1/2];
			}
			if(xk1 < BAM_MAX_READ_LEN)
				aln -> sequence[xk1] = read_int_char(0xf&(read_2_seq >> (xk1%2?0:4)));
		}
		aln -> sequence[min(BAM_MAX_READ_LEN-1,read_len)] = 0;
		if(read_len >= BAM_MAX_READ_LEN-1)
			SUBREADprintf("WARNING: read is too long to the buffer\n");

		
		for(xk1=0;xk1<read_len;xk1++)
		{
			read_2_seq = aln -> buff_for_seq[(read_len /2)+(read_len%2) + xk1] ;
			if(xk1 < BAM_MAX_READ_LEN)
				aln -> seq_quality[xk1] = 33+read_2_seq;
		}
		aln -> seq_quality[min(BAM_MAX_READ_LEN-1,read_len)] = 0;
		if(aln -> seq_quality[0]==' ')
			strcpy(aln -> seq_quality, "*");
	}
	else
	{
		aln -> sequence[0]='N';
		aln -> sequence[1]=0;
		aln -> seq_quality[0]='#';
		aln -> seq_quality[1]=0;
	}

	char * chro_name = "*";
	char * cigar = "*";
	unsigned int chro_offset = 0;

	if(aln -> chro_name){
		chro_name = aln -> chro_name;
		chro_offset = aln -> chro_offset+1;
		if(aln -> cigar[0])
			cigar = aln -> cigar;
	}

	char * mate_chro_name = "*";
	unsigned int mate_chro_offset = 0;
	if(aln -> mate_chro_name)
	{
		if(aln -> mate_chro_name == chro_name) mate_chro_name = "=";
		else
			mate_chro_name = aln -> mate_chro_name;
		mate_chro_offset = aln -> mate_chro_offset+1;
	}

	long long int templete_length = aln -> templete_length;


	//fprintf(stderr, "HN_TAG=%d\n", nh_val	);

	int plen = snprintf(buff, buff_len-1, "%s\t%u\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%s\t%s%s\n%c", aln -> read_name, aln -> flags , chro_name, chro_offset, aln -> mapping_quality, cigar, mate_chro_name, mate_chro_offset, templete_length, aln -> sequence , aln -> seq_quality, extra_tags, 0);

	//fprintf(stderr,"%s", buff);

	return plen;
}


int PBum_load_header(FILE * bam_fp, SamBam_Reference_Info** chro_tab, char * remainder_reads_data , int * remainder_reads_data_len)
{
	char * CDATA = malloc(80010);
	char * PDATA = malloc(1000000);

	int chro_tab_size = 0, chro_tab_items = 0, chro_tab_state = 0, header_remainder = 0, remainder_byte_len = 0, bam_is_broken = 0; 
	z_stream strm;
	while(1)
	{
		unsigned int real_len = 0;
		int rlen = PBam_get_next_zchunk(bam_fp,CDATA,80000, & real_len);
		if(rlen<0){
			bam_is_broken = (rlen == -2);
			if(bam_is_broken){
				SUBREADprintf("BAM file format error!\n");
				free(CDATA);
				free(PDATA);
				return -1;
			}
			break;
		}

		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		strm.avail_in = 0;
		strm.next_in = Z_NULL;
		int ret = inflateInit2(&strm, SAMBAM_GZIP_WINDOW_BITS);
		if (ret != Z_OK)
		{
			free(CDATA);
			free(PDATA);
			return -1;
		}
		strm.avail_in = (unsigned int)rlen;
		strm.next_in = (unsigned char *)CDATA;


		strm.avail_out = 1000000 - remainder_byte_len;
		strm.next_out = (unsigned char *)(PDATA + remainder_byte_len);
		ret = inflate(&strm, Z_FINISH);
		int have = 1000000 - strm.avail_out;
		int PDATA_ptr=0;

		inflateEnd(&strm);

		ret = PBam_chunk_headers(PDATA, &PDATA_ptr, have, chro_tab, &chro_tab_size, &chro_tab_items, &chro_tab_state, &header_remainder,&remainder_byte_len);
		memcpy(PDATA , PDATA + have - remainder_byte_len, remainder_byte_len);
		if(ret<0)
		{
			SUBREADprintf("Header error!\n");
			free(CDATA);
			free(PDATA);
			return -1;
		}
		else if(ret == 0)
		{
			//SUBREADprintf("Header loaded = %d\n", (chro_tab_items));
			remainder_byte_len=0;
		}
		if(chro_tab_state>3){
			if(remainder_reads_data && PDATA_ptr < have)
			{
				memcpy(remainder_reads_data , PDATA + PDATA_ptr, have - PDATA_ptr);
				(*remainder_reads_data_len) =  have - PDATA_ptr ;
			}
			break;
		}
	}
	free(CDATA);
	free(PDATA);
	return 0;
}


int test_pbam(char * fname)
{
	FILE * bam_fp = f_subr_open(fname, "rb");
	char * CDATA = malloc(80010);
	char * PDATA = malloc(1000000);

	z_stream strm;
	SamBam_Reference_Info * chro_tab;

	PBum_load_header(bam_fp, & chro_tab, NULL, NULL);

	while(1)
	{
		unsigned int real_len = 0;
		int rlen = PBam_get_next_zchunk(bam_fp,CDATA,80000, & real_len);
		if(rlen<0) break;

		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		strm.avail_in = 0;
		strm.next_in = Z_NULL;
		int ret = inflateInit2(&strm, SAMBAM_GZIP_WINDOW_BITS);
		if (ret != Z_OK)SUBREADprintf("Ohh!\n");

		strm.avail_in = (unsigned int)rlen;
		strm.next_in = (unsigned char *)CDATA;


		strm.avail_out = 1000000;
		strm.next_out = (unsigned char *)PDATA;
		ret = inflate(&strm, Z_FINISH);
		int have = 1000000 - strm.avail_out;
		inflateEnd(&strm);

		int PDATA_ptr=0;

		while(PDATA_ptr < have)
		{
			char * read_line = malloc(3000);
			SamBam_Alignment  aln;
			PBam_chunk_gets(PDATA, &PDATA_ptr, have, chro_tab, read_line , 2999, &aln, 0);
			SUBREADprintf("%s", read_line);
			free(read_line);
		}
	}
	free(CDATA);
	free(PDATA);
	fclose(bam_fp);

	return 0;
}
int test_bamview(int argc, char ** argv)
{
	if(argc>1)
	{
		SamBam_FILE * fp = SamBam_fopen(argv[1], SAMBAM_FILE_BAM);
		assert(fp);
		/*
		while(1)
		{
			char buf[3000];
			char * buf2 = SamBam_fgets(fp,buf, 3000);
			//SUBREADprintf(">>%s<<\n",buf);
			//if(buf2)
			//	fwrite(buf,strlen(buf), 1, stdout);
			//else break;
		}
		*/
		SamBam_fclose(fp);
	}
	return 0;
}

int SamBam_writer_create(SamBam_Writer * writer, char * BAM_fname)
{
	memset(writer, 0, sizeof(SamBam_Writer));

	if(BAM_fname)
	{
		writer -> bam_fp = f_subr_open(BAM_fname, "wb");
		if(!writer -> bam_fp) return -1;
	}
	#ifdef MAKE_STANDALONE
	else
		writer -> bam_fp = stdout;
	#endif
	writer -> chunk_buffer = malloc(70000); 
	writer -> compressed_chunk_buffer = malloc(70000); 
	writer -> chromosome_name_table = HashTableCreate(1603);
	writer -> chromosome_id_table = HashTableCreate(1603);
	writer -> chromosome_len_table = HashTableCreate(1603);
	writer -> header_plain_text_buffer = malloc(100000000);
	writer -> header_plain_text_buffer_max = 100000000;
	writer -> header_plain_text_buffer_used = 0;

	//memset(writer -> header_plain_text_buffer , 0 , 100000000);
	HashTableSetHashFunction(writer -> chromosome_name_table , fc_chro_hash);
	HashTableSetKeyComparisonFunction(writer -> chromosome_name_table , fc_strcmp_chro);
	HashTableSetDeallocationFunctions(writer -> chromosome_name_table , free, NULL);

	return 0;
}

void SamBam_writer_chunk_header(SamBam_Writer * writer, int compressed_size)
{

	// the four magic characters
	fputc(31,  writer -> bam_fp);
	fputc(139,  writer -> bam_fp);
	fputc(8,  writer -> bam_fp);
	fputc(4,  writer -> bam_fp);

	time_t time_now = 0;
	fwrite(&time_now,4,1, writer -> bam_fp);

	int tmp_i;
	// Extra flags and OS
	fputc(0,  writer -> bam_fp);
	fputc(0xff,  writer -> bam_fp); 

	// Extra length
	tmp_i = 6;
	fwrite(&tmp_i,2,1, writer -> bam_fp);


	// SI1 and SI2 magic numbers, and SLEN
	fputc(66,  writer -> bam_fp);
	fputc(67,  writer -> bam_fp);
	tmp_i = 2;
	fwrite(&tmp_i,2,1, writer -> bam_fp);
	tmp_i = compressed_size + 19 + 6;
	fwrite(&tmp_i,2,1, writer -> bam_fp);
}

unsigned int SamBam_CRC32(char * dat, int len)
{
	unsigned int crc0 = crc32(0, NULL, 0);
	unsigned int ret = crc32(crc0, (unsigned char *)dat, len);
	return ret;
}

void SamBam_writer_add_chunk(SamBam_Writer * writer)
{
	int compressed_size ; 
	unsigned int CRC32;
	writer -> output_stream.avail_out = 70000;
	writer -> output_stream.avail_in = writer ->chunk_buffer_used;
	CRC32 = SamBam_CRC32(writer -> chunk_buffer , writer ->chunk_buffer_used);

	//FILE * dfp = f_subr_open("my.xbin","ab");
	//fwrite( writer ->chunk_buffer,  writer ->chunk_buffer_used, 1, dfp);
	//fclose(dfp);

 	int Z_DEFAULT_MEM_LEVEL = 8;
	writer -> output_stream.zalloc = Z_NULL;
	writer -> output_stream.zfree = Z_NULL;
	writer -> output_stream.opaque = Z_NULL;

	deflateInit2(&writer -> output_stream, SAMBAM_COMPRESS_LEVEL, Z_DEFLATED,
		SAMBAM_GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
	
	writer -> output_stream.next_in = (unsigned char *)writer -> chunk_buffer;
	writer -> output_stream.next_out = (unsigned char *)writer -> compressed_chunk_buffer;

	deflate(&writer -> output_stream, Z_FINISH);
	deflateEnd(&writer -> output_stream);

	compressed_size = 70000 - writer -> output_stream.avail_out;
	//printf("ADDED BLOCK=%d; LEN=%d; S=%s\n", compressed_size, writer ->chunk_buffer_used,  writer ->chunk_buffer);
	SamBam_writer_chunk_header(writer, compressed_size);
	int chunk_write_size = fwrite(writer -> compressed_chunk_buffer, 1, compressed_size, writer -> bam_fp);

	fwrite(&CRC32 , 4, 1, writer -> bam_fp);
	fwrite(&writer ->chunk_buffer_used , 4, 1, writer -> bam_fp);
	if(chunk_write_size < compressed_size){
		if(!writer -> is_internal_error)SUBREADputs("ERROR: no space left in the output directory.");
		writer -> is_internal_error = 1;
	}

	writer ->chunk_buffer_used = 0;


}

double sambam_t1 = 0;

void SamBam_writer_write_header(SamBam_Writer * writer)
{
	int header_ptr=0, header_block_start = 0;
	while(header_ptr < writer->header_plain_text_buffer_used)
	{
		if(( header_ptr - header_block_start > 55000 || header_ptr >= writer->header_plain_text_buffer_used-1) && writer -> header_plain_text_buffer[header_ptr] == '\n')
		{
			writer -> chunk_buffer_used = 0;
			if(header_block_start == 0)	// the very first block
			{
				memcpy(writer -> chunk_buffer, "BAM\1",4);
				writer -> chunk_buffer_used  = 4;
				memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used, &writer -> header_plain_text_buffer_used, 4);
				writer -> chunk_buffer_used += 4;
		
			}

			memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , writer -> header_plain_text_buffer + header_block_start, header_ptr - header_block_start+1);
			writer -> chunk_buffer_used +=  header_ptr - header_block_start + 1;
			SamBam_writer_add_chunk(writer);
			header_block_start = header_ptr + 1;
		}
		header_ptr ++;
	}

	free(writer -> header_plain_text_buffer);
	writer -> header_plain_text_buffer = NULL;

	// reference sequences
	writer -> chunk_buffer_used = 0;
	memcpy(writer -> chunk_buffer, & writer -> chromosome_name_table -> numOfElements, 4);
	writer -> chunk_buffer_used = 4;

	for( header_ptr=0 ;  header_ptr < writer -> chromosome_name_table -> numOfElements ; header_ptr ++)
	{
		//printf("D=%d\n", writer -> chromosome_id_table -> numOfElements);
		char * chro_name = HashTableGet(writer -> chromosome_id_table, NULL + 1 + header_ptr);
		unsigned int chro_len = HashTableGet(writer -> chromosome_len_table, NULL + 1 + header_ptr) - NULL - 1;
		assert(chro_name);
		int chro_name_len = strlen(chro_name)+1;

		memcpy(writer -> chunk_buffer +  writer -> chunk_buffer_used , &chro_name_len, 4);
		writer -> chunk_buffer_used += 4;

		strcpy(writer -> chunk_buffer +  writer -> chunk_buffer_used , chro_name);
		writer -> chunk_buffer_used += chro_name_len;

		memcpy(writer -> chunk_buffer +  writer -> chunk_buffer_used , &chro_len, 4);
		writer -> chunk_buffer_used += 4;

		if(header_ptr ==  writer -> chromosome_name_table -> numOfElements - 1 || writer -> chunk_buffer_used > 55000)
		{
			SamBam_writer_add_chunk(writer);
			writer -> chunk_buffer_used = 0;
		}
	}

}

int SamBam_writer_close(SamBam_Writer * writer)
{
	if(writer -> writer_state == 0)	// no reads were added
	{
		if(writer -> header_plain_text_buffer)
			SamBam_writer_write_header(writer);
	}
	else if(writer -> chunk_buffer_used)
		SamBam_writer_add_chunk(writer);
	
	writer -> chunk_buffer_used = 0;
	SamBam_writer_add_chunk(writer);
//	fputc(0, writer -> bam_fp);

	writer -> output_stream.next_in= NULL;
	writer -> output_stream.avail_in= 0;
	writer -> output_stream.next_out= NULL;
	writer -> output_stream.avail_out= 0;

	free(writer -> chunk_buffer);
	free(writer -> compressed_chunk_buffer);
	HashTableDestroy(writer -> chromosome_name_table);
	HashTableDestroy(writer -> chromosome_id_table);
	HashTableDestroy(writer -> chromosome_len_table);
	#ifdef MAKE_STANDALONE
	if(stdout != writer -> bam_fp)
	#endif
	fclose(writer -> bam_fp);

	return 0;
}

int SamBam_writer_add_header(SamBam_Writer * writer, char * header_text, int add_chro)
{
	int new_text_len = strlen(header_text);

	if(writer -> header_plain_text_buffer_max <= writer -> header_plain_text_buffer_used + new_text_len + 1)
	{
		//return 0;
		writer -> header_plain_text_buffer_max *=2;
		writer -> header_plain_text_buffer = realloc(writer -> header_plain_text_buffer ,  writer -> header_plain_text_buffer_max);
		//printf("REAL: %d : %llX\n",writer -> header_plain_text_buffer_max, (long long ) writer -> header_plain_text_buffer);
	}

	strcpy(writer -> header_plain_text_buffer + writer -> header_plain_text_buffer_used, header_text);
	writer -> header_plain_text_buffer_used += new_text_len;
	strcpy(writer -> header_plain_text_buffer + writer -> header_plain_text_buffer_used, "\n");
	writer -> header_plain_text_buffer_used ++;
	if(add_chro && memcmp(header_text, "@SQ",3)==0)
	{
		char * chro = NULL;
		int chro_len = -1;
		char * toktmp = NULL;
		char * ret_tmp = strtok_r(header_text, "\t", &toktmp);

		while(1){
			if(!ret_tmp) break;

			if(memcmp(ret_tmp,"SN:", 3)==0) chro = ret_tmp + 3;
			else if(memcmp(ret_tmp,"LN:", 3)==0) chro_len = atoi(ret_tmp + 3);

			ret_tmp = strtok_r(NULL, "\t", &toktmp);
		}

		if(chro && (chro_len>0))
			SamBam_writer_add_chromosome(writer, chro, chro_len, 0);
		
	}

	//if(writer -> header_plain_text_buffer_used %97==0) printf("MV=%d\n",writer -> header_plain_text_buffer_used);

	return 0;
}

int SamBam_writer_add_chromosome(SamBam_Writer * writer, char * chro_name, unsigned int chro_length, int add_header)
{
	unsigned int chro_id = writer -> chromosome_name_table -> numOfElements;

	//assert(strlen(chro_name) < 30);

	char * chro_name_space = malloc(strlen(chro_name)+1);
	strcpy(chro_name_space , chro_name);
	HashTablePut(writer -> chromosome_name_table, chro_name_space, NULL+1+chro_id);
	HashTablePut(writer -> chromosome_id_table, NULL+1+chro_id, chro_name_space);
	HashTablePut(writer -> chromosome_len_table, NULL+1+chro_id, NULL + 1 + chro_length);

	if(add_header)
	{
		char * line_buf = malloc(1000);
		snprintf(line_buf,999, "@SQ\tSN:%s\tLN:%u", chro_name , chro_length);
		SamBam_writer_add_header(writer, line_buf, 0);
		free(line_buf);
	}

	return 0;
}


int SamBam_compress_cigar(char * cigar, int * cigar_int, int * ret_coverage, int max_secs)
{
	int tmp_int=0;
	int cigar_cursor = 0, num_opt = 0;
	int coverage_len = 0;
	(* ret_coverage) = 0;

	if(cigar[0]=='*') return 0;
	
	while(1)
	{
		char nch = cigar[cigar_cursor++];
		if(!nch)break;
		if(isdigit(nch))
		{
			tmp_int = tmp_int*10+(nch-'0');
		}
		else
		{
			int int_opt=0;
			if(nch == 'M' || nch == 'N' || nch == 'D') coverage_len += tmp_int;
			//if(nch == 'M' ||nch == 'D' || nch == '=' || nch == 'X') coverage_len += tmp_int;
			for(; int_opt<8; int_opt++) if("MIDNSHP=X"[int_opt] == nch)break;
			cigar_int[num_opt ++] = (tmp_int << 4) | int_opt; 
			tmp_int = 0;
			//SUBREADprintf("CIGARCOM: %d-th is %c\n", num_opt, nch);
			if(num_opt>=max_secs)break;
		}
	}

	(*ret_coverage) = coverage_len;
	return num_opt;
}

void SamBam_read2bin(char * read_txt, char * read_bin)
{
	int bin_cursor = 0, txt_cursor = 0;

	while(1)
	{
		char nch = read_txt[txt_cursor++];
		if(!nch)break;
		int fourbit;
		for(fourbit=0;fourbit<15;fourbit++) if("=ACMGRSVTWYHKDBN"[fourbit] == nch)break;

		if(bin_cursor %2 == 0)  read_bin[bin_cursor/2] =  fourbit<<4;
		else read_bin[bin_cursor/2] |=  fourbit;

		bin_cursor++;
	}
}

int SamBam_compress_additional(char * additional_columns, char * bin)
{
	int col_cursor = 0 , col_len = strlen(additional_columns);
	int bin_cursor = 0;

	while(col_cursor<col_len)
	{
		if(col_cursor==0 || additional_columns[col_cursor]=='\t')
		{
			if(additional_columns[col_cursor]=='\t') col_cursor++;

			bin[bin_cursor] = additional_columns[col_cursor];
			bin[bin_cursor+1] = additional_columns[col_cursor+1];

			char datatype = additional_columns[col_cursor+3];
			if(datatype=='i' || datatype == 'f')
			{
				int dig_len =0;
				while(additional_columns[dig_len+col_cursor+5] != '\t' && additional_columns[dig_len+col_cursor+5]) dig_len++;
				int val = 0;
				float fval = 0;
				if(datatype=='i') val = atoi(additional_columns+col_cursor+5);
				else val = atof(additional_columns+col_cursor+5);

				bin[bin_cursor+2]=datatype;
				memcpy(bin+bin_cursor+3, (datatype=='i')? ((void *)&val):((void *)&fval),4);
				bin_cursor += 3 + 4;
				col_cursor += 5 + dig_len;
			}
			else if(datatype=='Z' || datatype == 'H')
			{
				bin[bin_cursor+2]=datatype;
				bin_cursor +=3;
				int str_len = 0;
				col_cursor +=5;
				while(additional_columns[str_len+col_cursor] != '\t' && additional_columns[str_len+col_cursor])
				{
					bin[bin_cursor + str_len] = additional_columns[str_len+col_cursor];
					str_len++;
					if(bin_cursor + str_len > 280) break;
				}

				bin[bin_cursor + str_len] =0;

				bin_cursor += str_len + 1;
				col_cursor += str_len;
			}
			else if(datatype=='A')
			{
				bin[bin_cursor+2]='A';
				bin[bin_cursor+3]=additional_columns[col_cursor+5];
				col_cursor += 6;
				bin_cursor += 4;
			}
			else if(datatype=='B')
				//array
			{
				char celltype = additional_columns[col_cursor+5];
				int * items = (int *)(&bin[bin_cursor+4]);

				bin[bin_cursor+2]='B';
				bin[bin_cursor+3]=celltype;
				bin_cursor += 4 + 4;
				col_cursor += 7;

				(*items) = 0;

				int last_cursor = col_cursor;
				while(1){
					if(additional_columns[col_cursor] == ',' || additional_columns[col_cursor] == '\t' || additional_columns[col_cursor] == 0)
					{ // add new item 

						char cell_buff [30];
						if((col_cursor - last_cursor) < 29)
						{
							memcpy(cell_buff, additional_columns + last_cursor, (col_cursor - last_cursor));
							cell_buff[(col_cursor - last_cursor)] = 0;
							int intv = 0; float fltv = 0;
							if(celltype == 'i')intv = atoi(cell_buff);							
							else fltv = atof(cell_buff);
							if(bin_cursor < 280){
								memcpy(bin + bin_cursor, (celltype == 'i')?(void *)&intv:(void *)&fltv, 4);
								bin_cursor += 4;
								(*items) ++;
							}
						}
						last_cursor = col_cursor+1;
					}
					if(additional_columns[col_cursor] == '\t' || additional_columns[col_cursor] == 0)
						break;

					col_cursor++;
					
				}
				
			}

			if(bin_cursor>250) break;
			continue;
		}
		
		col_cursor++;
	}
	return bin_cursor;
}

int SamBam_reg2bin(int beg, int end)
{
	--end;
	if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
	if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
	if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
	if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
	if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
	return 0;
}

#define FC_MAX_CIGAR_SECTIONS 96

int SamBam_writer_add_read(SamBam_Writer * writer, char * read_name, unsigned int flags, char * chro_name, unsigned int chro_position, int mapping_quality, char * cigar, char * next_chro_name, unsigned int next_chro_position, int temp_len, int read_len, char * read_text, char * qual_text, char * additional_columns)
{
	if(writer -> writer_state == 0)	// no reads were added
	{
		if(writer -> header_plain_text_buffer)
			SamBam_writer_write_header(writer);
	}

	if(!qual_text || !read_text)	
	{
		SUBREADprintf("ERROR: sam file is incomplete.\n");
		return 1;
	}

	writer -> writer_state = 10;
	char additional_bin[300];
	int cigar_opts[FC_MAX_CIGAR_SECTIONS], xk1, cover_length = 0;
	int cigar_opt_len = SamBam_compress_cigar(cigar, cigar_opts, & cover_length, FC_MAX_CIGAR_SECTIONS);
	int read_name_len = 1+strlen(read_name) ;
	int additional_bin_len = SamBam_compress_additional(additional_columns, additional_bin);
	int record_length = 4 + 4 + 4 + 4 +  /* l_seq: */ 4 + 4 + 4 + 4 + /* read_name:*/ read_name_len + cigar_opt_len * 4 + (read_len + 1) /2 + read_len + additional_bin_len;

	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , & record_length , 4);
	writer -> chunk_buffer_used += 4;

	int bin = SamBam_reg2bin(chro_position -1, chro_position-1+cover_length);

	int refID = HashTableGet(writer -> chromosome_name_table, chro_name) - NULL - 1; 
	int bin_mq_nl = (bin<<16) | (mapping_quality << 8) | read_name_len ;
	int fag_nc = (flags<<16) | cigar_opt_len;
	int nextRefID = -1;

	if(next_chro_name[0] != '*' && next_chro_name[0]!='=')
		nextRefID = HashTableGet(writer -> chromosome_name_table, next_chro_name) - NULL - 1;
	else if(next_chro_name[0] == '=')
		nextRefID = refID;

	
	chro_position--;
	next_chro_position--;

	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , & refID , 4);
	writer -> chunk_buffer_used += 4;
	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , & chro_position , 4);
	writer -> chunk_buffer_used += 4;
	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , & bin_mq_nl , 4);
	writer -> chunk_buffer_used += 4;
	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , & fag_nc , 4);
	writer -> chunk_buffer_used += 4;
	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , & read_len , 4);
	writer -> chunk_buffer_used += 4;
	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , & nextRefID , 4);
	writer -> chunk_buffer_used += 4;
	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , & next_chro_position , 4);
	writer -> chunk_buffer_used += 4;
	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , & temp_len , 4);
	writer -> chunk_buffer_used += 4;
	strcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , read_name);
	writer -> chunk_buffer_used += read_name_len;
	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used , cigar_opts, 4*cigar_opt_len);
	writer -> chunk_buffer_used += 4*cigar_opt_len;
	SamBam_read2bin(read_text  , writer -> chunk_buffer + writer -> chunk_buffer_used);
	writer -> chunk_buffer_used += (read_len + 1) /2; 
	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used, qual_text, read_len);
	for(xk1=0; xk1<read_len; xk1++)
		writer -> chunk_buffer[writer -> chunk_buffer_used+xk1] -= 33;
	
	writer -> chunk_buffer_used += read_len; 
	memcpy(writer -> chunk_buffer + writer -> chunk_buffer_used, additional_bin, additional_bin_len);
	writer -> chunk_buffer_used += additional_bin_len;


	if(writer -> chunk_buffer_used>55000)
	{
		//	double t0 = miltime();
		SamBam_writer_add_chunk(writer);
		//	double t1 = miltime();
		//	if(sambam_t1 > 100)
		//		SUBREADprintf("Running = %.6f , Compress Time = %.6f\n", t0 - sambam_t1, t1 - t0);
		//	sambam_t1 = t1;
	



		writer -> chunk_buffer_used = 0;
	}	
	return 0;
}

int SamBam_unzip(char * out , char * in , int inlen)
{
	#define unzip_out_max_len 65537
	z_stream strm;
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	int ret = inflateInit2(&strm, SAMBAM_GZIP_WINDOW_BITS);
	if (ret != Z_OK)
		return -1;

	strm.avail_in = (unsigned int)inlen;
	strm.next_in = (unsigned char *)in;

	strm.avail_out = unzip_out_max_len;
	strm.next_out = (unsigned char *)out;
	ret = inflate(&strm, Z_FINISH);
	if(ret != Z_STREAM_END)
	{
		inflateEnd(&strm);
		SUBREADprintf("DATA ERROR! code=%d\n", ret);
		return -1;
	}
	int have = unzip_out_max_len - strm.avail_out;

	inflateEnd(&strm);
	//SUBREADprintf("DECOMPRESS GENERATED=%d\n", have);

	return have;
}




#ifdef MAKE_TEST_SAMBAM

int is_badBAM(char * fn)
{
	FILE * fp = f_subr_open(fn , "r");
	if(!fp) return -1;

	char * in_buff = malloc(70000);
	char * out_buff = malloc(170000);
	int blks=0;

	int state = 0;
	int chros = 0, all_chros;

	unsigned int data_ptr = 0;
	unsigned int chunk_start_ptr = 0;
	unsigned int head_text_len = 0;
	unsigned int tested_chunks = 0;
	unsigned int tested_reads = 0;

	int fret = 0;
	int last_len = 0, last_val = 0;

	while (!feof(fp))
	{

		int real_len = 0, BSIZE=0;
		int ID1=0, ID2=0, CM=0, FLG=0, XLEN=0;

		fread(&ID1, 1, 1, fp);
		if(feof(fp))
		{
			if(ID1!=0 || blks==0)fret = 2;
			break;
		}

		fread(&ID2, 1, 1, fp);
		fread(&CM, 1, 1, fp);
		fread(&FLG, 1, 1, fp);


		if(ID1!=31 || ID2!=139 || CM!=8 || FLG!=4)
		{
			fret = 2;
			break;
		}
		fseeko(fp, 6, SEEK_CUR);
		fread(&XLEN,1, 2, fp );

		int XLEN_READ = 0;
		while(1)
		{
			unsigned char SI1, SI2;
			unsigned short SLEN, BSIZE_MID;
			
			fread(&SI1, 1, 1, fp);
			fread(&SI2, 1, 1, fp);
			fread(&SLEN, 1, 2, fp);

			if(SI1==66 && SI2== 67 && SLEN == 2)
			{
				fread(&BSIZE_MID, 1,2 , fp);
				BSIZE = BSIZE_MID;
			}
			else	fseeko(fp, SLEN, SEEK_CUR);
			XLEN_READ += SLEN + 4;
			if(XLEN_READ>=XLEN) break;
		}

		if(BSIZE>19)
		{
			int CDATA_LEN = BSIZE - XLEN - 19;
			int CDATA_READING = CDATA_LEN;
			fread(in_buff, 1, CDATA_READING, fp);
			if(CDATA_READING<CDATA_LEN)
				fseeko(fp, CDATA_LEN-CDATA_READING, SEEK_CUR);
			fseeko(fp, 4, SEEK_CUR);
			fread(&real_len, 4, 1, fp);





			z_stream strm;


			strm.zalloc = Z_NULL;
			strm.zfree = Z_NULL;
			strm.opaque = Z_NULL;
			strm.avail_in = 0;
			strm.next_in = Z_NULL;
			int ret = inflateInit2(&strm, SAMBAM_GZIP_WINDOW_BITS);
			if (ret != Z_OK)
			{
				fret = 2;
				break;
			}
			strm.avail_in = (unsigned int)CDATA_READING;
			strm.next_in = (unsigned char *)in_buff;


			strm.avail_out = 70000;
			strm.next_out = (unsigned char *)out_buff;
			ret = inflate(&strm, Z_FINISH);

			if (ret != Z_STREAM_END)
			{
				fret = 2;
				break;
			}
		

			int have = 70000 - strm.avail_out;

			inflateEnd(&strm);

			if(state == 0)
			{
				data_ptr = 4;
				if(memcmp(out_buff, "BAM\1", 4)!=0)
				{
					fret=2;
					break;
				}
				memcpy(&head_text_len , out_buff+4 , 4);
				state = 1;

				//printf("header=%d\n", head_text_len);
			}

			if(state == 1)
			{
				//printf("chunk_end=%d\n", chunk_start_ptr + have );
				if(chunk_start_ptr + have >= head_text_len + 8)
				{
					data_ptr =  head_text_len + 8;
					state = 2;
				}
			}

			if(state == 2 && data_ptr <  chunk_start_ptr+have)
			{
				memcpy(& chros, out_buff + (data_ptr - chunk_start_ptr), 4);

				all_chros = chros;
				printf("chros=%d\n", chros);
				data_ptr +=4;
				state = 3;
			}

			if(state == 3 && data_ptr <  chunk_start_ptr+have)
			{

				while(data_ptr <= chunk_start_ptr + have - 4)
				{
					int ref_name_len ;
					memcpy(& ref_name_len ,  out_buff + (data_ptr - chunk_start_ptr), 4 - last_len);

					if(last_len)
					{
						ref_name_len = (ref_name_len<< (8*last_len)) + ref_name_len;
						last_len = 0;
					}

					{
					//	char chn[300];
					//	memcpy(chn, out_buff + (data_ptr - chunk_start_ptr) + 4, ref_name_len);
		//				printf("skipped=%d (%s)\n", ref_name_len+8, chn);
					}

					data_ptr += ref_name_len + 8 - last_len;
					if(chros ==1){
						state = 4;
						printf("header len = %d\n", data_ptr);
						break;
					}
					else
						chros --;
		//			printf("chros-=%d\n", chros);
				}

				if( data_ptr > chunk_start_ptr + have - 4 && data_ptr <  chunk_start_ptr + have)
				{
					last_len = chunk_start_ptr + have - data_ptr;
					last_val = 0;
					memcpy(&last_val, out_buff + (data_ptr - chunk_start_ptr) , last_len);
					data_ptr = chunk_start_ptr + have ;
				}
				else last_len = 0;
	
			}

			if(state == 4 && data_ptr <  chunk_start_ptr+have)
			{
				tested_chunks ++;
				if(tested_chunks > TEST_BAD_BAM_CHUNKS)
					break;
				printf("tested_chunks=%d; reads=%d; data_ptr=%d; start_ptr=%d; have=%d\n", tested_chunks, tested_reads, data_ptr, chunk_start_ptr,  have);
				if(data_ptr!= chunk_start_ptr && have > 9993000)
				{
					//printf("PTRS: %d != %d - %d\n", data_ptr, chunk_start_ptr, have);
					fret=1;
					break;
				}

				int my_r = 0;
				while(data_ptr <= chunk_start_ptr + have - 4)
				{
					int read_len =0, my_chro;

					memcpy(&read_len,  out_buff + (data_ptr - chunk_start_ptr) , (4 - last_len));
					if(last_len)
						read_len = (read_len << (8*last_len)) + last_val;
					

					data_ptr += (4-last_len);
					int bin_mq_nl, pos, flag_nc;
					memcpy(&pos, out_buff + data_ptr - chunk_start_ptr+ 4,4);
					memcpy(&bin_mq_nl, out_buff + data_ptr - chunk_start_ptr+ 8,4);
					memcpy(&flag_nc, out_buff + data_ptr - chunk_start_ptr+12,4);
					int read_name_ptr = data_ptr+32;
					if(memcmp(out_buff+read_name_ptr - chunk_start_ptr, "V0112_0155:7:1101:1818:190479", strlen("V0112_0155:7:1101:1818:190479")) == 0)
						SUBREADprintf("BIN_MG_NL = %08X ; POS=%d; FLAG=%04X\n", bin_mq_nl, pos, flag_nc);

					last_len = 0;
					memcpy(&my_chro,  out_buff + (data_ptr - chunk_start_ptr) , 4);
					data_ptr += read_len;

					printf("the %d-th read_len=%d\n", tested_reads+1, read_len);

					if(read_len>10000)
					{
						ret = 2;
						break;
					}
					my_r++;
					tested_reads ++;
				}

				if( data_ptr > chunk_start_ptr + have - 4 && data_ptr <  chunk_start_ptr + have)
				{
					last_len = chunk_start_ptr + have - data_ptr;
					last_val = 0;
					memcpy(&last_val, out_buff + (data_ptr - chunk_start_ptr) , last_len);
					data_ptr = chunk_start_ptr + have ;
				}
				else last_len = 0;
			}

			chunk_start_ptr += have;


			blks++;
			
		}
		else if(blks < 1) fret=2;

		if(fret) break;

	}

	fclose(fp);

	free(in_buff);
	free(out_buff);


	return fret;
}

int main(int argc , char ** argv)
{
	int bad = is_badBAM(argv[1]);
	printf("BAD BAM=%d\t\t%s\n", bad, argv[1]);

	SamBam_FILE * r2fp = SamBam_fopen(argv[1], SAMBAM_FILE_BAM);

	while(1)
	{
		char read_buff[3000];
		char * ret = SamBam_fgets(r2fp, read_buff, 2999, 0);
		if(!ret)break;
		SUBREADprintf("%s", ret);
		
	}

	return 0;
}

void test_bam_compress()
{
	SamBam_Writer writer;
	if(SamBam_writer_create(&writer , "my.bam")) printf("INIT ERROR\n");

	SamBam_writer_add_header(&writer, "@RG	ID:xxhxh", 0);
	SamBam_writer_add_chromosome(&writer, "chr1", 123123, 1);
	SamBam_writer_add_chromosome(&writer, "chr2", 223123, 1);
	SamBam_writer_add_header(&writer, "@PG	ID:subread	VN:1.4.0b2", 0);
	SamBam_writer_add_read(& writer, "Read1", 0, "chr1", 100000, 200, "50M", "*", 0, 0, 50, "ATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGA", "AAAAABBBBBAAAAABBBBBAAAAABBBBBAAAAABBBBBAAAAABBBBB", "XG:Z:OX	NM:i:2	RG:Z:MyGroup1");
	SamBam_writer_add_read(& writer, "Read2", 16, "chr2", 200000, 200, "50M", "*", 0, 0, 50, "ATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGA", "AAAAABBBBBAAAAABBBBBAAAAABBBBBAAAAABBBBBAAAAABBBBB", "NM:i:1	XX:i:8172736	RG:Z:nxnmn	XY:i:33999	XZ:Z:Zuzuzu");
	SamBam_writer_close(&writer);
}


#endif
