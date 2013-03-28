#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "subread.h"
#include "sambam-file.h"

BS_uint_16 gzread_B16(gzFile fp)
{
	BS_uint_16 ret;
	gzread(fp, &ret, 2);
	return ret;
}

BS_uint_32 gzread_B32(gzFile fp)
{
	BS_uint_32 ret;
	gzread(fp, &ret, 4);
	return ret;
}

BS_uint_8 gzread_B8(gzFile fp)
{
	BS_uint_8 ret;
	gzread(fp, &ret, 1);
	return ret;
}


SamBam_FILE * SamBam_fopen(const char * fname , int file_type)
{
	SamBam_FILE * ret = (SamBam_FILE *)malloc(sizeof(SamBam_FILE));
	ret -> file_type = file_type;

	if(file_type ==SAMBAM_FILE_SAM) 
	{
		ret -> os_file = fopen(fname, "rb");
		if(!ret -> os_file)
		{
			free(ret);
			return NULL;
		}
		fseek(ret -> os_file,0,SEEK_SET);
	}
	else
	{
		FILE * os_file = fopen(fname, "rb");
		if(os_file == NULL)
		{
			free(ret);
			return NULL;
		}
		unsigned char first_ch = fgetc(os_file);
		unsigned char second_ch = fgetc(os_file);

		fclose(os_file);
		if(first_ch!=31 || second_ch!=139)
		{
			free(ret);
			return NULL;
		}


		gzFile nf = gzopen(fname, "rb");
		if(!nf)
		{
			free(ret);
			return NULL;
		}

		ret -> gz_file = nf;
		ret -> bam_file_stage = BAM_FILE_STAGE_HEADER;
		
		BS_uint_32 magic_4 = gzread_B32(ret -> gz_file);
		if(magic_4 != 21840194) // this number is the four bytes of "BAM1"
		{
			free(ret);
			return NULL;
		}
		BS_uint_32 l_text = gzread_B32(ret -> gz_file);
		ret -> bam_file_next_section_start = gztell(ret -> gz_file) + l_text;
	}
	return ret;
}

char cigar_op_char(int ch)
{
	assert(ch<9);
	return "MIDNSHP=X"[ch];
}

char read_int_char(int ch)
{
	assert(ch<16);
	return "=ACMGRSVTWYHKDBN"[ch];
}

int SamBam_get_alignment(SamBam_FILE * fp, SamBam_Alignment * aln)
{
		if(gzeof(fp->gz_file)) return -1;
		unsigned long long head_pos = gztell(fp->gz_file);
		unsigned int block_size = gzread_B32(fp->gz_file);
		if(gzeof(fp->gz_file)) return -1;

		unsigned int ref_id = gzread_B32(fp->gz_file);
		if(gzeof(fp->gz_file)) return -1;

		assert(ref_id < fp->bam_chro_table_size|| ref_id == -1);


		if(ref_id == -1) aln -> chro_name = NULL;
		else aln -> chro_name = fp -> bam_chro_table[ref_id].chro_name; 
		aln -> chro_offset = gzread_B32(fp->gz_file);

		unsigned int comb1 = gzread_B32(fp->gz_file);
		aln -> mapping_quality = 0xff & (comb1 >> 8);

		unsigned int comb2 = gzread_B32(fp->gz_file);
		aln -> flags = 0xffff&(comb2 >> 16);

		unsigned int read_len = gzread_B32(fp->gz_file);

		unsigned int mate_ref_id = gzread_B32(fp->gz_file);
		if(gzeof(fp->gz_file)) return -1;

		assert(mate_ref_id < fp->bam_chro_table_size || mate_ref_id == -1);
		if(mate_ref_id == -1) aln -> mate_chro_name = NULL;
		else aln -> mate_chro_name = fp -> bam_chro_table[mate_ref_id].chro_name; 

		aln -> mate_chro_offset = gzread_B32(fp->gz_file);

		aln -> templete_length = (int)gzread_B32(fp->gz_file);

		int read_name_len = comb1 & 0xff;
		assert(read_name_len < BAM_MAX_READ_NAME_LEN);

		gzread(fp->gz_file, aln -> read_name, read_name_len);
		aln -> read_name[read_name_len] = 0;

		int cigar_ops = comb2 & 0xffff;
		int xk1;
		aln -> cigar[0]=0; 
		for(xk1=0; xk1<cigar_ops;xk1++)
		{
			char cigar_piece_buf[BAM_MAX_CIGAR_LEN];
			unsigned int cigar_piece = gzread_B32(fp->gz_file);

			sprintf(cigar_piece_buf, "%u%c", cigar_piece>>4, cigar_op_char(cigar_piece&0xf));
			if(strlen(cigar_piece_buf)+strlen(aln->cigar)<BAM_MAX_CIGAR_LEN-1)
				strcat(aln->cigar, cigar_piece_buf);
			else
				SUBREADprintf("WARNING: cigar string is too long to the buffer\n");
		}

		char read_2_seq;
		int seq_qual_bytes = read_len + (read_len /2)+(read_len%2);
		int gzread_len = gzread(fp->gz_file, aln-> buff_for_seq, seq_qual_bytes);
		if(gzread_len < seq_qual_bytes)
			return -1;

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

		unsigned long long tail_pos = gztell(fp->gz_file);

		unsigned long long skip_len = block_size - (tail_pos - head_pos - 4);

		long long int seek_ret= gzseek(fp->gz_file, skip_len, SEEK_CUR);
		if(seek_ret < 0) return -1;
		return 0;
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
		gzclose(fp->gz_file);
		free(fp);
	}
}

int SamBam_feof(SamBam_FILE * fp)
{
	if(fp->file_type==SAMBAM_FILE_SAM) return feof(fp->os_file);
	else{
			if(gzeof(fp->gz_file)) return 1;
			return 0;
	}
}

void SamBam_read_ref_info(SamBam_FILE * ret)
{
	unsigned int ref_info_size = gzread_B32(ret -> gz_file);

	int xk1;
	for(xk1=0;xk1<ref_info_size;xk1++)
	{
		int ref_name_len = gzread_B32(ret -> gz_file);
		int ref_readin_len = min(ref_name_len, BAM_MAX_CHROMOSOME_NAME_LEN-1);
		int ref_skip_len = ref_name_len - ref_readin_len;

		gzread(ret -> gz_file, ret -> bam_chro_table[xk1].chro_name , ref_readin_len);
		ret -> bam_chro_table[xk1].chro_name[ref_readin_len] = 0;
		if(ref_skip_len)gzseek(ret -> gz_file , ref_skip_len , SEEK_CUR);


		ret -> bam_chro_table[xk1].chro_length = gzread_B32(ret -> gz_file);


		//printf("CHRO[%d] : %s [%d]\n", xk1+1, ret -> bam_chro_table[xk1].chro_name , ret -> bam_chro_table[xk1].chro_length);
		if(xk1 >= BAM_MAX_CHROMOSOME_NUMBER)
			SUBREADprintf("WARNING: There are too many reference sequences in the BAM file!\n"); 
	}
	ret ->bam_chro_table_size = ref_info_size;
}

char * SamBam_fgets(SamBam_FILE * fp, char * buff , int buff_len)
{
	if(fp->file_type==SAMBAM_FILE_SAM){
		char * ret = fgets(buff, buff_len, fp->os_file);
		if(strlen(buff)<1) return NULL;
		else return ret;
	}
	else
	{
		int xk1;
		// decrypt the BAM mess.
		if(fp-> bam_file_stage == BAM_FILE_STAGE_HEADER)
		{
			char nch;
			xk1=0;
			while(1)
			{
				if(xk1 >= buff_len-2 || gztell(fp->gz_file) >= fp -> bam_file_next_section_start)
					break;

				nch = gzgetc(fp->gz_file);
				if(nch == '\r'||nch=='\n' || nch <0) break;
				buff[xk1]=nch;
				xk1++;
			}

			if(xk1<buff_len-1){
				buff[xk1]='\n';
				buff[xk1+1]=0;
			}
			if(gztell(fp->gz_file) >= fp -> bam_file_next_section_start)
			{
				SamBam_read_ref_info(fp);
				fp -> bam_file_stage = BAM_FILE_STAGE_ALIGNMENT;
			}
			return buff;
		}
		else
		{
			SamBam_Alignment *aln = &fp->aln_buff;
			int is_align_error =SamBam_get_alignment(fp, aln);
			if(is_align_error)return NULL;
			else
			{
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

					snprintf(buff, buff_len-1, "%s\t%u\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%s\t%s\n", aln -> read_name, aln -> flags , chro_name, chro_offset, aln -> mapping_quality, cigar, mate_chro_name, mate_chro_offset, templete_length, aln -> sequence , aln -> seq_quality);
			}
		
			return buff;
		}
	}
}






// test function
#ifdef MAKE_STANDALONE
int main(int argc , char ** argv)
#else
int test_bamview(int argc, char ** argv)
#endif
{
	if(argc>1)
	{
		SamBam_FILE * fp = SamBam_fopen(argv[1], SAMBAM_FILE_BAM);
		assert(fp);

		while(1)
		{
			char buf[3000];
			char * buf2 = SamBam_fgets(fp,buf, 3000);
			//printf(">>%s<<\n",buf);
			if(buf2)
				fwrite(buf,strlen(buf), 1, stdout);
			else break;
		}

		SamBam_fclose(fp);
	}
	return 0;
}

