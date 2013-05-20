#ifndef _SAMBAM_FILE_H_
#define _SAMBAM_FILE_H_

typedef unsigned char BS_uint_8;
typedef unsigned short BS_uint_16;
typedef unsigned int BS_uint_32;

#define BAM_MAX_CHROMOSOME_NAME_LEN 256 
#define BAM_MAX_CIGAR_LEN 64
#define BAM_MAX_READ_NAME_LEN 256
#define BAM_MAX_READ_LEN 3000

#define SAMBAM_FILE_SAM	10
#define SAMBAM_FILE_BAM 20

#define BAM_FILE_STAGE_HEADER 10
#define BAM_FILE_STAGE_ALIGNMENT 20


typedef struct
{
	char chro_name[BAM_MAX_CHROMOSOME_NAME_LEN];
	unsigned int chro_length;
} SamBam_Reference_Info;


typedef struct
{
	char read_name[BAM_MAX_READ_NAME_LEN];
	char * chro_name;
	unsigned int chro_offset;
	unsigned short flags;
	char * mate_chro_name;
	unsigned int mate_chro_offset;
	int templete_length;
	unsigned char mapping_quality;

	char cigar[BAM_MAX_CIGAR_LEN];
	char sequence[BAM_MAX_READ_LEN];
	char seq_quality[BAM_MAX_READ_LEN];

	char buff_for_seq[BAM_MAX_READ_LEN*2];

} SamBam_Alignment;


typedef struct
{
	union{
		FILE * os_file;
		gzFile gz_file;
	};
	int file_type;
	int bam_file_stage;
	unsigned long long bam_file_next_section_start;
	SamBam_Reference_Info * bam_chro_table;
	int bam_chro_table_size;
	SamBam_Alignment aln_buff;
} SamBam_FILE;


// This function opens a file, either SAM or BAM, in read-only mode.
// The "file_type" parameter specifies which type of file it is: SAMBAM_FILE_BAM or SAMBAM_FILE_SAM.
SamBam_FILE * SamBam_fopen(const char * fname , int file_type);

// This function closes any opened file and releases memory footprint. It works just like "fclose()".
void SamBam_fclose(SamBam_FILE * fp);

// This function tells if a file is exhausted.
// Note that a non-exhausted file can still contain no more alignment results.
// Hence, it is recommended to check the return value of SamBam_fgets() to tell if the file has reached its end.
int SamBam_feof(SamBam_FILE * fp);

// This function works like fgets except it decode the BAM file.
// If the buffer is not long enough to store the line, the remainder of this line is omitted and the next call will read the next alignment.
// A very important difference between fgets and SamBam_fgets is that this function returns NULL when there are no more lines.
// It is recommended to use the return value as the indicator of EOF like:
/**
 * SamBam_FILE * fp = SamBam_fopen("my.bam", SAMBAM_FILE_BAM);
 * while(1)
 * {
 *   char buf[3000];
 *   char * ret = SamBam_fgets(fp, buf, 3000);
 *   if(ret) puts(buf);
 *   else break;
 * }
 * SamBam_fclose(fp);
 */
char * SamBam_fgets(SamBam_FILE * fp , char * buff , int buff_len );

#endif
