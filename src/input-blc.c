#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>

#include "subread.h"
#include "HelperFunctions.h"
#include "seek-zlib.h"
#include "gene-algorithms.h"
#include "input-blc.h"


struct iBLC_scan_t{
	char out_format_string[MAX_FILE_NAME_LENGTH];
	char filter_format_string[MAX_FILE_NAME_LENGTH];
	int found_answer;
	int filter_is_gzipped;
	int bcl_is_gzipped;
	int reads_per_cluster;
	int read_lengths[INPUT_BLC_MAX_READS];
	int read_is_index[INPUT_BLC_MAX_READS];
};

int iBLC_guess_scan(struct iBLC_scan_t * scancon, char * data_dir ){
	DIR * this_level = opendir(data_dir);
	struct dirent *dp;
	int filter_found = 0, bcl_found = 0;
	char testfile_name[MAX_FILE_NAME_LENGTH];
	while ((dp = readdir (this_level)) != NULL) {
		if(dp -> d_type == DT_DIR && dp->d_name[0]!='.'){
			strcpy(testfile_name,data_dir);
			strcat(testfile_name, "/");
			strcat(testfile_name, dp->d_name);
			//SUBREADprintf("DIG: %s\n", testfile_name);
			if(iBLC_guess_scan( scancon, testfile_name))return -1;
		}else if(dp -> d_type == DT_REG){
			//SUBREADprintf( "%s  %s  %p  %p\n" , data_dir, dp->d_name , strstr( dp->d_name , "0001.bcl." ) , strstr( dp->d_name , ".bci") );
			if(0==strcmp(dp->d_name, "RunInfo.xml")){
				if(scancon->reads_per_cluster > 0){
					SUBREADprintf("ERROR: the root directory contains multiple scRNA data sets.\n");
					return -1;
				}

				strcpy(testfile_name, data_dir);    
				strcat(testfile_name, "/");
				strcat(testfile_name, dp->d_name);
				FILE *fp = fopen(testfile_name,"r");
				if(NULL == fp){
					SUBREADprintf("ERROR: cannot open the run info file: %s\n", testfile_name);
				}
				while(1){
					char inbuf[MAX_READ_LENGTH];
					if(!fgets( inbuf, MAX_READ_LENGTH-1, fp))break;
					if(strstr( inbuf, "<Read Number=\"" )){
						char * rbuf=NULL;
						int my_index = -1, is_idx = -1, rlen = -1, ii=0;
						strtok_r(inbuf, "\"", &rbuf);

						while(rbuf){
							char * sec = strtok_r(NULL, "\"", &rbuf);
							if(!sec) break;
							//printf("SEC %d : %s\n", ii, sec);
							if(ii == 0) my_index = atoi(sec);
							if(ii == 2) rlen = atoi(sec);
							if(ii == 4) is_idx = sec[0]=='Y';
							ii++;
						}
						assert(INPUT_BLC_MAX_READS>my_index);
						if(my_index >0 && is_idx >=0 && rlen>0){
							scancon -> read_lengths[my_index-1]=rlen;
							scancon -> read_is_index[my_index-1]=is_idx;
							scancon -> reads_per_cluster = max(scancon -> reads_per_cluster, my_index);
						}else assert( my_index >0 && is_idx >=0 && rlen>0 );
					}
				}
				fclose(fp);
				if(scancon -> reads_per_cluster <1){
					SUBREADprintf("ERROR: the format of RunInfo.xml is unknown\n");
					return -1;
				}
			}
			if(0==memcmp(data_dir+ strlen(data_dir)-5, "/L001",5 ) && strstr( dp->d_name , "s_1.filter")){
				autozip_fp tfp;
				strcpy(testfile_name, data_dir);    
				strcat(testfile_name, "/");
				strcat(testfile_name, dp->d_name);
				int resop = autozip_open(testfile_name, &tfp);
				if(0 <= resop){
					autozip_close(&tfp);
					char * gen_fmt = str_replace(dp->d_name , "s_1.filter", "s_%d.filter");
					char * gen_fmt2 = str_replace(data_dir , "/L001", "/L%03d");
					strcpy(scancon -> filter_format_string, gen_fmt2);
					strcat(scancon -> filter_format_string, "/");
					strcat(scancon -> filter_format_string, gen_fmt);
					free(gen_fmt2);
					free(gen_fmt);
					filter_found = resop + 1;
				}
			}
			if(0==memcmp(data_dir+ strlen(data_dir)-5, "/L001",5 ) && strstr( dp->d_name , "0001.bcl." ) && !strstr( dp->d_name , ".bci") ){
				int tti;
				bcl_found = 1;
				char * gen_fmt = str_replace(dp->d_name , "0001.bcl.", "%04d.bcl.");
				
				for(tti = 0; tti<22; tti++){
					strcpy(testfile_name, data_dir);	
					strcat(testfile_name, "/");
					sprintf(testfile_name+strlen(testfile_name), gen_fmt, 1, 2+tti);
					autozip_fp tfp;
					int resop = autozip_open(testfile_name, &tfp);
					//printf("%d === %s    %s\n", resop, gen_fmt, testfile_name);
					if(0<=resop){
						scancon -> bcl_is_gzipped = resop;
						autozip_close(&tfp);
					}else bcl_found=0;
				}
				if(bcl_found){
					char * gen_fmt2 = str_replace(data_dir , "/L001", "/L%03d");
					strcpy(scancon -> out_format_string, gen_fmt2);	
					free(gen_fmt2);
					strcat(scancon -> out_format_string, "/");
					strcat(scancon -> out_format_string, gen_fmt);
				}

				free(gen_fmt);
			}
		}
	}

	if(bcl_found && filter_found){
		scancon -> found_answer = 1;
		scancon -> filter_is_gzipped = filter_found - 1;
	}

	closedir(this_level);
	return 0;
}

int iBLC_guess_format_string(char * data_dir, int * cluster_bases, char * format_string, char * filter_format, int * bcl_is_gzipped, int * filter_is_gzipped, int * read_lens, int * is_index){
	struct iBLC_scan_t sct;
	memset(&sct, 0, sizeof(sct));
	int tii = iBLC_guess_scan(&sct, data_dir);

	if(tii || ! sct.found_answer) return -1;
	strcpy(format_string, sct.out_format_string);
	strcpy(filter_format, sct.filter_format_string);
	*filter_is_gzipped = sct.filter_is_gzipped;
	*bcl_is_gzipped = sct.bcl_is_gzipped;
	*cluster_bases=0;

	for(tii=0; tii<sct.reads_per_cluster; tii++){
		if(sct.read_lengths[tii]<1) return -1;
		read_lens[tii] = sct.read_lengths[tii];
		is_index[tii] = sct.read_is_index[tii];
		(*cluster_bases) += sct.read_lengths[tii];
		read_lens[tii+1]=0;
	}
		
	return 0;
}

void iBLC_close_batch(input_BLC_t * blc_input){
	int ii;
	if(blc_input->is_EOF) return;
	if(NULL == blc_input -> bcl_gzip_fps && blc_input -> bcl_is_gzipped)return;
	if(NULL == blc_input -> bcl_fps && !blc_input -> bcl_is_gzipped)return;
	for(ii=0; ii < blc_input->total_bases_in_each_cluster; ii++){
		if(blc_input -> bcl_is_gzipped){
			seekgz_close(blc_input -> bcl_gzip_fps[ii]);
			free(blc_input -> bcl_gzip_fps[ii]);
			blc_input -> bcl_gzip_fps[ii] = NULL;
		}else{
			fclose(blc_input -> bcl_fps[ii]);
			blc_input -> bcl_fps[ii] = NULL;
		}
	}
	if(blc_input -> filter_is_gzipped){
		seekgz_close(blc_input -> filter_gzip_fp);
		free(blc_input -> filter_gzip_fp);
		blc_input -> filter_gzip_fp = NULL;
	}else{
		fclose(blc_input -> filter_fp);
		blc_input -> filter_fp = NULL;
	}

	if( blc_input -> bcl_is_gzipped ){
		free(blc_input -> bcl_gzip_fps);
		blc_input -> bcl_gzip_fps = NULL;
	}else{
		free(blc_input -> bcl_fps);
		blc_input -> bcl_fps = NULL;
	}
}

int iBLC_open_batch(input_BLC_t * blc_input ){
	char fname[MAX_FILE_NAME_LENGTH];
	iBLC_close_batch(blc_input);
	int fii, xx;
	blc_input -> is_EOF=1;

	if(blc_input -> bcl_gzip_fps == NULL) blc_input -> bcl_gzip_fps = calloc( sizeof(void *), blc_input -> total_bases_in_each_cluster ); // for both FILE** and seekgz **
	for(fii = 0; fii < blc_input -> total_bases_in_each_cluster; fii++){
		sprintf(fname, blc_input -> bcl_format_string, blc_input -> current_lane, fii+1);
		if(blc_input -> bcl_is_gzipped){
			blc_input -> bcl_gzip_fps[fii] = calloc( sizeof(seekable_zfile_t), 1);
			int rv = seekgz_open(fname, blc_input -> bcl_gzip_fps[fii], NULL);
			if(rv){
				SUBREADprintf("ERROR: Unable to open %s\n", fname);
				return -1;
			}
			for(xx = 0; xx < 4; xx++) seekgz_next_int8(blc_input -> bcl_gzip_fps[fii]); // skip the first 32-b integer
		}else{
			blc_input -> bcl_fps[fii] = fopen(fname, "r");
			 if(NULL == blc_input -> bcl_fps[fii]){
				SUBREADprintf("ERROR: Unable to open %s\n", fname);
				return -1;
			}
			for(xx = 0; xx < 4; xx++) fgetc(blc_input -> bcl_fps[fii]); // skip the first 32-b integer
		}
	}

	sprintf(fname, blc_input -> filter_format_string, blc_input -> current_lane,blc_input -> current_lane);
	if(blc_input -> filter_is_gzipped){
		blc_input -> filter_gzip_fp = calloc( sizeof(seekable_zfile_t), 1);
		int rv = seekgz_open(fname, blc_input -> filter_gzip_fp, NULL);
		if(rv){
			SUBREADprintf("ERROR: Unable to open %s\n", fname);
			return -1;
		}
		for(xx = 0; xx < 12; xx++) seekgz_next_int8(blc_input -> filter_gzip_fp); // skip the 12-byte header
	}else{
		blc_input -> filter_fp = fopen(fname, "r");
		if(NULL == blc_input -> filter_fp){
			SUBREADprintf("ERROR: Unable to open %s\n", fname);
			return -1;
		}
		for(xx = 0; xx < 12; xx++) fgetc(blc_input -> filter_fp); // skip the 12-byte header
	}
	blc_input -> is_EOF=0;
	return 0;
}

int iCache_open_batch( cache_BCL_t * cache_input){
	cache_input -> bcl_gzip_fps = calloc(sizeof(autozip_fp), cache_input -> total_bases_in_each_cluster);
	return 0;
}

int cacheBCL_go_chunk_end( cache_BCL_t * cache_input ){
	cache_input -> read_no_in_chunk= cache_input -> reads_available_in_chunk;
	return 0;
}
int cacheBCL_go_chunk_start( cache_BCL_t * cache_input ){
	cache_input -> read_no_in_chunk=0;
	return 0;
}
void cacheBCL_close(cache_BCL_t * cache_input){
	int x1;
	for(x1=0;x1<cache_input -> total_bases_in_each_cluster;x1++){
		if( cache_input -> bcl_gzip_fps[x1]. plain_fp || cache_input -> bcl_gzip_fps[x1].gz_fp.gz_fp )autozip_close(&cache_input -> bcl_gzip_fps[x1]);
		free(cache_input -> bcl_bin_cache[x1]);
	}
	free(cache_input -> bcl_gzip_fps);
	if(cache_input -> filter_fp. plain_fp || cache_input -> filter_fp.gz_fp.gz_fp )autozip_close(&cache_input -> filter_fp);
	free(cache_input -> lane_no_in_chunk);
	free(cache_input -> flt_bin_cache);
}

int cacheBCL_init( cache_BCL_t * cache_input, char * data_dir, int reads_in_chunk, int all_threads ){
	memset(cache_input, 0, sizeof( cache_BCL_t));
	subread_init_lock(&cache_input -> read_lock);
	int rv = iBLC_guess_format_string(data_dir, &cache_input -> total_bases_in_each_cluster, cache_input -> bcl_format_string, cache_input -> filter_format_string, &cache_input -> bcl_is_gzipped, &cache_input -> filter_is_gzipped, cache_input -> single_read_lengths, cache_input -> single_read_is_index);
	if(rv) return -1;
	cache_input -> current_lane = 1;
	cache_input -> reads_per_chunk = reads_in_chunk;
	cache_input -> bcl_bin_cache = malloc(sizeof(char*) * cache_input -> total_bases_in_each_cluster);
	int x1;
	for(x1 = 0; x1 < cache_input -> total_bases_in_each_cluster; x1++)
		cache_input -> bcl_bin_cache[x1] = malloc(reads_in_chunk);
	cache_input -> flt_bin_cache = malloc(reads_in_chunk*2);
	cache_input -> flt_bin_cache_size = reads_in_chunk*2;
	cache_input -> lane_no_in_chunk = malloc(reads_in_chunk);
	cache_input -> chunk_end_lane = 1;	// no "0th lane"
	cache_input -> all_threads = all_threads;
	return iCache_open_batch(cache_input)?1:0;
}


void iCache_close_one_fp( cache_BCL_t * cache_input, int bcl_no){
	autozip_fp * tfp = bcl_no<0?&cache_input-> filter_fp :(cache_input -> bcl_gzip_fps + bcl_no);
//	SUBREADprintf("CLOSE_AUTO: %d\n", bcl_no);
	autozip_close(tfp);
	memset(tfp,0,sizeof(autozip_fp));
}

int iCache_open_one_fp( cache_BCL_t * cache_input, int bcl_no, int lane_no){
	autozip_fp * tfp = bcl_no<0?&cache_input-> filter_fp :(cache_input -> bcl_gzip_fps + bcl_no);
	char fname[MAX_FILE_NAME_LENGTH+1];
	if(bcl_no <0)
		sprintf(fname,  cache_input -> filter_format_string, lane_no, lane_no);
	else
		sprintf(fname,  cache_input -> bcl_format_string, lane_no, bcl_no+1);

	int rv = autozip_open(fname, tfp);
//	SUBREADprintf("OPEN_AUTO %s = %d\n", fname, rv);
	if(rv>=0){
		int sk = bcl_no<0?12:4;
		for(; sk>0; sk--){
			autozip_getch(tfp);
			//SUBREADprintf("SKEP=%d\n", nch);
		}
	}else{
		memset(tfp,0,sizeof(autozip_fp));
	}
	return rv < 0;
}

// it returns the number of bytes loaded. 0 if no bytes are available.
int iCache_continuous_read_lanes( cache_BCL_t * cache_input, int bcl_no){
	autozip_fp * tfp = bcl_no<0?&cache_input-> filter_fp :(cache_input -> bcl_gzip_fps + bcl_no);
	int my_lane = cache_input -> chunk_start_lane ,wptr=0;
	char * wpt =  bcl_no<0?cache_input-> flt_bin_cache:cache_input -> bcl_bin_cache[bcl_no];
	int total_valid_reads = 0;
	int raw_ptr = 0;
	while(1){
		if(!(tfp -> plain_fp || tfp -> gz_fp.gz_fp)){
			int fpr = iCache_open_one_fp( cache_input, bcl_no, my_lane );
			if(fpr){
				if(bcl_no<0)cache_input -> last_chunk_in_cache = 1;
				break;
			}
		}
		while(1){
			int nch = autozip_getch( tfp );
			//if(bcl_no>=0)SUBREADprintf("NCH=%d\n", nch);
			if(nch>=0){
				if( bcl_no < 0 || cache_input -> flt_bin_cache [raw_ptr] ){
					if(bcl_no <0){
						if(nch>0) cache_input -> lane_no_in_chunk[total_valid_reads ++] = my_lane;
						if(wptr == cache_input ->  flt_bin_cache_size){
							cache_input ->  flt_bin_cache_size *= 1.6;
							wpt = cache_input-> flt_bin_cache = realloc( wpt, cache_input ->  flt_bin_cache_size );
						}
					}else total_valid_reads++;

					//if(wptr == 1000000)SUBREADprintf("BaseNo=%d; Nch=%d\n", bcl_no, nch);
					wpt[wptr++] = nch&0xff;
					if(total_valid_reads == cache_input -> reads_per_chunk) break;
				}
			}else{
//				SUBREADprintf("CONTINUOUSLY [%d-th base] BREAK : %d\n", bcl_no, nch);
				break;
			}
			raw_ptr ++;
		}
		if(total_valid_reads == cache_input -> reads_per_chunk) break;
		iCache_close_one_fp(cache_input, bcl_no);
		my_lane++;
	}
	if(bcl_no <0){
		cache_input -> reads_available_in_chunk = total_valid_reads;
		cache_input -> chunk_end_lane = my_lane;
	}
	//if(bcl_no<0) SUBREADprintf("CONTINUOUSLY [%d-th base] READ %d reads into the %d-th lanel fpntr = %ld\n", bcl_no,  total_valid_reads, my_lane, tfp -> plain_fp ? ftello(tfp -> plain_fp ): -1);
	return total_valid_reads;
}

void * iCache_decompress_chunk_1T(void * arg){
	cache_BCL_t * cache_input = arg;
	while (1){
		int my_bcl_no ;
		subread_lock_occupy(&cache_input -> read_lock);
		for(my_bcl_no =0; my_bcl_no<cache_input -> total_bases_in_each_cluster; my_bcl_no++){
			if(!cache_input -> bcl_no_is_used[my_bcl_no]) {
				cache_input -> bcl_no_is_used[my_bcl_no]=1;
				break;
			}
		}
		subread_lock_release(&cache_input -> read_lock);
		if(my_bcl_no>= cache_input -> total_bases_in_each_cluster) return NULL;

		iCache_continuous_read_lanes( cache_input, my_bcl_no );
	}
	return NULL;
}

int cacheBCL_next_chunk(cache_BCL_t * cache_input){
	int x1;
	//SUBREADprintf("BCL: READ_CHUNK_NEXT ( %d )\n", cache_input -> chunk_no+1);
	cache_input -> chunk_start_lane = cache_input -> chunk_end_lane;
	memset(cache_input -> bcl_no_is_used, 0 , sizeof(int)* MAX_READ_LENGTH);
	pthread_t * threads = malloc(sizeof(pthread_t)*cache_input -> all_threads);
	iCache_continuous_read_lanes( cache_input, -1 ); // read filtering binary

	for(x1=0; x1<cache_input -> all_threads; x1++)
		pthread_create(threads+x1, NULL, iCache_decompress_chunk_1T, cache_input);

	for(x1=0; x1<cache_input -> all_threads; x1++)
		pthread_join(threads[x1],NULL);
	free(threads);
	cache_input -> read_no_in_chunk = 0;
	cache_input -> chunk_no ++;
	return 0;
}


int iCache_copy_read(cache_BCL_t * cache_input, char * read_name, char * seq, char * qual, long long rno){
	int bii, idx_offset, base_offset;
	int * srii = cache_input -> single_read_lengths;

	sprintf(read_name, "R%011llu:", rno);
	idx_offset  = srii[0];
	base_offset = srii[1] + idx_offset;

	read_name[13+idx_offset]='|';
	read_name[14+2*idx_offset]='|';
	read_name[15+base_offset+idx_offset]='|';
	sprintf(read_name +16 +2*base_offset, "|L%03d" , cache_input -> lane_no_in_chunk[cache_input -> read_no_in_chunk]);

	for(bii = 0; bii < cache_input -> total_bases_in_each_cluster; bii++){
		int nch = cache_input -> bcl_bin_cache[bii][cache_input -> read_no_in_chunk];
		//if(rno == 1000000){
		//	SUBREADprintf("COPY: %d base NCH=%d\n", bii, nch);
		//}
		if(nch<0) nch+=256;
		int nbase = 'N';
		int nqual = '#';
		if(nch > 0){
			nbase="ACGT"[nch%4];
			nqual=33+(nch>>2);
		}
		if(nqual >= '/' && bii < srii[0] + srii[1]) nqual++;
		if(bii < srii[0]){
			read_name[13+bii] = nbase;
			read_name[14+idx_offset+bii]= nqual;
		}else if(bii < srii[0] + srii[1]){
			read_name[15+idx_offset+bii] = nbase;
			read_name[16+base_offset+bii]= nqual;
		}else{
			seq[bii - srii[0]-srii[1] ] = nbase;
			qual[bii - srii[0]-srii[1] ] = nqual;
		}
	}

	cache_input -> read_no_in_chunk++;
	//SUBREADprintf("GOT READ #%d ; ret=%d\n" , cache_input -> read_no_in_chunk, srii[2]);
	return srii[2];
}



int cacheBCL_next_read(cache_BCL_t * cache_input, char * read_name, char * seq, char * qual, long long * read_number_in_all){
	if(cache_input -> read_no_in_chunk >= cache_input -> reads_available_in_chunk){
		if(cache_input -> last_chunk_in_cache) return 0;
		cacheBCL_next_chunk(cache_input);
		if(cache_input -> read_no_in_chunk >= cache_input -> reads_available_in_chunk) // no reads are loaded in the previous step
			return 0;
	}

	long long rnumb =(cache_input -> chunk_no -1)*1ll * cache_input -> reads_per_chunk +(cache_input -> read_no_in_chunk);
	if(read_number_in_all) *read_number_in_all = rnumb;
	return iCache_copy_read(cache_input, read_name, seq, qual, rnumb);
}

int input_BLC_init( input_BLC_t * blc_input , char * data_dir ){
	memset(blc_input, 0, sizeof(input_BLC_t));
	subread_init_lock(&blc_input -> read_lock);

	int rv = iBLC_guess_format_string(data_dir, &blc_input -> total_bases_in_each_cluster, blc_input -> bcl_format_string, blc_input -> filter_format_string, &blc_input -> bcl_is_gzipped, &blc_input -> filter_is_gzipped, blc_input -> single_read_lengths, blc_input -> single_read_is_index);
	if(rv) return -1;

	blc_input -> current_lane = 1;

	return iBLC_open_batch(blc_input)?1:0;
}
// load the next read W/O switch lane. 
// Return 0 if EOF, -1 if error or bases if loaded correctly.
int iBLC_current_lane_next_read(input_BLC_t * blc_input, char * readname , char * read, char * qual){
	int bii, idx_offset, base_offset;

	sprintf(readname, "R%011llu:", blc_input -> read_number +1);

	{
		idx_offset = blc_input -> single_read_lengths[0];
		base_offset = idx_offset + blc_input -> single_read_lengths[1];
	}

	readname[13+idx_offset]='|';
	readname[14+2*idx_offset]='|';
	readname[15+base_offset+idx_offset]='|';
	sprintf(readname +16 +2*base_offset, "|L%03d" , blc_input -> current_lane);

	while(1){
		int fch = blc_input -> filter_is_gzipped? seekgz_next_int8(blc_input -> filter_gzip_fp) :fgetc(blc_input -> filter_fp);
		if(fch < 0) return 0;
		int baseii =0;
		for(bii =0; bii< blc_input -> total_bases_in_each_cluster; bii++){
			int nch = blc_input -> bcl_is_gzipped?seekgz_next_int8(blc_input -> bcl_gzip_fps[bii]):fgetc(blc_input -> bcl_fps[bii]), bv, qv;
			assert(nch >=0 && nch <=255);

			if(0==nch){
				bv='N';
				qv='#';
			}else{
				bv="ACGT"[nch%4];
				qv=33+(nch>>2);
				if(qv >= '/') qv++;
			}
			if(bii < idx_offset){
				assert(bv !=0 && qv !=0);
				readname[13+ bii]=bv;
				readname[14+ bii+idx_offset]=qv;
			}else if(bii < base_offset){
				assert(bv !=0 && qv !=0);
				readname[15+ bii+ idx_offset]=bv;
				readname[16+ bii+ base_offset]=qv;
			}else{
				read[baseii] = bv;
				qual[baseii] = qv;
				baseii++;
			}
		}
		assert(fch==1||fch==0);
		if(fch==1){
		//	SUBREADprintf("ASS_RNAME=%s\n", readname);
			blc_input -> read_number ++;
			return baseii;
		}
	}
}

int iBLC_inc_lane(input_BLC_t * blc_input){
	blc_input -> current_lane ++;
	return iBLC_open_batch(blc_input); // this function automatically closes BCL fps and FILTER fp.
}

// return : -1: error, 0: end of all files and all lanes, >0: actual read is loaded (return the read len). The read name is the combination of the short-end and the index-end.
// NOTE: this only works with scRNA protocol!!
int input_BLC_next_read(input_BLC_t * blc_input , char * readname, char * read, char * qual){
	int nextlane, rrv=0;
	if(blc_input->is_EOF) return 0;

	subread_lock_occupy(&blc_input -> read_lock);
	for(nextlane = 0; nextlane <2; nextlane++){
		int rv = iBLC_current_lane_next_read(blc_input, readname, read, qual);
		if(rv >0 || rv <0){
			rrv = rv;
			break;
		}
		if(rv ==0 && nextlane){
			rrv = 0;
			break;
		}
		if(nextlane>0){
			rrv = -1;
			break;
		}
		
		rv = iBLC_inc_lane(blc_input);
		if(rv){
			rrv = 0;
			break;
		}
	}
	subread_lock_release(&blc_input -> read_lock);
	return rrv;
}

int input_BLC_tell ( input_BLC_t * blc_input , input_BLC_pos_t * pos ){
	int xx1;
	memset(pos,0, sizeof(*pos));
	pos -> lane_id = blc_input -> current_lane;
	pos -> read_number = blc_input -> read_number;
	pos -> is_EOF = blc_input -> is_EOF;
	if(pos->is_EOF) return 0;
	if(blc_input->bcl_is_gzipped){
		pos -> pos_of_bclgzs = calloc(sizeof(void *) , blc_input -> total_bases_in_each_cluster);
		for(xx1=0; xx1<blc_input -> total_bases_in_each_cluster; xx1++){
			pos -> pos_of_bclgzs[xx1] = malloc(sizeof(seekable_position_t));
			seekgz_tell(blc_input->bcl_gzip_fps[xx1], pos -> pos_of_bclgzs[xx1]);
		}
	}else{
		pos -> pos_of_bcls = calloc(sizeof(long long) , blc_input -> total_bases_in_each_cluster);
		for(xx1=0; xx1<blc_input -> total_bases_in_each_cluster; xx1++)
			pos -> pos_of_bcls[xx1] = ftello(blc_input->bcl_fps[xx1]);
	}


	if(blc_input->filter_is_gzipped){
		pos -> pos_of_filtergz = malloc(sizeof(seekable_position_t));
		seekgz_tell(blc_input->filter_gzip_fp, pos -> pos_of_filtergz);
	}else pos -> pos_of_filter = ftello(blc_input->filter_fp);

	return 0;
}

int input_BLC_seek( input_BLC_t * blc_input , input_BLC_pos_t * pos ){
	int xx1;
	blc_input -> read_number = pos -> read_number;
	if(pos -> is_EOF){
		iBLC_close_batch(blc_input);
		blc_input -> is_EOF = pos -> is_EOF;
		blc_input -> current_lane = pos -> lane_id;
		return 0;
	}

	if(pos -> lane_id != blc_input -> current_lane){
		blc_input -> current_lane = pos -> lane_id;
		iBLC_open_batch(blc_input);
	}

	for(xx1=0; xx1<blc_input -> total_bases_in_each_cluster; xx1++)
		if(blc_input->bcl_is_gzipped) seekgz_seek(blc_input->bcl_gzip_fps[xx1], pos -> pos_of_bclgzs[xx1]); 
		else fseeko(blc_input->bcl_fps[xx1], pos -> pos_of_bcls[xx1], SEEK_SET);

	if(blc_input->filter_is_gzipped) seekgz_seek(blc_input->filter_gzip_fp, pos -> pos_of_filtergz);
	else fseeko(blc_input->filter_fp, pos -> pos_of_filter, SEEK_SET);
	
	return 0;
}

void input_BLC_destroy_pos(input_BLC_t * blc_input , input_BLC_pos_t *pos){
	int xx1;
	for(xx1=0; xx1<blc_input -> total_bases_in_each_cluster; xx1++){
		if(blc_input->bcl_is_gzipped) free(pos -> pos_of_bclgzs[xx1]);
	}
	free((blc_input->bcl_is_gzipped? (void*)pos -> pos_of_bclgzs:(void*)pos -> pos_of_bcls));
}

void input_BLC_close(input_BLC_t * blc_input){
	iBLC_close_batch(blc_input);
	subread_destroy_lock(&blc_input -> read_lock);
}

void iBLC_free_sample_items(void * sample_arr){
	ArrayList * ar = (ArrayList * ) sample_arr;
	ArrayListDestroy(ar);
}

void iBLC_free_3tp(void * t){
	char **tt = (char**)t;
	free(tt[1]);
	free(tt);
}

HashTable * input_BLC_parse_SampleSheet(char * fname){
	HashTable * ret = StringTableCreate(30);
	HashTableSetDeallocationFunctions(ret, free, iBLC_free_sample_items);
	FILE * fp = fopen(fname, "r");
	if(fp==NULL) return NULL;
	char linebuf[MAX_FILE_NAME_LENGTH];
	int state = -1;
	while(!feof(fp)){
		fgets(linebuf, MAX_FILE_NAME_LENGTH-1, fp);
		if(strlen(linebuf)<5)continue;
		if(state < 0 && strstr(linebuf,"EMFileVersion,4")) state = 0;
		if(state == 1 && linebuf[0]=='[') state = 99999;
		if(state == 1){
			if(memcmp( linebuf, "Lane", 4 )==0)continue;

			char * tokp=NULL;
			int lane_no = atoi(strtok_r(linebuf, ",", &tokp));
			strtok_r(NULL, ",", &tokp);
			char * sample_name = strtok_r(NULL, ",", &tokp);
			char * sample_index = strdup(strtok_r(NULL, ",", &tokp));
			char ** entry = malloc(sizeof(void*)*3);
			entry[0] = NULL + lane_no;
			entry[1] = sample_index;

			ArrayList * arr = HashTableGet(ret, sample_name);
			if(NULL == arr){
				arr = ArrayListCreate(16);
				ArrayListSetDeallocationFunction(arr, iBLC_free_3tp);
				//SUBREADprintf("PUT_SAMPLE=%s\n", sample_name);
				HashTablePut( ret, strdup(sample_name), arr );
			}
			ArrayListPush(arr,entry);
		}
		if(state == 0 && strstr(linebuf,"ata]"))state = 1;
	}
	if(state <1){
		SUBREADprintf("ERROR: the sample sheet doesn't contain any sample.\n");
		return NULL;
	}
	return ret;
}

ArrayList * input_BLC_parse_CellBarcodes(char * fname){
	autozip_fp fp;
	int resop = autozip_open(fname, &fp);
	if(resop<0) return NULL;

	ArrayList * ret = ArrayListCreate(10000000);
	ArrayListSetDeallocationFunction(ret, free);

	while(1){
		char tmp_fl[MAX_BARCODE_LEN+1];
		int skr = autozip_gets(&fp, tmp_fl, MAX_BARCODE_LEN);
		if(skr<1) break;
		int x1;
		if(tmp_fl[skr-1]=='\n') tmp_fl[skr-1]=0;
		for(x1=0; tmp_fl[x1]; x1++) if(!isalpha(tmp_fl[x1])){
			tmp_fl[x1]=0;
			break;
		}
		ArrayListPush(ret, strdup(tmp_fl));
	}

	autozip_close(&fp);
	return ret;
}

#ifdef MAKE_TEST_SAMPLESHEET
int main(int argc, char ** argv){
	assert(argc>1);
	char * samplesheet = argv[1];
	HashTable * st = input_BLC_parse_SampleSheet(samplesheet);
	printf("P=%p\n", st);
}
#endif
#ifdef MAKE_TEST_ICACHE
int main(int argc, char ** argv){
	assert(argc>1);
	cache_BCL_t blc_input;
	int orv = cacheBCL_init(&blc_input, argv[1],5000000, 5), total_poses = 0;
	printf("orv=%d, bases=%d, filter_gzip=%d, data_gzip=%d\nBCL-pattern = %s\nFLT-pattern=%s\n", orv, blc_input.total_bases_in_each_cluster, blc_input.filter_is_gzipped, blc_input.bcl_is_gzipped, blc_input.bcl_format_string, blc_input.filter_format_string);

	while(1){
		char base[1000], qual[1000], rname[200];
		long long readno = 0;
		base[0]=qual[0]=rname[0]=0;
		orv = cacheBCL_next_read(&blc_input, rname, base, qual, &readno);
		assert(orv>=0);
		if(0==orv) break;
		if(1||readno%1000000==0)printf("%s %s %s\n", rname, base, qual);
	}

	cacheBCL_close(&blc_input);
}
#endif

#ifdef MAKE_TEST_IBLC
int main(int argc, char ** argv){
	assert(argc>1);
	input_BLC_t blc_input;
	input_BLC_pos_t poses[12];
	int orv = input_BLC_init(&blc_input, argv[1]), total_poses = 0;
	printf("orv=%d, bases=%d, filter_gzip=%d, data_gzip=%d\n", orv, blc_input.total_bases_in_each_cluster, blc_input.filter_is_gzipped, blc_input.bcl_is_gzipped);

	while(1){
		char base[1000], qual[1000], rname[200];
		base[0]=qual[0]=rname[0]=0;
		orv = input_BLC_next_read(&blc_input, rname, base, qual);
		assert(orv>=0);
		if(0==orv) break;
		if(blc_input.read_number%1000000==1)printf("%s %s %s\n", rname, base, qual);
		if(blc_input.read_number % 20000000 == 1){
			input_BLC_tell(&blc_input, poses + total_poses);
			total_poses++;
		}
		if(blc_input.read_number > 180000010) break;
	}

	int ii;
	for(ii = 0; ii < total_poses; ii++){
		long long jj;
		printf("\n\n========== SEEKING %d OF %d ================\n", ii,total_poses);
		input_BLC_seek( &blc_input, poses+ii );
		input_BLC_destroy_pos( &blc_input, poses+ii );
		for(jj = 0 ; jj < 3000011; jj++){
			char base[1000], qual[1000], rname[200];
			orv = input_BLC_next_read(&blc_input, rname, base, qual);
			assert(orv>=0);
			if(blc_input.read_number%1000000==2)printf("%s %s %s\n", rname, base, qual);
		}
	} 
	
	printf("END CORRECTLY  R=%llu\n", blc_input.read_number);
	input_BLC_close(&blc_input);
	return 0;
}
#endif
