//#define MAKE_CELLCOUNTS
#include <stdio.h>

#include <string.h>
#include <ctype.h>
#include <fcntl.h>
#include <dirent.h>
#include <limits.h>
#include <locale.h>
#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <sys/stat.h>
#include "subread.h"
#include "gene-algorithms.h"
#include "HelperFunctions.h"
#include "input-files.h"
#include "interval_merge.h"
#include "input-blc.h"
#include "core-indel.h"
#include "cell-counts.h"
#include "core-junction.h"
#include "softclip-test.h"

//#define DO_STARSOLO_THING
#define IMPOSSIBLE_MEMORY_SPACE 0x5CAFEBABE0000000llu

#define CCT_NIL_TXN_PLACEHOLDER "__DBPZ_com_nil_TXN" // SAF format: using placeholder for transcripts. DBPZ_com is my namespace.
#define CCT_GENE_ID_UMI_PREFIX "__DBPZ_com_g_ID_want_"
#define REVERSE_TABLE_BUCKET_LENGTH 131072

int has_reverse_table_reported =0;

void cellCounts_junckey_sort_exchange(void * inptr, int i, int j){

	char ** inp = (char **) inptr;
	char * tmpp = inp[j];
	inp[j]=inp[i];
	inp[i]=tmpp;
}


int cellCounts_junckey_sort_compare(void * inptr, int i, int j){
	char ** inp = (char **) inptr;
	int x1;

	int chrI=-1, chrJ=-1;

	if(atoi(inp[i])>0) chrI = atoi(inp[i]);
	if(atoi(inp[j])>0) chrJ = atoi(inp[j]);

	if(inp[i][0]=='X' && !isdigit(inp[i][1])&& !isalpha(inp[i][1])) chrI = 90;
	if(inp[i][0]=='Y' && !isdigit(inp[i][1])&& !isalpha(inp[i][1])) chrI = 91;
	if(inp[i][0]=='M' && !isdigit(inp[i][1])&& !isalpha(inp[i][1])) chrI = 99;
	if(inp[j][0]=='X' && !isdigit(inp[j][1])&& !isalpha(inp[j][1])) chrJ = 90;
	if(inp[j][0]=='Y' && !isdigit(inp[j][1])&& !isalpha(inp[j][1])) chrJ = 91;
	if(inp[j][0]=='M' && !isdigit(inp[j][1])&& !isalpha(inp[j][1])) chrJ = 99;



	if(memcmp(inp[i], "chr", 3)==0){
		chrI=atoi(inp[i]+3);
		if(0 == chrI && inp[i][3] == 'X') chrI = 90;
		if(0 == chrI && inp[i][3] == 'Y') chrI = 91;
		if(0 == chrI && inp[i][3] == 'M') chrI = 99;
	}
	if(memcmp(inp[j], "chr", 3)==0){
		chrJ=atoi(inp[j]+3);
		if(0 == chrJ && inp[j][3] == 'X') chrJ = 90;
		if(0 == chrJ && inp[j][3] == 'Y') chrJ = 91;
		if(0 == chrJ && inp[j][3] == 'M') chrJ = 99;
	}

	int len_I_long = 9;
	for(x1 = 0 ; x1 < FEATURE_NAME_LENGTH + 15 ; x1++){
		int c1 = inp[i][x1];
		int c2 = inp[j][x1];
		if(c1 == '\t' && c2 != '\t')
			len_I_long = -1;
		else if(c1 != '\t' && c2 == '\t')
			len_I_long = 1;
		else if(c1 == '\t' && c2 == '\t')
			len_I_long = 0;

		if(len_I_long != 9) break;
	}

	if(chrI != chrJ || len_I_long != 0){
		return (chrI * 100 + len_I_long) - (chrJ * 100);
	}

	for(x1 = 0 ; x1 < FEATURE_NAME_LENGTH + 15 ; x1++){
		int c1 = inp[i][x1];
		int c2 = inp[j][x1];
		if(c1 != c2){
			return c1 - c2;
		}else if(c1 == '\t' && c1 == c2){
			int pos1 = atoi(inp[i]+x1+1);
			int pos2 = atoi(inp[j]+x1+1);
			if( pos1 == pos2)
				return strcmp(inp[i], inp[j]);
			else
				return pos1 - pos2;
		}

		if(c1 == 0 || c2 == 0)return c1 - c2;
	}
	return 0;
}




void cellCounts_junckey_sort_merge(void * inptr, int start, int items1, int items2){
	char ** inp = (char **) inptr;
	char ** tmpp = malloc(sizeof(char *) * (items1+items2));
	int read_1_ptr = start, read_2_ptr = start+items1, outptr = 0;
	while(1){
		if(read_1_ptr == start+items1 && read_2_ptr == start+items1+items2) break;
		if((read_1_ptr == start+items1)||(read_2_ptr < start+items1+items2 && cellCounts_junckey_sort_compare(inptr, read_1_ptr, read_2_ptr) > 0 )) {
			// select 2
			tmpp[outptr++]=inp[read_2_ptr++];
		} else {
			// select 1
			tmpp[outptr++]=inp[read_1_ptr++];
		}
	}
	memcpy(inp + start, tmpp, sizeof(char *)*(items1+items2));
	free(tmpp);
}


int cellCounts_reduce_Cigar(char * cigar, char * cigarout, int * sum_ins){
	//LRMtest_move_buff( context, thread_context, iteration_context, thread_context -> dynamic_programming_indel_movement_buf, strlen(thread_context -> dynamic_programming_indel_movement_buf), iteration_context -> read_length);
	int tmpi = -1;
	int ci, nch, repeat_i = 0, old_opt = 0, wcur=0, rlen=0;
	for(ci = 0; ; ci++){
		nch = cigar[ci];
		if(!nch) break;
		if(isdigit(nch)){
			if(tmpi<0) tmpi = 0;
			tmpi = tmpi*10 + (nch-'0');
		}else{
			if(tmpi<0) tmpi = 1;
			if(old_opt != nch && repeat_i>0){
				if(old_opt=='M' || old_opt=='S' || old_opt=='I') rlen+=repeat_i;
				wcur += SUBreadSprintf( cigarout + wcur, 11, "%d%c", repeat_i, old_opt );
				repeat_i = 0;
			}
			repeat_i += tmpi;
			tmpi = -1;
			old_opt = nch;
		}
	}
	if(sum_ins) *(sum_ins)=0;
	if(repeat_i>0){
		SUBreadSprintf(cigarout + wcur, 11, "%d%c", repeat_i, old_opt);
		if(old_opt=='I'&&sum_ins) *(sum_ins)+= repeat_i;
		if(old_opt=='M' || old_opt=='S' || old_opt=='I') rlen+=repeat_i;
	}
	return rlen;
}

void cellCounts_cell_barcode_tabel_destroy(void *a){
	if(((a-NULL) & 0xfffffffff0000000llu ) ==IMPOSSIBLE_MEMORY_SPACE )return;
	ArrayListDestroy((ArrayList*)a);
}

int cellCounts_make_barcode_HT_table(cellcounts_global_t * cct_context){
	int xx1,xx2;
	cct_context -> cell_barcode_head_tail_table = StringTableCreate(600000);
	HashTableSetDeallocationFunctions(cct_context -> cell_barcode_head_tail_table, free, cellCounts_cell_barcode_tabel_destroy);

	for(xx1=0;xx1 < cct_context-> cell_barcodes_array -> numOfElements; xx1++){
		char * bc = ArrayListGet(cct_context-> cell_barcodes_array, xx1);
		int bcl =strlen(bc);
		if(cct_context->visium_hd_barcodes) cct_context -> known_cell_barcode_length= -1;
		else{
			if(cct_context -> known_cell_barcode_length==0) cct_context -> known_cell_barcode_length=bcl;
			if(bcl!=cct_context -> known_cell_barcode_length){
				SUBREADprintf("ERROR: the cell barcode list must contain equal-length strings!\n");
				return 1;
			}
		}

		char bctmp[20];
		HashTablePut(cct_context -> cell_barcode_head_tail_table, strdup(bc), NULL+xx1+IMPOSSIBLE_MEMORY_SPACE);
		for(xx2=0; xx2<2; xx2++){
			bctmp[0] = xx2?'S':'F';
			int xx3;
			int applied_bc_len = cct_context->visium_hd_barcodes?MIN_LEN_VISIUM_HD_CELLBC:bcl;
			for(xx3 = 0; xx3< applied_bc_len/2; xx3++)
				bctmp[xx3+1] = bc[ xx3*2+xx2 ];
			bctmp[applied_bc_len/2+1]=0;

			ArrayList * array_of_codes = HashTableGet(cct_context -> cell_barcode_head_tail_table, bctmp);
			if(!array_of_codes){
				array_of_codes = ArrayListCreate(4);
				HashTablePut(cct_context -> cell_barcode_head_tail_table, strdup(bctmp), array_of_codes);
			}
			ArrayListPush(array_of_codes, NULL+xx1);
		}
	}
	return 0;
}


void cellCounts_absoffset_to_posstr(cellcounts_global_t * cct_context, unsigned int pos, char * res);

#define SOFT_CLIPPING_WINDOW_SIZE 5
#define SOFT_CLIPPING_MAX_ERROR   1
#define gvindex_baseno2offset_m(base_number, index, offset_byte, offset_bit)    {offset_byte =  ((base_number) - index -> start_base_offset) >>2; offset_bit = (base_number) % 4 * 2;}
int cellCounts_get_index_int(gene_value_index_t * value_index, unsigned int pos){
	int offset_byte , offset_bit;
	gvindex_baseno2offset_m(pos, value_index,  offset_byte , offset_bit);
	return (value_index ->values [offset_byte] >> offset_bit)&3;
}

int cellCounts_get_read_int(char * read_bin, int read_offset){
	int read_byte , read_bit ;
	read_byte = read_offset/4;
	read_bit = read_offset%4 *2;
	return (read_bin[read_byte]>>read_bit)&3;
}

// it returns the number of bases to be clipped off.
int cellCounts_find_soft_clipping(cellcounts_global_t * cct_context, int thread_no, char * read_bin, int read_offset, unsigned int mapped_pos, int test_len,  int search_to_tail, int search_center) {
	int base_in_window = 0;
	int added_base_index = 0, removed_base_index = 0;
	int search_start = 0;
	int matched_in_window = SOFT_CLIPPING_WINDOW_SIZE;
	int last_matched_base_index = -1, delta;
	gene_value_index_t * current_value_index = cct_context->value_index;

	if(search_to_tail) {
		if(search_center < 0)
			search_start = 0;
		else if(search_center >= test_len)
			// SHOULD NOT HAPPEN!!!
			search_start = test_len - 1;
		else	search_start = search_center - 1;

		delta = 1;
	}else{
		if(search_center < 0)
			// SHOULD NOT HAPPEN!!!
			search_start = 0;
		else if(search_center >= test_len)
			search_start = test_len - 1;
		else	search_start = search_center + 1;

		delta = -1;
	}

//if(FIXLENstrcmp("GACTGACACATGAGCTGAGAATTTATTTTTTTAAGCAAGATGAAAGGGGCAGCTCTAAGACAGACAGGTCACGGGCTCCCTGATAAGTTTCTGAGCTC", read_text)==0)
//SUBREADprintf("SEARCH_SCLIP to_tail %d ; test_len %d ; mappos %u\n",search_to_tail , test_len, mapped_pos);
	for(added_base_index = search_start; added_base_index >= 0 && added_base_index < test_len; added_base_index += delta) {
		// add the new base
		char reference_base = cellCounts_get_index_int(current_value_index , added_base_index+mapped_pos);
		int added_is_matched = reference_base == cellCounts_get_read_int(read_bin, read_offset+added_base_index);

		matched_in_window += added_is_matched;
		if(added_is_matched)
			last_matched_base_index = added_base_index;

		base_in_window ++;

		if(base_in_window > SOFT_CLIPPING_WINDOW_SIZE){
			removed_base_index = added_base_index - delta * SOFT_CLIPPING_WINDOW_SIZE;
			char removing_ref_base = cellCounts_get_index_int(current_value_index, removed_base_index + mapped_pos );
			matched_in_window -= removing_ref_base == cellCounts_get_read_int(read_bin, removed_base_index + read_offset);
		}else{
			matched_in_window --;
		}

		if(matched_in_window < SOFT_CLIPPING_WINDOW_SIZE - SOFT_CLIPPING_MAX_ERROR){
			// clip, bondary is the last matched base.
			if(search_to_tail){
				if(last_matched_base_index < 0) return test_len - search_start;
				else return test_len - last_matched_base_index - 1;
			}else{
				if(last_matched_base_index >= 0) return last_matched_base_index;
				else return search_start - 1;
			}
		}
	}

	if(last_matched_base_index < 0) return test_len;

	if(search_to_tail){
		if(last_matched_base_index < 0) return test_len - search_start;
		else return test_len - last_matched_base_index - 1;
	}else{
		if(last_matched_base_index >= 0) return last_matched_base_index;
		else return search_start - 1;
	}
}

static struct option cellCounts_long_options[]={
	{"dataset", required_argument ,0,0},
	{"index", required_argument ,0,0},
	{"inputMode", required_argument ,0,0},
	{"output", required_argument ,0,0},
	{"threads", required_argument ,0,0},

	{"annotation", required_argument ,0,0},
	{"isGTFannotation", no_argument ,0,0},
	{"enableSoftClipping", no_argument ,0,0},
	{"geneIdColumn", required_argument ,0,0},
	{"annotationType", required_argument ,0,0},
	{"annotationChroAlias", required_argument ,0,0},
	{"reportExcludedBarcodes", required_argument ,0,0},

	{"readAssignmentFile",required_argument, 0,0},
	{"cellBarcodeFile",required_argument, 0,0},
	{"sampleSheetFile",required_argument, 0,0},
	{"reportMultiMappingReads", no_argument ,0,0},
	{"junctionDetection", no_argument ,0,0},
	{"binaryTempMemory", no_argument ,0,0},
	{"VisiumHD_barcode", required_argument ,0,0},
	{"cluster_junctions", required_argument ,0,0},
	{"cluster_map", required_argument ,0,0},

	{"maxDiffToTopVotes", required_argument ,0,0},
	{"maxMismatch", required_argument ,0,0},
	{"subreadsPerRead",required_argument,0,0},
	{"minVotesPerRead",required_argument,0,0},
	{"minMappedLength",required_argument,0,0},
	{"umiCutoff",required_argument,0,0},
	{"reportedAlignmentsPerRead",required_argument,0,0},

	{0,0,0,0}
};


void cellCounts_print_config(cellcounts_global_t * cct_context);
#define _gehash_hash(k) ((unsigned int)(k))
#define _gehash_get_bucket(tab, key)  ( (tab) -> buckets + _gehash_hash( key ) % (tab) -> buckets_number )

void prefill_votes(gehash_t * the_table, temp_votes_per_read_t * pnts, int applied_subreads, unsigned int subread, int offset, int subread_no, int is_negative_strand){
	struct gehash_bucket * current_bucket;
	current_bucket = _gehash_get_bucket (the_table, subread);
	int items = current_bucket -> current_items;
	int my_no = subread_no + applied_subreads * is_negative_strand;
	pnts -> votes[my_no] = 0;
	if(!items) return;

	short *current_keys = current_bucket -> new_item_keys;
	int imin=0, imax=items - 1;
	int last_accepted_index;
	short key = subread / the_table->buckets_number;
	while(1){
		last_accepted_index=(imin+imax)/2;
		short current_key = current_keys[last_accepted_index];
		if(current_key>key) imax = last_accepted_index - 1;
		else if(current_key<key) imin = last_accepted_index + 1;
		else break;

		if(imax<imin) return;
	}

	imax -= imin;
	int start_scan_idx = last_accepted_index, tl_last_accepted_index, prefill_big_search_step = imax/4;
	int stoploc;
	
	for(; prefill_big_search_step > 1 ; prefill_big_search_step /= 3)while(1){
		tl_last_accepted_index = last_accepted_index + prefill_big_search_step;
		if(tl_last_accepted_index >= items || current_keys[tl_last_accepted_index]!=key)break;
		last_accepted_index = tl_last_accepted_index;
	}

	while(1){
		last_accepted_index ++;
		if(last_accepted_index == items|| current_keys[last_accepted_index]!= key){
			stoploc = last_accepted_index ; 
			last_accepted_index = start_scan_idx;
			break;
		}
	}

	prefill_big_search_step = imax /4;
	for(; prefill_big_search_step > 1 ; prefill_big_search_step /= 3)while(1){
		tl_last_accepted_index = last_accepted_index - prefill_big_search_step;
		if(tl_last_accepted_index < imin|| current_keys[tl_last_accepted_index]!=key)break;
		last_accepted_index = tl_last_accepted_index;
	}

	while(1){
		if(last_accepted_index == imin || current_keys[last_accepted_index -1]!= key){
			pnts -> start_location_in_index[my_no] = current_bucket->item_values + last_accepted_index;
			pnts -> votes[my_no] = stoploc - last_accepted_index;
			pnts -> offsets[my_no] = offset;
			break;
		}
		last_accepted_index --;
	}
}

int cellCounts_args_context(cellcounts_global_t * cct_context, int argc, char** argv){
	int c , option_index=0;

	optind = 0;
	opterr = 1;
	optopt = 63;

	int cmd_rebuilt_size = 2000;
	char * cmd_rebuilt = malloc(cmd_rebuilt_size);

	cmd_rebuilt[0]=0;
	for(c = 0; c<argc;c++)
	{
		int needed_buff_len = strlen(cmd_rebuilt) + 100+strlen(argv[c]);
		if(needed_buff_len > cmd_rebuilt_size)
		{
			cmd_rebuilt_size = max(cmd_rebuilt_size *2, needed_buff_len);
			cmd_rebuilt = realloc(cmd_rebuilt, cmd_rebuilt_size);
		}
		SUBreadSprintf(cmd_rebuilt+strlen(cmd_rebuilt), cmd_rebuilt_size - strlen(cmd_rebuilt), "\"%s\" ", argv[c]);
	}

	cct_context -> input_mode = GENE_INPUT_BCL;
	cct_context -> total_threads = 10;
	cct_context -> features_annotation_file_type = FILE_TYPE_RSUBREAD;
	cct_context -> reads_per_chunk = 0x77000000;
	cct_context -> max_reported_alignments_per_read = 1;
//#warning "====== SHOULD REMOVE '* 4' BELOW ========="
	cct_context -> max_candidate_voteIJ_per_read = 3 /** 4*/;

//#warning "====== SHOULD REMOVE '+ 4' BELOW ========="
	cct_context -> max_differential_from_top_vote_number = 2 /*+ 4*/;


	cct_context -> chroEvent_lock_number = 241;	// arbitrary prime number.
	cct_context -> max_mismatching_bases_in_reads = 3;
	cct_context -> max_indel_length = 5;
	cct_context -> umi_cutoff = -1;
	cct_context -> min_votes_per_mapped_read = 3;
	cct_context -> total_subreads_per_read = 10;
	cct_context -> max_distinct_top_vote_numbers = cct_context -> total_subreads_per_read;
	cct_context -> is_BAM_and_FQ_out_generated = 1;
	cct_context -> current_dataset_no = 1;
	cct_context -> min_mapped_length_for_mapped_read = 40;
	cct_context -> cmd_rebuilt = cmd_rebuilt;
	cct_context -> need_check_strand = 1;
	cct_context -> do_cell_level_junction_detection = 0;

	if(0){
		SUBREADprintf("===== Strand is reversely checked for spatial =====\n");
		SUBREADprintf("===== Strand is reversely checked for spatial =====\n");
		SUBREADprintf("===== Strand is reversely checked for spatial =====\n");
		SUBREADprintf("===== Strand is reversely checked for spatial =====\n");
		cct_context -> need_check_strand = -1;
	}

	strcpy(cct_context -> temp_file_dir, "./");

	while (1){
		c = getopt_long(argc, argv, "", cellCounts_long_options, &option_index);
		if(c<0 || c==255)break;

		if(strcmp("maxMismatch", cellCounts_long_options[option_index].name)==0){
			cct_context -> max_mismatching_bases_in_reads = max(0, atoi(optarg));
		}
		if(strcmp("enableSoftClipping", cellCounts_long_options[option_index].name)==0){
			cct_context -> enable_soft_clipping=1;
		}
		if(strcmp("minMappedLength", cellCounts_long_options[option_index].name)==0){
			cct_context -> min_mapped_length_for_mapped_read = min(MAX_SCRNA_READ_LENGTH, max(-1, atoi(optarg)));
		}
		if(strcmp("minVotesPerRead", cellCounts_long_options[option_index].name)==0){
			cct_context -> min_votes_per_mapped_read = min(64, max(1, atoi(optarg)));
		}
		if(strcmp("subreadsPerRead", cellCounts_long_options[option_index].name)==0){
			cct_context -> total_subreads_per_read = min(SCRNA_SUBREADS_HARD_LIMIT, max(7, atoi(optarg)));
		}
		if(strcmp("reportExcludedBarcodes", cellCounts_long_options[option_index].name)==0){
			cct_context -> report_excluded_barcodes = atoi(optarg);
		}
		if(strcmp("dataset", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> input_dataset_name, optarg, MAX_FILE_NAME_LENGTH * MAX_SCRNA_FASTQ_FILES * 3 -1);
		}
		if(strcmp("maxDiffToTopVotes", cellCounts_long_options[option_index].name)==0){
			cct_context -> max_differential_from_top_vote_number = min(30, max(1, atoi(optarg)));
		}
		if(strcmp("index", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> index_prefix, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("inputMode", cellCounts_long_options[option_index].name)==0){
			if(strcmp("FASTQ", optarg)==0) cct_context -> input_mode = GENE_INPUT_SCRNA_FASTQ;
			if(strcmp("BAM", optarg)==0) cct_context -> input_mode = GENE_INPUT_SCRNA_BAM;
		}
		if(strcmp("output", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> output_prefix, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("reportedAlignmentsPerRead", cellCounts_long_options[option_index].name)==0){
			cct_context -> max_reported_alignments_per_read = min(SCRNA_HIGHEST_REPORTED_ALIGNMENTS, max(1, atoi(optarg)));
		}
		if(strcmp("threads", cellCounts_long_options[option_index].name)==0){
			cct_context -> total_threads = min(64, max(1, atoi(optarg)));
		}
		if(strcmp("annotation", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> features_annotation_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("annotationChroAlias", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> features_annotation_alias_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("annotationType", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> features_annotation_feature_type, optarg, MAX_READ_NAME_LEN-1);
		}
		if(strcmp("reportMultiMappingReads", cellCounts_long_options[option_index].name)==0){
			cct_context -> report_multi_mapping_reads = 1;
		}
		if(strcmp("geneIdColumn", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> features_annotation_gene_id_column, optarg, MAX_READ_NAME_LEN-1);
		}
		if(strcmp("isGTFannotation", cellCounts_long_options[option_index].name)==0){
			cct_context -> features_annotation_file_type = FILE_TYPE_GTF;
		}
		if(strcmp("readAssignmentFile", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> read_assignment_detail_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("binaryTempMemory", cellCounts_long_options[option_index].name)==0){
			cct_context -> cell_level_junction_memory_temp = 1;
		}
		if(strcmp("junctionDetection", cellCounts_long_options[option_index].name)==0){
			cct_context -> do_cell_level_junction_detection = 1;
			//#warning "====== Yang Liao added '+20' for better junction calling -- think about it in final release??? ======"
			cct_context -> total_subreads_per_read = (10 /* +20*/);
		}
		if(strcmp("cellBarcodeFile", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> cell_barcode_list_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}

		if(strcmp("cluster_junctions", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> cluster_junctions_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("cluster_map", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> cluster_map_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("VisiumHD_barcode", cellCounts_long_options[option_index].name)==0){
			strcpy(cct_context -> visium_hd_CellRanger_bam, optarg);
			cct_context -> visium_hd_barcodes = 1;
			cct_context -> UMI_length = 9; // observed from example data
		}

		if(strcmp("sampleSheetFile", cellCounts_long_options[option_index].name)==0){
			strncpy(cct_context -> bcl_sample_sheet_file, optarg, MAX_FILE_NAME_LENGTH -1);
		}
		if(strcmp("umiCutoff", cellCounts_long_options[option_index].name)==0){
			cct_context -> umi_cutoff = atof(optarg);
//			SUBREADprintf("UMI_CUT=%.2f\n", cct_context -> umi_cutoff);
		}
	}


	char * DBPZ_cellCounts_CHUNK_READS = getenv("DBPZ_cellCounts_CHUNK_READS");
	if(DBPZ_cellCounts_CHUNK_READS) {
		cct_context -> reads_per_chunk = atoi(DBPZ_cellCounts_CHUNK_READS);
		SUBREADprintf("SET total reads for short running: %d\n", cct_context -> reads_per_chunk);
	}

	char * DBPZ_var_candidates = getenv("DBPZ_CCNT_CANDIDATES");
	char * DBPZ_var_topdiff = getenv("DBPZ_CCNT_TOPDIFF");
	if( DBPZ_var_candidates && DBPZ_var_topdiff ){
		cct_context -> max_candidate_voteIJ_per_read = atoi(DBPZ_var_candidates);
		cct_context -> max_differential_from_top_vote_number = atoi(DBPZ_var_topdiff);
		SUBREADprintf("OverSET CAND and TOPDIFF %d and %d\n",  cct_context -> max_candidate_voteIJ_per_read,  cct_context -> max_differential_from_top_vote_number);
	}
	cct_context -> max_distinct_top_vote_numbers = min(cct_context -> max_distinct_top_vote_numbers, 1+cct_context -> max_differential_from_top_vote_number);
	cct_context -> max_candidate_voteIJ_per_read = max(cct_context -> max_candidate_voteIJ_per_read, cct_context -> max_reported_alignments_per_read);

	return 0;
}

void cellCounts_print_config(cellcounts_global_t * cct_context){

	SUBREADputs("        ==========     _____ _    _ ____  _____  ______          _____  ");
	SUBREADputs("        =====         / ____| |  | |  _ \\|  __ \\|  ____|   /\\   |  __ \\ ");
	SUBREADputs("          =====      | (___ | |  | | |_) | |__) | |__     /  \\  | |  | |");
	SUBREADputs("            ====      \\___ \\| |  | |  _ <|  _  /|  __|   / /\\ \\ | |  | |");
	SUBREADputs("              ====    ____) | |__| | |_) | | \\ \\| |____ / ____ \\| |__| |");
	SUBREADputs("        ==========   |_____/ \\____/|____/|_|  \\_\\______/_/    \\_\\_____/");
	SUBREADprintf("       %s\n",SUBREAD_VERSION);
	SUBREADputs("");


	print_in_box(80,1,PRINT_BOX_CENTER,"cellCounts settings");
	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"         Index : %s", cct_context -> index_prefix);
	print_in_box(80,0,0,"    Input mode : %s", cct_context -> input_mode == GENE_INPUT_SCRNA_FASTQ?"FASTQ files":(
	    cct_context -> input_mode == GENE_INPUT_SCRNA_BAM?"BAM files":"Raw BCL files"));
	print_in_box(80,0,0,"");
	print_in_box(80,2,PRINT_BOX_CENTER,"");
	SUBREADputs(  "");
}


int determine_total_index_blocks(cellcounts_global_t * cct_context){
	char tmp_fname[MAX_FILE_NAME_LENGTH+ 30];
	cct_context-> total_index_blocks = 0;
	while(1){
		SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+ 30, "%s.%02d.b.tab", cct_context->index_prefix, cct_context->total_index_blocks);
		if(!does_file_exist(tmp_fname))break;
		cct_context->total_index_blocks ++;
	}
	if(cct_context->total_index_blocks> 1){
		SUBREADprintf("ERROR: cellCounts can only run with one-block index. Please build the index with indexSplit=FALSE.\n");
		return 1;
	}
	return 0;
}

void sheet_convert_ss_to_arr( void * key, void * hashed_obj, HashTable * tab ){
	ArrayList * hashed_arr = hashed_obj ;
	cellcounts_global_t * cct_context = tab->appendix1;
	ArrayListPush(cct_context -> sample_id_to_name, key);
	hashed_arr -> appendix1 = NULL+ cct_context -> sample_id_to_name -> numOfElements; // One-based
					
	srInt_64 xx1;		   
	for(xx1 =0; xx1< hashed_arr -> numOfElements; xx1++){
		char ** push_arr = malloc(sizeof(char*)*4); 
		char ** sbc_lane_sample = ArrayListGet(hashed_arr, xx1);
		srInt_64 lane_sample_int = sbc_lane_sample[0]-(char*)NULL;
			
		ArrayListPush(cct_context -> sample_barcode_list, push_arr);
		push_arr[0] = NULL + lane_sample_int;
		push_arr[1] = NULL + cct_context -> sample_id_to_name -> numOfElements;
		push_arr[2] = sbc_lane_sample[1]; // Sample Barcode
		push_arr[3] = NULL + (sbc_lane_sample[1]!=NULL && strlen(sbc_lane_sample[1])>12);
		int line_no_in_sheet = sbc_lane_sample[2] - (char*)NULL;
		HashTablePut(cct_context -> lineno1B_to_sampleno1B_tab , NULL+line_no_in_sheet, NULL + cct_context -> sample_id_to_name -> numOfElements);
	}
}       


void cellCounts_close_sample_SamBam_writers(void *v){
	void ** vv = v;
	simple_bam_writer * wtr = vv[0];
	simple_bam_close(wtr);

	if(vv[1]){
		parallel_gzip_writer_t* gzfp = vv[1];
		parallel_gzip_writer_close(gzfp);

		gzfp = vv[2];
		parallel_gzip_writer_close(gzfp);

		gzfp = vv[3];
		if(gzfp)parallel_gzip_writer_close(gzfp);

		gzfp = vv[4];
		parallel_gzip_writer_close(gzfp);
	}

	cellCounts_lock_t * gz_lock = vv[5];
	cellCounts_destroy_lock(gz_lock);
	free(gz_lock);
	
	free(vv);
}

void cellCounts_add_simple_writer_header(cellcounts_global_t * cct_context, simple_bam_writer * wtr){
	int x1, outcapa, outsize=0;
	unsigned int last_end = 0;
	int rebuilt_len = strlen(cct_context->cmd_rebuilt) + 200;
	outcapa = rebuilt_len + MAX_CHROMOSOME_NAME_LEN *2 + 200;
	char * outbuff = malloc(outcapa);

	outsize += SUBreadSprintf(outbuff, outcapa, "@HD\tVN:1.0\tSO:%s\n", "coordinate");
	for(x1=0;x1<cct_context->chromosome_table.total_offsets; x1++){
		if(outsize + rebuilt_len + 200 + MAX_CHROMOSOME_NAME_LEN + 40 >= outcapa){
			outcapa *=2;
			outbuff = realloc(outbuff, outcapa);
		}
		char * chro_name = cct_context->chromosome_table.read_names+x1*MAX_CHROMOSOME_NAME_LEN;
		unsigned int this_end = cct_context->chromosome_table.read_offsets[x1];
		int this_size = this_end - last_end;
		last_end = this_end;
		outsize += SUBreadSprintf(outbuff+outsize, outcapa-outsize,"@SQ\tSN:%s\tLN:%d\n",chro_name,this_size);
	}

	// if the command line is longer than 16K the remaining part is omitted.
	outsize += snprintf(outbuff + outsize, outcapa - outsize , "@PG\tID:cellCounts\tPN:cellCounts\tVN:%s\tCL:%s", SUBREAD_VERSION, cct_context->cmd_rebuilt);
	outbuff[outsize++]='\n';
	outbuff[outsize++]='\0';
	simple_bam_write(&outsize,4,wtr,0);
	simple_bam_write(outbuff, outsize, wtr, 1);
	free(outbuff);

	simple_bam_write(&cct_context->chromosome_table.total_offsets,4,wtr,0);
	for(x1=0;x1<cct_context->chromosome_table.total_offsets; x1++){
		char * chro_name = cct_context->chromosome_table.read_names+x1*MAX_CHROMOSOME_NAME_LEN;
		unsigned int this_end = cct_context->chromosome_table.read_offsets[x1];
		int this_size = this_end - last_end;
		last_end = this_end;
		int sqname_len = strlen(chro_name)+1;
		simple_bam_write(&sqname_len, 4, wtr, 0);
		simple_bam_write(chro_name, sqname_len, wtr, 0);
		simple_bam_write(&this_size, 4, wtr, 0);
	}
	simple_bam_write("",0,wtr,1);
	wtr -> total_chromosomes = cct_context->chromosome_table.total_offsets;
}

void cellCounts_sample_SamBam_writers_new_files(void *k, void *v, HashTable * tab){
	HashTable * fp_tab = tab -> appendix1;
	cellcounts_global_t * cct_context = tab -> appendix2;
	ArrayList * scRNA_sample_id_to_name = tab -> appendix3;

	char * samplename = k;
	char fname [MAX_FILE_NAME_LENGTH+20], fnamet[MAX_FILE_NAME_LENGTH+20];
	SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "%s.bam", samplename);
	SUBreadSprintf(fnamet,MAX_FILE_NAME_LENGTH+20, "del4-cC-tmp0-%s.del", samplename);
	simple_bam_writer * wtr = simple_bam_create(fname);
	cellCounts_add_simple_writer_header(cct_context, wtr);
	parallel_gzip_writer_t * gzipR1fq=NULL, * gzipI1fq=NULL, * gzipR2fq=NULL, * gzipI2fq=NULL;

	if(cct_context -> input_mode == GENE_INPUT_BCL || cct_context -> input_mode == GENE_INPUT_SCRNA_BAM){
		gzipR1fq = calloc(sizeof(parallel_gzip_writer_t),1);
		gzipI1fq = calloc(sizeof(parallel_gzip_writer_t),1);
		if(cct_context -> is_dual_index)gzipI2fq = calloc(sizeof(parallel_gzip_writer_t),1);
		gzipR2fq = calloc(sizeof(parallel_gzip_writer_t),1);
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "%s_R1.fastq.gz", samplename);
		parallel_gzip_writer_init(gzipR1fq, fname, cct_context -> total_threads);
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "%s_I1.fastq.gz", samplename);
		parallel_gzip_writer_init(gzipI1fq, fname, cct_context -> total_threads);
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "%s_I2.fastq.gz", samplename);
		if(gzipI2fq)parallel_gzip_writer_init(gzipI2fq, fname, cct_context -> total_threads);
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "%s_R2.fastq.gz", samplename);
		parallel_gzip_writer_init(gzipR2fq, fname, cct_context -> total_threads);
	}

	cellCounts_lock_t * gzfp_lock = malloc(sizeof(cellCounts_lock_t));
	cellCounts_init_lock(gzfp_lock, 0);
	int x1;
	for(x1=0; x1<scRNA_sample_id_to_name -> numOfElements; x1++){
		char * sample_name = ArrayListGet( scRNA_sample_id_to_name, x1 );
		if(strcmp(sample_name, samplename)==0){
			void ** wtrptr = malloc(sizeof(void*)*6);
			wtrptr[0]=wtr;
			wtrptr[1]=gzipR1fq;
			wtrptr[2]=gzipI1fq;
			wtrptr[3]=gzipI2fq;
			wtrptr[4]=gzipR2fq;
			wtrptr[5]=gzfp_lock;
			HashTablePut(fp_tab, NULL+x1+1 , wtrptr);
			break;
		}
	}
}

void extract_sam_tags(const char *bambuff, char *BAM1R, char *BAM1Y, char *CB, char *UR, char *UY) {
	// Initialize pre-allocated buffers to empty strings in case tags aren't found
	if (BAM1R) BAM1R[0] = '\0';
	if (BAM1Y) BAM1Y[0] = '\0';
	if (CB) CB[0] = '\0';
	if (UR) UR[0] = '\0';
	if (UY) UY[0] = '\0';

	const char *p = bambuff;
	int tab_count = 0;

	// 1. Skip the 11 mandatory SAM fields
	while (*p && tab_count < 11) {
		if (*p == '\t') {
			tab_count++;
		}
		p++;
	}

	// If the line is malformed or doesn't have optional fields, exit early
	if (tab_count < 11) return;

	// 2. Parse the optional fields
	while (*p) {
		// Check if the current pointer matches any of our target tags followed by a colon
		if (strncmp(p, "1R:", 3) == 0 || strncmp(p, "1Y:", 3) == 0 || strncmp(p, "CB:", 3) == 0 || strncmp(p, "UR:", 3) == 0 || strncmp(p, "UY:", 3) == 0) {
			const char *tag_start = p;
			
			// SAM optional fields format is TAG:TYPE:VALUE (e.g., BAM1R:Z:ATGCA)
			// Ensure the structure has the second colon after the 1-character TYPE
			if (*(p + 2) == ':' && *(p + 3) != '\0' && *(p + 4) == ':') {
				const char *val_start = p + 5; // Value starts right after the second colon
				const char *val_end = val_start;
				
				// Find the end of the current field (tab or end of line)
				while (*val_end && *val_end != '\t' && *val_end != '\n' && *val_end != '\r') {
					val_end++;
				}
				
				size_t len = val_end - val_start;
				
				// Copy the value into the corresponding pre-allocated variable
				if (strncmp(tag_start, "1R", 2) == 0 && BAM1R) {
					strncpy(BAM1R, val_start, len);
					BAM1R[len] = '\0';
				} else if (strncmp(tag_start, "1Y", 2) == 0 && BAM1Y) {
					strncpy(BAM1Y, val_start, len);
					BAM1Y[len] = '\0';
				} else if (strncmp(tag_start, "UY", 2) == 0 && UY) {
					strncpy(UY, val_start, len);
					UY[len] = '\0';
				} else if (strncmp(tag_start, "UR", 2) == 0 && UR) {
					strncpy(UR, val_start, len);
					UR[len] = '\0';
				} else if (strncmp(tag_start, "CB", 2) == 0 && CB) {
					strncpy(CB, val_start, len);
					CB[len] = '\0';
				}
			}
		}

		// Advance to the next optional field (skip to the next tab)
		while (*p && *p != '\t' && *p != '\n' && *p != '\r') {
			p++;
		}
		
		if (*p == '\t') {
			p++; // Skip the tab character to point to the start of the next tag
		} else {
			break; // Reached end of the line (\n, \r, or \0)
		}
	}
}

static void copy_strings_fast(const char *readfl, char *BAM1R, char *BAM1Y, char *CB) {
	const char *start = readfl;
	const char *p = readfl;
	size_t len;

	while (*p && *p != '\t') {
		p++;
	}
	len = p - start;
	memcpy(BAM1R, start, len);
	BAM1R[len] = '\0'; // Manually null-terminate destination

	if (*p == '\t') p++; 
	
	start = p;
	while (*p && *p != '\t') {
		p++;
	}
	len = p - start;
	memcpy(BAM1Y, start, len);
	BAM1Y[len] = '\0';

	if (*p == '\t') p++;

	start = p;
	while (*p && *p != '\t' && *p != '\n' && *p != '\r') {
		p++;
	}
	len = p - start;
	memcpy(CB, start, len);
	CB[len] = '\0';
}

#define FILE_IS_BAM   1
#define FILE_IS_TEXT  0
#define FILE_ERROR   -1

static int check_file_type(const char *filename) {
	FILE *file = fopen(filename, "rb");
	if (file == NULL) {
		perror("Error opening file");
		return FILE_ERROR;
	}
	unsigned char header[2];
	size_t bytes_read = fread(header, 1, 2, file);
	fclose(file);

	if (bytes_read < 2) return FILE_IS_TEXT; 
	if (header[0] == 0x1f && header[1] == 0x8b) return FILE_IS_BAM; // gzipped
	return FILE_IS_TEXT;
}

int cellCounts_load_scRNA_tables(cellcounts_global_t * cct_context){
	int rv = 0;
	if(cct_context -> visium_hd_CellRanger_bam[0]){
		char* early_terminate = getenv("DBPZ_cellCounts_BARREF_TERMINAL");
		cct_context -> VisiumHD_barcode_to_best_mapping = StringTableCreate(32*1024*1024+39);
		HashTableSetDeallocationFunctions(cct_context -> VisiumHD_barcode_to_best_mapping,free,free); // key: BAM1R+CY; val: CB; all duplicated.
		int is_BAM_file = check_file_type(cct_context -> visium_hd_CellRanger_bam)==FILE_IS_BAM;
		void * fparby;
//fprintf(stderr,"BAMREF %d\n", is_BAM_file);

		if(is_BAM_file) fparby = SamBam_fopen(cct_context -> visium_hd_CellRanger_bam, SAMBAM_FILE_BAM);
		else fparby = fopen(cct_context -> visium_hd_CellRanger_bam, "r");
		int offset_of_1R = is_BAM_file? cct_context -> UMI_length :0;
		while(1){
			char BAM1R[100],BAM1Y[100],CB[100];
			char bambuff[5001];
			char * readfl;
			BAM1R[0]=BAM1Y[0]=CB[0]=0;

			if(is_BAM_file){
				readfl = SamBam_fgets((SamBam_FILE*)fparby,bambuff, 5000,0);
				if(!readfl)break;
				if(bambuff[0]=='@')continue;
//fprintf(stderr,"REFLINE %s\n", bambuff);
				extract_sam_tags(bambuff, BAM1R, BAM1Y, CB, NULL, NULL);
			} else {
				readfl = fgets(bambuff, 5000, (FILE*) fparby);
				if(!readfl)break;
				copy_strings_fast(readfl, BAM1R, BAM1Y, CB);
			}

			if(CB[0] && BAM1R[0] && BAM1Y[0]){
				char CKey[200];
				char CVal[200];
				int x1, bc1=-1, bc2=-1;
				for(x1=0; BAM1Y[x1]; x1++)if(BAM1Y[x1]>='/') BAM1Y[x1] ++; // when matching the quality string, the qualty char >= '/' is +1
				sprintf(CKey,"%s/%s",BAM1R + offset_of_1R ,BAM1Y + offset_of_1R); // UMI is before spot barcodes in Visium HD R1. And we don't need to index the UMI for having the barcodes.
				sscanf(CB, "s_%*[^_]_%d_%d", &bc1, &bc2);
				sprintf(CVal,"%05d_%05d",bc1,bc2);

				char * oldCB = HashTableGet(cct_context -> VisiumHD_barcode_to_best_mapping, CKey);
				if(oldCB){
					if(0)if(strcmp(CB, oldCB)!=0)SUBREADprintf("ERROR: the same BAM1R and BAM1Y are mapped to different CB: %s and %s have %s != %s\n", BAM1R, BAM1Y, CB, oldCB);
				}else HashTablePut(cct_context -> VisiumHD_barcode_to_best_mapping, strdup(CKey), strdup(CVal));
				if(cct_context -> VisiumHD_barcode_to_best_mapping->numOfElements % 3000000==0)//fprintf(stderr,"INSERT_FROM_BAN %s %s  OLD %p\n", CKey, CVal, oldCB);
					SUBREADprintf("Loaded the %lld-th barcode from Space Ranger reference.\n", cct_context -> VisiumHD_barcode_to_best_mapping->numOfElements);
			}
			if(early_terminate)  if(cct_context -> VisiumHD_barcode_to_best_mapping->numOfElements > 3654321)break;
		}
		if(is_BAM_file) SamBam_fclose((SamBam_FILE*)fparby);
		else fclose((FILE*)fparby);
	}else{
		cct_context -> cell_barcodes_array = input_BLC_parse_CellBarcodes( cct_context-> cell_barcode_list_file );

		if(NULL == cct_context-> cell_barcodes_array){
			SUBREADprintf("ERROR: cannot find valid cell barcodes from the cell barcode list. Please check the content and the accessibility of the file.\n");
			rv = 1;
		}
	}
	if(!rv){
		if(cct_context-> cell_barcodes_array)rv = cellCounts_make_barcode_HT_table( cct_context );
		if(!rv){
			cct_context-> sample_sheet_table = input_BLC_parse_SampleSheet( cct_context -> bcl_sample_sheet_file);
			if(NULL == cct_context-> sample_sheet_table) rv = 1;
			if(rv==0 && cct_context-> sample_sheet_table -> numOfElements > MAX_SCRNA_SAMPLE_NUMBER){
				SUBREADprintf("ERROR: too many samples in the sample sheet.\n");
				rv = 1;
			}
			if(!rv){
				cct_context -> sample_id_to_name = ArrayListCreate(64);
				cct_context -> lineno1B_to_sampleno1B_tab = HashTableCreate(40);

				cct_context -> sample_sheet_table -> appendix1 = cct_context;
				cct_context -> sample_barcode_list = ArrayListCreate(64);

				ArrayListSetDeallocationFunction(cct_context -> sample_barcode_list, free);
				HashTableIteration(cct_context-> sample_sheet_table, sheet_convert_ss_to_arr);

				if(cct_context -> is_BAM_and_FQ_out_generated){
					cct_context -> sample_BAM_writers = HashTableCreate(cct_context -> sample_sheet_table -> numOfElements);
					HashTableSetDeallocationFunctions(cct_context -> sample_BAM_writers, NULL, cellCounts_close_sample_SamBam_writers);
					cct_context -> sample_sheet_table ->appendix1 = cct_context -> sample_BAM_writers;
					cct_context -> sample_sheet_table ->appendix2 = cct_context;
					cct_context -> sample_sheet_table ->appendix3 = cct_context -> sample_id_to_name;
					HashTableIteration( cct_context -> sample_sheet_table, cellCounts_sample_SamBam_writers_new_files);
				}

			}
		}
	}
	return rv;
}

int cellCounts_load_base_value_indexes(cellcounts_global_t * cct_context){
	int rv=0;
	char tmp_fname[MAX_FILE_NAME_LENGTH+ 30];
	SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+30, "%s.%02d.b.array", cct_context ->index_prefix, 0);
	cct_context -> value_index = calloc(sizeof(gene_value_index_t),1);
	rv = rv || gvindex_load_largebuffer(cct_context -> value_index, tmp_fname, 1); // pre-decompress base values.
	return rv;
}

typedef struct {
	srInt_64 feature_name_pos;
	unsigned int start;
	unsigned int end;
	unsigned int sorted_order;

	unsigned short chro_name_pos_delta;
	char is_negative_strand;
	char * extra_columns;
} fc_feature_info_t;

typedef struct {
	unsigned int chro_number;
	unsigned int chro_features;
	unsigned int chro_feature_table_start;
	unsigned int chro_block_table_end;
	unsigned int chro_possible_length;

	unsigned short chro_reverse_table_current_size;
	unsigned int * reverse_table_start_index;
	int reverse_table_start_index_size;
	//unsigned int * reverse_table_end_index;
} fc_chromosome_index_info;

srInt_64 cellCounts_unistr_cpy(cellcounts_global_t * cct_context, char * str, int strl)
{
	srInt_64 ret;
	if(cct_context->unistr_buffer_used + strl >= cct_context->unistr_buffer_size-1)
	{
		if( cct_context->unistr_buffer_size < (1000llu*1000u*1000u*32)) // 32GB
		{
			cct_context -> unistr_buffer_size = cct_context->unistr_buffer_size /2 *3;
			cct_context -> unistr_buffer_space = realloc(cct_context -> unistr_buffer_space, cct_context->unistr_buffer_size);
		}
		else
		{
			SUBREADprintf("Error: exceed memory limit (32GB) for storing feature names.\n");
			return 0xffffffffu;
		}
	}

	strcpy(cct_context -> unistr_buffer_space + cct_context->unistr_buffer_used, str);
	ret = cct_context->unistr_buffer_used;

	cct_context->unistr_buffer_used += strl +1;

	return ret;
}

void cellCounts_write_one_chroEvent(void *k, void *v, HashTable * tab){
	chroEvent_t *ce = v;
	if(ce->event_type == chroEvent_t_TYPE_JUNCTION && ce -> step2_supported_reads>0){
		void ** params = tab -> appendix1;
		FILE * wfp = params[0];
		cellcounts_global_t* cct_context = params[1];

		char * chro=NULL;
		int pos1 = 0;
		locate_gene_position(ce -> left_edge +1, &cct_context -> chromosome_table, &chro, &pos1);
		int pos2 = pos1 + ( (void*) ce->length - NULL ) +1;
		
		fprintf(wfp,"%s\t%d\t%d\t%s\t%d\n", chro, pos1, pos2, ce -> from_truth?"ANNOTATED":"NEW",*ce -> step2_supported_reads);
	} 
}

void cellCounts_write_final_junctions(cellcounts_global_t * cct_context,  char * output_file_name);
void cellCounts_write_final_other_events(cellcounts_global_t * cct_context,  char * output_file_name, int type_to_write);
int cellCounts_write_junction_sumtable(cellcounts_global_t* cct_context){
	void * params [5];
#warning "===== CURRENT JUNCTION DETECTION DOESN'T NEED SAMPLE-LEVEL SUPPORT OUT ====="
	if(0)cellCounts_write_final_junctions(cct_context, cct_context -> output_prefix);
	// DO NOT DETELE: They are for indel and exon detection.
	//cellCounts_write_final_other_events(cct_context, "del4-indels.tab", chroEvent_t_TYPE_INDEL);
	//cellCounts_write_final_other_events(cct_context, "del4-exons.tab", chroEvent_t_TYPE_EXON);
	return 0;
}

// "tlen" is the number before "N" in cigar.
chroEvent_t * cellCounts_set_chroEvent_details(cellcounts_global_t* cct_context, int sample_i, int event_type, unsigned int linear_loc, unsigned int linear_loc2, int tlen, int from_truth){
	srUInt_64 envkey = (linear_loc*1LLU<<32)|(linear_loc2);

	chroEvent_t * old_or_new_env = HashTableGet( cct_context -> chroEvent_detail_table[sample_i], NULL+envkey);
	if(old_or_new_env){
		if(event_type == chroEvent_t_TYPE_INDEL && tlen <0){
			int known =0, x1;
			for(x1=0;x1 < old_or_new_env -> n_events;x1++){
				int oldlen;
				if( old_or_new_env -> n_events == 1 ) oldlen = (void *)( old_or_new_env -> length ) - NULL;
				else oldlen = old_or_new_env -> length[x1];
				if(oldlen == tlen){ 
					known=1;
					break;
				}
			}
			if(!known){
				if( old_or_new_env -> n_events == 1 ){
					int old_length = (void *)( old_or_new_env -> length ) - NULL;
					int * newlen = malloc(sizeof(int)*2);
					newlen[0]=old_length;
					newlen[1]=tlen;
					old_or_new_env -> length = newlen;
					old_or_new_env -> step2_supported_reads = calloc(sizeof(int),2);
					old_or_new_env -> step2_non_supported_reads = calloc(sizeof(int),2);

					old_or_new_env -> n_events =2;
				}else{
					old_or_new_env -> length = realloc(old_or_new_env -> length, sizeof(int)*(1+old_or_new_env -> n_events));
					old_or_new_env -> step2_supported_reads = realloc(old_or_new_env -> step2_supported_reads, sizeof(int)*(1+old_or_new_env -> n_events));
					old_or_new_env -> step2_non_supported_reads = realloc(old_or_new_env -> step2_non_supported_reads, sizeof(int)*(1+old_or_new_env -> n_events));

					old_or_new_env -> step2_supported_reads [ old_or_new_env -> n_events ] = old_or_new_env -> step2_non_supported_reads [ old_or_new_env -> n_events ] = 0;
					old_or_new_env -> length[old_or_new_env -> n_events++] = tlen;
				}
			}
		}
	}else{
		old_or_new_env = calloc(sizeof(chroEvent_t),1);
		old_or_new_env -> event_type = event_type;
		old_or_new_env -> left_edge = linear_loc;
		old_or_new_env -> n_events = 1;
		old_or_new_env -> from_truth = from_truth;
		old_or_new_env -> length = (void*)NULL+(tlen); // "tlen" is the number before "N" in cigar.
		if(event_type == chroEvent_t_TYPE_INDEL && tlen <0) old_or_new_env -> inserted_bases = calloc(1,1);

		HashTablePut(cct_context -> chroEvent_detail_table[sample_i], NULL+envkey, old_or_new_env);
	}
	return old_or_new_env;
}


int cellCounts_junction_in_table(cellcounts_global_t* cct_context, int sample_i, char * chro, int l, int r);

// l and r are 0-based.
int cellCounts_add_or_update_chroEvent_in_table(cellcounts_global_t* cct_context, int sample_i, int env_type, char * chro, int l, int r, int inslen_negative, int from_truth){
	cellCounts_lock_occupy(&cct_context -> chroEvent_entry_table_lock);
	int known = cellCounts_junction_in_table(cct_context, sample_i, chro, l, r), retv=0;

	char chro_strn_ky[MAX_CHROMOSOME_NAME_LEN+10];
	char negchar = '*';
	snprintf(chro_strn_ky, MAX_CHROMOSOME_NAME_LEN+10, "%s\t%c", chro, negchar);

	IVT_IntervalTreeNode **LR_roots = HashTableGet(cct_context -> chroEvent_entry_table[sample_i], chro_strn_ky);
	int is_new_key= !LR_roots;
	if(is_new_key) LR_roots = calloc(sizeof(void*),3);

	if(known) retv=1; else{
		LR_roots[0]=IVT_insert(LR_roots[0], l, l, NULL+r);
		LR_roots[1]=IVT_insert(LR_roots[1], r, r, NULL+l);
		if( env_type == chroEvent_t_TYPE_EXON )  LR_roots[2]=IVT_insert(LR_roots[2], l, r, NULL);
		if(is_new_key) HashTablePut(cct_context -> chroEvent_entry_table[sample_i], strdup( chro_strn_ky ), LR_roots);
	}

	int tlen=-9999990;
	if(0==known || (chroEvent_t_TYPE_INDEL == env_type && inslen_negative<0)){
		unsigned int linear_loc = linear_gene_position(&cct_context->chromosome_table , chro, l);
		unsigned int linear_loc2 = linear_loc - l + r;
		tlen = r-l-1;
		if( chroEvent_t_TYPE_INDEL == env_type && inslen_negative<0 ) tlen = inslen_negative;
		cellCounts_set_chroEvent_details(cct_context, sample_i, env_type, linear_loc, linear_loc2, tlen, from_truth);
	}

	cellCounts_lock_release(&cct_context -> chroEvent_entry_table_lock);
	return retv;
}


int cellCounts_junction_in_table(cellcounts_global_t* cct_context, int sample_i, char * chro, int l, int r){
	char chro_strn_ky[MAX_CHROMOSOME_NAME_LEN+10];
	int coloc_num = 0;
	char negchar = '*';
	snprintf(chro_strn_ky, MAX_CHROMOSOME_NAME_LEN+10, "%s\t%c", chro, negchar);
	IVT_IntervalTreeNode **LR_roots = HashTableGet(cct_context -> chroEvent_entry_table[sample_i], chro_strn_ky);
	if(!LR_roots)return 0;
	IVT_Interval* search_out [JUNCTION_MAX_COLOCATION];
	IVT_query_range( LR_roots[0], l,l, (IVT_Interval**)search_out, JUNCTION_MAX_COLOCATION, &coloc_num );

	if(coloc_num>=JUNCTION_MAX_COLOCATION-1)SUBREADprintf("WARNING: your annotation input contains very many exons ending at the same location. We only use %d of them.\n", JUNCTION_MAX_COLOCATION);
	int known = 0, x2;
	for(x2=0; x2<coloc_num; x2++) if( search_out[x2]->attr == NULL+r ){known=1;break; }
	return known;
}


void cellCounts_ins_IVT_for_a_junc(cellcounts_global_t* cct_context, char * chname_strand, int start, int end, void * attr, int is_negative){
	if(attr==NULL) SUBREADprintf("ERROR: NULL ATTR val\n");
	int x2;
	for(x2 = 0; x2 < 2;x2++){
		HashTable * mystrand_tab;
		if(x2==0) mystrand_tab = cct_context -> junction_ExonEdgeTree_table[0];
		else mystrand_tab = cct_context -> junction_ExonEdgeTree_table[1+is_negative];
		IVT_IntervalTreeNode * edge_tree_root = HashTableGet(mystrand_tab , chname_strand);
		edge_tree_root = IVT_insert(edge_tree_root, start, start, attr);
		edge_tree_root = IVT_insert(edge_tree_root, end,   end,   attr);
		HashTablePutReplaceEx(mystrand_tab, strdup(chname_strand), edge_tree_root,1,1,0);
	}
}


void cellCounts_copy_txn_to_juncs(void * ky, void * va, HashTable * me){
	int x1, sample_i;
	char * txn_chr_neg = ky;
	cellcounts_global_t* cct_context = me->appendix1;
	ArrayList * my_exons = va;
	ArrayListSort(my_exons, ArrayListLLUComparison);
	int last_end = -1, last_start = -1;
	char * chname_strand = strstr(txn_chr_neg,"\t")+1;
	char * chro_end = strstr(chname_strand,"\t");
	int is_negative = -1;
	// here: the strand info is only used for sorting transcript exons. After so, all events don't have strands in the chroEvent tables.
	if(*(chro_end+1)=='-') is_negative = 1;
	if(*(chro_end+1)=='+') is_negative = 0;

	*chro_end=0;
	*(chname_strand-1)=0;
	char * gene_name = HashTableGet(cct_context -> transcript_to_gene_name_table, txn_chr_neg);

//#warning "====== Not using input GTF for initing the event table ====="

	for(x1=0; x1 < my_exons->numOfElements; x1++){
		srInt_64 ve = ArrayListGet(my_exons, x1)-NULL;
		int start = ve >> 32;
		int end = ve & 0x7fffffff; 

		for(sample_i = 1; sample_i <=cct_context-> sample_sheet_table -> numOfElements ; sample_i ++){
			cellCounts_add_or_update_chroEvent_in_table(cct_context, sample_i, chroEvent_t_TYPE_EXON, chname_strand , start, end, 0, 1);

			if(last_end>0){
				if(start < last_end)SUBREADprintf("WARNING: Transcript '%s' contains overlapping exons (%d > %d). The junction between the overlapping exons is ignored.\n", ky, last_end, start);
				cellCounts_add_or_update_chroEvent_in_table(cct_context, sample_i, chroEvent_t_TYPE_JUNCTION, chname_strand , last_end, start, 0, 1);
			}
		}
		last_end = end;
		last_start = start;
		
		cellCounts_ins_IVT_for_a_junc(cct_context, chname_strand, start, end, gene_name, is_negative);
	}

	*chro_end='\t';
	*(chname_strand-1)='\t';
}

void cellCounts_sort_junc_feature_make_gaps(void *k, void *v, HashTable * tab){
	cct_junction_genebody_t * jg = v;
	int xk1, seek_strand, txx;
	ArrayList * transcripts = jg -> transcript_list; 
	for(seek_strand = 0; seek_strand < 2; seek_strand++){
		unsigned int min_start = 0xffffffffu;
		unsigned int max_stop = 0x0u;
		for(txx = 0; txx < transcripts -> numOfElements ; txx++){
			cct_junction_transcript_t * txl = ArrayListGet(transcripts , txx);
			ArrayList * exl = txl -> exons_in_transcript;
			for(xk1=0; xk1<exl->numOfElements; xk1++){
				cct_junction_exon_in_transcript_t * one_ex = ArrayListGet(exl,xk1);
				if(one_ex -> is_negative != seek_strand) continue;
				min_start = min(min_start, one_ex->chro_start);
				max_stop = max(max_stop, one_ex->chro_stop);
			}
		}
		if(max_stop == 0x0u) continue;
		cellcounts_global_t *cct_context = tab -> appendix1;
		IVT_IntervalTreeNode * IVT_rootnode = HashTableGet(cct_context -> junction_GenebodyTree_table, jg -> chromosome_name);
		IVT_rootnode = IVT_insert(IVT_rootnode, min_start, max_stop, jg); 
		HashTablePutReplaceEx(cct_context -> junction_GenebodyTree_table, jg -> chromosome_name, IVT_rootnode,0,0,0); // use old key ptr in table; don't free key and value.
	}
}

int cellCounts_sort_junc_feature_sort_exons_cmp(void * L_elem, void * R_elem, ArrayList * me){
	cct_junction_exon_in_transcript_t *l = L_elem, *r = R_elem;
	return l->chro_start - r->chro_start;
}

void cellCounts_sort_junc_feature_sort_exons(void *k, void *v, HashTable * tab){
	cct_junction_transcript_t * txnobj = v;
	ArrayListSort(txnobj -> exons_in_transcript, cellCounts_sort_junc_feature_sort_exons_cmp);
}

int cellCounts_extract_and_sort_juncs(cellcounts_global_t* cct_context){
	cct_context -> transcript_exon_table -> appendix1 = cct_context;
	HashTableIteration(cct_context -> transcript_exon_table, cellCounts_copy_txn_to_juncs);
//#warning "========= THE 2 TABLES ISN'T RELEASED FOR NOT CRASHING BUT IT SHOULD BE!!! ======"
	HashTableDestroy(cct_context -> transcript_exon_table);

	cct_context -> junction_genebody_table -> appendix1 = cct_context;
	HashTableIteration(cct_context -> junction_genebody_table, cellCounts_sort_junc_feature_make_gaps);
	HashTableIteration(cct_context -> junction_transcript_table, cellCounts_sort_junc_feature_sort_exons);
	return 0;
}

void cellCounts_register_junc_feature(cellcounts_global_t * cct_context, char * feature_name, char * transcript_id, char * chro, unsigned int start, unsigned int stop, int is_negative);

int features_load_one_line(char * gene_name, char * transcript_name, char * chro_name, unsigned int start, unsigned int end, int is_negative_strand, void * context){
	int txn_to_free = 0;
	if(transcript_name == NULL){
		transcript_name = malloc(FEATURE_NAME_LENGTH+10);
		SUBreadSprintf(transcript_name, FEATURE_NAME_LENGTH+10,"TXN_%s", gene_name);
		txn_to_free = 1;
	}
	cellcounts_global_t * cct_context = context;
	ArrayList * the_features = cct_context -> all_features_array;
	fc_feature_info_t * new_added = calloc(sizeof(fc_feature_info_t), 1);

	if(cct_context -> sam_chro_to_anno_chr_alias){
		char * sam_chro = get_sam_chro_name_from_alias(cct_context -> sam_chro_to_anno_chr_alias, chro_name);
		if(sam_chro!=NULL) chro_name = sam_chro;
	}

	if(!cct_context -> transcript_exon_table){
		cct_context ->transcript_exon_table = StringTableCreate(100000);
		HashTableSetDeallocationFunctions(cct_context -> transcript_exon_table, free, (void (*)(void *))ArrayListDestroy);

		cct_context -> transcript_to_gene_name_table = StringTableCreate(100000);
		HashTableSetDeallocationFunctions(cct_context -> transcript_to_gene_name_table, free, free);
	}

	char * Rgene_name = HashTableGet(cct_context -> transcript_to_gene_name_table, transcript_name);
	if(!Rgene_name) HashTablePut(cct_context -> transcript_to_gene_name_table, strdup(transcript_name), strdup(gene_name));

	char txn_tab_key[MAX_CHROMOSOME_NAME_LEN+FEATURE_NAME_LENGTH+10];
	snprintf(txn_tab_key, MAX_CHROMOSOME_NAME_LEN+FEATURE_NAME_LENGTH+10, "%s\t%s\t%c", transcript_name, chro_name, is_negative_strand?'-':'+');
	ArrayList * txn_exons = HashTableGet(cct_context ->transcript_exon_table,txn_tab_key);
	if(!txn_exons){
		txn_exons = ArrayListCreate(30);
		HashTablePut(cct_context ->transcript_exon_table, strdup(txn_tab_key), txn_exons);
	}
	ArrayListPush(txn_exons,NULL+((start-1)*1LLU<<32)+end-1);

	char tmp_chro_name[MAX_CHROMOSOME_NAME_LEN];
	int access_n = HashTableGet(cct_context -> chromosome_table.read_name_to_index, chro_name ) - NULL;
	if(access_n < 1){
		if(chro_name[0]=='c' && chro_name[1]=='h' && chro_name[2]=='r'){
			chro_name += 3;
		}else{
			strcpy(tmp_chro_name, "chr");
			strcat(tmp_chro_name, chro_name);
			chro_name = tmp_chro_name;
		}
	}

	// for the featureCounts part.
	cct_context -> longest_chro_name = max(cct_context -> longest_chro_name, strlen(chro_name));
	new_added -> feature_name_pos = cellCounts_unistr_cpy(cct_context, gene_name, strlen(gene_name));
	new_added -> chro_name_pos_delta = cellCounts_unistr_cpy(cct_context, chro_name, strlen(chro_name)) - new_added -> feature_name_pos;
	new_added -> start = start;
	new_added -> end = end;
	new_added -> sorted_order = the_features -> numOfElements;
	new_added -> is_negative_strand = is_negative_strand;
	ArrayListPush(the_features, new_added);

	fc_chromosome_index_info * chro_stub = HashTableGet(cct_context -> chromosome_exons_table, chro_name);
	if(chro_stub){
		if(chro_stub -> chro_possible_length < end+1){
			if(chro_stub -> reverse_table_start_index){
				SUBREADprintf("ERROR: chromosome '%s' in the index is shorter than the same chromosome in the annotations: %d > %d.\n", chro_name, end, chro_stub -> chro_possible_length);
				return -1;
			}
			else chro_stub -> chro_possible_length = end+1;
		}
	}else{
		chro_stub = calloc(sizeof(fc_chromosome_index_info),1);
		char * tmp_chro_name = malloc(CHROMOSOME_NAME_LENGTH);
		term_strncpy(tmp_chro_name, chro_name, CHROMOSOME_NAME_LENGTH);
		chro_stub -> chro_number = cct_context -> chromosome_exons_table -> numOfElements;
		chro_stub -> chro_possible_length = end+1;
		chro_stub -> reverse_table_start_index_size = 0;
		chro_stub -> reverse_table_start_index = NULL;
		HashTablePut(cct_context -> chromosome_exons_table, tmp_chro_name, chro_stub);
	}
	chro_stub -> chro_features ++;

	if(chro_stub -> reverse_table_start_index){
		if(start >=  chro_stub -> reverse_table_start_index_size){
			if(!has_reverse_table_reported)SUBREADprintf("WARNING: size of chromosome in the reference genome is smaller than it in the annotation!\n");
			has_reverse_table_reported = 1;
		}else{
			int bin_location = start / REVERSE_TABLE_BUCKET_LENGTH;
			chro_stub -> reverse_table_start_index[bin_location]++;
		}
	}

	// for the align part.
	unsigned int exonic_map_start = linear_gene_position(&cct_context->chromosome_table , chro_name, start);
	unsigned int exonic_map_stop = linear_gene_position(&cct_context->chromosome_table , chro_name, end), exonpos_i;
	if(exonic_map_start > 0xffffff00 || exonic_map_stop > 0xffffff00){
		if(txn_to_free) free(transcript_name);
		return -1;
	}

	//SUBREADprintf("LINE1ADD %s of %s:%d ~ %d\n", gene_name, chro_name, start, end );
	for(exonpos_i = exonic_map_start; exonpos_i <= exonic_map_stop; exonpos_i++){
		int exonic_map_byte = exonpos_i/ 8;
		int exonic_map_bit = exonpos_i % 8;

		cct_context ->exonic_region_bitmap[exonic_map_byte] |= (1<<exonic_map_bit);
	}

	for(exonpos_i = exonic_map_start -100; exonpos_i <= exonic_map_stop +100; exonpos_i++){
		int exonic_map_byte = exonpos_i / 8 + (4096 / 8)*1024*1024;
		int exonic_map_bit =  exonpos_i % 8;

		//if(strcmp("ENSMUSG00000024608", gene_name)==0)SUBREADprintf("LINE1ADD %s: %d > %d\n", gene_name, exonic_map_byte , exonic_map_bit );
		cct_context ->exonic_region_bitmap[exonic_map_byte] |= (1<<exonic_map_bit);
	}


	cellCounts_register_junc_feature(cct_context, gene_name, transcript_name,chro_name, start, end, is_negative_strand); // The cellCounts part uses 0-based coordinates. The featureCOunts part uses 1-based coordinates.
	if(txn_to_free) free(transcript_name);
	return 0;
}

void cellCounts_register_junc_feature(cellcounts_global_t * cct_context, char * feature_name, char * transcript_id, char * chro, unsigned int start, unsigned int stop, int is_negative){
	int xk1, edge_tab_i, sample_i;
	char transcript_id_place[FEATURE_NAME_LENGTH+30];
	if(transcript_id==NULL){
			SUBreadSprintf(transcript_id_place, FEATURE_NAME_LENGTH+30,"%s%s", feature_name, CCT_NIL_TXN_PLACEHOLDER);
			transcript_id=transcript_id_place;
	}

	cct_junction_exon_in_transcript_t * new_item = calloc(sizeof(cct_junction_exon_in_transcript_t),1);
	new_item -> transcript_id = HashTableGetKey(cct_context -> junction_transcript_table, transcript_id);
	char * edgetab_used_gene_name=NULL;
	if(NULL == new_item -> transcript_id){
			cct_junction_transcript_t * new_txp = calloc(sizeof(cct_junction_transcript_t),1);
			new_txp -> exons_in_transcript = ArrayListCreate(10);
			ArrayListSetDeallocationFunction(new_txp -> exons_in_transcript, free);
			ArrayListPush(new_txp -> exons_in_transcript, new_item);
			new_txp -> gene_name = strdup(feature_name);
			new_txp -> transcript_id = strdup(transcript_id);
			new_item -> transcript_id = new_txp -> transcript_id;
			HashTablePut(cct_context -> junction_transcript_table, new_item->transcript_id, new_txp);
			edgetab_used_gene_name=new_txp -> gene_name;
	}else{
			cct_junction_transcript_t * has_txp = HashTableGet(cct_context -> junction_transcript_table, transcript_id);
			ArrayListPush(has_txp -> exons_in_transcript, new_item);
			edgetab_used_gene_name=has_txp -> gene_name;
	}

	new_item -> chro_start = start;
	new_item -> chro_stop = stop;
	new_item -> is_negative = is_negative;

if(0){
/*
	IVT_IntervalTreeNode * IVT_rootnode = HashTableGet(cct_context -> junction_ExonTree_table, chro);
	if(NULL == IVT_rootnode){
			IVT_rootnode = IVT_insert(NULL, start, stop, new_item); // JCount module uses 1-based coordinates.
			char * new_name = strdup(chro);
			HashTablePut(cct_context -> junction_ExonTree_table, new_name, IVT_rootnode);
	}else{
			IVT_rootnode = IVT_insert(IVT_rootnode, start, stop, new_item); // JCount module uses 1-based coordinates.
			HashTablePutReplaceEx(cct_context -> junction_ExonTree_table, chro, IVT_rootnode,0,0,0); // use old key ptr in table; don't free key and value.
	}
*/

}

	IVT_IntervalTreeNode * IVT_edgenode = HashTableGet(cct_context -> junction_ExonEdgeTree_table[0], chro);
	IVT_edgenode = IVT_insert(IVT_edgenode, start, start, edgetab_used_gene_name);
	IVT_edgenode = IVT_insert(IVT_edgenode, stop, stop, edgetab_used_gene_name);
	HashTablePutReplaceEx(cct_context -> junction_ExonEdgeTree_table[0], strdup(chro), IVT_edgenode,1,1,0); // use old key ptr in table; don't free key and value.

	for(edge_tab_i=1; edge_tab_i<=2; edge_tab_i++){
			if(is_negative != ( edge_tab_i==2 )) continue;
			IVT_edgenode = HashTableGet(cct_context -> junction_ExonEdgeTree_table[edge_tab_i], chro);
			IVT_edgenode = IVT_insert(IVT_edgenode, start, start, edgetab_used_gene_name);
			IVT_edgenode = IVT_insert(IVT_edgenode, stop, stop, edgetab_used_gene_name);
			HashTablePutReplaceEx(cct_context -> junction_ExonEdgeTree_table[edge_tab_i], strdup(chro), IVT_edgenode,1,1,0);
	}

	// create gene body data structure for better junction reporting
	char gene_body_key[FEATURE_NAME_LENGTH + CHROMOSOME_NAME_LENGTH];
	SUBreadSprintf(gene_body_key, FEATURE_NAME_LENGTH + CHROMOSOME_NAME_LENGTH ,"%s\t%s%s", feature_name, chro, is_negative?"NEG":"POS");
	cct_junction_genebody_t * jgbody = HashTableGet(cct_context -> junction_genebody_table, gene_body_key);
	if(!jgbody){
		jgbody = malloc(sizeof(cct_junction_genebody_t));
		memset(jgbody,0,sizeof(cct_junction_genebody_t));
		strcpy(jgbody -> gene_name, feature_name);
		strcpy(jgbody -> chromosome_name, chro);
		jgbody -> transcript_list = ArrayListCreate(10); 
		// the items are destroyed in  global_context -> junction_transcript_table[txid] => exons_in_transcript
		HashTablePut(cct_context -> junction_genebody_table, strdup(gene_body_key), jgbody);
	}

	int transcript_had_in_gene = 0;
	for(xk1 = 0; xk1 < jgbody ->transcript_list -> numOfElements; xk1++){
		cct_junction_transcript_t * had_txp = ArrayListGet(jgbody -> transcript_list, xk1);
		if(strcmp(transcript_id, had_txp -> transcript_id)==0) transcript_had_in_gene = 1;
	}
	if(!transcript_had_in_gene) ArrayListPush(jgbody -> transcript_list, HashTableGet(cct_context -> junction_transcript_table, transcript_id));
}

int cellCounts_find_or_insert_gene_name(cellcounts_global_t * cct_context, unsigned char * feature_name) {
	HashTable * genetable = cct_context -> gene_name_table;

	srInt_64 gene_number = HashTableGet(genetable, feature_name) - NULL;
	if(gene_number>0)
		return gene_number-1;
	else {
		gene_number = genetable -> numOfElements; 
		HashTablePut(genetable, feature_name, NULL+gene_number+1);
		cct_context -> gene_name_array[gene_number] = feature_name;
			// real memory space of feature_name is in the "loaded_features" data structure.
			// now we only save its pointer.

		return gene_number;
	}
}

void cellCounts_register_reverse_table(int block_no, srInt_64 this_block_min_start, srInt_64 this_block_max_end, fc_chromosome_index_info * chro_inf) {
	unsigned int reversed_bucket_start = this_block_min_start /  REVERSE_TABLE_BUCKET_LENGTH;
	unsigned int reversed_bucket_end = this_block_max_end / REVERSE_TABLE_BUCKET_LENGTH;
	assert(this_block_min_start <= this_block_max_end);
	assert(reversed_bucket_end < chro_inf -> chro_possible_length);
	if(!chro_inf->reverse_table_start_index)return;
	int x1;
	for(x1 = reversed_bucket_start; x1 <= reversed_bucket_end; x1++)
		chro_inf->reverse_table_start_index[x1] = min(chro_inf->reverse_table_start_index[x1], block_no);
}

void cellCounts_feature_merge(void * arrv, int start, int items, int items2)
{

	void ** arr = (void **) arrv;

	srInt_64 * ret_start = (srInt_64 *) arr[0];
	srInt_64 * ret_end = (srInt_64 *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	int total_items = items+items2;
	srInt_64 * tmp_start = malloc(sizeof(srInt_64) * total_items);
	srInt_64 * tmp_end = malloc(sizeof(srInt_64) * total_items);
	unsigned char * tmp_strand = malloc(sizeof(char) * total_items);
	int * tmp_entyrez = malloc(sizeof(int) * total_items);
	fc_feature_info_t ** tmp_info_ptr = malloc(sizeof(fc_feature_info_t*) * total_items);

	int read_1_ptr = start;
	int read_2_ptr = start+items;
	int write_ptr;

	for(write_ptr=0; write_ptr<total_items; write_ptr++)
	{
		if((read_1_ptr >= start+items)||(read_2_ptr < start+total_items && ret_start[read_1_ptr] >= ret_start[read_2_ptr]))
		{
			tmp_start[write_ptr] = ret_start[read_2_ptr];
			tmp_end[write_ptr] = ret_end[read_2_ptr];
			tmp_strand[write_ptr] = ret_strand[read_2_ptr];
			tmp_entyrez[write_ptr] = ret_entyrez[read_2_ptr];
			tmp_info_ptr[write_ptr] = old_info_ptr[read_2_ptr];
			read_2_ptr++;
		}
		else
		{
			tmp_start[write_ptr] = ret_start[read_1_ptr];
			tmp_end[write_ptr] = ret_end[read_1_ptr];
			tmp_strand[write_ptr] = ret_strand[read_1_ptr];
			tmp_entyrez[write_ptr] = ret_entyrez[read_1_ptr];
			tmp_info_ptr[write_ptr] = old_info_ptr[read_1_ptr];
			read_1_ptr++;
		}
	}

	memcpy(ret_start+ start, tmp_start, sizeof(srInt_64) * total_items);
	memcpy(ret_end+ start, tmp_end, sizeof(srInt_64) * total_items);
	memcpy(ret_strand+ start, tmp_strand, sizeof(char) * total_items);
	memcpy(ret_entyrez+ start, tmp_entyrez, sizeof(int) * total_items);
	memcpy(old_info_ptr+ start, tmp_info_ptr, sizeof(fc_feature_info_t*) * total_items);

	free(tmp_start);
	free(tmp_end);
	free(tmp_strand);
	free(tmp_entyrez);
	free(tmp_info_ptr);
}


int cellCounts_feature_sort_compare(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	srInt_64 * ret_start = (srInt_64 *)arr[0];
	srInt_64 ll = ret_start[l];
	srInt_64 rl = ret_start[r];

	if(ll==rl) return 0;
	else if(ll>rl) return 1;
	else return -1;
}

void cellCounts_feature_sort_exchange(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	srInt_64 tmp;
	fc_feature_info_t * tmpptr;

	srInt_64 * ret_start = (srInt_64 *) arr[0];
	srInt_64 * ret_end = (srInt_64 *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	
	tmp = ret_start[r];
	ret_start[r]=ret_start[l];
	ret_start[l]=tmp;

	tmp = ret_end[r];
	ret_end[r]=ret_end[l];
	ret_end[l]=tmp;

	tmp = ret_strand[r];
	ret_strand[r]=ret_strand[l];
	ret_strand[l]=tmp;

	tmp = ret_entyrez[r];
	ret_entyrez[r]=ret_entyrez[l];
	ret_entyrez[l]=tmp;

	tmpptr = old_info_ptr[r];
	old_info_ptr[r]=old_info_ptr[l];
	old_info_ptr[l]=tmpptr;

}



void cellCounts_sort_feature_info(cellcounts_global_t * cct_context, unsigned int features, ArrayList * loaded_features, char *** sorted_chr_names, int ** sorted_entrezid, srInt_64 ** sorted_start, srInt_64 ** sorted_end, unsigned char ** sorted_strand, srInt_64 ** block_end_index, srInt_64 ** block_min_start_pos, srInt_64 ** block_max_end_pos) {
	unsigned int chro_pnt;
	unsigned int xk1,xk2;
	int * ret_entrez = malloc(sizeof(int) * features);
	srInt_64 * ret_start = malloc(sizeof(srInt_64) * features);
	srInt_64 * ret_end = malloc(sizeof(srInt_64) * features);
	int current_block_buffer_size = 2000;

	srInt_64 * ret_block_end_index = malloc(sizeof(srInt_64) * current_block_buffer_size);
	srInt_64 * ret_block_min_start = malloc(sizeof(srInt_64) * current_block_buffer_size);
	srInt_64 * ret_block_max_end = malloc(sizeof(srInt_64) * current_block_buffer_size);
	unsigned char * ret_strand = malloc(features);
	char ** ret_char_name = malloc(sizeof(void *) * features);
	fc_feature_info_t ** old_info_ptr = malloc(sizeof(void *) * features);
	unsigned int * chro_feature_ptr = calloc(sizeof(int) , cct_context -> chromosome_exons_table -> numOfElements);
	fc_chromosome_index_info ** tmp_chro_info_ptrs = calloc(sizeof(fc_chromosome_index_info *), cct_context -> chromosome_exons_table -> numOfElements);

	cct_context -> gene_name_array = malloc(sizeof(char *) * features);	// there should be much less identical names.
	cct_context -> gene_name_table = HashTableCreate(5000);
	HashTableSetHashFunction(cct_context -> gene_name_table, HashTableStringHashFunction);
	HashTableSetKeyComparisonFunction(cct_context -> gene_name_table, my_strcmp);

	// init start positions of each chromosome block.
	if(1) {
		KeyValuePair * cursor;
		int bucket;
		unsigned int sum_ptr = 0;
		for(bucket=0; bucket < cct_context -> chromosome_exons_table -> numOfBuckets; bucket++) {
			cursor = cct_context -> chromosome_exons_table -> bucketArray[bucket];
			while (1) {
				if (!cursor) break;
				fc_chromosome_index_info * tmp_chro_inf = cursor -> value;
				cursor = cursor->next;
				assert(tmp_chro_inf -> chro_number <  cct_context -> chromosome_exons_table -> numOfElements);
				chro_feature_ptr [tmp_chro_inf -> chro_number] = tmp_chro_inf -> chro_features;
				tmp_chro_info_ptrs[tmp_chro_inf -> chro_number] = tmp_chro_inf;
			}
		}

		for(xk1 = 0; xk1 < cct_context -> chromosome_exons_table -> numOfElements; xk1++) {
			unsigned int tmpv = chro_feature_ptr[xk1];
			chro_feature_ptr[xk1] = sum_ptr;
			tmp_chro_info_ptrs[xk1] -> chro_feature_table_start = sum_ptr;
			sum_ptr += tmpv;
		}
	}
	int current_block_id = 0, sort_i = 0;

	(*sorted_chr_names) = ret_char_name;
	(*sorted_entrezid) = ret_entrez;
	(*sorted_start) = ret_start;
	(*sorted_end) = ret_end;
	(*sorted_strand) = ret_strand;

	for(chro_pnt=0; chro_pnt < features; chro_pnt++) {
		fc_feature_info_t * cur_feature = ArrayListGet(loaded_features, chro_pnt);
		char * this_chro_name = cct_context -> unistr_buffer_space + cur_feature -> feature_name_pos + cur_feature -> chro_name_pos_delta;
		
		fc_chromosome_index_info * chro_stub = HashTableGet(cct_context -> chromosome_exons_table, this_chro_name);
		int this_chro_number = chro_stub -> chro_number;
		assert( this_chro_number < cct_context -> chromosome_exons_table -> numOfElements );
		unsigned int this_chro_table_ptr = chro_feature_ptr[this_chro_number];

		ret_char_name[this_chro_table_ptr] = this_chro_name;// (char *)cur_feature -> chro;
		ret_entrez[this_chro_table_ptr] = cellCounts_find_or_insert_gene_name(cct_context, (unsigned char *)(cct_context -> unistr_buffer_space + cur_feature -> feature_name_pos));
		ret_start[this_chro_table_ptr] = cur_feature -> start;
		ret_end[this_chro_table_ptr] = cur_feature -> end;
		ret_strand[this_chro_table_ptr] = cur_feature -> is_negative_strand;
		old_info_ptr[this_chro_table_ptr] = cur_feature;

		chro_feature_ptr[this_chro_number]++;
	}

	print_in_box(80,0,0,"Sort the %d genes...", cct_context -> gene_name_table -> numOfElements);
	for(xk1 = 0; xk1 < cct_context -> chromosome_exons_table -> numOfElements; xk1++) {
		fc_chromosome_index_info * tmp_chro_inf = tmp_chro_info_ptrs[xk1];
		int bins_in_chr = ( tmp_chro_inf->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2);
		short * features_per_block_bins = malloc(sizeof(short)*bins_in_chr);

		if(tmp_chro_inf -> reverse_table_start_index){
			for(xk2=0; xk2<bins_in_chr; xk2++)
				features_per_block_bins[xk2] = max(1,min(1000,(int)(0.9999999+sqrt(tmp_chro_inf -> reverse_table_start_index[xk2]))));
			memset(tmp_chro_inf -> reverse_table_start_index, 0xff, sizeof(int) *bins_in_chr);
		}else continue;

		unsigned int this_block_items = 0;
		srInt_64 this_block_min_start = 0x7fffffff, this_block_max_end = 0;
		unsigned int this_chro_tab_end =  tmp_chro_inf -> chro_features + tmp_chro_inf -> chro_feature_table_start;

		void * in_array[5];
		in_array[0] = ret_start + tmp_chro_inf -> chro_feature_table_start; 
		in_array[1] = ret_end + tmp_chro_inf -> chro_feature_table_start; 
		in_array[2] = ret_strand + tmp_chro_inf -> chro_feature_table_start; 
		in_array[3] = ret_entrez + tmp_chro_inf -> chro_feature_table_start; 
		in_array[4] = old_info_ptr + tmp_chro_inf -> chro_feature_table_start; 

		merge_sort(in_array, this_chro_tab_end - tmp_chro_inf -> chro_feature_table_start, cellCounts_feature_sort_compare, cellCounts_feature_sort_exchange, cellCounts_feature_merge);

		for(sort_i = tmp_chro_inf -> chro_feature_table_start; sort_i< this_chro_tab_end ; sort_i++) {
			old_info_ptr[sort_i]->sorted_order = sort_i;
			int feature_bin_location = ret_start[sort_i] / REVERSE_TABLE_BUCKET_LENGTH;
			int block_bin_location = this_block_min_start / REVERSE_TABLE_BUCKET_LENGTH;

			if(this_block_items && (this_block_items > features_per_block_bins[block_bin_location] || feature_bin_location != block_bin_location)){
				if(current_block_id >= current_block_buffer_size - 1) {
					current_block_buffer_size *= 1.3;
					ret_block_min_start = realloc(ret_block_min_start, sizeof(srInt_64)*current_block_buffer_size);
					ret_block_max_end = realloc(ret_block_max_end, sizeof(srInt_64)*current_block_buffer_size);
					ret_block_end_index = realloc(ret_block_end_index, sizeof(srInt_64)*current_block_buffer_size);
				}
				ret_block_end_index[current_block_id] = sort_i;	// FIRST UNWANTED ID
				ret_block_min_start[current_block_id] = this_block_min_start;
				ret_block_max_end[current_block_id] = this_block_max_end;
				cellCounts_register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
				current_block_id++;
				this_block_max_end = 0;
				this_block_items = 0;
				this_block_min_start = 0x7fffffff;
			}

			this_block_max_end = max(this_block_max_end, ret_end[sort_i]);
			this_block_min_start = min(this_block_min_start, ret_start[sort_i]);
			this_block_items ++;
		
		}

		if(this_block_items) {
			if(current_block_id >= current_block_buffer_size) {
				current_block_buffer_size *= 1.3;
				ret_block_min_start = realloc(ret_block_min_start, sizeof(srInt_64)*current_block_buffer_size);
				ret_block_max_end = realloc(ret_block_max_end, sizeof(srInt_64)*current_block_buffer_size);
				ret_block_end_index = realloc(ret_block_end_index, sizeof(srInt_64)*current_block_buffer_size);
			}

			ret_block_end_index[current_block_id] = this_chro_tab_end;	// FIRST UNWANTED ID
			ret_block_min_start[current_block_id] = this_block_min_start;
			ret_block_max_end[current_block_id] = this_block_max_end;
			cellCounts_register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
			current_block_id++;
		}

		tmp_chro_inf -> chro_block_table_end = current_block_id; 
		free(features_per_block_bins);
	}

	(*block_end_index) = ret_block_end_index;
	(*block_min_start_pos) = ret_block_min_start;
	(*block_max_end_pos) = ret_block_max_end;

	free(old_info_ptr);
	free(tmp_chro_info_ptrs);
	free(chro_feature_ptr);
}



int cellCounts_load_annotations(cellcounts_global_t * cct_context){
	int rv = 0;
	unsigned int last_loc = 0;
	has_reverse_table_reported = 0;
			
	if(cct_context -> features_annotation_alias_file[0])
		rv = (NULL != (cct_context -> sam_chro_to_anno_chr_alias = load_alias_table(cct_context->features_annotation_alias_file)));
	if(!rv){
		int x1;
		cct_context -> unistr_buffer_size = 1024*1024*2;
		cct_context -> unistr_buffer_space = malloc(cct_context -> unistr_buffer_size);
		cct_context -> chromosome_exons_table = HashTableCreate(163);
		HashTableSetHashFunction(cct_context -> chromosome_exons_table, HashTableStringHashFunction);
		HashTableSetKeyComparisonFunction(cct_context -> chromosome_exons_table, my_strcmp);

		for(x1 = 0; x1 < cct_context -> chromosome_table.total_offsets; x1++){
			fc_chromosome_index_info * chro_stub = calloc(sizeof(fc_chromosome_index_info),1);
			char * tmp_chro_name = malloc(CHROMOSOME_NAME_LENGTH);
			char * seq_name = cct_context -> chromosome_table.read_names + MAX_CHROMOSOME_NAME_LEN*x1;
			term_strncpy(tmp_chro_name, seq_name, CHROMOSOME_NAME_LENGTH);
			chro_stub -> chro_number = HashTableGet(cct_context -> chromosome_table.read_name_to_index, seq_name ) - NULL -1;
			chro_stub -> chro_possible_length = cct_context -> chromosome_table.read_offsets[x1] - last_loc;
			last_loc = cct_context -> chromosome_table.read_offsets[x1];
			chro_stub -> reverse_table_start_index_size = chro_stub -> chro_possible_length + 1024*1024;
			chro_stub -> reverse_table_start_index = calloc( chro_stub -> reverse_table_start_index_size  / REVERSE_TABLE_BUCKET_LENGTH +2, sizeof(int));
			HashTablePut(cct_context -> chromosome_exons_table, tmp_chro_name, chro_stub);
		}
		
		cct_context -> all_features_array = ArrayListCreate(350000);
		ArrayListSetDeallocationFunction(cct_context -> all_features_array, free);
		int loaded_features = load_features_annotation(cct_context->features_annotation_file, cct_context->features_annotation_file_type, cct_context->features_annotation_gene_id_column, "transcript_id" , cct_context-> features_annotation_feature_type, cct_context, features_load_one_line);
		if(loaded_features<1) rv = 1;

		if(!rv) rv = cellCounts_extract_and_sort_juncs(cct_context); // transcripts are copied to junctions and saved in the junction table.

		if(!rv){
			int anno_index_matched=0;
			ArrayList * annot_chros = HashTableKeys(cct_context -> chromosome_exons_table);
			for(x1=0; x1<annot_chros -> numOfElements; x1++){
				char * t1chro = (char*)ArrayListGet(annot_chros,x1);
				fc_chromosome_index_info * chro_stub = HashTableGet(cct_context -> chromosome_exons_table, t1chro);
				if(chro_stub -> chro_features<1) ArrayListSet(annot_chros,x1,NULL);
			}
			int all_chro_unmatched = warning_array_hash_numbers(annot_chros, cct_context-> chromosome_table.read_name_to_index, & anno_index_matched);
			rv=all_chro_unmatched;
			ArrayListDestroy(annot_chros);

			if(all_chro_unmatched) SUBREADprintf("ERROR: no matched chromosomes/contigs found between reference sequences and gene annotation.\n"); else{
				char tbuf[90];
				char_strftime(tbuf);
				SUBREADprintf("Number of chromosomes/contigs matched between reference sequences and gene annotation is %d.\n\n", anno_index_matched);
				cellCounts_print_config(cct_context);
				print_in_box(80,1,1,"Running (%s, pid=%d)", tbuf, getpid());
				print_in_box(80,0,0,"");
			}
			if(!rv) cellCounts_sort_feature_info(cct_context, loaded_features, cct_context -> all_features_array, &cct_context -> features_sorted_chr, &cct_context -> features_sorted_geneid, &cct_context -> features_sorted_start, &cct_context -> features_sorted_stop, &cct_context -> features_sorted_strand, &cct_context -> block_end_index, &cct_context -> block_min_start, &cct_context -> block_max_end);
		}
	}
	return rv;
}

int cellCounts_open_cellbc_batches(cellcounts_global_t * cct_context){
	int x1;
	for(x1=0;x1<CELLBC_BATCH_NUMBER+2; x1++){
		char fname[MAX_FILE_NAME_LENGTH+200];
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+200,"%s/temp-cellcounts-%06d-%03d.tmpbin",cct_context->temp_file_dir,getpid(), x1);
		cct_context -> batch_files[x1] = REP_fopen(fname,"wb");
		REP_setvbuf(cct_context -> batch_files[x1], cct_context -> cellbin_v_buffers[x1], _IOFBF, SCRNA_SMALLER_VBUFF_SIZE);
		cellCounts_init_lock(cct_context -> batch_file_locks+x1, 0);
	}
	int umfpi;
	if(cct_context->input_mode == GENE_INPUT_BCL)for(umfpi=1; umfpi<=4; umfpi++){
		char fname [MAX_FILE_NAME_LENGTH+20];
		char * ftype = "R1";
		if(3==umfpi && !cct_context->is_dual_index) continue;
		if(2==umfpi) ftype = "I1";
		if(3==umfpi) ftype = "I2";
		if(4==umfpi) ftype = "R2";
		SUBreadSprintf(fname, MAX_FILE_NAME_LENGTH+20, "UnassignedReads%03d_%s.fastq.gz", cct_context ->current_dataset_no, ftype);
		parallel_gzip_writer_init(cct_context-> fastq_unassigned_writer+(umfpi-1), fname, cct_context->total_threads);
	}
	cellCounts_init_lock(&cct_context->fastq_unassigned_lock, 0);
	return 0;
}

// we will never use spin-lock.
void cellCounts_init_lock(cellCounts_lock_t * lock, int is_spin){
	pthread_mutex_init(lock, NULL);
}

void cellCounts_destroy_lock(cellCounts_lock_t * lock){
	//if(lock -> is_spin_lock) pthread_spin_destroy(&lock->spinlock);
	//	else
	 pthread_mutex_destroy(lock);
}

/*
int cellCounts_lock_occupy(cellCounts_lock_t * lock){
	//if(lock -> is_spin_lock)
	//	return pthread_spin_lock(&lock->spinlock);
	//else
	return pthread_mutex_lock(lock);
}

int cellCounts_lock_release(cellCounts_lock_t * lock){
	//if(lock -> is_spin_lock) return pthread_spin_unlock(&lock->spinlock);
	//else
	return pthread_mutex_unlock(lock);
}
*/

#define CACHED_BCL_READ_NUMBER 0x1000000 // 16.78M

int cellCounts_open_input_fps(cellcounts_global_t * cct_context){
	int rv=0;
	if(cct_context -> input_mode == GENE_INPUT_BCL){
		rv = rv || geinput_open_bcl(cct_context -> input_dataset_name , & cct_context -> input_dataset , /*cct_context -> reads_per_chunk*/ CACHED_BCL_READ_NUMBER, cct_context -> total_threads);
		if(!rv)cct_context -> is_dual_index = cct_context -> input_dataset.bcl_input.is_dual_index;
	} else if(cct_context -> input_mode == GENE_INPUT_SCRNA_FASTQ)
		rv = rv || geinput_open_scRNA_fqs(cct_context -> input_dataset_name , & cct_context -> input_dataset , cct_context -> reads_per_chunk, cct_context -> total_threads);
	else if(cct_context -> input_mode == GENE_INPUT_SCRNA_BAM)
		rv = rv || geinput_open_scRNA_BAM(cct_context -> input_dataset_name , & cct_context -> input_dataset , cct_context -> reads_per_chunk, cct_context -> total_threads);
	else 	rv = -1;
	return rv;
}

void cellCounts_2IVT_freeTree(void * c2){
	void ** cc2 = c2;
	IVT_freeTree(cc2[0]);
	IVT_freeTree(cc2[1]);
	IVT_freeTree(cc2[2]);
}

void cellCounts_chroEvent_freebuff(void *c2){
	chroEvent_t * de = c2;
	if(de -> event_type== chroEvent_t_TYPE_INDEL){
		int inslen = (de -> n_events>1)? de -> length[0] :((void*)(de -> length) - NULL);
		if(inslen<0) free(de ->inserted_bases);
	}
	if(de -> n_events>1)free(de->length);
	free(de);
}

void cellCounts_junction_genebody_table_free(void * jcv){
	cct_junction_genebody_t * jc = jcv;
	ArrayListDestroy(jc -> transcript_list);
 	free(jc);
}

void cellCounts_junc_transcript_free(cct_junction_transcript_t * txp){
	free(txp -> gene_name);
	free(txp -> transcript_id);
	ArrayListDestroy(txp -> exons_in_transcript);
	free(txp);
}

int cellCounts_init_junction_related_context(cellcounts_global_t * cct_context){
	cellCounts_init_lock(&cct_context -> chroEvent_entry_table_lock, 0);
	int x1, sample_i;
	cct_context -> read_assignment_counter_locks = malloc(sizeof(cellCounts_lock_t)* cct_context -> chroEvent_lock_number);
	for(x1=0; x1 < cct_context -> chroEvent_lock_number; x1++) cellCounts_init_lock( cct_context -> read_assignment_counter_locks + x1,0);

	for(sample_i = 1; sample_i <=cct_context-> sample_sheet_table -> numOfElements ; sample_i ++){
		cct_context -> chroEvent_entry_table[sample_i] = StringTableCreate(200);
		HashTableSetDeallocationFunctions(cct_context -> chroEvent_entry_table[sample_i], free, cellCounts_2IVT_freeTree);

		cct_context -> chroEvent_detail_table[sample_i] = HashTableCreate(400000);
		HashTableSetDeallocationFunctions(cct_context -> chroEvent_detail_table[sample_i], NULL, cellCounts_chroEvent_freebuff);

		cct_context -> junction_to_cell_umi_table[sample_i] = HashTableCreate(100000);
		HashTableSetDeallocationFunctions(cct_context -> junction_to_cell_umi_table[sample_i] , NULL,  (void (*)(void *value))ArrayListDestroy);
	}

	cct_context -> junction_genebody_table = StringTableCreate(1603);
	HashTableSetDeallocationFunctions(cct_context -> junction_genebody_table, free, (void (*)(void *))cellCounts_junction_genebody_table_free); 

	cct_context -> junction_transcript_table = StringTableCreate(1603);
	HashTableSetDeallocationFunctions(cct_context -> junction_transcript_table, NULL, (void (*)(void *))cellCounts_junc_transcript_free); 


	for(x1=0; x1<3; x1++){
		cct_context -> junction_ExonEdgeTree_table[x1] = StringTableCreate(1603);
		HashTableSetDeallocationFunctions(cct_context -> junction_ExonEdgeTree_table[x1], free, (void (*)(void *))IVT_freeTree);
	}

	cct_context -> junction_GenebodyTree_table = StringTableCreate(1603);
	HashTableSetDeallocationFunctions(cct_context -> junction_GenebodyTree_table, NULL, (void (*)(void *))IVT_freeTree);
	return 0;
}

#define have_sample_i_calc \
                for(x1=0; x1< cct_context -> sample_id_to_name -> numOfElements; x1++)if(strcmp((char*) ArrayListGet(cct_context -> sample_id_to_name, x1) , sample_name ) ==0){ sample_i=x1+1;break;}\
                if(sample_i<=0){\
                        SUBREADprintf("ERROR: unknown sample name in junc table: %s\n", sample_name);\
                        return -1;\
                }

int cellCounts_load_cluster_related_junctions_n_cells(cellcounts_global_t * cct_context){
	int x1;
	if(cct_context -> cluster_junctions_file[0] || cct_context -> cluster_map_file[0]){
		if(!(cct_context -> cluster_junctions_file[0] && cct_context -> cluster_map_file[0] )){
			SUBREADprintf("ERROR: only gave one of the cluster-junction table or the cluster-cell mapping table.\n");
			return -1; 
		}
	}else return 0; // no junction neither cell map is given.

	cct_context -> cluster_spec_junction_table = HashTableCreate(613); // sample_no << 56 | cluster_no => table
	HashTableSetDeallocationFunctions(cct_context -> cluster_spec_junction_table, NULL, (void *)(void *)HashTableDestroy);
	cct_context -> cluster_cell_map_table = HashTableCreate(76543);

	FILE * junc_fp = fopen(cct_context -> cluster_junctions_file,"r");
	FILE * map_fp = fopen(cct_context -> cluster_map_file,"r");

	char linebuf[MAX_CHROMOSOME_NAME_LEN+32+FEATURE_NAME_LENGTH], chname_strand[MAX_CHROMOSOME_NAME_LEN+3];

	while(1){
		char * rv = fgets(linebuf, MAX_CHROMOSOME_NAME_LEN+31+FEATURE_NAME_LENGTH, junc_fp);
		if(!rv)break;

		char * rvtmp = NULL;
		char * cluster_name = strtok_r(linebuf, "\t", &rvtmp);
		char * junc_desc = strtok_r(NULL, "\t", &rvtmp);
		char * sample_name = strtok_r(NULL, "\t", &rvtmp);
		sample_name[strlen(sample_name)-1]=0;
		int sample_i = -1;

		rvtmp = NULL;
		char * chro_name = strtok_r(junc_desc,":", &rvtmp);
		int start_pos = atoi(strtok_r(NULL,":", &rvtmp)) - 1;
		int end_pos = atoi(strtok_r(NULL,":", &rvtmp)) - 1; // internal pos: 0-based, even for the in-chro pos.

		//for(x1=0; x1< cct_context -> sample_id_to_name -> numOfElements; x1++) SUBREADprintf("CHsample '%s' == '%s' -> %d\n", sample_name, ArrayListGet(cct_context -> sample_id_to_name, x1), x1 );

		have_sample_i_calc;

		unsigned int absstart_pos = linear_gene_position(&cct_context->chromosome_table, chro_name , start_pos);
		unsigned int absend_pos = linear_gene_position(&cct_context->chromosome_table, chro_name , end_pos);

		char Lb2[2], Rb2[2];
		Lb2[0] = gvindex_get(cct_context -> value_index, absstart_pos +1);
		Lb2[1] = gvindex_get(cct_context -> value_index, absstart_pos +2);
		Rb2[0] = gvindex_get(cct_context -> value_index, absend_pos -2);
		Rb2[1] = gvindex_get(cct_context -> value_index, absend_pos -1);

		int known = cellCounts_junction_in_table(cct_context,sample_i,chro_name,start_pos,end_pos);// cellCounts_junction_in_table's l and r are 0 based
		int is_negative = -1;
		if(Lb2[0]=='G' && Lb2[1]=='T' && Rb2[0]=='A' && Rb2[1]=='G') is_negative=0;
		else if(Lb2[0]=='C' && Lb2[1]=='T' && Rb2[0]=='A' && Rb2[1]=='C') is_negative=1;
		else if(!known) SUBREADprintf("Warning: junction has a non-canonical donor-acceptor pair.  %s : %d - %d.\n", chro_name, start_pos, end_pos);

if(0)if(is_negative>=0) SUBREADprintf("Normal junction accpted: %s : %d - %d.\n" , chro_name, start_pos, end_pos);

		snprintf(chname_strand, MAX_CHROMOSOME_NAME_LEN+3,"%s\t%c", chro_name, is_negative?'-':'+');

		//cellCounts_ins_IVT_for_a_junc(cct_context, chname_strand, start_pos, end_pos, NULL+IMPOSSIBLE_MEMORY_SPACE, is_negative);
		if(!known) cellCounts_add_or_update_chroEvent_in_table(cct_context, sample_i, chroEvent_t_TYPE_JUNCTION, chname_strand , absstart_pos, absend_pos , 0, 0);

		int cluster_no = atoi(cluster_name+7); // "Cluster01"
		HashTable * this_cluster_table = (HashTable*)HashTableGet(cct_context -> cluster_spec_junction_table, NULL+(sample_i*1LLU<<56)+ cluster_no);
		if(!this_cluster_table) {
			this_cluster_table = HashTableCreate(76543);
			HashTablePut(cct_context -> cluster_spec_junction_table, NULL+(sample_i*1LLU<<56)+ cluster_no , this_cluster_table);
		}
		HashTablePut(this_cluster_table, NULL+( absstart_pos*1LLU<<32 ) + absend_pos , NULL+1);
	}
	fclose(junc_fp);

	while(1){
		char * rv = fgets(linebuf, MAX_CHROMOSOME_NAME_LEN+31+FEATURE_NAME_LENGTH, map_fp);
		if(!rv) break;
		char * rvtmp = NULL;
		char * cluster_name = strtok_r(linebuf, "\t", &rvtmp);
		char * cellbc = strtok_r(NULL, "\t", &rvtmp);
		char * sample_name = strtok_r(NULL, "\t", &rvtmp);
		sample_name[strlen(sample_name)-1]=0;
		int sample_i=-1;

		have_sample_i_calc;
		int cluster_no = atoi(cluster_name+7); // "Cluster01"
		int cellbc_no = HashTableGet(cct_context -> cell_barcode_head_tail_table, cellbc) -NULL - IMPOSSIBLE_MEMORY_SPACE;
		HashTablePut(cct_context -> cluster_cell_map_table, NULL+(sample_i*1LLU<<56)+cellbc_no, NULL+cluster_no);
	}
	fclose(map_fp);
	return 0;
}

int cellCounts_load_context(cellcounts_global_t * cct_context){
	int rv = 0;
	cellCounts_init_lock(&cct_context -> input_dataset_lock, 1 || (cct_context -> input_mode == GENE_INPUT_BCL));
	if(cct_context -> read_assignment_detail_file[0]) {
		cellCounts_init_lock(&cct_context -> read_assignment_detail_lock,0);
		cct_context -> read_assignment_detail_fp = fopen(cct_context -> read_assignment_detail_file,"w");
	}

	rv = rv || cellCounts_open_input_fps(cct_context);
	rv = rv || load_offsets(& cct_context -> chromosome_table, cct_context -> index_prefix);
	rv = rv || determine_total_index_blocks(cct_context);
	int bitmap_size = (4096 / 8)*1024*1024 *2; // the last "*2" is for the extended exon regions
	rv = rv || ((cct_context -> exonic_region_bitmap = calloc(bitmap_size, 1))==NULL);
	rv = rv || cellCounts_load_base_value_indexes(cct_context);
	rv = rv || cellCounts_load_scRNA_tables(cct_context);
	rv = rv || cellCounts_init_junction_related_context(cct_context);
	rv = rv || cellCounts_load_annotations(cct_context);
	rv = rv || cellCounts_open_cellbc_batches(cct_context);
	rv = rv || cellCounts_load_cluster_related_junctions_n_cells(cct_context); // this step will add some new junctions into the junction table load from annotation.

	return rv;
}

void * delete_file_thread(void * arg);
int cellCounts_destroy_context(cellcounts_global_t * cct_context){
	int x1, sample_i;
	pthread_join(cct_context ->thread_delete_files,NULL);
	for(x1=0;x1<CELLBC_BATCH_NUMBER+2; x1++)
		cellCounts_destroy_lock(cct_context -> batch_file_locks+x1);
	cellCounts_destroy_lock(&cct_context -> input_dataset_lock);
	cellCounts_destroy_lock(&cct_context -> chroEvent_entry_table_lock);
	for(x1=0; x1 < cct_context -> chroEvent_lock_number; x1++) cellCounts_destroy_lock( cct_context -> read_assignment_counter_locks + x1);

	if(cct_context -> is_BAM_and_FQ_out_generated){
		HashTableDestroy(cct_context->sample_BAM_writers);
		cellCounts_destroy_lock(&cct_context->fastq_unassigned_lock);
		int umfpi;
		if(cct_context->input_mode == GENE_INPUT_BCL)
			for(umfpi = 0; umfpi < 4; umfpi++)
				if(umfpi!=2 || cct_context->is_dual_index)
					parallel_gzip_writer_close(cct_context->fastq_unassigned_writer+umfpi);
	}

	for(x1=0; x1<3; x1++) HashTableDestroy(cct_context -> junction_ExonEdgeTree_table[x1]);
	HashTableDestroy(cct_context->junction_transcript_table);
	HashTableDestroy(cct_context->junction_genebody_table);
	HashTableDestroy(cct_context->transcript_to_gene_name_table);

	geinput_close(&cct_context -> input_dataset);
	destroy_offsets(&cct_context->chromosome_table);

	for(sample_i = 1; sample_i <=cct_context-> sample_sheet_table -> numOfElements ; sample_i ++){
		HashTableDestroy(cct_context->chroEvent_entry_table[sample_i]);
		HashTableDestroy(cct_context->chroEvent_detail_table[sample_i]);
		if(cct_context -> junction_to_cell_umi_table[sample_i]) HashTableDestroy(cct_context -> junction_to_cell_umi_table[sample_i]);
	}
	HashTableDestroy(cct_context->sample_sheet_table);

	HashTableDestroy(cct_context->lineno1B_to_sampleno1B_tab);
	ArrayListDestroy(cct_context->sample_id_to_name);
	ArrayListDestroy(cct_context->sample_barcode_list);
	ArrayListDestroy(cct_context->all_features_array);
	HashTableDestroy(cct_context->gene_name_table);
	HashTableDestroy(cct_context->chromosome_exons_table);
	gvindex_destory(cct_context->value_index);

	// either having the barcode array or having the barcode map from Space Ranger BAM.
	if(cct_context -> VisiumHD_barcode_to_best_mapping) HashTableDestroy(cct_context -> VisiumHD_barcode_to_best_mapping);
	if(cct_context->cell_barcodes_array){
		ArrayListDestroy(cct_context->cell_barcodes_array);
		HashTableDestroy(cct_context->cell_barcode_head_tail_table);
	}

	if(cct_context -> cluster_spec_junction_table){
		HashTableDestroy(cct_context -> cluster_spec_junction_table);
		HashTableDestroy(cct_context -> cluster_cell_map_table);
	}

	free(cct_context -> cmd_rebuilt);
	free(cct_context -> value_index);
	free(cct_context -> exonic_region_bitmap);
	free(cct_context -> features_sorted_chr);
	free(cct_context -> features_sorted_geneid);
	free(cct_context -> features_sorted_start);
	free(cct_context -> features_sorted_stop);
	free(cct_context -> features_sorted_strand);
	free(cct_context -> block_end_index);
	free(cct_context -> block_min_start);
	free(cct_context -> block_max_end);
	free(cct_context -> gene_name_array);
	free(cct_context -> read_assignment_counter_locks);
	free(cct_context -> unistr_buffer_space);
	if(cct_context -> read_assignment_detail_file[0])fclose(cct_context -> read_assignment_detail_fp);

	print_in_box(80,0,0,"");
	print_in_box(80,2,0,"");
	SUBREADputs("");
	return 0;
}

void cellCounts_go_chunk_nextchunk(cellcounts_global_t * cct_context){
	cct_context -> running_processed_reads_in_chunk=0;
}

void cellCounts_clean_context_after_chunk(cellcounts_global_t * cct_context) {
	cct_context -> running_processed_reads_in_chunk = 0;
	cct_context -> processed_reads_in_chunk = 0;
}

void cellCounts_init_topKbuff(cellcounts_global_t * cct_context, int thread_no){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	topK_buffer_t * topKbuff = &thread_context -> topKbuff;
	topKbuff -> vote_simple_1_buffer = malloc(cct_context -> max_candidate_voteIJ_per_read * sizeof(simple_mapping_t));
}

void cellCounts_free_topKbuff(cellcounts_global_t * cct_context, int thread_no){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	topK_buffer_t * topKbuff = &thread_context -> topKbuff;
	free(topKbuff -> vote_simple_1_buffer);
}

typedef struct{
	int thread_id;
	int block_start;
	int block_end;
	HashTable * result_tab;
	int * small_side_ordered_event_ids, * large_side_ordered_event_ids;
	chromosome_event_t * event_space;
	cellcounts_global_t * cct_context;
} AT_context_t;


#define SCORING_MAX_QUALITY_MAPPING  10000000ll

srInt_64 cellCounts_calculate_pos_weight_1sec(cellcounts_global_t * cct_context, unsigned int pos, int len){
	unsigned int pos_i;
	srInt_64 ret = 10ll;
	for(pos_i = pos+1; pos_i <= pos+len ; pos_i ++){ //because exons are 1-based. 
		int exonic_map_byte= pos_i /8;
		int exonic_map_bit = pos_i %8;
		if (cct_context ->exonic_region_bitmap [exonic_map_byte] & (1<<exonic_map_bit))
			return SCORING_MAX_QUALITY_MAPPING;
		exonic_map_byte +=  (4096 / 8)*1024*1024;
		if (cct_context ->exonic_region_bitmap [exonic_map_byte] & (1<<exonic_map_bit))
			ret = 13ll;
	}
	return ret;
}

srInt_64 cellCounts_calculate_pos_weight(cellcounts_global_t * cct_context, unsigned int pos, char * cigar){ 
	int tmpi=0, nch;
	srInt_64 max_weight = 10;
	while(0!=(nch = *(cigar++))){
		if(isdigit(nch)){
			tmpi = tmpi*10+nch-'0';
		}else{
			int toadd_chro = 0;
			if(nch=='M'){
				toadd_chro = tmpi;
				max_weight = max(max_weight, cellCounts_calculate_pos_weight_1sec(cct_context, pos, toadd_chro));
				if(max_weight >= SCORING_MAX_QUALITY_MAPPING)return max_weight;
			} else if(nch == 'D' || nch == 'N' || nch == 'S'){ // At this step, the head-soft-clipped bases were not offset to the mapping position. 
				toadd_chro = tmpi;
			}
			pos += toadd_chro;
			tmpi = 0;
		}
	}

	return max_weight;
}

#define MAX_HIT_NUMBER (1000*1000)

void cellCounts_find_hits_for_mapped_section(cellcounts_global_t * cct_context, int thread_no, char * chro_name, int section_begin_pos, int section_end_pos, int is_fragment_negative_strand, int * nhits, char * read_name){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	int start_reverse_table_index = section_begin_pos / REVERSE_TABLE_BUCKET_LENGTH;
	int end_reverse_table_index = (1+section_end_pos) / REVERSE_TABLE_BUCKET_LENGTH;

	//SUBREADprintf("CHRO:%p(%s)\n", chro_name,chro_name);
	fc_chromosome_index_info * this_chro_info = HashTableGet(cct_context -> chromosome_exons_table, chro_name);
	if(this_chro_info == NULL) {
		if(cct_context -> sam_chro_to_anno_chr_alias) {
			char * anno_chro_name = HashTableGet(cct_context -> sam_chro_to_anno_chr_alias, chro_name);
			if(anno_chro_name) this_chro_info = HashTableGet(cct_context -> chromosome_exons_table, anno_chro_name);
		}
		if(this_chro_info == NULL && memcmp(chro_name, "chr", 3)==0) {
			this_chro_info = HashTableGet(cct_context -> chromosome_exons_table, chro_name+3);
		}
		if(this_chro_info == NULL && strlen(chro_name)<=2) {
			char chro_name_buff[MAX_CHROMOSOME_NAME_LEN+5];
			strcpy(chro_name_buff, "chr");
			strcpy(chro_name_buff+3, chro_name);
			this_chro_info = HashTableGet(cct_context -> chromosome_exons_table, chro_name_buff);
		}
	}

//	if( strcmp("chrX",chro_name)==0 && section_begin_pos >= 95077507-60 && section_begin_pos < 95077507+10 )SUBREADprintf("CHROINF %s=%p ; PossibleLen = %d\n", chro_name, this_chro_info, this_chro_info?this_chro_info-> chro_possible_length:-1);
	if(this_chro_info) {
		unsigned int search_start, search_end, search_block_id;
		assert(this_chro_info -> reverse_table_start_index);
		start_reverse_table_index = min(start_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH);
		end_reverse_table_index = min(end_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH+ 1);

		while(start_reverse_table_index<=end_reverse_table_index) {
			search_start = *(this_chro_info -> reverse_table_start_index +start_reverse_table_index);
			if(search_start<0xffffff00)break;
			start_reverse_table_index++;
		}
		//SUBREADprintf("REV_TAB RANGE: %d ~ %d for CHRO %s => SEARCH_START %d ; SEARCH_BLOCK_ID %d\n", start_reverse_table_index, end_reverse_table_index, chro_name, search_start, search_block_id);
		if(search_start>0xffffff00) return;
		search_end = this_chro_info -> chro_block_table_end;

		for(search_block_id=search_start; search_block_id<search_end; search_block_id++){
			//SUBREADprintf("TEST BLOCK %s : %d --  %d ~ %d ?? %d ~ %d\n", chro_name, search_block_id, cct_context -> block_min_start[search_block_id] , cct_context -> block_max_end[search_block_id], section_begin_pos, section_end_pos);
			if (cct_context -> block_min_start[search_block_id] > section_end_pos) break;
			if (cct_context -> block_max_end[search_block_id] < section_begin_pos) continue;

			int search_item_start = 0, search_item_end = cct_context -> block_end_index[search_block_id];
			if(search_block_id>0)search_item_start = cct_context -> block_end_index[search_block_id-1];
			//SUBREADprintf("REACHED SEARCH %s : %d : BLOCK=%d. ITEM START=%d ~ %d\n", chro_name, section_begin_pos, search_block_id, search_item_start, search_item_end);

			// search_item_id is the inner number of the exons.
			int search_item_id;
			for(search_item_id = search_item_start ; search_item_id < search_item_end; search_item_id++) {
				if (cct_context -> features_sorted_stop[search_item_id] >= section_begin_pos) {
					if (cct_context -> features_sorted_start[search_item_id] > section_end_pos) break;
					// there is an overlap >=1 between read and feature.
					// the overlap length is min(end_r, end_F) - max(start_r, start_F) + 1
					
					int is_strand_ok = 1;
					if(cct_context -> need_check_strand){
						is_strand_ok = (is_fragment_negative_strand == cct_context -> features_sorted_strand[search_item_id]);
						if(cct_context -> need_check_strand==-1) is_strand_ok =!is_strand_ok;
					}
//					if( strcmp("chrX",chro_name)==0 && section_begin_pos >= 95077507-60 && section_begin_pos < 95077507+10 )SUBREADprintf("Geneno=%d, StrandOK=%d\n", search_item_id, is_strand_ok); 

					if(is_strand_ok){
						if((*nhits) >= thread_context -> hits_number_capacity - 1) {
							thread_context -> hits_number_capacity = max(10, thread_context -> hits_number_capacity*2);
							thread_context -> hits_start_pos = realloc(thread_context -> hits_start_pos , sizeof(int) * thread_context -> hits_number_capacity);
							thread_context -> hits_length = realloc(thread_context -> hits_length, sizeof(short) * thread_context -> hits_number_capacity);
							thread_context -> hits_chro = realloc(thread_context -> hits_chro, sizeof(char *) * thread_context -> hits_number_capacity);
							thread_context -> hits_indices = realloc(thread_context -> hits_indices, sizeof(srInt_64) * thread_context -> hits_number_capacity);
						}

						if((*nhits) <= MAX_HIT_NUMBER - 1) {
							thread_context -> hits_indices[(*nhits)] = search_item_id;
							(*nhits)++;
						} else {
							SUBREADprintf("ERROR: the read overlapped with more than %d features.\n", (*nhits));
							return ;
						}
					}
				} 
			}
		}
	}
}

#define SCRNA_READ_NAME_SPLIT_CHAR '|'

int cellCounts_scan_read_name_str(cellcounts_global_t * cct_context, char * rbin, char * read_name, char ** sample_seq, char ** sample_qual, char ** seq_1R, char ** qual_1Y, char ** UMI_seq, char ** UMI_qual, char ** lane_str, char ** RG, int * rname_trimmed_len){
	char * testi;
	int field_i=0;
	if(NULL == read_name && rbin) read_name = rbin + 36;
	for(testi = read_name +1; * testi; testi ++){
		if((*testi)== SCRNA_READ_NAME_SPLIT_CHAR || ((*testi)== ':' && cct_context -> input_mode == GENE_INPUT_BCL )){
			field_i++;
			if(field_i == 1) {
				if(rname_trimmed_len) (*rname_trimmed_len)=testi-read_name;
				if(seq_1R)(*seq_1R) = testi+1;

				if(cct_context -> visium_hd_barcodes){
					if(UMI_seq)(*UMI_seq) = testi+1;// in VisiumHD: the R1 is UMI + cell_barcode1+cell_barcode2.
				}else{
					if(UMI_seq)(*UMI_seq) = testi+1+cct_context -> known_cell_barcode_length;
				}
			}else if(field_i == 2){
				if(qual_1Y)(*qual_1Y) = testi+1;

				if(cct_context -> visium_hd_barcodes){
					if(UMI_qual)(*UMI_qual) = testi+1;// in VisiumHD: the R1 is UMI + cell_barcode1+cell_barcode2.
				}else{
					if(UMI_qual)(*UMI_qual) = testi+1+cct_context -> known_cell_barcode_length;
				}
			}else if(field_i == 3){
				*sample_seq = testi + 1;
				if(RG)(*RG) = *sample_seq;
			}else if(field_i == 4){
				if(sample_qual)(*sample_qual) = testi + 1;
			}else if(field_i == 5){
				(*lane_str) = testi + 1;
				if(memcmp(*lane_str, "@RgLater@", 9)==0) (*lane_str) += 9;
				break;
			}
		}
	}
	if(cct_context -> UMI_length <1){ // no locking is needed because it can be safely done many times.
		int umi_end_pos=0,nch;
		for(umi_end_pos=0; 0!=(nch = (*UMI_seq) [umi_end_pos]); umi_end_pos++) if(!isalpha(nch))break;
		if(umi_end_pos > MAX_UMI_LEN){
			SUBREADprintf("ERROR: the UMI length is abnormally long (%d bases). This can be caused by an incorrect cell barcode file.\n", umi_end_pos);
		  	umi_end_pos = MAX_UMI_LEN;
			cct_context -> has_error = 1;
		}
		cct_context -> UMI_length = umi_end_pos; 
	}

	return field_i;
}


int cellCounts_get_cellbarcode_no(cellcounts_global_t * cct_context, int thread_no, char * seq_1R, char * qual_1Y){
	//return -1;
	char tmpc [MAX_READ_NAME_LEN];
	int xx1, xx2,tb1=-1;
	ArrayList * ret=NULL;

	if(cct_context->visium_hd_barcodes && NULL== cct_context->VisiumHD_barcode_to_best_mapping){
		int seq_1Rlen = strstr(seq_1R,"|")-seq_1R, bc2_end=0, bc1_end=0, bc1_start=0, xx3, bc1=-1, bc2=-1;
		int min_bc1_misma = 2, min_bc2_misma = 2;
		for(xx2 = cct_context -> UMI_length ; xx2 < cct_context -> UMI_length+2; xx2++){ // probe the start of bc1
			for(xx1=1;xx1<3;xx1++){
				tmpc[0] = (xx1==2)?'S':'F';
				for(xx3=0; xx3<MIN_LEN_VISIUM_HD_CELLBC/2 ; xx3++)
					tmpc[1+xx3] = seq_1R[2*xx3+xx2+xx1-1];
				tmpc[1+MIN_LEN_VISIUM_HD_CELLBC/2]=0;
				ArrayList *xrawarr = HashTableGet(cct_context -> cell_barcode_head_tail_table, tmpc);
				if(!xrawarr)continue;

				for(xx3=0;xx3<xrawarr->numOfElements;xx3++){
					int tbcn = ArrayListGet( xrawarr, xx3 )-NULL;
					char * known_cellbc = ArrayListGet(cct_context -> cell_barcodes_array, tbcn);
					int hc = hamming_dist_ATGC_max2( known_cellbc, seq_1R+xx2 );
					if(hc < min_bc1_misma || (hc==min_bc1_misma && bc1_end < xx2+strlen(known_cellbc))){
						bc1_end = xx2+strlen(known_cellbc);
						bc1_start = xx2;
						min_bc1_misma = hc;
						bc1 = tbcn;
					}
				}
			}
		}
		if(!bc1_end) return -1;
		seq_1R[ bc1_start ] +=0x20; // upper => lower
		for(xx2 = bc1_end -1; xx2 <bc1_end+2; xx2++){ // probe the start of bc2 . BC1 and BC2 may share a base!!!
			for(xx1=1;xx1<3;xx1++){
				tmpc[0] = (xx1==2)?'S':'F';
				for(xx3=0; xx3<MIN_LEN_VISIUM_HD_CELLBC/2 ; xx3++)
					tmpc[1+xx3] = seq_1R[2*xx3+xx2+xx1-1];
				tmpc[1+MIN_LEN_VISIUM_HD_CELLBC/2]=0;
				ArrayList *xrawarr = HashTableGet(cct_context -> cell_barcode_head_tail_table, tmpc);
				if(!xrawarr)continue;

				for(xx3=0;xx3<xrawarr->numOfElements;xx3++){
					int tbcn = ArrayListGet( xrawarr, xx3 )-NULL;
					char * known_cellbc = ArrayListGet(cct_context -> cell_barcodes_array, tbcn);
					int hc = hamming_dist_ATGC_max2( known_cellbc, seq_1R+xx2 );
					if(hc < min_bc2_misma || (hc==min_bc2_misma && bc2_end < xx2+strlen(known_cellbc))){
						bc2_end = xx2+strlen(known_cellbc); 
						min_bc2_misma = hc;
						bc2 = tbcn;
					}
				}
			}
		}
		seq_1R[ bc2_end ] +=0x20;
		tb1 = bc1 << 16 | bc2;
	}else if(cct_context->visium_hd_barcodes && cct_context->VisiumHD_barcode_to_best_mapping){
		int cbclen = strstr(seq_1R,"|")-seq_1R;
		char bcback = seq_1R[cbclen], bqback = qual_1Y[cbclen];
		seq_1R[cbclen] = 0;
		qual_1Y[cbclen] = 0;
		// space ranger CB is like s_002um_02768_00939-1
		char CKey[200];
		int bc1=-1, bc2=-1;
		sprintf(CKey,"%s/%s",seq_1R + cct_context -> UMI_length, qual_1Y + cct_context -> UMI_length); // UMI is before the two spot barcodes in Visium HD R1 reads. Hence we don't index the UMIs in the R1.
		char * spaceranger_CB = HashTableGet(cct_context->VisiumHD_barcode_to_best_mapping, CKey);
		if(spaceranger_CB){
			sscanf(spaceranger_CB, "%d_%d_", &bc1, &bc2);
			tb1 = bc1 << 16 | bc2;
//			fprintf(stderr,"CREATE_CELLID %s %s = %s\n",  seq_1R, qual_1Y, spaceranger_CB);
		}
		seq_1R[cbclen]=bcback;
		qual_1Y[cbclen]=bqback;
	}else{
		for(xx1=0;xx1<3;xx1++){
			if(xx1==1) ret = ArrayListCreate(100);
			if(xx1>0){
				tmpc[0] = (xx1==2)?'S':'F';
				for(xx2=0; xx2<cct_context -> known_cell_barcode_length/2 ; xx2++)
					tmpc[1+xx2] = seq_1R[2*xx2+xx1-1];
				tmpc[1+cct_context -> known_cell_barcode_length/2]=0;
			}else{
				memcpy(tmpc, seq_1R, cct_context -> known_cell_barcode_length);
				tmpc[cct_context -> known_cell_barcode_length]=0;
			}

			void *xrawarr = HashTableGet(cct_context -> cell_barcode_head_tail_table, tmpc);

			if(xx1 == 0){
				//if(xrawarr) SUBREADprintf("CAFE ? %p\n", xrawarr);
				srInt_64 xint = xrawarr - NULL;
				if(( xint & 0xFFFFFFFFF0000000llu)== IMPOSSIBLE_MEMORY_SPACE){
					int only_cell_id = xint - IMPOSSIBLE_MEMORY_SPACE;
					// no memory was allocated.
					return only_cell_id;
				}
			}else{
				ArrayList * rawarr = xrawarr;
				if(rawarr){
					int xx3,xx2, found;
					for(xx2=0; xx2<rawarr->numOfElements; xx2++){
						int bcno = ArrayListGet(rawarr, xx2)-NULL;
						found=0;
						for(xx3=0;xx3<ret -> numOfElements;xx3++){
							if(ArrayListGet(ret, xx3)==NULL+bcno){
								found=1;
								break;
							}
						}

						if(!found)ArrayListPush(ret, NULL+bcno);
					}
				}
			}
		}

		for(xx1=0; xx1<ret -> numOfElements; xx1++){
			int tbcn = ArrayListGet(ret,xx1)-NULL;
			char * known_cellbc = ArrayListGet(cct_context -> cell_barcodes_array, tbcn);
			int hc = hamming_dist_ATGC_max2( known_cellbc, seq_1R);

			if(hc==1){
				tb1 = tbcn;
				break;
			}
		}
		ArrayListDestroy(ret);
	}
	//SUBREADprintf("CANDIDATE CELL BARCODES=%ld ; hit = %d\n", ret->numOfElements, tb1);

	return tb1;
}

void cellCounts_build_read_bin(cellcounts_global_t * cct_context, int thread_no, char * rbin, char * read_name, int read_name_len, int trimmed_rname_len, int read_len, char * read_text, char * qual_text, char * chro_name, int chro_pos, int reporting_index, int multi_mapping_number, int this_multi_mapping_i, int editing_dist, int vote_for_aln, srInt_64 ma_misma_ins_Sclip){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	char * cigar = NULL;
	int mapping_quality = 255;

	int flags = 4;
	if(reporting_index >=0){
		flags = thread_context -> reporting_flags[reporting_index];
		cigar = thread_context -> reporting_cigars[reporting_index];
		mapping_quality = thread_context -> reporting_mapq[reporting_index];
		if(this_multi_mapping_i>1) flags += 256;
	}

	int cigar_opts[1+2*JUNCTION_REALIGNMENT_MAX_DEPTH *2], xk1, cover_length = 0;
	int cigar_opt_len = 0;
	if(cigar) cigar_opt_len = SamBam_compress_cigar(cigar, cigar_opts, & cover_length, 1+2*JUNCTION_REALIGNMENT_MAX_DEPTH *2);

	int record_length = 4 + 4 + 4 + 4 + 4 +  /* l_seq: */ 4 + 4 + 4 + 4 + /* read_name:*/ read_name_len + cigar_opt_len * 4 + (read_len + 1) /2 + read_len;

	int bin = SamBam_reg2bin(chro_pos -1, chro_pos-1+cover_length);
	int refID = -1;
	if(chro_name) refID = HashTableGet(cct_context -> chromosome_table.read_name_to_index, chro_name) - NULL - 1;
	int bin_mq_nl = (bin<<16) | (mapping_quality << 8) | read_name_len ;
	int fag_nc = (flags<<16) | cigar_opt_len;
	int nextRefID = -1, temp_len = 0, next_chro_pos=-1;
	chro_pos--;
	memcpy(rbin + 4 , & refID , 4);
	memcpy(rbin + 8 , & chro_pos , 4);
	memcpy(rbin +12 , & bin_mq_nl , 4);
	memcpy(rbin +16 , & fag_nc , 4);
	memcpy(rbin +20 , & read_len , 4);
	memcpy(rbin +24 , & nextRefID , 4);
	memcpy(rbin +28 , & next_chro_pos , 4);
	memcpy(rbin +32 , & temp_len , 4);
	memcpy(rbin +36 , read_name, read_name_len);
	memcpy(rbin +36 + read_name_len, cigar_opts, cigar_opt_len*4);
	SamBam_read2bin(read_text  , rbin +36 +read_name_len +cigar_opt_len*4);

	int basev = 36 +read_name_len +cigar_opt_len*4+(read_len + 1) /2;
	for(xk1=0; xk1<read_len; xk1++) *(rbin +basev+xk1) = qual_text[xk1]-33;
	
	if(multi_mapping_number>=0){
		rbin[record_length ++]='X';
		rbin[record_length ++]='V';
		rbin[record_length ++]='C';
		rbin[record_length ++]=vote_for_aln;

		rbin[record_length ++]='H';
		rbin[record_length ++]='I';
		rbin[record_length ++]='C';
		rbin[record_length ++]=this_multi_mapping_i;

		rbin[record_length ++]='N';
		rbin[record_length ++]='H';
		rbin[record_length ++]='C';
		rbin[record_length ++]=multi_mapping_number;
	}

	if(ma_misma_ins_Sclip){
		int ma = (ma_misma_ins_Sclip>>48)&0xffff;
		int misma = (ma_misma_ins_Sclip>>32)&0xffff;
		int ins = (ma_misma_ins_Sclip>>16)&0xffff;
		int clip = (ma_misma_ins_Sclip)&0xffff;
		rbin[record_length ++]='X';
		rbin[record_length ++]='B';
		rbin[record_length ++]='Z';
		char mamainf[50];
		int mamainf_len = snprintf(mamainf,50, "%d,%d,%d,%d", ma, misma, ins, clip);
		strcpy(rbin+record_length, mamainf);
		record_length += mamainf_len+1;
	}

	if(editing_dist>=0){
		rbin[record_length ++]='N';
		rbin[record_length ++]='M';
		rbin[record_length ++]='C';
		rbin[record_length ++]=editing_dist;
	}

	record_length -=4;
	memcpy(rbin     , & record_length , 4);
}

void cellCounts_write_one_read_bin(cellcounts_global_t * cct_context, int thread_no, REPFILE * binfp, int sample_no, int cellbarcode_no, char * umi_barcode, char * readbin, int nhits, int notmapped){
	int x1;
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;

	REP_fwrite(&sample_no,4,1,binfp);
	if(0==notmapped){
		REP_fwrite(&cellbarcode_no,4,1,binfp);
		if(nhits <1){
			srInt_64 zero_genes = 1LLU<<63;
			REP_fwrite(&zero_genes, 8,1,binfp);
		}else if(nhits<2){
			srInt_64 gene_no = thread_context -> hits_indices[0];
			REP_fwrite(&gene_no, 8,1,binfp);
		}else{
			srInt_64 total_genes = nhits + (1LLU<<63);
			REP_fwrite(&total_genes, 8,1,binfp);
			for(x1=0;x1<nhits;x1++){
				srInt_64 gene_no = thread_context -> hits_indices[x1];
				REP_fwrite(&gene_no, 8,1,binfp);
			}
		}
		REP_fwrite(umi_barcode, cct_context->UMI_length, 1, binfp);
	}
	memcpy(&x1, readbin, 4);
	x1+=4;
	REP_fwrite(readbin, x1, 1, binfp);

	if(cct_context -> read_assignment_detail_fp && nhits>0 && 0==notmapped){
		char * cellbc = NULL;
		if(cellbarcode_no>=0)cellbc = ArrayListGet(cct_context -> cell_barcodes_array, cellbarcode_no);
		if(cellbc){
			cellCounts_lock_occupy(&cct_context -> read_assignment_detail_lock);
			char * rname = readbin + 36;
			char tmpc = umi_barcode[cct_context->UMI_length];
			rname[12]=0;
			umi_barcode[cct_context->UMI_length] = 0;
			fprintf(cct_context -> read_assignment_detail_fp,"READ_TO_GENE\tSAMPLE%03d\t%s\t%s\t%s", sample_no, rname, cellbc, umi_barcode);
			for(x1=0; x1<nhits;x1++){
				srInt_64 entrez_no = thread_context -> hits_indices[x1];
				fprintf(cct_context -> read_assignment_detail_fp,"\t%s", cct_context ->gene_name_array[entrez_no]);
			}
			fprintf(cct_context -> read_assignment_detail_fp,"\n");
			rname[12]='|';
			umi_barcode[cct_context->UMI_length] = tmpc;
			cellCounts_lock_release(&cct_context -> read_assignment_detail_lock);
		}
	}
}

int cellCounts_get_sample_id(cellcounts_global_t * cct_context, char * sbc, int read_laneno){
	int x1;
	for(x1=0; x1 < cct_context -> sample_barcode_list -> numOfElements ; x1++ ){
		char ** lane_and_barcode = ArrayListGet(cct_context -> sample_barcode_list, x1);
		int sheet_lane_no = lane_and_barcode[0]-(char*)NULL;
		if(sheet_lane_no == LANE_FOR_ALL_LANES || read_laneno == sheet_lane_no){
			int sample_no = lane_and_barcode[1]-(char*)NULL;
			char * knownbar = lane_and_barcode[2];
			if(lane_and_barcode[3]){
				int hd = hamming_dist_ATGC_max1_2p( sbc, knownbar );
				if(hd<=2) return sample_no;
			}else{
				int hd = hamming_dist_ATGC_max1( sbc, knownbar );
				if(hd<=1) return sample_no;
			}
		}
	}
	return -1;
}

int cellCounts_parallel_gzip_writer_add_read_fqs_scRNA(parallel_gzip_writer_t**outfps, char * bambin, int thread_no, char * read_text_raw, char * qual_text_raw){ // the text-raw variables are the reads in their input form (not potentially reversed form)
	int reclen=0;
	parallel_gzip_writer_t * outR1fp = outfps[0];
	parallel_gzip_writer_t * outI1fp = outfps[1];
	parallel_gzip_writer_t * outI2fp = outfps[2];
	parallel_gzip_writer_t * outR2fp = outfps[3];

	memcpy(&reclen, bambin,4);
	int flag = 0, l_seq = 0, l_read_name = 0, n_cigar_ops = 0;
	memcpy(&l_read_name, bambin+12,1);
	memcpy(&n_cigar_ops, bambin+16,2);
	memcpy(&flag, bambin+18,2);
	memcpy(&l_seq, bambin+20,4);

	parallel_gzip_writer_add_text(outR2fp,"@",1,thread_no);
	parallel_gzip_writer_add_text(outR1fp,"@",1,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"@",1,thread_no);
	if(outI2fp) parallel_gzip_writer_add_text(outI2fp,"@",1,thread_no);
	char * readname = bambin+36;
	parallel_gzip_writer_add_text(outR1fp,readname, 12,thread_no);
	parallel_gzip_writer_add_text(outR2fp,readname, 12,thread_no);
	parallel_gzip_writer_add_text(outI1fp,readname, 12,thread_no);
	if(outI2fp) parallel_gzip_writer_add_text(outI2fp,readname, 12,thread_no);

	parallel_gzip_writer_add_text(outR1fp,"\n",1,thread_no);
	parallel_gzip_writer_add_text(outR2fp,"\n",1,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"\n",1,thread_no);
	if(outI2fp) parallel_gzip_writer_add_text(outI2fp,"\n",1,thread_no);

	//SUBREADprintf("WRITEFQ RNAME '%s'\n", bambin+36);
	char * R1seq = bambin+36+13;
	int R1len = 0;
	for(R1len=0; R1seq[R1len] && R1seq[R1len]!='|' ;R1len++);
	char * R1qual = R1seq + R1len + 1;
	parallel_gzip_writer_add_text(outR1fp,R1seq, R1len,thread_no);
	parallel_gzip_writer_add_text(outR1fp,"\n+\n",3,thread_no);
	parallel_gzip_writer_add_text_qual(outR1fp,R1qual, R1len,thread_no);
	parallel_gzip_writer_add_text(outR1fp,"\n",1,thread_no);

	char * I1seq = R1qual + R1len + 1;
	int I1I2len = 0;
	for(I1I2len=0; I1seq[I1I2len] && I1seq[I1I2len]!='|' ;I1I2len++);
	int I1len = I1I2len;
	if(outI2fp) I1len /=2;

	char * I2seq = I1seq + I1len;
	char * I1qual = I1seq + I1I2len + 1;
	char * I2qual = I1seq + I1I2len + I1len + 1;

	parallel_gzip_writer_add_text(outI1fp,I1seq, I1len,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"\n+\n",3,thread_no);
	parallel_gzip_writer_add_text_qual(outI1fp,I1qual, I1len,thread_no);
	parallel_gzip_writer_add_text(outI1fp,"\n",1,thread_no);

	if(outI2fp){
		parallel_gzip_writer_add_text(outI2fp,I2seq, I1len,thread_no);
		parallel_gzip_writer_add_text(outI2fp,"\n+\n",3,thread_no);
		parallel_gzip_writer_add_text_qual(outI2fp,I2qual, I1len,thread_no);
		parallel_gzip_writer_add_text(outI2fp,"\n",1,thread_no);
	}

	char * oseq = read_text_raw;
	parallel_gzip_writer_add_text(outR2fp, oseq, l_seq,thread_no);
	parallel_gzip_writer_add_text(outR2fp,"\n+\n",3,thread_no);

	oseq = qual_text_raw;
	parallel_gzip_writer_add_text(outR2fp, oseq, l_seq,thread_no);
	parallel_gzip_writer_add_text(outR2fp,"\n",1,thread_no);
	return 0;
}

void cellCounts_vote_and_add_count(cellcounts_global_t * cct_context, int thread_no, int sample_no, char * read_name, int rlen, char * read_text, char * qual_text, char * raw_text, char * raw_qual, char * chro_name, int chro_pos, int reporting_index, int nhits, int multi_mapping_number, int this_multi_mapping_i, int editing_dist, int vote_for_aln, srInt_64 ma_misma_ins_Sclip, char * BC_seq, char * UMI_seq, int rname_trimmed_len, int cell_barcode_no){
	int batch_no;
	
	if(cct_context->visium_hd_barcodes) UMI_seq = BC_seq;
	if(nhits > 1 && !cct_context -> allow_multi_overlapping_reads) nhits = 0;
	if(reporting_index >=0){
		if(cell_barcode_no>=0 && sample_no>0) batch_no = cell_barcode_no % CELLBC_BATCH_NUMBER;
		else if(sample_no>0) batch_no = CELLBC_BATCH_NUMBER;
	}else batch_no = CELLBC_BATCH_NUMBER+1;
	
	char readbin[READ_BIN_BUF_SIZE];
	cellCounts_build_read_bin(cct_context, thread_no, readbin, read_name, strlen(read_name), rname_trimmed_len, rlen, read_text, qual_text, chro_name, chro_pos, reporting_index, multi_mapping_number, this_multi_mapping_i, editing_dist, vote_for_aln, ma_misma_ins_Sclip);

	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	if(sample_no>0){
		cellCounts_lock_occupy(cct_context -> batch_file_locks + batch_no);
		REPFILE * binfp = cct_context -> batch_files [ batch_no ];

		cellCounts_write_one_read_bin(cct_context, thread_no, binfp, sample_no, cell_barcode_no, UMI_seq, readbin, nhits, batch_no == CELLBC_BATCH_NUMBER+1);
		cellCounts_lock_release(cct_context -> batch_file_locks + batch_no);

		if(reporting_index >=0 && this_multi_mapping_i==1){
			thread_context -> mapped_reads_per_sample[ sample_no -1 ]++;
			thread_context -> reads_per_sample[ sample_no -1 ]++;
			if(nhits >0) thread_context -> assigned_reads_per_sample[ sample_no -1 ]++;
		}else thread_context -> reads_per_sample[ sample_no -1 ]++;
	}else thread_context -> reads_per_sample[ cct_context-> sample_sheet_table -> numOfElements ]++; // unassigned reads in NO+1:

	void * pps[6];
	void ** sample_bam_2fps = (void**)pps;
	if(sample_no<=0){
		pps[0]=NULL;
		pps[1]=cct_context -> fastq_unassigned_writer;
		pps[2]=cct_context -> fastq_unassigned_writer+1;
		if(cct_context -> is_dual_index)pps[3]=cct_context -> fastq_unassigned_writer+2;
		else pps[3]=NULL;
		pps[4]=cct_context -> fastq_unassigned_writer+3;
		pps[5]=&cct_context -> fastq_unassigned_lock;
	}else sample_bam_2fps = HashTableGet(cct_context -> sample_BAM_writers, NULL+(sample_no-1) + 1); // sample_id-1: 0,1,2,...

	if(GENE_INPUT_SCRNA_FASTQ != cct_context -> input_mode){
		parallel_gzip_writer_t **gz3fps = (parallel_gzip_writer_t **)sample_bam_2fps+1;
		cellCounts_parallel_gzip_writer_add_read_fqs_scRNA(gz3fps, readbin, thread_no, raw_text, raw_qual);
		// I2 always has the same length as I1, hence no test is required.
		if( gz3fps[0]-> thread_objs[thread_no].in_buffer_used >= PARALLEL_GZIP_TXT_BUFFER_SIZE - PARALLEL_GZIP_TXT_BUFFER_MARGIN ||
		    gz3fps[1]-> thread_objs[thread_no].in_buffer_used >= PARALLEL_GZIP_TXT_BUFFER_SIZE - PARALLEL_GZIP_TXT_BUFFER_MARGIN ||
		    gz3fps[3]-> thread_objs[thread_no].in_buffer_used >= PARALLEL_GZIP_TXT_BUFFER_SIZE - PARALLEL_GZIP_TXT_BUFFER_MARGIN ){
			parallel_gzip_zip_texts(gz3fps[0], thread_no, 0);
			parallel_gzip_zip_texts(gz3fps[1], thread_no, 0);
			if(gz3fps[2]) parallel_gzip_zip_texts(gz3fps[2], thread_no, 0);
			parallel_gzip_zip_texts(gz3fps[3], thread_no, 0);
			cellCounts_lock_occupy(sample_bam_2fps[5]);
			parallel_gzip_writer_flush(gz3fps[0], thread_no);
			parallel_gzip_writer_flush(gz3fps[1], thread_no);
			if(gz3fps[2]) parallel_gzip_writer_flush(gz3fps[2], thread_no);
			parallel_gzip_writer_flush(gz3fps[3], thread_no);
			cellCounts_lock_release(sample_bam_2fps[5]);
		}
	}
}

void cellCounts_summarize_entrez_hits(cellcounts_global_t * cct_context, int thread_no, int * nhits){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	if((*nhits) == 0){
		return;
	}else if((*nhits) == 1){
		thread_context -> hits_indices [0] = cct_context -> features_sorted_geneid[thread_context -> hits_indices [0]];
		return;
	}else{
		int scaned_idx, wrt_ptr=0;
		for(scaned_idx = 0; scaned_idx < *nhits; scaned_idx++)
			thread_context -> hits_indices [scaned_idx] = cct_context -> features_sorted_geneid[thread_context -> hits_indices [scaned_idx]];

		for(scaned_idx = 0; scaned_idx < *nhits; scaned_idx++){
			srInt_64 geneid = thread_context -> hits_indices [scaned_idx];
			int k2, found=0;
			for(k2=0; k2<wrt_ptr; k2++) if(thread_context -> hits_indices [k2]==geneid) found=1;
			if(!found) thread_context -> hits_indices [wrt_ptr++] = geneid;
		}
		(*nhits) = wrt_ptr;
	}
}

int cellCounts_is_homopolymer_or_N(const char *s, int ulen) {
	if (!s || !*s) return 0;
	char first = s[0];
	if(first == 'N')return 1;
	int ret=1;
	for (int i = 1; i<ulen; i++) {
		if (s[i] != first) ret = 0;
		if (s[i] == 'N') return 1;
	}
	return ret;
}

void cellCounts_write_read_in_batch_bin(cellcounts_global_t * cct_context, int thread_no, int sample_i, int reporting_index, char * read_name, char * read_text, char * qual_text, char * raw_text, char * raw_qual, int rlen){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	char * chro_name = NULL;
	int chro_pos = 0;
	unsigned int linear_pos;

	int rname_trimmed_len=0;
	char * sample_seq=NULL, *sample_qual=NULL, *BC_qual=NULL, *BC_seq=NULL, *UMI_seq=NULL, *UMI_qual=NULL, *lane_str=NULL, *RG=NULL, *testi;
	cellCounts_scan_read_name_str(cct_context, NULL, read_name, &sample_seq, &sample_qual, &BC_seq, &BC_qual, &UMI_seq, &UMI_qual, &lane_str, &RG, &rname_trimmed_len);

	int cell_barcode_no = cellCounts_get_cellbarcode_no(cct_context, thread_no, BC_seq, BC_qual);
	//if(cell_barcode_no>=0&&cct_context->visium_hd_barcodes) fprintf(stderr,"    UMIseq=%s  UMIqual=%s  UMIlen=%d\n", UMI_seq, UMI_qual, cct_context->UMI_length);
	if(reporting_index>=0){
		linear_pos = thread_context -> reporting_positions[reporting_index];
		linear_pos += get_soft_clipping_length(thread_context -> reporting_cigars[reporting_index]);

		locate_gene_position(linear_pos+1, &cct_context -> chromosome_table, &chro_name, &chro_pos );
		if(!chro_name){
			reporting_index=-1;
			read_text = raw_text;
			qual_text = raw_qual;
		}

		int is_negative = (thread_context -> reporting_flags[reporting_index] & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0;

		int nch, nchi, add_curs, tmpi=0, chro_curs=chro_pos, nhits=0;
		for(nchi=0; 0!=(nch= thread_context -> reporting_cigars[reporting_index][nchi]); nchi++){
			if(isdigit(nch)){
				tmpi = 10*tmpi + (nch-'0');
			}else{
				add_curs = 0;
				if(nch=='M'||nch=='D'||nch=='N') //'S'  has been added to loc by get_soft_clipping_length above.
					add_curs =tmpi;

				if(nch=='M')
					cellCounts_find_hits_for_mapped_section(cct_context, thread_no, chro_name, chro_curs, chro_curs+add_curs -1, is_negative,&nhits, read_name); // "section end pos" is inclusive

				chro_curs += add_curs;
				tmpi=0;
			}
		}

		cellCounts_summarize_entrez_hits(cct_context, thread_no, &nhits);

		cellCounts_vote_and_add_count(cct_context, thread_no, sample_i, read_name, rlen, read_text, qual_text, raw_text, raw_qual, chro_name, chro_pos, reporting_index, nhits, thread_context -> total_voteIJs_to_write, thread_context -> writing_voteID_buf_index +1, thread_context -> reporting_editing_distance[reporting_index], thread_context -> reporting_vote_for_aln[reporting_index],thread_context -> reporting_ma_misma_ins_Sclip[reporting_index], BC_seq, UMI_seq, rname_trimmed_len, cell_barcode_no);
	}else //unmapped
		cellCounts_vote_and_add_count(cct_context, thread_no, sample_i, read_name, rlen, read_text, qual_text, raw_text, raw_qual, NULL, 0, -1, 0, 0, 0, -1, 0, 0, BC_seq, UMI_seq, rname_trimmed_len, cell_barcode_no);
}

int cellCounts_fetch_next_read_pair(cellcounts_global_t * cct_context, int thread_no,int *read_len, char * read_name, char * read_text, char * qual_text, subread_read_number_t * read_no_in_chunk) ;
int cellCounts_do_voting(cellcounts_global_t * cct_context, int thread_no) ;
int cellCounts_do_junctable(cellcounts_global_t * cct_context, int thread_no) ;
int cellCounts_do_realign(cellcounts_global_t * cct_context);
static int cellCounts_temp_realign_fp_finish_write(cellcounts_temp_file_point_t * temp_fp);

void * cellCounts_run_in_thread(void * params){
	void ** parameters = (void **)params;
	cellcounts_global_t * cct_context = (cellcounts_global_t * ) parameters[0];
	int thread_no = parameters[1]-NULL;
	int task = parameters[2]-NULL;
	int * ret_value_pointer = (int *)parameters[3];
	free(parameters);

	switch(task) {
		case STEP_JUNC_TABLE:
			*ret_value_pointer = cellCounts_do_junctable(cct_context, thread_no);
		break;

		case STEP_VOTING:
			*ret_value_pointer = cellCounts_do_voting(cct_context, thread_no);
		break;
	}
	return NULL;
}

int cellCounts_release_context_from_align(cellcounts_global_t * cct_context, int thread_no, int task) {
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	destroy_typical_dynamic_align((void***)thread_context -> dynamic_align_buffers, MAX_SCRNA_READ_LENGTH);
	if(thread_context -> temp_realign_record_buf)
		free(thread_context -> temp_realign_record_buf);
	if(thread_context -> temp_realign_work_buf)
		free(thread_context -> temp_realign_work_buf);
	thread_context -> temp_realign_work_buf = thread_context -> temp_realign_record_buf = NULL;

	if(cct_context -> do_cell_level_junction_detection){
		cellCounts_temp_realign_fp_finish_write(&thread_context -> realign_temp_fp);
		cellcounts_temp_file_fclose(&thread_context -> realign_temp_fp);

		if(cct_context -> cell_level_junction_memory_temp){
			cct_context -> all_thread_realign_fp_ptrs[thread_no] = thread_context -> realign_temp_fp.realign_temp_memspace;
			cct_context -> all_thread_realign_fp_ints[thread_no*2] = thread_context -> realign_temp_fp.realign_temp_usedmem;
			cct_context -> all_thread_realign_fp_ints[thread_no*2+1] = thread_context -> realign_temp_fp.realign_temp_capamem;
		}
	}

	return 0;
}

int cellCounts_prepare_context_for_align(cellcounts_global_t * cct_context, int thread_no, int task) {
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	init_typical_dynamic_align((void***)thread_context -> dynamic_align_buffers, thread_context -> dynamic_align_penalties, MAX_SCRNA_READ_LENGTH);

	return 0;
}

int cellCounts_sort_junction_entry_table(cellcounts_global_t * cct_context);


void cellCounts_join_thread_junc_cell_umi_table(void * key, void * val, HashTable * thtab){
	HashTable * cctab = thtab -> appendix1;
	ArrayList * thlist = val;
	ArrayList * cclist = HashTableGet(cctab, key);
	
	if(!cclist){
		cclist = ArrayListCreate(max(10,thlist->numOfElements));
		HashTablePut(cctab, key, cclist);
	}
	ArrayListExtend(cclist, thlist);
}


int cellCounts_run_maybe_threads(cellcounts_global_t * cct_context, int task){
	int ret_value =0;

	int current_thread_no ;
	cellcounts_align_thread_t * thread_contexts = calloc(sizeof(cellcounts_align_thread_t) , cct_context->total_threads);
	cct_context -> all_thread_contexts = thread_contexts;

	int ret_values[64];

	for(current_thread_no = 0 ; current_thread_no < cct_context->total_threads ; current_thread_no ++) {
		thread_contexts[current_thread_no].thread_no = current_thread_no;
		cellCounts_prepare_context_for_align(cct_context, current_thread_no, task);

		if(STEP_VOTING == task || STEP_JUNC_TABLE == task) cellCounts_init_topKbuff(cct_context, current_thread_no);
		void ** thr_parameters = malloc(sizeof(void*)*4);
		thr_parameters[0] = cct_context;
		thr_parameters[1] = NULL+current_thread_no;
		thr_parameters[2] = NULL+task;
		thr_parameters[3] = ret_values + current_thread_no;

		pthread_create(&thread_contexts[current_thread_no].thread, NULL, cellCounts_run_in_thread, thr_parameters);
	}

	for(current_thread_no = 0 ; current_thread_no < cct_context->total_threads ; current_thread_no ++) {
		pthread_join(thread_contexts[current_thread_no].thread, NULL);
		cct_context -> loconf_map  +=thread_contexts[current_thread_no].loconf_map;
		cct_context -> hiconf_map  +=thread_contexts[current_thread_no].hiconf_map;

		if(STEP_VOTING == task || STEP_JUNC_TABLE == task) cellCounts_free_topKbuff(cct_context, current_thread_no);
		ret_value += *(ret_values + current_thread_no);
		int smpno;
		if(STEP_VOTING == task){
			for(smpno = 0; smpno < cct_context-> sample_sheet_table -> numOfElements; smpno ++){ // mapped / assigned / all read counts use 0-based sample ids.
				cct_context -> mapped_reads_per_sample[smpno] += thread_contexts[current_thread_no].mapped_reads_per_sample[smpno];
				cct_context -> assigned_reads_per_sample[smpno] += thread_contexts[current_thread_no].assigned_reads_per_sample[smpno];
				cct_context -> reads_per_sample[smpno] += thread_contexts[current_thread_no].reads_per_sample[smpno];
//fprintf(stderr,"ADD_COUNT_READS : THR %d => Sample %d ; READS %d %d\n", current_thread_no, smpno, );
			}
			cct_context -> reads_per_sample[smpno] += thread_contexts[current_thread_no].reads_per_sample[smpno]; //  for non-assigned
		}
		cellCounts_release_context_from_align(cct_context, current_thread_no, task);
		if(ret_value)break;
	}

	free(thread_contexts);
	return ret_value;
}

typedef struct {
	unsigned int thread_bodytable_number;
	short thread_no;
} concatinating_events_record_t;

static int is_dual_index, idx_offset, base_offset, sread_len, total_bin_len, rname_tail_pos;
int cellCounts_copy_bin_to_textread(cellcounts_global_t * cct_context, int readlane, unsigned char* readbin, char * read_name, char * seq, char *qual, int * psread_lens, subread_read_number_t rno){
	int  bii;
	if(sread_len < 1){
		is_dual_index = psread_lens[3]>0;
		idx_offset  = psread_lens[0];
		total_bin_len = psread_lens[0]+psread_lens[1]+psread_lens[2]+psread_lens[3];
		base_offset = psread_lens[1] + ( is_dual_index?psread_lens[2]:0) + idx_offset;
		rname_tail_pos = 16 +2*base_offset;
		sread_len = psread_lens[2+is_dual_index];
	}

	//SUBREADprintf("RLENs=%d, idx=%d, base=%d\n", total_bin_len, idx_offset, base_offset);
	#ifdef __MINGW32__
	SUBreadSprintf(read_name, 15, "R%011" PRIu64 "|", rno);
	#else
	SUBreadSprintf(read_name, 15, "R%011llu|", rno);
	#endif

	read_name[13+idx_offset]='|';
	read_name[14+2*idx_offset]='|';
	read_name[15+base_offset+idx_offset]='|';
	SUBreadSprintf(read_name + rname_tail_pos, MAX_READ_NAME_LEN +1- rname_tail_pos , "|@RgLater@L%03d" , readlane);

	for(bii = 0; bii < total_bin_len; bii++){
		unsigned int nch = readbin[bii];
		char nbase, nqual;
		if(nch > 0){
			int nch1 = nch & 0x3;
			nbase = (1413956417) >> (nch1*8); // 1413956417 = 'TGCA'
			nqual=33+(nch >>2);
			if(nqual <='!'){
				nbase='N';
				nqual='#';
			}
		}else{
			nqual='#';
			nbase='N';
		}
		if(nqual >= '/' && bii < base_offset ) nqual++;
		if(bii < idx_offset){
			read_name[13+bii] = nbase;
			read_name[14+idx_offset+bii]= nqual;
		}else if(bii < base_offset ){
			read_name[15+idx_offset+bii] = nbase;
			read_name[16+base_offset+bii]= nqual;
		}else{
			seq[bii - base_offset] = nbase;
			qual[bii - base_offset] = nqual;
		}
	}
	qual[sread_len] =0;
	seq[sread_len] =0;
	return sread_len;
}

int cellCounts_fetch_next_read_pair(cellcounts_global_t * cct_context, int thread_no,int *read_len, char * read_name, char * read_text, char * qual_text, subread_read_number_t * read_no_in_chunk) {
	int rl1=0;
	subread_read_number_t this_number = -1;
	gene_input_t * ginp1 = &cct_context -> input_dataset;

	if( ginp1 -> file_type == GENE_INPUT_BCL ){
		int * read_lengths = ginp1 -> bcl_input.single_read_lengths;
		cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
		if(thread_context -> bcl_input_local_cached<1){
			cellCounts_lock_occupy(&cct_context -> input_dataset_lock); 
			int new_reads = cacheBCL_next_readbin(&ginp1 -> bcl_input, thread_context -> bcl_input_local_readlane , thread_context -> bcl_input_local_readbin, BCL_READBIN_ITEMS_LOCAL, &thread_context -> bcl_input_local_start_no);
			if(new_reads)
				thread_context -> bcl_input_local_cached = thread_context -> bcl_input_local_filled = new_reads;
			else if(!cct_context -> running_processed_reads_in_chunk)
				cct_context -> running_processed_reads_in_chunk = cacheBCL_get_readno_in_dataset(&ginp1 -> bcl_input);
			cellCounts_lock_release(&cct_context -> input_dataset_lock); 
		}
		if(thread_context -> bcl_input_local_cached >0 ){ // bcl_input_local_cached changed above.
			int posnumb = thread_context -> bcl_input_local_filled - thread_context -> bcl_input_local_cached;
			this_number = thread_context -> bcl_input_local_start_no + posnumb;
			thread_context -> bcl_input_local_cached --;
			rl1 = cellCounts_copy_bin_to_textread(cct_context, thread_context -> bcl_input_local_readlane [posnumb], (unsigned char*)thread_context -> bcl_input_local_readbin[posnumb],
				read_name, read_text , qual_text, read_lengths, this_number);
		}
		else rl1=0;

	}else{
		cellCounts_lock_occupy(&cct_context -> input_dataset_lock); 
		if(cct_context -> running_processed_reads_in_chunk < cct_context -> reads_per_chunk) {
			rl1 = geinput_next_read_with_lock(ginp1, read_name, read_text , qual_text);

			if(rl1 > 0){
				this_number = cct_context -> running_processed_reads_in_chunk;
				cct_context -> running_processed_reads_in_chunk ++;
			}
		}
		cellCounts_lock_release(&cct_context -> input_dataset_lock); 
	}

	if(rl1>0 && this_number>=0 && this_number <  1000llu*1000 *1000*1000) {
		*read_no_in_chunk = this_number;
		*read_len = rl1;
		read_text[rl1] = qual_text[rl1] = 0;
		return 0;
	} else {
		*read_no_in_chunk = -1;
		*read_len = -1;
		if(rl1 == -2) cct_context -> has_error=1;
		return 1;
	}
}

void cellCounts_update_top_three(cellcounts_global_t * cct_context, int * top_buffer_3i, int new_value){
	if(new_value > top_buffer_3i[cct_context -> max_distinct_top_vote_numbers - 1]){
		int x1;
		for(x1 = 0;x1 < cct_context -> max_distinct_top_vote_numbers ; x1++){
			if(new_value > top_buffer_3i[x1]){
				int x2;
				for(x2 = cct_context -> max_distinct_top_vote_numbers - 1 ; x2 > x1 ; x2 --){
					top_buffer_3i[x2] = top_buffer_3i[x2-1];
				}
				top_buffer_3i[x1] = new_value;
				break;
			}else if(new_value == top_buffer_3i[x1]) break;
		}
	}
}

void cellCounts_set_insertion_sequence(cellcounts_global_t * cct_context, int thread_no , char ** binary_bases , char * read_text , int insertions) {
	int xk1;

	(*binary_bases) = malloc((1+insertions)/4+2);
	//SUBREADprintf("ALLOC PTR=%p\n", (*binary_bases) );

	assert(insertions <= MAX_INSERTION_LENGTH);
	memset((*binary_bases),0, (1+insertions)/4+2);

	for(xk1=0; xk1<insertions;xk1++)
	{
		int byte_no = xk1/4;
		int bit_no = 2*(xk1%4);

		*((*binary_bases)+byte_no) |= (base2int(read_text[xk1]))<<bit_no;
	}
}
#define _test_record_size	if(current_record_number >= current_record_size - 2){\
		current_record_size *= 1.5;\
		records=realloc(records, sizeof(scanning_events_record_t)*current_record_size);\
		if(NULL == records) return -1;\
	}\

#define _add_record	records[current_record_number].scanning_positons = body -> event_small_side;\
	records[current_record_number].thread_bodytable_number = xx1;\
	current_record_number++;\
	records[current_record_number].scanning_positons = body -> event_large_side;\
	records[current_record_number].thread_bodytable_number = xx1;\
	current_record_number++;\

typedef struct {
	unsigned int scanning_positons;
	unsigned int thread_bodytable_number;
} scanning_events_record_t;

int cellCounts_indel_recorder_copy(cellcounts_vote_number_t * alnrec, cellcounts_vote_number_t * votrec, int indelrec_num, int applied_subreads_per_strand, int * first_base_offset_from_mapped_loc, int * span_chro, int * span_read, char * read_name, unsigned int absloc){
	if(indelrec_num<=3){
		alnrec[0]=votrec[0];
		alnrec[1]=votrec[1];
		alnrec[2]=0;
		alnrec[3]=0;
		return 0;
	}

	char offsets[applied_subreads_per_strand];
	int sri,toli;
	memset(offsets,0x77, sizeof(char)*applied_subreads_per_strand);
	for(toli=0; toli<indelrec_num; toli+=3)
		for(sri=votrec[toli]; sri<=votrec[toli+1]; sri++) offsets[sri -1] = votrec[toli+2]; // votrec[toli+0] and votrec[toli+1] are the start/end subread numbers + 1

	for(sri=0; sri < applied_subreads_per_strand ; sri++){
		if(offsets[sri]!=0x77){
			*first_base_offset_from_mapped_loc=offsets[sri];
			break;
		}
	}

	int offset_cur = 0, subread0 =-1, high_conf_index = 0;
	toli = 0;
	alnrec[3]=0;
	for(sri=0; ; sri++){
		if(offsets[sri]!=0x77){
			offsets[sri] -= (*first_base_offset_from_mapped_loc);
			if(subread0 < 0){
				subread0 = sri+1;
				offset_cur = offsets[sri];
			}
			high_conf_index = offsets[sri];
		}
		if(offset_cur != offsets[sri] && subread0>0){
			alnrec[toli]= subread0;
			alnrec[toli +1]= sri; // sri is the "next section", which is the subread_last + 1
			alnrec[toli +2]= offset_cur;

			offset_cur = offsets[sri];
			subread0 = sri+1;
			toli +=3;
			if(toli >= MAX_INDEL_TOLERANCE*3){
				//SUBREADprintf("TOLIbreak %d\n", toli);
				break;
			}else alnrec[toli+3]=0;
		}
		if(offsets[sri]==0x77) subread0 =-1;
		if(sri == applied_subreads_per_strand-1){
			if(subread0 >0){
				alnrec[toli]= subread0;
				alnrec[toli +1]= sri+1;
				alnrec[toli +2]= offset_cur;
				if( toli +3 < MAX_INDEL_TOLERANCE*3 )alnrec[toli+3]=0;
				toli+=3;
			}
			break;
		}
	}

//	int first_covered_base = find_subread_end(read_len, all_subreads , first_cov_subreadNo);
//	int last_covered_base = find_subread_end(read_len, all_subreads , last_cov_subreadNo) +16 -1;
	return high_conf_index;
}

int cellCounts_matchBin_chro(char * read_bin, int base_offset, gene_value_index_t * index, unsigned int pos, int test_len){
	int ret = 0;

	unsigned int idx_offset_byte, idx_offset_bit;
	unsigned int rbin_offset_byte, rbin_offset_bit;
	gvindex_baseno2offset_m(pos, index , idx_offset_byte, idx_offset_bit);
	if(idx_offset_byte >= index-> values_bytes)return 0;
	char idx_intv = index->values [idx_offset_byte];

	rbin_offset_byte = base_offset/4;
	rbin_offset_bit = (base_offset*2)%8;
	char read_intv = read_bin[rbin_offset_byte];
	int read_i;
	for(read_i = 0; read_i < test_len ; read_i ++){
		char tt = (idx_intv >> idx_offset_bit) & 3;
		char tv = (read_intv >> rbin_offset_bit) & 3;
		if(tt == tv)ret++;
		idx_offset_bit+=2;
		if(idx_offset_bit==8){
			idx_offset_byte++;
			if(idx_offset_byte == index-> values_bytes)return 0;
			idx_intv = index->values [idx_offset_byte];
			idx_offset_bit = 0;
		}
		rbin_offset_bit+=2;
		if(rbin_offset_bit==8){
			rbin_offset_byte ++;
			read_intv = read_bin[rbin_offset_byte];
			rbin_offset_bit =0;
		}
	}
	return ret;
}

#define MINM_INVALID_INDEL (-9999999)

int cellCounts_indel_meet_in_the_middle(cellcounts_global_t * cct_context, int thread_no, unsigned int first_half_abs_pos, char * read_bin, int read_bin_base, int gap_on_read_len, int expected_indel_len, char * read_name, int * gap_mismatch){
	int second_half_start_in_read; 
	unsigned short first_half_matched[MAX_SCRNA_READ_LENGTH], second_half_matched[MAX_SCRNA_READ_LENGTH];
	int max_matched_bases = -99999, best_second_start_in_read = MINM_INVALID_INDEL;
	int x1, summ1=0;
	gene_value_index_t * current_value_index = cct_context->value_index;

	for(x1=0; x1<gap_on_read_len; x1++){
		int idx_value = cellCounts_get_index_int(current_value_index, first_half_abs_pos+ x1);
		int read_value = cellCounts_get_read_int(read_bin, read_bin_base+x1);
		first_half_matched[x1] = summ1; // "if second half starts at x1 (included x1), then how many matched in the first half?""
		summ1 += idx_value == read_value;
	}

	summ1 = 0;
	int indel_offset_first = max(0,(-expected_indel_len));
	for(x1=gap_on_read_len-1; x1>=indel_offset_first; x1--){
		int idx_value = cellCounts_get_index_int(current_value_index, first_half_abs_pos+ x1 +expected_indel_len);
		int read_value = cellCounts_get_read_int(read_bin, read_bin_base+x1); // "if second half starts at x1 (included x1), then how many matched in the second half?"
		summ1 += idx_value == read_value;
		second_half_matched[x1] = summ1;
	}

	if(1) for(second_half_start_in_read = indel_offset_first; second_half_start_in_read < gap_on_read_len; second_half_start_in_read++){
		int sum_here = first_half_matched[ second_half_start_in_read - indel_offset_first ] + second_half_matched[second_half_start_in_read];
		if(sum_here > max_matched_bases){
			max_matched_bases = sum_here ;
			best_second_start_in_read = second_half_start_in_read ;
		}
	}else for(second_half_start_in_read = max(0,(-expected_indel_len)); second_half_start_in_read < gap_on_read_len; second_half_start_in_read++){
		int first_half_length = second_half_start_in_read - max(0,(-expected_indel_len));
		int second_half_length = gap_on_read_len - second_half_start_in_read;
		unsigned int second_half_abs_pos = first_half_abs_pos + second_half_start_in_read + expected_indel_len;
		int first_half_matched = cellCounts_matchBin_chro(read_bin, read_bin_base ,  cct_context -> value_index, first_half_abs_pos, first_half_length);
		int second_half_matched = cellCounts_matchBin_chro(read_bin, read_bin_base + second_half_start_in_read,  cct_context -> value_index, second_half_abs_pos, second_half_length);
		int both_matched = first_half_matched + second_half_matched;
		if(both_matched>max_matched_bases){
			max_matched_bases = both_matched;
			best_second_start_in_read = second_half_start_in_read;
		}
	}
	(* gap_mismatch) = gap_on_read_len - max_matched_bases + min(0, expected_indel_len);
//if(0)SUBREADprintf("FOUND %d indel MEET at %d in %d gap ; %d matched,\n", expected_indel_len, best_second_start_in_read, gap_on_read_len, max_matched_bases );
	return best_second_start_in_read + min(0, expected_indel_len);
}

srInt_64 cellCounts_test_score(cellcounts_global_t * cct_context, int thread_no, char * read_name, int read_len, unsigned int abs_pos, char * cigar, int head_soft_clipped, int tail_soft_clipped, int all_matched_bases, int all_mismatched_bases){
	if(all_mismatched_bases > cct_context -> max_mismatching_bases_in_reads) return 0;
	return all_matched_bases*1000000llu / (1llu+all_mismatched_bases);
}

void cellCounts_chroEvent_locks_opt(cellcounts_global_t * cct_context, int thread_no, void * entity, int tolock){
	srInt_64 vn = entity-NULL;
	vn = vn & 0xffffffffffffllu;
	int lock_no = vn % (1LLU*cct_context -> chroEvent_lock_number);
	if(tolock) cellCounts_lock_occupy(cct_context -> read_assignment_counter_locks + lock_no);
	else cellCounts_lock_release(cct_context -> read_assignment_counter_locks + lock_no);
}

unsigned int cellCounts_convert_stack_to_cigar_str(cellcounts_global_t * cct_context, int thread_no, realignment_event_stack_item_t * stack, char * cigar, int to3end){
	int stack_depth = 0,x1, cigar_ptr=0;

	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	while(stack_depth<JUNCTION_REALIGNMENT_MAX_DEPTH){
		if(stack[stack_depth].event_details == NULL+IMPOSSIBLE_MEMORY_SPACE)break;
		stack_depth++;
	}

	unsigned int read_base1_linear = 0;
	for(x1=0; x1<stack_depth;x1++){
		int stidx = to3end?x1:( stack_depth-1-x1 );
		realignment_event_stack_item_t * oneitem = stack + stidx;
		if(x1==0 && !to3end) read_base1_linear = oneitem -> linear_l - oneitem -> read_covered_last_base; 

		if(!to3end)cigar_ptr += SUBreadSprintf( cigar+cigar_ptr, 12, "%dM", oneitem -> read_covered_last_base - oneitem -> read_covered_first_base +1);

		chroEvent_t * event_details = NULL;
		if(stidx>0) event_details = oneitem -> event_details;

		if( event_details ){
/*
			chroEvent_t * event_details;
			if(to3end) event_details = oneitem -> event_details;
			else event_details = oneitem -> event_details;
*/
			int Nlen;
			if(event_details -> n_events == 1) Nlen = (void*) event_details -> length -NULL;
			else Nlen = event_details -> length[oneitem -> insertion_idx];

			char Nmode;
			if(event_details -> event_type == chroEvent_t_TYPE_JUNCTION) Nmode = 'N';
			else if( Nlen > 0 ) Nmode = 'D';
			else Nmode = 'I';
			cigar_ptr += SUBreadSprintf( cigar+cigar_ptr,11, "%d%c", abs(Nlen), Nmode); 
			//cellCounts_chroEvent_locks_opt(cct_context, thread_no, oneitem -> event_details, 1);
			// oneitem -> event_details -> step2_supported_reads++;
			//cellCounts_chroEvent_locks_opt(cct_context, thread_no, oneitem -> event_details, 0);
		}
		if(to3end)cigar_ptr += SUBreadSprintf( cigar+cigar_ptr, 11, "%dM", oneitem -> read_covered_last_base - oneitem -> read_covered_first_base +1);
	}

	return read_base1_linear; // if to 3 end: ignore the retured value.
}

char cellCounts_softclip_getbase(unsigned int pos, void * context){
	cellcounts_global_t * cct_context = context;
	gene_value_index_t * current_value_index = cct_context->value_index;
	return gvindex_get(current_value_index, pos);
}

unsigned int cellCounts_softclip_candidate(cellcounts_global_t * cct_context, int thread_no, unsigned int linear_mapped_pos, char * cigar, char * read_text, int noindel_cover_firstbase, int noindel_cover_lastbase, int allowd_mismatching_in_window, srInt_64 * ma_misma_ins_Sclip){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	SoftClipResult * clipres=calculate_soft_clipping(cct_context, linear_mapped_pos, cigar, read_text, noindel_cover_firstbase, noindel_cover_lastbase, allowd_mismatching_in_window, cellCounts_softclip_getbase); 

	strcpy(cigar, clipres -> new_cigar);
	(* ma_misma_ins_Sclip)=
		((1LL<<48)*clipres -> num_matched) |
		((1LL<<32)*clipres -> num_mismatched) |
		((1LL<<16)*clipres -> num_inserted) |
		( 1LL     *clipres -> num_clipped);

	unsigned int ret = clipres -> new_pos;
	free(clipres);
	return ret;
}

void cellCounts_end_build_candidature_from_stacks(cellcounts_global_t * cct_context, int thread_no, int noindel_coved_firstbase, int noindel_coved_lastbase, int is_reversed, int read_len, unsigned int default_mapped_linear, int vote_for_aln, char * read_text, int total_score){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;

	int x_3end, x_5end, x_candidate;
	for(x_5end = -1; x_5end < thread_context -> best_5end_stack_list -> numOfElements; x_5end++){
		if(thread_context -> populating_voteIJ_buf_index >= cct_context -> max_candidate_voteIJ_per_read)break;

		unsigned int mapped_1st_loc = 0;
		char cigar_being_built[11*(1+4*JUNCTION_REALIGNMENT_MAX_DEPTH)];

		if( x_5end < 0 && thread_context -> best_5end_stack_list -> numOfElements >0 )continue; // good 5' stack doesn't need full-body case.

		int stack_depth_5_3_ends=0;
		if(x_5end >=0){
			realignment_event_stack_item_t * stack_end5 = ArrayListGet(thread_context -> best_5end_stack_list , x_5end);
			mapped_1st_loc = cellCounts_convert_stack_to_cigar_str(cct_context, thread_no, stack_end5, cigar_being_built, 0);
			while(stack_depth_5_3_ends<JUNCTION_REALIGNMENT_MAX_DEPTH){
				if(stack_end5[stack_depth_5_3_ends].event_details == NULL+IMPOSSIBLE_MEMORY_SPACE)break;
				stack_depth_5_3_ends++;
			}
		} else cigar_being_built[0]=0;

		if(!mapped_1st_loc) mapped_1st_loc = default_mapped_linear;

		int sub1=0;
		if(!thread_context -> best_3end_stack_list -> numOfElements) sub1++;
		if(!thread_context -> best_5end_stack_list -> numOfElements) sub1++;

		int cigar_ptr = strlen(cigar_being_built);
		cigar_ptr += SUBreadSprintf( cigar_being_built + cigar_ptr , 11, "%dM", noindel_coved_lastbase - noindel_coved_firstbase -1+ sub1); // the first and last covered base in the non-event region are also included in the first item in the stack.

		for(x_3end = -1; x_3end < thread_context -> best_3end_stack_list -> numOfElements; x_3end++){
			unsigned int mapped_loc = mapped_1st_loc;
			if(thread_context -> populating_voteIJ_buf_index >= cct_context -> max_candidate_voteIJ_per_read)break;
			if(x_3end <0 && thread_context -> best_3end_stack_list -> numOfElements >0)continue;

			int stack_depth=0;
			if(x_3end >= 0){
				realignment_event_stack_item_t * stack_end3 = ArrayListGet(thread_context -> best_3end_stack_list , x_3end);
				cellCounts_convert_stack_to_cigar_str(cct_context, thread_no, stack_end3, cigar_being_built+ cigar_ptr, 1);
				while(stack_depth<JUNCTION_REALIGNMENT_MAX_DEPTH){
					if(stack_end3[stack_depth].event_details == NULL+IMPOSSIBLE_MEMORY_SPACE)break;
					stack_depth++;
				}
			}
			stack_depth_5_3_ends += stack_depth;

			char * final_cigar = thread_context -> reporting_cigars[thread_context -> populating_voteIJ_buf_index];
			int sum_ins=0;
			int cigar_rlen = cellCounts_reduce_Cigar(cigar_being_built,final_cigar,&sum_ins);

			srInt_64 do_ma_misma_ins_Sclip = 0;

			// repeated reporting check don't consider the soft clipping.
			srInt_64 hkey = HashTableStringHashFunction(final_cigar);
			hkey = (hkey<<24) ^ mapped_loc;

			// "S" will be removed later on from the mapping location. Here we don't move it!
			mapped_loc = cellCounts_softclip_candidate(cct_context, thread_no, mapped_loc, final_cigar, read_text, noindel_coved_firstbase, noindel_coved_lastbase, cct_context -> enable_soft_clipping?1:9999, &do_ma_misma_ins_Sclip);
			mapped_loc -= get_soft_clipping_length(final_cigar);

			int mismatching_bases = (int)((do_ma_misma_ins_Sclip>>32) & 0xffffllu);
			if(mismatching_bases > cct_context -> max_mismatching_bases_in_reads) continue;
			int matched_bases = (int)((do_ma_misma_ins_Sclip>>48) & 0xffffllu);

			void * found_reporting = HashTableGet(thread_context -> alignment_repating_table, NULL+hkey) ;
			if( found_reporting )continue;
			HashTablePut(thread_context -> alignment_repating_table, NULL+hkey, NULL+1);

			srInt_64 weight = cellCounts_calculate_pos_weight(cct_context, mapped_loc, final_cigar);
			srInt_64 score  = ((255ll-mismatching_bases)<<16) + (matched_bases<<8)  + (255ll - stack_depth_5_3_ends);  // #Mismatch first -- finally used

			thread_context -> reporting_scores[thread_context -> populating_voteIJ_buf_index] = score * weight;
			thread_context -> reporting_positions[thread_context -> populating_voteIJ_buf_index] = mapped_loc;
			thread_context -> reporting_flags[thread_context -> populating_voteIJ_buf_index] = is_reversed?SAM_FLAG_REVERSE_STRAND_MATCHED:0;
			thread_context -> reporting_mapq[thread_context -> populating_voteIJ_buf_index] = 40 ;//- all_mismatched_bases;
			thread_context -> reporting_editing_distance[thread_context -> populating_voteIJ_buf_index] = 0;// all_mismatched_bases + all_indel_length;
			thread_context -> reporting_vote_for_aln[thread_context -> populating_voteIJ_buf_index] = vote_for_aln;
			thread_context -> reporting_ma_misma_ins_Sclip[thread_context -> populating_voteIJ_buf_index] = do_ma_misma_ins_Sclip;

			thread_context -> populating_voteIJ_buf_index ++;
			thread_context -> total_voteIJs_to_write ++;
		}
	}

	ArrayListDestroy(thread_context -> best_5end_stack_list);
	ArrayListDestroy(thread_context -> best_3end_stack_list);
	thread_context -> best_5end_stack_list = thread_context -> best_3end_stack_list = NULL;
}

void cellCounts_reset_3end_5end_best_stacks(cellcounts_global_t * cct_context, int thread_no, int to3end){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	if(to3end){
		if(thread_context -> best_3end_stack_list) ArrayListDestroy(thread_context -> best_3end_stack_list);
		thread_context -> best_3end_stack_list = ArrayListCreate(5);
		ArrayListSetDeallocationFunction(thread_context -> best_3end_stack_list ,free);
	}else{
		if(thread_context -> best_5end_stack_list) ArrayListDestroy(thread_context -> best_5end_stack_list);
		thread_context -> best_5end_stack_list = ArrayListCreate(5);
		ArrayListSetDeallocationFunction(thread_context -> best_5end_stack_list ,free);
	}
}

void cellCounts_init_build_junctionread_reset_best_oneend(cellcounts_global_t * cct_context, int thread_no, int to3end, srInt_64 bsscore){
	cellCounts_reset_3end_5end_best_stacks(cct_context, thread_no, to3end);
}

void cellCounts_init_build_junctionread_context(cellcounts_global_t * cct_context, int thread_no, char * read_name){
	int to3end;

	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	thread_context -> realignment_event_read_name = read_name;
	for(to3end=0;to3end<2;to3end++) cellCounts_init_build_junctionread_reset_best_oneend(cct_context, thread_no, to3end, -1); // Reset left and right alignments. Left and right alignments are independent.
}

void cellCounts_build_junction_read_set_current_stack(cellcounts_global_t * cct_context, int thread_no,
	unsigned int linear_env_L, unsigned int linear_env_R, int insertion_length_index, short total_match, short total_mismatch,
	short covered_first_base_in_read, short covered_last_base_in_read, chroEvent_t * event_details){  

	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth -1] . matching_bases_in_alignment = total_match; 
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth -1] . mismatching_bases_in_alignment = total_mismatch; 
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth -1] . read_covered_first_base = covered_first_base_in_read;
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth -1] . read_covered_last_base = covered_last_base_in_read;

	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth]    . event_details = event_details;
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth]    . linear_l = linear_env_L;
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth]    . linear_r = linear_env_R;
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth]    . insertion_idx = insertion_length_index;

}


int cellCounts_build_junction_read_finalise_current_stack(cellcounts_global_t * cct_context, int thread_no,
	short total_match, short total_mismatch, short covered_first_base_in_read, short covered_last_base_in_read, int to3end){

	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth -1] . matching_bases_in_alignment = total_match; 
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth -1] . mismatching_bases_in_alignment = total_mismatch; 
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth -1] . read_covered_first_base = covered_first_base_in_read;
	thread_context -> realignment_event_current_stack[ thread_context -> realignment_event_stack_current_depth -1] . read_covered_last_base = covered_last_base_in_read;

	int all_match_in_stack = 0, all_mismatch_in_stack = 0,x1;

	for(x1=0; x1<= thread_context -> realignment_event_stack_current_depth -1; x1++){
		// NB: insertions in reads are nither match nor mismatch.
		all_match_in_stack += thread_context -> realignment_event_current_stack[x1].matching_bases_in_alignment;
		all_mismatch_in_stack += thread_context -> realignment_event_current_stack[x1].mismatching_bases_in_alignment;
	}

	int my_score = all_match_in_stack * 10 + (10- thread_context -> realignment_event_stack_current_depth);
	if(my_score > thread_context -> realignment_event_stack_best_score){
		thread_context -> realignment_event_stack_best_score = my_score;
		cellCounts_reset_3end_5end_best_stacks(cct_context, thread_no, to3end);
	}
	if(my_score == thread_context -> realignment_event_stack_best_score){
		realignment_event_stack_item_t* stack_copy_ptr = malloc(sizeof(realignment_event_stack_item_t) * JUNCTION_REALIGNMENT_MAX_DEPTH);
		memcpy( stack_copy_ptr , thread_context -> realignment_event_current_stack , sizeof(realignment_event_stack_item_t) *  thread_context -> realignment_event_stack_current_depth  );
		if( thread_context -> realignment_event_stack_current_depth < JUNCTION_REALIGNMENT_MAX_DEPTH ) stack_copy_ptr[thread_context -> realignment_event_stack_current_depth].event_details = NULL+IMPOSSIBLE_MEMORY_SPACE;
		ArrayListPush( to3end? thread_context -> best_3end_stack_list : thread_context -> best_5end_stack_list , stack_copy_ptr);
	}


//#warning "========= THIS IS FOR TESTING EFFICIENCY; NOT NEEDED IN RELASED VERSIOn  ============="
	thread_context -> realignment_event_stack_runcount++;

	return all_mismatch_in_stack;
}

void cellCounts_tree_iterative_search( cellcounts_global_t * cct_context, int thread_no, int sample_i, char * chro,  int this_end_last_correct_maiping_chro, int this_end_last_correct_mapping_read, char * read_name, char * read_text, int read_len, int to3end, int cellbc_no){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	if( thread_context -> realignment_event_stack_runcount > JUNCTION_REALIGNMENT_MAX_TRIES)return;
	thread_context -> realignment_event_stack_current_depth ++;
	if(thread_context -> realignment_event_stack_current_depth > JUNCTION_REALIGNMENT_MAX_DEPTH){
		SUBREADprintf("SHOULND'T REACH HERE: %d > %d\n", thread_context -> realignment_event_stack_current_depth,  JUNCTION_REALIGNMENT_MAX_DEPTH);
		return;
	}
	gene_value_index_t * current_value_index = cct_context->value_index;

	int edge_region_start = to3end?this_end_last_correct_maiping_chro:(this_end_last_correct_maiping_chro - this_end_last_correct_mapping_read);  // including the last base in read.

	int remaining_bases_to_scan = to3end?(read_len - 1 - this_end_last_correct_mapping_read):this_end_last_correct_mapping_read;
	int edge_region_end = edge_region_start + remaining_bases_to_scan; // inclusive, as defined by the IVT_query_range function.
	int scan_len = remaining_bases_to_scan + 1;

	unsigned short mismatch_prefix[scan_len + 1];
	int total_mismatching = 0, total_matching = 0;
	int read_pos = this_end_last_correct_mapping_read;
	unsigned int chro_linear_pos = linear_gene_position(&cct_context->chromosome_table, chro , this_end_last_correct_maiping_chro);
	int x1delta = to3end?1:-1;
	int scan_i;

	mismatch_prefix[0] = 0;
	for(scan_i = 0; scan_i < scan_len && read_pos >=0 && read_pos < read_len; scan_i++){ 	// When last_correct_mapping_read is 0, read_pos can be -1. 
							// Similarly, if last_correct_mapping_read is read_len -1, read_pos can be read_len. 
		char chr_base = gvindex_get(current_value_index, chro_linear_pos);
		char read_base = read_text[read_pos];
		int this_base_misma = read_base!=chr_base;
		mismatch_prefix[scan_i + 1] = mismatch_prefix[scan_i] + this_base_misma;
		chro_linear_pos += x1delta;
		read_pos += x1delta;
	}
	int scanned_len = scan_i;
	total_mismatching = mismatch_prefix[scanned_len];
	total_matching = scanned_len - total_mismatching;

	int cov_base0 =   to3end ?this_end_last_correct_mapping_read:0;
	int cov_base1 =   to3end ?read_len - 1: this_end_last_correct_mapping_read;
	// no matter if there are events or not, always test using all.
	int margin_for_misma = cellCounts_build_junction_read_finalise_current_stack(cct_context, thread_no, total_matching, total_mismatching, cov_base0 , cov_base1, to3end);

	int founditems = 0, x1;
	int eventbufsize = JUNCTION_MAX_COLOCATION;
	IVT_Interval * eventbuf[eventbufsize];

	char chro_strn_ky[MAX_CHROMOSOME_NAME_LEN+10];
	char negchar = '*';
	snprintf(chro_strn_ky, MAX_CHROMOSOME_NAME_LEN+10, "%s\t%c", chro, negchar);
	IVT_IntervalTreeNode **LR_roots = HashTableGet(cct_context -> chroEvent_entry_table[sample_i], chro_strn_ky); // LR_roots [0] : left-edge ; LR_roots [1] : right-edge ; LR_roots [2] : covered-range
	if(margin_for_misma >0 && LR_roots && thread_context -> realignment_event_stack_current_depth < JUNCTION_REALIGNMENT_MAX_DEPTH){
		IVT_IntervalTreeNode * myroot = NULL;
		myroot=LR_roots[!to3end];
		IVT_query_range(myroot , edge_region_start, edge_region_end, eventbuf, eventbufsize, &founditems);
	}

	if(founditems){
		int my_env_total_mismatch[founditems], my_env_total_match[founditems];

		memset(my_env_total_match, 0, sizeof(int)*founditems);
		memset(my_env_total_mismatch, 0, sizeof(int)*founditems);

		if(founditems >= eventbufsize - 1)SUBREADprintf("Warning: there are %d chromosomal events found in a read region. This is abnormally too many.\n", founditems);

		// for each event that don't have too many mismatching bases locally, go deeper. 
		for(x1=0; x1<founditems; x1++){
			int event_offset = to3end ? (eventbuf[x1] -> start - this_end_last_correct_maiping_chro) : (this_end_last_correct_maiping_chro - eventbuf[x1] -> start);
			int event_prefix_len;
			if(event_offset < 0 || event_offset >= scanned_len) continue;
			event_prefix_len = event_offset + 1;
			my_env_total_mismatch[x1] = mismatch_prefix[event_prefix_len];
			my_env_total_match[x1] = event_prefix_len - my_env_total_mismatch[x1];
			if(my_env_total_mismatch[x1]<=min(JUNCTION_MAX_MISMATCHING_BASES_IN_REALIGNMENT, margin_for_misma -1)){
				IVT_Interval * evb = eventbuf[x1];
				int evbposleft = evb -> start;
				int evbposright = evb -> attr - NULL;
				unsigned int linear_env_1 = linear_gene_position(&cct_context->chromosome_table, chro , evbposleft);
				unsigned int linear_env_2 = linear_gene_position(&cct_context->chromosome_table, chro , evbposright);
				unsigned int linear_env_L = min(linear_env_1, linear_env_2);
				unsigned int linear_env_R = max(linear_env_1, linear_env_2);

				if(cct_context -> cluster_spec_junction_table){
					int my_cluster = HashTableGet(cct_context -> cluster_cell_map_table, NULL+(sample_i *1LLU<<56)+cellbc_no)-NULL;
					void *event_known = NULL;
					if(my_cluster){
						HashTable * my_known_tab = HashTableGet(cct_context -> cluster_spec_junction_table, NULL+(sample_i *1LLU<<56)+my_cluster);
						if(my_known_tab){
							event_known = HashTableGet(my_known_tab, NULL+(linear_env_L*1LLU<<32)+linear_env_R);
						}
					}
					if(!event_known) continue;
				}

				srUInt_64 envkey = (linear_env_L*1LLU<<32)| linear_env_R;
				chroEvent_t * envdtl = HashTableGet( cct_context -> chroEvent_detail_table[sample_i], NULL+ envkey);
	
				if(envdtl -> event_type >= chroEvent_t_TYPE_EXON) continue;
				int remote_first_matching_base_chro = evb -> attr - NULL, x2;
				for(x2=0; x2< envdtl -> n_events; x2++){
					int inslen = 0;
					if(envdtl -> event_type == chroEvent_t_TYPE_INDEL){
						if( envdtl -> n_events >1 ) inslen = envdtl -> length[x2];
						else inslen = (void*)envdtl -> length - NULL;
						if(inslen>0)inslen=0; else inslen = -inslen;
					}
					int remote_first_matching_base_read = inslen * x1delta; 
					int evb_pos_in_read_after_last_correct_base = evb -> start - this_end_last_correct_maiping_chro; // can be negative if search is to 5'.
					remote_first_matching_base_read += x1delta + evb_pos_in_read_after_last_correct_base + this_end_last_correct_mapping_read; 

					if(remote_first_matching_base_read >0 && remote_first_matching_base_read < read_len){
						int cov_base0 =   to3end ?this_end_last_correct_mapping_read:(evb_pos_in_read_after_last_correct_base + this_end_last_correct_mapping_read);
						int cov_base1 = (!to3end)?this_end_last_correct_mapping_read:(evb_pos_in_read_after_last_correct_base + this_end_last_correct_mapping_read);



						cellCounts_build_junction_read_set_current_stack(cct_context, thread_no, linear_env_L, linear_env_R, x2,
							my_env_total_match[x1], my_env_total_mismatch[x1], cov_base0 , cov_base1, envdtl);
						cellCounts_tree_iterative_search( cct_context, thread_no, sample_i, chro, remote_first_matching_base_chro, remote_first_matching_base_read, read_name, read_text, read_len, to3end , cellbc_no);
					}
				}
			}
		}
	}
	thread_context -> realignment_event_stack_current_depth --;
}



void cellCounts_end_junctionread_one_end(cellcounts_global_t * cct_context, int thread_no){
}

void cellCounts_init_junctionread_one_end(cellcounts_global_t * cct_context, int thread_no){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;

	thread_context -> realignment_event_stack_current_depth = 0;
	thread_context -> realignment_event_stack_best_score = -1;
}

void cellCounts_explain_one_alignment(cellcounts_global_t * cct_context, int thread_no, int sample_i, char * read_name, char * read_text, int read_len, int noindel_coved_firstbase, int noindel_coved_lastbase, unsigned int linear_mapped_pos, int is_reversed, int votes_for_aln){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	cellCounts_init_build_junctionread_context(cct_context, thread_no, read_name);
	int to3end;

	char * sample_seq=NULL, *sample_qual=NULL, *BC_qual=NULL, *BC_seq=NULL, *UMI_seq=NULL, *UMI_qual=NULL, *lane_str=NULL, *RG=NULL, *testi;
	int rname_trimmed_len=0;
	cellCounts_scan_read_name_str(cct_context, NULL, read_name, &sample_seq, &sample_qual, &BC_seq, &BC_qual, &UMI_seq, &UMI_qual, &lane_str, &RG, &rname_trimmed_len);
	int cell_barcode_no = cellCounts_get_cellbarcode_no(cct_context, thread_no, BC_seq, BC_qual);

	noindel_coved_firstbase += JUNCTION_WIDDEN_GAP_LEN; // widden the gap to avoid same bases before/after event
	noindel_coved_lastbase -= JUNCTION_WIDDEN_GAP_LEN;
	int total_score = (noindel_coved_lastbase - noindel_coved_firstbase+1)*10; // my_score = all_match_in_stack * 10 + (10- thread_context -> realignment_event_stack_current_depth);

	for(to3end=0; to3end<2; to3end++){
		cellCounts_init_junctionread_one_end(cct_context, thread_no);

		unsigned int this_end_last_correct_mapping_chro = to3end?(linear_mapped_pos + noindel_coved_lastbase ): (linear_mapped_pos + noindel_coved_firstbase) ;
		char * chro_name = NULL;
		int this_end_last_correct_mapping_read = to3end?noindel_coved_lastbase:noindel_coved_firstbase, chro_pos = 0;
		if(to3end && this_end_last_correct_mapping_read == read_len-1) continue;
		if(0==to3end && this_end_last_correct_mapping_read == 0) continue;
		thread_context -> realignment_event_current_stack[0].linear_l  = thread_context -> realignment_event_current_stack[0].linear_r = this_end_last_correct_mapping_chro;

		locate_gene_position(this_end_last_correct_mapping_chro, &cct_context -> chromosome_table, &chro_name, &chro_pos);
		thread_context -> realignment_event_stack_runcount = 0;
		cellCounts_tree_iterative_search( cct_context, thread_no, sample_i, chro_name, chro_pos, this_end_last_correct_mapping_read, read_name, read_text, read_len, to3end, cell_barcode_no);
		cellCounts_end_junctionread_one_end(cct_context, thread_no);
		if(thread_context -> realignment_event_stack_best_score>0) total_score += thread_context -> realignment_event_stack_best_score;
	}
	cellCounts_end_build_candidature_from_stacks(cct_context, thread_no, noindel_coved_firstbase, noindel_coved_lastbase, is_reversed, read_len, linear_mapped_pos + noindel_coved_firstbase, votes_for_aln, read_text, total_score);
	return;
}

// do: 
//   1, indel detection (meet-in-the-middle or Smith-Waterman)
//   2, build CIGAR
//   3, calculate matched/mismatched
//   4, calculate and save scores in array
srInt_64 cellCounts_explain_PaperVersion_one_alignment(cellcounts_global_t * cct_context, int thread_no,char * read_name, char * read_bin, char * read_text, int read_len,  cellcounts_vote_number_t all_subreads, gene_sc_vote_t * votetab, int vote_i, int vote_j){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	cellcounts_vote_number_t indel_offsets [MAX_INDEL_TOLERANCE*3];
	int first_mapped_base_offset = 0, toli, in_cigar_readlen = 0, all_mismatched_bases = 0, all_matched_bases = 0, all_mapped_bases = 0;
	char newcigar[30];

	int rbin_offset_for_reversed = (votetab -> masks[vote_i][vote_j] & IS_NEGATIVE_STRAND)? REVERSED_READ_BIN_OFFSET :0;
	unsigned int abs_pos = votetab -> pos[vote_i][vote_j];
	int tolimax = votetab -> toli[vote_i][vote_j], all_indel_length=0;
	cellCounts_indel_recorder_copy(indel_offsets, votetab -> indel_recorder[vote_i][vote_j], tolimax, all_subreads, &first_mapped_base_offset, NULL, NULL, read_name, abs_pos);
	abs_pos += first_mapped_base_offset; // N.B., first_mapped_base_offset is the offset of the first mapped base, not the read-pos of the first mapped base.
					     // The extraction location of first base is added below ("meet_start"). 
	signed int last_indel = indel_offsets[2], last_section_subread_no = indel_offsets[1]-1, head_soft_clipped = -1, last_mapped_base_in_read = 0;
	thread_context -> reporting_cigars[thread_context -> populating_voteIJ_buf_index][0]=0;
	for(toli = 3; toli < tolimax; toli+=3){
		if(indel_offsets[toli]==0)break;
		int hiconf_vote_first = indel_offsets[toli]-1;
		int hiconf_vote_last = indel_offsets[toli+1]-1;
		signed int indel_offset = indel_offsets[toli+2];
		signed int indel_diff = indel_offset - last_indel;
		if(abs(indel_diff)>= cct_context -> max_indel_length) continue;
		int last_correct_base = find_subread_end(read_len, all_subreads , last_section_subread_no) - 8;
		int first_correct_base = find_subread_end(read_len, all_subreads , hiconf_vote_first) - 16 + 8;

		if(last_correct_base < in_cigar_readlen) last_correct_base= in_cigar_readlen;

		unsigned int meet_start = abs_pos + last_correct_base + last_indel;
		if(head_soft_clipped <0 && cct_context -> enable_soft_clipping ){
			int last_mapped_base_in_read = find_subread_end(read_len, all_subreads , indel_offsets[0]-1);
			head_soft_clipped = cellCounts_find_soft_clipping(cct_context, thread_no, read_bin+rbin_offset_for_reversed, 0, abs_pos /* this can only happen if no indel is in read */, last_mapped_base_in_read , 0, last_mapped_base_in_read);
			if(head_soft_clipped > 0)SUBreadSprintf(thread_context -> reporting_cigars[thread_context -> populating_voteIJ_buf_index], 11,"%dS", head_soft_clipped );
			if(meet_start < head_soft_clipped + abs_pos ) meet_start= head_soft_clipped + abs_pos;
			if(head_soft_clipped > last_correct_base) last_correct_base = head_soft_clipped;
			in_cigar_readlen = head_soft_clipped;
		}else if( !cct_context -> enable_soft_clipping ) head_soft_clipped = 0;

		int gap_mismatched = 0;
		int indel_pos = cellCounts_indel_meet_in_the_middle(cct_context, thread_no, meet_start, read_bin + rbin_offset_for_reversed, last_correct_base, first_correct_base - last_correct_base, indel_diff, read_name, &gap_mismatched);
		if(indel_pos == MINM_INVALID_INDEL){
			indel_pos = (first_correct_base - last_correct_base)/2;
			all_indel_length += abs(indel_diff);
		}


		int section_matched = cellCounts_matchBin_chro(read_bin +rbin_offset_for_reversed , in_cigar_readlen , cct_context -> value_index, abs_pos + in_cigar_readlen + last_indel , last_correct_base - in_cigar_readlen);

		all_mismatched_bases += gap_mismatched + (last_correct_base - in_cigar_readlen - section_matched);
		all_matched_bases += section_matched + first_correct_base - last_correct_base - gap_mismatched + min(0, indel_diff);

		SUBreadSprintf(newcigar, 45, "%dM%dM%d%c%dM", last_correct_base - in_cigar_readlen , indel_pos, abs(indel_diff), indel_diff>0?'D':'I', first_correct_base - last_correct_base - indel_pos + min(0, indel_diff) );
		all_mapped_bases += first_correct_base - in_cigar_readlen + min(0, indel_diff);
		strcat(thread_context -> reporting_cigars[thread_context -> populating_voteIJ_buf_index], newcigar);
		in_cigar_readlen = first_correct_base;
		last_mapped_base_in_read = find_subread_end(read_len, all_subreads , hiconf_vote_last) - 16 + 9;
		last_indel = indel_offset;
		last_section_subread_no = hiconf_vote_last;
	}

	if(head_soft_clipped <0 && cct_context -> enable_soft_clipping){
		int first_mapped_base_in_read = find_subread_end(read_len, all_subreads , indel_offsets[0]-1) - 8;
		head_soft_clipped = cellCounts_find_soft_clipping(cct_context, thread_no, read_bin+rbin_offset_for_reversed,0, abs_pos /* this can only happen if no indel is in read */, first_mapped_base_in_read , 0, first_mapped_base_in_read);

		if(head_soft_clipped > 0)SUBreadSprintf(thread_context -> reporting_cigars[thread_context -> populating_voteIJ_buf_index], 11,"%dS", head_soft_clipped );
		in_cigar_readlen = head_soft_clipped;
		last_mapped_base_in_read = find_subread_end(read_len, all_subreads, last_section_subread_no) - 16 + 8;
	}else if( !cct_context -> enable_soft_clipping ) head_soft_clipped = 0;

	int tail_soft_clipped = 0;
	if( cct_context -> enable_soft_clipping ) tail_soft_clipped = cellCounts_find_soft_clipping(cct_context, thread_no, read_bin+rbin_offset_for_reversed, last_mapped_base_in_read, abs_pos + last_mapped_base_in_read  + last_indel,  read_len - last_mapped_base_in_read , 1, 1);
	int section_matched = cellCounts_matchBin_chro(read_bin +rbin_offset_for_reversed , in_cigar_readlen , cct_context -> value_index, abs_pos + in_cigar_readlen + last_indel, read_len - in_cigar_readlen - tail_soft_clipped);
	all_mismatched_bases += (read_len - in_cigar_readlen - tail_soft_clipped  - section_matched);
	all_matched_bases += section_matched;

	SUBreadSprintf(newcigar, 11, "%dM", read_len - in_cigar_readlen - tail_soft_clipped);
	all_mapped_bases += read_len - in_cigar_readlen - tail_soft_clipped;
	strcat(thread_context -> reporting_cigars[thread_context -> populating_voteIJ_buf_index], newcigar);
	if(tail_soft_clipped){
		SUBreadSprintf(newcigar, 11, "%dS", tail_soft_clipped);
		strcat(thread_context -> reporting_cigars[thread_context -> populating_voteIJ_buf_index], newcigar);
	}

	char tmp_new_cigar[MAX_SCRNA_READ_LENGTH+20];
	int rebuilt_rlen = cellCounts_reduce_Cigar(thread_context -> reporting_cigars[thread_context -> populating_voteIJ_buf_index], tmp_new_cigar, NULL );

	strcpy(thread_context -> reporting_cigars[thread_context -> populating_voteIJ_buf_index], tmp_new_cigar);

	srInt_64 weight = cellCounts_calculate_pos_weight(cct_context, abs_pos, tmp_new_cigar);
	srInt_64 score = 0;
	if(rebuilt_rlen==read_len && all_mapped_bases >= cct_context -> min_mapped_length_for_mapped_read )score=cellCounts_test_score(cct_context, thread_no, read_name, read_len, abs_pos, thread_context -> reporting_cigars[thread_context -> populating_voteIJ_buf_index], head_soft_clipped, tail_soft_clipped, all_matched_bases, all_mismatched_bases)*weight;

	int votes_for_aln = votetab -> votes[vote_i][vote_j];
	thread_context -> reporting_scores[thread_context -> populating_voteIJ_buf_index] = score;
	thread_context -> reporting_positions[thread_context -> populating_voteIJ_buf_index] = abs_pos;
	thread_context -> reporting_flags[thread_context -> populating_voteIJ_buf_index] = rbin_offset_for_reversed?SAM_FLAG_REVERSE_STRAND_MATCHED:0;
	thread_context -> reporting_mapq[thread_context -> populating_voteIJ_buf_index] = 40 - all_mismatched_bases;
	thread_context -> reporting_editing_distance[thread_context -> populating_voteIJ_buf_index] = all_mismatched_bases + all_indel_length;
	thread_context -> reporting_vote_for_aln[thread_context -> populating_voteIJ_buf_index] = votes_for_aln;
	thread_context -> populating_voteIJ_buf_index++;
	if(score>0)thread_context -> total_voteIJs_to_write++;
	return score;
}

int sort_readscore_compare_LargeFirst(void * vp , int i , int j){
	void ** pp = vp;
	cellcounts_align_thread_t * thread_context = pp[0];
	int * sorting_index = pp[1];
	int idxI = sorting_index [i];
	int idxJ = sorting_index [j];
	if(thread_context -> reporting_scores[idxI] > thread_context -> reporting_scores[idxJ]) return -1; // large number first
	if(thread_context -> reporting_scores[idxI] < thread_context -> reporting_scores[idxJ]) return 1;
	return 0;
}

void sort_readscore_exchange(void * vp , int i , int j){
	void ** pp = vp;
	int * sorting_index = pp[1];
	int idxI = sorting_index [i];
	sorting_index [i] = sorting_index [j];
	sorting_index [j] = idxI ;
}

int cellCounts_junc_meet_in_the_middle(cellcounts_global_t * cct_context, cellcounts_align_thread_t * thread_context, char * read_text, int major_read_start, int major_read_end, int indels_in_major_cov, unsigned int major_chro_loc, int minor_read_start, int minor_read_end, int indels_in_minor_cov, unsigned int minor_chro_loc, int * total_misma, int * outgaplen, int * is_GTAG_negative){
	// "loc" = chro_location
	unsigned int left_last_matched_base_loc = 0;
	unsigned int right_last_matched_base_loc = 0;
	int gaplen = -1, gap_start_in_read=-1; // gap_start = first not mapped base pos in read.
	if(major_read_start >= minor_read_end){
		left_last_matched_base_loc = minor_read_end + minor_chro_loc -1 +indels_in_minor_cov;
		right_last_matched_base_loc = major_read_start + major_chro_loc;
		gaplen = major_read_start - minor_read_end;
		gap_start_in_read= minor_read_end;
	}
	if(minor_read_start >= major_read_end){
		left_last_matched_base_loc = major_read_end + major_chro_loc -1 +indels_in_major_cov;
		right_last_matched_base_loc = minor_read_start + minor_chro_loc;
		gaplen = minor_read_start - major_read_end;
		gap_start_in_read= major_read_end;
	}
	if(0==right_last_matched_base_loc){
		SUBREADprintf("ERROR: two halves overlapped.\n");
		return -1;
	}

	char left_match_tab [gaplen+1]; // gaplen = gap_length_on_read. NB: one-base gap has 2 potential spliting points
	char right_match_tab [gaplen+1];
	char split_can_be_here [gaplen+1]; // A split can be here, before X1? determined by AT/GC.
	char left_split_chro_base_tab [gaplen+2]; // left_last_matched_base_loc - 1 (inc) ~ ??
	char right_split_chro_base_tab [gaplen+2]; // ?? ~ right_last_matched_base_loc + 1 (inc)
	int x1;

	gene_value_index_t * current_value_index = cct_context->value_index;


	gvindex_get_range(current_value_index, 
                  left_last_matched_base_loc + 1, 
                  left_split_chro_base_tab, 
                  gaplen + 2);

	gvindex_get_range(current_value_index, 
                  right_last_matched_base_loc - gaplen - 2, 
                  right_split_chro_base_tab, 
                  gaplen + 2);

	if(0)for(x1=0; x1<gaplen+2;x1++){
		char left_charg = gvindex_get(current_value_index,left_last_matched_base_loc + x1 +1); 
		char right_charg = gvindex_get(current_value_index,right_last_matched_base_loc - gaplen + x1 -2); 

		left_split_chro_base_tab[x1] = left_charg;
		right_split_chro_base_tab[x1] = right_charg;
	}

	#define L_QBASE_1 left_split_chro_base_tab[ x1 ] 
	#define L_QBASE_2 left_split_chro_base_tab[ x1 +1] 

	#define R_QBASE_1 right_split_chro_base_tab[ x1 ] 
	#define R_QBASE_2 right_split_chro_base_tab[ x1 +1 ] 
	for(x1=0; x1<=gaplen;x1++){ // splice BEFORE x1?
		if(      L_QBASE_1 == 'G' && L_QBASE_2 == 'T' &&  R_QBASE_1 == 'A' && R_QBASE_2 == 'G' ) split_can_be_here[x1]=1;
		else if( L_QBASE_1 == 'C' && L_QBASE_2 == 'T' &&  R_QBASE_1 == 'A' && R_QBASE_2 == 'C' ) split_can_be_here[x1]=2;
		else split_can_be_here[x1]=0;
	}

	for(x1=0; x1<=gaplen;x1++){
		char rchar = read_text[x1+ gap_start_in_read];
		char leftchar_if_split_before_x1 = 0;
		char rightchar_if_split_before_x1 = 0;

		leftchar_if_split_before_x1  = left_split_chro_base_tab[x1]; 
		if(x1<gaplen)rightchar_if_split_before_x1 = right_split_chro_base_tab[x1 + 2]; // x1==gaplen: value is meaningless.
		left_match_tab [x1] =(rchar == leftchar_if_split_before_x1);
		right_match_tab [x1] =(rchar == rightchar_if_split_before_x1);
	}

	int best_matched = -1;
	int best_gap_loc = -1;
	for(x1=0; x1<=gaplen;x1++){
		if(!split_can_be_here[x1])continue;

		int my_matched = 0;
		int x2;
		for(x2 = 0; x2 < gaplen; x2++) my_matched += (x2 >= x1)? right_match_tab[ x2 ] : left_match_tab[ x2 ] ; 
		if(my_matched > best_matched){
			best_gap_loc = x1;
			best_matched = my_matched;
			(*is_GTAG_negative) = split_can_be_here[x1]==2;
		}
	}
	*total_misma = (gaplen - best_matched);
	*outgaplen = gaplen;
	return best_gap_loc;
}

int cellCounts_add_covered_indels_in_table_getval(unsigned int pos, void * context, char * sequence_space, int num_bases){
	cellcounts_global_t * cct_context = context;
	gene_value_index_t * current_value_index = cct_context->value_index;
	gvindex_get_range(current_value_index,pos, sequence_space, num_bases);
	return 0;
}

void cellCounts_add_covered_indels_in_table(cellcounts_global_t * cct_context, int thread_no, int sample_i, char * read, int cov_start, int cov_end, int expected_indel, unsigned int chro_loc, char * read_name){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	char *** dpbuf = thread_context -> dynamic_align_buffers;
	int * dpplt = thread_context -> dynamic_align_penalties;
	char result_buffer[3*MAX_SCRNA_READ_LENGTH+1];
	char result_cigar[3*MAX_SCRNA_READ_LENGTH+1];
	int moved_len = general_dynamic_align(read + cov_start, cov_end - cov_start, chro_loc + cov_start, result_buffer, expected_indel, cct_context ->max_indel_length, (void***)dpbuf, dpplt, cellCounts_add_covered_indels_in_table_getval, cct_context);
	int retlen = general_dynamic_align_moves_to_cigar(result_buffer, moved_len, result_cigar);
//if(strstr(read_name,"R00000000207")) fprintf(stderr,"DYNAMIC_INDEL for %s : %s\n", read_name, result_cigar);

	int nch, tmpi=0,x1;
	unsigned int chro_cursor = chro_loc + cov_start;
	for(x1=0; 0!=(nch=result_cigar[x1]) ; x1++){
		if(isdigit(nch)){
			tmpi = 10*tmpi+nch-'0';
		}else{
			if(nch=='M'){
				chro_cursor += tmpi;
			}else if(nch=='D' || nch == 'I'){
				char * chro_name = NULL;
				int chro_pos = 0, inslen_negative = (nch == 'I')? - tmpi:0;

				locate_gene_position(chro_cursor, &cct_context -> chromosome_table, &chro_name, &chro_pos);
				int rposv = chro_pos;
				if( nch=='D') rposv = chro_pos + tmpi;
				cellCounts_add_or_update_chroEvent_in_table(cct_context, sample_i, chroEvent_t_TYPE_INDEL , chro_name, chro_pos -1, rposv, inslen_negative, 0); // thread safe
				if(nch=='D') chro_cursor += tmpi;
			}
			tmpi=0;
		}
	}
}

int cellCounts_call_juncs_put_in_tab(cellcounts_global_t * cct_context, int thread_no, int sample_i, gene_sc_vote_t * votetab, char * read_name, char * read_text, char * read_bin, char * read_qual, int read_len, cellcounts_vote_number_t all_subreads) {
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	int i,j,reverse_text_offset, distinct_vote_number_i;
	int tstpos=0, mapos=0;
	ArrayList * jstub_potential_mainhalf_list = ArrayListCreate(200);

	if(votetab -> max_vote >= cct_context -> min_votes_per_mapped_read){
		HashTable * iijj_to_chro_ptr = HashTableCreate(1000);

		for (i=0; i<GENE_SCRNA_VOTE_TABLE_SIZE; i++){
			for (j=0; j< votetab->items[i]; j++){
				int vv = votetab -> votes[i][j];
				char * tstchro=NULL;
				locate_gene_position(votetab -> pos[i][j], & cct_context -> chromosome_table, & tstchro, & tstpos);
				srInt_64 iijj = i*1LLU << 24 | j;
				HashTablePut(iijj_to_chro_ptr , NULL+1+iijj, tstchro);
				
				if(vv >= JUNCTION_MINIMUM_MAIN_HALF){
					ArrayListPush( jstub_potential_mainhalf_list,NULL+iijj );
				}
			}
		}

		HashTable * indel_dp_exed = HashTableCreate(100);
		int mainhalf_i;
		for(mainhalf_i = 0; mainhalf_i < jstub_potential_mainhalf_list -> numOfElements; mainhalf_i++){
			int best_s1_aft_split = -1;
			int best_misma_no = -1;
			srInt_64 best_mate_score = -1;
			srInt_64 best_mate_iijj = -1;
			int best_GTAG_negative = -1;

			srInt_64 maiijj = ArrayListGet( jstub_potential_mainhalf_list , mainhalf_i ) - NULL;
			int mai = maiijj >> 24, maj = maiijj &0xffffffllu;
			char * machro=NULL;

			unsigned int maloc = votetab -> pos[mai][maj];
			locate_gene_position(maloc, & cct_context -> chromosome_table, & machro, & mapos);

			int manegative = votetab -> masks[mai][maj];
			int ma_cov_start = votetab -> coverage_start[mai][maj] + JUNCTION_WIDDEN_GAP_LEN;
			int ma_cov_end = votetab -> coverage_end[mai][maj] - JUNCTION_WIDDEN_GAP_LEN;
			int ma_toli = votetab -> toli[mai][maj];
			int mavotes = votetab -> votes[mai][maj];
			#define read_text_rev_DEFLEN MAX_SCRNA_READ_LENGTH+1
			char * read_maped_look = manegative?read_text+read_text_rev_DEFLEN:read_text;
			int ma_indels_in_coverage = votetab -> indel_recorder[mai][maj][ma_toli -3 +2];
//			int ma_indel              = votetab -> indel_recorder[mai][maj][ma_toli -3 +2];

			for (i=0; i<GENE_SCRNA_VOTE_TABLE_SIZE; i++){
				for (j=0; j< votetab->items[i]; j++){
					unsigned int tstloc = votetab -> pos[i][j];

					srInt_64 chro_dist = maloc;
					chro_dist -= tstloc;

					if(abs(chro_dist) > JUNCTION_MAX_CHRO_DISTANCE)continue;
					if(abs(chro_dist) <= cct_context -> max_indel_length)continue;  

					srInt_64 iijj = i*1LLU << 24 | j;
					char * tstchro = HashTableGet(iijj_to_chro_ptr, NULL+1+iijj);
					if(tstchro!=machro)continue ; // all chro names are the PTR in the chro table, hence comparible.

					int tstnegative = votetab -> masks[i][j];
					if(tstnegative!=manegative) continue;

					int tst_cov_start = votetab -> coverage_start[i][j] + JUNCTION_WIDDEN_GAP_LEN;
					int is_correct_chro_order = (tst_cov_start > ma_cov_start)==(tstloc > maloc);
					if(!is_correct_chro_order)continue;

					int tstvotes = votetab -> votes[i][j];
					if( tstvotes > mavotes || ( mavotes == tstvotes && tst_cov_start >=ma_cov_start ) ) continue ; // keep order of comparison

					int tst_cov_end = votetab -> coverage_end[i][j] - JUNCTION_WIDDEN_GAP_LEN;


					int is_overlap_in_read = ( tst_cov_start <= ma_cov_start && tst_cov_end > ma_cov_start ) ||  ( ma_cov_start <= tst_cov_start && ma_cov_end > tst_cov_start ) ;
					if(is_overlap_in_read)continue;

					int tst_toli = votetab -> toli[i][j];
					int tst_indel = votetab -> indel_recorder[i][j][tst_toli -3 +2];
					int dist_log2 = integer_log2_64(abs(chro_dist));

					int votes_multiplex = tstvotes * mavotes;
					int misma_in_GTAG_met = -1, gaplen=-1, is_negative_by_GTAG=-1;

					int offset_from_gap = cellCounts_junc_meet_in_the_middle(cct_context, thread_context, read_maped_look,
						ma_cov_start, ma_cov_end, ma_indels_in_coverage, maloc, tst_cov_start, tst_cov_end, tst_indel, tstloc, &misma_in_GTAG_met, &gaplen, &is_negative_by_GTAG);

					srInt_64 this_mate_score = -1llu;
					if(  offset_from_gap >=0 && misma_in_GTAG_met <= JUNCTION_MAX_MISMA_MEET)
						this_mate_score = (MAX_READ_LENGTH - misma_in_GTAG_met)* 30000llu * 10000llu +
							30000llu * votes_multiplex +
							( 30-dist_log2 ) * 1000llu;

					if(this_mate_score >0 && this_mate_score > best_mate_score){
						best_mate_score = this_mate_score;
						best_mate_iijj = (mainhalf_i*1LLU << 48)|iijj;
						best_s1_aft_split = offset_from_gap;
						best_misma_no = misma_in_GTAG_met;
						best_GTAG_negative = is_negative_by_GTAG;
					}
				}
			}
//			#warning "====  Think: could a read contribute multiple junctions? could a major half contribute multiple junctions? ===="
			if( best_s1_aft_split >=0 ){
				srInt_64 best_major_half = best_mate_iijj>>48;
				srInt_64 tstiijj= best_mate_iijj & 0xffffffffffffllu;

				if(best_major_half != mainhalf_i)SUBREADprintf("LOGIC HAS CHANGED HERE!!\n");

	//			srInt_64 maiijj= ArrayListGet(jstub_potential_mainhalf_list, best_major_half)-NULL;
	//			int mai = maiijj>>24, maj=maiijj & 0xffffff;

				int tsti = tstiijj>>24, tstj=tstiijj & 0xffffff;

				//unsigned int maloc = votetab -> pos[mai][maj];
				//int manegative = votetab -> masks[mai][maj];
				//int ma_cov_start = votetab -> coverage_start[mai][maj] + JUNCTION_WIDDEN_GAP_LEN;
				//int ma_cov_end = votetab -> coverage_end[mai][maj] - JUNCTION_WIDDEN_GAP_LEN;
				//int ma_toli = votetab -> toli[mai][maj];

				unsigned int tstloc = votetab -> pos[tsti][tstj];
				int tst_cov_start = votetab -> coverage_start[tsti][tstj] + JUNCTION_WIDDEN_GAP_LEN;
				int tst_cov_end = votetab -> coverage_end[tsti][tstj] - JUNCTION_WIDDEN_GAP_LEN;
				int tst_toli = votetab -> toli[tsti][tstj];
				int minor_indels_in_coverage = votetab -> indel_recorder[tsti][tstj][tst_toli -3 +2];

				unsigned int junction_left_last_exon_base, junction_right_first_exon_base;
				int rgap = -1;
				if(tst_cov_start < ma_cov_start) rgap = ma_cov_start - tst_cov_end; else  rgap = tst_cov_start - ma_cov_end;

				if(ma_cov_start > tst_cov_start){
					junction_left_last_exon_base = tstloc + minor_indels_in_coverage+ tst_cov_end + best_s1_aft_split -1;
					junction_right_first_exon_base = maloc+ tst_cov_end + best_s1_aft_split ;
				}else{
					junction_left_last_exon_base = maloc + ma_indels_in_coverage+ ma_cov_end + best_s1_aft_split -1;
					junction_right_first_exon_base = tstloc+ ma_cov_end + best_s1_aft_split ;
				}
				locate_gene_position(junction_left_last_exon_base, & cct_context -> chromosome_table, & machro, & mapos);
				locate_gene_position(junction_right_first_exon_base, & cct_context -> chromosome_table, & machro, & tstpos);
				int is_update = cellCounts_add_or_update_chroEvent_in_table(cct_context, sample_i, chroEvent_t_TYPE_JUNCTION, machro, mapos, tstpos, 0, 0);//mapos and tstpos are just borrowed variables. 

				if(NULL==HashTableGet(indel_dp_exed, NULL+1+(tstiijj & 0xffffffffffffllu))){
					HashTablePut(indel_dp_exed, NULL+1+(tstiijj & 0xffffffffffffllu), NULL+1);
					if(minor_indels_in_coverage)cellCounts_add_covered_indels_in_table(cct_context, thread_no, sample_i, read_maped_look, tst_cov_start, tst_cov_end, minor_indels_in_coverage, tstloc, read_name);
				}
			}
			if(NULL==HashTableGet(indel_dp_exed, NULL+1+(maiijj & 0xffffffffffffllu))){
				HashTablePut(indel_dp_exed, NULL+1+(maiijj & 0xffffffffffffllu), NULL+1);
				if(ma_indels_in_coverage)cellCounts_add_covered_indels_in_table(cct_context, thread_no, sample_i, read_maped_look, ma_cov_start, ma_cov_end, ma_indels_in_coverage, maloc, read_name);
			}
		}
		HashTableDestroy(indel_dp_exed);
		HashTableDestroy(iijj_to_chro_ptr);
	}
	ArrayListDestroy(jstub_potential_mainhalf_list);
	return 0;
}

#define ADD_USED_EVENT_DETAILS_PTR(pptr, ppx2) {\
                                added_event_array[added_event_array_no] = (pptr);  \
                                if((pptr) -> n_events >1)added_event_insertion_index[added_event_array_no] = (ppx2); else added_event_insertion_index[added_event_array_no] = -2; \
                                added_event_array_no ++; \
                                if(added_event_array_no >= JUNCTION_MAX_COLOCATION)SUBREADprintf("ERROR: the read overlaps with too many chromosome events.\n");\
				}


#define FIND_USED_EVENT_DETAILS(toadd_ptr, toadd_x2, is_used, x2_used) { int xed3; is_used=0; x2_used=0; \
						for(xed3=0; xed3<added_event_array_no; xed3++){\
							chroEvent_t * added_env = added_event_array[xed3];\
							if((toadd_ptr) == added_env){\
								is_used = 1; \
								if((toadd_ptr)-> n_events>1 && toadd_x2 == added_event_insertion_index[xed3]) x2_used = 1; \
							} }}

void cellCounts_add_supported_unsupported_reads_from_cigar( cellcounts_global_t * cct_context, int thread_no, int sample_i, int reporting_index){
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;

//#warning "===== BAD DEBUG ====="
//if(strstr(thread_context -> realignment_event_read_name, "ATTTTGATT|G;G,GG;G;;;GGGGG;,," ))return;
	char * cigar = thread_context -> reporting_cigars[reporting_index];

	// The added_event_array stores the events that are supported and/or non-supported.
	// A read can only give a count once to an event. 
	chroEvent_t * added_event_array[ JUNCTION_MAX_COLOCATION ];
	int added_event_insertion_index[ JUNCTION_MAX_COLOCATION ];
	int added_event_array_no = 0;

	unsigned int linear_cur = thread_context -> reporting_positions[reporting_index], linear0 = linear_cur;
	//int read_cur = 0;
	int nch, tmpi=0, x1=0, b1off=0, x3;
	char *ch = NULL;
	locate_gene_position(linear_cur, &cct_context -> chromosome_table, &  ch, &b1off);
	char chro_strn_ky[MAX_CHROMOSOME_NAME_LEN+10];
	char negchar = '*', LRi;
	snprintf(chro_strn_ky, MAX_CHROMOSOME_NAME_LEN+10, "%s\t%c", ch, negchar);
	IVT_IntervalTreeNode **LR_roots = HashTableGet(cct_context -> chroEvent_entry_table[sample_i], chro_strn_ky); // LR_roots [0] : left-edge ; LR_roots [1] : right-edge ; LR_roots [2] : covered_regions

	int eventbufsize = JUNCTION_MAX_COLOCATION, founditems = 0;
	IVT_Interval * eventbuf[eventbufsize];

	while(0!=(nch = cigar[x1++])){
		if(isdigit(nch)){
			tmpi = 10*tmpi+nch-'0';
		}else{
			unsigned int p2 = 0;
			chroEvent_t * oneitem = NULL;

			// add supp reads to events that EXACTLY match the two edge locations.
			if(nch == 'N' || nch == 'D') p2 = linear_cur + tmpi;
			else if(nch=='I') p2 = linear_cur;

			// Add junction or indel supp reads.
			// N.B. each junction is exactly the same as an intron in the transcript.
			// # supp for the junction is exactly the same as # non-supp for the intron.
			// # non-supp for the junction is exactly the same as # supp for the intron.
			if(p2){
				int last_added_insertion_index = -1;
				srUInt_64 envkey = ((linear_cur-1) *1LLU<<32) | p2;
				oneitem = HashTableGet( cct_context -> chroEvent_detail_table[sample_i], NULL+ envkey);
				if(!oneitem) SUBREADprintf("UNABLE to find a junction: %s : %d ~ %d ; CIGAR FOR %s is %s mapped to %s:%d --  CHR %c    ItemPTR=%p    ky=%p\n", chro_strn_ky, ( linear_cur-1 ), p2, thread_context -> realignment_event_read_name , cigar, ch, b1off, nch, oneitem, NULL+ envkey);
				cellCounts_chroEvent_locks_opt(cct_context, thread_no, oneitem, 1);
				if(oneitem -> n_events > 1){
					int x2;
					// You wont have two D/N events having exactly same splicing points. Only insertions with diff lengths can do so.
					assert(nch=='I');
					for(x2=0; x2< oneitem -> n_events;x2++)if( oneitem -> length[x2] == -tmpi ){
						last_added_insertion_index = x2;
						oneitem -> step2_supported_reads[ x2 ]++;
						break;
					}
				}else oneitem -> step2_supported_reads = (void*) oneitem -> step2_supported_reads+1;
				cellCounts_chroEvent_locks_opt(cct_context, thread_no, oneitem, 0);

				int on_chro_len = tmpi +1;
				if(nch=='I') on_chro_len = 1;
				ADD_USED_EVENT_DETAILS_PTR(oneitem, last_added_insertion_index);
			}

			if(nch=='M'){
				// add overlapping exons.
				if(LR_roots){
					int Mstart_base = linear_cur - linear0 + b1off, x2;
					int Mlast_base = Mstart_base + tmpi - 1;
					founditems=0;

					IVT_IntervalTreeNode * cover_root = LR_roots[2];
					IVT_query_range(cover_root , Mstart_base, Mlast_base, eventbuf, eventbufsize, &founditems);

					for(x2=0; x2 < founditems; x2++){
						IVT_Interval * evb = eventbuf [x2];
						int evbposleft = evb -> start;
						int evbposright = evb -> end;
						unsigned int linear_env_1 = linear_gene_position(&cct_context->chromosome_table, ch , evbposleft);
						unsigned int linear_env_2 = linear_gene_position(&cct_context->chromosome_table, ch , evbposright);
						unsigned int linear_env_L = min(linear_env_1, linear_env_2);
						unsigned int linear_env_R = max(linear_env_1, linear_env_2);

						srUInt_64 envkey = (linear_env_L*1LLU<<32)| linear_env_R;
						chroEvent_t * exed = HashTableGet(cct_context -> chroEvent_detail_table [sample_i], NULL+envkey);


						if(!exed) continue; // if it isn't an Exon event, it can be NULL. We only need EXONs.

						int is_the_last, x2_not_matter;
						if(exed -> event_type != chroEvent_t_TYPE_EXON) continue;

						FIND_USED_EVENT_DETAILS(exed, -1, is_the_last, x2_not_matter);
						if(is_the_last) continue;

						cellCounts_chroEvent_locks_opt(cct_context, thread_no, exed, 1);
						exed -> step2_supported_reads = (void*) exed -> step2_supported_reads+1;
						cellCounts_chroEvent_locks_opt(cct_context, thread_no, exed, 0);

						ADD_USED_EVENT_DETAILS_PTR(exed, -1); 
					}
				}
			}

			if(nch == 'D' || nch == 'S' /* mapping location of 'S' is shiftted right before writing */ || nch == 'M' || nch== 'N') linear_cur += tmpi; 
			tmpi = 0;
		}
	}

	x1=0;
	tmpi=0;
	linear_cur = linear0;
	while(0!=(nch = cigar[x1++])){
		if(isdigit(nch)){
			tmpi = 10*tmpi+nch-'0';
		}else{
			if(nch == 'M' || nch == 'N'){
				// Find all events in this region -- they aren't supported.
				// The events can be starting and/or ending in this region.
				// This includes the events at two edges. 
				if(LR_roots) for(LRi=0; LRi<2; LRi++){
					founditems=0;
					IVT_IntervalTreeNode * myroot = LR_roots[LRi];
					int Mstart_base = linear_cur - linear0 + b1off;
					int Mlast_base = Mstart_base + tmpi - 1;
					IVT_query_range(myroot , Mstart_base, Mlast_base, eventbuf, eventbufsize, &founditems);
					int evbi;
					for(evbi=0; evbi<founditems; evbi++){
						IVT_Interval * evb = eventbuf [evbi];
						int evbposleft = evb -> start;
						int evbposright = evb -> attr - NULL;
						unsigned int linear_env_1 = linear_gene_position(&cct_context->chromosome_table, ch , evbposleft);
						unsigned int linear_env_2 = linear_gene_position(&cct_context->chromosome_table, ch , evbposright);
						unsigned int linear_env_L = min(linear_env_1, linear_env_2);
						unsigned int linear_env_R = max(linear_env_1, linear_env_2);

						srUInt_64 envkey = (linear_env_L*1LLU<<32)| linear_env_R;
						chroEvent_t * uned = HashTableGet(cct_context -> chroEvent_detail_table[sample_i] , NULL+envkey);

						if(nch == 'N' && uned -> event_type != chroEvent_t_TYPE_EXON)continue; // "N" can only give non-supporting to the fly-over exons.
						if(nch == 'N' && !(uned -> left_edge >= linear_cur && uned -> left_edge + ((void*) uned -> length-NULL) <= linear_env_R )) continue; // fly-over exons must be entirely in the "N" part.
						int is_the_last = 0, is_same_x2 = 0;

						cellCounts_chroEvent_locks_opt(cct_context, thread_no, uned, 1);
						if(uned -> n_events == 1){
							FIND_USED_EVENT_DETAILS(uned, -1, is_the_last, is_same_x2);
							if(!is_the_last) uned -> step2_non_supported_reads = (void *) uned -> step2_non_supported_reads +1 ;
						}else{
							int x2;
							for(x2=0; x2<uned -> n_events; x2++){
								FIND_USED_EVENT_DETAILS(uned, x2, is_the_last, is_same_x2);
								if(!(is_the_last && is_same_x2)){
									ADD_USED_EVENT_DETAILS_PTR(uned, x2); 
									uned -> step2_non_supported_reads[x2]++;
								}
							}
						}
						cellCounts_chroEvent_locks_opt(cct_context, thread_no, uned, 0);
					}
				}
			}

			if(nch == 'D' || nch == 'M' || nch== 'N') linear_cur += tmpi; 
			tmpi = 0;
		}
	}
}

int cellCounts_select_and_write_alignments(cellcounts_global_t * cct_context, int thread_no, int sample_i, gene_sc_vote_t * votetab, char * read_name, char * read_text, char * read_bin, char * read_qual, int read_len, cellcounts_vote_number_t all_subreads) {
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	int i,j,reverse_text_offset, distinct_vote_number_i;

	int index_gap_width = cct_context -> current_index -> index_gap;

	thread_context -> total_voteIJs_to_write = 0;
	thread_context -> populating_voteIJ_buf_index=0;
	if(votetab && votetab -> max_vote >= cct_context -> min_votes_per_mapped_read){
		int top_distinct_vote_numbers[cct_context -> max_distinct_top_vote_numbers];
		memset(top_distinct_vote_numbers, 0 , cct_context -> max_distinct_top_vote_numbers * sizeof(int));

		for (i=0; i<GENE_SCRNA_VOTE_TABLE_SIZE; i++){
			for (j=0; j< votetab->items[i]; j++){
				int vv = votetab -> votes[i][j];
				if(vv>=cct_context -> min_votes_per_mapped_read)cellCounts_update_top_three(cct_context, top_distinct_vote_numbers, vv);
			}
		}

		thread_context -> alignment_repating_table = HashTableCreate(50);
		int doneit=0;
		for(distinct_vote_number_i = 0 ; distinct_vote_number_i < cct_context -> max_distinct_top_vote_numbers; distinct_vote_number_i ++){
			if(doneit||thread_context -> populating_voteIJ_buf_index >= cct_context -> max_candidate_voteIJ_per_read)break;
			int this_vote_N = top_distinct_vote_numbers[distinct_vote_number_i];

			if(this_vote_N < 1 || (top_distinct_vote_numbers[0] - this_vote_N > cct_context -> max_differential_from_top_vote_number )) break;

			for (i=0; i<GENE_SCRNA_VOTE_TABLE_SIZE; i++){
				if(doneit||thread_context -> populating_voteIJ_buf_index >= cct_context -> max_candidate_voteIJ_per_read)break;
				for (j=0; j< votetab->items[i]; j++){
					if(doneit||thread_context -> populating_voteIJ_buf_index >= cct_context -> max_candidate_voteIJ_per_read)break;

					int vv = votetab->votes[i][j];
					if(vv == this_vote_N){
						if(cct_context -> do_cell_level_junction_detection && sample_i >0){
							int CR15GLS = (read_len - 15 - index_gap_width)<<16;
							int subread_step =  CR15GLS /(cct_context -> total_subreads_per_read -1);
							if(subread_step<(index_gap_width<<16))subread_step = index_gap_width<<16;

							int perfect_align_srno = votetab->indel_recorder[i][j][0]-1; // indel_recorder is subread_no + 1
							int perfect_coved_firstbase = ((subread_step * perfect_align_srno) >> 16);
							perfect_align_srno = votetab->indel_recorder[i][j][1]-1;
							int perfect_coved_lastbase = ((subread_step * perfect_align_srno) >> 16) +15;

							int is_negative = votetab->masks[i][j];
							int reverse_text_offset = is_negative?MAX_SCRNA_READ_LENGTH+1:0;

							cellCounts_explain_one_alignment(cct_context, thread_no, sample_i, read_name, read_text + reverse_text_offset, read_len, perfect_coved_firstbase, perfect_coved_lastbase, votetab->pos[i][j], is_negative, vv);
						}else cellCounts_explain_PaperVersion_one_alignment (cct_context, thread_no, read_name, read_bin, read_text, read_len, all_subreads, votetab, i, j); // ( cct_context, thread_no, read_name, read_bin, read_text + reverse_text_offset, read_len, perfect_coved_firstbase, perfect_coved_lastbase, votetab,i,j);
					}
				}
			}
		}
		HashTableDestroy(thread_context -> alignment_repating_table);
		thread_context -> total_voteIJs_to_write = min(thread_context -> total_voteIJs_to_write, cct_context -> max_reported_alignments_per_read);
	}else thread_context -> total_voteIJs_to_write = 0;

//	SUBREADprintf("READQV %s : VOTE %d , ALN %d\n\n", read_name, votetab->max_vote, thread_context -> total_voteIJs_to_write);


	if(thread_context -> total_voteIJs_to_write) {
		int sorting_index [thread_context -> populating_voteIJ_buf_index];
		for(distinct_vote_number_i = 0 ; distinct_vote_number_i < thread_context -> populating_voteIJ_buf_index; distinct_vote_number_i ++) sorting_index [distinct_vote_number_i ] = distinct_vote_number_i ;
		void * sorting_ptr[2];
		sorting_ptr[0] = thread_context;
		sorting_ptr[1] = sorting_index;
		//sort : large number first
		quick_sort(sorting_ptr , thread_context -> populating_voteIJ_buf_index , sort_readscore_compare_LargeFirst, sort_readscore_exchange); // The last many records are 0-score records. Only "total_voteIJs_to_write" records ahead are worth writting (score > 0).

		for(thread_context -> writing_voteID_buf_index = 0 ; thread_context -> writing_voteID_buf_index < thread_context -> total_voteIJs_to_write; thread_context -> writing_voteID_buf_index ++){
			int myno = sorting_index[thread_context -> writing_voteID_buf_index ];
			if(thread_context -> reporting_scores[ myno ] < 1)continue;
			if(thread_context -> writing_voteID_buf_index >= cct_context -> max_reported_alignments_per_read) break;
			reverse_text_offset = (thread_context -> reporting_flags[myno] & SAM_FLAG_REVERSE_STRAND_MATCHED)?MAX_SCRNA_READ_LENGTH+1:0;
			if(reverse_text_offset >0 && 0==read_qual[reverse_text_offset]){
				strcpy(read_qual+reverse_text_offset, read_qual);
				reverse_quality(read_qual+reverse_text_offset, read_len);
			}
#warning "===== CURRENT JUNCTION DETECTION DOESN'T NEED ANTI_SUPPORT OUT ====="
			if(0)if(cct_context -> do_cell_level_junction_detection && sample_i>0)cellCounts_add_supported_unsupported_reads_from_cigar( cct_context, thread_no, sample_i, myno);
			cellCounts_write_read_in_batch_bin(cct_context, thread_no, sample_i, myno, read_name, read_text + reverse_text_offset, read_qual+reverse_text_offset, read_text , read_qual , read_len);
		}
	} else cellCounts_write_read_in_batch_bin(cct_context, thread_no, sample_i, -1, read_name, read_text, read_qual , read_text, read_qual, read_len);

	return 0;
}

int cellCounts_genekey2int(char *key) {
	int ret = 0;

	char * keyc = key +16;
	for (; key < keyc; key++) {
		char c1 = *key;
		ret = (ret << 2) | base2int(c1);
	}
	return ret;
}


int cellCounts_process_copy_ptrs_to_votes_compare(void * arrp, int i, int j){
	int * Xptrs = ((temp_votes_per_read_t**)arrp)[1] -> votes;
	int * trying_subread_no = ((int**)arrp)[0];
	int diffv = Xptrs[ trying_subread_no [i] ] - Xptrs[ trying_subread_no [j] ] ;
	return diffv;
}

void cellCounts_process_copy_ptrs_to_votes_exchange(void * arrp, int i, int j){
	int * trying_subread_no = ((int**)arrp)[0];
	int tmpi = trying_subread_no [i];
	trying_subread_no [i] = trying_subread_no [j];
	trying_subread_no [j] = tmpi;
}


#define INDEL_SEGMENT_SIZE 5
#define VOTING_PRIME_NUMBER 66889

//#define _index_vote(key) (((unsigned int)(key))%GENE_SCRNA_VOTE_TABLE_SIZE)
#define _index_vote_tol(key) (((unsigned int)(key)/INDEL_SEGMENT_SIZE)%GENE_SCRNA_VOTE_TABLE_SIZE)
#define cellCounts_voting_update_topK(vt, vn, ij)  { \
	int x3, x3done=0;\
	for(x3= 0; x3 < cct_context -> max_candidate_voteIJ_per_read +1; x3++){\
		if(ij == vt->topK_IJ[x3]){\
			if(x3 > 0 && vn > vt->topK_votes[x3-1]){\
				int x4;\
				for(x4 = x3; x4< cct_context -> max_candidate_voteIJ_per_read; x4++){\
					vt->topK_votes[x4] = vt->topK_votes[x4+1];\
					vt->topK_IJ[x4] = vt->topK_IJ[x4+1];\
				}\
				vt->topK_votes[cct_context -> max_candidate_voteIJ_per_read]=0;\
			}else{\
				vt->topK_votes[x3] = vn;\
				x3done=1;\
			}\
			break;\
		}\
	}\
	if(!x3done)for(x3= 0; x3 < cct_context -> max_candidate_voteIJ_per_read +1; x3++){\
		if(vn > vt->topK_votes[x3]){\
			int x4;\
			if(vt->topK_votes[x3]>0)for(x4 = cct_context -> max_candidate_voteIJ_per_read ; x4 > x3; x4--){\
				vt->topK_votes[x4] = vt->topK_votes[x4-1];\
				vt->topK_IJ[x4] = vt->topK_IJ[x4-1];\
			}\
			vt->topK_votes[x3]=vn;\
			vt->topK_IJ[x3]=ij;\
			break; \
		}\
	}\
}


int fix_indel_record_order_compare(void * arr, int i, int j){
	cellcounts_vote_number_t * rec = arr;
	return rec[i*3]-rec[j*3];
}

void fix_indel_record_order_exchange(void * arr, int i, int j){
	cellcounts_vote_number_t * rec = arr;
	cellcounts_vote_number_t ttv = rec[i*3];
	rec[i*3] = rec[j*3];
	rec[j*3] = ttv;

	ttv = rec[i*3+1];
        rec[i*3+1] = rec[j*3+1];
        rec[j*3+1] = ttv;

	ttv = rec[i*3+2];
        rec[i*3+2] = rec[j*3+2];
        rec[j*3+2] = ttv;
}


void cellCounts_process_copy_ptrs_to_votes(cellcounts_global_t * cct_context, int thread_no, temp_votes_per_read_t * ptrs, gene_sc_vote_t * vote, int applied_subreads_per_strand, char * read_name){
	int subreads = applied_subreads_per_strand*2;
	int x1, x2, trying_subread_no[subreads];
	for(x1=0; x1<subreads; x1++) trying_subread_no[x1]=x1;
	void * sort_arr[2];
	sort_arr [0] = trying_subread_no;
	sort_arr [1] = ptrs;
	quick_sort(sort_arr , subreads , cellCounts_process_copy_ptrs_to_votes_compare, cellCounts_process_copy_ptrs_to_votes_exchange);

	init_gene_vote(vote);
	int cct_indel_len = cct_context -> max_indel_length, cct_indel_neg = -cct_context -> max_indel_length;
	for(x1=0; x1<subreads; x1++){
		int myno = trying_subread_no[x1];
		int has_votes = ptrs->votes[myno];
		if(has_votes<1) continue;

		int mynoP1PStr=(myno%applied_subreads_per_strand)+1;
		int offset = ptrs->offsets[myno];
		int of_p_16 = offset + 16;
		int is_reversed = (myno >= applied_subreads_per_strand)?IS_NEGATIVE_STRAND:0;
		unsigned int * index_ptr = ptrs->start_location_in_index [myno];
		int vote_prime_sum = ptrs->votes[trying_subread_no[subreads-1]];

		for(x2 = 0 ; x2 < has_votes; x2++){
			unsigned int kv = index_ptr[ vote_prime_sum % has_votes ] - offset;
			vote_prime_sum += VOTING_PRIME_NUMBER;
			int iix, offsetX2, offsetX, datalen, datalen2, found = 0;
			offsetX = offsetX2 = _index_vote_tol(kv);
			datalen = datalen2 = vote -> items[offsetX];
			unsigned int * dat2, *dat;
			dat = dat2 = vote -> pos[offsetX];

			for(iix = 0; iix<=INDEL_SEGMENT_SIZE; iix = iix>0?-iix:(-iix+INDEL_SEGMENT_SIZE)){
				if(iix) {
					offsetX = _index_vote_tol(kv+iix);
					datalen = vote -> items[offsetX];
					if(!datalen)continue;
					dat = vote -> pos[offsetX];
				} else if(!datalen) continue;

				int itemidx;
				for (itemidx=0;itemidx<datalen;itemidx++){
					signed int dist0 = kv-dat[itemidx];
					if( dist0 >= cct_indel_neg  && dist0 <= cct_indel_len  && is_reversed == vote->masks[offsetX][itemidx]){
						int toli, tolimax=vote -> toli[offsetX][itemidx], known_indel=0;

						cellcounts_vote_number_t * indelrec = vote -> indel_recorder[offsetX][itemidx];
						for(toli = 0; toli < tolimax; toli +=3){
							if( indelrec [toli+2] == dist0 ){
								if(indelrec [toli] > mynoP1PStr) indelrec [toli]  = mynoP1PStr;
								if(indelrec [toli +1] < mynoP1PStr) indelrec [toli +1]  = mynoP1PStr;
								known_indel=1;
								break;
							}
						}

						if(tolimax < MAX_INDEL_TOLERANCE*3 && !known_indel){
							indelrec [tolimax] = mynoP1PStr;
							indelrec [tolimax+1] = mynoP1PStr;
							indelrec [tolimax+2] = dist0;
							vote -> toli[offsetX][itemidx]=tolimax+3;
						}

						cellcounts_vote_number_t test_max = (vote->votes[offsetX][itemidx]);
						test_max ++;
						vote -> votes[offsetX][itemidx] = test_max;
						if(vote->max_vote < test_max) vote->max_vote = test_max;

						if (offset < vote->coverage_start [offsetX][itemidx])
							vote->coverage_start [offsetX][itemidx] = offset;
						if (of_p_16 > vote->coverage_end [offsetX][itemidx])
							vote->coverage_end [offsetX][itemidx] = of_p_16;

						found = 1;
						break;
					}
				}
				if(found)break;
			}
			if((!found) && datalen2 < GENE_SCRNA_VOTE_SPACE){
				vote -> items[offsetX2] = datalen2+1;
				dat2[datalen2] = kv;
				vote -> masks[offsetX2][datalen2] = is_reversed;
				vote -> votes[offsetX2][datalen2] = 1;
				vote -> indel_recorder[offsetX2][datalen2][0] = mynoP1PStr;
				vote -> indel_recorder[offsetX2][datalen2][1] = mynoP1PStr;
				vote -> indel_recorder[offsetX2][datalen2][2] = 0;
				vote -> toli[offsetX2][datalen2] = 3;
				vote->coverage_start [offsetX2][datalen2] = offset;
				vote->coverage_end [offsetX2][datalen2] = of_p_16;

				if (vote->max_vote==0) vote->max_vote = 1;
			}
		}
	}

	// Addressing each stub: sorting indel_recorder then move pos
	// This only affects the junction-detection mode.
	// It also adds extra weights to exon-located locations.
	if(cct_context ->do_cell_level_junction_detection) for(x1 =0; x1 < GENE_SCRNA_VOTE_TABLE_SIZE; x1++){
		for(x2 = 0; x2 < vote -> items[x1]; x2++){
			srInt_64 weight_vv = cellCounts_calculate_pos_weight_1sec(cct_context, vote->pos[x1][x2] + vote->coverage_start[x1][x2],  vote->coverage_end[x1][x2] - vote->coverage_start[x1][x2] );

			if(weight_vv>10) vote -> votes[x1][x2] = 1.3*vote -> votes[x1][x2];
			if(vote -> toli[x1][x2]<=3) continue;

			basic_sort( vote -> indel_recorder[x1][x2], vote -> toli[x1][x2] /3, fix_indel_record_order_compare, fix_indel_record_order_exchange );
			int x3;
			int subtract_base_offset = vote -> indel_recorder[x1][x2][2];
			vote -> indel_recorder[x1][x2][2] = 0;
			for(x3=3; x3 <  vote -> toli[x1][x2]; x3+=3) vote -> indel_recorder[x1][x2][x3+2]-= subtract_base_offset;
			vote ->pos[x1][x2] += subtract_base_offset; 
		}
	}
}

int cellCounts_simple_mode_highconf(cellcounts_global_t * cct_context, int thread_no, int applied_subreads, gene_sc_vote_t * vote, char * read_name){
	int x1, vlast = vote -> max_vote;
	for(x1 = 1; x1 < cct_context -> max_candidate_voteIJ_per_read +1; x1++){
		int vdiff = vlast - vote->topK_votes[x1];
		if(vdiff >= 3) return 1;

		vlast = vote->topK_votes[x1];
	}
	return 0;

/*
	int max_I = vote -> max_vote_IJ >> 16;
	int max_J = vote -> max_vote_IJ & 0xffff;

	if(0){
		int tolimax = vote -> toli[max_I][max_J];
		if(tolimax >3)return 0;
	}

	int SRp1_start = vote -> indel_recorder[max_I][max_J][0];
	int SRp1_end = vote -> indel_recorder[max_I][max_J][1];
	return SRp1_start <= 4 && SRp1_end >= applied_subreads - 4 - ( applied_subreads>=12?1:0 );
*/
}

int cellCounts_build_simple_mode_subread_masks(cellcounts_global_t * cct_context, int thread_no, int applied_subreads){
	if(applied_subreads<9) return 0;
	int last_sr_0B = applied_subreads-1;
	int sr_step = (last_sr_0B -1) * 10000 / 4 + 1, sri, ret=0;
	for(sri = 0 ; sri < last_sr_0B*10000+100; sri+= sr_step) ret |= 1<<( sri / 10000 );
	return ret;
}

int cellCounts_get_sample_no_from_rname(cellcounts_global_t * cct_context, int thread_no, char * read_name){
	char * sample_seq=NULL, *sample_qual=NULL, *BC_qual=NULL, *BC_seq=NULL, *UMI_seq=NULL, *UMI_qual=NULL, *lane_str=NULL, *RG=NULL, *testi;
	int rname_trimmed_len=0;
	cellCounts_scan_read_name_str(cct_context, NULL, read_name, &sample_seq, &sample_qual, &BC_seq, &BC_qual, &UMI_seq, &UMI_qual, &lane_str, &RG, &rname_trimmed_len);

	int sample_no = -1;
	if(cct_context -> input_mode == GENE_INPUT_SCRNA_BAM){
		sample_no = 1;  // Only one sample in the BAM mode. A sample may have multiple BAM files but they have the same sample_no.
				// Multiple input samples are mapped/counted in multiple C_cellCounts calls. Each call only does one sample.
	}else if(lane_str){
		int laneno = 0;
		for(testi = lane_str+1; *testi; testi++){
			if(!isdigit(*testi))break;
			laneno = laneno*10 + (*testi)-'0';
		}
		sample_no = cellCounts_get_sample_id(cct_context, sample_seq, laneno); 
	}else if(memcmp("input#", sample_seq,6)==0){
		int lineno = (sample_seq[6]-'0')*1000+(sample_seq[7]-'0')*100+(sample_seq[8]-'0')*10+(sample_seq[9]-'0') +1;
		sample_no = HashTableGet(cct_context -> lineno1B_to_sampleno1B_tab, NULL+lineno)-NULL;
	}else SUBREADprintf("Wrong read name: %s\n", read_name);
	if(0==sample_no)SUBREADprintf("ERROR: ZERO sample_no UNEXPECTED!\n");
	return sample_no;
}

static unsigned char cellCounts_temp_realign_base2code(char base){
	switch(toupper((unsigned char)base)){
		case 'A': return 1;
		case 'C': return 2;
		case 'G': return 3;
		case 'T': return 4;
		default: return 0;
	}
}

static int cellCounts_temp_realign_pack_3bit(const char * bases, int len, unsigned char * out){
	int i, out_bytes = (len * 3 + 7) / 8;
	memset(out, 0, out_bytes);
	for(i = 0; i < len; i++){
		unsigned int code = cellCounts_temp_realign_base2code(bases[i]) & 7u;
		int bit_no = i * 3;
		int byte_no = bit_no >> 3;
		int bit_shift = bit_no & 7;
		out[byte_no] |= (unsigned char)(code << bit_shift);
		if(bit_shift > 5) out[byte_no + 1] |= (unsigned char)(code >> (8 - bit_shift));
	}
	return out_bytes;
}

static unsigned int cellCounts_temp_realign_pack_3bit_prefix(const char * bases, int len){
	unsigned int ret = 0;
	int i, use_len = min(len, 10);
	for(i = 0; i < use_len; i++) ret |= ((unsigned int)cellCounts_temp_realign_base2code(bases[i]) & 7u) << (i * 3);
	return ret;
}

static char cellCounts_temp_realign_code2base(unsigned int code){
	switch(code & 7u){
		case 1: return 'A';
		case 2: return 'C';
		case 3: return 'G';
		case 4: return 'T';
		default: return 'N';
	}
}

static void cellCounts_temp_realign_unpack_3bit(const unsigned char * in, int len, char * out){
	int i;
	for(i = 0; i < len; i++){
		int bit_no = i * 3;
		int byte_no = bit_no >> 3;
		int bit_shift = bit_no & 7;
		unsigned int code = (in[byte_no] >> bit_shift) & 7u;
		if(bit_shift > 5) code |= ((unsigned int)in[byte_no + 1] << (8 - bit_shift)) & 7u;
		out[i] = cellCounts_temp_realign_code2base(code);
	}
	out[len] = 0;
}

static int cellCounts_temp_realign_is_control_byte(unsigned char this_byte){
	return this_byte >= 0xD0 && this_byte <= 0xF0;
}

static void cellCounts_temp_realign_fp_put_byte(cellcounts_temp_file_point_t * temp_fp, unsigned char this_byte){
	if(temp_fp -> fp)putc((int)this_byte, temp_fp -> fp);
	else temp_fp -> realign_temp_memspace[temp_fp -> realign_temp_usedmem++] = this_byte;
}


static int cellCounts_temp_realign_fp_flush_run(cellcounts_temp_file_point_t * temp_fp){
	int repeats = 0;
	unsigned char this_byte = 0;

	if(!temp_fp -> rle_run_active) return 1;
	repeats = (int)temp_fp -> rle_run_repeats;
	this_byte = temp_fp -> rle_run_byte;

	if(repeats ==0){
		if(cellCounts_temp_realign_is_control_byte(this_byte)) cellCounts_temp_realign_fp_put_byte(temp_fp, 0xD0);
		cellCounts_temp_realign_fp_put_byte(temp_fp, this_byte);
	}else if(repeats ==1){
		if(cellCounts_temp_realign_is_control_byte(this_byte)) cellCounts_temp_realign_fp_put_byte(temp_fp, 0xD1);
		else cellCounts_temp_realign_fp_put_byte(temp_fp, this_byte);
		cellCounts_temp_realign_fp_put_byte(temp_fp, this_byte);
	} else {
		unsigned char marker = (unsigned char)(0xD0 + repeats);
		cellCounts_temp_realign_fp_put_byte(temp_fp, marker);
		cellCounts_temp_realign_fp_put_byte(temp_fp, this_byte);
	}

	temp_fp -> rle_run_active = 0;
	temp_fp -> rle_run_repeats = 0;
	temp_fp -> rle_buffer_used = 0;
	return 1;
}

static int cellCounts_temp_realign_fp_write_plain(cellcounts_temp_file_point_t * temp_fp, const unsigned char * plain, int plain_bytes){
	int i;

	for(i = 0; i < plain_bytes; i++){
		unsigned char this_byte = plain[i];
		if(!temp_fp -> rle_run_active){
			temp_fp -> rle_run_active = 1;
			temp_fp -> rle_run_byte = this_byte;
			temp_fp -> rle_run_repeats = 0;
			temp_fp -> rle_buffer_used = 0;
		}else if(this_byte == temp_fp -> rle_run_byte && temp_fp -> rle_run_repeats < 32){
			temp_fp -> rle_buffer[temp_fp -> rle_run_repeats] = this_byte;
			temp_fp -> rle_run_repeats++;
			temp_fp -> rle_buffer_used = temp_fp -> rle_run_repeats;
		}else{ // this byte != run_byte, OR repeats had been 32 (need to be flushed).
			if(!cellCounts_temp_realign_fp_flush_run(temp_fp)) return 0;
			temp_fp -> rle_run_active = 1;
			temp_fp -> rle_run_byte = this_byte;
			temp_fp -> rle_run_repeats = 0;
			temp_fp -> rle_buffer_used = 0;
		}
	}
	return 1;
}

int cellCounts_temp_realign_fp_fgetc(cellcounts_temp_file_point_t * temp_fp){
	if(temp_fp -> realign_temp_memspace){
		if(temp_fp -> realign_temp_usedmem == temp_fp -> realign_temp_capamem)return EOF;
		int rv = (int)temp_fp -> realign_temp_memspace[temp_fp -> realign_temp_usedmem++];
		return rv;
	}else return fgetc(temp_fp -> fp);
}

static int cellCounts_temp_realign_fp_read_plain(cellcounts_temp_file_point_t * temp_fp, unsigned char * plain, int plain_bytes){
	int out_used = 0;

	while(out_used < plain_bytes){
		if(temp_fp -> rle_buffer_used > 0){
			int copy_bytes = min((int)temp_fp -> rle_buffer_used, plain_bytes - out_used);
			memcpy(plain + out_used, temp_fp -> rle_buffer, copy_bytes);
			out_used += copy_bytes;
			temp_fp -> rle_buffer_used -= (unsigned char)copy_bytes;
			if(temp_fp -> rle_buffer_used > 0) memmove(temp_fp -> rle_buffer, temp_fp -> rle_buffer + copy_bytes, temp_fp -> rle_buffer_used);
			continue;
		}

		int marker = cellCounts_temp_realign_fp_fgetc(temp_fp);
		if(marker == EOF) return out_used;

		if(cellCounts_temp_realign_is_control_byte((unsigned char)marker)){
			int raw_byte = cellCounts_temp_realign_fp_fgetc(temp_fp);
			int repeats, emit_now, keep_now;
			if(raw_byte == EOF) return -1;
			repeats = marker - (0xD0-1);
			emit_now = min(repeats, plain_bytes - out_used);
			memset(plain + out_used, raw_byte, emit_now);
			out_used += emit_now;
			keep_now = repeats - emit_now;
			if(keep_now > 0){
				memset(temp_fp -> rle_buffer, raw_byte, keep_now);
				temp_fp -> rle_buffer_used = (unsigned char)keep_now;
			}
		}else{
			plain[out_used++] = (unsigned char)marker;
		}
	}

	return out_used;
}

static int cellCounts_temp_realign_fp_read_plain_exact(cellcounts_temp_file_point_t * temp_fp, unsigned char * plain, int plain_bytes, int allow_eof){
	int rlen = cellCounts_temp_realign_fp_read_plain(temp_fp, plain, plain_bytes);
	if(rlen < 0) return -1;
	if(rlen == plain_bytes) return 1;
	if(rlen == 0 && allow_eof) return 0;
	return -1;
}

static int cellCounts_temp_realign_fp_finish_write(cellcounts_temp_file_point_t * temp_fp){
	if(!cellCounts_temp_realign_fp_flush_run(temp_fp)) return 0;
	if(temp_fp -> fp)return 0 == fflush(temp_fp -> fp);
	else return 0;
}

void cellcounts_temp_file_destroy(cellcounts_global_t * cct_context, char * tmp_fname, cellcounts_temp_file_point_t *temp_fp){
	if(temp_fp -> realign_temp_memspace){
		free(temp_fp -> realign_temp_memspace);
		temp_fp -> realign_temp_memspace = NULL;
	}
	if(temp_fp -> fp){
		fclose(temp_fp -> fp);
		unlink(tmp_fname);
		temp_fp -> fp = NULL;
	}
}

void cellcounts_temp_file_fclose(cellcounts_temp_file_point_t * temp_fp){
	if(temp_fp -> fp) fclose(temp_fp -> fp);
	temp_fp -> fp = NULL;
	temp_fp -> rle_buffer_used = 0;
	temp_fp -> rle_run_byte = 0;
	temp_fp -> rle_run_repeats = 0;
	temp_fp -> rle_run_active = 0;
}

void cellcounts_temp_file_open(cellcounts_global_t * cct_context, char * tmp_fname, int for_writting, cellcounts_temp_file_point_t * temp_fp){
	if(cct_context -> cell_level_junction_memory_temp){
		if(for_writting){
			temp_fp -> realign_temp_memspace = malloc(TEMP_BINFILE_MEMORY_SIZE_INIT);
			if(!temp_fp -> realign_temp_memspace){
				SUBREADprintf("\nEEROR: cannot allocate memory for saving alignment results. Please disable the memory mode by specifying 'binaryTempMemory=FALSE'.\n");
				assert( temp_fp -> realign_temp_memspace );
				return;
			}
			temp_fp -> realign_temp_capamem = TEMP_BINFILE_MEMORY_SIZE_INIT;
		}else temp_fp -> realign_temp_capamem = temp_fp -> realign_temp_usedmem; // read mode: capa=current_available
		temp_fp -> realign_temp_usedmem = 0;
		temp_fp -> fp = NULL;
	}else{
		temp_fp -> fp = fopen(tmp_fname, for_writting?"wb":"rb");
		temp_fp -> realign_temp_memspace = NULL;
	}
	temp_fp -> for_writting = for_writting;
}

cellcounts_temp_file_point_t * cellCounts_select_and_write_temps_open_fp(cellcounts_global_t * cct_context, int thread_no){
	if(thread_no < 0 || thread_no >= 64) return NULL;

	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	cellcounts_temp_file_point_t * temp_fp = &thread_context -> realign_temp_fp;

	if(! (temp_fp -> fp || temp_fp -> realign_temp_memspace)){
		char tmp_fname[MAX_FILE_NAME_LENGTH + 120];
		SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH + 120, "%s/temp-cellcounts-realign-%06d-%03d.tmpbin", cct_context -> temp_file_dir, getpid(), thread_no);
		memset(temp_fp, 0, sizeof(cellcounts_temp_file_point_t));
		cellcounts_temp_file_open(cct_context, tmp_fname, 1, temp_fp);
		if(temp_fp -> fp)setvbuf(temp_fp -> fp, thread_context -> tempbin_v_buffer, _IOFBF , SCRNA_VBUFF_SIZE);
	}
	return temp_fp;
}

static unsigned char * cellCounts_temp_realign_ensure_buf(unsigned char ** buf, int * capacity, int required){
	unsigned char * new_buf;
	if(required < 1) required = 1;
	if(*buf && *capacity >= required) return *buf;
	new_buf = realloc(*buf, required+150);
	if(!new_buf) return NULL;
	*buf = new_buf;
	*capacity = (required+150);
	return *buf;
}

int cellCounts_select_and_write_temps(cellcounts_global_t * cct_context, int thread_no, int sample_i, gene_sc_vote_t * votetab, char * read_name, char * read_text, char * read_bin, char * read_qual, int read_len, cellcounts_vote_number_t all_subreads) {
	int i, j, distinct_vote_number_i;
	int saved_alignments = 0;
	int max_saved_alignments = min(cct_context -> max_candidate_voteIJ_per_read, SCRNA_HIGHEST_REPORTED_ALIGNMENTS);
	//int max_saved_alignments = min(cct_context -> max_reported_alignments_per_read, SCRNA_HIGHEST_REPORTED_ALIGNMENTS);
	int index_gap_width = cct_context -> current_index -> index_gap;
	char * sample_seq=NULL, *sample_qual=NULL, *BC_qual=NULL, *BC_seq=NULL, *UMI_seq=NULL, *UMI_qual=NULL, *lane_str=NULL, *RG=NULL;
	char * sample_seq_end = NULL, * sample_qual_end = NULL;
	int rname_trimmed_len=0;
	char tmp_bc[MAX_READ_NAME_LEN+1], temp_bq_buf[100];
	struct TempForRealign temprec;
	cellcounts_temp_file_point_t * temp_fp;

	(void)all_subreads;
	(void)index_gap_width;
	(void)sample_seq;
	(void)sample_qual;
	(void)lane_str;
	(void)RG;
	(void)rname_trimmed_len;

	memset(&temprec, 0, sizeof(temprec));
	cellCounts_scan_read_name_str(cct_context, NULL, read_name, &sample_seq, &sample_qual, &BC_seq, &BC_qual, &UMI_seq, &UMI_qual, &lane_str, &RG, &rname_trimmed_len);

	if(sample_seq){
		size_t sample_seq_len = 0;
		sample_seq_end = strchr(sample_seq, '|');
		sample_seq_len = sample_seq_end ? (size_t)(sample_seq_end - sample_seq) : strlen(sample_seq);
		if(sample_seq_len > USHRT_MAX) return 1;
		temprec.sample_seq_length = (unsigned short)sample_seq_len;
	}
	if(sample_qual){
		size_t sample_qual_len = 0;
		sample_qual_end = strchr(sample_qual, '|');
		sample_qual_len = sample_qual_end ? (size_t)(sample_qual_end - sample_qual) : strlen(sample_qual);
		if(sample_qual_len > USHRT_MAX) return 1;
		temprec.sample_qual_length = (unsigned short)sample_qual_len;
	}

	temprec.sample_number = sample_i > 0 ? (unsigned int)sample_i : 0u;
	temprec.read_number = (unsigned int)strtoul(read_name + 1, NULL, 10);
	temprec.read_length = read_len;
	temprec.raw_umi_sequence = cellCounts_temp_realign_pack_3bit_prefix(UMI_seq ? UMI_seq : "", UMI_seq ? cct_context -> UMI_length : 0);

	if(sample_i>=0){
		int top_distinct_vote_numbers[cct_context -> max_distinct_top_vote_numbers];
		memset(top_distinct_vote_numbers, 0, sizeof(top_distinct_vote_numbers));
		int CR15GLS = (read_len - 15 - index_gap_width)<<16;
		int subread_step =  CR15GLS /(cct_context -> total_subreads_per_read -1);
		if(subread_step<(index_gap_width<<16))subread_step = index_gap_width<<16;

		if(votetab -> max_vote >= cct_context -> min_votes_per_mapped_read){
			for(i=0; i<GENE_SCRNA_VOTE_TABLE_SIZE; i++){
				for(j=0; j< votetab->items[i]; j++){
					int vv = votetab -> votes[i][j];
					if(vv>=cct_context -> min_votes_per_mapped_read)cellCounts_update_top_three(cct_context, top_distinct_vote_numbers, vv);
				}
			}
	//for(i = 0; i < GENE_SCRNA_VOTE_TABLE_SIZE; i++)fprintf(stderr,"VOTE_N = %d # %d\n", top_distinct_vote_numbers[i],i);

			for(distinct_vote_number_i = 0 ; distinct_vote_number_i < cct_context -> max_distinct_top_vote_numbers; distinct_vote_number_i ++){
				int this_vote_N = top_distinct_vote_numbers[distinct_vote_number_i];
				if(this_vote_N < 1 || (top_distinct_vote_numbers[0] - this_vote_N > cct_context -> max_differential_from_top_vote_number )) break;
				for(i=0; i<GENE_SCRNA_VOTE_TABLE_SIZE; i++){
					for(j=0; j< votetab->items[i]; j++){
						int vv = votetab->votes[i][j];
						if(vv != this_vote_N) continue;
						if(saved_alignments >= max_saved_alignments) break;
						temprec.num_of_votes[saved_alignments] = (unsigned short)vv;
						temprec.voted_position[saved_alignments] = votetab->pos[i][j];

						int perfect_align_srno = votetab->indel_recorder[i][j][0]-1; // indel_recorder is subread_no + 1
						temprec.coverage_start[saved_alignments] = ((subread_step * perfect_align_srno) >> 16);
						perfect_align_srno = votetab->indel_recorder[i][j][1]-1;
						temprec.coverage_end[saved_alignments] = ((subread_step * perfect_align_srno) >> 16) +15;

						temprec.flags[saved_alignments] = votetab->masks[i][j] ? SAM_FLAG_REVERSE_STRAND_MATCHED : 0;
						saved_alignments++;
	//fprintf(stderr,"ADD_TEMP_REALIGN_FILE of ALN %d : vote = %d at VTAB %d_%d   TOP_N %d is %d / %d\n", saved_alignments, vv, i, j,  this_vote_N, distinct_vote_number_i , cct_context -> max_differential_from_top_vote_number );

					}
					if(saved_alignments >= max_saved_alignments) break;
				}
				if(saved_alignments >= max_saved_alignments) break;
			}
		}
	}

	temprec.saved_alignments = saved_alignments;
	temp_fp = cellCounts_select_and_write_temps_open_fp(cct_context, thread_no);
	if(temp_fp){
		int umi_len = cct_context -> UMI_length;
		int bc_len = 0;
		int bc_qual_len = 0;
		int packed_read_len = (read_len * 3 + 7) / 8;
		int packed_bcumi_len = 0;
		int total_bc_umi_len = 0;
		size_t record_size_sz;
		int record_size;
		size_t sample_seq_len_sz = temprec.sample_seq_length;
		size_t sample_qual_len_sz = temprec.sample_qual_length;
		unsigned char * record, * wp;
		cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;

		if(BC_seq){
			char * bc_end = strchr(BC_seq, '|');
			char * bc_qual_end = BC_qual ? strchr(BC_qual, '|') : NULL;
			bc_qual_len = bc_qual_end ? (int)(bc_qual_end - BC_qual) : (BC_qual ? (int)strlen(BC_qual) : 0);
			total_bc_umi_len = bc_end ? (int)(bc_end - BC_seq) : 0;
			bc_len = cct_context -> known_cell_barcode_length;
			if(bc_len < 0){
				bc_len = total_bc_umi_len - umi_len;
				if(bc_len < 0) bc_len = 0;
			}
			if(bc_end){
				char * temp_bq=NULL;
				int copy_len = total_bc_umi_len;
				if(copy_len > MAX_READ_NAME_LEN) copy_len = MAX_READ_NAME_LEN;
				memcpy(tmp_bc, BC_seq, copy_len);
				if(cct_context -> VisiumHD_barcode_to_best_mapping){
					temp_bq = temp_bq_buf;
					memcpy(temp_bq, BC_qual, copy_len);
				}
				tmp_bc[copy_len] = 0;

				{
					int cell_no = cellCounts_get_cellbarcode_no(cct_context, thread_no, tmp_bc, temp_bq);
					temprec.cell_number = cell_no >= 0 ? (unsigned int)(cell_no + 1) : 0u;
				}
			}
		}
		packed_bcumi_len = ((bc_len + umi_len) * 3 + 7) / 8;
		record_size_sz = sizeof(int) + sizeof(unsigned int) * 4 + sizeof(unsigned short) * 2 + sizeof(int) * 3
			+ sample_seq_len_sz + sample_qual_len_sz
			+ saved_alignments * (sizeof(unsigned short) + sizeof(unsigned int) + sizeof(unsigned short) + sizeof(unsigned short) + sizeof(unsigned char))
			+ read_len + packed_read_len + bc_qual_len + packed_bcumi_len;
		if(record_size_sz > (size_t)2147483647) return 1;
		record_size = (int)record_size_sz;
		record = cellCounts_temp_realign_ensure_buf(&thread_context -> temp_realign_record_buf, &thread_context -> temp_realign_record_capacity, record_size);
		if(!record) return 1;
		wp = record;

		memcpy(wp, &temprec.saved_alignments, sizeof(int)); wp += sizeof(int);
		memcpy(wp, &temprec.sample_number, sizeof(unsigned int)); wp += sizeof(unsigned int);
		memcpy(wp, &temprec.cell_number, sizeof(unsigned int)); wp += sizeof(unsigned int);
		memcpy(wp, &temprec.read_number, sizeof(unsigned int)); wp += sizeof(unsigned int);
		memcpy(wp, &temprec.raw_umi_sequence, sizeof(unsigned int)); wp += sizeof(unsigned int);
		memcpy(wp, &temprec.sample_seq_length, sizeof(unsigned short)); wp += sizeof(unsigned short);
		memcpy(wp, &temprec.sample_qual_length, sizeof(unsigned short)); wp += sizeof(unsigned short);
		memcpy(wp, &temprec.read_length, sizeof(int)); wp += sizeof(int);
		memcpy(wp, &bc_len, sizeof(int)); wp += sizeof(int);
		memcpy(wp, &umi_len, sizeof(int)); wp += sizeof(int);

		if(sample_seq_len_sz > 0){
			memcpy(wp, sample_seq, sample_seq_len_sz);
			wp += sample_seq_len_sz;
		}
		if(sample_qual_len_sz > 0){
			memcpy(wp, sample_qual, sample_qual_len_sz);
			wp += sample_qual_len_sz;
		}

		if(saved_alignments > 0){
			memcpy(wp, temprec.num_of_votes, sizeof(unsigned short) * saved_alignments);
			wp += sizeof(unsigned short) * saved_alignments;
			memcpy(wp, temprec.voted_position, sizeof(unsigned int) * saved_alignments);
			wp += sizeof(unsigned int) * saved_alignments;
			memcpy(wp, temprec.coverage_start, sizeof(unsigned short) * saved_alignments);
			wp += sizeof(unsigned short) * saved_alignments;
			memcpy(wp, temprec.coverage_end, sizeof(unsigned short) * saved_alignments);
			wp += sizeof(unsigned short) * saved_alignments;
			memcpy(wp, temprec.flags, sizeof(unsigned char) * saved_alignments);
			wp += sizeof(unsigned char) * saved_alignments;
		}

		memcpy(wp, read_qual, read_len);
		wp += read_len;
		if(packed_read_len > 0){
			cellCounts_temp_realign_pack_3bit(read_text, read_len, wp);
			wp += packed_read_len;
		}

		if(BC_qual && bc_qual_len > 0){
			memcpy(wp, BC_qual, bc_qual_len);
			wp += bc_qual_len;
		}
		if(packed_bcumi_len > 0 && BC_seq){
			cellCounts_temp_realign_pack_3bit(BC_seq, bc_len + umi_len, wp);
			wp += packed_bcumi_len;
		}

		if(temp_fp -> realign_temp_memspace && temp_fp -> realign_temp_usedmem >= temp_fp -> realign_temp_capamem - 100 * MAX_SCRNA_READ_LENGTH){
			temp_fp -> realign_temp_capamem *= 1.5;
			temp_fp -> realign_temp_memspace = realloc(temp_fp -> realign_temp_memspace, temp_fp -> realign_temp_capamem);
		}

		if(!cellCounts_temp_realign_fp_write_plain(temp_fp, (unsigned char *)&record_size, sizeof(int))) return 1;
		if(record_size > 0 && !cellCounts_temp_realign_fp_write_plain(temp_fp, record, record_size)) return 1;
	}
	return 0;
}
void * cellCounts_select_and_write_alignments_from_temp(void * pr);

int cellCounts_do_realign(cellcounts_global_t * cct_context){
	int input_thread_no, task=STEP_VOTING; // task for realignment equals voting (main step).
	// init thread contexts

	int current_thread_no, smpno;
	cellcounts_align_thread_t * thread_contexts = calloc(sizeof(cellcounts_align_thread_t) , cct_context->total_threads);
	cct_context -> all_thread_contexts = thread_contexts;

	int ret_values[64];

	for(current_thread_no = 0 ; current_thread_no < cct_context->total_threads ; current_thread_no ++) {
		thread_contexts[current_thread_no].thread_no = current_thread_no;
		cellCounts_prepare_context_for_align(cct_context, current_thread_no, task);
		cellCounts_init_topKbuff(cct_context, current_thread_no);
	}

	// For each input thread file, start all threads. Total runs: threads ^ 2.
	for(input_thread_no = 0; input_thread_no < cct_context->total_threads; input_thread_no++){
		char tmp_fname[MAX_FILE_NAME_LENGTH + 120];
		cellcounts_temp_file_point_t temp_fp, *ptr_temp_fp;
		ptr_temp_fp = &temp_fp;

		SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH + 120, "%s/temp-cellcounts-realign-%06d-%03d.tmpbin", cct_context -> temp_file_dir, getpid(), input_thread_no);

		if(!cct_context ->cell_level_junction_memory_temp){
			memset(&temp_fp, 0, sizeof(temp_fp));
			cellcounts_temp_file_open(cct_context, tmp_fname, 0, &temp_fp);
			if(temp_fp.fp)setvbuf(temp_fp.fp, thread_contexts[0].tempbin_v_buffer, _IOFBF , SCRNA_VBUFF_SIZE);
		}

		for(current_thread_no = 0 ; current_thread_no < cct_context->total_threads ; current_thread_no ++) {
			if(cct_context ->cell_level_junction_memory_temp){
				ptr_temp_fp = & thread_contexts[current_thread_no].realign_temp_fp;
				memset(ptr_temp_fp,0,sizeof(*ptr_temp_fp));
				ptr_temp_fp -> realign_temp_memspace = cct_context -> all_thread_realign_fp_ptrs[current_thread_no] ;
				ptr_temp_fp -> realign_temp_capamem  = cct_context -> all_thread_realign_fp_ints[current_thread_no*2] ; // this is the total written bytes
			}

			void ** thr_parameters = malloc(sizeof(void*)*4);
			thr_parameters[0] = cct_context;
			thr_parameters[1] = NULL+current_thread_no;
			thr_parameters[2] = ptr_temp_fp;
			pthread_create(&thread_contexts[current_thread_no].thread, NULL, cellCounts_select_and_write_alignments_from_temp, thr_parameters);
		}

		for(current_thread_no = 0 ; current_thread_no < cct_context->total_threads ; current_thread_no ++) {
			pthread_join(thread_contexts[current_thread_no].thread, NULL);
			for(smpno = 0; smpno < cct_context-> sample_sheet_table -> numOfElements; smpno ++){ 
				cct_context -> mapped_reads_per_sample[smpno] += thread_contexts[current_thread_no].mapped_reads_per_sample[smpno];
				cct_context -> assigned_reads_per_sample[smpno] += thread_contexts[current_thread_no].assigned_reads_per_sample[smpno];
				cct_context -> reads_per_sample[smpno] += thread_contexts[current_thread_no].reads_per_sample[smpno];
			}
			cct_context -> reads_per_sample[smpno] += thread_contexts[current_thread_no].reads_per_sample[smpno]; 
		}

		if(cct_context -> cell_level_junction_memory_temp){
			for(current_thread_no = 0 ; current_thread_no < cct_context->total_threads ; current_thread_no ++) 
				cellcounts_temp_file_destroy(cct_context, tmp_fname, &thread_contexts[current_thread_no].realign_temp_fp);
			break; // if it is run on the memry mode, each thread uses its own temp fp and only run once.
		} else cellcounts_temp_file_destroy(cct_context, tmp_fname, &temp_fp);
	}

	// release thread contexts
	for(current_thread_no = 0 ; current_thread_no < cct_context->total_threads ; current_thread_no ++) {
		cellCounts_free_topKbuff(cct_context, current_thread_no);
		cellCounts_release_context_from_align(cct_context, current_thread_no, task);
	}
	free(thread_contexts);
}

void * cellCounts_select_and_write_alignments_from_temp(void * pr){
	void ** my_parameters = pr;
	cellcounts_global_t * cct_context = my_parameters[0];
	int thread_no = my_parameters[1]-NULL;
	cellcounts_temp_file_point_t * temp_fp = my_parameters[2];

	free(pr);
	int rc = 0;
	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;
	if(!temp_fp || !(temp_fp -> fp|| temp_fp -> realign_temp_memspace )) return NULL+1;

	int processed_records =0;
	while(1){
		int record_size = 0;
		unsigned char * record = NULL;
		unsigned char * rp = NULL;
		unsigned char * record_end = NULL;
		unsigned char * work = NULL;
		int saved_alignments = 0, read_len = 0, bc_len = 0, umi_len = 0;
		unsigned short sample_seq_length = 0, sample_qual_length = 0;
		unsigned int sample_number = 0, cell_number = 0, read_number = 0, raw_umi_sequence = 0;
		unsigned short * votes_buf = NULL;
		unsigned int * pos_buf = NULL;
		unsigned short * cstart_buf = NULL;
		unsigned short * cend_buf = NULL;
		unsigned char * flags_buf = NULL;
		char * read_qual_fwd = NULL;
		char * read_qual_rev = NULL;
		char * read_text_fwd = NULL;
		char * read_text_rev = NULL;
		char * read_name = NULL;
		char * bcumi_seq_buf = NULL;
		char * bcumi_qual_buf = NULL;
		char * sample_seq = NULL;
		char * sample_qual = NULL;
		char * read_qual = NULL;
		char * bc_qual = NULL;
		int i;

		if(0==cct_context ->cell_level_junction_memory_temp)cellCounts_lock_occupy(&cct_context -> input_dataset_lock);

		rc = cellCounts_temp_realign_fp_read_plain_exact(temp_fp, (unsigned char *)&record_size, sizeof(int), 1);
//fprintf(stderr,"DO_GET %lld < %lld   ret %d   size %d\n",temp_fp->realign_temp_usedmem,temp_fp->realign_temp_capamem, rc, record_size);
		if(rc <= 0){
			if(0==cct_context ->cell_level_junction_memory_temp)cellCounts_lock_release(&cct_context -> input_dataset_lock);
			if(rc == 0) break;
			return NULL+1;
		}
		if(record_size < 0){
			if(0==cct_context ->cell_level_junction_memory_temp)cellCounts_lock_release(&cct_context -> input_dataset_lock);
			return NULL+1;
		}

		record = cellCounts_temp_realign_ensure_buf(&thread_context -> temp_realign_record_buf, &thread_context -> temp_realign_record_capacity, record_size);
		if(!record){
			if(0==cct_context ->cell_level_junction_memory_temp)cellCounts_lock_release(&cct_context -> input_dataset_lock);
			return NULL+1;
		}

		if(record_size > 0 && 1 != cellCounts_temp_realign_fp_read_plain_exact(temp_fp, record, record_size, 0)){
			if(0==cct_context ->cell_level_junction_memory_temp)cellCounts_lock_release(&cct_context -> input_dataset_lock);
			return NULL+1;
		}
//fprintf(stderr,"DO_P2 %lld < %lld   ret %p\n",temp_fp->realign_temp_usedmem,temp_fp->realign_temp_capamem, record);

		if(0==cct_context ->cell_level_junction_memory_temp)cellCounts_lock_release(&cct_context -> input_dataset_lock);

		rp = record;
		record_end = record + record_size;

		if((size_t)(record_end - rp) < sizeof(int)) return NULL+1;
		memcpy(&saved_alignments, rp, sizeof(int)); rp += sizeof(int);
		if((size_t)(record_end - rp) < sizeof(unsigned int) * 4 + sizeof(unsigned short) * 2 + sizeof(int) * 3) return NULL+1;
		memcpy(&sample_number, rp, sizeof(unsigned int)); rp += sizeof(unsigned int);
		memcpy(&cell_number, rp, sizeof(unsigned int)); rp += sizeof(unsigned int);
		memcpy(&read_number, rp, sizeof(unsigned int)); rp += sizeof(unsigned int);
		memcpy(&raw_umi_sequence, rp, sizeof(unsigned int)); rp += sizeof(unsigned int);
		memcpy(&sample_seq_length, rp, sizeof(unsigned short)); rp += sizeof(unsigned short);
		memcpy(&sample_qual_length, rp, sizeof(unsigned short)); rp += sizeof(unsigned short);
		memcpy(&read_len, rp, sizeof(int)); rp += sizeof(int);
		memcpy(&bc_len, rp, sizeof(int)); rp += sizeof(int);
		memcpy(&umi_len, rp, sizeof(int)); rp += sizeof(int);

		(void)cell_number;
		(void)raw_umi_sequence;

		if(saved_alignments < 0 || saved_alignments > SCRNA_HIGHEST_REPORTED_ALIGNMENTS) return NULL+1;
		if(read_len < 0 || read_len > MAX_SCRNA_READ_LENGTH) return NULL+1;
		if(bc_len < 0 || bc_len > MAX_CELLBC_LEN) return NULL+1;
		if(umi_len < 0 || umi_len > MAX_UMI_LEN) return NULL+1;

		size_t read_text_slot = (size_t)MAX_SCRNA_READ_LENGTH + 2;
		{
			size_t read_name_slot = (size_t)MAX_READ_NAME_LEN + 1;
			size_t bcumi_slot = (size_t)MAX_CELLBC_LEN + (size_t)MAX_UMI_LEN + 1;
			size_t votes_bytes = (size_t)saved_alignments * sizeof(unsigned short);
			size_t pos_bytes = (size_t)saved_alignments * sizeof(unsigned int);
			size_t cstart_bytes = (size_t)saved_alignments * sizeof(unsigned short);
			size_t cend_bytes = (size_t)saved_alignments * sizeof(unsigned short);
			size_t flags_bytes = (size_t)saved_alignments * sizeof(unsigned char);
			size_t read_qual_len = (size_t)read_len;
			size_t packed_read_len = (size_t)((read_len * 3 + 7) / 8);
			size_t bcqual_len = (size_t)(bc_len + umi_len);
			size_t packed_bcumi_len = (size_t)(((bc_len + umi_len) * 3 + 7) / 8);
			size_t need_bytes = (size_t)sample_seq_length + (size_t)sample_qual_length + votes_bytes + pos_bytes + cstart_bytes + cend_bytes + flags_bytes + read_qual_len + packed_read_len + bcqual_len + packed_bcumi_len;

			size_t work_required = read_text_slot * 4 /* RTEXT, REV_TEXT, RQUAL, REV_QUAL */ + read_name_slot + bcumi_slot * 2;

			if((size_t)(record_end - rp) < need_bytes) return NULL+1;
			if(work_required > (size_t)INT_MAX) return NULL+1;

			work = cellCounts_temp_realign_ensure_buf(&thread_context -> temp_realign_work_buf, &thread_context -> temp_realign_work_capacity, (int)work_required);
			if(!work) return NULL+1;
			read_text_fwd = (char *)work;
			read_text_rev = read_text_fwd + read_text_slot;
			read_qual_fwd = read_text_rev + read_text_slot;
			read_qual_rev = read_qual_fwd + read_text_slot;
			read_name = read_qual_rev + read_text_slot;
			bcumi_seq_buf = read_name + read_name_slot;
			bcumi_qual_buf = bcumi_seq_buf + bcumi_slot;

			if(read_len >= (int)read_text_slot || bc_len + umi_len >= (int)bcumi_slot){
				return NULL+1;
			}

			sample_seq = (char *)rp; rp += sample_seq_length;
			sample_qual = (char *)rp; rp += sample_qual_length;
			votes_buf = (unsigned short *)rp; rp += votes_bytes;
			pos_buf = (unsigned int *)rp; rp += pos_bytes;
			cstart_buf = (unsigned short *)rp; rp += cstart_bytes;
			cend_buf = (unsigned short *)rp; rp += cend_bytes;
			flags_buf = (unsigned char *)rp; rp += flags_bytes;
			read_qual = (char *)rp; rp += read_qual_len;
			memcpy(read_qual_fwd, read_qual, read_qual_len);
			{
				unsigned char * packed_read = rp;
				rp += packed_read_len;
				bc_qual = (char *)rp;
				rp += bcqual_len;
				{
					unsigned char * packed_bcumi = rp;
					rp += packed_bcumi_len;

					if(rp != record_end) return NULL+1;

					cellCounts_temp_realign_unpack_3bit(packed_read, read_len, read_text_fwd);
					cellCounts_temp_realign_unpack_3bit(packed_bcumi, bc_len + umi_len, bcumi_seq_buf);
					if(bcqual_len > 0){
						size_t copy_len = bcqual_len;
						if(copy_len > bcumi_slot - 1) copy_len = bcumi_slot - 1;
						memcpy(bcumi_qual_buf, bc_qual, copy_len);
						bcumi_qual_buf[copy_len] = 0;
					}else bcumi_qual_buf[0] = 0;
				}
			}
		}

		SUBreadSprintf(read_name, MAX_READ_NAME_LEN + 1, "R%011u|%s|%s|%.*s|%.*s", read_number, bcumi_seq_buf, bcumi_qual_buf, (int)sample_seq_length, sample_seq ? sample_seq : "", (int)sample_qual_length, sample_qual ? sample_qual : "");

		thread_context -> alignment_repating_table = HashTableCreate(50);

		thread_context -> total_voteIJs_to_write = 0;
		if(sample_number >0){
			thread_context -> populating_voteIJ_buf_index=0;
			for(i = 0; i < saved_alignments; i++){
				char * this_read_text = read_text_fwd;
				if(flags_buf[i] & SAM_FLAG_REVERSE_STRAND_MATCHED){
					memcpy(read_text_rev, read_text_fwd, read_len + 1);
					reverse_read(read_text_rev, read_len, GENE_SPACE_BASE);
					memcpy(read_qual_rev, read_qual_fwd, read_len + 1);
					reverse_quality(read_qual_rev, read_len);
					this_read_text = read_text_rev;
				}
				cellCounts_explain_one_alignment(cct_context, thread_no, sample_number ? (int)sample_number : -1, read_name, this_read_text, read_len, cstart_buf[i], cend_buf[i], pos_buf[i], (flags_buf[i] & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0, votes_buf[i]);
			}
			HashTableDestroy(thread_context -> alignment_repating_table);
			thread_context -> total_voteIJs_to_write = min(thread_context -> total_voteIJs_to_write, cct_context -> max_reported_alignments_per_read);
		}

		int distinct_vote_number_i;
		if(thread_context -> total_voteIJs_to_write) {
			int sorting_index [thread_context -> populating_voteIJ_buf_index];
			for(distinct_vote_number_i = 0 ; distinct_vote_number_i < thread_context -> populating_voteIJ_buf_index; distinct_vote_number_i ++) sorting_index [distinct_vote_number_i ] = distinct_vote_number_i ;
			void * sorting_ptr[2];
			sorting_ptr[0] = thread_context;
			sorting_ptr[1] = sorting_index;
			//sort : large number first
			basic_sort(sorting_ptr , thread_context -> populating_voteIJ_buf_index , sort_readscore_compare_LargeFirst, sort_readscore_exchange); // The last many records are 0-score records. Only "total_voteIJs_to_write" records ahead are worth writting (score > 0).

			for(thread_context -> writing_voteID_buf_index = 0 ; thread_context -> writing_voteID_buf_index < thread_context -> total_voteIJs_to_write; thread_context -> writing_voteID_buf_index ++){
				int myno = sorting_index[thread_context -> writing_voteID_buf_index ];
				if(thread_context -> reporting_scores[ myno ] < 1)break;
				int reverse_text_offset = (thread_context -> reporting_flags[myno] & SAM_FLAG_REVERSE_STRAND_MATCHED)?read_text_slot:0;
				if(0) /* DO NOT need to write read_qual. only need is read_qual_fwd or rev */if(reverse_text_offset >0 && 0==read_qual[reverse_text_offset]){
					strcpy(read_qual+reverse_text_offset, read_qual);
					reverse_quality(read_qual+reverse_text_offset, read_len);
				}
#warning "===== CURRENT JUNCTION DETECTION DOESN'T NEED ANTI_SUPPORT OUT ====="
				if(0)if(cct_context -> do_cell_level_junction_detection && sample_number>0)cellCounts_add_supported_unsupported_reads_from_cigar( cct_context, thread_no, sample_number, myno);
				cellCounts_write_read_in_batch_bin(cct_context, thread_no, sample_number, myno, read_name, read_text_fwd + reverse_text_offset, read_qual_fwd+reverse_text_offset, read_text_fwd , read_qual_fwd , read_len);
			}
		} else cellCounts_write_read_in_batch_bin(cct_context, thread_no, sample_number, -1, read_name, read_text_fwd, read_qual, read_text_fwd, read_qual, read_len);
		processed_records++;
	}
	return NULL;
}

int cellCounts_do_jtab_or_voting(cellcounts_global_t * cct_context, int thread_no, int task) {
	subread_read_number_t current_read_number=0;
	char * read_text, * qual_text;

	char read_name[MAX_READ_NAME_LEN+1];
	char read_bin[REVERSED_READ_BIN_OFFSET * 2];
	int read_len=0;

	read_text = malloc(MAX_SCRNA_READ_LENGTH * 2+2);
	qual_text = malloc(MAX_SCRNA_READ_LENGTH * 2+2);

	temp_votes_per_read_t prefill_ptrs;
	gene_sc_vote_t * vote_me = malloc(sizeof(gene_sc_vote_t));

	if(vote_me==NULL) {
		SUBREADprintf("Cannot allocate voting memory.\n");
		return -1;
	}

	int index_gap_width = cct_context -> current_index -> index_gap;

	while(!cct_context -> has_error) {
		int subread_no;
		int is_reversed, applied_subreads = 0;

		// Read name format: like R00000000059|TNACCCGCCTGCCTCGGCGCGGGGCGNG|D#DDDDDD-DDDDDDD-D-DDDDD<D#D|TNCCCGGGNNGTCGCNNCGN|D#DD-DDD##DDDDD##DD#|@RgLater@L001.
		// The index sequence / quality are in the 4th and 5th columns. They are only used for writing I1/I2 fastq.gz output.
		// For general read alignment, only the 2nd and 3rd columns for cell barcode and UMI are used.

		cellCounts_fetch_next_read_pair(cct_context, thread_no,  &read_len, read_name, read_text, qual_text, &current_read_number);
		if(current_read_number < 0) break;
		if(current_read_number >= cct_context-> reads_per_chunk) break;

		if(read_len< 16) continue;
		int sample_i = cellCounts_get_sample_no_from_rname(cct_context, thread_no, read_name); // Sample_i is 1-based. It can NEVER be 0 (see function). "Not found" = -1
		if(sample_i >=0){
			int CR15GLS = (read_len - 15 - index_gap_width)<<16;
			int subread_step =  CR15GLS /(cct_context -> total_subreads_per_read -1);
			if(subread_step<(index_gap_width<<16))subread_step = index_gap_width<<16;
			applied_subreads = 1 + CR15GLS / subread_step;

			int building_rbin_offset = 0, read_text_rev_offset =0;
			for(is_reversed = 0; is_reversed<2; is_reversed++) {
				gehash_key_t subread_integer = 0;
				int last_vote_rpos = -16;

				for(subread_no=0; subread_no < applied_subreads ; subread_no++) {
					int subread_offset = ((subread_step * subread_no) >> 16);
					#define SHIFT_SUBREAD_INT(ii, pp) { int nch = read_text [pp+read_text_rev_offset]; ii = (ii << 2) | base2int( nch );}
					#define BUILD_RBIN  {  int new2b = subread_integer & 3;\
						int rbin_byte = building_rbin_offset + (last_vote_rpos +16)/4;\
						int rbin_bit =(last_vote_rpos +16)%4 *2;\
						if(rbin_bit ==0) read_bin[rbin_byte]=0;\
						read_bin[rbin_byte] |= new2b<<rbin_bit;  }
					if(task==STEP_JUNC_TABLE){
						for(; last_vote_rpos  < subread_offset ; last_vote_rpos ++)
							SHIFT_SUBREAD_INT(subread_integer , last_vote_rpos  +16); // read_bin is not used in junction-detection mode.
					}else{
						for(; last_vote_rpos  < subread_offset ; last_vote_rpos ++){
							SHIFT_SUBREAD_INT(subread_integer , last_vote_rpos  +16);
							BUILD_RBIN;
						}
					}
					prefill_votes(cct_context->current_index, &prefill_ptrs, applied_subreads, subread_integer, subread_offset, subread_no, is_reversed);
				}

				if(last_vote_rpos > read_len - 16)SUBREADprintf("ERROR: exceeded offset %d > %d\n", last_vote_rpos , read_len - 16);

				if(task!=STEP_JUNC_TABLE) for(; last_vote_rpos  < read_len - 16 ; last_vote_rpos ++){
					SHIFT_SUBREAD_INT(subread_integer , last_vote_rpos  +16);
					BUILD_RBIN;
				} // SHIFT_SUBREAD_INT here is only for build read_bin, which isn't used on junction detection mode.

				if(is_reversed) {
					cellCounts_process_copy_ptrs_to_votes(cct_context, thread_no, &prefill_ptrs, vote_me, applied_subreads, read_name);
	#ifdef __MINGW32__
					if(current_read_number % 1000000 == 0 && current_read_number>0) print_in_box(80,0,0,"  Mapped : % 13" PRId64 " reads; time elapsed : % 5.1f mins\n", cct_context -> all_processed_reads_before_chunk + current_read_number, ( - cct_context -> program_start_time + miltime() ) / 60.);
	#else
					if(current_read_number % 1000000 == 0 && current_read_number>0) print_in_box(80,0,0,"  Mapped : % 13lld reads; time elapsed : % 5.1f mins\n", cct_context -> all_processed_reads_before_chunk + current_read_number, ( - cct_context -> program_start_time + miltime() ) / 60.);
	#endif
					
	//fprintf(stderr,"RNAME=%s\n", read_name);
					if(task==STEP_VOTING)
						cellCounts_select_and_write_alignments(cct_context, thread_no, sample_i, vote_me, read_name, read_text, read_bin, qual_text, read_len, applied_subreads);
					if(task==STEP_JUNC_TABLE){
						cellCounts_call_juncs_put_in_tab(cct_context, thread_no, sample_i, vote_me, read_name, read_text, /*read_bin -- not used for junction detection */ NULL, qual_text, read_len, applied_subreads);
						cellCounts_select_and_write_temps(cct_context, thread_no, sample_i, vote_me, read_name, read_text, /*read_bin -- not used for junction detection */ NULL, qual_text, read_len, applied_subreads);
					}
				} else {
					building_rbin_offset = REVERSED_READ_BIN_OFFSET;
					read_text_rev_offset = MAX_SCRNA_READ_LENGTH+1;
					strcpy(read_text+read_text_rev_offset, read_text);
					reverse_read(read_text+read_text_rev_offset, read_len, GENE_SPACE_BASE);
					qual_text[read_text_rev_offset] = 0;
				}
			}
		}else if(task==STEP_JUNC_TABLE){
			cellCounts_select_and_write_temps(cct_context, thread_no, -1, NULL, read_name, read_text, /* read_bin  -- read_bin is the 2-bit encoded read sequence and is not used for junction detection*/ NULL, qual_text, read_len, -1);
		}else if(task==STEP_VOTING){ // junction-detection mode doens't have the voting step.
			cellCounts_select_and_write_alignments(cct_context, thread_no, -1, NULL, read_name, read_text, read_bin, qual_text, read_len, -1);
		}
			// if the read belongs to "unassigned", then it doesn't need to be mapped at all (nowhere to write BAM).
	}

	free(vote_me);
	free(read_text);
	free(qual_text);

	return cct_context -> has_error;
}


int cellCounts_do_junctable(cellcounts_global_t * cct_context, int thread_no) {
	return cellCounts_do_jtab_or_voting(cct_context, thread_no, STEP_JUNC_TABLE);
}
int OLDcellCounts_do_junctable(cellcounts_global_t * cct_context, int thread_no) {
	subread_read_number_t current_read_number=0;
	char * read_text, * qual_text;
	char read_name[MAX_READ_NAME_LEN+1];
	int read_len=0;

	read_text = malloc(MAX_SCRNA_READ_LENGTH * 2+2);
	qual_text = malloc(MAX_SCRNA_READ_LENGTH * 2+2);

	int preads=0;
	while(!cct_context -> has_error) {
		cellCounts_fetch_next_read_pair(cct_context, thread_no,  &read_len, read_name, read_text, qual_text, &current_read_number);
		if(current_read_number < 0) break;
		if(read_len< 16) continue;
		preads++;
	}
	return cct_context -> has_error;
}

int cellCounts_do_voting(cellcounts_global_t * cct_context, int thread_no) {
	return cellCounts_do_jtab_or_voting(cct_context, thread_no, STEP_VOTING);
}

#define MAKE_SUBREAD_OFFSET	if(subread_no == applied_subreads -1) subread_offset= read_len-16; else subread_offset= ((subread_step * subread_no) >> 16);

#define MAKE_SUBREAD_INTVAL	subread_integer =0; for(xk1 = 0; xk1 < 16; xk1++){\
					int rbin_byte = (xk1+subread_offset)/4 + is_reversed * REVERSED_READ_BIN_OFFSET;\
					int rbin_bit = (xk1+subread_offset) %4 *2;\
					unsigned int vtmp = ( read_bin[ rbin_byte ] >> rbin_bit )&3;\
					subread_integer |= vtmp<<(2*(15-xk1)); }

void cellCounts_absoffset_to_posstr(cellcounts_global_t * cct_context, unsigned int pos, char * res){
	char * ch;
	int off;
	locate_gene_position(pos, &cct_context -> chromosome_table, &  ch, &off);
	SUBreadSprintf(res, 100 , "%s:%u", ch, off);
}

void known_pointer_strcat(char * targ, char * src, char ** buf){
	int srclen = strlen(src);
	if( (*buf) == NULL){
		(*buf) = targ;
	}
	memcpy((*buf), src, srclen);
	(*buf) += srclen;
	(**buf) = 0;
}

int cellCounts_write_gene_list(cellcounts_global_t * cct_context){
	int xk1;
	char ofname[MAX_FILE_NAME_LENGTH + 100];
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.Annot",cct_context->output_prefix);
	FILE * fp_out = fopen( ofname , "w" );
	fprintf(fp_out,"GeneID\tChr\tStart\tEnd\tStrand\tLength\n");

	unsigned int * gene_exons_number = calloc(sizeof(unsigned int) , cct_context -> gene_name_table -> numOfElements);
	unsigned int * gene_exons_pointer = calloc(sizeof(unsigned int) , cct_context -> gene_name_table -> numOfElements);
	unsigned int * gene_exons_start = malloc(sizeof(unsigned int) * cct_context -> all_features_array -> numOfElements);
	unsigned int * gene_exons_end = malloc(sizeof(unsigned int) * cct_context -> all_features_array -> numOfElements);
	char ** gene_exons_chr = malloc(sizeof(char *) * cct_context -> all_features_array -> numOfElements);
	char * gene_exons_strand = malloc(cct_context -> all_features_array -> numOfElements);

	for(xk1 = 0; xk1 < cct_context -> all_features_array -> numOfElements; xk1++) {
		int gene_id = cct_context -> features_sorted_geneid[xk1];
		gene_exons_number[gene_id]++;
	}

	unsigned int accumulative_no = 0;
	unsigned longest_gene_exons = 0;
	for(xk1 = 0 ; xk1 < cct_context -> gene_name_table -> numOfElements; xk1++) {
		unsigned int this_gene_exons = gene_exons_number[xk1];
		longest_gene_exons = max(longest_gene_exons, this_gene_exons);
		gene_exons_number[xk1] = accumulative_no;
		accumulative_no += this_gene_exons;
	}

	for(xk1 = 0; xk1 < cct_context -> all_features_array -> numOfElements; xk1++) {
		int gene_id = cct_context -> features_sorted_geneid[xk1];
		int gene_write_ptr = gene_exons_number[gene_id] + gene_exons_pointer[gene_id];

		gene_exons_chr[gene_write_ptr] = cct_context -> features_sorted_chr[xk1];
		gene_exons_start[gene_write_ptr] = cct_context -> features_sorted_start[xk1]; 
		gene_exons_end[gene_write_ptr] = cct_context -> features_sorted_stop[xk1]; 
		gene_exons_strand[gene_write_ptr] = cct_context -> features_sorted_strand[xk1]; 

		gene_exons_pointer[gene_id]++;
	}

	char *is_occupied = malloc(longest_gene_exons);
	unsigned int * input_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);
	unsigned int * output_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);
	int disk_is_full = 0;

	char * out_chr_list = malloc(longest_gene_exons * (1+cct_context -> longest_chro_name) + 1), * tmp_chr_list = NULL;
	char * out_start_list = malloc(11 * longest_gene_exons + 1), * tmp_start_list = NULL;
	char * out_end_list = malloc(11 * longest_gene_exons + 1), * tmp_end_list = NULL;
	char * out_strand_list = malloc(2 * longest_gene_exons + 1), * tmp_strand_list = NULL;

	for(xk1 = 0 ; xk1 < cct_context -> gene_name_table -> numOfElements; xk1++) {
		int xk2;
		
		memset(is_occupied,0,gene_exons_pointer[xk1]);
		tmp_chr_list = NULL;
		tmp_start_list = NULL;
		tmp_end_list = NULL;
		tmp_strand_list = NULL;
		out_chr_list[0]=0;
		out_start_list[0]=0;
		out_end_list[0]=0;
		out_strand_list[0]=0;
		int gene_nonoverlap_len =0;

		unsigned char * gene_symbol = cct_context -> gene_name_array [xk1];
		for(xk2=0; xk2<gene_exons_pointer[xk1]; xk2++) {
			if(!is_occupied[xk2]) {
				int xk3;
				char * matched_chr = gene_exons_chr[xk2 + gene_exons_number[xk1]];
				char matched_strand = gene_exons_strand[xk2 + gene_exons_number[xk1]];

				memset(input_start_stop_list, 0, gene_exons_pointer[xk1] * sizeof(int) * 2);
				int gap_merge_ptr = 1;
				input_start_stop_list[0] = gene_exons_start[xk2 + gene_exons_number[xk1]];
				input_start_stop_list[1] = gene_exons_end[xk2 + gene_exons_number[xk1]] + 1;

				for(xk3 = xk2; xk3 < gene_exons_pointer[xk1]; xk3++)
				{
					if(xk3==xk2)continue;

					if((!is_occupied[xk3]) && strcmp(matched_chr, gene_exons_chr[xk3+gene_exons_number[xk1]])==0 && matched_strand == gene_exons_strand[xk3 + gene_exons_number[xk1]])
					{
						is_occupied[xk3]=1;
						input_start_stop_list[gap_merge_ptr*2] = gene_exons_start[xk3+gene_exons_number[xk1]]; 
						input_start_stop_list[gap_merge_ptr*2+1] = gene_exons_end[xk3+gene_exons_number[xk1]]+1;

						gap_merge_ptr++;
					}
				}

				{
						int merged_gaps = mergeIntervals(input_start_stop_list, output_start_stop_list, gap_merge_ptr);

						for(xk3=0; xk3<gap_merge_ptr; xk3++)
						{
							char numbbuf[12];
							known_pointer_strcat(out_chr_list, matched_chr, &tmp_chr_list);
							known_pointer_strcat(out_chr_list, ";", &tmp_chr_list);

							SUBreadSprintf(numbbuf,12,"%u;", input_start_stop_list[xk3 * 2]);
							known_pointer_strcat(out_start_list, numbbuf, &tmp_start_list);
							SUBreadSprintf(numbbuf,12,"%u;", input_start_stop_list[xk3 * 2 + 1] - 1);
							known_pointer_strcat(out_end_list, numbbuf, &tmp_end_list);
							SUBreadSprintf(numbbuf,12,"%c;", (matched_strand==1)?'-':( ( matched_strand==0 )? '+':'.'));
							known_pointer_strcat(out_strand_list, numbbuf, &tmp_strand_list);

						}
						for(xk3=0; xk3<merged_gaps; xk3++)
							gene_nonoverlap_len += output_start_stop_list[xk3 * 2 + 1] - output_start_stop_list[xk3 * 2];
				}
			}
		}
		#define _cut_tail(x) (x)[strlen(x)-1]=0

		_cut_tail(out_chr_list);
		_cut_tail(out_start_list);
		_cut_tail(out_end_list);
		_cut_tail(out_strand_list);

		int wlen = fprintf(fp_out, "%s\t%s\t%s\t%s\t%s\t%d\n", gene_symbol, out_chr_list, out_start_list, out_end_list, out_strand_list, gene_nonoverlap_len);
		if(wlen < 6)disk_is_full = 1;
	}

	free(is_occupied);
	free(input_start_stop_list);
	free(output_start_stop_list);
	free(out_chr_list);
	free(out_strand_list);
	free(out_start_list);
	free(out_end_list);

	free(gene_exons_number);
	free(gene_exons_pointer);
	free(gene_exons_chr);
	free(gene_exons_start);
	free(gene_exons_end);
	free(gene_exons_strand);
	fclose(fp_out);

	if(disk_is_full){
		SUBREADprintf("ERROR: disk is full; the count file cannot be generated.\n");
		unlink(ofname);
		return -1;
	}
	return 0;
}



int cellCounts_run_mapping(cellcounts_global_t * cct_context){
	int chunk_no = 0;

	cct_context -> current_index = (gehash_t*) malloc(sizeof(gehash_t));
	sread_len = 0;

	if(1){ // Only load index once. No split index is supported.
		char tmp_fname[MAX_FILE_NAME_LENGTH+30];
		SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+30, "%s.%02d.b.tab", cct_context->index_prefix, cct_context->current_index_block_number);
		print_in_box(80,0,0, "Load the %d-%s index block...",1+ cct_context->current_index_block_number, cct_context->current_index_block_number==0?"st":(cct_context->current_index_block_number==1?"nd":"th"));
		if(gehash_load(cct_context -> current_index, tmp_fname)) return -1;
		print_in_box(80,0,0, "The index block has been loaded. Now map the reads...");
		print_in_box(80,0,0, "");
		//SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+30, "%s.%02d.b.array", cct_context->index_prefix, cct_context->current_index_block_number);
	}

	int main_step;
	for(main_step=0; main_step<2; main_step++){
		if(0==main_step && !cct_context -> do_cell_level_junction_detection)continue;
		cct_context -> all_processed_reads_before_chunk = 0;
		cct_context -> running_processed_reads_in_chunk=0;
		cct_context -> processed_reads_in_chunk=0;

		while(1) {
			int ret = 0;
			for(cct_context->current_index_block_number = 0; cct_context->current_index_block_number < cct_context->total_index_blocks; cct_context->current_index_block_number++) {
				if(cct_context->total_index_blocks == cct_context->current_index_block_number + 1)
					cct_context -> is_final_voting_run = 1;
				else	cct_context -> is_final_voting_run = 0;
				
				cct_context -> processed_reads_in_chunk = cct_context -> running_processed_reads_in_chunk;
				int is_last_chunk = cct_context -> processed_reads_in_chunk < cct_context-> reads_per_chunk;
				
				if(main_step && cct_context -> do_cell_level_junction_detection)cellCounts_do_realign(cct_context);
					 // because there are many input files, the thread control function is different
				else ret = cellCounts_run_maybe_threads(cct_context, main_step?STEP_VOTING:STEP_JUNC_TABLE);

//				if(cct_context->total_index_blocks > 1 || is_last_chunk)
				
				if(ret) break;
				if(!cct_context -> processed_reads_in_chunk) break;
			}

			cct_context -> all_processed_reads_before_chunk += cct_context -> processed_reads_in_chunk ;

			if(ret) return ret;

			if(1 || cct_context -> processed_reads_in_chunk < cct_context -> reads_per_chunk ||
			  (cct_context -> output_binfiles_are_full))
				// There will not be "chunks" for read processing. All reads are processed in one block in each pass.
				break;

			*(int*)(0x0) = 0x12345678;// this will NEVER be reached.

			cellCounts_go_chunk_nextchunk(cct_context);
			cellCounts_clean_context_after_chunk(cct_context);
			chunk_no++;
		}
		if(0==main_step){
			geinput_close(&cct_context -> input_dataset);
			cellCounts_open_input_fps(cct_context);
		}
	}
	gehash_destory_fast(cct_context -> current_index);

	free(cct_context -> current_index);
	return 0;
}

#define CELLCOUNTS_BAMBLOCK_SIZE 60000
#define CELLCOUNTS_BAMBLOCK_COMP_NUMBER 1
#define CELLRANGER_MERGER_WORKER_BINSIZE 62000

struct scRNA_merge_batches_worker_task{
	int sample_id;
	int inbin_len;
	int inbin_number;
	int inbin_batch_start_offsets [CELLCOUNTS_BAMBLOCK_COMP_NUMBER ];
	srInt_64 block_number;
	char inbin[(READ_BIN_BUF_SIZE+CELLCOUNTS_BAMBLOCK_SIZE )*CELLCOUNTS_BAMBLOCK_COMP_NUMBER ];
};

struct scRNA_merge_batches_worker_current{
	struct scRNA_merge_batches_worker_task * task;
	char outbin[CELLRANGER_MERGER_WORKER_BINSIZE * CELLCOUNTS_BAMBLOCK_COMP_NUMBER ];
	int outbin_len[CELLCOUNTS_BAMBLOCK_COMP_NUMBER ];
	unsigned int crc32[CELLCOUNTS_BAMBLOCK_COMP_NUMBER ];

	z_stream strm;
};

struct cell_gene_umi_supp{
	int cellbc;
	srInt_64 gene_no;
	char umi[MAX_UMI_LEN];
	int supp_reads;
};

int cellCounts_hamming_max2_fixlen(char * u1, char * u2, int ulen){
	int x, ret=0;
	for(x=0; x<ulen; x++){
		if(u1[x]!=u2[x]) ret++;
		if(ret>1)return ret;
	}
	return ret;
}

#define ADD_count_hash(bc,gn,no)  { HashTablePut(cellBCp0_genep0_P1_to_UMIs, NULL +1+(((1LLU*(bc))<<32)| (gn) ),  HashTableGet(   cellBCp0_genep0_P1_to_UMIs, NULL +1+(((1LLU*(bc))<<32)| (gn))) +(no) );\
    if( cct_context -> read_assignment_detail_fp ){\
        cellCounts_lock_occupy(&cct_context -> read_assignment_detail_lock);\
        fprintf( cct_context -> read_assignment_detail_fp,  "UMI_FINALLY_ASSIGN\tSAMPLE%03d\t%s\t%s\t%s\n",sample_no, (char*)ArrayListGet(cct_context -> cell_barcodes_array, bc),  str1 -> umi, cct_context ->gene_name_array[gn]);\
        cellCounts_lock_release(&cct_context -> read_assignment_detail_lock);\
       }\
    }
void cellCounts_do_one_batch_UMI_merge_one_cell(ArrayList* structs, int sec_start, int sec_end, int is_UMI_step2, HashTable * filtered_CGU_table, srInt_64 * remove_count, int sample_no){
	int x1;
	void ** app1 = structs -> appendix1;
	cellcounts_global_t * cct_context = app1[0];
	HashTable * cellBCp0_genep0_P1_to_UMIs = app1[2];
	int sample_id = app1[3]-NULL;

	if(is_UMI_step2){
		// NB: when this function is called, sec_end - sec_start MUST be >=2.
		for(x1 = sec_start; x1<sec_end; x1++) {
			struct cell_gene_umi_supp * str1 = ArrayListGet(structs, x1);
			if(x1 == sec_start){
				struct cell_gene_umi_supp * str2 = ArrayListGet(structs, sec_start+1);
				if(str1 -> supp_reads > str2 -> supp_reads){
					ADD_count_hash(str1->cellbc, str1->gene_no,1);
					continue;
				} else if(remove_count)(*remove_count)++;
			}


			char replaced_key[55+MAX_UMI_LEN];
#ifdef __MINGW32__
			//int keyptr = SUBreadSprintf(replaced_key, 40+MAX_UMI_LEN,"%d-%" PRId64 "-", str1 -> cellbc, str1 -> gene_no);
			int keyptr = SUBreadSprintf(replaced_key, 55+MAX_UMI_LEN,"%d-%d-%" PRId64 "-", sample_id, str1 -> cellbc, str1 -> gene_no);
#else
			//int keyptr = SUBreadSprintf(replaced_key, 40+MAX_UMI_LEN,"%d-%lld-", str1 -> cellbc, str1 -> gene_no);
			int keyptr = SUBreadSprintf(replaced_key, 55+MAX_UMI_LEN,"%d-%d-%lld-", sample_id, str1 -> cellbc, str1 -> gene_no);
#endif

			memcpy(replaced_key+keyptr, str1 -> umi, cct_context -> UMI_length);
			replaced_key[keyptr+cct_context -> UMI_length]=0;

			HashTablePut(filtered_CGU_table, strdup(replaced_key), NULL-1);

			str1 -> cellbc = -1;
		}
	}else{
		ArrayList * accepted_list =NULL;
		HashTable * looktable = NULL;
		int n_cutoff_looktab = 30;
#ifdef __DEBUG_NO_LOOK
		#warning "=====  DISABLED TABLE_BASED BATCHING !!! ====="
		n_cutoff_looktab = 0x3fffffff;
#endif
		if(sec_end - sec_start > n_cutoff_looktab){
			looktable = StringTableCreate((sec_end - sec_start)/5);
			HashTableSetDeallocationFunctions(looktable, free, (void (*)(void *value))ArrayListDestroy);
		}else accepted_list = ArrayListCreate(sec_end - sec_start);

		for(x1=sec_start; x1<sec_end; x1++){
			struct cell_gene_umi_supp * try_str = ArrayListGet(structs , x1);
			int x2, found = 0;
			ArrayList * test_accs;
			int hx;

			if(looktable){
				for(hx = 0; hx<2; hx++){
					char test_ky[MAX_UMI_LEN];
					test_ky[0] = hx?'S':'F';
					memcpy(test_ky +1, try_str -> umi + hx * cct_context -> UMI_length/2 , cct_context -> UMI_length/2);
					test_ky[1+cct_context -> UMI_length/2]=0;

					test_accs = HashTableGet(looktable, test_ky);
					if(!test_accs)continue;

					for(x2=0; x2<test_accs->numOfElements; x2++){
						struct cell_gene_umi_supp * acc_str = ArrayListGet(test_accs, x2);
						if(cellCounts_hamming_max2_fixlen(acc_str -> umi, try_str -> umi, cct_context -> UMI_length)<2){
							found=1;

							if(cct_context -> read_assignment_detail_fp){
								cellCounts_lock_occupy(&cct_context -> read_assignment_detail_lock);
								fprintf( cct_context -> read_assignment_detail_fp,  "FIX_UMI_SIMILAR\tSAMPLE%03d\t%s\t%s\t%s\t%s\n",sample_no, (char*)ArrayListGet(cct_context -> cell_barcodes_array, acc_str -> cellbc),  try_str -> umi, acc_str -> umi, cct_context ->gene_name_array[ acc_str -> gene_no ]);
								cellCounts_lock_release(&cct_context -> read_assignment_detail_lock);
							}
							acc_str -> supp_reads += try_str -> supp_reads;

							char replaced_key[55+MAX_UMI_LEN];
#ifdef __MINGW32__
							int keyptr = SUBreadSprintf(replaced_key, 55+MAX_UMI_LEN,"%d-%d-%" PRId64 "-", sample_id, try_str -> cellbc, try_str -> gene_no);
#else
							int keyptr = SUBreadSprintf(replaced_key, 55+MAX_UMI_LEN,"%d-%d-%lld-", sample_id, try_str -> cellbc, try_str -> gene_no);
#endif

							memcpy(replaced_key+keyptr, try_str -> umi, cct_context -> UMI_length);
							replaced_key[keyptr+cct_context -> UMI_length]=0;
							HashTablePut(filtered_CGU_table, strdup(replaced_key), acc_str -> umi);
							try_str -> cellbc = -1;
							break;
						}
					}
					if(found)break;
				}
			}else{
				test_accs = accepted_list;

				for(x2=0; x2<test_accs->numOfElements; x2++){
					struct cell_gene_umi_supp * acc_str = ArrayListGet(test_accs, x2);
					if(cellCounts_hamming_max2_fixlen(acc_str -> umi, try_str -> umi, cct_context -> UMI_length)<2){
						if(cct_context -> read_assignment_detail_fp){
							cellCounts_lock_occupy(&cct_context -> read_assignment_detail_lock);
							fprintf( cct_context -> read_assignment_detail_fp,  "FIX_UMI_SIMILAR\tSAMPLE%03d\t%s\t%s\t%s\t%s\n",sample_no, (char*)ArrayListGet(cct_context -> cell_barcodes_array, acc_str -> cellbc),  try_str -> umi, acc_str -> umi, cct_context ->gene_name_array[ acc_str -> gene_no ]);
							cellCounts_lock_release(&cct_context -> read_assignment_detail_lock);
						}
						found=1;
						acc_str -> supp_reads += try_str -> supp_reads;

						char replaced_key[55+MAX_UMI_LEN];
#ifdef __MINGW32__
						int keyptr = SUBreadSprintf(replaced_key, 55+MAX_UMI_LEN,"%d-%d-%" PRId64 "-", sample_id, try_str -> cellbc, try_str -> gene_no);
#else
						int keyptr = SUBreadSprintf(replaced_key, 55+MAX_UMI_LEN,"%d-%d-%lld-", sample_id, try_str -> cellbc, try_str -> gene_no);
#endif
						memcpy(replaced_key+keyptr, try_str -> umi, cct_context -> UMI_length);
						replaced_key[keyptr+cct_context -> UMI_length]=0;
						HashTablePut(filtered_CGU_table, strdup(replaced_key), acc_str -> umi);
						try_str -> cellbc = -1;
						break;
					}
				}
			}
			if(!found){
				if(looktable){
					for(hx = 0; hx<2; hx++){
						char test_ky[MAX_UMI_LEN];
						test_ky[0] = hx?'S':'F';
						memcpy(test_ky +1, try_str -> umi + hx * cct_context -> UMI_length/2 , cct_context -> UMI_length/2);
						test_ky[1+cct_context -> UMI_length/2]=0;
						test_accs = HashTableGet(looktable, test_ky);
						if(!test_accs){
							test_accs = ArrayListCreate(10);
							HashTablePut(looktable, strdup(test_ky), test_accs);
						}
						ArrayListPush(test_accs, try_str);
					}
				}else ArrayListPush(accepted_list, try_str);
			}
		}

		if(looktable)HashTableDestroy(looktable);
		else ArrayListDestroy(accepted_list);
	}
}

void cellCounts_do_one_batch_UMI_merge_one_step(ArrayList* structs, int is_UMI_step2, HashTable * filtered_CGU_table, srInt_64 * remove_count, int sample_no){
	void ** app1 = structs -> appendix1;
	cellcounts_global_t * cct_context = app1[0];
	HashTable * cellBCp0_genep0_P1_to_UMIs = app1[2];
	srInt_64 x1, sec_start = 0;
	srInt_64 old_sec_key = -1;

	for(x1=0; x1<=structs -> numOfElements; x1++){
		srInt_64 sec_key = -1;
		int is_umi_changed = 0;

		struct cell_gene_umi_supp * str1 =NULL;
		if(x1<structs -> numOfElements){
			str1 = ArrayListGet(structs, x1);
			if(str1 -> cellbc <0) continue;
			sec_key = str1 -> cellbc;
			sec_key = sec_key << 32;
			if(is_UMI_step2 && sec_key == old_sec_key){
				struct cell_gene_umi_supp * strold = ArrayListGet(structs, sec_start);
				is_umi_changed = memcmp(strold -> umi, str1-> umi, cct_context-> UMI_length);
			}else if(!is_UMI_step2) sec_key = sec_key | str1 -> gene_no;
				// gene_no itself is 64-bit, but it is nearly impossible to have two neighbouring
				// structures that have the same last 32-bit of gene_no.
		}


		if( (x1>sec_start && sec_key!=old_sec_key) || is_umi_changed){ // when x1 == numOfElements, sec_key is -1. If old_sec_key is also -1, no item is included in the list. If old_sec_key is >=0, the last sec is processed.
			struct cell_gene_umi_supp * str1 = ArrayListGet(structs, sec_start);

			if(x1 - sec_start>1 && str1->cellbc>=0) cellCounts_do_one_batch_UMI_merge_one_cell(structs, sec_start, x1, is_UMI_step2, filtered_CGU_table, remove_count, sample_no);
			else if(is_UMI_step2 && str1->cellbc>=0) ADD_count_hash(str1->cellbc,str1->gene_no,1);

			sec_start = x1;
		}
		old_sec_key = sec_key;
	}
}


int cellCounts_do_one_batch_sort_compare(void * ar, int l, int r){
	void ** arr = ar;
	void ** bin_ptrs = arr[0];
	cellcounts_global_t * cct_context = arr[1];

	char * Lptr = bin_ptrs[l];
	char * Rptr = bin_ptrs[r];
	srInt_64 Lgenes=0, Rgenes=0;
	memcpy(&Lgenes, Lptr+8, 8);
	memcpy(&Rgenes, Rptr+8, 8);
	if(Lgenes & (1LLU<<63))Lgenes=Lgenes & 0x7fffffffllu; else Lgenes=0;
	if(Rgenes & (1LLU<<63))Rgenes=Rgenes & 0x7fffffffllu; else Rgenes=0;
	srInt_64 Lpos= ((0LLU+*(int*)(Lptr+16+Lgenes*8+cct_context->UMI_length+4))<<32) | *(unsigned int*)(Lptr+16+Lgenes*8+cct_context->UMI_length+4+4);
	srInt_64 Rpos= ((0LLU+*(int*)(Rptr+16+Rgenes*8+cct_context->UMI_length+4))<<32) | *(unsigned int*)(Rptr+16+Rgenes*8+cct_context->UMI_length+4+4);
	if(Lpos>Rpos)return 1;
	if(Lpos<Rpos)return -1;
	return 0;
}

void cellCounts_do_one_batch_sort_exchange(void * ar, int l, int r){
	void ** arr = ar;
	void ** bin_ptrs = arr[0];
	void * tp = bin_ptrs[l];
	bin_ptrs[l]=bin_ptrs[r];
	bin_ptrs[r]=tp;
}

void cellCounts_do_one_batch_sort_merge(void * ar, int start, int items, int items2){
	void ** arr = ar;
	void ** bin_ptrs = arr[0];
	bin_ptrs +=start;

	void ** tmp = malloc(sizeof(void*)*(items2+items));
	int i1_cursor=0, i2_cursor=items, wptr=0;
	while(1){
		if(i1_cursor == items && i2_cursor == items + items2 )break;
		int select_items_1 = (i2_cursor == items + items2) || (i1_cursor < items && cellCounts_do_one_batch_sort_compare(ar, start+ i1_cursor,start + i2_cursor) <= 0);
		if(select_items_1) tmp[wptr++] = bin_ptrs[i1_cursor++];
		else tmp[wptr++] = bin_ptrs[i2_cursor++];
	}
	memcpy(bin_ptrs, tmp, sizeof(void*)*(items2+items));
	free(tmp);
}

int cellCounts_do_one_batch_tab_to_struct_list_compare(void * L_elem, void * R_elem, ArrayList * me){
	struct cell_gene_umi_supp *L = L_elem, *R = R_elem;
	void ** app1 = me -> appendix1;
	cellcounts_global_t * cct_context = app1[0];
	int sort_by_geneid_then_umi = app1[1] - NULL;

	if(L->cellbc > R->cellbc) return 1;
	if(L->cellbc < R->cellbc) return -1;

	if(sort_by_geneid_then_umi){
		if(L->gene_no>R->gene_no) return 1;
		if(L->gene_no<R->gene_no) return -1;
	}else{
		int umicmps = memcmp(L->umi, R->umi, cct_context -> UMI_length);
		if(umicmps) return umicmps;
	}

	if(L->supp_reads < R->supp_reads) return 1;
	if(L->supp_reads > R->supp_reads) return -1; // reversed by # supp reads

	if(sort_by_geneid_then_umi){
		int umicmps = memcmp(L->umi, R->umi, cct_context -> UMI_length);
		if(umicmps) return umicmps;
	}else{

		if(L->gene_no>R->gene_no) return 1;
		if(L->gene_no<R->gene_no) return -1;
	}
	return 0;
}

void cellCounts_do_one_batch_tab_to_struct_list(void *ky, void *val, HashTable * tab){
	int supp_reads = val-NULL;
	ArrayList ** cell_gene_umi_list = tab -> appendix1;
	int UMI_length = tab -> counter1;

	struct cell_gene_umi_supp * new_item = malloc(sizeof(struct cell_gene_umi_supp));
	char * kyptr = ky;
	int sample_id = atoi(kyptr); // one-based sample id
	for(; '-' != *kyptr; kyptr++);
	kyptr++;
	new_item -> cellbc = atoi(kyptr);
	for(; '-' != *kyptr; kyptr++);
	kyptr++;
	new_item -> gene_no = atoll(kyptr);
	for(; '-' != *kyptr; kyptr++);
	memcpy(new_item->umi, kyptr+1, UMI_length);
	new_item -> supp_reads = supp_reads;
	if(sample_id<1)SUBREADprintf("WRONG SAMPLE ID: %d from '%s'\n", sample_id, (char*)ky);
	ArrayListPush(cell_gene_umi_list[sample_id-1], new_item);
}

void cellCounts_do_one_batch_write_UMIs(void * vcell_gene, void * vumis, HashTable * me){
	REPFILE * fp = me->appendix1;
	vcell_gene --;
	REP_fwrite(&vcell_gene,1,8,fp);
	REP_fwrite(&vumis,1,8,fp);
}


void * cellCounts_merge_batches_worker(void * vp){
	void **vpp = vp;
	cellcounts_global_t * cct_context = vpp[0];
	worker_master_mutex_t * worker_mut  = vpp[1];
	int my_worker_id = vpp[2] - NULL;
	struct scRNA_merge_batches_worker_current * my_current_job = vpp[3];
	free(vp);

	int Z_DEFAULT_MEM_LEVEL = 8;
	worker_thread_start(worker_mut, my_worker_id);
	while(1){
		if(worker_wait_for_job(worker_mut, my_worker_id)) break;
		if(!cct_context -> is_BAM_and_FQ_out_generated) continue;


		struct scRNA_merge_batches_worker_task * current_input = my_current_job -> task;
		int current_blk;
		for(current_blk =0; current_blk < current_input -> inbin_number ; current_blk ++){
			char * inbin_blk = current_input -> inbin + current_input -> inbin_batch_start_offsets[current_blk ];
			int inblock_size = -1;
			if(current_blk == current_input -> inbin_number -1) inblock_size= current_input -> inbin_len - current_input -> inbin_batch_start_offsets[current_blk ];
			else if(CELLCOUNTS_BAMBLOCK_COMP_NUMBER>1) inblock_size= current_input  -> inbin_batch_start_offsets[current_blk +1] - current_input -> inbin_batch_start_offsets[current_blk ];

			deflateInit2(&my_current_job -> strm , SAMBAM_COMPRESS_LEVEL_NORMAL, Z_DEFLATED, SAMBAM_GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
			my_current_job -> strm.avail_in = inblock_size ;
			my_current_job -> strm.next_in = (unsigned char*)inbin_blk ;
			my_current_job -> strm.avail_out = CELLRANGER_MERGER_WORKER_BINSIZE;
			my_current_job -> strm.next_out = (unsigned char*)(my_current_job -> outbin + CELLRANGER_MERGER_WORKER_BINSIZE * current_blk);

			deflate(&my_current_job -> strm, Z_FINISH);
			my_current_job -> outbin_len [current_blk] = CELLRANGER_MERGER_WORKER_BINSIZE-my_current_job -> strm.avail_out;
			my_current_job -> crc32 [current_blk] = SamBam_CRC32(inbin_blk, inblock_size);
			deflateEnd(&my_current_job -> strm);
		}
	}
	return NULL;
}

int cellCounts_make_barcode_bam_bin(cellcounts_global_t * cct_context, char * rbin, char * new_rbin, int binlen, char * fixedbc_seq, char * fixedumi_seq, srInt_64 gene_no, srInt_64 * genes) {
	char * cellbc_seq=NULL,*umi_seq=NULL, * cellbc_qual=NULL,*umi_qual=NULL, *sample_seq=NULL, *sample_qual=NULL, *lane_str=NULL;
	int rname_trimmed_len=0;
	cellCounts_scan_read_name_str(cct_context, rbin, NULL, & sample_seq, & sample_qual, & cellbc_seq, & cellbc_qual, & umi_seq, & umi_qual, &lane_str, NULL, &rname_trimmed_len);
	if( cct_context->visium_hd_barcodes ) {
		umi_seq = cellbc_seq;
		umi_qual = cellbc_qual;
	}
	int new_rbin_len = 0, n_cigar_op =0, l_read_name=0, l_seq=0;

	memcpy(new_rbin, rbin, 36);
	new_rbin_len += 36;

	memcpy(&n_cigar_op, rbin+16,2);
	memcpy(&l_seq, rbin+20,4);
	l_read_name=((unsigned char*)rbin)[12];
	new_rbin[12] = rname_trimmed_len+1;
	memcpy(new_rbin+new_rbin_len, rbin+36, rname_trimmed_len);
	new_rbin[36+rname_trimmed_len]=0;
	new_rbin_len+= rname_trimmed_len+1;
	memcpy(new_rbin+new_rbin_len, rbin +36 + l_read_name, 4*n_cigar_op + l_seq + (l_seq+1)/2);
	new_rbin_len += 4*n_cigar_op + l_seq + (l_seq+1)/2;
	char * ext_bin_ptr = rbin + 36 + l_read_name +4*n_cigar_op + l_seq + (l_seq+1)/2;

#ifndef DO_STARSOLO_THING
	int CR_found=0, CB_found=0, CY_found=0, UR_found=0, UY_found=0, UB_found=0;
	int cellbc_content_len = cct_context -> known_cell_barcode_length;
	int cellbc_content_len_raw = cellbc_content_len;
	if(cct_context->visium_hd_barcodes){
		cellbc_content_len = 11; // 01234_01234
		int xx1, start_pos = -1;
		int end_pos = -1;
		for(xx1 = cct_context -> UMI_length; ; xx1++){
			int nch = cellbc_seq[xx1];
			if(nch >='a' && nch!='|'){
				nch -= 32;
				cellbc_seq[xx1] = nch;
				if(start_pos<0) start_pos = xx1;
				else end_pos = xx1+1;
			}
			if(nch=='|' || !nch) break;
		}
		cellbc_content_len_raw = end_pos - start_pos;
		cellbc_seq += start_pos;
		cellbc_qual += start_pos;
	}

	while(ext_bin_ptr < rbin+binlen+4){
		char * tagstr = NULL; 
		int taglen = 0;
		if(ext_bin_ptr[0]=='C' && ext_bin_ptr[1]=='R' && ext_bin_ptr[2]=='Z'){
			CR_found = 1;
			tagstr = cellbc_seq;
			taglen = cellbc_content_len;
		}else if(ext_bin_ptr[0]=='C' && ext_bin_ptr[1]=='B' && ext_bin_ptr[2]=='Z'){
			CB_found = 1;
			tagstr = fixedbc_seq;
			taglen = cellbc_content_len;
		}else if(ext_bin_ptr[0]=='C' && ext_bin_ptr[1]=='Y' && ext_bin_ptr[2]=='Z'){
			CY_found = 1;
			tagstr = cellbc_qual;
			taglen = cellbc_content_len;
		}else if(ext_bin_ptr[0]=='U' && ext_bin_ptr[1]=='R' && ext_bin_ptr[2]=='Z'){
			UR_found = 1;
			tagstr = umi_seq;
			taglen = cct_context -> UMI_length;
		}else if(ext_bin_ptr[0]=='U' && ext_bin_ptr[1]=='B' && ext_bin_ptr[2]=='Z'){
			UB_found = 1;
			tagstr = fixedumi_seq;
			taglen = cct_context -> UMI_length;
		}else if(ext_bin_ptr[0]=='U' && ext_bin_ptr[1]=='Y' && ext_bin_ptr[2]=='Z'){
			UY_found = 1;
			tagstr = umi_qual;
			taglen = cct_context -> UMI_length;
		}
	
		if(tagstr){
			new_rbin[new_rbin_len++]=*(ext_bin_ptr++);
			new_rbin[new_rbin_len++]=*(ext_bin_ptr++);
			new_rbin[new_rbin_len++]=*(ext_bin_ptr++);
			int taglenold = strlen(ext_bin_ptr);
			memcpy(new_rbin+new_rbin_len,tagstr, taglen);
			*(new_rbin+new_rbin_len+taglen)=0;
			ext_bin_ptr += taglenold+1;
			new_rbin_len += taglen+1;
		}else{
			int content_len = SAP_pairer_skip_tag_body_len(ext_bin_ptr);
			memcpy(new_rbin + new_rbin_len, ext_bin_ptr, content_len );
			new_rbin_len += content_len;
			ext_bin_ptr += content_len;
		}
	}
//	char * cellbc_out = cellbc_seq;

	int x2;
	if(cellbc_content_len_raw>0 && !CR_found){
		new_rbin[new_rbin_len++]='C';new_rbin[new_rbin_len++]='R';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, cellbc_seq, cellbc_content_len_raw);
		*(new_rbin+new_rbin_len+cellbc_content_len_raw)=0;
		new_rbin_len += cellbc_content_len_raw+1;
	}
	if(fixedbc_seq && !CB_found){
		new_rbin[new_rbin_len++]='C';new_rbin[new_rbin_len++]='B';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, fixedbc_seq, cellbc_content_len);
		*(new_rbin+new_rbin_len+cellbc_content_len)=0;
		new_rbin_len += cellbc_content_len+1;
	}
	if(cellbc_content_len_raw >0 && !CY_found){
		new_rbin[new_rbin_len++]='C';new_rbin[new_rbin_len++]='Y';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, cellbc_qual, cellbc_content_len_raw);
		for(x2=0; x2< cellbc_content_len_raw; x2++)if( new_rbin[ new_rbin_len+x2 ]>'/' ) new_rbin[ new_rbin_len+x2 ]--;
		*(new_rbin+new_rbin_len+cellbc_content_len_raw)=0;
		new_rbin_len += cellbc_content_len_raw+1;
	}

	if(!UR_found){
		new_rbin[new_rbin_len++]='U';new_rbin[new_rbin_len++]='R';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, umi_seq, cct_context -> UMI_length);
		*(new_rbin+new_rbin_len+cct_context -> UMI_length)=0;
		new_rbin_len += cct_context -> UMI_length+1;
	}
	if(fixedumi_seq && !UB_found){
		new_rbin[new_rbin_len++]='U';new_rbin[new_rbin_len++]='B';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, fixedumi_seq, cct_context -> UMI_length);
		*(new_rbin+new_rbin_len+cct_context -> UMI_length)=0;
		new_rbin_len += cct_context -> UMI_length+1;
	}
	if(!UY_found){
		new_rbin[new_rbin_len++]='U';new_rbin[new_rbin_len++]='Y';new_rbin[new_rbin_len++]='Z';
		memcpy(new_rbin+new_rbin_len, umi_qual, cct_context -> UMI_length);
		*(new_rbin+new_rbin_len+cct_context -> UMI_length)=0;

		for(x2=0; x2< cct_context -> UMI_length; x2++)if( new_rbin[ new_rbin_len+x2 ]>'/' ) new_rbin[ new_rbin_len+x2 ]--;
		new_rbin_len += cct_context -> UMI_length+1;
	}
	if(cct_context->visium_hd_barcodes){
		int x1,cbclen = strstr(umi_seq,"|")-umi_seq; // umi_seq had been set to start of the 1R.
		for(x1=0; x1<2; x1++){
			new_rbin[new_rbin_len++]='1';new_rbin[new_rbin_len++]=x1?'Y':'R';new_rbin[new_rbin_len++]='Z';
			memcpy(new_rbin+new_rbin_len, x1?umi_qual:umi_seq, cbclen);
			if(x1) for(x2=0; x2<cbclen; x2++)if( new_rbin[ new_rbin_len+x2 ]>'/' ) new_rbin[ new_rbin_len+x2 ]--;
			*(new_rbin+new_rbin_len+cbclen)=0;
			new_rbin_len += cbclen+1;
		}
	}

	new_rbin[new_rbin_len++]='X';new_rbin[new_rbin_len++]='Q';new_rbin[new_rbin_len++]='I';
	memcpy(new_rbin+new_rbin_len,&gene_no,4);
	new_rbin_len+=4;

	new_rbin[new_rbin_len++]='X';new_rbin[new_rbin_len++]='K';new_rbin[new_rbin_len++]='C';
	new_rbin[new_rbin_len++]=gene_no>>63;

#endif
	new_rbin_len-=4;
	memcpy(new_rbin, &new_rbin_len,4);
	return new_rbin_len;
}

void cellCounts_do_one_batch_write_extend_rbin(cellcounts_global_t * cct_context, char * rbin, int binlen, REPFILE * fp, char * fixedbc_seq, char * fixedumi_seq, srInt_64 gene_no, srInt_64 * genes){
	char new_rbin[ binlen + 150 ]; // removed barcodes/qual from read names, add them to extra fields if they weren't there. Gene names are not put here.
	int new_rbin_len = cellCounts_make_barcode_bam_bin( cct_context, rbin, new_rbin, binlen, fixedbc_seq, fixedumi_seq, gene_no, genes );
	REP_fwrite(new_rbin, 1, new_rbin_len+4, fp);
}

#ifdef __MINGW32__
#define ADD_key_FMT1 "%d-%d-%" PRId64 "-%s"
#else
#define ADD_key_FMT1 "%d-%d-%lld-%s"
#endif
#define ADD_key_struct { char my_key [50+MAX_UMI_LEN]; \
	SUBreadSprintf(my_key, 50+MAX_UMI_LEN,ADD_key_FMT1, sample_id, cell_no, gene_no, UMI_str); \
	srInt_64 supp_reads = HashTableGet(supp_reads_SCGU, my_key)-NULL; \
	if(1>supp_reads) HashTablePut(supp_reads_SCGU, strdup(my_key), NULL+1); \
	else HashTablePutReplaceEx(supp_reads_SCGU, my_key, NULL+supp_reads+1, 0,0,0); }

void * cellCounts_do_one_batch(void * paramsp1){
	srInt_64 x1;
	void ** params = paramsp1;
	cellcounts_global_t * cct_context = params[0];
	ArrayList * file_size_list = params[2];
	char *temp_dir = cct_context -> temp_file_dir;
	int thread_no = params[1]-NULL;
	free(paramsp1);

	cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + thread_no;

	int me_max_Rbin_len = 0;
	int me_max_genes = 0;
	char ** bin_ptrs = malloc(sizeof(char*) * 1500000), * batch_content=NULL;
	int bin_ptr_size = 1500000;
	srInt_64 removed_UMIs = 0;

	while(1){
		int this_batch_no = -1;
		cellCounts_lock_occupy(&cct_context -> input_dataset_lock);
		//cellCounts_lock_occupy(&cct_context -> input_dataset_lock);
		if(cct_context -> do_one_batch_runner_current < CELLBC_BATCH_NUMBER +1){
			int this_batch_sorted_idx = (cct_context -> do_one_batch_runner_current ++);
			srInt_64 this_batch_size_and_no = ArrayListGet(file_size_list, file_size_list->numOfElements-1 -this_batch_sorted_idx)-NULL;
			this_batch_no = (int)(this_batch_size_and_no&0xfffffllu);
		}
		if(me_max_genes > cct_context -> barcode_batched_max_genes) cct_context -> barcode_batched_max_genes = me_max_genes;
		if(me_max_Rbin_len > cct_context-> barcode_batched_max_Rbin_len) cct_context-> barcode_batched_max_Rbin_len = me_max_Rbin_len;
		cellCounts_lock_release(&cct_context -> input_dataset_lock);
		//cellCounts_lock_release(&cct_context -> input_dataset_lock);
		if(0>this_batch_no)break;
		char tmp_fname[MAX_FILE_NAME_LENGTH+80];
		SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+80, "%s/temp-cellcounts-%06d-%03d.tmpbin", temp_dir, getpid(), this_batch_no);
		REPFILE * fp = REP_fopen(tmp_fname, "rb");
		REP_setvbuf(fp, thread_context -> cellbin_v_buffer, _IOFBF , SCRNA_VBUFF_SIZE);

		srInt_64 batch_fsize = REP_filesize(fp);
		if(batch_content==NULL) batch_content = malloc(batch_fsize);
		srInt_64 batch_content_len = REP_fread(batch_content, 1, batch_fsize, fp);
		REP_fclose(fp);
		if(batch_content_len!=batch_fsize){
			SUBREADprintf("ERROR: Cannot load file at once: %d!\n", this_batch_no);
			return NULL;
		}

		HashTable * supp_reads_SCGU = StringTableCreate(500000);
		HashTableSetDeallocationFunctions(supp_reads_SCGU, free, NULL);
		srInt_64 scanptr = 0;
		int rbin_no = 0;
		char UMI_str[MAX_UMI_LEN+1];

		while(scanptr < batch_content_len-1){
			int cell_no=0, sample_id=0;
			srInt_64 gene_no=0;
			if(bin_ptr_size<=rbin_no){
				bin_ptr_size = bin_ptr_size*2;
				bin_ptrs = realloc(bin_ptrs, sizeof(char*)*bin_ptr_size);
			}
			bin_ptrs[rbin_no] = batch_content+scanptr;
			memcpy(&sample_id, batch_content+scanptr, 4);
			scanptr += 4; // sample_ID 
			memcpy(&cell_no, batch_content+scanptr, 4);
			scanptr += 4; // cellbarcode_NO
			memcpy(&gene_no, batch_content+scanptr, 8);
			scanptr += 8; // gene_id 
			if(gene_no & (1LLU<<63)){
				int genes = (int)(gene_no & 0x7fffffffllu);
				if(genes > me_max_genes)me_max_genes=genes;
				memcpy(UMI_str, batch_content+scanptr+8*genes, cct_context -> UMI_length);
				UMI_str[cct_context -> UMI_length]=0;

				for(x1=0; x1<genes; x1++){
					memcpy(&gene_no, batch_content+scanptr, 8);
					scanptr += 8;
					ADD_key_struct;
				}
			}else{
				UMI_str[cct_context -> UMI_length]=0;
				memcpy(UMI_str, batch_content+scanptr, cct_context -> UMI_length);
				ADD_key_struct;
			}
			scanptr += cct_context -> UMI_length ; // UMI str
			int rbinlen = 0;
			memcpy(&rbinlen, batch_content+scanptr, 4);
			if(me_max_Rbin_len < rbinlen) me_max_Rbin_len = rbinlen;
			scanptr += rbinlen +4; // read_bin
			rbin_no++;
		}
		ArrayList ** cell_gene_umi_list = malloc(sizeof(void*)*cct_context -> sample_sheet_table -> numOfElements);
		for(x1 =0; x1< cct_context -> sample_sheet_table -> numOfElements; x1++){
			cell_gene_umi_list[x1]=ArrayListCreate(2000000);
			ArrayListSetDeallocationFunction(cell_gene_umi_list[x1], free);
		}
		supp_reads_SCGU -> appendix1 = cell_gene_umi_list;
		supp_reads_SCGU -> appendix2 = cct_context;
		supp_reads_SCGU -> counter1 = cct_context -> UMI_length;
		HashTableIteration(supp_reads_SCGU, cellCounts_do_one_batch_tab_to_struct_list);
		HashTable * filtered_SCGU_table = StringTableCreate(max(10000,cell_gene_umi_list[0] -> numOfElements / 10));
		HashTableSetDeallocationFunctions(filtered_SCGU_table, free, NULL);

		fp = REP_fopen(tmp_fname, "wb");
		REP_setvbuf(fp, thread_context -> cellbin_v_buffer, _IOFBF , SCRNA_VBUFF_SIZE);
		for(x1 = 0; x1 < cct_context -> sample_sheet_table -> numOfElements; x1++){
			HashTable * cellbcP0_to_geneno0B_P1_to_UMIs = HashTableCreate(500000);

			void * app1[4];
			cell_gene_umi_list[x1] -> appendix1 = app1;
			app1[3] = NULL+x1+1; // sample_no
			app1[0] = cct_context;
			app1[1] = NULL+1;
						// 1 : sorted by cell_bc, then gene, then supported_reads, then UMIstr (this is for step1 UMI merging)
						// 0 : sorted by cell_bc, then UMIstr, then supported_reads, then gene (this is for step2 UMI merging)
						// supported_reads : large -> small; the other: small -> large
			ArrayListSort(cell_gene_umi_list[x1],  cellCounts_do_one_batch_tab_to_struct_list_compare);
			cellCounts_do_one_batch_UMI_merge_one_step(cell_gene_umi_list[x1], 0, filtered_SCGU_table, NULL, x1);

			app1[1] = NULL+0;
			app1[2] = cellbcP0_to_geneno0B_P1_to_UMIs;

			ArrayListSort(cell_gene_umi_list[x1], cellCounts_do_one_batch_tab_to_struct_list_compare);
			cellCounts_do_one_batch_UMI_merge_one_step(cell_gene_umi_list[x1], 1, filtered_SCGU_table, &removed_UMIs, x1);

			cellbcP0_to_geneno0B_P1_to_UMIs -> appendix1 = fp;

			REP_fwrite(&cellbcP0_to_geneno0B_P1_to_UMIs -> numOfElements,1,8,fp);
			HashTableIteration(cellbcP0_to_geneno0B_P1_to_UMIs, cellCounts_do_one_batch_write_UMIs);
			HashTableDestroy(cellbcP0_to_geneno0B_P1_to_UMIs);
		}

		void * sort_base[2];
		sort_base[0] = bin_ptrs;
		sort_base[1] = cct_context;
		merge_sort(sort_base, rbin_no, cellCounts_do_one_batch_sort_compare, cellCounts_do_one_batch_sort_exchange, cellCounts_do_one_batch_sort_merge);
		for(x1 = 0; x1 < rbin_no; x1++){
			char * binptr = bin_ptrs[x1];
			int cellid =0, sampleid = 0; // sampleid is 1-based
			srInt_64 gene_no =0, genes = 0, geneno_0 = 0;
			char * glist_ptr =NULL;
			memcpy(&sampleid, binptr, 4);
			memcpy(&cellid, binptr+4, 4);
			memcpy(&gene_no, binptr+8, 8);
			if(gene_no & (1LLU<<63)){
				glist_ptr =binptr + 16;
				genes = (int)(gene_no & 0x7fffffff);
				memcpy(&geneno_0, binptr+16, 8);
			}
			char * umi  = binptr + 16 + 8*genes;
			char * rbinptr = umi + cct_context ->UMI_length;

			int is_homopolymer_or_N_this_UR;
			if(cct_context->do_cell_level_junction_detection) is_homopolymer_or_N_this_UR = cellCounts_is_homopolymer_or_N(umi, cct_context->UMI_length);

			char SCGU_key [40+MAX_UMI_LEN];

			int remove_step;
#ifdef __MINGW32__
			int keyptr = SUBreadSprintf(SCGU_key, 40+MAX_UMI_LEN,"%d-%d-%" PRId64 "-", sampleid, cellid,  (gene_no & (1LLU<<63))? geneno_0: gene_no);
#else
			int keyptr = SUBreadSprintf(SCGU_key, 40+MAX_UMI_LEN,"%d-%d-%lld-", sampleid, cellid,  (gene_no & (1LLU<<63))? geneno_0: gene_no);
#endif
			SCGU_key[keyptr+ cct_context -> UMI_length] = 0;

			char * emptyUMI = "-----------------------------------------";
			for(remove_step =1; remove_step<=2; remove_step++){
				memcpy(SCGU_key+keyptr, umi, cct_context -> UMI_length);

				char * new_UMI = HashTableGet(filtered_SCGU_table, SCGU_key);
				if(new_UMI == NULL-1){
					umi = emptyUMI;
					break;
				}else if(new_UMI != NULL){
					umi = new_UMI;
				}else break;
			}


			if(cct_context->do_cell_level_junction_detection && umi[0]!='-' && !is_homopolymer_or_N_this_UR){
				int l_read_name, n_cigar_op = 0;
				memcpy(&n_cigar_op, rbinptr+16,2);
				l_read_name=((unsigned char*)rbinptr)[12];
				unsigned int * cigar_bin_ptr = ( unsigned int * )(rbinptr +36 + l_read_name);

				int ref_id = 0;
				unsigned int chro_pos = 0;

				memcpy(&ref_id, rbinptr + 4, 4);
				memcpy(&chro_pos, rbinptr + 8, 4);

				unsigned int linear_pos = cct_context -> chromosome_table . padding + chro_pos + 1; // offset: zero-based.
				if(ref_id>0) linear_pos += cct_context -> chromosome_table . read_offsets[ ref_id - 1 ]; // chr1: padding + offset; chr2: chr1_end + padding + offset, ... 

				int x1;
				for(x1=0; x1 < n_cigar_op; x1++){
					unsigned int cigar_oneop = ((unsigned int *)( rbinptr + 36 + l_read_name ))[x1];
					int cigar_op =cigar_oneop& 0xf; // cigar_op: MIDNSHP=X
					int op_len = (cigar_oneop>>4)& 0x0fffffff;

					int add_curs=0;
					if( cigar_op == 0 || cigar_op == 2 || cigar_op == 3 ) add_curs = op_len;
					if( cigar_op == 3 /* 'N' in cigar*/){
						unsigned int left_edge_last_exonbase = linear_pos -1;
						unsigned int right_edge_first_exonbase = linear_pos + add_curs;
						srUInt_64 junckey = (left_edge_last_exonbase *(1LLU<<32)) | right_edge_first_exonbase;
						ArrayList * mylist = HashTableGet(thread_context -> junction_to_cell_umi_table[sampleid], NULL+junckey);
						if(!mylist){
							mylist = ArrayListCreate(3); // [cell_and_umi_plus_one,...]
							HashTablePut(thread_context -> junction_to_cell_umi_table[sampleid], NULL+junckey, mylist);
						}
						srUInt_64 cell_and_umi =(cellid * 1LLU<<32)| convert_umi_to_2bit_int(umi, cct_context->UMI_length);
						ArrayListPush(mylist, NULL+cell_and_umi);
					}
					linear_pos += add_curs;
				}
			}

			REP_fwrite(&sampleid, 1, 4, fp);
			REP_fwrite(&cellid, 1, 4, fp);
			REP_fwrite(&gene_no, 1, 8, fp);
			if(gene_no & (1LLU<<63))REP_fwrite( glist_ptr, 1, 8*genes, fp );
			REP_fwrite(umi,1, cct_context -> UMI_length, fp);
			char * new_cellbc = NULL;
			char visiumHD_cellbc [12];// 01234_01234

			int binlen;
			memcpy(&binlen, binptr+16+8*genes+cct_context -> UMI_length,4 );
			if(cellid>=0){
				if(cct_context->visium_hd_barcodes){
					snprintf(visiumHD_cellbc,12,"%05d_%05d", (cellid&0xffff0000)>>16  , cellid&0xffff );
					//fprintf(stderr,"HAD_2D_BC %s\n", visiumHD_cellbc);
					new_cellbc = visiumHD_cellbc;
				}else new_cellbc = ArrayListGet(cct_context -> cell_barcodes_array, cellid);
			}
			cellCounts_do_one_batch_write_extend_rbin(cct_context, rbinptr, binlen, fp, new_cellbc, umi[0]=='-'?NULL:umi, gene_no, (srInt_64*)glist_ptr);
		}
		REP_fclose(fp);
		HashTableDestroy(supp_reads_SCGU);
		HashTableDestroy(filtered_SCGU_table);
		for(x1 =0; x1< cct_context -> sample_sheet_table -> numOfElements; x1++)ArrayListDestroy(cell_gene_umi_list[x1]);
		free(cell_gene_umi_list);
	}
	free(batch_content);
	free(bin_ptrs);
	return NULL + removed_UMIs;
}

#ifdef DO_STARSOLO_THING
#define DO_CREATE_BAI_FOR_BAM 0
#else
#define DO_CREATE_BAI_FOR_BAM 1
#endif

void cellCounts_save_BAM_result(cellcounts_global_t * cct_context, struct scRNA_merge_batches_worker_current * finished_job){
	if(!finished_job -> task)return;
	if(cct_context -> is_BAM_and_FQ_out_generated){
		int sample_id = finished_job -> task -> sample_id;
		void ** fps = HashTableGet(cct_context -> sample_BAM_writers, NULL+sample_id);
		simple_bam_writer * wtr = (*fps);
		int inbin_pos = 0, outblocki = 0;
		int nextoffset = -1;
		if(CELLCOUNTS_BAMBLOCK_COMP_NUMBER>1)nextoffset = finished_job -> task -> inbin_batch_start_offsets[1];
		int block_number_this = finished_job -> task -> block_number - finished_job -> task -> inbin_number +1;
		if(DO_CREATE_BAI_FOR_BAM) while(inbin_pos < finished_job -> task -> inbin_len){
			int binlen = 0;
			if(outblocki < finished_job -> task -> inbin_number -1 &&  inbin_pos == nextoffset ){
				outblocki ++;
				if(outblocki < finished_job -> task -> inbin_number -1 && CELLCOUNTS_BAMBLOCK_COMP_NUMBER>1)nextoffset = finished_job -> task -> inbin_batch_start_offsets[outblocki +1];
				block_number_this = finished_job -> task -> block_number - (finished_job -> task -> inbin_number -1 - outblocki);
			}
			binlen=*(int*)(finished_job -> task -> inbin+inbin_pos);
			simple_bam_writer_update_index(wtr, finished_job -> task -> inbin+inbin_pos, binlen, block_number_this, inbin_pos);
			inbin_pos += 4+binlen;
		}

		for(outblocki=0; outblocki < finished_job -> task -> inbin_number ; outblocki ++){
			int block_number_this = finished_job -> task -> block_number - (finished_job -> task -> inbin_number -1 - outblocki);
			int inblock_size = -1;
			if(outblocki  == finished_job -> task -> inbin_number -1) inblock_size = finished_job -> task -> inbin_len - finished_job -> task -> inbin_batch_start_offsets[outblocki];
			else if(CELLCOUNTS_BAMBLOCK_COMP_NUMBER>1)inblock_size = finished_job -> task -> inbin_batch_start_offsets[outblocki +1] -  finished_job -> task -> inbin_batch_start_offsets[outblocki];
			simple_bam_write_compressed_block(wtr, finished_job -> outbin + CELLRANGER_MERGER_WORKER_BINSIZE*outblocki , finished_job -> outbin_len[outblocki ], inblock_size, finished_job -> crc32[outblocki], block_number_this);
		}
	}
	finished_job -> task = NULL;
}

void cellCounts_merged_write_sparse_unique_genes(void * ky, void * va, HashTable * tab){
	HashTable * unique_geneno1B_tab = tab -> appendix1;
	HashTable * used_cellnoP1_tab = tab -> appendix2;

	int cellbcP1 = ky-NULL;
	if(used_cellnoP1_tab && !HashTableGet(used_cellnoP1_tab, NULL+cellbcP1))return;
	HashTable * g2u = va;
	ArrayList * g2ul = HashTableKeys(g2u);
	int x1;
	for(x1=0; x1<g2ul->numOfElements; x1++){
		void *geneno1B_ptr = ArrayListGet(g2ul,x1);
		if(!HashTableGet(unique_geneno1B_tab, ArrayListGet(g2ul,x1))) HashTablePut(unique_geneno1B_tab, geneno1B_ptr, NULL+1);
		tab -> counter1 += HashTableGet(g2u, geneno1B_ptr)?1:0;
	}
	ArrayListDestroy(g2ul);
}



int cellCounts_merged_write_sparse_matrix(cellcounts_global_t * cct_context, HashTable * cellP1_to_geneP1_to_umis_tab, ArrayList * used_cell_barcodes, int sample_index, char * tabtype, unsigned char ** feature_name_array){
	int x1,x2;

	char ofname[MAX_FILE_NAME_LENGTH + 100];
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.%s.summary",cct_context->output_prefix, sample_index+1,tabtype);
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.%s.BCtab",cct_context->output_prefix, sample_index+1,tabtype);
	FILE * ofp_bcs = fopen( ofname , "w" );
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.%s.GENEtab",cct_context->output_prefix, sample_index+1,tabtype);
	FILE * ofp_genes = fopen( ofname , "w" );
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.%s.spmtx",cct_context->output_prefix, sample_index+1,tabtype);
	FILE * ofp_mtx = fopen( ofname , "w" );
	fprintf(ofp_mtx,"%%%%MatrixMarket matrix coordinate integer general\n");

	HashTable * used_cellnoP1_tab = ArrayListToLookupTable_Int(used_cell_barcodes);
	HashTable * unique_NZ_geneno1B_table = HashTableCreate(50000);
	cellP1_to_geneP1_to_umis_tab -> counter1 = 0;
	cellP1_to_geneP1_to_umis_tab -> appendix1 = unique_NZ_geneno1B_table;
	cellP1_to_geneP1_to_umis_tab -> appendix2 = used_cellnoP1_tab;
	HashTableIteration(cellP1_to_geneP1_to_umis_tab, cellCounts_merged_write_sparse_unique_genes);
	srInt_64 total_lines = cellP1_to_geneP1_to_umis_tab -> counter1;
	ArrayList * unique_NZ_genenosP1_list = HashTableKeys(unique_NZ_geneno1B_table);
	HashTableDestroy(unique_NZ_geneno1B_table);
	HashTableDestroy(used_cellnoP1_tab);
	ArrayListSort(unique_NZ_genenosP1_list, NULL);

	#ifdef __MINGW32__
	fprintf(ofp_mtx, "%" PRId64 " %" PRId64 " %" PRId64 "\n", unique_NZ_genenosP1_list -> numOfElements , used_cell_barcodes -> numOfElements,  total_lines );
	#else
	fprintf(ofp_mtx, "%lld %lld %lld\n", unique_NZ_genenosP1_list -> numOfElements , used_cell_barcodes -> numOfElements,  total_lines );
	#endif

	for(x2=0; x2 < unique_NZ_genenosP1_list -> numOfElements; x2++){
		int gene_index_0B = ArrayListGet(unique_NZ_genenosP1_list, x2) - NULL-1;
		char* gene_name = (char*) feature_name_array[gene_index_0B]; // (char*)cct_context -> gene_name_array [gene_index_0B];
		fprintf(ofp_genes,"%s\n", gene_name);
	}

	for(x1 = 0; x1 < used_cell_barcodes -> numOfElements; x1++){
		srInt_64 cellno = ArrayListGet(used_cell_barcodes, x1)-NULL;
		if(cct_context -> visium_hd_barcodes){
			fprintf(ofp_bcs,"%05d_%05d\n", (cellno&0xffff0000)>>16, cellno&0xffff);
		}else{
			char * cellbc_seq = ArrayListGet(cct_context -> cell_barcodes_array, cellno);
			fprintf(ofp_bcs,"%s\n", cellbc_seq);
		}
	}

	for(x1 = 0; x1 < used_cell_barcodes -> numOfElements; x1++){
		srInt_64 cellno = ArrayListGet(used_cell_barcodes, x1)-NULL;
		HashTable * geneno1B_to_UMIs = HashTableGet(cellP1_to_geneP1_to_umis_tab, NULL+1+cellno);

		for(x2=0; x2 < unique_NZ_genenosP1_list -> numOfElements; x2++){
			int geneno1B = ArrayListGet(unique_NZ_genenosP1_list, x2)-NULL;
			int this_umis = HashTableGet(geneno1B_to_UMIs, NULL+geneno1B) -NULL;
			if(this_umis>0)fprintf(ofp_mtx,"%d %d %d\n", x2+1, x1+1, this_umis);
		}
	}
	ArrayListDestroy(unique_NZ_genenosP1_list);
	fclose(ofp_bcs);
	fclose(ofp_genes);
	fclose(ofp_mtx);

	return 0;
}

void cellCounts_merged_45K_to_90K_sum_SUM_Level2(void * GeneNo1B, void * vUMIs, HashTable * m2){
	HashTable * summed_gene_to_umis = m2 -> appendix1;
	HashTablePut(summed_gene_to_umis, GeneNo1B, vUMIs + (HashTableGet(summed_gene_to_umis, GeneNo1B)-NULL));
}

void cellCounts_merged_45K_to_90K_sum_SUM(void * keyCellNoP1, void * Vgno_umi_tab, HashTable * me){
	HashTable * summed_gene_to_umis  = me -> appendix1;
	HashTable * bcid_look_tab = me -> appendix2;
	HashTable * geneno1B_to_UMIs_tab = Vgno_umi_tab;
	if(!HashTableGet(bcid_look_tab, keyCellNoP1))return;
	geneno1B_to_UMIs_tab -> appendix1 = summed_gene_to_umis;
	HashTableIteration(geneno1B_to_UMIs_tab ,cellCounts_merged_45K_to_90K_sum_SUM_Level2 );
}

void cellCounts_merged_45K_to_90K_sum_WRT(void * kyGeneID, void * valUMIs, HashTable * me){
	cellcounts_global_t * cct_context = me -> appendix1;
	
	FILE * ofp = me -> appendix2;

	unsigned char * gene_name = cct_context -> gene_name_array[ kyGeneID - NULL-1 ];
	fprintf(ofp, "%s\t%u\n", gene_name, (unsigned int) (valUMIs-NULL));
}

void cellCounts_merged_45K_to_90K_sum(cellcounts_global_t * cct_context, HashTable * cellP1_geneP1_UMIs_tab, ArrayList * bcid_P0_arr, int sample_no, ArrayList * loaded_features, HashTable * sorted_index_p1_to_i_p1_tab){
	HashTable * summed_gene_to_umis = HashTableCreate( 3+cellP1_geneP1_UMIs_tab->numOfElements/6 );
	HashTable * bcid_look_tab = ArrayListToLookupTable_Int(bcid_P0_arr);
	cellP1_geneP1_UMIs_tab -> appendix1 = summed_gene_to_umis;
	cellP1_geneP1_UMIs_tab -> appendix2 = bcid_look_tab;
	cellP1_geneP1_UMIs_tab -> appendix3 = cct_context;
	HashTableIteration( cellP1_geneP1_UMIs_tab, cellCounts_merged_45K_to_90K_sum_SUM );

	char ofname[MAX_FILE_NAME_LENGTH + 100];
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.AmbSum",cct_context->output_prefix, sample_no+1);
	FILE * write_fp = fopen(ofname,"w");
	fprintf(write_fp,"GeneID\tUMIs\n");
	summed_gene_to_umis -> appendix1 = cct_context;
	summed_gene_to_umis -> appendix2 = write_fp;
	void * vp2[2];
	vp2[0]=loaded_features;
	vp2[1]=sorted_index_p1_to_i_p1_tab;
	summed_gene_to_umis -> appendix3 = vp2;
	summed_gene_to_umis -> counter1 = sample_no;
	HashTableIteration( summed_gene_to_umis, cellCounts_merged_45K_to_90K_sum_WRT );
	HashTableDestroy(bcid_look_tab);
	HashTableDestroy(summed_gene_to_umis);
	fclose(write_fp);
}

void cellCounts_merged_write_nozero_geneids_WRT(void *k, void *v, HashTable* me){
	FILE * fp = me->appendix1;
	cellcounts_global_t * cct_context = me->appendix2;
	unsigned char* gene_symbol = cct_context -> gene_name_array [k-NULL-1];
	fprintf(fp, "%s\n", gene_symbol);
}

void cellCounts_merged_write_nozero_geneids(cellcounts_global_t * cct_context, HashTable * no0genes, int samplenno, ArrayList * loaded_features, HashTable * sorted_order_p1_to_i_p1_tab){
	char ofname[MAX_FILE_NAME_LENGTH + 100];
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 100,"%s.scRNA.%03d.no0Genes", cct_context->output_prefix, samplenno+1);
	FILE * fp = fopen( ofname , "w" );
	no0genes -> appendix1 =fp;
	void * tv2[2];
	no0genes -> appendix2 = cct_context;
	tv2[0]=loaded_features;
	tv2[1]=sorted_order_p1_to_i_p1_tab;
	no0genes -> appendix3 =tv2;
	HashTableIteration(no0genes, cellCounts_merged_write_nozero_geneids_WRT);
	fclose(fp);
}

void cellCounts_merged_to_tables_write_build_UMIcount_in(void * ky, void * val, HashTable * tab){
	tab -> counter1 += (val-NULL);
}

void cellCounts_merged_to_tables_write_build_UMIcounts(void * ky, void * val, HashTable * tab){
	HashTable * cellbcP1_to_umis_tab = tab -> appendix1;
	int cell_no = ky-NULL-1;
	HashTable * geneP1_to_counts_tab = val;

	geneP1_to_counts_tab -> counter1 = 0;
	HashTableIteration(geneP1_to_counts_tab, cellCounts_merged_to_tables_write_build_UMIcount_in);
	HashTablePut(cellbcP1_to_umis_tab, NULL+1+cell_no, NULL+geneP1_to_counts_tab -> counter1);
}

#define MIN_UMIS_FOR_CANDIDATE_RESCUE 500
#define SCRNA_AMBIENT_RESCURE_MEDIAN_FRACTION 0.01
#define MAX_CANDIDATE_CELLS 20000
void cellCounts_merged_ambient_rescure(cellcounts_global_t * cct_context, HashTable * cellP1_to_geneP1_to_umis_tab, HashTable * cellnoP1_to_umis_tab, ArrayList * this_sample_45k_90k_barcode_no_P0, ArrayList * this_sample_ambient_rescure_candi, ArrayList * highconf_cellbc_list){
	ArrayList * sorted_bcno_p1 = HashTableSortedIndexes( cellnoP1_to_umis_tab, 1);
	HashTable * highconf_cellbc_list_tab = ArrayListToLookupTable_Int(highconf_cellbc_list);
	srInt_64 x1, high_conf_cells = 0;
	for(x1=0; x1 < sorted_bcno_p1 -> numOfElements; x1++){
		void * this_bc_pnt = ArrayListGet(sorted_bcno_p1 ,  x1);
		if(HashTableGet(highconf_cellbc_list_tab, this_bc_pnt)) high_conf_cells = x1+1;
		else break; // assuming that all high-umi barcodes are high-confident, this makes x1 being the # of total high-confidence barcodes.	
	}
	if(high_conf_cells >0){
		srInt_64 median_umis = HashTableGet(cellnoP1_to_umis_tab, ArrayListGet(sorted_bcno_p1 ,  (high_conf_cells-1)/2))-NULL;
		srInt_64 median_umis_001_cut = (srInt_64)(median_umis *1. *SCRNA_AMBIENT_RESCURE_MEDIAN_FRACTION +0.50000001);
		for(x1=0; x1 < sorted_bcno_p1 -> numOfElements; x1++){
			void * this_bc_pnt_p1 = ArrayListGet(sorted_bcno_p1 ,  x1);
			if(HashTableGet(highconf_cellbc_list_tab, this_bc_pnt_p1)){
				continue; // it is in high-conf list
			}
			srInt_64 this_bc_umis = HashTableGet(cellnoP1_to_umis_tab, this_bc_pnt_p1) - NULL;
			if(this_bc_umis < median_umis_001_cut) break;
			if(this_bc_umis < MIN_UMIS_FOR_CANDIDATE_RESCUE) break;
			if(x1 >= 45000) break;

			// The MAX_CANDIDATE_CELLS was added on 31 JAN 2023 when I realised that Cell Ranger only tests 20,000 candidates.
			if(this_sample_ambient_rescure_candi -> numOfElements < MAX_CANDIDATE_CELLS)
				ArrayListPush(this_sample_ambient_rescure_candi, this_bc_pnt_p1-1);
		}
	}
	for(x1=45000; x1 < sorted_bcno_p1 -> numOfElements; x1++){
		if(x1 >= 90000) break;
		ArrayListPush(this_sample_45k_90k_barcode_no_P0, ArrayListGet(sorted_bcno_p1 ,  x1)-1 );
	}
	ArrayListDestroy(sorted_bcno_p1);
	HashTableDestroy(highconf_cellbc_list_tab);
}



#define SCRNA_BOOTSTRAP_HIGH_INDEX 30
#define SCRNA_BOOTSTRAP_SAMPLING_TIMES 100 

// static is safe because only one thread;
static srUInt_64 bootstrap_seed = 1234567890123456789ULL, bootstrap_seed2 = 987654321000ULL; 

srUInt_64 bootstrap_rand_U64(srUInt_64 addvar) {
	srUInt_64 x = bootstrap_seed;
	srUInt_64 y = bootstrap_seed2;

	x += addvar;
	x ^= x >> 12; // a
	x ^= x << 25; // b
	x ^= x >> 27; // c

	bootstrap_seed = bootstrap_seed2 ^ (x*0x2545F4914F6CDD1DULL);
	bootstrap_seed2 = x;
	return y+bootstrap_seed2;
}

int cellCounts_merged_bootstrap_a_sample(cellcounts_global_t * cct_context, HashTable * cellP1_to_geneP1_to_umis_tab, HashTable * cellnoP1_to_umis_tab, ArrayList * highconf_cellbc_list){
	ArrayList * list_cellBCs_sorted_by_UMIs = HashTableSortedIndexes( cellnoP1_to_umis_tab, 1); // "1": large first
	srInt_64 x2, x1;
	float cellCounts_umi_cutoff = cct_context -> umi_cutoff;

	srInt_64 this_total = 0;
	bootstrap_seed = bootstrap_seed ^ list_cellBCs_sorted_by_UMIs -> numOfElements;
//#warning "XXXXXXXXXXXXXXXXX RANDOMIZE SEED XXXXXXXXXXXXXXXX"
//	double miltval = miltime()*1000;
//	bootstrap_seed ^= (srInt_64)miltval;

	int last_umi_no= -1;
	if(cellCounts_umi_cutoff >= 0.0){
		for(x1 = 0; x1 < list_cellBCs_sorted_by_UMIs -> numOfElements ; x1++){
			void * cellbc_p1_ptr = ArrayListGet(list_cellBCs_sorted_by_UMIs,x1);
			srInt_64 this_umis = HashTableGet(cellnoP1_to_umis_tab, cellbc_p1_ptr )-NULL;
			if(this_umis >= cellCounts_umi_cutoff-0.1){
				ArrayListPush(highconf_cellbc_list, ArrayListGet( list_cellBCs_sorted_by_UMIs, x1 ) - 1 );
				last_umi_no = this_umis;
			}else break;	// #UMI-sorted so no need to scan more
		}
	}else{
		for(x1 = 0; x1 < SCRNA_BOOTSTRAP_SAMPLING_TIMES; x1++){
			ArrayList * sorted_resampled_UMIs = ArrayListCreate( list_cellBCs_sorted_by_UMIs->numOfElements );
			for(x2 = 0; x2 < list_cellBCs_sorted_by_UMIs -> numOfElements ; x2++){
				srUInt_64 seed_rand = bootstrap_rand_U64( this_total ^ (cellP1_to_geneP1_to_umis_tab -> numOfElements<<12))  % (srUInt_64)list_cellBCs_sorted_by_UMIs -> numOfElements;
				void * cellbc_p1_ptr = ArrayListGet(list_cellBCs_sorted_by_UMIs, seed_rand);
				srInt_64 this_umis = HashTableGet( cellnoP1_to_umis_tab, cellbc_p1_ptr )-NULL;
				ArrayListPush(sorted_resampled_UMIs,NULL+this_umis);
			}
			ArrayListSort( sorted_resampled_UMIs, NULL );
			srInt_64 UMIs_30th_div10 = ArrayListGet(sorted_resampled_UMIs, sorted_resampled_UMIs -> numOfElements - SCRNA_BOOTSTRAP_HIGH_INDEX) -NULL;
			UMIs_30th_div10 = (srInt_64)(UMIs_30th_div10*1./10 + 0.500000001);

			for(x2 =0; x2< sorted_resampled_UMIs -> numOfElements; x2++){
				srInt_64 lli = sorted_resampled_UMIs -> numOfElements -1 -x2;
				srInt_64 this_umis = ArrayListGet(sorted_resampled_UMIs, lli)-NULL;
				if(this_umis >= UMIs_30th_div10) this_total ++;
				else break;
			}
			ArrayListDestroy(sorted_resampled_UMIs);
		}
		double total_f = this_total*1. / SCRNA_BOOTSTRAP_SAMPLING_TIMES;
		this_total = (int)(total_f + 0.500000001);

		void * last_ptr =NULL;
		for(x1 = 0; x1 < min(list_cellBCs_sorted_by_UMIs -> numOfElements, this_total) ; x1++){
			last_ptr = ArrayListGet( list_cellBCs_sorted_by_UMIs, x1 );
			ArrayListPush(highconf_cellbc_list, last_ptr - 1 );
		}
		last_umi_no = HashTableGet(cellnoP1_to_umis_tab ,last_ptr)-NULL;
	}
	ArrayListDestroy(list_cellBCs_sorted_by_UMIs);
	return last_umi_no;
}

void cellCounts_finalise_per_junction_cell_table(cellcounts_global_t * cct_context, ArrayList ** useable_cell_ids);
void cellCounts_merged_to_tables_write(cellcounts_global_t * cct_context, HashTable ** cellP1_to_geneP1_to_umis, ArrayList * loaded_features, srInt_64 nexons, ArrayList ** output_needed_cellids){
	char ofname[MAX_FILE_NAME_LENGTH + 20];
	SUBreadSprintf(ofname,MAX_FILE_NAME_LENGTH + 20,"%s.scRNA.SampleTable",cct_context->output_prefix);
	FILE * sample_tab_fp = fopen( ofname , "w" );
	int x1;

	fprintf(sample_tab_fp,"SampleName\tUMICutoff\tTotalReads\tMappedReads\tAssignedReads\tIndex\n");
	for(x1 = 0; x1 < cct_context -> sample_sheet_table -> numOfElements ; x1++){
		ArrayList * high_confid_barcode_index_list = ArrayListCreate(20000);
		output_needed_cellids[x1 +1] = high_confid_barcode_index_list; // output: 1-based sample; x1: 0-based.
		ArrayList * this_sample_ambient_rescure_candi = ArrayListCreate(10000);
		ArrayList * this_sample_45k_90k_barcode_no_P0 = ArrayListCreate(90000 - 45000 + 100);

		HashTable * cellbcP1_to_umis_tab = HashTableCreate(cellP1_to_geneP1_to_umis[x1] -> numOfElements);
		cellP1_to_geneP1_to_umis[x1] -> appendix1 = cellbcP1_to_umis_tab;
		HashTableIteration(cellP1_to_geneP1_to_umis[x1], cellCounts_merged_to_tables_write_build_UMIcounts);

		int applied_umi_cut = -1;
		if(cellbcP1_to_umis_tab -> numOfElements >0) applied_umi_cut = cellCounts_merged_bootstrap_a_sample(cct_context, cellP1_to_geneP1_to_umis[x1], cellbcP1_to_umis_tab, high_confid_barcode_index_list);
		cct_context -> applied_umi_cut[x1] = applied_umi_cut;
		cellCounts_merged_ambient_rescure(cct_context, cellP1_to_geneP1_to_umis[x1], cellbcP1_to_umis_tab, this_sample_45k_90k_barcode_no_P0, this_sample_ambient_rescure_candi, high_confid_barcode_index_list);

		int umi_cutoff = cct_context -> applied_umi_cut[x1];
		char * this_sample_name = ArrayListGet( cct_context -> sample_id_to_name, x1);
#ifdef __MINGW32__
		fprintf(sample_tab_fp,"%s\t%d\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t%d\n", this_sample_name, umi_cutoff,  cct_context -> reads_per_sample[x1],  cct_context -> mapped_reads_per_sample[x1],  cct_context -> assigned_reads_per_sample[x1] ,x1+1);
#else
		fprintf(sample_tab_fp,"%s\t%d\t%lld\t%lld\t%lld\t%d\n", this_sample_name, umi_cutoff, cct_context -> reads_per_sample[x1],  cct_context -> mapped_reads_per_sample[x1],  cct_context -> assigned_reads_per_sample[x1], x1+1);
#endif
		srInt_64 xk1;
		HashTable * sorted_order_p1_to_i_p1_tab = HashTableCreate(nexons/4);
		for(xk1 = 0; xk1 < nexons ; xk1++){
			fc_feature_info_t * feature1 = ArrayListGet(loaded_features, xk1);
			HashTablePut(sorted_order_p1_to_i_p1_tab, NULL+ feature1->sorted_order+1 , NULL+xk1+1 );
		}

		if(cct_context -> report_excluded_barcodes){
			ArrayList * all_barcode_index_list = HashTableKeys(cellbcP1_to_umis_tab);
			for(xk1=0; xk1<all_barcode_index_list->numOfElements; xk1++) all_barcode_index_list->elementList[xk1]--;
			cellCounts_merged_write_sparse_matrix(cct_context, cellP1_to_geneP1_to_umis[x1], all_barcode_index_list, x1, "RawOut", cct_context -> gene_name_array);
			ArrayListDestroy(all_barcode_index_list);
		}
		cellCounts_merged_write_sparse_matrix(cct_context, cellP1_to_geneP1_to_umis[x1], high_confid_barcode_index_list, x1, "HighConf",  cct_context -> gene_name_array);
		cellCounts_merged_write_sparse_matrix(cct_context, cellP1_to_geneP1_to_umis[x1], this_sample_ambient_rescure_candi, x1, "RescCand",  cct_context -> gene_name_array);
		cellCounts_merged_45K_to_90K_sum( cct_context, cellP1_to_geneP1_to_umis[x1], this_sample_45k_90k_barcode_no_P0, x1 , loaded_features, sorted_order_p1_to_i_p1_tab);
		HashTable * no0genes = HashTableCreate(50000);
		cellP1_to_geneP1_to_umis[x1] -> appendix1 = no0genes;
		cellP1_to_geneP1_to_umis[x1] -> appendix2 = NULL;
		HashTableIteration(cellP1_to_geneP1_to_umis[x1], cellCounts_merged_write_sparse_unique_genes);
		cellCounts_merged_write_nozero_geneids(cct_context, no0genes, x1, loaded_features, sorted_order_p1_to_i_p1_tab);

		HashTableDestroy(no0genes);
		ArrayListExtend(high_confid_barcode_index_list, this_sample_ambient_rescure_candi);
		ArrayListDestroy(this_sample_ambient_rescure_candi);
		ArrayListDestroy(this_sample_45k_90k_barcode_no_P0);
		HashTableDestroy(cellbcP1_to_umis_tab);
		HashTableDestroy(sorted_order_p1_to_i_p1_tab);
	}

	fclose(sample_tab_fp);

}

int bin_file_exists(char *filename) {
	struct stat   buffer;   
	return (stat (filename, &buffer) == 0);
}

void * delete_file_thread(void * arg){
	void ** ptrs =arg;
	srInt_64 *current_sorting_key = NULL;
	cellcounts_global_t * cct_context = NULL;
	current_sorting_key = ptrs[0];
	cct_context = ptrs[1];
	
	while(1){
		int all_closed = 1;
		int x1;
		for(x1=0; x1<CELLBC_BATCH_NUMBER +2; x1++){
			if(NULL== current_sorting_key || current_sorting_key[x1] == 0x7fffffffffffffffLLU){
				char tmp_fname[MAX_FILE_NAME_LENGTH+50];
				SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+50, "%s/temp-cellcounts-%06d-%03d.tmpbin", cct_context -> temp_file_dir, getpid(), x1);
				if(bin_file_exists(tmp_fname)) unlink(tmp_fname);
			}else all_closed = 0;
		}
		if(all_closed ) break;
		sleep(2);
	}
	if(current_sorting_key)free(current_sorting_key);
	return NULL;
}

int cellCounts_do_cellbc_batches(cellcounts_global_t * cct_context){
	pthread_t *threads = malloc(sizeof(pthread_t)* cct_context-> total_threads);
	int sample_i,xk2,xk1,compress_workers = max(1, cct_context-> total_threads-1);
	HashTable * cellnoP1_to_genenoP1_to_UMIs[cct_context -> sample_sheet_table -> numOfElements];
	struct scRNA_merge_batches_worker_task * task_buffers = malloc(sizeof(struct scRNA_merge_batches_worker_task) * (1+compress_workers)* cct_context->sample_sheet_table -> numOfElements);
	int current_filling_worker_per_sample [cct_context-> sample_sheet_table -> numOfElements];
	struct scRNA_merge_batches_worker_current * worker_current_jobs = calloc(sizeof(struct scRNA_merge_batches_worker_current), compress_workers);

	ArrayList * file_size_list = ArrayListCreate(CELLBC_BATCH_NUMBER +1);
	for(xk1=0; xk1<CELLBC_BATCH_NUMBER +2; xk1++){
		if(xk1<CELLBC_BATCH_NUMBER +1){
			srInt_64 batchsize = REP_filesize( cct_context -> batch_files[xk1]);
			ArrayListPush(file_size_list, NULL+( batchsize<<20 | xk1));
		}
		REP_fclose(cct_context -> batch_files[xk1]);
	}
	ArrayListSort(file_size_list, NULL);

	srInt_64 block_numbers_current [cct_context-> sample_sheet_table -> numOfElements];
	for(xk1=0; xk1< cct_context-> sample_sheet_table -> numOfElements; xk1++){
		cellnoP1_to_genenoP1_to_UMIs[xk1] = HashTableCreate(10000);
		HashTableSetDeallocationFunctions(cellnoP1_to_genenoP1_to_UMIs[xk1], NULL,(void (*) (void*))HashTableDestroy);
		current_filling_worker_per_sample[xk1] = 0;
		task_buffers[xk1].inbin_len = 0;
		task_buffers[xk1].inbin_number = 0;
		task_buffers[xk1].inbin_batch_start_offsets[0]=0;
		block_numbers_current[xk1] = 0;
	}

	for(xk1=0; xk1<compress_workers+1; xk1++)for(xk2 = 0; xk2 < cct_context-> sample_sheet_table -> numOfElements; xk2++) task_buffers[xk1*cct_context->sample_sheet_table -> numOfElements + xk2].sample_id = xk2+1;

	cellcounts_align_thread_t * thread_contexts = calloc(sizeof(cellcounts_align_thread_t) , cct_context->total_threads);
	cct_context -> all_thread_contexts = thread_contexts;

	for(xk1=0; xk1< cct_context-> total_threads; xk1++){
		void ** vpp = malloc(sizeof(void*)*3);
		vpp[0] = cct_context;
		vpp[1] = NULL + xk1;
		vpp[2] = file_size_list;

		cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + xk1;
		memset(thread_context -> junction_to_cell_umi_table, 0, sizeof(void*)*cct_context-> sample_sheet_table -> numOfElements);

		if(cct_context -> do_cell_level_junction_detection) for(sample_i = 1; sample_i <=cct_context-> sample_sheet_table -> numOfElements; sample_i++){
			thread_context -> junction_to_cell_umi_table[sample_i] = HashTableCreate(20000);
			HashTableSetDeallocationFunctions(thread_context -> junction_to_cell_umi_table[sample_i] , NULL,  (void (*)(void *value))ArrayListDestroy);
		}

		pthread_create(threads + xk1, NULL, cellCounts_do_one_batch, vpp);
	}

	srInt_64 removed_umis = 0;
	for(xk1=0; xk1<cct_context-> total_threads; xk1++){
		void * pret = NULL;
		pthread_join(threads[xk1], &pret);
		removed_umis += (pret - NULL);

		cellcounts_align_thread_t * thread_context = cct_context -> all_thread_contexts + xk1;

		if(cct_context -> do_cell_level_junction_detection) for(sample_i = 1; sample_i <=cct_context-> sample_sheet_table -> numOfElements; sample_i++){
			thread_context -> junction_to_cell_umi_table[sample_i]->appendix1 = cct_context -> junction_to_cell_umi_table[sample_i];
			HashTableIteration( thread_context -> junction_to_cell_umi_table[sample_i],  cellCounts_join_thread_junc_cell_umi_table);
			HashTableDestroy(thread_context -> junction_to_cell_umi_table[sample_i]);
		}
	}
	free(thread_contexts);
	cct_context -> all_thread_contexts = NULL;

	print_in_box(80,0,0,"");
	if(cct_context -> input_mode== GENE_INPUT_BCL){
		srInt_64 all_extracted_reads = 0;
		for(xk1 = 0; xk1 < cct_context-> sample_sheet_table -> numOfElements+1; xk1++) 
			all_extracted_reads += cct_context-> reads_per_sample[xk1];

		for(xk1 = 0; xk1 < cct_context-> sample_sheet_table -> numOfElements; xk1++) {
			srInt_64 extracted_reads = cct_context-> reads_per_sample[xk1];
			char * sample_name = ArrayListGet(cct_context-> sample_id_to_name, xk1);
#ifdef __MINGW32__
			print_in_box(81,0,0,"  % 13" PRId64 " (%4.1f%%%%) reads were assigned to %s.\n", extracted_reads, extracted_reads*100./all_extracted_reads, sample_name);
#else
			print_in_box(81,0,0,"  %'13lld (%4.1f%%%%) reads were assigned to %s.\n", extracted_reads, extracted_reads*100./all_extracted_reads, sample_name);
#endif
		}
		print_in_box(80,0,0,"");

#ifdef __MINGW32__
		if(cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements] < 0.005*all_extracted_reads) print_in_box(81,0,0,"  % 13" PRId64 "(%4.0f%%%%) reads were assigned to samples in total.", all_extracted_reads - cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements], 100.-cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements]*100./all_extracted_reads);
		else print_in_box(81,0,0,"  % 13" PRId64 " (%4.1f%%%%) reads were assigned to samples in total.", all_extracted_reads - cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements], 100.-cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements]*100./all_extracted_reads);
#else
		if(cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements] < 0.005*all_extracted_reads) print_in_box(81,0,0,"  %'13lld (%4.0f%%%%) reads were assigned to samples in total.", all_extracted_reads - cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements], 100.-cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements]*100./all_extracted_reads);
		else print_in_box(81,0,0,"  %'13lld (%4.1f%%%%) reads were assigned to samples in total.", all_extracted_reads - cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements], 100.-cct_context-> reads_per_sample[cct_context-> sample_sheet_table -> numOfElements]*100./all_extracted_reads);
#endif

	}else{
		for(xk1 = 0; xk1 < cct_context-> sample_sheet_table -> numOfElements; xk1++) {
			srInt_64 extracted_reads = cct_context-> reads_per_sample[xk1];
			char * sample_name = ArrayListGet(cct_context-> sample_id_to_name, xk1);
#ifdef __MINGW32__
			print_in_box(80,0,0,"  % 13" PRId64 " reads were processed for %s.\n", extracted_reads, sample_name);
#else
			print_in_box(80,0,0,"  %'13lld reads were processed for %s.\n", extracted_reads, sample_name);
#endif
		}

	}
	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"Generate UMI count tables...");
	ArrayListDestroy(file_size_list);

	worker_master_mutex_t worker_mut;
	worker_master_mutex_init(&worker_mut, max(1, cct_context-> total_threads-1));

	for(xk1=0; xk1<max(1, cct_context-> total_threads-1); xk1++){
		void ** vpp = malloc(sizeof(void*)*4);
		vpp[0] = cct_context;
		vpp[1] = &worker_mut;
		vpp[2] = NULL + xk1;
		vpp[3] = worker_current_jobs + xk1;
		pthread_create(threads + xk1, NULL, cellCounts_merge_batches_worker, vpp);
	}

	REPFILE * input_fps[CELLBC_BATCH_NUMBER+2];
	char * last_rbin_buffer[CELLBC_BATCH_NUMBER+1];
	srInt_64 * current_sorting_key = malloc(sizeof(srInt_64)*(CELLBC_BATCH_NUMBER+2));
	
	for(xk1=0; xk1< CELLBC_BATCH_NUMBER+2; xk1++){
		char tmp_fname[MAX_FILE_NAME_LENGTH+80];
		SUBreadSprintf(tmp_fname, MAX_FILE_NAME_LENGTH+80, "%s/temp-cellcounts-%06d-%03d.tmpbin", cct_context -> temp_file_dir, getpid(), xk1);
		input_fps[xk1] = REP_fopen(tmp_fname,"rb");
		REP_setvbuf(input_fps[xk1], cct_context -> cellbin_v_buffers[xk1], _IOFBF, SCRNA_SMALLER_VBUFF_SIZE);
		if(xk1 == CELLBC_BATCH_NUMBER+1)break;

		srInt_64 section1_items=0;
		for(sample_i = 0; sample_i < cct_context -> sample_sheet_table -> numOfElements; sample_i++){
			size_t frret = REP_fread(&section1_items,1, 8, input_fps[xk1]);
			for(xk2 = 0; xk2 < section1_items; xk2++){
				srInt_64 cellbcP0_geneno0B=0, umis=0;
				frret += REP_fread(&cellbcP0_geneno0B,1,8,input_fps[xk1]);
				frret += REP_fread(&umis,1,8,input_fps[xk1]);

				int cellbc_no = cellbcP0_geneno0B>>32;
				int gene_no0B = (int)(cellbcP0_geneno0B&0xffffffffu);
				HashTable *gene_tab = HashTableGet(cellnoP1_to_genenoP1_to_UMIs[sample_i], NULL+cellbc_no+1);
				if(gene_tab==NULL){
					gene_tab = HashTableCreate(300);
					HashTablePut(cellnoP1_to_genenoP1_to_UMIs[sample_i], NULL+cellbc_no+1, gene_tab); 
				}
				HashTablePut(gene_tab, NULL+gene_no0B+1 , NULL+umis);
			}
		}
		last_rbin_buffer[xk1] = malloc( cct_context -> barcode_batched_max_genes *8 + cct_context -> barcode_batched_max_Rbin_len + 4 + MAX_UMI_LEN + 16 + 10000);
		int rlen = REP_fread(last_rbin_buffer[xk1], 1, 16, input_fps[xk1]);
		if(rlen >0){
			int binlen = 0;
			srInt_64 genes = 0;
			memcpy(&genes, last_rbin_buffer[xk1]+8, 8);
			if(genes & (1LLU<<63))genes = genes & 0x7fffffff;
			else genes= 0;
			
			size_t frret = REP_fread(last_rbin_buffer[xk1]+16, 1, 8*genes+ cct_context -> UMI_length + 4, input_fps[xk1]);
			memcpy(&binlen, last_rbin_buffer[xk1] +16 +8*genes+ cct_context -> UMI_length  , 4);
			frret += REP_fread(last_rbin_buffer[xk1] + 16+ 8*genes+ cct_context -> UMI_length + 4, 1, binlen, input_fps[xk1]);

			srInt_64 sorting_key = *(int*)(last_rbin_buffer[xk1] + 16 +8*genes+cct_context -> UMI_length +4);
			sorting_key = sorting_key << 32;
			sorting_key |= *(int*)(last_rbin_buffer[xk1] + 16+ 8*genes+cct_context -> UMI_length +8);
			current_sorting_key[xk1] = sorting_key;
		}else{
			REP_fclose(input_fps[xk1]);
			free(last_rbin_buffer[xk1]);
			current_sorting_key[xk1] = 0x7fffffffffffffffLLU;
		}
	}

	current_sorting_key[CELLBC_BATCH_NUMBER+1]=0;
	void * del_ptrs[2];
	del_ptrs[0] = current_sorting_key;
	del_ptrs[1] = cct_context ;
	pthread_create(&cct_context->thread_delete_files, NULL, delete_file_thread, del_ptrs);
	int current_worker = 0;
	while(1){
		int selected_fp_no = 0;
		srInt_64 selected_fp_key = current_sorting_key[0];
		for(xk1=1; xk1<CELLBC_BATCH_NUMBER+1; xk1++){
			if(current_sorting_key[xk1] < selected_fp_key){
				selected_fp_key = current_sorting_key[xk1] ;
				selected_fp_no = xk1;
			}
		}
		if(selected_fp_key == 0x7fffffffffffffffLLU) break;

		int sample_id = 0, binlen = 0;
		srInt_64 genes = 0;
		memcpy(&sample_id, last_rbin_buffer[selected_fp_no], 4);
		memcpy(&genes, last_rbin_buffer[selected_fp_no]+8, 8);
		if(genes & (1LLU<<63)) genes = genes & 0x7fffffff;
		else genes = 0;
		memcpy(&binlen,last_rbin_buffer[selected_fp_no]+16+8*genes+cct_context -> UMI_length,4);
		struct scRNA_merge_batches_worker_task * tofill = task_buffers+(current_filling_worker_per_sample[sample_id-1] * cct_context->sample_sheet_table -> numOfElements +sample_id-1);
		memcpy(tofill->inbin + tofill-> inbin_len, last_rbin_buffer[selected_fp_no]+16+8*genes+cct_context -> UMI_length, binlen + 4);
		tofill -> inbin_len += (binlen + 4);
		if(tofill -> inbin_number ==0) tofill -> inbin_number =1;
		//SUBREADprintf("ADDING BLOCKKK = %d  WKR = %d  IT THINK IT'S %d ; GENES=%d\n", tofill -> inbin_len, current_worker, tofill -> sample_id, genes);
		if(tofill-> inbin_len > CELLCOUNTS_BAMBLOCK_SIZE * CELLCOUNTS_BAMBLOCK_COMP_NUMBER){
			master_wait_for_job_done(&worker_mut, current_worker);
			struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
			cellCounts_save_BAM_result(cct_context, my_finished_job);
			my_finished_job -> task = tofill;
			my_finished_job -> task -> block_number = (block_numbers_current[sample_id-1]++);
			master_notify_worker(&worker_mut, current_worker);

			current_filling_worker_per_sample[sample_id-1] ++;
			if(current_filling_worker_per_sample[sample_id-1] == compress_workers +1) current_filling_worker_per_sample[sample_id-1] = 0;
			tofill = task_buffers+(current_filling_worker_per_sample[sample_id-1] * cct_context->sample_sheet_table -> numOfElements +sample_id-1);
			tofill -> inbin_len = 0;
			tofill -> inbin_number = 0;
			tofill -> inbin_batch_start_offsets[0]=0;
			
			current_worker ++;
			if(current_worker == compress_workers) current_worker=0;
		}else if( tofill -> inbin_len - tofill -> inbin_batch_start_offsets[ tofill -> inbin_number - 1 ] > CELLCOUNTS_BAMBLOCK_SIZE ){
			tofill -> inbin_batch_start_offsets[tofill -> inbin_number]= tofill -> inbin_len;
			tofill -> inbin_number++;
			block_numbers_current[sample_id-1]++;
		}

		int rlen = REP_fread(last_rbin_buffer[selected_fp_no], 1, 16, input_fps[selected_fp_no]);
		if(rlen >0){
			int binlen = 0;
			srInt_64 genes = 0;
			memcpy(&genes, last_rbin_buffer[selected_fp_no]+8, 8);
			if(genes & (1LLU<<63))genes = genes & 0x7fffffff;
			else genes= 0;
			size_t frret = REP_fread(last_rbin_buffer[selected_fp_no]+16, 1, 8*genes+ cct_context -> UMI_length + 4, input_fps[selected_fp_no]);
			memcpy(&binlen, last_rbin_buffer[selected_fp_no] +16 +8*genes+ cct_context -> UMI_length  , 4);
			frret += REP_fread(last_rbin_buffer[selected_fp_no] + 16+ 8*genes+ cct_context -> UMI_length + 4, 1, binlen, input_fps[selected_fp_no]);
			srInt_64 sorting_key = *(int*)(last_rbin_buffer[selected_fp_no] + 16+8*genes +cct_context -> UMI_length +4);
			sorting_key = sorting_key << 32;
			sorting_key |= *(int*)(last_rbin_buffer[selected_fp_no] + 16 +8*genes+cct_context -> UMI_length +8);
			current_sorting_key[selected_fp_no] = sorting_key;
		} else {
			REP_fclose(input_fps[selected_fp_no]);
			free(last_rbin_buffer[selected_fp_no]);
			current_sorting_key[selected_fp_no] = 0x7fffffffffffffffLLU;
		}
	}

	for(xk1=0; xk1<cct_context -> sample_sheet_table -> numOfElements; xk1++){
		struct scRNA_merge_batches_worker_task * tofill = task_buffers+(current_filling_worker_per_sample[xk1] * cct_context-> sample_sheet_table -> numOfElements +xk1);
		if(tofill->inbin_len<1) continue;

		master_wait_for_job_done(&worker_mut, current_worker);
		struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
		cellCounts_save_BAM_result(cct_context, my_finished_job);
		my_finished_job -> task = tofill;
		my_finished_job -> task -> block_number = (block_numbers_current[xk1]++);
		master_notify_worker(&worker_mut, current_worker);
		current_worker ++;
		if(current_worker == compress_workers) current_worker=0;
	}
	for(xk1=0; xk1<compress_workers; xk1++){
		struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
		if(my_finished_job -> task)master_wait_for_job_done(&worker_mut, current_worker);
		cellCounts_save_BAM_result(cct_context, my_finished_job);

		current_worker ++;
		if(current_worker == compress_workers) current_worker=0;
	}

	for(xk1 = 0; xk1 < 1+compress_workers; xk1++) for(xk2 = 0; xk2 < cct_context->sample_sheet_table -> numOfElements;xk2++) {
		task_buffers[ xk1 * cct_context->sample_sheet_table -> numOfElements  + xk2 ].inbin_len = 0;
		task_buffers[ xk1 * cct_context->sample_sheet_table -> numOfElements  + xk2 ].inbin_number = 0;
	}
	current_worker = 0;
	REPFILE * notmapped_fp = input_fps[CELLBC_BATCH_NUMBER+1];
	while(1){
		int sample_id = 0, binlen = 0;
		int rlen = REP_fread(&sample_id, 1, 4, notmapped_fp);
		if(rlen < 4) break;
		struct scRNA_merge_batches_worker_task * tofill = task_buffers+(current_filling_worker_per_sample[sample_id -1] * cct_context->sample_sheet_table -> numOfElements +sample_id-1);
		size_t frret = REP_fread(&binlen, 1, 4, notmapped_fp);
		char old_bin[binlen+4];
		memcpy(old_bin, &binlen,4);
		frret += REP_fread(old_bin+4, 1, binlen, notmapped_fp);
		int new_binlen = cellCounts_make_barcode_bam_bin(cct_context, old_bin, tofill -> inbin + tofill -> inbin_len, binlen, NULL, NULL, -1, NULL);
		tofill -> inbin_len += 4+ new_binlen; // block size: Total length of the alignment record, excluding this field. Then, the alignment record.

		if(tofill -> inbin_number ==0) tofill -> inbin_number =1;
		if(tofill-> inbin_len > CELLCOUNTS_BAMBLOCK_SIZE * CELLCOUNTS_BAMBLOCK_COMP_NUMBER){
			struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
			if(my_finished_job -> task)master_wait_for_job_done(&worker_mut, current_worker);
			cellCounts_save_BAM_result(cct_context, my_finished_job);
			my_finished_job -> task = tofill;
			my_finished_job -> task -> block_number = (block_numbers_current[sample_id-1]++);
			master_notify_worker(&worker_mut, current_worker);

			current_filling_worker_per_sample[sample_id-1] ++;
			if(current_filling_worker_per_sample[sample_id-1] == compress_workers +1) current_filling_worker_per_sample[sample_id-1] = 0;
			tofill = task_buffers+(current_filling_worker_per_sample[sample_id-1] * cct_context-> sample_sheet_table -> numOfElements +sample_id-1);
			tofill -> inbin_len = 0;
			tofill -> inbin_number = 0;
			tofill -> inbin_batch_start_offsets[0]=0;
			
			current_worker ++;
			if(current_worker == compress_workers) current_worker=0;
		}else if( tofill -> inbin_len - tofill -> inbin_batch_start_offsets[ tofill -> inbin_number - 1 ] > CELLCOUNTS_BAMBLOCK_SIZE ){
			tofill -> inbin_batch_start_offsets[tofill -> inbin_number]= tofill -> inbin_len;
			block_numbers_current[sample_id-1]++;
			tofill -> inbin_number++;
		}
	}
	REP_fclose(notmapped_fp);
	current_sorting_key[CELLBC_BATCH_NUMBER+1] = 0x7fffffffffffffffLLU;

	for(xk1=0; xk1<cct_context -> sample_sheet_table -> numOfElements; xk1++){
		struct scRNA_merge_batches_worker_task * tofill = task_buffers+(current_filling_worker_per_sample[xk1] * cct_context-> sample_sheet_table -> numOfElements +xk1);
		if(tofill->inbin_len<1) continue;

		master_wait_for_job_done(&worker_mut, current_worker);
		struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
		cellCounts_save_BAM_result(cct_context, my_finished_job);
		my_finished_job -> task = tofill;
		my_finished_job -> task -> block_number = (block_numbers_current[xk1]++);
		master_notify_worker(&worker_mut, current_worker);
		current_worker ++;
		if(current_worker == compress_workers) current_worker=0;
	}

	for(xk1=0; xk1<compress_workers; xk1++){
		master_wait_for_job_done(&worker_mut, current_worker);
		struct scRNA_merge_batches_worker_current * my_finished_job = worker_current_jobs+current_worker;
		cellCounts_save_BAM_result(cct_context, my_finished_job);

		current_worker ++;
		if(current_worker == compress_workers) current_worker=0;
	}

	terminate_workers(&worker_mut);
	free(task_buffers);
	free(worker_current_jobs);

	for(xk1=0; xk1< compress_workers; xk1++) pthread_join(threads[xk1],NULL);

	worker_master_mutex_destroy(&worker_mut);

	ArrayList * to_keep_cell_ids_per_samples[cct_context -> sample_sheet_table -> numOfElements+1];
	memset(to_keep_cell_ids_per_samples, 0, sizeof(void*)*(cct_context -> sample_sheet_table -> numOfElements+1));
	if(!strstr(cct_context->index_prefix,".InternalTest.Rsr.step5@rand")){
		cellCounts_merged_to_tables_write(cct_context , cellnoP1_to_genenoP1_to_UMIs , cct_context -> all_features_array, cct_context -> all_features_array->numOfElements, to_keep_cell_ids_per_samples);
		if(cct_context -> do_cell_level_junction_detection) cellCounts_finalise_per_junction_cell_table(cct_context, to_keep_cell_ids_per_samples);
	}

	for(xk1=0; xk1< cct_context -> sample_sheet_table -> numOfElements; xk1++){ // cellnoP1_to_genenoP1_to_UMIs is 0-based sample ids; to_keep_cell_ids_per_samples is 1-based sample ids.
		if(to_keep_cell_ids_per_samples[xk1]) ArrayListDestroy(to_keep_cell_ids_per_samples[xk1 +1]);
		HashTableDestroy(cellnoP1_to_genenoP1_to_UMIs[xk1]);
	}
	return 0;
}

void cellCounts_finalise_per_junction_add_to_table(void * ky, void * val, HashTable * tab){
	srUInt_64 junc_LR_linear = tab -> counter1;
	HashTable * cellid_p1_to_junc_supp_tab = tab -> appendix1;
	HashTable * juncLR_to_juncid_p1_tab = tab -> appendix2;

	srInt_64 junc_id = HashTableGet(juncLR_to_juncid_p1_tab, NULL+junc_LR_linear)-NULL-1;
	int cellid = ky-NULL-1;
	srInt_64 count = val - NULL;
	HashTable * juncno_p1_to_supp = HashTableGet(cellid_p1_to_junc_supp_tab, NULL+1+cellid);
	if(!juncno_p1_to_supp){
		juncno_p1_to_supp = HashTableCreate(5000);
		HashTablePut(cellid_p1_to_junc_supp_tab, NULL+1+cellid, juncno_p1_to_supp);
	}
	if(count) HashTablePut(juncno_p1_to_supp, NULL+ junc_id+1, NULL+count);
}

void cellCounts_finalise_per_junc_sumcounts(void * ky, void * val, HashTable * tab){
	srUInt_64 junc_LR_linear = ky - NULL;

	ArrayList * cellid_umiseq_list = val;
	void ** sumparams = tab -> appendix1;
	HashTable * cellid_p1_to_junc_supp_tab = sumparams[0]; // cell_id +NULL+1 => JuncID_+1 => count
	HashTable * wanted_cellid_p1_tab = sumparams[1];
	HashTable * juncjunc_to_junc_id_p1_tab = sumparams[2];
	HashTable * cellid_umi_p1_seq_tab = HashTableCreate(1000);

	srInt_64 xx;
	for(xx=0; xx< cellid_umiseq_list -> numOfElements; xx++){
		srUInt_64 cellid_umiseq = ArrayListGet(cellid_umiseq_list, xx)-NULL;
		int cellid = (cellid_umiseq>>32)& 0x7fffffff;
		void * needed_cell = HashTableGet(wanted_cellid_p1_tab , NULL+cellid+1);
		if(!needed_cell) continue;
		HashTablePut(cellid_umi_p1_seq_tab, NULL+cellid_umiseq +1, NULL+1);
	}

	HashTable * cellid_p1_to_count_table = HashTableCreate(1000);
	ArrayList * uniq_cellid_umi_seq_p1_list = HashTableKeys(cellid_umi_p1_seq_tab);
	for(xx=0; xx < uniq_cellid_umi_seq_p1_list -> numOfElements; xx++){
		srUInt_64 cellid_umi_p1 = ArrayListGet(uniq_cellid_umi_seq_p1_list, xx)-NULL;
		int cellid =((cellid_umi_p1-1) >>32) & 0x7fffffff;
		srInt_64 count0 = HashTableGet(cellid_p1_to_count_table, NULL+cellid+1) - NULL;
		HashTablePut(cellid_p1_to_count_table, NULL+cellid+1, NULL+1+count0);
	}

	cellid_p1_to_count_table -> appendix1 = cellid_p1_to_junc_supp_tab;
	cellid_p1_to_count_table -> appendix2 = juncjunc_to_junc_id_p1_tab;
	cellid_p1_to_count_table -> counter1 = junc_LR_linear;

	ArrayListDestroy(uniq_cellid_umi_seq_p1_list);
	HashTableIteration(cellid_p1_to_count_table, cellCounts_finalise_per_junction_add_to_table);
	HashTableDestroy(cellid_p1_to_count_table);
	HashTableDestroy(cellid_umi_p1_seq_tab);
}

void cellCounts_finalise_per_junction_cell_table(cellcounts_global_t * cct_context, ArrayList** highconf_and_candidate_cell_ids_sps){
	int thread_no, sample_i;
	char *temp_dir = cct_context -> temp_file_dir;

	for(sample_i = 1; sample_i <=cct_context-> sample_sheet_table -> numOfElements ; sample_i ++){
		ArrayList * highconf_and_candidate_cell_ids = highconf_and_candidate_cell_ids_sps[sample_i];
		HashTable * needed_cellid_p1_tab = ArrayListToLookupTable_Int(highconf_and_candidate_cell_ids);
		HashTable * cct_junctab = cct_context -> junction_to_cell_umi_table[sample_i];
		ArrayList * cct_junc_LRlist = HashTableKeys(cct_junctab);
		HashTable * juncjunc_to_juncid_p1_table = HashTableCreate(20000);


		ArrayList * juncname_list = ArrayListCreate(juncjunc_to_juncid_p1_table -> numOfElements+1);
		ArrayListSetDeallocationFunction(juncname_list, free);
		HashTable * cellid_p1_output_table = HashTableCreate(15000);
		srInt_64 xx;
		for(xx=0; xx < cct_junc_LRlist -> numOfElements; xx++){
			srUInt_64 LRval = ArrayListGet( cct_junc_LRlist, xx )-NULL;
			HashTablePut(juncjunc_to_juncid_p1_table ,NULL+LRval,NULL+xx+1);
			char junction_name_tmp[MAX_CHROMOSOME_NAME_LEN + 12*2 +1], *chro_name=NULL, *chro_name_R=NULL;
			int chro_pos_L=0, chro_pos_R=0;

			locate_gene_position( (LRval>>32)&0xffffffffu , &cct_context -> chromosome_table, &chro_name, &chro_pos_L);
			locate_gene_position( (LRval)&0xffffffffu , &cct_context -> chromosome_table, &chro_name_R, &chro_pos_R);
			assert(chro_name==chro_name_R);

			SUBreadSprintf(junction_name_tmp,sizeof(junction_name_tmp),"%s:%d^%d", chro_name, chro_pos_L, chro_pos_R);
			ArrayListPush(juncname_list, strdup(junction_name_tmp));
		}

		void * sum_params[3];
		sum_params[0]= cellid_p1_output_table;
		sum_params[1]= needed_cellid_p1_tab;
		sum_params[2]= juncjunc_to_juncid_p1_table;
		cct_junctab -> appendix1 = sum_params;

		HashTableIteration(cct_junctab, cellCounts_finalise_per_junc_sumcounts );

		ArrayList * junc_highcand_cellids = ArrayList_Int_Hash_Intersect(highconf_and_candidate_cell_ids, cellid_p1_output_table);
		cellCounts_merged_write_sparse_matrix(cct_context, cellid_p1_output_table, junc_highcand_cellids, 
			sample_i -1, "cellJuncs", (unsigned char**)juncname_list -> elementList); // this function adds 1 to the sample no.
//		FILE * mtx_junc_fp???? // to write cellid_p1_output_table :  cell_id +NULL+1 => [L|R , count, L|R, count, ...]

		ArrayListDestroy(junc_highcand_cellids);
		ArrayListDestroy(juncname_list);
		ArrayListDestroy(cct_junc_LRlist);
		HashTableDestroy(needed_cellid_p1_tab);
		HashTableDestroy(juncjunc_to_juncid_p1_table);
	}
}


int cellCounts_run_counting(cellcounts_global_t * cct_context){
	int ret= 0;
	ret = ret || cellCounts_do_cellbc_batches(cct_context);
	ret = ret || cellCounts_write_gene_list(cct_context);
	if( cct_context -> do_cell_level_junction_detection ) ret = ret || cellCounts_write_junction_sumtable(cct_context);
	return ret;
}

void cellCounts_finalise_error_run(cellcounts_global_t * cct_context){
	void * argsv[2];
	argsv[0]=NULL;
	argsv[1]=cct_context;
	delete_file_thread(argsv);
}

#ifdef MAKE_STANDALONE
	#define cellCounts_main main
#endif
int cellCounts_main(int argc, char** argv){
	setlocale(LC_ALL,"");
	cellcounts_global_t * cct_context = calloc(sizeof(cellcounts_global_t),1);
	cct_context -> program_start_time = miltime();

	int ret = 0;
	ret = ret || cellCounts_args_context(cct_context, argc, argv);
	ret = ret || cellCounts_load_context(cct_context);
	ret = ret || cellCounts_run_mapping(cct_context);
	ret = ret || cellCounts_run_counting(cct_context);
	ret = ret || cellCounts_destroy_context(cct_context);
	if(cct_context -> has_error) cellCounts_finalise_error_run(cct_context);

	free(cct_context);
	if(ret) SUBREADprintf("cellCounts terminates with errors.\n");
	return ret;
}



#define JC_STATUS_NA 0
#define JC_STATUS_KNOWN 1
#define JC_STATUS_NOVEL 2
#define JC_STATUS_NOVEL_FUSION 3
#define MAX_OVERLAP_EDGE_NUMBER 1000
#define JC_OUT_GENE_COLUMNS_LENGTH (MAX_OVERLAP_EDGE_NUMBER * (4+CHROMOSOME_NAME_LENGTH)) 


void cellCounts_find_nearest_gene_dist(cellcounts_global_t * cct_context, int side_small, int side_large, char * dist_to_nearest_splice_side_str_SP1,  char * dist_to_nearest_splice_side_str_SP2,
				int junc_near_LLedge_no, int junc_near_LRedge_no, int junc_near_RLedge_no, int junc_near_RRedge_no,
				IVT_Interval ** junc_nearest_LLedges, IVT_Interval ** junc_nearest_LRedges, IVT_Interval ** junc_nearest_RLedges, IVT_Interval ** junc_nearest_RRedges, char * my_chro);

int cellCounts_determine_jcount_gene_transcript_report(cellcounts_global_t * cct_context, int side_small, int side_large, int junc_olay_genebody_left_result_no,int junc_olay_genebody_right_result_no, IVT_Interval ** junc_genebody_olayleft, IVT_Interval ** junc_genebody_olayright, char * gene_ids_str_SP1, char * gene_ids_str_SP2, char * transcript_ids_str_SP1, char * transcript_ids_str_SP2, char * dist_to_nearest_splice_side_str_SP1,  char * dist_to_nearest_splice_side_str_SP2, int strand_learnt_from_FASTA, int junc_near_LLedge_no, int junc_near_LRedge_no, int junc_near_RLedge_no, int junc_near_RRedge_no, IVT_Interval ** junc_nearest_LLedges, IVT_Interval ** junc_nearest_LRedges, IVT_Interval ** junc_nearest_RLedges, IVT_Interval ** junc_nearest_RRedges, char * my_chro){
	int xk1,xk2, txp_id, retv=0;

	gene_ids_str_SP1[0] = transcript_ids_str_SP1[0] = dist_to_nearest_splice_side_str_SP1[0] =
	gene_ids_str_SP2[0] = transcript_ids_str_SP2[0] = dist_to_nearest_splice_side_str_SP2[0] = 'N';

	gene_ids_str_SP1[1] = transcript_ids_str_SP1[1] = dist_to_nearest_splice_side_str_SP1[1] =
	gene_ids_str_SP2[1] = transcript_ids_str_SP2[1] = dist_to_nearest_splice_side_str_SP2[1] = 'A';

	gene_ids_str_SP1[2] = transcript_ids_str_SP1[2] = dist_to_nearest_splice_side_str_SP1[2] =
	gene_ids_str_SP2[2] = transcript_ids_str_SP2[2] = dist_to_nearest_splice_side_str_SP2[2] = '\0';

	HashTable * match1_txn_table = StringTableCreate(100);
	HashTable * edge1P1_table = StringTableCreate(100);
	HashTable * edge2P1_table = StringTableCreate(100);

	HashTable * edge1exonNo_table = StringTableCreate(100);
	HashTable * edge2exonNo_table = StringTableCreate(100);

	ArrayList * common_txn_list = ArrayListCreate(10);
	int small_site_exactly = 0;
	int large_site_exactly = 0;
	for(xk1=0; xk1<junc_olay_genebody_left_result_no; xk1++){
		cct_junction_genebody_t * jg_ptr = junc_genebody_olayleft[xk1]->attr;
		for(txp_id = 0; txp_id < jg_ptr -> transcript_list -> numOfElements; txp_id ++){
			srInt_64 edge1_dist = 0xffffffffu;
			cct_junction_transcript_t * txp_ptr = ArrayListGet(jg_ptr -> transcript_list, txp_id);
			ArrayList * exon_list = txp_ptr -> exons_in_transcript;
			int exactly_matched_exon_no = -1;
			for(xk2=0; xk2< exon_list->numOfElements; xk2++){
				cct_junction_exon_in_transcript_t * jte_ptr = ArrayListGet(exon_list,xk2);
				if(strand_learnt_from_FASTA >=0 && strand_learnt_from_FASTA !=  jte_ptr -> is_negative) continue;
				int exon_start_known = jte_ptr -> chro_start;
				int exon_end_known = jte_ptr -> chro_stop;
				int dist_to_any_side = min(abs(side_small - exon_end_known),abs(side_small - exon_start_known));
				if(dist_to_any_side < edge1_dist){
					edge1_dist = dist_to_any_side;
					if(0==edge1_dist) exactly_matched_exon_no = xk2;
				}
			}
			if(edge1_dist!=0xffffffffu){
				HashTablePut(edge1P1_table, txp_ptr -> transcript_id, NULL+1+edge1_dist);
				if(0==edge1_dist){
					small_site_exactly++;
					HashTablePut(edge1exonNo_table, txp_ptr -> transcript_id, NULL+1+exactly_matched_exon_no);
				}
			}
		}
	}
	
	for(xk1=0; xk1<junc_olay_genebody_right_result_no; xk1++){
		cct_junction_genebody_t * jg_ptr = junc_genebody_olayright[xk1]->attr;
		for(txp_id = 0; txp_id < jg_ptr -> transcript_list -> numOfElements; txp_id ++){
			srInt_64 edge2_dist = 0xffffffffu;
			cct_junction_transcript_t * txp_ptr = ArrayListGet(jg_ptr -> transcript_list, txp_id);
			ArrayList * exon_list = txp_ptr -> exons_in_transcript;
			int exactly_matched_exon_no = -1;
			for(xk2=0; xk2< exon_list->numOfElements; xk2++){
				cct_junction_exon_in_transcript_t * jte_ptr = ArrayListGet(exon_list,xk2);
				if(strand_learnt_from_FASTA >=0 && strand_learnt_from_FASTA !=  jte_ptr -> is_negative) continue;
				int exon_end_known = jte_ptr -> chro_stop; 
				int exon_start_known = jte_ptr -> chro_start; 
				int dist_to_any_side = min(abs(side_large - exon_end_known),abs(side_large - exon_start_known));
				if(dist_to_any_side < edge2_dist){
					edge2_dist = dist_to_any_side;
					if(0==edge2_dist) exactly_matched_exon_no = xk2;
				}
			}
			if(edge2_dist!=0xffffffffu){
				HashTablePut(edge2P1_table, txp_ptr -> transcript_id, NULL+1+edge2_dist);
				if(0==edge2_dist){
					large_site_exactly++;
					HashTablePut(edge2exonNo_table, txp_ptr -> transcript_id, NULL+1+exactly_matched_exon_no);
				}
			}
		}
	}

	if(small_site_exactly >0 && large_site_exactly >0){
		ArrayList * edge1_txnids = HashTableKeys(edge1P1_table);
		for(xk1=0; xk1<edge1_txnids->numOfElements; xk1++){
			char * edge1_txnid = ArrayListGet(edge1_txnids, xk1);
			int edge1_ptr = HashTableGet(edge1P1_table, edge1_txnid)-NULL;
			if(edge1_ptr!=1)continue;
			int edge2_ptr = HashTableGet(edge2P1_table, edge1_txnid)-NULL;
			if(edge2_ptr==1){
				int exonno_1 = HashTableGet(edge1exonNo_table, edge1_txnid)-NULL-1;
				int exonno_2 = HashTableGet(edge2exonNo_table, edge1_txnid)-NULL-1;
				if(abs(exonno_1 - exonno_2)==1 || strstr(edge1_txnid, CCT_NIL_TXN_PLACEHOLDER))ArrayListPush(common_txn_list , edge1_txnid);
			}
		}
	}

	if(common_txn_list -> numOfElements > 0){
		HashTable * out_gene_tab = StringTableCreate(10);
		HashTable * out_txn_tab = StringTableCreate(10);
		for(xk1=0; xk1<common_txn_list->numOfElements; xk1++){
			char *common_txn_id = ArrayListGet(common_txn_list,xk1);
			cct_junction_transcript_t *txnptr = HashTableGet(cct_context ->junction_transcript_table, common_txn_id);
			HashTablePut(out_gene_tab, txnptr -> gene_name, NULL+1);
			HashTablePut(out_txn_tab, txnptr -> transcript_id, NULL+1);
		}
		ArrayList * out_gene_list = HashTableKeys(out_gene_tab);
		ArrayList * out_txn_list = HashTableKeys(out_txn_tab);

		ArrayListSort(out_gene_list, ArrayListStringComparison);
		ArrayListSort(out_txn_list, ArrayListStringComparison);

		ArrayListStringJoin(out_gene_list, gene_ids_str_SP1,JC_OUT_GENE_COLUMNS_LENGTH);
		ArrayListStringJoin(out_txn_list, transcript_ids_str_SP1,JC_OUT_GENE_COLUMNS_LENGTH);

		// when they share the same genes and txns, the lists must be identical for two sides. 
		ArrayListStringJoin(out_gene_list, gene_ids_str_SP2,JC_OUT_GENE_COLUMNS_LENGTH);
		//ArrayListStringJoin(out_txn_list, transcript_ids_str_SP2,JC_OUT_GENE_COLUMNS_LENGTH);

		ArrayListDestroy(out_gene_list);
		ArrayListDestroy(out_txn_list);
		HashTableDestroy(out_txn_tab);
		HashTableDestroy(out_gene_tab);
		retv = JC_STATUS_KNOWN; 
	}
	if(retv==0){
		int side_i;
		for(side_i=0; side_i<2; side_i++){
			HashTable * this_side_reported_genes = StringTableCreate(10);
			int this_side_has_exactly = side_i?large_site_exactly:small_site_exactly;
			HashTable * this_side_olay_tab = side_i?edge2P1_table:edge1P1_table;
			ArrayList * edge_txnids = HashTableKeys(this_side_olay_tab);

			for(xk1 = 0; xk1 < edge_txnids -> numOfElements; xk1++){
				char * this_side_txn_id = ArrayListGet(edge_txnids,xk1);
				srInt_64 edge_dist = HashTableGet(this_side_olay_tab, this_side_txn_id)-NULL;
				if( this_side_has_exactly==0 || edge_dist == 1 ){
					cct_junction_transcript_t * txn_ptr = HashTableGet(cct_context -> junction_transcript_table,  this_side_txn_id);
					HashTablePut(this_side_reported_genes, txn_ptr -> gene_name, NULL+1);
				}
			}

			ArrayList * this_side_reported_list = HashTableKeys(this_side_reported_genes);
			ArrayListSort(this_side_reported_list, ArrayListStringComparison);
			ArrayListStringJoin(this_side_reported_list, side_i?gene_ids_str_SP2:gene_ids_str_SP1,JC_OUT_GENE_COLUMNS_LENGTH);
			ArrayListDestroy(this_side_reported_list);
			HashTableDestroy(this_side_reported_genes);
		}

		retv = JC_STATUS_NOVEL;
	}

	ArrayListDestroy(common_txn_list);
	HashTableDestroy(match1_txn_table);
	HashTableDestroy(edge1P1_table);
	HashTableDestroy(edge2P1_table);
	HashTableDestroy(edge1exonNo_table);
	HashTableDestroy(edge2exonNo_table);

	cellCounts_find_nearest_gene_dist(cct_context, side_small, side_large,
				dist_to_nearest_splice_side_str_SP1, dist_to_nearest_splice_side_str_SP2,
				junc_near_LLedge_no, junc_near_LRedge_no, junc_near_RLedge_no, junc_near_RRedge_no,
				junc_nearest_LLedges,  junc_nearest_LRedges,  junc_nearest_RLedges,  junc_nearest_RRedges, my_chro);

	if(strstr(transcript_ids_str_SP1,CCT_NIL_TXN_PLACEHOLDER)) strcpy(transcript_ids_str_SP1,"NA");
	if(cct_context -> ignore_transcript_junction_assignment) retv = JC_STATUS_NA; // i.e., if annotation input is from SAF, not GTF.
	
	return retv;
}

void cellCounts_find_nearest_gene_dist(cellcounts_global_t * cct_context, int side_small, int side_large, char * dist_to_nearest_splice_side_str_SP1,  char * dist_to_nearest_splice_side_str_SP2,
				int junc_near_LLedge_no, int junc_near_LRedge_no, int junc_near_RLedge_no, int junc_near_RRedge_no,
				IVT_Interval ** junc_nearest_LLedges, IVT_Interval ** junc_nearest_LRedges, IVT_Interval ** junc_nearest_RLedges, IVT_Interval ** junc_nearest_RRedges, char * exon_SE_chro){
	int side_i, xk1;
	for(side_i=0; side_i<2; side_i++){
		int Lscan_edge_no = side_i?junc_near_RLedge_no:junc_near_LLedge_no;
		int Rscan_edge_no = side_i?junc_near_RRedge_no:junc_near_LRedge_no;
		srInt_64 this_side = side_i?side_large:side_small;
		int L_scan_dist = -1, L_exon_SE_coord = -1;
		int R_scan_dist = -1, R_exon_SE_coord = -1;
		IVT_Interval ** L_scan_res = side_i?junc_nearest_RLedges:junc_nearest_LLedges;
		IVT_Interval ** R_scan_res = side_i?junc_nearest_RRedges:junc_nearest_LRedges;
		int show_genes_L = 0;
		int show_genes_R = 0;
		if(Lscan_edge_no >0){
			L_scan_dist = abs( this_side - L_scan_res[0]->start );
			L_exon_SE_coord = L_scan_res[0]->start;
		}
		if(Rscan_edge_no >0){
			R_scan_dist = abs( this_side - R_scan_res[0]->start );
			R_exon_SE_coord = R_scan_res[0]->start;
		}

		int final_dist = -1;
		if(Lscan_edge_no >0 && Rscan_edge_no<1){
			show_genes_L = 1;
			final_dist = L_scan_dist;
		}else if(Lscan_edge_no <1 && Rscan_edge_no>0){
			show_genes_R = 1;
			final_dist = R_scan_dist;
		}else if(Lscan_edge_no >0 && Rscan_edge_no>0){
			if(L_scan_dist < R_scan_dist) show_genes_L = 1;
			else if(L_scan_dist > R_scan_dist) show_genes_R = 1;
			else{
				show_genes_L = 1;
				show_genes_R = 1;
			}
			final_dist = min(R_scan_dist, L_scan_dist);
		}
		if(final_dist<1) show_genes_R = 0; //L and R search had the same results.
		if(show_genes_L || show_genes_R){
			char * outchrs = side_i?dist_to_nearest_splice_side_str_SP2:dist_to_nearest_splice_side_str_SP1;
			int lri, outchrs_ptr=0;
			for(lri=0; lri<2; lri++){
				HashTable * gene_name_tab = StringTableCreate(10);
				if(lri==0 && !show_genes_L) continue;
				if(lri==1 && !show_genes_R) continue;
				int this_scan_dir_items;
				IVT_Interval ** this_scan_dir_item_ptr;
				if(side_i){
					this_scan_dir_items = lri?junc_near_RRedge_no:junc_near_RLedge_no;
					this_scan_dir_item_ptr = lri?junc_nearest_RRedges:junc_nearest_RLedges;
				}else{
					this_scan_dir_items = lri?junc_near_LRedge_no:junc_near_LLedge_no;
					this_scan_dir_item_ptr = lri?junc_nearest_LRedges:junc_nearest_LLedges;
				}
				for(xk1 = 0; xk1 < this_scan_dir_items ; xk1++){
					char * gene_name = this_scan_dir_item_ptr[xk1] -> attr;
					if((void*)gene_name  != NULL+IMPOSSIBLE_MEMORY_SPACE) HashTablePut(gene_name_tab, gene_name, NULL+1);
				}

				ArrayList * gene_name_list = HashTableKeys(gene_name_tab);
				ArrayListSort(gene_name_list, ArrayListStringComparison);
				int me_exon_SE_coord = lri?R_exon_SE_coord:L_exon_SE_coord;
				outchrs_ptr += SUBreadSprintf(outchrs+outchrs_ptr, JC_OUT_GENE_COLUMNS_LENGTH - outchrs_ptr,"%s:%d,", exon_SE_chro, me_exon_SE_coord);
				outchrs_ptr += ArrayListStringJoin(gene_name_list, outchrs+outchrs_ptr, JC_OUT_GENE_COLUMNS_LENGTH - outchrs_ptr-12);
				ArrayListDestroy(gene_name_list);
				HashTableDestroy(gene_name_tab);

				char * lrstr = lri?",right":",left";
				if (final_dist < 1) lrstr = "";
				outchrs_ptr += SUBreadSprintf(outchrs+outchrs_ptr, JC_OUT_GENE_COLUMNS_LENGTH - outchrs_ptr,",%d%s", final_dist, lrstr);
				if(outchrs_ptr)outchrs[outchrs_ptr ++]=';';
			}
			if(outchrs[outchrs_ptr-1]==';') outchrs_ptr--;
			outchrs[outchrs_ptr] = 0;
		}
	}
}

void cellCounts_copy_junc_table_to_output_table(void * k, void * v, HashTable * oldtab){
	HashTable * newtab = oldtab -> appendix1;
	cellcounts_global_t * cct_context = oldtab -> appendix2;
	int sample_i = oldtab -> counter1;
	chroEvent_t * ed = v;
	
	if(ed -> event_type != chroEvent_t_TYPE_JUNCTION)return;
//#warning "========= ALL JUNCTIONS INCLUDE 0 JUNCS ARE REPORTED ======="
	int report_0_junc=0;
	if(!report_0_junc)if(!ed -> step2_supported_reads)return;

	char junckey[MAX_CHROMOSOME_NAME_LEN*2+24];
	char * chro=NULL;
	int pos_small=0, pos_large=0;
	locate_gene_position(ed -> left_edge, &cct_context -> chromosome_table, &chro, &pos_small);
	locate_gene_position(ed -> left_edge + ((void*)ed -> length-NULL) +1, &cct_context -> chromosome_table, &chro, &pos_large);
	
	// format: chro-small<TAB>pos-small<TAB>chro-large<TAB>pos-large
	SUBreadSprintf(junckey, MAX_CHROMOSOME_NAME_LEN*2+24, "%s\t%d\t%s\t%d", chro, pos_small+1, chro, pos_large+1);
	srInt_64 * srtab = HashTableGet(newtab, junckey);
	if(!srtab){
		srtab = calloc(sizeof(srInt_64), MAX_SCRNA_SAMPLE_NUMBER +1);
		HashTablePut(newtab, strdup(junckey), srtab);
	}
	srInt_64 cck = (void*)ed -> step2_supported_reads - NULL + report_0_junc;
	cck = cck << 32 |((void*)ed -> step2_non_supported_reads - NULL);
	srtab[sample_i] = cck;
}

void cellCounts_write_final_other_events(cellcounts_global_t * cct_context,  char * output_file_name, int type_mode){
//TODO: hasn't been tested!
	int sample_i = 0 ; // to be done

	ArrayList * small_large_keys = HashTableKeys(cct_context -> chroEvent_detail_table[sample_i]);
	ArrayListSort(small_large_keys, ArrayListLLUComparison);
	FILE * ofp = fopen(output_file_name,"w");
	int x1;
	for(x1=0; x1< small_large_keys-> numOfElements; x1++){
		chroEvent_t * ed = HashTableGet(cct_context -> chroEvent_detail_table[sample_i] , ArrayListGet(small_large_keys,x1));
		if(ed -> event_type != type_mode) continue;

		char * chro = NULL;
		int poss = 0;

		locate_gene_position(ed -> left_edge +1, &cct_context -> chromosome_table, &chro, &poss);
		int x2;
		for(x2=0; x2 < ed -> n_events; x2++){
			int envlen, nsupp, nonsupp;
			if(ed -> n_events >1) {
				envlen = ed -> length[x2];
				nsupp = ed -> step2_supported_reads[x2];
				nonsupp = ed -> step2_non_supported_reads[x2];
			}else{
				envlen = (void*)ed -> length - NULL;
				nsupp = (void*)ed -> step2_supported_reads - NULL;
				nonsupp = (void*)ed -> step2_non_supported_reads - NULL;
			}
			
			if(nsupp)fprintf(ofp,"%s\t%d\t%d%c\t%d\t%d\n", chro, poss, abs(envlen), envlen>0?'D':'I', nsupp, nonsupp);
		}
	}
	ArrayListDestroy(small_large_keys);
	fclose(ofp);
}

void cellCounts_write_final_junctions(cellcounts_global_t * cct_context,  char * output_file_name){
	int infile_i, disk_is_full = 0, sample_i;

	HashTable * junction_table = StringTableCreate(200000); //  str(chro, pos_small+1, chro, pos_large+1) => srInt_64[0..39] of nsup <<32 | nnonsup
	HashTableSetDeallocationFunctions(junction_table,free,free);
	for(sample_i=1; sample_i <=cct_context-> sample_sheet_table -> numOfElements ; sample_i ++){
		cct_context -> chroEvent_detail_table[sample_i] -> appendix1 = junction_table;
		cct_context -> chroEvent_detail_table[sample_i] -> appendix2 = cct_context;
		cct_context -> chroEvent_detail_table[sample_i] -> counter1 = sample_i;
		HashTableIteration(cct_context -> chroEvent_detail_table[sample_i],cellCounts_copy_junc_table_to_output_table);
	}

	gene_value_index_t * current_value_index = cct_context->value_index;

	char ** key_list;
	key_list = malloc(sizeof(char *) * junction_table -> numOfElements);

	KeyValuePair * cursor;
	int bucket, ky_i = 0;
	for(bucket=0; bucket < junction_table -> numOfBuckets; bucket++){
		cursor = junction_table -> bucketArray[bucket];
		while (cursor){
			char * ky = (char *)cursor -> key;

			key_list[ky_i ++] = ky;
			cursor = cursor -> next;
		}
	}

	merge_sort(key_list,  junction_table -> numOfElements , cellCounts_junckey_sort_compare, cellCounts_junckey_sort_exchange, cellCounts_junckey_sort_merge);

	char outfname[MAX_FILE_NAME_LENGTH];
	SUBreadSprintf(outfname, MAX_FILE_NAME_LENGTH, "%s.jcounts", output_file_name);

	int max_junction_genes = 3000;
	char * gene_names = malloc(max_junction_genes * FEATURE_NAME_LENGTH), * gene_name_tail;

	int ky_i1, ky_i2;
	FILE * ofp = fopen(outfname, "w");
	char * tmpp = NULL;

	fprintf(ofp, "Gene_SP1\tGene_SP2\tTranscript\t"
            "Status\tDonor\tAcceptor\t"
            "Chr_SP1\tLocation_SP1\tStrand_SP1\tChr_SP2\tLocation_SP2\tStrand_SP2\t"
            "NearestExonBoundary_SP1\tNearestExonBoundary_SP2" 
        );

	for(sample_i = 1; sample_i <=cct_context-> sample_sheet_table -> numOfElements ; sample_i ++)
		fprintf(ofp, "\tSupporting_Reads_%03d\tUnsupported_Reads_%03d", sample_i, sample_i);
	fprintf(ofp, "\n");

	IVT_Interval ** junc_genebody_olayleft = malloc(sizeof(void*) * MAX_OVERLAP_EDGE_NUMBER);
	IVT_Interval ** junc_genebody_olayright = malloc(sizeof(void*) * MAX_OVERLAP_EDGE_NUMBER);
	IVT_Interval ** junc_nearest_LLedges = malloc(sizeof(void*) * MAX_OVERLAP_EDGE_NUMBER);
	IVT_Interval ** junc_nearest_RLedges = malloc(sizeof(void*) * MAX_OVERLAP_EDGE_NUMBER);
	IVT_Interval ** junc_nearest_LRedges = malloc(sizeof(void*) * MAX_OVERLAP_EDGE_NUMBER);
	IVT_Interval ** junc_nearest_RRedges = malloc(sizeof(void*) * MAX_OVERLAP_EDGE_NUMBER);
	for(ky_i = 0; ky_i < junction_table -> numOfElements ; ky_i ++){
		int unique_junctions = 0;
		char * chro_small = strtok_r( key_list[ky_i] , "\t", &tmpp);
		char * pos_small_str = strtok_r( NULL, "\t", &tmpp);
		char * chro_large = strtok_r( NULL, "\t", &tmpp);
		char * pos_large_str = strtok_r( NULL, "\t", &tmpp);

		unsigned int pos_small = atoi(pos_small_str);
		unsigned int pos_large = atoi(pos_large_str);

		char * strand = "NA";
		if(1){
			unsigned int linear_small = 0, linear_large = 0;
			char donor[2], receptor[2];
			linear_small = linear_gene_position(&cct_context->chromosome_table , chro_small, pos_small-1);
			linear_large = linear_gene_position(&cct_context->chromosome_table , chro_small, pos_large-1);
			donor[0] = gvindex_get(current_value_index,linear_small+1);
			donor[1] = gvindex_get(current_value_index,linear_small+2);

			receptor[0] = gvindex_get(current_value_index,linear_large-2);
			receptor[1] = gvindex_get(current_value_index,linear_large-1);

			if(donor[0]=='G' && donor[1]=='T' && receptor[0]=='A' && receptor[1]=='G') strand = "+";
			else if(donor[0]=='C' && donor[1]=='T' && receptor[0]=='A' && receptor[1]=='C') strand = "-";
		}
		assert(0==strcmp(chro_small, chro_large));

		IVT_IntervalTreeNode * IVT_gbody_root = HashTableGet(cct_context -> junction_GenebodyTree_table, chro_small);

		int this_edge_tab_i = 0; // no strand info
		if(strand[0]=='+') this_edge_tab_i = 1;
		if(strand[0]=='-') this_edge_tab_i = 2;
		HashTable * this_edge_tab = cct_context -> junction_ExonEdgeTree_table[this_edge_tab_i];
		IVT_IntervalTreeNode * IVT_edge_root = HashTableGet(this_edge_tab, chro_small);

		char gene_ids_str_SP1[JC_OUT_GENE_COLUMNS_LENGTH], gene_ids_str_SP2[JC_OUT_GENE_COLUMNS_LENGTH],
		     transcript_ids_str_SP1[JC_OUT_GENE_COLUMNS_LENGTH], transcript_ids_str_SP2[JC_OUT_GENE_COLUMNS_LENGTH],
		     dist_to_nearest_splice_side_str_SP1[JC_OUT_GENE_COLUMNS_LENGTH], dist_to_nearest_splice_side_str_SP2[JC_OUT_GENE_COLUMNS_LENGTH];
		strcpy(gene_ids_str_SP1,"NA");
		strcpy(gene_ids_str_SP2,"NA");
		strcpy(transcript_ids_str_SP1,"NA");
		strcpy(transcript_ids_str_SP2,"NA");
		strcpy(dist_to_nearest_splice_side_str_SP1,"NA");
		strcpy(dist_to_nearest_splice_side_str_SP2,"NA");

		char * jc_retv_str = "NA";
		int jc_gene_status = -1;

		int strand_learned_from_FASTA = -1;
		if(strand[0]=='+') strand_learned_from_FASTA =0;
		if(strand[0]=='-') strand_learned_from_FASTA =1;
		if(IVT_gbody_root){
			int junc_olay_genebody_left_result_no, junc_olay_genebody_right_result_no, junc_near_LLedge_no, junc_near_LRedge_no, junc_near_RLedge_no, junc_near_RRedge_no;
			junc_olay_genebody_left_result_no = IVT_query(IVT_gbody_root, pos_small, junc_genebody_olayleft, MAX_OVERLAP_EDGE_NUMBER);
			junc_olay_genebody_right_result_no = IVT_query(IVT_gbody_root, pos_large, junc_genebody_olayright, MAX_OVERLAP_EDGE_NUMBER);

			junc_near_LLedge_no = IVT_edges_lr(IVT_edge_root,pos_small,junc_nearest_LLedges,MAX_OVERLAP_EDGE_NUMBER,1);
			junc_near_LRedge_no = IVT_edges_lr(IVT_edge_root,pos_small,junc_nearest_LRedges,MAX_OVERLAP_EDGE_NUMBER,0);

			junc_near_RLedge_no = IVT_edges_lr(IVT_edge_root,pos_large,junc_nearest_RLedges,MAX_OVERLAP_EDGE_NUMBER,1);
			junc_near_RRedge_no = IVT_edges_lr(IVT_edge_root,pos_large,junc_nearest_RRedges,MAX_OVERLAP_EDGE_NUMBER,0);

			if(junc_olay_genebody_left_result_no == MAX_OVERLAP_EDGE_NUMBER|| junc_olay_genebody_right_result_no==MAX_OVERLAP_EDGE_NUMBER ||
                             junc_near_LLedge_no == MAX_OVERLAP_EDGE_NUMBER || junc_near_RLedge_no == MAX_OVERLAP_EDGE_NUMBER||
                             junc_near_LRedge_no == MAX_OVERLAP_EDGE_NUMBER || junc_near_RRedge_no == MAX_OVERLAP_EDGE_NUMBER
                           ){
				SUBREADprintf("WARNING: Your annotation file contains very many exons that start/end at the same location. Consider to increase MAX_OVERLAP_EDGE_NUMBER in readSummary.c to accomodate these exons for junction counting.\n");
			}

			jc_gene_status = cellCounts_determine_jcount_gene_transcript_report(cct_context , pos_small, pos_large, junc_olay_genebody_left_result_no, junc_olay_genebody_right_result_no, junc_genebody_olayleft, junc_genebody_olayright, gene_ids_str_SP1,gene_ids_str_SP2, transcript_ids_str_SP1, transcript_ids_str_SP2, dist_to_nearest_splice_side_str_SP1, dist_to_nearest_splice_side_str_SP2, strand_learned_from_FASTA, junc_near_LLedge_no, junc_near_LRedge_no, junc_near_RLedge_no, junc_near_RRedge_no, junc_nearest_LLedges,junc_nearest_LRedges,junc_nearest_RLedges,junc_nearest_RRedges, chro_small);
		}
		if(jc_gene_status == JC_STATUS_KNOWN) jc_retv_str = "KNOWN";
		if(jc_gene_status == JC_STATUS_NOVEL) jc_retv_str = "NOVEL";
		if(jc_gene_status == JC_STATUS_NA) jc_retv_str = "NA";
		char donorside[20], acceptorside[20];
		if( strand_learned_from_FASTA == -1){
			strcpy(donorside,"NA");
			strcpy(acceptorside,"NA");
		}else if(strand_learned_from_FASTA == 0 ){
			strcpy(donorside,"SP1");
			strcpy(acceptorside,"SP2");
		}else{
			strcpy(donorside,"SP2");
			strcpy(acceptorside,"SP1");
		}

		fprintf(ofp, "%s\t%s\t%s\t"
				"%s\t%s\t%s\t"
				"%s\t%d\t%s\t%s\t%d\t%s\t"
				"%s\t%s",
				gene_ids_str_SP1,gene_ids_str_SP2,transcript_ids_str_SP1,
				jc_retv_str,donorside,acceptorside,
				chro_small,pos_small,strand,  chro_large,pos_large,strand,
				dist_to_nearest_splice_side_str_SP1,  dist_to_nearest_splice_side_str_SP2
		);

		*(pos_small_str-1)='\t';
		*(pos_large_str-1)='\t';
		chro_large[-1]='\t';

		srInt_64 *countarry = HashTableGet(junction_table, key_list[ky_i]);
		for(sample_i = 1; sample_i <=cct_context-> sample_sheet_table -> numOfElements ; sample_i ++){
			srInt_64 count = countarry[sample_i];
			srInt_64 unsupp = count & 0xffffffffLL;
			count = count >> 32;
			#ifdef __MINGW32__
			fprintf(ofp,"\t%" PRId64, count);
			#else
			fprintf(ofp,"\t%lld\t%lld", count, unsupp);
			#endif
		}

		int wlen = fprintf(ofp, "\n");
		if(wlen < 1) disk_is_full = 1;
	}
	fclose(ofp);
	free(gene_names);
	free(key_list);
	free(junc_genebody_olayleft);
	free(junc_genebody_olayright);
	free(junc_nearest_LLedges);
	free(junc_nearest_RLedges);
	free(junc_nearest_LRedges);
	free(junc_nearest_RRedges);
	HashTableDestroy(junction_table);

	//print_in_box(80,0,PRINT_BOX_CENTER,"Found %llu junctions in all the input files.", merged_junction_table -> numOfElements);
	//print_in_box(80,0,0,"");

	if(disk_is_full){
		unlink(outfname);
		SUBREADprintf("ERROR: disk is full; no junction counting table is generated.\n");
	}
}

