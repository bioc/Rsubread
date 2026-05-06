#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sambam-file.h"
#include "hashtable.h"

#ifdef _WIN32
static char *strsep(char **stringp, const char *delim) {
    char *start = *stringp;
    char *p;

    if (!start) return NULL;
    p = strpbrk(start, delim);
    if (p) {
        *p = '\0';
        *stringp = p + 1;
    } else {
        *stringp = NULL;
    }
    return start;
}
#endif

#define MAX_LINE (10*10*1024)
#define HASH_SIZE (2048*1024)

#ifdef MAKE_STANDALONE
#define EXTRACT_JUNCTION_MAIN main
#else
#define EXTRACT_JUNCTION_MAIN extract_junction_from_BAM_main
#endif

typedef struct {
    char *barcode;
    int index; 
} BarcodeEntry;

typedef struct UMIEntry {
    char *umi;
    struct UMIEntry *next;
} UMIEntry;

typedef struct JunctionEntry {
    char *key;
    int *counts;             // Only storing non-zero values here is better for RAM, 
                             // but for simplicity, we use the array and skip zeros during output.
    UMIEntry **seen_umis;    
    int row_id;
    struct JunctionEntry *next;
} JunctionEntry;

#define TOOL_CELLRANGER 100
#define TOOL_STARSOLO   200

int tool;
BarcodeEntry **barcode_hash;
int barcode_count = 0;
char **barcode_list;
JunctionEntry **junction_hash;
int junction_count = 0;

char *prefixed_filename(const char *prefix, const char *suffix) {
    size_t len = strlen(prefix) + strlen(suffix) + 1;
    char *filename = malloc(len);
    if (!filename) {
        perror("Filename allocation error");
        exit(1);
    }
    snprintf(filename, len, "%s%s", prefix, suffix);
    return filename;
}

// Check if all characters in the string are the same
int is_homopolymer(const char *s) {
    if (!s || !*s) return 0;
    char first = s[0];
    for (int i = 1; s[i] != '\0'; i++) {
        if (s[i] != first) return 0;
    }
    return 1;
}

unsigned int hash(const char *str) {
    unsigned int h = 5381;
    int c;
    while ((c = *str++)) h = ((h << 5) + h) + c;
    return h % HASH_SIZE;
}

int load_barcodes(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) return -1;
    barcode_hash = calloc(HASH_SIZE, sizeof(BarcodeEntry*));
    char line[128];
    while (fgets(line, sizeof(line), f)) {
        line[strcspn(line, "\n")] = 0;
        unsigned int h = hash(line);
        while (barcode_hash[h]) h = (h + 1) % HASH_SIZE;
        barcode_hash[h] = malloc(sizeof(BarcodeEntry));
        barcode_hash[h]->barcode = strdup(line);
        barcode_hash[h]->index = ++barcode_count; // 1-indexed for Matrix Market
    }
    fclose(f);
    return 0;
}

int get_barcode_index(const char *cb) {
    unsigned int h = hash(cb);
    while (barcode_hash[h]) {
        if (strcmp(barcode_hash[h]->barcode, cb) == 0) return barcode_hash[h]->index;
        h = (h + 1) % HASH_SIZE;
    }
    return -1;
}

void process_junction(const char *j_key, int cb_idx, const char *ub) {
    unsigned int h = hash(j_key);
    JunctionEntry *entry = junction_hash[h];
    while (entry && strcmp(entry->key, j_key) != 0) entry = entry->next;

    if (!entry) {
        entry = malloc(sizeof(JunctionEntry));
        entry->key = strdup(j_key);
        entry->counts = calloc(barcode_count + 1, sizeof(int));
        entry->seen_umis = calloc(barcode_count + 1, sizeof(UMIEntry*));
        entry->row_id = ++junction_count; // 1-indexed for Matrix Market
        entry->next = junction_hash[h];
        junction_hash[h] = entry;
    }

    // UMI de-duplication
    UMIEntry *curr = entry->seen_umis[cb_idx];
    while (curr) {
        if (strcmp(curr->umi, ub) == 0) return;
        curr = curr->next;
    }

    UMIEntry *new_umi = malloc(sizeof(UMIEntry));
    new_umi->umi = strdup(ub);
    new_umi->next = entry->seen_umis[cb_idx];
    entry->seen_umis[cb_idx] = new_umi;
    entry->counts[cb_idx]++;
}

int parse_sam(char * bamfilename) {
    char line[MAX_LINE];
    tool=-1;
    junction_hash = calloc(HASH_SIZE, sizeof(JunctionEntry*));
    SamBam_FILE *bamfp = SamBam_fopen(bamfilename, SAMBAM_FILE_BAM);
    if(!bamfp) return -1;

    while (1){
        void * eofnull = SamBam_fgets(bamfp, line, MAX_LINE, 0);
        if(eofnull==NULL)break;
        if (line[0] == '@'){
            if(tool<0){
                if(strstr(line,"ID:cellranger")) tool=TOOL_CELLRANGER;
                if(strstr(line,"--soloType")) tool=TOOL_STARSOLO;
            }
            continue;
        }
        
        if(tool<0){
            SUBREADprintf("ERROR: unable to determine the aligner (STARsolo or Cell Ranger) using the BAM header.");   
            return -1;
        }
//fprintf(stderr,"BAMLINE: %s\n", line);
        char *line_ptr = line;
        strsep(&line_ptr, "\t"); // qname
        strsep(&line_ptr, "\t"); // flag
        char *rname = strsep(&line_ptr, "\t");
        char *pos_s = strsep(&line_ptr, "\t");
        strsep(&line_ptr, "\t"); // mapq
        char *cigar = strsep(&line_ptr, "\t");
        for(int i=0; i<5; i++) strsep(&line_ptr, "\t"); // skip to tags
        
        char *cb = NULL, *ub = NULL, *ur = NULL;
        int gx_present = 0;
        int xf_val = 0;
        char *tag;
        while ((tag = strsep(&line_ptr, "\t"))) {
            if      (strncmp(tag, "CB:Z:", 5) == 0 && tag[5]!='-') cb = tag + 5;
            else if (strncmp(tag, "UB:Z:", 5) == 0 && tag[5]!='-') ub = tag + 5;
            else if (strncmp(tag, "UR:Z:", 5) == 0 && tag[5]!='-') ur = tag + 5;
            else if (strncmp(tag, "GX:Z:", 5) == 0 && tag[5]!='-') gx_present = 1;
            else if (strncmp(tag, "xf:i:", 5) == 0 && tag[5]!='-') xf_val = atoi(tag + 5);
        }
        if (!cb || !ub || !ur) continue;
        if (is_homopolymer(ur) || strstr(ur,"N")) continue;
        int cb_idx = get_barcode_index(cb);
        if (cb_idx == -1) continue;
        if (tool == TOOL_CELLRANGER && gx_present) {
            int condition = (xf_val & 1) && !(xf_val & 2) && !(xf_val & 4) && (xf_val & 16) && !(xf_val & 32);
            if (!condition) continue;
        }

        long curr_pos = atol(pos_s);
        char *c = cigar;
        while (*c) {
            char *end;
            long len = strtol(c, &end, 10);
            if (*end == 'N') {
                char j_key[256];
                snprintf(j_key, sizeof(j_key), "%s:%ld^%ld", rname, curr_pos-1, curr_pos + len);
                process_junction(j_key, cb_idx, ub);
            }
            if (strchr("MND=X", *end)) curr_pos += len;
            c = end + 1;
        }
    }
    SamBam_fclose(bamfp);
    return 0;
}

int EXTRACT_JUNCTION_MAIN(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "\nUsage: %s barcodes.txt output_prefix bam_file_name\nThe tool (Cell Ranger or STARsolo) will be inferenced from the BAM header.\nColumns are defined in the cell barcode input. Matrix output is written to output_prefix.mtx and junctions are written to output_prefix.junctions.tsv.\n\n", argv[0]); return 1; }

    int has_error = load_barcodes(argv[1]);
    has_error |= parse_sam(argv[3]);

    if(has_error){ // R will check input files. 
        return -1;
    }

    char *matrix_filename = prefixed_filename(argv[2], ".mtx");
    char *junction_filename = prefixed_filename(argv[2], ".junctions.tsv");

    // 1. Calculate non-zero entries (NNZ)
    long long nnz = 0;
    for (int i = 0; i < HASH_SIZE; i++) {
        for (JunctionEntry *entry = junction_hash[i]; entry; entry = entry->next) {
            for (int j = 1; j <= barcode_count; j++) {
                if (entry->counts[j] > 0) nnz++;
            }
        }
    }

    // 2. Output Matrix Market Header
    FILE *f_matrix = fopen(matrix_filename, "w");
    if (!f_matrix) {
        perror("Matrix output file error");
        exit(1);
    }
    fprintf(f_matrix, "%%%%MatrixMarket matrix coordinate integer general\n");
    fprintf(f_matrix, "%% Rows: Junctions, Columns: Barcodes\n");
    fprintf(f_matrix, "%d %d %lld\n", junction_count, barcode_count, nnz);

    // 3. Output Entries
    for (int i = 0; i < HASH_SIZE; i++) {
        for (JunctionEntry *entry = junction_hash[i]; entry; entry = entry->next) {
            for (int j = 1; j <= barcode_count; j++) {
                if (entry->counts[j] > 0) {
                    fprintf(f_matrix, "%d %d %d\n", entry->row_id, j, entry->counts[j]);
                }
            }
        }
    }
    fclose(f_matrix);

    // 4. Write the junction names to a separate file for reference
    FILE *f_genes = fopen(junction_filename, "w");
    if (!f_genes) {
        perror("Junction output file error");
        exit(1);
    }
    for (int i = 0; i < HASH_SIZE; i++) {
        for (JunctionEntry *entry = junction_hash[i]; entry; entry = entry->next) {
            fprintf(f_genes, "%d\t%s\n", entry->row_id, entry->key);
        }
    }
    fclose(f_genes);

    free(matrix_filename);
    free(junction_filename);

    return 0;
}

