#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "softclip-test.h"


/**
 * Calculates soft clipping based on mismatch density.
 * 
 * Returns a pointer to a SoftClipResult struct allocated on the heap.
 * The caller must free this memory.
 */
SoftClipResult* calculate_soft_clipping(
    void * context,
    unsigned int pos,
    const char* cigar,
    const char* read_seq,
    unsigned int perf_start,
    unsigned int perf_end,
    char (*get_ref_base)(unsigned int pos, void * context)
) {
    // Allocate result struct on the heap
    SoftClipResult* result = (SoftClipResult*)calloc(1, sizeof(SoftClipResult));
    if (!result) return NULL; // Allocation failure check

    int read_len = strlen(read_seq);
    if (read_len == 0) {
        // Empty read, return empty result
        return result; 
    }

    // Temporary arrays for per-base analysis
    unsigned int* ref_coords = (unsigned int*)malloc(read_len * sizeof(unsigned int));
    bool* is_mismatch = (bool*)malloc(read_len * sizeof(bool)); // True for SNP or Insertion
    bool* is_insertion = (bool*)malloc(read_len * sizeof(bool)); // Specific flag for I
    
    if (!ref_coords || !is_mismatch || !is_insertion) {
        free(result); if(ref_coords) free(ref_coords); 
        if(is_mismatch) free(is_mismatch); if(is_insertion) free(is_insertion);
        return NULL;
    }

    // 1. Parse CIGAR to map Read Indices -> Reference Coordinates
    unsigned int current_ref = pos;
    int current_read = 0;
    const char* p = cigar;
    
    while (*p) {
        char* end_ptr;
        long op_len = strtol(p, &end_ptr, 10);
        char op = *end_ptr;
        p = end_ptr + 1;

        if (op == 'M' || op == '=' || op == 'X') {
            for (int k = 0; k < op_len; k++) {
                if (current_read < read_len) {
                    ref_coords[current_read] = current_ref;
                    is_insertion[current_read] = false;
                    
                    char ref_base = get_ref_base(current_ref,context);
                    char read_base = read_seq[current_read];
                    
                    if (toupper(read_base) != toupper(ref_base)) {
                        is_mismatch[current_read] = true;
                    } else {
                        is_mismatch[current_read] = false;
                    }
                    current_read++;
                }
                current_ref++;
            }
        } else if (op == 'I' || op == 'S') {
            for (int k = 0; k < op_len; k++) {
                if (current_read < read_len) {
                    ref_coords[current_read] = current_ref; 
                    is_insertion[current_read] = true;
                    // Insertions count as mismatches for window density checks
                    is_mismatch[current_read] = true; 
                    current_read++;
                }
            }
        } else if (op == 'D' || op == 'N') {
            current_ref += op_len;
        }
    }

    // 2. Scan Left (from perfect segment start backwards)
    int new_read_start = 0; 
    
    for (int i = perf_start - 1; i >= 0; i--) {
        int start_win = i - SC_SOFT_CLIPPING_WINDOW_SIZE + 1;
        if (start_win < 0) start_win = 0;
        
        int mismatch_count = 0;
        for (int w = start_win; w <= i; w++) {
            if (is_mismatch[w]) mismatch_count++;
        }
        
        if (mismatch_count > SC_MAX_MISMATCHED_BASES) {
            // Window Failed. Find last matched base farthest from perf_start (smallest index)
            int kept_base_index = -1;
            for (int w = start_win; w <= i; w++) {
                if (!is_mismatch[w]) {
                    kept_base_index = w;
                    break; // Found leftmost match
                }
            }
            if (kept_base_index != -1) new_read_start = kept_base_index;
            else new_read_start = i + 1;
            break;
        }
    }

    // 3. Scan Right (from perfect segment end forwards)
    int new_read_end = read_len;
    
    for (int i = perf_end; i < read_len; i++) {
        int end_win = i + SC_SOFT_CLIPPING_WINDOW_SIZE - 1;
        if (end_win >= read_len) end_win = read_len - 1;
        
        int mismatch_count = 0;
        for (int w = i; w <= end_win; w++) {
            if (is_mismatch[w]) mismatch_count++;
        }
        
        if (mismatch_count > SC_MAX_MISMATCHED_BASES) {
            // Window Failed. Find last matched base farthest from perf_end (largest index)
            int kept_base_index = -1;
            for (int w = end_win; w >= i; w--) {
                if (!is_mismatch[w]) {
                    kept_base_index = w;
                    break; // Found rightmost match
                }
            }
            if (kept_base_index != -1) new_read_end = kept_base_index + 1;
            else new_read_end = i;
            break;
        }
    }

    // 4. Calculate Statistics
    int left_clip = new_read_start;
    int right_clip = read_len - new_read_end;
    if (right_clip < 0) right_clip = 0;
    
    result->num_clipped = left_clip + right_clip;
    
    if (new_read_start >= new_read_end) {
        // Edge case: All clipped
        sprintf(result->new_cigar, "%dS", read_len);
        result->new_pos = pos;
        result->num_clipped = read_len;
        result->num_matched = 0;
        result->num_mismatched = 0;
        result->num_inserted = 0;
        
        free(ref_coords); free(is_mismatch); free(is_insertion);
        return result;
    }

    // Tally stats in kept region
    for (int k = new_read_start; k < new_read_end; k++) {
        if (is_insertion[k]) {
            result->num_inserted++;
        } else if (is_mismatch[k]) {
            result->num_mismatched++;
        } else {
            result->num_matched++;
        }
    }

    // 5. Construct Result CIGAR and Pos
    result->new_pos = ref_coords[new_read_start];
    char* c_ptr = result->new_cigar;
    
    // Add Left Clipping
    if (left_clip > 0) {
        c_ptr += sprintf(c_ptr, "%dS", left_clip);
    }
    
    // Internal CIGAR
    current_read = 0;
    p = cigar;
    while (*p) {
        char* end_ptr;
        long op_len = strtol(p, &end_ptr, 10);
        char op = *end_ptr;
        p = end_ptr + 1;
        
        int op_width = (op == 'D' || op == 'N') ? 0 : (int)op_len;
        int op_start = current_read;
        int op_end = current_read + op_width;
        
        int keep_s = (op_start > new_read_start) ? op_start : new_read_start;
        int keep_e = (op_end < new_read_end) ? op_end : new_read_end;
        
        if (op_width > 0) {
            if (keep_e > keep_s) {
                c_ptr += sprintf(c_ptr, "%d%c", (keep_e - keep_s), op);
            }
        } else {
            // Keep deletion only if strictly inside the kept region
            if (current_read >= new_read_start && current_read < new_read_end) {
                c_ptr += sprintf(c_ptr, "%ld%c", op_len, op);
            }
        }
        current_read += op_width;
    }
    
    // Add Right Clipping
    if (right_clip > 0) {
        c_ptr += sprintf(c_ptr, "%dS", right_clip);
    }

    // Cleanup local arrays
    free(ref_coords);
    free(is_mismatch);
    free(is_insertion);
    
    return result;
}

char test_clipping_mock_genome_access(unsigned int pos, void * cont) {
    // Mock reference: A simple repeated pattern or specific bases
    // Let's assume the reference is all 'A's for simplicity of testing mismatches
    char * o="AAAAAAAAAACCAAAAAAAAAA";
    return o[pos-1]; 
}

int main_test_for_clipping() {
    // Scenario: Read is 20bp. Perfect middle [5, 15).
    // Ends are messy. 
    // Read: TTTTT (mismatch) AAAAAAAAAA (match) TTTTT (mismatch)
    // Ref:  AAAAA            AAAAAAAAAA         AAAAA
    
    unsigned int pos = 1;
    const char* cigar = "10M2D1M1I10M";
    const char* read = "CCAAAAAAAAATAAAAACCCC";
    unsigned int p_start = 6;
    unsigned int p_end = 7;
    
    SoftClipResult * res = calculate_soft_clipping(NULL,pos, cigar, read, p_start, p_end, test_clipping_mock_genome_access); 
    
    printf("Old Pos: %u, Old Cigar: %s\n", pos, cigar);
    printf("New Pos: %u, New Cigar: %s\n", res->new_pos, res->new_cigar);
    free(res);
    
    return 0;
}
