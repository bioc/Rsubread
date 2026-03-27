#ifndef __SOFTCLIP_TEST_H_
#define __SOFTCLIP_TEST_H_

// --- Configuration Constants ---
#define SC_SOFT_CLIPPING_WINDOW_SIZE 5
#define SC_MAX_CIGAR_LEN 1024

typedef struct {
    char new_cigar[SC_MAX_CIGAR_LEN];
    unsigned int new_pos;
    
    // Statistics for the final aligned (kept) segment + clipped amounts
    int num_matched;    // Bases matching reference
    int num_mismatched; // Bases mapping to reference but different (SNPs)
    int num_inserted;   // Bases that are insertions (I)
    int num_clipped;    // Bases that are soft clipped (S)
} SoftClipResult;

/**
 * Calculates soft clipping with "Keep Farthest Match" logic.
 *
 * @param context       Context data that is passed to the get_ref_base function.
 * @param pos           Original mapping position (0-based ref coordinate).
 * @param cigar         Original CIGAR string.
 * @param read_seq      Read sequence string.
 * @param perf_start    Start index of perfect segment (inclusive).
 * @param perf_end      End index of perfect segment (exclusive).
 * @param get_ref_base  Callback to get reference base (returns char).
 */
SoftClipResult * calculate_soft_clipping(
    void * context,
    unsigned int pos,
    const char* cigar,
    const char* read_seq,
    unsigned int perf_start,
    unsigned int perf_end,
    int max_mismatched_bases_in_window,
    char (*get_ref_base)(unsigned int pos, void * context)
);
#endif
