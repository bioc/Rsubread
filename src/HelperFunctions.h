#ifndef __HELPER_FUNCTIONS_H_
#define __HELPER_FUNCTIONS_H_


// This function parses CIGAR_Str and extract the relative starting points and lengths of all sections (i.e., the sections of read that are separated by 'N').
// CIGAR_Str is a CIGAR string containing 'S', 'M', 'I', 'D' and 'N' operations. Other operations are all ignored. The length of CIGAR_Str should be always less than 100 bytes or "-1" is returned.
// Staring_Points and Section_Length are empty arrays to write the sections. The minimum length of each array is 6 items.
// The length of a section is its length on the chromosome, namely 'I' is ignored but 'D' is added into the length.
// This function ignores all sections from the 7-th.

// This function returns the number of sections found in the CIGAR string. It returns -1 if the CIGAR string cannot be parsed.

int RSubread_parse_CIGAR_string(const char * CIGAR_Str, unsigned int * Staring_Points, unsigned short * Section_Length);

#endif
