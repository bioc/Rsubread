\name{Rsubreadnews}
\title{Rsubread News}
\encoding{UTF-8}

\section{Version 2.0.0}{\itemize{
\item
Rsubread package is ported to Windows OS.

\item
New function cellCounts(): generate UMI counts for Chromium 10X single-cell RNA-seq data.

\item
flattenGTF() function can merge or chop overlap features.

\item
Check and display the amount of memory available on the computer before starting read mapping.

\item
Optimize the data structure used in buildindex() function to reduce its memory use.

\item
qualityScores() function can optionally retrieve quality scores from all the reads.

\item
File paths included in column names of objects returned by featureCounts(), align(), subjunc() and propmapped() functions are removed or shortened where appropriate.

\item
featureCounts() will be terminated if both single-end and paired-end reads are found in the same input file.

\item
Limit on the length of input file names is increased to 1000 bytes for all functions.
}}


\section{Version 1.34.0}{\itemize{
\item
New functions: simReads() and scanFasta(). simReads() generates simulation RNA-seq reads for transcripts.

\item
align() and subjunc() estimate fragment length from mapped read pairs and use the estimated length to assist reporting the best alignment.

\item
align() and subjunc() prefer alignments with no indels included over indel-containing alignments when same number of matched bases are found.

\item
align() and subjunc() check if index files were successfully loaded before starting read mapping.

\item
align() and subjunc() detect indels arising from incorrect shifting of seed sequence when being mapped to low-complexity region and exclude such indels from read re-alignment and indel reporting.

\item
buildindex(), align() and subjunc() support gzipped FASTA format.

\item
featureCounts() allows mapped reads to have ‘*’ as their chromosome name.

\item
removeDupReads() supports BAM-format input and output.
}}


\section{Version 1.32.0}{\itemize{
\item
New function flattenGTF() that merges overlapping features into a single interval.

\item
New parameter for align() and subjunc(): sortReadsByCoordinates.

\item
New parameters for featureCounts(): readShiftType, readShiftSize and additionalAttributes.

\item
Specify strand protocol for each library individually in featureCounts().

\item
Much improved speed of align() and subjunc().

\item
align() and subjunc() return mapping statistics.

\item
Default setting of buildindex() is changed to building a one-block full index.
}}

