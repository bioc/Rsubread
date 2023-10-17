\name{Rsubreadnews}
\title{Rsubread News}
\encoding{UTF-8}
\section{Version 2.16.0}{\itemize{

\item
Add support for ARM64 platforms.

\item
Add support for CBCL format in cellCounts.
}}
\section{Version 2.12.0}{\itemize{
\item
Improve the data structure used by cellCounts to further improve its speed.

\item
More checks for input parameters to prevent cellCounts from crashing.

\item
Improve the screen output of cellCounts.
}}

\section{Version 2.10.0}{\itemize{
\item
Added inbuilt RefSeq annotation for mm39 (mouse genome Build 39).

\item
Streamlined the mapping and counting processes in cellCounts.

\item
Added support for processing dual-index 10x data in cellCounts.
}}


\section{Version 2.6.0}{\itemize{
\item
Improved the speed of cellCounts and also reduced its memory use.

\item
Added a parameter 'umi.cutoff' to cellCounts to call all the cells that had a total UMI count greater than the specified threshold.

\item
Added support for FASTQ-format read input in CellCounts.
}}


\section{Version 2.4.0}{\itemize{
\item
The 'isPairedEnd' parameter in featureCounts() is now used to check if the type of input reads (paired-end or single-end) is correctly specified.

\item
A new parameter 'countReadPairs' was added to featureCounts to specify if read pairs should be counted for paired-end read data.

\item
Changes to the input parameters and output of cellCounts() function. CellCounts will generate a Sample Sheet for samples included in the scRNA-seq data, based on the sample index set name provided by the user. Structure of the List object returned by cellCounts() is also simplified. cellCounts() now also outputs a BAM file and a gzipped FASTQ file including raw reads for each sample.
}}


\section{Version 2.2.0}{\itemize{
\item
Improve cellCounts() on the identification of cell barcodes arising from ambient RNAs.
}}


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

