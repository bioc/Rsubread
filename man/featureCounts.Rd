\name{featureCounts}
\alias{featureCounts}
\title{featureCounts: a general-purpose read summarization function}
\description{This function assigns mapped sequencing reads to genomic features}
\usage{
featureCounts(files,

	# annotation
	annot.inbuilt="mm10",
	annot.ext=NULL,
	isGTFAnnotationFile=FALSE,
	GTF.featureType="exon",
	GTF.attrType="gene_id",
    additionalAttributes=NULL,
	chrAliases=NULL,
	
	# level of summarization
	useMetaFeatures=TRUE,
	
	# overlap between reads and features
	allowMultiOverlap=FALSE,
	minOverlap=1,
	fracOverlap=0,
	fracOverlapFeature=0,
	largestOverlap=FALSE,
	nonOverlap=NULL,
	nonOverlapFeature=NULL,

	# Read extension and reduction
	readExtension5=0,
	readExtension3=0,
	read2pos=NULL,
	
	# multi-mapping reads
	countMultiMappingReads=FALSE,

	# fractional counting
	fraction=FALSE,

	# long reads
	isLongRead=FALSE,

	# read filtering
	minMQS=0,
	splitOnly=FALSE,
	nonSplitOnly=FALSE,
	primaryOnly=FALSE,
	ignoreDup=FALSE,
	
	# strandness
	strandSpecific=0,
	
	# exon-exon junctions
	juncCounts=FALSE,
	genome=NULL,
	
	# parameters specific to paired end reads
	isPairedEnd=FALSE,
	requireBothEndsMapped=FALSE,
	checkFragLength=FALSE,
	minFragLength=50,
	maxFragLength=600,
	countChimericFragments=TRUE,	
	autosort=TRUE,
	
	# number of CPU threads
	nthreads=1,

	# read group
	byReadGroup=FALSE,

	# report assignment result for each read
	reportReads=NULL,
	reportReadsPath=NULL,

	# miscellaneous
	maxMOp=10,
	tmpDir=".",
	verbose=FALSE)
}
\arguments{
  \item{files}{ a character vector giving names of input files containing read mapping results. The files can be in either SAM format or BAM format. The file format is automatically detected by the function.}
  \item{annot.inbuilt}{ a character string specifying an in-built annotation used for read summarization. It has four possible values including \code{"mm10"}, \code{"mm9"}, \code{"hg38"} and \code{"hg19"}, corresponding to the NCBI RefSeq annotations for genomes `mm10', `mm9', `hg38' and `hg19', respectively. \code{"mm10"} by default. The in-built annotation has a SAF format (see below).}
  \item{annot.ext}{ A character string giving name of a user-provided annotation file or a data frame including user-provided annotation data. If the annotation is in GTF format, it can only be provided as a file. If it is in SAF format, it can be provided as a file or a data frame. See below for more details about SAF format annotation. If an annotation file is provided, it can be uncompressed or gzip compressed. Note that \code{annot.ext} will override \code{annot.inbuilt} if both provided.}
  \item{isGTFAnnotationFile}{ logical indicating whether the annotation provided via the \code{annot.ext} argument is in GTF format or not. \code{FALSE} by default. This option is only applicable when \code{annot.ext} is not \code{NULL}.}
  \item{GTF.featureType}{ a character string giving the feature type used to select rows in the GTF annotation which will be used for read summarization. \code{"exon"} by default. This argument is only applicable when \code{isGTFAnnotationFile} is \code{TRUE}. Feature types can be found in the third column of a GTF annotation.}
  \item{GTF.attrType}{ a character string giving the attribute type in the GTF annotation which will be used to group features (eg. exons) into meta-features (eg. genes). \code{"gene_id"} by default. This argument is only applicable when \code{isGTFAnnotationFile} is \code{TRUE}. Attributes can be found in the ninth column of a GTF annotation.}
  \item{additionalAttributes}{ a character vector specifying extra GTF attribute types that will also be included in the counting output. These attribute types will not be used to group features. \code{NULL} by default.}
  \item{chrAliases}{ a character string giving the name of a chromosome name alias file. This should be a two-column comma-delimited text file. Chromosome name aliases included in this file are used to match chr names in annotation with those in the reads. First column in the file should include chr names in the annotation and second column should include chr names in the reads. Chr names are case sensitive. No column header should be included in the file.}
  \item{useMetaFeatures}{ logical indicating whether the read summarization should be performed at the feature level (eg. exons) or meta-feature level (eg. genes). If \code{TRUE}, features in the annotation (each row is a feature) will be grouped into meta-features, using their values in the \code{GeneID} column in the SAF-format annotation file or using the \code{GTF.attrType} attribute in the GTF-format annotation file, and then reads will be assiged to the meta-features instead of the features. See below for more details.}
  \item{allowMultiOverlap}{ logical indicating if a read is allowed to be assigned to more than one feature (or meta-feature) if it is found to overlap with more than one feature (or meta-feature). \code{FALSE} by default.}
  \item{minOverlap}{ \code{integer} giving the minimum number of overlapped bases required for assigning a read to a feature (or a meta-feature). For assignment of read pairs (fragments), number of overlapping bases from each read in the same pair will be summed. If a negative value is provided, then a gap of up to specified size will be allowed between read and the feature that the read is assigned to. \code{1} by default.}
  \item{fracOverlap}{ \code{numeric} giving minimum fraction of overlapping bases in a read that is required for read assignment. Value should be within range [0,1]. \code{0} by default. Number of overlapping bases is counted from both reads if paired end. Both this option and \code{minOverlap} option need to be satisfied before a read can be assigned.}
  \item{fracOverlapFeature}{ \code{numeric} giving minimum fraction of bases included in a feature that is required for overlapping with a read or a read pair. Value should be within range [0,1]. \code{0} by default.}
  \item{largestOverlap}{ If \code{TRUE}, a read (or read pair) will be assigned to the feature (or meta-feature) that has the largest number of overlapping bases, if the read (or read pair) overlaps with multiple features (or meta-features).}
  \item{nonOverlap}{ \code{integer} giving the maximum number of bases in a read (or a read pair) that are allowed not to overlap the assigned feature. \code{NULL} by default (ie. no limit is set).}
  \item{nonOverlapFeature}{ \code{integer} giving the maximum number of non-overlapping bases allowed in a feature during read assignment. \code{NULL} by default (ie. no limit is set).}
  \item{readExtension5}{ \code{integer} giving the number of bases extended upstream from 5' end of each read. \code{0} by default.}
  \item{readExtension3}{ \code{integer} giving the number of bases extended downstream from 3' end of each read. \code{0} by default.}
  \item{read2pos}{ Specifying whether each read should be reduced to its 5' most base or 3' most base. It has three possible values: \code{NULL}, \code{5} (denoting 5' most base) and \code{3} (denoting 3' most base). The default value is \code{NULL}. With the default value, the whole read is used for summarization. When \code{read2pos} is set to \code{5} (or \code{3}), read summarization will be performed based on the 5' (or 3') most base position. \code{read2pos} can be used together with \code{readExtension5} and \code{readExtension3} parameters to set any desired length for reads.}
  \item{countMultiMappingReads}{ logical indicating if multi-mapping reads/fragments should be counted, \code{FALSE} by default. If \code{TRUE}, a multi-mapping read will be counted up to N times if it has N reported mapping locations. This function uses the `NH' tag to find multi-mapping reads.}
  \item{fraction}{ logical indicating if fractional counts are produced for multi-mapping reads and/or multi-overlapping reads. \code{FALSE} by default. See below for more details.}
  \item{isLongRead}{ logical indicating if input data contain long reads. This option should be set to \code{TRUE} if counting Nanopore or PacBio long reads.}
  \item{minMQS}{ \code{integer} giving the minimum mapping quality score a read must satisfy in order to be counted. For paired-end reads, at least one end should satisfy this criteria. \code{0} by default.}
  \item{splitOnly}{ logical indicating whether only split alignments (their CIGAR strings contain letter 'N') should be included for summarization. \code{FALSE} by default. Example split alignments are exon-spanning reads from RNA-seq data. \code{useMetaFeatures} should be set to \code{FALSE} and \code{allowMultiOverlap} should be set to \code{TRUE}, if the purpose of summarization is to assign exon-spanning reads to all their overlapping exons.}
  \item{nonSplitOnly}{ logical indicating whether only non-split alignments (their CIGAR strings do not contain letter 'N') should be included for summarization. \code{FALSE} by default.}
  \item{primaryOnly}{ logical indicating if only primary alignments should be counted. Primary and secondary alignments are identified using bit 0x100 in the Flag field of SAM/BAM files. If \code{TRUE}, all primary alignments in a dataset will be counted no matter they are from multi-mapping reads or not (ie. \code{countMultiMappingReads} is ignored).}
  \item{ignoreDup}{ logical indicating whether reads marked as duplicates should be ignored. \code{FALSE} by default. Read duplicates are identified using bit Ox400 in the FLAG field in SAM/BAM files. The whole fragment (read pair) will be ignored if paired end.}
  \item{strandSpecific}{ an integer vector indicating if strand-specific read counting should be performed. Length of the vector should be either \code{1} (meaning that the value is applied to all input files), or equal to the total number of input files provided. Each vector element should have one of the following three values: \code{0} (unstranded), \code{1} (stranded) and \code{2} (reversely stranded). Default value of this parameter is \code{0} (ie. unstranded read counting is performed for all input files).}
  \item{juncCounts}{ logical indicating if number of reads supporting each exon-exon junction will be reported. Junctions are identified from those exon-spanning reads in input data. \code{FALSE} by default.}
  \item{genome}{ a character string giving the name of a FASTA-format file that includes the reference sequences used in read mapping that produced the provided SAM/BAM files. \code{NULL} by default. This argument should only be used when \code{juncCounts} is \code{TRUE}. Note that providing reference sequences is optional when \code{juncCounts} is set to \code{TRUE}.}
  \item{isPairedEnd}{ logical indicating if paired-end reads are used. If \code{TRUE}, fragments (templates or read pairs) will be counted instead of individual reads. \code{FALSE} by default.}
  \item{requireBothEndsMapped}{ logical indicating if both ends from the same fragment are required to be successfully aligned before the fragment can be assigned to a feature or meta-feature. This parameter is only appliable when \code{isPairedEnd} is \code{TRUE}.}
  \item{checkFragLength}{ logical indicating if the two ends from the same fragment are required to satisify the fragment length criteria before the fragment can be assigned to a feature or meta-feature. This parameter is only appliable when \code{isPairedEnd} is \code{TRUE}. The fragment length criteria are specified via \code{minFragLength} and \code{maxFragLength}.}
  \item{minFragLength}{ \code{integer} giving the minimum fragment length for paired-end reads. \code{50} by default.}
  \item{maxFragLength}{ \code{integer} giving the maximum fragment length for paired-end reads. \code{600} by default. \code{minFragLength} and \code{maxFragLength} are only applicable when \code{isPairedEnd} is \code{TRUE}. Note that when a fragment spans two or more exons, the observed fragment length might be much bigger than the nominal fragment length.}
  \item{countChimericFragments}{ logical indicating whether a chimeric fragment, which has its two reads mapped to different chromosomes, should be counted or not. \code{TRUE} by default.}
  \item{autosort}{ logical specifying if the automatic read sorting is enabled. This option is only applicable for paired-end reads. If \code{TRUE}, reads will be automatically sorted by their names if reads from the same pair are found not to be located next to each other in the input. No read sorting will be performed if there are no such reads found.}
  \item{nthreads}{ \code{integer} giving the number of threads used for running this function. \code{1} by default.}
  \item{byReadGroup}{ logical indicating if read counting will be performed for each individual read group. \code{FALSE} by default.}
  \item{reportReads}{ output detailed read assignment results for each read (or fragment if paired end). The detailed assignment results can be saved in three different formats including \code{CORE}, \code{SAM} and \code{BAM} (note that these values are case sensitive). Default value of this option is \code{NULL}, indicating not to report assignment results for reads. \code{CORE} format represents a tab-delimited text file that contains four columns including read name, status (assigned or the reason if not assigned), number of targets and target list. A target is a feature or a meta-feature. Items in the target lists is separated by comma. If a read is not assigned, its number of targets will be set as -1. For each input file, a text file is generated and its name is the input file name added with `.featureCounts'. When \code{SAM} or \code{BAM} format is specified, the detailed assignment results will be saved to SAM and BAM format files. Names of generated files are the input file names added with `.featureCounts.sam' or `.featureCounts.bam'. Three tags are used to describe read assignment results: XS, XN and XT. Tag XS gives the assignment status. Tag XN gives number of targets. Tag XT gives comma separated target list.}
  \item{reportReadsPath}{ a character string specifying the directory where files including detailed assignment results are saved. If \code{NULL}, the results will be saved to the current working directory.}
  \item{maxMOp}{ \code{integer} giving the maximum number of `M' operations (matches or mis-matches) allowed in a CIGAR string. \code{10} by default. Both `X' and `=' operations are treated as `M' and adjacent `M' operations are merged in the CIGAR string.}
  \item{tmpDir}{ a character string specifying the directory under which intermediate files are saved (later removed). By default, current working directory is used.}
  \item{verbose}{ logical indicating if verbose information for debugging will be generated. This may include information such as unmatched chromosomes/contigs between reads and annotation.}
  }
\details{
\code{featureCounts} is a general-purpose read summarization function, which assigns to the genomic features (or meta-features) the mapped reads that were generated from genomic DNA and RNA sequencing.

This function takes as input a set of files containing read mapping results output from a read aligner (e.g. \code{\link{align}} or \code{\link{subjunc}}), and then assigns mapped reads to genomic features.
Both SAM and BAM format input files are accepted.

\code{featureCounts} accepts two annotation formats: SAF (Simplified Annotation Format) and GTF/GFF formats.
Specification of GTF/GFF format can be found at \url{https://genome.ucsc.edu/FAQ/FAQformat.html}.
SAF is a simplified annotation format and below shows an example of this format:

\preformatted{
GeneID		Chr	Start	End	Strand
497097		chr1	3204563	3207049	-
497097		chr1	3411783	3411982	-
497097		chr1	3660633	3661579	-
100503874	chr1	3637390	3640590	-
100503874	chr1	3648928	3648985	-
100038431	chr1	3670236	3671869	-
...
}

SAF annotation includes the following five required columns: \code{GeneID}, \code{Chr}, \code{Start}, \code{End} and \code{Strand}.
The \code{GeneID} column includes identifiers of features.
These identifiers can be integer or character string.
The \code{Chr} column includes chromosome names of features and these names should match the chromosome names includes in the provided SAM/BAM files.
The \code{Start} and \code{End} columns include start and end coordinates of features, respectively.
Both start and end coordinates are inclusive.
The \code{Strand} column indicates the strand of features (\code{"+"} or \code{"-"}).
Column names in a SAF annotation should be the same as those shown in the above example (case insensitive).
Columns can be in any order.
Extra columns are allowed to be added into the annotation.

In-built annotations, which were generated based on NCBI RefSeq gene annotations, are provided to faciliate convenient read summarization.
We provide in-built annotations for the following genomes: \code{"hg38"}, \code{"hg19"}, \code{"mm10"} and \code{"mm9"}.
The content of in-built annotations can be accessed via the \code{\link{getInBuiltAnnotation}} function.
These annotations have a SAF format.

The in-built annotations are a modified version of NCBI RefSeq gene annotations.
We downloaded the RefSeq gene annotations from NCBI ftp server (eg. RefSeq annotation for \code{mm10} was downloaded from \url{ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/mapview/seq_gene.md.gz}).
We then used these annotations to create our in-built annotations.
For each gene, we used its CDS (coding DNA sequence) and UTR (untranslated) regions provided in the original annotation to construct a list of exons, by merging and concatenating all CDSs and UTRs belonging to the same gene.
Exons within each gene include all chromosomal bases included in the original CDS and UTR regions, but they include each base only once (they might be included multiple times in the original CDSs and UTRs).
Also, exons within the same gene do not overlap with each other.

Users may provide an external annotation for read summarization via the \code{annot.ext} argument.
If the external annotation is in SAF format, it can be provdied as either a \code{data.frame} or a tab-delimited text file with proper column names included.
If it is in GTF/GFF format, it should be provided as a file only (and \code{isGTFAnnotationFile} should be set to \code{TRUE}).

\code{featureCounts} function uses the \code{GTF.attrType} attribute in a GTF/GFF annotation to group features to form meta-features when performing read summarization at meta-feature level.

The argument \code{useMetaFeatures} specifies whether read summarization should be performed at feature level or at meta-feature level.
A feature represents a continuous genomic region and a meta-feature is a group of features.
For instance, an exon is a feature and a gene comprising one or more exons is a meta-feature.
To assign reads to meta-features, \code{featureCounts} firstly groups into meta-features the features that have the same gene identifiers.
\code{featureCounts} looks for gene identifiers in \code{GeneID} column of a SAF annotation or by using \code{GTF.attrType} attribute in a GTF/GFF annotation.
Then for each read \code{featureCounts} searches for meta-features that have at least one feature that overlaps with the read.
A read might be found to overlap with more than one feature within the same meta-feature (eg. an exon-spanning read overlaps with more than one exon from the same gene), however this read will still be counted only once for the meta-feature.

RNA-seq reads are often summarized to meta-features to produce read counts for genes.
Further downstream analysis can then be carried out to discover differentially expressed genes.
Feature-level summarization of RNA-seq data often yields exon-level read counts, which is useful for investigating alternative splicing of genes.

\code{featureCounts} provides multiple options to count multi-mapping reads (reads mapping to more than one location in the reference genome).
Users can choose to ignore such reads by setting \code{countMultiMappingReads} to \code{FALSE}, or fully count every alignment reported for a multi-mapping read by setting \code{countMultiMappingReads} to \code{TRUE} (each alignment carries \code{1} count), or count each alignment fractionally by setting both \code{countMultiMappingReads} and \code{fraction} to \code{TRUE} (each alignment carries \code{1/x} count where \code{x} is the total number of alignments reported for the read).

\code{featureCounts} also provides multiple options to count multi-overlapping reads (reads overlapping with more than one meta-feature when conducting meta-feature-level summarization or overlapping with more than one feature when conducting feature-level summarization).
Users can choose to ignore such reads by setting \code{allowMultiOverlap} to \code{FALSE}, or fully count them for each overlapping meta-feature/feature by setting \code{allowMultiOverlap} to \code{TRUE} (each overlapping meta-feature/feature receives a count of 1 from a read), or assign a fractional count to each overlapping meta-feature/feature by setting both \code{allowMultiOverlap} and \code{fraction} to \code{TRUE} (each overlapping meta-feature/feature receives a count of \code{1/y} from a read where \code{y} is the total number of meta-features/features overlapping with the read).

If a read is both multi-mapping and multi-overlapping, then each overlapping meta-feature/feature will receive a fractional count of \code{1/(x*y)} if \code{countMultiMappingReads}, \code{allowMultiOverlap} and \code{fraction} are all set to \code{TRUE}.

When \code{isPairedEnd} is \code{TRUE}, fragments (pairs of reads) instead of reads will be counted.
\code{featureCounts} function checks if reads from the same pair are adjacent to each other (this could happen when reads were for example sorted by their mapping locations), and it automatically reorders those reads that belong to the same pair but are not adjacent to each other in the input read file.
}
\value{ 
A list with the following components:
\item{counts}{ a data matrix containing read counts for each feature or meta-feature for each library.}
\item{counts_junction (optional)}{ a data frame including primary and secondary genes overlapping an exon junction, position information for the left splice site (`Site1') and the right splice site (`Site2') of a junction (including chromosome name, coordinate and strand) and number of supporting reads for each junction in each library. Both primary and secondary genes overlap at least one of the two splice sites of a junction. Secondary genes do not overlap more splice sites than the primary gene. When the primary and secondary genes overlap same number of splice sites, the gene with the smallest leftmost base position is selected as the primary gene. For each junction, no more than one primary gene is reported but there might be more than one secondary genes reported. If a junction does not overlap any genes, it has a missing value in the fields of primary gene and secondary gene. Note that this data frame is only generated when \code{juncCounts} is \code{TRUE}.}
\item{annotation}{ a data frame with six columns including \code{GeneID}, \code{Chr}, \code{Start}, \code{End} and \code{Length}. When read summarization was performed at feature level, each row in the data frame is a feature and columns in the data frame give the annotation information for the features. When read summarization was performed at meta-feature level, each row in the data frame is a meta-feature and columns in the data frame give the annotation information for the features included in each meta feature except the \code{Length} column. For each meta-feature, the \code{Length} column gives the total length of genomic regions covered by features included in that meta-feature.   Note that this length will be less than the sum of lengths of features included in the meta-feature when there are features overlapping with each other. Also note the \code{GeneID} column gives Entrez gene identifiers when the in-built annotations are used.}
\item{targets}{ a character vector giving sample information.}
\item{stat}{ a data frame giving numbers of unassigned reads and the reasons why they are not assigned (eg. ambiguity, multi-mapping, secondary alignment, mapping quality, fragment length, chimera, read duplicate, non-junction and so on), in addition to the number of successfully assigned reads for each library.}
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014.
}
\author{Wei Shi and Yang Liao}
%\note{}
\seealso{
\code{\link{getInBuiltAnnotation}}}
\examples{\dontrun{
library(Rsubread)

# Summarize SAM format single-end reads using built-in RefSeq annotation for mouse genome mm10:
featureCounts(files="mapping_results_SE.sam",annot.inbuilt="mm10")

# Summarize single-end reads using a user-provided GTF annotation file:
featureCounts(files="mapping_results_SE.sam",annot.ext="annotation.gtf",
isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id")

# Summarize single-end reads using 5 threads:
featureCounts(files="mapping_results_SE.sam",nthreads=5)

# Summarize BAM format single-end read data:
featureCounts(files="mapping_results_SE.bam")

# Perform strand-specific read counting (strandSpecific=2 if reversely stranded):
featureCounts(files="mapping_results_SE.bam",strandSpecific=1)

# Summarize paired-end reads and counting fragments (instead of reads):
featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE)

# Count fragments satisfying the fragment length criteria, eg. [50bp, 600bp]:
featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE,
checkFragLength=TRUE,minFragLength=50,maxFragLength=600)

# Count fragments that have both ends successfully aligned without checking the fragment length:
featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE,requireBothEndsMapped=TRUE)

# Exclude chimeric fragments from fragment counting:
featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE,countChimericFragments=FALSE)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
