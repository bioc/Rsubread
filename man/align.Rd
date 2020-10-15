\name{align}
\alias{align}
\alias{subjunc}
\title{Align Sequence Reads to a Reference Genome via Seed-and-Vote}
\description{The \code{align} function can align both DNA and RNA sequencing reads. \code{Subjunc} is an RNA-seq aligner and it reports full alignment of each read (\code{align} reports partial alignment for exon spanning reads).}
\usage{
align(

    # index for reference sequences
    index,
    
    # input reads and output
    readfile1,
    readfile2 = NULL,
    type = "rna",
    input_format = "gzFASTQ",
    output_format = "BAM",
    output_file = paste(readfile1,"subread",output_format,sep="."),
    
    # offset value added to Phred quality scores of read bases
    phredOffset = 33,
    
    # thresholds for mapping
    nsubreads = 10,
    TH1 = 3,
    TH2 = 1,
    maxMismatches = 3,
    
    # unique mapping and multi-mapping
    unique = FALSE,
    nBestLocations = 1,
    
    # indel detection
    indels = 5,
    complexIndels = FALSE,
    
    # read trimming
    nTrim5 = 0,
    nTrim3 = 0,
    
    # distance and orientation of paired end reads
    minFragLength = 50,
    maxFragLength = 600,
    PE_orientation = "fr",
    
    # number of CPU threads
    nthreads = 1,
    
    # read group
    readGroupID = NULL,
    readGroup = NULL,
    
    # read order
    keepReadOrder = FALSE,
    sortReadsByCoordinates = FALSE,
    
    # color space reads
    color2base = FALSE,
    
    # dynamic programming
    DP_GapOpenPenalty = -1,
    DP_GapExtPenalty = 0,
    DP_MismatchPenalty = 0,
    DP_MatchScore = 2,
    
    # detect structural variants
    detectSV = FALSE,
    
    # gene annotation
    useAnnotation = FALSE,
    annot.inbuilt = "mm10",
    annot.ext = NULL,
    isGTF = FALSE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    chrAliases = NULL)
    
subjunc(

    # index for reference sequences
    index,
    
    # input reads and output
    readfile1,
    readfile2 = NULL,
    input_format = "gzFASTQ",
    output_format = "BAM",
    output_file = paste(readfile1,"subjunc",output_format,sep = "."),
    
    # offset value added to Phred quality scores of read bases
    phredOffset = 33,
    
    # thresholds for mapping
    nsubreads = 14,
    TH1 = 1,
    TH2 = 1,
    maxMismatches = 3,
    
    # unique mapping and multi-mapping
    unique = FALSE,
    nBestLocations = 1,
    
    # indel detection
    indels = 5,
    complexIndels = FALSE,
    
    # read trimming
    nTrim5 = 0,
    nTrim3 = 0,
    
    # distance and orientation of paired end reads
    minFragLength = 50,
    maxFragLength = 600,
    PE_orientation = "fr",
    
    # number of CPU threads
    nthreads = 1,
    
    # read group
    readGroupID = NULL,
    readGroup = NULL,
    
    # read order
    keepReadOrder = FALSE,
    sortReadsByCoordinates = FALSE,
    
    # color space reads
    color2base = FALSE,
    
    # dynamic programming
    DP_GapOpenPenalty = -1,
    DP_GapExtPenalty = 0,
    DP_MismatchPenalty = 0,
    DP_MatchScore = 2,
    
    # detect all junctions including gene fusions
    reportAllJunctions = FALSE,
    
    # gene annotation
    useAnnotation = FALSE,
    annot.inbuilt = "mm10",
    annot.ext = NULL,
    isGTF = FALSE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    chrAliases = NULL)
}
\arguments{
  \item{index}{ character string giving the basename of index file. Index files should be located in the current directory.}
  \item{readfile1}{ a character vector including names of files that include sequence reads to be aligned. For paired-end reads, this gives the list of files including first reads in each library. File format is FASTQ/FASTA by default. See \code{input_format} option for more supported formats.}
  \item{readfile2}{ a character vector giving names of files that include second reads in paired-end read data. Files included in \code{readfile2} should be in the same order as their mate files included in \code{readfile1} .\code{NULL} by default.}
  \item{type}{ a character string or an integer giving the type of sequencing data. Possible values include \code{rna} (or \code{0}; RNA-seq data) and \code{dna} (or \code{1}; genomic DNA-seq data such as WGS, WES, ChIP-seq data etc.). Character strings are case insensitive.}
  \item{input_format}{ character string specifying format of read input files. \code{gzFASTQ} by default (this also includes FASTQ, FASTA, gzipped FASTQ and gzipped FASTA formats). Other supported formats include \code{SAM} and \code{BAM}. Character values are case insensitive.}
  \item{output_format}{ character string specifying format of output file. \code{BAM} by default. Acceptable formats include \code{SAM} and \code{BAM}.}
  \item{output_file}{ a character vector specifying names of output files. By default, names of output files are set as the file names provided in \code{readfile1} added with an suffix string.}
  \item{phredOffset}{ numeric value added to base-calling Phred scores to make quality scores (represented as ASCII letters). Possible values include \code{33} and \code{64}. By default, \code{33} is used.}
  \item{nsubreads}{ numeric value giving the number of subreads extracted from each read.}
  \item{TH1}{ numeric value giving the consensus threshold for reporting a hit. This is the threshold for the first reads if paired-end read data are provided.}
  \item{TH2}{ numeric value giving the consensus threhold for the second reads in paired-end data.}
  \item{maxMismatches}{ numeric value giving the maximum number of mis-matched bases allowed in the alignment. \code{3} by default. Mis-matches found in soft-clipped bases are not counted.}
  \item{unique}{ logical indicating if only uniquely mapped reads should be reported. A uniquely mapped read has one single mapping location that has less mis-matched bases than any other candidate locations. By default, multi-mapping reads will be reported in addition to uniquely mapped reads. Number of alignments reported for each multi-mapping read is determined by the \code{nBestLocations} parameter.}
  \item{nBestLocations}{ numeric value specifying the maximal number of equally-best mapping locations that will be reported for a multi-mapping read. 1 by default. The allowed value is between 1 to 16 (inclusive). In the mapping output, `NH' tag is used to indicate how many alignments are reported for the read and `HI' tag is used for numbering the alignments reported for the same read. This argument is only applicable when \code{unique} option is set to \code{FALSE}.}
  \item{indels}{ numeric value giving the maximum number of insertions/deletions allowed during the mapping. \code{5} by default.}
  \item{complexIndels}{ logical indicating if complex indels will be detected. If \code{TRUE}, the program will try to detect multiple short indels that occurs concurrently in a small genomic region (indels could be as close as 1bp apart).}
  \item{nTrim5}{ numeric value giving the number of bases trimmed off from 5' end of each read. \code{0} by default.}
  \item{nTrim3}{ numeric value giving the number of bases trimmed off from 3' end of each read. \code{0} by default.}
  \item{minFragLength}{ numeric value giving the minimum fragment length. 50 by default.}
  \item{maxFragLength}{ numeric value giving the maximum fragment length. 600 by default.}
  \item{PE_orientation}{ character string giving the orientation of the two reads from the same pair. It has three possible values including \code{fr}, \code{ff} and \code{rf}. Letter \code{f} denotes the forward strand and letter \code{r} the reverse strand. \code{fr} by default (ie. the first read in the pair is on the forward strand and the second read on the reverse strand).}
  \item{nthreads}{ numeric value giving the number of threads used for mapping. \code{1} by default.}
  \item{readGroupID}{ a character string giving the read group ID. The specified string is added to the read group header field and also be added to each read in the mapping output. \code{NULL} by default.}
  \item{readGroup}{ a character string in the format of \code{tag:value}. This string will be added to the read group (RG) header in the mapping output. \code{NULL} by default.}
  \item{keepReadOrder}{ logical indicating if the order of reads in the BAM output is kept the same as that in the input file. (Reads from the same pair are always placed next to each other regardless of this option).}
  \item{sortReadsByCoordinates}{ logical indicating if reads will be sorted by their mapping locations in the mapping output. This option is applicable for BAM output only. A BAI index file will also be generated for each BAM file so the generated BAM files can be directly loaded into a genome browser such as IGB or IGV. \code{FALSE} by default.}
  \item{color2base}{ logical. If \code{TRUE}, color-space read bases will be converted to base-space bases in the mapping output. Note that the mapping itself will still be performed at color-space. \code{FALSE} by default.}
  \item{DP_GapOpenPenalty}{ a numeric value giving the penalty for opening a gap when using the Smith-Waterman dynamic programming algorithm to detect insertions and deletions. The Smith-Waterman algorithm is only applied for those reads which are found to contain insertions or deletions. \code{-1} by default.}
  \item{DP_GapExtPenalty}{ a numeric value giving the penalty for extending the gap, used by the Smith-Waterman algorithm. \code{0} by default.}
  \item{DP_MismatchPenalty}{ a numeric value giving the penalty for mismatches, used by the Smith-Waterman algorithm. \code{0} by default.}
  \item{DP_MatchScore}{ a numeric value giving the score for matches used by the Smith-Waterman algorithm. \code{2} by default.}
  \item{detectSV}{ logical indicating if structural variants (SVs) will be detected during read mapping. See below for more details.}
  \item{reportAllJunctions}{ logical indicating if all discovered junctions will be reported. This includes exon-exon junctions and also gene fusions. Presence of donor/receptor sites is not required when \code{reportAllJunctions} is \code{TRUE}. This option should only be used for RNA-seq data.}
  \item{useAnnotation}{ logical indicating if gene annotation information will be used in the mapping. \code{FALSE} by default.}
  \item{annot.inbuilt}{ a character string specifying an in-built annotation used for read summarization. It has four possible values including \code{"mm10"}, \code{"mm9"}, \code{"hg38"} and \code{"hg19"}, corresponding to the NCBI RefSeq annotations for genomes `mm10', `mm9', `hg38' and `hg19', respectively. \code{"mm10"} by default. See \code{\link{featureCounts}} function for more details on the in-built annotations.}
  \item{annot.ext}{ A character string giving name of a user-provided annotation file or a data frame including user-provided annotation data. If the annotation is in GTF format, it can only be provided as a file. If it is in SAF format, it can be provided as a file or a data frame. The provided annotation file can be a uncompressed or gzip compressed file. If an annotation is provided via \code{annot.ext}, the \code{annot.inbuilt} parameter will be ignored. See \code{\link{featureCounts}} function for more details about this parameter.}
  \item{isGTF}{ logical indicating if the annotation provided via the \code{annot.ext} argument is in GTF format or not. \code{FALSE} by default. This option is only applicable when \code{annot.ext} is not \code{NULL}.}
  \item{GTF.featureType}{ a character string denoting the type of features that will be extracted from a GTF annotation. \code{"exon"} by default. This argument is only applicable when \code{isGTF} is \code{TRUE}.}
  \item{GTF.attrType}{ a character string denoting the type of attributes in a GTF annotation that will be used to group features. \code{"gene_id"} by default. The grouped features are called meta-features. For instance, a feature can be an exon and several exons can be grouped into a gene (meta-feature). This argument is only applicable when \code{isGTF} is \code{TRUE}.}
  \item{chrAliases}{ a character string providing the name of a comma-delimited text file that includes aliases of chromosome names. This file should contain two columns. First column contains names of chromosomes included in the SAF or GTF annotation and second column contains corresponding names of chromosomes in the reference genome. No column headers should be provided. Also note that chromosome names are case sensitive. This file can be used to match chromosome names between the annotation and the reference genome.}
}
\details{
\code{align} and \code{subjunc} are R versions of programs in the Subread suite of programs (\url{http://subread.sourceforge.net}) that implement fast, accurate sequence read alignment based on the "seed-and-vote" mapping paradigm (Liao et al. 2013; Liao et al. 2019).
\code{align} is the R version of the command-line program \code{subread-align} and \code{subjunc} is the R version of the command-line program \code{subjunc}.

\code{align} is a general-purpose aligner that can be used to align both genomic DNA-seq reads and RNA-seq reads.
\code{subjunc} is designed for mapping RNA-seq reads only.
The main difference between \code{subjunc} and \code{align} for RNA-seq alignment is that \code{subjunc} reports discovered exon-exon junctions and it also performs full alignments for every read including exon-spanning reads.
\code{align} instead reports the largest mappable regions for each read and soft-clips the remainder of the read.
\code{subjunc} requires the presence of donor and receptor sites when calling exon-exon junctions.
It can detect up to four junction locations in a junction read.

\code{align} and \code{subjunc} can be run with either a full index or a gapped index (see \code{\link{buildindex}} for more details).
Maximum mapping speed is achieved when the full index is used.
With a full index built for human/mouse genome, \code{align} maps ~14 million reads per minute with 10 threads.
\code{subjunc} is slightly slower than \code{align}.
Both \code{align} and \code{subjunc} require 17.8GB of memory for mapping.
With a gapped index built for human/mouse genome, \code{align} maps ~9 million reads per minute with 10 threads.
\code{align} requires 8.2GB of memory for mapping and \code{subjunc} requires 8.8GB of memory.

The \code{detectSV} argument of \code{align} can be used to detect structural variants (SVs) in genomic DNA sequencing data.
The \code{reportAllJunctions} argument of \code{subjunc} function can be used for the detection of SVs in RNA-seq data, as well as the detection of non-canonical junctions that do not have the canonical donor/receptor sites.
For each library, breakpoints detected from SV events are saved to a file with suffix name `.breakpoints.txt', which includes chromosomal coordinates of SV breakpoints and numbers of supporting reads.
The BAM/SAM output includes extra fields to describe the complete alignments of breakpoint-containing reads.
For a breakpoint-containing read, mapping of its major sequence segment is described in the main fields of BAM/SAM output whereas mapping of its minor sequence segment, which does not map along with the major segment due to the presence of a breakpoint, is described in the extra fields including 'CC'(Chr), 'CP'(Position),'CG'(CIGAR) and 'CT'(strand).
Note that each breakpoint-containing read occupies only one row in the BAM/SAM output.
}
\value{
For each library, mapping results including both aligned and unaligned reads are saved to a file.
Indels discovered during mapping are saved to a VCF format file (\code{.indel.vcf}).
\code{subjunc} also outputs a BED format file that contains list of discovered exon-exon junctions (\code{.junction.bed}).
If \code{detectSV} or \code{reportAllJunctions} options are set to \code{TRUE}, structural variants discovered during read mapping will be reported and saved to a file (\code{.breakpoints.txt}), which includes chromosomal coordinates of breakpoints and number of supporting reads for each breakpoint.

Other than outputting data to files, \code{align} and \code{subjunc} also return a \code{data.frame} that includes mapping statistics for each library, such as total number of reads, number of mapped reads, number of uniquely mapped reads, number of indels and number of junctions (\code{subjunc} only).
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi (2019).
The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads.
\emph{Nucleic Acids Research}, 47(8):e47.
\url{http://www.ncbi.nlm.nih.gov/pubmed/30783653}

Yang Liao, Gordon K Smyth and Wei Shi (2013).
The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.
\emph{Nucleic Acids Research}, 41(10):e108.
\url{http://www.ncbi.nlm.nih.gov/pubmed/23558742}
}
\author{Wei Shi and Yang Liao}
%\note{}
\seealso{
\code{\link{buildindex}}
}
\examples{
# Build an index for the artificial sequence included in file 'reference.fa'.
ref <- system.file("extdata","reference.fa",package="Rsubread")

buildindex(basename="./reference_index",reference=ref)

# align a sample read dataset ('reads.txt') to the sample reference
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")

align.stat <- align(index = "./reference_index", readfile1 = reads,
output_file = "./Rsubread_alignment.BAM", phredOffset = 64)

align.stat
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
