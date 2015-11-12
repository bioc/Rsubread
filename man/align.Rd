\name{align}
\alias{align}
\alias{subjunc}
\title{Read mapping for genomic DNA-seq and RNA-seq data via seed-and-vote (Subread and Subjunc)}
\description{Subread and Subjunc perform local and global alignments respectively. The seed-and-vote paradigm enables efficient and accurate alignments to be carried out.}
\usage{
align(index,readfile1,readfile2=NULL,type=0,input_format="gzFASTQ",output_format="BAM",
output_file=paste(readfile1,"subread",output_format,sep="."),nsubreads=10,
TH1=3,TH2=1,maxMismatches=3,nthreads=1,indels=5,complexIndels=FALSE,phredOffset=33,
unique=TRUE,nBestLocations=1,minFragLength=50,maxFragLength=600,PE_orientation="fr",
nTrim5=0,nTrim3=0,readGroupID=NULL,readGroup=NULL,color2base=FALSE,
DP_GapOpenPenalty=-1,DP_GapExtPenalty=0,DP_MismatchPenalty=0,DP_MatchScore=2,
detectSV=FALSE)

subjunc(index,readfile1,readfile2=NULL,input_format="gzFASTQ",output_format="BAM",
output_file=paste(readfile1,"subjunc",output_format,sep="."),nsubreads=14,
TH1=1,TH2=1,maxMismatches=3,nthreads=1,indels=5,complexIndels=FALSE,phredOffset=33,
unique=TRUE,nBestLocations=1,minFragLength=50,maxFragLength=600,PE_orientation="fr",
nTrim5=0,nTrim3=0,readGroupID=NULL,readGroup=NULL,color2base=FALSE,
reportAllJunctions=FALSE)
}
\arguments{
  \item{index}{ character string giving the basename of index file. Index files should be located in the current directory.}
  \item{readfile1}{ a character vector including names of files that include sequence reads to be aligned. For paired-end reads, this gives the list of files including first reads in each library. File format is FASTQ/FASTA by default. See \code{input_format} option for more supported formats.}
  \item{readfile2}{ a character vector giving names of files that include second reads in paired-end read data. Files included in \code{readfile2} should be in the same order as their mate files included in \code{readfile1} .\code{NULL} by default.}
  \item{type}{ an integer giving the type of sequencing data. Possible values include \code{0} (RNA-seq data) and \code{1} (genomic DNA-seq data such as WGS, WES, ChIP-seq data etc.).}
  \item{input_format}{ character string specifying format of read input files. \code{gzFASTQ} by default (this also includes FASTQ, FASTA and gzipped FASTA formats). Other supported formats include \code{SAM} and \code{BAM}. Character values are case insensitive.}
  \item{output_format}{ character string specifying format of output file. \code{BAM} by default. Acceptable formats include \code{SAM} and \code{BAM}.}
  \item{output_file}{ a character vector specifying names of output files. By default, names of output files are set as the file names provided in \code{readfile1} added with an suffix string.}
  \item{nsubreads}{ numeric value giving the number of subreads extracted from each read.} 
  \item{TH1}{ numeric value giving the consensus threshold for reporting a hit. This is the threshold for the first reads if paired-end read data are provided.}
  \item{TH2}{ numeric value giving the consensus threhold for the second reads in paired-end data.}
  \item{maxMismatches}{ numeric value giving the maximum number of mis-matched bases allowed in the alignment. \code{3} by default. Mis-matches found in soft-clipped bases are not counted.}
  \item{nthreads}{ numeric value giving the number of threads used for mapping. \code{1} by default.}
  \item{indels}{ numeric value giving the maximum number of insertions/deletions allowed during the mapping. \code{5} by default.}
  \item{complexIndels}{ logical indicating if complex indels will be detected. If \code{TRUE}, the program will try to detect multiple short indels that occurs concurrently in a small genomic region (indels could be as close as 1bp apart).}
  \item{phredOffset}{ numeric value added to base-calling Phred scores to make quality scores (represented as ASCII letters). Possible values include \code{33} and \code{64}. By default, \code{33} is used.}
  \item{unique}{ logical indicating if uniquely mapped reads should be reported only. \code{TRUE} by default. Uniquely mapped reads must have one (1) mapping location that has less mis-matched bases than other candidate locations.}
  \item{nBestLocations}{ numeric value giving the maximal number of equally-best mapping locations allowed to be reported for the read. 1 by default. The allowed value is between 1 to 16 (inclusive). `NH' tag is used to indicate how many alignments are reported for the read and `HI' tag is used for numbering the alignments reported for the same read, in the output. Note that the \code{unique} argument takes precedence over \code{nBestLocations} argument.}
  \item{minFragLength}{ numeric value giving the minimum fragment length. 50 by default.}
  \item{maxFragLength}{ numeric value giving the maximum fragment length. 600 by default.}
  \item{PE_orientation}{ character string giving the orientation of the two reads from the same pair. It has three possible values including \code{fr}, \code{ff} and \code{rf}. Letter \code{f} denotes the forward strand and letter \code{r} the reverse strand. \code{fr} by default (ie. the first read in the pair is on the forward strand and the second read on the reverse strand).}
  \item{nTrim5}{ numeric value giving the number of bases trimmed off from 5' end of each read. \code{0} by default.}
  \item{nTrim3}{ numeric value giving the number of bases trimmed off from 3' end of each read. \code{0} by default.}
  \item{readGroupID}{ a character string giving the read group ID. The specified string is added to the read group header field and also be added to each read in the mapping output. \code{NULL} by default.}
  \item{readGroup}{ a character string in the format of \code{tag:value}. This string will be added to the read group (RG) header in the mapping output. \code{NULL} by default.}
  \item{color2base}{ logical. If \code{TRUE}, color-space read bases will be converted to base-space bases in the mapping output. Note that the mapping itself will still be performed at color-space. \code{FALSE} by default.}
  \item{DP_GapOpenPenalty}{ a numeric value giving the penalty for opening a gap when using the Smith-Waterman dynamic programming algorithm to detect insertions and deletions. The Smith-Waterman algorithm is only applied for those reads which are found to contain insertions or deletions. \code{-1} by default.}
  \item{DP_GapExtPenalty}{ a numeric value giving the penalty for extending the gap, used by the Smith-Waterman algorithm. \code{0} by default.}
  \item{DP_MismatchPenalty}{ a numeric value giving the penalty for mismatches, used by the Smith-Waterman algorithm. \code{0} by default.}
  \item{DP_MatchScore}{ a numeric value giving the score for matches used by the Smith-Waterman algorithm. \code{2} by default.}
  \item{detectSV}{ logical indicating if structural variants (SVs) will be detected during read mapping. See below for more details.}
  \item{reportAllJunctions}{ logical indicating if all possible junctions and structural variants will be reported. Presence of donor/receptor sites is not required for junction calling. This argument should be used for RNA-seq data.}
  }
\details{
The \code{align} function implements the \code{Subread} aligner (Liao et al., 2013) that uses a new mapping paradigm called ``seed-and-vote".
\code{Subread} is general-purpose aligner that can be used to align both genomic DNA-seq reads and RNA-seq reads.

\code{Subjunc} is designed for mapping RNA-seq reads.
The major difference between \code{Subjunc} and \code{Subread} is that \code{Subjunc} reports discovered exon-exon junctions and it also performs full alignments for every read including exon-spanning reads.
\code{Subread} uses the largest mappable regions in the reads to find their mapping locations.
The seed-and-vote paradigm has been found to be not only more accurate than the conventional seed-and-extend (adopted by most aligners) in read mapping, but it is a lot more efficient (Liao et al., 2013).

Both \code{Subread} and \code{Subjunc} can be used to align reads generated from major sequencing platforms including Illumina GA/HiSeq, ABI SOLiD, Roche 454 and Ion Torrent sequencers.
Note that to map color-space reads (e.g. SOLiD reads), a color-space index should be built for the reference genome (see \code{\link{buildindex}} for details).

\code{Subread} and \code{Subjunc} have adjustable memory usage.
See \code{\link{buildindex}} function for more details.

Mapping performance is largely determined by the number of subreads extracted from each read \code{nsubreads} and the consensus threshold \code{TH1} (also \code{TH2} for paired-end read data).
Default settings are recommended for most of the read mapping tasks.

\code{Subjunc} requires donor/receptor sites to be present when detecting exon-exon junctions.
It can detect up to four junction locations in each exon-spanning read.

\code{detectSV} option should be used for SV detection in genomic DNA sequencing data.
For RNA-seq data, users may use \code{subjunc} with the \code{reportAllJunctions} option to detect SVs (and also junctions).
For each library, breakpoints detected from SV events will be saved to a file with suffix name `.breakpoints.txt', which includes chromosomal coordinates of SV breakpoints and numbers of supporting reads.
The BAM/SAM output includes extra fields to describe the complete alignments of breakpoint-containing reads.
For a breakpoint-containing read, mapping of its major sequence segment is described in the main fields of BAM/SAM ouptut whereas mapping of its minor sequence segment, which does not map along with the major segment due to the presence of a breakpoint, is described in the extra fields including `CC'(Chr), `CP'(Position),`CG'(CIGAR) and `CT'(strand).
Note that each breakpoint-containing read occupies only one row in BAM/SAM output. 
}
\value{ 
No value is produced but SAM or BAM format files are written to the current working directory. For \code{Subjunc}, BED files including discovered exon-exon junctions are also written to the current working directory.
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.  Nucleic Acids Research, 41(10):e108, 2013.
}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
\examples{
# Build an index for the artificial sequence included in file 'reference.fa'.
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)

# align a sample read dataset ('reads.txt') to the sample reference
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")
align(index="./reference_index",readfile1=reads,output_file="./Rsubread_alignment.BAM",phredOffset=64)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
