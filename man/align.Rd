\name{align}
\alias{align}
\title{Align next-gen sequencing reads to reference genome}
\description{This function performs read mapping using the seed-and-vote method. It can be used for mapping both gDNA and RNA sequencing data. It supports indel detection.}
\usage{
align(index,readfile1,readfile2=NULL,output_file,nsubreads=10,TH1=3,TH2=1,nthreads=1,indels=5,phredOffset=33,markJunctionReads=TRUE,tieBreakQS=FALSE,tieBreakHamming=FALSE,unique=FALSE,min_distance=50,max_distance=600,PE_orientation="fr",DP_GapOpenPenalty=-2,DP_GapExtPenalty=0,DP_MismatchPenalty=0,DP_MatchScore=2)
}
\arguments{
  \item{index}{ character string giving the basename of index file. Index files should be located in the current directory.}
  \item{readfile1}{ character string giving the name of file which includes sequencing reads. This will be the name of first file of paired-end data are provided. The read file should be in FASTQ or FASTA format.}
  \item{readfile2}{ character string giving the name of second file when paired-end read data are provided. \code{NULL} by default.}
  \item{output_file}{ character string giving the name of output file which includes the mapping results.}
  \item{nsubreads}{ numeric value giving the number of subreads extracted from each read. 10 by default.} 
  \item{TH1}{ numeric value giving the consensus threshold for reporting a hit. This threshold will be applied to the first read if paired-end read data are provided. 3 by default.}
  \item{TH2}{ numeric value giving the consensus threhold for the second read in a pair. 1 by default.}
  \item{nthreads}{ numeric value giving the number of threads used for mapping. 1 by default.} 
  \item{indels}{ numeric value giving the number of insertions/deletions allowed during the mapping. 5 by default.}
  \item{phredOffset}{ numeric value added to the phred score of each base in each read. \code{33} by default.}
  \item{markJunctionReads}{ a logical value indicating whether the unmapped bases of discovered junction reads should be marked by "S" operations in their CIGAR strings in the SAM output file. The unmapped bases should originate from an exon which is different from the one the corresponding read was mapped to. 'Subjunc' program can then be called to discover exon junction locations. This options is only applicable to RNAseq data. \code{TRUE} by default.}
  \item{tieBreakQS}{logical. If \code{TRUE}, use quality scores to break ties when more than one best locations were found. \code{FALSE} by default. }
  \item{tieBreakHamming}{ logical. If \code{TRUE}, use Hamming distance to break ties when more than one best locations were found (to use this option, the index must be built with reference sequence included).}
  \item{unique}{ logical. If \code{TRUE}, only uniquely mapped reads will be reported (reads mapped to multiple locations in the reference genome will not be reported). This option can be used together with option \code{tieBreakQS} or \code{tieBreakHamming}.}
  \item{min_distance}{ numeric value giving the minimal distance between the pair of reads. 50 by default.}
  \item{max_distance}{ numeric value giving the maximal distance between the pair of reads. 600 by default.}
  \item{PE_orientation}{ character string giving the orientation of the two reads in a pair. "fr" by default, which means the first read is on the forward strand and the second read is on the reverse strand.}
  \item{DP_GapOpenPenalty}{ a numeric value giving the penalty for opening a gap when using the Smith-Waterman dynamic programming algorithm to detect insertions and deletions. The Smith-Waterman algorithm is only applied for those reads which are found to contain insertions or deletions. \code{-2} by default.}
  \item{DP_GapExtPenalty}{ a numeric value giving the penalty for extending the gap, used by the Smith-Waterman algorithm. \code{0} by default.}
  \item{DP_MismatchPenalty}{ a numeric value giving the penalty for mismatches, used by the Smith-Waterman algorithm. \code{0} by default.}
  \item{DP_MatchScore}{ a numeric value giving the score for matches used by the Smith-Waterman algorithm. \code{2} by default.}
  }
\details{
This function implements a novel mapping strategy which uses a set of 16bp substrings (called subreads) extracted from each read to map them to the reference genome.
Different from "seed-and-extend" mapping strategy, this new strategy does not have an extension step therefore it is a lot faster than the competing aligners.
Due to the short length of the selected subsreads (16bp long), this strategy has a much high sensitivity than other aligners (seed length is usually around 30bp), i.e. it can align a lot more reads than competing aligners.
Our evaluation results (using both simulation dataset and real dataset) showed that the accuracy of the new strategy is comparable to or slightly better than other aligners.

Two key parameters used by this new strategy are the number of subreads selected \code{nsubreads} and the consensus threshold for determining mapping locations \code{TH1} (also \code{TH2} for paired-end read data).
We recommend using the default setting of these parameters to map reads of around 100bp long.
However, users can choose to use more subreads when mapping longer reads.
We recommend to set the value of consensus threshold to be 30 percent of the number of subreads used.

The C implemetation of this strategy can be found at http://subread.sourceforge.net/.
This R function calls the corresponding C function to perform the alignment.
Therefore, it has the mapping speed as the C program.

\code{\link{buildindex}} function should be called if the index has not been built for the reference genome.
The index can be re-used once it has been built.

If paired-end read data is provided, file \code{readfile1} will assumed to contain the first read from the read pair and \code{readfile2} the second read. 
}
\value{ 
No value is produced but a SAM format file is written to the current working directory.
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi, "The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote", Manuscript under review.
}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
\examples{
# build an index using the sample reference sequence provided in the package
# and save it to the current directory
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)

# align the sample read data provided in this packge to the sample reference
# and save the mapping results to the current directory
reads <- system.file("extdata","reads.txt",package="Rsubread")
align(index="./reference_index",readfile1=reads,output_file="./Rsubread_alignment.SAM")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
