\name{align}
\alias{align}
\alias{subjunc}
\title{Read mapping for genomic DNA-seq and RNA-seq data via seed-and-vote (Subread and Subjunc)}
\description{Subread and Subjunc are fast and accurate read aligners. They adopt a new mapping read mapping paradigm called ``seed-and-vote". Seed-and-vote enables highly scalable, efficient and accurate read mapping to be performed.}
\usage{
align(index,readfile1,readfile2=NULL,input_format="FASTQ",output_file,output_format="SAM",nsubreads=10,TH1=3,TH2=1,nthreads=1,indels=5,phredOffset=33,tieBreakQS=FALSE,tieBreakHamming=FALSE,unique=FALSE,nBestLocations=1,minFragLength=50,maxFragLength=600,PE_orientation="fr",nTrim5=0,nTrim3=0,readGroupID=NULL,readGroup=NULL,color2base=FALSE,DP_GapOpenPenalty=-2,DP_GapExtPenalty=0,DP_MismatchPenalty=0,DP_MatchScore=2)
subjunc(index,readfile1,readfile2=NULL,input_format="FASTQ",output_file,output_format="SAM",nsubreads=14,TH1=1,TH2=1,nthreads=1,indels=5,phredOffset=33,tieBreakQS=FALSE,tieBreakHamming=FALSE,unique=FALSE,minFragLength=50,maxFragLength=600,PE_orientation="fr",nTrim5=0,nTrim3=0,readGroupID=NULL,readGroup=NULL,color2base=FALSE,DNAseq=FALSE,reportAllJunctions=FALSE)
}
\arguments{
  \item{index}{ character string giving the basename of index file. Index files should be located in the current directory.}
  \item{readfile1}{ character string giving the name of file which includes sequencing reads. This will be the name of first file of paired-end data are provided. The read file should be in FASTQ or FASTA format.}
  \item{readfile2}{ character string giving the name of second file when paired-end read data are provided. \code{NULL} by default.}
  \item{input_format}{ character string specifying format of the read input file(s). \code{FASTQ} by default. Acceptable formats include \code{FASTQ}, \code{FASTA}, \code{SAM} and \code{BAM}. The character string is case insensitive.}
  \item{output_file}{ character string giving the name of output file which includes the mapping results.}
  \item{output_format}{ character string specifying format of the output file. \code{SAM} by default. Acceptable formats include \code{SAM} and \code{BAM}.}
  \item{nsubreads}{ numeric value giving the number of subreads extracted from each read.} 
  \item{TH1}{ numeric value giving the consensus threshold for reporting a hit. This threshold will be applied to the first read if paired-end read data are provided.}
  \item{TH2}{ numeric value giving the consensus threhold for the second read in a pair.}
  \item{nthreads}{ numeric value giving the number of threads used for mapping. \code{1} by default.} 
  \item{indels}{ numeric value giving the number of insertions/deletions allowed during the mapping. \code{5} by default.}
  \item{phredOffset}{ numeric value added to the phred score of each base in each read. \code{33} by default.}
  \item{tieBreakQS}{logical. If \code{TRUE}, use quality scores to break ties when more than one best locations were found. \code{FALSE} by default. }
  \item{tieBreakHamming}{ logical. If \code{TRUE}, use Hamming distance to break ties when more than one best locations were found (to use this option, the index must be built with reference sequence included).}
  \item{unique}{ logical. If \code{TRUE}, only uniquely mapped reads will be reported (reads mapped to multiple locations in the reference genome will not be reported). This option can be used together with option \code{tieBreakQS} or \code{tieBreakHamming}.}
  \item{nBestLocations}{ numeric value giving the maximal number of equally-best mapping locations allowed to be reported for the read. 1 by default. The value of \code{nBestLocations} has to be within the range of 1 to 16. The number of equally-best locations reported for a read will be less than or equal to this value. For example, if a read has two equally-best mapping locations, but the \code{nBestLocations} set to be 5, then only two locations will be reported for this read. Note that the \code{-u} takes precedence over \code{nBestLocations}, ie. no mapping locations will be reported for multi-mapping reads when \code{unique} is set to \code{TRUE}.}
  \item{minFragLength}{ numeric value giving the minimum fragment length. 50 by default.}
  \item{maxFragLength}{ numeric value giving the maximum fragment length. 600 by default.}
  \item{PE_orientation}{ character string giving the orientation of the two reads from the same pair. It has three possible values including \code{fr}, \code{ff} and \code{rf}. Letter \code{f} denotes the forward strand and letter \code{r} the reverse strand. \code{fr} by default (ie. the first read in the pair is on the forward strand and the second read on the reverse strand).}
  \item{nTrim5}{ numeric value giving the number of bases trimmed off from 5' end of each read. \code{0} by default.}
  \item{nTrim3}{ numeric value giving the number of bases trimmed off from 3' end of each read. \code{0} by default.}
  \item{readGroupID}{ a character string giving the read group ID. The specified string is added to the read group header field and also be added to each read in the mapping output. \code{NULL} by default.}
  \item{readGroup}{ a character string in the format of \code{tag:value}. This string will be added to the read group (RG) header in the mapping output. \code{NULL} by default.}
  \item{color2base}{ logical. If \code{TRUE}, color-space read bases will be converted to base-space bases in the mapping output. Note that the mapping itself will still be performed at color-space. \code{FALSE} by default.}
  \item{DP_GapOpenPenalty}{ a numeric value giving the penalty for opening a gap when using the Smith-Waterman dynamic programming algorithm to detect insertions and deletions. The Smith-Waterman algorithm is only applied for those reads which are found to contain insertions or deletions. \code{-2} by default.}
  \item{DP_GapExtPenalty}{ a numeric value giving the penalty for extending the gap, used by the Smith-Waterman algorithm. \code{0} by default.}
  \item{DP_MismatchPenalty}{ a numeric value giving the penalty for mismatches, used by the Smith-Waterman algorithm. \code{0} by default.}
  \item{DP_MatchScore}{ a numeric value giving the score for matches used by the Smith-Waterman algorithm. \code{2} by default.}
  \item{DNAseq}{ logical, indicating if the input read data are genomic DNA sequencing data. If \code{TRUE}, \code{subjunc} function will ignore the splicing signals (donor/receptor sites) when searching for junctions. Junctions that occur within the same chromosome or across different chromosomes (chimerism) will all be reported for DNA seqence data.}
  \item{reportAllJunctions}{ logical. This argument should be used for RNA-seq data. If \code{TRUE}, \code{subjunc} function will report those junctions that do not have the required donor/receptor sites (GT/AG), or cross different chromosomes or are located on different strands within the same chromosome, in addition to the canonical exon-exon junctions.}
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

\code{Subjunc} requires donor and receptor sites to be existent when detecting exon-exon junctions.
It can detect up to four junction locations in a single exon-spanning reads.
}
\value{ 
No value is produced but a SAM format file is written to the current working directory. For \code{Subjunc}, a BED file including discovered exon-exon junction will also be written to the current working directory.
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.  Nucleic Acids Research, 41(10):e108, 2013.
}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
\examples{
# build an index using a sample reference sequence ('reference.fa') that includes a chromosome called 'chr_dummy' and 90,000 bases
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)

# align a sample read dataset ('reads.txt') to the sample reference
reads <- system.file("extdata","reads.txt",package="Rsubread")
align(index="./reference_index",readfile1=reads,output_file="./Rsubread_alignment.SAM")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
