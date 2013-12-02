\name{align}
\alias{align}
\alias{subjunc}
\title{Read mapping for genomic DNA-seq and RNA-seq data via seed-and-vote (Subread and Subjunc)}
\description{Subread and Subjunc perform local and global alignments respectively. The seed-and-vote paradigm enables efficient and accurate alignments to be carried out.}
\usage{
align(index,readfile1,readfile2=NULL,input_format="FASTQ",output_format="SAM",
output_file=paste(readfile1,"subread",output_format,sep="."),nsubreads=10,TH1=3,TH2=1,
nthreads=1,indels=5,phredOffset=33,tieBreakQS=FALSE,tieBreakHamming=TRUE,unique=TRUE,nBestLocations=1,
minFragLength=50,maxFragLength=600,PE_orientation="fr",nTrim5=0,nTrim3=0,readGroupID=NULL,readGroup=NULL,color2base=FALSE,
DP_GapOpenPenalty=-1,DP_GapExtPenalty=0,DP_MismatchPenalty=0,DP_MatchScore=2,reportFusions=FALSE)

subjunc(index,readfile1,readfile2=NULL,input_format="FASTQ",output_format="SAM",
output_file=paste(readfile1,"subjunc",output_format,sep="."),nsubreads=14,TH1=1,TH2=1,
nthreads=1,indels=5,phredOffset=33,tieBreakQS=FALSE,tieBreakHamming=TRUE,unique=TRUE,
minFragLength=50,maxFragLength=600,PE_orientation="fr",nTrim5=0,nTrim3=0,readGroupID=NULL,readGroup=NULL,color2base=FALSE,
DNAseq=FALSE,reportAllJunctions=FALSE)
}
\arguments{
  \item{index}{ character string giving the basename of index file. Index files should be located in the current directory.}
  \item{readfile1}{ a character vector including names of files that include sequence reads to be aligned. For paired-end reads, this gives the list of files including first reads in each library. File format is FASTQ/FASTA by default. See \code{input_format} option for more supported formats.}
  \item{readfile2}{ a character vector giving names of files that include second reads in paired-end read data. Files included in \code{readfile2} should be in the same order as their mate files included in \code{readfile1} .\code{NULL} by default.}
  \item{input_format}{ character string specifying format of the read input file(s). \code{FASTQ} by default (this includes FASTA format as well). Acceptable formats include \code{FASTQ} (including FASTA), \code{gzFASTQ} (gzipped FASTQ or FASTA), \code{SAM} and \code{BAM}. The character string is case insensitive.}
  \item{output_format}{ character string specifying format of the output file. \code{SAM} by default. Acceptable formats include \code{SAM} and \code{BAM}.}
  \item{output_file}{ a character vector specifying names of output files. By default, names of output files are set as the file names provided in \code{readfile1} added with an suffix string.}
  \item{nsubreads}{ numeric value giving the number of subreads extracted from each read.} 
  \item{TH1}{ numeric value giving the consensus threshold for reporting a hit. This is the threshold for the first reads if paired-end read data are provided.}
  \item{TH2}{ numeric value giving the consensus threhold for the second reads in paired-end data.}
  \item{nthreads}{ numeric value giving the number of threads used for mapping. \code{1} by default.} 
  \item{indels}{ numeric value giving the maximum number of insertions/deletions allowed during the mapping. \code{5} by default.}
  \item{phredOffset}{ numeric value added to the Phred quality scores of read bases by the sequencer. \code{33} by default.}
  \item{tieBreakQS}{logical indicating if the mapping quality score should be to break the tie when more than one best location was found for a read. \code{FALSE} by default. Note that the mapping quality score used for tie breaking was calcuated only from the perfectly matched subreads (16mers) extracted from the read, whereas the mapping quality scores included in the mapping output for mapped reads were calcuated from using all the bases in the read. Also note that there may still be more than one best mapping location found after tie breaking using this option.}
  \item{tieBreakHamming}{ logical indicating if the Hamming distance should be used to break the tie when more than one best mapping location was found for a read. \code{TRUE} by default. The distance between the mapped read and the target region in the reference is calcuated using all the bases included in the read. Note that there may still be more than one best mapping location found after tie breaking using this option.}
  \item{unique}{ logical indicating if uniquely mapped reads should be reported only. \code{TRUE} by default. It is recommended that this option is used together with \code{tieBreakHamming} or \code{tieBreakQS}.}
  \item{nBestLocations}{ numeric value giving the maximal number of equally-best mapping locations allowed to be reported for the read. 1 by default. The allowed value is between 1 to 16 (inclusive). Note that the \code{unique} argument takes precedence over \code{nBestLocations} argument.}
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
  \item{reportFusions}{ logical indicating if genomic fusion events such as chimeras should be reported. If \code{TRUE}, \code{align} function will detect and report such events. This option should only be applied for genomic DNA sequencing data. For RNA-seq data, users should use \code{subjunc} with the \code{reportAllJunctions} option for detection of fusions. Discovered fusions will be saved to a file (*.fusions.txt). Detailed mapping results for fusion reads will also be saved to the SAM/BAM output file. Secondary alignments of fusion reads will be saved to the following optional fields: CC(Chr), CP(Position), CG(CIGAR) and CT(strand). Note that each fusion read occupies only one row in the SAM/BAM output file.}
  \item{DNAseq}{ logical indicating if the input read data are genomic DNA sequencing data. If \code{TRUE}, \code{subjunc} function will ignore the splicing signals (donor/receptor sites) when searching for junctions. Junctions that occur within the same chromosome or across different chromosomes (chimerism) will all be reported for DNA seqence data.}
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

Multiple input files (libraries) can be provided to \code{Subread} or \code{Subjunc} at the same time.

\code{Subjunc} requires donor/receptor sites to be present when detecting exon-exon junctions.
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
