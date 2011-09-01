\name{subjunc}
\alias{subjunc}
\title{Finding exon junctions}
\description{This function uses RNA-seq reads to discover exon junction locations}
\usage{
subjunc(index,readfile1,readfile2=NULL,output_file,nthreads=1,indels=0,min_distance=50,max_distance=600,PE_orientation="fr")
}
\arguments{
  \item{index}{ character string giving the basename of index file. Index files should be located in the current directory.}
  \item{readfile1}{ character string giving the name of file which includes sequencing reads. This will be the name of first file of paired-end data are provided. The read file should be in FASTQ or FASTA format.}
  \item{readfile2}{ character string giving the name of second file when paired-end read data are provided. \code{NULL} by default.}
  \item{output_file}{ character string giving the name of output file which includes the mapping results.}
  \item{nthreads}{ numeric value giving the number of threads used for mapping. 1 by default.} 
  \item{indels}{ numeric value giving the number of insertions/deletions allowed during the mapping. 0 by default.}
  \item{min_distance}{ numeric value giving the minimal distance between the pair of reads. 50 by default.}
  \item{max_distance}{ numeric value giving the maximal distance between the pair of reads. 600 by default.}
  \item{PE_orientation}{ character string giving the orientation of the two reads in a pair. "fr" by default, which means the first read is on the forward strand and the second read is on the reverse strand.}
  }
\details{
This function uses a number of overlapping 16 mers extracted from each read to find junction reads and uses the mapping locations of these 16 mers to determine the exon junction locations in the reference genome.
Exon junction locations with known donor and receptor sites are reported.
This function calls the underlying C program to perform the exon junction discovery.
The C program can be run on its own as well, which is included in the subread software package (http://subread.sourceforge.net/)

To use this function, the index has to built with \code{-b} option enabled (\code{\link{buildindex}}).
}
\value{ 
This function ouputs two files: one is a BED format file which inlcudes the discovered junction locations and numbers of supporting reads, and the other is a SAM format file which includes the mapping result for each read.
}
\references{
}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
