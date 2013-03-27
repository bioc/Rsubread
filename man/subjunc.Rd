\name{subjunc}
\alias{subjunc}
\title{Discovering exon-exon junctions from RNA-seq data}
\description{Discover exon-exon junction locations from RNA-seq data using seed-and-vote paradigm}
\usage{
subjunc(index,samfile,output_file,nsubreads=14,paired_end=FALSE,nthreads=1,indels=5,min_distance=50,max_distance=600,PE_orientation="fr")
}
\arguments{
  \item{index}{ character string giving the basename of index file. Index files should be located in the current directory.}
  \item{samfile}{ character string giving the name of SAM file generated from read alignment (eg. output from \code{align} function).}
  \item{output_file}{ character string giving the name of output file which includes the mapping results.}
  \item{nsubreads}{ numeric value giving the number of subreads extracted from each read. \code{14} by default.}
  \item{paired_end}{ logical, indicating whether the data is paired ended or single ended. \code{FALSE} by default.}
  \item{nthreads}{ numeric value giving the number of threads used for mapping. 1 by default.} 
  \item{indels}{ numeric value giving the number of insertions/deletions allowed during the mapping. 5 by default.}
  \item{min_distance}{ numeric value giving the minimal distance between the pair of reads. 50 by default.}
  \item{max_distance}{ numeric value giving the maximal distance between the pair of reads. 600 by default.}
  \item{PE_orientation}{ character string giving the orientation of the two reads in a pair. "fr" by default, which means the first read is on the forward strand and the second read is on the reverse strand.}
  }
\details{
This function takes as input the SAM format file generated from read alignment (eg. output from \code{align} function), and then tries to discover exon-exon junctions.
Donor and receptor sites are being examined when detecting junctions.
This function can detect junction locations located at any position of the read.
It can detect spliced exons which are as far as 500kb.
It supports INDEL detection during the junction detection as well.
\code{align} function should be run first before calling this function.

This function calls the underlying C program to perform the exon junction discovery.
The C program can be run on its own as well, which is included in the subread software package (http://subread.sourceforge.net/)

}
\value{ 
This function outputs two files: one is a SAM format file which includes the mapping results for all the reads (including both junction reads and exonic reads), and the other is a BED format file which inlcudes pairs of splicing points (exon-exon junctions) found from each junction read, number of supporting junction reads for each exon-exon junction and so on.
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi (2013), "The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote",  Nucleic Acids Research, Volume 41, Accepted on 8 March, In Press.
}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
