\name{subjunc}
\alias{subjunc}
\title{Discovering exon-exon junctions from using RNA-seq data}
\description{Use the seed-and-vote paradigm to accurately detect exon-exon junctions from RNA-seq data}
\usage{
subjunc(index,samfile,output_file,nsubreads=14,paired_end=FALSE,nthreads=1,indels=5,minFragLength=50,maxFragLength=600,PE_orientation="fr")
}
\arguments{
  \item{index}{ character string giving the basename of index file. Index files should be located in the current directory.}
  \item{samfile}{ character string giving the name of SAM file generated from read alignment (eg. output from \code{align} function).}
  \item{output_file}{ character string giving the name of output file which includes the mapping results.}
  \item{nsubreads}{ numeric value giving the number of subreads extracted from each read. \code{14} by default.}
  \item{paired_end}{ logical, indicating whether the data is paired ended or single ended. \code{FALSE} by default.}
  \item{nthreads}{ numeric value giving the number of threads used for mapping. 1 by default.} 
  \item{indels}{ numeric value giving the number of insertions/deletions allowed during the mapping. 5 by default.}
  \item{minFragLength}{ numeric value giving the minimum fragment length. 50 by default.}
  \item{maxFragLength}{ numeric value giving the maximum fragment length. 600 by default.}
  \item{PE_orientation}{ character string giving the orientation of the two reads in a pair. "fr" by default, which means the first read is on the forward strand and the second read is on the reverse strand.}
  }
\details{
The RNA-seq technology provides a unique opportunity to identify the alternative splicing events that occur during the gene transcription process.
This function makes use of a powerful read mapping paradigm called ``seed-and-vote" to accurately detect the exon-exon junctions (Liao et al. 2013).

This function takes as input the mapping results (in SAM format) from a read aligner, preferably the Subread aligner included in this package (\code{\link{align}}).
It extracts a number of subreads (16mers) from each read, maps them to the reference genome and then identifies the two best mapping locations for each read.
A verification step is then applied to determine whether each read should be mapped as an exon-spanning read or not and to report the discovered exon-exon junctions.
The donor and receptor sites, which are splicing signals, are required to be present when calling exon-exon junctions.
The indels identified in the mapped reads included in the provided SAM file will be kept in the output of this function.
This function also reports the indels discovered in the exon-spanning reads.

It has the capacity to detect junction locations that are present at any position of the reads.
The pair of exons that are spliced together are allowed to be up to 500kb apart in the genome.

This function can run on both single-ended and paired-end data.

It calls the \emph{Subjunc} program (written in C programming language) included in the Subread package (http://subread.sourceforge.net/) for this analysis.
}
\value{ 
This function outputs two files: one is a SAM format file which includes the mapping results for all the reads (including both junction reads and exonic reads), and the other is a BED format file which inlcudes pairs of splicing points (exon-exon junctions) found from each junction read, number of supporting junction reads for each exon-exon junction and so on.
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.  Nucleic Acids Research, 41(10):e108, 2013.
}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
