\name{qualityScores}
\alias{qualityScores}
\title{Extract quality score data in a sequencing read dataset}
\description{Extract quality strings and convert them to Phred scores
}
\usage{
qualityScores(filename, input_format="gzFASTQ", offset=33, nreads=10000)
}
\arguments{
  \item{filename}{ character string giving the name of an input file containing sequence reads.}
  \item{input_format}{ character string specifying format of the input file. \code{gzFASTQ} (gzipped FASTQ) by default. Acceptable formats include \code{gzFASTQ}, \code{FASTQ}, \code{SAM} and \code{BAM}. Character string is case insensitive.}
  \item{offset}{ numeric value giving the offset added to the Phred scores, \code{33} by default.}
  \item{nreads}{ numeric value giving the number of reads from which quality scores are extracted. \code{10000} by default.}
}
\details{
Quality scores of read bases are represented by ASCII characters in next-gen sequencing data.
This function extracts the quality characters from each base in each read and then converts them to Phred scores using the provided offset value (\code{offset}).

If the total number of reads in a dataset is \code{n}, then every \code{n}/\code{nreads} read is extracted from the input data.
}
\value{ 
A data matrix containing Phred scores for read bases. Rows in the matrix are reads and columns are base positions in each read.
}
%\references{}
\author{Wei Shi, Yang Liao and Zhiyin Dai}
%\note{}
%\seealso{}
\examples{
library(Rsubread)
reads <- system.file("extdata","reads.txt",package="Rsubread")
x <- qualityScores(filename=reads,input_format="FASTQ",offset=64,nreads=1000)
x[1:10,1:10]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

