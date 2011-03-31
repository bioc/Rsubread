\name{qualityScores}
\alias{qualityScores}
\title{Extract quality score infromation from a sequencing read dataset}
\description{Extract quality scores and convert them to ASCII code
}
\usage{
qualityScores(filename, offset=64, nreads=10000)
}
\arguments{
  \item{filename}{ character string giving the name of input FASTQ file.}
  \item{offset}{ numeric value giving the offset added to the original quality score, 64 by default.}
  \item{nreads}{ numeric value giving the number of reads from which quality scores are extracted}
}
\details{
Quality scores are given in the form of characters in datasets which contain sequencing reads.
This function extracts the quality scores and then convert them to the ASCII codes which encode these characters.
These ASCII codes are then subtracted by the \code{offset} to obtain the original quality scores.

If the total number of reads is \code{n}, then every \code{n}/\code{nreads} read will be used for quality score retrieval.
}
\value{ 
A data matrix containing the quality scores with rows being reads and columns being base positions in the read.
}
%\references{}
\author{Zhiyin Dai and Wei Shi}
%\note{}
%\seealso{}
\examples{
library(Rsubread)
reads <- system.file("extdata","reads.txt",package="Rsubread")
x <- qualityScores(filename=reads,nreads=1000)
boxplot(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

