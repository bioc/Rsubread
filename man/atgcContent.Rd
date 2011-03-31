\name{atgcContent}
\alias{atgcContent}
\title{Calculate percentages of nucletodies A, T, G and C in a sequencing read datafile}
\description{Calculate percentages of nucletodies A, T, G and C
}
\usage{
atgcContent(filename, basewise=FALSE)
}
\arguments{
  \item{filename}{ character string giving the name of input FASTQ/FASTA file}
  \item{basewise}{ logical. If \code{TRUE}, nucleotide percentages will be calculated for each base position in the read across all the reads. By default, percentages are calculated for the entire dataset.}
}
\details{
Sequencing reads could contain letter "N" besides "A", "T", "G" and "C".
Percentage of "N" in the read dataset is calcuated as well.

The \code{basewise} calculation is useful for examining the GC bias towards the base position in the read.
By default, the percentages of nucleotides in the entire dataset will be reported.
}

\value{ 
A named vector containing percentages for each nucleotide type if \code{basewise} is \code{FALSE}.
Otherwise, a data matrix containing nucleotide percentages for each base position of the reads. 
}
%\references{}
\author{Zhiyin Dai and Wei Shi}
%\note{}
%\seealso{}
\examples{
library(Rsubread)
reads <- system.file("extdata","reads.txt",package="Rsubread")
# Fraction of A,T,G and C in the entire dataset
x <- atgcContent(filename=reads,basewise=FALSE)
# Fraction of A,T,G and C at each base location across all the reads
xb <- atgcContent(filename=reads,basewise=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the 
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

