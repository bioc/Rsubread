\name{atgcContent}
\alias{atgcContent}
\title{Calculate Percentages of Nucleotides in Reads}
\description{Calculate percentages of nucleotides A, T, G and C in the input FASTQ or FASTQ file.  }
\usage{
atgcContent(

    filename,
    basewise = FALSE)
}
\arguments{
  \item{filename}{ a character string giving the name of input FASTQ/FASTA file}
  \item{basewise}{ logical specifying how the percentages are calculated. If \code{TRUE}, nucleotide percentages will be calculated for each base position in the read across all the reads. By default, percentages are calculated for the entire dataset.}
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
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the 
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

