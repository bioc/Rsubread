\name{removeDupReads}
\alias{removeDupReads}
\title{Remove sequencing reads which are mapped to identical locations}
\description{Remove reads which are mapped to identical locations, using mapping location of the first base of each read.}
\usage{
removeDupReads(SAMfile,threshold=50,nthreads=1,outputFile)
}
\arguments{
  \item{SAMfile}{ a character string giving the name of a SAM format input file.}
  \item{threshold}{ a numeric value giving the threshold for removing duplicated reads, 50 by default. Reads will be removed if they are found to be duplicated equal to or more than \code{threshold} times.}
  \item{nthreads}{ a numeric value giving the number of CPUs used.}
  \item{outputFile}{ a character string giving the base name of output files.}
  }
\details{
This function uses the mapping location of first base of each read to find duplicated reads.
Reads are removed if they are duplicated more than \code{threshold} times.\\
}
\value{
A SAM file which includes the remaining reads after duplicate removal.
}
%\references{
%}
\author{Yang Liao and Wei Shi}
%\note{}
%\seealso{}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
