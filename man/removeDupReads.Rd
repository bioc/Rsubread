\name{removeDupReads}
\alias{removeDupReads}
\title{Remove sequencing reads which are mapped to identical locations}
\description{Remove reads which are mapped to identical locations, using mapping location of the first base of each read.}
\usage{
removeDupReads(SAMfile,threshold=50,outputFile)
}
\arguments{
  \item{SAMfile}{ a character string giving the name of a SAM format input file.}
  \item{threshold}{ a numeric value giving the threshold for removing duplicated reads, 50 by default. Reads will be removed if they are found to be duplicated equal to or more than \code{threshold} times.}
  \item{outputFile}{ a character string giving the base name of output files.}
  }
\details{
This function uses the mapping location of the first base of each read to find duplicated reads (mapped to the same location).
Reads will be removed from the SAM file if they are found to be duplicated equal to or more than \code{threshold} times.\\

It requires the chromosome names included in the SAM file to be 'chr1', 'chr2', 'chr3', ..., 'chrX', and 'chrY'.
}
\value{
A SAM file which includes reads which are not duplicated (named \code{outputFile}.NoneDupReads), and a text file which includes removed duplicated reads (named \code{outputFile}.DupReads.txt).
}
%\references{
%}
\author{Wei Shi and Jenny Zhiyin Dai}
%\note{}
%\seealso{}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
