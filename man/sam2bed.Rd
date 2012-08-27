\name{sam2bed}
\alias{sam2bed}
\title{Convert SAM format file to BED format}
\description{SAM to BED conversion}
\usage{
sam2bed(samfile,bedfile,readlen)
}
\arguments{
  \item{samfile}{ character string giving the name of input file. Input format should be in SAM format.}
  \item{bedfile}{ character string giving the name of output file. Output file is in BED format.}
  \item{readlen}{ numeric value giving the length of reads included in the input file.}
  }
\details{
SAM format is the de facto standard format of output from read aligner.
This format not only includes the mapping coordinates of the reads but also includes other using information such as mapping quality, CIGAR information and so on.  
This function converts a SAM format file to a BED format file, which can then be displayed in a genome browser like UCSC genome browser, IGB, IGV etc.
}
\value{ 
A BED format file.
}
%\references{}
\author{Wei Shi}
%\note{}
%\seealso{}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
