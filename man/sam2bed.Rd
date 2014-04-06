\name{sam2bed}
\alias{sam2bed}
\title{Convert a SAM format file to a BED format file}
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
This function converts a SAM format file to a BED format file, which can then be displayed in a genome browser like UCSC genome browser, IGB, IGV.
}
\value{ 
No value is produced but a BED format file is written to the current working directory. This file contains six columns including chromosomal name, start position, end position, name(`.'), mapping quality score and strandness.
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
