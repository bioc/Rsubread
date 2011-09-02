\name{processExons}
\alias{processExons}
\title{Obtain chromosomal coordiates of each exon using NCBI annotation}
\description{This is for internal use.}
\usage{
processExons(filename="human_seq_gene.md", species="hg")
}
\arguments{
  \item{filename}{a character string giving the name of input .md file (NCBI annotation file)}
  \item{species}{a character string specifying the species}
}
\details{
The NCBI annotation file gives the chromosomal coordinates of UTR (Untranslated region) and CDS (Coding sequence).
This function uses these information to derive the chromosomal coordinates of exons.
The first and last exons of genes usually contain both UTR sequence and CDS sequence.
}

\value{ 
A text file containing chromosomal coordinates of each exon.  
}
%\references{}
\author{Zhiyin Dai and Wei Shi}
%\note{}
%\seealso{}
% Add one or more standard keywords, see file 'KEYWORDS' in the 
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

