\name{detectionCall}
\alias{detectionCall}
\title{Determine Detection P Values for Each Gene in an RNA-seq Dataset}
\description{Use GC content adjusted background read counts to determine the detection p values for each gene}

\usage{
detectionCall(

    dataset,
    species = "hg",
    plot = FALSE)
}
\arguments{
  \item{dataset}{ a character string giving the filename of a SAM format file, which is the output of read alignment.}
  \item{species}{ a character string specifing the species. Options are \code{hg} and \code{mm}.} 
  \item{plot}{ logical indicating whether a density plot of detection p values will be generated.}
}

\value{A data frame which includes detection p values and annotation information for each gene.}

%\references{}
\author{Zhiyin Dai and Wei Shi}
%\note{}
%\seealso{}
% Add one or more standard keywords, see file 'KEYWORDS' in the 
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

