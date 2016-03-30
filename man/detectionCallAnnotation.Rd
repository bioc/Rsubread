\name{detectionCallAnnotation}
\alias{detectionCallAnnotation}
\title{Generate annotation data used for calculating detection p values}
\description{This is for internal use only.}
\usage{
detectionCallAnnotation(species="hg", binsize=2000)
}
\arguments{
  \item{species}{character string specifying the species to analyase}
  \item{binsize}{binsize of integenic region}
}
\details{
This is an internal function and should not be called by users directly.

It takes as input the annotation files produced by \code{\link{processExons}} function, calculates GC percentages for each exon of genes and also for intergenic regions and add GC info into the annotations.
The new annotation data are then saved to files which can be used by \code{\link{detectionCall}} function for calling absolutely expressed genes.
}
\value{
Two annotation files, which contain GC content for exons of genes and for intergenic regions respectively, are written to the current working directory. This function returns a \code{NULL} object.
}
%\references{}
\author{Zhiyin Dai and Wei Shi}
%\note{}
%\seealso{}
% Add one or more standard keywords, see file 'KEYWORDS' in the 
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

