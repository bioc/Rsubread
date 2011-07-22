\name{featureCounts}
\alias{featureCounts}
\title{Count the number of mapped reads for each feature}
\description{Summarize read counts to features including genes and exons}
\usage{
featureCounts(SAMfiles,type="gene",species="mm")
}
\arguments{
  \item{SAMfiles}{ a character vector giving names of SAM format files.}
  \item{type}{ a character string giving the feature type. Its value could be \code{gene} or \code{exon}.}
  \item{species}{ a character string specifying the species. Only \code{mm} is supported at the moment. }
  }
\details{
This function takes as input a set of SAM format files and assigns reads included in them to the features.
The currently supported features include \code{gene} and \code{exon}.
Here, \code{gene} is the aggregation of all its exons, i.e. intron regions are excluded.

Currently only NCBI annotation for mouse is supported.
Work is under way to support more speices (e.g. human) and more annotations (e.g. Ensemble annotation).
}
\value{ 
A data frame containing read counts for each feature.
}
%\references{
%}
\author{Wei Shi}
%\note{}
%\seealso{}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
