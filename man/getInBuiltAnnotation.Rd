\name{getInBuiltAnnotation}
\alias{getInBuiltAnnotation}
\title{Retrieve in-built annotations provided by featureCounts function}
\description{Retrieve an in-built annotation and save it to a data frame}
\usage{
getInBuiltAnnotation(annotation="mm9")
}
\arguments{
  \item{annotation}{ a character string specifying the in-built annotation to be retrieved. It has three possible values including \code{mm9}, \code{mm10} and \code{hg19}, corresponding to the NCBI RefSeq annotations for genomes `mm9', `mm10' and `hg19', respectively. \code{mm9} by default.}
}
\details{
The \code{\link{featureCounts}} read summarization function provides in-built annotations for conveniently summarizing reads to genes or exons, and this function allows users to have access to those in-built annotations.

For more information about these annotations, please refer to the help page for \code{\link{featureCounts}} function.
}
\value{ 
A data frame with five columns including \code{GeneID}, \code{Chr}, \code{Start}, \code{End} and \code{Strand}.
}
%\references{
%}
\author{Wei Shi}
%\note{}
\seealso{
\code{\link{featureCounts}}}
\examples{
library(Rsubread)
x <- getInBuiltAnnotation("hg19")
x[1:5,]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
