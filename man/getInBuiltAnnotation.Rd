\name{getInBuiltAnnotation}
\alias{getInBuiltAnnotation}
\title{Retrieve In-Built Annotations}
\description{Retrieve an in-built annotation and save it to a data frame.}
\usage{
getInBuiltAnnotation(annotation = "mm10")
}
\arguments{
  \item{annotation}{ a character string specifying the in-built annotation to be retrieved. It has four possible values including \code{mm10}, \code{mm9}, \code{hg38} and \code{hg19}, corresponding to the NCBI RefSeq annotations for genomes `mm10', `mm9', `hg38' and `hg19', respectively. \code{mm10} by default.}
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
x <- getInBuiltAnnotation("hg38")
x[1:5,]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
