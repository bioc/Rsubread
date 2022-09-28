\name{getInBuiltAnnotation}
\alias{getInBuiltAnnotation}
\title{Retrieve inbuilt gene annotations}
\description{Retrieve an inbuilt gene annotation and save it to a data frame.}
\usage{
getInBuiltAnnotation(annotation = "mm39")
}
\arguments{
  \item{annotation}{ a character string specifying the inbuilt gene annotation to be retrieved. It has five possible values including \code{mm39}, \code{mm10}, \code{mm9}, \code{hg38} and \code{hg19}, which correspond to the NCBI RefSeq gene annotations for genomes `mm39', `mm10', `mm9', `hg38' and `hg19', respectively. \code{mm39} by default.}
}
\details{
The inbuilt gene annotations included in the Rsubread package were generated based on the NCBI RefSeq gene annotations.
The inbuilt annotations are much simpler than the original RefSeq annotations in that they only contain essential data required for the processing of sequencing data such as RNA-seq data.
Also in each inbuilt annotation, duplicate exons within the same gene are removed and overlapping exons within the same gene are merged.
For more information about the inbuilt annotations, please refer to the help page for \code{\link{featureCounts}} function.
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
head(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
