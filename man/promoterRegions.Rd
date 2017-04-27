\name{promoterRegions}
\alias{promoterRegions}
\title{Generate annotation for promoter region of genes}
\description{Create a SAF data.frame of genewise promoter regions}
\usage{
promoterRegions <- function(annotation="mm10", upstream=3000L, downstream=2000L)
}
\arguments{
\item{annotation}{ a character string giving the name of a genome with built-in annotation or a SAF format gene annotation. If using built-in annotation, the character string should be one of the following: \code{mm10}, \code{mm9}, \code{hg38} or \code{hg19}. They correspond to the NCBI RefSeq annotations for genomes `mm10', `mm9', `hg38' and `hg19', respectively. \code{mm10} by default.}
  \item{upstream}{ an integer giving the number of upstream bases that will be inclued in the promoter region generated for each gene. These bases are taken immediately upstream (5' end) from transcriptional start site of each gene.}
  \item{downstream}{ an integer giving the number of downstream bases that will be inclued in the promoter region generated for each gene. These bases are taken immediately downstream (3' end) from transcriptional start site of each gene.}
}
\details{
This function takes as input a SAF format gene annotation and produces a SAF format annotation that contains chromosomal coordinates of promoter region of each gene.
See \code{\link{featureCounts}} for definition of the SAF format.

Built-in annotations can be accessed via function \code{\link{getInBuiltAnnotation}}.
}
\value{ 
A SAF format \code{data.frame} with columns \code{GeneID}, \code{Chr}, \code{Start}, \code{End} and \code{Strand}.
}
%\references{
%}
\author{Gordon K Smyth}
%\note{}
\seealso{
\code{\link{featureCounts}}, \code{\link{getInBuiltAnnotation}}}
\examples{
library(Rsubread)
x <- promoterRegions()
head(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
