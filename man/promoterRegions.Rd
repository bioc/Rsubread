\name{promoterRegions}
\alias{promoterRegions}
\title{Generate Annotation for Promoter Regions of Genes}
\description{Create a SAF data-frame of genewise promoter regions.}
\usage{
promoterRegions(

    annotation = "mm10",
    upstream = 3000L,
    downstream = 2000L)
}
\arguments{
  \item{annotation}{ a data.frame containing gene annotation in SAF format or a character string giving the name of a genome with built-in annotation. If using built-in annotation, the character string should be one of the following: \code{mm10}, \code{mm9}, \code{hg38} or \code{hg19} corresponding to the NCBI RefSeq annotations for the genomes `mm10', `mm9', `hg38' and `hg19', respectively.}
  \item{upstream}{ an integer giving the number of upstream bases that will be inclued in the promoter region generated for each gene. These bases are taken immediately upstream (5' end) from transcriptional start site of each gene.}
  \item{downstream}{ an integer giving the number of downstream bases that will be inclued in the promoter region generated for each gene. These bases are taken immediately downstream (3' end) from transcriptional start site of each gene.}
}
\details{
This function takes as input a SAF format gene annotation and produces a SAF format data.frame containing the chromosomal coordinates of the specified promoter region for each gene.
See \code{\link{featureCounts}} for definition of the SAF format.

Regardless of the \code{upstream} or \code{downstream} values, the downstream end of the region never extends past the end of the gene and the upstream end never extends outside the relevant chromosome.
Setting \code{downstream} to an infinite or large value will cause the body of each gene to be included.
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
# To get whole gene bodies for the mouse genome:
x <- promoterRegions("mm10", upstream = 0, downstream = Inf)
head(x)
}
