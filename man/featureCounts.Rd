\name{featureCounts}
\alias{featureCounts}
\title{Count the number of mapped reads for each feature}
\description{Summarize read counts to features including genes and exons}
\usage{
featureCounts(SAMfiles,type="gene",species="mm",annot=NULL)
}
\arguments{
  \item{SAMfiles}{ a character vector giving names of SAM format files.}
  \item{type}{ a character string giving the feature type. Its value could be \code{gene} or \code{exon}.}
  \item{species}{ a character string specifying the species. It can be \code{mm} or \code{hg}. Values of this argument determines which in-built annotation file will be used, if \code{annot} is \code{NULL}.}
  \item{annot}{ a character string giving the name of an annotation file provided by the user, or a data frame containing annotation information. If \code{NULL} (by default), the in-built annotation file will be used. See below for more details.}
  }
\details{
This function takes as input a set of SAM format files and assigns reads to the features.
Currently, only feature types including \code{gene} and \code{exon} are supported.
\code{gene} is the aggregation of all the exons for each gene.

There are two in-built annotation files which are used by this function to summarize reads for genes or exons for mouse and human, respectively.
These annotation files include the exon annotation information downloaded from NCBI Build 37.2, including Entrez gene identifier and chromosomal coordinates for each exon.
The \code{species} argument specifies which annotation file should be used.
These in-built annotation will be used when \code{annot} is \code{NULL}.

Users can also provide their own annotation for read summarization, via the \code{annot} argument.
In this case, user provided annotation will override the in-built annotation.
Users can either provide the name of their own annotation file, or they can provide a data frame which contain annotation information.
If a data frame is provided, at least four columns, named \code{entrezid}, \code{chromosome}, \code{chr_start} and \code{chr_stop}, are required.
If users provide an annotation file, then it should be a tab delimited file and its first four columns should be named as \code{entrezid}, \code{chromosome}, \code{chr_start} and \code{chr_stop}, respectively.
The file has to be sorted by the \code{chromosome} column (either increasing or decreasing order is OK), and chromosome names should not contain more than 40 characters.
Below is an example for the annotation file provided by the user (sorted by the chromosome column):
\preformatted{
entrezid	chromosome	chr_start	chr_stop
497097	chr1	3204563	3207049
497097	chr1	3411783	3411982
497097	chr1	3660633	3661579
100503874	chr1	3637390	3640590
100503874	chr1	3648928	3648985
100038431	chr1	3670236	3671869
...
}

Although this function is designed for summarizing reads from RNA-seq experiments, it can be used to summarize reads from other next-gen sequencing experiments as well, for example ChIP-seq or other DNA sequencing experiments.
Simply by setting \code{type} to \code{exon} and providing an annotation, this function will yield numbers of mapped reads for each feature. 
}
\value{ 
A list with the following components:
\tabular{ll}{
  \code{counts}:\tab a data matrix containing read counts for each gene or exon.\cr
  \code{annotation}:\tab a data frame containing Entrez gene identifers, gene/exon length etc.\cr
  \code{targets}:\tab a character vector giving sample information.
}
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
