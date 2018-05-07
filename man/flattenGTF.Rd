\name{flattenGTF}
\alias{flattenGTF}
\title{Merge the overlapping exons of each gene in a GTF file}
\description{The overlapping exons of each gene were merged in the input GTF file. A data frame containing merged exons is returned.}
\usage{
flattenGTF(input.file.name, GTF.featureType = "exon", GTF.attrType = "gene_id")
}
\arguments{
  \item{input.file.name}{ a character string giving the name of a GTF file as input.}
  \item{GTF.featureType}{ a character string giving the feature type used to select rows in the GTF annotation which will be used for read summarization. \code{"exon"} by default. Feature types can be found in the third column of a GTF annotation.}
  \item{GTF.attrType}{ a character string giving the attribute type in the GTF annotation which will be used to group features (eg. exons) into meta-features (eg. genes). \code{"gene_id"} by default. Attributes can be found in the ninth column of a GTF annotation.}
}

\details{
\code{flattenGTF} merges the overlapping exons of each gene in a GTF file. Exons on different strands are not merged.

The return value from \code{flattenGTF} is a data frame. It can be directly given to \code{featureCounts} as the external annotation.
}
\author{Wei Shi and Yang Liao}
\seealso{
\code{\link{featureCounts}}}
\examples{\dontrun{
library(Rsubread)

flattenGTF("overlapping_GTF_annotation.GTF")
}}
