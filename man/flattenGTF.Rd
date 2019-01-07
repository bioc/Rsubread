\name{flattenGTF}
\alias{flattenGTF}
\title{Flatten features provided in a GTF annotation}
\description{Overlapping features in a GTF annotation, such as overlapping exons from the same gene, are collapsed into a single feature that covers all overlapping features.}
\usage{
flattenGTF(GTFfile, GTF.featureType = "exon", GTF.attrType = "gene_id")
}
\arguments{
  \item{GTFfile}{ a character string giving the name of a GTF file as input.}
  \item{GTF.featureType}{ a character string giving the feature type used to select rows in a GTF annotation. \code{"exon"} by default. Feature types can be found in the third column of a GTF annotation.}
  \item{GTF.attrType}{ a character string giving the attribute type in a GTF annotation which will be used to group features. \code{"gene_id"} by default. Attributes can be found in the ninth column of a GTF annotation.}
}

\details{
This function locates features in a GTF annotation via \code{GTF.featureType} and then groups them into meta-features via \code{GTF.attrType}.
Default values of these two parameters will make exons be identified and then grouped into genes.

If there is overlap found between features belonging to the same meta-feature, this function will merge all overlaping features to form a single larger feature that covers all overlapping features.
Output of this function is a SAF format annotation which can be fed to \code{\link{featureCounts}} function for read counting.
Description to SAF format annotation can also be found in \code{\link{featureCounts}}.
}
\value{
A \code{data.frame} including a SAF format annotation.
}
\author{Yang Liao and Wei Shi}
\seealso{\code{\link{featureCounts}}}
%\examples{}
