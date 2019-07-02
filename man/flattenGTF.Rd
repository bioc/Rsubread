\name{flattenGTF}
\alias{flattenGTF}
\title{Flatten features provided in a GTF annotation}
\description{Overlapping features in a GTF annotation, such as overlapping exons from the same gene, are collapsed into a single feature that covers all overlapping features.}
\usage{
flattenGTF(GTFfile, GTF.featureType = "exon", GTF.attrType = "gene_id", method = "merge")
}
\arguments{
  \item{GTFfile}{ a character string giving the name of a GTF file as input.}
  \item{GTF.featureType}{ a character string giving the feature type used to select rows in a GTF annotation. \code{"exon"} by default. Feature types can be found in the third column of a GTF annotation.}
  \item{GTF.attrType}{ a character string giving the attribute type in a GTF annotation which will be used to group features. \code{"gene_id"} by default. Attributes can be found in the ninth column of a GTF annotation.}
  \item{method}{ a character string specifying the method for flattening the GTF file. The method can be either \code{"merge"} that takes the union of all the overlapping features of the same meta-feature to generate a single continuous feature, or \code{"chop"} that splits the union of the overlapping exons by their edges. See the \code{details} section for more information. \code{"merge"} by default.}
}

\details{
This function locates features in a GTF annotation via \code{GTF.featureType} and then groups them into meta-features via \code{GTF.attrType}.
Default values of these two parameters will make exons be identified and then grouped into genes.

This function supports two methods for flattening the input features: \code{"merge"} and \code{"chop"}. When running on the \code{"merge"} mode, if there is overlap found between features belonging to the same meta-feature, this function will merge all overlaping features to form a single larger feature that covers all overlapping features. When running on the \code{"chop"} mode, the function may generate multiple non-overlapping features from the overlaping features belonging to the same meta-feature in the input, chopped by the edges (the start and end points) of the input exons.

Say, there are three exons belonging to the same gene, and the coordinates of the exons are [100, 200], [150, 500] and [250, 400]. When running on the \code{"merge"} mode, this function simply returns a single exon, [100, 500], which is the union of the input exons. On the other hand, when running on the \code{"chop"} mode, this function returns five non-overlapping exons: [100, 149], [150, 200], [201, 249], [250, 400] and [401, 500] for the three input exons. Namely, the union of the input exons is then chopped by all the edges in the input exons: 100, 150, 200, 250, 400 and 500.

Output of this function is a SAF format annotation which can be fed to \code{\link{featureCounts}} function for read counting.
Description to SAF format annotation can also be found in \code{\link{featureCounts}}.

}
\value{
A \code{data.frame} including a SAF format annotation.
}
\author{Yang Liao and Wei Shi}
\seealso{\code{\link{featureCounts}}}
%\examples{}
