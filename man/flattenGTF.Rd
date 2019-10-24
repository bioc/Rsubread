\name{flattenGTF}
\alias{flattenGTF}
\title{Flatten Features in GTF or GFF Annotation Files}
\description{Convert a GTF/GFF annotation to SAF format and flatten overlapping features.}
\usage{
flattenGTF(

    # basic input/output options
    GTFfile,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",

    # the option specifying the merging algorithm
    method = "merge")
}
\arguments{
  \item{GTFfile}{ a character string giving the name of a GTF file as input.}
  \item{GTF.featureType}{ a character string giving the feature type used to select rows in a GTF annotation. \code{"exon"} by default. Feature types can be found in the third column of a GTF annotation.}
  \item{GTF.attrType}{ a character string giving the attribute type in a GTF annotation which will be used to group features. \code{"gene_id"} by default. Attributes can be found in the ninth column of a GTF annotation.}
  \item{method}{ a character string specifying the method for how to flatten the GTF file. The method can be either \code{"merge"} or \code{"chop"}. \code{"merge"} by default. See the \code{details} section for more information.}
}

\details{
This function locates features in a GTF annotation via \code{GTF.featureType} and then groups them into meta-features via \code{GTF.attrType}.

When \code{method="merge"}, the overlapping features found in a meta-feature will be merged to form a single continuous feature encompassing all the overlapping features.
If \code{method="chop"}, overlapping features will be chopped into multiple non-overlapping bins.
Here is an example to illustrate the differences between the two methods.
Say there are three exons belonging to the same gene and the coordinates of these exons are [100, 200], [150, 500] and [250, 400].
When running on \code{"merge"} mode, a single exon will be returned for this gene and the coordinates of this exon are [100, 500], representing the union of the three original exons.
When running on \code{"chop"} mode, five non-overlapping bins will be returned, including [100, 149], [150, 200], [201, 249], [250, 400] and [401, 500].
Intervals of these bins are determined by start and end coordinates of the three original exons (100, 150, 200, 250, 400 and 500).

Output of this function is a SAF format annotation which can be fed to \code{\link{featureCounts}} function for read counting.
Description to SAF format annotation can also be found in \code{\link{featureCounts}}.

}
\value{
A \code{data.frame} including a SAF format annotation in which features included in each meta-feature are all distinct.
}
\author{Yang Liao and Wei Shi}
\seealso{\code{\link{featureCounts}}}
%\examples{}
