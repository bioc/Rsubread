\name{txUnique}
\alias{txUnique}
\title{txUnique: counting unique bases of transcripts}
\description{This function counts the unique bases and total bases of transcripts in a GTF annotation file
}
\usage{
txUnique(GTF_Annotation_File, Feature_Type = "exon", Gene_ID_Attribute = "gene_id", Transcript_ID_Attribute = "transcript_id")
}
\arguments{
  \item{GTF_Annotation_File}{ character string giving the name of a GTF file containing the transcript-level annotations.}
  \item{Feature_Type}{ character string specifying the type of annotations. Only the annotations of the given type will be loaded from the GTF file. \code{exon} by default. The third column in the GTF file contains the type of the annotation in each line}
  \item{Gene_ID_Attribute}{ character string specifying the attribute name in the 9-th column in the GTF file for the gene identifier. \code{gene_id} by default.}
  \item{Transcript_ID_Attribute}{ character string specifying the attribute name in the 9-th column in the GTF file for the transcript identifier. \code{transcript_id} by default.}
}
\details{
txUnique calculates the total number of bases in each transcript, as well as the number of bases that are only in this transcript, but not included in other transcript(s) belonging to the same gene.
}
\value{
A data matrix, where each row corresponds to a transcript in the GTP annotation file. There are four columns for each transcript: the identifier of the gene, the identifier of the transcript, the numbers of unique bases and total bases in the transcript.
}
%\references{}
\author{Yang Liao and Wei Shi}
%\note{}
%\seealso{}
\examples{
library(Rsubread)
base_counts <- txUnique("hg19-Ensembl.GTF")
base_counts[1:10,]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

