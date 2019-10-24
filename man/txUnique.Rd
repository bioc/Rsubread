\name{txUnique}
\alias{txUnique}
\title{Count Number of Bases Unique to Each Transcript}
\description{For each transcript, number of unique bases and total number of bases will be reported. }
\usage{
txUnique(

    # basic input/output options
    GTF_Annotation_File,
    Feature_Type = "exon",
    Gene_ID_Attribute = "gene_id",
    Transcript_ID_Attribute = "transcript_id")
}
\arguments{
  \item{GTF_Annotation_File}{ a character string giving the name of a GTF file containing the transcript-level annotations.}
  \item{Feature_Type}{ a character string specifying the type of annotations. Only the annotations of the given type will be loaded from the GTF file. \code{exon} by default. The third column in the GTF file contains the type of the annotation in each line}
  \item{Gene_ID_Attribute}{ a character string specifying the attribute name in the 9-th column in the GTF file for the gene identifier. \code{gene_id} by default.}
  \item{Transcript_ID_Attribute}{ a character string specifying the attribute name in the 9-th column in the GTF file for the transcript identifier. \code{transcript_id} by default.}
}
\details{
This function compares the transcripts belonging to the same gene and then counts the number of bases unique to each transcript. It also reports the total number of bases for each transcript. When a transcript is found to contain overlapping exons, the overlapping exons will be merged and each overlapping base will be counted only once.
}
\value{
A data matrix in which each row represents a transcript. The matrix includes four columns: gene identifier, transcript identifier, number of bases unique to each transcript and total number of bases each transcript has.
}
%\references{}
\author{Yang Liao and Wei Shi}
%\note{}
%\seealso{}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

