\name{processExons}
\alias{processExons}
\title{Obtain Chromosomal Coordiates of Each Exon Using NCBI Annotation}
\description{This is for internal use. It converts the NCBI exon annotations into data files for use in the \code{detectionCall} function.}
\usage{
processExons(

    filename = "human_seq_gene.md",
    species = "hg")
}
\arguments{
  \item{filename}{ a character string giving the name of input .md file (NCBI annotation file)}
  \item{species}{ a character string specifying the species}
}
\details{
This is an internal function and should not be called by users directly.

It processes annotation of genes in human genome GRCh37/hg19 or mouse genome GRCm37/mm9.
The annotation is available for download from the following links (these annotations include chromosomal coordinates of UTR and CDS regions of genes):

\url{ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.37.2/mapview/seq_gene.md.gz}

\url{ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.2/mapview/seq_gene.md.gz}

This function finds the chromosomal coordinates of intergenic regions (regions between neighbouring genes) and then outputs them to a file.
It also outputs to a file chromosomal coordinates of exons of genes by concatenating UTRs with CDSs and merging overlapping CDSs within each gene.
The generated annotation files will then be used by \code{\link{detectionCallAnnotation}} function to produce annotation data required by \code{\link{detectionCall}} function.
}
\value{ 
Two annotation files are written to the current working directory. This function returns a \code{NULL} object.
}
%\references{}
\author{Zhiyin Dai and Wei Shi}
%\note{}
%\seealso{}
% Add one or more standard keywords, see file 'KEYWORDS' in the 
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

