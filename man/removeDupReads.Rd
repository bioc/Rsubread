\name{removeDupReads}
\alias{removeDupReads}
\title{Remove Sequencing Reads Mapped to Identical Locations}
\description{Remove reads which are mapped to identical locations, using mapping location of the first base of each read.}
\usage{
removeDupReads(

    # basic input/output options
    inputFile,
    threshold = 50,
    outputFile,
    outputFormat = "BAM")
}
\arguments{
  \item{inputFile}{ a character string giving the name of a SAM or BAM format input file.}
  \item{threshold}{ a numeric value giving the threshold for removing duplicated reads, 50 by default. Reads will be removed if they are found to be duplicated equal to or more than \code{threshold} times.}
  \item{outputFile}{ a character string giving the base name of output files.}
  \item{outputFormat}{ a character string giving the format of the output file. \code{BAM} and \code{SAM} are accepted. \code{BAM} by default.}
}
\details{
This function uses the mapping location of first base of each read to find duplicated reads.
Reads are removed if they are duplicated more than \code{threshold} number of times.
}
\value{
This function generates a BAM or SAM file including the remaining reads after duplicate removal.
}
\author{Yang Liao and Wei Shi}
