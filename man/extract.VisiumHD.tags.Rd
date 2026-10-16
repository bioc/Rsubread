\name{extract.VisiumHD.tags}
\alias{extract.VisiumHD.tags}
\title{Extract Visium HD barcode tags from a Space Ranger BAM file}
\usage{
extract.VisiumHD.tags(SpaceRanger.bam.file, Output.text.file)
}
\arguments{
  \item{SpaceRanger.bam.file}{Space Ranger BAM file containing Visium HD read
  tags from Space Ranger mapping results.}
  \item{Output.text.file}{Output text file receiving the corresponding 1R, 1Y,
  and cell-barcode tags. An existing file is replaced.}
}
\description{Extracts Visium HD barcode tags from a Space Ranger BAM file and
writes a barcode reference file for Visium HD cell counting.}
\details{Requires \command{samtools} to be available on \code{PATH}. The
output file has three tab-separated columns: 1R tag, 1Y tag, and cell barcode.}
