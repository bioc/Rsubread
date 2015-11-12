\name{repair}
\alias{repair}
\title{Re-order paired-end reads to place reads from the same pair next to each other}
\description{Fast re-odering of paired-end reads using read names and mapping locations.}
\usage{
repair <- function(inFiles,inFormat="BAM",outFiles=paste(inFiles,"repair",sep="."),
addDummy=TRUE,fullData=TRUE,compress=FALSE,nthreads=8)
}
\arguments{
  \item{inFiles}{ a character vector giving names of input files. These files are typically location-sorted BAM files.}
  \item{inFormat}{ a character string specifying format of input files. Supported formats include \code{BAM} and \code{SAM}.}
  \item{outFiles}{ a character string giving names of output files. Re-ordered reads are saved to BAM-format files.}
  \item{addDummy}{ logical indicating if a dummy read will be added to each singleton read which has a missing pair in the input. \code{TRUE} by default.}
  \item{fullData}{ logical indicating if sequences and base-calling quality scores of reads will be included in the output. \code{TRUE} by default.}
  \item{compress}{ logical indicating if compression should be turned on when generating BAM output. \code{FALSE} by default.}
  \item{nthreads}{ a numeric value giving number of CPU threads. 8 by default.}
}
\details{
This function takes as input paired-end BAM or SAM files, re-orders reads to make reads from the same pair be adajacent to each other and then outputs the re-ordered reads into BAM files.

The function makes use of both read names and mapping information of reads (eg. mapping coordinates) to identify reads belonging to the same pair.
This makes sure that all paired-end reads are correctly re-ordered, especially those multi-mapping read pairs that include more than one reported alignment in the input.

The BAM files produced by this function are comptible with \code{\link{featureCounts}}, meaning that no read re-ordering will be performed when providing these files to \code{\link{featureCounts}}.
}
\value{
No value is produced but BAM files with re-ordered reads are written to the current working directory.
}
%\references{}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
