\name{qualityScores}
\alias{qualityScores}
\title{Extract Quality Score Data in a Sequencing Read Dataset}
\description{Extract quality strings and convert them to Phred scores for generating boxplots.
}
\usage{
qualityScores(

    # basic input/output options
    filename,
    input_format = "gzFASTQ",
    offset = 33,
    nreads = 10000)
}
\arguments{
  \item{filename}{ a character string giving the name of an input file containing sequence reads.}
  \item{input_format}{ a character string specifying format of the input file. \code{gzFASTQ} (gzipped FASTQ) by default. Acceptable formats include \code{gzFASTQ}, \code{FASTQ}, \code{SAM} and \code{BAM}. Character string is case insensitive.}
  \item{offset}{ a numeric value giving the offset added to the base-calling Phred scores. Possible values include \code{33} and \code{64}. By default, \code{33} is used.}
  \item{nreads}{ a numeric value giving the number of reads from which quality scores are extracted. \code{10000} by default. A value of \code{-1} indicates that quality scores will be extracted from all the reads.}
}
\details{
Quality scores of read bases are represented by ASCII characters in next-gen sequencing data.
This function extracts the quality characters from each base in each read and then converts them to Phred scores using the provided offset value (\code{offset}).

If the total number of reads in a dataset is \code{n}, then every \code{n}/\code{nreads} read is extracted from the input data.
}
\value{ 
A data matrix containing Phred scores for read bases. Rows in the matrix are reads and columns are base positions in each read.
}
%\references{}
\author{Wei Shi, Yang Liao and Zhiyin Dai}
%\note{}
%\seealso{}
\examples{
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")
x <- qualityScores(filename=reads,offset=64,nreads=1000)
x[1:10,1:10]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line

