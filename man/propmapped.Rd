\name{propmapped}
\alias{propmapped}
\title{Calculate the proportion of mapped reads/fragments in SAM/BAM files}
\description{Number of mapped reads/fragments will be counted and fraction of such reads/fragments will be calculated.}
\usage{
propmapped(files,countFragments=TRUE,properlyPaired=FALSE)
}
\arguments{
  \item{files}{ a character vector giving the names of SAM/BAM format files. Format of input files is automatically determined by the function.}
  \item{countFragments}{ logical, indicating if reads for fragments (read pairs) should be counted. If \code{TRUE}, fragments will be counted when paired-end read data are provided. The function automatically detects if the data are single end or paired end.}
  \item{properlyPaired}{ logical, indicating if only properly paired reads will be counted. This is only applicable for paired end data. \code{FALSE} by default.}
}
\details{
This function uses the \code{FLAG} field in the SAM/BAM to look for mapped reads and count them.
Reads/fragments, which have more than one reported location, will be reported only once.

When counting single end reads, counting reads has the same meaning with counting fragments (the results are identical). 
}
\value{ 
A data frame containing the total number of reads, number of mapped reads and proportion of mapped reads for each library.
}
%\references{
%}
\author{Wei Shi}
%\note{}
%\seealso{}
\examples{
# build an index using the sample reference sequence provided in the package
# and save it to the current directory
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)

# align the sample read data provided in this packge to the sample reference
# and save the mapping results to the current directory
reads <- system.file("extdata","reads.txt",package="Rsubread")
align(index="./reference_index",readfile1=reads,output_file="./Rsubread_alignment.SAM")

# get the percentage of successfully mapped reads
propmapped("./Rsubread_alignment.SAM")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
