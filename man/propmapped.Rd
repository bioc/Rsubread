\name{propmapped}
\alias{propmapped}
\title{Obtain the proportion of mapped reads}
\description{Use mapping information stored in a SAM format file to count the number of mapped reads}
\usage{
propmapped(samfiles)
}
\arguments{
  \item{samfiles}{ a character vector giving the names of SAM format files.}
  }
\details{
This function counts of number of mapped reads using the mapping informtion stored in SAM format files. 
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
