\name{propmapped}
\alias{propmapped}
\title{Calculate the Proportion of Mapped Reads or Fragments in SAM/BAM Files}
\description{Number of mapped reads/fragments will be counted and fraction of such reads/fragments will be calculated.}
\usage{
propmapped(

    files,
    countFragments = TRUE,
    properlyPaired = FALSE,
    verbose = FALSE)
}
\arguments{
  \item{files}{ a character vector giving the names of SAM/BAM format files. Format of input files is automatically determined by the function.}
  \item{countFragments}{ logical indicating whether reads or fragments (read pairs) should be counted. If \code{TRUE}, fragments will be counted when paired-end read data are provided. This function automatically detects if the data are single end or paired end. For single end data, each read is treated as a fragment and therefore the value of this parameter should be set to \code{TRUE}.}
  \item{properlyPaired}{ logical indicating if only properly paired reads will be counted. This is only applicable for paired end data. \code{FALSE} by default.}
  \item{verbose}{ logical indicating if verbose information should be displayed.}
}
\details{
This function uses the \code{FLAG} field in the SAM/BAM to look for mapped reads and count them.
Reads/fragments, which have more than one reported location, will be reported only once.

When counting single end reads, counting reads has the same meaning as counting fragments (the results are identical).

Note that \code{align} and \code{subjunc} return the same mapping statistics as \code{propmapped} as a by-product of the alignment process, so saving the output from \code{align} or \code{subjunc} will be faster and more informative than running \code{propmapped} separately.
}
\value{ 
A data frame containing the total number of reads, number of mapped reads and proportion of mapped reads for each library.
}
%\references{
%}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
\examples{
# build an index using the sample reference sequence provided in the package
# and save it to the current directory
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)

# align the sample read data provided in this package to the sample reference
# and save the mapping results to the current directory
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")
stat <- align(index="./reference_index",readfile1=reads,output_file="./Rsubread_alignment.BAM")
stat

# get the percentage of successfully mapped reads
propmapped("./Rsubread_alignment.BAM")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
