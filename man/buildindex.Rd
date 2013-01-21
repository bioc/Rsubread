\name{buildindex}
\alias{buildindex}
\title{Build index for a reference genome for read mapping}
\description{This function builds an index for the reference genome using a hash table. Memory usage can be adjusted when loading the index into the computer memory.}
\usage{
buildindex(basename,reference,colorspace=FALSE,memory=3700,TH_subread=24)
}
\arguments{
  \item{basename}{ character string giving the basename of created index files.}
  \item{reference}{ charater string giving the name of the file containing all the refernece sequences.}
  \item{colorspace}{ logical. If \code{TRUE}, a color space index will be built. Otherwise, a base space index will be built.}
  \item{memory}{ numeric value specifying the amount of memory to be requested in megabytes. 3700 MB by default.}
  \item{TH_subread}{ numeric value specifying the threshold for removing highly repetitive subreads (16 mers). 24 by default. Subreads will be excluded from the index if they occur more than threshold number of times in the genome.}
  }
\details{
A hash table will be built for the reference genome.
Keys in the hash table are the 16bp sequences and hash values are their chromosomal locations.
A 16bp sequence could have one or more than one chromosomal locations.
They are all recorded in the hash table.
Non-informative 16bp sequences, which are highly repetitive in the reference genome (occurring 24 time or more in the genome), are excluded from the hash table.

After the index is built, reads can then be mapped to the reference genome by using \code{\link{align}} function.

It takes around 1 hour to build an index for human genome. 
}
\value{ 
Index files with basename provided in \code{basename}.
These files are saved in the current directory.
}
%\references{}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
\examples{
# build an index using the sample reference sequence provided in the package
# and save it to the current directory
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
