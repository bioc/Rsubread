\name{buildindex}
\alias{buildindex}
\title{Build index for a reference genome for read mapping}
\description{This is an R wrapper function for building index for a reference genome. This function calls the underlying C function.}
\usage{
buildindex(basename,reference,colorspace=FALSE,memory=3700)
}
\arguments{
  \item{basename}{ character string giving the basename of created index files.}
  \item{reference}{ charater string giving the name of the file containing all the refernece sequences.}
  \item{colorspace}{ logical. If \code{TRUE}, a color space index will be built. Otherwise, a base space index will be built.}
  \item{memory}{ numeric value specifying the amount of memory to be requested in gigabytes}
  }
\details{
A hash table will be built for the reference genome.
Keys in the hash table are the 16bp sequences and hash values are their chromosomal locations.
A 16bp sequence could have one or more than one chromosomal locations.
They are all recorded in the hash table.
Non-informative 16bp sequences, which are highly repetitive in the reference genome, are not included in the hash table.

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
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
path <- system.file("extdata",package="Rsubread")
buildindex(basename=file.path(path,"reference_index"),reference=ref)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
