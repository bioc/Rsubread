\name{propmapped}
\alias{propmapped}
\title{Obtain the proportion of mapped reads}
\description{Use mapping information stored in a SAM format file to count the number of mapped reads}
\usage{
propmapped(samfile)
}
\arguments{
  \item{samfile}{ character string giving the name of a SAM format file.}
  }
\details{
This function uses the mapping information included in a SAM format file get the proportion mapped reads out of all the reads. 
}
\value{ 
Fraction of mapped reads is printed on the screen.
}
%\references{
%}
\author{Wei Shi}
%\note{}
%\seealso{}
\examples{
library(Rsubread)
results <- system.file("extdata","alignResults.SAM",package="Rsubread")
propmapped(results)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
