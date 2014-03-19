\name{buildindex}
\alias{buildindex}
\title{Build index for a reference genome}
\description{An index needs to be built before read mapping can be performed. This function creates a hash table for the reference genome, which can then be used by Subread and Subjunc aligners for read alignment.}
\usage{
buildindex(basename,reference,gappedIndex=TRUE,indexSplit=TRUE,memory=8000,TH_subread=24,colorspace=FALSE)
}
\arguments{
  \item{basename}{ character string giving the basename of created index files.}
  \item{reference}{ charater string giving the name of the file containing all the refernece sequences.}
  \item{gappedIndex}{ logical. If \code{FALSE}, 16mers (subreads) will be extracted from every chromosomal location of a reference genome and then they will be used to build a hash table index. By default(\code{TRUE}), subreads are extracted in every three bases from the genome.}
  \item{indexSplit}{ logical. If \code{TRUE}, the built index is allowed to be splitted into multiple segments. The number of such segments is determined by \code{memory} value, genome size and permitting of gaps between subreads(\code{gappedIndex}). If \code{indexSplit} is set to \code{FALSE}, a single-segment index (no splitting) will be generated regardless of what value is chosen for \code{memory}.}
  \item{memory}{ numeric value specifying the amount of memory to be requested in megabytes. 8000 MB by default.}
  \item{TH_subread}{ numeric value specifying the threshold for removing highly repetitive subreads (16 mers). 24 by default. Subreads will be excluded from the index if they occur more than threshold number of times in the genome.}
  \item{colorspace}{ logical. If \code{TRUE}, a color space index will be built. Otherwise, a base space index will be built.}
  }
\details{
This function generates a hash table (an index) for a reference genome, in which keys are subreads (16mers) and values are their chromosomal locations in the reference genome.
By default, subreads will be extracted in every three bases from a reference genome.
However, if \code{gappedIndex} is set to \code{FALSE}, then subreads will be extracted from every chromosomal location of genome for index building.
The built index can then be used by Subread (\code{\link{align}}) and \code{\link{subjunc}} aligners to map reads(Liao et al. 2013).

Highly repetitive subreads (or uninformative subreads) are excluded from the hash table so as to reduce mapping ambiguity.
\code{TH_subread} specifies the maximal number of times a subread is allowed to occur in the reference genome to be included in hash table.

The built index might be splitted into multiple segments if its size is greater than \code{memory} value.
The number of such segments is dependent on \code{memory} value, size of reference genome and whether gaps are allowed between subreads extracted from genome.
Only one segment is loaded into memory at any time when read alignment is being carried out.
The larger the \code{memory} value, the faster the read mapping will be.
If \code{indexSplit} is set to \code{FALSE}, the index will not be splitted and this will enable maximum mapping speed to be achieved.

The index needs to be built only once and it can then be re-used in the subsequent alignments.
}
\value{
No value is produced but index files are written to the current working directory.
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.  Nucleic Acids Research, 41(10):e108, 2013.
}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
\examples{
# Build an index for the artifical sequence included in file 'reference.fa'
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
