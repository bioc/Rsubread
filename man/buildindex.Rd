\name{buildindex}
\alias{buildindex}
\title{Build index for a reference genome}
\description{An index needs to be built before read mapping can be performed. This function creates a hash table for the reference genome, which can then be used by Subread and Subjunc aligners for read alignment.}
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
This function builds an index for a reference genome and the built index can then be used by Subread (\code{\link{align}}) and \code{\link{subjunc}} to map reads.
Subread is a general-purpose read aligner that can be used to map both genomic DNA sequencing (gDNA-seq) reads and RNA sequencing (RNA-seq) reads.
Subjunc is designed for detecting exon-exon junctions from using RNA-seq data.
Subread and Subjunc use a mapping paradigm called 'seed-and-vote', allowing more efficient and accurate read mapping (Liao et al. 2013).

This function generates a hash table for the reference genome, in which keys are subreads (16mers) and values are their locations in the reference genome.
Highly repetitive subreads (or uninformative subreads) are excluded from the hash table so as to reduce mapping ambiguity (mapping location of the read is directly determined from the mapping locations of subreads extracted from it).
The argument \code{TH_subread} specifies the maximal number of times a subread is allowed to occur in the reference genome to be included in the hash table.

The argument \code{memory} controls how many index segments will be generated after index building is completed.
Because only one index segment is present in the memory at any time when performing read alignments, this argument also controls the amount of computer memory (RAM) being used in read aligning.
The default value (\code{3700MB}) enables the index (ie. the hash table) built for the human or mouse genome to be partitioned into two segments, each about 4GB in size (including half of the hash table entries and half the actual reference sequences. Note that each reference base in encoded in 2 bits).
The running time of mapping 10 million 100bp reads to the human or mouse genome is typically 30 minutes.

The index for the human or mouse genome can be built into one whole segment by setting \code{memory} to \code{8000}.
This allows the maximal mapping speed to be achieved.
It takes less than 20 minutes to map 10 million reads to the human or mouse genome with this setting.

It takes typically less than one hour and half to build an index for the human or mouse genome.
The index only needs to be built once and can be re-used in the subsequent alignments.
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
