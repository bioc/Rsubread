\name{buildindex}
\alias{buildindex}
\title{Build index for a reference genome}
\description{An index needs to be built before read mapping can be performed. This function creates a hash table index for the reference genome, which can then be used by Subread and Subjunc aligners for read alignment.}
\usage{
buildindex(basename,reference,gappedIndex=FALSE,indexSplit=FALSE,memory=8000,
TH_subread=100,colorspace=FALSE)
}
\arguments{
  \item{basename}{ character string giving the basename of created index files.}
  \item{reference}{ charater string giving the name of a FASTA or gzipped FASTA file that includes sequences of all chromosomes and contigs.}
  \item{gappedIndex}{ logical indicating if a gapped index or a full index will be built. A gapped index contains 16mers (subreads) that are extracted every three bases from a reference genome, whereas a full index contains subreads extracted from every chromosomal location of a genome. The index contains a hash table, which includes sequences of subreads and their corresponding chromosomal locations. Default value of this argument is \code{FALSE} (ie. a full index is built). }
  \item{indexSplit}{ logical indicating if an index can be split into multiple blocks. The block size is determined by value of parameter \code{memory}. \code{FALSE} by default (ie. a single-block index is generated). }
  \item{memory}{ numeric value specifying the amount of memory (in megabytes) used for storing the index during read mapping. 8000 MB by default. Note that this option is ignored when \code{indexSplit} is \code{FALSE}. }
  \item{TH_subread}{ numeric value specifying the threshold for removing highly repetitive subreads (16bp mers). 100 by default. Subreads will be excluded from the index if they occur more than threshold number of times in the genome.}
  \item{colorspace}{ logical. If \code{TRUE}, a color space index will be built. Otherwise, a base space index will be built.}
}
\details{
This function generates a hash-table index for a reference genome, in which keys are subreads (16mers) and values are their chromosomal locations in the reference genome.
The built index can then be used by Subread (\code{\link{align}}) and \code{\link{subjunc}} aligners to map reads(Liao et al. 2013).

Highly repetitive subreads (or uninformative subreads) are excluded from the hash table so as to reduce mapping ambiguity.
\code{TH_subread} specifies the maximal number of times a subread is allowed to occur in the reference genome to be included in hash table.

When an index index is split into multiple blocks, only one block will be loaded into memory at any time during read mapping.
Therefore the size of memory used during mapping can be controlled.
The more memory is used in mapping, the faster the mapping speed.
Generating a one-block full index (eg. setting both \code{gappedIndex} and \code{indexSplit} to \code{FALSE}) will enable the maximum mapping speed to be achieved.

Once an index is built, it can be re-used in each mapping.

Sequences of reference genomes can be downloaded from public databases.
For instance, primary assembly of human genome GRCh38/hg38 or mouse genome GRCm38/mm10 can be downloaded from the GENCODE database via the following links:

\url{ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.primary_assembly.genome.fa.gz}

\url{ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/GRCm38.primary_assembly.genome.fa.gz}
}
\value{
No value is produced but index files are written to the current working directory.
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.  Nucleic Acids Research, 41(10):e108, 2013.

Yang Liao, Gordon K Smyth and Wei Shi. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. bioRxiv 26 July 2018. \url{https://doi.org/10.1101/377762}
}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
\examples{
# Build an index for the artifical sequence included in file 'reference.fa'
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
