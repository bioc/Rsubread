\name{buildindex}
\alias{buildindex}
\title{Build index for a reference genome}
\description{An index needs to be built before read mapping can be performed. This function creates a hash table index for the reference genome, which can then be used by Subread and Subjunc aligners for read alignment.}
\usage{
buildindex(basename, reference, gappedIndex = FALSE, indexSplit = FALSE,
           memory = 8000, TH_subread = 100, colorspace = FALSE)
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
The built index can then be used by Subread (\code{\link{align}}) and \code{\link{subjunc}} aligners to map reads (Liao et al. 2019; Liao et al. 2013).
Index building is an one-off operation.

Highly repetitive subreads (or uninformative subreads) are excluded from the hash table so as to reduce mapping ambiguity.
\code{TH_subread} specifies the maximal number of times a subread is allowed to occur in the reference genome to be included in hash table.

Maximum mapping speed can be achieved by building a full index for the reference genome.
By default \code{buildindex} builds a full index.
Building a gapped index will significantly reduce the memory use, at a modest cost to read mapping time.
It is recommended to use a gapped index on a personal computer due to the limited amount of computer memory available.
Memory use can be further reduced by splitting an index to multiple blocks.
The amount of memory to be used in read mapping is determined at the index building stage.

To build a full index for human/mouse genome, \code{buildindex} function requires 15GB memory.
When using a full index to map reads to human/mouse genome, \code{\link{align}} and \code{\link{subjunc}} requires 17.8GB memory.
To build a gapped index for human/mouse genome, \code{buildindex} function only requires 5.7GB memory.
When using a gapped index to map reads to human/mouse genome, \code{\link{align}} requires 8.2GB memory and \code{\link{subjunc}} requires 8.8GB memory.

Sequences of reference genomes can be downloaded from public databases.
For instance, primary assembly of human genome GRCh38/hg38 or mouse genome GRCm38/mm10 can be downloaded from the GENCODE database via the following links:

\url{ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh38.primary_assembly.genome.fa.gz}

\url{ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/GRCm38.primary_assembly.genome.fa.gz}
}
\value{
No value is produced but index files are written to the current working directory.
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi (2019).
The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads.
\emph{Nucleic Acids Research}, 47(8):e47.
\url{http://www.ncbi.nlm.nih.gov/pubmed/30783653}

Yang Liao, Gordon K Smyth and Wei Shi (2013).
The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.
\emph{Nucleic Acids Research}, 41(10):e108.
\url{http://www.ncbi.nlm.nih.gov/pubmed/23558742}
}
\author{Wei Shi and Yang Liao}
%\note{}
\seealso{\code{\link{align}}}
\examples{
# Build an index for the artifical sequence included in file 'reference.fa'
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
