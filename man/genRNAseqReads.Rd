\name{generateRNAseqReads}
\alias{generateRNAseqReads}
\alias{summarizeTranscripts}
\title{Generate simulative data from a contig database}
\description{\code{generateRNAseqReads} function can generate simulative reads from a contig database, such that each contig has a pre-specified expression level.}
\usage{
summarizeContigs(
  # the file containing the transcript database
  contig.file
)

generateRNAseqReads(
  # the contig database and the wanted expression levels
  contig.file,
  TPM,

  # the prefix of the names of the output file
  output.prefix,

  # output sample size and read length
  out.sample.size=1000000,
  read.length=75,
  simulate.sequencing.error=TRUE,
  quality.reference=NULL,

  # parameters for generating paired-end reads.
  isPairedEndOutput=F,
  Insertion.Length.Min=100,
  Insertion.Length.Max=500,
  Insertion.Length.Mean=150,
  Insertion.Length.Sigma=25
)
}
\arguments{
  \item{contig.file}{ character string giving the name of a file that contains the contig names and sequences. File format can be FASTA or gzipped-FASTA.}
  \item{TPM}{ a data frame giving the Transcripts Per Kilobase Million (TPM) values of the contig. The data frame should have two columns with the order of \code{ContigID} and \code{TPM}. The sum of the TPM values in the data frame must be one million. See https://www.biostars.org/p/273537/ for the difinition of TPM.}
  \item{output.prefix}{ character string giving the basename of output files.}
  \item{out.sample.size}{ numeric value giving the number of reads or read-pairs to be generated. One million by default.}
  \item{read.length}{ numeric value giving the length of each read in the output. Maximum length is 250bp. 75 by default.}
  \item{simulate.sequencing.error}{ logical indicating if sequencing errors should be simulated in the output reads. If \code{TRUE}, the \code{quality.reference} parameter must be specified unless the output read length is 100-bp or 75-bp. If the output read length is 100-bp or 75-bp, you can optionally omit the \code{quality.reference} parameter for simulating sequencing errors, and the function will use its inbuilt quality strings.}
  \item{quality.reference} { character string giving the name of a file that contains one or multiple sequencing quality strings in the Phred 33 format. The sequencing quality strings must have the same length as \code{read.length}.}
  \item{isPairedEndOutput}{ logical indicating if paired-end reads should be generated. FALSE by default.}
  \item{Insertion.Length.Min, Insertion.Length.Max, Insertion.Length.Mean and Insertion.Length.Sigma}{ numeric values giving the minimum, maximum, mean and standard error of the insertion lengths when paired-end reads are generated. The minimum insertion length must be greater than the output read length. Contig that are shorther than the minimum insertion length will \emph{not} be used for generating simulative reads.}
}
\details{
  The \code{generateRNAseqReads} function is for generating simulative reads from a set of contig sequences (eg all the transcript sequences of human) at wanted expression levels. The input to this function includes a contig file in the FASTA or gzipped-FASTA format, and a data frame that describes the wanted expression level for each contig. The output of this function is one or two gzipped FASTQ files that contain the simulative reads or read-pairs, and a truth file giving the actually counts of reads or read-pairs generated from each contig. The output reads usually do not have exactly the same expression levels as wanted due to round-off errors and randomness. It is also remarked that not all contigs can have reads generated from them. To have reads generated from it, a contig has to have a length equal or greater than the output read length, and also equal or greater than the minimum insertion length in case of paired-end reads. The contigs that are too short for read generation have \code{NA} read numbers in the output truth file. 

  When generating paired-end reads, the insertion lengths are drawn from a truncated normal distribution with the mean and standard error specified in \code{Insertion.Length.Mean} and \code{Insertion.Length.Sigma}; the minimum and maximum lengths are specified in \code{Insertion.Length.Min} and \code{Insertion.Length.Max}. The insertion lengths that are greater than the actually contig length are truncated to the contig length.

  Substitution sequencing errors can be simulated in the reads by emulating the sequencing quality of a real High-Throughput Sequencing sample. When a set of Phred+33 encoded quality strings are provided to \code{generateRNAseqReads}, it randomly chooses a quality string for each output read, and substitutes the read bases with random base values at the probabilities described in the quality string.

  The \code{summarizeContigs} function is provided in our package for facilitating the calculation of wanted TPM values. This function scans through a FASTA or gzipped-FASTA file that contains the contig sequences, and returns the list of contig names with their sequence lengths.
}
\value{
  The \code{summarizeContigs} function returns a \code{data.frame} object that has two columns: the name of the contigs (\code{ContigID}) and their lengths (\code{Length}). The \code{generateRNAseqReads} function returns a \code{data.frame} object that has three columns: the name of the contigs (\code{ContigID}), their lengths (\code{Length}) and the actual read counts (\code{Count}).
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.  Nucleic Acids Research, 41(10):e108, 2013.
}
\author{Wei Shi and Yang Liao}
\examples{
# Build an index for the artificial sequence included in file 'reference.fa'.
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)

# align a sample read dataset ('reads.txt') to the sample reference
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")
align.stat <- align(index="./reference_index",readfile1=reads,output_file="./Rsubread_alignment.BAM",phredOffset=64)
}
