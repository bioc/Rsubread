\name{generateSimulatedReads}
\alias{generateSimulatedReads}
\alias{summarizeContigs}
\title{Generate simulated reads from a set of contigs}
\description{The \code{generateSimulatedReads} function generates simulated reads from a set of contigs. Each contig has a predefined abundance in the output.}
\usage{
generateSimulatedReads(
  # the contig database and the wanted abundances
  contig.file,
  TPM,

  # the name of the output
  output.prefix,

  # options on the output
  out.sample.size=1000000,
  read.length=75,
  truth.in.read.names=FALSE,

  # simulating sequencing errors
  simulate.sequencing.error=TRUE,
  quality.reference=NULL,

  # parameters for generating paired-end reads.
  isPairedEndOutput=FALSE,
  Insertion.Length.Min=100,
  Insertion.Length.Max=500,
  Insertion.Length.Mean=150,
  Insertion.Length.Sigma=25,

  # manipulating contig names
  simplify.contig.names=FALSE
)

summarizeContigs(
  # the file containing the transcript database
  contig.file,

  # manipulating contig names
  simplify.contig.names=FALSE
)
}
\arguments{
  \item{contig.file}{ character string giving the name of a file that contains the contig names and sequences. The format can be FASTA or gzipped-FASTA.}
  \item{TPM}{ a data frame giving the Transcripts Per Kilobase Million (TPM) values of the contig. The data frame should have only two columns: \code{ContigID} and \code{TPM}. The sum of the TPM values in the data frame must be one million, and no negative TPM value is allowed. See https://www.biostars.org/p/273537/ for the difinition of TPM.}
  \item{output.prefix}{ character string giving the basename of all the output files.}
  \item{out.sample.size}{ numeric value giving the number of reads or read-pairs to be generated. One million by default.}
  \item{read.length}{ numeric value giving the length of each read in the output. Maximum length is 250bp. Contigs that are shorter than the read length will not be used for generating simulated reads. 75 by default.}
  \item{truth.in.read.names}{ logical indicating if the true mapping location of reads or read-pairs should be encoded into the read names. \code{FALSE} by default.}
  \item{simulate.sequencing.error}{ logical indicating if sequencing errors should be simulated in the output reads. If \code{TRUE}, the \code{quality.reference} parameter must be specified unless the output read length is 100-bp or 75-bp. If the output read length is 100-bp or 75-bp, the \code{quality.reference} parameter can be optionally omitted, and the function will use its inbuilt quality strings.}
  \item{quality.reference}{ character string giving the name of a file that contains one or multiple sequencing quality strings in the Phred+33 format. The sequencing quality strings must have the same length as \code{read.length}.}
  \item{isPairedEndOutput}{ logical indicating if paired-end reads should be generated. \code{FALSE} by default.}
  \item{Insertion.Length.Min}{ numeric values giving the minimum insertion length. The minimum insertion length must be equal or greater than the output read length. 100 by default.}
  \item{Insertion.Length.Max}{ numeric values giving the maximum insertion length. 500 by default.}
  \item{Insertion.Length.Mean}{ numeric values giving the mean of insertion lengths. 150 by default.}
  \item{Insertion.Length.Sigma}{ numeric values giving the standard deviation of insertion lengths. The insertion lengths are drawn for a truncated normal distribution that is defined by \code{Insertion.Length.Min}, \code{Insertion.Length.Max}, \code{Insertion.Length.Mean} and \code{Insertion.Length.Sigma}. 25 by default.}
  \item{simplify.contig.names}{ logical indicating if contig names should be simplified. If \code{TRUE}, the contig names are truncated to the first \code{|} or space. \code{FALSE} by default.}
}
\details{
  The \code{generateSimulatedReads} function is for generating simulated reads from a set of contig sequences (for example, all the transcript sequences of human) at predefined abundances. The input to this function includes a contig file in the FASTA or gzipped-FASTA format, and a data frame that describes the wanted abundance for each contig. The output of this function is one or two gzipped FASTQ files that contain the simulated reads or read-pairs, and a truth file giving the actually counts of reads or read-pairs generated from each contig. The output reads usually do not have exactly the same abundances as wanted due to round-off errors and randomness. It is also remarked that not all contigs can have reads generated from them. To have reads generated from it, a contig has to have a length equal or greater than the output read length, and also equal or greater than the minimum insertion length in case of paired-end reads. The contigs that are too short for read generation have \code{NA} read numbers in the output truth file. 

  When generating paired-end reads, the insertion lengths are drawn from a truncated normal distribution with the mean and standard deviation specified in \code{Insertion.Length.Mean} and \code{Insertion.Length.Sigma}; the minimum and maximum insertion lengths are specified in \code{Insertion.Length.Min} and \code{Insertion.Length.Max}. The insertion lengths that are greater than the actually contig length are truncated to the contig length.

  Substitution sequencing errors can be simulated in the reads by emulating the sequencing quality of a real High-Throughput Sequencing sample. When \code{simulate.sequencing.error} is \code{TRUE} and a set of Phred+33 encoded quality strings are provided to \code{generateSimulatedReads}, it randomly chooses a quality string for each output read, and substitutes the read bases with random base values at the probabilities described in the quality string. This function has inbuilt quality strings for generating 100-bp and 75-bp long reads, hence the \code{quality.reference} can be optionally omitted when \code{read.length} is 100 or 75.

  The \code{summarizeContigs} function is provided for facilitating the calculation of wanted TPM values. This function scans through a FASTA or gzipped-FASTA file that contains the contig sequences, and returns the list of contig names with their sequence lengths.
}
\value{
  The \code{summarizeContigs} function returns a \code{data.frame} object that has two columns: the name of the contigs (\code{ContigID}) and their lengths (\code{Length}). The \code{generateSimulatedReads} function returns a \code{data.frame} object that has three columns: the name of the contigs (\code{ContigID}), their lengths (\code{Length}) and the actual read counts (\code{Count}).
}
\references{
Yang Liao, Gordon K Smyth and Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.  Nucleic Acids Research, 41(10):e108, 2013.
}
\author{Wei Shi and Yang Liao}
\examples{
\dontrun{
# Build an index for the artificial sequence included in file 'reference.fa'.
library(Rsubread)

contigs <- summarizeContigs("GENCODE-Human-transcripts.fa.gz")

# Assigned a random TPM value to each contig. 
TPMs <- data.frame(ContigID=contigs$ContigID, TPM=rnorm(nrow(contigs), mean=2, sd=1))
TPMs$TPM <- ifelse(TPMs$TPM < 0, 0, TPMs$TPM) # no negative TPM is allowed
TPMs$TPM <- TPMs$TPM / sum(TPMs$TPM)

true.counts <- generateSimulatedReads("GENCODE-Human-transcripts.fa.gz", TPMs, "my-simulated-sample")
print(true.counts[1:10,])
}
}
