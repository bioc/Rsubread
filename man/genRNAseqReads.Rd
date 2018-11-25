\name{generateRNAseqReads}
\alias{generateRNAseqReads}
\alias{summarizeTranscripts}
\title{Generate simulation data from a transcript database}
\description{\code{generateRNAseqReads} function can generate simulation reads from a transcript database, such that each transcript has a pre-specified expression level.}
\usage{
summarizeTranscripts(
  # the file containing the transcript database
  transcript.file
)

generateRNAseqReads(
  # the transcript database and the wanted expression levels
  transcript.file,
  transcript.TPM,

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
  \item{transcript.file}{ character string giving the name of a file that contains the transcript names and sequences. File format can be FASTA or gzipped-FASTA.}
  \item{transcript.TPM}{ a data frame giving the Transcripts Per Kilobase Million (TPM) values of the transcripts. The data frame should have two columns with the order of TranscriptID and TPM. The sum of the TPM values in the data frame must be one million. TPMs are essencially the ratios of the RPKM values in a sample. See https://www.biostars.org/p/273537/ for more details about TPM and RPKM.}
  \item{output.prefix}{ character string giving the basename of output files.}
  \item{out.sample.size}{ numeric value giving the number of reads or read-pairs to be generated. One million by default.}
  \item{read.length}{ numeric value giving the length of each read in the output. Maximum length is 250bp. 75 by default.}
  \item{simulate.sequencing.error}{ logical indicating if sequencing errors should be simulated in the output reads. If \code{TRUE}, the \code{quality.reference} parameter must be specified unless the output read length is 100-bp or 75-bp. If the output read length is 100-bp or 75-bp, you can optionally omit the \code{quality.reference} parameter for simulating sequencing errors, and the function will use its inbuilt quality strings.}
  \item{quality.reference} { character string giving the name of a file that contains one or multiple sequencing quality strings in the Phred 33 format. The sequencing quality strings must have the same length as \code{read.length}.}
  \item{isPairedEndOutput}{ logical indicating if paired-end reads should be generated. FALSE by default.}
  \item{Insertion.Length.Min, Insertion.Length.Max, Insertion.Length.Mean and Insertion.Length.Sigma}{ numeric values giving the minimum, maximum, mean and standard error of the insertion lengths when paired-end reads are generated. The minimum insertion length must be higher than the output read length. Transcripts that are shorther than the minimum insertion length will \emph{not} be used for generating simulation reads.}
}
\details{
  The \code{summarizeTranscripts} and \code{generateRNAseqReads} functions are useful in generating simulative RNA-seq reads from transcript sequences. 
  The \code{summarizeTranscripts} function scans a FASTA or gzipped-FASTA file that contains transcript sequences, and reports the list of transcripts with their lengths. From the output of \code{summarizeTranscripts}, the user can create a data frame specifying the wanted Transcripts Per Kilobase Million (TPM) value for each transcript, and run \code{generateRNAseqReads} to generate reads that express the transcripts at the wanted TPM levels. 
  Substitution sequencing errors can be simulated in the reads by following real samples.
 when reference sequencing quality strings are provided. 
}
\value{
  The \code{summarizeTranscripts} function returns a \code{data.frame} object that includes
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
