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
  Fragment.Length.Min=100,
  Fragment.Length.Max=500,
  Fragment.Length.Mean=150,
  Fragment.Length.Sigma=25
)
}
\arguments{
  \item{transcript.file}{ character string giving the name of a file that contains the transcript names and sequences. File format can be FASTA or gzipped-FASTA.}
  \item{transcript.TPM}{ a data frame giving the Transcripts Per Kilobase Million (TPM) values of the transcripts. The data frame should have two columns: TranscriptID and TPM. The sum of the TPM values in the data frame must be one million.}
  \item{output.prefix}{ character string giving the basename of output files.}
  \item{out.sample.size}{ numeric value giving the number of reads or read-pairs to be generated. One million by default.}
  \item{read.length}{ numeric value giving the length of each read in the output. The current version can only generate 75- or 100-bp long reads. 75 by default.}
  \item{simulate.sequencing.error}{ logical indicating if sequencing errors should be simulated in the output reads. If \code{TRUE}, the }
  \item{quality.reference} { character string giving the name of a file that contains one or multiple sequencing quality strings in the Phred 33 format. }
}
