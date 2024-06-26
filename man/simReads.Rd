\name{simReads}
\alias{simReads}
\alias{scanFasta}
\title{Generate Simulated Reads from a Set of Transcripts}
\description{The \code{simReads} function generates simulated reads from a given set of transcripts with specified expression levels.}
\usage{
simReads(

    # the transcript sequences and relative expression levels
    transcript.file,
    expression.levels,
  
    # the name of the output files
    output.prefix,
  
    # options on the output
    library.size = 1e5,
    read.length = 75,
    truth.in.read.names = FALSE,
  
    # simulating sequencing errors
    simulate.sequencing.error = TRUE,
    quality.reference = NULL,
  
    # distribution of fragment lengths
    paired.end = FALSE,
    fragment.length.min = 100,
    fragment.length.max = 500,
    fragment.length.mean = 180,
    fragment.length.sd = 40,
  
    # manipulating transcript names
    simplify.transcript.names = FALSE)

scanFasta(
  
    # the file containing the transcript database
    transcript.file,
  
    # manipulating transcript names
    simplify.transcript.names = FALSE,
  
    # miscellaneous options
    quiet = FALSE)
}
\arguments{
  \item{transcript.file}{ a character string giving the name of a file that contains the transcript names and sequences. The format can be FASTA or gzipped-FASTA. All transcript names must be unique and any duplicate sequences (two or more transcripts with the same sequence) will trigger warning messages.}
  \item{expression.levels}{ a numeric vector giving the relative expression levels of the transcripts, in the order of transcripts in \code{transcript.file}. The sum of the values must be positive, and no negative value is allowed.}
  \item{output.prefix}{ a character string giving the basename of all the output files.}
  \item{library.size}{ a numeric value giving the number of reads or read-pairs to be generated. One million by default.}
  \item{read.length}{a anumeric value giving the length of each read in the output. Maximum length is 250bp. Transcripts that are shorter than the read length will not be used for generating simulated reads. 75 by default.}
  \item{truth.in.read.names}{ logical indicating if the true mapping location of reads or read-pairs should be encoded into the read names. \code{FALSE} by default.}
  \item{simulate.sequencing.error}{ logical indicating if sequencing errors should be simulated in the output reads. If \code{TRUE}, the \code{quality.reference} parameter must be specified unless the output read length is 100-bp or 75-bp. If the output read length is 100-bp or 75-bp, the \code{quality.reference} parameter can be optionally omitted, and the function will use its inbuilt quality strings.}
  \item{quality.reference}{ a character string giving the name of a file that contains one or multiple sequencing quality strings in the Phred+33 format. The file should have Unix rather than Windows line ends and the sequencing quality strings must have the same length as \code{read.length}.}
  \item{paired.end}{ logical indicating if paired-end reads should be generated. \code{FALSE} by default.}
  \item{fragment.length.min}{ a numeric value giving the minimum fragment length. The minimum fragment length must be equal to or greater than the output read length. 100 by default.}
  \item{fragment.length.max}{ a numeric value giving the maximum fragment length. 500 by default.}
  \item{fragment.length.mean}{ a numeric value giving the mean of fragment lengths. 180 by default.}
  \item{fragment.length.sd}{ a numeric value giving the standard deviation of fragment lengths. The fragment lengths are drawn from a truncated gamma distribution that is defined by \code{fragment.length.min}, \code{fragment.length.max}, \code{fragment.length.mean} and \code{fragment.length.sd}. 40 by default.}
  \item{simplify.transcript.names}{ logical indicating if transcript names should be simplified. If \code{TRUE}, the transcript names are truncated to the first \code{|} or space. \code{FALSE} by default.}
  \item{quiet}{ logical indicating if the warning message for repeated sequences should be suppressed in the \code{scanFasta} function. \code{FALSE} by default.}
}
\details{
  \code{simReads} generates simulated reads from a set of transcript sequences at specified abundances. The input includes a transcript file in FASTA format and a numeric vector specifying the relative expression levels the transcripts. The output consists of one or two gzipped FASTQ files that contain the simulated reads or read-pairs.

  The number of fragments for each transcript are drawn from a multinomial distribution, summing to the specified library.size.
  The RNA fragment lengths are drawn from a truncated gamma distribution with mean and standard deviation specified by \code{fragment.length.mean} and \code{fragment.length.sd} respectively and lower and upper limits specified by \code{fragment.length.min} and \code{fragment.length.max}.
  The sequence reads are then taken from one or both ends of the simulated fragments to simulate short-read sequencing.

  To have reads generated from it, a transcript must have length equal to or greater than the minimum fragment length.

  Substitution sequencing errors can be simulated in the reads by emulating the sequencing quality of a real High-Throughput Sequencing sample. When \code{simulate.sequencing.error = TRUE} and a set of Phred+33 encoded quality strings are provided to \code{simReads}, a quality string is randomly chosen for each output read and read bases are substituted with random base values according to the probabilities described in the quality string. Inbuilt quality strings are provided for generating 100-bp and 75-bp long reads, hence the \code{quality.reference} can be optionally omitted when \code{read.length} is 100 or 75.

  The \code{scanFasta} function checks and processes the FASTA or gzipped-FASTA file. It scans through the file that defines the transcript sequences and returns a data.frame of transcript names and sequence lengths. It additionally checks the transcript sequences for uniqueness.
  
  Baldoni et al (2024) uses \code{simReads} to generate realistic RNA-seq data with specified levels of biological and technical variation.
}
\value{
  \code{simReads} writes a FASTQ file, or a pair of FASTQ files if \code{paired.end}, containing the simulated reads.
  It also returns a \code{data.frame} with three columns:
  \tabular{lll}{
    \code{TranscriptID} \tab character \tab transcript name\cr
    \code{Length} \tab integer \tab length of transcript sequence\cr
    \code{Count} \tab integer \tab simulated read count
  }

  \code{scanFasta} returns a \code{data.frame} with six columns:
  \tabular{lll}{
    \code{TranscriptID} \tab character \tab transcript name\cr
    \code{Length} \tab integer \tab length of transcript sequence\cr
    \code{MD5} \tab character \tab MD5 digest of the transcript sequence\cr
    \code{Unique} \tab logical \tab is this transcript's sequence unique in the FASTA file?\cr
    \code{Occurrence} \tab integer \tab number of times this transcript's sequence was observed\cr
    \code{Duplicate} \tab logical \tab this transcript's sequence is a duplicate of a previous sequence.
  }
  Note that selecting transcripts with \code{Duplicate == FALSE} will ensure unique sequences, i.e., any sequence that was observed multiple times in the FASTQ file will be only included only once in the selection.
}
\references{
Baldoni PL, Chen Y, Hediyeh-zadeh S, Liao Y, Dong X, Ritchie ME, Shi W, Smyth GK (2024).
Dividing out quantification uncertainty allows efficient assessment of differential transcript expression with edgeR.
\emph{Nucleic Acids Research} 52(3), e13.
\doi{10.1093/nar/gkad1167}.
}
\author{Yang Liao, Gordon K Smyth and Wei Shi}
\examples{
\dontrun{
# Scan through the fasta file to get transcript names and lengths
transcripts <- scanFasta("GENCODE-Human-transcripts.fa.gz")
nsequences <- nrow(transcripts) - sum(transcripts$Duplicate)

# Assign a random TPM value to each non-duplicated transcript sequence
TPMs <- rep(0, nrow(transcripts))
TPMs[!transcripts$Duplicate] <- rexp(nsequences)

# Generate actual reads.
# The output read file is my-simulated-sample_R1.fastq.gz 
# The true read counts are returned.
true.counts <- simReads("GENCODE-Human-transcripts.fa.gz", TPMs, "my-simulated-sample")
print(true.counts[1:10,])
}
}
