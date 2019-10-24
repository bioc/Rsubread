\name{simReads}
\alias{simReads}
\alias{scanFasta}
\title{Generate simulated reads from a set of transcripts}
\description{The \code{simReads} function generates simulated reads from a set of transcripts. Each transcript has a predefined abundance in the output.}
\usage{
simReads(

  # the transcript database and the wanted abundances
  transcript.file,
  expression.levels,

  # the name of the output
  output.prefix,

  # options on the output
  library.size = 1000000,
  read.length = 75,
  truth.in.read.names = FALSE,

  # simulating sequencing errors
  simulate.sequencing.error = TRUE,
  quality.reference = NULL,

  # parameters for generating paired-end reads.
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
  \item{transcript.file}{ a character string giving the name of a file that contains the transcript names and sequences. The format can be FASTA or gzipped-FASTA. No duplicate sequence name is allowed in the file, and duplicate sequences will trigger warning messages.}
  \item{expression.levels}{ a numeric vector giving the relative expression levels of the transcripts, in the order of transcripts in \code{transcript.file}. The sum of the values must be positive, and no negative value is allowed.}
  \item{output.prefix}{ a character string giving the basename of all the output files.}
  \item{library.size}{ a numeric value giving the number of reads or read-pairs to be generated. One million by default.}
  \item{read.length}{a anumeric value giving the length of each read in the output. Maximum length is 250bp. Transcripts that are shorter than the read length will not be used for generating simulated reads. 75 by default.}
  \item{truth.in.read.names}{ logical indicating if the true mapping location of reads or read-pairs should be encoded into the read names. \code{FALSE} by default.}
  \item{simulate.sequencing.error}{ logical indicating if sequencing errors should be simulated in the output reads. If \code{TRUE}, the \code{quality.reference} parameter must be specified unless the output read length is 100-bp or 75-bp. If the output read length is 100-bp or 75-bp, the \code{quality.reference} parameter can be optionally omitted, and the function will use its inbuilt quality strings.}
  \item{quality.reference}{ a character string giving the name of a file that contains one or multiple sequencing quality strings in the Phred+33 format. The sequencing quality strings must have the same length as \code{read.length}.}
  \item{paired.end}{ logical indicating if paired-end reads should be generated. \code{FALSE} by default.}
  \item{fragment.length.min}{ a numeric value giving the minimum fragment length. The minimum fragment length must be equal to or greater than the output read length. 100 by default.}
  \item{fragment.length.max}{ a numeric value giving the maximum fragment length. 500 by default.}
  \item{fragment.length.mean}{ a numeric value giving the mean of fragment lengths. 180 by default.}
  \item{fragment.length.sd}{ a numeric value giving the standard deviation of fragment lengths. The fragment lengths are drawn for a truncated gamma distribution that is defined by \code{fragment.length.min}, \code{fragment.length.max}, \code{fragment.length.mean} and \code{fragment.length.sd}. 40 by default.}
  \item{simplify.transcript.names}{ logical indicating if transcript names should be simplified. If \code{TRUE}, the transcript names are truncated to the first \code{|} or space. \code{FALSE} by default.}
  \item{quiet}{ logical indicating if the warning message for repeated sequences should be suppressed in the \code{scanFasta} function. \code{FALSE} by default.}
}
\details{
  \code{simReads} generates simulated reads from a set of transcript sequences at specified abundances. The input includes a transcript file in FASTA or gzipped-FASTA format, and a numeric vector that describes the wanted abundance for each transcript. The output of this function is one or two gzipped FASTQ files that contain the simulated reads or read-pairs. To have reads generated from it, a transcript must have a length equal to or greater than the output read length, and also equal to or greater than the minimum fragment length in case of paired-end reads.

  When generating paired-end reads, the fragment lengths are drawn from a truncated normal distribution with the mean and standard deviation specified in \code{fragment.length.mean} and \code{fragment.length.sd}; the minimum and maximum fragment lengths are specified in \code{fragment.length.min} and \code{fragment.length.max}.

  Substitution sequencing errors can be simulated in the reads by emulating the sequencing quality of a real High-Throughput Sequencing sample. When \code{simulate.sequencing.error = TRUE} and a set of Phred+33 encoded quality strings are provided to \code{simReads}, it randomly chooses a quality string for each output read, and substitutes the read bases with random base values at the probabilities described in the quality string. This function has inbuilt quality strings for generating 100-bp and 75-bp long reads, hence the \code{quality.reference} can be optionally omitted when \code{read.length} is 100 or 75.

  The \code{scanFasta} function checks and processes the FASTA or gzipped-FASTA file. It scans through the file that defines the transcript sequences and returns a data.frame of transcript names and sequence lengths. It additionally checks the transcript sequences for uniqueness.
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
%\references{
%Yang Liao, Gordon K Smyth and Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.
%\emph{Nucleic Acids Research}, 41(10):e108, 2013.
%}
\author{Wei Shi and Yang Liao}
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
