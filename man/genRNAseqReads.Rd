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
  out.sample.size=1000000,
  read.length=75,
  truth.in.read.names=FALSE,

  # simulating sequencing errors
  simulate.sequencing.error=TRUE,
  quality.reference=NULL,

  # consider lowly expressed transcripts
  low.transcripts=TRUE,

  # parameters for generating paired-end reads.
  paired.end=FALSE,
  insertion.length.min=100,
  insertion.length.max=500,
  insertion.length.mean=150,
  insertion.length.sigma=25,

  # manipulating transcript names
  simplify.transcript.names=FALSE
)

scanFasta(
  # the file containing the transcript database
  transcript.file,

  # manipulating transcript names
  simplify.transcript.names=FALSE
)
}
\arguments{
  \item{transcript.file}{character string giving the name of a file that contains the transcript names and sequences. The format can be FASTA or gzipped-FASTA. No duplicate sequence name is allowed in the file, and duplicate sequences will trigger warning messages.}
  \item{expression.levels}{data frame with two columns specifying the transcript IDs (first column) and relative expression levels (second column). Transcript IDs must be unique and expression levels much not be negative.}
  \item{output.prefix}{character string giving the basename of all the output files.}
  \item{out.sample.size}{numeric value giving the number of reads or read-pairs to be generated. One million by default.}
  \item{read.length}{numeric value giving the length of each read in the output. Maximum length is 250bp. Transcripts that are shorter than the read length will not be used for generating simulated reads. 75 by default.}
  \item{truth.in.read.names}{logical indicating if the true mapping location of reads or read-pairs should be encoded into the read names. \code{FALSE} by default.}
  \item{simulate.sequencing.error}{logical indicating if sequencing errors should be simulated in the output reads. If \code{TRUE}, the \code{quality.reference} parameter must be specified unless the output read length is 100-bp or 75-bp. If the output read length is 100-bp or 75-bp, the \code{quality.reference} parameter can be optionally omitted, and the function will use its inbuilt quality strings.}
  \item{quality.reference}{character string giving the name of a file that contains one or multiple sequencing quality strings in the Phred+33 format. The sequencing quality strings must have the same length as \code{read.length}.}
  \item{low.transcripts}{logical indicating if the transcripts that have expected read numbers less than one should be used for generating reads. Giving \code{FALSE} to this option can result in less reads/read-pairs generated than the wanted read/read-pair number. \code{TRUE} by default.}
  \item{paired.end}{logical indicating if paired-end reads should be generated. \code{FALSE} by default.}
  \item{insertion.length.min}{numeric value giving the minimum insertion length. The minimum insertion length must be equal to or greater than the output read length. 100 by default.}
  \item{insertion.length.max}{numeric value giving the maximum insertion length. 500 by default.}
  \item{insertion.length.mean}{numeric value giving the mean of insertion lengths. 150 by default.}
  \item{insertion.length.sigma}{numeric value giving the standard deviation of insertion lengths. The insertion lengths are drawn for a truncated normal distribution that is defined by \code{insertion.length.min}, \code{insertion.length.max}, \code{insertion.length.mean} and \code{insertion.length.sigma}. 25 by default.}
  \item{simplify.transcript.names}{logical indicating if transcript names should be simplified. If \code{TRUE}, the transcript names are truncated to the first \code{|} or space. \code{FALSE} by default.}
}
\details{
  The \code{simReads} function is for generating simulated reads from a set of transcript sequences (for example, all the transcript sequences of human) at predefined abundances. The input to this function includes a transcript file in the FASTA or gzipped-FASTA format, and a data frame that describes the wanted abundance for each transcript. The output of this function is one or two gzipped FASTQ files that contain the simulated reads or read-pairs, and a truth file giving the actually counts of reads or read-pairs generated from each transcript. The output reads usually do not have exactly the same abundances as wanted due to round-off errors and randomness. It is also remarked that not all transcripts can have reads generated from them. To have reads generated from it, a transcript has to have a length equal to or greater than the output read length, and also equal to or greater than the minimum insertion length in case of paired-end reads. The transcripts that are too short for read generation have \code{NA} read numbers in the output truth file. 

  When generating paired-end reads, the insertion lengths are drawn from a truncated normal distribution with the mean and standard deviation specified in \code{insertion.length.mean} and \code{insertion.length.sigma}; the minimum and maximum insertion lengths are specified in \code{insertion.length.min} and \code{insertion.length.max}. The insertion lengths that are greater than the actually transcript length are truncated to the transcript length.

  Substitution sequencing errors can be simulated in the reads by emulating the sequencing quality of a real High-Throughput Sequencing sample. When \code{simulate.sequencing.error} is \code{TRUE} and a set of Phred+33 encoded quality strings are provided to \code{simReads}, it randomly chooses a quality string for each output read, and substitutes the read bases with random base values at the probabilities described in the quality string. This function has inbuilt quality strings for generating 100-bp and 75-bp long reads, hence the \code{quality.reference} can be optionally omitted when \code{read.length} is 100 or 75.

  The \code{scanFasta} function checks and processes the FASTA file.
  It scans through the FASTA file that defines the transcript sequences and returns a data.frame of transcript names and sequence lengths.
  It additionally checks the transcript sequences for uniqueness.
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
%Yang Liao, Gordon K Smyth and Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. #%%\emph{Nucleic Acids Research}, 41(10):e108, 2013.
%}
\author{Wei Shi and Yang Liao}
\examples{
\dontrun{
# Scan through the fasta file to get transcript names and lengths
transcripts <- scanFasta("GENCODE-Human-transcripts.fa.gz")
nsequences <- nrow(transcripts) - sum(transcripts$Duplicate)

# Assign a random TPM value to each non-duplicated transcript sequence
TPMs <- data.frame(TranscriptID = transcripts$TranscriptID, TPM = 0)
TPMs$TPM[!transcripts$Duplicate] <- rexp(nsequences)
TPMs$TPM <- TPMs$TPM / sum(TPMs$TPM) * 1e6

# Generate actual reads.
# The output read file is my-simulated-sample_R1.fastq.gz 
# The true read counts are returned, and also in file my-simulated-sample.trueCounts
true.counts <- simReads("GENCODE-Human-transcripts.fa.gz", TPMs, "my-simulated-sample")
print(true.counts[1:10,])
}
}
