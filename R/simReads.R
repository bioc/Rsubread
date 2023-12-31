.simFragments <- function(transcript.lengths, transcript.expressions=NULL, library.size=1000000L, fragment.length.min=100L, fragment.length.max=500L, fragment.length.mean=180, fragment.length.sd=40)
#  Randomly generate fragment lengths and starting positions for RNA-seq simulation.
#  Assume gamma distribution for fragment lengths.
#  Gordon Smyth
#  Created 13 March 2019. Last modified 15 March 2019.
{
# Default to equal expression levels 
  if(is.null(transcript.expressions)) transcript.expressions <- rep_len(1,length(transcript.lengths))

# To conserve memory, use integers
  transcript.lengths <- as.integer(transcript.lengths)
  fragment.length.min <- as.integer(fragment.length.min)
  fragment.length.max <- as.integer(fragment.length.max)

# First generate number of fragments from each transcript

# We will assume that any transcript with length >= fragment.length.min has a chance of being sequenced
  effective.transcript.lengths <- pmax(transcript.lengths - fragment.length.min + 1L, 0L)  
  prob <- effective.transcript.lengths * transcript.expressions
  prob <- prob / sum(prob)

# Randomly allocate total number of fragments (library size) to transcripts
  n.fragments <- drop(rmultinom(1L, size=library.size, prob=prob))

# Parameter values for gamma distribution
  alpha <- ( fragment.length.mean / fragment.length.sd )^2
  beta <- fragment.length.mean / alpha

# Probability distribution of fragment lengths
  frag.len <- seq.int(fragment.length.min,fragment.length.max)
  frag.len.p <- dgamma(frag.len,shape=alpha,scale=beta,log=TRUE)
  frag.len.p <- exp(frag.len.p - mean(frag.len.p))
  frag.len.n <- length(frag.len)

# Expand out transcript.lengths to library.size
  out <- matrix(0L,library.size,3)
  colnames(out) <- c("Transcript","FragmentLength","StartPosition")
  ntranscripts <- length(transcript.lengths)
  out[,"Transcript"] <- rep.int(seq_len(ntranscripts), n.fragments)
  TraLen <- rep.int(transcript.lengths, n.fragments)

# Randomly assign fragment lengths
  out[,"FragmentLength"] <- sample.int(frag.len.n, size=library.size, replace=TRUE, prob=frag.len.p)
  out[,"FragmentLength"] <- frag.len[out[,"FragmentLength"]]
# Rerun for short transcripts
  i <- which(out[,"FragmentLength"] > TraLen)
  if(length(i)) {
    FragLenShort <- TraLenShort <- TraLen[i]
    WhichFragMax <- max(TraLenShort) - fragment.length.min + 1L
    for (j in seq_len(WhichFragMax)) {
      k <- which(TraLenShort == frag.len[j])
      n <- length(k)
      if(n) FragLenShort[k] <- sample.int(j, size=n, replace=TRUE, prob=frag.len.p[1:j])
    }
    out[i,"FragmentLength"] <- frag.len[FragLenShort]
  }

# Generate start position
  out[,"StartPosition"] <- 1L + as.integer( runif(library.size) * (TraLen - out[,"FragmentLength"]) + 0.5 )

  list(n.fragments=n.fragments, read.positions=out)
}

simReads <- function(transcript.file, expression.levels, output.prefix, library.size=1e5, read.length=75L, truth.in.read.names=FALSE, simulate.sequencing.error=TRUE, quality.reference=NULL, paired.end=FALSE, fragment.length.min=100L, fragment.length.max=500L, fragment.length.mean=180, fragment.length.sd=40, simplify.transcript.names=FALSE)
# Simulate transcript reads and write FASTQ files
{
# Check expression.levels
  if(anyNA(expression.levels)) expression.levels[ is.na(expression.levels) ] <- 0
  if(min(expression.levels) < 0) stop("Negative expression values are not allowed.")
  if(max(expression.levels) == 0) stop("The specified transcript expression levels are all zero.")

# Check library.size
  if(!is.integer(library.size)) library.size <- as.integer(library.size+0.5)
  if(library.size < 1L) stop("At least one read should be generated.")
  if(library.size > 100L*1000L*1000L) stop("The current version cannot generate more than 100 million reads/fragments in a single run")

# Check read length
  if(!is.integer(read.length)) read.length <- as.integer(read.length+0.5)
  if(read.length < 1L || read.length > 250L) stop("The read length must be between 1 and 250.")

# Check fragment lengths
  if(fragment.length.min < read.length) stop("The minimum fragment length cannot be lower than the output read length.")
  if(fragment.length.min > fragment.length.max) stop("Maximum fragment length must be greater than or equal to the minimum fragment length.")
  if(fragment.length.mean < 1) stop("Mean fragment length must be >= 1")

# Check transcript.file
  transcript.file <- .check_and_NormPath(transcript.file, mustWork=TRUE, opt="transcript.file")
  output.prefix <- .check_and_NormPath(output.prefix, mustWork=FALSE, opt="output.prefix")
  fasta.meta <- scanFasta(transcript.file, simplify.transcript.names, quiet=TRUE)
  if(!identical(length(expression.levels),nrow(fasta.meta))) stop("Number of expression levels does not match the number of transcripts in the input fasta file.")

  if(simulate.sequencing.error){
    if(is.null(quality.reference)){
      if(read.length==75) quality.reference <- system.file("qualf","ref-quality-strings-20k-75bp-ERR1_59-SRR3649332.txt",package="Rsubread")
      if(read.length==100) quality.reference <- system.file("qualf","ref-quality-strings-20k-100bp-ERR2_70-SRR3045231.txt",package="Rsubread")
      if(is.null(quality.reference)) stop("To simulate sequencing errors in reads that are neither 100-bp nor 75-bp long, you need to provide a file containing reference quality strings of the same length as the output reads.")
    }else{
      quality.reference <- .check_and_NormPath(quality.reference,  mustWork=TRUE, opt="quality.reference")
    }
  }else{
    quality.reference <- NULL
  }

  sf <- .simFragments(fasta.meta$Length, expression.levels, library.size, fragment.length.min, fragment.length.max, fragment.length.mean, fragment.length.sd )
  C_args <- .C("R_genSimReads_at_poses", transcript.file, output.prefix, as.character(quality.reference), fasta.meta$TranscriptID, sf$read.positions[,'Transcript'], sf$read.positions[,'StartPosition'], sf$read.positions[,'FragmentLength'], as.integer(read.length), as.integer(library.size), nrow(fasta.meta), as.integer(simplify.transcript.names), as.integer(truth.in.read.names), as.integer(paired.end), PACKAGE="Rsubread")
  data.frame(fasta.meta[,1:2], NReads=sf$n.fragments)
}

scanFasta <- function(transcript.file, simplify.transcript.names=FALSE, quiet=FALSE){
  fout_sum <- file.path(".",paste(".Rsubread_sumfile_pid",Sys.getpid(),sep=""))
  transcript.file <- .check_and_NormPath(transcript.file, mustWork=T, opt="transcript.file")
  cmd <- paste("RscanFasta","--summarizeFasta","--transcriptFasta",transcript.file,"--outputPrefix", fout_sum, sep=.R_param_splitor)
  if(quiet) cmd<-paste(cmd, "--quiet", sep=.R_param_splitor)
  if(simplify.transcript.names) cmd<-paste(cmd, "--simpleTranscriptId", sep=.R_param_splitor)

  n <- length(unlist(strsplit(cmd,.R_param_splitor)))
  C_args <- .C("R_generate_random_RNAseq_reads",as.integer(n), as.character(cmd),PACKAGE="Rsubread")

  fout_sum <-paste0(fout_sum,".faSummary")
  summ <- read.delim(fout_sum, stringsAsFactors=FALSE, header=T)
  file.remove(fout_sum)
  summ
}
