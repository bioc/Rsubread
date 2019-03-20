simFragments <- function(transcript.lengths, transcript.expressions=NULL, library.size=1e6, fragment.length.min=100L, fragment.length.max=500L, fragment.length.mean=180, fragment.length.sd=40)
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
  effective.transcript.lengths <- pmax(transcript.lengths - fragment.length.min + 1, 0)  
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

  out
}

simReads <- function(transcript.file, expression.levels, output.prefix, library.size=1000000, read.length=75, truth.in.read.names=FALSE, simulate.sequencing.error=TRUE, quality.reference=NULL, paired.end=FALSE, fragment.length.min=100L, fragment.length.max=500L, fragment.length.mean=180, fragment.length.sd=40, simplify.transcript.names=FALSE){
  if(!paired.end) stop("current temp version does not support SE")
  if(library.size>100*1000*1000) stop("The current version cannot generate more than 100 million reads/fragments in a single run")
  transcript.file <- .check_and_NormPath(transcript.file, mustWork=TRUE, opt="transcript.file")
  output.prefix <- .check_and_NormPath(output.prefix, mustWork=FALSE, opt="output.prefix")

  fasta.meta <- scanFasta(transcript.file, simplify.transcript.names)
  expression.levels.MetaOrder <- expression.levels[ match( fasta.meta$TranscriptID, expression.levels[,1] ),2 ]

  read.positions <- simFragments(fasta.meta$Length, expression.levels.MetaOrder, library.size, fragment.length.min, fragment.length.max, fragment.length.mean, fragment.length.sd )
  if(simulate.sequencing.error){
    if(is.null(quality.reference)){
      if(read.length==75) quality.reference<- system.file("qualf","ref-quality-strings-20k-75bp-ERR1_59-SRR3649332.txt",package="Rsubread")
      if(read.length==100) quality.reference <- system.file("qualf","ref-quality-strings-20k-100bp-ERR2_70-SRR3045231.txt",package="Rsubread")
      if(is.null(quality.reference)) stop("To simulate sequencing errors in reads that are neither 100-bp nor 75-bp long, you need to provide a file containing reference quality strings of the same length as the output reads.")
    }else{
      quality.reference <- .check_and_NormPath(quality.reference,  mustWork=TRUE, opt="quality.reference")
    }
  }else{
    quality.reference <- NULL
  }
  
  C_args <- .C("R_genSimReads_at_poses", transcript.file, output.prefix, as.character(quality.reference), fasta.meta$TranscriptID , read.positions[,'Transcript'] , read.positions[,'StartPosition'],  read.positions[,'FragmentLength'], as.integer(read.length), as.integer(library.size), nrow(fasta.meta), as.integer(simplify.transcript.names), as.integer(truth.in.read.names), as.integer(paired.end), PACKAGE="Rsubread")
  rets <- table(fasta.meta$TranscriptID[read.positions[,"Transcript"]])
  rets <- data.frame(fasta.meta[,1:2], Count=as.vector(rets)[match( fasta.meta[,1] , names(rets))])
  rets[is.na(rets[,'Count']) ,'Count']<-0
  write.table(rets, paste0(output.prefix,".truthCounts"), quote=FALSE, sep="\t", row.names=FALSE)
  rets
}

scanFasta <- function(transcript.file, simplify.transcript.names=FALSE){
	fout_sum <- file.path(".",paste(".Rsubread_sumfile_pid",Sys.getpid(),sep=""))
	transcript.file <- .check_and_NormPath(transcript.file, mustWork=T, opt="transcript.file")
	cmd <- paste("RscanFasta,--summarizeFasta,--transcriptFasta",transcript.file,"--outputPrefix", fout_sum, sep=",")
    if(simplify.transcript.names) cmd<-paste(cmd, "--simpleTranscriptId", sep=",")

	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_generate_random_RNAseq_reads",as.integer(n), as.character(cmd),PACKAGE="Rsubread")

	fout_sum <-paste0(fout_sum,".faSummary")
	summ <- read.delim(fout_sum, stringsAsFactors=FALSE, header=T)
	file.remove(fout_sum)
	summ
}

fix_simReads <- function(transcript.file, expression.levels, output.prefix, out.sample.size=1000000, read.length=75, truth.in.read.names=FALSE, simulate.sequencing.error=TRUE, quality.reference=NULL, low.transcripts=TRUE, iterative.find.N=FALSE, paired.end=FALSE, fragment.length.min=100, fragment.length.max=500, fragment.length.mean=150, fragment.length.sigma=25, simplify.transcript.names=FALSE,gen.reads=TRUE){
	transcript.file <- .check_and_NormPath(transcript.file, mustWork=TRUE, opt="transcript.file")
	output.prefix <- .check_and_NormPath(output.prefix, mustWork=FALSE, opt="output.prefix")
	if(read.length > 250) stop("The current version can generate reads at most 250bp long.")
	if(read.length <1) stop("Read length must be a positive integer")
	if(iterative.find.N && low.transcripts) stop("When using the iterative algorithm to find the best sample size, the lowly expressed transcripts have to be excluded (ie low.transcripts=FALSE)")

	qualfile <- NULL
	if(simulate.sequencing.error){
		if(is.null(quality.reference)){
			if(read.length==75) qualfile<- system.file("qualf","ref-quality-strings-20k-75bp-ERR1_59-SRR3649332.txt",package="Rsubread")
			if(read.length==100) qualfile <- system.file("qualf","ref-quality-strings-20k-100bp-ERR2_70-SRR3045231.txt",package="Rsubread")
			if(is.null(qualfile)) stop("When you want to simulate sequencing errors in reads that are neither 100-bp nor 75-bp long, you need to provide a file containing reference quality strings that have the length as the output reads.")
		}else{
			qualfile <- .check_and_NormPath(quality.reference, mustWork=TRUE, opt="qualfile")
		}
	}

	if( paired.end ){
		if(fragment.length.min < read.length) stop("The minimum fragment length must be higher than the read length")
		if(fragment.length.min > fragment.length.max) stop("The minimum fragment length must be equal or lower than the maximum length")
		if(fragment.length.mean < fragment.length.min || fragment.length.mean >fragment.length.max) stop("Error: the mean fragment length must be between the minimum and maximum fragment lengths")
	}
	if(out.sample.size<1) stop("The output sample size must be a positive integer")
	if(out.sample.size>1000*1000*1000) stop("The current version cannot generate more than one billion reads/fragments in a single run")

	fin_TPMtab <- file.path(".",paste(".Rsubread_genReadTPM_pid",Sys.getpid(),sep=""))
	transcript.TPM <- as.data.frame(expression.levels)
	if( "TranscriptID" %in% colnames(transcript.TPM) && "TPM" %in% colnames(transcript.TPM) ){
		transcript.TPM <- data.frame(TranscriptID=transcript.TPM$TranscriptID, TPM=as.numeric(as.character(transcript.TPM$TPM)))
	}else{
		if(ncol(transcript.TPM)!=2) stop("The 'expression.levels' parameter must be a two-column data.frame. The first column contains the transcript names and the second column contains the relative expression levels")
		transcript.TPM[,2] <- as.numeric(as.character(transcript.TPM[,2]))
	}
	colnames( transcript.TPM ) <- c("TranscriptID","TPM")

	# we allow any ratio values for the expression values.
	# however, the C code only allows the sum of one million.
	# hence, there is a conversion in R.

	rangeTPM <- range(transcript.TPM$TPM)
	if( is.na(rangeTPM[1]) ) stop("NA expression levels not allowed")
	if( rangeTPM[1] < 0 ) stop("Negative expression levels not allowed")
	if( rangeTPM[2] <= 0 ) stop("At least some expression levels must be positive")
	transcript.TPM$TPM <- transcript.TPM$TPM / sum(transcript.TPM$TPM) * 1e6
	write.table(transcript.TPM, file=fin_TPMtab, sep="\t", row.names=FALSE, quote=FALSE)

	cmd <- paste("RgenerateRNAseqReads,--transcriptFasta",transcript.file,"--expressionLevels",fin_TPMtab,"--outputPrefix",output.prefix, "--totalReads",sprintf("%d",out.sample.size), "--readLen",read.length, sep=",")
	if(paired.end) cmd <- paste(cmd, "--pairedEnd,--fragmentLenMean",fragment.length.mean, "--fragmentLenMax",fragment.length.max,"--fragmentLenMin",fragment.length.min,"--fragmentLenSigma",fragment.length.sigma, sep=",")
	if(simplify.transcript.names) cmd <- paste(cmd, "--simpleTranscriptId", sep=",")
	if(truth.in.read.names) cmd <- paste(cmd, "--truthInReadNames", sep=",")

	if(iterative.find.N) cmd <- paste(cmd, "--floorStrategy","ITERATIVE", sep=",")
	if(!(low.transcripts || iterative.find.N )) cmd <- paste(cmd, "--floorStrategy","FLOOR", sep=",")

	if(!gen.reads) cmd <- paste(cmd, "--noActualReads" ,sep=",")
	if(!is.null(qualfile)) cmd <- paste(cmd, "--qualityRefFile", qualfile, sep=",")

	#print(substr(cmd,1,2000))
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_generate_random_RNAseq_reads",as.integer(n), as.character(cmd),PACKAGE="Rsubread")

	summ <- read.delim(paste0(output.prefix,".truthCounts"), stringsAsFactors=FALSE, header=T, comment.char="#", colClasses=c("character", "numeric", "numeric"))
	file.remove(fin_TPMtab)
	summ
}
