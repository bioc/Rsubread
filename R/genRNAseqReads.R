scanFasta <- function(transcript.file, simplify.transcript.names=FALSE){
	fout_sum <- file.path(".",paste(".Rsubread_sumfile_pid",Sys.getpid(),sep=""))
	transcript.file <- normalizePath(transcript.file, mustWork=T)
	cmd <- paste("RscanFasta,--summarizeFasta,--transcriptFasta",transcript.file,"--outputPrefix", fout_sum, sep=",")
    if(simplify.transcript.names) cmd<-paste(cmd, "--simpleTranscriptId", sep=",")

	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_generate_random_RNAseq_reads",as.integer(n), as.character(cmd),PACKAGE="Rsubread")

	fout_sum <-paste0(fout_sum,".faSummary")
	summ <- read.delim(fout_sum, stringsAsFactors=FALSE, header=T)
	file.remove(fout_sum)
	summ
}

simReads <- function(transcript.file, expression.levels, output.prefix, out.sample.size=1000000, read.length=75, truth.in.read.names=FALSE, simulate.sequencing.error=TRUE, quality.reference=NULL, low.transcripts=TRUE, isPairedEndOutput=FALSE, Insertion.Length.Min=100, Insertion.Length.Max=500, Insertion.Length.Mean=150, Insertion.Length.Sigma=25, simplify.transcript.names=FALSE){
	transcript.file <- normalizePath(transcript.file, mustWork=T)
	output.prefix <- normalizePath(output.prefix, mustWork=F)
	if(read.length > 250)stop("Error: the current version can generate reads at most 250bp long.")
	if(read.length <1) stop("Error: the read length must be positive.")

	qualfile <- NULL
	if(simulate.sequencing.error){
		if(is.null(quality.reference)){
			if(read.length==75) qualfile<- system.file("qualf","ref-quality-strings-20k-75bp-ERR1_59-SRR3649332.txt",package="Rsubread")
			if(read.length==100) qualfile <- system.file("qualf","ref-quality-strings-20k-100bp-ERR2_70-SRR3045231.txt",package="Rsubread")
			if(is.null(qualfile)) stop("When you want to simulate sequencing errors in reads that are neither 100-bp nor 75-bp long, you need to provide a file containing reference quality strings that have the length as the output reads.")
		}else{
			qualfile <- normalizePath(quality.reference, mustWork=T)
		}
	}

	if( isPairedEndOutput ){
		if(Insertion.Length.Min < read.length) stop("Error: the minimum insertion length must be higher than the read length")
		if(Insertion.Length.Min > Insertion.Length.Max) stop("Error: the minimum insertion length must be equal or lower than the maximum length")
		if(Insertion.Length.Mean < Insertion.Length.Min || Insertion.Length.Mean >Insertion.Length.Max) stop("Error: the mean insertion length must be between the minimum and maximum insertion lengths")
	}
	if(out.sample.size<1) stop("Error: the output sample size must be a positive number")
	if(out.sample.size>1000*1000*1000) stop("Error: the current version cannot generate more than one billion reads/insertions in a single run")

	fin_TPMtab <- file.path(".",paste(".Rsubread_genReadTPM_pid",Sys.getpid(),sep=""))
	transcript.TPM<-as.data.frame(expression.levels)
	if( "TranscriptID" %in% colnames(transcript.TPM) && "TPM" %in% colnames(transcript.TPM) ){
		transcript.TPM<-data.frame(TranscriptID=transcript.TPM$TranscriptID, TPM=as.numeric(as.character(transcript.TPM$TPM)))
	}else{
		if(ncol(transcript.TPM)!=2) stop("Error: the TPM parameter must be a two-column data.frame. The first column contains the transcript names and the second column contains the TPM values")
		transcript.TPM[,2] <- as.numeric(as.character(transcript.TPM[,2]))
	}
	colnames( transcript.TPM )<-c("TranscriptID","TPM")

	
	if(F){ # do not check one million.
		if(abs(sum( transcript.TPM$TPM ) - 1000000) > 100)
			stop(paste0("Error: the TPM parameter must be a two-column data.frame. The first column contains the transcript names and the second column contains the TPM values. The sum of all the TPM values must be 1,000,000 (but not ",sum( transcript.TPM$TPM ),"). See https://www.biostars.org/p/273537/ for the definition of TPM."))
	}

	# we allow any ratio values for the expression values.
	# however, the C code only allows the sum of one million.
	# hence, there is a conversion in R.

	if(any(transcript.TPM$TPM< -0.0000000001))stop("Error: no negative TPM is allowed")
	if( sum(transcript.TPM$TPM) <= 0 ) stop("Error: the sum of the expression levels is zero or negative.")
	transcript.TPM$TPM <- transcript.TPM$TPM / sum(transcript.TPM$TPM) * 1000000.
	write.table(transcript.TPM, fin_TPMtab, sep="\t", row.names=F, quote=F)

	cmd<-paste("RgenerateRNAseqReads,--transcriptFasta",transcript.file,"--expressionLevels",fin_TPMtab,"--outputPrefix",output.prefix, "--totalReads",sprintf("%d",out.sample.size), "--readLen",read.length, sep=",")
	if(isPairedEndOutput) cmd<-paste(cmd, "--pairedEnd,--insertionLenMean",Insertion.Length.Mean, "--insertionLenMax",Insertion.Length.Max,"--insertionLenMin",Insertion.Length.Min,"--insertionLenSigma",Insertion.Length.Sigma, sep=",")
	if(simplify.transcript.names) cmd<-paste(cmd, "--simpleTranscriptId", sep=",")
	if(truth.in.read.names) cmd<-paste(cmd, "--truthInReadNames", sep=",")
	if(!low.transcripts) cmd <- paste(cmd, "--noLowTranscripts", sep=",")
	if(!is.null(qualfile)) cmd<-paste(cmd, "--qualityRefFile", qualfile, sep=",")

	#print(substr(cmd,1,2000))
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_generate_random_RNAseq_reads",as.integer(n), as.character(cmd),PACKAGE="Rsubread")

	summ <- read.delim(paste0(output.prefix,".truthCounts"), stringsAsFactors=FALSE, header=T, comment.char="#", colClasses=c("character", "numeric", "numeric"))
	file.remove(fin_TPMtab)
	summ
}
