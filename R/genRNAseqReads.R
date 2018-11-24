summarizeTranscripts <- function(transcript.file){
	fout_sum <- file.path(".",paste(".Rsubread_sumfile_pid",Sys.getpid(),sep=""))
	transcript.file <- normalizePath(transcript.file, mustWork=T)
	cmd <- paste("RsummarizeTranscripts,-M,-t",transcript.file,"-o", fout_sum, sep=",")

	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_generate_random_RNAseq_reads",as.integer(n), as.character(cmd),PACKAGE="Rsubread")

	fout_sum <-paste0(fout_sum,".faSummary")
	summ <- read.delim(fout_sum, stringsAsFactors=FALSE, header=T)
	file.remove(fout_sum)
	summ
}

generateRNAseqReads <- function(transcript.file, transcript.TPM, output.prefix, out.sample.size=1000000, read.length=75, isPairedEndOutput=F, Fragment.Length.Min=100, Fragment.Length.Max=500, Fragment.Length.Mean=150, Fragment.Length.Sigma=25){
	transcript.file <- normalizePath(transcript.file, mustWork=T)
	output.prefix <- normalizePath(output.prefix, mustWork=F)
	if( !read.length %in% c(100,75) )
		stop("Error: the current version can only generate 75-bp or 100-bp reads")

	qualfile<- system.file("qualf","ref-quality-strings-20k-75bp-ERR1_59-SRR3649332.txt",package="Rsubread")
	if(read.length==100) qualfile <- system.file("qualf","ref-quality-strings-20k-100bp-ERR2_70-SRR3045231.txt",package="Rsubread")

	if( isPairedEndOutput ){
		if(Fragment.Length.Min < read.length) stop("Error: the minimum fragment length must be higher than the read length")
		if(Fragment.Length.Min < Fragment.Length.Max) stop("Error: the minimum fragment length must be equal or lower than the maximum length")
		if(Fragment.Length.Mean < Fragment.Length.Min || Fragment.Length.Mean >Fragment.Length.Max) stop("Error: the mean fragment length must be between the minimum and maximum fragment lengths")
	}
	if(out.sample.size<1) stop("Error: the output sample size must be a positive number")
	if(out.sample.size>1000*1000*1000) stop("Error: the current version cannot generate more than one billion reads/fragments in a single run")

	fin_TPMtab <- file.path(".",paste(".Rsubread_genReadTPM_pid",Sys.getpid(),sep=""))
	transcript.TPM<-as.data.frame(transcript.TPM)
	if( "TranscriptID" %in% colnames(transcript.TPM) && "TPM" %in% colnames(transcript.TPM) ){
		transcript.TPM<-data.frame(TranscriptID=transcript.TPM$TranscriptID, TPM=as.numeric(as.character(transcript.TPM$TPM)))
	}else{
		if(ncol(transcript.TPM)!=2) stop("Error: the transcript.TPM parameter must be a two-column data.frame. The first column contains the transcript names and the second column contains the TPM values")
		transcript.TPM[,2] <- as.numeric(transcript.TPM[,2])
	}
	colnames( transcript.TPM )<-c("TranscriptID","TPM")

	#print(transcript.TPM[1:6,])
	if(abs(sum( transcript.TPM$TPM ) - 1000000) > 100)
		stop(paste0("Error: the transcript.TPM parameter must be a two-column data.frame. The first column contains the transcript names and the second column contains the TPM values. The sum of all the TPM values must be 1,000,000 (but not ",sum( transcript.TPM$TPM ),"). See https://www.biostars.org/p/273537/ for the definition of TPM."))
	if(any(transcript.TPM$TPM< -0.0000000001))stop("Error: no negative TPM is allowed")

	write.table(transcript.TPM, fin_TPMtab, sep="\t", row.names=F, quote=F)

	cmd<-paste("RgenerateRNAseqReads,-t",transcript.file,"-e",fin_TPMtab,"-o",output.prefix,"-q",qualfile, "-r",out.sample.size, "-L",read.length, "-F",Fragment.Length.Mean, "-X",Fragment.Length.Max,"-N",Fragment.Length.Min,"-V",Fragment.Length.Sigma, sep=",")
	if(isPairedEndOutput) cmd<-paste0(cmd, ",-p")

	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_generate_random_RNAseq_reads",as.integer(n), as.character(cmd),PACKAGE="Rsubread")

	summ <- read.delim(paste0(output.prefix,".truthCounts"), stringsAsFactors=FALSE, header=T, comment.char="#")
	file.remove(fin_TPMtab)
	summ
}
