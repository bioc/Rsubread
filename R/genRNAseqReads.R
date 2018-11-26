summarizeContigs <- function(contig.file, simplify.contig.names=F){
	fout_sum <- file.path(".",paste(".Rsubread_sumfile_pid",Sys.getpid(),sep=""))
	contig.file <- normalizePath(contig.file, mustWork=T)
	cmd <- paste("RsummarizeContigs,--summarizeFasta,--contigFasta",contig.file,"--outputPrefix", fout_sum, sep=",")
    if(simplify.contig.names) cmd<-paste(cmd, "--simpleContigId", sep=",")

	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_generate_random_RNAseq_reads",as.integer(n), as.character(cmd),PACKAGE="Rsubread")

	fout_sum <-paste0(fout_sum,".faSummary")
	summ <- read.delim(fout_sum, stringsAsFactors=FALSE, header=T)
	file.remove(fout_sum)
	summ
}

generateSimulatedReads <- function(contig.file, TPM, output.prefix, out.sample.size=1000000, read.length=75, truth.in.read.names=F, simulate.sequencing.error=T, quality.reference=NULL, isPairedEndOutput=F, Insertion.Length.Min=100, Insertion.Length.Max=500, Insertion.Length.Mean=150, Insertion.Length.Sigma=25, simplify.contig.names=F){
	contig.file <- normalizePath(contig.file, mustWork=T)
	output.prefix <- normalizePath(output.prefix, mustWork=F)
	if( !read.length %in% c(100,75) )
		stop("Error: the current version can only generate 75-bp or 100-bp reads")

	qualfile <- NULL
	if(simulate.sequencing.error){
		if(quality.reference ==NULL){
			if(read.length==75) qualfile<- system.file("qualf","ref-quality-strings-20k-75bp-ERR1_59-SRR3649332.txt",package="Rsubread")
			if(read.length==100) qualfile <- system.file("qualf","ref-quality-strings-20k-100bp-ERR2_70-SRR3045231.txt",package="Rsubread")
			if(qualfile == NULL) stop("When you want to simulate sequencing error in the reads that are neither 100-bp nor 75-bp long, you need to provide a file containing reference quality strings.")
		}else{
			qualfile <- quality.reference
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
	contig.TPM<-as.data.frame(TPM)
	if( "ContigID" %in% colnames(contig.TPM) && "TPM" %in% colnames(contig.TPM) ){
		contig.TPM<-data.frame(ContigID=contig.TPM$ContigID, TPM=as.numeric(as.character(contig.TPM$TPM)))
	}else{
		if(ncol(contig.TPM)!=2) stop("Error: the TPM parameter must be a two-column data.frame. The first column contains the contig names and the second column contains the TPM values")
		contig.TPM[,2] <- as.numeric(as.character(contig.TPM[,2]))
	}
	colnames( contig.TPM )<-c("ContigID","TPM")

	#print(contig.TPM[1:6,])
	if(abs(sum( contig.TPM$TPM ) - 1000000) > 100)
		stop(paste0("Error: the TPM parameter must be a two-column data.frame. The first column contains the contig names and the second column contains the TPM values. The sum of all the TPM values must be 1,000,000 (but not ",sum( contig.TPM$TPM ),"). See https://www.biostars.org/p/273537/ for the definition of TPM."))
	if(any(contig.TPM$TPM< -0.0000000001))stop("Error: no negative TPM is allowed")

	write.table(contig.TPM, fin_TPMtab, sep="\t", row.names=F, quote=F)

	cmd<-paste("RgenerateRNAseqReads,--contigFasta",contig.file,"--expressionLevels",fin_TPMtab,"--outputPrefix",output.prefix, "--totalReads",sprintf("%d",out.sample.size), "--readLen",read.length, sep=",")
	if(isPairedEndOutput) cmd<-paste(cmd, "--pairedEnd,--insertionLenMean",Insertion.Length.Mean, "--insertionLenMax",Insertion.Length.Max,"--insertionLenMin",Insertion.Length.Min,"--insertionLenSigma",Insertion.Length.Sigma, sep=",")
	if(simplify.contig.names) cmd<-paste(cmd, "--simpleContigId", sep=",")
	if(truth.in.read.names) cmd<-paste(cmd, "--truthInReadNames", sep=",")
	if(qualfile != NULL) cmd<-paste(cmd, "--qualityRefFile", qualfile, sep=",")

	#print(substr(cmd,1,2000))
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_generate_random_RNAseq_reads",as.integer(n), as.character(cmd),PACKAGE="Rsubread")

	summ <- read.delim(paste0(output.prefix,".truthCounts"), stringsAsFactors=FALSE, header=T, comment.char="#", colClasses=c("character", "numeric", "numeric"))
	file.remove(fin_TPMtab)
	summ
}
