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

generateRNAseqReads <- function(transcript.file, transcript.TPM, output.prefix, read.length=75, isPairedEndOutput=F, Fragment.Length.Min=100, Fragment.Length.Max=500, Fragment.Length.Mean=150, Fragment.Length.Sigma=25){
}
