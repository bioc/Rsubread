txUnique <- function(GTF_Annotation_File, Feature_Type = "exon", Gene_ID_Attribute = "gene_id", Transcript_ID_Attribute = "transcript_id"){
	GTF_Annotation_File <- .check_and_NormPath(GTF_Annotation_File, mustWork=T, opt="GTF_Annotation_File")

	fout <- file.path(".",paste(".Rsubread_txUnique_pid",Sys.getpid(),sep=""))
	opt <- paste("-a",GTF_Annotation_File,sep="\t")
	opt <- paste(opt,"-o",fout,sep="\t")
	opt <- paste(opt,"-g",Gene_ID_Attribute,sep="\t")
	opt <- paste(opt,"-t",Transcript_ID_Attribute,sep="\t")
	opt <- paste(opt,"-f",Feature_Type,sep="\t")
	cmd <- paste("txUnique",opt,sep="\t")
	n <- length(unlist(strsplit(cmd,"\t")))
	C_args <- .C("R_txUnique_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
	if(file.exists(fout)){
		z <- read.delim(fout)
		z <- z[ order(z[,'Gene_ID'], z[,'Transcript_ID']) ,]
		file.remove(fout)
		z
	}else{
		print("No output was generated.")
		NA
	}
}
