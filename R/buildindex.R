buildindex <- function(basename,reference,gappedIndex=FALSE,indexSplit=FALSE,memory=8000,TH_subread=100,colorspace=FALSE)
{
	basename <- .check_and_NormPath(basename, mustWork=F, opt="basename")
	reference <- .check_and_NormPath(reference, mustWork=T, opt="reference")

	opt <- paste("-o",basename,"-f",TH_subread,"-M",memory,reference,sep=",")

	if(gappedIndex == FALSE) opt <- paste("-F",opt,sep=",")
	if(indexSplit == FALSE) opt <- paste("-B",opt,sep=",")
	if(colorspace) opt <- paste("-c",opt,sep=",")
		
	cmd <- paste("subread-buildindex",opt,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_buildindex_wrapper",argc=as.integer(n),argv=as.character(cmd),PACKAGE="Rsubread")
}
