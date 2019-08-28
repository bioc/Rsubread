buildindex <- function(basename,reference,gappedIndex=FALSE,indexSplit=FALSE,memory=8000,TH_subread=100,colorspace=FALSE)
{
    if(!.is.64bit.system()) cat("WARNING: your system seems to be 32-bit. Rsubread supports 32-bit sustems to a very limited extend.\nIt is highly recommended to run Rsubread on a 64-bit system to avoid errors.\n\n")

	basename <- .check_and_NormPath(basename, mustWork=F, opt="basename")
	reference <- .check_and_NormPath(reference, mustWork=T, opt="reference")

	opt <- paste("-o",basename,"-f",TH_subread,"-M",memory,reference,sep=.R_param_splitor)

	if(gappedIndex == FALSE) opt <- paste("-F",opt,sep=.R_param_splitor)
	if(indexSplit == FALSE) opt <- paste("-B",opt,sep=.R_param_splitor)
	if(colorspace) opt <- paste("-c",opt,sep=.R_param_splitor)
		
	cmd <- paste("subread-buildindex",opt,sep=.R_param_splitor)
	n <- length(unlist(strsplit(cmd,.R_param_splitor)))
	C_args <- .C("R_buildindex_wrapper",argc=as.integer(n),argv=as.character(cmd),PACKAGE="Rsubread")
}
