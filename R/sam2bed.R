sam2bed <- function(samfile,bedfile,readlen)
{
    .check_string_param(readlen,"readlen")
	samfile <- .check_and_NormPath(samfile, mustWork=T, opt="samfile")
	bedfile <- .check_and_NormPath(bedfile, mustWork=F, opt="bedfile")
	opt <- paste("-n",readlen,samfile,bedfile,sep=.R_param_splitor)
	cmd <- paste("sam2bed",opt,sep=.R_param_splitor)
	n <- length(unlist(strsplit(cmd,.R_param_splitor)))
	C_args <- .C("R_sam2bed_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
}
