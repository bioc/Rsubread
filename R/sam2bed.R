sam2bed <- function(samfile,bedfile,readlen)
{
	opt <- paste("-n",readlen,samfile,bedfile,sep=",")
	cmd <- paste("sam2bed",opt,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	.C("R_sam2bed_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
}
