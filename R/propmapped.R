propmapped <- function(samfile)
{
	opt <- samfile
	cmd <- paste("unmapped",opt,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	y <- .C("R_unmapped_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
}

