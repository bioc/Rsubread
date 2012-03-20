removeDupReads <- function(SAMfile,threshold=50)
{
  cmd <- paste("removeDupReads",SAMfile,threshold,sep=",")
  n <- length(unlist(strsplit(cmd,",")))
  C_args <- .C("R_removeDupReads_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

}

