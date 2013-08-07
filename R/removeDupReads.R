removeDupReads <- function(SAMfile,threshold=50,outputFile)
{
  cmd <- paste("removeDupReads","-i",SAMfile,"-r",threshold,"-o",outputFile,sep=",")
  n <- length(unlist(strsplit(cmd,",")))
  C_args <- .C("R_removeDupReads_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

}

