removeDupReads <- function(SAMfile,threshold=50,outputFile)
{
  SAMfile <- normalizePath(SAMfile, mustWork=T)
  outputFile <- normalizePath(outputFile, mustWork=F)
  cmd <- paste("removeDupReads","-i",SAMfile,"-r",threshold,"-o",outputFile,sep=",")
  n <- length(unlist(strsplit(cmd,",")))
  C_args <- .C("R_removeDupReads_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

}

