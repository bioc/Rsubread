removeDupReads <- function(SAMfile,threshold=50,outputFile)
{
  SAMfile <- .check_and_NormPath(SAMfile, mustWork=T, opt="SAMfile")
  outputFile <- .check_and_NormPath(outputFile, mustWork=F, opt="outputFile")
  cmd <- paste("removeDupReads","-i",SAMfile,"-r",threshold,"-o",outputFile,sep=",")
  n <- length(unlist(strsplit(cmd,",")))
  C_args <- .C("R_removeDupReads_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

}

