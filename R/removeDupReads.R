removeDupReads <- function(SAMfile,threshold=50,nthreads=1,outputFile)
{
  cmd <- paste("removeDupReads","-i",SAMfile,"-r",threshold,"-T",nthreads,"-o",outputFile,sep=",")
  n <- length(unlist(strsplit(cmd,",")))
  C_args <- .C("R_removeDupReads_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

}

