removeDupReads <- function(inputFile,threshold=50,outputFile,outputFormat="BAM")
{
  SAMfile <- .check_and_NormPath(inputFile, mustWork=T, opt="SAMfile")
  outputFile <- .check_and_NormPath(outputFile, mustWork=F, opt="outputFile")
  cmd <- paste("removeDupReads","-i",SAMfile,"-r",threshold,"-o",outputFile,sep=.R_param_splitor)

  if(outputFormat=="SAM") cmd <- paste(cmd, "-S", sep=.R_param_splitor)
  else if(outputFormat!="BAM") stop("ERROR: Unknown output format! Only 'BAM' or 'SAM' are accepted!")

  n <- length(unlist(strsplit(cmd, .R_param_splitor)))
  C_args <- .C("R_removeDupReads_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
}

