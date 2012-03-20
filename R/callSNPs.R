callSNPs <- function(SAMfile,readLength,refGenomeFile,outputFile,minReadCoverage=5,minAlleleFraction=0.5)
{

  cmd <- paste("SNPcalling",SAMfile,readLength,refGenomeFile,outputFile,minReadCoverage,minAlleleFraction,sep=",")
  n <- length(unlist(strsplit(cmd,",")))
  C_args <- .C("R_SNPcalling_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

}

