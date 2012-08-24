callSNPs <- function(SAMfile,refGenomeFile,outputFile,minBaseQuality=13,minReadCoverage=5,minAlleleFraction=0.5)
{

  cmd <- paste("SNPcalling","-i",SAMfile,"-g",refGenomeFile,"-o",outputFile,"-I",5,"-p",0.05,"-q",minBaseQuality,"-n",minReadCoverage,"-r",minAlleleFraction,sep=",")
  n <- length(unlist(strsplit(cmd,",")))
  C_args <- .C("R_SNPcalling_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

}

