exactSNP <- function(readFile,isBAM=FALSE,refGenomeFile,outputFile,qvalueCutoff=12,minAllelicFraction=0,minAllelicBases=1,minReads=1,maxReads=3000, minBaseQuality=13,nTrimmedBases=3,nthreads=1)
{
  opt <- paste("-i",readFile,sep=",")
  
  if(isBAM)
    opt <- paste(opt,"-b",sep=",")

  opt <- paste(opt,"-g",refGenomeFile,"-o",outputFile,"-Q",qvalueCutoff,"-f",minAllelicFraction,"-n",minAllelicBases,"-r",minReads,"-x",maxReads,"-s",minBaseQuality,"-t",nTrimmedBases,"-T",nthreads,sep=",")

  cmd <- paste("SNPcalling",opt,sep=",")
  n <- length(unlist(strsplit(cmd,",")))
  C_args <- .C("R_SNPcalling_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

}

