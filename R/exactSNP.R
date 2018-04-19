exactSNP <- function(readFile,isBAM=FALSE,refGenomeFile,SNPAnnotationFile=NULL,outputFile=paste(readFile,".exactSNP.VCF",sep=""),qvalueCutoff=12,minAllelicFraction=0,minAllelicBases=1,minReads=1,maxReads=1000000, minBaseQuality=13,nTrimmedBases=3,nthreads=1)
{

  if(length(readFile) > 1)
	stop("You are not allowed to provide more than one input file.")

  readFile <- normalizePath(readFile, mustWork=T)
  refGenomeFile <- normalizePath(refGenomeFile, mustWork=T);
  outputFile <- normalizePath(outputFile, mustWork=F)
  if(!is.null(SNPAnnotationFile))SNPAnnotationFile <- normalizePath(SNPAnnotationFile, mustWork=T)

  opt <- paste("-i",readFile,sep=",")
  
  if(isBAM)
    opt <- paste(opt,"-b",sep=",")
	
  if(!is.null(SNPAnnotationFile))
    opt <- paste(opt,"-a",SNPAnnotationFile,sep=",")

  opt <- paste(opt,"-g",refGenomeFile,"-o",outputFile,"-Q",qvalueCutoff,"-f",minAllelicFraction,"-n",minAllelicBases,"-r",minReads,"-x",maxReads,"-s",minBaseQuality,"-t",nTrimmedBases,"-T",nthreads,sep=",")

  cmd <- paste("SNPcalling",opt,sep=",")
  n <- length(unlist(strsplit(cmd,",")))
  C_args <- .C("R_SNPcalling_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

}

