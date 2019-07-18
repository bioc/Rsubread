exactSNP <- function(readFile,isBAM=FALSE,refGenomeFile,SNPAnnotationFile=NULL,outputFile=paste0(readFile,".exactSNP.VCF"),qvalueCutoff=12,minAllelicFraction=0,minAllelicBases=1,minReads=1,maxReads=1000000, minBaseQuality=13,nTrimmedBases=3,nthreads=1)
{
  if(length(readFile) > 1L)
	stop("only one input file is allowed.")

  readFile <- normalizePath(readFile, mustWork=T)
  refGenomeFile <- normalizePath(refGenomeFile, mustWork=T);
  outputFile <- normalizePath(outputFile, mustWork=F)
  if(!is.null(SNPAnnotationFile))SNPAnnotationFile <- normalizePath(SNPAnnotationFile, mustWork=T)

  opt <- paste("-i",readFile,sep=.R_param_splitor)
  
  if(isBAM)
    opt <- paste(opt,"-b",sep=.R_param_splitor)
	
  if(!is.null(SNPAnnotationFile))
    opt <- paste(opt,"-a",SNPAnnotationFile,sep=.R_param_splitor)

  opt <- paste(opt,"-g",refGenomeFile,"-o",outputFile,"-Q",qvalueCutoff,"-f",minAllelicFraction,"-n",minAllelicBases,"-r",minReads,"-x",maxReads,"-s",minBaseQuality,"-t",nTrimmedBases,"-T",nthreads,sep=.R_param_splitor)

  cmd <- paste("SNPcalling",opt,sep=.R_param_splitor)
  n <- length(unlist(strsplit(cmd,.R_param_splitor)))
  C_args <- .C("R_SNPcalling_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

}

