repair <- function(inFiles,inFormat="BAM",outFiles=paste(inFiles,"repair",sep="."),addDummy=TRUE,fullData=TRUE,compress=FALSE,nthreads=8)
{
	inFiles <- normalizePath(inFiles, mustWork=T)
	outFiles <- normalizePath(outFiles, mustWork=F)

	if(length(inFiles) != length(outFiles))
		stop("Number of input files is different from number of output files.")

	if(!all(file.exists(inFiles)))
		stop("One or more input files cannot be found.")

	if(!all(file.exists(dirname(outFiles))))
		stop("Invalid path was found in output file names.")
	
	for(i in 1:length(inFiles)){
		opt <- paste("-i",inFiles[i],sep=",")

		if(tolower(inFormat) == "sam")
			opt <- paste(opt,"-S",sep=",")	  

		if(!addDummy)
			opt <- paste(opt,"-d",sep=",")

		if(!fullData)
			opt <- paste(opt,"-t",sep=",")

		if(compress)
			opt <- paste(opt,"-c",sep=",")

		opt <- paste(opt,"-T",nthreads,"-o",outFiles[i],sep=",")

		cmd <- paste("repair",opt,sep=",")
		n <- length(unlist(strsplit(cmd,",")))
		C_args <- .C("R_repair_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
	}
}
