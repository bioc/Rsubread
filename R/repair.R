repair <- function(inFiles,inFormat="BAM",outFiles=paste(inFiles,"repair",sep="."),addDummy=TRUE,fullData=TRUE,compress=FALSE,nthreads=8)
{
    .check_string_param(inFormat,"inFormat")
	inFiles <- .check_and_NormPath(inFiles, mustWork=T, opt="inFiles")
	outFiles <- .check_and_NormPath(outFiles, mustWork=F, opt="outFiles")

	if(length(inFiles) != length(outFiles))
		stop("Number of input files is different from number of output files.")

	if(!all(file.exists(inFiles)))
		stop("One or more input files cannot be found.")

	if(!all(file.exists(dirname(outFiles))))
		stop("Invalid path was found in output file names.")
	
	for(i in 1:length(inFiles)){
		opt <- paste("-i",inFiles[i],sep=.R_param_splitor)

		if(tolower(inFormat) == "sam")
			opt <- paste(opt,"-S",sep=.R_param_splitor)	  

		if(!addDummy)
			opt <- paste(opt,"-d",sep=.R_param_splitor)

		if(!fullData)
			opt <- paste(opt,"-t",sep=.R_param_splitor)

		if(compress)
			opt <- paste(opt,"-c",sep=.R_param_splitor)

		opt <- paste(opt,"-T",nthreads,"-o",outFiles[i],sep=.R_param_splitor)

		cmd <- paste("repair",opt,sep=.R_param_splitor)
		n <- length(unlist(strsplit(cmd,.R_param_splitor)))
		C_args <- .C("R_repair_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
	}
}
