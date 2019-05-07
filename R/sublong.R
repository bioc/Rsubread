sublong <- function(index, readFiles, outputFiles, outputFormat="BAM" , nthreads=1) {
    .check_string_param(outputFormat,"outputFormat")
	readFiles <- .check_and_NormPath(as.character(readFiles), mustWork=T, opt="readFiles")
	outputFiles <- .check_and_NormPath(as.character(outputFiles), mustWork=F , opt="outputFiles")

	if(length(readFiles) != length(outputFiles))
		stop("The number of input file names is different from the number of output file names.")

	if(!all(file.exists(readFiles)))
		stop("One or more input files cannot be found.")

	if(!all(file.exists(dirname(outputFiles))))
		stop("Invalid path was found in output file name(s).")

	for(i in 1:length(readFiles)){
		opt <- paste("-i",index,sep=.R_param_splitor)
		if(tolower(outputFormat) == "sam")
			opt <- paste(opt,"--SAMoutput",sep=.R_param_splitor)
		opt <- paste(opt,"-r",readFiles[i],"-T",nthreads,"-o",outputFiles[i],sep=.R_param_splitor)
	
		cmd <- paste("Rsublong",opt,sep=.R_param_splitor)
		n <- length(unlist(strsplit(cmd,.R_param_splitor)))
		C_args <- .C("R_sublong_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
	}
}
