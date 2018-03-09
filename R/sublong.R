sublong <- function(index, readFile, outputFile, outputFormat="BAM" , nthreads=1) {
	opt <- paste("-i",index,sep=",")
	if(tolower(outputFormat) == "sam")
		opt <- paste(opt,"--SAMoutput",sep=",")	  
	opt <- paste(opt,"-r",readFile,"-T",nthreads,"-o",outputFile,sep=",")

	cmd <- paste("Rsublong",opt,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_sublong_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
}
