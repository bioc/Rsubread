qualityScores <- function(filename,input_format="gzFASTQ",offset=33,nreads=10000)
{
	if (file.exists(filename) == FALSE)
		stop("Can not find the input file. Pleaes check whether the file name is correct.")
	
	score_file = paste(".Rsubread_qualityScores_score_pid",Sys.getpid(),sep="")
	
	opt <- paste("-i",filename,"-o",score_file,sep=",")

	if(tolower(input_format) == "gzfastq")
		opt <- paste(opt,"--gzFASTQinput",sep=",")
	  	
	if(tolower(input_format) == "sam")
		opt <- paste(opt,"--SAMinput",sep=",")

	if(tolower(input_format) == "bam")
		opt <- paste(opt,"--BAMinput",sep=",")

	opt <- paste(opt,"--phred-offset",offset,"--counted-reads",nreads,sep=",")

	cmd <- paste("qualityScores",opt,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_qualityScores_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
	
	scores <- read.csv(score_file,header=FALSE,stringsAsFactors=FALSE)
	scores <- as.matrix(scores)
	
	file.remove(score_file)
	
	scores
}

