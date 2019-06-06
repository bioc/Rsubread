qualityScores <- function(filename,input_format="gzFASTQ",offset=33,nreads=10000)
{
    .check_string_param(input_format,"input_format")
	if(length(filename)>1) stop("The qualityScores function only allows one input file.")
	filename <- .check_and_NormPath(filename, mustWork=T, opt="filename")

	if (file.exists(filename) == FALSE)
		stop("Can not find the input file. Pleaes check whether the file name is correct.")
	
	score_file = paste(".Rsubread_qualityScores_score_pid",Sys.getpid(),sep="")
	
	opt <- paste("-i",filename,"-o",score_file,sep=.R_param_splitor)

	if(tolower(input_format) == "gzfastq")
		opt <- paste(opt,"--gzFASTQinput",sep=.R_param_splitor)
	  	
	if(tolower(input_format) == "sam")
		opt <- paste(opt,"--SAMinput",sep=.R_param_splitor)

	if(tolower(input_format) == "bam")
		opt <- paste(opt,"--BAMinput",sep=.R_param_splitor)

	opt <- paste(opt,"--phred-offset",offset,"--counted-reads",nreads,sep=.R_param_splitor)

	cmd <- paste("qualityScores",opt,sep=.R_param_splitor)
	n <- length(unlist(strsplit(cmd,.R_param_splitor)))
	C_args <- .C("R_qualityScores_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
	
	scores <- read.csv(score_file,header=FALSE,stringsAsFactors=FALSE)
    scores <- scores[, colSums( is.na(scores) ) < nrow(scores) ] # remove all-NA columns
	scores <- as.matrix(scores)
	colnames(scores) <- 1:ncol(scores)
	
	file.remove(score_file)
	
	scores
}

