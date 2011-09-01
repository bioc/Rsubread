qualityScores <- function(filename, offset=64, nreads = 10000)
{
	if (file.exists(filename) == FALSE){
		print("Source file specified doesn't exist.");
	} 
	else
	{
		score_file = paste("/tmp/.Rsubread_qualityScores_score_pid",Sys.getpid(),sep="") 
		temp_file = paste("/tmp/.Rsubread_qualityScores_temp_pid",Sys.getpid(),sep="")
		.C("retrieve_scores", as.character(filename), as.integer(offset), as.integer(nreads), as.character(temp_file), as.character(score_file), PACKAGE="Rsubread")
	if (file.exists(score_file)) {
		scores = as.matrix(read.csv(score_file,header=F));
		file.remove(score_file)
	}
	else
	{
		scores = as.matrix(read.csv(temp_file,header=F));
	}
	file.remove(temp_file)
	colnames(scores) <- 1:ncol(scores)
	return(scores)
	}
}
