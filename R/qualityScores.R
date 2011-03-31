qualityScores <- function(filename, offset=64, nreads = 10000)
{
	if (file.exists(filename) == FALSE){
		print("Source file specified doesn't exist.");
	} 
	else
	{
		score_file = ".__Rsubread_score_Rsubread__000000" 
		temp_file = ".__Rsubread_temp_Rsubread__000000"
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
