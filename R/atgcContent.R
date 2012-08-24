atgcContent <- function(filename, basewise=FALSE)
{
	if (file.exists(filename) == FALSE){
		print("Souce file specified doesn't exist!")
	} 
	else 
	{
	  sequence_file <- paste("./.Rsubread_atgcContent_sequence_pid",Sys.getpid(),sep="")
	  perc_file <- paste("./.Rsubread_atgcContent_percentage_pid",Sys.getpid(),sep="")
	  .C("retrieve_sequence", as.character(filename), as.character(sequence_file), PACKAGE="Rsubread")
	  .C("atgcContent", as.character(sequence_file), as.character(perc_file), as.integer(basewise))
	  data <- as.matrix(read.csv(perc_file, header=T))
	  file.remove(sequence_file)
	  file.remove(perc_file)
	  if (basewise){
		percentage <- data
		percentage <- as.matrix(percentage[-1,])
		percentage <- as.matrix(t(percentage))
		colnames(percentage) <- 1:ncol(percentage)
		return(percentage)
	  } 
	  else 
	  {
		return(data)
	  }
	}
}
