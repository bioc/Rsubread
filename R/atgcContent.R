atgcContent <- function(filename, basewise=FALSE)
{
	if (file.exists(filename) == FALSE){
		print("Souce file specified doesn't exist!")
	} 
	else 
	{
	  sequence_file <- ".__Rsubread_sequence_Rsubread__000000"
	  perc_file <- ".__Rsubread_percentage_Rsubread__000000"
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
