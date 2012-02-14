subjunc <- function(index,samfile,output_file,paired_end=FALSE,nthreads=1,indels=5,min_distance=50,max_distance=600,PE_orientation="fr")
{
	opt <- paste("-i",index,"-o",output_file,"-T",nthreads,"-I",indels,"-H",2,"--nofull","--extending",sep=",")

	if(paired_end) 
	  opt <- paste(opt,"--pairedSAM",samfile,"-d",min_distance,"-D",max_distance,sep=",")
	else
	  opt <- paste(opt,"--singleSAM",samfile,sep=",")

	cmd <- paste("subjunc",opt,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_junction_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
}
