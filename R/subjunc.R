subjunc <- function(index,samfile,output_file,nsubreads=14,paired_end=FALSE,nthreads=1,indels=5,minFragLength=50,maxFragLength=600,PE_orientation="fr")
{
	opt <- paste("-i",index,"-o",output_file,"-n",nsubreads,"-T",nthreads,"-I",indels,sep=",")

	if(paired_end) 
	  opt <- paste(opt,"-2",samfile,"-d",minFragLength,"-D",maxFragLength,sep=",")
	else
	  opt <- paste(opt,"-1",samfile,sep=",")

	cmd <- paste("subjunc",opt,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_junction_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
}
