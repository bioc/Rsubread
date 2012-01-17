subjunc <- function(index,readfile1,readfile2=NULL,output_file,nthreads=1,indels=5,min_distance=50,max_distance=600,PE_orientation="fr")
{
	opt <- paste("-i",index,"-r",readfile1,"-o",output_file,"-T",nthreads,"-I",indels,sep=",")
	if(!is.null(readfile2)) opt <- paste(opt,"-R",readfile2,"-d",min_distance,"-D",max_distance,"-S",PE_orientation,sep=",")
	cmd <- paste("subjunc",opt,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_junction_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
}
