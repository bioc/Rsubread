align <- function(index,readfile1,readfile2=NULL,output_file,nsubreads=10,TH1=3,TH2=3,nthreads=1,indels=0,min_distance=50,max_distance=500,PE_orientation="fr")
{
	opt <- paste("-i",index,"-r",readfile1,"-o",output_file,"-n",nsubreads,"-m",TH1,"-T",nthreads,"-I",indels,sep=",")
	if(!is.null(readfile2)) opt <- paste(opt,"-R",readfile2,"-p",TH2,"-d",min_distance,"-D",max_distance,"-S",PE_orientation,sep=",")
	cmd <- paste("subread-align",opt,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_align_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
}
