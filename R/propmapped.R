propmapped <- function(samfiles)
{
    fout <- file.path("/tmp",paste(".Rsubread_propmapped_pid",Sys.getpid(),sep=""))

    for(i in 1:length(samfiles)){
	cmd <- paste("propmapped",samfiles[i],fout,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	y <- .C("R_propmapped_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
    }

    x1 <- read.delim(fout,header=FALSE,stringsAsFactors=FALSE,sep=",")
    file.remove(fout)

    colnames(x1) <- c("Samples","NumTotal","NumMapped","PropMapped")
    x1
}

