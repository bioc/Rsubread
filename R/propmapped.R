propmapped <- function(files,countFragments=TRUE,properlyPaired=FALSE,verbose=FALSE)
{
  files <- .check_and_NormPath(files, mustWork=T, opt="files")
  fout <- file.path(".",paste(".Rsubread_propmapped_pid",Sys.getpid(),sep=""))

  for(i in 1:length(files)){
    opt <- paste("-i",files[i],sep=",")

    if(countFragments)
    opt <- paste(opt,"-f",sep=",")
    if(properlyPaired)
    opt <- paste(opt,"-p",sep=",")

    if(verbose)
    opt <- paste(opt,"-V",sep=",")

    opt <- paste(opt,"-o",fout,sep=",")
    cmd <- paste("propmapped",opt,sep=",")
    n <- length(unlist(strsplit(cmd,",")))
    C_args <- .C("R_propmapped_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
  }

  x1 <- read.delim(fout,header=FALSE,stringsAsFactors=FALSE,sep=",")
  file.remove(fout)

  colnames(x1) <- c("Samples","NumTotal","NumMapped","PropMapped")
  x1
}
