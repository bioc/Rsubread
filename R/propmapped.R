propmapped <- function(files,countFragments=TRUE,properlyPaired=FALSE,verbose=FALSE)
{
  out.base.names <- basename(files)
  if(any(duplicated(out.base.names))){
    out.table.rows <- files
  }else{
    out.table.rows <- out.base.names
  }

  files <- .check_and_NormPath(files, mustWork=T, opt="files")
  if(any(duplicated(files)))stop("Error: duplicated input files are provided. No result is generated.")
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

  x1 <- read.csv(fout,header=FALSE,row.names=1,stringsAsFactors=FALSE)
  rownames(x1) <- out.table.rows
  file.remove(fout)

  colnames(x1) <- c("NumTotal","NumMapped","PropMapped")
  x1
}
