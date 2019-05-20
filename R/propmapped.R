propmapped <- function(files,countFragments=TRUE,properlyPaired=FALSE,verbose=FALSE)
{


  out.base.names <- basename(files)
  if(any(duplicated(out.base.names))){
    out.table.rows <- files
  }else{
    out.table.rows <- out.base.names
  }

  out.table.rows <-gsub("[[:punct:]]+", ".", out.table.rows)
  out.table.rows <-gsub(" ", ".", out.table.rows)

  files <- .check_and_NormPath(files, mustWork=T, opt="files")
  fout <- file.path(".",paste(".Rsubread_propmapped_pid",Sys.getpid(),sep=""))

  for(i in 1:length(files)){
    opt <- paste("-i",files[i],sep=.R_param_splitor)

    if(countFragments)
    opt <- paste(opt,"-f",sep=.R_param_splitor)
    if(properlyPaired)
    opt <- paste(opt,"-p",sep=.R_param_splitor)

    if(verbose)
    opt <- paste(opt,"-V",sep=.R_param_splitor)

    opt <- paste(opt,"-o",fout,sep=.R_param_splitor)
    cmd <- paste("propmapped",opt,sep=.R_param_splitor)
    n <- length(unlist(strsplit(cmd,.R_param_splitor)))
    C_args <- .C("R_propmapped_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
  }

  x1 <- read.csv(fout,header=FALSE,row.names=1,stringsAsFactors=FALSE)
  rownames(x1) <- out.table.rows
  file.remove(fout)

  colnames(x1) <- c("NumTotal","NumMapped","PropMapped")
  x1
}
