buildindex <- function(basename,reference,colorspace=FALSE,memory=3700)
{
opt <- paste("-o",basename,"-M",memory,reference,sep=",")
if(colorspace) opt <- paste("-c",opt,sep=",")
cmd <- paste("subread-buildindex",opt,sep=",")
n <- length(unlist(strsplit(cmd,",")))
C_args <- .C("R_buildindex_wrapper",argc=as.integer(n),argv=as.character(cmd),PACKAGE="Rsubread")
}
