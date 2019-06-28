flattenGTF <- function(GTFfile, GTF.featureType="exon",GTF.attrType="gene_id"){
  .check_string_param(GTF.featureType,"GTF.featureType")
  .check_string_param(GTF.attrType,"GTF.attrType")

  GTFfile <- .check_and_NormPath(GTFfile, mustWork=T, opt="GTFfile")

  fout <- file.path(".",paste(".Rsubread_flattenGTF_pid",Sys.getpid(),sep=""))
  cmd <- paste("RflattenGTF","-a",GTFfile,sep=.R_param_splitor)
  cmd <- paste(cmd,"-g",GTF.attrType,sep=.R_param_splitor)
  cmd <- paste(cmd,"-t",GTF.featureType,sep=.R_param_splitor)
  cmd <- paste(cmd,"-o",fout,sep=.R_param_splitor)

  n <- length(unlist(strsplit(cmd,.R_param_splitor)))
  C_args <- .C("R_flattenGTF_wrapper",n,cmd,PACKAGE="Rsubread")
  if(file.exists(fout)){
	z <- read.delim(fout,stringsAsFactors=FALSE,colClasses=c("character","character","integer","integer","character"))
	file.remove(fout)
	z
  }else{
	warning("No output was generated.")
	data.frame()
  }
}
