flattenGTF <- function(GTFfile, GTF.featureType="exon",GTF.attrType="gene_id"){
  .check_string_param(GTF.featureType,"GTF.featureType")
  .check_string_param(GTF.attrType,"GTF.attrType")

  GTFfile <- .check_and_NormPath(GTFfile, mustWork=T, opt="GTFfile")

  fout <- file.path(".",paste(".Rsubread_flattenGTF_pid",Sys.getpid(),sep=""))
  cmd <- paste("RflattenGTF","-a",GTFfile,sep=.R_param_splitor);
  cmd <- paste(cmd,"-g",GTF.attrType,sep=.R_param_splitor);
  cmd <- paste(cmd,"-t",GTF.featureType,sep=.R_param_splitor);
  cmd <- paste(cmd,"-o",fout,sep=.R_param_splitor);

  n <- length(unlist(strsplit(cmd,.R_param_splitor)))
  C_args <- .C("R_flattenGTF_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
  if(file.exists(fout)){
	z <- read.delim(fout)
	file.remove(fout)
	z
  }else{
	print("No output was generated.")
	NA
  }
}
