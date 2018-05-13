flattenGTF <- function(GTFfile, GTF.featureType="exon",GTF.attrType="gene_id"){

  GTFfile <- normalizePath(GTFfile, mustWork=T)

  fout <- file.path(".",paste(".Rsubread_flattenGTF_pid",Sys.getpid(),sep=""))
  cmd <- paste("RflattenGTF","-a",GTFfile,sep=",");
  cmd <- paste(cmd,"-g",GTF.attrType,sep=",");
  cmd <- paste(cmd,"-t",GTF.featureType,sep=",");
  cmd <- paste(cmd,"-o",fout,sep=",");

  n <- length(unlist(strsplit(cmd,",")))
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
