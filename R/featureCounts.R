featureCounts <- function(SAMfiles,type="gene",species="mm")
{
	if(species != "mm")
	  stop("Only mouse data are supported at this stage.") 
	ann <- system.file("annot","mm9_NCBI_exon.txt",package="Rsubread")

	fout <- system.file("extdata",package="Rsubread")
	fout <- file.path(fout,paste(".Rsubread_pid",Sys.getpid(),sep=""))

	for(i in 1:length(SAMfiles)){
	  cat("Processing", SAMfiles[i], " ...\n")
	  cmd <- paste("readSummary",ann,SAMfiles[i],fout,sep=",")
	  n <- length(unlist(strsplit(cmd,",")))
	  C_args <- .C("R_readSummary_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
	  x1 <- read.delim(fout,stringsAsFactors=FALSE)  
	  if(i == 1)
	    x <- x1
	  else
	    x <- cbind(x,x1$nreads)
	}
	colnames(x)[1:4] <- c("EntrezID","Chr","Start","End")
	colnames(x)[-c(1:4)] <- SAMfiles
	file.remove(fout)
	if(type == "gene"){
	  y <- rowsum(x[,-c(1:4)],x$EntrezID)
	  y <- data.frame(EntrezID=rownames(y),y,stringsAsFactors=FALSE)
	  colnames(y)[-1] <- SAMfiles
	  return(y)
	}
	else
	  return(x)
}

