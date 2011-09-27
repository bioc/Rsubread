featureCounts <- function(SAMfiles,type="gene",species="mm",annot=NULL)
{
	if(is.null(annot)){
	  if(species == "mm"){
	    ann <- system.file("annot","mm9_NCBI_exon.txt",package="Rsubread")
	    cat("Mouse annotation from NCBI Build 37.2 is used.\n")
	  }
	  else{
	    if(species == "hg"){
	      ann <- system.file("annot","hg19_NCBI_exon.txt",package="Rsubread")
	      cat("Human annotation from NCBI Build 37.2 is used.\n")
	    }
	    else
	      stop("In-built annotation for species", species, "is not available. Please make sure the species name is correct or provide an annotation file for it.\n")
	  }
	}
	else{
	  ann <- annot
	}

	fout <- file.path("/tmp",paste(".Rsubread_featureCounts_pid",Sys.getpid(),sep=""))

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
	  gene_length <- rowsum(x$End-x$Start+1,x$EntrezID)
	  z <- list(counts=as.matrix(y),annotation=data.frame(EntrezID=rownames(y),GeneLength=gene_length,stringsAsFactors=FALSE),targets=SAMfiles)
	}
	else{
	  y <- as.matrix(x[,-c(1:4)]) 
	  rownames(y) <- x$EntrezID
	  z <- list(counts=y,annotation=data.frame(x[,1:4],ExonLength=x$End-x$Start+1,stringsAsFactors=FALSE),targets=SAMfiles)
	}
	z	
}

