getInBuiltAnnotation <- function(annotation="mm9")
{

	 switch(tolower(as.character(annotation)),
	    mm9={
	      ann <- system.file("annot","mm9_RefSeq_exon.txt",package="Rsubread")
	      cat("NCBI RefSeq annotation for mm9 (build 37.2).\n")
		},
	    mm10={
	      ann <- system.file("annot","mm10_RefSeq_exon.txt",package="Rsubread")
	      cat("NCBI RefSeq annotation for mm10 (build 38.1).\n")
		 },
	    hg19={
	      ann <- system.file("annot","hg19_RefSeq_exon.txt",package="Rsubread")
	      cat("NCBI RefSeq annotation for hg19 (build 37.2).\n")
	       },
	       {
		stop("In-built annotation for ", annot.inbuilt, " is not available.\n")
	       }
	  ) # end switch
	  
	  x <- read.delim(ann, stringsAsFactors=FALSE)
	  x
}

