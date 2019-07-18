getInBuiltAnnotation <- function(annotation="mm10")
{
     .check_string_param(annotation,"annotation")

	 switch(tolower(as.character(annotation)),
	    mm9={
	      ann <- system.file("annot","mm9_RefSeq_exon.txt",package="Rsubread")
	      message("NCBI RefSeq annotation for mm9 (build 37.2).")
		},
	    mm10={
	      ann <- system.file("annot","mm10_RefSeq_exon.txt",package="Rsubread")
	      message("NCBI RefSeq annotation for mm10 (build 38.1).")
		 },
	    hg19={
	      ann <- system.file("annot","hg19_RefSeq_exon.txt",package="Rsubread")
	      message("NCBI RefSeq annotation for hg19 (build 37.2).")
	       },
	    hg38={
	      ann <- system.file("annot","hg38_RefSeq_exon.txt",package="Rsubread")
	      message("NCBI RefSeq annotation for hg38 (build 38.2).")
	       },
	       {
		stop("In-built annotation for ", annotation, " is not available.\n")
	       }
	  ) # end switch
	  
	  x <- read.delim(ann, sep="\t", quote="", colClasses=c("character","character","integer","integer","character"), comment.char="")
	  x
}

