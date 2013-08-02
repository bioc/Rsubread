featureCounts <- function(files,file.type="SAM",annot.inbuilt="mm9",annot.ext=NULL,isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",useMetaFeatures=TRUE,allowMultiOverlap=FALSE,nthreads=1,strandSpecific=0,countMultiMappingReads=FALSE,minMQS=0,isPairedEnd=FALSE,requireBothEndsMapped=FALSE,checkFragLength=FALSE,minFragLength=50,maxFragLength=600,countChimericFragments=TRUE)
{
	flag <- FALSE

	if(is.null(annot.ext)){
	  switch(tolower(as.character(annot.inbuilt)),
	    mm9={
	      ann <- system.file("annot","mm9_RefSeq_exon.txt",package="Rsubread")
	      cat("Mouse annotation NCBI Build 37.2 is used.\n")
		},
	    mm10={
	      ann <- system.file("annot","mm10_RefSeq_exon.txt",package="Rsubread")
	      cat("Mouse annotation NCBI Build 38.1 is used.\n")
		 },
	    hg19={
	      ann <- system.file("annot","hg19_RefSeq_exon.txt",package="Rsubread")
	      cat("Human annotation NCBI Build 37.2 is used.\n")
	       },
	       {
		stop("In-built annotation for ", annot.inbuilt, " is not available.\n")
	       }
	  ) # end switch
	}
	else{
	  if(is.character(annot.ext)){
	    ann <- annot.ext
	  }
	  else{
	    annot_df <- as.data.frame(annot.ext,stringsAsFactors=FALSE)
	    if(sum(c("geneid","chr","start","end", "strand") %in% tolower(colnames(annot_df))) != 5)
	      stop("One or more required columns are missing in the provided annotation data. Please refer to help page for annotation format.\n")
		colnames(annot_df) <- tolower(colnames(annot_df))
		annot_df <- data.frame(geneid=annot_df$geneid,chr=annot_df$chr,start=annot_df$start,end=annot_df$end,strand=annot_df$strand,stringsAsFactors=FALSE)		
	    annot_df$chr <- as.character(annot_df$chr)
	    fout_annot <- file.path(".",paste(".Rsubread_UserProvidedAnnotation_pid",Sys.getpid(),sep=""))
	    write.table(x=annot_df,file=fout_annot,sep="\t",row.names=FALSE,quote=FALSE)
	    ann <- fout_annot
	    flag <- TRUE
	  }
	}

	fout <- file.path(".",paste(".Rsubread_featureCounts_pid",Sys.getpid(),sep=""))

	for(i in 1:length(files)){
	  cat("Processing", files[i], " ...\n")
	  cmd <- paste("readSummary",ann,files[i],fout,as.numeric(isPairedEnd),minFragLength,maxFragLength,as.numeric(tolower(file.type)=="sam"),as.numeric(allowMultiOverlap),as.numeric(useMetaFeatures),nthreads,as.numeric(isGTFAnnotationFile),strandSpecific,as.numeric(FALSE),as.numeric(requireBothEndsMapped),as.numeric(!countChimericFragments),as.numeric(checkFragLength),GTF.featureType,GTF.attrType,minMQS,as.numeric(countMultiMappingReads),sep=",")
	  n <- length(unlist(strsplit(cmd,",")))
	  C_args <- .C("R_readSummary_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
	  x1 <- read.delim(fout,stringsAsFactors=FALSE)  
	  if(i == 1)
	    x <- x1
	  else
	    x <- cbind(x,x1$nreads)
	}
	colnames(x)[1:5] <- c("GeneID","Chr","Start","End","Strand")
	colnames(x)[-c(1:5)] <- files

	file.remove(fout)
	if(flag) 
	  file.remove(fout_annot)

	if(useMetaFeatures){
	  y <- rowsum(x[,-c(1:5)],x$GeneID)
	  gene_length <- rowsum(x$End-x$Start+1,x$GeneID)
	  z <- list(counts=as.matrix(y),annotation=data.frame(GeneID=rownames(y),Length=gene_length,stringsAsFactors=FALSE),targets=files)
	}
	else{
	  y <- as.matrix(x[,-c(1:5)]) 
	  rownames(y) <- x$GeneID
	  z <- list(counts=y,annotation=data.frame(x[,1:5],Length=x$End-x$Start+1,stringsAsFactors=FALSE),targets=files)
	}
	z	
}

