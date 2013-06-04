featureCounts <- function(files,file.type="SAM",annot=NULL,genome="mm9",isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",useMetaFeatures=TRUE,allowMultiOverlap=FALSE,nthreads=1,strandSpecific=0,countMultiMappingReads=FALSE,minMQS=0,isPairedEnd=FALSE,requireBothEndsMapped=FALSE,checkPairedEndDistance=FALSE,min.distance=50,max.distance=600,countChimericFragments=TRUE)
{
	flag <- FALSE

	if(is.null(annot)){
	  switch(tolower(as.character(genome)),
	    mm9={
	      ann <- system.file("annot","mm9_RefSeq_exon.txt",package="Rsubread")
	      cat("Mouse annotation NCBI Build 37.2 is used.\n")
		},
	    mm10={
	      ann <- system.file("annot","mm10_RefSeq_exon.txt",package="Rsubread")
	      cat("Mouse annotation NCBI Build 38.1 is used.\n")
		 },
	    hg={
	      ann <- system.file("annot","hg19_RefSeq_exon.txt",package="Rsubread")
	      cat("Human annotation NCBI Build 37.2 is used.\n")
	       },
	       {
		stop("In-built annotation for ", genome, " is not available.\n")
	       }
	  ) # end switch
	}
	else{
	  if(is.character(annot)){
	    ann <- annot
	  }
	  else{
	    annot_df <- as.data.frame(annot,stringsAsFactors=FALSE)
	    if(sum(c("geneid","chr","start","end", "strand") %in% tolower(colnames(annot_df))) != 5)
	      stop("The format of the provided annotation data is incorrect. Please refer to the help page for the correct data format.\n")
		colnames(annot_df) <- tolower(colnames(annot_df))
		annot_df <- data.frame(geneid=annot_df$geneid,chr=annot_df$chr,start=annot_df$start,end=annot_df$end,strand=annot_df$strand,stringsAsFactors=FALSE)		
	    annot_df$chr <- as.character(annot_df$chr)
	    if(any(nchar(annot_df$chr) > 40))
	      annot_df$chr <- substring(annot_df$chr,1,40)
	    fout_annot <- file.path(".",paste(".Rsubread_UserProvidedAnnotation_pid",Sys.getpid(),sep=""))
	    write.table(x=annot_df,file=fout_annot,sep="\t",row.names=FALSE,quote=FALSE)
	    ann <- fout_annot
	    flag <- TRUE
	  }
	}

	fout <- file.path(".",paste(".Rsubread_featureCounts_pid",Sys.getpid(),sep=""))

	for(i in 1:length(files)){
	  cat("Processing", files[i], " ...\n")
	  cmd <- paste("readSummary",ann,files[i],fout,as.numeric(isPairedEnd),min.distance,max.distance,as.numeric(tolower(file.type)=="sam"),as.numeric(allowMultiOverlap),as.numeric(useMetaFeatures),nthreads,as.numeric(isGTFAnnotationFile),strandSpecific,as.numeric(FALSE),as.numeric(requireBothEndsMapped),as.numeric(!countChimericFragments),as.numeric(checkPairedEndDistance),GTF.featureType,GTF.attrType,minMQS,as.numeric(countMultiMappingReads),sep=",")
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

