featureCounts <- function(SAMfiles,type="gene",species="mm9",annot=NULL,isPairedEnd=FALSE,min_distance=10,max_distance=2000)
{
	flag <- FALSE

	if(is.null(annot)){
	  switch(tolower(as.character(species)),
	    mm9={
	      ann <- system.file("annot","mm9_NCBI_exon.txt",package="Rsubread")
	      cat("Mouse annotation NCBI Build 37.2 is used.\n")
		},
	    mm10={
	      ann <- system.file("annot","mm10_NCBI_annotation_exon_sorted.txt",package="Rsubread")
	      cat("Mouse annotation NCBI Build 38.1 is used.\n")
		 },
	    hg={
	      ann <- system.file("annot","hg19_NCBI_exon.txt",package="Rsubread")
	      cat("Human annotation NCBI Build 37.2 is used.\n")
	       },
	       {
		stop("In-built annotation for ", species, " is not available.\n")
	       }
	  ) # end switch
	}
	else{
	  if(is.character(annot)){
	    cat("Please make sure the annotation file you provided has been SORTED by chromosome names.\n")
	    ann <- annot
	  }
	  else{
	    annot_df <- as.data.frame(annot,stringsAsFactors=FALSE)
	    if(sum(c("entrezid","chromosome","chr_start","chr_stop") %in% colnames(annot_df)) != 4)
	      stop("annot is not a data frame or required columns can not be found in annot.\n")
	    annot_df$chromosome <- as.character(annot_df$chromosome)
	    annot_df <- annot_df[order(annot_df$chromosome),]
	    if(any(nchar(annot_df$chromosome) > 40))
	      annot_df$chromosome <- substring(annot_df$chromosome,1,40)
	    fout_annot <- file.path(".",paste(".Rsubread_UserProvidedAnnotation_pid",Sys.getpid(),sep=""))
	    write.table(x=annot_df,file=fout_annot,sep="\t",row.names=FALSE,quote=FALSE)
	    ann <- fout_annot
	    flag <- TRUE
	  }
	}

	fout <- file.path(".",paste(".Rsubread_featureCounts_pid",Sys.getpid(),sep=""))

	for(i in 1:length(SAMfiles)){
	  cat("Processing", SAMfiles[i], " ...\n")
	  cmd <- paste("readSummary",ann,SAMfiles[i],fout,as.numeric(isPairedEnd),min_distance,max_distance,sep=",")
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
	if(flag) 
	  file.remove(fout_annot)

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

