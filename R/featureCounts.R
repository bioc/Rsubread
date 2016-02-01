featureCounts <- function(files,annot.inbuilt="mm9",annot.ext=NULL,isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",chrAliases=NULL,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,minOverlap=1,largestOverlap=FALSE,readExtension5=0,readExtension3=0,read2pos=NULL,countMultiMappingReads=FALSE,fraction=FALSE,minMQS=0,splitOnly=FALSE,nonSplitOnly=FALSE,primaryOnly=FALSE,ignoreDup=FALSE,strandSpecific=0,juncCounts=FALSE,genome=NULL,isPairedEnd=FALSE,requireBothEndsMapped=FALSE,checkFragLength=FALSE,minFragLength=50,maxFragLength=600,PE_orientation="fr",countChimericFragments=TRUE,autosort=TRUE,nthreads=1,maxMOp=10,reportReads=FALSE)
{
	flag <- FALSE

	if(is.null(annot.ext)){
	  switch(tolower(as.character(annot.inbuilt)),
	    mm9={
	      ann <- system.file("annot","mm9_RefSeq_exon.txt",package="Rsubread")
	      cat("NCBI RefSeq annotation for mm9 (build 37.2) is used.\n")
		},
	    mm10={
	      ann <- system.file("annot","mm10_RefSeq_exon.txt",package="Rsubread")
	      cat("NCBI RefSeq annotation for mm10 (build 38.1) is used.\n")
		 },
	    hg19={
	      ann <- system.file("annot","hg19_RefSeq_exon.txt",package="Rsubread")
	      cat("NCBI RefSeq annotation for hg19 (build 37.2) is used.\n")
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
		oldScipen <- options(scipen=999)
	    write.table(x=annot_df,file=fout_annot,sep="\t",row.names=FALSE,quote=FALSE)
		options(oldScipen)
	    ann <- fout_annot
	    flag <- TRUE
	  }
	}

	fout <- file.path(".",paste(".Rsubread_featureCounts_pid",Sys.getpid(),sep=""))

	files_C <- paste(files,collapse=";")
	
	if(nchar(files_C) == 0) stop("No read files provided!")

	genome_C <- genome
	if(is.null(genome))
	  genome_C <- " "
	
	chrAliases_C <- chrAliases
	if(is.null(chrAliases))
	  chrAliases_C <- " "
	  
	countMultiMappingReads_C <- countMultiMappingReads
	if(primaryOnly) countMultiMappingReads_C <- 2

	read2pos_C <- read2pos
	if(is.null(read2pos)) read2pos_C <- 0

	split_C <- 0
	if(splitOnly) split_C <- 1
	if(nonSplitOnly) split_C <- 2
	  
	cmd <- paste("readSummary",ann,files_C,fout,as.numeric(isPairedEnd),minFragLength,maxFragLength,0,as.numeric(allowMultiOverlap),as.numeric(useMetaFeatures),nthreads,as.numeric(isGTFAnnotationFile),strandSpecific,as.numeric(reportReads),as.numeric(requireBothEndsMapped),as.numeric(!countChimericFragments),as.numeric(checkFragLength),GTF.featureType,GTF.attrType,minMQS,as.numeric(countMultiMappingReads_C),chrAliases_C," ",as.numeric(FALSE),14,readExtension5,readExtension3,minOverlap,split_C,read2pos_C," ",as.numeric(ignoreDup),as.numeric(!autosort),as.numeric(fraction),as.numeric(largestOverlap),PE_orientation,as.numeric(juncCounts),genome_C,maxMOp,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	C_args <- .C("R_readSummary_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

	x <- read.delim(fout,stringsAsFactors=FALSE)
	colnames(x)[1:6] <- c("GeneID","Chr","Start","End","Strand","Length")

	x_summary <- read.delim(paste(fout,".summary",sep=""), stringsAsFactors=FALSE)

	if(juncCounts)
		x_jcounts <- read.delim(paste(fout,".jcounts",sep=""), stringsAsFactors=FALSE)

	file.remove(fout)
	file.remove(paste(fout,".summary",sep=""))

	if(juncCounts)
		file.remove(paste(fout,".jcounts",sep=""))
	
	if(flag) 
	  file.remove(fout_annot)
	
	if(ncol(x) == 6){
	  stop("No count data were generated.")
	}
	
	y <- as.matrix(x[,-c(1:6)])
	colnames(y) <- colnames(x)[-c(1:6)]
	rownames(y) <- x$GeneID

	if(juncCounts)
		z <- list(counts=y,counts_junction=x_jcounts,annotation=x[,1:6],targets=colnames(y),stat=x_summary)
	else
		z <- list(counts=y,annotation=x[,1:6],targets=colnames(y),stat=x_summary)	

	z	
}

