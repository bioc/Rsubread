.cellCounts_try_cellbarcode <- function( input.directory, sample.sheet, cell.barcode.list, nreads.testing ){ # the three parameters can only be one string, not strings!
    input.directory <- .check_and_NormPath(input.directory)
	sample.sheet <- .check_and_NormPath(sample.sheet)
	cell.barcode.list <- .check_and_NormPath(cell.barcode.list)
	cmd <- paste0(c(input.directory, sample.sheet, cell.barcode.list, as.character(nreads.testing)), collapse=.R_param_splitor)
	rvs <- as.integer(rep(0,5))
	C_args <- .C("R_try_cell_barcode_wrapper",nargs=as.integer(4),argv=as.character(cmd),retv=rvs ,PACKAGE="Rsubread")
	return_val <- ifelse(C_args$retv[1]==0,"FINISHED","ERROR")
	tested_reads <- C_args$retv[2]
	good_sample <- C_args$retv[3]
	good_cell <- C_args$retv[4]
	if(return_val=="ERROR"){
		return(NA)
	}
	return(c(tested_reads, good_sample, good_cell))
}

.find_best_cellbarcode <- function( input.directory, sample.sheet){
}

cellCounts <- function(index, input.directory, output.BAM, sample.sheet, cell.barcode.list, input.mode="BCL", nthreads=16, annot.inbuilt="mm10",annot.ext=NULL,isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",GTF.attrType.extra=NULL,chrAliases=NULL,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,countMultiMappingReads=FALSE){

  fc <- list()

  if(length(input.directory) != length(output.BAM) || length(input.directory) != length( sample.sheet ))stop("The arguments to the input.directory, output.BAM and sample.sheet options must have the same length.")

  for(ii in 1:length(input.directory)){
	  input.1 <- input.directory[ii]
	  output.1 <- output.BAM[ii]
	  sample.1 <- sample.sheet[ii]
	  align(index, input.1, output_file=output.1, nthreads=nthreads, isBCLinput=TRUE)
	  fc[[paste0("counts.", output.1)]]<-featureCounts(output.1, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTFAnnotationFile=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, GTF.attrType.extra=GTF.attrType.extra, chrAliases=chrAliases, useMetaFeatures=useMetaFeatures, allowMultiOverlap=allowMultiOverlap, countMultiMappingReads=countMultiMappingReads, sampleSheet=sample.1, cellBarcodeList=cell.barcode.list, nthreads=nthreads)
      fc[[paste0("scRNA.table.", output.1)]] <- read.delim(paste0(output.1,".scRNA.table"), header=T, stringsAsFactors=F)
  }

  fc
}
