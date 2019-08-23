.cellCounts_try_cellbarcode <- function( input.directory, sample.sheet, cell.barcode.list, nreads.testing ){ # the three parameters can only be one string, not strings!
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

.dataURL <- function(uri){
	sprintf("http://bioinf.wehi.edu.au/Rsubread/%s",uri)
}

.find_best_cellbarcode <- function( input.directory, sample.sheet){
	barcode.database.file <- path.expand("~/.Rsubread/cellCounts/known_barcode_sets.txt")
	if(!file.exists(barcode.database.file)){
		dir.create("~/.Rsubread/cellCounts", recursive=TRUE)
		rr <- download.file(.dataURL("cellCounts/known_barcode_sets.txt"), barcode.database.file)
		if(rr!=0)stop("ERROR: the barcode database cannot be retrieved from the Internet. You may still run cellCounts by specifying a local barcode list file to the `cell.barcode.list` option.")
	}
	bcb.sets <- read.delim(barcode.database.file, stringsAsFactors=F, header=T)
	max.cell.good <- -1
	for(libf in bcb.sets$File){
		listfile <- path.expand(paste0("~/.Rsubread/cellCounts/",libf))
		if(!file.exists(listfile)){
			rr <- download.file(.dataURL(paste0("cellCounts/", libf)), listfile)
			if(rr!=0)stop("ERROR: the barcode list cannot be retrieved from the Internet. You may still run cellCounts by specifying a local barcode list file to the `cell.barcode.list` option.")
		}
		barcode_res <- .cellCounts_try_cellbarcode(input.directory, sample.sheet, listfile, 30000)
		if(length(barcode_res)<3)stop("ERROR: the input sample cannot be processed.")
		sample.good.rate <- barcode_res[2]/barcode_res[1]
		cell.good.rate <- barcode_res[3]/barcode_res[1]
		max.cell.good <- max(max.cell.good, cell.good.rate)
		if(sample.good.rate < 0.5)cat(sprintf("WARNING: there are only %.1f%% reads having known sample indices. Please check if the sample sheet is correct.\n", sample.good.rate*100.))
		if(cell.good.rate > 0.6){
			cat(sprintf("Found cell-barcode list '%s' for the input data: supported by %.1f%% reads.\n", libf, cell.good.rate*100.))
			return(listfile)
		}
	}

	stop(sprintf("ERROR: no known cell barcode set was found for the data set. The highest percentage of cell-barcode matched reads is %.1f%%\n", max.cell.good*100.))
}

cellCounts <- function(index, input.directory, output.BAM, sample.sheet, cell.barcode.list=NULL, input.mode="BCL", nthreads=16, annot.inbuilt="mm10",annot.ext=NULL,isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",GTF.attrType.extra=NULL,chrAliases=NULL,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,countMultiMappingReads=FALSE){
  input.directory <- .check_and_NormPath(input.directory,  mustWork=T, opt="input.directory")
  sample.sheet <- .check_and_NormPath(sample.sheet,  mustWork=T, opt="sample.sheet")
  output.BAM <- .check_and_NormPath(output.BAM,  mustWork=F, opt="output.BAM")

  fc <- list()

  if(length(input.directory) != length(output.BAM) || length(input.directory) != length( sample.sheet ))stop("The arguments to the input.directory, output.BAM and sample.sheet options must have the same length.")
  if(is.null(cell.barcode.list)){
    cell.barcode.list <- .find_best_cellbarcode(input.directory, sample.sheet)
  }else{
    cell.barcode.list <- .check_and_NormPath(cell.barcode.list, mustWork=T, opt="cell.barcode.list")
  }

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
