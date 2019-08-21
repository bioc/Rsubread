cellCounts <- function(input.directory, index, output.BAM, sample.sheet, cell.barcode.list, input.mode="BCL", nthreads=16, annot.inbuilt="mm10",annot.ext=NULL,isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",GTF.attrType.extra=NULL,chrAliases=NULL,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,countMultiMappingReads=FALSE){

  fc <- list()

  for(ii in 1:length(input.directory)){
	  input.1 <- input.directory[ii]
	  output.1 <- output.BAM[ii]
	  sample.1 <- sample.sheet[ii]
	  align(index, input.1, output_file=output.1, nthread=nthreads, isBCLinput=TRUE)
	  fc[[paste0("counts.", ob)]]<-featureCounts(output.1, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTFAnnotationFile=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, GTF.attrType.extra=GTF.attrType.extra, chrAliases=chrAliases, useMetaFeatures=useMetaFeatures, allowMultiOverlap=allowMultiOverlap, countMultiMappingReads=countMultiMappingReads, sampleSheet=sample.1, cellBarcodeList=cell.barcode.list)
  }
  for(ob in output.BAM){
    fc[[paste0("scRNA.table.", ob)]] <- read.delim(paste0(ob,".scRNA.table"), header=T, stringsAsFactors=F)
  }
  fc
}
