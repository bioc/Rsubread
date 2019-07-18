promoterRegions <- function(annotation="mm10", upstream=3000L, downstream=2000L)
#	Create a SAF data.frame of genewise promoter regions
#	Gordon Smyth
#	24 April 2017
{
    .check_string_param(annotation,'annotation')
#	annotation can be a SAF format data.frame or can be the name of a genome with built-in annotation
	if(is.character(annotation)) annotation <- getInBuiltAnnotation(annotation)

#	Remove unassembled contigs
	N <- grep("^N",annotation$Chr)
	annotation <- annotation[-N,]

#	Check upstream and downstream limits
	upstream <- as.integer(upstream[1])
	if(upstream < 0L) upstream <- 0L
	downstream <- as.integer(downstream[1])
	if(downstream < 0L) downstream <- 0L

#	Combine Chr and GeneID
	annotation$ChrGeneID <- paste(annotation$Chr,annotation$GeneID,sep=".")

#	Get start of each gene
	o <- order(annotation$ChrGeneID,annotation$Start)
	anno.start <- annotation[o,]
	isdup <- duplicated(anno.start$ChrGeneID)
	exon.first <- which(!isdup)
	anno.start <- anno.start[exon.first,]

#	Get end of each gene
	o <- order(annotation$ChrGeneID,annotation$End)
	anno.end <- annotation[o,]
	exon.last <- c(exon.first[-1]-1L,nrow(annotation))
	anno.end <- anno.end[exon.last,]

#	Assemble start and end into one data.frame
	ann.gene <- anno.start
	ann.gene$End <- anno.end$End

#	Compute promoter region for + strand
	prom.start <- pmax(ann.gene$Start-upstream,1L)
	prom.end <- pmin(ann.gene$Start+downstream,ann.gene$End)

#	Compute promoter region for - strand
	neg <- ann.gene$Strand=="-"
	prom.start[neg] <- pmax(ann.gene$End[neg]-downstream,ann.gene$Start[neg])
	prom.end[neg] <- ann.gene$End[neg]+upstream

#	Output
	ann.gene$Start <- prom.start
	ann.gene$End <- prom.end
	ann.gene$ChrGeneID <- NULL
	o <- order(ann.gene$Chr,ann.gene$Start)
	ann.gene[o,]
}
