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
	sprintf("http://bioinf.wehi.edu.au/cellCounts/cell-barcodes/%s",uri)
}

.find_best_cellbarcode <- function( input.directory, sample.sheet){
	barcode.database.file <- path.expand("~/.Rsubread/cellCounts/known_barcode_sets.txt")
	if(!file.exists(barcode.database.file)){
		dir.create("~/.Rsubread/cellCounts", recursive=TRUE)
		rr <- download.file(.dataURL("known_barcode_sets.txt"), barcode.database.file)
		if(rr!=0)stop("ERROR: the barcode database cannot be retrieved from the Internet. You may still run cellCounts by specifying a local barcode list file to the `cell.barcode.list` option.")
	}
	bcb.sets <- read.delim(barcode.database.file, stringsAsFactors=F, header=T)
	cat(sprintf("Found %d known cell barcode sets.\n", nrow(bcb.sets)))
	max.cell.good <- -1
	for(libf in bcb.sets$File){
		listfile <- path.expand(paste0("~/.Rsubread/cellCounts/",libf))
		if(!file.exists(listfile)){
			rr <- download.file(.dataURL(libf), listfile)
			if(rr!=0)stop("ERROR: the barcode list cannot be retrieved from the Internet. You may still run cellCounts by specifying a local barcode list file to the `cell.barcode.list` option.")
		}
		cat(sprintf("Testing the cell barcodes in %s.\n", libf))
		barcode_res <- .cellCounts_try_cellbarcode(input.directory[1], sample.sheet[1], listfile, 30000)
		if(length(barcode_res)<3)stop("ERROR: the input sample cannot be processed.")
		sample.good.rate <- barcode_res[2]/barcode_res[1]
		cell.good.rate <- barcode_res[3]/barcode_res[1]
		max.cell.good <- max(max.cell.good, cell.good.rate)
		#cat(sprintf("Sample supporting rate : %.1f%% ; cell supporting rate : %.1f%%.\n", sample.good.rate*100., cell.good.rate*100.))
		if(sample.good.rate < 0.5)cat(sprintf("WARNING: there are only %.1f%% reads having known sample indices. Please check if the sample sheet is correct.\n", sample.good.rate*100.))
		if(cell.good.rate > 0.6){
			cat(sprintf("Found cell-barcode list '%s' for the input data: supported by %.1f%% reads.\n", libf, cell.good.rate*100.))
			return(listfile)
		}
	}

	stop(sprintf("ERROR: no known cell barcode set was found for the data set. The highest percentage of cell-barcode matched reads is %.1f%%\n", max.cell.good*100.))
}

.read.sparse.mat <- function (fn){
  library(Matrix)
  mtx <- readMM(paste0(fn, ".spmtx"))
  coln <- read.delim(paste0(fn, ".BCtab"), stringsAsFactors=F, header=F)$V1
  rown <- read.delim(paste0(fn, ".GENEtab"), stringsAsFactors=F, header=F)$V1
  colnames(mtx) <- coln
  rownames(mtx) <- rown

  mtx
}

.simu.multinomial <- function(candi.mat, gene.profile.freq , times=10000){
  bcsizes <- sort(unique( colSums(candi.mat) ))
  ret.nUMI.LLH.tab <- list()
  N10000.LLH <- NA
  Old_bs_one <- NA
  log_E_GTE <- log(gene.profile.freq)
  for(bcsize_one in bcsizes){
  UMI_step_diff <- NA
  if(!any(is.na(Old_bs_one))) UMI_step_diff <- bcsize_one - Old_bs_one

  if(any(is.na(N10000.LLH))){
      N10000 <- rmultinom(n=times, size=bcsize_one, prob=gene.profile.freq )
    N10000.LLH <- apply(N10000, 2, function(x) dmultinom(x, prob=gene.profile.freq, log =T ))
  }else if(UMI_step_diff >= 1000){
    UMI_step_diff_N10000 <- rmultinom(n=times, size=UMI_step_diff, prob=gene.profile.freq )
    N10000 <- N10000 + UMI_step_diff_N10000
    N10000.LLH <- apply(N10000, 2, function(x) dmultinom(x, prob=gene.profile.freq, log =T ))
  }else{
    for(curi in (Old_bs_one+1):bcsize_one){
    UMI_step_indices_N10000 <- sample(1:nrow(candi.mat), size=times, prob=gene.profile.freq, replace=T)
    for(idi in 1:times){
      N10000[ UMI_step_indices_N10000[idi] ,idi ] <- N10000[ UMI_step_indices_N10000[idi] ,idi ]+1
    }
    UMI_step_values_N10000 <- rep(0, times)
    for(idi in 1:times){
      UMI_step_values_N10000[idi] <- N10000[UMI_step_indices_N10000[idi] ,idi ]
    }
    #print(UMI_step_values_N10000)
    N10000.LLH <- N10000.LLH + log_E_GTE[ UMI_step_indices_N10000 ] + log((curi +1.) / UMI_step_values_N10000)
    }
  }

  # cat("\n\n ++++++++++ ", bcsize_one," +++++++++ \n")
  # print(summary(N10000.LLH))
  do_CTRL <- F
  if(do_CTRL && 0== bcsize_one %% 15){
    CTRL <- rmultinom(n=times, size=bcsize_one, prob=gene.profile.freq )
    CTRL.LLH <- apply(CTRL, 2, function(x) dmultinom(x, prob=gene.profile.freq, log =T ))
    print(summary(CTRL.LLH))
  }
  ret.nUMI.LLH.tab[[ bcsize_one ]] <- N10000.LLH
  Old_bs_one <- bcsize_one
  }
  ret.nUMI.LLH.tab
}

.cellCounts.rescue <- function( BAM.name, FC.gene.ids, sample.no ){
  fname <- sprintf("%s.scRNA.%03d", BAM.name, sample.no)
  nozero.anywhere.genes <- read.delim(paste0(fname,".no0Genes"), stringsAsFactors=F, header=F)$V1
  ambient.accumulate <- read.delim(paste0(fname,".AmbSum"), stringsAsFactors=F)
  ambient.accumulate <- ambient.accumulate[ match(FC.gene.ids , ambient.accumulate$GeneID), ]
  ambient.accumulate$UMIs[is.na(ambient.accumulate$UMIs)] <- 0
  ambient.accumulate <- ambient.accumulate$UMIs
  names(ambient.accumulate) <- FC.gene.ids 

  library(edgeR)
  ambient.accumulate <- ambient.accumulate[ names(ambient.accumulate) %in%  nozero.anywhere.genes]
  gte <- goodTuringProportions (ambient.accumulate)


  # This function returns "times" log-likelihoods.

  rescue.candidates <- as.matrix(.read.sparse.mat(paste0(fname,".RescCand")))
  rescue.candidates <- rescue.candidates[ match(names(ambient.accumulate), rownames(rescue.candidates)), ]
  rownames(rescue.candidates) <- names(ambient.accumulate)
  rescue.candidates [is.na(rescue.candidates )] <- 0
  # Compute observed log-likelihood of barcodes being generated from ambient RNA
  # Compute the multinomial log PMF for many barcodes -- log-likelihoods
  # Only use the "non-zero anywhere" genes.
  log.like.cands <- apply(rescue.candidates, 2, function(x) dmultinom(x, prob=gte, log =T ))

  # Simulate log likelihoods
  # This step is to build like 10000 sets of N_GENES vectors from the multinomial distribution of "gte"
  # See what are their log-likelihoods against "gte" -- most should be very small.

  simu.pvalues <- .simu.multinomial( rescue.candidates, gte )

  # Compute p-values : the p-value is the chance of a barcode is actually from ambient RNA (ie smaller the p-value, more likely this barcode is for a real cell)
  actual.pvalues <- rep(0, ncol(rescue.candidates) )
  names(actual.pvalues) <- colnames(rescue.candidates)
  for(candi in names(log.like.cands)){
    cand_umis <- sum(rescue.candidates[,candi])
    cand.simu.pvs <- simu.pvalues[[ cand_umis ]]
    cand.actual.pv <- log.like.cands[candi]
    actual.pvalues[candi] <- sum(cand.simu.pvs > cand.actual.pv) / length(cand.simu.pvs)
  }

  # p-value => FDR
  actual.FDR <- p.adjust(actual.pvalues, method='BH')

  # select cells that has FDR < cutoff
  FDR.Cutoff <- 0.01
  Rescured.Barcodes <- names(actual.FDR)[actual.FDR <= FDR.Cutoff]
  rescue.candidates[,Rescured.Barcodes]
}

.load.one.scSample <- function( BAM.name, FC.gene.ids, sample.no ){
  fname <- sprintf("%s.scRNA.%03d", BAM.name, sample.no)
  highconf <- as.matrix(.read.sparse.mat(paste0(fname,".HighConf")))
  rescued <- .cellCounts.rescue(BAM.name, FC.gene.ids, sample.no)
  cbind( highconf,rescued )
}

.load.all.scSamples <- function( BAM.name, FC.gene.ids){
  sum.tab <- read.delim(paste0(BAM.name,".SampleTable"), stringsAsFactors=F)
  ret <- list()
  for(roiw in 1:nrow(sum.tab)){
    sname <- sum.tab$SampleName
    sid <- sum.tab$Index
    count.tab <- .load.one.scSample(BAM.name, FC.gene.ids, sid)
	ret[[sname]] <- count.tab
  }

  ret
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
	  fc[[paste0("scRNA.", output.1)]] <- .load.all.scSamples(output.1, fc[[paste0("counts.", output.1)]]$annotation$GeneID)
  }

  fc
}
