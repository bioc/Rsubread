
.rstest <- function(r, coef) r * (1 + 1/r)^(1 + coef)

.averaging_transform <- function(r, nr){
    d = c(1, r[ 2:length(r) ] - r[1:(length(r)-1)])
    dr = c( 0.5 * (d[2:length(d)] + d[1:(length(d)-1)]),d[ length(d) ])
       nr/dr
}


.simple.Good.Turing.Freq<-function( raw.freq ){
  n_zero <- sum(raw.freq == 0)
  times <- table(raw.freq[raw.freq>0])
  obs <- sort(as.numeric(names(times)))
  times <- times[as.character(obs)]

  nobs = length(obs)
  if(nobs < .min_mySGT_input_size){
    cat(sprintf("Warning: not enough element in the observation array : %d < %d.\nThe `Rescued' element in the returned object will be set to NA.\n",  nobs, .min_mySGT_input_size))
    return(NA)
  }else{
      if(nobs != length(times)) stop("The oservation array doesn't match the frequency array.")
      total_observed <- sum(obs * times)
      tval <- matrix(0, ncol=2, nrow=nobs)
      ks <- c(obs[2:nobs], 2*obs[nobs] - obs[nobs-1])
      tval[,1] <- log(obs)
      tval[,2] <- log( .averaging_transform(obs,times))
      meanx <- mean(tval[,1])
      meany <- mean(tval[,2])
      slope <- sum( (tval[,1]-meanx) * (tval[,2]-meany) ) / sum((tval[,1]-meanx) ^2)
      intercept <- meany - slope * meanx
  
      obs.in.next <- (obs+1) %in% obs
      SD <- ( 1+(1:nobs) ) / times *(c(times[2:nobs], 0)* (1+c(times[2:nobs],0) /times))^.5
      SD [!obs.in.next] <- 1
  
      GT <- (obs+1)/obs * c(times[2:nobs],0)/times
      GT [!obs.in.next] <- 0
      y <- .rstest(obs,slope)
      LGT <- y / obs
  
      xy.sim <- ( abs(GT-LGT) * (1:nobs) / SD ) > 1.65
      min.no.turing = min(which(!xy.sim))
      r <- GT
      if(min.no.turing>0)r[min.no.turing:nobs] <- LGT[min.no.turing:nobs]
      p_zero <- times[1] / total_observed
      rs <- (1-p_zero) * obs * r / sum(r*times*obs/total_observed)
      total_rs <- sum(rs * times)
      puniq <- (1-p_zero) * rs / total_rs
      pres <- rep(p_zero / n_zero, length(raw.freq))
      pres [raw.freq>0] <- puniq[match( raw.freq[raw.freq>0], obs )]
      ret <- list(p0=p_zero, p=pres)
      ret
  }
}

.min_mySGT_input_size <- 5
.smoothed <- function(i,intercept,slope) 2.718281828^(intercept + slope * log(i))

# Simple Good Turing
.mySGT <- function(obs.per.spe){
    obs.tab <- table(obs.per.spe[obs.per.spe>0])
    obs <- sort( as.numeric(names(obs.tab)))
    obs.tab <- obs.tab[ as.character(obs) ]
    #saveRDS(obs.tab,"del4-obs.tab.RDS")
    #saveRDS(obs.per.spe, "del4-obs.per.spe.RDS")
    sgtr <- .mySGTsorted(as.numeric(names(obs.tab)), obs.tab)
    #saveRDS(sgtr, "del4-sgtr.RDS")
    res <- obs.per.spe
    res[obs.per.spe!=0] <- sgtr$p[match(as.character(obs.per.spe[obs.per.spe!=0]), names(obs.tab))]
    n_zero <- sum(obs.per.spe==0)
    res[obs.per.spe==0] <- sgtr$p0/n_zero
    res
}

.mySGTsorted<-function( obs, times ){
    if(any(obs != sort(obs)))stop("Observations have to be sorted.")
    nobs = length(obs)
    if(nobs < .min_mySGT_input_size) stop("Not enough element in the observation array.")
    if(nobs != length(times)) stop("The oservation array doesn't match the frequency array.")
    total_observed <- sum(obs * times)
    p_zero <- 0
    if(1 %in% obs) p_zero <- times[which(1==obs)[1]] / total_observed
    tval <- matrix(0, ncol=2, nrow=nobs)
    ks <- c(obs[2:nobs], 2*obs[nobs] - obs[nobs-1])
    tval[,1] <- log(obs)
    tval[,2] <- log(2*times / (ks - c(0, obs[1:(nobs-1)])))
    meanX <- mean(tval[,1])
    meanY <- mean(tval[,2])
    slope <- sum( (tval[,1]-meanX) * (tval[,2]-meanY) ) / sum((tval[,1]-meanX) ^2)
    intercept <- meanY - slope * meanX

    y <- (1+obs) * .smoothed(1+obs, intercept, slope ) / .smoothed(obs, intercept, slope)
    x <- (1+obs) * c((times[2:nobs] / times[1:(nobs-1)]), y[nobs])
    xy.sim <- abs(x-y) <= 1.96*( (obs+1)^2 * c(times[2:nobs], 0)/ times^2 * (1+ c(times[2:nobs], 0) / times ))^0.5
    obs.in.next <- (obs+1) %in% obs
    min.inv <- min(which( xy.sim | !obs.in.next ))
    r<-x
    if(min.inv>0)r[min.inv:nobs] <- y[min.inv:nobs]
    list(p0=p_zero, p=(1-p_zero) * r / sum(r*times))
}


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

library(Matrix)
.read.sparse.mat <- function (fn){
#  print(fn)
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
      #print(summary(CTRL.LLH))
    }
    ret.nUMI.LLH.tab[[ bcsize_one ]] <- N10000.LLH
    Old_bs_one <- bcsize_one
  }
  ret.nUMI.LLH.tab
}

.is.no.candBC <- function(fn){
  info = file.info(fn)
  return(info$size < 2)
}

.cellCounts.rescue <- function( BAM.name, FC.gene.ids, sample.no ){
  fname <- sprintf("%s.scRNA.%03d", BAM.name, sample.no)
  nozero.anywhere.genes <- read.delim(paste0(fname,".no0Genes"), stringsAsFactors=F, header=F)$V1
  ambient.accumulate <- read.delim(paste0(fname,".AmbSum"), stringsAsFactors=F)
  #saveRDS(list(BAM.name, FC.gene.ids, sample.no), "del4-debug.RDS")
  ambient.accumulate <- ambient.accumulate[ match(FC.gene.ids , ambient.accumulate$GeneID), ]
  ambient.accumulate$UMIs[is.na(ambient.accumulate$UMIs)] <- 0
  ambient.accumulate <- ambient.accumulate$UMIs
  names(ambient.accumulate) <- FC.gene.ids 

  ambient.accumulate <- ambient.accumulate[ names(ambient.accumulate) %in%  nozero.anywhere.genes]

  resc.bc.fn <- paste0(fname,".RescCand")
  rstfq <- .simple.Good.Turing.Freq(ambient.accumulate)
  if(any(is.na(rstfq)) || .is.no.candBC(paste0(resc.bc.fn,".BCtab"))){
    return(NA)
  }else{
    gte <- rstfq$p
    print("Summary of the ambient RNA profile (proportions)")
    print(summary(gte))
    # This function returns "times" log-likelihoods.
  
    rescue.candidates <- as.matrix(.read.sparse.mat(resc.bc.fn))
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
    rescue.candidates[,Rescured.Barcodes, drop=FALSE]
  }
}

.load.one.scSample <- function( BAM.name, FC.gene.ids, sample.no ){
  fname <- sprintf("%s.scRNA.%03d", BAM.name, sample.no)
  highconf <- as.matrix(.read.sparse.mat(paste0(fname,".HighConf")))
  rescued <- .cellCounts.rescue(BAM.name, FC.gene.ids, sample.no)
  if(!any(is.na(rescued))) rescued <- rescued[rowSums(rescued)>0,]
  list( HighConf=highconf, Rescued=rescued )
}

.load.all.scSamples <- function( BAM.name, FC.gene.ids){
  sum.tab <- read.delim(paste0(BAM.name,".scRNA.SampleTable"), stringsAsFactors=F)
  ret <- list()
  for(roiw in 1:nrow(sum.tab)){
    sname <- sum.tab$SampleName[roiw]
    sid <- sum.tab$Index[roiw]
    count.tab <- .load.one.scSample(BAM.name, FC.gene.ids, sid)
    ret[[sprintf("Sample.%03d",sid)]] <- count.tab
  }
  ret[["Sample.Table"]] <- sum.tab

  ret
}

cellCounts <- function(index, input.directory, output.BAM, sample.sheet, cell.barcode.list=NULL, input.mode="BCL", aligner="align", nthreads=16, annot.inbuilt="mm10",annot.ext=NULL,isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",GTF.attrType.extra=NULL,chrAliases=NULL,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,countMultiMappingReads=TRUE){
  input.directory <- .check_and_NormPath(input.directory,  mustWork=T, opt="input.directory")
  sample.sheet <- .check_and_NormPath(sample.sheet,  mustWork=T, opt="sample.sheet")
  output.BAM <- .check_and_NormPath(output.BAM,  mustWork=F, opt="output.BAM")
  if(!is.null(aligner)) aligner <- match.arg(aligner,c("subjunc","align")) 
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
      if(is.null(aligner)){
        if(!file.exists(output.1)) stop("No aligner is specified but the BAM file does not exist. Please specify 'align' or 'subjunc' as the aligner.")
      }else if(aligner=="align"){
        align(index, input.1, output_file=output.1, nthreads=nthreads, isBCLinput=TRUE)
      }else if(aligner=="subjunc"){
        stop("Using subjunc as the aligner is not yet supported.")
      }
      fc[[paste0("counts.", ii)]]<-featureCounts(output.1, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTFAnnotationFile=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, GTF.attrType.extra=GTF.attrType.extra, chrAliases=chrAliases, useMetaFeatures=useMetaFeatures, allowMultiOverlap=allowMultiOverlap, countMultiMappingReads=countMultiMappingReads, sampleSheet=sample.1, cellBarcodeList=cell.barcode.list, nthreads=nthreads)
      fc[[paste0("scRNA.", ii)]] <- .load.all.scSamples(output.1, as.character(fc[[paste0("counts.", ii)]]$annotation$GeneID))
  }
  fc[["Input.Files"]] <- input.directory

  fc
}

