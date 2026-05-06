
bam2junctions <- function(BAMfile, cell.barcodes){
   temp_prefix <- file.path(".",paste(".Rsubread_featureCounts_pid",Sys.getpid(),sep=""))

   BAMfile <- .check_and_NormPath(BAMfile, mustWork=T, opt="BAM file")
   cellBC.file <- paste0(temp_prefix,".wanted.cellbc")
   write.table(data.frame(cell.barcodes)[,1,drop=F], cellBC.file , row.names=F, col.names=F, quote=F)

   n <- 3
   cmd <- paste(c(cellBC.file , temp_prefix, BAMfile), collapse=.R_param_splitor) 
   .C("R_extract_junction_from_BAM",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
   mtx.file <- paste0(temp_prefix,".mtx")
   junc.file <- paste0(temp_prefix,".junctions.tsv")
   if(!(file.exists( mtx.file  )&&file.exists( junc.file ))) stop("Unable to parse the BAM file.")
   mtx <- Matrix::readMM(mtx.file)
   junc.list <- read.delim(junc.file, header=F)
   rownames(mtx) <- junc.list [order(as.numeric(junc.list[,1])),2]
   colnames(mtx) <- cell.barcodes
   file.remove(c( mtx.file , junc.file,  cellBC.file ))
   mtx
}
