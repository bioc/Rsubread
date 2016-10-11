detectionCall <- function(dataset, species="hg", plot=FALSE)
{
	
	if (file.exists(dataset) == FALSE) {
		stop("File for the given dataset doesn't exist!\n")
	} 
	
	if (species == "mm"){
		exon_file <- system.file("annot","mm9_NCBI_exon_GC.txt",package="Rsubread")
        if(exon_file == ""){
            ret_df <- download.file("http://bioinf.wehi.edu.au/software/Rsubread/mm9_NCBI_exon_GC.txt",file.path(system.file("annot",package="Rsubread"),"mm9_NCBI_exon_GC.txt"))
            if(ret_df != 0)
                stop("Failed to download file http://bioinf.wehi.edu.au/software/Rsubread/mm9_NCBI_exon_GC.txt")
            exon_file <- system.file("annot","mm9_NCBI_exon_GC.txt",package="Rsubread")
        }
        
		ir_file <- system.file("annot","mm9_NCBI_binned_integenic_region_GC.txt",package="Rsubread")
        if(ir_file == ""){
            ret_df <- download.file("http://bioinf.wehi.edu.au/software/Rsubread/mm9_NCBI_binned_integenic_region_GC.txt",file.path(system.file("annot",package="Rsubread"),"mm9_NCBI_binned_integenic_region_GC.txt"))
            if(ret_df != 0)
                stop("Failed to download file http://bioinf.wehi.edu.au/software/Rsubread/mm9_NCBI_binned_integenic_region_GC.txt")
            ir_file <- system.file("annot","mm9_NCBI_binned_integenic_region_GC.txt",package="Rsubread")
        }

	} 
	
	if(species == "hg"){
		exon_file <- system.file("annot","hg19_NCBI_exon_GC.txt",package="Rsubread")
        if(exon_file == ""){
            ret_df <- download.file("http://bioinf.wehi.edu.au/software/Rsubread/hg19_NCBI_exon_GC.txt",file.path(system.file("annot",package="Rsubread"),"hg19_NCBI_exon_GC.txt"))
            if(ret_df != 0)
                stop("Failed to download file http://bioinf.wehi.edu.au/software/Rsubread/hg19_NCBI_exon_GC.txt")
            exon_file <- system.file("annot","hg19_NCBI_exon_GC.txt",package="Rsubread")
        }
        
		ir_file <- system.file("annot","hg19_NCBI_binned_integenic_region_GC.txt",package="Rsubread")
        if(ir_file == ""){
            ret_df <- download.file("http://bioinf.wehi.edu.au/software/Rsubread/hg19_NCBI_binned_integenic_region_GC.txt",file.path(system.file("annot",package="Rsubread"),"hg19_NCBI_binned_integenic_region_GC.txt"))
            if(ret_df != 0)
                stop("Failed to download file http://bioinf.wehi.edu.au/software/Rsubread/hg19_NCBI_binned_integenic_region_GC.txt")
            ir_file <- system.file("annot","hg19_NCBI_binned_integenic_region_GC.txt",package="Rsubread")
        }

	}

	temp_header <- 	"."

	temp_header <- paste(temp_header,paste("/.Rsubread_pid",Sys.getpid(),sep=""), sep="")
	
	temp_header <- paste(temp_header, "_", sep="")
	
	temp_header <- paste(temp_header, dataset,sep="")
	
	##temp_header <- dataset
	
	
	.C("detectionCall", as.character(dataset), as.character(exon_file), as.character(ir_file), as.character(temp_header), PACKAGE="Rsubread")
	
	##########################################################
	### After processing SAM file and Annotation file in C####
	###  The result will be processed here for statistics ####
	##########################################################

	### ALL CHROMOSOMES
	chrs <- c("chr1", "chr10", "chr11","chr12", "chr13","chr14","chr15","chr16", "chr17","chr18","chr19", "chr2","chr20","chr21","chr22","chr3","chr4", "chr5","chr6", "chr7", "chr8","chr9", "chrX","chrY")


	bg_file <- paste(temp_header, "_mapping_binned_integenic_region_GC.txt", sep="")
	signal_file <- paste(temp_header, "_mapping_exon_GC.txt", sep="")

	if ((file.exists(bg_file) == FALSE) || (file.exists(signal_file) == FALSE)) {

		print("Mapping results for this dataset are not available.\n")
		
	} else {

		######################################################
		########## INTEGENIC REGION EXPRESSION DATA ########## 
		######################################################

		bg <- read.delim(bg_file, header=TRUE)
		
		bg$percent <- bg$gcnum / (bg$chr_stop - bg$chr_start + 1 - bg$nnum)
		
		bg$sqrtnreads <- sqrt(bg$nreads)
		
		lowessfit <- lowess(x=bg$percent, y=bg$sqrtnreads)
		
		lowessfun <- approxfun(x=lowessfit$x, y=lowessfit$y)

		##########################################
		########## GENE EXPRESSION DATA ########## 
		##########################################

		exon <- read.delim(signal_file, header=TRUE)
		exon <- exon[exon$chr %in% chrs, ]
		exon$length <- exon$chr_stop - exon$chr_start + 1
		
		# filter out exons which are covered by other exons
		exon <- exon[exon$nnum + exon$gcnum + exon$atnum >= 0.9*exon$length, ]
		
		# group data into gene level
		gene_id <- sort(unique(exon$entrezid))
		
		gene_nreads <- tapply(exon$nreads, factor(exon$entrezid), sum)
		
		gene_length <- tapply(exon$length, factor(exon$entrezid), sum)
		
		gene_gcnum <- tapply(exon$gcnum, factor(exon$entrezid), sum)
		gene_nnum <- tapply(exon$nnum, factor(exon$entrezid), sum)
		
		gene_percent <- gene_gcnum / (gene_length - gene_nnum)

		gene <- cbind(matrix(gene_id), matrix(gene_length),matrix(gene_percent), matrix(gene_nreads))
		gene <- as.data.frame(gene)
		
		colnames(gene) <- c("entrezid", "length","GC","nreads")

		### For each gene, get the read intensity of same GC content, then calculate its p_values

		# predict the background intensity of a gene using previously fitted lowess function
		# and translate it back from sqrt scale to normal scale
		gene$background <- (lowessfun(gene$GC))^2
		
		gene$background[is.na(gene$background)] <- 0
		
		# calculate the expected nreads under null hypothesis for each gene
		gene$expected_nreads <- gene$background*gene$length / 2000
		
		gene$p_value <- ppois(gene$nreads, gene$expected_nreads, lower.tail=FALSE) + runif(nrow(gene))*dpois(gene$nreads, gene$expected_nreads)
	
		###########################################################
		########## DENSITY PLOT OF P_VALUE FOR EACH GENE ########## 
		###########################################################

		if (plot){
			plot_output_file <- paste(dataset, "_p_value_density_plot.pdf", sep="")
			pdf(plot_output_file)
			plot(density(gene$p_value), xlab = "p_value", ylab = "density", main="Density Plot of P_values for Each Gene(Binsize=2000)", col=1)
			dev.off()
		}


		##############################################
		########## CLEAN UP AND DATA EXPORT ########## 
		##############################################
		
		unlink(bg_file)
		unlink(signal_file)
		return(gene)
		
	}
}

