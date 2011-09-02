detectionCall <- function(dataset, species="hg", plot=FALSE)

{
	if (tolower(substring(dataset, nchar(dataset)-3, nchar(dataset))) == ".sam"){
		dataset <- substring(dataset, 1, nchar(dataset)-4)
	}
	
	samfile <- paste(dataset,".sam",sep="")
	
	if (file.exists(samfile) == FALSE) {

		stop("SAM file for the given dataset doesn't exist!\n")

	} 
	
	if (species == "mm"){
		
		### mouse annotation file not available yet
		exon_file <- system.file("annot","mm9_NCBI_exon_GC.txt",package="Rsubread")
		ir_file <- system.file("annot","mm9_NCBI_binned_integenic_region_GC.txt",package="Rsubread") 
			
	} 
	
	if(species == "hg"){
			
		exon_file <- system.file("annot","hg19_NCBI_exon_GC.txt",package="Rsubread")
		ir_file <- system.file("annot","hg19_NCBI_binned_integenic_region_GC.txt",package="Rsubread") 
			
	}

	temp_header <- 	"/tmp"

	temp_header <- paste(temp_header,paste("/.Rsubread_pid",Sys.getpid(),sep=""), sep="")
	
	temp_header <- paste(temp_header, "_", sep="")
	
	temp_header <- paste(temp_header, dataset,sep="")
	
	.C("detectionCall", as.character(dataset), as.character(exon_file), as.character(ir_file), as.character(temp_header), PACKAGE="Rsubread")

	##########################################################
	### After processing SAM file and Annotation file in C####
	###  The result will be processed here for statistics ####
	##########################################################

	### ALL CHROMOSOMES
	chrs <- c("chr1", "chr10", "chr11","chr12", "chr13","chr14","chr15","chr16", "chr17","chr18","chr19", "chr2","chr20","chr21","chr22","chr3","chr4", "chr5","chr6", "chr7", "chr8","chr9", "chrX","chrY")


	bg_file <- paste(temp_header, "_mapping_binned_integenic_region_GC.txt", sep="")
	signal_file <- paste(temp_header, "_mapping_exon_GC.txt", sep="")

	plot_output_file <- paste(dataset, "_p_value_density_plot.pdf", sep="")


	if ((file.exists(bg_file) == FALSE) || (file.exists(signal_file) == FALSE)) {

		print("Mapping results for this dataset are not available.\n")
		print("Please execute detectionCall(dataset) first. \n")
		
	} else {


		######################################################
		########## INTEGENIC REGION EXPRESSION DATA ########## 
		######################################################

		bg <- read.delim(bg_file, header=TRUE)
		bg$percent <- bg$gcnum / (bg$chr_stop - bg$chr_start + 1 - bg$nnum)
		bg <- bg[((bg$chr_stop - bg$chr_start +1) == 2000),]
		bg <- bg[(bg$nnum + bg$gcnum + bg$atnum) >= 1990,]
		bg$gccontent = round(bg$percent*100)


		######	MEAN VS GC CONTENT CALCULATION

		mean_nreads <- rep(0,101)
		corrected_mean_nreads <- rep(0,101)
		sd <- rep(0,101)
		num_entries <- rep(0,101)

		for (i in 0:100){
			bg_percent <- bg[bg$gccontent == i, ]
			mean_nreads[i+1] <- mean(bg_percent$nreads)
			sd[i+1] <- sd(bg_percent$nreads)	
			selected_bg_percent <- bg_percent[bg_percent$nreads >= (mean_nreads[i] - 2*sd[i]),]	
			selected_bg_percent <- selected_bg_percent[selected_bg_percent$nreads <= (mean_nreads[i] + 2*sd[i]),]
			corrected_mean_nreads[i+1] <- mean(selected_bg_percent$nreads)
			num_entries[i+1] <- nrow(bg_percent)
		}	

		bg_intensities <- corrected_mean_nreads
		### bg_intensities is an array of 101 elements, recording the mean n reads for each GC content(outliers removed)

		##########################################
		########## GENE EXPRESSION DATA ########## 
		##########################################

		exon <- read.delim(signal_file, header=TRUE)
		exon <- exon[exon$chr %in% chrs, ]
		exon$length <- exon$chr_stop - exon$chr_start + 1


		gene_id <- sort(unique(exon$entrezid))
		gene_nreads <- tapply(exon$nreads, factor(exon$entrezid), sum)
		gene_gcnum <- tapply(exon$gcnum, factor(exon$entrezid), sum)
		gene_length <- tapply(exon$length, factor(exon$entrezid), sum)
		gene_nnum <- tapply(exon$nnum, factor(exon$entrezid), sum)
		gene_percent <- gene_gcnum / (gene_length - gene_nnum)
		### gene_signal is the adjusted gene nreads per 2000 base
		gene_signal <- gene_nreads / gene_length * 2000

		gene <- cbind(matrix(gene_id), matrix(gene_length),matrix(gene_percent),  matrix(gene_signal))
		gene <- as.data.frame(gene)
		colnames(gene) <- c("entrezid", "length","GC", "signal")

		### For each gene, get the read intensity of same GC content, then calculate its p_values

		gene$gccontent <- round(gene$GC*100)
		gene_bg <- rep(0, nrow(gene))

		for (i in 1:nrow(gene)){
			this_gccontent <- gene$gccontent[i]
			this_bg_intensity <- bg_intensities[this_gccontent+1]
			if (is.na(this_bg_intensity)) {this_bg_intensity <- 0}
			gene_bg[i] <- this_bg_intensity
		}
		gene$gccontent <- NULL
		gene$background <- gene_bg

		gene$p_value <- ppois(gene$background, gene$signal)

		###########################################################
		########## DENSITY PLOT OF P_VALUE FOR EACH GENE ########## 
		###########################################################

		if (plot){
			pdf(plot_output_file )
			plot(density(gene$p_value), xlab = "p_value", ylab = "density", main="Density Plot of P_values for Each Gene(Binsize=2000)", col=1)
			dev.off()
		}

		##############################################
		########## CLEAN UP AND DATA EXPORT ########## 
		##############################################
		
		file.remove(bg_file)
		file.remove(signal_file)
		return(gene)
		
	}
}

