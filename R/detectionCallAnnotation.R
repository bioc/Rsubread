detectionCallAnnotation <- function(species="hg", binsize=2000) {
	if (species=="hg"){
		exonfile <- "hg19_exon.txt"
		irfile <- "hg19_integenic_region.txt"
	}
	if (species == "mm"){
		exonfile <- "mm9_exon.txt"
		irfile <- "mm9_integenic_region.txt"
	}
	
	if ((file.exists(exonfile) == FALSE) || (file.exists(irfile) == FALSE)){

		print("Source files specified doesn't exist!")

	} 

	else 

	{
	
		.C("detectionCallAnnotation", as.character(exonfile), as.character(irfile), as.character(species), as.integer(binsize), PACKAGE="Rsubread")
	
	}
	
	NULL

}

