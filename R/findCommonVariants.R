findCommonVariants <- function(VCF_files)
{
	VCF_files <- .check_and_NormPath(VCF_files, mustWork=T, opt="VCF_files")
	fout <- file.path(".",paste(".Rsubread_featureCounts_pid",Sys.getpid(),sep=""))

	files_C <- paste(VCF_files,collapse=";")

	cmd <- paste('-o', fout, files_C, sep=';')

	n <- 2 + length(VCF_files)

	C_args <- .C("R_mergeVCF",as.integer(n),as.character(cmd),PACKAGE="Rsubread")

	x <- try(read.delim(fout,stringsAsFactors=FALSE, comment.char="#"),silent=TRUE)
	file.remove(fout)
	if(class(x)!='try-error') {
		colnames(x)[1:8] <- c("Chr","Pos","Reference","Identity","Alternative","Quality", "Filter","Info")
		x
	}
	else{
		NA
	}

}

.findCommonVariants_inner <- function(VCF_files)
{
	
	for(vcf_file in VCF_files)
	{
		print(paste("Doing", vcf_file))
		vcf_body <- read.delim(vcf_file, comment.char='#', header=FALSE)
		
		nky <- 0
		for(row_no in 1:nrow(vcf_body))
		{
			alts <- strsplit(as.vector(vcf_body$V5[row_no]), ',')
			alts <- alts[[1]]
			nky <- nky + length(alts)
		}

		print(paste("NROW =", nky))
		
		data_mat <- matrix(ncol = 3, nrow = nky)

		n2 <- 1
		for(row_no in 1:nrow(vcf_body))
		{
			vcf_row <- vcf_body[row_no,]
			key_values <- strsplit(as.vector(vcf_row$V8), ';')
			key_values <- key_values[[1]]
			key_type <- ifelse('INDEL' %in% key_values, 'INDEL','SNP')
			alts <- strsplit(as.vector(vcf_row$V5), ',')
			alts <- alts[[1]]

			chro <- vcf_row$V1
			pos <- vcf_row$V2
			qual <- vcf_row$V6
			for(alt in alts)
			{
				ky = paste(key_type, chro, pos, alt)
				data_mat[n2,1] <- ky
				data_mat[n2,2] <- qual
				data_mat[n2,3] <- as.vector(vcf_row$V8)
				n2<-n2+1
			}
			print(paste("CUR_ROW" , n2))
				n2<-n2+1
		}

	}
}
