processExons <- function(filename = "human_seq_gene.md", species="hg")

{
	if (species == "mm"){
			filename <- "mouse_seq_gene.md"
	}
		
	if (file.exists(filename) == FALSE){

		print("Souce file specified doesn't exist!")

	} 

	else 

	{
		x <- read.delim(filename, stringsAsFactors=FALSE)
		if (species == "hg"){
			y <- x[x$group_label=="GRCh37.p2-Primary Assembly" | x$group_label=="GRCh37.p2-non-nuclear", ]
		} else {
			if (species == "mm"){
				y <- x[x$group_label=="MGSCv37-C57BL/6J",]
			}
		}
		y <- y[y$feature_type=="CDS" | y$feature_type=="UTR",]
		y <- y[,c("chromosome","chr_start","chr_stop","chr_orient","feature_id")]

		if (species == "hg"){
			chr_max <- 22
		} else {
			chr_max <- 19
		}
		for(i in 1:chr_max){
			index <- y$chromosome %in% i
			y$chromosome[index] <- rep(paste("chr",i,sep=""),sum(index))
		}
		
		index <- y$chromosome %in% "X"
		y$chromosome[index] <- rep("chrX",sum(index))
		index <- y$chromosome %in% "Y"
		y$chromosome[index] <- rep("chrY",sum(index))


		###NT <- unique(y$chromosome[grep("NT_",y$chromosome)])
		###for(i in NT){
		###	index <- y$chromosome %in% i
		###	i1 <- unlist(strsplit(i,"|",fixed=TRUE))
		###	i1 <- i1[2]
		###	y$chromosome[index] <- rep(i1,sum(index))
		###}
		z <- y[substr(y$chromosome,1,3) == "chr",]
		z$feature_id <- substring(z$feature_id, 8, 1000000L)
		
		### sort utr and cds according to their chromosome and chr_start 
		
		sort.z <- z[order(z$chromosome, z$chr_start),]
		
		### specify file header according to species
		if (species == "mm"){
			header <- "mm9_"
		} else {
			header <- "hg19" 
		}
		utr_cds_file <- paste(header, "utr_cds.txt", sep="")
		exon_file <- paste(header,"exon.txt",sep="")
		gene_file <- paste(header,"gene.txt",sep="")
		ir_file <- paste(header,"integenic_region.txt",sep="")

		if (length(unique(z$feature_id)) > 1000000){
			print("The input file contains too many genes.\n")
		} else {
			write.table(sort.z,file=utr_cds_file, quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
			
			.C("processExons", as.character(utr_cds_file), as.character(exon_file), as.character(gene_file), PACKAGE="Rsubread")
			
			gene <- read.delim(gene_file, stringsAsFactors=FALSE)
			colnames(gene) <- c("entrezid","chromosome", "chr_start","chr_stop")
			sort.gene <- gene[order(gene$chromosome, gene$chr_start),]
			file.remove(gene_file)
			write.table(sort.gene,file=gene_file, quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
			
			exon <- read.delim(exon_file, stringsAsFactors=FALSE)
			colnames(exon) <- c("entrezid","chromosome", "chr_start","chr_stop")
			sort.exon <- exon[order(exon$chromosome, exon$chr_start),]
			file.remove(exon_file)
			write.table(sort.exon,file=exon_file, quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
			
			.C("generate_ir_file", as.character(gene_file), as.character(ir_file), PACKAGE="Rsubread")
			
			file.remove(gene_file)
			file.remove(utr_cds_file)
		}
	
	}

}

