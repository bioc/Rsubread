.load.delete.summary <- function(bam.name){
  sumfile <- paste0(bam.name,".summary")
  if(!file.exists(sumfile)){
    stop(paste("ERROR: Summary file",sumfile,"was not generated! The program terminated wrongly!"))
  } 

  tmp.frame <- read.delim(sumfile, header=FALSE, row.names=1, stringsAsFactors=FALSE)
#  file.remove(sumfile)
  return(tmp.frame)
}

align <- function(index,readfile1,readfile2=NULL,type="rna",input_format="gzFASTQ",output_format="BAM",output_file=paste(readfile1,"subread",output_format,sep="."),phredOffset=33,nsubreads=10,TH1=3,TH2=1,maxMismatches=3,unique=FALSE,nBestLocations=1,indels=5,complexIndels=FALSE,nTrim5=0,nTrim3=0,minFragLength=50,maxFragLength=600,PE_orientation="fr",nthreads=1,readGroupID=NULL,readGroup=NULL,keepReadOrder=FALSE,sortReadsByCoordinates=FALSE,color2base=FALSE,DP_GapOpenPenalty=-1,DP_GapExtPenalty=0,DP_MismatchPenalty=0,DP_MatchScore=2,detectSV=FALSE,useAnnotation=FALSE,annot.inbuilt="mm10",annot.ext=NULL,isGTF=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",chrAliases=NULL)
{
  readfile1 <- as.character(readfile1)
  readfile1 <- .check_and_NormPath(readfile1, mustWork=TRUE, opt="readfile1")
  output_file <- .check_and_NormPath(output_file, mustWork=FALSE, opt="output_file")
  index <- .check_and_NormPath(index, mustWork=FALSE, opt="index")
  if((!is.null(annot.ext)) && is.character(annot.ext)) annot.ext = .check_and_NormPath(annot.ext, mustWork=TRUE, opt="annot.ext")
  if(!is.null(chrAliases)) chrAliases = .check_and_NormPath(chrAliases, mustWork=TRUE, opt="chrAliases")

  if(length(readfile1) != length(output_file))
    stop("The number of input file names is different from the number of output file names.")

  if(!is.null(readfile2)){
    readfile2 <- as.character(readfile2)
    readfile2 <- .check_and_NormPath(readfile2, mustWork=TRUE, "readfile2")
    if(length(readfile1) != length(readfile2))
      stop("The number of file names for the first reads is different from the number of file names for the second reads.")
  }

  if(!all(file.exists(readfile1)))
    stop("One or more input files cannot be found (readfile1).")

  if(!is.null(readfile2)){
    if(!all(file.exists(readfile2)))
      stop("One or more input files cannot be found (readfile2).")
  }

  if(!all(file.exists(dirname(output_file))))
    stop("Invalid path was found in output file name(s).")

  opt <- paste("-i",index,sep=",")
  
  type_C <- -1
  if(is.character(type)){
    if(tolower(type) == "rna") type_C <- 0
    if(tolower(type) == "dna") type_C <- 1
  }
  if(is.numeric(type)){
    if (type == 0 || type == 1) type_C <- type
  }
  if(type_C == -1) stop("Invalid value was provided for data type.")
  
  opt <- paste(opt,"--type",type_C,sep=",")

  if(tolower(input_format) == "sam")
    opt <- paste(opt,"--SAMinput",sep=",")
  if(tolower(input_format) == "bam")
    opt <- paste(opt,"--BAMinput",sep=",")

  if(tolower(output_format) == "sam")
    opt <- paste(opt,"--SAMoutput",sep=",")

  opt <- paste(opt,"-n",nsubreads,"-m",TH1,"-p",TH2,"-M",maxMismatches,"-T",nthreads,"-I",indels,sep=",")

  if(complexIndels)
    opt <- paste(opt,"--complexIndels",sep=",")

  if(!unique)
    opt <- paste(opt,"--multiMapping",sep=",")

  opt <- paste(opt,"-B",nBestLocations,"-d",minFragLength,"-D",maxFragLength,"-S",PE_orientation,"--trim5",nTrim5,"--trim3",nTrim3,sep=",")

  if(!is.null(readGroupID))
    opt <- paste(opt,"--rg-id",readGroupID,sep=",")
  if(!is.null(readGroup))
    opt <- paste(opt,"--rg",readGroup,sep=",")
  if(color2base)
    opt <- paste(opt,"-b",sep=",")
  if(keepReadOrder)
    opt <- paste(opt,"--keepReadOrder",sep=",")
  if(sortReadsByCoordinates)
    opt <- paste(opt,"--sortReadsByCoordinates",sep=",")

  opt <- paste(opt,"-G",DP_GapOpenPenalty,"-E",DP_GapExtPenalty,"-X",DP_MismatchPenalty,"-Y",DP_MatchScore,sep=",")

  if(detectSV)
    opt <- paste(opt,"--sv",sep=",")

  flag <- FALSE
  if(useAnnotation){
    annot.display <- "R data.frame"
    if(is.null(annot.ext)){
      switch(tolower(as.character(annot.inbuilt)),
      mm9={
        ann <- system.file("annot","mm9_RefSeq_exon.txt",package="Rsubread")
        annot.display <- "inbuilt (mm9)"
      },
      mm10={
        ann <- system.file("annot","mm10_RefSeq_exon.txt",package="Rsubread")
        annot.display <- "inbuilt (mm10)"
      },
      hg19={
        ann <- system.file("annot","hg19_RefSeq_exon.txt",package="Rsubread")
        annot.display <- "inbuilt (hg19)"
      },
      hg38={
        ann <- system.file("annot","hg38_RefSeq_exon.txt",package="Rsubread")
        annot.display <- "inbuilt (hg38)"
      },
      {
        stop("In-built annotation for ", annot.inbuilt, " is not available.\n")
      }
      ) # end switch
    }
    else{
      if(is.character(annot.ext)){
        ann <- annot.ext
        annot.display <- paste0(basename(ann), " (", ifelse(isGTF,"GTF","SAF"),")")
      }
      else{
        annot_df <- as.data.frame(annot.ext,stringsAsFactors=FALSE)
        if(sum(c("geneid","chr","start","end", "strand") %in% tolower(colnames(annot_df))) != 5)
          stop("One or more required columns are missing in the provided annotation data. Please refer to help page for annotation format.\n")
        colnames(annot_df) <- tolower(colnames(annot_df))
        annot_df <- data.frame(geneid=annot_df$geneid,chr=annot_df$chr,start=annot_df$start,end=annot_df$end,strand=annot_df$strand,stringsAsFactors=FALSE)
        annot_df$chr <- as.character(annot_df$chr)
        fout_annot <- file.path(".",paste(".Rsubread_align_UserProvidedAnnotation_pid",Sys.getpid(),sep=""))
        oldScipen <- options(scipen=999)
        write.table(x=annot_df,file=fout_annot,sep="\t",row.names=FALSE,quote=FALSE)
        options(oldScipen)
        ann <- fout_annot
        flag <- TRUE
      }
    }

    opt <- paste(opt,"-a",ann,sep=",")

    if(isGTF)
      opt <- paste(opt,"-F","GTF","--gtfFeature",GTF.featureType,"--gtfAttr",GTF.attrType,sep=",")
    else
      opt <- paste(opt,"-F","SAF",sep=",")

    if(!is.null(chrAliases))
      opt <- paste(opt,"-A",chrAliases,sep=",")

    opt <- paste(opt, "--exonAnnotationScreenOut", annot.display, sep=",")
  }

  if(phredOffset == 33)
    opt <- paste(opt,"-P",3,sep=",")
  else
    opt <- paste(opt,"-P",6,sep=",")

  for(i in 1:length(readfile1)){
    opt_files <- paste("-r",readfile1[i],sep=",")
    if(!is.null(readfile2)) 
      opt_files <- paste(opt_files,"-R",readfile2[i],sep=",")
    opt_files <- paste(opt_files,"-o",output_file[i],sep=",")

    cmd <- paste("subread-align",opt_files,opt,sep=",")
    n <- length(unlist(strsplit(cmd,",")))
    C_args <- .C("R_align_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
    summary.data <- .load.delete.summary(output_file[i])
    if(i == 1){
      return.summary <- summary.data
    }else{
      return.summary <- data.frame(return.summary, summary.data)
    }
  }

  if(flag)
    file.remove(fout_annot)

  names(return.summary) <- basename(output_file)
  return.summary
}
