subjunc <- function(index,readfile1,readfile2=NULL,input_format="gzFASTQ",output_format="BAM",output_file=paste(readfile1,"subjunc",output_format,sep="."),phredOffset=33,nsubreads=14,TH1=1,TH2=1,maxMismatches=3,unique=FALSE,nBestLocations=1,indels=5,complexIndels=FALSE,nTrim5=0,nTrim3=0,minFragLength=50,maxFragLength=600,PE_orientation="fr",nthreads=1,readGroupID=NULL,readGroup=NULL,keepReadOrder=FALSE,sortReadsByCoordinates=FALSE,color2base=FALSE,DP_GapOpenPenalty=-1,DP_GapExtPenalty=0,DP_MismatchPenalty=0,DP_MatchScore=2,reportAllJunctions=FALSE,useAnnotation=FALSE,annot.inbuilt="mm10",annot.ext=NULL,isGTF=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",chrAliases=NULL)
{
  if(!.is.64bit.system()) cat("WARNING: your system seems to be 32-bit. Rsubread supports 32-bit systems to a very limited level.\nIt is highly recommended to run Rsubread on a 64-bit system to avoid errors.\n\n")

  extra.params <- .retrieve.tmp.parameters()
  isBCLinput <- FALSE
  isScRNAFastqinput <- FALSE
  if("isBCLinput" %in% names(extra.params))isBCLinput <- extra.params[["isBCLinput"]]
  if("isScRNAFastqinput" %in% names(extra.params))isScRNAFastqinput <- extra.params[["isScRNAFastqinput"]]

  .check_string_param(input_format,"input_format")
  .check_string_param(output_format,"output_format")
  .check_string_param(PE_orientation,"PE_orientation")
  .check_string_param(annot.inbuilt,"annot.inbuilt")
  .check_string_param(GTF.featureType,"GTF.featureType")
  .check_string_param(GTF.attrType,"GTF.attrType")
  .check_string_param(readGroupID,"readGroupID")
  .check_string_param(readGroup,"readGroup")

  out.base.names <- basename(output_file)
  if(any(duplicated(out.base.names))){
    out.table.cols <- output_file
  }else{
    out.table.cols <- out.base.names
  }
  #out.table.cols <-gsub("[[:punct:]]+", ".", out.table.cols)
  #out.table.cols <-gsub(" ", ".", out.table.cols)

  readfile1 <- as.character(readfile1)
  if(!isScRNAFastqinput)readfile1 <- .check_and_NormPath(readfile1, mustWork=TRUE, opt="readfile1")
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

  if((!isScRNAFastqinput) && !all(file.exists(readfile1)))
    stop("One or more input files cannot be found (readfile1).")

  if(!is.null(readfile2)){
    if(!all(file.exists(readfile2)))
      stop("One or more input files cannot be found (readfile2).")
  }

  if(!all(file.exists(dirname(output_file))))
    stop("Invalid path was found in output file name(s).")

  opt <- paste("-i",index,sep=.R_param_splitor)

  if(tolower(input_format) == "sam")
    opt <- paste(opt,"--SAMinput",sep=.R_param_splitor)
  if(tolower(input_format) == "bam")
    opt <- paste(opt,"--BAMinput",sep=.R_param_splitor)

  if(tolower(output_format) == "sam")
    opt <- paste(opt,"--SAMoutput",sep=.R_param_splitor)

  if(isBCLinput)
    opt <- paste(opt,"--BCLinput",sep=.R_param_splitor)
	
  if(isScRNAFastqinput)
    opt <- paste(opt,"--scRNA_FQinput",sep=.R_param_splitor)

  opt <- paste(opt,"-n",nsubreads,"-m",TH1,"-p",TH2,"-M",maxMismatches,"-T",nthreads,"-I",indels,sep=.R_param_splitor)

  if(complexIndels)
    opt <- paste(opt,"--complexIndels",sep=.R_param_splitor)

  if(!unique)
    opt <- paste(opt,"--multiMapping",sep=.R_param_splitor)

  opt <- paste(opt,"-B",nBestLocations,"-d",minFragLength,"-D",maxFragLength,"-S",PE_orientation,"--trim5",nTrim5,"--trim3",nTrim3,sep=.R_param_splitor)

  if(!is.null(readGroupID))
    opt <- paste(opt,"--rg-id",readGroupID,sep=.R_param_splitor)
  if(!is.null(readGroup))
    opt <- paste(opt,"--rg",readGroup,sep=.R_param_splitor)
  if(color2base)
    opt <- paste(opt,"-b",sep=.R_param_splitor)
  if(keepReadOrder)
    opt <- paste(opt,"--keepReadOrder",sep=.R_param_splitor)
  if(sortReadsByCoordinates)
    opt <- paste(opt,"--sortReadsByCoordinates",sep=.R_param_splitor)

  opt <- paste(opt,"-G",DP_GapOpenPenalty,"-E",DP_GapExtPenalty,"-X",DP_MismatchPenalty,"-Y",DP_MatchScore,sep=.R_param_splitor)

  if(reportAllJunctions)
    opt <- paste(opt,"--allJunctions",sep=.R_param_splitor)

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
      mm39={
        ann <- system.file("annot","mm39_RefSeq_exon.txt",package="Rsubread")
        annot.display <- "inbuilt (mm39)"
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
        fout_annot <- file.path(".",paste(".Rsubread_subjunc_UserProvidedAnnotation_pid",Sys.getpid(),sep=""))
        oldScipen <- options(scipen=999)
        write.table(x=annot_df,file=fout_annot,sep="\t",row.names=FALSE,quote=FALSE)
        options(oldScipen)
        ann <- fout_annot
        flag <- TRUE
      }
    }

    opt <- paste(opt,"-a",ann,sep=.R_param_splitor)

    if(isGTF)
      opt <- paste(opt,"-F","GTF","--gtfFeature",GTF.featureType,"--gtfAttr",GTF.attrType,sep=.R_param_splitor)
    else
      opt <- paste(opt,"-F","SAF",sep=.R_param_splitor)

    if(!is.null(chrAliases))
      opt <- paste(opt,"-A",chrAliases,sep=.R_param_splitor)

    opt <- paste(opt, "--exonAnnotationScreenOut", annot.display, sep=.R_param_splitor)
  }

  if(phredOffset == 33)
    opt <- paste(opt,"-P",3,sep=.R_param_splitor)
  else
    opt <- paste(opt,"-P",6,sep=.R_param_splitor)

  for(i in 1:length(readfile1)){
    opt_files <- paste("-r",readfile1[i],sep=.R_param_splitor)
    if(!is.null(readfile2)) 
      opt_files <- paste(opt_files,"-R",readfile2[i],sep=.R_param_splitor)
    opt_files <- paste(opt_files,"-o",output_file[i],sep=.R_param_splitor)

    cmd <- paste("subjunc",opt_files,opt,sep=.R_param_splitor)
    n <- length(unlist(strsplit(cmd,.R_param_splitor)))
    C_args <- .C("R_junction_wrapper",as.integer(n),as.character(cmd),PACKAGE="Rsubread")
    summary.data <- .load.delete.summary(output_file[i])
    if(i ==1){
      return.summary <- summary.data
    }else{
      return.summary <- data.frame(return.summary, summary.data)
    }
  }

  if(flag)
    file.remove(fout_annot)

  names(return.summary) <- out.table.cols # basename(output_file)
  return.summary
}
