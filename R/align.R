align <- function(index,readfile1,readfile2=NULL,input_format="gzFASTQ",output_format="BAM",output_file=paste(readfile1,"subread",output_format,sep="."),nsubreads=10,TH1=3,TH2=1,maxMismatches=3,nthreads=1,indels=5,phredOffset=33,tieBreakQS=FALSE,tieBreakHamming=TRUE,unique=TRUE,nBestLocations=1,minFragLength=50,maxFragLength=600,PE_orientation="fr",nTrim5=0,nTrim3=0,readGroupID=NULL,readGroup=NULL,color2base=FALSE,DP_GapOpenPenalty=-1,DP_GapExtPenalty=0,DP_MismatchPenalty=0,DP_MatchScore=2,reportFusions=FALSE)
{
	if(length(readfile1) != length(output_file))
		stop("The number of input file names is different from the number of output file names.")

	if(!is.null(readfile2)){
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

	if(tolower(input_format) == "sam")
		opt <- paste(opt,"--SAMinput",sep=",")
	if(tolower(input_format) == "bam")
		opt <- paste(opt,"--BAMinput",sep=",")
	if(tolower(input_format) == "gzfastq")
		opt <- paste(opt,"--gzFASTQinput",sep=",")

	if(tolower(output_format) == "bam")
		opt <- paste(opt,"--BAMoutput",sep=",")	  

	opt <- paste(opt,"-n",nsubreads,"-m",TH1,"-p",TH2,"-M",maxMismatches,"-T",nthreads,"-I",indels,sep=",")

	if(tieBreakQS)
		opt <- paste(opt,"-Q",sep=",")
	if(tieBreakHamming)
		opt <- paste(opt,"-H",sep=",")
	if(unique)
		opt <- paste(opt,"-u",sep=",")

	opt <- paste(opt,"-B",nBestLocations,"-d",minFragLength,"-D",maxFragLength,"-S",PE_orientation,"--trim5",nTrim5,"--trim3",nTrim3,sep=",")

	if(!is.null(readGroupID))
	  opt <- paste(opt,"--rg-id",readGroupID,sep=",")
	if(!is.null(readGroup))
	  opt <- paste(opt,"--rg",readGroup,sep=",")
	if(color2base)
		opt <- paste(opt,"-b",sep=",")

	opt <- paste(opt,"-G",DP_GapOpenPenalty,"-E",DP_GapExtPenalty,"-X",DP_MismatchPenalty,"-Y",DP_MatchScore,sep=",")

	if(reportFusions)
		opt <- paste(opt,"--reportFusions",sep=",")
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
	}
}
