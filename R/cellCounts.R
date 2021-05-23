.write.tmp.parameters <- function(param){
  temp.file.prefix <- file.path(".",paste(".Rsubread_PARAM_cellCounts_",Sys.getpid(),sep=""))
  param$PMR_CREATE_TIME <- Sys.time()
  saveRDS(param, temp.file.prefix)
}

.retrieve.tmp.parameters <- function(){
  temp.file.prefix <- file.path(".",paste(".Rsubread_PARAM_cellCounts_",Sys.getpid(),sep=""))
  if(!file.exists(temp.file.prefix)) return (NA)
  rv <- readRDS(temp.file.prefix)
  file.remove(temp.file.prefix)
  if( Sys.time() - rv$PMR_CREATE_TIME > 5) return (NA) # the file should be created no longer than 5 seconds ago
  rv
}

.index.names.to.sheet.BAM.mode<-function(fqtab, fname){
  lines <- paste(fqtab$BAMFile, fqtab$SampleName , sep=",")
  fileConn<-file(fname)
  writeLines(c("EMFileVersion,4","[data]","BAMFile,Sample_Name",lines), fileConn)
  close(fileConn)
}

.index.names.to.sheet.FASTQ.mode<-function(fqtab, fname){
  lines <- paste(fqtab$BarcodeUMIFile, fqtab$ReadFile, fqtab$SampleName , sep=",")
  fileConn<-file(fname)
  writeLines(c("EMFileVersion,4","[data]","BarcodeUMIFile,ReadFile,Sample_Name",lines), fileConn)
  close(fileConn)
}

.index.names.to.sheet.raw.dir.mode<-function(dirname, nametab, fname, sample.name=NA){
  lanes <- c()
  index.names <- c()
  sample.names <- c()
  index.seq <- c()
  for(cli in 1:nrow(nametab)){
    if(nametab$InputDirectory[cli]!=dirname)next
    if((!is.na(sample.name)) && as.character(nametab$SampleName[cli])!=sample.name)next
    seqs <- .convert.sample_index.id.to.seq(nametab$IndexSetName[cli])
    for(seq in seqs){
      lanes <-c(lanes, nametab$Lane[cli])
      index.names <-c(index.names, nametab$IndexSetName[cli])
      sample.names <- c(sample.names, as.character(nametab$SampleName[cli]))
      index.seq <- c(index.seq, seq)
    }
  }
  sdf<-data.frame(Lane=lanes, Sample_ID=index.names, Sample_Name=sample.names, index=index.seq, Sample_Project=rep("cellCounts", length(lanes)), stringsAsFactors=F)
  fileConn<-file(fname)
  lines <- paste( sdf$Lane, sdf$Sample_ID, sdf$Sample_Name, sdf$index, sdf$Sample_Project ,sep=",")
  writeLines(c("EMFileVersion,4","[data]","Lane,Sample_ID,Sample_Name,index,Sample_Project",lines), fileConn)
  close(fileConn)
}

.convert.sample_index.id.to.seq <- function(spid){
  seqs <- .one.convert.sample_index.id.to.seq(spid)
  if(any(is.na(seqs))){
    spid <- gsub("-","_", spid)
    seqs <- .one.convert.sample_index.id.to.seq(spid)
  }
  if(any(is.na(seqs)))stop(paste("Unable to find sample index for", spid))
  seqs
}
.one.convert.sample_index.id.to.seq <- function(spid){
  spseqs <- NA
  if(spid == "SI_001") spseqs <- c('TCGCCATA', 'GTATACAC', 'AATGGTGG', 'CGCATGCT')
  if(spid == "SI_002") spseqs <- c('TATCCTCG', 'GCGAGGTC', 'CGCTTCAA', 'ATAGAAGT')
  if(spid == "SI_003") spseqs <- c('TGACGTCG', 'CTTGTGTA', 'ACGACCGT', 'GACTAAAC')
  if(spid == "SI_004") spseqs <- c('ATCTAGCT', 'GAGCGTAC', 'TCAGCCTG', 'CGTATAGA')
  if(spid == "SI_005") spseqs <- c('CCGTTCCC', 'ATACAGTT', 'TGTAGTAA', 'GACGCAGG')
  if(spid == "SI_006") spseqs <- c('TCAATTGG', 'AGTTAGAA', 'GAGCGCTT', 'CTCGCACC')
  if(spid == "SI_007") spseqs <- c('CTGCCTTG', 'ACTAGCCC', 'GGCTAGAT', 'TAAGTAGA')
  if(spid == "SI_008") spseqs <- c('GGCAGAAA', 'ACGGTTCT', 'CATTCGTC', 'TTACACGG')
  if(spid == "SI_P01_A1") spseqs <- c('TTGTAAGA', 'GGCGTTTC', 'CCTACCAT', 'AAACGGCG')
  if(spid == "SI_P01_B1") spseqs <- c('CTAGCTGT', 'GCCAACAA', 'AGGCTACC', 'TATTGGTG')
  if(spid == "SI_P01_C1") spseqs <- c('GATGCAGT', 'AGACTTTC', 'TTCTACAG', 'CCGAGGCA')
  if(spid == "SI_P01_D1") spseqs <- c('AGCTGCGT', 'GTGGAGCA', 'TCTATTAG', 'CAACCATC')
  if(spid == "SI_P01_E1") spseqs <- c('CGCCCGTA', 'GTTTGCCT', 'TAAGTTAG', 'ACGAAAGC')
  if(spid == "SI_P01_F1") spseqs <- c('TGACTAGT', 'GATAGGTA', 'CCGTACAG', 'ATCGCTCC')
  if(spid == "SI_P01_G1") spseqs <- c('CATATGCG', 'ATGCGATT', 'TCCGCTAC', 'GGATACGA')
  if(spid == "SI_P01_H1") spseqs <- c('TCTTGTCC', 'CGGGAGTA', 'GTCACAGG', 'AAACTCAT')
  if(spid == "SI_P01_A2") spseqs <- c('AGCCCTTT', 'TCTTAGGC', 'GTGAGAAG', 'CAAGTCCA')
  if(spid == "SI_P01_B2") spseqs <- c('GGTCGAGC', 'TTAGATTG', 'CCCACCCA', 'AAGTTGAT')
  if(spid == "SI_P01_C2") spseqs <- c('CCGAGAAC', 'TGCTCTGT', 'GTAGTGCG', 'AATCACTA')
  if(spid == "SI_P01_D2") spseqs <- c('ACATTCCG', 'GACACAAT', 'CTGCGGTA', 'TGTGATGC')
  if(spid == "SI_P01_E2") spseqs <- c('TCATCAAG', 'GTTGGTCC', 'AGGCTGGT', 'CACAACTA')
  if(spid == "SI_P01_F2") spseqs <- c('TGCCCGCT', 'GCAAACGC', 'CATTGATA', 'ATGGTTAG')
  if(spid == "SI_P01_G2") spseqs <- c('CGGTGAGC', 'ATAACCTA', 'TCCGAGCG', 'GATCTTAT')
  if(spid == "SI_P01_H2") spseqs <- c('CCGAACTC', 'AACGGTCA', 'TTATTGGT', 'GGTCCAAG')
  if(spid == "SI_P01_A3") spseqs <- c('AAAGCATA', 'GCCTTTAT', 'CTGCAGCC', 'TGTAGCGG')
  if(spid == "SI_P01_B3") spseqs <- c('TCATCCTT', 'ATTGGACG', 'CAGCTTAC', 'GGCAAGGA')
  if(spid == "SI_P01_C3") spseqs <- c('ACGTTACA', 'TTACCTAC', 'GACGACGG', 'CGTAGGTT')
  if(spid == "SI_P01_D3") spseqs <- c('GAGCACGC', 'CGAAGTTG', 'TTCGTGAA', 'ACTTCACT')
  if(spid == "SI_P01_E3") spseqs <- c('TCTGCAGG', 'CGGCTCCA', 'AACAAGTC', 'GTATGTAT')
  if(spid == "SI_P01_F3") spseqs <- c('TTAGGACC', 'AGTCTGTA', 'GCCTCCGT', 'CAGAATAG')
  if(spid == "SI_P01_G3") spseqs <- c('TACGAGTT', 'ATGTCCAG', 'GCTATAGC', 'CGACGTCA')
  if(spid == "SI_P01_H3") spseqs <- c('TTGGGCTT', 'GACAAACC', 'ACACCTAA', 'CGTTTGGG')
  if(spid == "SI_P01_A4") spseqs <- c('CATGGCAG', 'AGAACGCC', 'GTCTTTGA', 'TCGCAATT')
  if(spid == "SI_P01_B4") spseqs <- c('GACAGGCT', 'CCTCTAAC', 'AGGGACTG', 'TTATCTGA')
  if(spid == "SI_P01_C4") spseqs <- c('ACATTGGC', 'GAGCCCAT', 'CTTAGTCA', 'TGCGAATG')
  if(spid == "SI_P01_D4") spseqs <- c('AAATCGTC', 'GCGATCGG', 'CTTCGAAT', 'TGCGATCA')
  if(spid == "SI_P01_E4") spseqs <- c('GTAAACAT', 'TATCCTGA', 'AGCTGACG', 'CCGGTGTC')
  if(spid == "SI_P01_F4") spseqs <- c('GCATGATA', 'CGTCCTCT', 'AACGACAC', 'TTGATGGG')
  if(spid == "SI_P01_G4") spseqs <- c('CCTGCGGT', 'GTACAACG', 'AGCTTCTC', 'TAGAGTAA')
  if(spid == "SI_P01_H4") spseqs <- c('TTCCATCT', 'ACTGGAGC', 'CGGTCGTG', 'GAAATCAA')
  if(spid == "SI_P01_A5") spseqs <- c('CAGTCTGG', 'TCACACTC', 'ATTGGGAA', 'GGCATACT')
  if(spid == "SI_P01_B5") spseqs <- c('TCATGCGA', 'ATCGTACT', 'CATCAGTG', 'GGGACTAC')
  if(spid == "SI_P01_C5") spseqs <- c('TCAGTCAA', 'CACTGACT', 'ATGCATTC', 'GGTACGGG')
  if(spid == "SI_P01_D5") spseqs <- c('GTCGACTC', 'AGACGGAT', 'CCTTTAGA', 'TAGACTCG')
  if(spid == "SI_P01_E5") spseqs <- c('CCGTTGAA', 'TATGCTCT', 'ATCCAAGG', 'GGAAGCTC')
  if(spid == "SI_P01_F5") spseqs <- c('TCTGACTA', 'GTACGGCT', 'CGGTTTAG', 'AACACAGC')
  if(spid == "SI_P01_G5") spseqs <- c('ATGAAGTA', 'GAAGCTCG', 'TCTTTCGT', 'CGCCGAAC')
  if(spid == "SI_P01_H5") spseqs <- c('ATAGTATG', 'TATAAGGA', 'GGCTCCAC', 'CCGCGTCT')
  if(spid == "SI_P01_A6") spseqs <- c('CTTTCGAC', 'ACGGGACT', 'TGCATCTG', 'GAACATGA')
  if(spid == "SI_P01_B6") spseqs <- c('GCGCACCT', 'AACGCGAA', 'CTATTTGG', 'TGTAGATC')
  if(spid == "SI_P01_C6") spseqs <- c('CGCTCAGG', 'GAGGTTTA', 'ACTCAGAC', 'TTAAGCCT')
  if(spid == "SI_P01_D6") spseqs <- c('GAAGTCTT', 'TGCAGGGC', 'ATGCCAAA', 'CCTTATCG')
  if(spid == "SI_P01_E6") spseqs <- c('TGATGGCT', 'GCCATTTG', 'ATTGAAAC', 'CAGCCCGA')
  if(spid == "SI_P01_F6") spseqs <- c('GACTTCCT', 'TGAGGAAG', 'ATGCCGGC', 'CCTAATTA')
  if(spid == "SI_P01_G6") spseqs <- c('GTGCGACA', 'TCAGTGTT', 'AGCACTGG', 'CATTACAC')
  if(spid == "SI_P01_H6") spseqs <- c('AAGCATAA', 'CCCATCGC', 'TTAGCGCT', 'GGTTGATG')
  if(spid == "SI_P01_A7") spseqs <- c('CTCATCAT', 'TAACGTCC', 'AGGTCATA', 'GCTGAGGG')
  if(spid == "SI_P01_B7") spseqs <- c('TCCACACG', 'CTTCTGTT', 'GAATGCAC', 'AGGGATGA')
  if(spid == "SI_P01_C7") spseqs <- c('GGCGGAAT', 'ACACCGGG', 'CATAATCC', 'TTGTTCTA')
  if(spid == "SI_P01_D7") spseqs <- c('CCGGATCC', 'GGTCGCAT', 'TTAACGTG', 'AACTTAGA')
  if(spid == "SI_P01_E7") spseqs <- c('AAGACGTG', 'CCATGTGT', 'GTTCACAA', 'TGCGTACC')
  if(spid == "SI_P01_F7") spseqs <- c('GGTTAGAC', 'CAAACTTT', 'ACCCGAGA', 'TTGGTCCG')
  if(spid == "SI_P01_G7") spseqs <- c('GCCGGTAA', 'TGACTGCC', 'ATTACCGG', 'CAGTAATT')
  if(spid == "SI_P01_H7") spseqs <- c('TGGCACGA', 'AACGGGTG', 'CTAATTCT', 'GCTTCAAC')
  if(spid == "SI_P01_A8") spseqs <- c('GACTGTTC', 'ATGATACG', 'CCACAGAA', 'TGTGCCGT')
  if(spid == "SI_P01_B8") spseqs <- c('ACGTTCAC', 'TGCCCAGA', 'CAAGGTCT', 'GTTAAGTG')
  if(spid == "SI_P01_C8") spseqs <- c('TTTATCCC', 'GCACGTTT', 'CAGGAAGA', 'AGCTCGAG')
  if(spid == "SI_P01_D8") spseqs <- c('AATCTTTG', 'GGATGAGT', 'CTCAAGAC', 'TCGGCCCA')
  if(spid == "SI_P01_E8") spseqs <- c('GCTTACAT', 'TAGGGTGC', 'AGCCTATG', 'CTAACGCA')
  if(spid == "SI_P01_F8") spseqs <- c('AGTTGGGA', 'TACATTCT', 'CCAGAAAG', 'GTGCCCTC')
  if(spid == "SI_P01_G8") spseqs <- c('AAGTACTC', 'GGAACTCT', 'TCCCTGAG', 'CTTGGAGA')
  if(spid == "SI_P01_H8") spseqs <- c('AAGAGCGG', 'TCATAGCA', 'GGCCCATC', 'CTTGTTAT')
  if(spid == "SI_P01_A9") spseqs <- c('GAGTGCGT', 'CTCCAACA', 'ACAACTTG', 'TGTGTGAC')
  if(spid == "SI_P01_B9") spseqs <- c('AAGCGTGT', 'CTTGACCG', 'TGATTAAC', 'GCCACGTA')
  if(spid == "SI_P01_C9") spseqs <- c('AGATCGGT', 'CATCGTCG', 'GTCATATA', 'TCGGACAC')
  if(spid == "SI_P01_D9") spseqs <- c('CAAGGGAC', 'ACCTACTG', 'GGGACACA', 'TTTCTTGT')
  if(spid == "SI_P01_E9") spseqs <- c('AGTAAGCA', 'TACCGCGG', 'CCGGTAAT', 'GTATCTTC')
  if(spid == "SI_P01_F9") spseqs <- c('AGTTAGTT', 'GTACTTAA', 'CACGCACG', 'TCGAGCGC')
  if(spid == "SI_P01_G9") spseqs <- c('TTGACTTC', 'GCCGAAGT', 'CAATGGCA', 'AGTCTCAG')
  if(spid == "SI_P01_H9") spseqs <- c('GGAATATG', 'ACCTGCCA', 'CTTCATAC', 'TAGGCGGT')
  if(spid == "SI_P01_A10") spseqs <- c('ACAGCAAC', 'TTTCGCGA', 'CGCAATTT', 'GAGTTGCG')
  if(spid == "SI_P01_B10") spseqs <- c('ACCATTAA', 'CTGGACGT', 'GAACGGTC', 'TGTTCACG')
  if(spid == "SI_P01_C10") spseqs <- c('CGTGCTAA', 'TCACTCCT', 'ATCTGATC', 'GAGAAGGG')
  if(spid == "SI_P01_D10") spseqs <- c('CTTATTTG', 'GCGGGCAT', 'AGATAACA', 'TACCCGGC')
  if(spid == "SI_P01_E10") spseqs <- c('GCACCAGT', 'CGCAGGAG', 'TAGTACCA', 'ATTGTTTC')
  if(spid == "SI_P01_F10") spseqs <- c('TCGACAAT', 'GAATACTG', 'ATTCGTGC', 'CGCGTGCA')
  if(spid == "SI_P01_G10") spseqs <- c('CGGAGACT', 'TCCTATGA', 'ATACTGAG', 'GATGCCTC')
  if(spid == "SI_P01_H10") spseqs <- c('GACCGCCA', 'TCGAATTG', 'ATTTCAGC', 'CGAGTGAT')
  if(spid == "SI_P01_A11") spseqs <- c('CTTTCCTT', 'TAGGTAAA', 'ACCAGTCC', 'GGACAGGG')
  if(spid == "SI_P01_B11") spseqs <- c('TCCAGATA', 'GATTCGCT', 'CGACATAG', 'ATGGTCGC')
  if(spid == "SI_P01_C11") spseqs <- c('GTTTGTGG', 'ACCGAACA', 'TAGACGAC', 'CGACTCTT')
  if(spid == "SI_P01_D11") spseqs <- c('GCTACTTC', 'CACCTCAG', 'ATATGAGA', 'TGGGAGCT')
  if(spid == "SI_P01_E11") spseqs <- c('ATCGCCAT', 'TCACGGTA', 'GGGTTTCC', 'CATAAAGG')
  if(spid == "SI_P01_F11") spseqs <- c('GAACCCGG', 'AGCAGTTA', 'TCGTAGAT', 'CTTGTACC')
  if(spid == "SI_P01_G11") spseqs <- c('AGGGCGTT', 'CTATACGC', 'TACATAAG', 'GCTCGTCA')
  if(spid == "SI_P01_H11") spseqs <- c('TCTCGACT', 'AGGATCGA', 'CACGATTC', 'GTATCGAG')
  if(spid == "SI_P01_A12") spseqs <- c('TTATGGAA', 'ACTACTGT', 'CGGGAACG', 'GACCTCTC')
  if(spid == "SI_P01_B12") spseqs <- c('GAAAGACA', 'CGCTACAT', 'ACGCTTGG', 'TTTGCGTC')
  if(spid == "SI_P01_C12") spseqs <- c('TAAGCCAC', 'CCGTTATG', 'GGTAATGT', 'ATCCGGCA')
  if(spid == "SI_P01_D12") spseqs <- c('GCTGTGTA', 'AGAAACGT', 'CACTCAAC', 'TTGCGTCG')
  if(spid == "SI_P01_E12") spseqs <- c('CGCTATCC', 'ACGCGGAA', 'TAAATCGT', 'GTTGCATG')
  if(spid == "SI_P01_F12") spseqs <- c('AATTGAAC', 'TGGACCCT', 'CCAGTGGA', 'GTCCATTG')
  if(spid == "SI_P01_G12") spseqs <- c('CATGCGTA', 'ACCCGCAC', 'TGATATCG', 'GTGATAGT')
  if(spid == "SI_P01_H12") spseqs <- c('TGTGTATA', 'GTCCAGGC', 'CAATCCCT', 'ACGAGTAG')
  
  if(spid == "SI_3A_A1" || spid == "SI_P03_A1") spseqs <- c('AAACGGCG', 'CCTACCAT', 'GGCGTTTC', 'TTGTAAGA')
  if(spid == "SI_3A_B1" || spid == "SI_P03_B1") spseqs <- c('AGGCTACC', 'CTAGCTGT', 'GCCAACAA', 'TATTGGTG')
  if(spid == "SI_3A_C1" || spid == "SI_P03_C1") spseqs <- c('AGACTTTC', 'CCGAGGCA', 'GATGCAGT', 'TTCTACAG')
  if(spid == "SI_3A_D1" || spid == "SI_P03_D1") spseqs <- c('AGCTGCGT', 'CAACCATC', 'GTGGAGCA', 'TCTATTAG')
  if(spid == "SI_3A_E1" || spid == "SI_P03_E1") spseqs <- c('ACGAAAGC', 'CGCCCGTA', 'GTTTGCCT', 'TAAGTTAG')
  if(spid == "SI_3A_F1" || spid == "SI_P03_F1") spseqs <- c('ATCGCTCC', 'CCGTACAG', 'GATAGGTA', 'TGACTAGT')
  if(spid == "SI_3A_G1" || spid == "SI_P03_G1") spseqs <- c('ATGCGATT', 'CATATGCG', 'GGATACGA', 'TCCGCTAC')
  if(spid == "SI_3A_H1" || spid == "SI_P03_H1") spseqs <- c('AAACTCAT', 'CGGGAGTA', 'GTCACAGG', 'TCTTGTCC')
  if(spid == "SI_3A_A2" || spid == "SI_P03_A2") spseqs <- c('AGCCCTTT', 'CAAGTCCA', 'GTGAGAAG', 'TCTTAGGC')
  if(spid == "SI_3A_B2" || spid == "SI_P03_B2") spseqs <- c('AAGTTGAT', 'CCCACCCA', 'GGTCGAGC', 'TTAGATTG')
  if(spid == "SI_3A_C2" || spid == "SI_P03_C2") spseqs <- c('AATCACTA', 'CCGAGAAC', 'GTAGTGCG', 'TGCTCTGT')
  if(spid == "SI_3A_D2" || spid == "SI_P03_D2") spseqs <- c('ACATTCCG', 'CTGCGGTA', 'GACACAAT', 'TGTGATGC')
  if(spid == "SI_3A_E2" || spid == "SI_P03_E2") spseqs <- c('AGGCTGGT', 'CACAACTA', 'GTTGGTCC', 'TCATCAAG')
  if(spid == "SI_3A_F2" || spid == "SI_P03_F2") spseqs <- c('ATGGTTAG', 'CATTGATA', 'GCAAACGC', 'TGCCCGCT')
  if(spid == "SI_3A_G2" || spid == "SI_P03_G2") spseqs <- c('ATAACCTA', 'CGGTGAGC', 'GATCTTAT', 'TCCGAGCG')
  if(spid == "SI_3A_H2" || spid == "SI_P03_H2") spseqs <- c('AACGGTCA', 'CCGAACTC', 'GGTCCAAG', 'TTATTGGT')
  if(spid == "SI_3A_A3" || spid == "SI_P03_A3") spseqs <- c('AAAGCATA', 'CTGCAGCC', 'GCCTTTAT', 'TGTAGCGG')
  if(spid == "SI_3A_B3" || spid == "SI_P03_B3") spseqs <- c('ATTGGACG', 'CAGCTTAC', 'GGCAAGGA', 'TCATCCTT')
  if(spid == "SI_3A_C3" || spid == "SI_P03_C3") spseqs <- c('ACGTTACA', 'CGTAGGTT', 'GACGACGG', 'TTACCTAC')
  if(spid == "SI_3A_D3" || spid == "SI_P03_D3") spseqs <- c('ACTTCACT', 'CGAAGTTG', 'GAGCACGC', 'TTCGTGAA')
  if(spid == "SI_3A_E3" || spid == "SI_P03_E3") spseqs <- c('AACAAGTC', 'CGGCTCCA', 'GTATGTAT', 'TCTGCAGG')
  if(spid == "SI_3A_F3" || spid == "SI_P03_F3") spseqs <- c('AGTCTGTA', 'CAGAATAG', 'GCCTCCGT', 'TTAGGACC')
  if(spid == "SI_3A_G3" || spid == "SI_P03_G3") spseqs <- c('ATGTCCAG', 'CGACGTCA', 'GCTATAGC', 'TACGAGTT')
  if(spid == "SI_3A_H3" || spid == "SI_P03_H3") spseqs <- c('ACACCTAA', 'CGTTTGGG', 'GACAAACC', 'TTGGGCTT')
  if(spid == "SI_3A_A4" || spid == "SI_P03_A4") spseqs <- c('AGAACGCC', 'CATGGCAG', 'GTCTTTGA', 'TCGCAATT')
  if(spid == "SI_3A_B4" || spid == "SI_P03_B4") spseqs <- c('AGGGACTG', 'CCTCTAAC', 'GACAGGCT', 'TTATCTGA')
  if(spid == "SI_3A_C4" || spid == "SI_P03_C4") spseqs <- c('ACATTGGC', 'CTTAGTCA', 'GAGCCCAT', 'TGCGAATG')
  if(spid == "SI_3A_D4" || spid == "SI_P03_D4") spseqs <- c('AAATCGTC', 'CTTCGAAT', 'GCGATCGG', 'TGCGATCA')
  if(spid == "SI_3A_E4" || spid == "SI_P03_E4") spseqs <- c('AGCTGACG', 'CCGGTGTC', 'GTAAACAT', 'TATCCTGA')
  if(spid == "SI_3A_F4" || spid == "SI_P03_F4") spseqs <- c('AACGACAC', 'CGTCCTCT', 'GCATGATA', 'TTGATGGG')
  if(spid == "SI_3A_G4" || spid == "SI_P03_G4") spseqs <- c('AGCTTCTC', 'CCTGCGGT', 'GTACAACG', 'TAGAGTAA')
  if(spid == "SI_3A_H4" || spid == "SI_P03_H4") spseqs <- c('ACTGGAGC', 'CGGTCGTG', 'GAAATCAA', 'TTCCATCT')
  if(spid == "SI_3A_A5" || spid == "SI_P03_A5") spseqs <- c('ATTGGGAA', 'CAGTCTGG', 'GGCATACT', 'TCACACTC')
  if(spid == "SI_3A_B5" || spid == "SI_P03_B5") spseqs <- c('ATCGTACT', 'CATCAGTG', 'GGGACTAC', 'TCATGCGA')
  if(spid == "SI_3A_C5" || spid == "SI_P03_C5") spseqs <- c('ATGCATTC', 'CACTGACT', 'GGTACGGG', 'TCAGTCAA')
  if(spid == "SI_3A_D5" || spid == "SI_P03_D5") spseqs <- c('AGACGGAT', 'CCTTTAGA', 'GTCGACTC', 'TAGACTCG')
  if(spid == "SI_3A_E5" || spid == "SI_P03_E5") spseqs <- c('ATCCAAGG', 'CCGTTGAA', 'GGAAGCTC', 'TATGCTCT')
  if(spid == "SI_3A_F5" || spid == "SI_P03_F5") spseqs <- c('AACACAGC', 'CGGTTTAG', 'GTACGGCT', 'TCTGACTA')
  if(spid == "SI_3A_G5" || spid == "SI_P03_G5") spseqs <- c('ATGAAGTA', 'CGCCGAAC', 'GAAGCTCG', 'TCTTTCGT')
  if(spid == "SI_3A_H5" || spid == "SI_P03_H5") spseqs <- c('ATAGTATG', 'CCGCGTCT', 'GGCTCCAC', 'TATAAGGA')
  if(spid == "SI_3A_A6" || spid == "SI_P03_A6") spseqs <- c('ACGGGACT', 'CTTTCGAC', 'GAACATGA', 'TGCATCTG')
  if(spid == "SI_3A_B6" || spid == "SI_P03_B6") spseqs <- c('AACGCGAA', 'CTATTTGG', 'GCGCACCT', 'TGTAGATC')
  if(spid == "SI_3A_C6" || spid == "SI_P03_C6") spseqs <- c('ACTCAGAC', 'CGCTCAGG', 'GAGGTTTA', 'TTAAGCCT')
  if(spid == "SI_3A_D6" || spid == "SI_P03_D6") spseqs <- c('ATGCCAAA', 'CCTTATCG', 'GAAGTCTT', 'TGCAGGGC')
  if(spid == "SI_3A_E6" || spid == "SI_P03_E6") spseqs <- c('ATTGAAAC', 'CAGCCCGA', 'GCCATTTG', 'TGATGGCT')
  if(spid == "SI_3A_F6" || spid == "SI_P03_F6") spseqs <- c('ATGCCGGC', 'CCTAATTA', 'GACTTCCT', 'TGAGGAAG')
  if(spid == "SI_3A_G6" || spid == "SI_P03_G6") spseqs <- c('AGCACTGG', 'CATTACAC', 'GTGCGACA', 'TCAGTGTT')
  if(spid == "SI_3A_H6" || spid == "SI_P03_H6") spseqs <- c('AAGCATAA', 'CCCATCGC', 'GGTTGATG', 'TTAGCGCT')
  if(spid == "SI_3A_A7" || spid == "SI_P03_A7") spseqs <- c('AGGTCATA', 'CTCATCAT', 'GCTGAGGG', 'TAACGTCC')
  if(spid == "SI_3A_B7" || spid == "SI_P03_B7") spseqs <- c('AGGGATGA', 'CTTCTGTT', 'GAATGCAC', 'TCCACACG')
  if(spid == "SI_3A_C7" || spid == "SI_P03_C7") spseqs <- c('ACACCGGG', 'CATAATCC', 'GGCGGAAT', 'TTGTTCTA')
  if(spid == "SI_3A_D7" || spid == "SI_P03_D7") spseqs <- c('AACTTAGA', 'CCGGATCC', 'GGTCGCAT', 'TTAACGTG')
  if(spid == "SI_3A_E7" || spid == "SI_P03_E7") spseqs <- c('AAGACGTG', 'CCATGTGT', 'GTTCACAA', 'TGCGTACC')
  if(spid == "SI_3A_F7" || spid == "SI_P03_F7") spseqs <- c('ACCCGAGA', 'CAAACTTT', 'GGTTAGAC', 'TTGGTCCG')
  if(spid == "SI_3A_G7" || spid == "SI_P03_G7") spseqs <- c('ATTACCGG', 'CAGTAATT', 'GCCGGTAA', 'TGACTGCC')
  if(spid == "SI_3A_H7" || spid == "SI_P03_H7") spseqs <- c('AACGGGTG', 'CTAATTCT', 'GCTTCAAC', 'TGGCACGA')
  if(spid == "SI_3A_A8" || spid == "SI_P03_A8") spseqs <- c('ATGATACG', 'CCACAGAA', 'GACTGTTC', 'TGTGCCGT')
  if(spid == "SI_3A_B8" || spid == "SI_P03_B8") spseqs <- c('ACGTTCAC', 'CAAGGTCT', 'GTTAAGTG', 'TGCCCAGA')
  if(spid == "SI_3A_C8" || spid == "SI_P03_C8") spseqs <- c('AGCTCGAG', 'CAGGAAGA', 'GCACGTTT', 'TTTATCCC')
  if(spid == "SI_3A_D8" || spid == "SI_P03_D8") spseqs <- c('AATCTTTG', 'CTCAAGAC', 'GGATGAGT', 'TCGGCCCA')
  if(spid == "SI_3A_E8" || spid == "SI_P03_E8") spseqs <- c('AGCCTATG', 'CTAACGCA', 'GCTTACAT', 'TAGGGTGC')
  if(spid == "SI_3A_F8" || spid == "SI_P03_F8") spseqs <- c('AGTTGGGA', 'CCAGAAAG', 'GTGCCCTC', 'TACATTCT')
  if(spid == "SI_3A_G8" || spid == "SI_P03_G8") spseqs <- c('AAGTACTC', 'CTTGGAGA', 'GGAACTCT', 'TCCCTGAG')
  if(spid == "SI_3A_H8" || spid == "SI_P03_H8") spseqs <- c('AAGAGCGG', 'CTTGTTAT', 'GGCCCATC', 'TCATAGCA')
  if(spid == "SI_3A_A9" || spid == "SI_P03_A9") spseqs <- c('ACAACTTG', 'CTCCAACA', 'GAGTGCGT', 'TGTGTGAC')
  if(spid == "SI_3A_B9" || spid == "SI_P03_B9") spseqs <- c('AAGCGTGT', 'CTTGACCG', 'GCCACGTA', 'TGATTAAC')
  if(spid == "SI_3A_C9" || spid == "SI_P03_C9") spseqs <- c('AGATCGGT', 'CATCGTCG', 'GTCATATA', 'TCGGACAC')
  if(spid == "SI_3A_D9" || spid == "SI_P03_D9") spseqs <- c('ACCTACTG', 'CAAGGGAC', 'GGGACACA', 'TTTCTTGT')
  if(spid == "SI_3A_E9" || spid == "SI_P03_E9") spseqs <- c('AGTAAGCA', 'CCGGTAAT', 'GTATCTTC', 'TACCGCGG')
  if(spid == "SI_3A_F9" || spid == "SI_P03_F9") spseqs <- c('AGTTAGTT', 'CACGCACG', 'GTACTTAA', 'TCGAGCGC')
  if(spid == "SI_3A_G9" || spid == "SI_P03_G9") spseqs <- c('AGTCTCAG', 'CAATGGCA', 'GCCGAAGT', 'TTGACTTC')
  if(spid == "SI_3A_H9" || spid == "SI_P03_H9") spseqs <- c('ACCTGCCA', 'CTTCATAC', 'GGAATATG', 'TAGGCGGT')
  if(spid == "SI_3A_A10" || spid == "SI_P03_A10") spseqs <- c('ACAGCAAC', 'CGCAATTT', 'GAGTTGCG', 'TTTCGCGA')
  if(spid == "SI_3A_B10" || spid == "SI_P03_B10") spseqs <- c('ACCATTAA', 'CTGGACGT', 'GAACGGTC', 'TGTTCACG')
  if(spid == "SI_3A_C10" || spid == "SI_P03_C10") spseqs <- c('ATCTGATC', 'CGTGCTAA', 'GAGAAGGG', 'TCACTCCT')
  if(spid == "SI_3A_D10" || spid == "SI_P03_D10") spseqs <- c('AGATAACA', 'CTTATTTG', 'GCGGGCAT', 'TACCCGGC')
  if(spid == "SI_3A_E10" || spid == "SI_P03_E10") spseqs <- c('ATTGTTTC', 'CGCAGGAG', 'GCACCAGT', 'TAGTACCA')
  if(spid == "SI_3A_F10" || spid == "SI_P03_F10") spseqs <- c('ATTCGTGC', 'CGCGTGCA', 'GAATACTG', 'TCGACAAT')
  if(spid == "SI_3A_G10" || spid == "SI_P03_G10") spseqs <- c('ATACTGAG', 'CGGAGACT', 'GATGCCTC', 'TCCTATGA')
  if(spid == "SI_3A_H10" || spid == "SI_P03_H10") spseqs <- c('ATTTCAGC', 'CGAGTGAT', 'GACCGCCA', 'TCGAATTG')
  if(spid == "SI_3A_A11" || spid == "SI_P03_A11") spseqs <- c('ACCAGTCC', 'CTTTCCTT', 'GGACAGGG', 'TAGGTAAA')
  if(spid == "SI_3A_B11" || spid == "SI_P03_B11") spseqs <- c('ATGGTCGC', 'CGACATAG', 'GATTCGCT', 'TCCAGATA')
  if(spid == "SI_3A_C11" || spid == "SI_P03_C11") spseqs <- c('ACCGAACA', 'CGACTCTT', 'GTTTGTGG', 'TAGACGAC')
  if(spid == "SI_3A_D11" || spid == "SI_P03_D11") spseqs <- c('ATATGAGA', 'CACCTCAG', 'GCTACTTC', 'TGGGAGCT')
  if(spid == "SI_3A_E11" || spid == "SI_P03_E11") spseqs <- c('ATCGCCAT', 'CATAAAGG', 'GGGTTTCC', 'TCACGGTA')
  if(spid == "SI_3A_F11" || spid == "SI_P03_F11") spseqs <- c('AGCAGTTA', 'CTTGTACC', 'GAACCCGG', 'TCGTAGAT')
  if(spid == "SI_3A_G11" || spid == "SI_P03_G11") spseqs <- c('AGGGCGTT', 'CTATACGC', 'GCTCGTCA', 'TACATAAG')
  if(spid == "SI_3A_H11" || spid == "SI_P03_H11") spseqs <- c('AGGATCGA', 'CACGATTC', 'GTATCGAG', 'TCTCGACT')
  if(spid == "SI_3A_A12" || spid == "SI_P03_A12") spseqs <- c('ACTACTGT', 'CGGGAACG', 'GACCTCTC', 'TTATGGAA')
  if(spid == "SI_3A_B12" || spid == "SI_P03_B12") spseqs <- c('ACGCTTGG', 'CGCTACAT', 'GAAAGACA', 'TTTGCGTC')
  if(spid == "SI_3A_C12" || spid == "SI_P03_C12") spseqs <- c('ATCCGGCA', 'CCGTTATG', 'GGTAATGT', 'TAAGCCAC')
  if(spid == "SI_3A_D12" || spid == "SI_P03_D12") spseqs <- c('AGAAACGT', 'CACTCAAC', 'GCTGTGTA', 'TTGCGTCG')
  if(spid == "SI_3A_E12" || spid == "SI_P03_E12") spseqs <- c('ACGCGGAA', 'CGCTATCC', 'GTTGCATG', 'TAAATCGT')
  if(spid == "SI_3A_F12" || spid == "SI_P03_F12") spseqs <- c('AATTGAAC', 'CCAGTGGA', 'GTCCATTG', 'TGGACCCT')
  if(spid == "SI_3A_G12" || spid == "SI_P03_G12") spseqs <- c('ACCCGCAC', 'CATGCGTA', 'GTGATAGT', 'TGATATCG')
  if(spid == "SI_3A_H12" || spid == "SI_P03_H12") spseqs <- c('ACGAGTAG', 'CAATCCCT', 'GTCCAGGC', 'TGTGTATA')
  
  
  if(spid == "SI_T2_1") spseqs <- c('GGGTGATC', 'TTACCGAT', 'AATGACGA', 'CCCATTCG')
  if(spid == "SI_T2_2") spseqs <- c('GGGTCGAA', 'ATCCGCCC', 'TCTATAGT', 'CAAGATTG')
  if(spid == "SI_T2_3") spseqs <- c('GCTGATAT', 'TGCCGAGC', 'AAATTGCG', 'CTGACCTA')
  if(spid == "SI_T2_4") spseqs <- c('ACTTCTGA', 'TTCATCTT', 'CGACGACG', 'GAGGAGAC')
  if(spid == "SI_T2_5") spseqs <- c('GAATACAA', 'AGCATACC', 'TCGGGTTT', 'CTTCCGGG')
  if(spid == "SI_T2_6") spseqs <- c('TATTGAGA', 'GTAGTCAG', 'CGCCATTC', 'ACGACGCT')
  if(spid == "SI_T2_7") spseqs <- c('AAATCTGT', 'GTCCAACC', 'TCTGGCTG', 'CGGATGAA')
  if(spid == "SI_T2_8") spseqs <- c('CCTTGAAC', 'GAAATCGG', 'TGGCCTCT', 'ATCGAGTA')
  
  
  if(spid == "SI_GA_A1" || spid == "SI_P2_A1") spseqs <- c('GGTTTACT', 'CTAAACGG', 'TCGGCGTC', 'AACCGTAA')
  if(spid == "SI_GA_A2" || spid == "SI_P2_A2") spseqs <- c('TTTCATGA', 'ACGTCCCT', 'CGCATGTG', 'GAAGGAAC')
  if(spid == "SI_GA_A3" || spid == "SI_P2_A3") spseqs <- c('CAGTACTG', 'AGTAGTCT', 'GCAGTAGA', 'TTCCCGAC')
  if(spid == "SI_GA_A4" || spid == "SI_P2_A4") spseqs <- c('TATGATTC', 'CCCACAGT', 'ATGCTGAA', 'GGATGCCG')
  if(spid == "SI_GA_A5" || spid == "SI_P2_A5") spseqs <- c('CTAGGTGA', 'TCGTTCAG', 'AGCCAATT', 'GATACGCC')
  if(spid == "SI_GA_A6" || spid == "SI_P2_A6") spseqs <- c('CGCTATGT', 'GCTGTCCA', 'TTGAGATC', 'AAACCGAG')
  if(spid == "SI_GA_A7" || spid == "SI_P2_A7") spseqs <- c('ACAGAGGT', 'TATAGTTG', 'CGGTCCCA', 'GTCCTAAC')
  if(spid == "SI_GA_A8" || spid == "SI_P2_A8") spseqs <- c('GCATCTCC', 'TGTAAGGT', 'CTGCGATG', 'AACGTCAA')
  if(spid == "SI_GA_A9" || spid == "SI_P2_A9") spseqs <- c('TCTTAAAG', 'CGAGGCTC', 'GTCCTTCT', 'AAGACGGA')
  if(spid == "SI_GA_A10" || spid == "SI_P2_A10") spseqs <- c('GAAACCCT', 'TTTCTGTC', 'CCGTGTGA', 'AGCGAAAG')
  if(spid == "SI_GA_A11" || spid == "SI_P2_A11") spseqs <- c('GTCCGGTC', 'AAGATCAT', 'CCTGAAGG', 'TGATCTCA')
  if(spid == "SI_GA_A12" || spid == "SI_P2_A12") spseqs <- c('AGTGGAAC', 'GTCTCCTT', 'TCACATCA', 'CAGATGGG')
  if(spid == "SI_GA_B1" || spid == "SI_P2_B1") spseqs <- c('GTAATCTT', 'TCCGGAAG', 'AGTTCGGC', 'CAGCATCA')
  if(spid == "SI_GA_B2" || spid == "SI_P2_B2") spseqs <- c('TACTCTTC', 'CCTGTGCG', 'GGACACGT', 'ATGAGAAA')
  if(spid == "SI_GA_B3" || spid == "SI_P2_B3") spseqs <- c('GTGTATTA', 'TGTGCGGG', 'ACCATAAC', 'CAACGCCT')
  if(spid == "SI_GA_B4" || spid == "SI_P2_B4") spseqs <- c('ACTTCATA', 'GAGATGAC', 'TGCCGTGG', 'CTAGACCT')
  if(spid == "SI_GA_B5" || spid == "SI_P2_B5") spseqs <- c('AATAATGG', 'CCAGGGCA', 'TGCCTCAT', 'GTGTCATC')
  if(spid == "SI_GA_B6" || spid == "SI_P2_B6") spseqs <- c('CGTTAATC', 'GCCACGCT', 'TTACTCAG', 'AAGGGTGA')
  if(spid == "SI_GA_B7" || spid == "SI_P2_B7") spseqs <- c('AAACCTCA', 'GCCTTGGT', 'CTGGACTC', 'TGTAGAAG')
  if(spid == "SI_GA_B8" || spid == "SI_P2_B8") spseqs <- c('AAAGTGCT', 'GCTACCTG', 'TGCTGTAA', 'CTGCAAGC')
  if(spid == "SI_GA_B9" || spid == "SI_P2_B9") spseqs <- c('CTGTAACT', 'TCTAGCGA', 'AGAGTGTG', 'GACCCTAC')
  if(spid == "SI_GA_B10" || spid == "SI_P2_B10") spseqs <- c('ACCGTATG', 'GATTAGAT', 'CTGACTGA', 'TGACGCCC')
  if(spid == "SI_GA_B11" || spid == "SI_P2_B11") spseqs <- c('GTTCCTCA', 'AGGTACGC', 'TAAGTATG', 'CCCAGGAT')
  if(spid == "SI_GA_B12" || spid == "SI_P2_B12") spseqs <- c('TACCACCA', 'CTAAGTTT', 'GGGTCAAG', 'ACTGTGGC')
  if(spid == "SI_GA_C1" || spid == "SI_P2_C1") spseqs <- c('CCACTTAT', 'AACTGGCG', 'TTGGCATA', 'GGTAACGC')
  if(spid == "SI_GA_C2" || spid == "SI_P2_C2") spseqs <- c('CCTAGACC', 'ATCTCTGT', 'TAGCTCTA', 'GGAGAGAG')
  if(spid == "SI_GA_C3" || spid == "SI_P2_C3") spseqs <- c('TCAGCCGT', 'CAGAGGCC', 'GGTCAATA', 'ATCTTTAG')
  if(spid == "SI_GA_C4" || spid == "SI_P2_C4") spseqs <- c('ACAATTCA', 'TGCGCAGC', 'CATCACTT', 'GTGTGGAG')
  if(spid == "SI_GA_C5" || spid == "SI_P2_C5") spseqs <- c('CGACTTGA', 'TACAGACT', 'ATTGCGTG', 'GCGTACAC')
  if(spid == "SI_GA_C6" || spid == "SI_P2_C6") spseqs <- c('ATTACTTC', 'TGCGAACT', 'GCATTCGG', 'CAGCGGAA')
  if(spid == "SI_GA_C7" || spid == "SI_P2_C7") spseqs <- c('GTCTCTCG', 'AATCTCTC', 'CGGAGGGA', 'TCAGAAAT')
  if(spid == "SI_GA_C8" || spid == "SI_P2_C8") spseqs <- c('GTTGAGAA', 'AGATCTGG', 'TCGATACT', 'CACCGCTC')
  if(spid == "SI_GA_C9" || spid == "SI_P2_C9") spseqs <- c('GCGCAGAA', 'ATCTTACC', 'TATGGTGT', 'CGAACCTG')
  if(spid == "SI_GA_C10" || spid == "SI_P2_C10") spseqs <- c('TCTCAGTG', 'GAGACTAT', 'CGCTTAGC', 'ATAGGCCA')
  if(spid == "SI_GA_C11" || spid == "SI_P2_C11") spseqs <- c('GAGGATCT', 'AGACCATA', 'TCCTGCGC', 'CTTATGAG')
  if(spid == "SI_GA_C12" || spid == "SI_P2_C12") spseqs <- c('TCTCGTTT', 'GGCTAGCG', 'ATGACCGC', 'CAAGTAAA')
  if(spid == "SI_GA_D1" || spid == "SI_P2_D1") spseqs <- c('CACTCGGA', 'GCTGAATT', 'TGAAGTAC', 'ATGCTCCG')
  if(spid == "SI_GA_D2" || spid == "SI_P2_D2") spseqs <- c('TAACAAGG', 'GGTTCCTC', 'ATCATGCA', 'CCGGGTAT')
  if(spid == "SI_GA_D3" || spid == "SI_P2_D3") spseqs <- c('ACATTACT', 'TTTGGGTA', 'CAGCCCAC', 'GGCAATGG')
  if(spid == "SI_GA_D4" || spid == "SI_P2_D4") spseqs <- c('CCCTAACA', 'ATTCCGAT', 'TGGATTGC', 'GAAGGCTG')
  if(spid == "SI_GA_D5" || spid == "SI_P2_D5") spseqs <- c('CTCGTCAC', 'GATCAGCA', 'ACAACAGG', 'TGGTGTTT')
  if(spid == "SI_GA_D6" || spid == "SI_P2_D6") spseqs <- c('CATGCGAT', 'TGATATTC', 'GTGATCGA', 'ACCCGACG')
  if(spid == "SI_GA_D7" || spid == "SI_P2_D7") spseqs <- c('ATTTGCTA', 'TAGACACC', 'CCACAGGG', 'GGCGTTAT')
  if(spid == "SI_GA_D8" || spid == "SI_P2_D8") spseqs <- c('GCAACAAA', 'TAGTTGTC', 'CGCCATCG', 'ATTGGCGT')
  if(spid == "SI_GA_D9" || spid == "SI_P2_D9") spseqs <- c('AGGAGATG', 'GATGTGGT', 'CTACATCC', 'TCCTCCAA')
  if(spid == "SI_GA_D10" || spid == "SI_P2_D10") spseqs <- c('CAATACCC', 'TGTCTATG', 'ACCACGAA', 'GTGGGTGT')
  if(spid == "SI_GA_D11" || spid == "SI_P2_D11") spseqs <- c('CTTTGCGG', 'TGCACAAA', 'AAGCAGTC', 'GCAGTTCT')
  if(spid == "SI_GA_D12" || spid == "SI_P2_D12") spseqs <- c('GCACAATG', 'CTTGGTAC', 'TGCACCGT', 'AAGTTGCA')
  if(spid == "SI_GA_E1" || spid == "SI_P2_E1") spseqs <- c('TGGTAAAC', 'GAAAGGGT', 'ACTGCTCG', 'CTCCTCTA')
  if(spid == "SI_GA_E2" || spid == "SI_P2_E2") spseqs <- c('GTGGTACC', 'TACTATAG', 'ACAAGGTA', 'CGTCCCGT')
  if(spid == "SI_GA_E3" || spid == "SI_P2_E3") spseqs <- c('AGGTATTG', 'CTCCTAGT', 'TCAAGGCC', 'GATGCCAA')
  if(spid == "SI_GA_E4" || spid == "SI_P2_E4") spseqs <- c('TTCGCCCT', 'GGATGGGC', 'AATCAATG', 'CCGATTAA')
  if(spid == "SI_GA_E5" || spid == "SI_P2_E5") spseqs <- c('CATTAGCG', 'TTCGCTGA', 'ACAAGAAT', 'GGGCTCTC')
  if(spid == "SI_GA_E6" || spid == "SI_P2_E6") spseqs <- c('CTGCGGCT', 'GACTCAAA', 'AGAAACTC', 'TCTGTTGG')
  if(spid == "SI_GA_E7" || spid == "SI_P2_E7") spseqs <- c('CACGCCTT', 'GTATATAG', 'TCTCGGGC', 'AGGATACA')
  if(spid == "SI_GA_E8" || spid == "SI_P2_E8") spseqs <- c('ATAGTTAC', 'TGCTGAGT', 'CCTACGTA', 'GAGCACCG')
  if(spid == "SI_GA_E9" || spid == "SI_P2_E9") spseqs <- c('TTGTTTCC', 'GGAGGAGG', 'CCTAACAA', 'AACCCGTT')
  if(spid == "SI_GA_E10" || spid == "SI_P2_E10") spseqs <- c('AAATGTGC', 'GGGCAAAT', 'TCTATCCG', 'CTCGCGTA')
  if(spid == "SI_GA_E11" || spid == "SI_P2_E11") spseqs <- c('AAGCGCTG', 'CGTTTGAT', 'GTAGCACA', 'TCCAATGC')
  if(spid == "SI_GA_E12" || spid == "SI_P2_E12") spseqs <- c('ACCGGCTC', 'GAGTTAGT', 'CGTCCTAG', 'TTAAAGCA')
  if(spid == "SI_GA_F1" || spid == "SI_P2_F1") spseqs <- c('GTTGCAGC', 'TGGAATTA', 'CAATGGAG', 'ACCCTCCT')
  if(spid == "SI_GA_F2" || spid == "SI_P2_F2") spseqs <- c('TTTACATG', 'CGCGATAC', 'ACGCGGGT', 'GAATTCCA')
  if(spid == "SI_GA_F3" || spid == "SI_P2_F3") spseqs <- c('TTCAGGTG', 'ACGGACAT', 'GATCTTGA', 'CGATCACC')
  if(spid == "SI_GA_F4" || spid == "SI_P2_F4") spseqs <- c('CCCAATAG', 'GTGTCGCT', 'AGAGTCGC', 'TATCGATA')
  if(spid == "SI_GA_F5" || spid == "SI_P2_F5") spseqs <- c('GACTACGT', 'CTAGCGAG', 'TCTATATC', 'AGGCGTCA')
  if(spid == "SI_GA_F6" || spid == "SI_P2_F6") spseqs <- c('CGGAGCAC', 'GACCTATT', 'ACTTAGGA', 'TTAGCTCG')
  if(spid == "SI_GA_F7" || spid == "SI_P2_F7") spseqs <- c('CGTGCAGA', 'AACAAGAT', 'TCGCTTCG', 'GTATGCTC')
  if(spid == "SI_GA_F8" || spid == "SI_P2_F8") spseqs <- c('CATGAACA', 'TCACTCGC', 'AGCTGGAT', 'GTGACTTG')
  if(spid == "SI_GA_F9" || spid == "SI_P2_F9") spseqs <- c('CAAGCTCC', 'GTTCACTG', 'TCGTGAAA', 'AGCATGGT')
  if(spid == "SI_GA_F10" || spid == "SI_P2_F10") spseqs <- c('GCTTGGCT', 'AAACAAAC', 'CGGGCTTA', 'TTCATCGG')
  if(spid == "SI_GA_F11" || spid == "SI_P2_F11") spseqs <- c('GCGAGAGT', 'TACGTTCA', 'AGTCCCAC', 'CTATAGTG')
  if(spid == "SI_GA_F12" || spid == "SI_P2_F12") spseqs <- c('TGATGCAT', 'GCTACTGA', 'CACCTGCC', 'ATGGAATG')
  if(spid == "SI_GA_G1" || spid == "SI_P2_G1") spseqs <- c('ATGAATCT', 'GATCTCAG', 'CCAGGAGC', 'TGCTCGTA')
  if(spid == "SI_GA_G2" || spid == "SI_P2_G2") spseqs <- c('TGATTCTA', 'ACTAGGAG', 'CAGCCACT', 'GTCGATGC')
  if(spid == "SI_GA_G3" || spid == "SI_P2_G3") spseqs <- c('CCTCATTC', 'AGCATCCG', 'GTGGCAAT', 'TAATGGGA')
  if(spid == "SI_GA_G4" || spid == "SI_P2_G4") spseqs <- c('GCGATGTG', 'AGATACAA', 'TTTCCACT', 'CACGGTGC')
  if(spid == "SI_GA_G5" || spid == "SI_P2_G5") spseqs <- c('GAGCAAGA', 'TCTGTGAT', 'CGCAGTTC', 'ATATCCCG')
  if(spid == "SI_GA_G6" || spid == "SI_P2_G6") spseqs <- c('CTGACGCG', 'GGTCGTAC', 'TCCTTCTT', 'AAAGAAGA')
  if(spid == "SI_GA_G7" || spid == "SI_P2_G7") spseqs <- c('GGTATGCA', 'CTCGAAAT', 'ACACCTTC', 'TAGTGCGG')
  if(spid == "SI_GA_G8" || spid == "SI_P2_G8") spseqs <- c('TATGAGCT', 'CCGATAGC', 'ATACCCAA', 'GGCTGTTG')
  if(spid == "SI_GA_G9" || spid == "SI_P2_G9") spseqs <- c('TAGGACGT', 'ATCCCACA', 'GGAATGTC', 'CCTTGTAG')
  if(spid == "SI_GA_G10" || spid == "SI_P2_G10") spseqs <- c('TCGCCAGC', 'AATGTTAG', 'CGATAGCT', 'GTCAGCTA')
  if(spid == "SI_GA_G11" || spid == "SI_P2_G11") spseqs <- c('TTATCGTT', 'AGCAGAGC', 'CATCTCCA', 'GCGGATAG')
  if(spid == "SI_GA_G12" || spid == "SI_P2_G12") spseqs <- c('ATTCTAAG', 'CCCGATTA', 'TGGAGGCT', 'GAATCCGC')
  if(spid == "SI_GA_H1" || spid == "SI_P2_H1") spseqs <- c('GTATGTCA', 'TGTCAGAC', 'CACGTCGG', 'ACGACATT')
  if(spid == "SI_GA_H2" || spid == "SI_P2_H2") spseqs <- c('TAATGACC', 'ATGCCTTA', 'GCCGAGAT', 'CGTATCGG')
  if(spid == "SI_GA_H3" || spid == "SI_P2_H3") spseqs <- c('CCAAGATG', 'AGGCCCGA', 'TACGTGAC', 'GTTTATCT')
  if(spid == "SI_GA_H4" || spid == "SI_P2_H4") spseqs <- c('GCCATTCC', 'CAAGAATT', 'TTGCCGGA', 'AGTTGCAG')
  if(spid == "SI_GA_H5" || spid == "SI_P2_H5") spseqs <- c('CCACTACA', 'GATTCTGG', 'TGCGGCTT', 'ATGAAGAC')
  if(spid == "SI_GA_H6" || spid == "SI_P2_H6") spseqs <- c('TAGGATAA', 'CCTTTGTC', 'GTACGCGG', 'AGCACACT')
  if(spid == "SI_GA_H7" || spid == "SI_P2_H7") spseqs <- c('AGCTATCA', 'CATATAAC', 'TCAGGGTG', 'GTGCCCGT')
  if(spid == "SI_GA_H8" || spid == "SI_P2_H8") spseqs <- c('TTGTTGAT', 'GCTCAACC', 'CAAAGTGG', 'AGCGCCTA')
  if(spid == "SI_GA_H9" || spid == "SI_P2_H9") spseqs <- c('ACACTGTT', 'CAGGATGG', 'GGCTGAAC', 'TTTACCCA')
  if(spid == "SI_GA_H10" || spid == "SI_P2_H10") spseqs <- c('GTAATTGC', 'AGTCGCTT', 'CACGAGAA', 'TCGTCACG')
  if(spid == "SI_GA_H11" || spid == "SI_P2_H11") spseqs <- c('GGCGAGTA', 'ACTTCTAT', 'CAAATACG', 'TTGCGCGC')
  if(spid == "SI_GA_H12" || spid == "SI_P2_H12") spseqs <- c('GACAGCAT', 'TTTGTACA', 'AGGCCGTG', 'CCATATGC')

  if(spid == "SI-NN-C2") spseqs <- c("GACGTACGGA", "ACCTCGCTTT")
  if(spid == "SI-NN-C3") spseqs <- c("CTTAGATGCA", "CTATTCGACT")
  if(spid == "SI-NN-C4") spseqs <- c("ACTCCTCAAC", "ATTCACGGGA")
  if(spid == "SI-NN-C5") spseqs <- c("TCCATACTGT", "GAGCGACTTC")
  if(spid == "SI-NN-C6") spseqs <- c("GCCAGAATGG", "GAATTGTATC")
  if(spid == "SI-NN-C7") spseqs <- c("ATTGCAAGAC", "CGGTGACCAT")
  if(spid == "SI-NN-C8") spseqs <- c("CAACATATCG", "AACTAAAGGC")
  if(spid == "SI-NN-C9") spseqs <- c("CAGGGACCTC", "AATAGCAGCA")
  if(spid == "SI-NN-C10") spseqs <- c("ATTTGGGAAT", "TCAGCCGAGT")
  if(spid == "SI-NN-C11") spseqs <- c("TGGATAATTG", "TTTGTGATCA")
  if(spid == "SI-NN-C12") spseqs <- c("CAGATTATGA", "GTGCCTGTAC")
  if(spid == "SI-NN-D1") spseqs <- c("AATGGTAAGC", "CTGTGAAGTT")
  if(spid == "SI-NN-D2") spseqs <- c("CATCCCTCGG", "AGGGTAGGCG")
  if(spid == "SI-NN-D3") spseqs <- c("CTGGTCCGTC", "ACGAATAAGG")
  if(spid == "SI-NN-D4") spseqs <- c("GTCGTAAATT", "GTAGTAGGCT")
  if(spid == "SI-NN-D5") spseqs <- c("CAACACTGAT", "CAAACCAGGG")
  if(spid == "SI-NN-D6") spseqs <- c("CAGACTGAAT", "TGCCTAACCG")
  if(spid == "SI-NN-D7") spseqs <- c("CGGACACGCT", "TCACGATATA")
  if(spid == "SI-NN-D8") spseqs <- c("CTCGGCTCTC", "TATACTTCCA")
  if(spid == "SI-NN-D9") spseqs <- c("AAGCCTATCA", "GTAGTGAACA")
  if(spid == "SI-NN-D10") spseqs <- c("ACTGTGCCAA", "TAACCGTTAT")
  if(spid == "SI-NN-D11") spseqs <- c("CTACTCCCAT", "TACGTTGCGG")
  if(spid == "SI-NN-D12") spseqs <- c("AGACTCTAAC", "GTGTCGATGT")
  if(spid == "SI-NN-E1") spseqs <- c("TGAGCCGGCA", "TGCTGCATGT")
  if(spid == "SI-NN-E2") spseqs <- c("TAGCTTGCGT", "GCTGCTAAAG")
  if(spid == "SI-NN-E3") spseqs <- c("GTACACCGGG", "AACAGTCGCG")
  if(spid == "SI-NN-E4") spseqs <- c("CGCGTAGACG", "AAAGTAGCTG")
  if(spid == "SI-NN-E5") spseqs <- c("CCGATATATT", "AACTATCCGA")
  if(spid == "SI-NN-E6") spseqs <- c("TTGGGCGGGA", "ATGTATCTGT")
  if(spid == "SI-NN-E7") spseqs <- c("GTTTATGAGG", "GCCAAACAGC")
  if(spid == "SI-NN-E8") spseqs <- c("TTGGATACTC", "CGTACTGAAG")
  if(spid == "SI-NN-E9") spseqs <- c("GCTATCTATC", "AAAGATGGAT")
  if(spid == "SI-NN-E10") spseqs <- c("ATGTGATCGC", "TGACCTACCT")
  if(spid == "SI-NN-E11") spseqs <- c("TCCTCACATG", "AACAGATTCA")
  if(spid == "SI-NN-E12") spseqs <- c("GCCTCCTAAT", "CTCCTCCTGT")
  if(spid == "SI-NN-F1") spseqs <- c("CGAGGCTGAT", "TCGACTTTCT")
  if(spid == "SI-NN-F2") spseqs <- c("CTCATGACAT", "TGTTCCGCAC")
  if(spid == "SI-NN-F3") spseqs <- c("AACGCGGCAA", "GCCCGTTTAC")
  if(spid == "SI-NN-F4") spseqs <- c("AACACCTCTC", "GATTGTATGA")
  if(spid == "SI-NN-F5") spseqs <- c("AAATAACGCG", "ATAGAGGAGC")
  if(spid == "SI-NN-F6") spseqs <- c("GAAATCTTGT", "TAGACGCCAC")
  if(spid == "SI-NN-F7") spseqs <- c("TAGATGGACT", "CCAGGGATAT")
  if(spid == "SI-NN-F8") spseqs <- c("CAATCTGTAT", "TCCCATCGTT")
  if(spid == "SI-NN-F9") spseqs <- c("GCGAAAGTAT", "GTGTCACAGA")
  if(spid == "SI-NN-F10") spseqs <- c("CCCAGAGCTT", "ATGCTGTAAT")
  if(spid == "SI-NN-F11") spseqs <- c("GACTACAGGA", "AGTGACAAAG")
  if(spid == "SI-NN-F12") spseqs <- c("CAATTCATGG", "ACGTAGTTCG")
  if(spid == "SI-NN-G1") spseqs <- c("TGGGAAGGGC", "ACTAACGATG")
  if(spid == "SI-NN-G2") spseqs <- c("TACTTCGCAT", "AACGGTGTCG")
  if(spid == "SI-NN-G3") spseqs <- c("TAAACTGGCA", "CCGACCTCTT")
  if(spid == "SI-NN-G4") spseqs <- c("TTGAAGAACT", "CAACACCAGG")
  if(spid == "SI-NN-G5") spseqs <- c("TCCCTCGTCA", "ATCAGGTGTG")
  if(spid == "SI-NN-G6") spseqs <- c("ACTGGCCTAA", "CGCTCCCTTG")
  if(spid == "SI-NN-G7") spseqs <- c("TCATAATGCA", "GATGAGGCGA")
  if(spid == "SI-NN-G8") spseqs <- c("CGCCACGTTA", "CTTTCATAGG")
  if(spid == "SI-NN-G9") spseqs <- c("TGGACTAAGC", "CTTAACATCG")
  if(spid == "SI-NN-G10") spseqs <- c("AAGCCCGGCA", "GAGTTTCCGT")
  if(spid == "SI-NN-G11") spseqs <- c("ATAGTGGCGG", "TACTGTTGTT")
  if(spid == "SI-NN-G12") spseqs <- c("CCAACCGGAC", "GTTCTATCTA")
  if(spid == "SI-NN-H1") spseqs <- c("TTACAGAGGG", "GTTTATGGCA")
  if(spid == "SI-NN-H2") spseqs <- c("TCCAACAACG", "GTCCGGGTTT")
  if(spid == "SI-NN-H3") spseqs <- c("AAAGCATGTC", "AGATACTACT")
  if(spid == "SI-NN-H4") spseqs <- c("AGAAGTCGTT", "AGGACATTTA")
  if(spid == "SI-NN-H5") spseqs <- c("GCTAGCGTTC", "TCGCGTGGTG")
  if(spid == "SI-NN-H6") spseqs <- c("CAACTAGCTT", "TCAACTAATC")
  if(spid == "SI-NN-H7") spseqs <- c("AGGGCCAGTT", "TAGGGCAAAT")
  if(spid == "SI-NN-H8") spseqs <- c("CTCTATGCGT", "CCTCTGATGC")
  if(spid == "SI-NN-H9") spseqs <- c("GTCCATACAA", "AGCCACTCAC")
  if(spid == "SI-NN-H10") spseqs <- c("AGGGATGCAA", "CTCATATCCT")
  if(spid == "SI-NN-H11") spseqs <- c("TTCAAGTCCT", "CAGCAGCTAA")
  if(spid == "SI-NN-H12") spseqs <- c("CTCCAACCTA", "GCACCGCATC")
  if(spid == "SI-NT-E4") spseqs <- c("CAGCTCCAAT", "TAGGAGGGTG")
  if(spid == "SI-NT-E5") spseqs <- c("TTTGGAGAAA", "CAGCCCTTGA")
  if(spid == "SI-NT-E6") spseqs <- c("TGTCGCTTTC", "CCTAAAGGGT")
  if(spid == "SI-NT-E7") spseqs <- c("GTGTTGCGAG", "CACTTTGAGC")
  if(spid == "SI-NT-E8") spseqs <- c("GTCTGTTCTG", "TCTGATGCTT")
  if(spid == "SI-NT-E9") spseqs <- c("TCGCTCAATC", "GCTTACGTGA")
  if(spid == "SI-NT-E10") spseqs <- c("CGTTACTAGC", "AATGGCCATT")
  if(spid == "SI-NT-E11") spseqs <- c("GTGCAGATTT", "TAAAGGCGAT")
  if(spid == "SI-NT-E12") spseqs <- c("GTTACCGAAA", "CGTTTGAATC")
  if(spid == "SI-NT-F1") spseqs <- c("CGGTATGGAA", "TGCAACTGGT")
  if(spid == "SI-NT-F2") spseqs <- c("TCGGTTTAGT", "TAAAGTGTTC")
  if(spid == "SI-NT-F3") spseqs <- c("CTTTGCTTTG", "ACATTCATTC")
  if(spid == "SI-NT-F4") spseqs <- c("ATGATCGCAC", "GCGGTATGTA")
  if(spid == "SI-NT-F5") spseqs <- c("CGGCGTGTTA", "TGAGACGCGG")
  if(spid == "SI-NT-F6") spseqs <- c("TCATTTCCGC", "ACAGCTGTTG")
  if(spid == "SI-NT-F7") spseqs <- c("TTGTGCCCAG", "TTGCGAGTTG")
  if(spid == "SI-NT-F8") spseqs <- c("AGCCAACGTG", "TCAAAGGAAC")
  if(spid == "SI-NT-F9") spseqs <- c("GATCCATACT", "ATCGAGAGAA")
  if(spid == "SI-NT-F10") spseqs <- c("TTTCCATGAA", "CGGTTCTAAC")
  if(spid == "SI-NT-F11") spseqs <- c("TCTTTATCGG", "CCCAATAGCG")
  if(spid == "SI-NT-F12") spseqs <- c("GAACCCGATG", "TCTATCCGGG")
  if(spid == "SI-NT-G1") spseqs <- c("GATGTCTGTG", "CGTGCTACCG")
  if(spid == "SI-NT-G2") spseqs <- c("TGTTAGACTA", "TCATGACGAC")
  if(spid == "SI-NT-G3") spseqs <- c("TAGAACGCTT", "GTCGCAAACG")
  if(spid == "SI-NT-G4") spseqs <- c("AATTTGACGT", "AGAAGTGAGG")
  if(spid == "SI-NT-G5") spseqs <- c("GAGTGGCCAC", "GATGCACAGT")
  if(spid == "SI-NT-G6") spseqs <- c("ATCAGTTACG", "TACAATCTCT")
  if(spid == "SI-NT-G7") spseqs <- c("TGTAAGCTTT", "CGGACTAGTC")
  if(spid == "SI-NT-G8") spseqs <- c("GCTTCGGTCT", "CGCGGGTTAG")
  if(spid == "SI-NT-G9") spseqs <- c("AGTCTAAGAC", "CTGAAGTGAA")
  if(spid == "SI-NT-G10") spseqs <- c("TGAAACATTC", "CCTCCTGCCT")
  if(spid == "SI-NT-G11") spseqs <- c("CCTCACGAAA", "TGCCCTGTGC")
  if(spid == "SI-NT-G12") spseqs <- c("CTACGGGCTT", "AGGCGTCCCT")
  if(spid == "SI-NT-H1") spseqs <- c("GCGGCGTAAA", "CACGCGATCA")
  if(spid == "SI-NT-H2") spseqs <- c("TCTCCTTCGG", "TGTCCGTCGC")
  if(spid == "SI-NT-H3") spseqs <- c("CCTTACCCAT", "TCCCGGCAAC")
  if(spid == "SI-NT-H4") spseqs <- c("AAGGAACATC", "TATCTGTGGG")
  if(spid == "SI-NT-H5") spseqs <- c("GTGGGAACTT", "GCACTTACAG")
  if(spid == "SI-NT-H6") spseqs <- c("ACCGCACCAC", "TCGGATTATA")
  if(spid == "SI-NT-H7") spseqs <- c("ACCTTTATCT", "CGAGTGCCAA")
  if(spid == "SI-NT-H8") spseqs <- c("AGTAGCCCGT", "TTCGGGACCT")
  if(spid == "SI-NT-H9") spseqs <- c("AGATAGCATA", "CATTTCCGGA")
  if(spid == "SI-NT-H10") spseqs <- c("TGGCGTTAAA", "TGTTAAGATG")
  if(spid == "SI-NT-H11") spseqs <- c("TCAGGCGAAA", "TCTCTGGAAG")
  if(spid == "SI-NT-H12") spseqs <- c("AAGACATAGC", "ATGTGGACGT")
  if(spid == "SI-NN-A1") spseqs <- c("GCCTTGTCAA", "GATCCAATCA")
  if(spid == "SI-NN-A2") spseqs <- c("AACTTCTTTG", "ATAAGCGGTA")
  if(spid == "SI-NN-A3") spseqs <- c("TGGTTATAGA", "GAGCTCATAA")
  if(spid == "SI-NN-A4") spseqs <- c("GTAGTGCATC", "ATTTGACGGT")
  if(spid == "SI-NN-A5") spseqs <- c("GCTTGATTTA", "CGGGTGTGAC")
  if(spid == "SI-NN-A6") spseqs <- c("GAAGTTTCGC", "TTTCAAAGCA")
  if(spid == "SI-NN-A7") spseqs <- c("GCAGTTGTTT", "ACACTACTTT")
  if(spid == "SI-NN-A8") spseqs <- c("GCCGCAAGTA", "AAAGGTCCGC")
  if(spid == "SI-NN-A9") spseqs <- c("AAGTTGTATG", "TGTACAATGC")
  if(spid == "SI-NN-A10") spseqs <- c("AGCAACCTTT", "TTCAGCGCTA")
  if(spid == "SI-NN-A11") spseqs <- c("GACGAGGGTT", "TTACTTTGAC")
  if(spid == "SI-NN-A12") spseqs <- c("TGAACAGTCG", "ATTTCTCAGA")
  if(spid == "SI-NN-B1") spseqs <- c("CGAAAGGAGT", "ACACCAAACA")
  if(spid == "SI-NN-B2") spseqs <- c("GAATGACTTG", "TGAGCTCGCT")
  if(spid == "SI-NN-B3") spseqs <- c("GCCAGGGATT", "CAAGTTTAGG")
  if(spid == "SI-NN-B4") spseqs <- c("CTTTGTCCGC", "GCTGTCCCGT")
  if(spid == "SI-NN-B5") spseqs <- c("AGCTTGGGTG", "AGCGTATAGT")
  if(spid == "SI-NN-B6") spseqs <- c("TCGCACCGGT", "TAGCCTCAGG")
  if(spid == "SI-NN-B7") spseqs <- c("TTTAGGTAGG", "ACTGAGGGAC")
  if(spid == "SI-NN-B8") spseqs <- c("GTTCGAGCGG", "CTGTGGGTTT")
  if(spid == "SI-NN-B9") spseqs <- c("AGGGTATGGC", "GAGTCCGACT")
  if(spid == "SI-NN-B10") spseqs <- c("TACCCACTAT", "TGGAACATAT")
  if(spid == "SI-NN-B11") spseqs <- c("TTCAGGCTAC", "CTGGGTATTG")
  if(spid == "SI-NN-B12") spseqs <- c("AGGAATAGAT", "TGCGGTCACG")
  if(spid == "SI-NN-C1") spseqs <- c("ACAGCGCTGT", "GAGGGTCATC")
  if(spid == "SI-TT-G10") spseqs <- c("ACTTTACGTG", "TGAACGCCCT")
  if(spid == "SI-TT-G11") spseqs <- c("GATAACCTGC", "CATTAGAAAC")
  if(spid == "SI-TT-G12") spseqs <- c("CTTGCATAAA", "ATCAGGGCTT")
  if(spid == "SI-TT-H1") spseqs <- c("ACAATGTGAA", "CGTACCGTTA")
  if(spid == "SI-TT-H2") spseqs <- c("TAGCATAGTG", "CGGCTCTGTC")
  if(spid == "SI-TT-H3") spseqs <- c("CCCGTTCTCG", "GACGGATTGG")
  if(spid == "SI-TT-H4") spseqs <- c("AGTTTCCTGG", "TGCCACACAG")
  if(spid == "SI-TT-H5") spseqs <- c("AGCAAGAAGC", "TTGTGTTTCT")
  if(spid == "SI-TT-H6") spseqs <- c("CCTATCCTCG", "GAATACTAAC")
  if(spid == "SI-TT-H7") spseqs <- c("ACCTCGAGCT", "TGTGTTCGAT")
  if(spid == "SI-TT-H8") spseqs <- c("ATAAGGATAC", "ATAGATAGGG")
  if(spid == "SI-TT-H9") spseqs <- c("AGAACTTAGA", "CGAGTCCTTT")
  if(spid == "SI-TT-H10") spseqs <- c("TTATCTAGGG", "AAAGGCTCTA")
  if(spid == "SI-TT-H11") spseqs <- c("ACAATCGATC", "TGACGGAATG")
  if(spid == "SI-TT-H12") spseqs <- c("TGATGATTCA", "GTAGGAGTCG")
  if(spid == "SI-NT-A1") spseqs <- c("ATTTACCGCA", "GACAATAAAG")
  if(spid == "SI-NT-A2") spseqs <- c("TTGTCGTAGA", "CAATGTAGCA")
  if(spid == "SI-NT-A3") spseqs <- c("AGTCCTGCGG", "TGACACAAGT")
  if(spid == "SI-NT-A4") spseqs <- c("TTGTTCATGT", "TGACAGCTGA")
  if(spid == "SI-NT-A5") spseqs <- c("TCAGGAAGGA", "GAACGTGCTT")
  if(spid == "SI-NT-A6") spseqs <- c("CTGTTAGAGG", "CCACGCTTCG")
  if(spid == "SI-NT-A7") spseqs <- c("AGATGAGAAT", "GTCGACGGGT")
  if(spid == "SI-NT-A8") spseqs <- c("CCAAAGCCGG", "ACCGTGCACA")
  if(spid == "SI-NT-A9") spseqs <- c("AGTCATAATG", "AGGCTTGAAA")
  if(spid == "SI-NT-A10") spseqs <- c("TTCATCAGAG", "TTGTCGTCTC")
  if(spid == "SI-NT-A11") spseqs <- c("GAAGCGCGAA", "CAGCGAAATT")
  if(spid == "SI-NT-A12") spseqs <- c("ATCCGCCGAA", "GCTACAGAAT")
  if(spid == "SI-NT-B1") spseqs <- c("CTTATTGTGG", "AGCTGTGGGT")
  if(spid == "SI-NT-B2") spseqs <- c("ATTCGTTGGG", "CATCAGGAGC")
  if(spid == "SI-NT-B3") spseqs <- c("GTGGCCTCAT", "TCGAAAGTGA")
  if(spid == "SI-NT-B4") spseqs <- c("CTATGGCATC", "CTCTGAGCGC")
  if(spid == "SI-NT-B5") spseqs <- c("AAACCACAGT", "CCGCAAATGG")
  if(spid == "SI-NT-B6") spseqs <- c("ATGTATCCAC", "CAGGCTGAGG")
  if(spid == "SI-NT-B7") spseqs <- c("GATGAGTCTG", "TATGTAACCG")
  if(spid == "SI-NT-B8") spseqs <- c("CTCACAATAA", "CCGTGTTTAA")
  if(spid == "SI-NT-B9") spseqs <- c("TCTGTCGCAA", "ATCCGAACTG")
  if(spid == "SI-NT-B10") spseqs <- c("CATCGAGAAG", "GTGCATCCGC")
  if(spid == "SI-NT-B11") spseqs <- c("ACGTATTGGG", "CAATCGCAAA")
  if(spid == "SI-NT-B12") spseqs <- c("CTTACGCGAC", "CCTGTCGGAT")
  if(spid == "SI-NT-C1") spseqs <- c("GATAGTGAAG", "CGTTCGCCAG")
  if(spid == "SI-NT-C2") spseqs <- c("CTGTGCAATG", "ATGGTGCTTT")
  if(spid == "SI-NT-C3") spseqs <- c("AAAGCTGAGC", "CAGGAACGAG")
  if(spid == "SI-NT-C4") spseqs <- c("GTCCAAGTCG", "GTCATGGCAC")
  if(spid == "SI-NT-C5") spseqs <- c("GATTTCCATC", "CTTATGGTCA")
  if(spid == "SI-NT-C6") spseqs <- c("CTAAATGGAT", "AATCACATAC")
  if(spid == "SI-NT-C7") spseqs <- c("ATTTCAGTTC", "GCAGCATTAA")
  if(spid == "SI-NT-C8") spseqs <- c("TTGCGACGTC", "TTAATCCACA")
  if(spid == "SI-NT-C9") spseqs <- c("GTACTGTCCA", "TTTGTTGGAA")
  if(spid == "SI-NT-C10") spseqs <- c("AAACGGTTTA", "TCTTCGTTAC")
  if(spid == "SI-NT-C11") spseqs <- c("CGATAGTTCT", "GCGATTTCCT")
  if(spid == "SI-NT-C12") spseqs <- c("TACTTTAGCT", "GAGTTATTTG")
  if(spid == "SI-NT-D1") spseqs <- c("CTGGGATTAA", "AGTAAAGTTC")
  if(spid == "SI-NT-D2") spseqs <- c("TACACTATTC", "ATTGCAACGT")
  if(spid == "SI-NT-D3") spseqs <- c("AGAGCTATTT", "CTGGCACCCA")
  if(spid == "SI-NT-D4") spseqs <- c("AGGTCGTTAT", "AGGTTATCCA")
  if(spid == "SI-NT-D5") spseqs <- c("GCTGGACAGG", "TAATCTCTGT")
  if(spid == "SI-NT-D6") spseqs <- c("ACTTAGATGG", "AATGTTGAGT")
  if(spid == "SI-NT-D7") spseqs <- c("AGGCGCGAGT", "TGTGTCTATA")
  if(spid == "SI-NT-D8") spseqs <- c("GAGGTCGTAC", "CAACTGCGGG")
  if(spid == "SI-NT-D9") spseqs <- c("CCCAATGAGC", "AGTTTCTGAT")
  if(spid == "SI-NT-D10") spseqs <- c("ACATTCAGGG", "GTCCTGAGAG")
  if(spid == "SI-NT-D11") spseqs <- c("AGAAACGGTG", "TCACCCAGCA")
  if(spid == "SI-NT-D12") spseqs <- c("AGTTGCAAGT", "ATACTTCATG")
  if(spid == "SI-NT-E1") spseqs <- c("GCCTAGAAAT", "ACTACTCCGG")
  if(spid == "SI-NT-E2") spseqs <- c("TACAAAGATG", "CCGATGACGT")
  if(spid == "SI-NT-E3") spseqs <- c("GCAGGCTTAG", "CGTCTATTGC")
  if(spid == "SI-NT-E4") spseqs <- c("CAGCTCCAAT", "TAGGAGGGTG")
  if(spid == "SI-TT-B5") spseqs <- c("TCGGCTCTAC", "CCGATGGTCT")
  if(spid == "SI-TT-B6") spseqs <- c("AATGCCATGA", "TACGTAATGC")
  if(spid == "SI-TT-B7") spseqs <- c("GCCTTCGGTA", "CCAACGATTT")
  if(spid == "SI-TT-B8") spseqs <- c("GCACTGAGAA", "TATGCGTGAA")
  if(spid == "SI-TT-B9") spseqs <- c("TATTGAGGCA", "CAGGTAAGTG")
  if(spid == "SI-TT-B10") spseqs <- c("GCCCGATGGA", "AATCGTCTAG")
  if(spid == "SI-TT-B11") spseqs <- c("TCTTACTTGC", "TGACCTCTAG")
  if(spid == "SI-TT-B12") spseqs <- c("CGTCAAGGGC", "TAGGTCACTC")
  if(spid == "SI-TT-C1") spseqs <- c("TGCGCGGTTT", "CAAGGATAAA")
  if(spid == "SI-TT-C2") spseqs <- c("CAATCCCGAC", "CCGAGTAGTA")
  if(spid == "SI-TT-C3") spseqs <- c("ATGGCTTGTG", "GAATGTTGTG")
  if(spid == "SI-TT-C4") spseqs <- c("TTCTCGATGA", "TGTCGGGCAC")
  if(spid == "SI-TT-C5") spseqs <- c("TCCGTTGGAT", "ACGTTCTCGC")
  if(spid == "SI-TT-C6") spseqs <- c("ACGACTACCA", "ACGACCCTAA")
  if(spid == "SI-TT-C7") spseqs <- c("CGCGCACTTA", "CCTGTATTCT")
  if(spid == "SI-TT-C8") spseqs <- c("GCTACAAAGC", "CACGTGCCCT")
  if(spid == "SI-TT-C9") spseqs <- c("TATCAGCCTA", "GTTTCGTCCT")
  if(spid == "SI-TT-C10") spseqs <- c("AGAATGGTTT", "GAGGGTGGGA")
  if(spid == "SI-TT-C11") spseqs <- c("ATGGGTGAAA", "CTTGGGAATT")
  if(spid == "SI-TT-C12") spseqs <- c("TCGTCAAGAT", "GCAACTCAGG")
  if(spid == "SI-TT-D1") spseqs <- c("TGCAATGTTC", "GCTTGTCGAA")
  if(spid == "SI-TT-D2") spseqs <- c("TTAATACGCG", "CACCTCGGGT")
  if(spid == "SI-TT-D3") spseqs <- c("CCTTCTAGAG", "AATACAACGA")
  if(spid == "SI-TT-D4") spseqs <- c("GCAGTATAGG", "TTCCGTGCAC")
  if(spid == "SI-TT-D5") spseqs <- c("TGGTTCGGGT", "GTGGCAGGAG")
  if(spid == "SI-TT-D6") spseqs <- c("CCCAGCTTCT", "GACACCAAAC")
  if(spid == "SI-TT-D7") spseqs <- c("CCTGTCAGGG", "AGCCCGTAAC")
  if(spid == "SI-TT-D8") spseqs <- c("CGCTGAAATC", "AGGTGTCTGC")
  if(spid == "SI-TT-D9") spseqs <- c("TGGTCCCAAG", "CCTCTGGCGT")
  if(spid == "SI-TT-D10") spseqs <- c("ATGCGAATGG", "ACAAGTGTCG")
  if(spid == "SI-TT-D11") spseqs <- c("CGAATATTCG", "CTGGAAGCAA")
  if(spid == "SI-TT-D12") spseqs <- c("GAATTGGTTA", "ACTCTAGTAG")
  if(spid == "SI-TT-E1") spseqs <- c("TTATTCGAGG", "CTGTCCTGCT")
  if(spid == "SI-TT-E2") spseqs <- c("ATGGAGGGAG", "ATAACCCATT")
  if(spid == "SI-TT-E3") spseqs <- c("ACCAGACAAC", "AGGAACTAGG")
  if(spid == "SI-TT-E4") spseqs <- c("AACCACGCAT", "ATTCAGGTTA")
  if(spid == "SI-TT-E5") spseqs <- c("CGCGGTAGGT", "CAGGATGTTG")
  if(spid == "SI-TT-E6") spseqs <- c("TTGAGAGTCA", "AACCTGGTAG")
  if(spid == "SI-TT-E7") spseqs <- c("GTCCTTCGGC", "TCATGCACAG")
  if(spid == "SI-TT-E8") spseqs <- c("GAGCAAGGGC", "ATTGACTTGG")
  if(spid == "SI-TT-E9") spseqs <- c("TGTCCCAACG", "TCGATGTCCA")
  if(spid == "SI-TT-E10") spseqs <- c("CACAATCCCA", "ATATCCACAA")
  if(spid == "SI-TT-E11") spseqs <- c("TCCGGGACAA", "GTGAATGCCA")
  if(spid == "SI-TT-E12") spseqs <- c("CGTCCACCTG", "CATTCATGAC")
  if(spid == "SI-TT-F1") spseqs <- c("AAGATTGGAT", "AGCGGGATTT")
  if(spid == "SI-TT-F2") spseqs <- c("AAGGGCCGCA", "CTGATTCCTC")
  if(spid == "SI-TT-F3") spseqs <- c("GAGAGGATAT", "TTGAAATGGG")
  if(spid == "SI-TT-F4") spseqs <- c("CCCACCACAA", "ACCTCCGCTT")
  if(spid == "SI-TT-F5") spseqs <- c("CGGCTGGATG", "TGATAAGCAC")
  if(spid == "SI-TT-F6") spseqs <- c("TTGCCCGTGC", "GCGTGAGATT")
  if(spid == "SI-TT-F7") spseqs <- c("AATGTATCCA", "AATGAGCTTA")
  if(spid == "SI-TT-F8") spseqs <- c("CTCCTTTAGA", "GACATAGCTC")
  if(spid == "SI-TT-F9") spseqs <- c("GTCCCATCAA", "CGAACGTGAC")
  if(spid == "SI-TT-F10") spseqs <- c("CCGGCAACTG", "CGGTTTAACA")
  if(spid == "SI-TT-F11") spseqs <- c("TTCACACCTT", "TAGTGTACAC")
  if(spid == "SI-TT-F12") spseqs <- c("GAGACGCACG", "CTATGAACAT")
  if(spid == "SI-TT-G1") spseqs <- c("TGTAGTCATT", "CTTGATCGTA")
  if(spid == "SI-TT-G2") spseqs <- c("CATGTGGGTT", "GATTCCTTTA")
  if(spid == "SI-TT-G3") spseqs <- c("ATGACGTCGC", "AGGTCAGGAT")
  if(spid == "SI-TT-G4") spseqs <- c("GCGCTTATGG", "GCCTGGCTAG")
  if(spid == "SI-TT-G5") spseqs <- c("ATAGGGCGAG", "TGCATCGAGT")
  if(spid == "SI-TT-G6") spseqs <- c("GCGGGTAAGT", "TAGCACTAAG")
  if(spid == "SI-TT-G7") spseqs <- c("GTTTCACGAT", "TTCGGCCAAA")
  if(spid == "SI-TT-G8") spseqs <- c("TAAGCAACTG", "CTATACTCAA")
  if(spid == "SI-TT-G9") spseqs <- c("CCGGAGGAAG", "TGCGGATGTT")
  if(spid == "SI-TT-A1") spseqs <- c("GTAACATGCG", "AGTGTTACCT")
  if(spid == "SI-TT-A2") spseqs <- c("GTGGATCAAA", "GCCAACCCTG")
  if(spid == "SI-TT-A3") spseqs <- c("CACTACGAAA", "TTAGACTGAT")
  if(spid == "SI-TT-A4") spseqs <- c("CTCTAGCGAG", "TATCTTCATC")
  if(spid == "SI-TT-A5") spseqs <- c("GTAGCCCTGT", "GAGCATCTAT")
  if(spid == "SI-TT-A6") spseqs <- c("TAACGCGTGA", "CCCTAACTTC")
  if(spid == "SI-TT-A7") spseqs <- c("TCCCAAGGGT", "TACTACCTTT")
  if(spid == "SI-TT-A8") spseqs <- c("CGAAGTATAC", "GAACTTGGAG")
  if(spid == "SI-TT-A9") spseqs <- c("AAGTGGAGAG", "TTCCTGTTAC")
  if(spid == "SI-TT-A10") spseqs <- c("CGTGACATGC", "ATGGTCTAAA")
  if(spid == "SI-TT-A11") spseqs <- c("CGGAACCCAA", "GATTCGAGGA")
  if(spid == "SI-TT-A12") spseqs <- c("CACCGCACCA", "GACTGTCAAT")
  if(spid == "SI-TT-B1") spseqs <- c("ACAGTAACTA", "ACAGTTCGTT")
  if(spid == "SI-TT-B2") spseqs <- c("TCTACCATTT", "CGGGAGAGTC")
  if(spid == "SI-TT-B3") spseqs <- c("CACGGTGAAT", "GTTCGTCACA")
  if(spid == "SI-TT-B4") spseqs <- c("GTAGACGAAA", "CTAGTGTGGT")

  spseqs
}

.rstest <- function(r, coef) r * (1 + 1/r)^(1 + coef)

.averaging_transform <- function(r, nr){
    d <- c(1, r[ 2:length(r) ] - r[1:(length(r)-1)])
    dr <- c( 0.5 * (d[2:length(d)] + d[1:(length(d)-1)]),d[ length(d) ])
       nr/dr
}


.do.EXACT.DEBUG <- FALSE 
#.do.EXACT.DEBUG <- TRUE 
.simple.Good.Turing.Freq<-function( raw.freq ){
  n_zero <- sum(raw.freq == 0)
  times <- table(raw.freq[raw.freq>0])
  obs <- sort(as.numeric(names(times)))
  times <- times[as.character(obs)]

  nobs <- length(obs)
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
      min.no.turing <- min(which(!xy.sim))
      r <- GT
      if(min.no.turing>0)r[min.no.turing:nobs] <- LGT[min.no.turing:nobs]
      p_zero <- times[1] / total_observed
      rs <- (1-p_zero) * obs * r / sum(r*times*obs/total_observed)
      total_rs <- sum(rs * times)
      puniq <- (1-p_zero) * rs / total_rs
      pres <- rep(p_zero / n_zero, length(raw.freq))
      pres [raw.freq>0] <- puniq[match( raw.freq[raw.freq>0], obs )]
      ret <- list(p0=p_zero, p=pres)

      if(.do.EXACT.DEBUG){
        print("DEBUG_EXACT")
        skf <- "/tmp/del4-YangLiao-GoodTuring-Pr-vals.txt"
        if(file.exists(skf))file.remove(skf)
        sink(skf)
        for(rr in pres){
          cat(sprintf("%.20g\n",rr))
        }
        sink()
      }

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
    sgtr <- .mySGTsorted(as.numeric(names(obs.tab)), obs.tab)
    res <- obs.per.spe
    res[obs.per.spe!=0] <- sgtr$p[match(as.character(obs.per.spe[obs.per.spe!=0]), names(obs.tab))]
    n_zero <- sum(obs.per.spe==0)
    res[obs.per.spe==0] <- sgtr$p0/n_zero
    res
}

.mySGTsorted<-function( obs, times ){
    if(any(obs != sort(obs)))stop("Observations have to be sorted.")
    nobs <- length(obs)
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


.cellCounts_try_cellbarcode <- function( input.directory, sample.sheet, cell.barcode.list, nreads.testing, input.mode ){ # the three parameters can only be one string, not strings!
    if(is.null(sample.sheet)) sample.sheet<-"."
    cmd <- paste0(c(input.directory, sample.sheet, cell.barcode.list, as.character(nreads.testing), input.mode), collapse=.R_param_splitor)
    rvs <- as.integer(rep(0,5))
    C_args <- .C("R_try_cell_barcode_wrapper",nargs=as.integer(5),argv=as.character(cmd),retv=rvs ,PACKAGE="Rsubread")
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

.find_best_cellbarcode <- function( input.directory, sample.sheet=NULL, input.mode="bcl"){
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
        barcode_res <- .cellCounts_try_cellbarcode(input.directory[1], sample.sheet[1], listfile, 30000, input.mode)
        if(length(barcode_res)<3)stop("ERROR: the input sample cannot be processed.")
        sample.good.rate <- barcode_res[2]/barcode_res[1]
        cell.good.rate <- barcode_res[3]/barcode_res[1]
        max.cell.good <- max(max.cell.good, cell.good.rate)
        if(input.mode == "fastq" || input.mode == "bam"){
          cat(sprintf("Cell supporting rate : %.1f%%.\n", cell.good.rate*100.))
        }else{
          cat(sprintf("Sample supporting rate : %.1f%% ; cell supporting rate : %.1f%%.\n", sample.good.rate*100., cell.good.rate*100.))
        }
        if(input.mode=="bcl" && sample.good.rate < 0.5)cat(sprintf("WARNING: there are only %.1f%% reads having known sample indices. Please check if the sample sheet is correct.\n", sample.good.rate*100.))
        if(cell.good.rate > 0.6){
            cat(sprintf("Found cell-barcode list '%s' for the input data: supported by %.1f%% reads.\n", libf, cell.good.rate*100.))
            return(listfile)
        }
    }

    stop(sprintf("ERROR: no known cell barcode set was found for the data set. The highest percentage of cell-barcode matched reads is %.1f%%\n", max.cell.good*100.))
}

library(Matrix)
.read.sparse.mat <- function (fn){
  cat("Loading matrix from",fn,"\n")
  mtx <- readMM(paste0(fn, ".spmtx"))
  coln <- read.delim(paste0(fn, ".BCtab"), stringsAsFactors=F, header=F)$V1
  rown <- read.delim(paste0(fn, ".GENEtab"), stringsAsFactors=F, header=F)$V1
  colnames(mtx) <- coln
  rownames(mtx) <- rown

  mtx
}

.use.DEBUG.data <- NA
# .use.DEBUG.i <- 1

.get.multi.nomial <- function(n, size, prob){
  if(.do.EXACT.DEBUG){
    .use.DEBUG.i <- 1
    ret <- read.table(paste0("/tmp/del4-YangLiao-multi-nomial-SIZE",size,"-",.use.DEBUG.i,".txt"), header=F)
    ret <- t(ret)
   #.use.DEBUG.i <<-  .use.DEBUG.i+1
   #print(ret[1:5,1:5])
    ret
  }else{
    rmultinom(n, size, prob )
  }
}

.simu.multinomial <- function(candi.mat, gene.profile.freq , times=10000){
  bcsizes <- sort(unique( colSums(candi.mat) ))
  ret.nUMI.LLH.tab <- list()
  N10000.LLH <- NA
  Old_bs_one <- NA
  log_E_GTE <- log(gene.profile.freq)

  if(.do.EXACT.DEBUG){
    tff<-"/tmp/del4-YangLiao-for-multi-nomial-steps.txt"
    if(file.exists(tff))file.remove(tff)
    sink(tff)
    for(bcsize_one in bcsizes){
      UMI_step_diff <- NA
      if(!any(is.na(Old_bs_one))) UMI_step_diff <- bcsize_one - Old_bs_one
      st1 <- NA
      if(is.na(Old_bs_one)) {
        st1 <- bcsize_one
      }else if(UMI_step_diff >= 1000){
        st1 <- UMI_step_diff
      }
    
      if(!is.na(st1))cat(st1,"\n")
      Old_bs_one <- bcsize_one
    }
    sink()
    system("python /usr/local/work/liao/subread/scripts/Cellranger-replicate/CrepPY-multi-nomial.py")
  }

  Old_bs_one <- NA
  N1000.FstLLH <- c()

  if(.do.EXACT.DEBUG) M50.gene.ids <- read.table("/tmp/del4-YangLiao-50M-gene-nos.txt", header=F)$V1
  all.steps.len <- 0

  for(bcsize_one in bcsizes){
    UMI_step_diff <- NA
    if(!is.na(Old_bs_one)) UMI_step_diff <- bcsize_one - Old_bs_one

    if(is.na(Old_bs_one)){
    }else if(UMI_step_diff >= 1000){
    }else{
      step_len <- bcsize_one - Old_bs_one
      all.steps.len <- step_len + all.steps.len
    }
    Old_bs_one <- bcsize_one
  }
  if(5000 <= all.steps.len)print("WARNING: SMALL STEPS >= 5000!")

  Old_bs_one <- NA
  M50.in.loop.K <- 0
  for(bcsize_one in bcsizes){
    UMI_step_diff <- NA
    if(!any(is.na(Old_bs_one))) UMI_step_diff <- bcsize_one - Old_bs_one
  
    if(any(is.na(N10000.LLH))){
      N10000 <- .get.multi.nomial(n=times, size=bcsize_one, prob=gene.profile.freq )
      print("======== N10000 and PROF DIM =======")
      print(dim(N10000))
      print(length(gene.profile.freq))
      N10000.LLH <- apply(N10000, 2, function(x) dmultinom(x, prob=gene.profile.freq, log =T ))
    }else if(UMI_step_diff >= 1000){
      UMI_step_diff_N10000 <- .get.multi.nomial(n=times, size=UMI_step_diff, prob=gene.profile.freq )
      N10000 <- N10000 + UMI_step_diff_N10000
      N10000.LLH <- apply(N10000, 2, function(x) dmultinom(x, prob=gene.profile.freq, log =T ))
    }else{
      step_len <- bcsize_one - Old_bs_one
      for(curi in (Old_bs_one+1):bcsize_one){
        if(.do.EXACT.DEBUG){
          UMI_step_indices_N10000 <- M50.gene.ids[1+M50.in.loop.K + (0:9999) * all.steps.len ]
          if(curi == 588){
            tff<-"/tmp/del4-YangLiao-DEBUG-588-10KJ.txt"
            if(file.exists(tff))file.remove(tff)
            sink(tff)
            for(pi in 1:10000){
                cat(pi," ", N10000[UMI_step_indices_N10000[pi],pi],"\n")
            }
            sink()
          }
          M50.in.loop.K <- 1+M50.in.loop.K
        }else{
          UMI_step_indices_N10000 <- sample(1:nrow(candi.mat), size=times, prob=gene.profile.freq, replace=T)
        }
        for(idi in 1:times){
          N10000[ UMI_step_indices_N10000[idi] ,idi ] <- N10000[ UMI_step_indices_N10000[idi] ,idi ]+1
        }
        UMI_step_values_N10000 <- rep(0, times)
        for(idi in 1:times){
          UMI_step_values_N10000[idi] <- N10000[UMI_step_indices_N10000[idi] ,idi ]
        }
        #print(UMI_step_values_N10000)
        N10000.LLH <- N10000.LLH + log_E_GTE[ UMI_step_indices_N10000 ] + log((curi) / UMI_step_values_N10000)

        if(.do.EXACT.DEBUG && curi == 588)for(uii in 1:30){
            tff<-"/tmp/del4-YangLiao-DEBUG-1-30-E_GTE.txt"
            if(file.exists(tff))file.remove(tff)
            sink(tff)
            for(i in 1:30){
                cat(sprintf("%.5g %.5g %.5g\n", N10000.LLH[i], log_E_GTE[ UMI_step_indices_N10000 [i]],  log((curi) / UMI_step_values_N10000[i])))
            }
            sink()
        }
      }
    }
  
    ret.nUMI.LLH.tab[[ bcsize_one ]] <- N10000.LLH
    N1000.FstLLH <- c(  N1000.FstLLH, N10000.LLH[1] )
    Old_bs_one <- bcsize_one
  }
  if(.do.EXACT.DEBUG){
    tff<-"/tmp/del4-YangLiao-10K-LLH-from-R.txt"
    if(file.exists(tff))file.remove(tff)
    sink(tff)
    for(flh in N1000.FstLLH){
      cat(sprintf("%.6g\n", flh))
    }
    sink()
  }
  ret.nUMI.LLH.tab
}

.is.no.candBC <- function(fn){
  info <- file.info(fn)
  return(info$size < 2)
}

.cellCounts.rescue <- function( BAM.name, FC.gene.ids, sample.no ){
  fname <- sprintf("%s.scRNA.%03d", BAM.name, sample.no)
  nozero.anywhere.genes <- read.delim(paste0(fname,".no0Genes"), stringsAsFactors=F, header=F)$V1
  ambient.accumulate <- read.delim(paste0(fname,".AmbSum"), stringsAsFactors=F)
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
      actual.pvalues[candi] <- (1+sum(cand.simu.pvs < cand.actual.pv))/(length(cand.simu.pvs)+1)
    }
  
    # p-value => FDR
    actual.FDR <- p.adjust(actual.pvalues, method='BH')
  
    if(.do.EXACT.DEBUG){
      tff<-"/tmp/del4-YangLiao-final-FDR.txt"
      if(file.exists(tff))file.remove(tff)
      sink(tff)
      for(candi in names(log.like.cands)) cat(sprintf("%s %.10g %.10g\n", candi, actual.pvalues[candi], actual.FDR[candi]))
      sink()
    }
    # select cells that has FDR < cutoff
    FDR.Cutoff <- 0.01
    head(actual.FDR)
    Rescured.Barcodes <- names(actual.FDR)[actual.FDR <= FDR.Cutoff]
    head(Rescured.Barcodes)
    rescue.candidates[,Rescured.Barcodes, drop=FALSE]
  }
}

.match.exons <- function(annot.tab, counts){
  annot.str <- paste(annot.tab$GeneID, annot.tab$Chr, annot.tab$Start,
    annot.tab$End, ifelse(annot.tab$Strand == "+", "P","N"), sep=":fc@R@Spl:")
  ret <- matrix(0, ncol=ncol(counts), nrow=nrow(annot.tab))
  colnames(ret) <- colnames(counts)
  rownames(ret) <- annot.str
  ret[ rownames(counts),] <- counts
  ret
}

.load.one.scSample <- function( BAM.name, FC.gene.ids, sample.no, use.meta.features, annot.tab, umi.cutoff){
  set.seed(0)
  fname <- sprintf("%s.scRNA.%03d", BAM.name, sample.no)
  cat("Loading high-conf matrix from '",fname,"'\n")
  highconf <- as.matrix(.read.sparse.mat(paste0(fname,".HighConf")))
  rescued <- NA
  if(is.null(umi.cutoff)) rescued <- .cellCounts.rescue(BAM.name, FC.gene.ids, sample.no)

  if(use.meta.features){
    if(!any(is.na(rescued))) rescued <- rescued[rowSums(rescued)>0,]
    ncolRescued <- 0
    if(any(is.na(rescued))){
      ret <- matrix(0,ncol=ncol(highconf), nrow=length(FC.gene.ids))
      colnames(ret) <- colnames(highconf)
    }else{
      ret <- matrix(0,ncol=ncol(highconf)+ncol(rescued), nrow=length(FC.gene.ids))
      colnames(ret) <- c( colnames(highconf), colnames(rescued) )
    }
    rownames(ret) <- FC.gene.ids 
    ret[rownames(highconf), colnames(highconf) ] <- highconf
    if(!any(is.na(rescued)))ret[rownames(rescued), colnames(rescued) ] <- rescued

    list(Counts=ret, HighConfidneceCell=colnames(ret) %in% colnames(highconf))
  }else{
    fcmat <- .match.exons(annot.tab,highconf)
    rownames(fcmat)<-NULL
    list(Counts=fcmat, HighConfidneceCell=rep(T, ncol(highconf))) 
  }
}

.load.all.scSamples <- function( BAM.name, FC.gene.ids, use.meta.features, annot.tab, umi.cutoff){
  sum.tab <- read.delim(paste0(BAM.name,".scRNA.SampleTable"), stringsAsFactors=F)
  ret <- list()
  for(roiw in 1:nrow(sum.tab)){
    sname <- as.character(sum.tab$SampleName[roiw])
    sid <- sum.tab$Index[roiw]
    count.tab <- .load.one.scSample(BAM.name, FC.gene.ids, sid, use.meta.features, annot.tab, umi.cutoff)
    ret[[sprintf("Sample.%d",sid)]] <- count.tab
  }
  ret[["Sample.Table"]] <- sum.tab

  ret
}

.del.temp.files  <- function(pfx){
  if(nchar(pfx)<16) stop("Prefix should be longer.")
  flist <- dir(".", "^[.]Rsubread", all.files=TRUE)
  to.delete <- flist[grepl(pfx, flist)]
  cat(paste(flist,";"),"Delete temporary files :", paste(to.delete, collapse=","),"\nthat have a prefix of",pfx,"\n")
  for(df in to.delete){
    if(grepl(pfx, df)){
      file.remove(df)
#cat ("Goinng to remove",df,"\n")
    }else{
      stop("Wrong selection of files to delete")
    }
  }
}

.extract.sample.table.cols <- function(rdir, smr, input.mode="bcl", umi.cutoff=NULL){
  total.cells <- c()
  hiconf.cells <- c()
  res.cells <- c()
  umis <- c()
  umi.statistics <- data.frame(MinUMI=NULL, MedianUMI=NULL, MaxUMI=NULL, MeanUMI=NULL,stringsAsFactors=F)[,0]
  for(spi in 1:nrow(smr[["Sample.Table"]])){
    samplename <- as.character(smr[["Sample.Table"]]$SampleName[spi])
    sampleno <- sprintf("Sample.%d", spi)
    total.cells <- c(total.cells,length(smr[[sampleno]][["HighConfidneceCell"]]))
    hiconf.cells <- c(hiconf.cells,sum(smr[[sampleno]][["HighConfidneceCell"]]))
    res.cells <- c(res.cells, sum(!(smr[[sampleno]][["HighConfidneceCell"]])))
    cell.umis <- colSums(smr[[sampleno]][["Counts"]])
    umis <- c(umis,sum(cell.umis))
    umi.statistics <- rbind(umi.statistics, list(MinUMI=min(cell.umis), MedianUMI=median(cell.umis), MaxUMI=max(cell.umis), MeanUMI=mean(cell.umis)))
  }
  ret <- NULL
  if(is.null(umi.cutoff)){
    if(input.mode=="fastq" || input.mode=="bam")
      ret <- cbind( SampleName=as.character(smr[["Sample.Table"]]$SampleName), TotalCells=total.cells, HighConfidenceCells=hiconf.cells, RescuedCells=res.cells, TotalUMI=umis, umi.statistics, smr[["Sample.Table"]][,c("TotalReads","MappedReads","AssignedReads")],stringsAsFactors=F )
    else
      ret <- cbind( SampleName=as.character(smr[["Sample.Table"]]$SampleName), InputDirectory=rep(rdir, nrow(smr[["Sample.Table"]])), TotalCells=total.cells, HighConfidenceCells=hiconf.cells, RescuedCells=res.cells, TotalUMI=umis, umi.statistics, smr[["Sample.Table"]][,c("TotalReads","MappedReads","AssignedReads")],stringsAsFactors=F )
  }else{
    if(input.mode=="fastq" || input.mode=="bam")
      ret <- cbind( SampleName=as.character(smr[["Sample.Table"]]$SampleName), TotalCells=total.cells, TotalUMI=umis, umi.statistics, smr[["Sample.Table"]][,c("TotalReads","MappedReads","AssignedReads")], stringsAsFactors=F)
    else
      ret <- cbind( SampleName=as.character(smr[["Sample.Table"]]$SampleName), InputDirectory=rep(rdir, nrow(smr[["Sample.Table"]])), TotalCells=total.cells, TotalUMI=umis, umi.statistics, smr[["Sample.Table"]][,c("TotalReads","MappedReads","AssignedReads")], stringsAsFactors=F)
  }
}

.scan.fastq.dir <- function(dirname){
  ddf <- list.files(dirname , recursive=TRUE)
  ddf <- ddf[grepl("L[0-9][0-9][1-9]_R1_[0-9][0-9][1-9].fastq.gz$",ddf) &(!grepl("Undetermined_S",ddf)) &!grepl("/fork0/", ddf)]
  sample.sheet <- data.frame(stringsAsFactors=F)
  for(df1 in sort(ddf)){
    df1 <- file.path(dirname,df1)
    r1_part <- substr(df1, nchar(df1)-12, 9999)
    r1_pref <- substr(df1, 1, nchar(df1)-15)
    df2 <- paste0(r1_pref, "R2", r1_part)
    sample.name <- unlist(strsplit(df2, "\\/"))
    sample.name0<- sample.name
    sample.name <- sample.name[length(sample.name)]
    if(substr(sample.name, nchar(sample.name)-22, nchar(sample.name)-22)=='S') sample.name <- substr(sample.name, 1, nchar(sample.name)-24)
    else if(substr(sample.name, nchar(sample.name)-23, nchar(sample.name)-23)=='S') sample.name <- substr(sample.name, 1, nchar(sample.name)-25)
    else stop(paste("ERROR: Unable to parse the file name,",sample.name))

    if(substr(sample.name,1,1)=="_") sample.name <- substr(sample.name,2,999)
    if(substr(sample.name,1,1)=="_") sample.name <- substr(sample.name,2,999)
    #if(substr(sample.name,1,1)!="S")stop(sprintf("Unable to parse the FASTQ file name : '%s'", sample.name0))

    sample.sheet <- rbind(sample.sheet, data.frame( SampleName=as.character(sample.name), BarcodeUMIFile=df1, ReadFile=df2, stringsAsFactors=F))
  }
  sample.sheet
}

.SCRNA_FASTA_SPLIT1 <- "|Rsd:cCounts:mFQs|"
.SCRNA_FASTA_SPLIT2 <- "|Rsd:cCounts:1mFQ|"

cellCounts <- function(index, sample,input.mode="BCL", cell.barcode=NULL, aligner="align", annot.inbuilt="mm10",annot.ext=NULL,isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",useMetaFeatures=TRUE, umi.cutoff=NULL, nthreads=10, nBestLocations =1, unique.mapping=FALSE, ...){
  .remove.temp.files <- T
  if(!is.null(aligner)) aligner <- match.arg(aligner,c("subjunc","align")) 
  if(!is.null(umi.cutoff)){
    umi.cutoff <- as.numeric(umi.cutoff)
    if(umi.cutoff < 0.0) stop("UMI cutoff must be a positive number.")
  }
  BAM_is_Rerun_Persample <- is.null(aligner)
  sample.info.idx <- sample
  fc <- list()

  temp.file.prefix <- file.path(".",paste(".Rsubread_cCounts_Tmp_for_Pid_",Sys.getpid(),"_Rproc",sep=""))
  raw.fc.annot <- NA
  df.sample.info <- data.frame(stringsAsFactors=F)
  sample.1 <- paste0(temp.file.prefix,".samplesheet")

  if(input.mode=="BCL"){
    sample.info.idx$SampleName <- as.character(sample.info.idx$SampleName)
    if('IndexSetName' %in% colnames(sample.info.idx))sample.info.idx$IndexSetName <- as.character(sample.info.idx$IndexSetName)

    dirs <- unique( sample.info.idx$InputDirectory )
    fc[["counts"]] <- list()
#    if(is.null(umi.cutoff))fc[["cell.confidence"]] <- list()
  
    for(dirname in dirs) .check_and_NormPath(dirname, mustWork=TRUE, "InputDirectory in sample.info.idx") # check files before the slow mapping/counting step.

    for(dirname in dirs){
      unique.samples <- unique( as.character(sample.info.idx$SampleName[ sample.info.idx$InputDirectory == dirname ] ))
  
      if(is.null(cell.barcode)){
        .index.names.to.sheet.raw.dir.mode(dirname, sample.info.idx, sample.1)
        cell.barcode <- .find_best_cellbarcode(dirname, sample.1)
      }else{
        cell.barcode <- .check_and_NormPath(cell.barcode, mustWork=T, opt="cell.barcode")
      }
  
      full_dirname <- .check_and_NormPath(dirname, mustWork=TRUE, "InputDirectory in sample.info.idx")
      if(is.null(aligner)){
        generate.scRNA.BAM <- FALSE 
        if(!all(file.exists(paste0( unique.samples,".bam" )))) stop("No aligner is specified but the BAM file does not exist. Please specify 'align' or 'subjunc' as the aligner.")
        for(samplename in unique.samples){
          .index.names.to.sheet.raw.dir.mode(dirname, sample.info.idx, sample.1, samplename)
          one.bam.name <- paste0(samplename, ".bam")
          .write.tmp.parameters(list(sampleSheet=sample.1, umi.cutoff=umi.cutoff, cellBarcodeList=cell.barcode, generate.scRNA.BAM=generate.scRNA.BAM, BAM_is_Rerun_Persample=BAM_is_Rerun_Persample))
          one.raw.fc <- featureCounts(one.bam.name, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTFAnnotationFile=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, useMetaFeatures=useMetaFeatures, nthreads=nthreads, strandSpecific=1, ...)
          if(any(is.na(raw.fc.annot))) raw.fc.annot<- one.raw.fc$annotation
          cat("Processing sample '",samplename,"'\n")
          one.result <- .load.all.scSamples(paste0(samplename,".bam"), as.character(raw.fc.annot$GeneID), useMetaFeatures, raw.fc.annot, umi.cutoff)
          fc[["counts"]][[samplename]] <- one.result[["Sample.1"]][["Counts"]] # only one sample.
          if(is.null(umi.cutoff))fc[["cell.confidence"]][[samplename]] <- one.result[["Sample.1"]][["HighConfidneceCell"]]
          stt <- .extract.sample.table.cols(full_dirname,one.result, umi.cutoff=umi.cutoff)
          df.sample.info <- rbind(df.sample.info, stt)
        }
      }else{
        .index.names.to.sheet.raw.dir.mode(dirname, sample.info.idx, sample.1)
        generate.scRNA.BAM <- TRUE
        .write.tmp.parameters(list(isBCLinput=TRUE))
        if(aligner=="align"){
          align(index, full_dirname, output_file=temp.file.prefix, nthreads=nthreads, useAnnotation=TRUE, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTF=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, ...)
        }else if(aligner=="subjunc"){
          subjunc(index, full_dirname, output_file=temp.file.prefix, nthreads=nthreads, ...)
        }
        .write.tmp.parameters(list(sampleSheet=sample.1, umi.cutoff=umi.cutoff, cellBarcodeList=cell.barcode, generate.scRNA.BAM=generate.scRNA.BAM,BAM_is_Rerun_Persample=BAM_is_Rerun_Persample))
        raw.fc<-featureCounts(temp.file.prefix, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTFAnnotationFile=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, useMetaFeatures=useMetaFeatures,nthreads=nthreads,strandSpecific=1, ...)
        if(any(is.na(raw.fc.annot))) raw.fc.annot<-raw.fc$annotation
        some.results <- .load.all.scSamples(temp.file.prefix, as.character(raw.fc.annot$GeneID), useMetaFeatures, raw.fc.annot, umi.cutoff)
        for(spi in 1:nrow(some.results[["Sample.Table"]])){
          samplename <- as.character(some.results[["Sample.Table"]][["SampleName"]][spi])
          fc[["counts"]][[samplename]] <- some.results[[sprintf("Sample.%d", spi)]][["Counts"]] # only one sample.
          if(is.null(umi.cutoff))fc[["cell.confidence"]][[samplename]] <- some.results[[sprintf("Sample.%d", spi)]][["HighConfidneceCell"]]
        }
        stt <- .extract.sample.table.cols(full_dirname,some.results, umi.cutoff=umi.cutoff)
        df.sample.info <- rbind(df.sample.info, stt)
      }
    }
  } else if(input.mode == "BAM"){
    if(is.null(cell.barcode)){
      cell.barcode <- .find_best_cellbarcode(sample$BAMFile, "N/A", input.mode="bam")
    }else{
      cell.barcode <- .check_and_NormPath(cell.barcode, mustWork=T, opt="cell.barcode")
    }

    unique.samples <- unique(as.character(sample.info.idx$SampleName))
    for(this.sample in unique.samples){
      .index.names.to.sheet.BAM.mode(data.frame(BAMFile="COMBINED.INPUT",SampleName=as.character(this.sample),stringsAsFactors=F), sample.1)
      BAM.names <- paste(sample.info.idx$BAMFile[ as.character(sample.info.idx$SampleName) == this.sample ], collapse=.SCRNA_FASTA_SPLIT1)
      .write.tmp.parameters(list(isScRNABAMinput=TRUE))
      align(index, BAM.names, output_file=temp.file.prefix, nthreads=nthreads, useAnnotation =TRUE,annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTF=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType,...)

      generate.scRNA.BAM <- TRUE
      .write.tmp.parameters(list(BAM_is_ScRNA_BAM=TRUE, sampleSheet=sample.1, umi.cutoff=umi.cutoff, cellBarcodeList=cell.barcode, generate.scRNA.BAM=generate.scRNA.BAM,BAM_is_Rerun_Persample=BAM_is_Rerun_Persample))
      raw.fc<-featureCounts(temp.file.prefix, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTFAnnotationFile=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, useMetaFeatures=useMetaFeatures,nthreads=nthreads,strandSpecific=1, ...)

      if(any(is.na(raw.fc.annot))) raw.fc.annot<-raw.fc$annotation
      some.results <- .load.all.scSamples(temp.file.prefix, as.character(raw.fc.annot$GeneID), useMetaFeatures, raw.fc.annot, umi.cutoff)
      for(spi in 1:nrow(some.results[["Sample.Table"]])){
        samplename <- as.character(some.results[["Sample.Table"]][["SampleName"]][spi])
        fc[["counts"]][[samplename]] <- some.results[[sprintf("Sample.%d", spi)]][["Counts"]] # only one sample.
        if(is.null(umi.cutoff))fc[["cell.confidence"]][[samplename]] <- some.results[[sprintf("Sample.%d", spi)]][["HighConfidneceCell"]]
      }
    }
    df.sample.info <- .extract.sample.table.cols(NA,some.results, input.mode="bam", umi.cutoff=umi.cutoff)
  } else if(input.mode == "FASTQ" || input.mode == "FASTQ-dir"){
    if(input.mode == "FASTQ-dir") sample.info.idx <- .scan.fastq.dir(sample)
    print(sample.info.idx)
    if(!("BarcodeUMIFile" %in%  colnames(sample.info.idx) && "ReadFile" %in%  colnames(sample.info.idx) )) stop("You need to provide BC+UMI and Genomic sequence files")
 
    combined.fastq.names <- ""
    for(rowi in 1:nrow(sample.info.idx)){
      df1 <- as.character(sample.info.idx[rowi,"BarcodeUMIFile"])
      df2 <- as.character(sample.info.idx[rowi,"ReadFile"])
      R1.file.name <- .check_and_NormPath(df1, mustWork=TRUE, "The barcode FASTQ file in sample.info.idx")
      R2.file.name <- .check_and_NormPath(df2, mustWork=TRUE, "The genomic read FASTQ file in sample.info.idx")
      combined.fastq.names <- paste0( combined.fastq.names,.SCRNA_FASTA_SPLIT1, R1.file.name, .SCRNA_FASTA_SPLIT2,".",.SCRNA_FASTA_SPLIT2, R2.file.name)
    }
    combined.fastq.names <- substr(combined.fastq.names, nchar(.SCRNA_FASTA_SPLIT1)+1, 9999999)
    if(is.null(cell.barcode)){
      cell.barcode <- .find_best_cellbarcode(combined.fastq.names, "N/A", input.mode="fastq")
    }else{
      cell.barcode <- .check_and_NormPath(cell.barcode, mustWork=T, opt="cell.barcode")
    }
 
    bam.for.FC <- c()
    if(is.null(aligner)){
      unique.samples <- unique(as.character(sample.info.idx$SampleName))
      bam.for.FC <- paste0(sample.info.idx$SampleName,".bam")
      generate.scRNA.BAM <- FALSE
      
      for(samplename in unique.samples){
        leftover.bam <- paste0(samplename, ".bam")

        .index.names.to.sheet.FASTQ.mode(sample.info.idx[ as.character(sample.info.idx$SampleName) == samplename, ][1,], sample.1)
        .write.tmp.parameters(list(BAM_is_ScRNA_Fastq=TRUE, sampleSheet=sample.1, umi.cutoff=umi.cutoff, cellBarcodeList=cell.barcode, generate.scRNA.BAM=generate.scRNA.BAM,BAM_is_Rerun_Persample=BAM_is_Rerun_Persample))
        raw.fc<-featureCounts(leftover.bam, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTFAnnotationFile=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, useMetaFeatures=useMetaFeatures,nthreads=nthreads,strandSpecific=1, ...)

        if(any(is.na(raw.fc.annot))) raw.fc.annot<-raw.fc$annotation
        some.results <- .load.all.scSamples(leftover.bam, as.character(raw.fc.annot$GeneID), useMetaFeatures, raw.fc.annot, umi.cutoff)
        fc[["counts"]][[samplename]] <- some.results[[sprintf("Sample.%d", 1)]][["Counts"]] # featureCounts treats each left-over BAM file as a single run and has only one sample 
        if(is.null(umi.cutoff))fc[["cell.confidence"]][[samplename]] <- some.results[[sprintf("Sample.%d", 1)]][["HighConfidneceCell"]]

        newrow<- .extract.sample.table.cols(NA,some.results,input.mode="fastq", umi.cutoff=umi.cutoff)
        df.sample.info <- rbind(df.sample.info,newrow)
      }
    }else{
      .write.tmp.parameters(list(isScRNAFastqinput=TRUE))
      align(index, combined.fastq.names, output_file=temp.file.prefix, nthreads=nthreads, useAnnotation =TRUE, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTF=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, nBestLocations = nBestLocations, unique=unique.mapping)
      bam.for.FC <- temp.file.prefix
      generate.scRNA.BAM <- TRUE
      .index.names.to.sheet.FASTQ.mode(sample.info.idx, sample.1)
      .write.tmp.parameters(list(BAM_is_ScRNA_Fastq=TRUE, sampleSheet=sample.1, umi.cutoff=umi.cutoff, cellBarcodeList=cell.barcode, generate.scRNA.BAM=generate.scRNA.BAM,BAM_is_Rerun_Persample=BAM_is_Rerun_Persample))
      raw.fc<-featureCounts(bam.for.FC, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTFAnnotationFile=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, useMetaFeatures=useMetaFeatures,nthreads=nthreads,strandSpecific=1, ...)

     if(any(is.na(raw.fc.annot))) raw.fc.annot<-raw.fc$annotation
      some.results <- .load.all.scSamples(temp.file.prefix, as.character(raw.fc.annot$GeneID), useMetaFeatures, raw.fc.annot, umi.cutoff)

      for(spi in 1:nrow(some.results[["Sample.Table"]])){
        samplename <- as.character(some.results[["Sample.Table"]][["SampleName"]][spi])
        fc[["counts"]][[samplename]] <- some.results[[sprintf("Sample.%d", spi)]][["Counts"]] # only one sample.
        if(is.null(umi.cutoff))fc[["cell.confidence"]][[samplename]] <- some.results[[sprintf("Sample.%d", spi)]][["HighConfidneceCell"]]
      }
      df.sample.info <- .extract.sample.table.cols(NA,some.results,input.mode="fastq", umi.cutoff=umi.cutoff)
    }
  }

  fc[["annotation"]] <- raw.fc.annot
  fc[["sample.info"]] <- df.sample.info
  if(.remove.temp.files).del.temp.files(substr(temp.file.prefix,4,99))
  fc
}
