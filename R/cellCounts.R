.convert.sample_index.id.to.seq <- function(spid){
  spseqs = NA
  if(spid == "SI_001") spseqs = c('TCGCCATA', 'GTATACAC', 'AATGGTGG', 'CGCATGCT')
  if(spid == "SI_002") spseqs = c('TATCCTCG', 'GCGAGGTC', 'CGCTTCAA', 'ATAGAAGT')
  if(spid == "SI_003") spseqs = c('TGACGTCG', 'CTTGTGTA', 'ACGACCGT', 'GACTAAAC')
  if(spid == "SI_004") spseqs = c('ATCTAGCT', 'GAGCGTAC', 'TCAGCCTG', 'CGTATAGA')
  if(spid == "SI_005") spseqs = c('CCGTTCCC', 'ATACAGTT', 'TGTAGTAA', 'GACGCAGG')
  if(spid == "SI_006") spseqs = c('TCAATTGG', 'AGTTAGAA', 'GAGCGCTT', 'CTCGCACC')
  if(spid == "SI_007") spseqs = c('CTGCCTTG', 'ACTAGCCC', 'GGCTAGAT', 'TAAGTAGA')
  if(spid == "SI_008") spseqs = c('GGCAGAAA', 'ACGGTTCT', 'CATTCGTC', 'TTACACGG')
  if(spid == "SI_P01_A1") spseqs = c('TTGTAAGA', 'GGCGTTTC', 'CCTACCAT', 'AAACGGCG')
  if(spid == "SI_P01_B1") spseqs = c('CTAGCTGT', 'GCCAACAA', 'AGGCTACC', 'TATTGGTG')
  if(spid == "SI_P01_C1") spseqs = c('GATGCAGT', 'AGACTTTC', 'TTCTACAG', 'CCGAGGCA')
  if(spid == "SI_P01_D1") spseqs = c('AGCTGCGT', 'GTGGAGCA', 'TCTATTAG', 'CAACCATC')
  if(spid == "SI_P01_E1") spseqs = c('CGCCCGTA', 'GTTTGCCT', 'TAAGTTAG', 'ACGAAAGC')
  if(spid == "SI_P01_F1") spseqs = c('TGACTAGT', 'GATAGGTA', 'CCGTACAG', 'ATCGCTCC')
  if(spid == "SI_P01_G1") spseqs = c('CATATGCG', 'ATGCGATT', 'TCCGCTAC', 'GGATACGA')
  if(spid == "SI_P01_H1") spseqs = c('TCTTGTCC', 'CGGGAGTA', 'GTCACAGG', 'AAACTCAT')
  if(spid == "SI_P01_A2") spseqs = c('AGCCCTTT', 'TCTTAGGC', 'GTGAGAAG', 'CAAGTCCA')
  if(spid == "SI_P01_B2") spseqs = c('GGTCGAGC', 'TTAGATTG', 'CCCACCCA', 'AAGTTGAT')
  if(spid == "SI_P01_C2") spseqs = c('CCGAGAAC', 'TGCTCTGT', 'GTAGTGCG', 'AATCACTA')
  if(spid == "SI_P01_D2") spseqs = c('ACATTCCG', 'GACACAAT', 'CTGCGGTA', 'TGTGATGC')
  if(spid == "SI_P01_E2") spseqs = c('TCATCAAG', 'GTTGGTCC', 'AGGCTGGT', 'CACAACTA')
  if(spid == "SI_P01_F2") spseqs = c('TGCCCGCT', 'GCAAACGC', 'CATTGATA', 'ATGGTTAG')
  if(spid == "SI_P01_G2") spseqs = c('CGGTGAGC', 'ATAACCTA', 'TCCGAGCG', 'GATCTTAT')
  if(spid == "SI_P01_H2") spseqs = c('CCGAACTC', 'AACGGTCA', 'TTATTGGT', 'GGTCCAAG')
  if(spid == "SI_P01_A3") spseqs = c('AAAGCATA', 'GCCTTTAT', 'CTGCAGCC', 'TGTAGCGG')
  if(spid == "SI_P01_B3") spseqs = c('TCATCCTT', 'ATTGGACG', 'CAGCTTAC', 'GGCAAGGA')
  if(spid == "SI_P01_C3") spseqs = c('ACGTTACA', 'TTACCTAC', 'GACGACGG', 'CGTAGGTT')
  if(spid == "SI_P01_D3") spseqs = c('GAGCACGC', 'CGAAGTTG', 'TTCGTGAA', 'ACTTCACT')
  if(spid == "SI_P01_E3") spseqs = c('TCTGCAGG', 'CGGCTCCA', 'AACAAGTC', 'GTATGTAT')
  if(spid == "SI_P01_F3") spseqs = c('TTAGGACC', 'AGTCTGTA', 'GCCTCCGT', 'CAGAATAG')
  if(spid == "SI_P01_G3") spseqs = c('TACGAGTT', 'ATGTCCAG', 'GCTATAGC', 'CGACGTCA')
  if(spid == "SI_P01_H3") spseqs = c('TTGGGCTT', 'GACAAACC', 'ACACCTAA', 'CGTTTGGG')
  if(spid == "SI_P01_A4") spseqs = c('CATGGCAG', 'AGAACGCC', 'GTCTTTGA', 'TCGCAATT')
  if(spid == "SI_P01_B4") spseqs = c('GACAGGCT', 'CCTCTAAC', 'AGGGACTG', 'TTATCTGA')
  if(spid == "SI_P01_C4") spseqs = c('ACATTGGC', 'GAGCCCAT', 'CTTAGTCA', 'TGCGAATG')
  if(spid == "SI_P01_D4") spseqs = c('AAATCGTC', 'GCGATCGG', 'CTTCGAAT', 'TGCGATCA')
  if(spid == "SI_P01_E4") spseqs = c('GTAAACAT', 'TATCCTGA', 'AGCTGACG', 'CCGGTGTC')
  if(spid == "SI_P01_F4") spseqs = c('GCATGATA', 'CGTCCTCT', 'AACGACAC', 'TTGATGGG')
  if(spid == "SI_P01_G4") spseqs = c('CCTGCGGT', 'GTACAACG', 'AGCTTCTC', 'TAGAGTAA')
  if(spid == "SI_P01_H4") spseqs = c('TTCCATCT', 'ACTGGAGC', 'CGGTCGTG', 'GAAATCAA')
  if(spid == "SI_P01_A5") spseqs = c('CAGTCTGG', 'TCACACTC', 'ATTGGGAA', 'GGCATACT')
  if(spid == "SI_P01_B5") spseqs = c('TCATGCGA', 'ATCGTACT', 'CATCAGTG', 'GGGACTAC')
  if(spid == "SI_P01_C5") spseqs = c('TCAGTCAA', 'CACTGACT', 'ATGCATTC', 'GGTACGGG')
  if(spid == "SI_P01_D5") spseqs = c('GTCGACTC', 'AGACGGAT', 'CCTTTAGA', 'TAGACTCG')
  if(spid == "SI_P01_E5") spseqs = c('CCGTTGAA', 'TATGCTCT', 'ATCCAAGG', 'GGAAGCTC')
  if(spid == "SI_P01_F5") spseqs = c('TCTGACTA', 'GTACGGCT', 'CGGTTTAG', 'AACACAGC')
  if(spid == "SI_P01_G5") spseqs = c('ATGAAGTA', 'GAAGCTCG', 'TCTTTCGT', 'CGCCGAAC')
  if(spid == "SI_P01_H5") spseqs = c('ATAGTATG', 'TATAAGGA', 'GGCTCCAC', 'CCGCGTCT')
  if(spid == "SI_P01_A6") spseqs = c('CTTTCGAC', 'ACGGGACT', 'TGCATCTG', 'GAACATGA')
  if(spid == "SI_P01_B6") spseqs = c('GCGCACCT', 'AACGCGAA', 'CTATTTGG', 'TGTAGATC')
  if(spid == "SI_P01_C6") spseqs = c('CGCTCAGG', 'GAGGTTTA', 'ACTCAGAC', 'TTAAGCCT')
  if(spid == "SI_P01_D6") spseqs = c('GAAGTCTT', 'TGCAGGGC', 'ATGCCAAA', 'CCTTATCG')
  if(spid == "SI_P01_E6") spseqs = c('TGATGGCT', 'GCCATTTG', 'ATTGAAAC', 'CAGCCCGA')
  if(spid == "SI_P01_F6") spseqs = c('GACTTCCT', 'TGAGGAAG', 'ATGCCGGC', 'CCTAATTA')
  if(spid == "SI_P01_G6") spseqs = c('GTGCGACA', 'TCAGTGTT', 'AGCACTGG', 'CATTACAC')
  if(spid == "SI_P01_H6") spseqs = c('AAGCATAA', 'CCCATCGC', 'TTAGCGCT', 'GGTTGATG')
  if(spid == "SI_P01_A7") spseqs = c('CTCATCAT', 'TAACGTCC', 'AGGTCATA', 'GCTGAGGG')
  if(spid == "SI_P01_B7") spseqs = c('TCCACACG', 'CTTCTGTT', 'GAATGCAC', 'AGGGATGA')
  if(spid == "SI_P01_C7") spseqs = c('GGCGGAAT', 'ACACCGGG', 'CATAATCC', 'TTGTTCTA')
  if(spid == "SI_P01_D7") spseqs = c('CCGGATCC', 'GGTCGCAT', 'TTAACGTG', 'AACTTAGA')
  if(spid == "SI_P01_E7") spseqs = c('AAGACGTG', 'CCATGTGT', 'GTTCACAA', 'TGCGTACC')
  if(spid == "SI_P01_F7") spseqs = c('GGTTAGAC', 'CAAACTTT', 'ACCCGAGA', 'TTGGTCCG')
  if(spid == "SI_P01_G7") spseqs = c('GCCGGTAA', 'TGACTGCC', 'ATTACCGG', 'CAGTAATT')
  if(spid == "SI_P01_H7") spseqs = c('TGGCACGA', 'AACGGGTG', 'CTAATTCT', 'GCTTCAAC')
  if(spid == "SI_P01_A8") spseqs = c('GACTGTTC', 'ATGATACG', 'CCACAGAA', 'TGTGCCGT')
  if(spid == "SI_P01_B8") spseqs = c('ACGTTCAC', 'TGCCCAGA', 'CAAGGTCT', 'GTTAAGTG')
  if(spid == "SI_P01_C8") spseqs = c('TTTATCCC', 'GCACGTTT', 'CAGGAAGA', 'AGCTCGAG')
  if(spid == "SI_P01_D8") spseqs = c('AATCTTTG', 'GGATGAGT', 'CTCAAGAC', 'TCGGCCCA')
  if(spid == "SI_P01_E8") spseqs = c('GCTTACAT', 'TAGGGTGC', 'AGCCTATG', 'CTAACGCA')
  if(spid == "SI_P01_F8") spseqs = c('AGTTGGGA', 'TACATTCT', 'CCAGAAAG', 'GTGCCCTC')
  if(spid == "SI_P01_G8") spseqs = c('AAGTACTC', 'GGAACTCT', 'TCCCTGAG', 'CTTGGAGA')
  if(spid == "SI_P01_H8") spseqs = c('AAGAGCGG', 'TCATAGCA', 'GGCCCATC', 'CTTGTTAT')
  if(spid == "SI_P01_A9") spseqs = c('GAGTGCGT', 'CTCCAACA', 'ACAACTTG', 'TGTGTGAC')
  if(spid == "SI_P01_B9") spseqs = c('AAGCGTGT', 'CTTGACCG', 'TGATTAAC', 'GCCACGTA')
  if(spid == "SI_P01_C9") spseqs = c('AGATCGGT', 'CATCGTCG', 'GTCATATA', 'TCGGACAC')
  if(spid == "SI_P01_D9") spseqs = c('CAAGGGAC', 'ACCTACTG', 'GGGACACA', 'TTTCTTGT')
  if(spid == "SI_P01_E9") spseqs = c('AGTAAGCA', 'TACCGCGG', 'CCGGTAAT', 'GTATCTTC')
  if(spid == "SI_P01_F9") spseqs = c('AGTTAGTT', 'GTACTTAA', 'CACGCACG', 'TCGAGCGC')
  if(spid == "SI_P01_G9") spseqs = c('TTGACTTC', 'GCCGAAGT', 'CAATGGCA', 'AGTCTCAG')
  if(spid == "SI_P01_H9") spseqs = c('GGAATATG', 'ACCTGCCA', 'CTTCATAC', 'TAGGCGGT')
  if(spid == "SI_P01_A10") spseqs = c('ACAGCAAC', 'TTTCGCGA', 'CGCAATTT', 'GAGTTGCG')
  if(spid == "SI_P01_B10") spseqs = c('ACCATTAA', 'CTGGACGT', 'GAACGGTC', 'TGTTCACG')
  if(spid == "SI_P01_C10") spseqs = c('CGTGCTAA', 'TCACTCCT', 'ATCTGATC', 'GAGAAGGG')
  if(spid == "SI_P01_D10") spseqs = c('CTTATTTG', 'GCGGGCAT', 'AGATAACA', 'TACCCGGC')
  if(spid == "SI_P01_E10") spseqs = c('GCACCAGT', 'CGCAGGAG', 'TAGTACCA', 'ATTGTTTC')
  if(spid == "SI_P01_F10") spseqs = c('TCGACAAT', 'GAATACTG', 'ATTCGTGC', 'CGCGTGCA')
  if(spid == "SI_P01_G10") spseqs = c('CGGAGACT', 'TCCTATGA', 'ATACTGAG', 'GATGCCTC')
  if(spid == "SI_P01_H10") spseqs = c('GACCGCCA', 'TCGAATTG', 'ATTTCAGC', 'CGAGTGAT')
  if(spid == "SI_P01_A11") spseqs = c('CTTTCCTT', 'TAGGTAAA', 'ACCAGTCC', 'GGACAGGG')
  if(spid == "SI_P01_B11") spseqs = c('TCCAGATA', 'GATTCGCT', 'CGACATAG', 'ATGGTCGC')
  if(spid == "SI_P01_C11") spseqs = c('GTTTGTGG', 'ACCGAACA', 'TAGACGAC', 'CGACTCTT')
  if(spid == "SI_P01_D11") spseqs = c('GCTACTTC', 'CACCTCAG', 'ATATGAGA', 'TGGGAGCT')
  if(spid == "SI_P01_E11") spseqs = c('ATCGCCAT', 'TCACGGTA', 'GGGTTTCC', 'CATAAAGG')
  if(spid == "SI_P01_F11") spseqs = c('GAACCCGG', 'AGCAGTTA', 'TCGTAGAT', 'CTTGTACC')
  if(spid == "SI_P01_G11") spseqs = c('AGGGCGTT', 'CTATACGC', 'TACATAAG', 'GCTCGTCA')
  if(spid == "SI_P01_H11") spseqs = c('TCTCGACT', 'AGGATCGA', 'CACGATTC', 'GTATCGAG')
  if(spid == "SI_P01_A12") spseqs = c('TTATGGAA', 'ACTACTGT', 'CGGGAACG', 'GACCTCTC')
  if(spid == "SI_P01_B12") spseqs = c('GAAAGACA', 'CGCTACAT', 'ACGCTTGG', 'TTTGCGTC')
  if(spid == "SI_P01_C12") spseqs = c('TAAGCCAC', 'CCGTTATG', 'GGTAATGT', 'ATCCGGCA')
  if(spid == "SI_P01_D12") spseqs = c('GCTGTGTA', 'AGAAACGT', 'CACTCAAC', 'TTGCGTCG')
  if(spid == "SI_P01_E12") spseqs = c('CGCTATCC', 'ACGCGGAA', 'TAAATCGT', 'GTTGCATG')
  if(spid == "SI_P01_F12") spseqs = c('AATTGAAC', 'TGGACCCT', 'CCAGTGGA', 'GTCCATTG')
  if(spid == "SI_P01_G12") spseqs = c('CATGCGTA', 'ACCCGCAC', 'TGATATCG', 'GTGATAGT')
  if(spid == "SI_P01_H12") spseqs = c('TGTGTATA', 'GTCCAGGC', 'CAATCCCT', 'ACGAGTAG')
  
  if(spid == "SI_3A_A1" || spid == "SI_P03_A1") spseqs = c('AAACGGCG', 'CCTACCAT', 'GGCGTTTC', 'TTGTAAGA')
  if(spid == "SI_3A_B1" || spid == "SI_P03_B1") spseqs = c('AGGCTACC', 'CTAGCTGT', 'GCCAACAA', 'TATTGGTG')
  if(spid == "SI_3A_C1" || spid == "SI_P03_C1") spseqs = c('AGACTTTC', 'CCGAGGCA', 'GATGCAGT', 'TTCTACAG')
  if(spid == "SI_3A_D1" || spid == "SI_P03_D1") spseqs = c('AGCTGCGT', 'CAACCATC', 'GTGGAGCA', 'TCTATTAG')
  if(spid == "SI_3A_E1" || spid == "SI_P03_E1") spseqs = c('ACGAAAGC', 'CGCCCGTA', 'GTTTGCCT', 'TAAGTTAG')
  if(spid == "SI_3A_F1" || spid == "SI_P03_F1") spseqs = c('ATCGCTCC', 'CCGTACAG', 'GATAGGTA', 'TGACTAGT')
  if(spid == "SI_3A_G1" || spid == "SI_P03_G1") spseqs = c('ATGCGATT', 'CATATGCG', 'GGATACGA', 'TCCGCTAC')
  if(spid == "SI_3A_H1" || spid == "SI_P03_H1") spseqs = c('AAACTCAT', 'CGGGAGTA', 'GTCACAGG', 'TCTTGTCC')
  if(spid == "SI_3A_A2" || spid == "SI_P03_A2") spseqs = c('AGCCCTTT', 'CAAGTCCA', 'GTGAGAAG', 'TCTTAGGC')
  if(spid == "SI_3A_B2" || spid == "SI_P03_B2") spseqs = c('AAGTTGAT', 'CCCACCCA', 'GGTCGAGC', 'TTAGATTG')
  if(spid == "SI_3A_C2" || spid == "SI_P03_C2") spseqs = c('AATCACTA', 'CCGAGAAC', 'GTAGTGCG', 'TGCTCTGT')
  if(spid == "SI_3A_D2" || spid == "SI_P03_D2") spseqs = c('ACATTCCG', 'CTGCGGTA', 'GACACAAT', 'TGTGATGC')
  if(spid == "SI_3A_E2" || spid == "SI_P03_E2") spseqs = c('AGGCTGGT', 'CACAACTA', 'GTTGGTCC', 'TCATCAAG')
  if(spid == "SI_3A_F2" || spid == "SI_P03_F2") spseqs = c('ATGGTTAG', 'CATTGATA', 'GCAAACGC', 'TGCCCGCT')
  if(spid == "SI_3A_G2" || spid == "SI_P03_G2") spseqs = c('ATAACCTA', 'CGGTGAGC', 'GATCTTAT', 'TCCGAGCG')
  if(spid == "SI_3A_H2" || spid == "SI_P03_H2") spseqs = c('AACGGTCA', 'CCGAACTC', 'GGTCCAAG', 'TTATTGGT')
  if(spid == "SI_3A_A3" || spid == "SI_P03_A3") spseqs = c('AAAGCATA', 'CTGCAGCC', 'GCCTTTAT', 'TGTAGCGG')
  if(spid == "SI_3A_B3" || spid == "SI_P03_B3") spseqs = c('ATTGGACG', 'CAGCTTAC', 'GGCAAGGA', 'TCATCCTT')
  if(spid == "SI_3A_C3" || spid == "SI_P03_C3") spseqs = c('ACGTTACA', 'CGTAGGTT', 'GACGACGG', 'TTACCTAC')
  if(spid == "SI_3A_D3" || spid == "SI_P03_D3") spseqs = c('ACTTCACT', 'CGAAGTTG', 'GAGCACGC', 'TTCGTGAA')
  if(spid == "SI_3A_E3" || spid == "SI_P03_E3") spseqs = c('AACAAGTC', 'CGGCTCCA', 'GTATGTAT', 'TCTGCAGG')
  if(spid == "SI_3A_F3" || spid == "SI_P03_F3") spseqs = c('AGTCTGTA', 'CAGAATAG', 'GCCTCCGT', 'TTAGGACC')
  if(spid == "SI_3A_G3" || spid == "SI_P03_G3") spseqs = c('ATGTCCAG', 'CGACGTCA', 'GCTATAGC', 'TACGAGTT')
  if(spid == "SI_3A_H3" || spid == "SI_P03_H3") spseqs = c('ACACCTAA', 'CGTTTGGG', 'GACAAACC', 'TTGGGCTT')
  if(spid == "SI_3A_A4" || spid == "SI_P03_A4") spseqs = c('AGAACGCC', 'CATGGCAG', 'GTCTTTGA', 'TCGCAATT')
  if(spid == "SI_3A_B4" || spid == "SI_P03_B4") spseqs = c('AGGGACTG', 'CCTCTAAC', 'GACAGGCT', 'TTATCTGA')
  if(spid == "SI_3A_C4" || spid == "SI_P03_C4") spseqs = c('ACATTGGC', 'CTTAGTCA', 'GAGCCCAT', 'TGCGAATG')
  if(spid == "SI_3A_D4" || spid == "SI_P03_D4") spseqs = c('AAATCGTC', 'CTTCGAAT', 'GCGATCGG', 'TGCGATCA')
  if(spid == "SI_3A_E4" || spid == "SI_P03_E4") spseqs = c('AGCTGACG', 'CCGGTGTC', 'GTAAACAT', 'TATCCTGA')
  if(spid == "SI_3A_F4" || spid == "SI_P03_F4") spseqs = c('AACGACAC', 'CGTCCTCT', 'GCATGATA', 'TTGATGGG')
  if(spid == "SI_3A_G4" || spid == "SI_P03_G4") spseqs = c('AGCTTCTC', 'CCTGCGGT', 'GTACAACG', 'TAGAGTAA')
  if(spid == "SI_3A_H4" || spid == "SI_P03_H4") spseqs = c('ACTGGAGC', 'CGGTCGTG', 'GAAATCAA', 'TTCCATCT')
  if(spid == "SI_3A_A5" || spid == "SI_P03_A5") spseqs = c('ATTGGGAA', 'CAGTCTGG', 'GGCATACT', 'TCACACTC')
  if(spid == "SI_3A_B5" || spid == "SI_P03_B5") spseqs = c('ATCGTACT', 'CATCAGTG', 'GGGACTAC', 'TCATGCGA')
  if(spid == "SI_3A_C5" || spid == "SI_P03_C5") spseqs = c('ATGCATTC', 'CACTGACT', 'GGTACGGG', 'TCAGTCAA')
  if(spid == "SI_3A_D5" || spid == "SI_P03_D5") spseqs = c('AGACGGAT', 'CCTTTAGA', 'GTCGACTC', 'TAGACTCG')
  if(spid == "SI_3A_E5" || spid == "SI_P03_E5") spseqs = c('ATCCAAGG', 'CCGTTGAA', 'GGAAGCTC', 'TATGCTCT')
  if(spid == "SI_3A_F5" || spid == "SI_P03_F5") spseqs = c('AACACAGC', 'CGGTTTAG', 'GTACGGCT', 'TCTGACTA')
  if(spid == "SI_3A_G5" || spid == "SI_P03_G5") spseqs = c('ATGAAGTA', 'CGCCGAAC', 'GAAGCTCG', 'TCTTTCGT')
  if(spid == "SI_3A_H5" || spid == "SI_P03_H5") spseqs = c('ATAGTATG', 'CCGCGTCT', 'GGCTCCAC', 'TATAAGGA')
  if(spid == "SI_3A_A6" || spid == "SI_P03_A6") spseqs = c('ACGGGACT', 'CTTTCGAC', 'GAACATGA', 'TGCATCTG')
  if(spid == "SI_3A_B6" || spid == "SI_P03_B6") spseqs = c('AACGCGAA', 'CTATTTGG', 'GCGCACCT', 'TGTAGATC')
  if(spid == "SI_3A_C6" || spid == "SI_P03_C6") spseqs = c('ACTCAGAC', 'CGCTCAGG', 'GAGGTTTA', 'TTAAGCCT')
  if(spid == "SI_3A_D6" || spid == "SI_P03_D6") spseqs = c('ATGCCAAA', 'CCTTATCG', 'GAAGTCTT', 'TGCAGGGC')
  if(spid == "SI_3A_E6" || spid == "SI_P03_E6") spseqs = c('ATTGAAAC', 'CAGCCCGA', 'GCCATTTG', 'TGATGGCT')
  if(spid == "SI_3A_F6" || spid == "SI_P03_F6") spseqs = c('ATGCCGGC', 'CCTAATTA', 'GACTTCCT', 'TGAGGAAG')
  if(spid == "SI_3A_G6" || spid == "SI_P03_G6") spseqs = c('AGCACTGG', 'CATTACAC', 'GTGCGACA', 'TCAGTGTT')
  if(spid == "SI_3A_H6" || spid == "SI_P03_H6") spseqs = c('AAGCATAA', 'CCCATCGC', 'GGTTGATG', 'TTAGCGCT')
  if(spid == "SI_3A_A7" || spid == "SI_P03_A7") spseqs = c('AGGTCATA', 'CTCATCAT', 'GCTGAGGG', 'TAACGTCC')
  if(spid == "SI_3A_B7" || spid == "SI_P03_B7") spseqs = c('AGGGATGA', 'CTTCTGTT', 'GAATGCAC', 'TCCACACG')
  if(spid == "SI_3A_C7" || spid == "SI_P03_C7") spseqs = c('ACACCGGG', 'CATAATCC', 'GGCGGAAT', 'TTGTTCTA')
  if(spid == "SI_3A_D7" || spid == "SI_P03_D7") spseqs = c('AACTTAGA', 'CCGGATCC', 'GGTCGCAT', 'TTAACGTG')
  if(spid == "SI_3A_E7" || spid == "SI_P03_E7") spseqs = c('AAGACGTG', 'CCATGTGT', 'GTTCACAA', 'TGCGTACC')
  if(spid == "SI_3A_F7" || spid == "SI_P03_F7") spseqs = c('ACCCGAGA', 'CAAACTTT', 'GGTTAGAC', 'TTGGTCCG')
  if(spid == "SI_3A_G7" || spid == "SI_P03_G7") spseqs = c('ATTACCGG', 'CAGTAATT', 'GCCGGTAA', 'TGACTGCC')
  if(spid == "SI_3A_H7" || spid == "SI_P03_H7") spseqs = c('AACGGGTG', 'CTAATTCT', 'GCTTCAAC', 'TGGCACGA')
  if(spid == "SI_3A_A8" || spid == "SI_P03_A8") spseqs = c('ATGATACG', 'CCACAGAA', 'GACTGTTC', 'TGTGCCGT')
  if(spid == "SI_3A_B8" || spid == "SI_P03_B8") spseqs = c('ACGTTCAC', 'CAAGGTCT', 'GTTAAGTG', 'TGCCCAGA')
  if(spid == "SI_3A_C8" || spid == "SI_P03_C8") spseqs = c('AGCTCGAG', 'CAGGAAGA', 'GCACGTTT', 'TTTATCCC')
  if(spid == "SI_3A_D8" || spid == "SI_P03_D8") spseqs = c('AATCTTTG', 'CTCAAGAC', 'GGATGAGT', 'TCGGCCCA')
  if(spid == "SI_3A_E8" || spid == "SI_P03_E8") spseqs = c('AGCCTATG', 'CTAACGCA', 'GCTTACAT', 'TAGGGTGC')
  if(spid == "SI_3A_F8" || spid == "SI_P03_F8") spseqs = c('AGTTGGGA', 'CCAGAAAG', 'GTGCCCTC', 'TACATTCT')
  if(spid == "SI_3A_G8" || spid == "SI_P03_G8") spseqs = c('AAGTACTC', 'CTTGGAGA', 'GGAACTCT', 'TCCCTGAG')
  if(spid == "SI_3A_H8" || spid == "SI_P03_H8") spseqs = c('AAGAGCGG', 'CTTGTTAT', 'GGCCCATC', 'TCATAGCA')
  if(spid == "SI_3A_A9" || spid == "SI_P03_A9") spseqs = c('ACAACTTG', 'CTCCAACA', 'GAGTGCGT', 'TGTGTGAC')
  if(spid == "SI_3A_B9" || spid == "SI_P03_B9") spseqs = c('AAGCGTGT', 'CTTGACCG', 'GCCACGTA', 'TGATTAAC')
  if(spid == "SI_3A_C9" || spid == "SI_P03_C9") spseqs = c('AGATCGGT', 'CATCGTCG', 'GTCATATA', 'TCGGACAC')
  if(spid == "SI_3A_D9" || spid == "SI_P03_D9") spseqs = c('ACCTACTG', 'CAAGGGAC', 'GGGACACA', 'TTTCTTGT')
  if(spid == "SI_3A_E9" || spid == "SI_P03_E9") spseqs = c('AGTAAGCA', 'CCGGTAAT', 'GTATCTTC', 'TACCGCGG')
  if(spid == "SI_3A_F9" || spid == "SI_P03_F9") spseqs = c('AGTTAGTT', 'CACGCACG', 'GTACTTAA', 'TCGAGCGC')
  if(spid == "SI_3A_G9" || spid == "SI_P03_G9") spseqs = c('AGTCTCAG', 'CAATGGCA', 'GCCGAAGT', 'TTGACTTC')
  if(spid == "SI_3A_H9" || spid == "SI_P03_H9") spseqs = c('ACCTGCCA', 'CTTCATAC', 'GGAATATG', 'TAGGCGGT')
  if(spid == "SI_3A_A10" || spid == "SI_P03_A10") spseqs = c('ACAGCAAC', 'CGCAATTT', 'GAGTTGCG', 'TTTCGCGA')
  if(spid == "SI_3A_B10" || spid == "SI_P03_B10") spseqs = c('ACCATTAA', 'CTGGACGT', 'GAACGGTC', 'TGTTCACG')
  if(spid == "SI_3A_C10" || spid == "SI_P03_C10") spseqs = c('ATCTGATC', 'CGTGCTAA', 'GAGAAGGG', 'TCACTCCT')
  if(spid == "SI_3A_D10" || spid == "SI_P03_D10") spseqs = c('AGATAACA', 'CTTATTTG', 'GCGGGCAT', 'TACCCGGC')
  if(spid == "SI_3A_E10" || spid == "SI_P03_E10") spseqs = c('ATTGTTTC', 'CGCAGGAG', 'GCACCAGT', 'TAGTACCA')
  if(spid == "SI_3A_F10" || spid == "SI_P03_F10") spseqs = c('ATTCGTGC', 'CGCGTGCA', 'GAATACTG', 'TCGACAAT')
  if(spid == "SI_3A_G10" || spid == "SI_P03_G10") spseqs = c('ATACTGAG', 'CGGAGACT', 'GATGCCTC', 'TCCTATGA')
  if(spid == "SI_3A_H10" || spid == "SI_P03_H10") spseqs = c('ATTTCAGC', 'CGAGTGAT', 'GACCGCCA', 'TCGAATTG')
  if(spid == "SI_3A_A11" || spid == "SI_P03_A11") spseqs = c('ACCAGTCC', 'CTTTCCTT', 'GGACAGGG', 'TAGGTAAA')
  if(spid == "SI_3A_B11" || spid == "SI_P03_B11") spseqs = c('ATGGTCGC', 'CGACATAG', 'GATTCGCT', 'TCCAGATA')
  if(spid == "SI_3A_C11" || spid == "SI_P03_C11") spseqs = c('ACCGAACA', 'CGACTCTT', 'GTTTGTGG', 'TAGACGAC')
  if(spid == "SI_3A_D11" || spid == "SI_P03_D11") spseqs = c('ATATGAGA', 'CACCTCAG', 'GCTACTTC', 'TGGGAGCT')
  if(spid == "SI_3A_E11" || spid == "SI_P03_E11") spseqs = c('ATCGCCAT', 'CATAAAGG', 'GGGTTTCC', 'TCACGGTA')
  if(spid == "SI_3A_F11" || spid == "SI_P03_F11") spseqs = c('AGCAGTTA', 'CTTGTACC', 'GAACCCGG', 'TCGTAGAT')
  if(spid == "SI_3A_G11" || spid == "SI_P03_G11") spseqs = c('AGGGCGTT', 'CTATACGC', 'GCTCGTCA', 'TACATAAG')
  if(spid == "SI_3A_H11" || spid == "SI_P03_H11") spseqs = c('AGGATCGA', 'CACGATTC', 'GTATCGAG', 'TCTCGACT')
  if(spid == "SI_3A_A12" || spid == "SI_P03_A12") spseqs = c('ACTACTGT', 'CGGGAACG', 'GACCTCTC', 'TTATGGAA')
  if(spid == "SI_3A_B12" || spid == "SI_P03_B12") spseqs = c('ACGCTTGG', 'CGCTACAT', 'GAAAGACA', 'TTTGCGTC')
  if(spid == "SI_3A_C12" || spid == "SI_P03_C12") spseqs = c('ATCCGGCA', 'CCGTTATG', 'GGTAATGT', 'TAAGCCAC')
  if(spid == "SI_3A_D12" || spid == "SI_P03_D12") spseqs = c('AGAAACGT', 'CACTCAAC', 'GCTGTGTA', 'TTGCGTCG')
  if(spid == "SI_3A_E12" || spid == "SI_P03_E12") spseqs = c('ACGCGGAA', 'CGCTATCC', 'GTTGCATG', 'TAAATCGT')
  if(spid == "SI_3A_F12" || spid == "SI_P03_F12") spseqs = c('AATTGAAC', 'CCAGTGGA', 'GTCCATTG', 'TGGACCCT')
  if(spid == "SI_3A_G12" || spid == "SI_P03_G12") spseqs = c('ACCCGCAC', 'CATGCGTA', 'GTGATAGT', 'TGATATCG')
  if(spid == "SI_3A_H12" || spid == "SI_P03_H12") spseqs = c('ACGAGTAG', 'CAATCCCT', 'GTCCAGGC', 'TGTGTATA')
  
  
  if(spid == "SI_T2_1") spseqs = c('GGGTGATC', 'TTACCGAT', 'AATGACGA', 'CCCATTCG')
  if(spid == "SI_T2_2") spseqs = c('GGGTCGAA', 'ATCCGCCC', 'TCTATAGT', 'CAAGATTG')
  if(spid == "SI_T2_3") spseqs = c('GCTGATAT', 'TGCCGAGC', 'AAATTGCG', 'CTGACCTA')
  if(spid == "SI_T2_4") spseqs = c('ACTTCTGA', 'TTCATCTT', 'CGACGACG', 'GAGGAGAC')
  if(spid == "SI_T2_5") spseqs = c('GAATACAA', 'AGCATACC', 'TCGGGTTT', 'CTTCCGGG')
  if(spid == "SI_T2_6") spseqs = c('TATTGAGA', 'GTAGTCAG', 'CGCCATTC', 'ACGACGCT')
  if(spid == "SI_T2_7") spseqs = c('AAATCTGT', 'GTCCAACC', 'TCTGGCTG', 'CGGATGAA')
  if(spid == "SI_T2_8") spseqs = c('CCTTGAAC', 'GAAATCGG', 'TGGCCTCT', 'ATCGAGTA')
  
  
  if(spid == "SI_GA_A1" || spid == "SI_P2_A1") spseqs = c('GGTTTACT', 'CTAAACGG', 'TCGGCGTC', 'AACCGTAA')
  if(spid == "SI_GA_A2" || spid == "SI_P2_A2") spseqs = c('TTTCATGA', 'ACGTCCCT', 'CGCATGTG', 'GAAGGAAC')
  if(spid == "SI_GA_A3" || spid == "SI_P2_A3") spseqs = c('CAGTACTG', 'AGTAGTCT', 'GCAGTAGA', 'TTCCCGAC')
  if(spid == "SI_GA_A4" || spid == "SI_P2_A4") spseqs = c('TATGATTC', 'CCCACAGT', 'ATGCTGAA', 'GGATGCCG')
  if(spid == "SI_GA_A5" || spid == "SI_P2_A5") spseqs = c('CTAGGTGA', 'TCGTTCAG', 'AGCCAATT', 'GATACGCC')
  if(spid == "SI_GA_A6" || spid == "SI_P2_A6") spseqs = c('CGCTATGT', 'GCTGTCCA', 'TTGAGATC', 'AAACCGAG')
  if(spid == "SI_GA_A7" || spid == "SI_P2_A7") spseqs = c('ACAGAGGT', 'TATAGTTG', 'CGGTCCCA', 'GTCCTAAC')
  if(spid == "SI_GA_A8" || spid == "SI_P2_A8") spseqs = c('GCATCTCC', 'TGTAAGGT', 'CTGCGATG', 'AACGTCAA')
  if(spid == "SI_GA_A9" || spid == "SI_P2_A9") spseqs = c('TCTTAAAG', 'CGAGGCTC', 'GTCCTTCT', 'AAGACGGA')
  if(spid == "SI_GA_A10" || spid == "SI_P2_A10") spseqs = c('GAAACCCT', 'TTTCTGTC', 'CCGTGTGA', 'AGCGAAAG')
  if(spid == "SI_GA_A11" || spid == "SI_P2_A11") spseqs = c('GTCCGGTC', 'AAGATCAT', 'CCTGAAGG', 'TGATCTCA')
  if(spid == "SI_GA_A12" || spid == "SI_P2_A12") spseqs = c('AGTGGAAC', 'GTCTCCTT', 'TCACATCA', 'CAGATGGG')
  if(spid == "SI_GA_B1" || spid == "SI_P2_B1") spseqs = c('GTAATCTT', 'TCCGGAAG', 'AGTTCGGC', 'CAGCATCA')
  if(spid == "SI_GA_B2" || spid == "SI_P2_B2") spseqs = c('TACTCTTC', 'CCTGTGCG', 'GGACACGT', 'ATGAGAAA')
  if(spid == "SI_GA_B3" || spid == "SI_P2_B3") spseqs = c('GTGTATTA', 'TGTGCGGG', 'ACCATAAC', 'CAACGCCT')
  if(spid == "SI_GA_B4" || spid == "SI_P2_B4") spseqs = c('ACTTCATA', 'GAGATGAC', 'TGCCGTGG', 'CTAGACCT')
  if(spid == "SI_GA_B5" || spid == "SI_P2_B5") spseqs = c('AATAATGG', 'CCAGGGCA', 'TGCCTCAT', 'GTGTCATC')
  if(spid == "SI_GA_B6" || spid == "SI_P2_B6") spseqs = c('CGTTAATC', 'GCCACGCT', 'TTACTCAG', 'AAGGGTGA')
  if(spid == "SI_GA_B7" || spid == "SI_P2_B7") spseqs = c('AAACCTCA', 'GCCTTGGT', 'CTGGACTC', 'TGTAGAAG')
  if(spid == "SI_GA_B8" || spid == "SI_P2_B8") spseqs = c('AAAGTGCT', 'GCTACCTG', 'TGCTGTAA', 'CTGCAAGC')
  if(spid == "SI_GA_B9" || spid == "SI_P2_B9") spseqs = c('CTGTAACT', 'TCTAGCGA', 'AGAGTGTG', 'GACCCTAC')
  if(spid == "SI_GA_B10" || spid == "SI_P2_B10") spseqs = c('ACCGTATG', 'GATTAGAT', 'CTGACTGA', 'TGACGCCC')
  if(spid == "SI_GA_B11" || spid == "SI_P2_B11") spseqs = c('GTTCCTCA', 'AGGTACGC', 'TAAGTATG', 'CCCAGGAT')
  if(spid == "SI_GA_B12" || spid == "SI_P2_B12") spseqs = c('TACCACCA', 'CTAAGTTT', 'GGGTCAAG', 'ACTGTGGC')
  if(spid == "SI_GA_C1" || spid == "SI_P2_C1") spseqs = c('CCACTTAT', 'AACTGGCG', 'TTGGCATA', 'GGTAACGC')
  if(spid == "SI_GA_C2" || spid == "SI_P2_C2") spseqs = c('CCTAGACC', 'ATCTCTGT', 'TAGCTCTA', 'GGAGAGAG')
  if(spid == "SI_GA_C3" || spid == "SI_P2_C3") spseqs = c('TCAGCCGT', 'CAGAGGCC', 'GGTCAATA', 'ATCTTTAG')
  if(spid == "SI_GA_C4" || spid == "SI_P2_C4") spseqs = c('ACAATTCA', 'TGCGCAGC', 'CATCACTT', 'GTGTGGAG')
  if(spid == "SI_GA_C5" || spid == "SI_P2_C5") spseqs = c('CGACTTGA', 'TACAGACT', 'ATTGCGTG', 'GCGTACAC')
  if(spid == "SI_GA_C6" || spid == "SI_P2_C6") spseqs = c('ATTACTTC', 'TGCGAACT', 'GCATTCGG', 'CAGCGGAA')
  if(spid == "SI_GA_C7" || spid == "SI_P2_C7") spseqs = c('GTCTCTCG', 'AATCTCTC', 'CGGAGGGA', 'TCAGAAAT')
  if(spid == "SI_GA_C8" || spid == "SI_P2_C8") spseqs = c('GTTGAGAA', 'AGATCTGG', 'TCGATACT', 'CACCGCTC')
  if(spid == "SI_GA_C9" || spid == "SI_P2_C9") spseqs = c('GCGCAGAA', 'ATCTTACC', 'TATGGTGT', 'CGAACCTG')
  if(spid == "SI_GA_C10" || spid == "SI_P2_C10") spseqs = c('TCTCAGTG', 'GAGACTAT', 'CGCTTAGC', 'ATAGGCCA')
  if(spid == "SI_GA_C11" || spid == "SI_P2_C11") spseqs = c('GAGGATCT', 'AGACCATA', 'TCCTGCGC', 'CTTATGAG')
  if(spid == "SI_GA_C12" || spid == "SI_P2_C12") spseqs = c('TCTCGTTT', 'GGCTAGCG', 'ATGACCGC', 'CAAGTAAA')
  if(spid == "SI_GA_D1" || spid == "SI_P2_D1") spseqs = c('CACTCGGA', 'GCTGAATT', 'TGAAGTAC', 'ATGCTCCG')
  if(spid == "SI_GA_D2" || spid == "SI_P2_D2") spseqs = c('TAACAAGG', 'GGTTCCTC', 'ATCATGCA', 'CCGGGTAT')
  if(spid == "SI_GA_D3" || spid == "SI_P2_D3") spseqs = c('ACATTACT', 'TTTGGGTA', 'CAGCCCAC', 'GGCAATGG')
  if(spid == "SI_GA_D4" || spid == "SI_P2_D4") spseqs = c('CCCTAACA', 'ATTCCGAT', 'TGGATTGC', 'GAAGGCTG')
  if(spid == "SI_GA_D5" || spid == "SI_P2_D5") spseqs = c('CTCGTCAC', 'GATCAGCA', 'ACAACAGG', 'TGGTGTTT')
  if(spid == "SI_GA_D6" || spid == "SI_P2_D6") spseqs = c('CATGCGAT', 'TGATATTC', 'GTGATCGA', 'ACCCGACG')
  if(spid == "SI_GA_D7" || spid == "SI_P2_D7") spseqs = c('ATTTGCTA', 'TAGACACC', 'CCACAGGG', 'GGCGTTAT')
  if(spid == "SI_GA_D8" || spid == "SI_P2_D8") spseqs = c('GCAACAAA', 'TAGTTGTC', 'CGCCATCG', 'ATTGGCGT')
  if(spid == "SI_GA_D9" || spid == "SI_P2_D9") spseqs = c('AGGAGATG', 'GATGTGGT', 'CTACATCC', 'TCCTCCAA')
  if(spid == "SI_GA_D10" || spid == "SI_P2_D10") spseqs = c('CAATACCC', 'TGTCTATG', 'ACCACGAA', 'GTGGGTGT')
  if(spid == "SI_GA_D11" || spid == "SI_P2_D11") spseqs = c('CTTTGCGG', 'TGCACAAA', 'AAGCAGTC', 'GCAGTTCT')
  if(spid == "SI_GA_D12" || spid == "SI_P2_D12") spseqs = c('GCACAATG', 'CTTGGTAC', 'TGCACCGT', 'AAGTTGCA')
  if(spid == "SI_GA_E1" || spid == "SI_P2_E1") spseqs = c('TGGTAAAC', 'GAAAGGGT', 'ACTGCTCG', 'CTCCTCTA')
  if(spid == "SI_GA_E2" || spid == "SI_P2_E2") spseqs = c('GTGGTACC', 'TACTATAG', 'ACAAGGTA', 'CGTCCCGT')
  if(spid == "SI_GA_E3" || spid == "SI_P2_E3") spseqs = c('AGGTATTG', 'CTCCTAGT', 'TCAAGGCC', 'GATGCCAA')
  if(spid == "SI_GA_E4" || spid == "SI_P2_E4") spseqs = c('TTCGCCCT', 'GGATGGGC', 'AATCAATG', 'CCGATTAA')
  if(spid == "SI_GA_E5" || spid == "SI_P2_E5") spseqs = c('CATTAGCG', 'TTCGCTGA', 'ACAAGAAT', 'GGGCTCTC')
  if(spid == "SI_GA_E6" || spid == "SI_P2_E6") spseqs = c('CTGCGGCT', 'GACTCAAA', 'AGAAACTC', 'TCTGTTGG')
  if(spid == "SI_GA_E7" || spid == "SI_P2_E7") spseqs = c('CACGCCTT', 'GTATATAG', 'TCTCGGGC', 'AGGATACA')
  if(spid == "SI_GA_E8" || spid == "SI_P2_E8") spseqs = c('ATAGTTAC', 'TGCTGAGT', 'CCTACGTA', 'GAGCACCG')
  if(spid == "SI_GA_E9" || spid == "SI_P2_E9") spseqs = c('TTGTTTCC', 'GGAGGAGG', 'CCTAACAA', 'AACCCGTT')
  if(spid == "SI_GA_E10" || spid == "SI_P2_E10") spseqs = c('AAATGTGC', 'GGGCAAAT', 'TCTATCCG', 'CTCGCGTA')
  if(spid == "SI_GA_E11" || spid == "SI_P2_E11") spseqs = c('AAGCGCTG', 'CGTTTGAT', 'GTAGCACA', 'TCCAATGC')
  if(spid == "SI_GA_E12" || spid == "SI_P2_E12") spseqs = c('ACCGGCTC', 'GAGTTAGT', 'CGTCCTAG', 'TTAAAGCA')
  if(spid == "SI_GA_F1" || spid == "SI_P2_F1") spseqs = c('GTTGCAGC', 'TGGAATTA', 'CAATGGAG', 'ACCCTCCT')
  if(spid == "SI_GA_F2" || spid == "SI_P2_F2") spseqs = c('TTTACATG', 'CGCGATAC', 'ACGCGGGT', 'GAATTCCA')
  if(spid == "SI_GA_F3" || spid == "SI_P2_F3") spseqs = c('TTCAGGTG', 'ACGGACAT', 'GATCTTGA', 'CGATCACC')
  if(spid == "SI_GA_F4" || spid == "SI_P2_F4") spseqs = c('CCCAATAG', 'GTGTCGCT', 'AGAGTCGC', 'TATCGATA')
  if(spid == "SI_GA_F5" || spid == "SI_P2_F5") spseqs = c('GACTACGT', 'CTAGCGAG', 'TCTATATC', 'AGGCGTCA')
  if(spid == "SI_GA_F6" || spid == "SI_P2_F6") spseqs = c('CGGAGCAC', 'GACCTATT', 'ACTTAGGA', 'TTAGCTCG')
  if(spid == "SI_GA_F7" || spid == "SI_P2_F7") spseqs = c('CGTGCAGA', 'AACAAGAT', 'TCGCTTCG', 'GTATGCTC')
  if(spid == "SI_GA_F8" || spid == "SI_P2_F8") spseqs = c('CATGAACA', 'TCACTCGC', 'AGCTGGAT', 'GTGACTTG')
  if(spid == "SI_GA_F9" || spid == "SI_P2_F9") spseqs = c('CAAGCTCC', 'GTTCACTG', 'TCGTGAAA', 'AGCATGGT')
  if(spid == "SI_GA_F10" || spid == "SI_P2_F10") spseqs = c('GCTTGGCT', 'AAACAAAC', 'CGGGCTTA', 'TTCATCGG')
  if(spid == "SI_GA_F11" || spid == "SI_P2_F11") spseqs = c('GCGAGAGT', 'TACGTTCA', 'AGTCCCAC', 'CTATAGTG')
  if(spid == "SI_GA_F12" || spid == "SI_P2_F12") spseqs = c('TGATGCAT', 'GCTACTGA', 'CACCTGCC', 'ATGGAATG')
  if(spid == "SI_GA_G1" || spid == "SI_P2_G1") spseqs = c('ATGAATCT', 'GATCTCAG', 'CCAGGAGC', 'TGCTCGTA')
  if(spid == "SI_GA_G2" || spid == "SI_P2_G2") spseqs = c('TGATTCTA', 'ACTAGGAG', 'CAGCCACT', 'GTCGATGC')
  if(spid == "SI_GA_G3" || spid == "SI_P2_G3") spseqs = c('CCTCATTC', 'AGCATCCG', 'GTGGCAAT', 'TAATGGGA')
  if(spid == "SI_GA_G4" || spid == "SI_P2_G4") spseqs = c('GCGATGTG', 'AGATACAA', 'TTTCCACT', 'CACGGTGC')
  if(spid == "SI_GA_G5" || spid == "SI_P2_G5") spseqs = c('GAGCAAGA', 'TCTGTGAT', 'CGCAGTTC', 'ATATCCCG')
  if(spid == "SI_GA_G6" || spid == "SI_P2_G6") spseqs = c('CTGACGCG', 'GGTCGTAC', 'TCCTTCTT', 'AAAGAAGA')
  if(spid == "SI_GA_G7" || spid == "SI_P2_G7") spseqs = c('GGTATGCA', 'CTCGAAAT', 'ACACCTTC', 'TAGTGCGG')
  if(spid == "SI_GA_G8" || spid == "SI_P2_G8") spseqs = c('TATGAGCT', 'CCGATAGC', 'ATACCCAA', 'GGCTGTTG')
  if(spid == "SI_GA_G9" || spid == "SI_P2_G9") spseqs = c('TAGGACGT', 'ATCCCACA', 'GGAATGTC', 'CCTTGTAG')
  if(spid == "SI_GA_G10" || spid == "SI_P2_G10") spseqs = c('TCGCCAGC', 'AATGTTAG', 'CGATAGCT', 'GTCAGCTA')
  if(spid == "SI_GA_G11" || spid == "SI_P2_G11") spseqs = c('TTATCGTT', 'AGCAGAGC', 'CATCTCCA', 'GCGGATAG')
  if(spid == "SI_GA_G12" || spid == "SI_P2_G12") spseqs = c('ATTCTAAG', 'CCCGATTA', 'TGGAGGCT', 'GAATCCGC')
  if(spid == "SI_GA_H1" || spid == "SI_P2_H1") spseqs = c('GTATGTCA', 'TGTCAGAC', 'CACGTCGG', 'ACGACATT')
  if(spid == "SI_GA_H2" || spid == "SI_P2_H2") spseqs = c('TAATGACC', 'ATGCCTTA', 'GCCGAGAT', 'CGTATCGG')
  if(spid == "SI_GA_H3" || spid == "SI_P2_H3") spseqs = c('CCAAGATG', 'AGGCCCGA', 'TACGTGAC', 'GTTTATCT')
  if(spid == "SI_GA_H4" || spid == "SI_P2_H4") spseqs = c('GCCATTCC', 'CAAGAATT', 'TTGCCGGA', 'AGTTGCAG')
  if(spid == "SI_GA_H5" || spid == "SI_P2_H5") spseqs = c('CCACTACA', 'GATTCTGG', 'TGCGGCTT', 'ATGAAGAC')
  if(spid == "SI_GA_H6" || spid == "SI_P2_H6") spseqs = c('TAGGATAA', 'CCTTTGTC', 'GTACGCGG', 'AGCACACT')
  if(spid == "SI_GA_H7" || spid == "SI_P2_H7") spseqs = c('AGCTATCA', 'CATATAAC', 'TCAGGGTG', 'GTGCCCGT')
  if(spid == "SI_GA_H8" || spid == "SI_P2_H8") spseqs = c('TTGTTGAT', 'GCTCAACC', 'CAAAGTGG', 'AGCGCCTA')
  if(spid == "SI_GA_H9" || spid == "SI_P2_H9") spseqs = c('ACACTGTT', 'CAGGATGG', 'GGCTGAAC', 'TTTACCCA')
  if(spid == "SI_GA_H10" || spid == "SI_P2_H10") spseqs = c('GTAATTGC', 'AGTCGCTT', 'CACGAGAA', 'TCGTCACG')
  if(spid == "SI_GA_H11" || spid == "SI_P2_H11") spseqs = c('GGCGAGTA', 'ACTTCTAT', 'CAAATACG', 'TTGCGCGC')
  if(spid == "SI_GA_H12" || spid == "SI_P2_H12") spseqs = c('GACAGCAT', 'TTTGTACA', 'AGGCCGTG', 'CCATATGC')
  spseqs
}

.rstest <- function(r, coef) r * (1 + 1/r)^(1 + coef)

.averaging_transform <- function(r, nr){
    d = c(1, r[ 2:length(r) ] - r[1:(length(r)-1)])
    dr = c( 0.5 * (d[2:length(d)] + d[1:(length(d)-1)]),d[ length(d) ])
   	nr/dr
}


.simple.Good.Turing.Freq<-function( raw.freq ){
  n_zero <- sum(raw.freq == 0)
  times <- table(raw.freq[raw.freq>0])
  obs <- sort(as.numeric(names(times)))
  times <- times[as.character(obs)]

  nobs = length(obs)
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
  	min.no.turing = min(which(!xy.sim))
  	r <- GT
  	if(min.no.turing>0)r[min.no.turing:nobs] <- LGT[min.no.turing:nobs]
  	p_zero <- times[1] / total_observed
  	rs <- (1-p_zero) * obs * r / sum(r*times*obs/total_observed)
  	total_rs <- sum(rs * times)
  	puniq <- (1-p_zero) * rs / total_rs
  	pres <- rep(p_zero / n_zero, length(raw.freq))
  	pres [raw.freq>0] <- puniq[match( raw.freq[raw.freq>0], obs )]
  	ret <- list(p0=p_zero, p=pres)
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
	#saveRDS(obs.tab,"del4-obs.tab.RDS")
	#saveRDS(obs.per.spe, "del4-obs.per.spe.RDS")
	sgtr <- .mySGTsorted(as.numeric(names(obs.tab)), obs.tab)
	#saveRDS(sgtr, "del4-sgtr.RDS")
	res <- obs.per.spe
	res[obs.per.spe!=0] <- sgtr$p[match(as.character(obs.per.spe[obs.per.spe!=0]), names(obs.tab))]
	n_zero <- sum(obs.per.spe==0)
	res[obs.per.spe==0] <- sgtr$p0/n_zero
	res
}

.mySGTsorted<-function( obs, times ){
	if(any(obs != sort(obs)))stop("Observations have to be sorted.")
	nobs = length(obs)
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


.cellCounts_try_cellbarcode <- function( input.directory, sample.sheet, cell.barcode.list, nreads.testing ){ # the three parameters can only be one string, not strings!

	cmd <- paste0(c(input.directory, sample.sheet, cell.barcode.list, as.character(nreads.testing)), collapse=.R_param_splitor)
	rvs <- as.integer(rep(0,5))
	C_args <- .C("R_try_cell_barcode_wrapper",nargs=as.integer(4),argv=as.character(cmd),retv=rvs ,PACKAGE="Rsubread")
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

.find_best_cellbarcode <- function( input.directory, sample.sheet){
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
		barcode_res <- .cellCounts_try_cellbarcode(input.directory[1], sample.sheet[1], listfile, 30000)
		if(length(barcode_res)<3)stop("ERROR: the input sample cannot be processed.")
		sample.good.rate <- barcode_res[2]/barcode_res[1]
		cell.good.rate <- barcode_res[3]/barcode_res[1]
		max.cell.good <- max(max.cell.good, cell.good.rate)
		#cat(sprintf("Sample supporting rate : %.1f%% ; cell supporting rate : %.1f%%.\n", sample.good.rate*100., cell.good.rate*100.))
		if(sample.good.rate < 0.5)cat(sprintf("WARNING: there are only %.1f%% reads having known sample indices. Please check if the sample sheet is correct.\n", sample.good.rate*100.))
		if(cell.good.rate > 0.6){
			cat(sprintf("Found cell-barcode list '%s' for the input data: supported by %.1f%% reads.\n", libf, cell.good.rate*100.))
			return(listfile)
		}
	}

	stop(sprintf("ERROR: no known cell barcode set was found for the data set. The highest percentage of cell-barcode matched reads is %.1f%%\n", max.cell.good*100.))
}

library(Matrix)
.read.sparse.mat <- function (fn){
#  print(fn)
  mtx <- readMM(paste0(fn, ".spmtx"))
  coln <- read.delim(paste0(fn, ".BCtab"), stringsAsFactors=F, header=F)$V1
  rown <- read.delim(paste0(fn, ".GENEtab"), stringsAsFactors=F, header=F)$V1
  colnames(mtx) <- coln
  rownames(mtx) <- rown

  mtx
}

.simu.multinomial <- function(candi.mat, gene.profile.freq , times=10000){
  bcsizes <- sort(unique( colSums(candi.mat) ))
  ret.nUMI.LLH.tab <- list()
  N10000.LLH <- NA
  Old_bs_one <- NA
  log_E_GTE <- log(gene.profile.freq)
  for(bcsize_one in bcsizes){
    UMI_step_diff <- NA
    if(!any(is.na(Old_bs_one))) UMI_step_diff <- bcsize_one - Old_bs_one
  
    if(any(is.na(N10000.LLH))){
        N10000 <- rmultinom(n=times, size=bcsize_one, prob=gene.profile.freq )
      N10000.LLH <- apply(N10000, 2, function(x) dmultinom(x, prob=gene.profile.freq, log =T ))
    }else if(UMI_step_diff >= 1000){
      UMI_step_diff_N10000 <- rmultinom(n=times, size=UMI_step_diff, prob=gene.profile.freq )
      N10000 <- N10000 + UMI_step_diff_N10000
      N10000.LLH <- apply(N10000, 2, function(x) dmultinom(x, prob=gene.profile.freq, log =T ))
    }else{
      for(curi in (Old_bs_one+1):bcsize_one){
      UMI_step_indices_N10000 <- sample(1:nrow(candi.mat), size=times, prob=gene.profile.freq, replace=T)
      for(idi in 1:times){
        N10000[ UMI_step_indices_N10000[idi] ,idi ] <- N10000[ UMI_step_indices_N10000[idi] ,idi ]+1
      }
      UMI_step_values_N10000 <- rep(0, times)
      for(idi in 1:times){
        UMI_step_values_N10000[idi] <- N10000[UMI_step_indices_N10000[idi] ,idi ]
      }
      #print(UMI_step_values_N10000)
      N10000.LLH <- N10000.LLH + log_E_GTE[ UMI_step_indices_N10000 ] + log((curi +1.) / UMI_step_values_N10000)
      }
    }
  
    # cat("\n\n ++++++++++ ", bcsize_one," +++++++++ \n")
    # print(summary(N10000.LLH))
    do_CTRL <- F
    if(do_CTRL && 0== bcsize_one %% 15){
      CTRL <- rmultinom(n=times, size=bcsize_one, prob=gene.profile.freq )
      CTRL.LLH <- apply(CTRL, 2, function(x) dmultinom(x, prob=gene.profile.freq, log =T ))
      #print(summary(CTRL.LLH))
    }
    ret.nUMI.LLH.tab[[ bcsize_one ]] <- N10000.LLH
    Old_bs_one <- bcsize_one
  }
  ret.nUMI.LLH.tab
}

.is.no.candBC <- function(fn){
  info = file.info(fn)
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
    print("Summary of the ambient RNA profile (proportions)")
    print(summary(gte))
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
      actual.pvalues[candi] <- (1+sum(cand.simu.pvs > cand.actual.pv))/(length(cand.simu.pvs)+1)
    }
  
    # p-value => FDR
    actual.FDR <- p.adjust(actual.pvalues, method='BH')
  
    # select cells that has FDR < cutoff
    FDR.Cutoff <- 0.01
    Rescured.Barcodes <- names(actual.FDR)[actual.FDR <= FDR.Cutoff]
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

.load.one.scSample <- function( BAM.name, FC.gene.ids, sample.no, use.meta.features, annot.tab){
  fname <- sprintf("%s.scRNA.%03d", BAM.name, sample.no)
  highconf <- as.matrix(.read.sparse.mat(paste0(fname,".HighConf")))
  rescued <- .cellCounts.rescue(BAM.name, FC.gene.ids, sample.no)

  if(use.meta.features){
    if(!any(is.na(rescued))) rescued <- rescued[rowSums(rescued)>0,]
    ncolRescued <- 0
    if(is.na(rescued)){
      ret <- matrix(0,ncol=ncol(highconf), nrow=length(FC.gene.ids))
      colnames(ret) <- colnames(highconf)
    }else{
      ret <- matrix(0,ncol=ncol(highconf)+ncol(rescued), nrow=length(FC.gene.ids))
      colnames(ret) <- c( colnames(highconf), colnames(rescued) )
    }
    rownames(ret) <- FC.gene.ids 
    ret[rownames(highconf), colnames(highconf) ] <- highconf
    if(!is.na(rescued)) ret[rownames(rescued), colnames(rescued) ] <- rescued

    list(Counts=ret, HighConfidneceCell=colnames(ret) %in% colnames(highconf))
  }else{
    fcmat <- .match.exons(annot.tab,highconf)
    rownames(fcmat)<-NULL
    list(Counts=fcmat, HighConfidneceCell=rep(T, ncol(highconf))) 
  }
}

.load.all.scSamples <- function( BAM.name, FC.gene.ids, use.meta.features, annot.tab){
  sum.tab <- read.delim(paste0(BAM.name,".scRNA.SampleTable"), stringsAsFactors=F)
  ret <- list()
  for(roiw in 1:nrow(sum.tab)){
    sname <- sum.tab$SampleName[roiw]
    sid <- sum.tab$Index[roiw]
    count.tab <- .load.one.scSample(BAM.name, FC.gene.ids, sid, use.meta.features, annot.tab)
	ret[[sprintf("Sample.%d",sid)]] <- count.tab
  }
  ret[["Sample.Table"]] <- sum.tab

  ret
}

cellCounts <- function(index, input.directory, output.BAM, sample.sheet, cell.barcode.list=NULL, input.mode="BCL", aligner="align", nthreads=16, annot.inbuilt="mm10",annot.ext=NULL,isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",GTF.attrType.extra=NULL,chrAliases=NULL,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,countMultiMappingReads=TRUE,reportAllJunctions=FALSE){
  set.seed(0)
  input.directory <- .check_and_NormPath(input.directory,  mustWork=T, opt="input.directory")
  sample.sheet <- .check_and_NormPath(sample.sheet,  mustWork=T, opt="sample.sheet")
  output.BAM <- .check_and_NormPath(output.BAM,  mustWork=F, opt="output.BAM")
  if(!is.null(aligner)) aligner <- match.arg(aligner,c("subjunc","align")) 
  fc <- list()

  if(length(input.directory) != length(output.BAM) || length(input.directory) != length( sample.sheet ))stop("The arguments to the input.directory, output.BAM and sample.sheet options must have the same length.")
  if(is.null(cell.barcode.list)){
    cell.barcode.list <- .find_best_cellbarcode(input.directory, sample.sheet)
  }else{
    cell.barcode.list <- .check_and_NormPath(cell.barcode.list, mustWork=T, opt="cell.barcode.list")
  }

  raw.fc.annot <- NA
  for(ii in 1:length(input.directory)){
	  input.1 <- input.directory[ii]
	  output.1 <- output.BAM[ii]
	  sample.1 <- sample.sheet[ii]
	  if(is.null(aligner)){
		if(!file.exists(output.1)) stop("No aligner is specified but the BAM file does not exist. Please specify 'align' or 'subjunc' as the aligner.")
	  }else if(aligner=="align"){
		align(index, input.1, output_file=output.1, nthreads=nthreads, isBCLinput=TRUE)
	  }else if(aligner=="subjunc"){
		subjunc(index, input.1, output_file=output.1, nthreads=nthreads, isBCLinput=TRUE, reportAllJunctions=reportAllJunctions)
	  }
      raw.fc<-featureCounts(output.1, annot.inbuilt=annot.inbuilt, annot.ext=annot.ext, isGTFAnnotationFile=isGTFAnnotationFile, GTF.featureType=GTF.featureType, GTF.attrType=GTF.attrType, GTF.attrType.extra=GTF.attrType.extra, chrAliases=chrAliases, useMetaFeatures=useMetaFeatures, allowMultiOverlap=allowMultiOverlap, countMultiMappingReads=countMultiMappingReads, sampleSheet=sample.1, cellBarcodeList=cell.barcode.list, nthreads=nthreads)
      if(is.na(raw.fc.annot)) raw.fc.annot<-raw.fc$annotation
	  fc[[paste0("Dataset.", ii)]] <- .load.all.scSamples(output.1, as.character(raw.fc.annot$GeneID), useMetaFeatures, raw.fc.annot)
  }
  fc[["Input.Files"]] <- input.directory
  fc[["Annotation"]] <- raw.fc.annot

  fc
}

