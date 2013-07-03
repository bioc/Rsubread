createAnnotationFile <- function(GR) {

  if(!(is(GR, 'GRanges'))) {
    message('To create Rsubread annotations from a TranscriptDb or FeatureDb,')
    message('extract a GRanges using exons(), transcripts(), or features().')
    stop('Aborting.')
  }
  
  if(!('id' %in% names(values(GR)))) {
    message('Your GRanges elementMetadata must contain a column called "id".')
    stop('Aborting.')
  }

  annot <- data.frame(GeneID=as.character(unlist(values(GR)[,'id'])),
                         Chr=as.character(seqnames(GR)), 
                         Start=start(GR),
                         End=end(GR),
						 Strand=as.character(strand(GR)),
						 stringsAsFactors=FALSE)
  return(annot)
}

# an alias for backward compatibility:
 write.Rsubread <- createAnnotationFile
