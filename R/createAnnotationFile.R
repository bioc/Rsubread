createAnnotationFile <- function(GR, fname=NULL) {

  if(!(is(GR, 'GRanges'))) { # {{{
    message('To create Rsubread annotations from a TranscriptDb or FeatureDb,')
    message('extract a GRanges using exons(), transcripts(), or features().')
    stop('Aborting.')
  } # }}}
  if(!('id' %in% names(values(GR)))) { # {{{
    message('Your GRanges elementMetadata must contain a column called "id".')
    stop('Aborting.')
  } # }}}
  if(is.null(fname)) { # {{{
    fname = paste(as.character(match.call()["GR"]), 'annot', 'txt', sep='.')
  } # }}}

  values(GR)[,'id'] <- as.numeric(as.factor(unlist(values(GR)[,'id'])))
  write.table(sep="\t", quote=FALSE, file=fname, row.names=FALSE,
              data.frame(entrezid=values(GR)[,'id'],
                         chromosome=as(seqnames(GR), 'character'), 
                         chr_start=start(GR),
                         chr_end=end(GR))
              ) 
  message(paste('Annotations were sucessfully written to', fname))
}

# an alias for backward compatibility:
write.Rsubread <- createAnnotationFile
