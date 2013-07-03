\name{createAnnotationFile}
\alias{createAnnotationFile}
\alias{write.Rsubread}
\title{
  Create an annotation file from a GRanges object, suitable for featureCounts() 
}
\description{
  Any of rtracklayer::import.bed('samplesubjunc.bed', asRangedData=FALSE),
  unlist(spliceGraph(TxDb)), transcripts(TxDb), exons(TxDb), or features(FDB)
  will produce a GRanges object containing usable features for read counting.

  This function converts a suitably streamlined GRanges object into annotations
  which can then be used by featureCounts() to quickly count aligned reads.

  The GRanges object must contain an elementMetadata column named 'id'.
}
\usage{
createAnnotationFile(GR)
write.Rsubread(GR)
}
\arguments{
  \item{GR}{
    The GRanges object to convert to an Rsubread annotation file
  }
}
\value{
  A data frame with five columns named \code{GeneID}, \code{Chr}, \code{Start}, \code{End} and \code{Strand}.
}
\author{
  Tim Triche, Jr. and Wei Shi
}
\examples{
\dontrun{
  library(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
  hg19LincRNAs <- transcripts(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
  names(values(hg19LincRNAs)) <- gsub('tx_id','id',names(values(hg19LincRNAs)))
  annot_for_featureCounts <- createAnnotationFile(hg19LincRNAs)
}
}
