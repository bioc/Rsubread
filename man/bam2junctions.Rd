\name{bam2junctions}
\alias{bam2junctions}
\title{Generate junction-level UMI counts in single-cell RNA-seq data}

\description{
\code{bam2junctions} constructs a UMI count matrix of splice junctions across cells
from a single-cell RNA-seq BAM file. Splice junctions are inferred directly from
read alignments, and UMI deduplication is performed at the junction level within each cell.
}

\usage{
bam2junctions(
  BAMfile,
  cell.barcodes
)
}

\arguments{

\item{BAMfile}{
Character string specifying the input BAM file containing read alignments.
The BAM file must include cell and UMI barcodes stored in standard SAM tags.
Supported aligners include Cell Ranger and STARsolo.
}

\item{cell.barcodes}{
Character vector of cell barcodes to include in the analysis.
Only reads assigned to these barcodes are considered.
}

}

\details{

\strong{Overview}

This function extracts splice junctions from aligned reads and constructs
a sparse matrix of UMI counts, where rows correspond to junctions and
columns correspond to cells.

\strong{BAM Tags Used}

The function relies on the following BAM tags:

\itemize{

\item \strong{CB tag (Corrected Cell Barcode):}
The \code{CB} tag stores the error-corrected cell barcode assigned to a read.
This barcode represents the inferred cell of origin after processing by the aligner
(e.g., Cell Ranger or STARsolo).

\item \strong{UR tag (Raw UMI Barcode):}
The \code{UR} tag stores the raw (uncorrected) UMI sequence as observed in the read.
It is used to identify low-complexity UMIs (e.g., homopolymers such as \code{AAAA...})
or UMIs containing ambiguous bases (e.g., \code{N}).

\item \strong{UB tag (Corrected UMI Barcode):}
The \code{UB} tag stores the corrected UMI sequence after error correction.
This tag is used for molecule identification and UMI deduplication.

\item \strong{GX tag (Gene Annotation):}
The \code{GX} tag indicates the gene(s) to which a read is assigned.
Reads lacking a \code{GX} tag or with multiple assignments may be handled differently
depending on the aligner.

\item \strong{xf tag (Feature Flags; Cell Ranger-specific):}
The \code{xf} tag is a bitwise flag used by Cell Ranger to encode read-level annotations,
including mapping confidence and transcript compatibility. Individual bits correspond
to specific properties and are used for filtering (see below).
}

\strong{Read Filtering}

Only reads satisfying all of the following conditions are retained:

\itemize{

\item The read contains a \code{CB} tag matching one of the specified
\code{cell.barcodes}.

\item The read has a valid raw UMI in the \code{UR} tag, excluding low-complexity
UMIs (e.g., homopolymers) and UMIs containing ambiguous bases (e.g., \code{N}).

\item The read contains a valid corrected UMI in the \code{UB} tag.
}

Additional filtering is applied depending on the aligner:

\itemize{

\item \strong{STARsolo:}
All reads passing the above criteria are retained without further filtering.

\item \strong{Cell Ranger:}
\itemize{
  \item Reads without gene assignment (i.e., missing \code{GX} tag) are retained.
  \item For reads assigned to one or more genes, the \code{xf} tag is evaluated.
  Reads are retained only if:

\itemize{
  \item bits 1 and 5 are set, indicating that the read is confidently mapped
  to the transcriptome and overlaps exonic regions; and

  \item bits 2, 3, and 6 are not set, excluding reads that are genomic-only,
  antisense, or intronic.
}

These bitwise conditions match the criteria used by Cell Ranger to select reads
for UMI counting in its gene expression quantification. 
}
}

\strong{Junction Extraction}

Splice junctions are identified directly from the CIGAR string of each read,
specifically from \code{N} (skipped region) operations representing introns.

Each junction is recorded as genomic coordinates in the format
\code{chr:start^end}.

\strong{UMI Counting}

For each cell, each unique combination of splice junction and corrected UMI
(\code{UB}) is counted once. This prevents PCR duplicates from inflating counts.

}

\value{

A sparse matrix in Matrix Market format, where:

\itemize{
\item Rows correspond to splice junctions, labeled by genomic coordinates
(e.g., \code{chr1:10000123^10000456}).
\item Columns correspond to cell barcodes.
\item Entries represent the number of unique UMIs supporting each junction
in each cell.
}

}

\author{
Yang Liao and Wei Shi
}

\examples{

\dontrun{
junction_matrix <- bam2junctions(
  BAMfile = "sample.bam",
  cell.barcodes = c("AAACCTGAG", "AAACCTGTC")
)
}

}