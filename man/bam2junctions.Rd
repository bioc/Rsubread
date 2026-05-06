\name{bam2junctions}
\alias{bam2junctions}
\title{Generate junction-level UMI counts in single-cell RNA-seq data}

\description{
\code{bam2junctions} constructs a UMI count matrix of splice junctions across cells
from a single-cell RNA-seq BAM file. Splice junctions are inferred directly from
read alignments, and UMI deduplication is performed at the junction level per cell.
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
The BAM file must include cell barcode and UMI annotations in standard tags.
Supported aligners include Cell Ranger and STARsolo.
}

\item{cell.barcodes}{
Character vector of cell barcodes to include in the analysis.
Only reads assigned to these barcodes are retained in the output.
}

}

\details{

\strong{Overview}

This function extracts splice junctions from aligned reads and constructs
a sparse matrix of UMI counts, where each row corresponds to a junction
and each column corresponds to a cell.

\strong{BAM Tags Used}

The function relies on the following BAM tags:

\itemize{

\item \strong{CB tag (Cell Barcode):}
The \code{CB} tag stores the corrected cell barcode assigned to a read.
This barcode typically represents the droplet or cell of origin after
error correction by the aligner (e.g., Cell Ranger or STARsolo).

\item \strong{UB tag (UMI Barcode):}
The \code{UB} tag contains the corrected Unique Molecular Identifier (UMI),
which is used to identify unique transcript molecules and remove PCR duplicates.

\item \strong{GX tag (Gene Annotation):}
The \code{GX} tag indicates the gene(s) to which a read is assigned.
Reads without a \code{GX} tag or with ambiguous assignments may be treated
differently depending on the aligner.

\item \strong{xf tag (Feature Flags, Cell Ranger-specific):}
The \code{xf} tag is a bitwise flag used by Cell Ranger to annotate read properties.
It encodes information such as mapping confidence, exon/intron assignment,
and transcript compatibility. Specific bits are used for filtering (see below).
}

\strong{Read Filtering}

Only reads satisfying all of the following conditions are used:

\itemize{

\item The read contains a valid cell barcode in the \code{CB} tag that matches
one of the provided \code{cell.barcodes}.

\item The read contains a valid UMI in the \code{UB} tag.
}

Additional filtering depends on the aligner:

\itemize{

\item \strong{STARsolo:}
All reads passing the above criteria are retained without additional filtering.

\item \strong{Cell Ranger:}
\itemize{
  \item Reads without gene assignment (no \code{GX} tag) are included.
  \item For reads assigned to one or more genes, the \code{xf} tag is evaluated.
  Reads are retained only if:
  \itemize{
    \item Bit 1 and bit 5 are set (indicating confident transcript/exon assignment),
    \item Bits 2, 3, and 6 are not set (filtering out ambiguous or low-confidence reads).
  }
}
}

\strong{Junction Extraction}

Splice junctions are identified directly from the CIGAR string of each read,
specifically from \code{N} (skipped region) operations representing introns.

Each junction is recorded as genomic coordinates in the format:
\code{chr:start^end}.

\strong{UMI Counting}

For each cell, unique combinations of splice junction and UMI are counted once.
This ensures that PCR duplicates do not inflate junction counts.

}

\value{

A sparse matrix in Matrix Market format, where:

\itemize{
\item Rows correspond to splice junctions, labeled by genomic coordinates
(e.g., \code{chr1:10000123^10000456}).
\item Columns correspond to cell barcodes.
\item Matrix entries represent the number of unique UMIs supporting each
junction in each cell.
}

}

\author{
Yang Liao and Wei Shi
}

\seealso{
\code{\link{readMM}} for reading Matrix Market files.
}

\examples{

\dontrun{
# Example usage
junction_matrix <- bam2junctions(
  BAMfile = "sample.bam",
  cell.barcodes = c("AAACCTGAG", "AAACCTGTC")
)
}

}
