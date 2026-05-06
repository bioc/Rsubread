\name{bam2junctions}
\alias{bam2junctions}
\title{Generate junction-level UMI counts in single-cell RNA-seq data}
\description{
\code{bam2junctions} generates a UMI count matrix of splice junctions across cells from a single-cell RNA-seq BAM file.
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
    Supported aligners include \code{Cell Ranger} and \code{STARsolo}.
  }

  \item{cell.barcodes}{
    Character vector of target cell barcodes. Only reads from these specified cells are included in the output.
  }
}



\details{
\code{bam2junctions} constructs a UMI count matrix of splice junctions from a single-cell RNA-seq BAM file.

Only reads satisfying both of the following conditions are used:
\itemize{
  \item The final cell barcode is stored in the \code{CB} tag and matches one of the provided \code{cell.barcodes}.
  \item The final UMI is stored in the \code{UB} tag.
}

Additional read filtering rules are applied depending on the aligner:
\itemize{
  \item \code{STARsolo}: All reads passing the above filtering criteria are used.
  \item \code{Cell Ranger}: Reads not assigned to any gene (as indicated by the \code{GX} tag) are included. For reads assigned to one or more genes, the \code{xf} tag is further evaluated. Reads are retained only if the 1st and 5th bits are set; and the 2nd, 3rd and 6th bits are not set.
}

Splice junctions are extracted directly from the CIGAR strings of the selected reads.

Each unique combination of junction and UMI is counted once per cell to generate UMI-level junction counts.
}


\value{
A sparse matrix in Matrix Market format is returned.

\itemize{
  \item Rows correspond to splice junctions, labeled as genomic coordinates in the format \code{chr:start^end} (e.g., \code{chr1:10000123^10000456}).
  \item Columns correspond to cell barcodes.
  \item Matrix values represent unique UMI counts supporting each junction in each cell.
}
}

\author{Yang Liao and Wei Shi}
