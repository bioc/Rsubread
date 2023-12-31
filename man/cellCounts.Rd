\name{cellCounts}
\alias{cellCounts}
\title{Map and quantify single-cell RNA-seq data generated by 10X Genomics}
\description{Process raw 10X scRNA-seq data and generate UMI counts for each gene in each cell.}
\usage{
cellCounts(

    # input data
    index,
    sample,
    input.mode = "BCL",
    cell.barcode = NULL,
  
    # parameters used for mapping reads
    nsubreads = 15,
    minVotes = 1,
    maxMismatches = 10,
    minMappedLength = 1,
  
    # parameters used for assigning and counting UMIs
    annot.inbuilt = "mm39",
    annot.ext = NULL,
    isGTFAnnotationFile = FALSE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    useMetaFeatures = TRUE,
    
    # user provided UMI cutoff for cell calling
    umi.cutoff = NULL,

    # number of threads
    nthreads = 10,

    # parameters for alignment of multi-mapping reads
    nBestLocations = 1,
    uniqueMapping = FALSE,

    # report UMI counts for cell barcodes that were filtered out
    reportExcludedBarcodes = FALSE
    )
}
\arguments{
  \item{index}{ A character string providing the base name of index files that were created for a reference genome by using \code{\link{buildindex}} function.}
  \item{sample}{ A data frame or a character string providing sample-related information. The value provided to this parameter depends on the \code{input.mode} parameter. See Details section below for more details.}
  \item{input.mode}{ A character string specifying the input mode. The supported input modes include \code{BCL}, \code{FASTQ}, \code{FASTQ-dir} and \code{BAM}. The \code{BCL} mode includes \code{BCL} and \code{CBCL} formats, which are used by Illumina for storing the raw reads directly generated from their sequencers. When \code{BCL} mode is specified, \code{cellCounts} will automatically identify whether the input data are in \code{BCL} format or \code{CBCL} format. \code{FASTQ} is the FASTQ format of sequencing reads. \code{FASTQ-dir} is a directory where FASTQ-format reads are saved. \code{FASTQ-dir} is useful for providing \code{cellCounts} the FASTQ data generated by \code{bcl2fastq} program or \code{bamtofastq} program (developed by 10X). \code{BAM} is the BAM format of mapped read data with cell barcodes and UMI sequences included (See the Details section below for more details). \code{BCL} by default.}
  \item{cell.barcode}{ A character string giving the name of a text file (can be gzipped) that contains the set of cell barcodes used in sample preparation. If \code{NULL}, a cell barcode set will be determined for the input data by \code{cellCounts} based on the matching of cell barcodes sequences of the first 100,000 reads in the data with the three cell barcode sets used by 10X Genomics. \code{NULL} by default.}
  \item{nsubreads}{ A numeric value giving the number of subreads (seeds) extracted from each read. \code{15} by default.}
  \item{minVotes}{ A numeric value giving the minimum number of votes required for reporting a hit. \code{1} by default.}
  \item{maxMismatches}{ A numeric value giving the maximum number of mis-matched bases allowed in the mapping of a read. \code{10} by default. Mis-matches found in soft-clipped bases are not counted.}
  \item{minMappedLength}{ A numeric value giving the minimum number of mapped bases in a read required for reporting a hit. \code{1} by default.}
  \item{annot.inbuilt}{ Specify an inbuilt annotation for UMI counting. See \code{\link{featureCounts}} for more details. \code{mm39} by default.}
  \item{annot.ext}{ Specify an external annotation for UMI counting. See \code{\link{featureCounts}} for more details. \code{NULL} by default.}
  \item{isGTFAnnotationFile}{ See \code{\link{featureCounts}} for more details. \code{FALSE} by default.}
  \item{GTF.featureType}{ See \code{\link{featureCounts}} for more details. \code{exon} by default.}
  \item{GTF.attrType}{ See \code{\link{featureCounts}} for more details. \code{gene_id} by default.} 
  \item{useMetaFeatures}{ Specify if UMI counting should be carried out at the meta-feature level (eg. gene level). See \code{\link{featureCounts}} for more details. \code{TRUE} by default.}
  \item{umi.cutoff}{ Specify a UMI count cutoff for cell calling. All the cells with a total UMI count greater than this cutoff will be called. If \code{NULL}, a bootstrapping procedure will be performed to determine this cutoff. \code{NULL} by default. }
  \item{nthreads}{ A numeric value giving the number of threads used for read mapping and counting. \code{10} by default.}
  \item{nBestLocations}{A numeric value giving the maximum number of reported alignments for each multi-mapping read. \code{1} by default.}
  \item{uniqueMapping}{A logical value specifying if the multi-mapping reads should not be reported as mapped (i.e. reporting uniquely mapped reads only). \code{FALSE} by default.}
  \item{reportExcludedBarcodes}{If \code{TRUE}, report UMI counts for those cell barcodes that were filtered out during cell calling. \code{FALSE} by default.}
}

\details{
The \code{cellCounts} function takes as input scRNA-seq reads generated by the 10X platform, maps them to the reference genome and then produces UMI (Unique Molecular Identifier) counts for each gene in each cell.
It also performs sample demultiplexing, cell barcode demultiplexing and read deduplication, before generating UMI counts.
It can process multiple datasets at the same time.

The sample information can be provided to \code{cellCounts} via the \code{sample} parameter.
If the input format is \code{BCL} or \code{CBCL}, the provided sample information should include the location where the read data are stored, flowcell lanes used for sequencing, sample names and names of index sets used for indexing samples. 
If a sample was sequenced in all lanes, then its lane number can be set as \code{*}.
Alternatively, all the lane numbers can be listed for this sample. 
The sample information should be saved to a \code{data.frame} object and then provided to the \code{sample} parameter.
Below shows an example of this data frame:


\preformatted{
InputDirectory		Lane		SampleName	IndexSetName
/path/to/dataset1	1		Sample1		SI-GA-E1
/path/to/dataset1	1		Sample2		SI-GA-E2
/path/to/dataset1	2		Sample1		SI-GA-E1
/path/to/dataset1	2		Sample2		SI-GA-E2
/path/to/dataset2	1		Sample3		SI-GA-E3
/path/to/dataset2	1		Sample4		SI-GA-E4
/path/to/dataset2	2		Sample3		SI-GA-E3
/path/to/dataset2	2		Sample4		SI-GA-E4
...
}


It is compulsory to have the four column headers shown in the example above when generating this data frame for a 10X dataset.
If more than one datasets are provided for analysis, the \code{InputDirectory} column should include more than one distinct directory.
Note that this data frame is different from the Sample Sheet generated by the Illumina sequencer.
The \code{cellCounts} function uses the index set names included in this data frame to generate an Illumina Sample Sheet and then uses it to demultiplex all the samples.  

If the input format is \code{FASTQ}, a \code{data.frame} object containing the following three columns, \code{BarcodeUMIFile}, \code{ReadFile} and \code{SampleName}, should be provided to the \code{sample} parameter.
Each row in the data frame represents a sample.
The \code{BarcodeUMIFile} column should include names of FASTQ files containing cell barcode and UMI sequences for each read, and the \code{ReadFile} column should include names of FASTQ files containing genomic sequences for corresponding reads.

If the input format is \code{FASTQ-dir}, a character string, which includes the path to the directory where the FASTQ-format read data are stored, should be provided to the \code{sample} parameter.
The data in this directory are expected to be generated by the \code{bcl2fastq} program (developed by Illumina), or by the \code{cellranger mkfastq} command (developed by 10x), or by the \code{bamtofastq} program (developed by 10x).

Finally, if the input format is \code{BAM}, a \code{data.frame} object containing the following two columns, \code{BAMFile} and \code{SampleName}, should be provided to the \code{sample} parameter.
Each row in the data frame represents a sample.
The \code{BAMFile} column should include names of BAM files containing genomic sequences of reads with associated cell barcode and UMI sequences. 
The cell barcode and UMI sequences should be provided in the "CR" and "UR" tags included in the BAM files. 
The Phred-encoded quality strings of the cell barcode and UMI sequences should be provided in the "CY" and "UY" tags. 
These tags were originally defined in CellRanger BAM files.
}

\value{
The \code{cellCounts} function returns a \code{List} object to R.
It also outputs one BAM file for each sample. The BAM file includes location-sorted read mapping results.
If the input mode is neither \code{FASTQ} nor \code{FASTQ-dir}, it also outputs three or four gzipped FASTQ files, which include cell barcode and UMI sequences (R1), sample index sequences (I1, and optionally I2 for dual-indexed samples) and the actual genomic sequences of the reads (R2), respectively. 

The returned \code{List} object contains the following components:
\item{counts}{a \code{List} object including UMI counts for each sample. Each component in this object is a matrix that contains UMI counts for a sample. Rows in the matrix are genes and columns are cells.}

\item{annotation}{a \code{data.frame} object containing a gene annotation. This is the annotation that was used for the assignment of UMIs to genes during quantification. Rows in the annotation are genes. Columns of the annotation include \code{GeneID}, \code{Chr}, \code{Start}, \code{End} and \code{Length}.}

\item{sample.info}{a \code{data.frame} object containing sample information and quantification statistics.
It includes the following columns: \code{SampleName}, \code{InputDirectory} (if the input format is \code{BCL}), \code{TotalCells}, \code{HighConfidenceCells} (if \code{umi.cutoff} is \code{NULL}), \code{RescuedCells} (if \code{umi.cutoff} is \code{NULL}), \code{TotalUMI}, \code{MinUMI}, \code{MedianUMI}, \code{MaxUMI}, \code{MeanUMI}, \code{TotalReads}, \code{MappedReads} and \code{AssignedReads}.
Each row in the data frame is a sample.}

\item{cell.confidence}{a \code{List} object indicating if a cell is a high-confidence cell or a rescued cell (low confidence). Each component in the object is a logical vector indicating which cells in a sample are high-confidence cells. \code{cell.confidence} is included in the output only if \code{umi.cutoff} is \code{NULL}.}

\item{counts.excluded.barcodes}{a \code{List} object including UMI counts for excluded cell barcodes for each sample.}
}

\author{Yang Liao and Wei Shi}

\seealso{
\code{\link{buildindex}}, \code{\link{align}}, \code{\link{featureCounts}}
}

\examples{
# Data used in this example (71MB) can be downloaded from https://shilab-bioinformatics.github.io/cellCounts-Example/cellCounts-Example.zip
\dontrun{
library(Rsubread)
buildindex("chr1", "hg38_chr1.fa.gz")
sample.sheet <- data.frame(BarcodeUMIFile="reads_R1.fastq.gz", ReadFile="reads_R2.fastq.gz", SampleName="Example")
counts <- cellCounts("chr1", sample.sheet, input.mode="FASTQ", annot.inbuilt="hg38")
}
}
