\name{featureCounts}
\alias{featureCounts}
\title{featureCounts: a general-purpose read summarization function}
\description{This function assigns mapped sequencing reads to genomic features}
\usage{
featureCounts(files,file.type="SAM",annot.inbuilt="mm9",annot.ext=NULL,isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",useMetaFeatures=TRUE,allowMultiOverlap=FALSE,nthreads=1,strandSpecific=0,countMultiMappingReads=FALSE,minMQS=0,isPairedEnd=FALSE,requireBothEndsMapped=FALSE,checkFragLength=FALSE,minFragLength=50,maxFragLength=600,countChimericFragments=TRUE)
}
\arguments{
  \item{files}{ a character vector giving names of SAM/BAM format files.}
  \item{file.type}{ a character string giving the type of files provided. \code{SAM} and \code{BAM} formats are supported. \code{SAM} by default.}
  \item{annot.inbuilt}{ a character string specifying an in-built annotation used for read summarization. It has three possible values including \code{mm9}, \code{mm10} and \code{hg19}, corresponding to the NCBI RefSeq annotations for genomes `mm9', `mm10' and `hg19', respectively. \code{mm9} by default. The in-built annotation has a SAF format (see below).}
  \item{annot.ext}{ A character string giving name of a user-provided annotation file or a data frame including user-provided annotation data. If the annotation is in GTF format, it can only be provided as a file. If it is in SAF format, it can be provided as a file or a data frame. See below for more details about SAF format annotation. \code{annot.ext} will override \code{annot.inbuilt} if they are both provided.}
  \item{isGTFAnnotationFile}{ logical indicating whether the annotation provided via the \code{annot.ext} argument is in GTF format or not. \code{FALSE} by default. This option is only applicable when \code{annot.ext} is not \code{NULL}.}
  \item{GTF.featureType}{ a character string giving the feature type used to select rows in the GTF annotation which will be used for read summarization. \code{exon} by default. This argument is only applicable when \code{isGTFAnnotationFile} is \code{TRUE}.}
  \item{GTF.attrType}{ a character string giving the attribute type in the GTF annotation which will be used to group features (eg. exons) into meta-features (eg. genes). \code{gene_id} by default. This argument is only applicable when \code{isGTFAnnotationFile} is \code{TRUE}.}
  \item{useMetaFeatures}{ logical indicating whether the read summarization should be performed at the feature level (eg. exons) or meta-feature level (eg genes). If \code{TRUE}, features in the annotation (each row is a feature) will be grouped into meta-features using their values in the ``GeneID" column or using the ``gene_id" attribute in the GTF-format annotation file, and reads will assiged to the meta-features instead of the features. See below for more details.}
  \item{allowMultiOverlap}{ logical indicating if a read is allowed to be assigned to more than one feature (or meta-feature) if it is found to overlap with more than one feature (or meta-feature). \code{FALSE} by default.}
  \item{nthreads}{ \code{integer} giving the number of threads used for running this function. \code{1} by default.}
  \item{strandSpecific}{ \code{integer} indicating if strand-specific read counting should be performed. It has three possible values: \code{0} (unstranded), \code{1} (stranded) and \code{2} (reversely stranded). \code{0} by default.}
  \item{countMultiMappingReads}{ logical indicating if multi-mapping reads/fragments should be counted, \code{FALSE} by default. If \code{TRUE}, a multi-mapping read will be counted up to N times if it has N reported mapping locations. This function uses the `NH' tag to find multi-mapping reads.}
  \item{minMQS}{ \code{integer} giving the minimum mapping quality score a read must have so as to be counted. For paired-end reads, at least one end should satisfy this criteria. \code{0} by default.}
  \item{isPairedEnd}{ logical indicating if paired-end reads are used. If \code{TRUE}, fragments (two reads are generated from each fragment in paired-end data) will be counted rather than individual reads. \code{FALSE} by default.}
  \item{requireBothEndsMapped}{ logical indicating if both ends from the same fragment are required to be successfully aligned before the fragment can be assigned to a feature or meta-feature. This parameter is only appliable when \code{isPairedEnd} is \code{TRUE}.}
  \item{checkFragLength}{ logical indicating if the two ends from the same fragment are required to satisify the fragment length criteria before the fragment can be assigned to a feature or meta-feature. This parameter is only appliable when \code{isPairedEnd} is \code{TRUE}. The fragment length criteria are specified via \code{minFragLength} and \code{maxFragLength}.}
  \item{minFragLength}{ \code{integer} giving the minimum fragment length for paired-end reads. \code{50} by default.}
  \item{maxFragLength}{ \code{integer} giving the maximum fragment length for paired-end reads. \code{600} by default. \code{minFragLength} and \code{maxFragLength} are only applicable when \code{isPairedEnd} is \code{TRUE}. Note that when a fragment spans two or more exons, the observed fragment length might be much bigger than the nominal fragment length.}
  \item{countChimericFragments}{ logical indicating whether a chimeric fragment, which has its two reads mapped to different chromosomes, should be counted or not. If this fragment overlaps with only one feature (or one meta-feature), typically by one of its two read, this fragment will be assigned to that feature (or meta-feature). If it is found to overlap more than one feature (or meta-feature), for example one of its two reads overlaps meta-feature A and the other overlaps meta-feature B, and \code{allowMultiOverlap} is \code{FALSE}, then this fragment will not be counted. This parameter is only appliable when \code{isPairedEnd} is \code{TRUE}.}
  }
\details{
\code{featureCounts} is a general-purpose read summarization function, which assigns to the genomic features (or meta-features) the mapped reads that were generated from genomic DNA and RNA sequencing.

This function takes as input a set of files containing read mapping results output from a read aligner (e.g. \code{\link{align}}), and then assigns mapped reads to genomic features.
Both SAM and BAM format input files are accepted.

The argument \code{useMetaFeatures} specifies the read summarization should be performed at the feature level or at the meta-feature level.
Each entry in the annotation data is a feature, which for example could be an exon.
When \code{useMetaFeatures} is \code{TRUE}, the \code{featureCounts} function creates meta-features by grouping features using the gene identifiers included in the ``GeneID" column in the annotation data (or in the ``gene_id" attribute in the GTF format annotation file) and then assigns reads to meta-features instead of features.
The \code{useMetaFeatures} is particularly useful for gene-level expression analysis, because it instructs this function to count reads for genes (meta-features) instead of exons (features).
Note that when meta-features are used in the read summarization, if a read is found to overlap two or more features belong to the same meta-feature it will be only counted once for that meta-feature.

The argument \code{allowMultiOverlap} specifies how those reads, which are found to overlap with more than one feature (or meta-feature), should be assigned.
When \code{allowMultiOverlap} is \code{FALSE}, a read overlapping multiple features (or meta-features) will not be assigned to any of them (not counted).
Otherwise, it will be assigned to all of them.
A read overlaps a meta-feature if it overlaps at least one of the features belonging to this meta-feature.

\code{gene} and \code{exon} are typically used when summarizing RNA-seq read data, which will yield read counts for genes and exons, respectively.

The in-built annotations for human and mouse genomes (\code{hg19}, \code{mm9} and \code{mm10}) provided in this function can be conveniently used for read summarization.
These annotations were downloaded from the NCBI ftp server (\url{ftp://ftp.ncbi.nlm.nih.gov/genomes/}) and were then postprocessed by removing redundant chromosomal regions within each gene and combining adjacent CDS and UTR sequences. 
The in-built annotations use the SAF annotation format (see below) and their content can be retrieved using the \code{\link{getInBuiltAnnotation}} function.

Users may also choose to provide their own annotation for summarization.
If users provide a SAF (Simplified Annotation Format) annotation, the annotation should have the following format:

\preformatted{
GeneID		Chr	Start	End	Strand
497097		chr1	3204563	3207049	-
497097		chr1	3411783	3411982	-
497097		chr1	3660633	3661579	-
100503874	chr1	3637390	3640590	-
100503874	chr1	3648928	3648985	-
100038431	chr1	3670236	3671869	-
...
}

The SAF annotation format has five required columns, including \code{GeneID}, \code{Chr}, \code{Start}, \code{End} and \code{Strand}.
These columns can be in any order.
More columns can be included in the annotation.
Columns are tab-delimited.
Column names are case insensitive.
\code{GeneID} column may contain integers or character strings.
Chromosomal names included in the \code{Chr} column must match those used inclued in the mapping results, otherwise reads will fail to be assigned. 
Users may provide a SAF annotation in the form of a data frame or a file using the \code{annot.ext} argument.\cr\cr

Users may also provide a GTF/GFF format annotation via \code{annot.ext} argument.
But GTF/GFF annotation should only be provided as a file, and \code{isGTFAnnotationFile} should be set to \code{TRUE} when such a annotation is provided.
\code{featureCounts} function uses the `gene_id' attribute in a GTF/GFF annotation to group features to form meta-features when performing read summarization at meta-feature level.\cr\cr
 
When \code{isPairedEnd} is \code{TRUE}, fragments (pairs of reads) instead of reads will be counted.
}
\value{ 
A list with the following components:
\tabular{ll}{
  \code{counts}:\tab a data matrix containing read counts for each feature in each library.\cr
  \code{annotation}:\tab a data frame including annotation information such as gene identifiers (\code{GeneID}) and gene/exon length (\code{Length}). The \code{GeneID} column includes Entrez gene identifiers when the in-built annotations are used.\cr
  \code{targets}:\tab a character vector giving sample information.
}
}
%\references{
%}
\author{Wei Shi and Yang Liao}
%\note{}
\seealso{
\code{\link{getInBuiltAnnotation}}}
\examples{\dontrun{
library(Rsubread)

# Summarizing single-end reads using built-in RefSeq annotation for mouse genome mm9:
featureCounts(files="mapping_results_SE.sam",annot.inbuilt="mm9")

# Summarizing single-end reads using a user-provided GTF annotation file:
featureCounts(files="mapping_results_SE.sam",annot.ext="annotation.gtf",isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id")

# Summarizing single-end reads using 5 threads:
featureCounts(files="mapping_results_SE.sam",nthreads=5)

# Summarizing BAM format single-end read data:
featureCounts(files="mapping_results_SE.bam",file.type="BAM")

# Performing strand-specific read counting (strandSpecific=2 if reversely stranded):
featureCounts(files="mapping_results_SE.sam",strandSpecific=1)

# Summarizing paired-end reads and counting fragments (instead of reads):
featureCounts(files="mapping_results_PE.sam",isPairedEnd=TRUE)

# Counting fragments satisfying the fragment length criteria, eg. [50bp, 600bp]:
featureCounts(files="mapping_results_PE.sam",isPairedEnd=TRUE,checkFragLength=TRUE,minFragLength=50,maxFragLength=600)

# Counting fragments that have both ends successfully aligned without checking the fragment length:
featureCounts(files="mapping_results_PE.sam",isPairedEnd=TRUE,requireBothEndsMapped=TRUE)

# Excluding chimeric fragments from fragment counting:
featureCounts(files="mapping_results_PE.sam",isPairedEnd=TRUE,countChimericFragments=FALSE)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
