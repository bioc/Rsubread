\name{featureCounts}
\alias{featureCounts}
\title{featureCounts: a general-purpose read summarization function}
\description{This function assigns mapped sequencing reads to genomic features}
\usage{
featureCounts(files,file.type="SAM",annot=NULL,genome="mm9",isGTFAnnotationFile=FALSE,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,nthreads=1,isStrandSpecific=FALSE,isPairedEnd=FALSE,min.distance=10,max.distance=2000)
}
\arguments{
  \item{files}{ a character vector giving names of SAM format files.}
  \item{file.type}{ a character string giving the type of files provided. \code{SAM} and \code{BAM} formats are supported. \code{SAM} by default.}
  \item{annot}{ a data frame containing annotation information. When \code{annot} is \code{NULL}, an in-built annotation file will be used. \code{NULL} by default.}
  \item{genome}{ a character string specifying the reference genome to which the sequencing reads were mapped. It has the values of \code{mm9}, \code{mm10} and \code{hg}, corresponding to mouse genome build `mm9', mouse genome build `mm10' and human genome build `hg19', respectively. \code{mm9} by default. This argument is only applicable when \code{annot} is \code{NULL}.}
  \item{isGTFAnnotationFile}{ logical indicating if the annotation file provided via \code{annot} is in GTF format. \code{FALSE} by default. If \code{isGTFAnnotationFile} is \code{TRUE}, a GTF format annotation file must be provided via \code{annot} argument.}
  \item{useMetaFeatures}{ logical indicating whether the read summarization should be performed at the feature level (eg. exons) or meta-feature level (eg genes). If \code{TRUE}, features in the annotation (each row is a feature) will be grouped into meta-features using their values in the ``GeneID" column or using the ``gene_id" attribute in the GTF-format annotation file, and reads will assiged to the meta-features instead of the features. See below for more details.}
  \item{allowMultiOverlap}{ logical indicating if a read is allowed to be assigned to more than one feature (or meta-feature) if it is found to overlap with more than one feature (or meta-feature). \code{FALSE} by default.}
  \item{nthreads}{ \code{integer} giving the number of threads used for running this function. \code{1} by default.}
  \item{isStrandSpecific}{ logical indicating whether the strandness should be considered when assigning reads to features. \code{FALSE} by default.}
  \item{isPairedEnd}{ logical indicating if paired-end reads are used. \code{FALSE} by default.}
  \item{min.distance}{ \code{integer} giving the minimal allowable fragment/template length. \code{10} by default.}
  \item{max.distance}{ \code{integer} giving the maximal allowable fragment/template length. \code{2000} by default. \code{min.distance} and \code{max.distance} are only applicable when \code{isPairedEnd} is \code{TRUE}. Note that when a template spans two or more exons, the mapping locations of the pair of reads sequenced from this template may have a distance much larger than the nominal template length.}
  }
\details{
\code{featureCounts} is a general-purpose read summarization function, which assigns to the genomic features (or meta-features) the mapped reads that were generated from genomic DNA and RNA sequencing.

This function takes as input a set of files containing read mapping results output from a read aligner (e.g. \code{\link{align}}), and then assigns the mapped reads to the genomic features.
The acceptable formats of the input files are SAM and BAM.

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

In-built annotation files for a few model organisms are included in the \code{Rsubread} package and these annotations can be easily used with this function to summarize reads for genes and exons.
These annotations were downloaded from the NCBI ftp server (\url{ftp://ftp.ncbi.nlm.nih.gov/genomes/}) and then processed to remove redundant chromosomal regions within each gene.
The resulting annotations contain non-overlapping exon regions for each gene.
Users can use these in-built annotation files simply by setting the parameter \code{genome} to a proper value (\code{mm9}, \code{mm10} or \code{hg}).
The format of in-built annotation files is the same as that of an example format shown below (a data frame provided by users).

This function supports the use of GTF format annotation file.
Users should be provide the name of the annotation file via the \code{annot} argument and set the argument \code{isGTFAnnotationFile} to be \code{TRUE}.
This function uses the `gene_id' attribute to group features to form meta-features if meta-feature level read summarization is needed.

Alternatively, users may create a data frame including the annotation information and then provide that data frame to this function via the \code{annot} parameter for read summarization.
The data frame must have the following format:

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

The \code{GeneID} column may include numbers or character strings.
Note that the chromosomal names included in the \code{Chr} column must match the chromosomal names used in the sequences of the reference genome to which the reads were mapped, otherwise reads cannot be assigned to the features.
 
When paired-end read data are provided and the argument \code{isPairedEnd} is set to \code{TRUE}, the fragments (not reads) will be counted for each feature.
For RNA-seq data, the distance between the mapping locations of the pair of reads sequenced from the same fragment could be much greater than the fragment length when the fragment spans two or more exons.
Such fragments can be counted when a large value of the argument \code{max.distance} is used (\code{2000} by default).
Note that the Subread (\code{\link{align}}) and Subjunc (\code{\link{subjunc}}) aligners can report the mapping coordinates of the paired-end reads even when the two reads have a mapping distance much larger than the nominal fragment length.
 
}
\value{ 
A list with the following components:
\tabular{ll}{
  \code{counts}:\tab a data matrix containing read counts for each feature in each library\cr
  \code{annotation}:\tab a data frame including columns named \code{GeneID}, \code{Chr}, \code{Start} and \code{End}. The \code{GeneID} column includes Entrez gene identifiers when the in-built annotations are used.\cr
  \code{targets}:\tab a character vector giving sample information.
}
}
%\references{
%}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
\examples{
# build an index using a sample reference sequence ('reference.fa') that includes a chromosome called 'chr_dummy' and 90,000 bases
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="./reference_index",reference=ref)

# align a sample read dataset ('reads.txt') to the sample reference
reads <- system.file("extdata","reads.txt",package="Rsubread")
align(index="./reference_index",readfile1=reads,output_file="./Rsubread_alignment.SAM")

# assign mapped reads to dummy features
dummy_features <- data.frame(GeneID=1:5,Chr="chr_dummy",Start=c(1,2000,5000,20000,50000),End=c(1,2000,5000,20000,50000)+1000,stringsAsFactors=FALSE)
fc <- featureCounts("Rsubread_alignment.SAM",feature.type="other",annot=dummy_features)
fc
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
