\name{callSNPs}
\alias{callSNPs}
\title{Finding Single Nucleotide Polymorphisms (SNPs) and INDELs (insertions and deletions) from next-gen sequencing data}
\description{Use read coverage and allele frequency to call SNPs and INDELs}
\usage{
callSNPs(SAMfile,readLength,refGenomeFile,outputFile,minReadCoverage=5,minAlleleFraction=0.5)
}
\arguments{
  \item{SAMfile}{ a character string giving the name of a SAM format file.}
  \item{readLength}{ a numeric value giving the read length.}
  \item{refGenomeFile}{ a character string giving the name of a FASTA format file which contains reference sequences.}
  \item{outputFile}{ a character string giving the base name of the output files. The file which contains discovered SNPs will have a suffix '.SNP' and the file which contains discovered INDELS will have suffix '.INDEL'.}
  \item{minReadCoverage}{ a numeric value giving the minimal number of reads which cover a chromosomal location, when considering if it contains a SNP or INDEL. 5 by default.}
  \item{minAlleleFraction}{ a numeric value giving the minimal fraction of reads which supports the called SNP or INDEL, out of all the reads which cover that chromosomal location. 0.5 by default.}
  }
\details{
This function takes as input a SAM format file, which includes mapping results for a set of reads, and then calls SNPs and INDELs using read coverage and allele frequency for each chromosomal location.\\

It requires the chromosome names included in the SAM file to be 'chr1', 'chr2', 'chr3', ..., 'chrX', and 'chrY'.\\ 

To report a SNP or INDEL, a chromosomal location must be covered by at least \code{minReadCoverage} reads and the fraction of reads, which include bases different from the reference base at that location, must be equal to or greater than \code{minAlleleFraction}.

}
\value{
Two files are generated: one including called SNPs(named \code{outputFile}.SNP), and the other including disovered insertions and deletions (named \code{outputFile}.INDEL).
These files include detailed SNP or INDEL information such as chromosomal locations, read coverages, allele frequencies and so on. 
}
%\references{
%}
\author{Wei Shi and Jenny Zhiyin Dai}
%\note{}
%\seealso{}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
