\name{callSNPs}
\alias{callSNPs}
\title{Finding Single Nucleotide Polymorphisms (SNPs) from next-gen sequencing data}
\description{Fisher's exact tests are performed to call SNPs.}
\usage{
callSNPs(SAMfile,refGenomeFile,outputFile,minBaseQuality=13,minReadCoverage=5,minAlleleFraction=0.5)
}
\arguments{
  \item{SAMfile}{ a character string giving the name of a SAM format file.}
  \item{refGenomeFile}{ a character string giving the name of the reference sequence file (FASTA format).}
  \item{outputFile}{ a character string giving the name of the output file including called SNPs (BED format).}
  \item{minBaseQuality}{ a numeric value giving the base quality cutoff for selecting the bases included in SNP calling. 13 by default(corresponding to base calling p value of 0.05). Bases with quality scores less than 13 will be excluded from reads when using them to discover SNPs.}
  \item{minReadCoverage}{ a numeric value giving the minimal number of reads required for considering a chromosomal location to include a SNP. 5 by default.}
  \item{minAlleleFraction}{ a numeric value giving the minimal required fraction of reads supporting a called SNP. 0.5 by default.}
  }
\details{
This function takes as input a SAM format file, which includes mapping results for a set of reads, and then calls SNPs by performing Fisher's exact tests.
Bases of low sequencing quality or low mapping quality were removed so as to use high confidence bases to call SNPs at each chromosomal loci.
}
\value{
A BED format file including detailed information about called SNPs, such as chromosomal location, reference base, alternative base, read coverage, allele frequency, p value etc. 
}
%\references{
%}
\author{Yang Liao and Wei Shi}
%\note{}
%\seealso{}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
