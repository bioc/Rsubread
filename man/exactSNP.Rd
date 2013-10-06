\name{exactSNP}
\alias{exactSNP}
\alias{callSNPs}
\title{exactSNP - an accurate and efficient SNP caller}
\description{Measure background noises and perform Fisher's Exact tests to detect SNPs.}
\usage{
exactSNP(readFile,isBAM=FALSE,refGenomeFile,outputFile,qvalueCutoff=12,minAllelicFraction=0,minAllelicBases=1,minReads=1,maxReads=3000, minBaseQuality=13,nTrimmedBases=3,nthreads=1)
callSNPs(readFile,isBAM=FALSE,refGenomeFile,outputFile,qvalueCutoff=12,minAllelicFraction=0,minAllelicBases=1,minReads=1,maxReads=3000, minBaseQuality=13,nTrimmedBases=3,nthreads=1)
}
\arguments{
  \item{readFile}{ a character string giving the name of a file including read mapping results. This function takes as input a SAM file by default. If a BAM file is provided, the \code{isBAM} argument should be set to \code{TRUE}. }
  \item{isBAM}{ logical indicating if the file provided via \code{readFile} is a BAM file. \code{FALSE} by default.}
  \item{refGenomeFile}{ a character string giving the name of a file containing reference sequences (FASTA format).}
  \item{outputFile}{ a character string giving the name of the output file including the discovered SNPs (VCF format).}
  \item{qvalueCutoff}{ a numeric value giving the q-value cutoff for SNP calling at sequencing depth of 50X. \code{12} by default. The q-value is calcuated as -log10(p), where p is the p-value yielded from the Fisher's Exact test. Note that this function automatically adjusts the q-value cutoff for each chromosomal location according to its sequencing depth, based on this cutoff.}
  \item{minAllelicFraction}{ a numeric value giving the minimum fraction of allelic bases out of all read bases included at a chromosomal location required for SNP calling. Its value must be within \code{0} and \code{1}. \code{0} by default.}
  \item{minAllelicBases}{ a numeric value giving the minimum number of allelic (mis-matched) bases a SNP must have at a chromosomal location. \code{1} by default.}
  \item{minReads}{ a numeric value giving the minimum number of mapped reads a SNP-containing location must have (ie. the minimum coverage). \code{1} by default.}
  \item{maxReads}{ Specify the maximum number of mapped reads a SNP-containing location can have. \code{3000} by default. Any location having more than this threshold number of reads will not be considered for SNP calling. This option is useful for removing PCR artefacts.}
  \item{minBaseQuality}{ a numeric value giving the minimum base quality score (Phred score) read bases should satisfy before being used for SNP calling. \code{13} by default(corresponding to base calling p value of 0.05). Read bases with quality scores less than 13 will be excluded from analysis.}
  \item{nTrimmedBases}{ a numeric value giving the number of bases trimmed off from each end of the read. \code{3} by default.}
  \item{nthreads}{ a numeric value giving the number of threads/CPUs used. \code{1} by default.}
}
\details{
This function takes as input a SAM/BAM format file, measures local background noise for each chromosomal location and then performs Fisher's exact tests to find statistically significant SNPs .

This function implements a novel algorithm for discovering SNPs.
This algorithm is comparable with or better than existing SNP callers, but it is fast more efficient.
It can be used to call SNPs for individual samples (ie. no control samples are required).
Detail of the algorithm is described in a manuscript which is currently under preparation.
}
\value{
A VCF format file including detailed information for reported SNPs, such as chromosomal location, reference base, alternative base, read coverage, allele frequency, p value etc. 
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
