\name{align}
\alias{align}
\title{Align next-gen sequencing reads to reference genome}
\description{This is a R wrapper function for aligning reads. This function calls the underlying C function}
\usage{
align(index,readfile1,readfile2=NULL,output_file,nsubreads=10,TH1=3,TH2=3,nthreads=1,indels=0,min_distance=50,max_distance=500,PE_orientation="fr")
}
\arguments{
  \item{index}{ character string giving the basename of index file. Index files should be located in the current directory.}
  \item{readfile1}{ character string giving the name of file which includes sequencing reads. This will be the name of first file of paired-end data are provided. The read file should be in FASTQ or FASTA format.}
  \item{readfile2}{ character string giving the name of second file when paired-end read data are provided. \code{NULL} by default.}
  \item{output_file}{ character string giving the name of output file which includes the mapping results.}
  \item{nsubreads}{ numeric value giving the number of subreads extracted from each read. 10 by default.} 
  \item{TH1}{ numeric value giving the consensus threshold for reporting a hit. This threshold will be applied to the first read if paired-end read data are provided. 3 by default.}
  \item{TH2}{ numeric value giving the consensus threhold for the second read in a pair. 3 by default.}
  \item{nthreads}{ numeric value giving the number of threads used for mapping. 1 by default.} 
  \item{indels}{ numeric value giving the number of insertions/deletions allowed during the mapping. 0 by default.}
  \item{min_distance}{ numeric value giving the minimal distance between the pair of reads. 50 by default.}
  \item{max_distance}{ numeric value giving the maximal distance between the pair of reads. 500 by default.}
  \item{PE_orientation}{ character string giving the orientation of the two reads in a pair. "fr" by default, which means the first read is on the forward strand and the second read is on the reverse strand.}
  }
\details{
This function implements a novel mapping strategy which uses a set of 16bp substrings (called subreads) extracted from each read to map them to the reference genome.
Different from "seed-and-extend" mapping strategy, this new strategy does not have an extension step therefore it is a lot faster than the competing aligners.
Due to the short length of the selected subsreads (16bp long), this strategy has a much high sensitivity than other aligners (seed length is usually around 30bp), i.e. it can align a lot more reads than competing aligners.
Our evaluation results (using both simulation dataset and real dataset) showed that the accuracy of the new strategy is comparable to or slightly better than other aligners.

Two key parameters used by this new strategy are the number of subreads selected \code{nsubreads} and the consensus threshold for determining mapping locations \code{TH1} (also \code{TH2} for paired-end read data).
We recommend using the default setting of these parameters to map reads of around 100bp long.
However, users can choose to use more subreads when mapping longer reads.
We recommend to set the value of consensus threshold to be 30 percent of the number of subreads used.

The C implemetation of this strategy can be found at http://sourceforge.net/projects/subread/.
This R function calls the corresponding C function to perform the alignment.
Therefore, it has the mapping speed as the C program.

\code{\link{buildindex}} function should be called if the index has not been built for the reference genome.
The index can be re-used once it has been built.

If paired-end read data is provided, file \code{readfile1} will assumed to contain the first read from the read pair and \code{readfile2} the second read. 
}
\value{ 
A file of SAM format which includes the mapping results.
}
\references{
Yang Liao and Wei Shi, "Subread: a superfast read aligner with high sensitivity and accuracy", In preparation.
}
\author{Wei Shi and Yang Liao}
%\note{}
%\seealso{}
\examples{
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
path <- system.file("extdata",package="Rsubread")
buildindex(basename=file.path(path,"reference_index"),reference=ref)
reads <- system.file("extdata","reads.txt",package="Rsubread")
align(index=file.path(path,"reference_index"),readfile1=reads,output_file=file.path(path,"alignResults.SAM"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
