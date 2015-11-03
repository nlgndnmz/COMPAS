\name{compasCore}
\alias{compasCore}
\title{
Core method of the compas package
}
\description{
This is the core function of the compas package, an R package to detect alternative splicing using RNA-Seq data. The best way to use this (or any other) function in the compas package is through the R script that is provided at https://github.com/nlgndnmz/COMPAS.git .
}
\usage{
compasCore(numSamples, numReads, readLength, inputFile,
	outputPrefix, baseReads, histFiles, enableLoess, 
	ignoranceCutoff, sigLevel)
}
\arguments{
  \item{numSamples}{
Number of samples to be compared. Currently only 1 or 2 samples are allowed.
}
  \item{numReads}{
A vector of size \code{numSamples} containing the number of reads per each sample.
}
  \item{readLength}{
A vector of size \code{readLength} containing the maximum read length per each sample.
}
  \item{inputFile}{
The path to an input file containing the coverage values for exons and junctions in each gene.
}
  \item{outputPrefix}{
A prefix for output files to be written on disk.
}
  \item{baseReads}{
This is the base number of reads used to normalize read coverage in the given samples.
}
  \item{histFiles}{
A vector of size \code{numSamples} containing the paths for optional read histogram files per each sample. If provided, these files are used to estimate the 5' to 3' bias of the datasets.
}
  \item{enableLoess}{
A logical parameter that enables the Loess regression on the histogram files given in \code{histFiles}.
}
  \item{ignoranceCutoff}{
A threshold value that is used to decide if the coverage at a junction is sufficient to make calls.
}
  \item{sigLevel}{
A threshold value that is used to decide if a differential splicing event is significant enough to report. 
}
}
\value{
Does not return a value.
}
\author{
Nilgun Donmez
}
\note{
Part of compas package. More information can be found at:

https://github.com/nlgndnmz/COMPAS.git
}
\seealso{
See also \code{\link{gpsrBasic}}
}
