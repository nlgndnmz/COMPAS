#!/usr/bin/env Rscript

##-------------------------------------------------------------------
# Author: Nilgun Donmez
# Date: July 17, 2013
#
# This is a wrapper script that calls the methods in the "compas" R package.
#
# Input:	.counts file (Output of processSam.pl)
#           .hist file (Output of processSam.pl)
#
# Output:	.gtf (Contains the predicted isoforms)
# 			.stat (Contains statistics about the processed genes)
#			.out (Comparison of predicted isoforms in multi-sample mode)
#
# Part of COMPAS package. See README.md for more information.
##-------------------------------------------------------------------

library("compas")

argline <- commandArgs(TRUE)
myargs <- (strsplit(argline, " "))

if(length(myargs) < 7)
{
	warning("Usage: runCOMPAS.R <#samples> -i <input file> -o <output prefix> -b <baseline reads> -c <ignorance cutoff> -s <SI cutoff> -l[2] <read length> -n[2] <#reads> -h[2] <hist file>")		
	quit()
}

baseReads <- 10*1000*1000

numSamples <- as.numeric(myargs[1][[1]])
numReads <- rep(baseReads, numSamples)
readLength <- rep(100, numSamples)
histFiles <- rep("none", numSamples)

inputFile <- "none"
outputFile <- "none"

ignoranceCutoff <- 50		# ignore an alternative splicing event if it's supported by less reads 
sigLevel <- 0.3		# significance threshold to report an alternative splicing event

enableLoess <- FALSE

cmnd <- 2
while(cmnd < length(myargs))
{		
	if(myargs[cmnd][[1]] == "-i")
		inputFile = myargs[cmnd+1][[1]]
		
	if(myargs[cmnd][[1]] == "-o")
		outputFile = myargs[cmnd+1][[1]]
		
	if(myargs[cmnd][[1]] == "-b")
		baseReads = as.numeric(myargs[cmnd+1][[1]])
		
	if(myargs[cmnd][[1]] == "-c")
		ignoranceCutoff = as.numeric(myargs[cmnd+1][[1]])
		
	if(myargs[cmnd][[1]] == "-s")
		sigLevel = as.numeric(myargs[cmnd+1][[1]])
		
	if(myargs[cmnd][[1]] == "-h")
	{
		histFiles[1] = myargs[cmnd+1][[1]]
		enableLoess <- TRUE
	}
	
	if(myargs[cmnd][[1]] == "-l")
		readLength[1] = as.numeric(myargs[cmnd+1][[1]])
	
	if(myargs[cmnd][[1]] == "-n")
		numReads[1] = as.numeric(myargs[cmnd+1][[1]])
			
	if(myargs[cmnd][[1]] == "-h2" && numSamples > 1)
		histFiles[2] = myargs[cmnd+1][[1]]
		
	if(myargs[cmnd][[1]] == "-l2" && numSamples > 1)
		readLength[2] = as.numeric(myargs[cmnd+1][[1]])
		
	if(myargs[cmnd][[1]] == "-n2" && numSamples > 1)
		numReads[2] = as.numeric(myargs[cmnd+1][[1]])
		
	cmnd <- cmnd + 2
}

compasCore(numSamples, numReads, readLength, inputFile, outputFile, baseReads, histFiles, enableLoess, ignoranceCutoff, sigLevel)

proc.time()
warnings()
