
compasCore <- function(numSamples, numReads, readLength, inputFile, outputPrefix, baseReads, histFiles, enableLoess, ignoranceCutoff, sigLevel)
{
	conGtf <- list()
	
	loessFit <- list()
	xrange <- 1:5000 * (1/5000)		
	plotLoess <- FALSE
	maskNoise <- TRUE
	
	if(numSamples == 1)
	{	
		plotLoess <- TRUE	
		maskNoise <- FALSE
	}
		
	cat("Processing: ", outputPrefix, "\n")		
	baseline <- (1000/readLength)*(baseReads/numReads)		# this is basically like RPKM, but adjusted with readlength
	ignorance <- ignoranceCutoff*baseline
	
	for(i in 1:numSamples)
	{
		conGtf[[i]] <- file(paste(outputPrefix, i, "gtf", sep="."), open = "w")
		s <- sprintf('#!CompreS: v1.0beta\n#!InputFile: %s\n#!Baseline: %.2f\n#!Ignorance: %.2f\n#!SignificanceLevel: %.2f', inputFile, baseline[i], ignorance[i], sigLevel);
		writeLines(s, conGtf[[i]])
		if(enableLoess)
		{				
			lofit <- getLoessFit(histFiles[i], 1000, 5000, plotLoess)		# if numSamples is 1, then it also plots the figures						
			loessFit[[i]] <- predict(lofit, xrange)	
		}
	}

	con <- file(inputFile, open = "r")
	conOut <- file(paste(outputPrefix, "out", sep="."), open = "w")
	conStat <- file(paste(outputPrefix, "stat", sep="."), open = "w")	

	chr <- "-"					# some dummy values
	numIso <- c(1, 1)

	genes <- matrix()
	transIDs <- matrix()
	known <- matrix()
	counts <- matrix()
	pairings <- matrix()

	totIso <- 0
	maxExon <- 0
	smallestIsoform <- 0

	paired <- FALSE
	isoforms <- FALSE
	introns <- FALSE

	maxIntron <- rep(0, numSamples)

	headcount <- 0
	while(length(oneLine <- readLines(con, n=1, warn=FALSE)) > 0) 
	{
		myVector <- strsplit(oneLine, " ")
		
		if(myVector[[1]][1] == "PAIRED")	# paired counts
		{		
			pairings <- matrix(0, 0, 2+numSamples)		
			
			paired <- TRUE
			isoforms <- FALSE
			introns <- FALSE
			next
		}
		
		if(myVector[[1]][1] == "ISOFORMS")	# known transcripts
		{		
			numGenes <- (length(myVector[[1]])-1)/4
			genes <- matrix("-", numGenes, 3)		
			colnames(genes) <- c("gene", "HRname", "strand")		
			numIso <- rep(1, numGenes)
		
			j <- 1
			i <- 2
			totIso <- 0
			while(i < length(myVector[[1]]))
			{
				genes[j, "gene"] <- myVector[[1]][i]
				genes[j, "HRname"] <- myVector[[1]][i+1]
				genes[j, "strand"] <- myVector[[1]][i+2]
				numIso[j] <- as.numeric(myVector[[1]][i+3])
				totIso <- totIso + numIso[j]
				i <- i + 4
				j <- j + 1
			}
		
			rownames(genes) <- genes[,1]
			
			known <- matrix(0, 0, maxExon)
			transIDs <- matrix("-", totIso, 4)
			colnames(transIDs) <- c("transcript", "gene", "HRname", "strand")		
			smallestIsoform <- maxExon
			
			i <- 1
			for(j in 1:length(numIso))
			{
				for(k in 1:numIso[j])
				{
					transIDs[i, "gene"] <- genes[j, "gene"]
					transIDs[i, "HRname"] <- genes[j, "HRname"]				
					transIDs[i, "strand"] <- genes[j, "strand"]
					i <- i + 1
				}
			}		
			
			paired <- FALSE
			isoforms <- TRUE
			introns <- FALSE	
			next
		}
		
		if(myVector[[1]][1] == "INTRONS")
		{
			maxIntron <- rep(0, numSamples)
			allIntrons <- matrix(0, 0, 4)
			
			paired <- FALSE
			isoforms <- FALSE
			introns <- TRUE
			next
		}
		
		if(myVector[[1]][1] == "CHR")	# gene change
		{	
			chr <- myVector[[1]][2]						
			counts <- matrix(0, 0, 4+numSamples)
							
			maxExon <- 0

			paired <- FALSE
			isoforms <- FALSE
			introns <- FALSE

			headcount <- headcount + 1
			next
		}
		
		if(introns)
		{
			vec <- matrix(as.numeric(myVector[[1]]), 1, 4)
			allIntrons <- rbind(allIntrons, vec)		
			
			smpl <- vec[1,1]
			value <- vec[1,2]
			
			if(smpl > numSamples)
				smpl <- 1
			
			if(value > maxIntron[smpl])
				maxIntron[smpl] <- value
			
			next	
		}
		
		if(paired)		
		{	
			vec <- matrix(as.numeric(myVector[[1]]), 1, 2+numSamples)			
			pairings <- rbind(pairings, vec)			
			next
		}
		
		if(isoforms)
		{
			vec <- matrix(as.numeric(myVector[[1]][1:maxExon]), 1, maxExon)			
			known <- rbind(known, vec)					
			transIDs[nrow(known), "transcript"] <- myVector[[1]][maxExon+1]
				
			if(sum(vec > 0) < smallestIsoform)
				smallestIsoform <- sum(vec > 0)		# the number of exons in the smallest isoform
			
			if(nrow(known) < totIso)	# gene information is not yet complete
				next
				
			# it's time to process the gene now!
			if(maxExon == 1)						# a single exon gene
			{
				for(i in 1:numSamples)
				{				
					s <- sprintf('%s\tCompreS\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; exon_number "1"; gene_name "%s"; coverage "%.2f"; novel "no"; success "1.00";', chr, counts[1,3+numSamples], counts[1,4+numSamples], transIDs[1, "strand"], transIDs[1, "gene"], transIDs[1, "transcript"], transIDs[1, "HRname"], counts[1,2+i]*baseline[i])			
					writeLines(s, conGtf[[i]])
					
					s <- sprintf('%s %s %d %s %f PASS', genes[,"gene"], genes[, "HRname"], i, "success", counts[1,2+i]*baseline[i])							
					writeLines(s, conStat)									
				}
			}
			if(maxExon > 1)
			{				
				noisy <- rep("PASS", numSamples)		# first we will decide if noise level dominates the gene coverage
				for(i in 1:numSamples)
				{								
					if(sum(counts[,2+i] > maxIntron[i]) < 0.3*smallestIsoform)
					{
						noisy[i] <- "FAIL"					
						if(maskNoise)
						{
							counts[,2+i] <- 0			# then mask out the necessary vectors
							if(nrow(pairings) > 0)
								pairings[,2+i] <- 0							
							next	# nothing more to do here
						}
					}
					
					if(!maskNoise || maxIntron[i] < 0)		# either no intron is found or the level is too low, doesn't worth masking
						next
					
					myIntrons <- matrix(c(i, 0 ,0 , 1), 1, 4)
					myIntrons <- rbind(myIntrons, allIntrons[allIntrons[,1]==i,])					
					myIntrons <- rbind(myIntrons, c(i, 0, max(counts[,4+numSamples])+1, max(counts[,4+numSamples])+2))		# add a fake intron to the very end
					for(j in 2:nrow(myIntrons))
					{														
						if((myIntrons[j-1,4]+1) == myIntrons[j,3])	# adjacent introns, nothing to do
							next
												
						# this is 80% of the average noise from these introns
						noise <- 0.4*(myIntrons[j-1,2] + myIntrons[j,2])															
						indices <- which(counts[,3+numSamples] > myIntrons[j-1,4] & counts[,4+numSamples] < myIntrons[j,3])		# get the exons that are enclosed by these introns
							
						if(length(indices) > 0)
						{
							for(m in 1:length(indices))
							{
								k <- indices[m]
								if(counts[k,1] == counts[k,2])		# an exon, ok to normalize				
									counts[k,2+i] <- max(0, counts[k,2+i] - noise)
								
								if(counts[k,1] == (counts[k,2]-1) && counts[k,3+numSamples] == (counts[k,4+numSamples]-1))	# a fake junction, ok to normalize
									counts[k,2+i] <- max(0, counts[k,2+i] - noise)
							}
						}
					}				
				}
						
				for(i in 1:numSamples)
				{
					counts[,2+i] <- counts[,2+i]*baseline[i]		# normalize counts		
					if(nrow(pairings) > 0)		# otherwise it means no pairing was found
						pairings[,2+i] <- pairings[,2+i]*baseline[i]
				}							
																		
				res <- compasDuo(conOut, conGtf, genes, chr, counts, pairings, known, transIDs, max(readLength), numSamples, loessFit, enableLoess, ignorance, sigLevel)		
				for(sampleNo in 1:numSamples)
				{					
					s <- sprintf('%s %s %d %s %f %s', genes[,"gene"], genes[, "HRname"], sampleNo, res$outcome, res$expLevels[,sampleNo], noisy[sampleNo])						
					writeLines(s, conStat)	
				}								
			}
			next
		}
		
		vec <- matrix(as.numeric(myVector[[1]]), 1, 4+numSamples)	
		if(vec[1] > maxExon)
			maxExon <- vec[1]
			
		counts <- rbind(counts, vec)	
	}

	for(i in 1:numSamples)
		close(conGtf[[i]])
		
	close(con)
	close(conOut)
	close(conStat)
}

