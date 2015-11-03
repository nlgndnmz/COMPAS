
checkCulprit <- function(conOut, exonMatrix, exonCoords, chr, known, transIDs, ignorance, sigLevel)
{
	nExons <- nrow(exonMatrix[[1]])
	significance <- 0.0
	geneID <- ""
	if(length(exonMatrix) != 2)		# this function does not make sense unless there are two samples
		return(list(significance = significance, geneID = geneID))
				
	culprit <- 0
	culpritEnd <- 0		
	for(i in 1:nExons)
	{		
		junA <- 0		
		junB <- 0		
		j <- i+1
		while(j <= nExons)
		{
			junA <- junA + exonMatrix[[1]][i,j]
			junB <- junB + exonMatrix[[2]][i,j]	
			j <- j+1
		}		
		
		if(junA < ignorance[1] || junB < ignorance[2])	# expression is too small, not worth it
			next
			
		maxDiff <- c(0, 0)		
		j <- i+1
		k <- c(0, 0)
		while(j <= nExons)
		{
			jDiff <- abs(exonMatrix[[1]][i,j]/junA - exonMatrix[[2]][i,j]/junB)		
			if(jDiff > maxDiff[1])		
			{			
				maxDiff[2] <- maxDiff[1]
				k[2] <- k[1]
				maxDiff[1] <- jDiff											
				k[1] <- j
			}
			else
			{
				if(jDiff > maxDiff[2])
				{	
					maxDiff[2] <- jDiff
					k[2] <- j
				}
			}									
			j <- j+1
		}				
		
		if(maxDiff[1] > sigLevel && i > culpritEnd)
		{
			if(maxDiff[1] > significance)
				significance <- maxDiff[1]
			
			culprit <- i
			culpritEnd <- max(k)
			
			unresolved <- FALSE
			geneID <- ""
			geneName <- ""
			strand <- ""
			for(ktrans in 1:nrow(known))
			{
				if(known[ktrans, culprit] == 1 && known[ktrans, culpritEnd] == 1)
				{
					candidate <- transIDs[ktrans, "gene"]
					if(geneID == "")
					{
						geneID <- candidate
						geneName <- transIDs[ktrans, "HRname"]
						strand <- transIDs[ktrans, "strand"]
					}
					else
					{
						if(candidate != geneID) 	# multiple genes containing the region, could be fishy
							unresolved <- TRUE				
					}			
				}
			}					
			if(geneID != "" && unresolved == FALSE)
			{
				aPercent <- exonMatrix[[1]][culprit,culpritEnd]/junA
				bPercent <- exonMatrix[[2]][culprit,culpritEnd]/junB				
				s <- sprintf('Gene %s Name %s Chr %s Diff: %.2f St: %d End: %d ( %s ) %.2f %.2f Support: %.2f %.2f Range: %d %d', geneID, geneName, chr, maxDiff[1], exonCoords[culprit, 2], exonCoords[culpritEnd, 1], strand, aPercent, bPercent, junA, junB, exonCoords[1, 1],  exonCoords[nExons, 2])
				writeLines(s, conOut)			
			}		
		}
	}							
	return(list(significance = significance, geneID = geneID))
}

adjustWeights <- function(B, biasWeight, junctionIDs, annotation, exonCoords, loessFit)
{			
	xrange <- length(loessFit)
	numofrows <- nrow(exonCoords) + nrow(junctionIDs)
	A <- B
	for(j in 1:ncol(B))	
	{
		totLength <- 0		
		for(i in 1:nrow(exonCoords))
		{
			if(B[i,j] == 1)			
				totLength <- totLength + (exonCoords[i,2] - exonCoords[i,1]) + 1			
		}		
					
		relativeLength <- c(0, 0)
		exonWeights <- matrix(0, nrow(exonCoords), 2)		
		for(i in 1:nrow(exonCoords))
		{
			if(B[i,j] == 1)			
			{
				relativeLength[1] <- relativeLength[2]
				relativeLength[2] <- relativeLength[1] + (exonCoords[i,2] - exonCoords[i,1]) + 1
					
				if(annotation[j, "strand"] == "-")
				{
					exonWeights[i,1] <- relativeLength[1]/totLength
					exonWeights[i,2] <- relativeLength[2]/totLength
				}
				else
				{
					exonWeights[i,1] <- (totLength-relativeLength[1])/totLength
					exonWeights[i,2] <- (totLength-relativeLength[2])/totLength
				}				
				exonWeights[i, exonWeights[i,]<0.05] <- 0.05		# otherwise predict may fail
				exonWeights[i, exonWeights[i,]>0.95] <- 0.95
				x <- c(loessFit[ceiling(xrange * exonWeights[i,1])], loessFit[ceiling(xrange * exonWeights[i,2])])			
				A[i,j] <- (biasWeight*mean(x) + (1.0-biasWeight))	# this is just to dampen the effect of unbiasing
			}
		}	
		if(nrow(junctionIDs) < 1)	# then nothing to do here
			next						
		for(k in 1:nrow(junctionIDs))
		{
			i <- k + nrow(exonCoords)
			if(B[i,j] == 1)
			{				
				x <- loessFit[ceiling(xrange * exonWeights[junctionIDs[k,1],2])]
				A[i,j] <- (biasWeight*x + (1.0-biasWeight))						
			}
		}		
	}	
	return(A)
}

removeLowJunctions <- function(exonMatrix, numSamples, threshold)
{
	nExons <- nrow(exonMatrix[[1]])
	for(i in 1:nExons)
	{	
		cutoffs <- rep(0, numSamples)
		for(k in 1:numSamples)
		{
			if(i < nExons)	
			{
				tot <- sum(exonMatrix[[k]][i, (i+1):nExons])
				cutoffs[k] <- threshold * tot		
			}
		}
				
		for(j in i:nExons)
		{
			e1 <- i
			e2 <- j			
			if(e1==e2)
				next										
			
			ok2remove <- 0
			for(k in 1:numSamples)
			{
				if(exonMatrix[[k]][e1, e2] < cutoffs[k])	# if the junction is below threshold in all samples
					ok2remove <- ok2remove + 1
			}
			
			if(ok2remove == numSamples)
			{
				for(k in 1:numSamples)
					exonMatrix[[k]][e1, e2] <- 0
			}											
		}
	}	
}

getPaths <- function(nExons, exonMatrix, sampleNo, geneNo, canonicals)
{		
	numPaths <- 0
	paths <- list()	
	
	numEdges <- 0
	edges <- matrix(0, nExons, nExons)
	for(i in 1:nExons)
	{
		for(j in i:nExons)
		{
			if(canonicals[geneNo,i] == FALSE || canonicals[geneNo,j] == FALSE)	# at least one exon does not belong to this gene
				next
		
			if(i == j || !(exonMatrix[[sampleNo]][i,j] > 0))
				next		
				
			edges[i,j] <- exonMatrix[[sampleNo]][i,j]
			numEdges <- numEdges + 1
		}
	}
	
	while(numEdges > 0)
	{
		longdist <- matrix(0, nExons, 2)	# all nodes start with a distance of 0
		for(i in 1:nExons)
		{					
			for(j in 1:i)
			{
				if(i == j || edges[j,i] == 0)
					next
				
				if(longdist[i, 1] < longdist[j, 1] + edges[j,i])
				{
					longdist[i, 1] <- longdist[j, 1] + edges[j,i]
					longdist[i, 2] <- j
				}
			}
		}	

		indices <- which(longdist[,1] == max(longdist[,1]))
		nextNode <- indices[1]		
		y <- c()
		while(nextNode > 0)
		{
			y <- c(nextNode, y)
			x <- nextNode
			nextNode <- longdist[x, 2]
			if(nextNode == 0)
				break
			edges[nextNode, x] <- 0		# remove the edge from the graph		
			numEdges <- numEdges - 1
		}		
		if(length(y) > 0)
		{
			numPaths <- numPaths + 1
			paths[[numPaths]] <- y
		}
	}
	return(paths)
}

getRawExp <- function(canonicals, exonMatrix, genes, exonCoords, numSamples)
{
	nExons <- nrow(exonMatrix[[1]])		# there is always at least one sample
	nGenes <- nrow(genes)
	levels <- matrix(0, nrow(genes), numSamples)

	for(geneNo in 1:nGenes)
	{
		uniqExons <- canonicals[geneNo,]
		if(nGenes > 1)
			uniqExons <- 2*canonicals[geneNo,] - colSums(canonicals)		
		
		for(sampleNo in 1:numSamples)
		{
			tot <- 0.0
			len <- 1	# just to avoid diving by zero
			for(i in 1:nExons)
			{
				if(uniqExons[i] < 1)
					next
					
				exonLen <- (exonCoords[i,2] - exonCoords[i,1] + 1)
				tot <- tot + exonLen*exonMatrix[[sampleNo]][i,i]
				len <- len + exonLen 				
			}						
			levels[geneNo, sampleNo] <- tot/len
		}		
	}
	return(levels)
}

getCanonicals <- function(known, transIDs, genes)
{	
	nExons <- ncol(known)
	num <- 0
	curr <- ""
	cans <- matrix(0, nrow(genes), nExons)		
	for(k in 1:nrow(known))
	{
		if(transIDs[k, "gene"] != curr)
		{
			num <- num + 1
			curr <- transIDs[k, "gene"]
		}					
		cans[num,] <- cans[num,] + known[k,]
	}			
	canonicals <- 1*(cans > 0)	
	return(canonicals)
}

getCandidates <- function(exonMatrix, genes, canonicals, sampleNo)
{						
	nExons <- nrow(exonMatrix[[1]])
	A <- matrix(0, 0, nExons)
	for(geneNo in 1:nrow(genes))
	{		
		paths <- getPaths(nExons, exonMatrix, sampleNo, geneNo, canonicals)			
		if(length(paths) == 0)			
			next
			
		candidates <- matrix(0, length(paths), nExons)
		for(k in 1:length(paths))
		{		
			if(k == 1)
				candidates [1,paths[[k]]] <- 1
			else
			{
				x1 <- paths[[k]][1]
				x2 <- paths[[k]][length(paths[[k]])]		

				if(x2 > nExons)				
					warning("Path exceeds number of exons: ", x2, " ", nExons, "\n")									
				
				candidates[k,paths[[k]]] <- 1
				if(candidates[1,x1] == 1)		# otherwise the first exons should stay as zero
					candidates[k,1:x1] <- candidates [1,1:x1]
					
				if(candidates[1,x2] == 1)		# otherwise the last exons should stay as zero
					candidates[k,x2:nExons] <- candidates [1,x2:nExons]									
			}
		}			
		A <- unique(rbind(A, candidates))
	}	
	return(A)	
}

# assign the closest gene for each isoform model
getAnnotation <- function(A, genes, canonicals, offset)
{
	annotation <- matrix("unknown", nrow(A), 4)
	colnames(annotation) <- c("transcript", "gene", "HRname", "strand")
	for(j in 1:nrow(A))
	{
		maxsim <- 0
		winner <- 1
		for(k in 1:nrow(genes))
		{
			x <- sum(canonicals[k,] == A[j,])
			if(x > maxsim)
			{
				maxsim <- x
				winner <- k
			}				
		}			
		trID <- paste(genes[winner, "gene"], j+offset, sep="_")
		annotation[j,] <- c(trID, genes[winner, "gene"], genes[winner, "HRname"], genes[winner, "strand"])				
	}	
	return(annotation)
}

# determine which set of isoforms explain the data best
selectBest <- function(matX, y, B, eps0)
{
	best <- 1
	success <- 0	
	B <- 1*(B>0)
	for(k in 1:ncol(matX))	
	{
		z <- matX[,k]	
		
		covered <- rep(0, length(y))
		for(j in 1:length(z))
		{
			if(z[j] > eps0)		
				covered <- covered + B[,j]
		}
		covered <- y*(covered>0)
		score <- sum(covered)/sum(y)
		if(score > success)
		{	
			success <- score
			best <- k
		}		
	}	
	return(list(best = best, success = success))
}

assignExp <- function(z, minExp, maxEx, eps0, eps1)
{	
	totExp <- 0
	for(i in 1:length(z))
	{
		# if(z[i] > eps0)		
		totExp <- totExp + minExp[i]		
	}
	newZ <- rep(0, length(z))
	if(totExp > eps0)
	{
		for(i in 1:length(z))
		{
			# if(z[i] > eps0 && (minExp[i]/totExp) > eps1)				
			if((minExp[i]/totExp) > eps1)
				newZ[i] <- (minExp[i]/totExp)*maxEx		
		}
	}
	return(newZ)
}

