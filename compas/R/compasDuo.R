
compasDuo <- function(conOut, conGtf, genes, chr, counts, pairings, known, transIDs, readLength, numSamples, loessFit, enableLoess, ignorance, sigLevel, minJunction=0.3, eps0=0.001, eps1=0.1)
{					
	nExons <- 1
	for(i in 1:nrow(counts))
	{
		if(counts[i,2] > nExons)
			nExons <- counts[i,2]
	}	
	
	exonMatrix <- list()
	for(i in 1:numSamples)
		exonMatrix[[i]] <- matrix(0, nExons, nExons)		
			
	exonCoords <- matrix(0, nExons, 2)
	exonLengths <- matrix(0, nExons, 1)
	exSt <- 3+numSamples
	exEnd <- 4+numSamples	
	
	for(i in 1:nrow(counts))
	{
		for(j in 1:numSamples)		
			exonMatrix[[j]][counts[i,1], counts[i,2]] <- counts[i,2+j]											
			
		if(counts[i,1] == counts[i,2])
		{
			exonCoords[counts[i,1], 1] <- counts[i, exSt]
			exonCoords[counts[i,1], 2] <- counts[i, exEnd]
			exonLengths[counts[i,1]] <- counts[i, exEnd] - counts[i, exSt] + 1
		}
	}	
		
	cc <- checkCulprit(conOut, exonMatrix, exonCoords, chr, known, transIDs, ignorance, sigLevel)
	removeLowJunctions(exonMatrix, numSamples, minJunction)		
	
	canonicals <- getCanonicals(known, transIDs, genes)		
	levels <- getRawExp(canonicals, exonMatrix, genes, exonCoords, numSamples)
	
	BB <- list()
	AA <- list()
	yy <- list()	
	anno <- list()
	maxEx <- list()
	minExp <- list()
	
	offset <- 0
	for(i in 1:numSamples)
	{						
		A <- getCandidates(exonMatrix, genes, canonicals, i)		
		if(nrow(A) == 0)			
			return(list(outcome="none", expLevels=levels))
			
		annotation <- getAnnotation(A, genes, canonicals, offset)
		offset <- nrow(annotation)
		
		y <- diag(exonMatrix[[i]])		# exon counts
		B <- t(A)		# here each column becomes an isoform
		usedJunctions <- matrix(FALSE, nExons, nExons)		
		for(j in 1:ncol(B))
		{			
			prevExon <- 0
			for(k in 1:nExons)		# now we will mark all junctions that are in use
			{
				if(B[k,j] == 1)	
				{
					if(prevExon > 0)
						usedJunctions[prevExon, k] <- TRUE
					prevExon <- k
				}
			}					
		}
				
		nJunctions <- 0		
		junctionIDs <- matrix(0, sum(usedJunctions[,] == TRUE), 2)		
		for(e1 in 1:nExons)
		{			
			for(e2 in e1:nExons)
			{	
				if(usedJunctions[e1, e2] == FALSE)	# this includes the case e1==e2
					next
					
				juncRow <- rep(0, ncol(B))
				indices <- which(B[e1,] + B[e2,] > 1)
				for(j in 1:length(indices))		# there should at least be one
				{
					x <- indices[j]
					if(sum(B[(e1+1):e2,x]) == 1)		# then all in between exons are zero					
						juncRow[x] <- 1					
				}
				B <- rbind(B, juncRow)
				y <- c(y, exonMatrix[[i]][e1, e2])

				nJunctions <- nJunctions + 1
				junctionIDs[nJunctions, ] <- c(e1, e2)
			}
		}			
		
		anno[[i]] <- annotation				
		maxEx[[i]] <- max(exonMatrix[[i]])
		AA[[i]] <- t(A)		
		
		z <- rep(0, ncol(B))		
		for(j in 1:ncol(B))
		{			
			minJunc <- maxEx[[i]]
			for(k in (nExons+1):(nExons+nJunctions))
			{
				if(B[k,j] > 0 && y[k] < minJunc)					
					minJunc <- y[k]				
			}
			z[j] <- minJunc
		}
		
		if(enableLoess)
			B <- adjustWeights(B, 0.2, junctionIDs, annotation, exonCoords, loessFit[[i]])							
						
		BB[[i]] <- B		
		yy[[i]] <- y/maxEx[[i]]						
		minExp[[i]] <- z
		
	}
	
	B <- BB[[1]]		
	y <- yy[[1]]		
		
	sharedTrans <- matrix(0, 0, 2)
	if(numSamples == 2)		# then adjust a few things
	{
		B1 <- rbind(BB[[1]], matrix(0, length(yy[[2]]), ncol(BB[[1]])))
		B2 <- rbind(matrix(0, length(yy[[1]]), ncol(BB[[2]])), BB[[2]])
		B <- cbind(B1, B2)
				
		for(i in 1:ncol(AA[[1]]))
		{
			trans1 <- AA[[1]][,i] 						
			for(j in 1:ncol(AA[[2]]))
			{
				trans2 <- AA[[2]][,j]
				if(identical(trans1, trans2))
				{
					B <- cbind(B, rbind(BB[[1]][,i,drop=FALSE], BB[[2]][,j,drop=FALSE]))					
					anno[[2]][j,] <- anno[[1]][i,]					
					sharedTrans <- rbind(sharedTrans, c(i, j))
				}
			}
		}	
		y <- c(yy[[1]], yy[[2]])
	}
		
	tauVec <- max(t(B) %*% y)
	result <- gpsrBasic(y, B, ncol(B), tauVec)	
	if(result$solFound == FALSE)
		return(list(outcome="failed", expLevels=levels))	# solution failed for whatever reason
		
	matX <- result$solMatrix
	if(sum(matX > eps0) < 1)
		return(list(outcome="none", expLevels=levels))
		
	sb <- selectBest(matX, y, B, eps0)			
	x <- matX[, sb$best]			
	
	z <- list()	
	z[[1]] <- x
	if(numSamples == 2)
	{				
		z[[1]] <- x[1:ncol(BB[[1]])]
		offset <- ncol(BB[[1]])
		z[[2]] <- x[(1+offset):(offset+ncol(BB[[2]]))]				
		offset <- offset + ncol(BB[[2]])		
		
		k <- 1
		while(k <= nrow(sharedTrans))
		{
			z[[1]][sharedTrans[k,1]] <- z[[1]][sharedTrans[k,1]] + x[offset+k]
			z[[2]][sharedTrans[k,2]] <- z[[2]][sharedTrans[k,2]] + x[offset+k]
			k <- k + 1
		}		
	}	
	
	for(i in 1:numSamples)
	{		
		newZ <- assignExp(z[[i]], minExp[[i]], maxEx[[i]], eps0, eps1)
		printGTF(conGtf[[i]], genes, chr, t(AA[[i]]), newZ, exonCoords, anno[[i]], sb$success)
	}

	return(list(outcome="success", expLevels=levels))
}

