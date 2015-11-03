
debiasplotter <- function(inputFile, subdir, suffix, minLen, maxLen, maxIso=1, minCount=1000, minAveCount=15, smoothing=0.5)
{
	con  <- file(inputFile, open = "r")
	while(length(oneLine <- readLines(con, n=1, warn=FALSE)) > 0) 
	{
		myVector <- (strsplit(oneLine, " "))
		gene <- myVector[[1]][1]		
		bins <- as.numeric(myVector[[1]][4])
		len <- as.numeric(myVector[[1]][5])
		isoforms <- as.numeric(myVector[[1]][6])
		vals <- matrix(0, 0, 2)
		for(i in 1:bins)
		{
			vec <- (strsplit(readLines(con, n=1, warn=FALSE), " "))
			vals <- rbind(vals, as.numeric(vec[[1]]))
		}
		tot <- sum(vals)
		averageCov <- tot/bins		
				
		if(isoforms > maxIso || tot < minCount || averageCov < minAveCount)		# not reliable
			next
			
		if(len < minLen || len > maxLen)
			next
		
		strandSum <- vals[,1] + vals[,2]
		
		coords <- 1:bins * (1/bins)		# x coordinates for the plot
		
		info <- paste("Gene:", gene, "Length:", len, "Total:", tot, "(", isoforms, ")", sep=' ')		
		outputFile <- paste(gene, suffix, "ori.jpg", sep='.')		
		jpeg(paste(subdir, outputFile, sep=''), quality=90)
		plot(coords, strandSum, type="s", xlab="Relative dist. from 3' end", ylab="Read count", main=info)
		dev.off()
		
		quarter <- floor(bins/4)					
		for(i in 1:bins)
		{
			scaler <- 2
			if(i < quarter)		# this is the closest to 3' end
				scaler <- 1.5
			if(i > 3*quarter)	# this is farthest from 3' end
				scaler <- 3
			
			maxi <- max(vals[i,])
			mini <- min(vals[i,])
			if(maxi > averageCov)
				strandSum[i] <- scaler*mini
			else
				strandSum[i] <- scaler*maxi
		}
		
		info <- paste("Gene:", gene, "Length:", len, "Total:", sum(strandSum), "(", isoforms, ")", sep=' ')
		outputFile <- paste(gene, suffix, "raw.jpg", sep='.')		
		jpeg(paste(subdir, outputFile, sep=''), quality=90)
		plot(coords, strandSum, type="s", xlab="Relative dist. from 3' end", ylab="Read count", main=info)
		dev.off()
	}
	close(con)
}

getLoessFit <- function(inputFile, minLen, maxLen, doPlot, maxIso=1, minCount=100, minAveCount=5, smoothing=0.5)
{		
	con <- file(inputFile, open = "r")		
	aggregative <- matrix(c(1.0, 0.0), 1, 2)	
		
	count <- 0
	while(length(oneLine <- readLines(con, n=1, warn=FALSE)) > 0) 
	{
		myVector <- (strsplit(oneLine, " "))
		gene <- myVector[[1]][1]		
		bins <- as.numeric(myVector[[1]][4])
		len <- as.numeric(myVector[[1]][5])
		isoforms <- as.numeric(myVector[[1]][6])		
		tvals <- matrix(0, 0, 2)
		for(i in 1:bins)
		{
			vec <- (strsplit(readLines(con, n=1, warn=FALSE), " "))
			tvals <- rbind(tvals, as.numeric(vec[[1]]))
		}
		vals <- tvals[,1] + tvals[,2]
		tot <- sum(vals)
		averageCov <- tot/bins		
		
		if(isoforms > maxIso)		
			next
		
		if(tot < minCount || averageCov < minAveCount)		
			next
									
		coords <- 1:bins * (1/bins)		# x coordinates for the plot					
		if(len >= minLen && len <= maxLen)		
		{			
			smoothed <- lowess(y=vals, x=coords, f=smoothing)		
			meany <- mean(smoothed$y)
			if(meany > 0.001)
			{
				count <- count + 1
				smoothed$y <- smoothed$y / meany			# normalize to have a mean of 1						
				aggregative <- rbind(aggregative, cbind(matrix(smoothed$x, bins, 1), matrix(smoothed$y, bins, 1)))																				
			}
		}	
	}	
	close(con)

	x <- aggregative[,1]
	y <- aggregative[,2]	
	loessFit <- loess(y ~ x)		
		
	x1 <- 5:995
	x1 <- x1/1000.0
	y1 <- predict(loessFit, x1)
		
	if(doPlot)
	{
		s <- paste(minLen, maxLen, sep="-")	
		s2 <- paste(inputFile, s, sep="_")
		jpeg(paste(s2, "jpg", sep="."), quality=90)		
		plot(x, y, pch='.', col="gray", xlab="Relative dist. from 3' end", ylab="Relative read count", ylim=c(-1.0, 3.0), main=paste("Genes of length:", s, "(", count, ")", sep=" "))			
		lines(x1, y1, col="red", lwd=2)		
		dev.off()			
		
	}	
	return(loessFit)
}

