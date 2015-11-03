
printGTF <- function(outcon, genes, chr, transcripts, abundances, exons, annotation, success, cutoff=0.0001)
{	
	reported <- 0
	for(i in 1:nrow(transcripts))
	{
		if(abundances[i] < cutoff)
			next
		
		novel <- "no"
		trID <- annotation[i, "transcript"]
		geneID <- annotation[i, "gene"]
		geneName <- annotation[i, "HRname"]
		direction <- annotation[i, "strand"]								
			
		exonBegin <- exons[transcripts[i,]>0,1]
		exonEnd <- exons[transcripts[i,]>0,2]
			
		numExons <- sum(transcripts[i,]>0)
		exBegin <- exonBegin[1]
		exEnd <- exonEnd[1]
		normCov <- abundances[i]  
		s <- sprintf('%s\tCompreS\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; exon_number "1"; gene_name "%s"; coverage "%.2f"; novel "%s"; success "%.2f";', chr, exBegin, exEnd, direction, geneID, trID, geneName, normCov, novel, success)
		exonum <- 1
		if(numExons > 1)
		{
			for(j in 2:numExons)
			{
				if(exonBegin[j] == (1 + exonEnd[j-1]))	# then merge with the old one
					exEnd <- exonEnd[j]
				else
				{
					writeLines(s, outcon)	# flush the string we already have
					exonum <- exonum + 1
					exBegin <- exonBegin[j]
					exEnd <- exonEnd[j]
				}
				s <- sprintf('%s\tCompreS\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s"; exon_number "%d"; gene_name "%s"; coverage "%.2f"; novel "%s"; success "%.2f";', chr, exBegin, exEnd, direction, geneID, trID, exonum, geneName, normCov, novel, success)
			}
		}
		writeLines(s, outcon)	
		reported <- reported + 1
	}
	return(reported)
}
