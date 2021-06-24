#!/usr/bin/env perl

##-------------------------------------------------------------------
# Author: Nilgun Donmez
# Date: October 29, 2015
#
# This script processes the read mappings in .sam format and 
# converts them to an internal format to be used by the R scripts.
#
# Input:  *.sam (read mappings, read from standard input)
#         *.gtf (Gene annotation file that has been used for the mappings)
#
# Output: *.counts (Internal file to be used by runCOMPAS.R)
#         *.hist (Internal file containing statistics about read coverage)
#         *.skipped (Internal file containing the IDs of overlapping genes)
#
# Part of COMPAS package. See README.md for more information.
##-------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;

sub usage
{
	print STDERR "Usage: cat <tophat_results.sam> | processSam -g <gene_annotations.gtf> -o <output_prefix> -r <read_length> -s [0|1|2]\n\n";
	print STDERR "Set the -s option to 1 for the forward strand genes and 2 for the reverse strand genes in a strand-specific dataset; set it to 0 for non-strand-specific datasets.\n";
	print STDERR "For strand-specific datasets, it is assumed that the SECOND read in the pair matches the strand of the gene. Otherwise, use the toggle '-t' option along with the commands above.\n";
	print STDERR "This script normally expects the SAM file to be sorted by genomic coordinates. For unsorted SAM files, use the '-u' option.\n\n";
	print STDERR "Part of COMPAS package. See README.md for more information.\n\n";
	exit(0);
}
if($#ARGV < 7) { usage(); }	

my $gtfFile = "";
my $prefix = "compas_output";
my $readlen = 100;
my $ssRNAseq = 0;
my $hlp = 0;
my $tgl = 0;
my $unsorted = 0;

my $result = GetOptions("genes|g=s" => \$gtfFile,
						"output|o=s" => \$prefix,
						"read|r=i" => \$readlen,						
						"strand|s=i" => \$ssRNAseq,
						"toggle|t" => \$tgl,
						"unsorted|u" => \$unsorted,
						"help|h" => \$hlp);		

if($hlp or !$result) { usage(); }		

my $sortedSam = ($unsorted == 0) ? 1 : 0;

if($ssRNAseq == 0)	
{	$tgl = 0; }		# this should not have been set!

open(GF, "<$gtfFile") || die "Can not open $gtfFile\n";

my $mIntron = int(0.9 * $readlen);
my $hist = int(0.4 * $readlen);
my $mOverlap = 5;
my $step = 500;		# window size for the bins

my $numGenes = 1;
my %geneHash;		# hashes the gene names to an integer
my @geneInfo;		# holds the chromosome, strand and ensemble id for each gene
my @transcripts;	# holds the transcript id for each gene
my @geneModels;		# holds the start and end of all exons
my %isoforms;		# holds the known isoform models

my $geneNo = 1;
my $strand = "+";

my %geneCoords;		# this is to figure out which genes have overlapping exons

my %exons;
my %introns;
my %bins;

my $geneID;
my $transID;

while(<GF>)
{
	my $ln = $_;	
	if($ln =~ /gene_id \"(.+?)\";/)
	{
		$geneID = $1;		
		if($ln =~ /transcript_id \"(.+?)\";/)		# this should work for both Ensembl38.77 and earlier versions
		{	$transID = $1;	}
		else
		{	next;	}

		my @parts = split(/\t/, $ln);
		unless($parts[2] eq "exon")
		{	next; }
		
		$strand = $parts[6];		
		if(($ssRNAseq eq 2 and $strand eq "+") || ($ssRNAseq eq 1 and $strand eq "-"))		# decide whether to skip this gene
		{	next;	}

		my $chr = $parts[0];		
		if(!defined $bins{$chr})
		{
			@{$bins{$chr}} = ();
			@{$exons{$chr}} = ();
			@{$introns{$chr}} = ();
		}
		
		my $geneName = "NOTAVAILABLE";
		if($ln =~ /gene_name \"(.+?)\";/)
		{
			$geneName = $1;
		}
		
		if(defined $geneHash{$geneID})
		{
			$geneNo = $geneHash{$geneID};
		}
		else
		{
			$geneNo = $numGenes;
			$numGenes++;
			$geneHash{$geneID} = $geneNo;

			push @{$geneInfo[$geneNo]}, $chr;
			push @{$geneInfo[$geneNo]}, $strand;		
			push @{$geneInfo[$geneNo]}, $geneID;
			push @{$geneInfo[$geneNo]}, $parts[1];		# RNA type (protein_coding, lincRNA, etc)		
			push @{$geneInfo[$geneNo]}, $geneName;
		}

		my @t1 = ($parts[3], $parts[4] + 1);			# make the ends exclusive, we will change that later
		push @{$geneModels[$geneNo]}, \@t1;		
				
		my @t2 = ($parts[3], $parts[4] + 1, $geneNo);
		push @{$geneCoords{$chr}}, \@t2;

		unless(defined $isoforms{$transID})
		{	push @{$transcripts[$geneNo]}, $transID; }

		my @t3 = ($parts[3], $parts[4]);
		push @{$isoforms{$transID}}, \@t3;
	}
	else
	{
		# don't do anything
	}
}
close(GF);

print STDERR "Finished reading the GTF file ($numGenes genes)\n";

my %coupled;
my %hostgenes;
while(my ($key, $values) = each %geneCoords)
{
	my @sorted = sort {@{$a}[0] <=> @{$b}[0] || @{$a}[1] <=> @{$b}[1]} @{$values};	# sorted exon coords
	
	for(my $i=0; $i<=$#sorted; $i++)
	{
		my $st = $sorted[$i]->[0];
		my $end = $sorted[$i]->[1];
		my $g1 = $sorted[$i]->[2];				
			
		for(my $j=$i+1; $j<=$#sorted; $j++)
		{
			if($sorted[$j]->[0] >= $end)		# since the end is exclusive, we allow an equality
			{	last; }
			
			my $g2 = $sorted[$j]->[2];
			if($g1 ne $g2)		# the exons overlap but the genes are not the same
			{				
				if(defined $coupled{$g1} and defined $coupled{$g2})
				{
					if($coupled{$g1} != $coupled{$g2})		# otherwise nothing to do
					{
						my $g3 = $coupled{$g2};
						foreach my $k (@{$hostgenes{$g3}})
						{
							$coupled{$k} = $coupled{$g1};
						}
						$coupled{$g3} = $coupled{$g1};
						push @{$hostgenes{$coupled{$g1}}}, $g3;
						push @{$hostgenes{$coupled{$g1}}}, @{$hostgenes{$g3}};						
						delete $hostgenes{$g3};
					}					
				}
				elsif(defined $coupled{$g1})
				{
					$coupled{$g2} = $coupled{$g1};
					push @{$hostgenes{$coupled{$g1}}}, $g2;				
				}
				elsif(defined $coupled{$g2})
				{
					$coupled{$g1} = $coupled{$g2};
					push @{$hostgenes{$coupled{$g2}}}, $g1;			
				}
				else	# neither is defined so pick one as the host
				{
					$coupled{$g1} = $g1;
					$coupled{$g2} = $g1;					
					push @{$hostgenes{$g1}}, $g2;			
				}
			}
		}
	}		
}

my $skipFile = $prefix.".skipped";
open(SF, ">$skipFile") || die "Can not open $skipFile\n";
while(my ($key, $values) = each %hostgenes)
{		
	print SF $geneInfo[$key][0], " ", $geneInfo[$key][1], " ", $geneInfo[$key][2];
	foreach my $val (@{$values})
	{
		print SF " ", $geneInfo[$val][1], " ", $geneInfo[$val][2];
	}
	print SF "\n";		
}
close(SF);

my @exonSt;
my @exonEnd;
my @exonCovPos;		# this is actually for cDNA bias estimation and should only be used for single isoform genes
my @exonCovNeg;		# this is actually for cDNA bias estimation and should only be used for single isoform genes

my @intronSt;
my @intronEnd;
my @intronCov;

my @geneMatrices;
$#geneMatrices = $numGenes;

for(my $g=1; $g<$numGenes; $g++)
{
	my $chr = $geneInfo[$g][0];
	my $strand = $geneInfo[$g][1];
		
	my @sorted;		
	if(defined $coupled{$g}) 
	{	
		if($coupled{$g} != $g)		# gene is coupled but it is not the host
		{	next; }
		
		my @temparr = @{$geneModels[$g]};
		foreach my $g2 (@{$hostgenes{$g}})
		{		
			@temparr = (@temparr, @{$geneModels[$g2]});
		}		
		@sorted = sort {@{$a}[0] <=> @{$b}[0] || @{$a}[1] <=> @{$b}[1]} (@temparr);
	}
	else
	{
		@sorted = sort {@{$a}[0] <=> @{$b}[0] || @{$a}[1] <=> @{$b}[1]} @{$geneModels[$g]};
	}
	
	my %arrhash;
	my $stlabel = 0;
	my $endlabel = 1;
	
	for(my $i=0; $i<=$#sorted; $i++)
	{
		my $st = @{$sorted[$i]}[0];
		my $end = @{$sorted[$i]}[1];
		
		$arrhash{$st} = $stlabel;
		unless(defined $arrhash{$end} and $arrhash{$end} eq $stlabel)
		{
			$arrhash{$end} = $endlabel;			
		}
		
		for(my $j=0; $j<=$#sorted; $j++)
		{
			if($i == $j)
			{	next; }
			
			if(@{$sorted[$j]}[1] > $st and @{$sorted[$j]}[1] < $end)	# this exon's end is overlapping with our current exon
			{
				$arrhash{ @{$sorted[$j]}[1] } = $stlabel;				
			}	
			
			if(@{$sorted[$j]}[0] >= $end)	# no need to go further
			{	last; }
		}
	}
	
	my $pos = 0;
	my @arr;			
	foreach my $k (sort {$a <=> $b} keys %arrhash)
	{
		$arr[$pos][0] = $k;
		$arr[$pos][1] = $arrhash{$k};
		$pos++;		
	}
	
	my $prevB = -2;
	my $numExons = 0;
	my $totGeneLen = 0;

	for(my $i=0; $i<$pos-1; $i++)
	{
		if($arr[$i][1] eq $endlabel and $arr[$i+1][1] eq $stlabel)			# this is an intron
		{	next; }
		elsif($arr[$i][1] eq $endlabel and $arr[$i+1][1] eq $endlabel)		# this shouldn't happen
		{
			for(my $j=0; $j<$pos; $j++)
			{	
				print STDERR $arr[$j][0], " ", $arr[$j][1], "\n";
			}			
			die "Fatal error: exon boundaries do not make sense\n"; 
		}
		else
		{
			my $a = $arr[$i][0];
			my $b = $arr[$i+1][0] - 1;		# make the end inclusive, so that the exons will never overlap					

			push @{$exonSt[$g]}, $a;
			push @{$exonEnd[$g]}, $b;
			
			$totGeneLen += ($b - $a + 1);

			my $st = int($a/$step);
			my $end = int($b/$step);

			for(my $j=$st; $j<=$end; $j++)
			{
				my @t = ($g, $numExons);
				push @{$bins{$chr}->[$j]}, \@t;
			}			
			$numExons += 1;
		}
	}
	
	if(!defined $coupled{$g} and $totGeneLen > 2*$readlen)		
	{		
		my $num = 0;	
		my $leftover = 0;
		
		if($strand eq "-")	# gene is on the reverse strand of the reference	
		{
			for(my $i=0; $i<$numExons; $i++)	# this will already start from 3' end
			{		
				my $a = $exonSt[$g][$i];
				my $b = $exonEnd[$g][$i];
				
				my $pos = $a + ($hist - $leftover);
				while($pos < $b)
				{
					my $ind = int($pos/$step);					
					
					my @t = ($g, $pos, $num);
					push @{$exons{$chr}->[$ind]}, \@t;
					
					$exonCovPos[$g][$num] = 0;
					$exonCovNeg[$g][$num] = 0;
					
					$num += 1;
					$pos += $hist;					
				}
				$leftover = $pos - $b;
			}
		}
		else	# gene is on the forward strand of the reference
		{
			for(my $i=$numExons-1; $i>=0; $i--)		# this will start from the 3' end
			{		
				my $a = $exonEnd[$g][$i];
				my $b = $exonSt[$g][$i];
				
				my $pos = $a - ($hist - $leftover);
				while($pos > $b)
				{
					my $ind = int($pos/$step);

					my @t = ($g, $pos, $num);
					push @{$exons{$chr}->[$ind]}, \@t;
					
					$exonCovPos[$g][$num] = 0;
					$exonCovNeg[$g][$num] = 0;
					
					$num += 1;
					$pos -= $hist;
				}
				$leftover = $b - $pos;
			}		
		}
	}

	my $numIntrons = 0;
	for(my $i=0; $i<$numExons-1; $i++)		# fill the introns
	{
		if($exonSt[$g]->[$i+1] > (2 + $exonEnd[$g]->[$i]))		# use only long enough introns
		{			
			my $a = (1 + $exonEnd[$g]->[$i]);
			my $intEnd = ($exonSt[$g]->[$i+1] - 1);
			
			my $b = $a + $mIntron;
			
			while($a < $intEnd - 1)
			{
				if($b > $intEnd)
				{	$b = $intEnd; }
				
				push @{$intronSt[$g]}, $a;
				push @{$intronEnd[$g]}, $b;

				my $st = int($a/$step);
				my $end = int($b/$step);

				for(my $j=$st; $j<=$end; $j++)
				{
					my @t = ($g, $numIntrons);
					push @{$introns{$chr}->[$j]}, \@t;
				}
				$intronCov[$g][$numIntrons] = 0;	# this is to initialize the array
				$numIntrons += 1;
			
				$a = $b+1;
				$b = $a + $mIntron;
			}
		}				
	}
	
	for(my $i=0; $i<$numExons; $i++)		
	{
		for(my $j=0; $j<$numExons; $j++)
		{
			$geneMatrices[$g][$i][$j] = 0; 	# this is to initialize the 3D array
		}
	}
}

print STDERR "Finished processing the gene models\n";
my $countFile = $prefix.".counts";
my $histFile = $prefix.".hist";

open(CF, ">$countFile") || die "Can not open $countFile\n";	
open(HF, ">$histFile") || die "Can not open $histFile\n";

sub printCounts
{	
	my $prevchr = $_[0];		
	if($sortedSam && !defined $bins{$prevchr})
	{	return; }
		
	for(my $g=1; $g<$numGenes; $g++)
	{
		my $chr = $geneInfo[$g][0];
		my $strand = $geneInfo[$g][1];
		my $key = $geneInfo[$g][2];
		my $type = $geneInfo[$g][3];
	
		if($sortedSam && $chr ne $prevchr)
		{	next; }				
		
		if(defined $coupled{$g} and $coupled{$g} != $g)
		{	next; }
		
		my $numEx = scalar(@{$exonSt[$g]});

		# first decide whether there is any expressed exon
		my $ok = 0;
		for(my $j=0; $j<$numEx; $j++)
		{
			if( $geneMatrices[$g][$j][$j] > 0)
			{
				$ok = 1;
				last;
			}
		}

		unless($ok)			# do not even bother reporting this gene
		{	next;	}

		print CF "CHR $chr GENE $key STRAND $strand TYPE $type\n";

		my @before;
		my @after;	

		for(my $j=0; $j<$numEx; $j++)		# this is a sneak peek to figure out first and last exons
		{		
			for(my $k=$j+1; $k<$numEx; $k++)
			{
				if($geneMatrices[$g][$j][$k] > 0)
				{
					$before[$k] += 1;		# exon j comes before exon k
					$after[$j] += 1;		# exon k comes after exon j
				}
			}
		}

		for(my $j=0; $j<$numEx; $j++)
		{
			if($geneMatrices[$g][$j][$j] eq 0)		# zero exon coverage
			{
				print CF $j+1, " ", $j+1, " 0 ", $exonSt[$g][$j], " ", $exonEnd[$g][$j], "\n";
				next;
			}

			my $cov = $geneMatrices[$g][$j][$j] / ($exonEnd[$g][$j] - $exonSt[$g][$j] +1);
			print CF $j+1, " ", $j+1, " ", $cov, " ", $exonSt[$g][$j], " ", $exonEnd[$g][$j], "\n";

			for(my $k=$j+1; $k<$numEx; $k++)		# by design it is a square matrix
			{
				if($geneMatrices[$g][$j][$k] > 0)
				{
					print CF $j+1, " ", $k+1, " ", $geneMatrices[$g][$j][$k], " ", $exonEnd[$g][$j], " ", $exonSt[$g][$k], "\n";
				}
			}
		}

		print CF "INTRONS $key\n";
		if(defined $intronCov[$g])
		{
			my $numInt = scalar(@{$intronCov[$g]});			
			for(my $j=0; $j<$numInt; $j++)
			{
				if($intronCov[$g][$j] > 0)
				{
					my $cov = $intronCov[$g][$j] / ($intronEnd[$g][$j] - $intronSt[$g][$j] + 1);
					print CF $j+1, " ", $cov, " ", $intronSt[$g][$j] , " ", $intronEnd[$g][$j], "\n";
				}
			}
		}

		my $totGeneLen = 0;
		print CF "PAIRED $key\n";
		for(my $j=0; $j<$numEx; $j++)
		{
			$totGeneLen += ($exonEnd[$g][$j] - $exonSt[$g][$j] + 1);
			for(my $k=$j-1; $k>=0; $k--)
			{
				if($geneMatrices[$g][$j][$k] > 0)
				{
					print CF $k+1, " ", $j+1, " ", $geneMatrices[$g][$j][$k], "\n";
				}
			}
		}		
		
		my @genes;
		push @genes, $g;
		if(defined $coupled{$g})		# no need to check if it is a host because if it is not, we would not arrive here
		{			
			push @genes, @{$hostgenes{$g}};			
		}
	
		print CF "ISOFORMS";		
		foreach my $g2 (@genes)
		{		
			my $numIso = scalar(@{$transcripts[$g2]});											
			print CF " ", $geneInfo[$g2][2], " ", $geneInfo[$g2][4], " ", $geneInfo[$g2][1], " ", $numIso;
		}
		print CF "\n";
		
		foreach my $g2 (@genes)
		{		
			my $numIso = scalar(@{$transcripts[$g2]});														
			for(my $i=0; $i<$numIso; $i++)
			{
				my $transID = $transcripts[$g2][$i];
				my @sorted = sort {@{$a}[0] <=> @{$b}[0]} @{$isoforms{$transID}};	# exons in the same isoform should never be overlapping
				my $k = 0;
				foreach my $item (@sorted)
				{
					while($k<$numEx)
					{
						if($exonSt[$g][$k] > @{$item}[1])	# proceed to the next exon in isoform
						{	last; }								
						
						if(@{$item}[0] <= $exonSt[$g][$k] and $exonEnd[$g][$k] <= @{$item}[1])		# the exon contains the partial exon
						{	print CF "1 "; }
						else
						{	print CF "0 "; }

						$k += 1;
					}
				}
				while($k<$numEx)
				{
					print CF "0 ";
					$k += 1;
				}
				print CF "$transID\n";
			}
		}								
	
		if(!defined $coupled{$g} and $totGeneLen > 2*$readlen)	# only report histograms for uncoupled and long enough genes
		{
			my $numBins = scalar(@{$exonCovPos[$g]});
			my $numIso = scalar(@{$transcripts[$g]});	
			print HF "$key $chr $strand $numBins $totGeneLen ", $numIso, "\n";

			for(my $i=0; $i<$numBins; $i++)
			{
				print HF $exonCovPos[$g][$i], " ", $exonCovNeg[$g][$i], "\n";
			}
		}
	}	
	print STDERR "Finished writing the counts for chromosome $prevchr\n";
}

my $prevChromosome = "-1";
my $readCounter = 0;
my $badReads = 0;
my @lenDist;
$#lenDist = $readlen;
for(my $i=0; $i<=$readlen; $i++)
{	$lenDist[$i] = 0; }

if($tgl && $ssRNAseq != 0)		# toggle the desired read orientation in strand-specific
{
	$ssRNAseq = ($ssRNAseq == 1) ? 2 : 1;
}

my $milestone = 1000000;

LOOP1:
while(<STDIN>)
{
	my $s = $_;
	if($s =~ /^@/)	# get rid of the headers
	{	next; }

	my @parts = split(' ', $s, 12);		# this behaves as awk

	my $chr = $parts[2];
	#if($chr =~ /^chr/)   
	#{	$chr = substr($chr, 3); }
	
	my $flg = $parts[1];	
	
	if($ssRNAseq eq 2)
	{	
		unless( (($flg & 32) and ($flg & 64)) || (($flg & 16) and ($flg & 128)) )
		{	next;	}
	}
	if($ssRNAseq eq 1)
	{
		unless( (($flg & 16) and ($flg &64)) || (($flg & 32) and ($flg & 128)) )
		{	next;	}
	}	
	
	if($sortedSam && $chr ne $prevChromosome)
	{
		printCounts($prevChromosome);		
		$prevChromosome = $chr;
	}
	
	if(!defined $bins{$chr})	# could be one of those additional chromosomes
	{	next; }

	my $st = $parts[3];
	my $cigar = $parts[5];
	if(($cigar =~ /[\*SHPX=]/))
	{	next; }
	
	my $refStrand = 1;
	if($flg & 16)		# then it means the sequence is reverse complemented; meaning it maps to the reverse strand
	{	$refStrand = 0; }

	my $gene1 = undef;
	my $exon1 = undef;
	my $gene2 = undef;
	my $exon2 = undef;
	my $ok = 1;
	
	$readCounter += 1;
	if(!$sortedSam && $readCounter > $milestone)
	{
		my $datestring = localtime();
		print STDERR "Processed $milestone reads $datestring\n";
		$milestone = 2 * $milestone;
	}
	
	my $nomLength = length($parts[9]);	
	if($nomLength > $readlen)
	{	die("Read exceeds given length!"); }	
	$lenDist[$nomLength] += 1;

	if($parts[6] eq "=")	# the other end has also been mapped
	{
		my $st2 = $parts[7];
		my $num = int($st2/$step);
		foreach my $item (@{$bins{$chr}->[$num]})
		{
			my $gene = @{$item}[0];
			my $exon = @{$item}[1];

			if( $exonEnd[$gene][$exon] >= $st2 and $exonSt[$gene][$exon] <= $st2 )	# check if the start of the read overlaps with the exon
			{
				if(!defined $gene2)
				{
					$gene2 = $gene;
					$exon2 = $exon;
				}			
			}
		}
	}

	my $end = 0;
	my @fields = split(/(M|I|N|D)/, $parts[5]);

	for(my $i=0; $i<$#fields; $i+=2)
	{
		if($fields[$i+1] eq "N")
		{
			$st = $end + $fields[$i] + 1;
		}
		elsif($fields[$i+1] eq "M")
		{
			$end = $st + $fields[$i] - 1;		# the end points are always inclusive!!
			my $num = int($st/$step);

			my $d = 1;
			if(int($end/$step) > $num)	# decide if we need to look at the next bin or not
			{	$d = 2; }				# reads are much shorter than the window size, we at most have to look at two bins

			if(($i==0 || ($i+1)==$#fields) && $fields[$i] < $mOverlap)		# if the overlap is at the beginning or the end of the read, then it's unreliable if it's too short
			{
				$d = -1;	# this is simply to skip the for loop below
			}

			for(my $j=0; $j<$d; $j++)
			{
				$num += $j;

				foreach my $item (@{$exons{$chr}->[$num]})		
				{
					my $g = @{$item}[0];		# gene no
					my $pos = @{$item}[1];		# absolute position of this bin on the chromosome
					my $n = @{$item}[2];		# bin no				

					if($pos >= $st and $pos <= $end)	# check if the read overlaps with the point
					{
						if($refStrand)
						{
							$exonCovPos[$g][$n] += 1;
						}
						else
						{						
							$exonCovNeg[$g][$n] += 1;
						}
					}
				}
								
				foreach my $item (@{$introns{$chr}->[$num]})
				{
					my $gene = @{$item}[0];
					my $intron = @{$item}[1];									

					my $intSt = $intronSt[$gene][$intron];
					my $intEnd = $intronEnd[$gene][$intron];

					if($st < $intEnd and $end > $intSt)		# then the read overlaps with the intron
					{						
						my $intronLen = $intEnd - $intSt + 1;
						my $overlap = (($end-$intSt) < ($intEnd-$st)) ? ($end-$intSt+1) : ($intEnd-$st+1);
						$overlap = ($nomLength < $overlap) ? $nomLength : $overlap;
						$overlap = ($intronLen < $overlap) ? $intronLen : $overlap;
														
						$intronCov[$gene][$intron] += $overlap;		# increment the intron coverage
					}
				}

				foreach my $item (@{$bins{$chr}->[$num]})
				{
					my $gene = @{$item}[0];
					my $exon = @{$item}[1];

					my $exSt = $exonSt[$gene][$exon];
					my $exEnd = $exonEnd[$gene][$exon];

					unless($exEnd >= $st and $exSt <= $end)	# check if the read overlaps with the exon
					{	next; }
					
					if(defined $gene1 and $gene1 ne $gene)		# this should not happen unless there is a fusion or readthrough
					{
						$badReads += 1;
						next LOOP1;				# terminate the processing of the read
					}
					
					if(defined $gene2 and $gene2 eq $gene)
					{
						if($exon > $exon2)
						{	$geneMatrices[$gene][$exon][$exon2] += 1; }	# fill the lower triangle of the matrix
						elsif($exon < $exon2)
						{	$geneMatrices[$gene][$exon2][$exon] += 1; }	# if exons are the same, don't do anything
					}

					my $overlap = 1;
					my $exonlen = $exEnd - $exSt + 1;
					if($exonlen > 1)
					{
						$overlap = (($end-$exSt) < ($exEnd-$st)) ? ($end-$exSt+1) : ($exEnd-$st+1);
						$overlap = ($nomLength < $overlap) ? $nomLength : $overlap;
						$overlap = ($exonlen < $overlap) ? $exonlen : $overlap;
					}										

					$geneMatrices[$gene][$exon][$exon] += $overlap;		# increment the exon coverage
					if(defined $exon1 and $exon ne $exon1)
					{
						($exon > $exon1) || die "($d , $num) Bad exon order at read $parts[0] on $gene: $exon is not greater than $exon1. Start: $st End: $end\n";
						$geneMatrices[$gene][$exon1][$exon] += 1;		# increment the junction count by one
					}			
					$gene1 = $gene;
					$exon1 = $exon;
				}
			}
		}
		elsif($fields[$i+1] eq "I" || $fields[$i+1] eq "D")		# insertion or deletion in the read
		{
			$st = $end + 1;		# just skip that part
		}
		else
		{	die "Cigar string fail\n"; }	# should not reach here
	}	
}

printCounts($prevChromosome);		# the last chromosome (or if the SAM file is unsorted, everything)

close(CF);
close(HF);

print STDERR "Total reads mapped: $readCounter \n";
print STDERR "Number of bad reads: $badReads \n\n";

for(my $i=0; $i<=$readlen; $i++)
{	print STDERR $i, " ", $lenDist[$i], "\n"; }

