#!/usr/bin/env perl

##-------------------------------------------------------------------
# Author: Nilgun Donmez
# Date: October 29, 2015
#
# This script combines two .counts files belonging to different samples
# and creates a single .counts file to be used by the R scripts in
# multi-sample mode.
#
# Input:  *.counts (counts file belonging to sample 1)
#         *.counts (counts file belonging to sample 2)
#
# Output: *.counts (written to standard output)
#
# Part of COMPAS package. See README.md for more information.
##-------------------------------------------------------------------

use strict;
use warnings;
use Getopt::Long;

sub usage
{
	print STDERR "Usage: mergeCounts.pl -q <sample1.counts> -t <sample2.counts> > <combined.counts>\n";
	print STDERR "Part of COMPAS package. See README.md for more information.\n\n"
	exit(0);
}
if($#ARGV < 3) { usage(); }

my $fname1;
my $fname2;
my $hlp = 0;

my $result = GetOptions("query|q=s" => \$fname1,
						"target|t=s" => \$fname2,						
						"help|h" => \$hlp);		

if($hlp or !$result) { usage(); }	

open(FH1, "<$fname1") || die "can not open $fname1\n";
open(FH2, "<$fname2") || die "can not open $fname2\n";

my @linesA = <FH1>;
my @linesB = <FH2>;

chomp(@linesA);
chomp(@linesB);

close(FH1);
close(FH2);

my $i = 0;
my $j = 0;
my %genes;
while($i<=$#linesA)		# figure out which genes are common
{
	my $ln = $linesA[$i++];
	if($ln =~ /GENE/)
	{	
		my @parts = split(' ', $ln);		
		my @t = ($i, -1);
		$genes{$parts[3]} = \@t; 
	}						
}				
while($j<=$#linesB)
{
	my $ln = $linesB[$j++];
	if($ln =~ /GENE/)
	{			
		my @parts = split(' ', $ln);
		if(defined $genes{$parts[3]})
		{	
			$genes{$parts[3]}->[1] = $j;
		}
	}			
}			

my @exonsA;
my @exonsB;
my @exonBegin;
my @exonEnd;

my @pairsA;
my @pairsB;	
	
my $numEx = 0;	
foreach my $key (keys %genes)
{
	$i = $genes{$key}->[0];
	$j = $genes{$key}->[1];	
	if($j < 0) { next; }
									
	$numEx = 0;	
	print $linesB[$j-1], "\n";	# the first line
	
	# --- EXONS ---	
	my $s = $i;		# first perform a sneak peek
	while($s<=$#linesA)
	{
		my $ln = $linesA[$s++];
		if($ln =~ /INTRONS/)
		{	last; }	
		my @parts = split(' ', $ln);		
		if($parts[0] eq $parts[1])
		{
			$exonBegin[$parts[0]] = $parts[3];
			$exonEnd[$parts[0]] = $parts[4];			
		}		
		if($parts[0] > $numEx)
		{	$numEx = $parts[0]; }
	}	
	for(my $k=0; $k<=$numEx; $k++)		# initialize
	{			
		for(my $t=0; $t<=$numEx; $t++)
		{
			$exonsA[$k][$t] = 0;
			$exonsB[$k][$t] = 0;
			$pairsA[$k][$t] = 0;
			$pairsB[$k][$t] = 0;				
		}
	}			
	while($i<=$#linesA)		# then fill in the matrices
	{
		my $ln = $linesA[$i++];
		if($ln =~ /INTRONS/)
		{	last; }	
		
		my @parts = split(' ', $ln);
		$exonsA[$parts[0]][$parts[1]] = $parts[2];				
	}				
	while($j<=$#linesB)
	{
		my $ln = $linesB[$j++];
		if($ln =~ /INTRONS/)
		{	last; }			
		
		my @parts = split(' ', $ln);
		$exonsB[$parts[0]][$parts[1]] = $parts[2];	
	}				
	for(my $k=1; $k<=$numEx; $k++)
	{
		print "$k $k ", $exonsA[$k][$k], " ", $exonsB[$k][$k], " ", $exonBegin[$k], " ", $exonEnd[$k], "\n";
		for(my $t=$k+1; $t<=$numEx; $t++)
		{
			if($exonsA[$k][$t] > 0 || $exonsB[$k][$t] > 0)	# if neither is positive, do not print
			{
				print "$k $t ", $exonsA[$k][$t], " ", $exonsB[$k][$t], " ", $exonEnd[$k], " ", $exonBegin[$t], "\n";
			}
		}		
	}				
	print $linesB[$j-1], "\n";
	
	# --- INTRONS ---	
	while($i<=$#linesA)
	{
		my $ln = $linesA[$i++];
		if($ln =~ /PAIRED/)
		{	last; }			
		my @parts = split(' ', $ln);			
		print "1 ", $parts[1], " ", $parts[2], " ", $parts[3], "\n";			
	}				
	while($j<=$#linesB)
	{
		my $ln = $linesB[$j++];
		if($ln =~ /PAIRED/)
		{
			print $ln, "\n";							
			last;
		}		
		my @parts = split(' ', $ln);			
		print "2 ", $parts[1], " ", $parts[2], " ", $parts[3], "\n";				
	}		
	
	# --- PAIRS ---	
	while($i<=$#linesA)		
	{
		my $ln = $linesA[$i++];
		if($ln =~ /ISOFORMS/)
		{	last; }				
		my @parts = split(' ', $ln);
		$pairsA[$parts[0]][$parts[1]] = $parts[2];		
	}				
	while($j<=$#linesB)
	{
		my $ln = $linesB[$j++];
		if($ln =~ /ISOFORMS/)
		{	last; }								
		my @parts = split(' ', $ln);			
		$pairsB[$parts[0]][$parts[1]] = $parts[2];	
	}						
	for(my $k=1; $k<=$numEx; $k++)
	{			
		for(my $t=$k+1; $t<=$numEx; $t++)
		{				
			if($pairsA[$k][$t] > 0 || $pairsB[$k][$t] > 0)
			{
				print "$k $t ", $pairsA[$k][$t], " ", $pairsB[$k][$t], "\n";
			}
		}		
	}							
	print $linesB[$j-1], "\n";
		
	# --- ISOFORMS ---
	while($i<=$#linesA)		
	{
		my $ln = $linesA[$i++];
		if($ln =~ /GENE/)
		{	last; }			
	}		
	while($j<=$#linesB)
	{
		my $ln = $linesB[$j++];
		if($ln =~ /GENE/)
		{	last; }		
		print $ln, "\n"; 
	}					
}

