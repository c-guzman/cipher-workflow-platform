#!/usr/bin/perl -w
if (@ARGV < 3) {
	print STDERR "\t<prefix name i.e. human> <name-mRNA.seq> <refGene.txt>\n";
	print STDERR "\t\twill output name-mRNA-3UTR.seq name-mRNA-5UTR.seq\n";
	exit;
}
my $fname = $ARGV[0];
my $seqfile = $ARGV[1];
my $reffile = $ARGV[2];
my %seq = ();
open IN, $seqfile;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$seq{$line[0]} = $line[1];
}
close IN;

open UTR5, ">$fname-mRNA-5UTR.seq";
open UTR3, ">$fname-mRNA-3UTR.seq";

open IN ,$reffile;
while(<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $id = $line[1];
	if (!exists($seq{$id})) {
		next;
	}
	my $dir = 0;
	$dir = 1 if ($line[3] eq '-');

	my $cdsStart = $line[6]+1;
	my $cdsEnd = $line[7]+1;

	my $startPos = 0;
	my $endPos = 0;
	
	my @starts = split /\,/,$line[9];
	foreach(@starts) {
		$_++;
	}
	my @ends = split /\,/,$line[10];
	my $curLen = 0;
	my $flag = 0;
	for (my $i=0;$i<@starts;$i++) {
		#print STDERR "exon $i\t$curLen\t$flag\n";
		if ($flag == 0) {
			if ($starts[$i] <= $cdsStart
			 	&& $ends[$i] >= $cdsStart) {
				$startPos = $curLen+($cdsStart-$starts[$i]);
				$flag = 1;
			}
		}
		if ($flag == 1) {
			if ($starts[$i] <= $cdsEnd 
				&& $ends[$i] >= $cdsEnd) {
				$endPos = $curLen+($cdsEnd-$starts[$i]);
				$flag = 2;
			}
		}
		$curLen += $ends[$i]-$starts[$i]+1;
	}

	#print STDERR "$id\t$startPos\t$endPos\t$curLen\n";
	if ($dir == 0) {
		if ($startPos == $endPos && $startPos == $curLen) {
			# noncoding RNA - assign to both 3' and 5' UTRs for now
			$endPos = 0;
		}
		if ($startPos > 0) {
			my $utr5 = substr($seq{$id},0,$startPos);
			print UTR5 "$id\t$utr5\n";
		}
		if ($endPos < $curLen) {
			my $utr3 = substr($seq{$id},$endPos);
			print UTR3 "$id\t$utr3\n";
		}
	} else {
		if ($startPos == $endPos && $startPos == $curLen) {
			# noncoding RNA - assign to both 3' and 5' UTRs for now
			$endPos = 0;
		}
		if ($endPos < $curLen) {
			my $utr5 = substr($seq{$id},0,$curLen-$endPos);
			print UTR5 "$id\t$utr5\n";
		}
		if ($startPos > 0) {
			my $utr3 = substr($seq{$id},$curLen-$startPos);
			print UTR3 "$id\t$utr3\n";
		}
	}
		

}
close IN;
