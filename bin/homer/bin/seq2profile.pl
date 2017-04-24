#!/usr/bin/env perl
use warnings;

if (@ARGV < 1) {
	print STDERR "\n\tUsage: seq2profile.pl <consensus sequence> [# mismatches | default=0] [name | default=consensus]\n";
	print STDERR "\n\tExample: seq2profile.pl ACNNGGTGT 1 motif1 > output.motif\n";
	print STDERR "\n\tConsensus sequence must be composed of IUPAC characters:  ACGT RYMKWS BDHV N\n";
	print STDERR "\tNote: # of mismatches is a tricky conversion for degenerate characters [RYMKWS/BDHV]\n";
	print STDERR "\t\tand should be checked in the results - mimatches are based on ACGT values\n";
	print STDERR "\n";
	exit;
}

my $minP = 0.001;
my $eps = 0.01;

my $seq = $ARGV[0];
my $mismatches = 0;
if (@ARGV > 1) {
	$mismatches = $ARGV[1];
}
my $name = $seq;
if (@ARGV > 2) {
	$name = $ARGV[2];
}

my $len = length($seq);
my $nSeq = $seq;
$nSeq =~ s/N//g;
my $nlen = length($nSeq);
my $numNs = $len - $nlen;

$nSeq = $seq;
$nSeq =~ s/[RYMKWS]//g;
my $n2len = length($nSeq);
my $num2s = $len - $n2len;

$nSeq = $seq;
$nSeq =~ s/[BDHV]//g;
my $n3len = length($nSeq);
my $num3s = $len - $n3len;

my $num1s = $len - $numNs - $num2s - $num3s;

my $maxP = 1-$minP*3;
my $maxP2 = (1-2*$minP)/2;
my $maxP3 = (1-$minP)/3;

my $unitScore = log($maxP/0.25);
my $unitBadScore = log($minP/0.25);
my $unitNScore = log(0.25/0.25);
my $unit2Score = log($maxP2/0.25);
my $unit3Score = log($maxP3/0.25);


my $expectedScore = $unitScore*($num1s - $mismatches) + $unitBadScore*$mismatches + $numNs*$unitNScore 
							+ $num2s*$unit2Score + $num3s*$unit3Score;

$expectedScore -= $eps;


print ">$seq\t$name\t$expectedScore\n";
for (my $i=0;$i<$len;$i++) {
	my $c = substr($seq, $i,1);
	if ($c eq 'A') {
		print "$maxP\t$minP\t$minP\t$minP\n";
	} elsif ($c eq 'C') {
		print "$minP\t$maxP\t$minP\t$minP\n";
	} elsif ($c eq 'G') {
		print "$minP\t$minP\t$maxP\t$minP\n";
	} elsif ($c eq 'T') {
		print "$minP\t$minP\t$minP\t$maxP\n";
	} elsif ($c eq 'N') {
		print "0.25\t0.25\t0.25\t0.25\n";
	} elsif ($c eq 'R') {
		print "$maxP2\t$minP\t$maxP2\t$minP\n";
	} elsif ($c eq 'Y') {
		print "$minP\t$maxP2\t$minP\t$maxP2\n";
	} elsif ($c eq 'M') {
		print "$maxP2\t$maxP2\t$minP\t$minP\n";
	} elsif ($c eq 'K') {
		print "$minP\t$minP\t$maxP2\t$maxP2\n";
	} elsif ($c eq 'W') {
		print "$maxP2\t$minP\t$minP\t$maxP2\n";
	} elsif ($c eq 'S') {
		print "$minP\t$maxP2\t$maxP2\t$minP\n";
	} elsif ($c eq 'B') {
		print "$minP\t$maxP3\t$maxP3\t$maxP3\n";
	} elsif ($c eq 'D') {
		print "$maxP3\t$minP\t$maxP3\t$maxP3\n";
	} elsif ($c eq 'H') {
		print "$maxP3\t$maxP3\t$minP\t$maxP3\n";
	} elsif ($c eq 'V') {
		print "$maxP3\t$maxP3\t$maxP3\t$minP\n";
	}
}

