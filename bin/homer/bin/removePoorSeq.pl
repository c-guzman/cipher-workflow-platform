#!/usr/bin/env perl
use warnings;

my $maxPercentN = 0.70;
my $minLength = 5;
my $maxLength = 100000;
my $errorRate = 0.25;

if (@ARGV < 1) {
	print STDERR "\n\tUsage: removePoorSeq.pl <seq file> [max % of Ns | default $maxPercentN] [minlength | $minLength] [maxlength | $maxLength]\n";
	print STDERR "\tRemoves sequences that are less than $minLength bp, more than $maxLength bp, or $maxPercentN N\n";
	print STDERR "\n";
	exit;
}

if (@ARGV > 1) {
	$maxPercentN = $ARGV[1];
}
if (@ARGV > 2) {
	$minLength = $ARGV[2];
}
if (@ARGV > 3) {
	$maxLength = $ARGV[3];
}

my $total = 0;
my $tooMuchN= 0;
my $tooShort= 0;
my $tooLong= 0;

open IN, $ARGV[0];
while (<IN>) {
	my $og = $_;
	chomp;
	my @line = split /\t/;
	my $L = length($line[1]);
	$total++;
	if ($L < $minLength) {
		$tooShort++;
		next;
	}
	if ($L > $maxLength) {
		$tooLong++;
		next;
	}
	next if ($L > $maxLength);
	$line[1] =~ s/[^ACGTN]/N/g;
	my $s = $line[1];
	$line[1] =~ s/N//g;
	my $l = length($line[1]);
	if (1-($l/$L) > $maxPercentN) {
		$tooMuchN++;
		next;
	}
	print "$line[0]\t$s\n";
}

my $numBad = $tooMuchN+$tooShort+$tooLong;
my $fraction = 0;
exit if ($total < 1);
if ($numBad/$total > $errorRate) {
	print STDERR "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	print STDERR "\tPotential Problem with Sequences (ignore if preparsing genome):\n";
}

if ($tooShort > 0) {
	print STDERR "\t$tooShort of $total removed because they were shorter than $minLength bp\n";
}
if ($tooLong > 0) {
	print STDERR "\t$tooLong of $total removed because they were longer than $maxLength bp\n";
}
if ($tooMuchN > 0) {
	my $pp = sprintf("%.2f",$maxPercentN*100) . "%";
	print STDERR "\t\t$tooMuchN of $total removed because they had >$pp Ns (i.e. masked repeats)\n";
	if ($tooMuchN/$total > $errorRate) {
		print STDERR "\tDue to 	high repeat content, consider using a non-repeat masked genome\n";
	}
}
if ($numBad/$total > $errorRate) {
	print STDERR "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
}

