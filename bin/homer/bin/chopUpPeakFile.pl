#!/usr/bin/env perl
use warnings;


use POSIX;

if (@ARGV < 2) {
	print STDERR "\n\tusage: chopUpPeakFile.pl <target peak file> <background peak file> > choppedBackground.peaks.tsv\n\n";
	exit;
}

my $avg = 0;
my $N = 0;
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (@line < 5);
	my $size = $line[3]-$line[2];
	$avg += $size;
	$N++;
}
close IN;

if ($N < 1) {
	print STDERR "!!! No target peaks !!!\n";
	exit;
}
$avg /= $N;
$avg = floor($avg);
print STDERR "\tAverage target sequence length = $avg bp\n";

my $oldN = 0;
my $newN = 0;
open IN, $ARGV[1];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (@line < 5);
	next if ($line[1] eq '');
	my $size = $line[3]-$line[2];

	$oldN++;	
	my $x = $size/$avg;

	for (my $i=0;$i<$x-0.5;$i++) {
		my $s = $line[2]+$i*$avg;
		my $e = $s+$avg;
		my $v = $i+1;
		print "$line[0]-chopped$v\t$line[1]\t$s\t$e\t$line[4]\n";
		$newN++;
	}
}
close IN;

print STDERR "\tCreated $newN background regions from $oldN (valid) original regions\n";

