#!/usr/bin/env perl
use warnings;


use POSIX;

if (@ARGV < 2) {
	print STDERR "\n\tusage: chopUpBackground.pl <target seq file> <background seq file> > choppedBack.seq.tsv\n\n";
	exit;
}

my $avg = 0;
my $N = 0;
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (@line < 2);
	next if ($line[1] eq '');
	
	$avg += length($line[1]);
	$N++;
}
close IN;

if ($N < 1) {
	print STDERR "!!! No target sequences !!!\n";
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
	next if (@line < 2);
	next if ($line[1] eq '');
	my $seq = $line[1];

	$oldN++;	
	my $L = length($line[1]);
	my $x = $L/$avg;

	for (my $i=0;$i<$x-0.5;$i++) {
		my $s = substr($seq,$i*$avg,$avg);
		my $v = $i+1;
		print "$line[0]-chopped$v\t$s\n";
		$newN++;
	}
}
close IN;

print STDERR "\tCreated $newN background sequences from $oldN (valid) original sequences\n";

