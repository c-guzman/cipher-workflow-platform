#!/usr/bin/env perl
use warnings;
#


use POSIX;

$defaultFraction=5;
if (@ARGV < 1) {
	printCMD();
}

sub printCMD {
	print STDERR "\n\tUsage: scrambleFasta.pl <input.fasta> [options]\n";
	print STDERR "\tOutput FASTA file sent to stdout\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-n <#> (Number of scrambled output sequences, defulat: 5x input)\n";
	print STDERR "\n";
	exit;
}

my $inputFile = $ARGV[0];
my $numOutputSeqs = -1;
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-n') {
		$numOutputSeqs = $ARGV[++$i];
	} else {
		printCMD();
	}
}

my $rand = rand();
my $tmpfile = $rand . '.tmp';
`fasta2tab.pl "$inputFile" > "$tmpfile"`;

my $numInputSeqs = 0;
open IN, "$tmpfile";
while (<IN>) {
	$numInputSeqs++;
}
close IN;

if ($numInputSeqs < 1) {
	print STDERR "!!! $inputFile appears to be empty!?\n";
	exit;
}

my $fraction = $defaultFraction;
if ($numOutputSeqs>0) {
	$fraction = $numOutputSeqs/$numInputSeqs;
}

my $floor = floor($fraction);
my $leftover = $fraction-$floor;

my $id = 1;
open IN, "$tmpfile";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my @seq = split //,$line[1];
	
	for (my $i=0;$i<$floor;$i++) {
		my $s = scrambleSeq(\@seq);
		print ">$id\n$s\n";
		$id++;
	}
	my $r = rand();
	if ($r < $leftover) {
		my $s = scrambleSeq(\@seq);
		print ">$id\n$s\n";
		$id++;
	}
}
close IN;

`rm "$tmpfile"`;

sub scrambleSeq {
	my ($seq) = @_;
	my $N = @$seq;
	my $s = '';
	for (my $i=0;$i<$N;$i++) {
		$s .= $seq->[floor(rand()*$N)];
	}
	return $s;
}
