#!/usr/bin/env perl
use warnings;



if (@ARGV < 1) {
	print STDERR "\n";
	print STDERR "\tcleanUpSequences.pl <tsv sequence file> [-min #] [-max #]\n";
	print STDERR "\tRemoves bad sequence characters\n";
	print STDERR "\t-max # and -min # removes sequences that are longer than max and shorter than min(optional)\n";
	print STDERR "\n";
	exit;
}
my $min = 0;
my $max = 1e10;
my $badSeq = 0;
my $tooLong= 0;
my $tooShort = 0;
my $totalSeq = 0;
my $seqfile = "";
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-max') {
		$max = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$min = $ARGV[++$i];
	} else {
		$seqfile = $ARGV[$i];
	}
}
open IN, $seqfile;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $len = length($line[1]);
	if ($len < $min) {
		$tooShort++;
		next;
	}
	if ($len > $max) {
		$tooLong++;
		next;
	}
	$line[1] =~ s/a/A/g;
	$line[1] =~ s/c/C/g;
	$line[1] =~ s/g/G/g;
	$line[1] =~ s/t/T/g;
	$line[1] =~ s/n/N/g;
	$totalSeq++;
	if ($line[1] =~ /[^ACGTN]/) {
		$badSeq++;
		$line[1] =~ s/[^ACGTN]/N/g;
	} else {
	}
	print "$line[0]\t$line[1]\n";
}
close IN;
if ($badSeq > 0) {
	print STDERR "\t!! $badSeq of $totalSeq contained bad nucleotide characters [not ACGTN], replaced with N\n";
}
if ($tooShort > 0) {
	print STDERR "\t\t$tooShort sequences < $min bp were removed\n"; 
}
if ($tooLong > 0) {
	print STDERR "\t\t$tooLong sequences > $max bp were removed\n"; 
}
