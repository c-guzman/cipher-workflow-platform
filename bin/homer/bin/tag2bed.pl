#!/usr/bin/env perl
use warnings;

if (@ARGV < 1) {
	print STDERR "<tag file i.e. chr1.tags.tsv> [Tag Length | 30]\n";
	exit;
}
my $len = 30;
if (@ARGV > 1) {
	$len = $ARGV[1];
}
my $id = 1;

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line =split /\t/;
	my $name = $line[0];
	my $chr= $line[1];
	my $start = $line[2]-1;
	my $end = $start+$len;
	my $dir= "+";
	my $v = $line[4];
	$v =~ s/\.0+$//;
	if ($line[3] == 1) {
		$dir = "-";
		$end = $line[2];
		$start = $line[2]-$len;
	}
	$name = "none". $id++ if ($name eq '');
	print "$chr\t$start\t$end\t$name\t$v\t$dir\n";
}
close IN;	
