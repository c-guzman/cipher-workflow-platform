#!/usr/bin/env perl
use warnings;



if (@ARGV < 1) {
	print STDERR "\n\tcleanUpPeakFile.pl <peak file> [peak ID prefix]\n";
	print STDERR "\n\tFor making sure peak IDs are unique\n";
	print STDERR "\n";
	exit;
}
my $prefix = "";
if (@ARGV > 1) {
	$prefix = $ARGV[1];
}

my %ids = ();	
my %counts = ();	
open IN, $ARGV[0];
my $numDuplicates = 0;
my $total = 0;
while (<IN>) {
	chomp;	
	s/\r//g;
	my @line = split /\t/;
	my $id = $prefix . $line[0];
	$nid = $id;
	if (exists($ids{$id})) {
		my $count = 1;
		$nid = $id . "-Dup" . $counts{$id};
		$numDuplicates++;
	}
	$ids{$nid} = 1;
	$counts{$id}++;
	print "$nid";
	for (my $i=1;$i<@line;$i++) {
		print "\t$line[$i]";
	}
	print "\n";
	$total++;
}
if ($numDuplicates > 0) {
	print STDERR "\t$numDuplicates duplicate peak IDs out of $total total peaks\n";
}
