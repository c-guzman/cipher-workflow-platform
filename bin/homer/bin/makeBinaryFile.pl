#!/usr/bin/env perl
use warnings;

if (@ARGV < 2) {
	print STDERR "<base> <special> [1 - force special to be in base]\n";
	exit;
}

open IN, $ARGV[1];
my %conv = ();
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$conv{$line[0]} = 1;
}
close IN;
open IN, $ARGV[0];
my %done = ();
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (exists($done{$line[0]}));
	if (exists($conv{$line[0]})) {
		print "$line[0]\t1\n";
		$conv{$line[0]} = 0;
	} else {
		print "$line[0]\t0\n";
	}
	$done{$line[0]} = 1;
}
close IN;
if (@ARGV < 3) {
foreach(keys %conv) {
	next if ($conv{$_}==0);
	print "$_\t1\n";
}
}
