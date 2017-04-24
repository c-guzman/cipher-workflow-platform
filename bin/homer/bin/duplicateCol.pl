#!/usr/bin/env perl
use warnings;



if (@ARGV < 1) {
	print STDERR "<file> [col:starting at 0]\n";
	exit;
}
my $col = 0;
if (@ARGV > 1) {
	$col= $ARGV[1];
}
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $c=0;
	for (my $i=0;$i<@line;$i++) {
		if ($c > 0) {
			print "\t";
		}
		$c++;
		print "$line[$i]";
		if ($i==$col) {
			print "\t$line[$i]";
		}
	}
	print "\n";
		
}
close IN;
