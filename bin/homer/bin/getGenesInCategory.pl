#!/usr/bin/env perl
use warnings;

if (@ARGV < 2) {
	print STDERR "<ontology.genes> <GO ID>\n";
	exit;
}
my $found = 0;
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if ($line[0] ne $ARGV[1] && $line[1] ne $ARGV[1]);
	$found = 1;
	my @ids = split /\,/, $line[2];
	my $N = 0;
	foreach(@ids) {
		print "$_\n";
		$N++;
	}
	print STDERR "\n\tFound term $line[0]/$line[1] ($N gene ids)\n\n";
}
close IN;
if ($found == 0) {
	print STDERR "\n\t!!! Couldn't find term $ARGV[1] !!!\n\n";
}
