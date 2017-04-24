#!/usr/bin/env perl
use warnings;


if (@ARGV < 3) {
	print STDERR "<promoter file> <start> <end> <ref | -2000>\n";
	exit;
}
my $ref = -2000;
if (@ARGV > 3) {
	$ref = $ARGV[3];
}

open IN, $ARGV[0] or die "Could not open sequence file $ARGV[0]\n";
while (<IN>) {
	chomp;
	my @line =split /\t/;
	my $s = substr($line[1],$ARGV[1]-$ref, $ARGV[2]-$ARGV[1]);
	print "$line[0]\t$s\n";
}
close IN;
