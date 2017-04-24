#!/usr/bin/env perl
use warnings;


if (@ARGV < 1) {
	print STDERR "<BED file>\n";
	exit;
}
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $chr = $line[0];
	my $pos = $line[1]+1;
	my $name = $line[3];
	my $v = $line[4];
	my $dir = 0;
	if ($line[5] eq '-') {
		$dir = 1;
		$pos = $line[2];
	}
	print "$name\t$chr\t$pos\t$dir\t$v\n";

}
close IN;
