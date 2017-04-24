#!/usr/bin/env perl
use warnings;

my $lineLength = 100;

if (@ARGV < 1) {
	print STDERR "<tab seq file>\n";
	exit;
}
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $len = length($line[1]);
	next if ($len < 1);
	print ">$line[0]\n";
	for (my $i=0;$i<$len;$i+= $lineLength) {
		my $s = substr($line[1],$i,$lineLength);
		print "$s\n";
	}

}
close IN;
