#!/usr/bin/env perl
use warnings;


if (@ARGV  < 1) {
	print STDERR "<group file> [1=get negatives]\n";
	exit;
}
my $flag = 0;
if (@ARGV > 1) {
	if ($ARGV[1] == 1) {
		$flag = 1;
	}
}
open IN, $ARGV[0];
while (<IN>) {
	$og = $_;
	chomp;
	s/\r//g;
	my @line= split /\t/;
	next if (@line < 2);
	if ($flag) {
		if ($line[1] eq '0') {
			print $og;
		}
	} else {
		if ($line[1] eq '1') {
			print $og;
		}
	}
}
close IN;


