#!/usr/bin/env perl
#
#
#
#
if (@ARGV < 1) {
	print STDERR "<cons file>\n";
	exit;
}
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	my @line = split /\t/;
	my @values = split undef, $line[1];
	my $avg = 0;
	my $n = 0;
	foreach(@values) {
		$n++;
		$avg += $_;
	}	
	$avg /= $n if ($n > 0);
	print "$line[0]\t$avg\n";
}

