#!/usr/bin/perl -w
#
if (@ARGV < 2) {
	print STDERR "<basic annotation file> <chrom.sizes>\n";
	exit;
}
my %max = ();
open IN, $ARGV[1];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$max{$line[0]} = $line[1];
}
close IN;

open IN, $ARGV[0];
while (<IN>) {
	if (/^Intergenic/) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if (exists($max{$line[1]})) {
			if ($line[3] > $max{$line[1]}) {
				$line[3] = $max{$line[1]};
			}
		}
		print "$line[0]";
		for (my $i=1;$i<@line;$i++) {
			print "\t$line[$i]";
		}
		print "\n";
	}
}
close IN;
