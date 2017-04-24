#!/usr/bin/env perl
use warnings;


if (@ARGV < 1) {
	print STDERR "\n\tUsage: renamePeaks.pl <peak file> [prefix | RandID]\n";
	print STDERR "\tWill name peaks RandID1, RandID2, ...\n";
	print STDERR "\tOutputs to stdout\n";
	print STDERR "\n\ti.e. renamePeaks.pl ERa.peaks > ERa.new.peaks\n";
	print STDERR "\n";
	exit;
}
my $prefix = "RandID";
if (@ARGV > 1) {
	$prefix = $ARGV[1];
}

my $id = 0;
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $name = $prefix . $id++;
	print "$name";
	for (my $i=1;$i<@line;$i++) {
		print "\t$line[$i]";
	}
	print "\n";
}
