#!/usr/bin/env perl
use warnings;


if (@ARGV < 2) {
	print STDERR "<file> <data file to add> [place holder value | '']\n";
	exit;
}
my $NULL = '';
if (@ARGV > 2) {
	$NULL = $ARGV[2];
}
%data = ();
open IN, $ARGV[1];
$totalColumn = 0;
$count = 0;
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($totalColumn < @line) {
		$totalColumn = @line;
	}
	my $name = shift @line;
	$data{$name} = \@line;
}
close IN;


open IN, $ARGV[0];
$count = 0;
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $name = $line[0];
	print $name;
	for (my $i=1;$i<@line;$i++) {
		print "\t$line[$i]";
	}
	if (exists($data{$name})) {
		my $j=0;
		foreach(@{$data{$name}}) {
			print "\t$_";
			$j++;
		}
		for (my $k=$j;$k<$totalColumn;$k++) {
			print "\t$NULL";
		}
	} else {
		for (my $k=0;$k<$totalColumn;$k++) {
			print "\t$NULL";
		}
	}
	print "\n";
}
close IN;


