#!/usr/bin/env perl
use warnings;
use POSIX;

if (@ARGV < 1) {
	print STDERR "<group file> - this program will randomize labels and output to stdout\n";
	exit;
}

my %group = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$group{$line[0]} = $line[1];
}
close IN;


my @genes = keys %group;
my $N = scalar(@genes);
my $P = 0;
my %rand = ();
foreach(@genes) {
	$rand{$_} = rand();
}

foreach(values %group) {
	if ($_ == 1) {
		$P++;
	}
}

@genes = sort {$rand{$a} <=> $rand{$b}} @genes;

for (my $i=0;$i<@genes;$i++) { 
	my $id = $genes[$i];
	if ($i < $P) {
		print "$id\t1\n";
	} else {
		print "$id\t0\n";
	}
}

