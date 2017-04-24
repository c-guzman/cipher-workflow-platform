#!/usr/bin/env perl
use warnings;



if (@ARGV < 4) {
	print STDERR "<file 1> <col for 1> <file 2> <col for 2> <1= header>\n";
	print STDERR "\tColumn counting starts with 1\n";
	exit;
}
my %data = ();
my $header = 0;
if (@ARGV > 4) {
	$header = $ARGV[4];
}
my $col1 = $ARGV[1];
my $col2 = $ARGV[3];
$col1--;
$col2--;
my @header = ();
my $count = 0;
my $nheader1 = 0;
my $nheader2 = 0;

open IN, $ARGV[2];
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($header && $count == 1) {
		@header = @line;
		$nheader2 = scalar(@header);
		next;
	}
	my $id = $line[$col2];
	$data{$id} = \@line;
}
close IN;

open IN, $ARGV[0];
$count = 0;
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($count == 1) {
		$nheader1 = scalar(@line);
	}

	print "$line[0]";
	my $printed = 1;
	for (my $i=1;$i<@line;$i++){ 
		print "\t$line[$i]";
		$printed++;
	}
	while ($printed < $nheader1) {
		print "\t";
		$printed++;
	}
	if ($count == 1 && $header) {
		foreach(@header) {
			print "\t$_";
		}
		print "\n";
		next;
	}
	if (exists($data{$line[$col1]})) {
		foreach(@{$data{$line[$col1]}}) {
			print "\t$_";
		}
		print "\n";
	} else {
		foreach(@header) {
			print "\t";
		}
		print "\n";
	}
	
}
close IN;
