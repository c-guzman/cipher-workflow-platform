#!/usr/bin/env perl
use warnings;



if (@ARGV < 2) {
	print STDERR "<file> <data file to add> [place holder | '']\n";
	exit;
}

my $placeHolder = '';
if (@ARGV > 2) {
	$placeHolder = $ARGV[2];
}

my @head = ();
%data = ();
open IN, $ARGV[1];
$count = 0;
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($count == 1) {
		for (my $i=1;$i<@line;$i++) {
			push(@head, $line[$i]);
		}
		next;
	}
	my $name = shift @line;
	$data{$name} = \@line;
}
close IN;
#print STDERR "@head\n";

open IN, $ARGV[0];
$numBefore = 0;
$count = 0;
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if ($count==1) {
		print "$line[0]";
		$numBefore = @line;
		for (my $i=1;$i<@line;$i++) {
			print "\t$line[$i]";
		}
		foreach(@head) {
			print "\t$_";
		}
		print "\n";
		next;
	}
	my $name = $line[0];
	print $name;
	for (my $i=1;$i<@line;$i++) {
		print "\t$line[$i]";
	}
	if (@line < $numBefore) { 
		for (my $i=@line;$i<$numBefore;$i++) {
			print "\t$placeHolder";
		}
	}
	if (exists($data{$name})) {
		foreach(@{$data{$name}}) {
			print "\t$_";
		}
	}
	else {
		foreach(@head) {
			print "\t$placeHolder";
		}
	}
	print "\n";
}
close IN;


