#!/usr/bin/env perl
use warnings;


if (@ARGV < 2) {
	print STDERR "<data> <list of ids to remove> <# lines of header> <set to 1 if want to keep IDs>\n";
	exit;
}

my $header = 0;
if (@ARGV > 2) {
	$header = $ARGV[2];
}
my $oppFlag = 0;
if (@ARGV > 3) {
	$oppFlag = $ARGV[3];
}

my %list = ();
open IN, $ARGV[1];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $id = shift @line;
	$list{$id} = 1;
}
close IN;

my $c= 0;
open IN, $ARGV[0];
while (<IN>) {
	$c++;
	chomp;
	s/\r//g;
	if ($c <= $header) {
		print "$_\n";
		next;
	}
	my @line = split /\t/;
	my $id = shift @line;
	my $pflag = 0;
	if ($oppFlag) {
		if (exists($list{$id})) {
			$pflag = 1;
		}
	} else {
		if (!exists($list{$id})) {
			$pflag = 1;
		}
	}
	if ($pflag) {
		print "$id";
		foreach(@line) {
			print "\t$_";
		}
		print "\n";
	}

}
close IN;
	



