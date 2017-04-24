#!/usr/bin/env perl
use warnings;

if (@ARGV < 2) {
	print STDERR "<file> <convert file> <This Col> <to this> <#header|-1 for ontology> <keepIDFlag> <keepall flag>\n";
	print STDERR "column starts with 0\n";
	exit;
}

my $col1 = 0;
if (@ARGV > 2) {
	$col1 = $ARGV[2];
}
my $col2 = 1;
if (@ARGV > 3) {
	$col2 = $ARGV[3];
}
my $header = 0;
if (@ARGV > 4) {
	$header = $ARGV[4];
}
my $keepFlag = 0;
if (@ARGV > 5) {
	if ($ARGV[5] ne '0') {
		$keepFlag = 1;
	}
}
my $keepAllFlag = 0;
if (@ARGV > 6) {
	if ($ARGV[6] ne '0') {
		$keepAllFlag = 1;
	}
}

my %conv = ();
my $maxN = $col1;
if ($col2 > $col1) {
	$maxN = $col2;
}
open IN, $ARGV[1];
while (<IN>){ 
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (@line <= $maxN);
	next if ($line[$col2] eq '');
	next if ($line[$col1] eq '');
	$conv{$line[$col1]} = $line[$col2];
}
close IN;
 
my $count = 0;
open IN, $ARGV[0];
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	if ($count <= $header) {
		print "ID\t" if ($keepFlag == 1);
		print "$_\n";
		next;
	}
	next if ($_ eq '');
	my @line = split /\t/;

	if ($header < 0) { #i.e. ontology conversion
		my @ids = split /\,/, $line[1];
		my %con = ();
		foreach(@ids) {
			if (!exists($conv{$_})) {
				next;
			}
			$con{$conv{$_}} = 1;
		}
		my @a = keys %con;
		next if (@a < 1);
		print "$line[0]\t";
		my $c = 0;
		foreach(@a) {
			print "," if ($c > 0);
			$c++;
			print "$_";
		}
		print "\n";
	} else { #normal
		$line[0] =~ s/\s+$//;
		my $noDup = $line[0];
		$noDup =~ s/\-Dup\d+//;
		$noDup =~ s/\-HOMER\d+//;
		if (!exists($conv{$line[0]}) && !exists($conv{$noDup})) {
			#print STDERR "|$line[0]|\n";
			#next;
			if ($keepAllFlag) {
				print "$line[0]";
			} else {
				next;
			}
		} elsif (!exists($conv{$line[0]})) {
			print "$conv{$noDup}";
		} else {
			print "$conv{$line[0]}";
		}
		print "\t$line[0]" if ($keepFlag);
		for (my $i=1;$i<@line;$i++) {
			print "\t$line[$i]";
		}
		print "\n";
	}
}
close IN;
