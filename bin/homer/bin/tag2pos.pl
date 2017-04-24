#!/usr/bin/env perl
use warnings;

my $start = -50;
my $end = 50;
my $oppDir = 0;
my $ID = 1;

if (@ARGV < 1) {
	print STDERR "\n\tUsage: tag2pos.pl <tag file> [start offset | \"given\" | 50] [end offset | \"given\" | 50] [-3p]\n";
	print STDERR "\n\tCreates a position/peak file out of a tag file\n";
	print STDERR "\tCurrent Start: $start\n";
	print STDERR "\tCurrent End: $end\n";
	print STDERR "\tUsing \"given\" will give on the dementions of the tag (i.e. uses the length field)\n";
	print STDERR "\tUsing -3p as the last argument will make offsets relative to the 3' end (uses length field)\n";
	print STDERR "\n";
	exit;
}
if (@ARGV > 1) {
	$start = $ARGV[1];
	if ($start eq 'given') {
		$start = 0;
	}
}
if (@ARGV > 2) {
	$end = $ARGV[2];
}

my $p3Flag = 0;
if (@ARGV > 3) {
	if ($ARGV[3] eq '-3p') {
		$p3Flag = 1;
	}
}
	
print STDERR "\tStart: $start relative to tag 5' end\n";
print STDERR "\tEnd:   $end relative to tag 5' end\n";

my $defLen = 1;

open IN, $ARGV[0];
my $id = 0;
while (<IN>) {
	chomp;
	my @line = split /\t/;
	my $name = $line[0];
	my $d = 0;
	my $len = $defLen;
	my $pos = $line[2];
	if (@line > 5) {
		$len = $line[5];
		$len = $defLen if ($len < $defLen);
	}
	if ($line[3] == 1) {
		$d = 1;
	}

	if ($p3Flag) {
		if ($d == 0) {
			$pos += $len;
		} else {
			$pos -= $len;
		}
	}

	my $s = $pos+$start;
	my $e = $pos;

	if ($end eq 'given') {
		$e = $pos+$len;
	} else {
		$e = $pos+$end;
	}

	if ($d == 1) {
		$e = $pos-$start;
		$s = $pos;
		if ($end eq 'given') {
			$s = $pos-$len;
		} else {
			$s = $pos-$end;
		}
	}
	if ($oppDir) {
		if ($d == 0) {
			$d=1;
		} else {
			$d=0;
		}
	}
	if ($name eq '') {
		$name = $ID++;
	}
	print "$name\t$line[1]\t$s\t$e\t$d\t$line[4]\n";

	$id++;
}
close IN;
