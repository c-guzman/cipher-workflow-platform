#!/usr/bin/env perl
use warnings;
use POSIX;

if (@ARGV < 1) {
	print STDERR "\n\tresizePosFile.pl <peak/position file> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t2nd argument: # - new size of peaks\n";
	print STDERR "\t\t3rd argument: # - adjust peak position\n";
	print STDERR "\t\t4th argument: yes - for a header line\n";
	#print STDERR "\t\t-size <#> (resize peak [around center] to this size, supports -size <#,#>)\n";
	#print STDERR "\t\t-rsize <#> (resize peak by this relative size, supports -size <#,#>)\n";
	#print STDERR "\t\t-move <#> (move peak [relative to strand], default: 0)\n";
	#print STDERR "\t\t-5p (recenter peak on 5' end of region)\n";
	#print STDERR "\t\t-3p (recenter peak on 3' end of region)\n";
	print STDERR "\n";
	exit;
}
my $size = 0;
if (@ARGV > 1) {
	$size = $ARGV[1];
}
my $move = 0;
if (@ARGV > 2) {
	$move = $ARGV[2];
}
my $headerFlag = 0;
if (@ARGV > 3) {
	$headerFlag = 1;
}
my $relativeFlag = 0;
if ($size =~ s/^\+//) {
	$relativeFlag = 1;
} elsif ($size =~ s/^\-//) {
	$relativeFlag = 1;
	$size = -1*$size;
}

my $half = $size/2;

my $removed = 0;	

open IN, $ARGV[0];
my $count = 0;
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	next if (/^#/);
	if ($count == 1 && $headerFlag) {
		print "$_\n";
		next;
	}
	my @line = split /\t/;

	next if (@line < 5);
	if (!($line[2] =~ /^\d+$/)) {
		next;
	}
	foreach(@line) {
		s/^\s*//;
		s/\s*$//;
	}

	my $center = (($line[2]+$line[3])/2);
	my $ogHalf = $line[3]-$center;
	if ($line[4] eq '-' || $line[4] eq '1') {
		$center -= $move;
	} else {
		$center += $move;
	}
	my $s = floor($center - $half);
	my $e = floor($center + $half);
	if ($relativeFlag) {
		$s = $center - $ogHalf - $half;
		$e = $center + $ogHalf + $half;
		if ($e < $s) {
			$removed++;
			next;
		}
	}

	print "$line[0]\t$line[1]\t$s\t$e\t$line[4]";
	for (my $i=5;$i<@line;$i++) {
		print "\t$line[$i]";
	}
	print "\n";
}
close IN;


if ($removed > 0) {
	print STDERR " $removed peaks removed because there size was less than 0\n";
}
