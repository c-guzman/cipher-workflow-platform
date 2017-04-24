#!/usr/bin/env perl
use warnings;



if (@ARGV < 1) {
	print STDERR "\n\tusage: checkPeakFile.pl <peak/pos file>\n";
	print STDERR "\n\tChecks if there are any formatting/data entry issues with peak file\n";
	
	print STDERR "\n";
	exit;
}
my $totalPeakLines = 0;
my $totalCommentLines = 0;
my $totalRedun= 0;
my $totalShort = 0;
my $badNumbers = 0;
my $badDirection = 0;
open IN, $ARGV[0];
my %data = ();
while (<IN>) {
	chomp;
	s/\r//g;
	if (/^#/) {
		$totalCommentLines++;
		next;
	}
	my @line = split /\t/;
	$totalPeakLines++;
	if (exists($data{$line[0]})) {
		$totalRedun++;
	}
	$data{$line[0]}++;
	if (@line < 5) {
		$totalShort++;
		next;
	}
	if ($line[2] !~ /^[\d\-\.]+$/ || $line[3] !~ /^[\d\-\.]+$/) {
		$badNumbers++;
	}
	if ($line[4] !~ /[\+\-01]/) {
		$badDirection++;
	}
}
close IN;

print STDERR "\n";
print STDERR "\tPeak File Statistics:\n";
print STDERR "\t\tTotal Peaks: $totalPeakLines\n";
print STDERR "\t\tRedundant Peak IDs: $totalRedun\n";
print STDERR "\t\tPeaks lacking information: $totalShort (need at least 5 columns per peak)\n";
print STDERR "\t\tPeaks with misformatted coordinates: $badNumbers (should be integer)\n";
print STDERR "\t\tPeaks with misformatted strand: $badDirection (should be either +/- or 0/1)\n";
print STDERR "\n";
my $good = 1;
if ($totalPeakLines < 2) {
	print STDERR "\tOnly $totalPeakLines peaks found! If many more expected, maybe the file format is for Macs only\n";
	print STDERR "\tTry running the command: changeNewLine.pl \"$ARGV[0]\"\n";
	print STDERR "\n";
	exit;
}
if ($totalRedun > 0) {
	print STDERR "\tRedunant Peaks found: Remove or rename these or some programs may have trouble...\n";
	print STDERR "\n";
	$good = 0;
}
if ($good == 1) {
	print STDERR "\tPeak file looks good!\n";
	print STDERR "\n";
}
