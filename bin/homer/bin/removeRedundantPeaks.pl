#!/usr/bin/env perl
use warnings;


use POSIX;

sub printCMD() {

	print STDERR "\n\tUsage: removeRedundantPeaks.pl <peakfile> [options]\n";
	print STDERR "\n\tBy default, peaks with redundant positions [exact] will be removed\n";
	print STDERR "\tOutput is a new peak file sent to stdout. Additionally, peaks can be removed if they\n";
	print STDERR "\thave similar sequence composition (checked with BLAT).\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-blat <#> (Where % is the percentage match to be considered redundant, default=0.33)\n";
	print STDERR "\t\t-size <#> (size of peaks to be used for sequence similarity, default=given)\n";
	print STDERR "\t\t-genome <directory of FASTA files> (genome for extracting sequence)\n";
	print STDERR "\t\t-mask (use repeat masked sequence)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 1) {
	printCMD();
}
my $peakFile = $ARGV[0];
my $minBLATscore = 0.33;
my $size = '';
my $genome = "";
my $maskFlag = "";

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-blat') {
		$minBLATscore = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-mask') {
		$maskFlag = " -mask ";
	} elsif ($ARGV[$i] eq '-size') {
		$size = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
	} else {
		printCMD();
	}
}

my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";

`bed2pos.pl "$peakFile" -check > "$tmpFile"`;

my $initialPromoters = 0;
my %data = ();
my $positionFiltered = 0;
open IN, $tmpFile;
while (<IN>) {
	my $og = $_;
	chomp;
	s/\r//g;
	my @line  = split /\t/;

	next unless ($line[2] =~ /\d/);
	my $mid = floor($line[2]+$line[3])/2;
	my $dir = '+';
	if ($line[4] eq '-' || $line[4] eq '1') {
		$dir = '-';
	}
	my $loc = $line[1] . '_' . $mid . '_' . $dir;
	if (!exists($data{$loc})) {
		my %a = ();
		$data{$loc} = \%a;
	}
	my $v = 0;
	if (@line > 6) {
		$v =$line[6];
	}
	my $d = {og=>$og, v=>$v};
	$data{$loc}->{$line[0]}=$d;
	$initialPromoters++;
}
close IN;
`rm "$tmpFile"`;

my %newData = ();
my %pos = ();
foreach(keys %data) {	
	my $loc = $_;
	my @a = sort {$data{$loc}->{$b}->{'v'} <=> $data{$loc}->{$a}->{'v'} || $a cmp $b} keys %{$data{$loc}};
	my $best = $a[0];
	my $og = $data{$loc}->{$best}->{'og'};
	$newData{$loc} = $og;
	$pos{$best} = $og;
	$positionFiltered++;
}
print STDERR "\n\tKept $positionFiltered of $initialPromoters based on genome position\n";

if ($genome eq '') {
	foreach(values %newData) {
		print $_;
	}
	exit;
}

open OUT, ">$tmpFile";
foreach(values %newData) {
	print OUT $_;
}
close OUT;
if ($size ne '') {
	print STDERR "\tResizing...\n";
	`resizePosFile.pl "$tmpFile" $size > "$tmpFile2"`;
	`mv "$tmpFile2" "$tmpFile"`;
}
`makeBinaryFile.pl "$tmpFile" "$tmpFile" > "$tmpFile3"`;

`homerTools extract "$tmpFile" "$genome" $maskFlag > "$tmpFile2"`;
`findRedundantBLAT.pl "$tmpFile2" $minBLATscore > "$tmpFile"`;
`adjustRedunGroupFile.pl "$tmpFile3" "$tmpFile" > "$tmpFile2"`;
open IN, $tmpFile2;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	if (exists($pos{$line[0]})) {
		print "$pos{$line[0]}";
	}
}
close IN;

`rm "$tmpFile" "$tmpFile2" "$tmpFile3"`;
