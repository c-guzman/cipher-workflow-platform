#!/usr/bin/env perl
use warnings;
use lib "/genomics/homer/.//bin";
my $homeDir = "/genomics/homer/./";

# Copyright 2009-2014 Christopher Benner <cbenner@salk.edu>
# 
# This file is part of HOMER
#
# HOMER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HOMER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


use POSIX;
use HomerConfig;

my $minFilter = 0.20;
my $stdFilter = 4;
my $res = 2000;
my $superRes = 10000;
my $minDist = -1234;
my $maxDist = 10000000;
my $pvalue = 0.01;
my $zscore = 1.5;
my $maxCPUs = 1;
my $ped = "";

my $config = HomerConfig::loadConfigFile();

sub printCMD() {
    print STDERR "\n\tUsage makeHiCWashUfile.pl <output prefix> <HiC directory> [options]\n";
	print STDERR "\n\tCreates a bed.gz & tbi file to visualize with WashU Epigenome Browser\n";
	print STDERR "\t\tRequires bgzip and tabix to be installed!\n";
    print STDERR "\n\tOptions:\n";
    print STDERR "\t\t-res <#res> <min distance> <max distance> (resolution & max distance in bp\n";
	print STDERR "\t\t\tcan be specified multiple times for different resolutions.  The defaults are:\n";
	print STDERR "\t\t\t\t-res 10000   0        200000\n";
	print STDERR "\t\t\t\t-res 25000   200000   5000000\n";
	print STDERR "\t\t\t\t-res 100000  5000000  20000000\n";
	print STDERR "\t\t\t\t-res 500000  20000000 80000000\n";
	print STDERR "\t\t\t\t-res 1000000 80000000 300000000\n";
	print STDERR "\t\t-cpu <#> (number of CPUs to use, default: $maxCPUs)\n";
	#print STDERR "\t\t-ped <HiC directory 2> (Background HiC directory)\n";
    #print STDERR "\t\t-std <#> (exclude regions with sequencing depth exceeding # std deviations, default: $stdFilter)\n";
    #print STDERR "\t\t-min <#> (exclude regions with sequencing depth less than this fraction of mean, default: $minFilter)\n";
    print STDERR "\n\tOutput files:\n";
    print STDERR "\t\t<outputPrefix>.gz\n";
    print STDERR "\t\t<outputPrefix>.gz.tbi\n";
    print STDERR "\n";
    exit;
}



if (@ARGV < 2) {
	printCMD();
}

my %oglevels = ();
$oglevels{0} = {        s=>0,       e=>200000,    r=>10000};
$oglevels{200000} = {   s=>200000,  e=>600000,   r=>20000};
$oglevels{600000} = {   s=>600000,  e=>1800000,  r=>40000};
$oglevels{1800000} = {  s=>1800000, e=>5400000,  r=>80000};
$oglevels{5400000} = {  s=>5400000, e=>20000000, r=>160000};
%levels = ();

my $superResSet = 0;
my $resSet = 0;
my $outputPrefix = $ARGV[0];
my $directory = $ARGV[1];
for (my $i=2;$i<@ARGV;$i++) {
    if ($ARGV[$i] eq '-res') {
        my $res = $ARGV[++$i];
		my $start = $ARGV[++$i];
		my $end = $ARGV[++$i];
		$levels{$res} = {s=>$start,e=>$end,r=>$res};
    } elsif ($ARGV[$i] eq '-superRes') {
        $superRes = $ARGV[++$i];
		$superResSet = 1;
    } elsif ($ARGV[$i] eq '-ped') {
        $ped = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-std') {
        $stdFilter = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-minDist') {
        $minDist = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-maxDist') {
        $maxDist = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-cpu') {
        $maxCPUs = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-pvalue') {
        $pvalue = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-zscore') {
        $zscore = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-min') {
        $minFilter = $ARGV[++$i];
    } else {
        print STDERR "!!! Didn't recognize option \"$ARGV[$i]\"!!!!\n";
        printCMD();
    }
}

if (scalar(keys %levels) < 1) {
	%levels = %oglevels;
}

print STDERR "\n\tVerifying background models are present for each resolution level:\n";

@starts = sort {$a <=> $b} keys %levels;
foreach(@starts) {
	my $res = $levels{$_}->{'r'};
	my $minDist = $levels{$_}->{'s'};
	my $maxDist = $levels{$_}->{'e'};
	print STDERR "\t\tResolution: $res (range: $minDist - $maxDist)\n";
	my $possibleRes = HomerConfig::getHiCBgRes($directory,$res,$maxCPUs);
}

my $rand = rand();
my $tmpFile = $rand . ".tmp";

print STDERR "\n\tWill analyze chrs:";
my @chrs = ();
`ls -1 "$directory"/*.tags.tsv > "$tmpFile"`;
open IN, $tmpFile;
while (<IN>) {
    chomp;
    s/\r//g;
    s/\.tags\.tsv//;
	s/^.*\///;
    push(@chrs,$_);
    print STDERR " $_";
}
close IN;
`rm "$tmpFile"`;
print STDERR "\n";

#@chrs = ("chr19");

print STDERR "\tWill analyze chrs: @chrs\n";

my @resultFiles = ();
my $resultFileStr = "";

my $cpus=0;
foreach(@chrs) {
	my $chr = $_;
	my $outFile = "$rand.$chr.bed";
	push(@resultFiles, $outFile);
	$resultFileStr .= " \"$outFile\"";

	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		#child process
		my $fileStr = "";
		@starts = sort {$a <=> $b} keys %levels;
		for (my $i=0;$i<@starts;$i++) {
			$s = $starts[$i];
			my $res = $levels{$s}->{'r'};
			my $minDist = $levels{$s}->{'s'};
			my $maxDist = $levels{$s}->{'e'};
			print STDERR "\tResolution: $res (range: $minDist - $maxDist)\n";
			my $out = $outFile . ".$res";
			`analyzeHiC "$directory" -res $res -chr $chr -std $stdFilter -min $minFilter -minDist $minDist -maxDist $maxDist -washu -relative > \"$out\"`;
			$fileStr .= " \"$out\"";
		}
		`cat $fileStr > "$outFile.cat"`;
		`sort -k1,1 -k2,2g "$outFile.cat" > "$outFile"`;
		`rm $fileStr "$outFile.cat"`;
		exit(0);
	}
	if ($cpus >= $maxCPUs) {
		my $id = wait();
		$cpus--;
	}
}
my $id = 0;
while ($id >=0) {
	$id = wait();
}

`cat $resultFileStr > "$outputPrefix.bed"`;
`rm $resultFileStr`;
`bgzip -f "$outputPrefix.bed"`;
`tabix -p bed "$outputPrefix.bed.gz"`;

