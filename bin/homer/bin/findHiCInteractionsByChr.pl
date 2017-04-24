#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";

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
    print STDERR "\n\tUsage findHiCInteractionsByChr.pl <HiC directory> [options]\n";
	print STDERR "\n\tPurpose of this program is to automate the running of analyzeHiC for high-res interactions\n";
	print STDERR "\tby running it chromsome by chromosome and combining the results [will NOT find interchr]\n";
    print STDERR "\n\tOptions:\n";
    print STDERR "\t\t-res <#> (resolution in bp, default: $res)\n";
    print STDERR "\t\t-superRes <#> (super resolution in bp, i.e. window size, default: $superRes)\n";
	print STDERR "\t\t-minDist <#> (minimum interaction distance to search, default: superRes distance)\n";
	print STDERR "\t\t-maxDist <#> (minimum interaction distance to search, default: $maxDist)\n";
	print STDERR "\t\t-pvalue <#> (pvalue cutoff, default: $pvalue)\n";
	print STDERR "\t\t-zscore <#> (z-score cutoff, default: $zscore)\n";
	print STDERR "\t\t-cpu <#> (number of CPUs to use, default: $maxCPUs)\n";
	print STDERR "\t\t-ped <HiC directory 2> (Background HiC directory)\n";
    print STDERR "\t\t-std <#> (exclude regions with sequencing depth exceeding # std deviations, default: $stdFilter)\n";
    print STDERR "\t\t-min <#> (exclude regions with sequencing depth less than this fraction of mean, default: $minFilter)\n";
    print STDERR "\n\tOutput files:\n";
    print STDERR "\t\tinteractions sent to stdout\n";
    print STDERR "\n";
    exit;
}



if (@ARGV < 1) {
	printCMD();
}


my $superResSet = 0;
my $resSet = 0;
my $directory = $ARGV[0];
for (my $i=1;$i<@ARGV;$i++) {
    if ($ARGV[$i] eq '-res') {
        $res = $ARGV[++$i];
		$resSet = 1;
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

if ($minDist == -1234) {
	$minDist = $superRes;
}
if ($resSet && !$superResSet) {
	$superRes = $res;
}
if (!$resSet && $superResSet) {
	$res = $superRes;
}
print STDERR "\n\tres set to $res\n";
print STDERR "\tsuperRes set to $superRes\n\n";


my $possibleRes = HomerConfig::getHiCBgRes($directory,$superRes,$maxCPUs);
if ($ped ne '') {
	$possibleRes = HomerConfig::getHiCBgRes($ped,$superRes,$maxCPUs);
    $ped = "-ped \"$ped\"";
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

my $cpus=0;
foreach(@chrs) {
	my $chr = $_;
	my $outFile = "$rand.$chr.interactions";
	push(@resultFiles, $outFile);

	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		#child process
		print STDERR "`analyzeHiC $directory -res $res -superRes $superRes -chr $chr -std $stdFilter -min $minFilter -nomatrix -pvalue $pvalue -zscore $zscore -minDist $minDist -maxDist $maxDist -interactions $outFile $ped -center`;\n";
		`analyzeHiC "$directory" -res $res -superRes $superRes -chr $chr -std $stdFilter -min $minFilter -nomatrix -pvalue $pvalue -zscore $zscore -minDist $minDist -maxDist $maxDist -interactions "$outFile" $ped -center`;
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
my $z = 0;
foreach(@resultFiles) {
	my $file = $_;
	$z++;
	open IN , $file;
	my $c = 0;
	while (<IN>) {
		$c++;
		next if ($c == 1 && $z > 1);
		print $_;
	}
	close IN;
	`rm "$file"`;
}
