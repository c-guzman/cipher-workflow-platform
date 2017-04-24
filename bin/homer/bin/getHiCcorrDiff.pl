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



my $genome = "";
my $config = HomerConfig::loadConfigFile();


my $minFilter = 0.10;
my $stdFilter = 4;
my $res = 50000;
my $superRes = 100000;
my $maxDist = "";
my $corrDepth = 3;
my $maxCPUs = 1;


sub printCMD() {
    print STDERR "\n\tUsage getHiCcorrDiff.pl <output prefix> <HiC directory1> <HiC directory2> [options]\n";
    print STDERR "\n\tOptions:\n";
    print STDERR "\t\t-res <#> (resolution in bp, default: $res)\n";
    print STDERR "\t\t-superRes <#> (super resolution in bp, i.e. window size, default: $superRes)\n";
    print STDERR "\t\t-corrDepth <#> (number of expected reads needed per data point when calculating correlation, default: 3)\n";
    print STDERR "\t\t-std <#> (exclude regions with sequencing depth exceeding # std deviations, default: $stdFilter)\n";
    print STDERR "\t\t-min <#> (exclude regions with sequencing depth less than this fraction of mean, default: $minFilter)\n";
    print STDERR "\t\t-maxDist <#> (maximum distance around regions to calculate similarity metrics, default: none)\n";
	print STDERR "\t\t-cpu <#> (number of CPUs to use, default: 1)\n";
    print STDERR "\n\tOutput files:\n";
    print STDERR "\t\t<prefix>.corrDiff.txt - peak file containing correlation values for each region\n";
    print STDERR "\t\t<prefix>.corrDiff.bedGraph - UCSC upload file showing correlation values across the genome\n";
    print STDERR "\n";
    exit;
}


if (@ARGV < 3) {
	printCMD();
}


my $prefix = $ARGV[0];
my $dir1 = $ARGV[1];
my $dir2 = $ARGV[2];

for (my $i=3;$i<@ARGV;$i++) {
    if ($ARGV[$i] eq '-res') {
        $res = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-superRes') {
        $superRes = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-corrDepth') {
        $corrDepth = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-std') {
        $stdFilter = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-min') {
        $minFilter = $ARGV[++$i];
    } elsif ($ARGV[$i] eq '-maxDist') {
        $maxDist = " -maxDist $ARGV[++$i] ";
    } elsif ($ARGV[$i] eq '-cpu') {
        $maxCPUs = $ARGV[++$i];
    } else {
        print STDERR "!!! Didn't recognize option \"$ARGV[$i]\"!!!!\n";
        printCMD();
    }
}

print STDERR "\tOutput files will start with: $prefix\n";
print STDERR "\tAnalyzing HiC directories: $dir1 vs $dir2\n\n";


if ($res > $superRes) {
	print STDERR "!!! Warning -superRes ($superRes) is smaller than -res ($res)...\n\t\tPausing for 5s...";
	`sleep 5`;
}
my $possibleRes = HomerConfig::getHiCBgRes($dir1,$superRes,$maxCPUs);
$possibleRes = HomerConfig::getHiCBgRes($dir2,$superRes,$maxCPUs);


my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";

print STDERR "\n\tWill analyze chrs:";
my @chrs = ();
`ls -1 "$dir1"/*.tags.tsv > "$tmpFile"`;
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
print STDERR "\n";
`rm "$tmpFile"`;
`sleep 2`;
#@chrs = ();
#for (my $i=1;$i<=$numChrs;$i++) {
#	push(@chrs, "chr" . $i);
#}
#push(@chrs, "chrX", "chrY", "chrM");
#@chrs = ('chr10','chr11');



my %chrFiles = ();

my $cpus = 0;
foreach(@chrs) {
	my $chr = $_;

	my $chrTmpFile = $rand . ".$chr.tmp";
	$chrFiles{$chr}= $chrTmpFile;

	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		#child process

		`analyzeHiC "$dir1" -ped "$dir2" -res $res -superRes $superRes -corr -corrDepth $corrDepth -chr $chr -min $minFilter -std $stdFilter $maxDist > "$chrTmpFile"`;
		
		exit(0);
	}
	if ($cpus >= $maxCPUs) {
		my $id = wait();
		$cpus--;
	}
}
my $id = 0;
while ($id >= 0) {
    $id = wait();
}
		
my %results = ();
foreach(keys %chrFiles) {
	my $chr = $_;
	my $file = $chrFiles{$chr};

	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $regionID = $line[0];
		my $v = $line[1];
		next if ($v < -1);
		$regionID =~ /^(.+?)\-(\d+)$/;
		my $start = $2;
		my $end = $start+$res;

		$results{$regionID} = {diff=>$v,c=>$chr,s=>$start,e=>$end};
	}
	close IN;
	`rm "$file"`;

}
my @all = sort {$results{$a}->{'c'} cmp $results{$b}->{'c'} || $results{$a}->{'s'} <=> $results{$b}->{'s'}} keys %results;
open OUT, ">$prefix.corrDiff.bedGraph";
open OUT2, ">$prefix.corrDiff.txt";
print OUT "track name=\"$prefix Correlation difference ($dir1 vs. $dir2)\" type=bedGraph  yLineMark=\"0.0\" alwaysZero=on maxHeightPixels=100:75:11 visibility=full viewLimits=-1:1 autoScale=off\n";
for (my $i=0;$i<@all-1;$i++) {
	my $id = $all[$i];
	my $c = $results{$id}->{'c'};
	my $c2 = $results{$all[$i+1]}->{'c'};
	next if ($c ne $c2);
	my $s = $results{$id}->{'s'};
	my $e = $results{$id}->{'e'};
	my $v = $results{$id}->{'diff'};
	print OUT "$c\t$s\t$e\t$v\n";
	print OUT2 "$id\t$c\t$s\t$e\t+\t$v\n";
}
close OUT;
close OUT2;
