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


my $stdFilter = 4;
my $minFilter = 0.15;
my $res = 50000;
my $superRes = 100000;
my $corrDepth = 3;
my $activeZones = "";
my $inActiveZones = "";
my $maxCPUs = 1;

my $genome = "";

my $config = HomerConfig::loadConfigFile();


sub printCMD() {
	print STDERR "\n\tUsage runHiCpca.pl <output prefix> <HiC directory> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-res <#> (resolution in bp, default: $res)\n";
	print STDERR "\t\t-superRes <#> (super resolution in bp, i.e. window size, default: $superRes)\n";
	print STDERR "\t\t-active <K4me+ regions> (Regions to use to help decide sign on principal component [active=+])\n";
	print STDERR "\t\t-inactive <K4me- regions> (Regions to use to help decide sign on principal component [inactive=-])\n";
	print STDERR "\t\t\t-genome <genome> (If you don't have seed regions, this will use the TSS file as active seeds)\n";
	print STDERR "\t\t-corrDepth <#> (number of expected reads needed per data point when calculating correlation, default: 3)\n";
	print STDERR "\t\t-std <#> (exclude regions with sequencing depth exceeding # std deviations, default: $stdFilter)\n";
	print STDERR "\t\t-min <#> (exclude regions with sequencing depth less than this fraction of mean, default: $minFilter)\n";
	print STDERR "\t\t-rpath <path to R executable> (If R is not accessible via the \$PATH variable)\n";
	print STDERR "\t\t-pc <#> (number of principal components to report, default: 1)\n";
	print STDERR "\t\t-cpu <#> (number of CPUs to use, default: 1)\n";
	print STDERR "\n\tOutput files:\n";
	print STDERR "\t\t<prefix>.PC1.txt - peak file containing coordinates along the first 2 principal components\n";
	print STDERR "\t\t<prefix>.PC1.bedGraph - UCSC upload file showing PC1 values across the genome\n";
	print STDERR "\n";
	exit;
}


if (@ARGV < 2) {
	printCMD();
}

my $prefix = $ARGV[0];
my $directory = $ARGV[1];
my $rpath = "R";
my $numPC=1;
for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-res') {
		$res = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-superRes') {
		$superRes = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rpath') {
		$rpath = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-active') {
		$activeZones = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-inactive') {
		$inActiveZones = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-corrDepth') {
		$corrDepth = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pc') {
		$numPC = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-std') {
		$stdFilter = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$minFilter = $ARGV[++$i];
	} else {
		print STDERR "!!! Didn't recognize option \"$ARGV[$i]\"!!!!\n";
		printCMD();
	}
}

if ($res > $superRes) {
	print STDERR "!!! Warning -superRes ($superRes) is smaller than -res ($res)...\n\t\tPausing for 5s...";
	`sleep 5`;
}

my $check = `which "$rpath"`;
if ($check eq "") {
    print STDERR "\n\t!!! Could not detect the program R in the executable path\n";
    print STDERR "\t\tIf R is install, but not on the executable path, specify it's location with -rpath\n";
    print STDERR "\t\tOtherwise, download and install it from http://www.r-project.org/\n\n";
    exit;
}


$genome =~ s/r$//;
my $genomeDir = "";
if ($genome eq 'none') {
    print STDERR "\tNo genome selected (\"none\") - ignored\n";
} elsif (!exists($config->{'GENOMES'}->{$genome})) {
	if ($genome ne '') {
		print STDERR "\t\"$genome\" not recognized - ignored\n";
	}
} else {
    $genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'};
}

if ($activeZones eq '' && $inActiveZones eq '') {
	if ($genomeDir ne '') {
		$activeZones = $genomeDir . "/" . $genome . ".tss";
		print STDERR "\tUsing TSS from $genome to assign sign to active(+) vs. inactive(-) PC1 values\n"
	} else {
		print STDERR "\t!Warning: No seed regions used - sign on PC1 value may not be consistent across chromosomes...\n";
		print STDERR "\t(pausing for 5 seconds...)\n";
		`sleep 5`;
	}
}


print STDERR "\tOutput files will start with: $prefix\n";
print STDERR "\tAnalyzing HiC directory: $directory\n\n";

my $force = 1;
my $possibleRes = HomerConfig::getHiCBgRes($directory,$superRes,$maxCPUs);


my $rand = rand();
my $tmpFile = $rand . ".tmp";

my @chrs = ();
`ls -1 "$directory"/*.tags.tsv > "$tmpFile"`;
open IN, $tmpFile;
while (<IN>) {
	chomp;
	s/\r//g;
	s/\.tags\.tsv//;
	s/^.*\///;
	push(@chrs,$_);
}
close IN;
`rm "$tmpFile"`;
#for (my $i=1;$i<=$numChrs;$i++) {
#	push(@chrs, "chr" . $i);
#}
#push(@chrs, "chrX", "chrY", "chrM");
#@chrs = ("chr19");

print STDERR "\tWill analyze chrs: @chrs\n";

my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";
my $tmpFile4 = $rand . ".4.tmp";
my $tmpFile5 = $rand . ".5.tmp";
my $tmpFile6 = $rand . ".6.tmp";

my $rInputFile = $tmpFile5;
my $rOutputFile = $tmpFile6;

open EIGEN, ">$prefix.PC1.txt";
print EIGEN "#peakID\tchr\tstart\tend\tstrand\tPC1";
open ALL, ">$prefix.all.txt";
print ALL "Cluster\tPC1";
for (my $i=1;$i<$numPC;$i++){ 
	my $j = $i+1;
	print EIGEN "\tPC$j";
	print ALL "\tPC$j";
}
print EIGEN "\n";
print ALL "\n";

my %results = ();

my %rOutputFiles = ();

my $cpus = 0;
foreach(@chrs) {
	my $chr = $_;

	my $chrTmpFile = $rand . "$chr.tmp";
	my $chrTmpFile2 = $rand . "$chr.2.tmp";
	my $rInputFile = $rand . "$chr.PCA.input.txt";
	my $rOutputFile = $rand . "$chr.PCA.txt";
	$rOutputFiles{$chr} = $rOutputFile;

	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		#child process

		`analyzeHiC "$directory" -res $res -superRes $superRes -corr -corrDepth $corrDepth -chr $chr -std $stdFilter -min $minFilter > "$chrTmpFile"`;

		`prepForR.pl "$chrTmpFile" > "$chrTmpFile2"`;

		open RSCRIPT, ">$rInputFile";
		print RSCRIPT "data <- as.matrix(read.table(\"$chrTmpFile2\"))\n";
		print RSCRIPT "pcdat = princomp(data,cor=TRUE)\n";
		my $end = 1+$numPC;
		print RSCRIPT "write.table(as.matrix(pcdat\$scores[,1:$end]),file=\"$rOutputFile\",sep=\"\t\")\n";
		close RSCRIPT;

		`"$rpath" --no-save < "$rInputFile"`;
		`rm -f "$rInputFile" "$chrTmpFile" "$chrTmpFile2"`;
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

foreach(@chrs) {
	my $chr = $_;
	my $rOutputFile = $rOutputFiles{$chr};
	my $coeff = getDirection($rOutputFile, $activeZones, $inActiveZones);


	open IN, $rOutputFile;
	my $count = 0;
	while (<IN>) {
		$count++;
		next if ($count < 2);
		s/\"//g;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $regionID = $line[0];
		my $v1 = $line[1]*$coeff;
		my @v = ();
		push(@v, $v1);
		for (my $i=1;$i<$numPC;$i++){ 
			push(@v, $line[1+$i]);
		}
		$regionID =~ /^(.+?)\-(\d+)$/;
		my $start = $2;
		my $end = $start+$res;

		print EIGEN "$line[0]\t$chr\t$start\t$end\t+\t$v1";
		print ALL "$line[0]\t$v1";
		for (my $i=1;$i<$numPC;$i++){ 
			print EIGEN "\t$v[$i]";
			print ALL "\t$v[$i]";
		}
		print EIGEN "\n";
		print ALL "\n";
	}
	close IN;
	`rm "$rOutputFile"`;

}
close ALL;
close EIGEN;
`cluster2bedgraph.pl "$prefix.all.txt" $res $prefix > "$prefix.PC1.bedGraph"`;
`rm "$prefix.all.txt"`;
#print STDERR "\tOuput Files:\n\t\t$prefix.active.txt\n\t\t$prefix.inactive.txt\n\t\t$prefix.ActiveVsInactive.bedGraph\n\n";

sub getDirection {
	my ($rfile,$activeZones, $inActiveZones) = @_;
	return 1 if ($activeZones eq '' && $inActiveZones eq '');
	open POS, ">$tmpFile3";
	open NEG, ">$tmpFile4";
	open IN, $rfile;
	my $count = 0;
	my $posTotal=0;
	my $negTotal=0;
	while (<IN>) {
		$count++;
		next if ($count < 2);
		chomp;
		s/\"//g;
		s/\r//g;
		my @line = split /\t/;
		my $regionID = $line[0];
		my $v1 = $line[1];
		my $v2 = $line[2];
		$regionID =~ /^(.+?)\-(\d+)$/;
		my $chr = $1;
		my $start = $2;
		my $end = $start+$res;
		if ($v1 > 0) {
			print POS "$line[0]\t$chr\t$start\t$end\t+\t$v1\n";
			$posTotal++;
		} else {
			print NEG "$line[0]\t$chr\t$start\t$end\t+\t$v1\n";
			$negTotal++;
		}
	}
	close IN;
	close POS;
	close NEG;
	my $posActive = 1;
	my $negActive = 1;
	if ($activeZones ne '') {
		`mergePeaks -d given "$tmpFile3" "$activeZones" -cobound 1`;
		$r = `wc -l coBoundBy1.txt`;
		$r =~ s/^\s+//;
		$r =~ s/\s.*//;
		$posActive = $r;
		`mergePeaks -d given "$tmpFile4" "$activeZones" -cobound 1`;
		$r = `wc -l coBoundBy1.txt`;
		$r =~ s/^\s+//;
		$r =~ s/\s.*//;
		$negActive = $r;
	}
	my $posInActive = 1;
	my $negInActive = 1;
	if ($inActiveZones ne '') {
		`mergePeaks -d given "$tmpFile3" "$inActiveZones" -cobound 1`;
		$r = `wc -l coBoundBy1.txt`;
		$r =~ s/^\s+//;
		$r =~ s/\s.*//;
		$posInActive = $r;
		`mergePeaks -d given "$tmpFile4" "$inActiveZones" -cobound 1`;
		$r = `wc -l coBoundBy1.txt`;
		$r =~ s/^\s+//;
		$r =~ s/\s.*//;
		$negInActive = $r;
	}

	`rm "$tmpFile3" "$tmpFile4" "coBoundBy0.txt" "coBoundBy1.txt"`;

	my $diffPos = 1;
	my $diffNeg = 1;
	if ($activeZones ne '' && $inActiveZones ne '') {
		$diffPos = (1+$posActive)/(1+$posInActive);
		$diffNeg = (1+$negActive)/(1+$negInActive);
	} elsif ($activeZones ne '') {
		$diffPos = (1+$posActive)/($posTotal+1);
		$diffNeg = (1+$negActive)/($negTotal+1);
	} elsif ($inActiveZones ne '') {
		$diffNeg = (1+$posInActive)/($posTotal+1);
		$diffPos = (1+$negInActive)/($negTotal+1);
	}
	print STDERR "\tRatio of pos in active to inactive: $diffPos\n";
	print STDERR "\tRatio of neg in active to inactive: $diffNeg\n";
	if ($diffPos > $diffNeg) {
		return 1;
	} else {
		return -1;
	}
}
