#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";


# Copyright 2009, 2010 Christopher Benner <cbenner@ucsd.edu>
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

my $knownMotifs = $homeDir . "/data/knownTFs/JASPAR.matrix";

my $config = HomerConfig::loadConfigFile();
$outputFileName = "tagFreq.txt";
$outputGCFileName = "tagGCcontent.txt";
$outputCpGFileName = "tagCpGcontent.txt";

print STDERR "!!! It is highly recommend that you use 'makeTagDirectory' with the '-update' option instead !!!!\n";
if (@ARGV < 2) {
	cmdLineOptions();
}
print STDERR "\tFor example: makeTagDirectory \"$ARGV[0]\" -update  -checkGC -genome \"$ARGV[1]\"\n";
print STDERR "\twill continue in 5 seconds...\n";
`sleep 5`;
sub cmdLineOptions {
	print STDERR "\n\tProgram will check for sequence bias in your mapped sequencing tags\n";
	print STDERR "\n\tUsage: checkTagBias.pl <Tag Directory> <genome> [additional options]\n";
	print STDERR "\tOutput files (NOTE: located in the <Tag Directory>)\n";
	print STDERR "\t\t$outputFileName\n";
	print STDERR "\t\t$outputGCFileName\n";
	print STDERR "\t\t$outputCpGFileName\n";
	print STDERR "\t\tchr(N).tags.tsv.seq (if -keep is used)\n";

	print STDERR "\n\tPossible Genomes:\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}\n";
	}

	print STDERR "\n\tBasic options:\n";
	print STDERR "\t\t-start <#> (offset to start frequency calculation, default=-50)\n";
	print STDERR "\t\t-end <#> (offset to end frequency calculation, default=200)\n";
	print STDERR "\t\t-gcstart <#> (offset to start GC% calculation, default=0)\n";
	print STDERR "\t\t-gcend <#> (offset to end GC% calculation, default=100)\n";
	print STDERR "\t\t-keep (do not delete sequence files for each tag position)\n";
	print STDERR "\t\t-3p (if tags have lengths, align them at the 3' end)\n";
	print STDERR "\t\t-pos or -neg (only check positive or negative strands)\n";
	print STDERR "\t\t-skipGC (skip GC% calculation)\n";
	print STDERR "\t\t-skipFreq (skip nucleotide frequency calculation)\n";
	print STDERR "\t\t-prefix <filename> (output files will start with prefix name, default: directory name)\n";
	print STDERR "\t\t-mask (use repeat-masked genome)\n";
	print STDERR "\n";
	exit;
}


my $directory = $ARGV[0];
my $genome = $ARGV[1];
my $maskFlag = "";
if ($genome =~ s/r$//) {
	$maskFlag = " -mask ";
}
if (!exists($config->{'GENOMES'}->{$genome})) {
	print STDERR "!!!!Genome $genome not found in $homeDir/config.txt\n\n";
    print STDERR "\tTo check if is available, run \"perl $homeDir/configureHomer.pl -list\"\n";
    print STDERR "\tIf so, add it by typing \"perl $homeDir/configureHomer.pl -install $genome\"\n";
    print STDERR "\n";
	exit;
}
my $genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'} . "/";
my $start = -50;
my $end = 200;
my $gcStart = 0;
my $gcEnd = 150;
my $gcFlag = 1;
my $freqFlag = 1;
my $keepFlag = 0;
my $p3Flag = 0;
my $dFlag = '';
my $prefix = "";

for (my $i=2;$i<@ARGV; $i++) {
	if ($ARGV[$i] eq '-start') {
		$start = $ARGV[++$i];
		print STDERR "\tNew Start Position is $start\n";
	} elsif ($ARGV[$i] eq '-end') {
		$end = $ARGV[++$i];
		print STDERR "\tNew End Position is $end\n";
	} elsif ($ARGV[$i] eq '-prefix') {
		$prefix = $ARGV[++$i];
		print STDERR "\tUsing $prefix as a prefix file name\n";
	} elsif ($ARGV[$i] eq '-gcstart') {
		$gcStart = $ARGV[++$i];
		if ($gcStart < $start) {
			$start = $gcStart;
		}
		print STDERR "\tNew GC% Start Position is $gcStart\n";
	} elsif ($ARGV[$i] eq '-gcend') {
		$gcEnd = $ARGV[++$i];
		print STDERR "\tNew GC% End Position is $gcEnd\n";
		if ($gcEnd > $end) {
			$end = $gcEnd;
		}
	} elsif ($ARGV[$i] eq '-skipGC') {
		$gcFlag = 0;
		print STDERR "\tSkipping GC% calculation\n";
	} elsif ($ARGV[$i] eq '-skipFreq') {
		$freqFlag = 0;
		print STDERR "\tSkipping nucleotide frequency calculation\n";
	} elsif ($ARGV[$i] eq '-3p') {
		$p3Flag = 1;
		print STDERR "\tAlignment will be at the 3' end of the tags\n";
	} elsif ($ARGV[$i] eq '-keep') {
		$keepFlag = 1;
		print STDERR "\tSequence files will not be deleted\n";
	} elsif ($ARGV[$i] eq '-pos') {
		$dFlag = '+';
		print STDERR "\tOnly using + strand tags\n";
	} elsif ($ARGV[$i] eq '-neg') {
		$dFlag = '-';
		print STDERR "\tOnly using - strand tags\n";
	} else {
		print STDERR "!!! $ARGV[$i] not recognized!\n\n";
		cmdLineOptions();
		exit;
	}
		
}

my $tmpID = rand();
my $tmpFile = $tmpID . ".tmp";


my @tagFiles = ();
#print STDERR "\tWill Process the following tag files:\n";
`ls -1 "$directory"/*.tags.tsv > $tmpFile`;
open IN, $tmpFile;
while (<IN>) {
	chomp;
	my @line = split /\t/;
	push(@tagFiles, $_);
	#print STDERR "\t\t$_\n";
}
close IN;
`rm $tmpFile`;
@tagFiles = sort {$a cmp $b} @tagFiles;

my $seqFileStr = '';
my $posFileStr = '';
foreach(@tagFiles) {
	my $tagFile = $_;
	print STDERR "\tProcessing $tagFile\n";
	my $posFile = $tagFile . ".pos";
	my $seqFile = $tagFile . ".seq";
	my $p3opt = "";
	if ($p3Flag) {
		$p3opt = " -3p";
	}
	`tag2pos.pl "$tagFile" $start $end $p3opt > "$posFile"`;
	if ($dFlag ne '') {
		open IN, $posFile;
		open OUT, ">$tmpFile";
		while (<IN>) {
			my $og = $_;
			chomp;
			my @line = split /\t/;
			if ($line[4] eq '0' && $dFlag eq '+') {
				print OUT $og;
			} elsif ($line[4] eq '1' && $dFlag eq '-') {
				print OUT $og;
			}
		}
		close IN;
		close OUT;
		`mv $tmpFile $posFile`;
	}
			
	`homerTools extract "$posFile" "$genomeDir" $maskFlag > "$seqFile"`;
	$seqFileStr .= " \"$seqFile\"";
	$posFileStr .= " \"$posFile\"";
}

my $allSeqFile = $directory . "/allSeq.tsv";
my $freqFile = $directory . "/" . $outputFileName;
my $gcFile = $directory . "/" . $outputGCFileName;
my $cpgFile = $directory . "/" . $outputCpGFileName;
if ($prefix ne '') {
	print STDERR "\tUsing prefix \"$prefix\" for output filenames\n";
	$freqFile = $directory . $outputFileName;
	$gcFile = $directory . $outputGCFileName;
	$cpgFile = $directory . $outputCpGFileName;
}

`cat $seqFileStr > \"$allSeqFile\"`;

if ($freqFlag) {
	`homerTools freq \"$allSeqFile\" -offset $start > \"$freqFile\"`;
	print STDERR "\tNucleotide Frequency Results can be found in $directory/$outputFileName\n";
}
if ($gcFlag) {
	my $startoffset = $gcStart - $start;
	my $endoffset = $gcEnd - $start;
	`homerTools freq \"$allSeqFile\" -gc "$tmpFile" > /dev/null`;
	#`homerTools freq \"$allSeqFile\" -offset $startoffset $endoffset > $tmpFile`;
	my @gcBins = ();
	my @cpgBins = ();
	my $gcBinSize = 1/($gcEnd-$gcStart);
	my $cpgBinSize = 1/($gcEnd-$gcStart);
	my $halfGCBinSize = $gcBinSize/3;
	my $halfCpGBinSize = $cpgBinSize/3;
	my $gcBinMAX = ($gcEnd-$gcStart);
	my $cpgBinMAX = ($gcEnd-$gcStart);
	for (my $i=0;$i<=$gcBinMAX;$i++) {
		push(@gcBins, 0);
	}
	for (my $i=0;$i<=$cpgBinMAX;$i++) {
		push(@cpgBins, 0);
	}
	open IN, $tmpFile;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $cpgIndex = floor($line[1]/$cpgBinSize+$halfCpGBinSize);
		my $gcIndex = floor($line[2]/$gcBinSize+$halfGCBinSize);
		$cpgIndex = 0 if ($cpgIndex < 0);
		$gcIndex = 0 if ($gcIndex < 0);
		$cpgIndex = $cpgBinMAX if ($cpgIndex > $cpgBinMAX);
		$gcIndex = $gcBinMAX if ($gcIndex > $gcBinMAX);
		$cpgBins[$cpgIndex]++;
		$gcBins[$gcIndex]++;
	}
	close IN;
	`rm $tmpFile`;

	my $avgGC = 0;
	my $avgGCN = 0;
	open OUT, ">$gcFile";
	for (my $i=0;$i<=$gcBinMAX;$i++) {
		my $value = $i*$gcBinSize;
		print OUT "$value\t$gcBins[$i]\n";
		$avgGC+=$gcBins[$i]*$value;
		$avgGCN+=$gcBins[$i];
	}
	close OUT;
	$avgGC /= $avgGCN if ($avgGCN > 0);	
	print STDERR "\tGC fragment distribution can be found at $gcFile (Avg GC%=$avgGC)\n";
	my $avgCpG = 0;
	my $avgCpGN = 0;
	open OUT, ">$cpgFile";
	for (my $i=0;$i<=$cpgBinMAX;$i++) {
		my $value = $i*$cpgBinSize;
		print OUT "$value\t$cpgBins[$i]\n";
		$avgCpG+=$cpgBins[$i]*$value;
		$avgCpGN+=$cpgBins[$i];
	}
	close OUT;
	$avgCpG /= $avgCpGN if ($avgCpGN > 0);	
	print STDERR "\tCpG fragment distribution can be found at $cpgFile (Avg CpG%=$avgCpG)\n";
}

`rm $posFileStr`;
`rm "$allSeqFile"`;
if ($keepFlag == 0) {
	`rm $seqFileStr`;
}

exit;
