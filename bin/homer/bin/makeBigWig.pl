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


use HomerConfig;
my $config = HomerConfig::loadConfigFile();

my $wwwDir = "";
my $httpDir = "";
if (exists($config->{'SETTINGS'}->{'bigWigUrl'})) {
	$httpDir = $config->{'SETTINGS'}->{'bigWigUrl'};
}
if (exists($config->{'SETTINGS'}->{'bigWigDir'})) {
	$wwwDir = $config->{'SETTINGS'}->{'bigWigDir'};
}


my $check = `which bedGraphToBigWig`; 
if ($check eq "") {
	print STDERR "\n\t!!! Could not detect bedGraphToBigWig program in the executable path !!!\n";
	print STDERR "\t\tDownload it from http://hgdownload.cse.ucsc.edu/admin/exe/\n\n";
	exit;
}


if (@ARGV < 2) {
	if (@ARGV == 1) {
		print STDERR "!! makeBigWig.pl now requires you enter the genome !!\n";
	}
	printCMD();
}
	
sub printCMD {
	print STDERR "\n\tScript for automating the process of creating bigWigs\n";
	print STDERR "\n\tUsage: makeBigWig.pl <tag directory> <genome> [special options] [options]\n";
	print STDERR "\n\tSpecial Options for bigWigs [choose one, don't combine]:\n";
	print STDERR "\t\t-normal (ChIP-Seq style, default)\n";
	print STDERR "\t\t-strand (Strand specific, for RNA-Seq and GRO-Seq)\n";
	print STDERR "\t\t-dnase (Special options for Crawford-lab style DNase-Seq)\n";
	print STDERR "\t\t-cage (Special options for CAGE/TSS-Seq)\n";
	print STDERR "\t\t-cpg (Special options for mCpG/CpG)\n";
	print STDERR "\n\tOther options:\n";
	print STDERR "\t\tWhatever options you want to pass to makeUCSCfile\n";
	print STDERR "\t\t!!Warning!!: do not try to specify \"-strand separate\" - use the special option above.\n";
	print STDERR "\n\tFile options:\n";
	print STDERR "\t\t-fsize <#> (Use to limit the size of the bigwig files)\n";
	print STDERR "\t\t-chromSizes <chrom.size file> (specify the chromosome sizes, default: automatic)\n";
	print STDERR "\t\t\t-forceSizeCalc (force the creation of a new chromsome size file)\n";
	print STDERR "\t\t-url <URL> (URL directory -no filename- to tell UCSC where to look)\n";
	print STDERR "\t\t-webdir <directory> (name of directory to place resulting bigWig file)\n";
	print STDERR "\t\t-update or -force (overwrite bigwigs in the webDir directory, otherwise random numbers are\n";
	print STDERR "\t\t\tadded to make the file unique.\n";
	print STDERR "\n\tCurrent url target (-url):         $httpDir\n";
	print STDERR "\tCurrent web directory (-webDir):   $wwwDir\n";
	print STDERR "\n\tYou're going to want to modify the \$wwwDir and \$httpDir variables at the top of\n";
	print STDERR "\tthe makeBigWig.pl program file to accomidate your system so you don't have to\n";
	print STDERR "\tspecify -url and -webdir all the time.\n\n";
	exit;
}

my $dir = $ARGV[0];
my $genome = $ARGV[1];
$dir =~ s/\/$//;

my $updateFlag = 0;
my $opt = "";
my $strandFlag = 0;
my $cpgFlag = 0;
my $fsize = 1e20;
my $chromSizeFile = "";
my $forceCalc = 0;
my @washU = ();

for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-strand') {
		$strandFlag = 1;
	} elsif ($ARGV[$i] eq '-dnase') {
		$opt .= ' -fragLength 80 -adjust -40';
	} elsif ($ARGV[$i] eq '-cpg') {
		$cpgFlag = 1;
		$opt .= " -fragLength 1 -normLength 0";
	} elsif ($ARGV[$i] eq '-cage') {
		$strandFlag = 1;
		$opt .= " -fragLength 1";
	} elsif ($ARGV[$i] eq '-normal') {
	} elsif ($ARGV[$i] eq '-fsize') {
		$fsize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-forceSizeCalc') {
		$forceCalc = 1;
	} elsif ($ARGV[$i] eq '-chromSizes') {
		$chromSizeFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-update' || $ARGV[$i] eq '-force') {
		$updateFlag = 1;
	} elsif ($ARGV[$i] eq '-webdir') {
		$wwwDir = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-url') {
		$httpDir = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-help') {
		printCMD();
	} elsif ($ARGV[$i] eq '--help') {
		printCMD();
	} else {
		$opt .= " $ARGV[$i]";
	}
}

if ($httpDir eq '') {
	print STDERR "!!! \"-url\" must be set to something (can also modify config.txt file)\n";
	exit;
}
if ($wwwDir eq '') {
	print STDERR "!!! \"-webdir\" must be set to something (can also modify config.txt file)\n";
	exit;
}
my $rand = sprintf("%d",rand()*1e7);
my $tmpFile = $rand . ".tmp";

if ($chromSizeFile eq '') {
	$chromSizeFile = $homeDir . "/data/genomes/$genome/chrom.sizes";
	if (-d $chromSizeFile && $forceCalc == 0) {
	} else {
		my $customGenomeFlag = 0;
		if (-f $genome) {
			$chromSizeFile = "$genome.chrom.sizes";
			$customGenomeFlag = 1;
		}
		if (-e $chromSizeFile && $forceCalc == 0) {
		} else {

			print STDERR "\t!! Could not find chromosome size file ($chromSizeFile)\n";
			print STDERR "\t   Generating one from genome sequence files...\n";
			if ($customGenomeFlag) {
				#print STDERR "`homerTools extract stats \"$genome\" > $tmpFile`;\n";
				`homerTools extract stats \"$genome\" > "$tmpFile"`;
			} else {
				`homerTools extract stats \"$homeDir/data/genomes/$genome/\" > "$tmpFile"`;
			}
			open IN, $tmpFile;
			open OUT, ">$chromSizeFile";
			while (<IN>) {
				chomp;
				my @line = split /\t/;
				next if ($line[0] eq 'ChrName');
				next if ($line[0] eq 'genome');
				print OUT "$line[0]\t$line[1]\n";
			}
			close IN;
			close OUT;
			`rm "$tmpFile"`;
		}
	}
}

my $topDir = $dir;
$topDir =~ s/\/+$//;
$topDir =~ s/^.*?\///;
#print STDERR "$topDir, $dir\n";

if ($strandFlag == 0 && $cpgFlag == 0) {


	`makeUCSCfile "$dir" -o auto -fsize $fsize -bigWig "$chromSizeFile" $opt > $tmpFile`;


	my $bigWigFile = $topDir . ".ucsc.bigWig";
	my $bigWigFileTrack = $bigWigFile;
	my $targetfile = $wwwDir . "/" . $bigWigFile;
	if (-e $targetfile && $updateFlag == 0) {
		print STDERR "Found file $targetfile - adding random number ($rand)\n";
		$targetfile = $targetfile . "." . $rand;
		$bigWigFileTrack = $bigWigFile . "." . $rand;
	}

	open IN, "$tmpFile";
	open OUT, ">$dir/$topDir.ucsc.bigWig.track.txt";
	while (<IN>) {
		chomp;
		my $start = $_;
		my $end = $_;
		$start =~ s/bigDataUrl=.*$//;
		$end =~ s/^.*?in\)//;
		print OUT $start . "bigDataUrl=$httpDir" . $bigWigFileTrack . $end . "\n";
		push(@washU, $httpDir . $bigWigFileTrack);
		last;
	}
	close IN;
	close OUT;

	`rm $tmpFile`;
	`mv $dir/$bigWigFile $targetfile`;
} else {
	my $fwdopt = $opt;
	if ($cpgFlag == 1) {
		$fwdopt .= " -style unmethylated";
	} else {
		$fwdopt .= " -strand +";
	}

	`makeUCSCfile "$dir" -o auto -fsize $fsize -bigWig "$chromSizeFile" $fwdopt > $tmpFile`;
	$bigWigFile = $topDir . ".ucsc.bigWig";
	$bigWigPosFile = $topDir . "pos.ucsc.bigWig";

	my $bigWigFileTrackPos = $bigWigPosFile;
	my $targetfilePos = $wwwDir . "/" . $bigWigPosFile;
	if (-e $targetfilePos && $updateFlag == 0) {
		print STDERR "Found file $targetfilePos - adding random number ($rand)\n";
		$targetfilePos = $targetfilePos . "." . $rand;
		$bigWigFileTrackPos = $bigWigPosFile . "." . $rand;
	}

	#print STDERR "`mv $dir/$bigWigFile $targetfilePos`;\n";
	`mv "$dir/$bigWigFile" "$targetfilePos"`;

	open IN, "$tmpFile";
	open OUT, ">$dir/$topDir.ucsc.bigWig.track.txt";
	while (<IN>) {
		chomp;
		my $start = $_;
		my $end = $_;
		$start =~ s/bigDataUrl=.*$//;
		$end =~ s/^.*?in\)//;
		print OUT $start . "bigDataUrl=$httpDir" . $bigWigFileTrackPos . $end . "\n";
		push(@washU, $httpDir . $bigWigFileTrackPos);
		last;
	}
	close IN;

	my $negopt = $opt;
	if ($cpgFlag == 1) {
		$negopt .= " -style methylated";
	} else {
		$negopt .= " -strand -";
	}


	`makeUCSCfile "$dir" -o auto -fsize $fsize -bigWig "$chromSizeFile" $negopt > $tmpFile`;
	$bigWigFile = $topDir . ".ucsc.bigWig";
	$bigWigNegFile = $topDir . "neg.ucsc.bigWig";

	my $bigWigFileTrackNeg = $bigWigNegFile;
	my $targetfileNeg = $wwwDir . "/" . $bigWigNegFile;
	if (-e $targetfileNeg && $updateFlag == 0) {
		print STDERR "Found file $targetfileNeg - adding random number ($rand)\n";
		$targetfileNeg = $targetfileNeg . "." . $rand;
		$bigWigFileTrackNeg = $bigWigNegFile . "." . $rand;
	}

	`mv "$dir/$bigWigFile" "$targetfileNeg"`;
	#print STDERR "`mv $dir/$bigWigFile $targetfileNeg`;\n";

	open IN, "$tmpFile";
	while (<IN>) {
		chomp;
		my $start = $_;
		my $end = $_;
		$start =~ s/bigDataUrl=.*$//;
		$end =~ s/^.*?in\)//;
		print OUT $start . "bigDataUrl=$httpDir" . $bigWigFileTrackNeg . $end . "\n";
		push(@washU, $httpDir . $bigWigFileTrackNeg);
		last;
	}
	close IN;
	close OUT;

	`rm $tmpFile`;
}

print STDERR "\n\tIf you want to upload these tracks to the WashU Epigenome Browser, use:\n";
foreach(@washU) {
	print STDERR "\t\t$_\n";
}
print STDERR "\n";


