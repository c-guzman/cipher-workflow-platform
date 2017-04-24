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

my $config = HomerConfig::loadConfigFile();

my $wwwDir = "";
my $httpDir = "";
if (exists($config->{'SETTINGS'}->{'hubsUrl'})) {
    $httpDir = $config->{'SETTINGS'}->{'hubsUrl'};
}
if (exists($config->{'SETTINGS'}->{'hubsDir'})) {
    $wwwDir = $config->{'SETTINGS'}->{'hubsDir'};
}



my $check = `which bedGraphToBigWig`; 
if ($check eq "") {
	print STDERR "\n\t!!! Could not detect bedGraphToBigWig program in the executable path !!!\n";
	print STDERR "\t\tDownload it from http://hgdownload.cse.ucsc.edu/admin/exe/\n\n";
	exit;
}
	


if (@ARGV < 1) {
	printCMD();
}
	
sub printCMD {
	print STDERR "\n\tScript for automating the process of creating multiWig tracks\n";
	print STDERR "\n\tUsage: makeMultiWigHub.pl <hubname> <genome> [options] -d <tag directory1> [tag directory2]...\n";
	print STDERR "\n\tSpecial Options for bigWigs [choose one, don't combine]:\n";
	print STDERR "\t\t-normal (ChIP-Seq style, default)\n";
	print STDERR "\t\t-strand (Strand specific, for RNA-Seq and GRO-Seq)\n";
	print STDERR "\t\t-dnase (Special options for Crawford-lab style DNase-Seq)\n";
	print STDERR "\t\t-cage (Special options for CAGE/TSS-Seq)\n";
	print STDERR "\t\t-cpg (Special options for mCpG/CpG)\n";
	print STDERR "\n\tOther options:\n";
	print STDERR "\t\tWhatever options you want to pass to makeUCSCfile\n";
	print STDERR "\t\t!!Warning!!: do not try to specify \"-strand separate\" - use the special option above.\n";
	print STDERR "\t\tAlso, for the genome, do NOT use repeat version (mm9r) - use mm9 instead\n";
	print STDERR "\n\tFile options:\n";
	print STDERR "\t\t-force (overwrite existing hub)\n";
	print STDERR "\t\t-fsize <#> (limit the file size of the bigwig files to this value)\n";
	print STDERR "\t\t-chromSizes <chrom.size file> (specify the chromosome sizes, default: automatic)\n";
	print STDERR "\t\t\t-forceSizeCalc (force the creation of a new chromsome size file)\n";
	print STDERR "\t\t-url <URL> (URL directory -no filename- to tell UCSC where to look)\n";
	print STDERR "\t\t-webdir <directory> (name of directory to place resulting hub directory)\n";
	print STDERR "\n\tCurrent url target (-url):		 $httpDir\n";
	print STDERR "\tCurrent web directory (-webDir):   $wwwDir\n";
	print STDERR "\n\tYou're going to want to modify the \$wwwDir and \$httpDir variables at the top of\n";
	print STDERR "\tthe makeMultiWigHub.pl program file to accomidate your system so you don't have to\n";
	print STDERR "\tspecify -url and -webdir all the time.\n\n";
	exit;
}

if (@ARGV < 4) {
	print STDERR "!!! Too few arguments !!!\n";
	printCMD();
}

my $dir = $ARGV[0];
$dir =~ s/\/$//;

my $genome = $ARGV[1];

my $updateFlag = 0;
my $opt = "";
my $strandFlag = 0;
my $cpgFlag = 0;
my $forceFlag = 0;
my $fsize = 1e20;
my $chromSizeFile = "";
my $forceCalc = 0;

my @tagDirs = ();

for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-strand') {
		$strandFlag = 1;
	} elsif ($ARGV[$i] eq '-dnase') {
		$opt .= ' -fragLength 80 -adjust -40';
	} elsif ($ARGV[$i] eq '-force') {
		$forceFlag = 1;
	} elsif ($ARGV[$i] eq '-fsize') {
		$fsize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpg') {
		$cpgFlag = 1;
		$opt .= " -fragLength 1 -normLength 0";
	} elsif ($ARGV[$i] eq '-cage') {
		$strandFlag = 1;
		$opt .= " -fragLength 1";
	} elsif ($ARGV[$i] eq '-update') {
		$updateFlag = 1;
	} elsif ($ARGV[$i] eq '-forceSizeCalc') {
		$forceCalc = 1;
	} elsif ($ARGV[$i] eq '-chromSizes') {
		$chromSizeFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-webdir') {
		$wwwDir = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-d') {
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@tagDirs, $ARGV[$i]);
			#print STDERR "\t\t$ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
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

my $customGenomeFlag = 0;
if (-d $homeDir . "/data/genomes/$genome/") {
} elsif (-f $genome) {
	$customGenomeFlag = 1;
	print STDERR "Using custom genome (i.e. fasta file)\n";
}
	

my $rand = sprintf("%d",rand()*1e7);
my $tmpFile = $rand . ".tmp";

if ($chromSizeFile eq '') {
    $chromSizeFile = $homeDir . "/data/genomes/$genome/chrom.sizes";
	if ($customGenomeFlag) {
		$chromSizeFile = "$genome.chrom.sizes";
	}
    if (-f $chromSizeFile && $forceCalc == 0) {
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




my $colorIndex = 0;
my @color = ();
my @rcolor = ();
$color[0] = "255,150,150";
$color[1] = "150,150,255";
$color[2] = "150,255,150";
$color[3] = "255,200,150";
$color[4] = "150,255,220";
$color[5] = "200,150,255";
$color[6] = "200,200,150";
$color[7] = "150,200,200";
$color[8] = "200,150,200";
$rcolor[0] = "255,180,180";
$rcolor[1] = "180,180,255";
$rcolor[2] = "180,255,180";
$rcolor[3] = "255,210,180";
$rcolor[4] = "180,255,240";
$rcolor[5] = "210,180,255";
$rcolor[6] = "210,210,170";
$rcolor[7] = "170,210,210";
$rcolor[8] = "210,170,210";
foreach(@tagDirs) {
	push(@color, floor(rand()*255). "," . floor(rand()*255) . "," . floor(rand()*255));
	push(@rcolor, floor(rand()*255). "," . floor(rand()*255) . "," . floor(rand()*255));
}

my $hubGenome = $genome;
if ($customGenomeFlag) {
	$hubGenome =~ s/^.*\///g;
	$hubGenome =~ s/\.fasta//g;
	$hubGenome =~ s/\.fa//g;
	print STDERR "Hub genome set to $hubGenome\n";
}

my $upload = $httpDir . "/$dir/hub.txt";
my $uploadW = $httpDir . "/$dir/washU.hub.txt";
my $washuBase = $httpDir . "/$dir/$hubGenome/";

print STDERR "\n\tOnce finished, you will want to upload the following hub URL:\n";
print STDERR "\t\t$upload\n\n";
print STDERR "\tIf loading to the Wash U Epigenome Browser, use:\n";
print STDERR "\t\t$uploadW\n\n";



my $hubDir = $wwwDir . "/" .  $dir;
if (-d $hubDir) {
	if ($forceFlag) {
		print STDERR "\tOverwriting contents in $hubDir\n";
	} else {
		print STDERR "!!! Warning - $hubDir exists !!!\n";
		print STDERR "!!! Rerun command with -force to overwrite !!!\n";
		exit;
	}
}

`mkdir -p "$hubDir"`;

$user = `whoami`;
$user =~ s/\n//g;

open HUB, ">$hubDir/hub.txt";
print HUB "hub $dir\n";
print HUB "shortLabel $dir\n";
print HUB "longLabel $dir\n";
print HUB "genomesFile genomes.txt\n";
print HUB "email $user" . '@ucsd.edu' . "\n";
close HUB;

open WASHU, ">$hubDir/washU.hub.txt";
print WASHU "# This is a hub file for the Wash U Epigenome Browser\n";


open GENOMES, ">$hubDir/genomes.txt";
print GENOMES "genome $hubGenome\n";
print GENOMES "trackDb $hubGenome/trackDb.txt\n";
close GENOMES;


`mkdir -p "$hubDir/$hubGenome"`;
open TRACKDB, ">$hubDir/$hubGenome/trackDb.txt";
print TRACKDB "track $dir\n";
print TRACKDB "container multiWig\n";
print TRACKDB "noInherit on\n";
print TRACKDB "shortLabel $dir\n";
print TRACKDB "longLabel $dir\n";
print TRACKDB "type bigWig\n";
print TRACKDB "configurable on\n";
print TRACKDB "visibility full\n";
print TRACKDB "aggregate transparentOverlay\n";
print TRACKDB "showSubtrackColorOnUi on\n";
print TRACKDB "autoScale on\n";
if (0) { #if ($cpgFlag) {
	print TRACKDB "windowingFunction average\n";
} else {
	print TRACKDB "windowingFunction maximum\n";
}
print TRACKDB "priority 1.4\n";
print TRACKDB "alwaysZero on\n";
print TRACKDB "yLineMark 0\n";
print TRACKDB "yLineOnOff on\n";

if ($strandFlag == 1) {
	print TRACKDB "maxHeightPixels 125:125:11\n";
} else {
	print TRACKDB "maxHeightPixels 100:75:11\n";
}
print TRACKDB "\n";

foreach(@tagDirs) {
	my $tagDir = $_;
	$tagDir =~ s/\/$//;

	my $tagDirName = $tagDir;
	$tagDirName =~ s/\/+$//g;
	$tagDirName =~ s/.*\///;
	my $trackName = $tagDirName;
	$trackName =~ s/\///g;
	#print STDERR "trackName=$trackName\n";

	my $hubListFile = $tagDir . "/availableHubs.txt";
	if (-e $hubListFile) {
		open HUBLIST, ">>$hubListFile";
	} else {
		open HUBLIST, ">$hubListFile";
		print HUBLIST "Hubs containing data from $trackName\tLocation of Hub on server\n";
	}
	print HUBLIST "$upload\t$hubDir\t\n";
	close HUBLIST;

	if ($strandFlag == 0 && $cpgFlag == 0) {
		`makeUCSCfile "$tagDir" -o auto -fsize $fsize -bigWig "$chromSizeFile" $opt > /dev/null`;

		my $bigWigFile = $tagDirName . ".ucsc.bigWig";
		my $targetfile = "$hubDir/$hubGenome/$bigWigFile";
		#print STDERR "`mv $tagDir/$bigWigFile $targetfile`;\n";
		`mv "$tagDir/$bigWigFile" "$targetfile"`;


		print TRACKDB "track $trackName\n";
		print TRACKDB "bigDataUrl $bigWigFile\n";
		print TRACKDB "shortLabel $trackName\n";
		print TRACKDB "longLabel $trackName\n";
		print TRACKDB "type bigWig\n";
		print TRACKDB "parent $dir\n";
		print TRACKDB "color $color[$colorIndex++]\n";
		print TRACKDB "\n";

		print WASHU "bigwig\t$washuBase" . $bigWigFile 
					. "\t$trackName\tshow\tcolorpositive:#990000,colornegative:#004080\n";

	} else {

		my $fwdopt = $opt;
		if ($cpgFlag == 1) {
			$fwdopt .= " -style methylated";
		} else {
			$fwdopt .= " -strand +";
		}

		`makeUCSCfile "$tagDir" -o auto -fsize $fsize -bigWig "$chromSizeFile" $fwdopt > /dev/null`;

		$bigWigFile = $tagDirName . ".ucsc.bigWig";
		$bigWigPosFile = $tagDirName . "pos.ucsc.bigWig";

		my $targetfilePos = "$hubDir/$hubGenome/$bigWigPosFile";
		#print STDERR "`mv $tagDir/$bigWigFile $targetfilePos`;\n";
		`mv "$tagDir/$bigWigFile" "$targetfilePos"`;

		my $negopt = $opt;
		if ($cpgFlag == 1) {
			#$negopt .= " -style methylated -mintbp 1";
			$negopt .= " -style unmethylated";
		} else {
			$negopt .= " -strand - -neg";
		}

		`makeUCSCfile "$tagDir" -o auto -fsize $fsize -bigWig "$chromSizeFile" $negopt > /dev/null`;
		$bigWigFile = $tagDirName . ".ucsc.bigWig";
		$bigWigNegFile = $tagDirName . "neg.ucsc.bigWig";

		my $targetfileNeg = "$hubDir/$hubGenome/$bigWigNegFile";
		#print STDERR "`mv $tagDir/$bigWigFile $targetfileNeg`;\n";
		`mv "$tagDir/$bigWigFile" "$targetfileNeg"`;

		$tagDirName = $tagDir;
		#$tagDirName =~ s/\s//g;

		my $pname = $trackName . "+";
		print TRACKDB "track $pname\n";
		print TRACKDB "bigDataUrl $bigWigPosFile\n";
		print TRACKDB "shortLabel $pname\n";
		print TRACKDB "longLabel $pname\n";
		print TRACKDB "type bigWig\n";
		print TRACKDB "parent $dir\n";
		if (0) { #if ($cpgFlag) {
			print TRACKDB "color 0,0,50\n";
		} else {
			print TRACKDB "color $color[$colorIndex]\n";
		}
		print TRACKDB "\n";

		print WASHU "bigwig\t$washuBase" . $bigWigPosFile 
					. "\t$pname\tshow\tcolorpositive:#990000,colornegative:#004080\n";

		$pname = $trackName . "-";
		print TRACKDB "track $pname\n";
		print TRACKDB "bigDataUrl $bigWigNegFile\n";
		print TRACKDB "shortLabel $pname\n";
		print TRACKDB "longLabel $pname\n";
		print TRACKDB "type bigWig\n";
		print TRACKDB "parent $dir\n";
		if (0) { #if ($cpgFlag) {
			print TRACKDB "color 220,220,220\n";
		} else {
			print TRACKDB "color $rcolor[$colorIndex++]\n";
		}
		print TRACKDB "\n";

		print WASHU "bigwig\t$washuBase" . $bigWigNegFile 
					. "\t$pname\tshow\tcolorpositive:#990000,colornegative:#004080\n";

	}
}


close TRACKDB;

print STDERR "\n\tAll finished: you will want to upload the following hub URL:\n";
print STDERR "\t\t$upload\n\n";
print STDERR "\tIf loading to the Wash U Epigenome Browser, use:\n";
print STDERR "\t\t$uploadW\n\n";
print STDERR "\n\tTo edit track settings, edit files in:\n\t\t$hubDir\n\n";

exit;

