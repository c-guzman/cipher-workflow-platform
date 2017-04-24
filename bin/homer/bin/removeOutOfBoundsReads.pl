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




if (@ARGV < 2) {
	printCMD();
}
	
sub printCMD {
	print STDERR "\n\tScript to automatically remove reads mapping outside of chromosome coordinates\n";
	print STDERR "\n\tUsage: removeOutOfBoundsReads.pl <tag directory> <genome> [options]\n";
	print STDERR "\n\tAltUsage: removeOutOfBoundsReads.pl <peak/BED file> <genome> [options] > newPeakFile.txt\n";
	print STDERR "\n\tOtions:\n";
	print STDERR "\t\t-chromSizes <chrom.size file> (specify the chromosome sizes, default: automatic)\n";
	print STDERR "\t\t-force (force calculation of chromsome sizes for genome FASTA files)\n";
	print STDERR "\n\t\tchrom.size file is tab delimited: chr<tab>size\n";
	print STDERR "\n\t\tFor custom genomes, use removeOutOfBoundsReads.pl <directory> none -chromSizes <chrom.size file>\n";
	print STDERR "\n";
	exit;
}

my $dir = $ARGV[0];
my $genome = $ARGV[1];
$dir =~ s/\/$//;

my $fileFlag = 0;
if (-d $dir) {
	$fileFlag = 0;
} elsif (-e $dir) {
	$fileFlag = 1;
	print STDERR "\tTreating input as a peak/BED file...\n";
}

my $chromSizeFile = "";
my $forceCalc = 0;

for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-chromSizes') {
		$chromSizeFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-force') {
		$forceCalc = 1;
	} elsif ($ARGV[$i] eq '-help') {
		printCMD();
	} elsif ($ARGV[$i] eq '--help') {
		printCMD();
	} else {
		printCMD();
	}
}
my $rand = sprintf("%d",rand()*1e7);
my $tmpFile = $rand . ".tmp";

if ($chromSizeFile eq '') {
	$chromSizeFile = $homeDir . "/data/genomes/$genome/chrom.sizes";
	if (-e $chromSizeFile && $forceCalc == 0) {
	} else {
		print STDERR "\t!! Could not find chromosome size file ($chromSizeFile)\n";
		print STDERR "\t   Generating one from genome sequence files...\n";
		`homerTools extract stats \"$homeDir/data/genomes/$genome/\" > "$tmpFile"`;
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

my %chrSizes = ();
open IN, $chromSizeFile or die "!!! Could not open a valid chromosome sizes file !!!\ntried: $chromSizeFile)\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$chrSizes{$line[0]}=$line[1];
}
close IN;

if ($fileFlag) {
	my $totalLines = 0;
	my $totalBad = 0;
	open IN, $dir;
	while (<IN>) {
		chomp;
		if (/^#/) {
			print "$_\n";
			next;
		}
		my @line = split /\t/;
		my $cool = 0;
		if (exists($chrSizes{$line[0]}) && $line[1]=~/\d+/ && $line[2] =~ /\d+/) {
			#bed file
			$totalLines++;
			next if ($line[1] > $chrSizes{$line[0]});
			if ($line[2] > $chrSizes{$line[0]}) {
				$totalBad++;
				$line[2] = $chrSizes{$line[0]};
			}
			$cool = 1;
		} elsif (exists($chrSizes{$line[1]}) && $line[2]=~/\d+/ && $line[3] =~ /\d+/) {
			#peak file
			$totalLines++;
			next if ($line[2] > $chrSizes{$line[1]});
			if ($line[3] > $chrSizes{$line[1]}) {
				$totalBad++;
				$line[3] = $chrSizes{$line[1]};
			}
			$cool = 1;
		}
		if ($cool) {
			print "$line[0]";
			for (my $i=1;$i<@line;$i++){ 
				print "\t$line[$i]";
			}
			print "\n";
		}
	}
	print STDERR "\tTotal of $totalBad of $totalLines were outside of range\n\n";

	exit;
}





`ls -1 "$dir/"*tags.tsv > "$tmpFile"`;
open IN, $tmpFile;
my %chrFiles = ();
while (<IN>) {
	chomp;
	my $f = $_;
	s/^.*\///;
	s/\.tags\.tsv//;
	my $c = $_;
	$chrFiles{$c} = $f;
}
close IN;
`rm "$tmpFile"`;


my @chr = sort {$a cmp $b} keys %chrFiles;
foreach(@chr) {
	my $chr = $_;
	my $file = $chrFiles{$chr};
	print STDERR "\tProcessing $chr:";

	if (!exists($chrSizes{$chr})) {
		print STDERR " !! No limit for $chr !! - removing...\n";
		`rm "$file"`;
		next;
	}
	my $limit = $chrSizes{$chr};

	my $totalPos = 0;
	my $badPos = 0;
	open OUT, ">$tmpFile";
	open IN, $file or die "!!! Could not open file ($file) for $chr !!!\n";
	while (<IN>) {
		my $og = $_;
		chomp;
		my @line = split /\t/;
		$totalPos++;
		if ($line[2] > $limit || $line[2] < 1) {
			$badPos++;
		} else {
			print OUT $og;
		}
	}
	close OUT;
	close IN;
	`mv "$tmpFile" "$file"`;

	my $DD = $totalPos;
	$DD = 1 if ($DD < 1);
	my $percentBad = $badPos/$DD*100;
	my $v = sprintf("%.3lf",$percentBad);
	print STDERR " $badPos of $totalPos positions ($v%) out of bounds\n";
}
