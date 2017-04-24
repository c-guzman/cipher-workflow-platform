#!/usr/bin/env perl
use warnings;
use lib "/sanger/nfs_data/chris/genomics/homer/.//bin";
my $homeDir = "/sanger/nfs_data/chris/genomics/homer/./";


# Copyright 2009 - 2014 Christopher Benner <cbenner@salk.edu>
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
use Statistics;

my $gbin = 50;
my $bin = 100;
my $size = 5000;
my $minGeneSize = 3000;
my $maxGeneSize = -1;

my $config = HomerConfig::loadConfigFile();

sub printCMD {

	print STDERR "\n\tUsage: makeMetaGeneProfile.pl <peak file> <genome> [options]\n";
	print STDERR   "\t       makeMetaGeneProfile.pl rna <genome> [options]\n";
	print STDERR "\n\tAvailable Genomes (required argument): (name,org,directory,default promoter set)\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}"
				. "\t$config->{'GENOMES'}->{$_}->{'promoters'}\n";
	}
	print STDERR "\t\t\t-- or --\n";
	print STDERR "\t\tCustom: provide the path to genome FASTA files (directory or single file)\n";
	print STDERR "\t\tIf no genome is available, specify 'none'.\n";

	print STDERR "\n\tOptions controlling meta-gene creation:\n";
	print STDERR "\t\t-min <#> (minimum size of genes/regions to use, default: $minGeneSize)\n";
	print STDERR "\t\t\t(This prevents extremely small regions from being used)\n";
	print STDERR "\t\t-max <#> (max size of genes/regions to use, default: no max)\n";
	print STDERR "\t\t-gbin <#> (Number of bins in gene body, default: $gbin)\n";
	print STDERR "\t\t-gRatio <#> (Ratio of gene region to flanks, default: 2)\n";
	print STDERR "\t\t-bin <#> (Bin size for 5' and 3' flanks, default: $bin)\n";
	print STDERR "\t\t-size <#> (Size of 5' and 3' flanks, default: $size)\n";
	print STDERR "\n\tAll other options are passed to annotatePeaks.pl.  For example, to see the read density\n";
	print STDERR "\tfrom a tag directory, add \"-d <tagDir>\", or for peak density, use \"-p <peakfile>\"\n";
	print STDERR "\t\n";
	print STDERR "\n";
	exit;
}


if (@ARGV < 2) {
	printCMD();
}
my $options = "";
my $peakFile = $ARGV[0];
my $givenGenome = $ARGV[1];
my $gratio = 2;
for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-gbin') {
		$gbin = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bin') {
		$bin = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-size') {
		$size = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-gRatio') {
		$gratio = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-max') {
		$maxGeneSize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$minGeneSize = $ARGV[++$i];
	} else {
		$options .= " $ARGV[$i]";
	}
	
}

my $genome = $givenGenome;
my $maskFlag = 1;
if ($genome =~ s/\r$//) {
    $maskFlag = 1;
}

my $genomeDir = "";
my $organism = "unknown";
if (!exists($config->{'GENOMES'}->{$genome})) {
    if ($genome eq 'none') {
        $organism = 'unknown';
    } else {
        print STDERR "\tCan't find genome \"$genome\"...\n";
        print STDERR "\t\tAssuming $genome is the name of the orgaism (or use \"none\")\n";
        $organism = $genome;
        $genome = 'none';
    }  
} else {
    $genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'};
    $organism = $config->{'GENOMES'}->{$genome}->{'org'};
}

if (-f $peakFile) {
} else {
	my $input = $peakFile;
    if ($input eq 'rna') {
        if ($genomeDir eq '') {
            print STDERR "!!! Can't use \"repeats\", \"rna\", or \"miRNA\" as input file if you don't use a valid homer genome\n";  
            exit;
        }
        $peakFile = $genomeDir . "/" . $genome . "." . $input;
        if (! -f $peakFile) {
            print STDERR "!!! Cannot find file \"$peakFile\" - likely a configuration error\n";
        }  
    } else {
        print STDERR "!!! What is input \"$input\"?? Not recognized as a file (check the path).\n";
        exit;
    }
}




my $offset = $size*$gratio;
$size *= 2;


my $rand = rand();
my $tmpFile1 = $rand . ".1.tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";
my $tmpFile4 = $rand . ".4.tmp";
my $tmpFile5 = $rand . ".5.tmp";
my $tmpFile6 = $rand . ".6.tmp";

my $opt = "";
if ($maxGeneSize > 0) {
	$opt =  " -max $maxGeneSize";
}
`adjustPeakFile.pl "$peakFile" -min $minGeneSize $opt > "$tmpFile1"`;
`adjustPeakFile.pl "$tmpFile1" -5p > "$tmpFile2"`;
`adjustPeakFile.pl "$tmpFile1" -3p > "$tmpFile3"`;

`annotatePeaks.pl "$tmpFile1" $givenGenome -size given -hist $gbin $options > "$tmpFile4"`;
`annotatePeaks.pl "$tmpFile2" $givenGenome -size $size -hist $bin $options > "$tmpFile5"`;
`annotatePeaks.pl "$tmpFile3" $givenGenome -size $size -hist $bin $options > "$tmpFile6"`;

open IN, $tmpFile5;
my $c = 0;
while (<IN>) {
	chomp;
	s/\r//g;
	my $og = $_;
	$c++;
	if ($c == 1) {
		print "$og\n";
		next;
	}
	my @line = split /\t/;
	if ($line[0] <= 0) {
		print "$og\n";
	} else {
		last;
	}
}
close IN;


open IN, $tmpFile4;
$c = 0;
while (<IN>) {
	chomp;
	s/\r//g;
	my $og = $_;
	$c++;
	if ($c == 1) {
		next;
	}
	my @line = split /\t/;
	$line[0] *= $offset;
	print "$line[0]";
	for (my $i=1;$i<@line;$i++) {
		print "\t$line[$i]";
	}
	print "\n";
}
close IN;

open IN, $tmpFile6;
$c = 0;
while (<IN>) {
	chomp;
	s/\r//g;
	my $og = $_;
	$c++;
	if ($c == 1) {
		next;
	}
	my @line = split /\t/;
	if ($line[0] <= 0) {
		next;
	}
	$line[0] += $offset;
	print "$line[0]";
	for (my $i=1;$i<@line;$i++) {
		print "\t$line[$i]";
	}
	print "\n";
}
close IN;
`rm "$tmpFile1" "$tmpFile2" "$tmpFile3" "$tmpFile4" "$tmpFile5" "$tmpFile6"`;
