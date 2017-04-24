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

if (@ARGV < 1) {
	cmdLineOptions();
}
sub cmdLineOptions {
	print STDERR "\n\tProgram preparse the genome for use with motif finding.\n";
	print STDERR "\n\tUsage: preparseGenome.pl <genome> -size <#> [additional options]\n";
	print STDERR "\tOutput files will be placed in $homeDir/data/genomes/[genome]/preparsed/\n";
	print STDERR "\t\t*.###.seq *.###.pos *.###.cgfreq *.###.cgbins *.###.gcbins\n";
	print STDERR "\tFirst argument must be <genome>\n";

	print STDERR "\n\tPossible Genomes:\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}\n";
	}
	print STDERR "\t\t\t-- or --\n";
	print STDERR "\t\tCustom: provide the path to genome FASTA files (directory or single file)\n";
	print STDERR "\t\t\tHeads up: will create the directory \"preparsed/\" in same location.\n";

	print STDERR "\n\tBasic options:\n";
	print STDERR "\t\t-size <#> (size of fragments to use for preparsing the genome)\n";
	print STDERR "\t\t-mask (mask repeats - i.e. lower-case bases)\n";
	print STDERR "\t\t-ref <peak file> (reference position file, default: [genome].tss)\n";
	print STDERR "\t\t\tIf no reference file is given or found, random regions will be used\n";
	print STDERR "\t\t\tTo force the use of random regions, use \"-ref random\"\n";
	print STDERR "\t\t-window <#> (size of window around ref positions to prepare, default=50000)\n";
	print STDERR "\t\t-max <#> (maximum number of preparesed fragments to create, default=2e6)\n";
	print STDERR "\t\t-minInc <#> (minimum size of region near ref pos to include, default=1000)\n";
	print STDERR "\t\t-preparsedDir <directory> (alternative directory to place the preparsed output files)\n";
	print STDERR "\n";
	exit;
}

my $genome = $ARGV[0];
my $customGenome = 0;
my $genomeDir = "";
my $genomeParseDir = "";
my $forceRandom = 0;
my $maskFlag = 0;

if ($genome =~ s/r$//) {
	$maskFlag = 1;
}

if (!exists($config->{'GENOMES'}->{$genome})) {
	($genome,$genomeDir,$genomeParseDir) = HomerConfig::parseCustomGenome($genome);
} else {
	$genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'} . "/";
	$genomeParseDir = $genomeDir . "preparsed/";
}
#print STDERR "\tPreparsed File will be put in: $genomeParseDir\n";
#print STDERR "\tgenome name: $genome\n";

my $defFile = $genomeDir . "$genome.tss";
if ($customGenome == 0) {
	if (-f $defFile) {
	} else {
		#print STDERR "!!! You should update the $genome genome using configureHomer.pl !!!\n";
		#$defFile = $genomeDir . "$genome.promoters";
	}
}


my $refFile = $defFile;
my $regionSize = 200;
my $window = 50000;
my $minInclude = 1000;
my $maxRegions = 2e6;
my $noRefFileFlag = 0;


for (my $i=1;$i<@ARGV; $i++) {
	if ($ARGV[$i] eq '-size') {
		$regionSize = $ARGV[++$i];
		print STDERR "\tpreparse size set to $regionSize\n";
	} elsif ($ARGV[$i] eq '-ref') {
		$refFile = $ARGV[++$i];
		if ($refFile eq 'random') {
			$forceRandom = 1;
			$noRefFileFlag = 1;
			print STDERR "\tUsing random positions as reference points\n";
		} else {
			print STDERR "\tPeak file $refFile will be used as reference positions\n";
		}
	} elsif ($ARGV[$i] eq '-preparsedDir') {
		$genomeParseDir = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-mask') {
		$maskFlag = 1;
	} elsif ($ARGV[$i] eq '-window') {
		$window = $ARGV[++$i];
		print STDERR "\tWindow size set to $window\n";
	} elsif ($ARGV[$i] eq '-max') {
		$maxRegions = $ARGV[++$i];
		print STDERR "\tMax number of preparsed regions set to $maxRegions\n";
	} elsif ($ARGV[$i] eq '-minInc') {
		$minInclude = $ARGV[++$i];
		print STDERR "\tA minimum of $minInclude will be used around reference positions\n";
	} else {
		print STDERR "!!! $ARGV[$i] not recognized as a valid option!\n";
		cmdLineOptions();
	}
}
$genomeParseDir .= "/";

my $regionSizeHalf = floor($regionSize/2);

print STDERR "\n";
if ($refFile eq $defFile) {
	unless (-e $refFile) {
		$noRefFileFlag = 1;
		print STDERR "\tNo reference file found - using random positions...\n";
	} else {
		print STDERR "\tBy default, using $defFile for reference positions\n";
	}
}
my $outputGenome = $genome;
if ($maskFlag) {
	$outputGenome .= "r";
}
print STDERR "\tOutput files will be placed in $genomeParseDir". "$outputGenome.$regionSize.* ...\n";
`mkdir -p "$genomeParseDir"`;

if (-w "$genomeParseDir") {
} else {
	print STDERR "!!! Warning: Looks like you do not have permission to write the preparsed\n";
	print STDERR "!!! genome files to the following directory:\n";
	print STDERR "!!!   $genomeParseDir\n";
	print STDERR "!!! Consider one of the two options:\n";
	print STDERR "!!!   1.) Get your system admin to set the directory to be writeable to you and/or\n";
	print STDERR "!!!       your group.\n";
	print STDERR "!!!          -or-\n";
	print STDERR "!!!   2.) Use the \"-preparsedDir <directory>\" option to specify a directory that\n";
	print STDERR "!!!       you do have permission to write files to.\n";
	print STDERR "\n";
	exit;
}

print STDERR "\tBy default, HOMER will set the the preparsed Directory to be group read/write\n";
`chmod g+rw "$genomeParseDir"`;

my $tmpID = rand();
my $tmpFile = $tmpID . ".tmp";
my $tmpFile2 = $tmpID . ".2.tmp";

my $totalTSS = 0;
my %tss = ();

if ($noRefFileFlag) {

	my $numRefRegions = 10000;

	`homerTools extract stats "$genomeDir" > "$tmpFile"`;
	my %chr = ();
	open IN, $tmpFile;
	my $lineNum=0;
	my $gtotal = 0;
	while (<IN>) {
		$lineNum++;
		next if ($lineNum < 2);
		chomp;
		my @line = split /\t/;
		if ($line[0] eq 'genome') {
			$gtotal = $line[1];
		} else {
			$chr{$line[0]}=$line[1];
		}
	}
	close IN;
	foreach(keys %chr) {
		my $chr = $_;
		my $chrSize = $chr{$chr};
		my %c = ();
		$tss{$chr} = \%c;
		
		my $numToGet = floor($numRefRegions*$chrSize/$gtotal);
		for (my $i=0;$i<$numToGet;$i++) {
			my $pos = floor(rand()*($chrSize));
			my $dir = 0;
			if (rand() > 0.5) {
				$dir = 1;
			}
			my $bestID = "$chr-$pos-$dir";
			$tss{$chr}->{$bestID} = {p=>$pos,d=>$dir};
			$totalTSS++;
		}
	}

} else {

	open IN, $refFile or die "Could not open reference file $refFile\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $line[0];
		my $chr = $line[1];
		my $dir = $line[4];
		if ($dir eq '-' || $dir eq '1') {
			$dir = 1;
		} else {
			$dir = 0;
		}
		my $pos = floor(($line[2]+$line[3])/2);
		my $bestID = $id . "-$chr-$pos-$dir";
		if (!exists($tss{$chr})) {
			my %c = ();
			$tss{$chr} = \%c;
		}
		$tss{$chr}->{$bestID} = {p=>$pos,d=>$dir};
		$totalTSS++;
	}
	close IN;

}


my $posFile = $genomeParseDir . "$outputGenome.$regionSize.pos";

open POS, ">$posFile";

my $totalPossible = $totalTSS*($window/$regionSize);
my $successRate = $maxRegions/$totalPossible;

foreach(keys %tss) {
	my $chr= $_;

	my @tss = sort {$tss{$chr}->{$a}->{'p'} <=> $tss{$chr}->{$b}->{'p'}} keys %{$tss{$chr}};
	
	for (my $i=0;$i<@tss;$i++) {
		my $gid = $tss[$i];
		my $tss = $tss{$chr}->{$tss[$i]}->{'p'};
		my $dir = $tss{$chr}->{$tss[$i]}->{'d'};

		my $start = $tss-$window;
		if ($start < 1) {
			$start = 1;
		}
		if ($i>0) {
			my $ptss = $tss{$chr}->{$tss[$i-1]}->{'p'};
			my $midpoint = floor(($ptss+$tss)/2);
			if ($start < $midpoint) {
				$start = $midpoint;
			}
		}

		my $end = $tss+$window;
		if ($i<@tss-1) {
			my $ntss = $tss{$chr}->{$tss[$i+1]}->{'p'};
			my $midpoint = floor(($ntss+$tss)/2);
			if ($end > $midpoint) {
				$end = $midpoint;
			}
		}
		for (my $p=$tss;$p>$start;$p-=$regionSize) {
			my $diff= $p-$tss;
			if ($dir == 1) {
				$diff *= -1;
			}
			if (abs($diff) > $minInclude) {
				next if (rand() > $successRate);
			}
			my $id = $gid . '+' . $diff;
			if ($diff < 0) {
				$id = $gid . $diff;
			}
			my $s = $p-$regionSizeHalf;
			my $e = $p+$regionSizeHalf;
			my $len = $e-$s;
			#print TSS "$id\t$diff\n";
			print POS "$id\t$chr\t$s\t$e\t$dir\n";
		}
		for (my $p=$tss+$regionSize;$p<$end;$p+=$regionSize) {
			my $diff= $p-$tss;
			if ($dir == 1) {
				$diff *= -1;
			}
			if (abs($diff) > $minInclude) {
				next if (rand() > $successRate);
			}
			my $id = $gid . '+' . $diff;
			if ($diff < 0) {
				$id = $gid . $diff;
			}
			my $s = $p-$regionSizeHalf;
			my $e = $p+$regionSizeHalf;
			my $len = $e-$s;
			#print TSS "$id\t$diff\n";
			print POS "$id\t$chr\t$s\t$e\t$dir\n";
		}
	}
}

close POS;


my $seqFile = $genomeParseDir . "$outputGenome.$regionSize.seq";
my $freqFile = $genomeParseDir . "$outputGenome.$regionSize.cgfreq";
my $cgBinsFile = $genomeParseDir . "$outputGenome.$regionSize.cgbins";
my $gcBinsFile = $genomeParseDir . "$outputGenome.$regionSize.gcbins";
#my $tssBinsFile = "$genome.$regionSize.tssbins";
#my $cgtssBinsFile = "$genome.$regionSize.cgtssbins";
#my $gctssBinsFile = "$genome.$regionSize.gctssbins";


print STDERR "\tExtracting sequences\n";
my $mflag = '';
if ($maskFlag) {
	$mflag = " -mask ";
}
`homerTools extract "$posFile" "$genomeDir" $mflag > "$tmpFile"`;
`cleanUpSequences.pl "$tmpFile" > "$tmpFile2"`;
`removePoorSeq.pl "$tmpFile2" > "$seqFile"`;
`homerTools freq -gc "$freqFile" -format tsv "$seqFile" > "$tmpFile2"`;
`freq2group.pl "$freqFile" > "$cgBinsFile"`;
`freq2group.pl "$freqFile" 2 > "$gcBinsFile"`;
#`getTSSBins.pl $tssFile > $tssBinsFile`;
#`combineBinFiles.pl $cgBinsFile $tssBinsFile > $cgtssBinsFile`;
#`combineBinFiles.pl $gcBinsFile $tssBinsFile > $gctssBinsFile`;
`rm "$tmpFile" "$tmpFile2"`;
print STDERR "\n";
