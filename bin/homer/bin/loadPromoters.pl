#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";

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

my $defaultdirectory = $homeDir . "/data/promoters/";
my $defaultSize = 4000;
my $fastaFlag = 0;
my $defaultStart = -300;
my $defaultEnd = 50;
my $redunDist = 500;
my $config = HomerConfig::loadConfigFile();
if (@ARGV < 1) {
	printCMD();
}
sub printCMD {

	print STDERR "\n\tProgram will prepare a custom promoter set for motif analysis\n\n";
	print STDERR "\tUsage: loadPromoters.pl <required parameters ...> [options]\n";
	print STDERR "\n\tProgram with generate promoter analysis files and update homer configuration\n";
	print STDERR "\tBy default this will create files *.base and *.base.gene using all IDs.  These files are\n";
	print STDERR "\tused as default background sets for motif finding and GO analysis respectively\n";
	print STDERR "\n\tGeneral Parameters:\n";
	print STDERR "\t\tRequired parameters:\n";
	print STDERR "\t\t-name <promoter set name> (used to refer to it later with findMotifs.pl)\n";
	print STDERR "\t\t-org <organism> (Name of organism, ok to set to 'null' if not in HOMER)\n";
	print STDERR "\t\t-id <id type> (specify one: gene, refseq, unigene, ensembl or custom)\n";
	print STDERR "\n\t\tOptional parameters:\n";
	print STDERR "\t\t-force (Overwrite existing promoter set definition)\n";
	print STDERR "\t\t-keepAccVer (By default, promoter IDs with an accession version number will be\n";
	print STDERR "\t\t\ttrimmed off - i.e. NM_012345.2 -> NM_012345 - use this flag keep the .#)\n";
	print STDERR "\t\t-version  <version id> (Assign version, versions starting with 'v' are managed\n";
    print STDERR "\t\t\tby the configureHomer.pl script and my be overwritten, default: custom)\n";
	#print STDERR "\t\t-d <output directory> (default: $defaultdirectory)\n";
	print STDERR "\t\t-as <start> (Redundant/CpG analysis start, default: $defaultStart)\n";
	print STDERR "\t\t-ae <end> (Redundant/CpG analysis end, default: $defaultEnd)\n";
	print STDERR "\n\tPromoters may be specified by either FASTA or genomic position:\n";
	print STDERR "\n\tOption 1: Genome+TSS positions:\n";
	print STDERR "\t\tRequired parameters:\n";
	print STDERR "\t\t-genome <genome> (homer genome version -OR- genome FASTA file)\n";
	print STDERR "\t\t-tss <TSS peak file> (peak file centered on TSS positions)\n";
	print STDERR "\n\t\tOptional parameters:\n";
	my $half = $defaultSize/2;
	print STDERR "\t\t-size <size> (default: $defaultSize, i.e. +/- $half from the TSS)\n";
	print STDERR "\t\t-dist <#> (Distance between promoters to consider redundant, default: 500)\n";
	print STDERR "\n\tOption 2: Promoter FASTA file:\n";
	print STDERR "\t\tRequired parameters:\n";
	print STDERR "\t\t-fasta <promoter FASTA file> (FASTA file of promoter regions)\n";
	print STDERR "\t\t\tEach promoter should be the same length with only the ID after each \">\"\n";
	print STDERR "\t\t\tFASTA Files will be considered \"masked\" so that it will be the default option\n";
	print STDERR "\t\t-offset <#> (offset of the first base, i.e. -1000 for 1kb upstream)\n";


	print STDERR "\n\tAvailable Genomes:\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}\n";
	}
	print STDERR "\n";
	exit;
}

my $astart = $defaultStart;
my $aend = $defaultEnd;
my $fastafile = '';
my $promoterName = "";
my $tssPeakFile = "";
my $organism = "";
my $genome = "";
my $genomeDir = "";
my $placeholder = "";
my $version = 'custom';
my $idtype = "";
my $fastaOffset = "";
my $forceFlag = 0;
my $keepAccVer = 0;

my $directory = $defaultdirectory;
my $size = $defaultSize;

for (my $i=0;$i<@ARGV;$i++){ 
	if ($ARGV[$i] eq '-d') {
		$directory = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-name') {
		$promoterName = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-org') {
		$organism = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-id') {
		$idtype = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-version') {
		$version = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-tss') {
		$tssPeakFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-fasta') {
		$fastafile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-offset') {
		$fastaOffset = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-dist') {
		$redunDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-size') {
		$size = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-keepAccVer') {
		$keepAccVer = 1;
	} elsif ($ARGV[$i] eq '-force') {
		$forceFlag = 1;
	} elsif ($ARGV[$i] eq '-as') {
		$astart = $ARGV[++$i];
		$redunDist = -1;
	} elsif ($ARGV[$i] eq '-ae') {
		$aend = $ARGV[++$i];
		$redunDist = -1;
	} else {
		print STDERR "!!! $ARGV[$i] not recognized\n";
		printCMD();
	}
}


if ($promoterName eq '' || $organism eq '' || $idtype eq '') {
	print STDERR "!!! -name, -org, and -id parameters are all required !!!\n";
	printCMD();
}
if ($genome ne '' && $tssPeakFile ne '') {
	print STDERR "\tWill configure promoters from genomic positions ($genome, $tssPeakFile)\n";
	$fastaFlag = 0;
	if (!exists($config->{'GENOMES'}->{$genome}) && $genome ne 'null') {
		($genome,$genomeDir,$placeholder) = HomerConfig::parseCustomGenome($genome);
	} else {
	    $genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'} . "/";
	}
} elsif ($fastafile ne '' && $fastaOffset ne '') {
	print STDERR "\tWill configure promoters from FASTA file ($fastafile, offset=$fastaOffset)\n";
	if ($genome eq '') {
		$genome = 'null';
	}
	$fastaFlag = 1;
} else {
	print STDERR "!!! Either -genome, -tss OR -fasta, -offset are required !!!\n";
	printCMD();
}

if (!exists($config->{'ORGANISMS'}->{$organism})) {
	print STDERR "!!! Warning: Organism \"$organism\" not found in config file !!!\n";
}

if (exists($config->{'PROMOTERS'}->{$promoterName})) {
	print STDERR "!!! Warning: Promoter set $promoterName exists !!!\n";
	if ($forceFlag == 0) {
		print STDERR "!!! Add -force to the command to overwrite !!!\n";
		exit;
	}
}
print STDERR "\n\tWating 10 seconds in case you want to review the changes (hit ctrl+C to cancel)\n\t\t";
for (my $i=10;$i>0;$i--) {
	print STDERR " $i";
	`sleep 1`;
}
print STDERR "\n";



my $posFile = $directory . "/$promoterName.pos";
my $seqFile = $directory . "/$promoterName.seq";
my $maskFile = $directory . "/$promoterName.mask";
my $consFile = $directory . "/$promoterName.cons";
my $redunFile = $directory . "/$promoterName.redun";
my $cgfreqFile = $directory . "/$promoterName.cgfreq";
my $cgbinsFile = $directory . "/$promoterName.cgbins";
my $gcbinsFile = $directory . "/$promoterName.gcbins";
my $baseFile = $directory . "/$promoterName.base";
my $geneBaseFile = $directory . "/$promoterName.base.gene";
my $rand = rand();
my $tmpSeqFile = $rand . ".seq.tmp";
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $consDone = 0;


my $analysisSource = "";
my $offset = 0;

if ($fastaFlag == 0) {
	$offset = -1*floor($size/2);
	if (abs($astart) > abs($offset)) {
		print STDERR "-as must be set within +/$offset (currently = $astart)\n";
		exit;
	}
	if (abs($aend) > abs($offset)) {
		print STDERR "-ae must be set within +/$offset (currently = $aend)\n";
		exit;
	}

	print STDERR "\tCreating TSS position file...\n";
	`resizePosFile.pl "$tssPeakFile" $size > "$posFile"`;
	if ($keepAccVer == 0) {
		`removeAccVersion.pl "$posFile" > "$tmpFile"`;
		`mv "$tmpFile" "$posFile"`;
	}
	
	print STDERR "\tExtracting Unmasked Promoter Sequences...\n";
	`homerTools extract "$posFile" "$genomeDir" > "$tmpSeqFile"`;
	`cleanUpSequences.pl "$tmpSeqFile" > "$seqFile"`;
	
	$analysisSource = $seqFile;
	if (!$consDone && open IN, $genomeDir . "/conservation/chr1.fa") {
		close IN;
		`homerTools extract "$posFile" "$genomeDir/conservation/" > "$consFile"`;
		$consDone = 1;
	}

	print STDERR "\tExtracting Masked Promoter Sequences...\n";
	`homerTools extract "$posFile" "$genomeDir" -mask > "$tmpSeqFile"`;
	`cleanUpSequences.pl "$tmpSeqFile" > "$maskFile"`;
	$analysisSource = $maskFile;
	if (!$consDone && open IN, $genomeDir . "/conservation/chr1.fa") {
		close IN;
		`homerTools extract "$posFile" "$genomeDir/conservation/" > "$consFile"`;
		$consDone = 1;
	}

} else {
	`fasta2tab.pl "$fastafile" > "$tmpSeqFile"`;
	`cleanUpSequences.pl "$tmpSeqFile" > "$maskFile"`;
	if ($keepAccVer == 0) {
		`removeAccVersion.pl "$maskFile" > "$tmpFile"`;
		`mv "$tmpFile" "$maskFile"`;
	}
	$offset = $fastaOffset;
	$analysisSource = $maskFile;
}

`getPartOfPromoter.pl "$analysisSource" $astart $aend $offset > "$tmpSeqFile"`;
print STDERR "\tFinding Redundant Promoters\n";

if ($fastaFlag == 0 && $redunDist >= 0) {
	`mergePeaks "$posFile" -d $redunDist > "$tmpFile"`;
	`annotateRelativePosition.pl "$posFile," "$tmpFile," > "$tmpFile2"`;
	open IN, $tmpFile2;
	my %meta = ();
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		if (!exists($meta{$line[1]})) {
			my @a = ();
			$meta{$line[1]} = \@a;
		}
		push(@{$meta{$line[1]}},$line[0]);
	}
	close IN;
	open OUT, ">$redunFile";
	foreach(keys %meta) {
		my @promoters = @{$meta{$_}};
		if (@promoters > 1) {
			for (my $i=0;$i<@promoters;$i++) {
				print OUT "$promoters[$i]\t";
				my $c = 0;
				for (my $j=0;$j<@promoters;$j++) {
					next if ($i==$j);
					print OUT "," if ($c > 0);
					print OUT "$promoters[$j]";
				}
				print OUT "\n";
			}
		}
	}
	`rm "$tmpFile" "$tmpFile2"`;
} elsif ($fastaFlag == 1) {
	`findRedundantBLAT.pl "$tmpSeqFile" > "$redunFile"`;
}

print STDERR "\tCalculating CpG Content\n";
`homerTools freq "$tmpSeqFile" -gc "$cgfreqFile" >  /dev/null`;
print STDERR "\tCalculating CpG Bins\n";
`freq2group.pl "$cgfreqFile" > "$cgbinsFile"`;
`freq2group.pl "$cgfreqFile" 2 > "$gcbinsFile"`;

print STDERR "\tCreating base files (default background)\n";
`cut -f1 "$analysisSource" > "$baseFile"`;
`convertIDs.pl "$baseFile" $organism gene > "$geneBaseFile"`;
`rm "$tmpSeqFile"`;


print STDERR "\n\tUpdating Configuration File...\n";
$genome = 'null' if ($genome eq '');
my @params = ($organism,$genome,$idtype,$offset,$offset*-1);
my $g = {org=>$organism,version=>$version,location=>"data/promoters/",
			directory=>$directory,url=>"LocalUpdate",
			desc=>"$promoterName promoters ($organism)",
            parameters=>\@params};

$config->{'PROMOTERS'}->{$promoterName} = $g;
HomerConfig::printConfigFile($config);
print STDERR "\n";
