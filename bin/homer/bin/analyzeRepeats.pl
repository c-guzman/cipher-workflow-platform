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
use Statistics;

my $config = HomerConfig::loadConfigFile();

sub printCMD {
	print STDERR "\n\tUsage: analyzeRepeats.pl <input file> <genome> [options] ...\n";
	print STDERR "\n\tWhere <input file> can be one of:\n";
	print STDERR "\t\trna (quantify gene expression of RefSeq annotated genes)\n";
	print STDERR "\t\trepeats (just type \"repeats\" to load the default repeats for genome)\n";
	print STDERR "\t\t<gtf file> (quantify gene expression from GTF file)\n";
	#print STDERR "\t\t<rmsk file> (UCSC formated repeatmasker file from annotation database)\n";
	#print STDERR "\t\t<rmsk directory> (Directory containing multiple rmsk files)\n";
	print STDERR "\t\t<transcript file> (Homer formatted peak/transcript file i.e. hg19.repeats)\n";
	print STDERR "\n\tAvailable Genomes (required argument): (name,org,directory,default promoter set)\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}"
					. "\t$config->{'GENOMES'}->{$_}->{'promoters'}\n";
	}
	print STDERR "\t\tGenome can also specify a FASTA file, or just specify \"none\"\n";

	print STDERR "\n\tSpecifying distinct classes of repeats or filtering based on parameters:\n";
	print STDERR "\t\t-L1 <repeat name> (level one 'repeat name', i.e. AluSx3)\n";
	print STDERR "\t\t-L2 <repeat name> (level two 'repeat family', i.e. SINE)\n";
	print STDERR "\t\t-L3 <repeat name> (level three 'repeat class', i.e. Alu)\n";
	print STDERR "\t\t-maxdiv (max divergence, i.e. -div 0.10, default: 1.0)\n";
	print STDERR "\t\t-mindiv (min divergence, default: 0)\n";
	print STDERR "\t\t-minLength <#> (only return repeats at lest this length, default: 0)\n";
	print STDERR "\t\t-maxLength <#> (only return repeats less than % of full length, default: no max)\n";
	print STDERR "\t\t-minLengthP <#> (only return repeats at lest % of full length, default: 0%)\n";
	print STDERR "\t\t-maxLengthP <#> (only return repeats less than % of full length, default: 100%)\n";
	print STDERR "\t\t-noexon (do not consider repeats found within coding exons)\n";
	print STDERR "\t\t-condenseL2, -condenseL3 (combine read counts for repeats for same L2 or L3 annotation)\n";
	print STDERR "\t\t-rmsk (force filtering of repeats in case it thinks your file is for RNA)\n";

	print STDERR "\n\tAdjusting coordinates returned:\n";
	print STDERR "\t\t-5p (return peak files centered on 5' end of repeats)\n";
	print STDERR "\t\t-3p (return peak files centered on 3' end of repeats)\n";
	print STDERR "\t\t-og (return positions relative to full length repeats)\n";

	print STDERR "\n\tGTF file options if specifying a GTF file:\n";
	print STDERR "\t\t-gff/-gff3 (for GFF or GFF3 formated files - ideally use a GTF formated file, default)\n";
	print STDERR "\t\t-gid (use gene_id instead of transcript_id when parsing GTF file)\n";

	print STDERR "\n\tExpression/Read Coverage Reporting Options:\n";
	print STDERR "\t\t-d <tag directory 1> [tag directory 2] ... (list of experiment directories to show\n";
	print STDERR "\t\t\ttag counts for) NOTE: -dfile <file> where file is a list of directories in first column\n";
	print STDERR "\t\t-count <genes|exons|introns|5utr|3utr|cds|pausing> (regions to count reads in, default: genes)\n";
	print STDERR "\t\t\tSpecify \"-pausingOptions\" to get more options for pausing analysis for GRO-Seq data\n";
	
	print STDERR "\t\t-strand <+|-|both> (count tags on indicated strand, default: both)\n";
	print STDERR "\t\t-pc <#> or -tbp <#> (maximum tags to count per position, default: 0=no limit)\n";
	print STDERR "\t\t-log (output tag counts as randomized log2 values - for scatter plots)\n";
	print STDERR "\t\t-sqrt (output tag counts as randomized sqrt values - for scatter plots)\n";
	print STDERR "\t\t-condenseGenes (Only report one transcript from each gene locus with highest expression)\n";
	print STDERR "\t\t-noCondensing (do not condense counts from entries will same ID, default: do condense)\n";
	print STDERR "\t\t-noCondensingParts (i.e. report exons separately)\n";
	print STDERR "\t\t-min <#> (minimum expression value to print, default: none)\n";
	print STDERR "\t\tNormalization:\n";
	print STDERR "\t\t\t-rpkm/-fpkm (Report results as fragments per kb per million mapped)\n";
	print STDERR "\t\t\t-norm <#> (Normalize to total mapped tags: default 1e7)\n";
	print STDERR "\t\t\t-normMatrix <#> (Normalize to total tags in gene expression matrix: not used)\n";
	print STDERR "\t\t\t-noadj (Don't normalize, -raw works too)\n";
	print STDERR "\t\tAdvanced Normalization: (will adjust normalization from above)\n";
	print STDERR "\t\t\t-quantile (quantile normalization on reported genes)\n";
	print STDERR "\n\tChecking read-through expression: (example: \"-upstream -1000 -L 2\")\n";
	print STDERR "\t\t-upstream <#> (Distance upstream of each repeat to check for reads, default: 0 [off])\n";
	print STDERR "\t\t-downstream <#> (Distance downstream to each repeat to check for reads, default: 0 [off])\n";
	print STDERR "\t\t-L <#> (Only keep repeats with local enrichment greater than #, default: keep all)\n";

	print STDERR "\n\tAnnotation:\n";
	print STDERR "\t\t-gwasCatalog <gwasCatalog.txt> (Find overlapping GWAS Catalog file from UCSC)\n";

	print STDERR "\n\tOther options:\n";
	#print STDERR "\t\t-cpu <#> (number of CPUs - spawns separate tag counting processes)\n";

	#print STDERR "\n\tHistogram options:\n";
	#print STDERR "\t\t-hist <bin size in bp> (i.e 1, 2, 5, 10, 20, 50, 100 etc.)\n";

	print STDERR "\n";
	exit;
}

sub printPausingCMD {
	print STDERR "\n\tPausing Options (use with \"-count pausing\"):\n";
	print STDERR "\t\t-pPromoterStart <#> (start of promoter region relative to TSS, default: -50)\n";
	print STDERR "\t\t-pPromoterEnd <#> (end of promoter region relative to TSS, default: 200)\n";
	print STDERR "\t\t-pBodyStart <#> (start of gene body region relative to TSS, default: 200)\n";
	print STDERR "\t\t-pBodyEnd <#> (end of gene body region relative to TSS, default: 5000)\n";
	print STDERR "\t\t-pPsuedoRatio <#> (Reads per bp [pseudo count] for ratio calculation, default: 0.005)\n";
	exit;
}

if (@ARGV < 2) { 
	printCMD();
}

my $maskFlag = 0;
my $upstream = 0;
my $downstream = 0;
my $localFold = 0;
my $condenseFlag = 1;
$condenseGenesFlag = 0;
my $condensePartsFlag = 1;
my $partsFlag = 0;
$gwasCatalog = "";
my $rmskFlag = 0;
my $gtfFlag = '';
my $gidFlag = '';
$quantileNorm = 0;

my $logFlag = 0;
my $sqrtFlag = 0;
my $numCPUs = 1;
$pausePromoterStart = -50;
$pausePromoterEnd = 200;
$pauseGeneBodyStart = 200;
$pauseGeneBodyEnd = 5000;
$pausePseudoRatio = 0.005; #one read per 200 bp

my $fragLength = 150;
my $tbp = 0;

my $countStrand = "both";

my $histBinSize = 0;
my $size = 4000;

my $countMethod = 'genes';

my $level1 = "";
my $level2 = "";
my $level3 = "";

my $maxDiv = 1e20;
my $minDiv = -1e20;
$minLengthPercent = -1e10;
$maxLengthPercent = 1e10;
$minLength = -1e20;
$maxLength = 1e20;

my $p5Flag = 0;
my $p3Flag = 0;
my $ogFlag = 0;
my $defHalfSize= 100;
my $chosenFormat="";

my @tagDirs = ();

my $rpkmFlag = 0;
my $rpkmLength= 1000;
my $normTotal = 1e7;
my $normMatrix = 0;
my $noNormFlag = 0;

my $minExpValue = -1e20;

my %possibleCounts = ();
$possibleCounts{'genes'} = 0;
$possibleCounts{'exons'} = 1;
$possibleCounts{'introns'} = 1;
$possibleCounts{'3utr'} = 1;
$possibleCounts{'5utr'} = 1;
$possibleCounts{'cds'} = 1;
$possibleCounts{'pausing'} = 1;

my $input = $ARGV[0];
my $genome = $ARGV[1];
my $cmd = 'analyzeRepeats.pl';
foreach(@ARGV) {
	$cmd .= " $_";
}

for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-L1') {
		$level1 = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pausingOptions') {
		printPausingCMD();
	} elsif ($ARGV[$i] eq '-gid') {
		$gidFlag = " -gid";
	} elsif ($ARGV[$i] eq '-gff') {
		$gtfFlag = " -gff";
	} elsif ($ARGV[$i] eq '-gff3') {
		$gtfFlag = " -gff3";
	} elsif ($ARGV[$i] eq '-L2') {
		$level2 = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-L3') {
		$level3 = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-L3') {
		$level3 = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-5p') {
		$p5Flag = 1;
	} elsif ($ARGV[$i] eq '-3p') {
		$p3Flag = 1;
	} elsif ($ARGV[$i] eq '-og') {
		$ogFlag = 1;
	} elsif ($ARGV[$i] eq '-pPromoterStart') {
		$pausePromoterStart = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pPromoterEnd') {
		$pausePromoterEnd = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pBodyStart') {
		$pauseGeneBodyStart = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pBodyEnd') {
		$pauseGeneBodyEnd = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$numCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rmsk') {
		$rmskFlag = 1;
	} elsif ($ARGV[$i] eq '-gwasCatalog') {
		$gwasCatalog = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minLength') {
		$minLength = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-condenseGenes') {
		$condenseGenesFlag = 1;
	} elsif ($ARGV[$i] eq '-quantile') {
		$quantileNorm = 1;
	} elsif ($ARGV[$i] eq '-maxLength') {
		$maxLength = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minLengthP') {
		$minLengthPercent = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxLengthP') {
		$maxLengthPercent = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxdiv') {
		$maxDiv = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-mindiv') {
		$minDiv = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-format') {
		$chosenFormat = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-count') {
		$countMethod = $ARGV[++$i];
		if (!exists($possibleCounts{$countMethod})) {
			print STDERR "!!! -count $countMethod is not a valid option!\n";
			exit;
		}
		$partsFlag = $possibleCounts{$countMethod};
	} elsif ($ARGV[$i] eq '-upstream') {
		$upstream = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-downstream') {
		$downstream = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-noCondensing') {
		$condenseFlag = 0;
	} elsif ($ARGV[$i] eq '-noCondensingParts') {
		$condensePartsFlag = 0;
	} elsif ($ARGV[$i] eq '-condenseL1') {
		$condenseFlag = 1;
	} elsif ($ARGV[$i] eq '-condenseL2') {
		$condenseFlag = 2;
	} elsif ($ARGV[$i] eq '-condenseL3') {
		$condenseFlag = 3;
	} elsif ($ARGV[$i] eq '-L') {
		$localFold = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rpkm' || $ARGV[$i] eq '-fpkm') {
		$rpkmFlag = 1;
		$normTotal = 1e6;
	} elsif ($ARGV[$i] eq '-norm') {
		$normTotal = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-normMatrix') {
		$normMatrix = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$minExpValue = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-noadj' || $ARGV[$i] eq '-raw') {
		$noNormFlag = 1;
	} elsif ($ARGV[$i] eq '-log' || $ARGV[$i] eq '-log2')  {
		$logFlag = 1;
	} elsif ($ARGV[$i] eq '-sqrt')  {
		$sqrtFlag = 1;
	} elsif ($ARGV[$i] eq '-tbp' || $ARGV[$i] eq '-pc') {
		$tbp = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-strand') {
		$countStrand = $ARGV[++$i];
		if ($countStrand ne '-' && $countStrand ne '+' && $countStrand ne 'both') {
			print STDERR "!!! Invalid option for -strand <+|-|both>: \"$countStrand\"\n";
			exit;
		}
	} elsif ($ARGV[$i] eq '-dfile') {
		open IN, $ARGV[++$i];
		print STDERR "\tAdding Tag Directories:\n";
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			push(@tagDirs, $line[0]);
			print STDERR "\t\t$line[0]\n";
		}
		close IN;
	} elsif ($ARGV[$i] eq '-d') {
		print STDERR "\tTag Directories:\n";
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@tagDirs, $ARGV[$i]);
			print STDERR "\t\t$ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} else {
		print STDERR "!!! Option \"$ARGV[$i]\" not recognized...\n";
		printCMD();
	}
}

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

my $rand = rand();
$tmpFile = $rand . ".tmp";
$tmpFile2 = $rand . ".2.tmp";

my @inputFiles = ();
if (-d $input) {
	`ls -1 "$input/"*rmsk* > "$tmpFile"`;
	open IN, $tmpFile;
	while (<IN>) {
		chomp;
		push(@inputFiles, $_);
	}
	close IN;
	$rmskFlag = 1;
	`rm "$tmpFile"`;
} elsif (-f $input) {
	push(@inputFiles, $input);
} else {
	$input = 'repeats' if ($input eq 'repeat');
	if ($input eq 'repeats') {
		$rmskFlag = 1;
	}
	if ($input eq 'repeats' || $input eq 'rna' || $input eq 'miRNA') {
		if ($genomeDir eq '') {
			print STDERR "!!! Can't use \"repeats\", \"rna\", or \"miRNA\" as input file if you don't use a valid homer genome\n";
			exit;
		}
		$input = $genomeDir . "/" . $genome . "." . $input;
		if (! -f $input) {
			print STDERR "!!! Cannot find file \"$input\" - likely a configuration error\n";
		}
	} else {
		print STDERR "!!! What is input \"$input\"?? Not recognized...\n";
		exit;
	}
	push(@inputFiles, $input);
}

my %repeatsHash = ();
my $repeats = \%repeatsHash;
my $readIntoMemory = 1;
if ($gwasCatalog eq '' && (@tagDirs < 1 && $histBinSize == 0)) {
	$readIntoMemory = 0;
}

%repeatIDs = ();
$totalRepeats = 0;
$totalRepeatsKept = 0;
open TMP, ">$tmpFile";
my $filehandle = *TMP;
foreach(@inputFiles) {
	my $filename = $_;
	my $tmpInputFile = "";
	my $format = checkFormat($filename);
	print STDERR "\tInput file format: $format\n";
	if ($chosenFormat ne '') {
		$format = $chosenFormat;
	}
	if ($format eq 'bed') {
		$tmpInputFile = rand() . ".tmp";
		print STDERR "\tInitial format set to BED\n";
		`bed2pos.pl "$filename" > "$tmpInputFile"`;
		$filename = $tmpInputFile;
		$format = 'homerRmsk';
	}
	if ($format eq 'gtf') {
		$tmpInputFile = rand() . ".tmp";
		`parseGTF.pl "$filename" rna $gtfFlag $gidFlag > "$tmpInputFile"`;
		$filename = $tmpInputFile;
		$format = 'homerRmsk';
	}
	if ($format eq 'rmsk') {
		readRMSKfile($filename, $filehandle, $level1, $level2, $level3, $ogFlag, $p5flag, $p3flag, $defHalfSize, 
							$maxDiv, $minDiv, $repeats, $readIntoMemory);
	} elsif ($format eq 'homerRmsk') {
		readHomerRmsk($filename, $filehandle, $level1, $level2, $level3, $ogFlag, $p5flag, $p3flag, $defHalfSize, 
							$maxDiv, $minDiv, $repeats, $readIntoMemory,$countMethod,$rmskFlag);
	}
	if ($tmpInputFile ne '') {
		`rm "$tmpInputFile"`;
	}
}
close $filehandle;

print STDERR "\tFiltering based on repeat parameters: kept $totalRepeatsKept of $totalRepeats\n";
if ($gwasCatalog ne '') {
	annotateGWAS($repeats,$gwasCatalog);
}

if (@tagDirs < 1 && $histBinSize == 0) {
	print "RepeatID (cmd=$cmd)\tchr\tstart\tend\tstrand\tCode\tDivergence\tFullStart\tFullEnd";
	if ($gwasCatalog ne '') {
		print "\tGWAS Risk SNPs";
	}
	print "\n";
	open IN, $tmpFile;
	while (<IN>) { 
		chomp;
		s/\r//g;
		print $_;
		if ($gwasCatalog ne '') {
			print "\t";
			my @line = split /\t/;
			if (exists($repeats->{$line[0]}->{'gwas'})) {
				print $repeats->{$line[0]}->{'gwas'};
			}
		}
		print "\n";
	}
	close IN;
	`rm "$tmpFile"`;
	exit;
}
my %exp = ();

if (@tagDirs > 0 && $histBinSize == 0) {

	my $tmpUpstream = "";
	my $tmpDownstream = "";

	my $indexSize = 1;
	if ($upstream ne 0) {
		$tmpUpstream = $rand . ".upstream.txt";
		repositionFile($tmpFile,$tmpUpstream,$upstream,0);
		$indexSize++;
	}
	if ($downstream ne 0) {
		$tmpDownstream = $rand . ".downstream.txt";
		repositionFile($tmpFile,$tmpDownstream,0,$downstream);
		$indexSize++;
	}



	my @tagTotals = ();
	my @matrixTotals = ();
	my $curIndex = 0;
	for (my $i=0;$i<@tagDirs;$i++) {
		print STDERR "\tCalculating read coverage for $tagDirs[$i]\n";
		my $directory = $tagDirs[$i];
		my ($totalDirTags, $totalDirPos, $fragEst, $peakEst) = HomerConfig::readTagInfo($directory,$tbp);
		push(@tagTotals, $totalDirTags);

		my $mtotal = 0;
		`getPeakTags "$tmpFile" "$directory" -fixed -count -strand $countStrand -tagAdjust 0 -tbp $tbp > "$tmpFile2"`;
		#print STDERR "`getPeakTags $tmpFile $directory -fixed -count -strand $countStrand -tagAdjust 0 -tbp $tbp > $tmpFile2`;\n";
		
		open IN, $tmpFile2;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			if (!exists($exp{$line[0]})) {
				my @a = ();
				$exp{$line[0]} = \@a;
			}
			$mtotal += $line[1];
			push(@{$exp{$line[0]}}, $line[1]);
		}
		close IN;
		`rm "$tmpFile2"`;	
		my $localIndex = $curIndex;
		$curIndex++;

		if ($upstream ne 0) {
			print STDERR "\t\tChecking upstream read coverage\n";
			`getPeakTags "$tmpUpstream" "$directory" -fixed -count -strand $countStrand -tagAdjust 0 -tbp $tbp > "$tmpFile2"`;
			my $localLength = abs($upstream);
			open IN, $tmpFile2;
			while (<IN>) {
				chomp;
				s/\r//g;
				my @line = split /\t/;
				if (!exists($exp{$line[0]})) {
					print STDERR "!!! Something is wrong - can't find $line[0] during upstream calculation\n";
					next;
				}
				my $len = $repeats->{$line[0]}->{'len'};
				my $v= $exp{$line[0]}->[$localIndex];
				$line[1] = 1 if ($line[1] < 1);
				my $fold = ($v/$len) / ($line[1]/$localLength);
				push(@{$exp{$line[0]}}, $fold);
			}
			close IN;
			`rm "$tmpFile2"`;	
			$curIndex++;
		}
		if ($downstream ne 0) {
			print STDERR "\t\tChecking downstream read coverage\n";
			`getPeakTags "$tmpDownstream" "$directory" -fixed -count -strand $countStrand -tagAdjust 0 -tbp $tbp > "$tmpFile2"`;
			my $localLength = abs($downstream);
			open IN, $tmpFile2;
			while (<IN>) {
				chomp;
				s/\r//g;
				my @line = split /\t/;
				if (!exists($exp{$line[0]})) {
					print STDERR "!!! Something is wrong - can't find $line[0] during downstream calculation\n";
					next;
				}
				my $len = $repeats->{$line[0]}->{'len'};
				my $v= $exp{$line[0]}->[$localIndex];
				$line[1] = 1 if ($line[1] < 1);
				my $fold = ($v/$len) / ($line[1]/$localLength);
				push(@{$exp{$line[0]}}, $fold);
			}
			close IN;
			`rm "$tmpFile2"`;	
			$curIndex++;
		}

		push(@matrixTotals, $mtotal);
	}

	if ($countMethod eq 'pausing') {
		my $numDataPoints = scalar(@tagDirs)*$indexSize;

		my $plen = $pausePromoterEnd-$pausePromoterStart;
		my $glen = $pauseGeneBodyEnd-$pauseGeneBodyStart;
		foreach(keys %$repeats) {
			my $oid = $_;
			my $id = $oid;
			if ($id =~ s/\--Part(\d*)$//) {
			} else {
				print STDERR "Problem with IDs: \'$id\'\n";
			}
			my $partNumber = $1;
			next if ($partNumber ne '1');
			my @pausing =();
			for (my $i=0;$i<$numDataPoints;$i++) {
				my $cPromoter = $exp{$id . "--Part1"}->[$i];
				my $cGeneBody = $exp{$id . "--Part2"}->[$i];
				my $ratio = (($cPromoter/$plen) + $pausePseudoRatio)/(($cGeneBody/$glen)+$pausePseudoRatio);
				push(@pausing, $ratio, $cPromoter,$cGeneBody);
			}
			$repeats->{$id} = $repeats->{$oid};
			$repeats->{$id}->{'s'} = $repeats->{$id}->{'ogs'};
			$repeats->{$id}->{'e'} = $repeats->{$id}->{'oge'};
			$exp{$id} = \@pausing;
			delete $repeats->{$id . "--Part1"};
			delete $repeats->{$id . "--Part2"};
			delete $exp{$id . "--Part1"};
			delete $exp{$id . "--Part2"};
		}
		$indexSize+=2; #this is becuase we added
	} elsif ($partsFlag && $condensePartsFlag) {
		my $numDataPoints = scalar(@tagDirs)*$indexSize;
		foreach(keys %$repeats) {
			my $oid = $_;
			my $id = $oid;
			$id =~ s/\--Part\d*$//;
			if ($oid ne $id) {
				if (exists($repeats->{$id})) {
					$repeats->{$id}->{'nexons'}+=$repeats->{$oid}->{'nexons'};
					$repeats->{$id}->{'len'}+=$repeats->{$oid}->{'len'};
					$repeats->{$id}->{'div'}+=$repeats->{$oid}->{'div'};
					my $preTotal=0;
					my $newTotal = 0;
					for (my $i=0;$i<$numDataPoints;$i++) {
						$preTotal += $exp{$id}->[$i];
						$newTotal += $exp{$oid}->[$i];
						$exp{$id}->[$i] += $exp{$oid}->[$i];
					}
				} else {
					$repeats->{$id} = $repeats->{$oid};
					$repeats->{$id}->{'s'} = $repeats->{$id}->{'ogs'};
					$repeats->{$id}->{'e'} = $repeats->{$id}->{'oge'};
					$exp{$id} = $exp{$oid};
				}
				delete $repeats->{$oid};
				delete $exp{$oid};
			}
		}
		foreach(keys %$repeats) {
		#	$repeats->{$_}->{'len'} /= $repeats->{$_}->{'n'};
			$repeats->{$_}->{'div'} /= $repeats->{$_}->{'nexons'};
		}
	}

	if ($condenseFlag) {

		my $numDataPoints = scalar(@tagDirs)*$indexSize;
		foreach(keys %$repeats) {
			my $oid = $_;
			my $id = $oid;
			$id =~ s/\-HOMER\d*$//;
			if ($condenseFlag > 1) {
				$id =~ s/^.+?\|//;
			}
			if ($condenseFlag == 2) {
				$id =~ s/\|.+?$//;
			}
			if ($oid ne $id) {
				if (exists($repeats->{$id})) {
					$repeats->{$id}->{'n'}+=$repeats->{$oid}->{'n'};
					$repeats->{$id}->{'len'}+=$repeats->{$oid}->{'len'};
					$repeats->{$id}->{'div'}+=$repeats->{$oid}->{'div'};
					my $preTotal=0;
					my $newTotal = 0;
					for (my $i=0;$i<$numDataPoints;$i++) {
						$preTotal += $exp{$id}->[$i];
						$newTotal += $exp{$oid}->[$i];
						$exp{$id}->[$i] += $exp{$oid}->[$i];
					}
					if ($newTotal > $preTotal) {
						$repeats->{$id}->{'s'} = $repeats->{$oid}->{'s'};
						$repeats->{$id}->{'e'} = $repeats->{$oid}->{'e'};
						$repeats->{$id}->{'c'} = $repeats->{$oid}->{'c'};
					}
				} else {
					$repeats->{$id} = $repeats->{$oid};
					$exp{$id} = $exp{$oid};
				}
				delete $repeats->{$oid};
				delete $exp{$oid};
			}
		}
		foreach(keys %$repeats) {
			$repeats->{$_}->{'len'} /= $repeats->{$_}->{'n'};
			if (isNumber($repeats->{$_}->{'div'})) {
				$repeats->{$_}->{'div'} /= $repeats->{$_}->{'n'};
			}

			if ($countMethod eq 'pausing') {
				my $plen = $pausePromoterEnd-$pausePromoterStart;
				my $glen = $pauseGeneBodyEnd-$pauseGeneBodyStart;
				for (my $i=0;$i<$numDataPoints;$i+=3) {
					my $cPromoter = $exp{$_}->[$i+1];
					my $cGeneBody += $exp{$_}->[$i+2];
					my $ratio = (($cPromoter/$plen) + $pausePseudoRatio)/(($cGeneBody/$glen)+$pausePseudoRatio);
					$exp{$_}->[$i] = $ratio;
				}
			}
		}
	}


	my @headers = ();
	#normalization
	for (my $i=0;$i<@tagDirs;$i++) {
		my $scaleFactor = 1;
		print STDERR "\t\tNormalizing...\n" if ($noNormFlag==0);

		if ($noNormFlag) {
		} elsif ($normMatrix > 0) {
			if ($matrixTotals[$i] > 0) {
				$scaleFactor = $normMatrix/$matrixTotals[$i];
			} else {
				print STDERR "!!! Warning: Total of $matrixTotals[$i] reads found for $tagDirs[$i]\n";
			}
		} else {
			if ($tagTotals[$i] > 0) {
				$scaleFactor = $normTotal/$tagTotals[$i];
			} else {
				print STDERR "!!! Warning: Total of $tagTotals[$i] reads found for $tagDirs[$i]\n";
			}
		}

		my $hname = $tagDirs[$i];
		if ($noNormFlag) {
		} elsif ($rpkmFlag) {
			$hname .= " FPKM";
		} elsif ($normMatrix > 0) {
			$hname .= " reads normalized to $normMatrix per column (" . sprintf("%.2lf",$scaleFactor) . ")";
		} else {
			$hname .= " reads normalized to $normTotal mapped reads (" . sprintf("%.2lf",$scaleFactor) . ")";
		}
			
		my $plen = 1;	
		my $glen = 1;
		my $index = $i*$indexSize;
		my @indices = ($index);
		if ($countMethod eq 'pausing') {
			$plen = $pausePromoterEnd-$pausePromoterStart;
			$glen = $pauseGeneBodyEnd-$pauseGeneBodyStart;
			@indices = ($i*$indexSize+1, $i*$indexSize+2);
		}
		for (my $z=0;$z<@indices;$z++) {
			my $index = $indices[$z];
			foreach(keys %$repeats) {
				if (!exists($exp{$_})) {
					print STDERR "missing $_...\n";
					next;
				}
				if ($noNormFlag == 0) {
					$exp{$_}->[$index] *= $scaleFactor;
	
					if ($rpkmFlag) {
						my $len = $repeats->{$_}->{'len'};
						if ($countMethod eq 'pausing') {
							if ($z==0) {
								$len = $plen;
							} else {
								$len = $glen;
							}
						}
						if ($len < 1) {
							print STDERR "Transcript $_ has length less than 1!!!\n";
						} else {
							$exp{$_}->[$index] *= $rpkmLength/$len;
						}
					} else {
						if ($sqrtFlag) {
							$exp{$_}->[$index] += rand()*$scaleFactor;
						}
						if ($logFlag) {
							$exp{$_}->[$index] += 1 + rand()*$scaleFactor;
						}
					}
				}
				if ($sqrtFlag) {
					$exp{$_}->[$index] = sqrt($exp{$_}->[$index]);
					}
				if ($logFlag) {
					$exp{$_}->[$index] = log($exp{$_}->[$index])/log(2.0);
				}
			}
		}
		if ($countMethod eq 'pausing') {
			push(@headers, "$tagDirs[$i] Pausing Ratio");
			push(@headers, "$tagDirs[$i] Promoter Reads ($pausePromoterStart,$pausePromoterEnd)");
			push(@headers, "$tagDirs[$i] GeneBody Reads ($pauseGeneBodyStart,$pauseGeneBodyEnd)");
		} else {
			push(@headers, $hname);
		}
		if ($upstream ne 0) {
			push(@headers, "$tagDirs[$i] Readthrough Enrichment (vs. upstream $upstream)");
		}
		if ($downstream ne 0) {
			push(@headers, "$tagDirs[$i] Readthrough Enrichment (vs. downstream $downstream)");
		}
	}	

	
	if ($tmpUpstream ne "") {
		`rm "$tmpUpstream"`;
	}
	if ($tmpDownstream ne "") {
		`rm "$tmpDownstream"`;
	}

	print STDERR "\tPrinting output\n";
	printExpression($repeats,$organism,\@headers,\%exp,$minExpValue,$localFold,$indexSize);
}

`rm -f "$tmpFile"`;
exit;


sub annotateGWAS {
	my ($repeats, $gwasCatalog) = @_;

	print STDERR "\tChecking for GWAS risk SNP overlap\n";
	open IN, $gwasCatalog;
	my %gwas = ();
	my $idCount = 1;
	my $rand = rand();
	my $tmpFile = $rand . ".tmp";
	my $tmpFile2 = $rand . ".2.tmp";
	my $tmpFile3 = $rand . ".3.tmp";
	open OUT, ">$tmpFile";
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = 'gwas' . $idCount++;
		my $chr = $line[1];
		my $start = $line[2];
		my $end = $line[3];
		my $strand = "+";
		my $snp = $line[4];
		my $study = $line[10];
		$gwas{$id} = {c=>$chr,s=>$start,e=>$end,d=>$strand,snp=>"$snp|$study"};
		print OUT "$id\t$chr\t$start\t$end\t$strand\n";
	}
	close OUT;
	close IN;

	open OUT, ">$tmpFile2";
	foreach(keys %$repeats) {
		print OUT "$_\t$repeats->{$_}->{'c'}\t$repeats->{$_}->{'s'}\t$repeats->{$_}->{'e'}\t$repeats->{$_}->{'d'}\n";
	}
	close OUT;

	`mergePeaks -d given "$tmpFile2" "$tmpFile" > "$tmpFile3" 2> /dev/null`;

	open IN, $tmpFile3;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if (@line < 10);
		my @peaks = split /\,/, $line[8];
		my @snps = split /\,/, $line[9];

		for (my $i=0;$i<@peaks;$i++) {
			my $pid = $peaks[$i];
			next if (!exists($repeats->{$pid}));

			for (my $j=0;$j<@snps;$j++) {
				my $snpID = $snps[$j];
				next if (!exists($gwas{$snpID}));
				my $str = $gwas{$snpID}->{'snp'};
				if (!exists($repeats->{$pid}->{'gwas'})) {
					$repeats->{$pid}->{'gwas'}= $str;
				} else {
					$repeats->{$pid}->{'gwas'} .= "," . $str;
				}
			}
		}
	}
	close IN;

	`rm -f "$tmpFile" "$tmpFile2" "$tmpFile3"`;


}

sub printExpression {
	my ($repeats, $org,$headers, $exp, $minExpValue,$localFold,$indexSize) = @_;

	my %ann = ();
	my %valid = ();
	if ($org ne '' && $org ne 'null') {
		my %map = ();
		my $rand = rand();
		my $tmpFile3 = $rand . ".tmp";
		my $tmpFile4 = $rand . ".2.tmp";
		open OUT, ">$tmpFile3";
		foreach(keys %$repeats) {
			$_ =~ s/-HOMER\d+//;
			print OUT "$_\n";
		}
		if ($org ne 'unknown') {
			`addGeneAnnotation.pl "$tmpFile3" $org no > "$tmpFile4"`;
			open IN, $tmpFile4;
			while (<IN>) {
				chomp;
				s/\r//g;
				my @line = split /\t/;
				if (@line < 10 || $line[5] eq '') {
					$ann{$line[1]} = '';
				} else {
					$ann{$line[1]} = "$line[5]|$line[6]|$line[7]|$line[8]|$line[10]";
					if (!exists($map{$line[0]})) {
						my @a = ();
						$map{$line[0]} = \@a;
					}
					push(@{$map{$line[0]}},$line[1]);
				}
			}
			close IN;
		} else {
			open IN, "$tmpFile3";
			while (<IN>) {
				chomp;
				s/\r//g;
				my @line = split /\t/;
				my @a = ($line[0]);
				$map{$line[0]} = \@a;
			}
			close IN;

		}
		`rm "$tmpFile3" "$tmpFile4"`;

		if ($condenseGenesFlag) {
			foreach(keys %map) {
				my $gid=$_;
				
				my $bestTID = '';
				my $bestExp = -1e20;
				foreach(@{$map{$gid}}) {
					my $tid = $_;
					if (!exists($exp->{$tid})) {
						next;
					}
					my $sum = 0;
					for (my $i=0;$i<@{$exp->{$tid}};$i+=$indexSize) {
						$sum += $exp->{$tid}->[$i];
					}
					if ($sum > $bestExp) {
						$bestExp = $sum;
						$bestTID = $tid;
					}
				}
				if ($bestTID ne '') {
					$valid{$bestTID}=1;
				}
			}
		}
	}

	if ($quantileNorm) {
		print STDERR "\tPerforming Quantile Normalization\n";
		my @matrix = ();
		my @ids = ();
		foreach (keys %$repeats) {
			my $id = $_;
			next if (!exists($exp->{$id}));
			if ($condenseGenesFlag) {
				next if (!exists($valid{$id}));
			}
			my @a = ();
			push(@ids,$id);
			for (my $i=0;$i<@{$exp->{$id}};$i+=$indexSize) {
				push(@a, $exp->{$id}->[$i]);
			}
			push(@matrix,\@a);
		}
		my $qmatrix = Statistics::quantileNorm(\@matrix);
		my $j = -1;
		foreach (@ids) {
			my $id = $_;
			$j++;
			for (my $i=0;$i<@{$exp->{$id}};$i+=$indexSize) {
				my $ii = $i/$indexSize;
				$exp->{$id}->[$i] = $qmatrix->[$j][$ii];
			}
		}
	}


	print "Transcript/RepeatID (cmd=$cmd)\tchr\tstart\tend\tstrand\tLength\tCopies\tAnnotation/Divergence";
	if ($gwasCatalog ne '') {
		print "\tGWAS Risk SNPs";
	}
	foreach(@$headers) {
		print "\t$_";
	}
	print "\n";
	my $totalPossible = 0;
	my $passedFilter = 0;
	foreach (keys %$repeats) {
		$totalPossible++;
		my $id = $_;
		if (!exists($exp->{$id})) {
			next;
		}
		if ($condenseGenesFlag) {
			next if (!exists($valid{$id}));
		}
		my $good = 0;
		my $localGood = 0;
		$localGood = 1 if ($localFold < 0.000001);
		my $expStr = "";
		for (my $i=0;$i<@{$exp->{$id}};$i+=$indexSize) {
			my $v = $exp->{$id}->[$i];
			$expStr .= "\t" . sprintf("%.3lf",$v);
			$good = 1 if ($v >= $minExpValue);
			for (my $j=1;$j<$indexSize;$j++) {
				my $f = $exp->{$id}->[$i+$j];
				$localGood = 1 if ($f >= $localFold);
				$expStr .= "\t" . sprintf("%.3lf",$f);
			}
		}
		next if ($good ==0 || $localGood == 0);
		$passedFilter++;
		my $len = sprintf("%.1lf",$repeats->{$id}->{'len'});
		my $ogDiv = $repeats->{$id}->{'div'};
		my $div = $ogDiv;
		if (isNumber($ogDiv)) {
			$div = sprintf("%.3lf",$repeats->{$id}->{'div'});
		}
		if (exists($ann{$id}) && $ann{$id} ne '') {
			$div = $ann{$id};
		}
		my $gwasStr = "";
		if ($gwasCatalog ne '') {
			$gwasStr = "\t";
			if (exists($repeats->{$id}->{'gwas'})) {
				$gwasStr .= $repeats->{$id}->{'gwas'};
			}
		}
		print "$id\t$repeats->{$id}->{'c'}\t$repeats->{$id}->{'s'}\t$repeats->{$id}->{'e'}\t$repeats->{$id}->{'d'}"
				. "\t$len\t$repeats->{$id}->{'n'}\t$div" . $gwasStr . $expStr . "\n";

	}
	print STDERR "\t\tPrinted $passedFilter of $totalPossible repeats (expression >= $minExpValue)\n";
}
sub isNumber {

	my ($str) = @_;
	if ($str =~ /^[-.\de\+]+$/) {
		return 1;
	} else {
		return 0;
	}
}
	
sub repositionFile {
	my ($ogFile,$newFile,$upstream,$downstream) = @_;
	open OUT, ">$newFile";
	open IN, $ogFile;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $start = 0;
		my $end = 0;
		if ($upstream != 0) {
			$start = $line[2];
			$end = $start + $upstream;
			if ($line[4] eq '-') {
				$start = $line[3];
				$end = $start - $upstream;
			}
		} elsif ($downstream != 0) {
			$start = $line[3];
			$end = $start + $downstream;
			if ($line[4] eq '-') {
				$start = $line[2];
				$end = $start - $downstream;
			}
		}
		if ($start > $end) {
			my $t = $end;
			$end = $start;
			$start = $t;
		}
		print OUT "$line[0]\t$line[1]\t$start\t$end\t$line[4]\n";
	}
	close IN;
	close OUT;
}
sub checkFormat {
	my ($filename) = @_;
	open IN, $filename;
	my $count = 0;
	my $homerCount=0;
	my $rmskCount=0;
	my $gtfCount=0;
	while (<IN>) {
		$count++;
		chomp;
		s/\r//g;
		next if (/^\s*\#/);
		my @line = split /\t/;
		if (@line > 4) {
			if ($line[1] =~ /^chr/) {
				$homerCount++;
			}
		}
		if (@line > 5) {
			if ($line[5] =~ /^chr/) {
				$rmskCount++;
			}
		}
		if (@line > 8) {
			if ($line[8] =~ /gene_id/) {
				$gtfCount++;
				last;
			}
		}
		last if ($count > 20000);
	}
	if ($homerCount  >= $rmskCount && $homerCount >= $gtfCount) {
		return 'homerRmsk';
	} elsif ($rmskCount  > $homerCount && $rmskCount > $gtfCount) {
		return 'rmsk';
	} else {
		return 'gtf';
	}
}
sub readHomerRmsk {
	my ($file, $tmpFileHandle, $level1, $level2, $level3, $ogFlag, $p5flag, $p3flag, $defHalfSize, 
							$maxDiv, $minDiv, $repeats,$readIntoMemory,$countMethod,$rmskFlag) = @_;
	open IN, $file;

	if ($rmskFlag) {
		print STDERR "\t\tWill filter repeats based on parameters...\n";
	}

	my $numL1 = 0;
	my $numL2 = 0;
	my $numL3 = 0;
	my $div = 0;

	while (<IN>) {
		chomp;
		next if (/^\s*#/);
		my @line = split /\t/;
		next if (@line < 5);
		next unless ($line[2] =~ /^[\d\.\-\e\+]+$/);
		next unless ($line[3] =~ /^[\d\.\-\e\+]+$/);
		my @rname = split /\|/, $line[0];
		$totalRepeats++;
		if (@rname > 2) {
			$rname[2]=~s/\-HOMER.*$//;
		}
		if (@rname > 2 && $level1 ne '' && $level1 ne $rname[0]) {
			$numL1++;
			next;
		}
		if (@rname > 2 && $level2 ne '' && $level2 ne $rname[1]) {
			$numL2++;
			next;
		}
		if (@rname > 2 && $level3 ne '' && $level3 ne $rname[2]) {
			$numL3++;
			next;
		}
		my $div = 0;
		if (@line > 6) {
			$div = $line[6];
			if ($rmskFlag && $div =~ /^[\d\.\-\e\+]+$/) {
				next if ($div > $maxDiv);
				next if ($div < $minDiv);
			}
		}
		my $chr = $line[1];
		my $start = $line[2];
		my $end = $line[3];
		my $strand = $line[4];
		if ($strand eq '1') {
			$strand = '-';
		} elsif ($strand eq '0') {
			$strand = '+';
		}
	
		my $ogStart = $start;	
		my $ogEnd = $end;
		my $ogLength  = $ogEnd-$ogStart;
		if (@line > 8 && $rmskFlag) {
			$ogStart = $line[7];
			$ogEnd = $line[8];
			if ($ogStart =~ /^[\d\.\-\e\+]+$/ && $ogEnd =~ /^[\d\.\-\e\+]+$/) {
				$ogLength  = $ogEnd-$ogStart;
			}
		}
		$ogLength = 1 if ($ogLength < 1);
		my $fracLength = ($end-$start)/$ogLength;
		my $length = $end - $start;
		if ($rmskFlag) {
			next if ($fracLength < $minLengthPercent);
			next if ($fracLength > $maxLengthPercent);
		}
		next if ($length < $minLength);
		next if ($length > $maxLength);

		if ($ogFlag) {
			$start = $ogStart;
			$end = $ogEnd;
		}

		if ($p5Flag || $p3Flag) {
			my $p = $start;
			if ($strand eq '-') {
				$p = $end;
			}
			if ($p3Flag) {
				$p = $end;
				if ($strand eq '-') {
					$p = $start;
				}
			}
			$start = $p-$defHalfSize;
			$end = $p+$defHalfSize;
		}
		my $name = $line[0];
		if (exists($repeatIDs{$name})) {
			$name .= "-HOMER" . $repeatIDs{$name}++;
		} else {
			$repeatIDs{$name} = 1;
		}
		$totalRepeatsKept++;

		if ($countMethod ne 'genes' && $countMethod ne 'pausing' && @line > 5) {
			$ogStart = $start;
			$ogEnd = $end;
			my $features = getFeatureArray($strand,$end,$line[5]);
			my $numKept = 0;
			foreach(@$features) {
				my $keep = 0;
				if ($countMethod eq 'exons') {
					if ($_->{'n'} =~ /^E/) {
						$keep = 1;
					}
				} elsif ($countMethod eq '5utr') {
					if ($_->{'n'} =~ /5UTR/) {
						$keep = 1;
					}
				} elsif ($countMethod eq '3utr') {
					if ($_->{'n'} =~ /3UTR/) {
						$keep = 1;
					}
				} elsif ($countMethod eq 'cds') {
					if ($_->{'n'} =~ /^E\d+$/) {
						$keep = 1;
					}
				} elsif ($countMethod eq 'introns') {
					if ($_->{'n'} =~ /^I/) {
						$keep = 1;
					}
				}
				next if ($keep == 0);
				$numKept++;
				my $partName = "$name--Part" . $numKept;
				$start = $_->{'s'};
				$end = $_->{'e'};
				print $tmpFileHandle "$partName\t$chr\t$start\t$end\t$strand\tE:$start\t$div\t$ogStart\t$ogEnd\n";
				if ($readIntoMemory) {
					$repeats->{$partName} = {c=>$chr,s=>$start,e=>$end,d=>$strand,div=>$div,ogs=>$ogStart,
										oge=>$ogEnd,len=>$end-$start,n=>1,nexons=>1};
				}
				
			}
		} elsif ($countMethod eq 'pausing') {
			my $ps = $start + $pausePromoterStart;
			my $pe = $start + $pausePromoterEnd;
			my $gs = $start + $pauseGeneBodyStart;
			my $ge = $start + $pauseGeneBodyEnd;
			$ogStart = $ps;
			$ogEnd = $ge;
			if ($strand eq '-') {
				$pe = $end - $pausePromoterStart;
				$ps = $end - $pausePromoterEnd;
				$ge = $end - $pauseGeneBodyStart;
				$gs = $end - $pauseGeneBodyEnd;
				$ogStart = $gs;
				$ogEnd = $pe;
			}
			print $tmpFileHandle "$name--Part1\t$chr\t$ps\t$pe\t$strand\tE:$ps\t$div\t$ogStart\t$ogEnd\n";
			print $tmpFileHandle "$name--Part2\t$chr\t$gs\t$ge\t$strand\tE:$gs\t$div\t$ogStart\t$ogEnd\n";
			if ($readIntoMemory) {
				$repeats->{"$name--Part1"} = {c=>$chr,s=>$ps,e=>$pe,d=>$strand,div=>$div,ogs=>$ogStart,
									oge=>$ogEnd,len=>$pe-$ps,n=>1,nexons=>1};
				$repeats->{"$name--Part2"} = {c=>$chr,s=>$gs,e=>$ge,d=>$strand,div=>$div,ogs=>$ogStart,
									oge=>$ogEnd,len=>$ge-$gs,n=>1,nexons=>1};
			}
		} else {
			print $tmpFileHandle "$name\t$chr\t$start\t$end\t$strand\tE:$start\t$div\t$ogStart\t$ogEnd\n";
			if ($readIntoMemory) {
				$repeats->{$name} = {c=>$chr,s=>$start,e=>$end,d=>$strand,div=>$div,ogs=>$ogStart,
									oge=>$ogEnd,len=>$end-$start,n=>1,nexons=>1};
			}
		}

	}
	close IN;
}
sub getFeatureArray {
	my ($strand, $end, $featureStr) = @_;
	my @rawfeatures = split /\,/, $featureStr;
	my @features = ();
	foreach(@rawfeatures) {
		my @info = split /\:/;
		my $f = {n=>$info[0],s=>$info[1],e=>''};
		push(@features, $f);
	}
	@features = sort {$a->{'s'} <=> $b->{'s'}} @features;
	for (my $i=0;$i<@features-1;$i++) {
		$features[$i]->{'e'} = $features[$i+1]->{'s'}-1;
	}
	$features[@features-1]->{'e'} = $end;
	return \@features;
}
sub readRMSKfile {
	my ($file, $tmpFileHandle, $level1, $level2, $level3, $ogFlag, $p5flag, $p3flag, $defHalfSize, 
							$maxDiv, $minDiv, $repeats,$readIntoMemory) = @_;
	open IN, $file;
	while (<IN>) {
		chomp;
		next if (/^\s*#/);
		$totalRepeats++;
		my @line = split /\t/;
		next if ($level1 ne '' && $level1 ne $line[10]);
		next if ($level2 ne '' && $level2 ne $line[11]);
		next if ($level3 ne '' && $level3 ne $line[12]);
		my $div = $line[2]/1000.0;
		next if ($div > $maxDiv);
		next if ($div < $minDiv);
		my $chr = $line[5];
		my $start = $line[6]+1;
		my $end = $line[7];
		my $strand = $line[9];

		my $ogStart = $start-$line[13];
		my $ogEnd = $end+$line[15];
		if ($strand eq '-') {
			$ogStart = $start+$line[13];
			$ogEnd = $end+$line[15];
		}
		my $ogLength  = $ogEnd-$ogStart;
		$ogLength = 1 if ($ogLength < 1);
		my $fracLength = ($end-$start)/$ogLength;
		my $length = $end - $start;
		next if ($fracLength < $minLengthPercent);
		next if ($fracLength > $maxLengthPercent);
		next if ($length < $minLength);
		next if ($length > $maxLength);

		if ($ogFlag) {
			$start = $ogStart;
			$end = $ogEnd;
		}

		if ($p5Flag || $p3Flag) {
			my $p = $start;
			if ($strand eq '-') {
				$p = $end;
			}
			if ($p3Flag) {
				$p = $end;
				if ($strand eq '-') {
					$p = $start;
				}
			}
			$start = $p-$defHalfSize;
			$end = $p+$defHalfSize;
		}

		my $name = "$line[10]|$line[11]|$line[12]";
		if (exists($repeatIDs{$name})) {
			$name .= "-HOMER" . $repeatIDs{$name}++;
		} else {
			$repeatIDs{$name} = 1;
		}
		$totalRepeatsKept++;
		print $tmpFileHandle "$name\t$chr\t$start\t$end\t$strand\tE:$start\t$div\t$ogStart\t$ogEnd\n";
		if ($readIntoMemory) {
			$repeats->{$name} = {c=>$chr,s=>$start,e=>$end,d=>$strand,div=>$div,ogs=>$ogStart,oge=>$ogEnd,
								len=>$end-$start,n=>1,nexons=>1};
		}

	}
	close IN;
}
	
