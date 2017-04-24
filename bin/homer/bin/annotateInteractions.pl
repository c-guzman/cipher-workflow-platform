#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";

# Copyright 2009-2016 Christopher Benner <cbenner@ucsd.edu>
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
my $pseudo = 0.5;

sub printCMD {
	print STDERR "\n\tUsage: annotateInteractions.pl <interaction file> <genome version> <output directory>\n";
	print STDERR "\t\t\t[additional options]\n";
	print STDERR "\n\tGeneral Options:\n";
	print STDERR "\t\t-res <#> (Resolution of analysis, default: auto detect)\n";
	print STDERR "\t\t-hubCount <#> (Minimum number of interactions to define a hub, default: 5)\n";
	print STDERR "\t\t-skipann (skip all endpoint annotations)\n";
	print STDERR "\t\t-noann (skip detailed annotation of endpoints)\n";
	print STDERR "\t\t-washu (create Wash U Epigenome Browser output files: outDirName.bed.gz and tbi)\n";
	print STDERR "\t\t-cpu # (number of cores to use)\n";

	print STDERR "\n\tFiltering Options:\n";
	print STDERR "\t\t-minDist <#> (filter out interactions spaced less than # bp - set high for only interchr)\n";
	print STDERR "\t\t-maxDist <#> (filter out interactions spaced more than # bp, will remove interchr)\n";
	print STDERR "\t\t-pvalue <#> (filter out interactions with p-value greater than #)\n";
	print STDERR "\t\t\t-dpvalue <#> (filter out interactions with p-value (vs bg) greater than #)\n";
	print STDERR "\t\t-zscore <#> (filter out interactions with zscore less than #)\n";
	print STDERR "\t\t\t-dzscore <#> (filter out interactions with zscore (vs bg) less than #)\n";
	print STDERR "\t\t-filter <peakfile> (only look at interactions with endpoints in peakfile)\n";
	print STDERR "\t\t\t-filter2 <peakfile2> (only look at interactions connecting -filter and -filter2 peak files)\n";
	#print STDERR "\t\t-rmpreproB (remove interactions between chr7 and chr12)\n";
	#
	print STDERR "\n\tEnrichment Options:\n";
	print STDERR "\t\t-p <peak file 1> [peak file 2] ... (Check overlap with peak files)\n";
#	print STDERR "\t\t-cytoscape (Creates cytoscape.* files to load into Cytoscape, requires -p <...> too)\n";
#
	print STDERR "\n\tAssessing Interactions across Hi-C Experiments:\n";
	print STDERR "\t\t-d <Hi-C TagDirectory> [2nd Hi-C TagDirectory] ...\n";

	print STDERR "\n\tSpecial Operations:\n";
	print STDERR "\t\t-circos (Convert interactions to circos interactions format - stdout)\n";
	print STDERR "\t\t-i <interaction file2> [interaction file3] ... (Compare 1st file interactions to these)\n";
	print STDERR "\t\t-connect <peakFile1> <peakFile2> (returns association table between sets of peaks)\n";
	print STDERR "\t\t-pout (Convert interactions to a non-redundant peak file, sent to stdout)\n";

	print STDERR "\n\tSpecifying Background (i.e. regions used to find interactions - default: whole genome)\n";
	print STDERR "\t\t-gsize <#> (size of genome, default: 2e9)\n";
	print STDERR "\t\t-pos chrN:XXX-YYY (specific, continuous region)\n";
	print STDERR "\t\t-bgp <peak file> (peak file)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 3) {
	printCMD();
}

$interactionFile = $ARGV[0];
my $genome = $ARGV[1];
$washuFlag = 0;
$outputDirectory = $ARGV[2];
my $res = -1;
my $annFlag = 1;
@peakFiles = ();
@hicDirs = ();
@interactionFiles = ();
my $poutFlag = 0;
$gsizeOpt = "";
$minDist = -1;
$maxDist = -1;
$gsize = 'default';
my $circosFlag = 0;
my $hubCount = 5;
my $cytoscapePrefix = "";
$cytoscapePrefix = "cytoscape";
my $rmpreproB = 0;
$pvalueCutoff = 1e10;
$zscoreCutoff = -1e10;
$dpvalueCutoff = 1e10;
$dzscoreCutoff = -1e10;
$interactionHeader = '';
my $filterPeakFile1 = "";
my $filterPeakFile2 = "";
$numLinesInteractionFile = 0;
$bgChr = '';
$bgStart = -1e10;
$bgEnd = 1e10;
$bgPeakFile = '';
$annOptions = "";
my $tolerance = -1;
my $connectFlag= 0;
my $connectPeakFile1 = "";
my $connectPeakFile2 = "";
my $numCPUs = 1;


for (my $i=3;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-res') {
		$res = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-noann') {
		$annOptions .= " -noann ";
	} elsif ($ARGV[$i] eq '-skipann') {
		$annFlag = 0;
	} elsif ($ARGV[$i] eq '-pout') {
		$poutFlag=1;
	} elsif ($ARGV[$i] eq '-pvalue') {
		$pvalueCutoff = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$numCPUs = $ARGV[++$i]
	} elsif ($ARGV[$i] eq '-washu') {
		$washuFlag = 1;
	} elsif ($ARGV[$i] eq '-dpvalue') {
		$dpvalueCutoff = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-hubCount') {
		$hubCount = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-zscore') {
		$zscoreCutoff = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-dzscore') {
		$dzscoreCutoff = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cytoscape') {
		#$cytoscapePrefix = $ARGV[++$i];
		$cytoscapePrefix = "cytoscape";
	} elsif ($ARGV[$i] eq '-rmpreproB') {
		$rmpreproB = 1;
	} elsif ($ARGV[$i] eq '-bgp') {
		$bgPeakFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pos') {
		($bgChr,$bgStart,$bgEnd) = HomerConfig::parseUCSCStr($ARGV[++$i]);
	} elsif ($ARGV[$i] eq '-circos') {
		$circosFlag=1;
	} elsif ($ARGV[$i] eq '-gsize') {
		$gsize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-filter' || $ARGV[$i] eq '-filter1') {
		$filterPeakFile1=$ARGV[++$i];
	} elsif ($ARGV[$i] eq '-filter2') {
		$filterPeakFile2=$ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxDist') {
		$maxDist=$ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minDist') {
		$minDist=$ARGV[++$i];
	} elsif ($ARGV[$i] eq '-connect') {
		$connectFlag= 1;
		$connectPeakFile1 = $ARGV[++$i];
		$connectPeakFile2 = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-d') {
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@hicDirs, $ARGV[$i]);
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-p') {
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@peakFiles, $ARGV[$i]);
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-i') {
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@interactionFiles, $ARGV[$i]);
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} else {
		print STDERR "!! What is \"$ARGV[$i]\"??\n";
		printCMD();
	}
}



my $rand = rand();
my $tmpFile1 = $rand . ".1.tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";
my $tmpFile4 = $rand . ".4.tmp";
my $tmpFile5 = $rand . ".5.tmp";
my $intPeakFile = $rand . ".intpeaks.tmp";
my $leftPeakFile = $rand . ".left.tmp";
my $rightPeakFile = $rand . ".right.tmp";
$bgChoppedFile = $rand . ".bgChopped.tmp";

#check if bgzip and tabix are in PATH
if ($washuFlag) {
	my $whichOutput = `which bgzip 2> /dev/null`;
	if ($whichOutput eq '') {
		print STDERR "!!! Could not execute program bgzip needed for -washu option!!! \n";
		print STDERR "!!! Please download and install it (google for 'bgzip tabix download')!!! \n";
		print STDERR "!!! Also make sure the programs are included in your PATH !!!\n";
		exit;
	}
	$whichOutput = `which tabix 2> /dev/null`;
	if ($whichOutput eq '') {
		print STDERR "!!! Could not execute program tabix needed for -washu option!!! \n";
		print STDERR "!!! Please download and install it (google for 'bgzip tabix download')!!! \n";
		print STDERR "!!! Also make sure the programs are included in your PATH !!! \n";
		exit;
	}
}

if ($outputDirectory ne '') {
	`mkdir -p "$outputDirectory"`;
	$intPeakFile = "$outputDirectory/interactionPeaks.txt";
	$interactionsOutputFile = "$outputDirectory/interactions.txt";
	$interactionsOutputBed = "$outputDirectory/interactions.bed";
	$bedGraphFilename = "$outputDirectory/endpoint.bedGraph";

	$washuFilename = "$outputDirectory";
	$washuFilename =~ s/\///g;
	$washuFilename = "$outputDirectory/$washuFilename.bed";
}

$trackname="$outputDirectory interaction endpoints";

my @resPeakFiles = ();

my $organism = "none";
my $customGenome = 0;
my $genomeDir = "";
my $genomeParseDir = "";
my $consDir = "";
my $promoter = "";

if ($genome eq 'none') {
	print STDERR "\tNo genome selected (\"none\") - limits what you can do\n";
	$customGenome = -1;
	$genomeDir = "none";
} elsif (!exists($config->{'GENOMES'}->{$genome})) {
	$customGenome = 1;
	($genome,$genomeDir,$genomeParseDir) = HomerConfig::parseCustomGenome($genome);
} else {
	$genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'};
	$organism = $config->{'GENOMES'}->{$genome}->{'org'};
	$promoter = $config->{'GENOMES'}->{$genome}->{'promoters'};
	$consDir = $config->{'GENOMES'}->{$genome}->{'directory'} . "/conservation/";
}

print STDERR "\tInteraction file = $interactionFile\n";
print STDERR "\tGenome = $genome\n";
print STDERR "\tOrganism = $organism\n";


my $ogRes = $res;
my $interactions = readInteractions($interactionFile);
if ($ogRes == -1) {
	print STDERR "\tResolution set to $res\n";
} else {
	print STDERR "\tResolution set to $ogRes (auto detect guesses $res)\n";
	$res = $ogRes;
}
if (@hicDirs > 0) {
	my $dir1 = $hicDirs[0];
	my $options = '';
	$options .= " -maxDist $maxDist" if ($maxDist > 0);
	$options .= " -minDist $minDist" if ($minDist > 0);
	$options .= " -pvalue $pvalueCutoff" if ($pvalueCutoff < 2);
	$options .= " -zscore $zscoreCutoff" if ($zscoreCutoff > 0);
	`analyzeHiC "$hicDirs[0]" -res $res -i "$interactionFile" -cpu $numCPUs $options -interactions "$tmpFile1" > /dev/null`;
	#change the interactions to reflect what was just analyzed...
	$interactions = readInteractions($tmpFile1);
	if ($ogRes == -1) {
	} else {
		$res = $ogRes;
	}

	if (@hicDirs > 1) {
		readInteractions($tmpFile1,$interactions,$hicDirs[0]);
		for (my $i=1;$i<@hicDirs;$i++) {
			`analyzeHiC "$hicDirs[$i]" -res $res -i "$interactionFile" -cpu $numCPUs $options -interactions "$tmpFile1" > /dev/null`;
			readInteractions($tmpFile1,$interactions,$hicDirs[$i]);
		}
	}
	`rm "$tmpFile1"`;
}



my $halfRes = $res/2;
if ($tolerance < 0) {
	$tolerance = $halfRes;
}
my %annIndex = ();
my $annNum = 0;

printInteractionPeaks($interactions,$leftPeakFile,$rightPeakFile,$intPeakFile);	
if ($bgPeakFile ne '') {
	$interactions = filterInteractions($interactions, $bgPeakFile);
	printInteractionPeaks($interactions,$leftPeakFile,$rightPeakFile,$intPeakFile);	
}
if ($filterPeakFile1 ne '') {
	$interactions = filterPeakInteractions($interactions, $filterPeakFile1, $filterPeakFile2);
	printInteractionPeaks($interactions,$leftPeakFile,$rightPeakFile,$intPeakFile);	
}
`cp "$intPeakFile" "$outputDirectory/peaks.txt"`;

if (@interactionFiles > 0) {
	my $common = $interactions;
	print STDERR "\n\tChecking for common Interactions:\n";
	for (my $i=0;$i<@interactionFiles;$i++) {
		my $ogRes = $res;
		my $inters = readInteractions($interactionFiles[$i]);
		$res = $ogRes;
		my $outSuffix = $interactionFiles[$i];
		$outSuffix =~ s/\//\_/g;
		my $outputFile = $outputDirectory . "/commonWith_" . $outSuffix;
		my $cInteractions = getIntersectingInteractions($interactions, $inters, $tolerance,$outputFile);
		$common = getIntersectingInteractions($common, $inters, $tolerance,"");
	}
	$interactions = $common;
}

if ($gsize eq 'default') {
	if ($bgChr ne '') {
		$gsize = $bgEnd - $bgEnd;
		$gsizeOpt = "-gsize " . $gsize;
	} elsif ($bgPeakFile) {
		`mergePeaks -d given $bgChoppedFile -coverage "$tmpFile1" > /dev/null 2> /dev/null`;
		open IN, $tmpFile1;
		my $count = 0;
		while (<IN>) {
			$count++;
			if ($count == 2) {
				my @line = split /\t/;
				$gsize = $line[2];
			}
		}
		close IN;
		`rm "$tmpFile1"`;
		$gsizeOpt = "-gsize " . $gsize;
	} else {
	}
} else {
	$gsizeOpt = "-gsize " . $gsize;
}
print STDERR "\tEffective genome size = $gsize\n";
printInteractions($interactions, $interactionsOutputFile, \@hicDirs);
printInteractionsBed($interactions, $interactionsOutputBed);


if ($washuFlag) {
	createWashUFiles($interactions,$washuFilename);
}


makePeakBedGraph($interactions, $bedGraphFilename, $trackname,$hubCount);

if ($circosFlag) {
	printCircos($interactions);
	`rm -f "$leftPeakFile" "$rightPeakFile" "$intPeakFile"`;
	exit;
}
if ($connectFlag) {
	connectPeaks($interactions,$connectPeakFile1,$connectPeakFile2);
	`rm -f "$leftPeakFile" "$rightPeakFile" "$intPeakFile"`;
	exit;
}
if ($outputDirectory ne '') {
	getLengthDistribution($interactions);
}


if ($poutFlag) {
	open IN, "$intPeakFile";
	while (<IN>) {
		print $_;
	}
	close IN;
	`rm -f "$leftPeakFile" "$rightPeakFile" "$intPeakFile"`;
	exit;
}

if ($hubCount > 0 && $outputDirectory ne '') {
	printHUBs($intPeakFile, $hubCount);
}


if (@peakFiles > 0) {
	getPeakFileEnrichment($interactions);
}

if ($annFlag) {
	annotatePeakFiles($interactions,$leftPeakFile,$rightPeakFile);
	printInteractionAnnotation($interactions, "", \@hicDirs);
}

#print STDERR "-=========================== $outputDirectory/interactionPeaks.txt \n";
#`annotatePeaks.pl "$intPeakFile" $genome > "$outputDirectory/interactionPeaks.txt"`;

`rm -f "$leftPeakFile" "$rightPeakFile" "$intPeakFile"`;

exit;




sub createWashUFiles {
	my ($interactions,$washuFilename) = @_;

	print STDERR "\tCreating WashU files...\n";
	my %washu = ();
	my $wid = 1;
	foreach(values %$interactions) {
		my $wid1 = $wid++;
		my $wid2 = $wid++;
		my $v = $_->{'lp'};
		my $str1 = "$_->{'c1'}\t$_->{'s1'}\t$_->{'e1'}\t$_->{'c2'}:$_->{'s2'}"."-"."$_->{'e2'},$v\t$wid1\t.";
		my $str2 = "$_->{'c2'}\t$_->{'s2'}\t$_->{'e2'}\t$_->{'c1'}:$_->{'s1'}"."-"."$_->{'e1'},$v\t$wid2\t.";

		$washu{$wid1} = {c=>$_->{'c1'},s=>$_->{'s1'},v=>$str1};
		$washu{$wid2} = {c=>$_->{'c2'},s=>$_->{'s2'},v=>$str2};
	}
	my @ids = sort {$washu{$a}->{'c'} cmp $washu{$b}->{'c'} ||
						$washu{$a}->{'s'} <=> $washu{$b}->{'s'}} keys %washu;
	
	open OUT, ">$washuFilename";
	foreach(@ids) {
		print OUT "$washu{$_}->{'v'}\n";
	}
	close OUT;
	
	`bgzip -f "$washuFilename"`;
	`tabix -f -p bed "$washuFilename.gz"`;
	print STDERR "\t\tFile: $washuFilename.gz\n";
		
}


sub makePeakBedGraph {
	my ($interactions, $filename, $trackname,$hubThreshold) = @_;


	my $v = 1;

	my %lastPos = ();
	my %diffGraph = ();
	foreach(keys %$interactions) {
		my $c1 = $interactions->{$_}->{'c1'};
		my $s1 = $interactions->{$_}->{'s1'};
		my $e1 = $interactions->{$_}->{'e1'};
		if (!exists($diffGraph{$c1})) {
			my %a = ();
			$diffGraph{$c1}=\%a;
			$lastPos{$c1}=0;
		}		
		if (!exists($diffGraph{$c1}->{$s1})) {
			$diffGraph{$c1}->{$s1}=0;
		}
		$diffGraph{$c1}->{$s1}+=$v;
		$lastPos{$c1} = $s1 if ($s1 > $lastPos{$c1});
		if (!exists($diffGraph{$c1}->{$e1})) {
			$diffGraph{$c1}->{$e1}=0;
		}
		$diffGraph{$c1}->{$e1}-=$v;

		my $c2 = $interactions->{$_}->{'c2'};
		my $s2 = $interactions->{$_}->{'s2'};
		my $e2 = $interactions->{$_}->{'e2'};
		if (!exists($diffGraph{$c2})) {
			my %a = ();
			$diffGraph{$c2}=\%a;
		}		
		if (!exists($diffGraph{$c2}->{$s2})) {
			$diffGraph{$c2}->{$s2}=0;
		}
		if (!exists($lastPos{$c2})) {
			$lastPos{$c2}=$s2;
		} else {
			$lastPos{$c2} = $s2 if ($s2 > $lastPos{$c2});
		}
		$diffGraph{$c2}->{$s2}+=$v;
		if (!exists($diffGraph{$c2}->{$e2})) {
			$diffGraph{$c2}->{$e2}=0;
		}
		$diffGraph{$c2}->{$e2}-=$v;
	}
	open OUT, ">$filename";
	open PEAKS, ">$filename.peaks";
	print OUT "track name=\"$trackname\" type=bedGraph\n";
	foreach(keys %diffGraph) {
		my $chr = $_;
		my @pos = sort {$a <=> $b} keys %{$diffGraph{$chr}};
		my $value = 0;
		my $aboveFlag = 0;
		my $max = 0;
		my $str = '';
		for (my $i=0;$i<@pos-1;$i++) {
			$value += $diffGraph{$chr}->{$pos[$i]};
			last if ($pos[$i] > $lastPos{$chr});	
			print OUT "$chr\t$pos[$i]\t$pos[$i+1]\t$value\n";

			if ($value > $hubThreshold) {
				$aboveFlag = 1;
				if ($value > $max) {
					$max = $value;
					$str = "$chr-$pos[$i]\t$chr\t$pos[$i]\t$pos[$i+1]\t+\t$value\n";
				}
			} else {
				if ($aboveFlag) {
					print PEAKS $str;
				}
				$aboveFlag = 0;
				$max = 0;
				$str= "";
			}
		}
	}
	close OUT;
	close PEAKS;

}

sub printInteractions {

	my ($interactions, $filename, $hicDirs) = @_;

	if ($filename eq '') {
		$filename = "$outputDirectory/interactionAnnotation.txt";
	}

	open OUT, ">$filename";
	if ($interactionHeader ne '') {
		print OUT "$interactionHeader";
	} else {
		print OUT "Interaction\tPeakID(1)\tchr(1)\tstart(1)\tend(1)\tstrand(1)\tTotal Reads(1)";
		print OUT "\tPeakID(2)\tchr(2)\tstart(2)\tend(2)\tstrand(2)\tTotal Reads(2)";
		print OUT "\tDistance\tInteraction Reads\tExpected Reads\tZ-score\tLogP\tFDR";
		for (my $i=17;$i<$numLinesInteractionFile;$i++) {
			print OUT "\tdata";
		}
	}
	foreach(@$hicDirs) {
		my $name = $_;
		print OUT "\t$name Total Reads(1)\t$name Total Reads(2)\t$name Interaction Reads"
				. "\t$name Expected Reads\t$name Log2 Ratio\t$name Z-score\t$name LogP";
	}
	print OUT "\n";


	my @names = sort {$interactions->{$a}->{'lp'} <=> $interactions->{$b}->{'lp'} ||
						$interactions->{$b}->{'v'} <=> $interactions->{$a}->{'v'} } keys %$interactions;
	foreach(@names) {
		my $n = $_;
		my $og = $interactions->{$n}->{'og'};
		print OUT "$og";
		foreach(@$hicDirs) {
			my $name = $_;
			if (exists($interactions->{$n}->{'hic'}->{$name})) {
				my $dd = $interactions->{$n}->{'hic'}->{$name};
				my $lr = log(($pseudo+$dd->{'v'})/($dd->{'ex'}+$pseudo))/log(2);
				print OUT "\t$dd->{'v1'}\t$dd->{'v2'}\t$dd->{'v'}\t$dd->{'ex'}\t$lr\t$dd->{'z'}"
						. "\t$dd->{'lp'}";
			} else {
				print OUT "\t\t\t\t\t\t";
			}
		}
		print OUT "\n";
	}
	close OUT;
}

sub printInteractionsBed {

	my ($interactions, $filename) = @_;

	if ($filename eq '') {
		print STDERR "!! no filename specified in printInteractionsBed !!\n";
		return;
	}

	open OUT, ">$filename";

	print OUT "track name=\"$interactionFile\"\n";

	my @names = sort {$interactions->{$a}->{'lp'} <=> $interactions->{$b}->{'lp'} ||
						$interactions->{$b}->{'v'} <=> $interactions->{$a}->{'v'} } keys %$interactions;
	foreach(@names) {
		my $n = $_;
		my $int = $interactions->{$n};
		my $score = $int->{'lp'}*-10;
		my $intname = $n;

		if ($int->{'c1'} eq $int->{'c2'}) {
			my $low = 0;
			if ($int->{'s1'} <= $int->{'s2'}) {
				print OUT "$int->{'c1'}\t$int->{'s1'}\t$int->{'e2'}\t$intname\t$score\t$int->{'d1'}"
						. "\t$int->{'s1'}\t$int->{'s1'}\t0";
				if ($int->{'e1'} < $int->{'s2'}) {
					my $size1 = $int->{'e1'}-$int->{'s1'};
					my $size2 = $int->{'e2'}-$int->{'s2'};
					my $start1 = 0;
					my $start2 = $int->{'s2'}-$int->{'s1'};
					print OUT "\t2\t$size1,$size2\t$start1,$start2\n";
				} else {
					my $size1 = $int->{'e2'}-$int->{'s1'};
					print OUT "\t1\t$size1,\t0\n";
				}
			} else {
				print OUT "$int->{'c1'}\t$int->{'s2'}\t$int->{'e1'}\t$intname\t$score\t$int->{'d2'}"
						. "\t$int->{'s2'}\t$int->{'s2'}\t0";
				if ($int->{'e2'} < $int->{'s1'}) {
					my $size1 = $int->{'e2'}-$int->{'s2'};
					my $size2 = $int->{'e1'}-$int->{'s1'};
					my $start1 = 0;
					my $start2 = $int->{'s1'}-$int->{'s2'};
					print OUT "\t2\t$size1,$size2\t$start1,$start2\n";
				} else {
					my $size1 = $int->{'e1'}-$int->{'s2'};
					print OUT "\t1\t$size1,\t0\n";
				}
			}
		} else {
			my $size1 = $int->{'e1'}-$int->{'s1'};
			$intName = $n . '->' . "$int->{'c2'}:$int->{'s2'}" . "-" . $int->{'e2'};
			print OUT "$int->{'c1'}\t$int->{'s1'}\t$int->{'e1'}\t$intname\t$score\t$int->{'d1'}"
					. "\t$int->{'s1'}\t$int->{'s1'}\t0\t1\t$size1\t0\n";

			$size1 = $int->{'e2'}-$int->{'s2'};
			$intName = $n . '->' . "$int->{'c1'}:$int->{'s1'}" . "-" . $int->{'e1'};
			print OUT "$int->{'c2'}\t$int->{'s2'}\t$int->{'e2'}\t$intname\t$score\t$int->{'d2'}"
					. "\t$int->{'s2'}\t$int->{'s2'}\t0\t1\t$size1\t0\n";
			
		}
	}
	close OUT;
}


sub printInteractionAnnotation {

	my ($interactions, $filename, $hicDirs) = @_;

	if ($filename eq '') {
		$filename = "$outputDirectory/interactionAnnotation.txt";
	}

	open OUT, ">$filename";
	if ($interactionHeader ne '') {
		print OUT "$interactionHeader";
	} else {
		print OUT "Interaction\tPeakID(1)\tchr(1)\tstart(1)\tend(1)\tstrand(1)\tTotal Reads(1)";
		print OUT "\tPeakID(2)\tchr(2)\tstart(2)\tend(2)\tstrand(2)\tTotal Reads(2)";
		print OUT "\tDistance\tInteraction Reads\tExpected Reads\tZ-score\tLogP\tFDR";
		for (my $i=17;$i<$numLinesInteractionFile;$i++) {
			print OUT "\tdata";
		}
	}
	print OUT "\tTotal Number of Significant Interactions at region(1)";
	print OUT "\tTotal Number of Significant Interactions at region(2)";

	print OUT "\tAnnotation(1)\tDetailed Annotation(1)\tDistance to TSS(1)\tNearest PromoterID(1)";
	print OUT "\tGene Name(1)\tGene Alias(1)\tGene Description(1)";
	print OUT "\tAnnotation(2)\tDetailed Annotation(2)\tDistance to TSS(2)\tNearest PromoterID(2)";
	print OUT "\tGene Name(2)\tGene Alias(2)\tGene Description(2)";
	if (@peakFiles > 0) {
		print OUT "\tPeak Links";
	}
	foreach(@$hicDirs) {
		my $name = $_;
		print OUT "\t$name Total Reads(1)\t$name Total Reads(2)\t$name Interaction Reads"
				. "\t$name Expected Reads\t$name Log2 Ratio\t$name Z-score\t$name LogP";
	}
	print OUT "\n";


	my @names = sort {$interactions->{$a}->{'lp'} <=> $interactions->{$b}->{'lp'} ||
						$interactions->{$b}->{'v'} <=> $interactions->{$a}->{'v'} } keys %$interactions;
	foreach(@names) {
		my $n = $_;
		my $og = $interactions->{$n}->{'og'};
		print OUT "$og";
		if (exists($interactions->{$n}->{'hub1'})) {
			print OUT "\t$interactions->{$n}->{'hub1'}";
		} else {
			print OUT "\t1";
		}
		if (exists($interactions->{$n}->{'hub2'})) {
			print OUT "\t$interactions->{$n}->{'hub2'}";
		} else {
			print OUT "\t1";
		}
		if (exists($interactions->{$n}->{'annstr1'})) {
			print OUT "\t$interactions->{$n}->{'annstr1'}";
		} else {
			print OUT "\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
		}
		if (exists($interactions->{$n}->{'annstr2'})) {
			print OUT "\t$interactions->{$n}->{'annstr2'}";
		} else {
			print OUT "\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
		}
		if (@peakFiles > 0) {
			print OUT "\t$interactions->{$n}->{'linkStr'}";
		}
		foreach(@$hicDirs) {
			my $name = $_;
			if (exists($interactions->{$n}->{'hic'}->{$name})) {
				my $dd = $interactions->{$n}->{'hic'}->{$name};
				my $lr = log(($pseudo+$dd->{'v'})/($dd->{'ex'}+$pseudo))/log(2);
				print OUT "\t$dd->{'v1'}\t$dd->{'v2'}\t$dd->{'v'}\t$dd->{'ex'}\t$lr\t$dd->{'z'}"
						. "\t$dd->{'lp'}";
			} else {
				print OUT "\t\t\t\t\t\t";
			}
		}


		print OUT "\n";
	}
	close OUT;
}

sub annotatePeakFiles {
	my ($interactions, $peakFile1, $peakFile2) = @_;

	print STDERR "\tAnnotating interaction endpoints...\n";
	`annotatePeaks.pl "$peakFile1" $genome $annOptions > $tmpFile3 2> /dev/null`;
	open IN, $tmpFile3;
	my $count = 0;
	while (<IN>) {
		$count++;
		next if ($count == 1);
		chomp;
		my @line = split /\t/;
		my $id = $line[0];
		while (@line < 18) {
			push(@line, "NA");
		}
		my $annstr = "$line[7]\t$line[8]\t$line[9]\t$line[10]\t$line[15]\t$line[16]\t$line[17]";
		if (exists($interactions->{$id})) {
			$interactions->{$id}->{'annstr1'} = $annstr;
		}
	}
	close IN;

	`annotatePeaks.pl "$peakFile2" $genome $annOptions > $tmpFile3 2> /dev/null`;
	open IN, $tmpFile3;
	$count = 0;
	while (<IN>) {
		$count++;
		next if ($count == 1);
		chomp;
		my @line = split /\t/;
		my $id = $line[0];
		while (@line < 18) {
			push(@line, "NA");
		}
		my $annstr = "$line[7]\t$line[8]\t$line[9]\t$line[10]\t$line[15]\t$line[16]\t$line[17]";
		if (exists($interactions->{$id})) {
			$interactions->{$id}->{'annstr2'} = $annstr;
		}
	}
	close IN;
	`rm -f $tmpFile3`;
}


sub getPeakFileEnrichment {

	my ($interactions) = @_;	


	for (my $i=0;$i<@peakFiles;$i++) {
		my $tmpPeakFile = "$rand.$i.peak.tmp";
		my $res1 = $res+1;
		`chopify.pl "$peakFiles[$i]" $res1 > "$tmpFile5" 2> /dev/null`;
		#print STDERR "`mergePeaks -d $res $tmpFile5 > $tmpPeakFile 2> /dev/null`;\n";
		`mergePeaks -d $res "$tmpFile5" > "$tmpPeakFile" 2> /dev/null`;
		`rm -f "$tmpFile5"`;
		push(@resPeakFiles,$tmpPeakFile);
	}

	for (my $i=0;$i<@peakFiles;$i++) {
		my $sname = getShortName($peakFiles[$i]);
		$annIndex{$annNum} = {name=>$peakFiles[$i],id=>$annNum,total=>"",file=>$resPeakFiles[$i],sname=>$sname};
		$annIndex{$annNum}->{'p'} = readPeakFile($annIndex{$annNum},$resPeakFiles[$i]);
		`mergePeaks -d $halfRes -cobound 1 "$leftPeakFile" "$resPeakFiles[$i]" -prefix "$tmpFile3" 2> /dev/null`;
		`mergePeaks -d $halfRes -cobound 1 "$rightPeakFile" "$resPeakFiles[$i]" -prefix "$tmpFile4" 2> /dev/null`;
		`mergePeaks -d $halfRes "$intPeakFile" "$resPeakFiles[$i]" -matrix "$tmpFile5" $gsizeOpt > /dev/null 2> /dev/null`;
		readInteractionOverlap($interactions,$annIndex{$annNum},"$tmpFile3.coBoundBy1.txt", 
								"$tmpFile4.coBoundBy1.txt","$tmpFile5.logPvalue.matrix.txt",
								"$tmpFile5.logRatio.matrix.txt");
		`rm -f "$tmpFile3.coBoundBy"* "$tmpFile4.coBoundBy"* "$tmpFile5"*`;
		printPeakInteractionPartners($interactions, $annIndex{$annNum}, "$outputDirectory/$sname.partners.txt");
		$annNum++;
	}

	my @annIndexes = sort {$a <=> $b} keys %annIndex;
	my $totalInteractions = 0;
	foreach (keys %$interactions) {
		$totalInteractions++;
		my $name = $_;
		my %sets = ();
		foreach(keys %{$interactions->{$name}->{'a1'}}) {
			my $a1 = $_;
			$interactions->{$name}->{'sets'} = \%sets;
			foreach(keys %{$interactions->{$name}->{'a2'}}) {
				my $a2 = $_;
				my $str = $a2 . "x" . $a1;
				if ($a1 <= $a2) {
					$str = $a1 . "x" . $a2;
				}
				$sets{$str} = 1;
			}
		}
		$interactions->{$name}->{'sets'} = \%sets;
	
		my $pos1 = $interactions->{$name}->{'c1'} . ":" . $interactions->{$name}->{'s1'} 
												. "-" . $interactions->{$name}->{'e1'};	
		my $pos2 = $interactions->{$name}->{'c2'} . ":" . $interactions->{$name}->{'s2'} 
												. "-" . $interactions->{$name}->{'e2'};	
		
		my $linkStr = "";
		foreach(keys %sets) {
			$linkStr .= "$_,";
		}
		$interactions->{$name}->{'linkStr'} = $linkStr;
	}
	
	
	if ($cytoscapePrefix ne '') {
		open CYTONET, ">$outputDirectory/$cytoscapePrefix.network.sif.txt";
		open CYTOEDGE, ">$outputDirectory/$cytoscapePrefix.edge.logp.attributes.txt";
		open CYTOEDGERATIO, ">$outputDirectory/$cytoscapePrefix.edge.ratio.attributes.txt";
		open CYTONODESIZE, ">$outputDirectory/$cytoscapePrefix.node.size.txt";
		open CYTONODELOGP, ">$outputDirectory/$cytoscapePrefix.node.logp.txt";
		open CYTONODERATIO, ">$outputDirectory/$cytoscapePrefix.node.ratio.txt";

		print CYTOEDGE "InteractionLogP\n";
		print CYTOEDGERATIO "InteractionRatio\n";
		print CYTONODESIZE "NodeSize\n";
		print CYTONODELOGP "NodeEnrichment\n";
		print CYTONODERATIO "NodeRatio\n";
	}
	
	open OUT, ">$outputDirectory/pairwiseFeatureEnrichment.txt";
		
	print OUT "FeatureSetID\tFeature1\tFeature2\tInteractions($totalInteractions)\tInteractions Expected\tEnrichment Ratio (log2)";
	print OUT "\tEnrichment logP (+values NOT enriched)\tP-value (- values NOT enriched)\tTotal Feature 1\tTotal Feature 2\tFeature Overlap";
	print OUT "\tFeature 1 General Enrichment\tFeature 2 General Enrichment\tTotal Interactions";
	print OUT "\n";

	open OUT2, ">$outputDirectory/featureEnrichment.txt";
	print OUT2 "Feature Name\tNumber of overlapping Features\tEnrichment Ratio (log2)\tEnrichment logP, +values NOT enriched)\tfile\n";

	for (my $i=0;$i<@annIndexes;$i++) {
	
		my $index1 = $annIndexes[$i];
		my $shortName1 = $annIndex{$index1}->{'sname'};
		my $longName1 = $annIndex{$index1}->{'name'};
		my $a1t1 = $annIndex{$index1}->{'t1'};
		my $a1t2 = $annIndex{$index1}->{'t2'};
		my $total1 = $a1t1+$a1t2;
		print OUT2 "$shortName1\t$total1\t$annIndex{$index1}->{'logRatio'}\t$annIndex{$index1}->{'enrichment'}\t$longName1\n";
		if ($cytoscapePrefix ne '') {
			print CYTONODESIZE "$shortName1 = $total1\n";
			print CYTONODELOGP "$shortName1 = $annIndex{$index1}->{'enrichment'}\n";
			print CYTONODERATIO "$shortName1 = $annIndex{$index1}->{'logRatio'}\n";
		}
		for (my $j=$i;$j<@annIndexes;$j++) {
			
			my $index1 = $annIndexes[$i];
			my $index2 = $annIndexes[$j];
			my $id = $index1 . "x" . $index2;
	
	
			my $name1 = $annIndex{$index1}->{'name'};
			my $name2 = $annIndex{$index2}->{'name'};
			my $shortName1 = $annIndex{$index1}->{'sname'};
			my $shortName2 = $annIndex{$index2}->{'sname'};
	
			print OUT "$id\t$shortName1\t$shortName2";
	
			my $featureLogP = 0.0;
			if ($i!=$j) {
				#$featureLogP = getPeakOverlapEnrichment($annIndex{$index1}->{'file'},$annIndex{$index2}->{'file'});
			}
	
			my $x = 0;
			my $N = 0;
			my $overlap = 0;
			foreach(values %$interactions) {
				if (exists($_->{'a1'}->{$index1}) && exists($_->{'a1'}->{$index2})) {
					$overlap++;
				}
				if (exists($_->{'a2'}->{$index1}) && exists($_->{'a2'}->{$index2})) {
					$overlap++;
				}
				if (exists($_->{'sets'}->{$id})) {
					$x++;
				}
				$N++;
			}
	
			my $a1t1 = $annIndex{$index1}->{'t1'};
			my $a1t2 = $annIndex{$index1}->{'t2'};
			my $a2t1 = $annIndex{$index2}->{'t1'};
			my $a2t2 = $annIndex{$index2}->{'t2'};
	
			my $rate1 = ($a1t1+$a1t2)/(2*$N);
			my $rate2 = ($a2t1+$a2t2)/(2*$N);
	
			my $total1 = $a1t1+$a1t2;
			my $total2 = $a2t1+$a2t2;
	
			my $expected = $rate1*$rate2*$N;
			if ($index1 ne $index2) {
				$expected = 0;
				my $overRate = $overlap/(2*$N);
				my $r1 = ($a1t1+$a1t2-$overlap)/(2*$N);
				my $r2 = ($a2t1+$a2t2-$overlap)/(2*$N);
	
				my $e1 = $overRate*($overRate+$r1+$r2)*$N;
				my $e2 = $r1*($overRate+$r2)*$N;
				my $e3 = $r2*($overRate+$r1)*$N;
	
				my $nn1 = ($a1t1+$a1t2-$overlap);
				my $nn2 = ($a2t1+$a2t2-$overlap);
				#print STDERR "E1=$e1\nE2=$e2\nE3=$e3\n";
				$expected = $e1+$e2+$e3;
				#print STDERR "Obs: $x vs. $expected ($nn1 $overlap $nn2) ($total1 vs. $total2 : $N [$index1 vs. $index2])\n";
			}
	
			my $logp = Statistics::logbinomial($N,$x,$expected/$N,$N);
			my $pvalue = exp($logp);
			if ($logp > log(0.5)) {
				#$logp = Statistics::ilogbinomial($N,$x,$expected/$N,$N);
				#print STDERR "logp=$logp\n";
				$logp = Statistics::ilogbinomial($N,$x,$expected/$N,$N);
				$pvalue = exp($logp);
				$logp *= -1;
				$pvalue *= -1;
			}
	
			my $LRatio = 0;	
			my $X1 = $x;
			my $X2 = $expected;
			$X1 = 1 if ($X1 == 0);
			$X2 = 1 if ($X2 == 0);
			$LRatio = log($X1/$X2);

			print OUT "\t$x\t$expected\t$LRatio";
			print OUT "\t$logp\t$pvalue";
			print OUT "\t$total1\t$total2\t$overlap";
			print OUT "\t" . $annIndex{$index1}->{'enrichment'};
			print OUT "\t" . $annIndex{$index2}->{'enrichment'};
			#print OUT "\t$featureLogP";
			print OUT "\t$N";
			print OUT "\n";
	
			if ($cytoscapePrefix ne '') {
				print CYTONET "$shortName1\tpp\t$shortName2\n";
				#print CYTONET "$shortName2\tpp\t$shortName1\n";
				print CYTOEDGE "$shortName1 (pp) $shortName2 = $logp\n";
				print CYTOEDGERATIO "$shortName1 (pp) $shortName2 = $LRatio\n";
				#print CYTOEDGE "$shortName2 (pp) $shortName1\n";
			}
	
		}
	}
	close OUT;
	close OUT2;

	if ($cytoscapePrefix ne '') {
		close CYTONET;
		close CYTOEDGE;
		close CYTOEDGERATIO;
		close CYTONODESIZE;
		close CYTONODELOGP;
		close CYTONODERATIO;
	}

	foreach(@resPeakFiles) {
		`rm -f "$_"`;
	}

	return $interactions;
}


########################################################################################

sub getShortName {
	my ($name) = @_;

	my $og = $name;	
	$name =~ s/\.ann.txt//;
	$name =~ s/peaks.txt//;
	$name =~ s/regions.txt//;
	$name =~ s/\.txt$//;
	$name =~ s/\.tsv$//;
	$name =~ s/^.*\/(.*\/.*)$/$1/;
	$name =~ s/\/+//;
	#print STDERR "\t\t$og -> $name\n";
	return $name;
}

sub getPeakOverlapEnrichment {
	my ($pfile1, $pfile2) = @_;
	`mergePeaks -d $halfRes "$pfile1" "$pfile2" -matrix "$tmpFile5" $gsizeOpt > /dev/null 2> /dev/null`;
	open IN, "$tmpFile5.logPvalue.matrix.txt";
	my $count = 0;
	my $logpvalue = 0;
	while (<IN>) {
		$count++;
		chomp;
		next if ($count < 2);
		s/\r//g;
		my @line = split /\t/;
		$logpvalue = $line[2];
		last;
	}
	close IN;
	`rm -f "$tmpFile5"*`;
	return $logpvalue;
}

sub readPeakFile {

	my ($ann, $peakFile) = @_;
	#$ann->{'name'} = $peakFile;
	my %data = ();
	open IN, $peakFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^\s*#/);
		my @line = split /\t/;
		my $name = $line[0];
		my $chr = $line[1];
		my $start = $line[2];
		my $end = $line[3];
		my $strand = $line[4];
		if (checkInt($start) || checkInt($end)) {
			next;
		}
		my $mid = floor(($start+$end)/2);
		$strand = '+' if ($strand eq '0');
		$strand = '-' if ($strand eq '1');
		$data{$name} = {s=>$start,e=>$end,d=>$strand,c=>$chr,m=>$mid};
	}
	$ann->{'peaks'} = \%data;
}

sub printPeakInteractionPartners {
	my ($ints, $ann, $outfile) = @_;
	my $annIndex = $ann->{'id'};
	open POUT, ">$outfile" or die "!!! Could not open file $outfile for writing!!!\n";
	print POUT "Interaction\tchr\tstart\tend\tstrand\tTotalReads\tDistance\tInteractions\tExpected Interactions";
	print POUT "\tZ-score\tLogP\n";
	foreach(keys %$ints) {
		my $iname = $_;
		my $left = 0;
		my $right = 0;
		if (exists($ints->{$iname}->{'a1'}->{$annIndex})) {
			$left = 1;
		}
		if (exists($ints->{$iname}->{'a2'}->{$annIndex})) {
			$right = 1;
		}
		next if ($left == $right);

		if ($right == 1) {
			print POUT "$iname\t$ints->{$iname}->{'c1'}\t$ints->{$iname}->{'s1'}\t$ints->{$iname}->{'e1'}"
					. "\t$ints->{$iname}->{'d1'}\t$ints->{$iname}->{'v1'}";
		} elsif ($left == 1) {
			print POUT "$iname\t$ints->{$iname}->{'c2'}\t$ints->{$iname}->{'s2'}\t$ints->{$iname}->{'e2'}"
					. "\t$ints->{$iname}->{'d2'}\t$ints->{$iname}->{'v2'}";

			#$interactions{$name} = {c1=>$chr1,s1=>$start1,e1=>$end1,d1=>$strand1,v1=>$v1,m1=>$m1,
			#			c2=>$chr2,s2=>$start2,e2=>$end2,d2=>$strand2,v2=>$v2,m2=>$m2,
			#			dist=>$dist,v=>$v,ex=>$expected,a1=>\%ann1,a2=>\%ann2,z=>$zscore,lp=>$logp,og=>$og,
			#			dz=>$dzscore,dlp=>$dlogp,th=>$thickness,hic=>\%hicExps};
		}
		print POUT "\t$ints->{$iname}->{'dist'}"
					. "\t$ints->{$iname}->{'v'}\t$ints->{$iname}->{'ex'}\t$ints->{$iname}->{'z'}"
					. "\t$ints->{$iname}->{'lp'}\n";
	}
	close POUT;
}


sub readInteractionOverlap {
	my ($ints,$ann,$file1, $file2, $matrixFile, $ratioFile) = @_;
	my $annIndex = $ann->{'id'};
	my $t1 = 0;
	open IN1, "$file1" or die "!!! Couldn't open file: $file1 for reading!!!\n";
	while (<IN1>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if (exists($ints->{$line[0]})) {
			$ints->{$line[0]}->{'a1'}->{$annIndex}=1;
			$t1++;
		}
	}
	close IN1;
	my $t2 = 0;
	open IN2, "$file2" or die "!!! Couldn't open file: $file2 for reading!!!\n";
	while (<IN2>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if (exists($ints->{$line[0]})) {
			$ints->{$line[0]}->{'a2'}->{$annIndex}=1;
			$t2++;
		}
	}
	close IN2;
	open IN, $matrixFile;
	my $count = 0;
	while (<IN>) {
		$count++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($count == 2) {
			$ann->{'enrichment'} = $line[2];
			print STDERR "\t$ann->{'name'} enrichment at interaction ends (logp): $line[2]\n";
			last;
		}
	}
	close IN;
	open IN, $ratioFile;
	$count = 0;
	while (<IN>) {
		$count++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($count == 2) {
			$ann->{'logRatio'} = $line[2];
			last;
		}
	}
	close IN;
	$ann->{'t1'} = $t1;
	$ann->{'t2'} = $t2;
}

sub printInteractionPeaks {
	my ($interactions, $leftPeakFile, $rightPeakFile, $intPeakFile) = @_;
	open OUT1, ">$leftPeakFile" or die "!!! Couldn't open file: $leftPeakFile for writing!!!\n";
	open OUT2, ">$rightPeakFile" or die "!!! Couldn't open file: $rightPeakFile for writing!!!\n";
	foreach(keys %$interactions) {
		print OUT1 "$_\t$interactions->{$_}->{'c1'}\t$interactions->{$_}->{'s1'}\t$interactions->{$_}->{'e1'}\t$interactions->{$_}->{'d1'}\n";
		print OUT2 "$_\t$interactions->{$_}->{'c2'}\t$interactions->{$_}->{'s2'}\t$interactions->{$_}->{'e2'}\t$interactions->{$_}->{'d2'}\n";
	}
	close OUT1;
	close OUT2;
	`mergePeaks -d $res "$leftPeakFile" "$rightPeakFile" > "$intPeakFile" 2> /dev/null`;

	my %hubs = ();
	open IN, $intPeakFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $v = 1;
		if (@line > 7) {
			$v = $line[7];
		}
		$hubs{$line[0]} = $v;
	}
	close IN;

	`annotateRelativePosition.pl "$leftPeakFile", "$intPeakFile", > "$tmpFile3"`;
	open IN, $tmpFile3;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $v = 1;
		if (exists($hubs{$line[1]})) {
			$v = $hubs{$line[1]};
		}
		if (exists($interactions->{$line[0]})) {
			$interactions->{$line[0]}->{'hub1'} = $v;
		}
	}
	close IN;

	`annotateRelativePosition.pl "$rightPeakFile", "$intPeakFile", > "$tmpFile3"`;
	open IN, $tmpFile3;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $v = 1;
		if (exists($hubs{$line[1]})) {
			$v = $hubs{$line[1]};
		}
		if (exists($interactions->{$line[0]})) {
			$interactions->{$line[0]}->{'hub2'} = $v;
		}
	}
	close IN;
	`rm -f "$tmpFile3"`;

}


sub readInteractions {
	my ($file,$ints,$hicDir) = @_;

	my $addFlag = 0;
	if (defined($ints) && defined($hicDir)) {
		$addFlag = 1;
		print STDERR "\tAdding interaction information for $hicDir...\n\n";
	}
	my %interactions = ();
	my %uniq = ();
	my $total = 0;
	my $N = 0;
	my $totalInteractions = 0;
	my $locationFiltered=0;;
	open IN, $file or die "Could not open interaction file ($file)!!!\n";
	my $count = 0;
	my $preproBInteractions = 0;
	while (<IN>) {
		$count++;
		chomp;
		s/\r//g;
		my $og = $_;
		my @line = split /\t/;

		my $name = $line[0];
		next if ($addFlag == 1 && !exists($ints->{$name}));
		my $name1 = $line[1];
		my $chr1 = $line[2];
		my $start1 = $line[3];
		my $end1 = $line[4];
		my $strand1 = $line[5];
		my $v1 = $line[6];
		my $name2 = $line[7];
		my $chr2 = $line[8];
		my $start2 = $line[9];
		my $end2 = $line[10];
		my $strand2 = $line[11];
		my $v2 = 0;
		$v2 = $line[12] if (@line > 12);
		my $v = 1;
		$v = $line[14] if (@line > 14);
		my $expected = 1;
		$expected = $line[15] if (@line > 15);

		if ($count == 1) {
			$numLinesInteractionFile = scalar(@line);
			if (checkInt($start1) || checkInt($end1)
					|| checkInt($start2) || checkInt($end2)) {
				$interactionHeader = $og;
				next;
			} else {
				$interactionHeader = '';
			}
		}
		if ($bgChr ne '') {
			if ($bgChr ne $chr1 || $bgChr ne $chr2
						|| $end1 < $bgStart || $end2 < $bgStart
						|| $start1 > $bgEnd || $start2 > $bgEnd) {
				$locationFiltered++;
				next;
			}
		}

		if (checkInt($start1) || checkInt($end1)
					|| checkInt($start2) || checkInt($end2)) {
			#print STDERR "$start1\t$end1\t$end2\t$start2\n";
			next;
		}
		$total += ($end2-$start2) + ($end1-$start1);
		$N+=2;

		$strand1 = '+' if ($strand1 eq '0');
		$strand1 = '-' if ($strand1 eq '1');
		$strand2 = '+' if ($strand2 eq '0');
		$strand2 = '-' if ($strand2 eq '1');
		my $m1 = floor(($start1+$end1)/2);
		my $m2 = floor(($start2+$end2)/2);
		my $dist = 'interchromosomal';
		if ($chr1 eq $chr2) {
			$dist = abs($m1-$m2);
			next if ($dist < $minDist);
			next if ($maxDist > 0 && $dist > $maxDist);
		} else {
			next if ($maxDist > 0);
		}
		

		if (exists($interactions{$name})) {
			$uniq{$name}++;
			$name = $name . "-" . $uniq{$name};
		}
		my %ann1 = ();
		my %ann2 = ();
		#my $zscore = 'NA';
		my $zscore = 0;
		$zscore = $line[16] if (@line > 16);
		#my $logp = 'NA';
		my $logp = 0;
		$logp = $line[17] if (@line > 17);
		my $dlogp = 'NA';
		$dlogp = $line[25] if (@line > 25);
		my $dzscore = 'NA';
		$dzscore = $line[26] if (@line > 26);
		my $thickness = $line[@line-1];

		
		if ($rmpreproB == 1) {
			if (($chr1 eq 'chr7' && $chr2 eq 'chr12')
					|| ($chr2 eq 'chr7' && $chr1 eq 'chr12')) {
				$preproBInteractions++;
				next;
			}
		}
		
		if ($addFlag == 1) {
			my $statline = {v1=>$v1,v2=>$v2,v=>$v,z=>$zscore,lp=>$logp,dz=>$dzscore,
						ex=>$expected,dlp=>$dlogp,th=>$thickness};
			$ints->{$name}->{'hic'}->{$hicDir} = $statline;

		} else {

			next if ($logp ne 'NA' && ($logp > log($pvalueCutoff)));
			next if ($zscore ne 'NA' && ($zscore < $zscoreCutoff));
			next if ($dlogp ne 'NA' && ($dlogp > log($dpvalueCutoff)));
			next if ($dzscore ne 'NA' && ($dzscore < $dzscoreCutoff));

			my %hicExps = ();
			$interactions{$name} = {c1=>$chr1,s1=>$start1,e1=>$end1,d1=>$strand1,v1=>$v1,m1=>$m1,
						c2=>$chr2,s2=>$start2,e2=>$end2,d2=>$strand2,v2=>$v2,m2=>$m2,
						dist=>$dist,v=>$v,ex=>$expected,a1=>\%ann1,a2=>\%ann2,z=>$zscore,lp=>$logp,og=>$og,
						dz=>$dzscore,dlp=>$dlogp,th=>$thickness,hic=>\%hicExps};
			$totalInteractions++;
		}
	}
	if ($addFlag == 1) {
		return;
	}
	$res = $total/$N if ($N > 0);
	if ($rmpreproB == 1) {
		print STDERR "\tpreproB interactions: $preproBInteractions\n";
	}
	if ($bgChr ne '') {
		print STDERR "\t$locationFiltered interactions filtered for not being within $bgChr:$bgStart-$bgEnd\n";
	}
	print STDERR "\tFound $totalInteractions (estimated resolution: $res)\n";
	return \%interactions;
}
sub getLengthDistribution {
	my ($interactions) = @_;

	my %bins = ();
	my $inter = 0;
	my $max = 0;
	foreach(keys %$interactions) {
		my $d  = $interactions->{$_}->{'dist'};
		if ($d eq 'interchromosomal') {
			$inter++;
		} else {
			my $index = floor(($d/$res+0.5));
			$bins{$index}++;
			$max = $index if ($max < $index);
		}
	}	
	if ($outputDirectory ne '') {
		open OUT, ">$outputDirectory/lengthDist.txt";
		print OUT "Distance Between Interactions\tCount\n";	
		print OUT "Intchr\t$inter\n";
		for (my $i=0;$i<=$max;$i++) {
			my $v = $res*$i;
			my $v2 = 0;
			$v2 = $bins{$i} if (exists($bins{$i}));
			print OUT "$v\t$v2\n";
		}
		close OUT;
	}
	
}
sub checkInt {
	my ($num) = @_;
	if ($num =~ /^[-+\d\.]+$/) {
		return 0;
	} else {
		return 1;
	}
}

sub printHUBs {
	my ($peakIntFile, $hubCount) = @_;

	my %hubDist = ();
	
	open OUT, ">$outputDirectory/hubs.gt$hubCount.interactons.txt";
	print OUT "#Hubs >= $hubCount interactions\n";
	my $totalHubs = 0;
	my $totalPeaks = 0;
	open IN, $peakIntFile;
	while (<IN>) {
		my $og = $_;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($line[0] =~ /^#/) {
			print OUT $og;
			next;
		}
		$totalPeaks++;
		if (@line > 7) {
			$hubDist{$line[7]}++;
			if ($line[7] >= $hubCount) {
				print OUT $og;
				$totalHubs++;
				next;
			}
		}
	}
	close IN;
	close OUT;
	print STDERR "\tHubs (>=$hubCount interactions): $totalHubs/$totalPeaks\n";
	open OUT, ">$outputDirectory/hubs.distribution.txt";
	print OUT "Number of Interactions per Peak\tNumber of Peaks\n";	
	my @bins = sort {$a <=> $b} keys %hubDist;
	foreach(@bins) {
		print OUT "$_\t$hubDist{$_}\n";
	}
	close OUT;

}
sub filterPeakInteractions {
	my ($interactions, $filterPeakFile1,$filterPeakFile2) = @_;

	my %score = ();
	foreach(keys %$interactions) {
		$score{$_}={L1=>"",L2=>"",R1=>"",R2=>""};
	}
	`chopify.pl "$filterPeakFile1" $res > "$tmpFile3" 2> /dev/null`;
	`mergePeaks -d $res -cobound 1 "$leftPeakFile" "$tmpFile3" -prefix "$tmpFile2" 2> /dev/null`;
	#print STDERR "`chopify.pl $filterPeakFile1 $res > $tmpFile3 2> /dev/null`\n";
	#print STDERR "`mergePeaks -d $res -cobound 1 $leftPeakFile $tmpFile3 -prefix $tmpFile2 2> /dev/null`\n";
	open IN, "$tmpFile2.coBoundBy1.txt";
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		$score{$line[0]}->{'L1'} = 1;
	}
	close IN;
	`mergePeaks -d $res -cobound 1 "$rightPeakFile" "$tmpFile3" -prefix "$tmpFile2" 2> /dev/null`;
	open IN, "$tmpFile2.coBoundBy1.txt";
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		$score{$line[0]}->{'R1'} = 1;
		$present{$line[0]}++;
	}
	close IN;
	`rm -f "$tmpFile2.coBoundBy0.txt" "$tmpFile2.coBoundBy1.txt" "$tmpFile3"`;
	if ($filterPeakFile2 eq '') {
		foreach(keys %$interactions) {
			$score{$_}->{'L2'} = 1;
			$score{$_}->{'R2'} = 1;
		}
	} else {
		`chopify.pl "$filterPeakFile2" $res > "$tmpFile3" 2> /dev/null`;
		`mergePeaks -d $res -cobound 1 "$leftPeakFile" "$tmpFile3" -prefix "$tmpFile2" 2> /dev/null`;
		open IN, "$tmpFile2.coBoundBy1.txt";
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			$score{$line[0]}->{'L2'} = 1;
		}
		close IN;
		`mergePeaks -d $res -cobound 1 "$rightPeakFile" "$tmpFile3" -prefix "$tmpFile2" 2> /dev/null`;
		open IN, "$tmpFile2.coBoundBy1.txt";
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			$score{$line[0]}->{'R2'} = 1;
			$present{$line[0]}++;
		}
		close IN;
		`rm -f "$tmpFile2.coBoundBy0.txt" "$tmpFile2.coBoundBy1.txt" "$tmpFile3"`;
	}
	my %newInteractions = ();
	my $totalInters = 0;
	my $keptInters = 0;
	foreach (keys %$interactions) {
		$totalInters++;
		if ($score{$_}->{'L1'} eq '1' && $score{$_}->{'R2'} eq '1') {
			$newInteractions{$_} = $interactions->{$_};
			$keptInters++;
		} elsif ($score{$_}->{'L2'} eq '1' && $score{$_}->{'R1'} eq '1') {
			$newInteractions{$_} = $interactions->{$_};
			$keptInters++;
		}
	}
	my $n = $filterPeakFile2;
	$n = "genome" if ($n eq '');
	print STDERR "\tKept $keptInters/$totalInters that connect from $filterPeakFile1 to $n\n";
	$interactions = \%newInteractions;
}
sub filterInteractions {
	my ($interactions, $bgPeakFile) = @_;

	my %newInteractions = ();
	`chopify.pl "$bgPeakFile" $res > "$bgChoppedFile" 2> /dev/null`;
	`mergePeaks -d $res -cobound 1 "$leftPeakFile" "$bgChoppedFile" -prefix "$tmpFile2" 2> /dev/null`;
	my %present = ();
	open IN, "$tmpFile2.coBoundBy1.txt";
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		$present{$line[0]} = 1;
	}
	close IN;
	`mergePeaks -d $res -cobound 1 "$rightPeakFile" "$bgChoppedFile" -prefix "$tmpFile2" 2> /dev/null`;
	open IN, "$tmpFile2.coBoundBy1.txt";
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		$present{$line[0]}++;
	}
	close IN;
	`rm -f "$tmpFile2.coBoundBy0.txt" "$tmpFile2.coBoundBy1.txt"`;
	my $totalInters = 0;
	my $keptInters = 0;
	foreach (keys %$interactions) {
		if (exists($present{$_})) {
			if ($present{$_} == 2) {
				$newInteractions{$_} = $interactions->{$_};
			}
		}
	}
	print STDERR "\tKept $keptInters/$totalInters that overlapped with peaks in $bgPeakFile\n";
	$interactions = \%newInteractions;
}
sub getIntersectingInteractions {
	my ($inters1, $inters2, $tolerance, $outputFile) = @_;

	my @set1 = ();
	my @set2 = ();
	foreach(keys %$inters1) {
		my $id = $_;
		my $chr1 = $inters1->{$_}->{'c1'};
		my $p1 = ($inters1->{$_}->{'s1'}+$inters1->{$_}->{'e1'})/2;
		my $chr2 = $inters1->{$_}->{'c2'};
		my $p2 = ($inters1->{$_}->{'s2'}+$inters1->{$_}->{'e2'})/2;
		my $switch = 0;
		if (($chr1 cmp $chr2) > 0) {
			$switch = 1;
		} elsif ($chr1 eq $chr2) {
			if ($p1 > $p2) {
				$switch=1;
			}
		}
		if ($switch) {
			my $tmp = $chr1;
			$chr1 = $chr2;
			$chr2 = $tmp;
			$tmp = $p1;
			$p1 = $p2;
			$p2 = $tmp;
		}
		my $d = {id=>$id,c1=>$chr1,p1=>$p1,c2=>$chr2,p2=>$p2};
		push(@set1, $d);
	}
	foreach(keys %$inters2) {
		my $id = $_;
		my $chr1 = $inters2->{$_}->{'c1'};
		my $p1 = ($inters2->{$_}->{'s1'}+$inters2->{$_}->{'e1'})/2;
		my $chr2 = $inters2->{$_}->{'c2'};
		my $p2 = ($inters2->{$_}->{'s2'}+$inters2->{$_}->{'e2'})/2;
		my $switch = 0;
		if (($chr1 cmp $chr2) > 0) {
			$switch = 1;
		} elsif ($chr1 eq $chr2) {
			if ($p1 > $p2) {
				$switch=1;
			}
		}
		if ($switch) {
			my $tmp = $chr1;
			$chr1 = $chr2;
			$chr2 = $tmp;
			$tmp = $p1;
			$p1 = $p2;
			$p2 = $tmp;
		}
		my $d = {id=>$id,c1=>$chr1,p1=>$p1,c2=>$chr2,p2=>$p2};
		push(@set2, $d);
	}

	#iterate through lists and find common interactions
	@set1 = sort {$a->{'c1'} cmp $b->{'c1'} || $a->{'p1'} <=> $b->{'p1'}} @set1;
	@set2 = sort {$a->{'c1'} cmp $b->{'c1'} || $a->{'p1'} <=> $b->{'p1'}} @set2;

	my $i2 = 0;
	my %commonInter1 = ();
	my %commonInter2 = ();
	for (my $i=0;$i<@set1;$i++) {
		my $s1c1 = $set1[$i]->{'c1'};
		my $s1p1 = $set1[$i]->{'p1'};
		my $s1c2 = $set1[$i]->{'c2'};
		my $s1p2 = $set1[$i]->{'p2'};
		#print STDERR "$i: $s1c1 $s1p1 $s1c2 $s1p2 | $i2: $set2[$i2]->{'c1'} $set2[$i2]->{'p1'} $tolerance\n";
		while ($i2 < @set2 && 
					(  (($set2[$i2]->{'c1'} cmp $s1c1) < 0)
						|| (($set2[$i2]->{'c1'} eq $s1c1) && $set2[$i2]->{'p1'} < ($s1p1 - $tolerance)) ) ) {
			$i2++;
		}
		last if ($i2 == @set2);
	
		for (my $j=$i2;$j<@set2;$j++) {
			my $diff1 = abs($set2[$j]->{'p1'} - $s1p1);
			last if ($set2[$j]->{'c1'} ne $s1c1 || $diff1 > $tolerance);

			my $diff2 = abs($set2[$j]->{'p2'} - $s1p2);
			if ($set2[$j]->{'c2'} eq $s1c2 && $diff2 < $tolerance) {
				# HA! found one!
				my $id1 = $set1[$i]->{'id'};
				my $id2 = $set2[$j]->{'id'};
				$commonInter1{$id1} = $inters1->{$id1};
				$commonInter2{$id2} = $inters2->{$id2};
			}
				
		}	
	}
	my $total1 = scalar(@set1);
	my $total2 = scalar(@set2);
	my $common1 = scalar(keys %commonInter1);
	my $common2 = scalar(keys %commonInter2);

	if ($outputFile ne '') {
		printInteractionAnnotation($interactions, "");
		#print STDERR "\tCommon Interactions (saved to $outputFile)\n";
		print STDERR "\t\tPrimary set: $common1 of $total1\n";
		print STDERR "\t\tSecndry set: $common2 of $total2\n";
	}
	return \%commonInter1;
}

sub printCircos {

	my ($inter) = @_;
	foreach(keys %$inter) {
		my $thickness = $inter->{$_}->{'th'};
		#$thickness = floor(abs($inter->{$_}->{'lp'}/8))+1;
		#$thickness = sprintf("%.0f",(-1*$inter->{$_}->{'lp'})/10+1);
		print "$_ $inter->{$_}->{'c1'} $inter->{$_}->{'s1'} $inter->{$_}->{'e1'} thickness=$thickness\n";
		print "$_ $inter->{$_}->{'c2'} $inter->{$_}->{'s2'} $inter->{$_}->{'e2'} thickness=$thickness\n";
	}
}
sub connectPeaks {
	my ($interactions, $peakFile1, $peakFile2) = @_;
	my %mapping = ();
	my $sameFlag = 0;
	if ($peakFile1 eq $peakFile2) {
		print STDERR "\tConnecting peaks to themselves\n";
		$sameFlag = 1;
	}
	foreach(keys %$interactions) {
		my %a = (); my %b = (); my %c = (); my %d = ();
		$mapping{$_} = {left1=>\%a,right1=>\%b,left2=>\%c,right2=>\%d};
	}
	`annotateRelativePosition.pl "$leftPeakFile," "$peakFile1," > "$tmpFile2"`;
	open IN, $tmpFile2;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		if (abs($line[2]) < $halfRes) {
			$mapping{$line[0]}->{'left1'}->{$line[1]}=$line[2];
		}
	}
	close IN;
	`annotateRelativePosition.pl "$leftPeakFile," "$peakFile2," > "$tmpFile2"`;
	open IN, $tmpFile2;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		if (abs($line[2]) < $halfRes) {
			$mapping{$line[0]}->{'left2'}->{$line[1]}=$line[2];
		}
	}
	close IN;
	`annotateRelativePosition.pl "$rightPeakFile," "$peakFile1," > "$tmpFile2"`;
	open IN, $tmpFile2;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		if (abs($line[2]) < $halfRes) {
			$mapping{$line[0]}->{'right1'}->{$line[1]}=$line[2];
		}
	}
	close IN;
	`annotateRelativePosition.pl "$rightPeakFile," "$peakFile2," > "$tmpFile2"`;
	open IN, $tmpFile2;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		if (abs($line[2]) < $halfRes) {
			$mapping{$line[0]}->{'right2'}->{$line[1]}=$line[2];
		}
	}
	close IN;
	`rm "$tmpFile2"`;

	foreach(keys %mapping) {
		my $intID = $_;
		foreach(keys %{$mapping{$intID}->{'left1'}}) {
			my $p1ID = $_;
			foreach(keys %{$mapping{$intID}->{'right2'}}) {
				my $p2ID = $_;
				print "$p1ID\t$p2ID\t$interactions->{$intID}->{'dist'}\n";
			}
		}
		foreach(keys %{$mapping{$intID}->{'left2'}}) {
			my $p2ID = $_;
			foreach(keys %{$mapping{$intID}->{'right1'}}) {
				my $p1ID = $_;
				print "$p1ID\t$p2ID\t$interactions->{$intID}->{'dist'}\n";
			}
		}
	}
	

}
