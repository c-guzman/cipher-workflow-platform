#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";


# Copyright 2009, 2010, 2011, 2012 Christopher Benner <cbenner@ucsd.edu>
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

	print STDERR "\n\t*** For standard analysis of rna or repeats, analyzeRepeats.pl is faster and recommended\n";
	print STDERR "\n\tUsage: analyzeRNA.pl <rna | repeats> <genome version>  [additional options...]\n";
	print STDERR "\t -or-: analyzeRNA.pl <custom RNA/GTF file> <organism|none>  [additional options...]\n";
	print STDERR "\n\tProgram for quantifying RNA tag counts.  The first argument can be \"rna\" (refseq genes),\n";
	print STDERR "\t\"repeats\" (repeat classes), or a custom RNA definition file (GTF format).\n";
	print STDERR "\n\tAvailable Genomes (required argument): (name,org,directory,default promoter set)\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}"
					. "\t$config->{'GENOMES'}->{$_}->{'promoters'}\n";
	}
	print STDERR "\n\tPrimary Annotation Options:\n";
	print STDERR "\t\t-d <tag directory 1> [tag directory 2] ... (list of experiment directories to show\n";
	print STDERR "\t\t\ttag counts for) NOTE: -dfile <file> where file is a list of directories in first column\n";
	print STDERR "\t\t-rpkm (Report results as reads per kb per million mapped)\n";
	print STDERR "\t\t-norm <#> (Normalize to total mapped tags: default 1e7)\n";
	print STDERR "\t\t-normMatrix <#> (Normalize to total tags in gene expression matrix: not used)\n";
	print STDERR "\t\t-noadj (Don't normalize, -raw works too)\n";
	print STDERR "\t\t-count <exons|introns|genes|5utr|3utr|cds> (Count tags in introns, exons, etc., default: genes)\n";
	print STDERR "\t\t-noCondensing (do not condense counts from entries will same ID, default: do condense)\n";
	print STDERR "\t\t-pc <#> (maximum tags to count per position, default: 0=no limit)\n";
	print STDERR "\t\t-strand <+|-|both> (count tags on indicated strand, default: +)\n";
	print STDERR "\t\t-gene <data file> ... (Adds additional data to result based on the closest gene.\n";
	print STDERR "\t\t\tThis is useful for adding gene expression data.  The file must have a header,\n";
	print STDERR "\t\t\tand the first column must be a GeneID, Accession number, etc.  If the peak\n";
	print STDERR "\t\t\tcannot be mapped to data in the file then the entry will be left empty.\n";
	print STDERR "\t\t-log (output tag counts as randomized log2 values - for scatter plots)\n";
	print STDERR "\t\t-sqrt (output tag counts as randomized sqrt values - for scatter plots)\n";
	print STDERR "\t\t-tss (estimate actual TSS in 1st exon and report as the centered position in columns 3 & 4)\n";
	print STDERR "\t\t-start <#> (start counting tags relative # offset of beginning of gene)\n";
	print STDERR "\t\t-end <#> (finish counting tags relative # offset to end of the gene)\n";
	print STDERR "\t\t-maxLength <#> (Don't count tags past # bp from the TSS, useful for GroSeq)\n";
	#print STDERR "\tOutlier filtering:\n";
	#print STDERR "\t\t-quant <#low> <#high> (calculate expression using only quantiles between low and high)\n";
	#print STDERR "\t\t-quantLength <#> (length of bins for quantile calculation, default: 400)\n";
	print STDERR "\t\t-pausing <#> (calculate ratio of pausing first [# bp of transcript] to gene body)\n";
	print STDERR "\t\t\tProduces 3 columns - promoter rpk, body rpk, and ratio (add -log for log versions)\n";
	print STDERR "\t\t\tAlso sets \"-count genes\".  Use \"-strand both\" when analyzing Pol II ChIP-Seq\n";
	print STDERR "\t\t\trpk is reads per kb - set -norm 1e6 or -normMatrix 1e6 to get rpkm\n";
	#print STDERR "\t\t-hist <bin size in bp> (i.e 1, 2, 5, 10, 20, 50, 100 etc.)\n";

	print STDERR "\n";
	exit;
}

if (@ARGV < 2) { 
	printCMD();
}
$cmd = "analyzeRNA.pl";
for (my $i=0;$i<@ARGV;$i++) {
	$cmd .= " $ARGV[$i]";
}

$upstream = -1000;
my $downstream = 1000;

my %possibleCountMethods = ();
$possibleCountMethods{"exons"}=1;
$possibleCountMethods{"introns"}=1;
$possibleCountMethods{"genes"}=1;
$possibleCountMethods{"5utr"}=1;
$possibleCountMethods{"3utr"}=1;
$possibleCountMethods{"cds"}=1;
$possibleCountMethods{"all"}=1;

my @countMethods = ("genes");

my $logFlag = 0;
my $sqrtFlag = 0;


my $fragLength = 150;
my $revoppFlag = 1;
my $init2One = 0;
my $maxLength = 1e10;
my $maskFlag =0;


my $size = 300;

my $strand = "+";
my $pauseStart = "";
my $offsetStart = 0;
my $offsetEnd = 0;
my $annFlag = 1;
my $local = 0;
my $adjustFlag = 1;
my $normValue = 1e7;
my $normMatrix = 0;
my $tssFlag = 0;
my $geneListFile = '';
my $ugFile = '';
my @geneDataFiles = ();
my $histBinSize = 0;
my $diFlag = 0;
my $nucFlag = 0;
my $centerMotif = '';
my $cpromoter = '';
my $histNorm = 0;
my $ghistFlag = 0;
my $goDir = '';
my $pDistFlag = 0;
my $noAnnFlag = 0;
my $minHistogramThreshold = 20;
my $exoncalledYet = 0;
my $rpkmFlag = 0;
my $condenseFlag = 1;
my $repeatFlag = 0;
my $quantStart = -1;
my $quantEnd = -1;
my $quantLength = 400;

my $rnaDefFile = $ARGV[0];
my $genome = $ARGV[1];


my @motifFiles = ();
my @tagDirs = ();
my @peakFiles = ();
my %filterMotifFiles = ();

for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-d') {
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
	} elsif ($ARGV[$i] eq '-gene') {
		print STDERR "\tGene Data Files:\n";
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@geneDataFiles, $ARGV[$i]);
			print STDERR "\t\t$ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-hist') {
		$histBinSize = $ARGV[++$i];
		print STDERR "\t-----------------------------------------------------\n";
		print STDERR "\tHistogram mode activated (bin size = $histBinSize bp)\n";
		print STDERR "\t-----------------------------------------------------\n";
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
	} elsif ($ARGV[$i] eq '-count') {
		$i++;
		if ($exoncalledYet == 0) {
			@countMethods = ();
			$exoncalledYet = 1;
		}
		if (!exists($possibleCountMethods{$ARGV[$i]})) {
			print STDERR "!!! \"-count $ARGV[$i]\" is not a valid option!\n";
			exit;
		}
		if ($ARGV[$i] eq 'all') {
			push(@countMethods, "exons", "introns", "genes");
		} else {
			push(@countMethods, $ARGV[$i]);
		}
	} elsif ($ARGV[$i] eq '-pausing') {
		print STDERR "\tWill output ratio of promoter/gene body tag counts\n";
		@countMethods = ('genes');
		$pauseStart = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pc') {
		$init2One = $ARGV[++$i];
		print STDERR "\tMaximum count per bp will be set to $init2One\n";
	} elsif ($ARGV[$i] eq '-log') {
		$logFlag = 1;
		print STDERR "\tWill output log(1+rand()+x) for tag counts\n";
	} elsif ($ARGV[$i] eq '-sqrt') {
		$sqrtFlag = 1;
		print STDERR "\tWill output sqrt(rand()+x) for tag counts\n";
	} elsif ($ARGV[$i] eq '-quantLength') {
		$quantLength = $ARGV[++$i];
		print STDERR "\tWill calculate quantile bins using a size of $quantLength\n";
		if ($quantStart == -1 && $quantEnd == -1) {
			$quantStart = 0.2;
			$quantEnd = 0.8;
			print STDERR "\tBy default using quantiles between $quantStart and $quantEnd for expression\n";
		}
	} elsif ($ARGV[$i] eq '-tss') {
		$tssFlag = 1;
		print STDERR "\tWill estimate TSS position within first exon (reported as midpoint of columns 3 & 4)\n";
	} elsif ($ARGV[$i] eq '-maxLength') {
		$maxLength = $ARGV[++$i];
		print STDERR "\tWill count reads within $maxLength bp from the TSS\n";
	} elsif ($ARGV[$i] eq '-quant') {
		$quantStart = $ARGV[++$i];
		$quantEnd = $ARGV[++$i];
		print STDERR "\tWill only use quantiles between $quantStart and $quantEnd to calculated expression\n";
	} elsif ($ARGV[$i] eq '-noCondensing') {
		$condenseFlag = 0;
		print STDERR "\tWill NOT condense counts from genes will the same accession numbers\n";
	} elsif ($ARGV[$i] eq '-ghist') {
		$ghistFlag = 1;
		print STDERR "\tWill create histogram for each gene\n";
	} elsif ($ARGV[$i] eq '-noadj' || $ARGV[$i] eq '-raw') {
		$adjustFlag = 0;
		print STDERR "\tWill NOT normalize tag counts\n";
	} elsif ($ARGV[$i] eq '-list') {
		$geneListFile = $ARGV[++$i];
		print STDERR "\tAnalysis will be restricted to genes found in file: $geneListFile\n";
	} elsif ($ARGV[$i] eq '-rpkm') {
		$normValue = 1e6;
		$rpkmFlag = 1;
		print STDERR "\tReporting tag counts as RPKM (reads per kb per million reads mapped)\n";
	} elsif ($ARGV[$i] eq '-norm') {
		$normValue = $ARGV[++$i];
		print STDERR "\tWill normalize tag counts to $normValue per experiment\n";
	} elsif ($ARGV[$i] eq '-normMatrix') {
		$normValue = $ARGV[++$i];
		$normMatrix=1;
		print STDERR "\tWill normalize gene expression martrix tag counts to $normValue per experiment\n";
	} elsif ($ARGV[$i] eq '-strand') {
		$strand = $ARGV[++$i];
		print STDERR "\tWill counts tags on $strand strand(s)\n";
	} elsif ($ARGV[$i] eq '-len') {
		$fragLength = $ARGV[++$i];
		print STDERR "\tFragment Length set to $fragLength\n";
	} elsif ($ARGV[$i] eq '-start') {
		$offsetStart = $ARGV[++$i];
		if ($offsetStart < $upstream) {
			$upstream = $offsetStart;
		}
	} elsif ($ARGV[$i] eq '-end') {
		$offsetEnd = $ARGV[++$i];
		if ($offsetEnd > $downstream) {
			$downstream = $offsetEnd;
		}
	} elsif ($ARGV[$i] eq '-histNorm') {
		$histNorm = $ARGV[++$i];
		print STDERR "\tWill normalize Tag histograms with minimum total of $histNorm\n";
	} else {
		print STDERR "$ARGV[$i] not recognized\n\n";
		printCMD();
	
	}
}

if (@tagDirs < 1) {
	print STDERR "!!! No tag directories selected - need to count something...\n";
	exit;
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
	if ($rnaDefFile eq 'rna' || $rnaDefFile eq 'repeats') {
		print STDERR "!!! $rnaDefFile can only be used with a valid genome !!!\n";
		exit;
	}
} else {
	$genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'};	
	$organism = $config->{'GENOMES'}->{$genome}->{'org'};	
	if ($rnaDefFile eq 'rna') {
		$rnaDefFile = $genomeDir . $genome . ".rna";
	} elsif ($rnaDefFile eq 'repeats') {
		$repeatFlag = 1;
		$rnaDefFile = $genomeDir . $genome . ".repeats";
		$annFlag = 0;
	}
}

print STDERR "\n\tRNA Definition File = $rnaDefFile\n";
print STDERR "\tGenome = $genome\n";
print STDERR "\tOrganism = $organism\n";

#tmp files
my $rand = rand();
my $tmpfile = $rand . ".tmp";
my $tmpfile2 = $rand . ".tmp";
my $tmpRnaDefFile = $rand . ".rna";
my $tmpRnaDefGTFFile = $rand . ".rna";
my %toDelete = ();


#check RNA def file
my $gtfFlag = 0;
open IN, $rnaDefFile;
my $Z = 0;
while (<IN> ){
	$Z++;
	chomp;
	s/\r//g;
	next if (/^\s*\#/);
	my @line = split /\t/;
	if (@line > 8) {
		if ($line[8] =~ /gene_id/) {
			$gtfFlag = 1;
			last;
		}
	}
	if ($Z > 100) {
		last;
	}
}
close IN;

if ($gtfFlag) {
	print STDERR "\tParsing GTF file...\n";
	`parseGTF.pl "$rnaDefFile" rna > "$tmpRnaDefGTFFile"`;
	$rnaDefFile = $tmpRnaDefGTFFile;
}



my %geneList = ();
if ($geneListFile ne '') {
	my $idType = "refseq";
	`convertIDs.pl "$geneListFile" $organism $idType | cut -f1 > "$tmpfile"`;
	`cat "$tmpfile" "$geneListFile" | cut -f1 | sort | uniq > "$tmpfile2"`;
	`mv "$tmpfile2" "$tmpfile"`;
	`mergeData.pl "$tmpfile" "$rnaDefFile" > "$tmpRnaDefFile"`;
	$rnaDefFile = $tmpRnaDefFile;
	`rm "$tmpfile"`;
	$toDelete{$tmpRnaDefFile}=1;
} 

my $rnaDef = readRNADefinition($rnaDefFile,$repeatFlag);

my $percentResolution = $histBinSize;


my @exp = ();
my @norms = ();
my @dists = ();
my @columnNames = ();
for (my $i=0;$i<@tagDirs;$i++) {

	my $directory = $tagDirs[$i];
	#`getRelativeTagPositions.pl "$rnaDefFile",0 "$directory" size,$upstream,$downstream > $tmpfile`;
	`getPeakTags "$rnaDefFile" "$directory" -peaktags -fixed -start $upstream -end $downstream > $tmpfile`;
	my $rnaTags = readTagPositions($tmpfile,$strand);
	`rm $tmpfile`;

	my ($totalDirTags, $totalDirPos, $fragEst, $peakEst) = HomerConfig::readTagInfo($directory,$init2One);

	if ($histBinSize > 0) {
		my ($p5dist, $p3dist) = geneHistogram($rnaDef,$rnaTags, 30);
		push(@dists, $p5dist, $p3dist);
		push(@columnNames, "5p", "3p");
	} else {

		foreach(@countMethods) {
			my $countMethod = $_;

			my $expression = getGeneExpression($rnaDef,$rnaTags, $countMethod,$quantStart,$quantEnd,
												$quantLength,$offsetStart,$offsetEnd,$strand,$maxLength);
			if ($condenseFlag == 1) {
				$expression = condenseCopies($expression);
			}

			my $totalRNA = 1;
			if ($normMatrix) {
				foreach(values %$expression) {
					$totalRNA += $_;
				}
			} else {
				$totalRNA = $totalDirTags;
			}
			$totalRNA = 1 if ($totalRNA < 1);

			my $normRatio = $normValue/$totalRNA;
			$normRatio = 1 if ($adjustFlag == 0);
			
			if ($pauseStart eq '') {
				my $colName = $directory . " " . $countMethod;
				my $normRatioStr = sprintf("%.2f",$normRatio);
				$colName .= " (Total: $totalDirTags)";
				if ($rpkmFlag) {
					$colName .= " RPKM";
					$normRatio .= "L$countMethod";
				} else {
					$colName .= " normFactor " . $normRatioStr;
				}
				if ($sqrtFlag) {
					$colName .= " (Sqrt)";
				}
				if ($logFlag) {
					$colName .= " (log2)";
				}

				push(@norms, $normRatio);
				push(@exp, $expression);
				push(@columnNames, $colName);
			} else {
				my $minimumReadDensity = $totalRNA/4e9;
				my %lengths = ();
				my %dups = ();
				foreach(keys %$rnaDef) {
					my $Lg = $rnaDef->{$_}->{'Lgenes'};
					my $id = $_;
					$id =~ s/\-HOMER\d+// if ($condenseFlag);
					$lengths{$id} = $Lg;
					$dups{$id}++;
				}
				my $bodyExpression = getGeneExpression($rnaDef,$rnaTags, $countMethod,$quantStart,$quantEnd,
												$quantLength,$pauseStart,$offsetEnd,$strand,$maxLength);
				if ($condenseFlag == 1) {
					$bodyExpression = condenseCopies($bodyExpression);
				}
				my %ratios = ();
				foreach(keys %$expression) {
					my $id = $_;
					my $promExp = $expression->{$id}-$bodyExpression->{$id};
					my $pLength = $pauseStart-$offsetStart;
					my $gLength = floor($lengths{$id}/$dups{$id});
					$gLength -= $offsetStart;
					$gLength += $offsetEnd;
					my $bLength = $gLength - $pLength;

					my $pRPK = 0;
					$pLength = 1 if ($pLength < 1);
					$pRPK = $promExp/$pLength*1000;

					my $bRPK = 0;
					$bLength = 1 if ($bLength < 1);
					$bRPK = $bodyExpression->{$id}/$bLength*1000;

					my $unitPDensity = 1/($pLength)*1000;
					my $unitBDensity = 1/($bLength)*1000;

					my $den = $bRPK;
					$den = $unitBDensity if ($bRPK < $unitBDensity);
					my $num = $pRPK;
					$num = $unitBDensity if ($pRPK < $unitBDensity);
					#print STDERR "BEFORE $num $den\n";
					$ratios{$id} = $num/$den;
					#print STDERR "$ratios{$id} $num $den\n";
					$expression->{$id} = $pRPK*$normRatio;
					$bodyExpression->{$id} = $bRPK*$normRatio;
					if ($logFlag) {
						$ratios{$id} = log($ratios{$id})/log(2);	
						my $rd = $expression->{$id};
						$unitPDensity*=$normRatio;
						if ($rd < $unitPDensity) {
							$rd = log($rd/2+rand()*$unitPDensity/2)/log(2);
						} else {
							$rd = log($rd+rand()*$unitPDensity)/log(2);
						}
						$expression->{$id} = $rd;
						$rd = $bodyExpression->{$id};
						$unitBDensity*=$normRatio;
						if ($rd < $unitBDensity) {
							$rd = log($rd/2+rand()*$unitBDensity/2)/log(2);
						} else {
							$rd = log($rd+rand()*$unitBDensity)/log(2);
						}
						$bodyExpression->{$id} = $rd;
					}
				}

				my $colName1 = $directory . " Pause Ratio";
				my $colName2 = $directory . " Promoter RPK";
				my $colName3 = $directory . " Gene Body RPK";
				my $normRatioStr = sprintf("%.2f",$normRatio);
				my $colName = " (Total: $totalDirTags)";
				$colName .= " normFactor " . $normRatioStr;
				if ($sqrtFlag) {
					$colName .= " (Sqrt)";
				}
				if ($logFlag) {
					$colName .= " (log2)";
				}
				$colName1 .= $colName;
				$colName2 .= $colName;
				$colName3 .= $colName;

				push(@norms, 1.0);
				if ($logFlag) {
					push(@norms, 1.0);
					push(@norms, 1.0);
				} else {
					push(@norms, $normRatio);
					push(@norms, $normRatio);
				}
				push(@exp, \%ratios);
				push(@exp, $expression);
				push(@exp, $bodyExpression);
				push(@columnNames, $colName1);
				push(@columnNames, $colName2);
				push(@columnNames, $colName3);
			}
		}
	}
}

if ($histBinSize > 0) {
	printHistogram(\@dists, \@columnNames);
} else {
	addGeneAnnotation($rnaDef);
	$logFlag = 0 if ($pauseStart ne '');
	$sqrtFlag = 0 if ($pauseStart ne '');
	printGeneExpMatrix($rnaDef, \@exp, \@columnNames, \@norms,$rpkmFlag, $annFlag,$tssFlag);
}

foreach(keys %toDelete) {
	`rm \"$_\"`;
}

if ($rnaDefFile eq $tmpRnaDefGTFFile) {
	`rm -f "$tmpRnaDefGTFFile"`;
}

print STDERR "\n";

exit;

sub printHistogram {
	my ($dists, $columnNames) = @_;

	print "Relative Position";
	foreach(@$columnNames) {
		print "\t$_";
	}
	print "\n";
		
	for (my $i=0;$i<=$percentResolution;$i++) {
		my $v = $i/$percentResolution;
		print "$v";
		foreach(@$dists) {
			if (exists($_->{$i})) {
				print "\t$_->{$i}";
			} else {
				print "\t0";
			}
		}
		print "\n";
	}
}

sub geneHistogram {
	my ($rnaDef, $rnaTags, $len) = @_;


	my %total5p = ();
	my %total3p = ();
	my %totalCoverage = ();

	my %exp = ();
	foreach(keys %$rnaDef) {
		my $id = $_;
		next if (!exists($rnaTags->{$id}));

		my $dir = $rnaDef->{$id}->{'d'};
		if ($dir eq '-' || $dir eq 1) {
			$dir = 1;
		} else {
			$dir = 0;
		}

		my %p5tags = ();
		my %p3tags = ();
		my %coverageDiff = ();
		my %coverage = ();

		my @defPos = sort {$a <=> $b} keys %{$rnaDef->{$id}->{'v'}};
		my @tagPos = sort {$a <=> $b} keys %{$rnaTags->{$id}};
		if (@defPos < 1) {
			print STDERR "$id - no definition\n";
		}


		my %histInfo = ();
		my $totalLength = 0;
		for (my $i=0;$i<@defPos-1;$i++) {
			my $p = $defPos[$i];
			my $def = $rnaDef->{$id}->{'v'}->{$p};
			if ($def =~ /^E/) {
				$histInfo{$def}=$totalLength;
				my $end = $defPos[$i+1];
				my $len = $end-$p;
				$totalLength+=$len;
			}
		}


		my $defIndex = 0;
		my $defPosStart = $defPos[$defIndex];
		my $defPosEnd = 1e12;
		my $def = $rnaDef->{$id}->{'v'}->{$defPos[$defIndex]};
		if (@defPos > $defIndex+1) {
			$defPosEnd = $defPos[$defIndex+1];
		}

		my $total5Tags = 0;
		my $total3Tags = 0;

		for (my $i=0;$i<@tagPos;$i++) {
			if ($tagPos[$i] < $defPosStart) {
				next;
			}
			my $stop = 0;
			while ($tagPos[$i] > $defPosEnd) {
				$defIndex++;
				if ($defIndex >= @defPos) {
					$stop = 1;
					last;
				}	
				$defPosStart = $defPos[$defIndex];
				$defPosEnd = 1e12;
				$def = $rnaDef->{$id}->{'v'}->{$defPos[$defIndex]};
				if (@defPos > $defIndex+1) {
					$defPosEnd = $defPos[$defIndex+1];
				}
			}
			last if ($stop);

			if (exists($histInfo{$def})) {
				my $p = $histInfo{$def} + ($tagPos[$i]-$defPosStart);
				my $val = $rnaTags->{$id}->{$tagPos[$i]};
				if ($val ne 'NA') {
					$p5tags{$p} += $val;
					$coverageDiff{$p} += $val;
					$coverageDiff{$p+$len} -= $val;
					$total5Tags += $val;
				}
				$val = $rnaTags->{$id}->{$tagPos[$i]};
				if ($val ne 'NA') {
					$p3tags{$p} += $val;
					$coverageDiff{$p} -= $val;
					$coverageDiff{$p-$len} += $val;
					$total3Tags += $val;
				}

			}

		}

		# create gene specific histogram
		my %norm5p = ();
		my %norm3p = ();
		my %normCoverage = ();
		my @histPos = sort {$a <=> $b} keys %coverageDiff;

		if ($totalLength < 1) {
			$totalLength=1;
		}
	
		my $current = 0;
		foreach(@histPos) {
			my $p = $_;
			my $presult = floor($p/$totalLength*$percentResolution + 0.5);
			$current += $coverageDiff{$p};
			$normCoverage{$presult} = $current;
			$norm5p{$presult}=0;
			$norm3p{$presult}=0;
			$norm5p{$presult}=$p5tags{$p} if (exists($p5tags{$p}));
			$norm3p{$presult}=$p3tags{$p} if (exists($p3tags{$p}));
		}

		if ($total5Tags+$total3Tags > $minHistogramThreshold) {
			my @pos = sort {$a <=> $b} keys %normCoverage;
			foreach(@pos) {
				my $v = $normCoverage{$_}/($total5Tags+$total3Tags);
				$totalCoverage{$_} += $v;
			}
		}
		if ($total5Tags > $minHistogramThreshold) {
			my @pos = sort {$a <=> $b} keys %norm5p;
			foreach(@pos) {
				my $v = $norm5p{$_}/$total5Tags;
				$total5p{$_} += $v;
			}
		}
		if ($total3Tags > $minHistogramThreshold) {
			my @pos = sort {$a <=> $b} keys %norm3p;
			foreach(@pos) {
				my $v = $norm3p{$_}/$total3Tags;
				$total3p{$_} += $v;
			}
		}
	}
	return (\%totalCoverage,\%total5p,\%total3p);
}




sub addGeneAnnotation {
	my ($rnaDef) = @_;

	if ($organism eq 'unknown') {
		return;
	}

	my $convFile = $homeDir . "/data/accession/$organism" . "2gene.tsv";
	unless (-f $convFile) {
		return;
	}

	my %altIDs = ();
	foreach(keys %$rnaDef) {
		my $id = $_;
		my $altID =$id;
		$altID =~ s/\.\d+?$//;
		if (!exists($altIDs{$altID})) {
			my @a = ();
			$altIDs{$altID} = \@a;
		}
		push(@{$altIDs{$altID}}, $id);
	}

	open IN, $convFile;
	while (<IN>) {
		chomp;
		my @line= split /\t/;
		if (exists($rnaDef->{$line[0]})) {
			$rnaDef->{$line[0]}->{'gid'} = $line[1];
		} elsif (exists($altIDs{$line[0]})) {
			foreach(@{$altIDs{$line[0]}}) {
				$rnaDef->{$_}->{'gid'} = $line[1];
			}
		}
	}
	close IN;

	my $descriptionFile = $homeDir . "/data/accession/$organism.description";
	open IN, $descriptionFile;
	my %desc = ();
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line= split /\t/;
		my $gid = $line[0];
		$desc{$gid}->{'gid'} = $line[0];
		$desc{$gid}->{'ug'} = $line[1];
		$desc{$gid}->{'refseq'} = $line[2];
		$desc{$gid}->{'ensembl'} = $line[3];
		$desc{$gid}->{'name'} = $line[4];
		$desc{$gid}->{'alias'} = $line[5];
		$desc{$gid}->{'chr'} = $line[7];
		$desc{$gid}->{'desc'} = $line[8];
		$desc{$gid}->{'ttype'} = "NA";
		if (@line > 9) {
			$desc{$gid}->{'ttype'} = $line[8];
		}
	}
	close IN;

	foreach(keys %$rnaDef) {
		my $id = $_;
		next if (!exists($rnaDef->{$id}->{'gid'}));
		my $gid = $rnaDef->{$id}->{'gid'};
		next if (!exists($desc{$gid}));
		$rnaDef->{$id}->{'desc'}=$desc{$gid};
	}
}


sub printGeneExpMatrix {
	my ($rnaDef,$exp, $tagDirs, $norms,$rpkmFlag,$annFlag,$tssFlag) = @_;

	print "Gene ID (cmd=$cmd)\tchr\tstart\tend\tstrand\tmRNA Length(avg)\tGene Length(avg)\tCopies in Genome";
	if ($annFlag) {
		print "\tSymbol\tAlias\tDescription\ttype\tUnigene\tGeneID\tEnsembl";
	}
	foreach(@$tagDirs) {
		print "\t$_";
	}
	print "\n";
	my %done = ();

	my %counts = ();
	my %lengths = ();
	my %totalLengths = ();
	foreach(keys %$rnaDef) {
		my $Le = $rnaDef->{$_}->{'Lexons'};
		my $Lg = $rnaDef->{$_}->{'Lgenes'};
		my $id = $_;
		if ($condenseFlag) {
			$id =~ s/\-HOMER\d+//;
		}
		$counts{$id}++;
		$lengths{$id}+=$Le;
		$totalLengths{$id}+=$Lg;
	}

	foreach(keys %$rnaDef) {
		my $gene = $_;
		my $egene = $gene;
		if ($condenseFlag) {
			$egene =~ s/\-HOMER\d+//;
			next if (exists($done{$egene}));
			$done{$egene}=1;
		}

		my $c = $rnaDef->{$gene}->{'c'};
		my $s = $rnaDef->{$gene}->{'s'};
		my $e = $rnaDef->{$gene}->{'e'};
		my $d = $rnaDef->{$gene}->{'d'};

		if ($tssFlag) {
			my $est = estimateTSS($rnaDef->{$gene}->{'tssReads'});
			if ($est =~ /NA/) {
				$s = $est;
				$e = $est;
			} elsif ($d eq '+' || $d eq '0') {
				my $tssPos = $s + $est;
				$s = $tssPos-100;
				$e = $tssPos+100;
			} else {
				my $tssPos = $e - $est;
				$s = $tssPos-100;
				$e = $tssPos+100;
			}
		}

		my $count = $counts{$egene};
		my $len = floor($lengths{$egene}/$count);
		my $glen = floor($totalLengths{$egene}/$count);
		print "$egene\t$c\t$s\t$e\t$d\t$len\t$glen\t$count";
		if ($annFlag) {
			if (exists($rnaDef->{$egene}->{'desc'})) {
				my $x = $rnaDef->{$egene}->{'desc'};
				print "\t$x->{'name'}\t$x->{'alias'}\t$x->{'desc'}\t$x->{'ttype'}\t$x->{'ug'}\t$x->{'gid'}\t$x->{'ensembl'}";
			} else {
				print "\t\t\t\t\t\t\t";
			}
		}
		for (my $i=0;$i<@$exp;$i++) {
			my $v = 0;
			my $norm = 1;
			if (exists($exp->[$i]->{$egene})) {
				$v = $exp->[$i]->{$egene};
				if ($tagDirs->[$i] =~ /Pause Ratio/) {
				} else {
					if ($adjustFlag) {
						$norm = $norms->[$i];
						if ($norm =~ s/(L.*)$//) {
							my $code = $1;
							my $lenn = 10;
							if (exists($rnaDef->{$gene}->{$code})) {
								$lenn = $rnaDef->{$gene}->{$code};
								if ($lenn < 1) {
									$lenn = 1 
								}
							}
							$norm *= 1000/$lenn;
							$v *= $norm;
							$v = sqrt($v) if ($sqrtFlag);
							$v = log($v)/log(2) if ($logFlag);
						} else {
							$v *= $norm;
							if ($sqrtFlag) {
								$v = sqrt($v+$norm*rand());
							}
							if ($logFlag) {
								$v = log($v+1+$norm*rand())/log(2);
							}
						}
					} elsif ($logFlag) {
						$v = log($v+1+rand())/log(2);
					} elsif ($sqrtFlag) {
						$v = sqrt($v+rand());
					}
				} 
			}
			print "\t$v";
		}
		print "\n";
	}
}


sub condenseCopies {
	my ($exp) = @_;
	my %ids = ();
	foreach(keys %$exp) {
		my $g = $_;
		my $v = $exp->{$g};
		my $new = $g;
		$new =~ s/\-HOMER\d+$//;
		$ids{$new}+=$v;
	}
	return \%ids;
}

sub getGeneExpression {
	my ($rnaDef,$rnaTags,$method,$quantStart,$quantEnd,$quantLength,$offsetStart,
											$offsetEnd,$strand,$maxLength) = @_;

	my %exp = ();
	foreach(keys %$rnaDef) {
		my $id = $_;
		next if (!exists($rnaTags->{$id}));

		my $dir = $rnaDef->{$id}->{'d'};
		if ($dir eq '-' || $dir eq 1) {
			$dir = 1;
		} else {
			$dir = 0;
		}
		my $v = 0;

		my @defPos = sort {$a <=> $b} keys %{$rnaDef->{$id}->{'v'}};
		my @tagPos = sort {$a <=> $b} keys %{$rnaTags->{$id}};
		if (@defPos < 1) {
			print STDERR "$id - no definition\n";
		}
		my $defIndex = 0;
		my $defPosStart = $defPos[$defIndex];
		my $defPosAbsStart = 0;
		my $GENELENGTH = $rnaDef->{$id}->{'Lgenes'};
		my $ASDF = $GENELENGTH/$quantLength;
		my $defPosEnd = 1e12;
		my $def = $rnaDef->{$id}->{'v'}->{$defPos[$defIndex]};
		if (@defPos > $defIndex+1) {
			$defPosEnd = $defPos[$defIndex+1];
		}

		my @quantBins=();
		my $curQuantBin=0;
		my $curQuantLen=0;


		my $hardStart = 0;
		my $hardEnd = abs($rnaDef->{$id}->{'e'} - $rnaDef->{$id}->{'s'}+1);
		$hardStart += $offsetStart;
		$hardEnd += $offsetEnd;
		$hardEnd = $maxLength if ($hardEnd > $maxLength);

		for (my $i=0;$i<@tagPos;$i++) {
			#my $val = $rnaTags->{$id}->{$tagPos[$i]};
			#print STDERR "\t$i\t$tagPos[$i]\t$val\n";
				
			my $stop = 0;
			while ($tagPos[$i] >= $defPosEnd) {
				$defIndex++;
				if ($defIndex >= @defPos) {
					$stop = 1;
					last;
				}	
				$defPosStart = $defPos[$defIndex];
				$defPosEnd = 1e12;
				$def = $rnaDef->{$id}->{'v'}->{$defPos[$defIndex]};
				if (@defPos > $defIndex+1) {
					$defPosEnd = $defPos[$defIndex+1];
				}
			}
			last if ($stop);
			my $val = $rnaTags->{$id}->{$tagPos[$i]};

			if ($tssFlag) {
				if ($def =~ /E1[_:]/ || $def eq 'Upstream') {
					$rnaDef->{$id}->{'tssReads'}->{$tagPos[$i]} += $val;
				}
			}
			if ($tagPos[$i] < $defPosStart) {
				next;
			}
			next if ($tagPos[$i] < $hardStart);
			last if ($tagPos[$i] > $hardEnd);

			my $willAdd = 0;
			
			if ($method eq 'exons') {
				if ($def =~ /^E/ || $def =~ /^Upstream/ || $def =~ /^Downstream/) {
					if ($val ne 'NA') {
						$willAdd=1;
					}
				}
			} elsif ($method eq '5utr') {
				if ($def =~ /5UTR/ || $def =~ /^Upstream/) {
					if ($val ne 'NA') {
						$willAdd=1;
					}
				}
			} elsif ($method eq '3utr') {
				if ($def =~ /3UTR/ || $def =~ /^Downstream/) {
					if ($val ne 'NA') {
						$willAdd=1;
					}
				}
			} elsif ($method eq 'cds') {
				if ($def =~ /^E\d+$/) {
					if ($val ne 'NA') {
						$willAdd=1;
					}
				}
			} elsif ($method eq 'introns') {
				if ($def =~ /^I/) {
					if ($val ne 'NA') {
						$willAdd=1;
					}
				}
			} elsif ($method eq 'genes') {
				if ($def =~ /^I/ || $def =~ /^E/ || $def =~ /^Upstream/ || $def =~ /^Downstream/) {
					if ($val ne 'NA') {
						$willAdd=1;
					}
				}
			}
			$v += $val if ($willAdd);
			if ($willAdd && $quantEnd > 0) {
				my $curOffsetPos = abs($tagPos[$i]-$defPosAbsStart);
				my $quantIndex = floor($curOffsetPos/$quantLength);
				my $val = $rnaTags->{$id}->{$tagPos[$i]};
				$quantBins[$quantIndex]+=$val;
			}
			#if ($def =~ /^E/ || $def =~ /^I/) {
		}

		my $NB = @quantBins;
		if ($NB > 2 && $quantEnd > 0) {
			for (my $i=0;$i<$NB;$i++) {
				if (!defined($quantBins[$i])) {
					$quantBins[$i]=0;
				}
			}
			$v = 0;
			@quantBins = sort {$a <=> $b} @quantBins;
			for (my $i=0;$i<$NB;$i++) {
				my $frac = ($i+0.5)/$NB;
				#print STDERR "$i\t$frac\t$quantBins[$i]\n";
				if (defined($quantBins[$i]) && $frac >= $quantStart && $frac <= $quantEnd) {
					$v += $quantBins[$i];
				}
			}
		}
		#print STDERR "$id\t$NB ($ASDF)\t$GENELENGTH\t$v\n";
		$exp{$id} = $v;
	}
	return \%exp;
}

sub estimateTSS {

	my ($reads) = @_;
	my @pos = sort {$a <=> $b} keys %$reads;

	return 0 if (@pos < 1);

	my $fragLength = 50;
	my $EstimateSize = 100;

	for (my $i=@pos-1;$i>=0;$i--) {
		my $pi = $pos[$i];
		#print STDERR "$i\t$pi\t$reads->{$pi}\n";
		for (my $j=$i+1;$j<@pos;$j++) {
			my $pj = $pos[$j];
			last if ($pj - $pi > $fragLength);
			$reads->{$pj} += $reads->{$pi};
		}
	}
	my $lastPos = $pos[@pos-1];
	my $level = 0;
	my $levelN = 0;
	for (my $i=@pos-1;$i>=0;$i--) {
		last if ($pos[$i] < $lastPos-$EstimateSize);
		$level += $reads->{$pos[$i]};
		$levelN++;
	}
	if ($levelN > 0) {
		$level /= $levelN;
	}
	my $tssEst = 'NA';
	my $startLevel = 0;
	$levelN = 0;
	for (my $i=0;$i<@pos;$i++) {
		if ($pos[$i] < $upstream+$EstimateSize) {
			$startLevel+=$reads->{$pos[$i]};
			$levelN++;
		} elsif ($tssEst ne 'NA') {
			last;
		}
		if ($tssEst eq 'NA' && $reads->{$pos[$i]} > $level/2) {
			$tssEst = $pos[$i];
		}
	}
	if ($levelN > 0) {
		$startLevel /= $levelN;
	}
	if ($startLevel > $level/4) {
		$tssEst = "NA|$startLevel|$level";
	}
	return $tssEst;
}

sub readRNADefinition {
	my ($file,$repeatFlag) = @_;
	my %rna = ();
	open IN, $file or die "!!! Cannot open RNA definition file: $file\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if (@line < 5);
		my $id = $line[0];
		my $chr = $line[1];
		my $start = $line[2];
		my $end = $line[3];
		my $dir = $line[4];
		next unless ($start =~ /^[\d\-]/);
		my @a= ();
		if (@line > 5) {
			@a = split /\,/,$line[5];
			unless ($a[0] =~ /\:/) {
				@a = ();
				$a[0] = "E1:$start";
			}
		} else {
			$a[0] = "E1:$start";
		}

		my %d = ();
		foreach(@a) {
			my @b = split /\:/;
			my $p = $b[1];
			$d{$p} = $b[0];
		}
		my %e = ();
		if ($dir eq '-' || $dir eq '1') {
			my @z = sort {$a <=> $b} keys %d;
			for (my $i=0;$i<@z;$i++) {
				if ($i+1 < @z) {
					$e{$end-$z[$i+1]} = $d{$z[$i]};
				} else {
					$e{0} = $d{$z[$i]};
				}
			}
		} else {
			foreach(keys %d) {
				my $p = $_ - $start;
				$e{$p} = $d{$_};
			}
		}
		my $lengthExon = 0;
		my $lengthIntron = 0;
		my @bits = sort {$a <=> $b} keys %e;
		for (my $i=0;$i<@bits;$i++) {
			my $secSize = 0;
			if ($i<@bits-1) {
				$secSize = $bits[$i+1]-$bits[$i];
			} else {
				$secSize = ($end-$start)-$bits[$i];
			}
			if ($e{$bits[$i]} =~ /^E/) {
				$lengthExon += $secSize;
			} elsif ($e{$bits[$i]} =~ /^I/) {
				$lengthIntron += $secSize;
			}
		}

		$e{-1e100} = "Upstream";
		$e{$end-$start} = "Downstream";
		my $lengthTotal = $end-$start;
		my %tssReads = ();
		$rna{$id} = {id=>$id,s=>$start,e=>$end,d=>$dir,c=>$chr,v=>\%e,
					Lexons=>$lengthExon,Lintrons=>$lengthIntron,Lgenes=>$lengthTotal,tssReads=>\%tssReads};
	}
	return \%rna;
}

sub readTagPositions {
	my ($file,$strand) = @_;
	my %gene = ();
	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $line[0];
		my %d = ();
		if (@line > 1) {
			my @a = split /\,/, $line[1];
			foreach(@a) {
				my @b = split /\=/;
				my $p = $b[0];
				my @c = split /\|/,$b[1];
				my $v = 0;
				if ($strand ne '-' && $c[0] ne 'NA') {
					if ($init2One > 0 && $c[0] > $init2One) {
						$v += $init2One;
					} else {
						$v += $c[0];
					}
				}
				if ($strand ne '+' && $c[1] ne 'NA') {
					if ($init2One > 0 && $c[1] > $init2One) {
						$v += $init2One;
					} else {
						$v += $c[1];
					}
				}
				$d{$p} = $v if ($v > 0);
			}
		}
		$gene{$id}=\%d;
	}
	close IN;
	return \%gene;
}





