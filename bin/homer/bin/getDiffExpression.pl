#!/usr/bin/env perl
use warnings;



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

sub printCMD {
	print STDERR "\n\tThis program uses DESeq/edgeR to find differential expression between sets of genes\n";
	print STDERR "\t(R must be installed in the executable path, and the edgeR package must be installed)\n";
	print STDERR "\t   Step 1:  Run analyzeRepeats.pl, but use -raw (or analyzeRNA.pl or annotatePeaks.pl)\n";
	print STDERR "\t   Step 2:  Run this program using that file (use -repeats/-rna/-peaks to match program)\n";
	print STDERR "\tThe output is sent to stdout - appends columns to original file containing data\n";
	print STDERR "\tGroup/Batch names below correspond to order of experiments given as tag directories\n";
	print STDERR "\n\tUsage: getDiffExpression.pl <data file> <group code...> [options]\n";
	print STDERR "\n\t\ti.e. getDiffExpression.pl inputFile.txt groupName1 groupName1 groupName2 > output.txt\n";
	print STDERR "\t\t\t(for 2 replicates in first group and 1 in second group)\n";
	print STDERR "\n\tTo normalize batch effects, add \"-batch <batch code 1> [batch code2] ...\" to the command\n";
	print STDERR "\t\tFor now batch normalization is only available with edgeR and DESeq2\n";
	print STDERR "\t\ti.e. getDiffExpression.pl inputFile.txt gName1 gName1 gName2 gName2 -batch 1 2 1 2 > output.txt\n";
	print STDERR "\t\t\t(for 2 replicates in each group, where the 1st in each and 2nd in are paired)\n";
	print STDERR "\n\tInput File format options (will try to autodetect):\n";
	print STDERR "\t\t-rna (for analyzeRNA.pl formatted input, default)\n";
	print STDERR "\t\t-repeats (for analyzeRepeats.pl formatted input file)\n";
	print STDERR "\t\t-peaks (for annotatePeaks.pl formatted input file)\n";
	print STDERR "\t\t-basic (for simple file with one column of gene identifiers and then the count data)\n";
	print STDERR "\n\tAdditional Options:\n";
	print STDERR "\t\t-dispersion <#> (edgeR common dispersion to use if no replicates, default: 0.05)\n";
	print STDERR "\t\t-norm2total (normalize using tag directory totals, default: normalize to gene totals(i.`e.table)\n";
	print STDERR "\t\t-AvsA (compare each group vs. each group, default: compare 1st group vs. others)\n";
	print STDERR "\t\t-showR (Show R status updates, command output)\n";
	print STDERR "\n\tDifferential Expression program selection:\n";
	print STDERR "\t\t-edgeR (use edgeR, - doesn't support -norm2total for normalization to total mapped reads)\n";
	print STDERR "\t\t-DESeq (use DESeq instead of edgeR - doesn't support batch mode)\n";
	print STDERR "\t\t-DESeq2 (use DESeq2 instead of edgeR)\n";
	print STDERR "\n\tDE List Reporting:\n";
	print STDERR "\t\t-export <prefix> (output differential expression gene lists with filename prefix)\n";
	print STDERR "\t\t-fdr <#> (FDR threshold for diff. expression reporting, default: 0.05)\n";
	print STDERR "\t\t-log2fold <#> (Log2 Fold threshold for diff. expression reporting, default: 1.0)\n";
	print STDERR "\n";
	print STDERR "\tCurrent versions that work:\n";
	print STDERR "\t\tR: v3.1.1\n";
	print STDERR "\t\tBioconductor: v2.13\n";
	print STDERR "\t\tDESeq v1.16.0\n";
	print STDERR "\t\tDESeq2 v1.4.5\n";
	print STDERR "\t\tedgeR v3.6.8\n";
	print STDERR "\t\tlimma v3.20.9\n";
	print STDERR "\n";
	exit;
}
if (@ARGV < 3) {
	printCMD();
}


my $cmd = "getDiffExpression.pl";
for (my $i=0;$i<@ARGV;$i++) {
	$cmd .= " $ARGV[$i]";
}

#annotatePeaks
#my $col = 18;
my $normFlag = 1;
my $showRFlag = 0;
#annotateRNA.pl
my $col = -1;
my $commonDispersion = 0.05;
my $tagDirNormFlag = 0;
my $AvsAflag = 0;

my %og = ();
my @header = ();
my @names = ();
my %data = ();

my $mode = 'edgeR';
my %groups= ();
my %batchs= ();
my @groups = ();
my @batchs = ();
my @argvGroups = ();
my @argvBatchs = ();
my $batchID = 1;
my @tagDirTotals = ();
my $fdrCutoff=0.05;
my $log2FoldCutoff=1.0;
my $export = "";


my $rand = rand();

my $batchFlag = 0;
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-batch') {
		$batchFlag = 1;
		next;
	}
	if ($ARGV[$i] eq '-export') {
		$export = $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-AvsA') {
		$AvsAflag = 1;
		next;
	}
	if ($ARGV[$i] eq '-edgeR') {
		$mode = 'edgeR';
		next;
	}
	if ($ARGV[$i] eq '-showR') {
		$showRFlag = 1;
		next;
	}
	if ($ARGV[$i] eq '-DESeq') {
		$mode = 'DESeq';
		next;
	}
	if ($ARGV[$i] eq '-DESeq2') {
		$mode = 'DESeq2';
		next;
	}
	if ($ARGV[$i] eq '-rna') {
		$col = 15;
		next;
	}
	if ($ARGV[$i] eq '-repeats') {
		$col = 8;
		next;
	}
	if ($ARGV[$i] eq '-fdr') {
		$fdrCutoff= $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-log2fold') {
		$log2FoldCutoff= $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] eq '-peaks') {
		$col = 19;
		next;
	}
	if ($ARGV[$i] eq '-basic') {
		$col = 1;
		next;
	}
	if ($ARGV[$i] eq '-norm2total') {
		$tagDirNormFlag = 1;
		next;
	}
	if ($ARGV[$i] eq '-dispersion') {
		$commonDispersion = $ARGV[++$i];
		next;
	}
	if ($ARGV[$i] =~ /^\-/) {
		print STDERR "!!! $ARGV[$i] not recognized!\n";
		printCMD();
	}
	if ($batchFlag == 1) {
		if (!exists($batchs{$ARGV[$i]})) {
			$batchs{$ARGV[$i]} = $batchID++;
			push(@batchs, $ARGV[$i]);
		}
		push(@argvBatchs, $ARGV[$i]);
	} else {
		if (!exists($groups{$ARGV[$i]})) {
			push(@groups, $ARGV[$i]);
		}
		$groups{$ARGV[$i]} = 1;
		push(@argvGroups, $ARGV[$i]);
	}
}
print STDERR "\n\tDifferential Expression Program: $mode\n";
my $numGroups = scalar(@groups);
my $numBatchs = 0;
my $numComp = $numGroups*($numGroups-1);

if ($batchFlag) {
	if ($mode eq 'DESeq') {
		print STDERR "!!! Sorry - this script will only run -batch mode with edgeR or DESeq2\n";
		exit;
	}
}
if ($tagDirNormFlag) {
	if ($mode eq 'edgeR') {
		print STDERR "!!! Unfortunately -norm2total will only work with DESeq or DESeq2\n";
		exit;
	}
}



if ($col == 15) {
	print STDERR "\tTreating input as file generated by analyzeRNA.pl (-rna)\n";
} elsif ($col == 8) {
	print STDERR "\tTreating input as file generated by analyzeRepeats.pl (-repeats)\n";
} elsif ($col == 19) {
	print STDERR "\tTreating input as file generated by annotatePeaks.pl (-peaks)\n";
} elsif ($col == -1) {
	print STDERR "\tAutodetecting input file format...\n";
}


if ($batchFlag) {
	print STDERR "\tWill attempt to remove batch effects\n";
	if (scalar(@argvGroups) != scalar(@argvBatchs)) {
		print STDERR "!!!! Entered different number of groups and batches on command line !!!!\n";
		exit;
	}
	$numBatchs = scalar(@groups);
	if ($numBatchs < 2) {
		print STDERR "!!!! Warning - you need at least 2 batch codes... otherwise omit \"-batch\" !!!\n";
		exit;
	}
}

print STDERR "\tUsing $mode to calculate differential expression/enrichment...\n";

my @colTotals = ();

my $defHeader = '';

open IN, $ARGV[0];
my $count = 0;
while (<IN>) {
	$count++;
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $str = $line[0];

	if ($count == 1 && $col == -1) {
		if ($line[0] =~ /analyzeRepeats\.pl/) {
			$col = 8;
			print STDERR "\tAutodetected analyzeRepeats.pl file\n";
		} elsif ($line[0] =~ /analyzeRNA\.pl/) {
			$col = 15;
			print STDERR "\tAutodetected analyzeRNA.pl file\n";
		} elsif ($line[0] =~ /annotatePeaks\.pl/) {
			$col = 19;
			print STDERR "\tAutodetected annotatePeaks.pl file\n";
		} else {
			print STDERR "!!! Warning - could not autodetect the file format !!!\n";
			print STDERR "!!! Assuming it might be a simple/basic format with one column of geneIDs... !!\n";
			$col = 1;
		}
	}
	if ($count == 1) {
		$str .= " (cmd=$cmd)";
	}
			


	for (my $i=1;$i<@line;$i++) {
		$str .= "\t$line[$i]";
	}

	if ($count == 1) {
		print "$str";
		$defHeader = $str;
		for (my $i=$col;$i<@line;$i++) {
			$line[$i] =~ s/\"//g;
			$line[$i] =~ /^(.*?)\s/;
			push(@names, $1);
			push(@colTotals, 0);

			my $tagDirTotal=0;
			if ($line[$i] =~ /\(([\d\.]+) Total,/) {
				$tagDirTotal = $1;
				#print STDERR "\t$line[$i] => $tagDirTotal\n";
			}
			push(@tagDirTotals, $tagDirTotal);

		}
		next;
	}	
	$og{$line[0]} = $str;

	my @d = ();
	for (my $i=$col;$i<@line;$i++) {
		push(@d,$line[$i]);
		$colTotals[$i-$col] += $line[$i];
	}
	$data{$line[0]} = \@d;
}
close IN;

my $avgTotal = 0;
foreach(@colTotals) {
	$avgTotal+=$_;
}
if ($tagDirNormFlag == 1) {
	$avgTotal = 0;
	foreach(@tagDirTotals) {
		$avgTotal+=$_;
	}
}

$avgTotal /= scalar(@colTotals);
if ($normFlag) {
	foreach(keys %og) {
		my $id = $_;
		my $v = $og{$id};
		my @a = split /\t/,$v;
		for (my $i=$col;$i<@a;$i++) {
			if ($tagDirNormFlag == 1) {
				$a[$i] = sprintf("%.2lf",$a[$i]*$avgTotal/$tagDirTotals[$i-$col]);
			} else {
				$a[$i] = sprintf("%.2lf",$a[$i]*$avgTotal/$colTotals[$i-$col]);
			}
		}
		$v = $a[0];
		for (my $i=1;$i<@a;$i++) {
			$v .= "\t$a[$i]";
		}
		$og{$id} = $v;
	}
}

my $min = 1e10;
my $extra = '';
my @comps = ();
for (my $i=0;$i<$numGroups;$i++) {
	last if ($i>0 && $AvsAflag == 0);
	for (my $j=$i+1;$j<$numGroups;$j++) {


		my $replicateFlag = 0;

		my $inputFile = $rand . ".edgeR.in.data.txt";
		my $gFile = $rand . ".edgeR.groups.data.txt";
		my $bFile = $rand . ".edgeR.batch.data.txt";
		my $normFile = $rand . ".edgeR.norm.data.txt";
		my $libSizeFile = $rand . ".edgeR.libSize.data.txt";
		my $deseq2File = $rand . ".deseq2.groups.data.txt";
		my $rOutputFile = $rand . ".edgeR.out.data.txt";
		my $rStderr = $rand . ".edgeR.out.stderr.txt";
		my $rScript = $rand . ".edgeR.R";
		open OUT, ">$inputFile";
		open DESEQ2, ">$deseq2File";
		open GROUPS, ">$gFile";
		open BATCH, ">$bFile";
		open NORM, ">$normFile";
		open LIBSIZE, ">$libSizeFile";
		my $gline = "";
		my $bline = "";
		my @dsamples = ();
		my @dgline = ();
		my @dbline = ();
		my %sampleCounts = ();
		my @totals = ();

		for (my $k=0;$k<@argvGroups;$k++) {
			if ($argvGroups[$k] eq $groups[$i]) {
				#print OUT "\t$i";
				print OUT "\t$k";
				push(@dsamples, $k);
				$sampleCounts{$i}++;
				$gline .= "\t$i";
				push(@totals, $tagDirTotals[$k]);
				push(@dgline, $i);
				if ($batchFlag) {
					my $b = $batchs{$argvBatchs[$k]};
					$bline .= "\t$b";
					push(@dbline, $b);
				}
			} elsif ($argvGroups[$k] eq $groups[$j]) {
				#print OUT "\t$j";
				print OUT "\t$k";
				push(@dsamples, $k);
				$sampleCounts{$j}++;
				$gline .= "\t$j";
				push(@totals, $tagDirTotals[$k]);
				push(@dgline, $j);
				if ($batchFlag) {
					my $b = $batchs{$argvBatchs[$k]};
					$bline .= "\t$b";
					push(@dbline, $b);
				}
			}
		}
		print OUT "\n";

		foreach(values %sampleCounts) {
			if ($_ > 1) {
				$replicateFlag = 1;
			}
		}

		if ($tagDirNormFlag == 1) {
			my $T = 0;
			$T+=$_ foreach(@totals);
			my $avg = $T/@totals if (@totals > 0);
			my $normLine="";
			my $sizeLine="";
			foreach(@totals) {
				$sizeLine .= "\t$_";
				my $v = 1;
				if ($_ > 0) {
					$v = $_/$avg;
				}
				$normLine .= "\t$v";
			}
			print NORM "$gline\n";
			print NORM "X" . "$normLine\n";
			print LIBSIZE "$gline\n";
			print LIBSIZE "X" . "$sizeLine\n";
		}
		close NORM;
		close LIBSIZE;
		print GROUPS "$gline\n";
		print GROUPS "X" . "$gline\n";
		close GROUPS;
		print BATCH "$bline\n";
		print BATCH "X" . "$bline\n";
		close BATCH;
		print DESEQ2 "\tgroup";
		print DESEQ2 "\tbatch" if ($batchFlag);
		print DESEQ2 "\n";
		for (my $z=0;$z<@dsamples;$z++) {
			print DESEQ2 "$dsamples[$z]\t$dgline[$z]";
			print DESEQ2 "\t$dbline[$z]" if ($batchFlag);
			print DESEQ2 "\n";
		}
		close DESEQ2;


		foreach(keys %data) {
			my $gene = $_;
			print OUT "$gene";
			for (my $k=0;$k<@argvGroups;$k++) {
				if ($argvGroups[$k] eq $groups[$i] || $argvGroups[$k] eq $groups[$j]) {
					my $v = $data{$gene}->[$k];
					$v = floor($v+0.5);
					print OUT "\t$v";
				}
			}
			print OUT "\n";
		}
		close OUT;

		open SCRIPT, ">$rScript";
		if ($mode eq 'edgeR') {
			print SCRIPT "####### basic script for running edgeR (generated by HOMER) ########\n";
			print SCRIPT "#### input is matrix data - first column is gene id, others are tag un-normalized tag counts\n";
			print SCRIPT "#### <blank> Exp1 Exp2\n#### NM_123   1233  444\n##### NM_234    424  198\n";
			print SCRIPT "library(edgeR)\n";
			print SCRIPT "#Read Data in\n";
			print SCRIPT "data <- as.matrix(read.table(\"$inputFile\"))\n";
			print SCRIPT "gin <- as.matrix(read.table(\"$gFile\"))\n\n";
			print SCRIPT "#set groups\n";
			print SCRIPT "g <- gin[1:length(gin)];\n";
			print SCRIPT "group=factor(g)\n";
			if ($tagDirNormFlag == 1) {
				print SCRIPT "normFactors <- as.vector(as.matrix(read.table(\"$normFile\")))\n";
				print SCRIPT "libSizes <- as.vector(as.matrix(read.table(\"$libSizeFile\")))\n";
				print SCRIPT "d <- DGEList(counts=data,group=g,lib.size=libSizes,norm.factors=normFactors)\n\n";
			} else {
				print SCRIPT "libSizes <- as.vector(colSums(data))\n";
				print SCRIPT "d <- DGEList(counts=data,group=g,lib.size=libSizes)\n\n";
				print SCRIPT "d <- calcNormFactors(d)\n\n";
			}
	
			if ($batchFlag) {
				print SCRIPT "bin <- as.matrix(read.table(\"$bFile\"))\n\n";
				print SCRIPT "b <- bin[1:length(bin)]\n";
				print SCRIPT "batch=factor(b)\n";
	
				print SCRIPT "design <-model.matrix(~batch+d\$samples\$group)\n";
				print SCRIPT "d <- estimateGLMCommonDisp(d, design)\n";
				print SCRIPT "d <- estimateGLMTrendedDisp(d,design)\n";
					print SCRIPT "d <- estimateGLMTagwiseDisp(d, design)\n";
				print SCRIPT "glmfit <- glmFit(d, design, dispersion = d\$tagwise.dispersion)\n";
				print SCRIPT "lrt <- glmLRT(glmfit, coef=dim(design)[[2]] )\n";
				print SCRIPT "results <- topTags(lrt,n = length(data))\n";
			} else {
				print SCRIPT "#Estimate Error Model\n";
				if ($replicateFlag) {
					print SCRIPT "d <- estimateCommonDisp(d)\n";
					print SCRIPT "d <- estimateTagwiseDisp(d)\n";
				} else {
					print SCRIPT "d\$common.dispersion <- $commonDispersion\n";
				}
				print SCRIPT "#compute p-values, then output\n";
				print SCRIPT "de.com <- exactTest(d)\n";
				print SCRIPT "results <- topTags(de.com,n = length(data))\n";
			}
			print SCRIPT "write.table(as.matrix(results\$table),file=\"$rOutputFile\",sep=\"\\t\")\n\n";
			print SCRIPT "write(\"Dispersion = \",stderr())\n";
			print SCRIPT "write(d\$common.dispersion,stderr())\n";
		} elsif ($mode eq 'DESeq') {
			print SCRIPT "####### basic script for running DESeq (generated by HOMER) ########\n";
			print SCRIPT "library(DESeq)\n";
			print SCRIPT "#Read Data in\n";
			print SCRIPT "data <- as.matrix(read.table(\"$inputFile\"))\n";
			print SCRIPT "gin <- as.matrix(read.table(\"$gFile\"))\n\n";
			print SCRIPT "#set groups\n";
			print SCRIPT "g <- gin[1:length(gin)];\n";

			if ($batchFlag) {
				print SCRIPT "bin <- as.matrix(read.table(\"$bFile\"))\n\n";
				print SCRIPT "b <- bin[1:length(bin)]\n";
				print SCRIPT "batch=factor(b)\n";
				#print SCRIPT "dds <- DESeqDataSetFromMatrix(countData = data, colData = gin,design=~batch+group)\n";
				print SCRIPT "cds = newCountDataSet(data,~batch+g)\n";
			} else {
				print SCRIPT "cds = newCountDataSet(data,g)\n";
			}
			if ($tagDirNormFlag == 1) {
				print SCRIPT "normFactors <- as.vector(as.matrix(read.table(\"$normFile\")))\n";
				print SCRIPT "pData(cds)\$sizeFactor <- normFactors\n";
			
			} else {
				print SCRIPT "cds = estimateSizeFactors(cds)\n";
			}
			if ($replicateFlag) {
				print SCRIPT "cds = estimateDispersions(cds)\n";
			} else {
				print SCRIPT "cds = estimateDispersions(cds,method=\"blind\", sharingMode=\"fit-only\")\n";
			}
			print SCRIPT "res = nbinomTest(cds,\"0\",\"1\")\n";
			print SCRIPT "write.table(as.matrix(res),file=\"$rOutputFile\",sep=\"\\t\")\n\n";
		} elsif ($mode eq 'DESeq2') {
			print SCRIPT "####### basic script for running DESeq (generated by HOMER) ########\n";
			print SCRIPT "library(DESeq2)\n";
			print SCRIPT "#Read Data in\n";
			print SCRIPT "data <- as.matrix(read.table(\"$inputFile\"))\n";
			print SCRIPT "gin <- read.table(\"$deseq2File\")\n\n";

			if ($batchFlag) {
				#print SCRIPT "design <-model.matrix(~batch+d\$samples\$group)\n";
				print SCRIPT "bin <- as.matrix(read.table(\"$bFile\"))\n\n";
				print SCRIPT "b <- bin[1:length(bin)]\n";
				print SCRIPT "batch=factor(b)\n";
				print SCRIPT "dds <- DESeqDataSetFromMatrix(countData = data, colData = gin,design=~batch+group)\n";
			} else {
				print SCRIPT "dds <- DESeqDataSetFromMatrix(countData = data, colData = gin,design=~group)\n";
			}

			if ($tagDirNormFlag == 1) {
				print SCRIPT "normFactors <- as.vector(as.matrix(read.table(\"$normFile\")))\n";
				print SCRIPT "colData(dds)\$sizeFactor <- normFactors\n";
				print SCRIPT "dds <- estimateDispersions(dds)\n";
				print SCRIPT "dds <- nbinomWaldTest(dds)\n";

			} else {
				print SCRIPT "dds <- DESeq(dds)\n";
			}
			print SCRIPT "res <- results(dds)\n";
			print SCRIPT "write.table(as.matrix(res),file=\"$rOutputFile\",sep=\"\\t\")\n\n";
		}
		close SCRIPT;

		#print STDERR "`R --no-save < $rScript`;\n";
		if ($showRFlag) {
			`R --no-save < "$rScript"`;
		} else {
			`R --no-save < "$rScript" 2> "$rStderr"`;
		}

		my $numUp = 0;
		my $numDown = 0;
		my $numTotal = 0;

		if ($export ne '') {
			open UP, ">$export.Up_$groups[$j]" . "_vs_" . $groups[$i] . ".txt";
			open DOWN, ">$export.Down_$groups[$j]" . "_vs_" . $groups[$i] . ".txt";
		}

		my %comp = ();
		open IN, $rOutputFile;
		$count = 0;
		while (<IN>) {
			$count++;
			chomp;
			s/\r//g;
			s/\"//g;
			my @line = split /\t/;
			foreach(@line) {
				s/^\s+//;
				s/\s+$//;
			}


			if ($count == 1) {
				print "\t$groups[$i] vs. $groups[$j] Log2 Fold Change";
				print "\t$groups[$i] vs. $groups[$j] p-value";
				print "\t$groups[$i] vs. $groups[$j] adj. p-value";
				if ($export ne '') {
					print UP "$defHeader";
					print UP "\t$groups[$i] vs. $groups[$j] Log2 Fold Change";
					print UP "\t$groups[$i] vs. $groups[$j] p-value";
					print UP "\t$groups[$i] vs. $groups[$j] adj. p-value";
					print UP "\n";
					print DOWN "$defHeader";
					print DOWN "\t$groups[$i] vs. $groups[$j] Log2 Fold Change";
					print DOWN "\t$groups[$i] vs. $groups[$j] p-value";
					print DOWN "\t$groups[$i] vs. $groups[$j] adj. p-value";
					print DOWN "\n";
				}
				next;
			}


			my $id = '';
			my $f = 0.0;
			my $p = 1.0;
			my $fdr = 1.0;
			my $str = "";

			if ($mode eq 'edgeR') {
				$id = $line[0];
				$f = $line[1];
				$p = $line[3];
				$fdr = $line[4];
				if ($batchFlag) {
					$p = $line[4];
					$fdr = $line[5];
				}
			} elsif ($mode eq 'DESeq') {
				$id = $line[1];
				$f = $line[6];
				$p = $line[7];
				$fdr = $line[8];
			} elsif ($mode eq 'DESeq2') {
				$id = $line[0];
				$f = $line[2];
				$p = $line[5];
				$fdr = $line[6];
			}

			if ($f eq '-Inf' || $f eq '-inf') {
				$f = -10;
			}
			if ($f eq 'Inf' || $f eq 'inf') {
				$f = 10;
			}
			if ($fdr eq 'NA' || $fdr > 1) {
				$fdr = 1;
			}
			if ($p eq 'NA' || $p > 1) {
				$p = 1;
			}

			if (!exists($og{$id})) {
				next;
			}

			$str = "\t$f\t$p\t$fdr";

			$numTotal++;
			if ($fdr < $fdrCutoff) {
				if ($f > $log2FoldCutoff) {
					$numUp++;
					print UP "$og{$id}" . "$str\n" if ($export ne '');
				} elsif ($f < -1*$log2FoldCutoff) {
					$numDown++;
					print DOWN "$og{$id}" . "$str\n" if ($export ne '');
				}
			}

			#$min = $line[$sIndex+1] if ($line[$sIndex+1] ne 'NA' && $line[$sIndex+1] < $min);
			$comp{$id}=$str;

		}
		close IN;

		if ($export ne '') {
			close UP;
			close DOWN;
		}

		push(@comps, \%comp);


		if ($numTotal == 0) {
			print STDERR "!!! Warnining - something likely failed during R execution.\n";
			print STDERR "!!! R script that was used:\n";
			open IN1, "$rScript";	
			while (<IN1>) {
				print STDERR "\t$_";
			}
			close IN1;
			if ($showRFlag==0) {
				print STDERR "!!! R execution details:\n";
				open IN1, "$rStderr";	
				while (<IN1>) {
					print STDERR "\t$_";
				}
				close IN1;
			}
		}
		if ($showRFlag==1) {
			open IN1, "$rScript";	
			while (<IN1>) {
				print STDERR "\t$_";
			}
			close IN1;

		}

		`rm "$inputFile" "$gFile" "$bFile" "$rOutputFile" "$rScript" "$deseq2File" "$normFile" "$libSizeFile"`;
		if ($showRFlag==0) {
			`rm "$rStderr"`;
		}

		if ($numTotal == 0) {
			exit;
		}

		print STDERR "\tOutput Stats $groups[$i] vs. $groups[$j]:\n";
		print STDERR "\t\tTotal Genes: $numTotal\n";
		my $r = sprintf("%.3f",$numUp/$numTotal*100);
		print STDERR "\t\tTotal Up-regulated in $groups[$j] vs. $groups[$i]: $numUp ($r%) [log2fold>$log2FoldCutoff, FDR<$fdrCutoff]\n";
		$r = sprintf("%.3f",$numDown/$numTotal*100);
		print STDERR "\t\tTotal Dn-regulated in $groups[$j] vs. $groups[$i]: $numDown ($r%) [log2fold<-$log2FoldCutoff, FDR<$fdrCutoff]\n";
	}
}

print "\n";

foreach(keys %og) {
	print "$og{$_}";
	my $gene = $_;
	foreach(@comps) {
		if (exists($_->{$gene})) {
			print $_->{$gene};
		} else {
			print "\t0\t1\t1";
		}
	}
	print "\n";
}
print STDERR "\n";
