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

my $reduceThresh = 0.6;
my $percentSimilar = 0.20;
my $matchThresh = "T10";
my $knownPvalueThresh = 0.01;
my $motifInfoFileName = "motifFindingParameters.txt";
my $deleteFlag = 1;
my %toDelete = ();

my $config = HomerConfig::loadConfigFile();

if (@ARGV < 3) {
	printCMD();
}

sub printCMD {
	print STDERR "\n\tProgram will find de novo and known motifs in regions in the genome\n";
	print STDERR "\n\tUsage: findMotifsGenome.pl <pos file> <genome> <output directory> [additional options]\n";
	print STDERR "\tExample: findMotifsGenome.pl peaks.txt mm8r peakAnalysis -size 200 -len 8\n";

	print STDERR "\n\tPossible Genomes:\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\n";
	}
	print STDERR "\t\t\t-- or --\n";
    print STDERR "\t\tCustom: provide the path to genome FASTA files (directory or single file)\n";
    print STDERR "\t\t\tHeads up: will create the directory \"preparsed/\" in same location.\n";

	print STDERR "\n\tBasic options:\n";
	print STDERR "\t\t-mask (mask repeats/lower case sequence, can also add 'r' to genome, i.e. mm9r)\n";
	print STDERR "\t\t-bg <background position file> (genomic positions to be used as background, default=automatic)\n";
	print STDERR "\t\t\tremoves background positions overlapping with target positions\n";
	print STDERR "\t\t\t-chopify (chop up large background regions to the avg size of target regions)\n";
	print STDERR "\t\t-len <#>[,<#>,<#>...] (motif length, default=8,10,12) [NOTE: values greater 12 may cause the program\n";
	print STDERR "\t\t\tto run out of memory - in these cases decrease the number of sequences analyzed (-N),\n";
	print STDERR "\t\t\tor try analyzing shorter sequence regions (i.e. -size 100)]\n";
	print STDERR "\t\t-size <#> (fragment size to use for motif finding, default=200)\n";
	print STDERR "\t\t\t-size <#,#> (i.e. -size -100,50 will get sequences from -100 to +50 relative from center)\n";
	print STDERR "\t\t\t-size given (uses the exact regions you give it)\n";
	print STDERR "\t\t-S <#> (Number of motifs to optimize, default: 25)\n";
	print STDERR "\t\t-mis <#> (global optimization: searches for strings with # mismatches, default: 2)\n";
	print STDERR "\t\t-norevopp (don't search reverse strand for motifs)\n";
	print STDERR "\t\t-nomotif (don't search for de novo motif enrichment)\n";
	print STDERR "\t\t-rna (output RNA motif logos and compare to RNA motif database, automatically sets -norevopp)\n";

	print STDERR "\n\tScanning sequence for motifs\n";
	print STDERR "\t\t-find <motif file> (This will cause the program to only scan for motifs)\n";

	print STDERR "\n\tKnown Motif Options/Visualization\n";
	print STDERR "\t\t-mset <vertebrates|insects|worms|plants|yeast|all> (check against motif collects, default: auto)\n";
	print STDERR "\t\t-basic (just visualize de novo motifs, don't check similarity with known motifs)\n";
	print STDERR "\t\t-bits (scale sequence logos by information content, default: doesn't scale)\n";
	print STDERR "\t\t-nocheck (don't search for de novo vs. known motif similarity)\n";
	print STDERR "\t\t-mcheck <motif file> (known motifs to check against de novo motifs,\n";
	#print STDERR "\t\t-float (allow adjustment of the degeneracy threshold for known motifs to improve p-value[dangerous])\n";
	print STDERR "\t\t-noknown (don't search for known motif enrichment, default: -known)\n";
	print STDERR "\t\t-mknown <motif file> (known motifs to check for enrichment,\n";
	print STDERR "\t\t-nofacts (omit humor)\n";
	
	print STDERR "\n\tSequence normalization options:\n";
	#print STDERR "\t\t-tss (normalize based on distance from TSS)\n";
	#print STDERR "\t\t-cgtss (normalize based on CpG content and distance from TSS)\n";
	#print STDERR "\t\t-cg DEFAULT (normalize based on CpG content)\n";
	print STDERR "\t\t-gc (use GC% for sequence content normalization, now the default)\n";
	print STDERR "\t\t-cpg (use CpG% instead of GC% for sequence content normalization)\n";
	print STDERR "\t\t-noweight (no CG correction)\n";
	print STDERR "\t\tAlso -nlen <#>, -olen <#>, see homer2 section below.\n";

	print STDERR "\n\tAdvanced options:\n";
	print STDERR "\t\t-h (use hypergeometric for p-values, binomial is default)\n";
	print STDERR "\t\t-N <#> (Number of sequences to use for motif finding, default=max(50k, 2x input)\n";
	#print STDERR "\t\t-noforce (will attempt to reuse sequence files etc. that are already in output directory)\n";
	print STDERR "\t\t-local <#> (use local background, # of equal size regions around peaks to use i.e. 2)\n";
	print STDERR "\t\t-redundant <#> (Remove redundant sequences matching greater than # percent, i.e. -redundant 0.5)\n";
	print STDERR "\t\t-maxN <#> (maximum percentage of N's in sequence to consider for motif finding, default: 0.7)\n";
	print STDERR "\t\t-maskMotif <motif file1> [motif file 2]... (motifs to mask before motif finding)\n";
	print STDERR "\t\t-opt <motif file1> [motif file 2]... (motifs to optimize or change length of)\n";
	#print STDERR "\t\t-refine <motif file1> (motif to optimize)\n";
	print STDERR "\t\t-rand (randomize target and background sequences labels)\n";
	print STDERR "\t\t-ref <peak file> (use file for target and background - first argument is list of peak ids for targets)\n";
	print STDERR "\t\t-oligo (perform analysis of individual oligo enrichment)\n";
	print STDERR "\t\t-dumpFasta (Dump fasta files for target and background sequences for use with other programs)\n";
	print STDERR "\t\t-preparse (force new background files to be created)\n";
	print STDERR "\t\t-preparsedDir <directory> (location to search for preparsed file and/or place new files)\n";
	print STDERR "\t\t-keepFiles (keep temporary files)\n";
	print STDERR "\t\t-fdr <#> (Calculate empirical FDR for de novo discovery #=number of randomizations)\n";

	print STDERR "\n";
	print STDERR "\thomer2 specific options:\n";
	print STDERR "\t\t-homer2 (use homer2 instead of original homer, default)\n";
	print STDERR "\t\t-nlen <#> (length of lower-order oligos to normalize in background, default: -nlen 3)\n";
	print STDERR "\t\t\t-nmax <#> (Max normalization iterations, default: 160)\n";
	print STDERR "\t\t\t-neutral (weight sequences to neutral frequencies, i.e. 25%, 6.25%, etc.)\n";
	print STDERR "\t\t-olen <#> (lower-order oligo normalization for oligo table, use if -nlen isn't working well)\n";
	print STDERR "\t\t-p <#> (Number of processors to use, default: 1)\n";
	print STDERR "\t\t-e <#> (Maximum expected motif instance per bp in random sequence, default: 0.01)\n";
	#print STDERR "\t\t-ps <#> (remove motifs sharing > # sites, default: 0.2 [i.e. 20%])\n";
	print STDERR "\t\t-cache <#> (size in MB for statistics cache, default: 500)\n";
	print STDERR "\t\t-quickMask (skip full masking after finding motifs, similar to original homer)\n";
	print STDERR "\t\t-minlp <#> (stop looking for motifs when seed logp score gets above #, default: -10)\n";
	print STDERR "\n\tOriginal homer specific options:\n";
	print STDERR "\t\t-homer1 (to force the use of the original homer)\n";
	print STDERR "\t\t-depth [low|med|high|allnight] (time spent on local optimization default: med)\n";
	print STDERR "\n";
	exit;
}


my $cmd = parseCMDLine(\@ARGV);

sub parseCMDLine {
	my ($argv) = @_;

	my @len = ();
	my @motifFiles = ();
	my @motifFiles2 = ();
	my $cmd = {posfile=>$ARGV[0], fg=>'', bg=>'', 
				genome=>$ARGV[1], output=>$ARGV[2], 
				size=>200, redundant=>2, gc=>1,
				cgweight=>1,tssweight=>0,motifMask=>\@motifFiles,motifOpt=>\@motifFiles2,
				nomotif=>0,noknown=>0, groupFlag=>0, norevopp=>0,local=>0,
				motif=>'', nomask=>0, len=>\@len, mis=>2, numSeq=>50000,rnaMode=>0,mask=>0,
				S=>25, reuse=>0, find=>'', alg=>' -alg binomial ', maxN=>0.7,
				force=>1, rand=>0,checkFlag=>1, depth=>0.5,float=>0,
				mknown=>'', mcheck=>'', refFlag=>0,
				oligoFlag=>0,dumpFasta=>0,preparseForce=>0,keepFiles=>0,
				homer2=>1,nlen=>3,olen=>0,cpus=>1,expect=>0,percentSimilar=>0.20,
				cache=>500,bits=>"",nmax=>160,quickMask=>0,sizeMove=>0,chopify=>0,
				nofacts=>"",fdr=>0,neutral=>"",minlp=>-10,mset=>'auto',org=>'',
				preparsedDir=>'/'
			};
	print STDERR "\n";

	my $genome = $ARGV[1];
	if ($genome =~ s/r$//) {
		$cmd->{'mask'} = 1;
	}
	$cmd->{'genome'} = $genome;
	print STDERR "\tPosition file = $cmd->{'posfile'}\n";
	print STDERR "\tGenome = $cmd->{'genome'}\n";
	print STDERR "\tOutput Directory = $cmd->{'output'}\n";
	
	for (my $i=3;$i<@$argv;$i++) { 
		if ($ARGV[$i] eq '-bg') {
			$cmd->{'bg'} = $ARGV[++$i];
			print STDERR "\tbackground position file: $cmd->{'bg'}\n";
		} elsif ($ARGV[$i] eq '-chopify') {
			$cmd->{'chopify'} = 1;
			print STDERR "\tChopping up large background regions to the avg size of target regions\n";
		} elsif ($ARGV[$i] eq '-len' ) {
			my @lengths = split /\,/, $ARGV[++$i];
			$cmd->{'len'} = \@lengths;
			my $str = '';
			print STDERR "\tMotif length set at ";
			foreach(@lengths) {
				print STDERR "$_,";
				if ($_ > 12) {
					print STDERR "*** might want to hit ctrl+C to quit if it takes too long!\n";
					print STDERR "*** If you run out of memmory try reducing the number background sequences\n";
				}
			}
			print STDERR "\n";
		} elsif ($ARGV[$i] eq '-fdr' ) {
			$cmd->{'fdr'} = $ARGV[++$i];
			print STDERR "\tWill randomize and repeat motif finding $cmd->{'fdr'} times to estimate FDR\n";
		} elsif ($ARGV[$i] eq '-minlp' ) {
			$cmd->{'minlp'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-preparsedDir' ) {
			$cmd->{'preparsedDir'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-nofacts' ) {
			$cmd->{'nofacts'} .= " -nofacts ";
		} elsif ($ARGV[$i] eq '-mset' ) {
			$cmd->{'mset'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-maxN' ) {
			$cmd->{'maxN'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-basic' ) {
			$cmd->{'nofacts'} .= " -basic ";
		} elsif ($ARGV[$i] eq '-neutral' ) {
			$cmd->{'neutral'} = " -neutral ";
		} elsif ($ARGV[$i] eq '-oligo' ) {
			$cmd->{'oligoFlag'} = 1; 
			print STDERR "\tWill perform analysis of enriched oligos\n";
		} elsif ($ARGV[$i] eq '-mask' ) {
			print STDERR "\tWill use repeat masked sequences\n";
			$cmd->{'mask'} = 1;
		} elsif ($ARGV[$i] eq '-norevopp' ) {
			$cmd->{'norevopp'} = 1; 
			print STDERR "\tWill not search the reverse strand\n";
		} elsif ($ARGV[$i] eq '-keepFiles' ) {
			$cmd->{'keepFiles'} = 1; 
			print STDERR "\tWill keeep temp files\n";
		} elsif ($ARGV[$i] eq '-opt') {
			my $bail = 0;
			print STDERR "\tWill optimize motifs in the following files:\n";
			print STDERR "\t\tskipping known enrichment...\n";
			while ($ARGV[++$i] !~ /^\-/) {
				print STDERR "\t\t$ARGV[$i]\n";
				push(@{$cmd->{'motifOpt'}}, $ARGV[$i]);
				if ($i>=@ARGV-1) {
					$bail=1;
					last;
				}
			}
			$cmd->{'noknown'} = 1;
			last if ($bail==1);
			$i--;
		} elsif ($ARGV[$i] eq '-maskMotif') {
			my $bail = 0;
			print STDERR "\tWill mask motifs in the following files:\n";
			while ($ARGV[++$i] !~ /^\-/) {
				print STDERR "\t\t$ARGV[$i]\n";
				push(@{$cmd->{'motifMask'}}, $ARGV[$i]);
				if ($i>=@ARGV-1) {
					$bail=1;
					last;
				}
			}
			last if ($bail==1);
			$i--;
		} elsif ($ARGV[$i] eq '-depth' ) {
			if ($ARGV[$i+1] eq 'low') {
				$cmd->{'depth'} = 1;
			} elsif ($ARGV[$i+1] eq 'med') {
				$cmd->{'depth'} = 0.5;
			} elsif ($ARGV[$i+1] eq 'high') {
				$cmd->{'depth'} = 0.1;
			} elsif ($ARGV[$i+1] eq 'allnight') {
				$cmd->{'depth'} = 0.01;
			} else {
				print STDERR "Don't understand $ARGV[$i+1] as an optimization depth!\n";
				exit;
			}
			print STDERR "\tLocal optimization set on $ARGV[$i+1] depth\n";
			$i++;
		} elsif ($ARGV[$i] eq '-rna' ) {
			$cmd->{'noknown'} = 1; 
			$cmd->{'rnaMode'} = 1; 
			$cmd->{'norevopp'} = 1; 
			$cmd->{'mknown'} = $homeDir . "/data/knownTFs/known.rna.motifs";
			$cmd->{'mcheck'} = $homeDir . "/data/knownTFs/all.rna.motifs";
			print STDERR "\tOperating in RNA mode\n";
		} elsif ($ARGV[$i] eq '-reuse' ) {
			$cmd->{'reuse'} = 1; 
			print STDERR "\tOld files will be used to repeat homer analysis\n";
		} elsif ($ARGV[$i] eq '-ref' ) {
			$cmd->{'refFlag'} = 1; 
			$cmd->{'bg'} = $ARGV[++$i];
			print STDERR "\tUsing a reference peak position file\n";
		} elsif ($ARGV[$i] eq '-preparse' ) {
			$cmd->{'preparseForce'} = 1;
		} elsif ($ARGV[$i] eq '-dumpFasta' ) {
			$cmd->{'dumpFasta'} = 1;
			print STDERR "\tOutputing target.fa and background.fa\n";
		} elsif ($ARGV[$i] eq '-float' ) {
			$cmd->{'float'} = 1;
			print STDERR "\tAllowing known motifs to be optimized for degeneracy\n";
		} elsif ($ARGV[$i] eq '-h' ) {
			$cmd->{'alg'} = ''; 
			print STDERR "\tUsing hypergeometric distribution for p-values\n";
		} elsif ($ARGV[$i] eq '-gc' ) {
			$cmd->{'gc'} = 1; 
			print STDERR "\tNormalizing sequences using GC% instead of CpG%\n";
		} elsif ($ARGV[$i] eq '-cpg' ) {
			$cmd->{'gc'} = 0;
			print STDERR "\tNormalizing sequences using CpG% instead of overall GC%\n";
		} elsif ($ARGV[$i] eq '-quickMask' ) {
			$cmd->{'quickMask'} = 1; 
			print STDERR "\tWill employ quick masking...\n";
		} elsif ($ARGV[$i] eq '-noforce' ) {
			$cmd->{'force'} = 0; 
			print STDERR "\tWill for the creation of intermediate files\n";
		} elsif ($ARGV[$i] eq '-homer1' ) {
			$cmd->{'homer2'} = 0;
			print STDERR "\tForcing the use of the original homer\n";
		} elsif ($ARGV[$i] eq '-homer2' ) {
			$cmd->{'homer2'} = 1;
			print STDERR "\tUsing homer2 (warning, unpredictably awesome results may, or may not, ensue)\n";
		} elsif ($ARGV[$i] eq '-p' || $ARGV[$i] eq '-cpu') {
			$cmd->{'cpus'} = $ARGV[++$i];
			print STDERR "\tUsing $cmd->{'cpus'} CPUs\n";
			print STDERR "\t\tWill only work with homer2 (i.e. use -homer2)\n" if ($cmd->{'homer2'} == 0);
		} elsif ($ARGV[$i] eq '-cache' ) {
			$cmd->{'cache'} = $ARGV[++$i];
			print STDERR "\tUsing $cmd->{'cache'} MB for statistics cache\n";
			print STDERR "\t\tWill only work with homer2 (i.e. use -homer2)\n" if ($cmd->{'homer2'} == 0);
		} elsif ($ARGV[$i] eq '-olen' ) {
			$cmd->{'olen'} = $ARGV[++$i];
			print STDERR "\tWill normalize background oligos for oligos of length $cmd->{'olen'} bp\n";
			print STDERR "\t\tWill only work with homer2 (i.e. use -homer2)\n" if ($cmd->{'homer2'} == 0);
		} elsif ($ARGV[$i] eq '-nmax' ) {
			$cmd->{'nmax'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-nlen' ) {
			$cmd->{'nlen'} = $ARGV[++$i];
			print STDERR "\tWill normalize background sequences for oligos of length $cmd->{'nlen'} bp\n";
			print STDERR "\t\tWill only work with homer2 (i.e. use -homer2)\n" if ($cmd->{'homer2'} == 0);
		} elsif ($ARGV[$i] eq '-pc' ) {
			$cmd->{'percentSimilar'} = $ARGV[++$i];
			print STDERR "\tWill remove motifs sharing more than $cmd->{'percentSimilar'} of their sites\n";
			print STDERR "\t\tWill only work with homer2 (i.e. use -homer2)\n" if ($cmd->{'homer2'} == 0);
		} elsif ($ARGV[$i] eq '-e' ) {
			$cmd->{'expect'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-bits' ) {
			$cmd->{'bits'} = "-bits";
		} elsif ($ARGV[$i] eq '-noweight' ) {
			$cmd->{'cgweight'} = 0; 
			print STDERR "\tWill not adjust sequences for CG content\n";
		} elsif ($ARGV[$i] eq '-tss' ) {
			$cmd->{'tssweight'} = 1; 
			$cmd->{'cgweight'} = 0; 
			print STDERR "\tWill normalize sequences for distance to TSS\n";
		} elsif ($ARGV[$i] eq '-cgtss' ) {
			$cmd->{'tssweight'} = 1; 
			print STDERR "\tWill normalize sequences for distance to TSS and CpG content\n";
		} elsif ($ARGV[$i] eq '-redundant' ) {
			$cmd->{'redundant'} = $ARGV[++$i]; 
			my $pp = 100* $cmd->{'redundant'};
			print STDERR "\tWill remove >$pp% redundant sequences\n";
		} elsif ($ARGV[$i] eq '-nomotif' ) {
			$cmd->{'nomotif'} = 1; 
			print STDERR "\tWill not run homer for de novo motifs\n";
		} elsif ($ARGV[$i] eq '-known' ) {
			$cmd->{'noknown'} = 0; 
		} elsif ($ARGV[$i] eq '-noknown' ) {
			$cmd->{'noknown'} = 1; 
			print STDERR "\tWill not search for known motifs\n";
		} elsif ($ARGV[$i] eq '-nomask' ) {
			$cmd->{'nomask'} = 1; 
			$cmd->{'mask'} = 0;
			print STDERR "\tUsing non-repeat masked sequences\n";
		} elsif ($ARGV[$i] eq '-rand' ) {
			$cmd->{'rand'} = 1;
			print STDERR "\tWill randomize target and background labels\n";
		} elsif ($ARGV[$i] eq '-mknown' ) {
			$cmd->{'mknown'} = $ARGV[++$i];
			print STDERR "\tWill search for known motifs in file: $cmd->{'mknown'}\n";
		} elsif ($ARGV[$i] eq '-mcheck' ) {
			$cmd->{'mcheck'} = $ARGV[++$i];
			print STDERR "\tWill search for known motifs (for de novo checking) in file: $cmd->{'mcheck'}\n";
		} elsif ($ARGV[$i] eq '-mis' ) {
			$cmd->{'mis'} = $ARGV[++$i];
			print STDERR "\tWill search for strings with $cmd->{'mis'} mismatches\n";
		} elsif ($ARGV[$i] eq '-S' ) {
			$cmd->{'S'} = $ARGV[++$i];
			print STDERR "\tWill optimize $cmd->{'S'} putative motifs\n";
		} elsif ($ARGV[$i] eq '-N' ) {
			$cmd->{'numSeq'} = $ARGV[++$i];
			print STDERR "\tTotal number of sequences (including background) to use: $cmd->{'numSeq'}\n";
		} elsif ($ARGV[$i] eq '-local' ) {
			$cmd->{'local'} = $ARGV[++$i];
			print STDERR "\tUsing local background ($cmd->{'local'} bp)\n";
		} elsif ($ARGV[$i] eq '-size' ) {
			$cmd->{'size'} = $ARGV[++$i];
			my $size = $cmd->{'size'};
			if ($size eq 'given') {
				print STDERR "\tUsing actual sizes of regions (-size given)\n";
			} elsif ($size =~ /\,/) {
				my @a = split /\,/, $size;
				my $sizeStart = $a[0];
				my $sizeEnd = $a[1];
				if ($sizeEnd <= $sizeStart) {
					print STDERR "!!! not using -size correctly, $sizeStart > $sizeEnd\n";
					exit;
				}
				print STDERR "\tUsing regions from $sizeStart to $sizeEnd relative to peak centers\n";
				$cmd->{'sizeMove'} = floor(($sizeStart+$sizeEnd)/2);
				$cmd->{'size'} = $sizeEnd-$sizeStart;
			}
			print STDERR "\tFragment size set to $cmd->{'size'}\n";
		} elsif ($ARGV[$i] eq '-nocheck' ) {
			$cmd->{'checkFlag'} = 0;
			print STDERR "\tWill not compare known motifs and de novo motifs\n";
		} elsif ($ARGV[$i] eq '-find' ) {
			$cmd->{'find'} = $ARGV[++$i];
			print STDERR "\tWill find motif(s) in $cmd->{'find'}\n";
		} else {
			print STDERR "!!! $ARGV[$i] not recognized!!\n";
			printCMD();
		}
	}

	if ($cmd->{'size'} eq 'NA') {
		print STDERR "!!! findMotifsGenome.pl now REQUIRES that you specify the -size option!!!\n";
		print STDERR "!!! Examples:\n";
		print STDERR "!!!   -size <#>       -size 200\n";
		print STDERR "!!!   -size <#>,<#>   -size -150,50\n";
		print STDERR "!!!   -size given     -size given (use actual regions in peak file)\n";
		print STDERR "!!! Try again with a valid option for -size... (this is to avoid confusion)\n";
		exit;
	}
	if ($cmd->{'find'}) {
		$cmd->{'maxN'} = 10;
	}


	if (@{$cmd->{'len'}} < 1) {
		push(@{$cmd->{'len'}},8);
		push(@{$cmd->{'len'}},10);
		push(@{$cmd->{'len'}},12);
	}

	return $cmd;

}	


my $tmpID = rand();

my $tmpFile = $cmd->{'output'} . "/" . $tmpID . ".tmp";
my $tmpFile2 = $cmd->{'output'} . "/" . $tmpID . ".2.tmp";
my $tmpFile3 = $cmd->{'output'} . "/" . $tmpID . ".3.tmp";
my $tmpFile4 = $cmd->{'output'} . "/" . $tmpID . ".4.tmp";
my $findFile = $cmd->{'output'} . "/" . $tmpID . ".find.tmp";
my $posFileSize = $cmd->{'output'} . "/target" . $cmd->{'size'} . ".pos";
my $targetSeq = $cmd->{'output'} . "/target" . $cmd->{'size'} . ".seq";
my $targetRedun = $cmd->{'output'} . "/target" . $cmd->{'size'} . ".redun";
my $targetNonRedun = $cmd->{'output'} . "/target" . $cmd->{'size'} . ".nonredun.seq";
my $targetCGFreq = $cmd->{'output'} . "/target" . $cmd->{'size'} . ".cgfreq";
my $targetCGBins = $cmd->{'output'} . "/target" . $cmd->{'size'} . ".cgbins";
my $targetTSSBins = $cmd->{'output'} . "/target" . $cmd->{'size'} . ".tssbins";
my $targetCGTSSBins = $cmd->{'output'} . "/target" . $cmd->{'size'} . ".cgtssbins";
my $knownMotifHTML = $cmd->{'output'} . "/knownResults.html";
my $groupCGAll = $cmd->{'output'} . "/group.cgbins.all";
my $groupCGbins = $cmd->{'output'} . "/group.cgbins";
my $groupFile = $cmd->{'output'} . "/group.adj";
my $seqFile = $cmd->{'output'} . "/seq.tsv";
my $noutFile = $cmd->{'output'} . "/seq.autonorm.tsv";
my $nooutFile = $cmd->{'output'} . "/oligo.autonorm.tsv";
my $oligoFilePrefix = $cmd->{'output'} . "/oligos";
my $localFile = $cmd->{'output'} . "/localBackground.txt";
my $motifInfoFile =  $cmd->{'output'} . "/" . $motifInfoFileName;

my $cleanPosFile = $cmd->{'output'} . "/target.clean..pos";
my $cleanBgFile = $cmd->{'output'} . "/bg.clean.pos";

my $dumpFastaTarget = $cmd->{'output'} . "/target.fa";
my $dumpFastaBackground = $cmd->{'output'} . "/background.fa";

my $mflag = "";
$mflag = " -mask " if ($cmd->{'mask'});

my $genomeDir = "";
my $preparsedDirFromConfig = "";
my $customGenome = "";
$cmd->{'org'} = 'null';
if (!exists($config->{'GENOMES'}->{$cmd->{'genome'}})) {
	$customGenome = $cmd->{'genome'};
	($cmd->{'genome'},$genomeDir,$preparsedDirFromConfig) = HomerConfig::parseCustomGenome($cmd->{'genome'});
} else {
	$genomeDir = $config->{'GENOMES'}->{$cmd->{'genome'}}->{'directory'} . "/";
	$cmd->{'org'} =  $config->{'GENOMES'}->{$cmd->{'genome'}}->{'org'};
	$preparsedDirFromConfig = $genomeDir . "preparsed/";
}
my $preparsedDir = $cmd->{'preparsedDir'};
if ($preparsedDir eq '/') {
	$preparsedDir = $preparsedDirFromConfig;
}

open IN, $cmd->{'posfile'} or die "!!! Could not open peak/position file $cmd->{'posfile'} !!!\n";
close IN;

#print STDERR "mset=$cmd->{'mset'}, org = $cmd->{'org'};\n";
($msetCheck, $msetKnown) = HomerConfig::checkMSet($cmd->{'mset'}, $cmd->{'org'});
$cmd->{'mcheck'} = $msetCheck if ($cmd->{'mcheck'} eq '');
$cmd->{'mknown'} = $msetKnown if ($cmd->{'mknown'} eq '');



`mkdir -p "$cmd->{'output'}"`;

open INFO, ">$motifInfoFile";
print INFO "cmd =";
foreach(@ARGV) {
	print INFO " $_";
}
print INFO "\n";
close INFO;

if ($cmd->{'refFlag'} == 1) {
	`cp "$cmd->{'posfile'}" "$cleanPosFile"`;
} else {
	`bed2pos.pl "$cmd->{'posfile'}" -check > "$cleanPosFile"`;
	`checkPeakFile.pl "$cleanPosFile"`;
	`cleanUpPeakFile.pl "$cleanPosFile" > "$tmpFile"`;
	`mv "$tmpFile" "$cleanPosFile"`;
}
$toDelete{$cleanPosFile}=1;
	


if ($cmd->{'bg'} ne '') {
	`bed2pos.pl "$cmd->{'bg'}" -check > "$cleanBgFile"`;
	$toDelete{$cleanBgFile}=1;

	if ($cmd->{'refFlag'} == 1) {
		`cleanUpPeakFile.pl "$cleanBgFile" > "$tmpFile"`;
	} else {
		if ($cmd->{'chopify'}) {
			if ($cmd->{'size'} eq 'given') {
				`cp "$cleanPosFile" "$tmpFile2"`;
			} else {
				`resizePosFile.pl "$cleanPosFile" $cmd->{'size'} > "$tmpFile2"`;
			}
			`chopUpPeakFile.pl "$tmpFile2" "$cleanBgFile" > "$tmpFile"`;
			`mv "$tmpFile" "$cleanBgFile"`;
			`rm "$tmpFile2"`;
		}
		`cleanUpPeakFile.pl "$cleanBgFile" BG > "$tmpFile"`;
	}
	`mv "$tmpFile" "$cleanBgFile"`;
}

$bestFragSize = $cmd->{'size'};

if ($cmd->{'size'} eq 'given' && $cmd->{'bg'} eq '') {
	my $N = 0;
	my $avg = 0;
	open IN, $cleanPosFile;
	while (<IN> ){
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next unless ($line[2] =~ /\d+/);
		next unless ($line[3] =~ /\d+/);
		$avg += abs($line[3]-$line[2]);
		$N++;
	}
	$avg /= $N if ($N > 0);
	$avg = floor($avg);
	$bestFragSize = $avg;
	print STDERR "\tBackground fragment size set to $bestFragSize (avg size of targets)\n";
}
		
my $bgCGBins = "";
my $bgSeq = "";

##############################################

if ($cmd->{'bg'} eq '') {
	$preparseFound = 0;
	my $prefix = $cmd->{'genome'};
	if ($cmd->{'mask'}) {
		$prefix .= "r";
	}
	`ls "$preparsedDir"/$prefix.*.cgbins > "$tmpFile"`;
	open IN, $tmpFile;
	my @availableSizes = ();
	while (<IN>) {
		chomp;
		/\.(\d+?)\.cgbins/;
		my $psize = $1;
		push(@availableSizes, $psize);
		if ($psize == $bestFragSize) {
			$preparseFound = 1;
		}
	}
	close IN;
	`rm "$tmpFile"`;
	
	if ($preparseFound == 1 && $cmd->{'preparseForce'} == 0) {
		print STDERR "\tBackground files for $bestFragSize bp fragments found.\n";
	} else {
		if ($cmd->{'preparseForce'} == 1) {
			print STDERR "\tGenome preparsing was FORCED.\n";
		} else {
			print STDERR "\tCould not find background files for $bestFragSize bp fragments\n";
			print STDERR "\tBelow are the sizes that are already available prepared.\n";
			@availableSizes = sort {$a <=> $b} @availableSizes;
			foreach(@availableSizes) {
				print STDERR "\t\t$_\n";
			}
			print STDERR "\tHOMER will now create background files for $bestFragSize bp fragments\n";
			print STDERR "\tTo CANCEL and rerun with a differet \"-size <#>\", hit <CTRL+C> now!\n";
			for (my $i=5;$i>=1;$i--) {
				print STDERR "\t\t$i\n";
				`sleep 1`;
			}
		}
		print STDERR "\tPreparsing genome for $bestFragSize bp fragments...(will probably take 1-5 min)\n"; 
		my $gg = $cmd->{'genome'};
		if ($customGenome ne '') {
			$gg = $customGenome;
		}
		# check if we can write to the preparsed directory
		if (-w "$preparsedDir") {
		} else {
			`mkdir -p "$preparsedDir"`;
			if (-w "$preparsedDir") {
			} else {
				print STDERR "!!! Warning: Looks like you do not have permission to write the preparsed\n";
				print STDERR "!!! genome files to the following directory:\n";
				print STDERR "!!!   $preparsedDir\n";
				print STDERR "!!! Consider one of the two options:\n";
				print STDERR "!!!   1.) Get your system admin to set the directory to be writeable to you and/or\n";
				print STDERR "!!!       your group.\n";
				print STDERR "!!!          -or-\n";
				print STDERR "!!!   2.) Use the \"-preparsedDir <directory>\" option to specify a directory that\n";
				print STDERR "!!!       you do have permission to write files to.\n";
				print STDERR "\n";
				exit;
			}
		}
	
		`preparseGenome.pl "$gg" $mflag -size $bestFragSize -preparsedDir "$preparsedDir"`;
	}

	#$bestFragSize = $cmd->{'size'};
		
	$bgCGBins = $preparsedDir . "/$prefix.$bestFragSize.cgbins";
	if ($cmd->{'gc'}==1) {
		$bgCGBins = $preparsedDir . "/$prefix.$bestFragSize.gcbins";
	}
	$bgSeq = $preparsedDir . "/$prefix.$bestFragSize.seq";

	my $wcOutput = `wc -l "$bgCGBins"`;
	$wcOutput =~ /^\s*(\d+)\s/;
	my $lineCount = $1;
	if ($lineCount < 1000) {
		print STDERR "!!!! Might have something wrong with preparsed files\n";
		print STDERR "!!!! Rerun and add \"-preparse\" to the command line to force recreation of the files\n";
		`sleep 10`;
	}
} else {
	# remove overlapping target and background positions
	`mergePeaks "$cleanBgFile" "$cleanPosFile" -d $cmd->{'size'} -cobound 1 -prefix "$tmpFile2"`;
	`mv "$tmpFile2".coBoundBy0.txt "$cleanBgFile"`;
}

if ($cmd->{'homer2'}==1) {
	if ($cmd->{'alg'} eq ' -alg binomial ') {
		$cmd->{'alg'} = " -stat binomial ";
	} elsif ($cmd->{'alg'} eq '') {
		$cmd->{'alg'} = " -stat hypergeo ";
	}
}


##############################################
`ls "$cmd->{'output'}"/* > $tmpFile`;
open IN, $tmpFile;
my %files = ();
while (<IN>) {
	chomp;
	$files{$_} = 1;
}
close IN;
`rm $tmpFile`;

my $executeFlag = 0;
#print STDERR "Checking for extracted target sequences...";
if (exists($files{$targetSeq}) && $cmd->{'force'} == 0) {
	#print STDERR "yes\n";
} else {
	#print STDERR "no\n";
	$executeFlag = 1;
}
if ($executeFlag==1) {
	if ($cmd->{'refFlag'} == 1) {
		if ($cmd->{'size'} ne 'given') {
			`resizePosFile.pl "$cleanBgFile" $cmd->{'size'} $cmd->{'sizeMove'} > "$posFileSize"`;
		} else {
			`cp "$cleanBgFile" "$posFileSize"`;
		}
		`homerTools extract "$posFileSize" "$genomeDir" $mflag > "$tmpFile"`;
		`cleanUpSequences.pl "$tmpFile" > "$tmpFile2"`;
		`removePoorSeq.pl "$tmpFile2" $cmd->{'maxN'} > "$targetSeq"`;
		`rm "$tmpFile" "$tmpFile2"`;
	} elsif ($cmd->{'bg'} ne '') {
		`cat "$cleanPosFile" "$cleanBgFile" > "$tmpFile2"`;
		if ($cmd->{'size'} ne 'given') {
			`resizePosFile.pl "$tmpFile2" $cmd->{'size'} $cmd->{'sizeMove'} > "$posFileSize"`;
		} else {
			`cp "$tmpFile2" "$posFileSize"`;
		}
		`homerTools extract "$posFileSize" "$genomeDir" $mflag > "$tmpFile"`;
		`cleanUpSequences.pl "$tmpFile" > "$tmpFile2"`;
		`removePoorSeq.pl "$tmpFile2" $cmd->{'maxN'} > "$targetSeq"`;
		`rm "$tmpFile" "$tmpFile2"`;
	} elsif ($cmd->{'local'} > 0) {
		if ($cmd->{'size'} ne 'given') {
			`resizePosFile.pl "$cleanPosFile" $cmd->{'size'} $cmd->{'sizeMove'}  > "$tmpFile"`;
		} else {
			`cp "$cleanPosFile" "$tmpFile"`;
		}
		my $local = $cmd->{'local'};
		open OUT, ">$posFileSize";
		open OUT2, ">$localFile";
		open IN, $tmpFile;
		while (<IN>) {
			chomp;
			s/\r//g;
			next if (/^#/);
			print OUT "$_\n";
			my @line = split /\t/;
			next unless ($line[2] =~ /^[\d]/);
			for (my $i=1;$i<=$local;$i++) {
				my $s = $line[2]+$bestFragSize*($i);
				my $e = $line[3]+$bestFragSize*($i);
				my $id = $line[0] . "-local+$i";
				print OUT "$id\t$line[1]\t$s\t$e\t$line[4]\n";
				print OUT2 "$id\t0\n";
				$s = $line[2]-$bestFragSize*($i);
				$e = $line[3]-$bestFragSize*($i);
				$id = $line[0] . "-local-$i";
				print OUT "$id\t$line[1]\t$s\t$e\t$line[4]\n";
				print OUT2 "$id\t0\n";
			}
		}
		close IN;
		close OUT;
		close OUT2;
		`homerTools extract "$posFileSize" "$genomeDir" $mflag > "$tmpFile"`;
		`cleanUpSequences.pl "$tmpFile" > "$tmpFile2"`;
		`removePoorSeq.pl "$tmpFile2" $cmd->{'maxN'} > "$targetSeq"`;
		`rm "$tmpFile" "$tmpFile2"`;
		$toDelete{$localFile}=1;
	} else {
		if ($cmd->{'size'} ne 'given') {
			`resizePosFile.pl "$cleanPosFile" $cmd->{'size'} $cmd->{'sizeMove'} > "$posFileSize"`;
		} else {
			`cp "$cleanPosFile" "$posFileSize"`;
		}
		`homerTools extract "$posFileSize" "$genomeDir" $mflag > "$tmpFile"`;
		`cleanUpSequences.pl "$tmpFile" > "$tmpFile2"`;
		`removePoorSeq.pl "$tmpFile2" $cmd->{'maxN'} > "$targetSeq"`;
		`rm "$tmpFile" "$tmpFile2"`;
	}
	$toDelete{$posFileSize}=1;
	$toDelete{$targetSeq}=1;
	print STDERR "\n";
}

if ($cmd->{'find'} eq '') {
	#print STDERR "Checking for non-redundant target sequences...";
	if (exists($files{$targetNonRedun}) && $cmd->{'force'} == 0) {
		#print STDERR "yes\n";
	} else {
		#print STDERR "no\n";
		$executeFlag = 1;
	}
	if ($executeFlag == 1) {
		if ($cmd->{'redundant'} > 1.99) { 
			print STDERR "\tNot removing redundant sequences\n";
			`cp "$targetSeq" "$targetNonRedun"`;
		} else {
			print STDERR "\tRemoving redundant sequences\n";
			`findRedundantBLAT.pl "$targetSeq" $cmd->{'redundant'} > "$targetRedun"`;
			if ($cmd->{'bg'} ne '') {
				`makeBinaryFile.pl "$cleanBgFile" "$cleanPosFile" > "$tmpFile"`;
			} elsif ($cmd->{'local'} > 0) {
				`makeBinaryFile.pl "$localFile" "$targetSeq" > "$tmpFile"`;
			} else {
				`makeBinaryFile.pl "$targetSeq" "$targetSeq" > "$tmpFile"`;
			}
			`adjustRedunGroupFile.pl "$tmpFile" "$targetRedun" | cut -f1 > "$tmpFile2"\n`;
			`addData.pl "$tmpFile2" "$targetSeq" > "$targetNonRedun"`;
			`rm "$tmpFile2" "$tmpFile"`;
			$toDelete{$targetRedun}=1;
		}
		$toDelete{$targetNonRedun}=1;
	}
	print STDERR "\n";
	
	#print STDERR "Checking for target CpG/TSS bins...";
	if (exists($files{$targetCGBins}) && $cmd->{'force'} == 0) {
		#print STDERR "yes\n";
	} else {
		#print STDERR "no\n";
		$executeFlag = 1;
	}
	if ($executeFlag == 1) {	
		`homerTools freq -gc "$targetCGFreq" "$targetNonRedun" > "$tmpFile"`;
		`rm "$tmpFile"`;
		my $col = 1;
		$col = 2 if ($cmd->{'gc'} == 1);
		`freq2group.pl "$targetCGFreq" $col > "$targetCGBins"`;
		#`annotatePeaks.pl "$posFileSize" $cmd->{'genome'} -noann > "$tmpFile"`;
		#`getTSSBins.pl "$tmpFile" > "$targetTSSBins"`;
		#`combineBinFiles.pl "$targetCGBins" "$targetTSSBins" > "$targetCGTSSBins"`;
		$toDelete{$targetCGFreq}=1;
		$toDelete{$targetCGBins}=1;
	}
	
	#print STDERR "Checking for group and sequence Homer files...";
	if (exists($files{$seqFile}) && exists($files{$groupFile}) && $cmd->{'force'} == 0) {
		#print STDERR "yes\n";
	} else {
		#print STDERR "no\n";
		$executeFlag = 1;
	}
	if ($executeFlag == 1) {
		if ($cmd->{'bg'} ne '' || $cmd->{'local'} > 0) {
			`cp "$targetNonRedun" "$seqFile"`;
			`makeBinaryFile.pl "$targetCGBins" "$cleanPosFile" > "$tmpFile"`;
			if ($cmd->{'redundant'} > 1.9) {
				`cp "$tmpFile" "$tmpFile2"`;
			}  else {
				`adjustRedunGroupFile.pl "$tmpFile" "$targetRedun" > "$tmpFile2"`;
			}
			
			if ($cmd->{'cgweight'} == 1 && $cmd->{'tssweight'} == 1) { 
				#`assignGeneWeights.pl "$tmpFile2" "$targetCGTSSBins" > "$groupFile"`;
			} elsif ($cmd->{'cgweight'} == 1) {
				`assignGeneWeights.pl "$tmpFile2" "$targetCGBins" > "$groupFile"`;
			} elsif ($cmd->{'tssweight'} == 1) {
				#`assignGeneWeights.pl "$tmpFile2" "$targetTSSBins" > "$groupFile"`;
			} else {
				`cp "$tmpFile2" "$groupFile"`;
			}

			`rm "$tmpFile" "$tmpFile2"`;
			#print STDERR "rm $tmpFile $tmpFile2\n";
		} else {
			
			my $wcOutput = `wc -l "$targetCGBins"`;
			$wcOutput =~ /^\s*(\d+)\s/;
			$numTargets = $1;
			if ($numTargets*2 > $cmd->{'numSeq'}) {
				$cmd->{'numSeq'} = $numTargets*2;
			}
			print STDERR "\n\tTotal sequences set to $cmd->{'numSeq'}\n";

			print STDERR "\n\tChoosing background that matches in CpG/GC Content...\n";
			if ($cmd->{'cgweight'} == 1 && $cmd->{'tssweight'} == 1) { 
				#`cat "$targetCGTSSBins" "$bgCGTSSBins" > "$tmpFile"`;
				#`makeBinaryFile.pl "$tmpFile" "$targetCGTSSBins" > "$tmpFile2"`;
				#`randRemoveBackground.pl "$tmpFile2" "$cmd->{'numSeq'}" "$tmpFile" > "$tmpFile3"`;
				#`assignGeneWeights.pl "$tmpFile3" "$tmpFile" > "$groupFile"`;
			} elsif ($cmd->{'cgweight'} == 1) { 
				`wc -l "$targetCGBins" "$bgCGBins"`;
				`cat "$targetCGBins" "$bgCGBins" > "$tmpFile"`;
				#`wc -l "$tmpFile" "$targetCGBins"`;
				`makeBinaryFile.pl "$tmpFile" "$targetCGBins" > "$tmpFile2"`;
				#`wc -l "$tmpFile2" "$tmpFile"`;
				`randRemoveBackground.pl "$tmpFile2" "$cmd->{'numSeq'}" "$tmpFile" > "$tmpFile3"`;
				#`wc -l "$tmpFile3" "$tmpFile"`;
				`assignGeneWeights.pl "$tmpFile3" "$tmpFile" > "$groupFile"`;
			} elsif ($cmd->{'tssweight'} == 1) { 
				#`cat "$targetTSSBins" "$bgTSSBins" > "$tmpFile"`;
				#`makeBinaryFile.pl "$tmpFile" "$targetTSSBins" > "$tmpFile2"`;
				#`randRemoveBackground.pl "$tmpFile2" "$cmd->{'numSeq'}" "$tmpFile" > "$tmpFile3"`;
				#`assignGeneWeights.pl "$tmpFile3" "$tmpFile" > "$groupFile"`;
			} else {
				print STDERR "\tRandomly choosing background...\n";
				`cat "$targetCGBins" "$bgCGBins" > "$tmpFile"`;
				`makeBinaryFile.pl "$tmpFile" "$targetCGBins" > "$tmpFile2"`;
				`randRemoveBackground.pl "$tmpFile2" "$cmd->{'numSeq'}" null > "$tmpFile3"`;
				`cp "$tmpFile3" "$groupFile"`;
			}
		
			print STDERR "\tAssembling sequence file...\n";
			`filterListBy.pl  "$bgSeq" "$groupFile" 0 1 > "$tmpFile"`;
			if ($cmd->{'dumpFasta'}==1) {
				`tab2fasta.pl "$tmpFile" > "$dumpFastaBackground"`;
				`tab2fasta.pl "$targetNonRedun" > "$dumpFastaTarget"`;
			}
			`cat "$targetNonRedun" "$tmpFile" > "$seqFile"`;
			`rm "$tmpFile" "$tmpFile2" "$tmpFile3"`;
		}
		$toDelete{$groupFile}=1;
		$toDelete{$seqFile}=1;
	}

	if ($cmd->{'rand'} == 1) {
		print STDERR "Randomizing Target and Background sequences\n";
		`randomizeGroupFile.pl "$groupFile" > "$tmpFile"`;
		`mv "$tmpFile" "$groupFile"`;
	}

	if ($cmd->{'homer2'} && $cmd->{'nlen'} > 0) {
		print STDERR "\tNormalizing lower order oligos using homer2\n";
		my $options = "";
		$options .= " -strand + " if ($cmd->{'norevopp'} == 1);
		$options .= " -nmax $cmd->{'nmax'} ";
		$options .= " -nlen $cmd->{'nlen'} ";
		$options .= " -nout \"$noutFile\" ";
		$options .= $cmd->{'neutral'};
		`homer2 norm -g "$groupFile" -s "$seqFile" $options  > "$tmpFile"`;
		`mv "$tmpFile" "$groupFile"`;
	}


	if (scalar(@{$cmd->{'motifMask'}}) > 0) {
		print STDERR "\tMasking given motifs...\n";
		my $files = '';
		foreach(@{$cmd->{'motifMask'}}) {
			$files .= " \"$_\"";
		}
		`cat $files > "$tmpFile"`;
		if ($cmd->{'homer2'}) {
			my $options = "";
			$options .= " -strand + " if ($cmd->{'norevopp'} == 1);
			`homer2 mask -s "$seqFile" -m "$tmpFile" $options > "$tmpFile2"`;
		} else {
			`homer -s "$seqFile" -m "$tmpFile" -a REMOVE > "$tmpFile2"`;
		}
		`mv "$tmpFile2" "$seqFile"`;
		`rm "$tmpFile"`;
	}

	print STDERR "\tFinished preparing sequence/group files\n\n";

	##########################################################
	print STDERR "\t----------------------------------------------------------\n";
	print STDERR "\tKnown motif enrichment\n";
	my $floatAction = "GETPVALUE";
	$floatAction = "OPTPVALUE" if ($cmd->{'float'} == 1);

	if ($cmd->{'noknown'}==1) { #|| ($cmd->{'force'} == 0 && exists($files{$knownMotifHTML}))) {
		#print STDERR "\tSkipping ($knownMotifHTML might exist - delete it if you want to repeat knownMotif analysis)...\n";
	} else {
		my $options = $cmd->{'bits'};
		$options = " -optimize" if ($floatAction eq 'OPTPVALUE');
		if ($cmd->{'homer2'}) {
			$options .= " -homer2 -p $cmd->{'cpus'}";
			$options .= " $cmd->{'alg'}";
			$options .= " -cache $cmd->{'cache'}";
		}
		#`findKnownMotifs.pl "$seqFile" "$groupFile" "$cmd->{'output'}" $knownPvalueThresh "$cmd->{'mknown'}" $floatAction`;
		`findKnownMotifs.pl -s "$seqFile" -g "$groupFile" -o "$cmd->{'output'}" -pvalue $knownPvalueThresh -m "$cmd->{'mknown'}" $options`;
	}

	if ($cmd->{'nomotif'} == 1) {
		print STDERR "\tSkipping...\n";
	} else {
		print STDERR "\t----------------------------------------------------------\n";
		print STDERR "\tDe novo motif finding (HOMER)\n";
		my $options = " -S $cmd->{'S'} ";
		my $coptions = $cmd->{'alg'} . " ";
		$options .= $cmd->{'alg'};
		$options .= " -mis $cmd->{'mis'}";
		my $cpuOptions = "";

		if ($cmd->{'homer2'}) {


			if (scalar(@{$cmd->{'motifOpt'}}) > 0) {
				print STDERR "\tOptimizing/Changing length of given motifs...\n";
				my $files = '';
				foreach(@{$cmd->{'motifOpt'}}) {
					$files .= " \"$_\"";
				}
				`cat $files > "$tmpFile"`;
				$options .= " -opt \"$tmpFile\" ";
			}

			if ($cmd->{'norevopp'} == 1) {
				$options .= ' -strand + ';
				$coptions .= ' -strand + ';
			}
			$cpuOptions = " -p $cmd->{'cpus'} ";

			$options .= " -cache $cmd->{'cache'}";
			$coptions .= " -cache $cmd->{'cache'}";
			if ($cmd->{'quickMask'} > 0) {
				$options .= " -quickMask";
			}
			if ($cmd->{'expect'} > 0) {
				$options .= " -e $cmd->{'expect'} ";
			}
			if ($cmd->{'olen'} > 0) {
				$options .= " -olen $cmd->{'olen'} ";
				$options .= " -omax $cmd->{'nmax'} ";
			}
			$options .= " -minlp $cmd->{'minlp'} ";
		} else {
			$options .= " -o \"$cmd->{'output'}/homerMotifs\" ";
			$options .= " -branch $cmd->{'depth'}";
			if ($cmd->{'cgweight'} == 1 || $cmd->{'tssweight'} == 1) {
				$options .= ' -w ';
			}
			if ($cmd->{'norevopp'} == 1) {
				$options .= ' -norevopp';
			}
		}

		foreach(@{$cmd->{'len'}}) {
			my $len = $_;
			if ($cmd->{'homer2'}) {

				if ($cmd->{'oligoFlag'} == 1) {
					`homer2 denovo -s "$seqFile" -g "$groupFile" $options $cpuOptions -len $_ -oligos "$oligoFilePrefix.$_.txt"`;
				}
				if ($cmd->{'nomotif'} == 0) {
					my $output = " -o \"$cmd->{'output'}/homerMotifs.motifs$len\" ";
					#print STDERR "`homer2 denovo -s $seqFile -g $groupFile $options $output -len $_`;\n";
					`homer2 denovo -s "$seqFile" -g "$groupFile" $options $cpuOptions $output -len $len`;

					if ($cmd->{'fdr'} > 0) {
						`mkdir -p \"$cmd->{'output'}/randomizations/\"`;
						print STDERR "\tPerforming empirical FDR calculation for length $len (n=$cmd->{'fdr'})\n";

						my $realPvalues = readPvalues("$cmd->{'output'}/homerMotifs.motifs$len");
						my @randPvalues = ();
						my @randFiles = ();

						my $cpus = 0;
						for (my $i=0;$i<$cmd->{'fdr'};$i++) {
							my $ii = $i+1;
							print STDERR "\t\t$ii of $cmd->{'fdr'}\n";
							my $outputfile = "$cmd->{'output'}/randomizations/homerMotifs.r$i.motifs$len";
							push(@randFiles,$outputfile);
							$pid = fork();
							$cpus++;
							if ($pid == 0) {
								#child process
								my $randGroupFile = "$cmd->{'output'}/randomizations/rand$i.group";
								`randomizeGroupFile.pl "$groupFile" > "$randGroupFile"`;
								my $output = " -o \"$outputfile\" ";
								my $cmdStr = "homer2 denovo -s \"$seqFile\" -g \"$randGroupFile\" $options -len $len $output";
								`$cmdStr 2> /dev/null`;
								`rm -f "$randGroupFile"`;
								exit(0);
							}
							if ($cpus >= $cmd->{'cpus'}) {
								wait();
								$cpus--;
							}
						}
						my $id = 0;
						while ($id >= 0) {
							$id = wait();
						}
						foreach(@randFiles) {
							my $file = $_;
							my $rPvalues = readPvalues($file);
							foreach(@$rPvalues) {
								push(@randPvalues, $_);
							}
						}
						my $fdrs = Statistics::empiricalFDR2($realPvalues,\@randPvalues, $cmd->{'fdr'});
						addFDR("$cmd->{'output'}/homerMotifs.motifs$len",$fdrs);
					}
				}

			} else {
				if ($cmd->{'oligoFlag'} == 1) {
					`homer -s "$seqFile" -g "$groupFile" $options -len $_ -a MERS > "$oligoFilePrefix.$_.txt"`;
				}
				if ($cmd->{'nomotif'} == 0) {
					#print STDERR "`homer -s $seqFile -g $groupFile $options -len $_ -a MOTIFS`;\n";
					`homer -s "$seqFile" -g "$groupFile" $options -len $_ -a MOTIFS `;
				}
			}
		}
		$toDelete{".tmp.motifs"}=1;

		if ($cmd->{'checkFlag'} == 1 && $cmd->{'nomotif'} == 0) {
			########visualization of homer... need to work on
			if ($cmd->{'homer2'}) {
				`cat "$cmd->{'output'}"/homerMotifs.motifs* > "$cmd->{'output'}/homerMotifs.all.motifs"`;
				#`homer2 known -s "$seqFile" -g "$groupFile" $coptions -m "$cmd->{'output'}/homerMotifs.all.approx.motifs" -siteReduce $cmd->{'percentSimilar'} -mout "$cmd->{'output'}/homerMotifs.all.motifs" > /dev/null`;
				#`rm "$cmd->{'output'}/homerMotifs.all.approx.motifs"`;
			} else {
				`cat "$cmd->{'output'}"/homerMotifs.motifs* > "$cmd->{'output'}/homerMotifs.all.motifs"`;
			}
			my $rnaFlag = "";
			$rnaFlag = " -rna -norevopp " if ($cmd->{'rnaMode'});
			`compareMotifs.pl "$cmd->{'output'}/homerMotifs.all.motifs" "$cmd->{'output'}/" -reduceThresh $reduceThresh -matchThresh $matchThresh -known "$cmd->{'mcheck'}" $cmd->{'bits'} $cmd->{'nofacts'} -cpu $cmd->{'cpus'} $rnaFlag`;
		}
		`rm -f "$tmpFile"`;
	}
	print STDERR "\tJob finished - if results look good, please send beer to ..\n\n";

} else {

	my $options = '';
	my $offset = 0;
	if ($cmd->{'size'} ne 'given') {
		$offset = -1*floor($cmd->{'size'}/2);
	}

	if ($cmd->{'homer2'}) {
		$options .= " -strand + " if ($cmd->{'norevopp'} == 1);
		`homer2 find -s "$targetSeq" $options -m "$cmd->{'find'}" -offset $offset > "$findFile"`;
		print "PositionID\tOffset\tSequence\tMotif Name\tStrand\tMotifScore\n";
	} else {
		if ($cmd->{'norevopp'} == 1) {
			$options .= " -norevopp";
		}
		`homer -s "$targetSeq" $options -a FIND -m "$cmd->{'find'}" -offset $offset > "$findFile"`;
		print "PositionID\tOffset\tSequence\tConservation\tStrand\tMotif Name\tMotifScore\n";
		#\tRefseq\tEnsembl\tName\tAlias\tOrf\tChr\tDescription\n";
	}

	open IN, $findFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		for (my $i=0;$i<@line;$i++) {
			#next if ($i==2 || $i==4);
			print "\t" if ($i > 0);
			print "$line[$i]";
		}
		print "\n";
	}
	close IN;
	`rm "$findFile"`;
}

if ($cmd->{'keepFiles'} == 0) {
	deleteFiles(\%toDelete);
}
print STDERR "\n";
exit;

sub deleteFiles {
	my ($files) = @_;
	print STDERR "\tCleaning up tmp files...\n";
	foreach(keys %$files) {
		my $file = $_;
		`rm -f \"$file\"`;
	}
}
sub readPvalues {
	my ($motifFile) = @_;
	open IN, $motifFile;
	my @values = ();
	while (<IN>) {
		chomp;
		if (/^>/) {
			my @line = split /\t/;
			push(@values, $line[3]);
		}
	}
	close IN;
	return \@values;
}
sub addFDR {
	my ($mfile, $fdrs) = @_;
	my $index = 0;
	my $tmp = rand() . ".tmp";
	open OUT, ">$tmp";
	open IN, $mfile;
	while (<IN>) {
		my $og = $_;
		if (/^>/) {
			chomp;
			my @line = split /\t/;
			print OUT "$line[0]";
			for (my $i=1;$i<@line;$i++) {
				print OUT "\t$line[$i]";
				if ($i==5) {
					print OUT ",FDR:$fdrs->[$index]";
				}
			}
			print OUT "\n";
			$index++;
		} else {
			print OUT $og;
		}
	}
	close IN;
	close OUT;
	`mv "$tmp" "$mfile"`;
}
