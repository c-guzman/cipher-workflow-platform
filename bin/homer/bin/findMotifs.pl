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

#default options...
my $accDir = $homeDir . "/data/accession/";
my $seqDir = $homeDir . "/data/promoters/";
my $reduceThresh = 0.6;
my $matchThresh = "T10";
my $knownPvalueThresh = 0.01;
my $motifInfoFileName = "motifFindingParameters.txt";
my $config = HomerConfig::loadConfigFile();
my $idtype = "";
our %toDelete = ();


if (@ARGV < 3) {
	printCMD();
}

sub printCMD {

	print STDERR "\n\tProgram will find de novo and known motifs in a gene list\n\n";
	print STDERR "\t\tUsage:  findMotifs.pl <input list> <promoter set> <output directory> [additoinal options]\n";
	print STDERR "\n\t\texample: findMotifs.pl genelist.txt mouse motifResults/ -len 10\n";
	print STDERR "\n\t\tFASTA example: findMotifs.pl targets.fa fasta motifResults/ -fasta background.fa\n";

	print STDERR "\n\tAvailable Promoter Sets: Add custom promoters sets with loadPromoters.pl\n";
	my $z = 0;
	foreach (keys %{$config->{'PROMOTERS'}}) {
		print STDERR "\t\t$_\t$config->{'PROMOTERS'}->{$_}->{'org'}\t$config->{'PROMOTERS'}->{$_}->{'directory'}"
						. "\t$config->{'PROMOTERS'}->{$_}->{'start'}\t$config->{'PROMOTERS'}->{$_}->{'end'}"
						. "\t$config->{'PROMOTERS'}->{$_}->{'idtype'}\n";
		$z++;
	}
	print STDERR "\n\t\tTry typing \"perl $homeDir/configureHomer.pl -list\" to see available promoter sets\n";
	print STDERR "\t\tTyping \"perl $homeDir/configureHomer.pl -install NNN\" to install promoter set NNN\n";

	print STDERR "\n\tBasic options:\n";
	print STDERR "\t\t-len <#>[,<#>,<#>...] (motif length, default=8,10,12) [NOTE: values greater 12 may cause the program\n";
	print STDERR "\t\t\tto run out of memmory - in these cases decrease the number of sequences analyzed]\n";
	print STDERR "\t\t-bg <background file> (ids to use as background, default: all genes)\n";
	print STDERR "\t\t-start <#> (offset from TSS, default=-300) [max=based on Promoter Set]\n";	
	print STDERR "\t\t-end <#> (offset from TSS, default=50) [max=based on Promoter Set]\n";	
	print STDERR "\t\t-rna (output RNA motif logos and compare to RNA motif database, automatically sets -norevopp)\n";
	print STDERR "\t\t-mask/-nomask (use/don't use repeatmasked files, default: -mask)\n";
	print STDERR "\t\t-S <#> (Number of motifs to optimize, default: 25)\n";
	print STDERR "\t\t-mis <#> (global optimization: searches for strings with # mismatches, default: 1)\n";
	print STDERR "\t\t-noconvert (will not worry about converting input files into unigene ids)\n";
	print STDERR "\t\t-norevopp (do not search the reverse strand for motifs)\n";
	print STDERR "\t\t-nomotif (don't search for de novo motif enrichment)\n";

	print STDERR "\n\tScanning sequence for motifs\n";
	print STDERR "\t\t-find <motif file> (This will cause the program to only scan for motifs)\n";
	print STDERR "\n\tIncluding Enhancers - peak files of enhancer location, peak ID should be gene ID\n";
	print STDERR "\t\t-enhancers <peak file> <genome verion>\n";
	print STDERR "\t\t\t(enhancers to include in search space, peaks/sequences should be named with a gene ID\n";
	print STDERR "\t\t\tIf multiple enhancers per gene, use the same gene ID, and all will be included)\n";
	print STDERR "\t\t-enhancersOnly (do not include promoter sequence in motif search)\n";

	print STDERR "\n\tFASTA files: If you prefer to use your own fasta files, place target sequences and \n";
	print STDERR "\t\tbackground sequences in two separate FASTA formated files (must have unique identifiers)\n";
	print STDERR "\t\tTarget File - use in place of <input list> (i.e. the first argument)\n";
	print STDERR "\t\tBackground File - after output directory (with additional options) use the argument:\n";
	print STDERR "\t\t\t-fastaBg <background fasta file> (This is recommended for fasta based analysis)\n";
	print STDERR "\t\tIn place of the promoter set use \"fasta\", or any valid set (this parameter is ignored)\n";
	print STDERR "\t\tWhen finding motifs [-find], only the target file with be searched)\n";
	print STDERR "\t\t\t-chopify (chops up background regions to match size of target regions)\n";
	print STDERR "\t\t\t\ti.e. if background is a full genome or all mRNAs\n";

	print STDERR "\n\tKnown Motif Options/Visualization:\n";
	print STDERR "\t\t-mset <vertebrates|insects|worms|plants|yeast|all> (check against motif collects, default: auto)\n";
	print STDERR "\t\t-basic (don't check de novo motifs for similarity to known motifs)\n";
	print STDERR "\t\t-bits (scale sequence logos by information content, default: doesn't scale)\n";
	print STDERR "\t\t-nocheck (don't check for similarity between novo motif motifs and known motifs)\n";
	print STDERR "\t\t-mcheck <motif file> (known motifs to check against de novo motifs,\n";
	print STDERR "\t\t-noknown (don't search for known motif enrichment, default: -known)\n";
	print STDERR "\t\t-mknown <motif file> (known motifs to check for enrichment,\n";
	print STDERR "\t\t-nofacts (omit humor)\n";

	print STDERR "\n\tAdvanced options:\n";
	print STDERR "\t\t-b (use binomial distribution to calculate p-values, hypergeometric is default)\n";
	print STDERR "\t\t-nogo (don't search for gene ontology enrichment)\n";
	print STDERR "\t\t-humanGO (Convert IDs to human for GO analysis)\n";
	print STDERR "\t\t-noweight (no CG correction)\n";
	print STDERR "\t\t-noredun (Don't remove predetermined redundant promoters/sequences)\n";
	print STDERR "\t\t-g (input file is a group file, i.e. 1st column = id, 2nd = 0 or 1 [1=target,0=back])\n";
	print STDERR "\t\t-cpg (use CpG% instead of GC% for sequence normalization)\n";
	print STDERR "\t\t-rand (randomize labels for target and backgound sequences)\n";
	print STDERR "\t\t-maskMotif <motif file 1> [motif file 2] ... (motifs to mask before motif finding)\n";
	print STDERR "\t\t-opt <motif file 1> [motif file 2] ... (motifs to optimize/change length)\n";
	print STDERR "\t\t-peaks (will produce peak file of promoters to use with findMotifsGenome.pl)\n";
	print STDERR "\t\t-nowarn (no warnings)\n";
	print STDERR "\t\t-keepFiles (don't delete temporary files)\n";
	print STDERR "\t\t-dumpFasta (create target.fa and background.fa files)\n";
	print STDERR "\t\t-min <#> (remove sequences shorter than #, default: 0)\n";
	print STDERR "\t\t-max <#> (remove sequences longer than #, default: 1e10)\n";
	print STDERR "\t\t-reuse (rerun homer using old seq files etc. with new options\n";
	print STDERR "\t\t\t  and ignores input list, organism)\n";
	print STDERR "\t\t-fdr <#> (Calculate empirical FDR for de novo discovery #=number of randomizations)\n";

	print STDERR "\n";
	print STDERR "\thomer2 specific options:\n";
	print STDERR "\t\t-homer2 (use homer2 instead of original homer, default)\n";
	print STDERR "\t\t-nlen <#> (length of lower-order oligos to normalize - general sequences, default: 3)\n";
	print STDERR "\t\t\t-nmax <#> (Max normalization iterations, default: 160)\n";
    print STDERR "\t\t\t-neutral (weight sequences to neutral frequencies, i.e. 25%, 6.25%, etc.)\n";
	print STDERR "\t\t-olen <#> (lower-order oligo normalization for oligo table, use if -nlen isn't working well)\n";
	print STDERR "\t\t-p <#> (Number of processors to use, default: 1)\n";
	print STDERR "\t\t-e <#> (Maximum expected motif instance per bp in random sequence, default: 0.01)\n";
	print STDERR "\t\t-cache <#> (size in MB for statistics cache, default: 500)\n";
	print STDERR "\t\t-quickMask (skip full masking after finding motifs, similar to original homer)\n";
	print STDERR "\t\t-homer1 (to force the use of the original homer)\n";
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
	my $cmd = {fg=>$ARGV[0], bg=>'', 
				promoters=>$ARGV[1], output=>$ARGV[2], 
				start=>-300, end=>50, rand=>0,maxlen=>10, gc=>1,
				noconv=>0,noweight=>0,redundant=>1.5, nogo=>0,motifMask=>\@motifFiles,motifOpt=>\@motifFiles2,
				nomotif=>0,noknown=>0, groupFlag=>0, norevopp=>0,float=>0,
				motif=>'', nomask=>0, len=>\@len, mis=>1,fasta=>'',mask=>1,
				S=>25, reuse=>0, g=>0, find=>'', alg=>'', peaks=>0,
				mcheck=>'', mknown=>'', checkFlag=>1,depth=>0.5,
				nowarn=>0,keepFiles=>0,homer2=>1,nlen=>3,olen=>0,expect=>0,cpus=>1,
				percentSimilar=>0.20,cache=>500,bits=>"",quickMask=>0,rnaMode=>0,
				minLen=>0,maxLen=>1e10,chopify=>0,nofacts=>"",fdr=>0,nmax=>160,neutral=>"",
				enhancersOnly=>0,enhancers=>"",genome=>"",dumpFasta=>0,minlp=>-10,mset=>'auto',
				humanGO=>0
				};
	print STDERR "\nSelected Options:\n";
	print STDERR "\tInput file = $cmd->{'fg'}\n";
	print STDERR "\tPromoter Set = $cmd->{'promoters'}\n";
	print STDERR "\tOutput Directory = $cmd->{'output'}\n";

	if ($cmd->{'promoters'} =~ /\-mRNA/) {
		print STDERR "\tRunning in mRNA mode (-norevopp -min 200 -max 10000 -noknown)\n";
		$cmd->{'rnaMode'} = 1;
		$cmd->{'noknown'} = 1;
		$cmd->{'norevopp'} = 1;
		$cmd->{'minLen'} = 200;
		$cmd->{'maxLen'} = 10000;
	}
	
	for (my $i=3;$i<@$argv;$i++) { 
		if ($ARGV[$i] eq '-bg') {
			$cmd->{'bg'} = $ARGV[++$i];
			print STDERR "\tbackground file: $cmd->{'bg'}\n";
		} elsif ($ARGV[$i] eq '-len' ) {
			my @lengths = split /\,/, $ARGV[++$i];
			$cmd->{'len'} = \@lengths;
			my $str = '';
			print STDERR "\tMotif length set at ";
			$cmd->{'maxlen'} = 0;
			foreach(@lengths) {
				print STDERR "$_, ";
				$cmd->{'maxlen'} = $_ if ($_ > $cmd->{'maxlen'});
				if ($_ > 12) {
					print STDERR "*** might want to hit ctrl+C to quit if it takes too long!\n";
					print STDERR "*** If you run out of memmory try reducing the number background sequences\n";
				}
			}
			print STDERR "\n";
		} elsif ($ARGV[$i] eq '-enhancersOnly' ) {
			$cmd->{'enhancersOnly'} = 1;
		} elsif ($ARGV[$i] eq '-enhancers' ) {
			if ($i+2 >= @ARGV ) {
				print STDERR "Something's wrong with your command line argument: -enhancers <peak/fasta file> <genome/fasta>\n";
				exit(1);
			}
			$cmd->{'enhancers'} = $ARGV[++$i];
			$cmd->{'genome'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-neutral' ) {
			$cmd->{'neutral'} = " -neutral";
		} elsif ($ARGV[$i] eq '-humanGO' ) {
			$cmd->{'humanGO'} = 1;
		} elsif ($ARGV[$i] eq '-minlp' ) {
			$cmd->{'minlp'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-mset' ) {
			$cmd->{'mset'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-dumpFasta' ) {
			$cmd->{'dumpFasta'} = 1;
		} elsif ($ARGV[$i] eq '-nmax' ) {
			$cmd->{'nmax'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-fdr' ) {
			$cmd->{'fdr'} = $ARGV[++$i];
			print STDERR "\tWill randomize and repeat motif finding $cmd->{'fdr'} times to estimate FDR\n";
		} elsif ($ARGV[$i] eq '-g' ) {
			$cmd->{'groupFlag'} = 1; 
			print STDERR "\t$ARGV[0] will be treated as a group file\n";
		} elsif ($ARGV[$i] eq '-basic' ) {
			$cmd->{'nofacts'} = " -basic ";
		} elsif ($ARGV[$i] eq '-nofacts' ) {
			$cmd->{'nofacts'} = " -nofacts ";
		} elsif ($ARGV[$i] eq '-keepFiles' ) {
			$cmd->{'keepFiles'} = 1; 
			print STDERR "\tWill keep temporary files\n";
		} elsif ($ARGV[$i] eq '-chopify' ) {
			$cmd->{'chopify'} = 1; 
			print STDERR "\tWill chopify the background FASTA files into the right size chunks\n";
		} elsif ($ARGV[$i] eq '-norevopp' ) {
			$cmd->{'norevopp'} = 1; 
			print STDERR "\tWill not search the reverse strand\n";
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
		} elsif ($ARGV[$i] eq '-reuse' ) {
			$cmd->{'reuse'} = 1; 
			print STDERR "\tOld files will be used to repeat homer analysis\n";
		} elsif ($ARGV[$i] eq '-rna' ) {
			print STDERR "\tModifying analysis for RNA analysis\n";
			$cmd->{'rnaMode'} = 1;
			$cmd->{'norevopp'} = 1;
			$cmd->{'noknown'} = 1;
			$cmd->{'mknown'} = $homeDir . "/data/knownTFs/known.rna.motifs";
			$cmd->{'mcheck'} = $homeDir . "/data/knownTFs/all.rna.motifs";
		} elsif ($ARGV[$i] eq '-opt') {
			my $bail = 0;
			print STDERR "\tWill optimize motifs in the following files:\n";
			print STDERR "\t\tSkipping known enrichment check\n";
			$cmd->{'noknown'}=1;
			while ($ARGV[++$i] !~ /^\-/) {
				print STDERR "\t\t$ARGV[$i]\n";
				push(@{$cmd->{'motifOpt'}}, $ARGV[$i]);
				if ($i>=@ARGV-1) {
					$bail=1;
					last;
				}
			}
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
		} elsif ($ARGV[$i] eq '-cpg' || $ARGV[$i] eq '-CpG') {
			$cmd->{'gc'} = 0;
			print STDERR "\tUsing GC% instead of CpG%\n";
		} elsif ($ARGV[$i] eq '-noweight' ) {
			$cmd->{'noweight'} = 1; 
			print STDERR "\tWill not adjust sequences for CG content\n";
		} elsif ($ARGV[$i] eq '-noredun' ) {
			$cmd->{'redundant'} = 2;
			#print STDERR "\tRemoving redundant sequences > $cmd->{'redundant'} % similar\n";
		} elsif ($ARGV[$i] eq '-nomotif' ) {
			$cmd->{'nomotif'} = 1; 
			print STDERR "\tWill not run homer for de novo motifs\n";
		} elsif ($ARGV[$i] eq '-peaks' ) {
			$cmd->{'peaks'} = 1; 
			print STDERR "\tWill create peak file of target promoters (will print to stdout)\n";
		} elsif ($ARGV[$i] eq '-known' ) {
			$cmd->{'noknown'} = 0;
		} elsif ($ARGV[$i] eq '-noknown' ) {
			$cmd->{'noknown'} = 1; 
			print STDERR "\tWill not search for known motifs\n";
		} elsif ($ARGV[$i] eq '-bits' ) {
			$cmd->{'bits'} = "-bits";
		} elsif ($ARGV[$i] eq '-homer1' ) {
			$cmd->{'homer2'} = 0;
			print STDERR "\tUsing original homer\n";
		} elsif ($ARGV[$i] eq '-homer2' ) {
			$cmd->{'homer2'} = 1;
			print STDERR "\tUsing homer2 (warning, unpredictably awesome results may, or may not, ensue)\n";
		} elsif ($ARGV[$i] eq '-cache' ) {
			$cmd->{'cache'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-ps' ) {
			$cmd->{'percentSimilar'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-quickMask' ) {
			$cmd->{'quickMask'} = 1;
		} elsif ($ARGV[$i] eq '-p' || $ARGV[$i] eq '-cpu') {
			$cmd->{'cpus'} = $ARGV[++$i];
			print STDERR "\tUsing $cmd->{'cpus'} CPUs\n";
			print STDERR "\t\tWill only work with homer2 (i.e. use -homer2)\n" if ($cmd->{'homer2'} == 0);
		} elsif ($ARGV[$i] eq '-olen' ) {
			$cmd->{'olen'} = $ARGV[++$i];
			print STDERR "\tWill normalize background oligos for oligos of length $cmd->{'olen'} bp\n";
			print STDERR "\t\tWill only work with homer2 (i.e. use -homer2)\n" if ($cmd->{'homer2'} == 0);
		} elsif ($ARGV[$i] eq '-nlen' ) {
			$cmd->{'nlen'} = $ARGV[++$i];
			print STDERR "\tWill normalize background sequences for oligos of length $cmd->{'nlen'} bp\n";
			print STDERR "\t\tWill only work with homer2 (i.e. use -homer2)\n" if ($cmd->{'homer2'} == 0);
		} elsif ($ARGV[$i] eq '-e' ) {
			$cmd->{'expect'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-min' ) {
			$cmd->{'minLen'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-max' ) {
			$cmd->{'maxLen'} = $ARGV[++$i];
		} elsif ($ARGV[$i] eq '-nowarn' ) {
			$cmd->{'nowarn'} = 1; 
			print STDERR "\tWill not stop program on warnings\n";
		} elsif ($ARGV[$i] eq '-rand' ) {
			$cmd->{'rand'} = 1; 
			print STDERR "\tWill randomize sequence labels (i.e. find motifs in random dataset)\n";
		} elsif ($ARGV[$i] eq '-nogo' ) {
			$cmd->{'nogo'} = 1; 
			print STDERR "\tWill Skip Gene Ontology Analysis\n";
		} elsif ($ARGV[$i] eq '-mask' ) {
			$cmd->{'nomask'} = 0; 
			$cmd->{'mask'} = 1;
			print STDERR "\tUsing repeat masked sequences\n";
		} elsif ($ARGV[$i] eq '-nomask' ) {
			$cmd->{'nomask'} = 1; 
			print STDERR "\tUsing non-repeat masked sequences\n";
		} elsif ($ARGV[$i] eq '-noconvert' ) {
			$cmd->{'noconv'} = 1; 
			print STDERR "\tWill not convert ids (Input files better be Unigene...)\n";
		} elsif ($ARGV[$i] eq '-b' ) {
			$cmd->{'alg'} = ' -alg binomial '; 
			print STDERR "\tUsing Binomial Distribution for p-values (instead of hypergeometric)\n";
		} elsif ($ARGV[$i] eq '-mis' ) {
			$cmd->{'mis'} = $ARGV[++$i];
			print STDERR "\tWill search for strings with $cmd->{'mis'} mismatches\n";
		} elsif ($ARGV[$i] eq '-S' ) {
			$cmd->{'S'} = $ARGV[++$i];
			print STDERR "\tWill optimize $cmd->{'S'} putative motifs\n";
		} elsif ($ARGV[$i] eq '-fasta' || $ARGV[$i] eq '-fastaBg') {
			$cmd->{'fasta'} = $ARGV[++$i];
			print STDERR "\tWill use FASTA files for motif finding\n";
			print STDERR "\t\tTarget Sequences = $cmd->{'fg'}\n";
			print STDERR "\t\tBackground Sequences = $cmd->{'fasta'}\n";
		} elsif ($ARGV[$i] eq '-start' ) {
			$cmd->{'start'} = $ARGV[++$i];
			print STDERR "\tNew start is $cmd->{'start'} relative to the TSS\n";
		} elsif ($ARGV[$i] eq '-find' ) {
			$cmd->{'find'} = $ARGV[++$i];
			print STDERR "\tWill find motif(s) in $cmd->{'find'}\n";
		} elsif ($ARGV[$i] eq '-end' ) {
			$cmd->{'end'} = $ARGV[++$i];
			print STDERR "\tNew end is $cmd->{'end'} relative to the TSS\n";
		} elsif ($ARGV[$i] eq '-nocheck' ) {
			$cmd->{'checkFlag'} = 0;
			print STDERR "\tWill not check for similarity between de novo and known motifs\n";
		} elsif ($ARGV[$i] eq '-mknown' ) {
			$cmd->{'mknown'} = $ARGV[++$i];
			print STDERR "\tKnown motif file set to $cmd->{'mknown'} (for known motif enrichment)\n";
		} elsif ($ARGV[$i] eq '-mcheck' ) {
			$cmd->{'mcheck'} = $ARGV[++$i];
			print STDERR "\tKnown motif file set to $cmd->{'mcheck'} (for checking de novo motifs)\n";
		} else {
			print STDERR "!! $ARGV[$i] is not recognized option!\n";
			printCMD();
		}
	}	
	if ($cmd->{'genome'} =~ s/r$//) {
		$cmd->{'mask'} = 1;
	}

	if (@{$cmd->{'len'}} < 1) {
		push(@{$cmd->{'len'}}, 8);
		push(@{$cmd->{'len'}}, 10);
		push(@{$cmd->{'len'}}, 12);
	}

	return $cmd;

}

if ($cmd->{'homer2'}==1) {
	if ($cmd->{'alg'} eq ' -alg binomial ') {
		$cmd->{'alg'} = " -stat binomial ";
	} elsif ($cmd->{'alg'} eq '') {
		$cmd->{'alg'} = " -stat hypergeo ";
	}
}


my $tmpID = rand();
my $ugFgFile = $cmd->{'output'} . '/' . "targetIDs.ug.txt";
my $targetGroupFile = $cmd->{'output'} . '/' . "targetIDs.group.txt";
my $ugBgFile = $cmd->{'output'} . '/' . "backgroundIDs.ug.txt";
my $ugGroupFile = $cmd->{'output'} . '/' . "group.ug.txt";
my $randGroupFile = $cmd->{'output'} . '/' . "group.rand.txt";
my $redunFile = $cmd->{'output'} . '/' . $tmpID . ".redun.tmp";
my $findFile = $cmd->{'output'} . '/' . $tmpID . ".find.tmp";
my $findFile2 = $cmd->{'output'} . '/' . $tmpID . ".find2.tmp";
my $findFile3 = $cmd->{'output'} . '/' . $tmpID . ".find3.tmp";
my $tmpSeq1 = $cmd->{'output'} . '/' . $tmpID . ".seq1.tmp";
my $tmpSeq2 = $cmd->{'output'} . '/' . $tmpID . ".seq2.tmp";
my $tmpSeq3 = $cmd->{'output'} . '/' . $tmpID . ".seq3.tmp";
my $tmpRedunDatFile = $cmd->{'output'}  . '/' . $tmpID . "redun.tmp";
my $tmpCGFreq = $cmd->{'output'}  . '/' . $tmpID . "cgfreq.tmp";
my $tmpCGBins = $cmd->{'output'}  . '/' . $tmpID . "cgbins.tmp";
my $adjFile = $cmd->{'output'} . "/group.adj";
my $seqFile = $cmd->{'output'} . "/seq.tsv";
my $noutFile = $cmd->{'output'} . "/seq.autonorm.tsv";
my $knownFile = $cmd->{'output'} . "/knownResults.html";
my $motifInfoFile =  $cmd->{'output'} . "/" . $motifInfoFileName;
my $enhancerSeq = $cmd->{'output'} . '/' . "enhancer.sequence.txt";
my $enhancerConvFile = $cmd->{'output'} . '/' . "enhancer.conv.txt";
my $enhancerCGBins = $cmd->{'output'} . '/' . "enhancer.cgbins";
my $targetFastaDump = $cmd->{'output'} . '/' . "target.fa";
my $backgroundFastaDump = $cmd->{'output'} . '/' . "background.fa";
my $tmpFile1 = $cmd->{'output'} . '/' . "$tmpID.1.tmp";
my $tmpFile2 = $cmd->{'output'} . '/' . "$tmpID.2.tmp";
my $tmpFile3 = $cmd->{'output'} . '/' . "$tmpID.3.tmp";
my $scrambledFasta = $cmd->{'output'} . '/' . "scrambleBg.fasta";


my $numLines = checkFile($cmd->{'fg'});
if ($numLines < 1) {
	print STDERR "There is no data in your input file ($cmd->{'fg'})\n";
	exit;
}
`mkdir -p "$cmd->{'output'}"`;

# check validity of Promoter Set
if ($cmd->{'promoters'} eq 'FASTA' || $cmd->{'promoters'} eq 'fasta') {
	if ($cmd->{'fasta'} eq '') {
		if ($cmd->{'find'}) {
			$cmd->{'fasta'} = 'placeholder';
		} else {
			print STDERR "\n\t!Warning - no background FASTA file specified (Highly recommended)\n";
			print STDERR "\t!Your input sequences will be randomized to serve as a background instead.\n\n";
			`scrambleFasta.pl "$cmd->{'fg'}" > "$scrambledFasta"`;
			$cmd->{'fasta'} = $scrambledFasta;
		}
	}
} elsif (!exists($config->{'PROMOTERS'}->{$cmd->{'promoters'}}) && $cmd->{'fasta'} eq '') {
	print STDERR "\n!!! $cmd->{'promoters'} not found in $homeDir/config.txt\n";
	print STDERR "\tTry typing \"perl $homeDir/configureHomer.pl -list\" to see available promoter sets\n";
	print STDERR "\tIf avaliable, type \"perl $homeDir/configureHomer.pl -install $cmd->{'promoters'}\" to install\n";
	exit;
}
my $promoterSeqOffset = 0;
if (exists($config->{'PROMOTERS'}->{$cmd->{'promoters'}})) {
	$cmd->{'org'} = $config->{'PROMOTERS'}->{$cmd->{'promoters'}}->{'org'};
	$seqDir = $config->{'PROMOTERS'}->{$cmd->{'promoters'}}->{'directory'} . "/";
	$maxTSSDist = $config->{'PROMOTERS'}->{$cmd->{'promoters'}}->{'end'};
	$idtype = $config->{'PROMOTERS'}->{$cmd->{'promoters'}}->{'idtype'};
	$promoterSeqOffset = $config->{'PROMOTERS'}->{$cmd->{'promoters'}}->{'start'};
} else {
	$cmd->{'org'} = 'null';
}
my $customGenome = "";
if ($cmd->{'enhancers'} ne '') {
	if (!exists($config->{'GENOMES'}->{$cmd->{'genome'}})) {
		$customGenome = $cmd->{'genome'};
		($cmd->{'genome'},$genomeDir,$preparsedDir) = HomerConfig::parseCustomGenome($cmd->{'genome'});
	} else {
		$genomeDir = $config->{'GENOMES'}->{$cmd->{'genome'}}->{'directory'} . "/";
		$preparsedDir = $genomeDir . "preparsed/";
	}
}

#check mset
my ($msetCheck, $msetKnown) = HomerConfig::checkMSet($cmd->{'mset'}, $cmd->{'org'});
$cmd->{'mcheck'} = $msetCheck if ($cmd->{'mcheck'} eq '');
$cmd->{'mknown'} = $msetKnown if ($cmd->{'mknown'} eq '');



open INFO, ">$motifInfoFile";
print INFO "cmd =";
foreach(@ARGV) {
	print INFO " $_";
}
print INFO "\n";
close INFO;

if ($cmd->{'enhancers'} ne '') {
	if ($cmd->{'genome'} ne 'fasta' && $cmd->{'genome'} ne 'FASTA') {
		my $mflag = "";
		if ($cmd->{'mask'}) {
			$mflag = " -mask ";
		}
		`bed2pos.pl "$cmd->{'enhancers'}" -check > "$tmpSeq1"`;
		`checkPeakFile.pl "$tmpSeq1"`;
		`cleanUpPeakFile.pl "$tmpSeq1" > "$tmpSeq2"`;
		`homerTools extract "$tmpSeq2" "$genomeDir" $mflag > "$tmpSeq1"`;
		`cleanUpSequences.pl "$tmpSeq1" > "$tmpSeq2"`;
		`removePoorSeq.pl "$tmpSeq2" > "$enhancerSeq"`;
		`homerTools freq -gc "$tmpSeq2" "$enhancerSeq" > "$tmpSeq1"`;
		my $col = 1;
		$col = 2 if ($cmd->{'gc'} == 1);
		`freq2group.pl "$tmpSeq2" $col > "$enhancerCGBins"`;
		
		`cut -f1 "$enhancerCGBins" > "$tmpSeq1"`; 	
		if ($cmd->{'noconv'} == 1 || $idtype eq 'null' || $idtype eq 'custom') {
			`duplicateCol.pl "$tmpSeq1" 0 | cut -f1,2 > "$enhancerConvFile"`;
		} else {
			`convertIDs.pl "$tmpSeq1" "$cmd->{'org'}" $idtype no yes | cut -f1,2 > "$enhancerConvFile"`;
		}
		`rm -f "$tmpSeq1" "$tmpSeq2"`;
		$toDelete{$enhancerSeq} = 1;
		$toDelete{$enhancerConvFile} = 1;
		$toDelete{$enhancerCGBins} = 1;
	}
}

###########################################################
# if input is group file, separate it into targets and background
if ($cmd->{'groupFlag'}) {
	`getPos.pl "$cmd->{'fg'}" > "$targetGroupFile"`;
	$cmd->{'bg'} = $cmd->{'fg'};
	$cmd->{'fg'} = $targetGroupFile;
	$numLines = checkFile($cmd->{'fg'});
}

if ($cmd->{'reuse'} == 0 && $cmd->{'fasta'} eq '') {
	print STDERR "\n\tProgress: Step1 - Convert input file to $idtype IDs\n";
	if ($cmd->{'noconv'} == 1) {
		print STDERR "\tskipping - file already converted\n";
		`cp "$cmd->{'fg'}" "$ugFgFile"`;
		if ($cmd->{'bg'} ne '') {
			`cp "$cmd->{'bg'}" "$ugBgFile"`;
			$toDelete{$ugBgFile}=1;
		}
	} else {
		if ($idtype eq 'null' || $idtype eq 'custom') {
			`cp "$cmd->{'fg'}" "$ugFgFile"`;
			if ($cmd->{'bg'} ne '') {
				`cp "$cmd->{'bg'}" "$ugBgFile"`;
				$toDelete{$ugBgFile}=1;
			}
		} else {
			#print STDERR "`convertIDs.pl $cmd->{'fg'} $cmd->{'org'} $idtype > $ugFgFile`;\n";
			`convertIDs.pl "$cmd->{'fg'}" "$cmd->{'org'}" $idtype > "$ugFgFile"`;
			if ($cmd->{'bg'} ne '') {
				`convertIDs.pl "$cmd->{'bg'}" "$cmd->{'org'}" $idtype > "$ugBgFile"`;
				$toDelete{$ugBgFile}=1;
			}
		}
	}
	$toDelete{$ugFgFile}=1;

	my $newNumLines = checkFile($ugFgFile);
	my $convPercent = $newNumLines/$numLines;
	my $zzz = sprintf("%.1f",$convPercent*100);
	print STDERR "\tPercentage of IDs converted into $idtype: $zzz"."% ($newNumLines out of $numLines)\n";
	if ($convPercent < 0.05) {
		print STDERR "!!!! Homer converted less than 5% of the rows in your file.\n";
		print STDERR "!!!! Check to be sure the input file has valid gene identifiers\n";
		print STDERR "!!!! Check to be sure the input file has valid gene identifiers\n";
		cleanUpAndExit() if ($cmd->{'nowarn'} == 0);
	}

	if ($cmd->{'peaks'} == 1) {
		my $promoterFile = $seqDir . $cmd->{'promoters'} . ".pos";
		my %ids = ();
		open IN, $ugFgFile;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			$ids{$line[0]} =1;
		}
		close IN;
		open IN, $promoterFile;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line= split /\t/;
			next if (!exists($ids{$line[0]}));
			my $start = $line[2]-$promoterSeqOffset-150;
			my $end = $line[2]-$promoterSeqOffset+50;
			if ($line[4] == 1) {
				$start = $line[3]+$promoterSeqOffset-50;
				$end = $line[3]+$promoterSeqOffset+150;
			}
			print "$line[0]\t$line[1]\t$start\t$end\t$line[4]\t$line[5]\n";
		}
		cleanUpAndExit();
	}

	############################################################
	print STDERR "\n\tProgress: Step2 - prepare sequence files\n";

	my $startFile = $seqDir . $cmd->{'promoters'} . ".mask";
	if ($cmd->{'nomask'}) {
		$startFile = $seqDir . $cmd->{'promoters'} . ".seq";
	}
	if ($cmd->{'rnaMode'} == 0) {
		my $start = $cmd->{'start'};
		my $end = $cmd->{'end'};
		if ($start < -1*$maxTSSDist || $start > $maxTSSDist) {
			print STDERR "Sequence start = $start is out of range (+/-$maxTSSDist)\n";
			exit;
		} 
		if ($end < -1*$maxTSSDist || $end > $maxTSSDist) {
			print STDERR "Sequence end = $end is out of range (+/-$maxTSSDist)\n";
			exit;
		}
		if ($end < $start || $end-$start < $cmd->{'maxlen'}) {
			print STDERR "Start and End values are too close!\n"; 
			exit;
		}
		`getPartOfPromoter.pl "$startFile" $start $end $promoterSeqOffset > "$tmpSeq1"`;
		`cleanUpSequences.pl "$tmpSeq1" -min $cmd->{'minLen'} -max $cmd->{'maxLen'} > "$seqFile"`;
		`rm "$tmpSeq1"`;
	} else {
		`cleanUpSequences.pl "$startFile" -min $cmd->{'minLen'} -max $cmd->{'maxLen'} > "$seqFile"`;
	}

	$toDelete{$seqFile}=1;

} elsif ($cmd->{'fasta'} ne '') {

	#############################################################
	# prepare and parse fasta files
	print STDERR "\tParsing FASTA format files...\n";
	`fasta2tab.pl "$cmd->{'fg'}" > "$tmpSeq1"`;
	`cleanUpSequences.pl "$tmpSeq1" -min $cmd->{'minLen'} -max $cmd->{'maxLen'} > "$tmpSeq2"`;
	`cleanUpPeakFile.pl "$tmpSeq2" > "$tmpSeq1"`;
	#`mv "$tmpSeq2" "$tmpSeq1"`;
	if ($cmd->{'find'} ne '') {
		my $options = '';

		if ($cmd->{'homer2'}) {
			if ($cmd->{'norevopp'} == 1) {
				$options .= " -strand +";
			}
			`homer2 find -s "$tmpSeq1" $options -m "$cmd->{'find'}" > "$findFile"`;
			print "FASTA ID\tOffset\tSequence\tMotif Name\tStrand\tMotifScore\n";
		} else {
			if ($cmd->{'norevopp'} == 1) {
				$options .= " -norevopp";
			}
			`homer -s "$tmpSeq1" $options -a FIND -m "$cmd->{'find'}" > "$findFile"`;
			print "FASTA ID\tOffset\tSequence\tConservation\tStrand\tMotif Name\tMotifScore\n";
		}

		open IN, $findFile;
		while (<IN>) {
			print $_;
		}
		close IN;
		`rm "$tmpSeq1" "$findFile"`;
		exit;
	}
	`fasta2tab.pl "$cmd->{'fasta'}" > "$tmpSeq2"`;
	if ($cmd->{'chopify'}==1) {
		`chopUpBackground.pl "$tmpSeq1" "$tmpSeq2" > "$tmpSeq3"`;
		`mv "$tmpSeq3" "$tmpSeq2"`;
	}
	`cleanUpSequences.pl "$tmpSeq2" > "$tmpSeq3"`;
	`cleanUpPeakFile.pl "$tmpSeq3" BG > "$tmpSeq2"`;

	`cat "$tmpSeq1" "$tmpSeq2" > "$seqFile"`;
	`makeBinaryFile.pl "$seqFile" "$tmpSeq1" > "$ugGroupFile"`;
	`rm -f "$tmpSeq1" "$tmpSeq2" "$tmpSeq3"`;

	$toDelete{$seqFile}=1;
	$toDelete{$ugGroupFile}=1;
}

#################################################
if ($cmd->{'find'} eq '' && $cmd->{'reuse'} == 0 && $cmd->{'fasta'} eq '') {
	print STDERR "\n\tProgress: Step3 - creating foreground/background file\n";
	if ($cmd->{'g'} == 1) {
		`cp "$ugFgFile" "$ugGroupFile"`;
	} else {
		if ($cmd->{'bg'} ne '') {
			`makeBinaryFile.pl "$ugBgFile" "$ugFgFile" > "$ugGroupFile"`;
		} else {
			my $baseFile = $seqDir .  $cmd->{'promoters'} . '.base';
			`makeBinaryFile.pl "$baseFile" "$ugFgFile" > "$ugGroupFile"`;
		} 
	}
	$toDelete{$ugGroupFile}=1;
}
##############################################################
if ($cmd->{'rand'} ne '0') {
	print STDERR "Randomizing Group File!!\n";
	`randomizeGroupFile.pl "$ugGroupFile" > "$randGroupFile"`;
	`mv "$randGroupFile" "$ugGroupFile"`;
}

if ($cmd->{'find'} eq '' && $cmd->{'reuse'} == 0) {
	###########################################################
	print STDERR "\n\tProgress: Step4 - removing redundant promoters\n";
	if ($cmd->{'redundant'} > 1.99) {
		print STDERR "\tskipping...\n";
		`cp "$ugGroupFile" "$redunFile"`;
	} else {
		if ($cmd->{'fasta'} ne '') {
			`cp "$ugGroupFile" "$redunFile"`;
			# need to calculate redundant sequences
			#`findRedundantBLAT.pl "$seqFile" $cmd->{'redundant'} > "$tmpRedunDatFile"`;
			#`adjustRedunGroupFile.pl "$ugGroupFile" "$tmpRedunDatFile" > "$redunFile"`;
			#`rm "$tmpRedunDatFile"`;
		} else {
			my $redunDatFile = $seqDir . $cmd->{'promoters'} . ".redun";
			#print STDERR "`adjustRedunGroupFile.pl $ugGroupFile $redunDatFile > $redunFile`\n";
			`adjustRedunGroupFile.pl "$ugGroupFile" "$redunDatFile" > "$redunFile"`;
		}
	}
	$toDelete{$redunFile}=1;
	if ($cmd->{'dumpFasta'} == 1) {
		print STDERR "\tDumping FASTA files of target and background sequences...\n";
		`getPos.pl "$redunFile" > "$tmpFile1"`;
		`filterListBy.pl "$seqFile" "$tmpFile1" 0 1 > "$tmpFile2"`;
		`tab2fasta.pl "$tmpFile2" > "$targetFastaDump"`;

		`getPos.pl "$redunFile" 1 > "$tmpFile1"`;
		`filterListBy.pl "$seqFile" "$tmpFile1" 0 1 > "$tmpFile2"`;
		`tab2fasta.pl "$tmpFile2" > "$backgroundFastaDump"`;

		`rm "$tmpFile1" "$tmpFile2"`;
	}


	if ($cmd->{'enhancers'} ne '') {
		`cat "$seqFile" "$enhancerSeq" > "$tmpSeq1"`;
		`mv "$tmpSeq1" "$seqFile"`;
		my %groups = ();
		open IN, "$ugGroupFile";
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			$groups{$line[0]} = $line[1];
		}
		close IN;
		if ($cmd->{'enhancersOnly'} ) {
			open REDUN, ">$redunFile";
		} else {
			open REDUN, ">>$redunFile";
		}
		open IN, "$enhancerConvFile";
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			if (exists($groups{$line[0]})) {
				print REDUN "$line[1]\t$groups{$line[0]}\n";
			}
		}
		close IN;
	}

	############################################################
	print STDERR "\n\tProgress: Step5 - adjusting background sequences for GC/CpG content...\n";
	if ($cmd->{'noweight'} == 1) {
		print STDERR "\tskipping...\n";
		`cp "$redunFile" "$adjFile"`;
	} else {
		if ($cmd->{'rnaMode'} == 0) {
			if ($cmd->{'fasta'} eq '') {
				my $weightFile = $seqDir . $cmd->{'promoters'} . ".cgbins";
				if ($cmd->{'gc'} == 1) {
					$weightFile = $seqDir . $cmd->{'promoters'} . ".gcbins";
				}
				if ($cmd->{'enhancers'} ne '') {
					`cat "$weightFile" "$enhancerCGBins" > "$tmpSeq1"`;
					`assignGeneWeights.pl "$redunFile" "$tmpSeq1" > "$adjFile"`;
					`rm "$tmpSeq1";`
				} else {
					`assignGeneWeights.pl "$redunFile" "$weightFile" > "$adjFile"`;
				}
			} else {
				`homerTools freq "$seqFile" -gc "$tmpCGFreq" > /dev/null`;
				my $col = 1;
				if ($cmd->{'gc'} == 1) {
					$col = 2;
				}
				`freq2group.pl "$tmpCGFreq" $col > "$tmpCGBins"`;
				`assignGeneWeights.pl "$redunFile" "$tmpCGBins" > "$adjFile"`;
				`rm "$tmpCGFreq" "$tmpCGBins"`;
			}
		} else {
			`cp "$redunFile" "$adjFile"`;
		}

		if ($cmd->{'nlen'} > 0) {
			print STDERR "\n\tNormalizing lower order oligos using homer2\n";
			my $options = "";
			$options .= " -strand + " if ($cmd->{'norevopp'} == 1);
			$options .= " -nmax $cmd->{'nmax'}";
			$options .= $cmd->{'neutral'};
			#print STDERR "`homer2 norm -g $adjFile -s $seqFile -nlen $cmd->{'nlen'} $options -nout $noutFile > $tmpCGFreq`;\n";
			`homer2 norm -g "$adjFile" -s "$seqFile" -nlen $cmd->{'nlen'} $options -nout "$noutFile" > "$tmpCGFreq"`;
			`mv "$tmpCGFreq" "$adjFile"`;
		}
	}
	$toDelete{$adjFile}=1;
}
if ($cmd->{'find'} eq '') {

	##########################################################
	print STDERR "\n\tProgress: Step6 - Gene Ontology Enrichment Analysis\n";
	if ($cmd->{'nogo'} == 1 || $cmd->{'fasta'} ne '') {
		print STDERR "\tSkipping...\n";
	} else {
		my $options = " -cpu $cmd->{'cpus'}";
		if ($cmd->{'bg'} ne '') {
			$options .= " -bg \"$cmd->{'bg'}\"";
		}
		if ($cmd->{'humanGO'} == 1) {
			$options .= " -human ";
		}
			
		`findGO.pl "$cmd->{'fg'}" "$cmd->{'org'}" "$cmd->{'output'}" $options`;
	}

	####################################################
	if (scalar(@{$cmd->{'motifMask'}}) > 0) {
		print STDERR "Masking given motifs...\n";
		my $files = '';
		foreach(@{$cmd->{'motifMask'}}) {
			$files .= " \"$_\"";
		}
		`cat $files > "$tmpSeq1"`;

		if ($cmd->{'homer2'}) {
			my $options = "";
			$options .= " -strand + " if ($cmd->{'norevopp'} == 1);
			`homer2 mask -s "$seqFile" -m "$tmpSeq1" $options > "$tmpSeq2"`;
		} else {
			`homer -s "$seqFile" -m "$tmpSeq1" -a REMOVE > "$tmpSeq2"`;
		}
		`mv "$tmpSeq2" "$seqFile"`;
		`rm "$tmpSeq1"`;
	}

	##########################################################
	print STDERR "\n\tProgress: Step7 - Known motif enrichment\n";
	if ($cmd->{'noknown'} == 1 || $cmd->{'reuse'}==1 ) {
		print STDERR "\tSkipping...\n";
	} else {
		my $floatAction = "GETPVALUE";
		$floatAction = "OPTPVALUE" if ($cmd->{'float'} == 1);
		my $options = "";
		$options = " -optimize" if ($floatAction eq 'OPTPVALUE');
		if ($cmd->{'homer2'}) {
			$options .= " -homer2 -p $cmd->{'cpus'}";
			$options .= " $cmd->{'alg'}";
		}
		`findKnownMotifs.pl -s "$seqFile" -g "$adjFile" -o "$cmd->{'output'}" -pvalue $knownPvalueThresh -m "$cmd->{'mknown'}" $options`;
	}

	##########################################################
	print STDERR "\n\tProgress: Step8 - De novo motif finding (HOMER)\n";
	if ($cmd->{'nomotif'} == 1) {
		print STDERR "\tSkipping...\n";
	} else {
		my $options = " -S $cmd->{'S'} ";
		$options .= $cmd->{'alg'};
		$options .= " -mis $cmd->{'mis'}";
		my $coptions = $cmd->{'alg'};
		my $cpuOptions = "";

		if ($cmd->{'homer2'}) {

			if (scalar(@{$cmd->{'motifOpt'}}) > 0) {
				print STDERR "\tOptimizing given motifs...\n";
				my $files = '';
				foreach(@{$cmd->{'motifOpt'}}) {
					$files .= " \"$_\"";
				}
				`cat $files > "$tmpSeq1"`;
				$options .= " -opt \"$tmpSeq1\" ";
			}

			if ($cmd->{'norevopp'} == 1) {
				$options .= " -strand + ";
				$coptions .= " -strand + ";
			}
			$cpuOptions .= " -p $cmd->{'cpus'} ";

			$options .= " -cache $cmd->{'cache'} ";
			$coptions .= " -cache $cmd->{'cache'} ";
			if ($cmd->{'quickMask'} > 0) {
				$options .= " -quickMask";
			}
			if ($cmd->{'expect'} > 0) {
				$options .= " -e $cmd->{'expect'} ";
			}
			$options .= " -minlp $cmd->{'minlp'}  ";
			$options .= " -o $cmd->{'olen'} " if ($cmd->{'olen'} > 0);
			foreach(@{$cmd->{'len'}}) {
				my $len = $_;
				my $outfile .= " -o \"$cmd->{'output'}/homerMotifs.motifs$len\" ";
				`homer2 denovo -s "$seqFile" -g "$adjFile" $options $cpuOptions -len $len $outfile`;

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
							#child proces
							my $randGroupFile = "$cmd->{'output'}/randomizations/rand$i.group";
							`randomizeGroupFile.pl "$adjFile" > "$randGroupFile"`;
							my $output = " -o \"$outputfile\" ";
							my $cmdStr = "homer2 denovo -s \"$seqFile\" -g \"$randGroupFile\" $options -len $len";
							$cmdStr .= $output;
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
			if ($cmd->{'noweight'} == 0) {
				$options .= ' -w ';
			}
			if ($cmd->{'norevopp'} == 1) {
				$options .= ' -norevopp';
			}
			$options .= " -o $cmd->{'output'}/homerMotifs ";
			$options .= " -branch $cmd->{'depth'}";
			foreach(@{$cmd->{'len'}}) {
				`homer -s "$seqFile" -g "$adjFile" $options -len $_ -a MOTIFS`;
			}
			$toDelete{".tmp.motifs"}=1;
			$toDelete{".mer.motifs"}=1;
		}


		if ($cmd->{'checkFlag'} == 1) {
			my $outdir = $cmd->{'output'};
			$outdir =~ s/ /\\ /g;
			if ($cmd->{'homer2'}) {
				`cat "$outdir/homerMotifs.motifs"* > "$cmd->{'output'}/homerMotifs.all.motifs"`;
				#`homer2 known -s "$seqFile" -g "$adjFile" $coptions -m "$cmd->{'output'}/homerMotifs.all.approx.motifs" -siteReduce $cmd->{'percentSimilar'} -mout "$cmd->{'output'}/homerMotifs.all.motifs" -offset $cmd->{'start'} > /dev/null`;
			} else {
				`cat "$outdir/homerMotifs.motifs"* > "$cmd->{'output'}/homerMotifs.all.motifs"`;
			}
			my $rnaOpt = "";
			$rnaOpt = " -rna " if ($cmd->{'rnaMode'});
			#print STDERR "`compareMotifs.pl $cmd->{'output'}/homerMotifs.all.motifs $cmd->{'output'}/ -reduceThresh $reduceThresh -matchThresh $matchThresh -known $cmd->{'mcheck'} BITS: $cmd->{'bits'} FACTS: $cmd->{'nofacts'} $rnaOpt `;\n";
			`compareMotifs.pl "$cmd->{'output'}/homerMotifs.all.motifs" "$cmd->{'output'}/" -reduceThresh $reduceThresh -matchThresh $matchThresh -known "$cmd->{'mcheck'}" $cmd->{'bits'} $cmd->{'nofacts'} -cpu $cmd->{'cpus'} $rnaOpt `;
		}
	}
	

	print STDERR "\tJob finished\n\n";

} else {
	my $options = '';
	`makeBinaryFile.pl "$ugFgFile" "$ugFgFile" > "$findFile"`;
	`mv "$findFile" "$ugFgFile"`;
	if ($cmd->{'homer2'}) {
		$options .= " -strand +" if ($cmd->{'norevopp'} == 1);
		`homer2 find -s "$seqFile" -g "$ugFgFile" $options -m "$cmd->{'find'}" -offset $cmd->{'start'} > "$findFile"`;
		print "GeneID\tPromoterID\tOffset\tSequence\tMotif Name\tStrand\tMotifScore\tUnigene\tRefseq\tEnsembl\tName\tAlias\tOrf\tChr\tDescription\tType\n";
	} else {
		$options .= " -norevopp" if ($cmd->{'norevopp'} == 1);
		`homer -s "$seqFile" -g "$ugFgFile" $options -a FIND -m "$cmd->{'find'}" -offset $cmd->{'start'} > "$findFile"`;
		print "GeneID\tPromoterID\tOffset\tSequence\tConservation\tStrand\tMotif Name\tMotifScore\tUnigene\tRefseq\tEnsembl\tName\tAlias\tOrf\tChr\tDescription\tType\n";
	}
	`convertIDs.pl "$findFile" $cmd->{'org'} gene no yes > "$findFile2"`;
	`addData.pl "$findFile2" "$accDir/$cmd->{'org'}.description" > "$findFile"`;
	open IN, $findFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		for (my $i=0;$i<@line;$i++) {
			print "\t" if ($i > 0);
			print "$line[$i]";
		}
		print "\n";
	}
	close IN;
	`rm "$findFile2" "$findFile"`;
}
cleanUpAndExit();
exit;

sub checkFile {
	my ($file) = @_;
	open IN, $file or die "!!!!\nCould not open file $file\n!!!!\n";
	my $c = 0;
	while (<IN>) {
		$c++;
	}
	return $c;
}

sub cleanUpAndExit {
	if ($cmd->{'keepFiles'} == 0) {
		foreach(keys %toDelete) {
			`rm "$_"`;
		}
	}
	exit;
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

