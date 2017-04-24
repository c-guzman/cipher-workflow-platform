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
use Statistics;


my $config = HomerConfig::loadConfigFile();

sub printCMD {

	print STDERR "\n\tUsage: annotatePeaks.pl <peak file | tss> <genome version>  [additional options...]\n";
	print STDERR "\n\tAvailable Genomes (required argument): (name,org,directory,default promoter set)\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}"
					. "\t$config->{'GENOMES'}->{$_}->{'promoters'}\n";
	}
	print STDERR "\t\t\t-- or --\n";
    print STDERR "\t\tCustom: provide the path to genome FASTA files (directory or single file)\n";
    print STDERR "\t\tIf no genome is available, specify 'none'.\n";
	print STDERR "\t\tIf using FASTA file or none, may want to specify '-organism <...>'\n";

	print STDERR "\n\tUser defined annotation files (default is UCSC refGene annotation):\n";
	print STDERR "\t\tannotatePeaks.pl accepts GTF (gene transfer formatted) files to annotate positions relative\n";
	print STDERR "\t\tto custom annotations, such as those from de novo transcript discovery or Gencode.\n";
	print STDERR "\t\t-gtf <gtf format file> (Use -gff and -gff3 if appropriate, but GTF is better)\n";
	print STDERR "\t\t-gid (by default the GTF file is processed by transcript_id, use this option for gene_id)\n";
	print STDERR "\t\t-ann <custom homer annotation file> (created by assignGenomeAnnotation, see website)\n";
	print STDERR "\n\tPeak vs. tss/tts/rna mode (works with custom GTF file):\n";
	print STDERR "\t\tIf the first argument is \"tss\" (i.e. annotatePeaks.pl tss hg18 ...) then a TSS centric\n";
	print STDERR "\t\tanalysis will be carried out.  Tag counts and motifs will be found relative to the TSS.\n";
	print STDERR "\t\t(no position file needed) [\"tts\" now works too - e.g. 3' end of gene]\n";
	print STDERR "\t\t[\"rna\" specifies gene bodies, will automaticall set \"-size given\"]\n";
	print STDERR "\t\tNOTE: The default TSS peak size is 4000 bp, i.e. +/- 2kb (change with -size option)\n";
	print STDERR "\t\t-list <gene id list> (subset of genes to perform analysis [unigene, gene id, accession,\n";
	print STDERR "\t\t\t probe, etc.], default = all promoters)\n";

	#print STDERR "\t\t-TSS <promoter set> (promoter definitions, default=genome default)\n";
	print STDERR "\t\t-cTSS <promoter position file i.e. peak file> (should be centered on TSS)\n";

	#print STDERR "\n\tAvailable Promoter Sets (use with -TSS): (name,org,directory,genome,masked genome)\n";
	#foreach(keys %{$config->{'PROMOTERS'}}) {
	#	print STDERR "\t\t$_\t$config->{'PROMOTERS'}->{$_}->{'org'}\t$config->{'PROMOTERS'}->{$_}->{'directory'}"
	#				. "\t$config->{'PROMOTERS'}->{$_}->{'genome'}\t$config->{'PROMOTERS'}->{$_}->{'mgenome'}\n";
	#}
	
	print STDERR "\n\tPrimary Annotation Options:\n";
    print STDERR "\t\t-mask (Masked repeats, can also add 'r' to end of genome name)\n";
	print STDERR "\t\t-m <motif file 1> [motif file 2] ... (list of motifs to find in peaks)\n";
	print STDERR "\t\t\t-mscore (reports the highest log-odds score within the peak)\n";
	print STDERR "\t\t\t-nmotifs (reports the number of motifs per peak)\n";
	print STDERR "\t\t\t-mdist (reports distance to closest motif)\n";
	print STDERR "\t\t\t-mfasta <filename> (reports sites in a fasta file - for building new motifs)\n";
	print STDERR "\t\t\t-fm <motif file 1> [motif file 2] (list of motifs to filter from above)\n";
	print STDERR "\t\t\t-rmrevopp <#> (only count sites found within <#> on both strands once, i.e. palindromic)\n";
	print STDERR "\t\t\t-matrix <prefix> (outputs a motif co-occurrence files:\n";
	print STDERR "\t\t\t\tprefix.count.matrix.txt - number of peaks with motif co-occurrence\n";
	print STDERR "\t\t\t\tprefix.ratio.matrix.txt - ratio of observed vs. expected  co-occurrence\n";
	print STDERR "\t\t\t\tprefix.logPvalue.matrix.txt - co-occurrence enrichment\n";
	print STDERR "\t\t\t\tprefix.stats.txt - table of pair-wise motif co-occurrence statistics\n";
	print STDERR "\t\t\t\tadditional options:\n";
	print STDERR "\t\t\t\t-matrixMinDist <#> (minimum distance between motif pairs - to avoid overlap, default: 4)\n";
	print STDERR "\t\t\t\t-matrixMaxDist <#> (maximum distance between motif pairs)\n";
	print STDERR "\t\t\t-mbed <filename> (Output motif positions to a BED file to load at UCSC (or -mpeak))\n";
	print STDERR "\t\t\t-mlogic <filename> (will output stats on common motif orientations)\n";
	print STDERR "\t\t-d <tag directory 1> [tag directory 2] ... (list of experiment directories to show\n";
	print STDERR "\t\t\ttag counts for) NOTE: -dfile <file> where file is a list of directories in first column\n";
	print STDERR "\t\t-bedGraph <bedGraph file 1> [bedGraph file 2] ... (read coverage counts from bedGraph files)\n";
	print STDERR "\t\t-wig <wiggle file 1> [wiggle file 2] ... (read coverage counts from wiggle files)\n";
	print STDERR "\t\t-p <peak file> [peak file 2] ... (to find nearest peaks)\n";
	print STDERR "\t\t\t-pdist to report only distance (-pdist2 gives directional distance)\n";
	print STDERR "\t\t\t-pcount to report number of peaks within region\n";
	print STDERR "\t\t-vcf <VCF file> (annotate peaks with genetic variation infomation, one col per individual)\n";
	print STDERR "\t\t\t-editDistance (Computes the # bp changes relative to reference)\n";
	print STDERR "\t\t\t-individuals <name1> [name2] ... (restrict analysis to these individuals)\n";
	print STDERR "\t\t\t-editDistance (Computes the # bp changes relative to reference)\n";
	print STDERR "\t\t\t-individuals <name1> [name2] ... (restrict analysis to these individuals)\n";
	print STDERR "\t\t-gene <data file> ... (Adds additional data to result based on the closest gene.\n";
	print STDERR "\t\t\tThis is useful for adding gene expression data.  The file must have a header,\n";
	print STDERR "\t\t\tand the first column must be a GeneID, Accession number, etc.  If the peak\n";
	print STDERR "\t\t\tcannot be mapped to data in the file then the entry will be left empty.\n";
	print STDERR "\t\t-go <output directory> (perform GO analysis using genes near peaks)\n";
	print STDERR "\t\t-genomeOntology <output directory> (perform genomeOntology analysis on peaks)\n";
	print STDERR "\t\t\t-gsize <#> (Genome size for genomeOntology analysis, default: 2e9)\n";

	print STDERR "\n\tAnnotation vs. Histogram mode:\n";
	print STDERR "\t\t-hist <bin size in bp> (i.e 1, 2, 5, 10, 20, 50, 100 etc.)\n";
	print STDERR "\t\tThe -hist option can be used to generate histograms of position dependent features relative\n";
	print STDERR "\t\tto the center of peaks.  This is primarily meant to be used with -d and -m options to map\n";
	print STDERR "\t\tdistribution of motifs and ChIP-Seq tags.  For ChIP-Seq peaks for a Transcription factor\n";
	print STDERR "\t\tyou might want to use the -center option (below) to center peaks on the known motif\n";
	print STDERR "\t\t** If using \"-size given\", histogram will be scaled to each region (i.e. 0-100%), with\n";
	print STDERR "\t\tthe -hist parameter being the number of bins to divide each region into.\n";
	print STDERR "\t\t\tHistogram Mode specific Options:\n";
	print STDERR "\t\t\t-nuc (calculated mononucleotide frequencies at each position,\n";
	print STDERR "\t\t\t\tWill report by default if extracting sequence for other purposes like motifs)\n";
	print STDERR "\t\t\t-di (calculated dinucleotide frequencies at each position)\n";
	print STDERR "\t\t\t-histNorm <#> (normalize the total tag count for each region to 1, where <#> is the\n";
	print STDERR "\t\t\t\tminimum tag total per region - use to avoid tag spikes from low coverage\n";
	print STDERR "\t\t\t-ghist (outputs profiles for each gene, for peak shape clustering)\n";
	print STDERR "\t\t\t-rm <#> (remove occurrences of same motif that occur within # bp)\n";

	print STDERR "\n\tPeak Centering: (other options are ignored)\n";
	print STDERR "\t\t-center <motif file> (This will re-center peaks on the specified motif, or remove peak\n";
	print STDERR "\t\t\tif there is no motif in the peak.  ONLY recentering will be performed, and all other\n";
	print STDERR "\t\t\toptions will be ignored.  This will output a new peak file that can then be reanalyzed\n";
	print STDERR "\t\t\tto reveal fine-grain structure in peaks (It is advised to use -size < 200) with this\n";
	print STDERR "\t\t\tto keep peaks from moving too far (-mirror flips the position)\n";
	print STDERR "\t\t-multi (returns genomic positions of all sites instead of just the closest to center)\n";

	print STDERR "\n\tGenome comparisons (need genome & liftOver)\n";
	print STDERR "\t\t-cmpGenome <genome1> [genome2] (Genomes to compare for sequence/motifs)\n";
	print STDERR "\t\t-cmpLiftover <liftover1> [genome2] (Genomes to compare for sequence/motifs)\n";
	#print STDERR "\t\t-revLiftover <liftover1> [genome2] (Genomes to compare for sequence/motifs)\n";
	#
	print STDERR "\n\tNormalization options:\n";
	print STDERR "\t\t-fpkm (normalize read counts to million reads or fragments per kilobase mapped)\n";
	print STDERR "\t\t-raw (do not adjust the tag counts based on total tags sequenced, -noadj works too)\n";
	print STDERR "\t\t-norm <#> (normalize tags to this tag count, default=1e7, 0=average tag count in all directories)\n";
	print STDERR "\t\t-normLength <#> (Fragment length to normlize to for experiments with different lens, def: 100)\n";
	print STDERR "\t\t-log (output tag counts as log2(x+1+rand) values - for scatter plots)\n";
	print STDERR "\t\t-sqrt (output tag counts as sqrt(x+rand) values - for scatter plots)\n";
	print STDERR "\t\t-ratio (process tag values as ratios - i.e. chip-seq, or mCpG/CpG)\n";

	print STDERR "\n\tAdvanced Options:\n";
	print STDERR "\t\t-len <#> / -fragLength <#> (Fragment length, default=auto, might want to set to 1 for 5'RNA)\n";
	print STDERR "\t\t-size <#> (Peak size[from center of peak], default=inferred from peak file)\n";
	print STDERR "\t\t\t-size #,# (i.e. -size -10,50 count tags from -10 bp to +50 bp from center)\n";
	print STDERR "\t\t\t-size \"given\" (count tags etc. using the actual regions - for variable length regions)\n";
	print STDERR "\t\t-strand <+|-|both> (Count tags on specific strands relative to peak, default: both)\n";
	#print STDERR "\t\t-local # (size in bp to count tags as a local background\n";
	print STDERR "\t\t-pc <#> (maximum number of tags to count per bp, default=0 [no maximum], -tbp <#> works too)\n";
	#print STDERR "\t\t-cons (Retrieve conservation information for peaks/sites)\n";
	print STDERR "\t\t-CpG (Calculate CpG/GC content)\n";
	print STDERR "\t\t-nfr (report nuclesome free region scores instead of tag counts, also -nfrSize <#>)\n";
	print STDERR "\t\t-norevopp (do not search for motifs on the opposite strand [works with -center too])\n";
	print STDERR "\t\t-gwasCatalog <gwasCatalog file from UCSC> (list overlapping GWAS risk SNPs)\n";
	#print STDERR "\t\t-snp <file> <id> (genotype file)\n";
	print STDERR "\t\t-pdist (only report distance to nearest peak using -p, not peak name)\n";
	print STDERR "\t\t-map <mapping file> (mapping between peak IDs and promoter IDs, overrides closest assignment)\n";
	print STDERR "\t\t-noann, -nogene (skip genome annotation step, skip TSS annotation)\n";
	print STDERR "\t\t-homer1/-homer2 (by default, the new version of homer [-homer2] is used for finding motifs)\n";
	print STDERR "\t\t-cpu <#> (Number of processors to use when possible - only some parts utilize multiple cores)\n";
	print STDERR "\t\t-noblanks (remove peaks/rows with missing data)\n";
	print STDERR "\n";
	exit;
	

}

if (@ARGV < 2) { 
	printCMD();
}

my $cmd = "annotatePeaks.pl";
for (my $i=0;$i<@ARGV;$i++) {
	$cmd .= " $ARGV[$i]";
}

print STDERR "\n";
my %toDelete = ();

$maxHomer2SeqLength = 1e7;

my $maxCPUs = 1;
my $seqFlag = 0;
my $consFlag = 0;
my $cpgFlag = 0;
my $skipBlastn = 1;
my $maskFlag = 0;
my $noblanksFlag = 0;

my $annStatFile = '';
my $normLength = 100;
my $gwasCatalog = "";
my $customAnnotationFile = "";
my $newPeakHistogramFlag=1;
my $homer2Flag = 1;
my $nfrFlag = 0;
my $nfrSize = 100;
my $mfastaFile = "";
my $fpkmFlag = 0;
my $gsize=2e9;
my $removeCloseMotifs=0;
my $strandFlag = "both";
my $logFlag = 0;
my $mdistFlag= 0;
my $mlogicFile = '';
my $sqrtFlag = 0;
my $pCountFlag = 0;
my $ratioFlag = 0;
my $snpFile = '';
my $snpID = '';
my $mscoreFlag = 0;
my $nscoreFlag = 0;
my $normCpGFlag = 0;
my $mirrorFlag = 0;
my $fragLength = 'auto';
my $revoppFlag = 1;
my $init2One = 0;
my $size = 300;
my $sizeMove = 0;
my $updateSize = 0;
my $local = 0;
my $adjustFlag = 1;
my $normValue = 1e7;
my $tssFlag = 0;
my $geneListFile = '';
my $ugFile = '';
my @geneDataFiles = ();
my $histBinSize = 0;
my $diFlag = 0;
my $nucFlag = 0;
my $mapFile = "";
my $centerMotif = '';
my $cpromoter = '';
my $histNorm = 0;
my $ghistFlag = 0;
my $goDir = '';
my $genomeOntologyDir = '';
my $pDistFlag = 0;
my $multiFlag = 0;
my $matrixPrefix = '';
my $matrixMinDist = 4;
my $matrixMaxDist = 1e10;
my $mbedFile = '';
my $noAnnFlag = 0;
my $noGeneFlag = 0;
my $mpeak = 0;
my $gtfFile = "";
my $rmMotifThresh = 10;
my $removeRevoppMotifs = 0;
my $vcfFile = "";
my @individuals = ();
my $editDistanceFlag = 0;
my $gffFlag = '';
my $bestScoreFlag = 0;
my $gidFlag = '';
my @cmpGenomes = ();
my @cmpLiftover = ();
my @revLiftover = ();

my $peakFile = $ARGV[0];
my $genome = $ARGV[1];
if ($genome =~ s/r$//) {
	$maskFlag = 1;
}

my $organism = "unknown";
my $promoter = "default";
my $consDir = "";
my $genomeDir = "";
my $genomeParseDir = "";
my $customGenome = 0;
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
if ($ARGV[0] eq 'rna') {
	$size = 'given';
	$updateSize = 1;
}


print STDERR "\tPeak file = $peakFile\n";
print STDERR "\tGenome = $genome\n";
print STDERR "\tOrganism = $organism\n";


my @motifFiles = ();
my @tagDirs = ();
my @wigFiles = ();
my @bedGraphFiles = ();
my @peakFiles = ();
my %filterMotifFiles = ();

for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-m' || $ARGV[$i] eq '-fm') {
		my $code = $ARGV[$i];
		$seqFlag =1;
		my $bail = 0;
		print STDERR "\tMotif files:\n";
		while ($ARGV[++$i] !~ /^\-/) {
			push(@motifFiles, $ARGV[$i]);
			$filterMotifFiles{$ARGV[$i]}=$code;
			print STDERR "\t\t$ARGV[$i]\t$code\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-p') {
		print STDERR "\tPeak Files:\n";
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@peakFiles, $ARGV[$i]);
			print STDERR "\t\t$ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-bedGraph') {
		print STDERR "\tbedGraph Files:\n";
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@bedGraphFiles, $ARGV[$i]);
			print STDERR "\t\t$ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-wig') {
		print STDERR "\tWiggle Files:\n";
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@wigFiles, $ARGV[$i]);
			print STDERR "\t\t$ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-cmpGenome') {
		print STDERR "\tGenomes for comparison:\n";
		$seqFlag =1;
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@cmpGenomes, $ARGV[$i]);
			print STDERR "\t\t$ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-cmpLiftover') {
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@cmpLiftover, $ARGV[$i]);
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-revLiftover') {
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@revLiftover, $ARGV[$i]);
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
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
	} elsif ($ARGV[$i] eq '-individuals') {
		print STDERR "\tVCF file individuals to analyze:\n";
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@individuals, $ARGV[$i]);
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
		if ($size eq 'given') {
			print STDERR "\tHistogram mode activated (bin size = 1/$histBinSize)\n";
		} else {
			print STDERR "\tHistogram mode activated (bin size = $histBinSize bp)\n";
		}
		#print STDERR "\tHistogram mode activated (bin size = $histBinSize bp)\n";
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
	} elsif ($ARGV[$i] eq '-size') {
		$size = $ARGV[++$i];
		if ($size eq 'given') {
			print STDERR "\tUsing actual sizes of regions\n";
		} elsif ($size =~ /\,/) {
			my @a = split /\,/, $size;
			my $sizeStart= $a[0];
			my $sizeEnd = $a[1];
			if ($sizeEnd < $sizeStart) {
				print STDERR "!!! Size end must be less than the size start range in -size $sizeStart,$sizeEnd\n";
				exit;
			}
			$sizeMove = floor(($sizeStart+$sizeEnd)/2);
			$size = $sizeEnd - $sizeStart;
		}
		$updateSize = 1;
		print STDERR "\tPeak Region set to $size\n";
	} elsif ($ARGV[$i] eq '-matrixMaxDist') {
		$matrixMaxDist = $ARGV[++$i];
		print STDERR "\tWhen producing a motif co-occurence matrix, will only consider co-bound if < $matrixMaxDist bp away\n";
	} elsif ($ARGV[$i] eq '-map') {
		print STDERR "\tWill map peaks to promoters using map file: $mapFile\n";
		$mapFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-matrixMinDist') {
		$matrixMinDist = $ARGV[++$i];
		print STDERR "\tWhen producing a motif co-occurence matrix, will only consider co-bound if > $matrixMinDist bp away\n";
	} elsif ($ARGV[$i] eq '-noblanks' || $ARGV[$i] eq '-noBlanks') {
		print STDERR "\tWill remove rows with data values of '' or 'NA'\n";
		$noblanksFlag = 1;
	} elsif ($ARGV[$i] eq '-matrix') {
		$matrixPrefix = $ARGV[++$i];
		print STDERR "\tWill produce a motif co-occurence of motifs analysis files, prefix: $matrixPrefix\n";
	} elsif ($ARGV[$i] eq '-rmrevopp') {
		$removeCloseMotifs = -1;
		$removeRevoppMotifs = 1;
		$rmMotifThresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rm') {
		$removeCloseMotifs = 1;
		$rmMotifThresh = $ARGV[++$i];
		if ($rmMotifThresh < 0) {
			$removeCloseMotifs = -1;
			$rmMotifThresh = abs($rmMotifThresh);
		}
	} elsif ($ARGV[$i] eq '-organism') {
		$organism = $ARGV[++$i];
		print STDERR "\tSetting organism to: $organism\n";
	} elsif ($ARGV[$i] eq '-mlogic') {
		$mlogicFile = $ARGV[++$i];
		print STDERR "\tWill outptu motif orientation stats to file: $mlogicFile\n";
	} elsif ($ARGV[$i] eq '-normLength') {
		$normLength = $ARGV[++$i];
		print STDERR "\tSetting normalization length to $normLength (set to 0 to disable)\n";
	} elsif ($ARGV[$i] eq '-vcf') {
		$vcfFile = $ARGV[++$i];
		print STDERR "\tWill get SNP info from VCF file: $vcfFile\n";
	} elsif ($ARGV[$i] eq '-editDistance') {
		$editDistanceFlag = 1;
		print STDERR "\tWill calculate total variation (edit disance) from reference sequence\n";
	} elsif ($ARGV[$i] eq '-gwasCatalog') {
		$gwasCatalog = $ARGV[++$i];
		print STDERR "\tWill annotate transcripts using GWAS catalog (from file: $gwasCatalog)\n";
	} elsif ($ARGV[$i] eq '-mbed') {
		$mbedFile = $ARGV[++$i];
		print STDERR "\tWill produce a motif bed file: $mbedFile\n";
	} elsif ($ARGV[$i] eq '-mpeak') {
		$mbedFile = $ARGV[++$i];
		$mpeak = 1;
		print STDERR "\tWill produce a motif peak file: $mbedFile\n";
	} elsif ($ARGV[$i] eq '-mask') {
		$maskFlag = 1;
		print STDERR "\tWill use repeat-masked sequences\n";
	} elsif ($ARGV[$i] eq '-TSS') {
		$promoter = $ARGV[++$i];
		print STDERR "\tPromoter Set will be $promoter\n";
	} elsif ($ARGV[$i] eq '-annStats') {
		$annStatFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-snp') {
		$snpFile = $ARGV[++$i];
		$snpID = $ARGV[++$i];
		print STDERR "\tWill use $snpFile with ID: $snpID for genotype\n";
	} elsif ($ARGV[$i] eq '-strand') {
		$strandFlag = $ARGV[++$i];
		print STDERR "\tWill count tags on strand: $strandFlag\n";
	} elsif ($ARGV[$i] eq '-mfasta') {
		$mfastaFile = $ARGV[++$i];
		print STDERR "\tWill output motif site sequences to FASTA file: $mfastaFile\n";
	} elsif ($ARGV[$i] eq '-cTSS') {
		$cpromoter = $ARGV[++$i];
		print STDERR "\tCustom promoter set will be $cpromoter\n";
	} elsif ($ARGV[$i] eq '-ann') {
		$customAnnotationFile = $ARGV[++$i];
		print STDERR "\tCustom annotation file: $customAnnotationFile\n";
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
		print STDERR "\tWill use up to $maxCPUs where possible\n";
	} elsif ($ARGV[$i] eq '-pc' || $ARGV[$i] eq '-tbp') {
		$init2One = $ARGV[++$i];
		print STDERR "\tMaximum count per bp will be set to $init2One\n";
	} elsif ($ARGV[$i] eq '-nmotifs') {
		$nscoreFlag = 1;
		print STDERR "\tWill report the number of motifs in each peak\n";
	} elsif ($ARGV[$i] eq '-mscore') {
		$mscoreFlag = 1;
		print STDERR "\tWill report max log-odds score for motif in each peak\n";
	} elsif ($ARGV[$i] eq '-log') {
		$logFlag = 1;
		print STDERR "\tWill output log(1+rand()+x) for tag counts\n";
	} elsif ($ARGV[$i] eq '-sqrt') {
		$sqrtFlag = 1;
		print STDERR "\tWill output sqrt(rand()+x) for tag counts\n";
	} elsif ($ARGV[$i] eq '-normCpG') {
		$normCpGFlag = 1;
		print STDERR "\tTreating tag files like CpG methylation ratios\n";
	} elsif ($ARGV[$i] eq '-gsize') {
		$gsize = $ARGV[++$i];
		print STDERR "\tEffective Genome size set to $gsize\n";
	} elsif ($ARGV[$i] eq '-nfr') {
		$nfrFlag = 1;
		print STDERR "\tReporting tag directory counts as Nucleosome Free Region scores ($nfrSize bp size)\n";
	} elsif ($ARGV[$i] eq '-nfrSize') {
		$nfrSize = $ARGV[++$i];
		print STDERR "\tUsing Nucleosome Free Region size of $nfrSize bp\n";
	} elsif ($ARGV[$i] eq '-ratio') {
		$ratioFlag = 1;
		print STDERR "\tTreating tag values like ratios\n";
	} elsif ($ARGV[$i] eq '-pdist2') {
		$pDistFlag = 2;
		print STDERR "\tWill only report distance for nearest peaks\n";
	} elsif ($ARGV[$i] eq '-mdist') {
		$mdistFlag = 1;
		print STDERR "\tReports distance to nearest motif\n";
	} elsif ($ARGV[$i] eq '-pdist') {
		$pDistFlag = 1;
		print STDERR "\tWill only report absolute distance for nearest peaks\n";
	} elsif ($ARGV[$i] eq '-pcount') {
		$pCountFlag = 1;
		print STDERR "\tWill report number of peaks\n";
	} elsif ($ARGV[$i] eq '-ghist') {
		$ghistFlag = 1;
		print STDERR "\tWill create histogram for each gene\n";
	} elsif ($ARGV[$i] eq '-noadj' || $ARGV[$i] eq '-raw') {
		$adjustFlag = 0;
		print STDERR "\tWill NOT normalize tag counts\n";
	} elsif ($ARGV[$i] eq '-mirror') {
		$mirrorFlag =1;
		print STDERR "\tReturning mirrored positions\n";
	} elsif ($ARGV[$i] eq '-list') {
		$geneListFile = $ARGV[++$i];
		print STDERR "\tPromoters for genes found in $geneListFile will be analyzed\n";
	} elsif ($ARGV[$i] eq '-nogene') {
		$noGeneFlag = 1;
		print STDERR "\tWill Skip closest gene annotation\n";
	} elsif ($ARGV[$i] eq '-noann') {
		$noAnnFlag = 1;
		print STDERR "\tWill Skip peak annotation\n";
	} elsif ($ARGV[$i] eq '-bestScore') {
		$bestScoreFlag = 1;
	} elsif ($ARGV[$i] eq '-center') {
		$centerMotif = $ARGV[++$i];
		$seqFlag = 1;
		print STDERR "\tPeaks/Regions will be centered on motif in file $centerMotif\n";
	} elsif ($ARGV[$i] eq '-norm') {
		$normValue = $ARGV[++$i];
		if ($normValue == 0) {
			print STDERR "\tNormalzing tags to the Average Tag totals\n";
		} else {
			print STDERR "\tWill normalize tag counts to $normValue per experiment\n";
		}
	} elsif ($ARGV[$i] eq '-norevopp') {
		$revoppFlag = 0;
		print STDERR "\tWill not search for motifs on the opposite strand\n";
	} elsif ($ARGV[$i] eq '-go') {
		$goDir = $ARGV[++$i];
		print STDERR "\tWill perform Gene Ontology analysis - output to directory = $goDir\n";
	} elsif ($ARGV[$i] eq '-genomeOntology') {
		$genomeOntologyDir = $ARGV[++$i];
		print STDERR "\tWill perform Genome Ontology analysis - output to directory = $genomeOntologyDir\n";
		print STDERR "\t\tWarning - might want to set the genome size with -gsize (currently $gsize)\n";
	} elsif ($ARGV[$i] eq '-multi') {
		$multiFlag = 1;
		print STDERR "\tWill return all motif positions when centering...\n";
	} elsif ($ARGV[$i] eq '-local') {
		$local = $ARGV[++$i];
		print STDERR "\tWill count tags in local backgound in $local bp around peak\n";
	} elsif ($ARGV[$i] eq '-blastn') {
		$skipBlastn = 0;
	} elsif ($ARGV[$i] eq '-cons') {
		$consFlag = 1;
		print STDERR "\tWill extract conservation information\n";
	} elsif ($ARGV[$i] eq '-len' || $ARGV[$i] eq '-fragLength') {
		$fragLength = $ARGV[++$i];
		print STDERR "\tFragment Length set to $fragLength\n";
	} elsif ($ARGV[$i] eq '-di') {
		$seqFlag =1;
		$diFlag = 1;
		print STDERR "\tWill report dinucleotide frequencies\n";
	} elsif ($ARGV[$i] eq '-gtf') {
		$gtfFile = $ARGV[++$i];
		print STDERR "\tCustom annotation GTF file: $gtfFile (using transcript_id)\n";
	} elsif ($ARGV[$i] eq '-gid') {
		$gidFlag = " -gid";
		print STDERR "\tUsing gene_ids for GTF file\n";
	} elsif ($ARGV[$i] eq '-gff') {
		$gtfFile = $ARGV[++$i];
		$gffFlag = ' -gff';
		print STDERR "\tCustom annotation GFF file: $gtfFile (better to get GTF file)\n";
	} elsif ($ARGV[$i] eq '-gff3') {
		$gtfFile = $ARGV[++$i];
		$gffFlag = ' -gff3';
		print STDERR "\tCustom annotation GFF3 file: $gtfFile (better to get GTF file)\n";
	} elsif ($ARGV[$i] eq '-nuc') {
		$seqFlag =1;
		print STDERR "\tWill report nucleotide frequencies\n";
	} elsif ($ARGV[$i] eq '-histNorm') {
		$histNorm = $ARGV[++$i];
		print STDERR "\tWill normalize Tag histograms with minimum total of $histNorm\n";
	} elsif ($ARGV[$i] eq '-homer2') {
		$homer2Flag = 1;
		print STDERR "\tUsing homer2...\n";
	} elsif ($ARGV[$i] eq '-fpkm' || $ARGV[$i] eq '-rpkm') {
		$fpkmFlag = 1;
		$normValue = 1e6;
		print STDERR "\tWill normalized reads to FPKM\n";
	} elsif ($ARGV[$i] eq '-homer1') {
		$homer2Flag = 0;
		print STDERR "\tUsing original homer...\n";
	} elsif ($ARGV[$i] eq '-CpG') {
		$cpgFlag = 1;
		$seqFlag = 1;
		print STDERR "\tWill calculate CpG/GC content\n";
	} else {
		print STDERR "$ARGV[$i] not recognized\n\n";
		printCMD();
	
	}
}

my %alpha2index = ();
$alpha2index{'A'} = 0;
$alpha2index{'a'} = 0;
$alpha2index{'C'} = 1;
$alpha2index{'c'} = 1;
$alpha2index{'G'} = 2;
$alpha2index{'g'} = 2;
$alpha2index{'T'} = 3;
$alpha2index{'t'} = 3;


if ($updateSize == 0 && $histBinSize == 0) {
	#print STDERR "size set to given\n";
	$size = "given";
	$updateSize = 1;
}

my $cmpGenomeFlag = 0;
if (scalar(@cmpGenomes) > 0 || scalar(@cmpLiftover) > 0) {
	$cmpGenomeFlag = 1;
	if (scalar(@cmpGenomes) != scalar(@cmpLiftover)) {
		print STDERR "!!! Error - Each -cmpGenome genome entry needs a matching -cmpLiftover\n";
		exit;
	}
	if (scalar(@revLiftover) > 0) {
		if (scalar(@revLiftover) != scalar(@cmpLiftover)) {
			print STDERR "!!! Error - Each -cmpLiftover should match a -revLiftover\n";
			exit;
		}
	}
	if ($skipBlastn) {
		print STDERR "\tAdd '-blastn' to command to calculate % identities (can take a while...)\n";
	} else {
		print STDERR "\tWill use blastn to calculate % identities\n";
	}
}

my %cmpGenomeInfo = ();
for (my $i=0;$i<@cmpGenomes;$i++) {
	my $cgenome = $cmpGenomes[$i];
	my $cmaskFlag = 0;
	if ($cgenome =~ s/r$//) {
		$cmaskFlag = 1;
	}

	my $cgenomeDir = "";
	my $ccustomGenome = 0;

	if ($cgenome eq 'none') {
		print STDERR "\tA -cmpGenome cannot be 'none'\n";
		exit;
		$ccustomGenome = -1;
		$cgenomeDir = "none";
	} elsif (!exists($config->{'GENOMES'}->{$cgenome})) {
		$ccustomGenome = 1;
		($cgenome,$cgenomeDir,$cgenomeParseDir) = HomerConfig::parseCustomGenome($cgenome);
		$cgenomeParseDir = '';
	} else {
		$cgenomeDir = $config->{'GENOMES'}->{$cgenome}->{'directory'};	
	}
	$cmpGenomeInfo{$cmpGenomes[$i]} = {genome=>$cgenome,dir=>$cgenomeDir,custom=>$ccustomGenome,mask=>$cmaskFlag};
}

my $mflag = '';
if ($maskFlag) {
	$mflag = " -mask ";
}

$promoterIDtype = 'refseq';
if ($promoter ne 'default') {
	if (!exists($config->{'PROMOTERS'}->{$promoter})) {
		print STDERR "!! Promoter Set $promoter not recognized\n";
		exit;
	}
	$promoterIDtype = $config->{'PROMOTERS'}->{$promoter}->{'idtype'};
}
if ($gtfFile) {
	$promoterIDtype = 'custom';
}


my $halfLocal = floor($local/2);

my %peaks = ();
my %gene = ();

#tmp files
my $rand = rand();
my $tmpfile = $rand . ".tmp";
my $tmpfile2 = $rand . ".2.tmp";
my $tmpfile3 = $rand . ".3.tmp";
my $tmpfile4 = $rand . ".4.tmp";
my $tmpfile5 = $rand . ".5.tmp";
my $tmpfile6 = $rand . ".6.tmp";
my $tmpfile7 = $rand . ".7.tmp";
my $seqFile = $rand . ".seq";
my $seqFaFile = $rand . ".seq.fa";
my $consFile = $rand . ".cons";
my $posFile = $rand . ".pos";
my $gtfTSSFile = $rand . ".gtf.tss";
my $tmpPeakFile = $rand . ".peak";
my $cleanPosFile = $rand . ".clean.pos";

my $ogPeakFile = $peakFile;

my %geneList = ();
my $specialFile = $genomeDir . "/" . $genome . "." . $peakFile;
if ($peakFile eq 'tss' || $peakFile eq 'tts' || $peakFile eq 'rna' || -e $specialFile) {

	print STDERR "\tFound special file: $specialFile\n";
	my $postfix = $peakFile;

	$tssFlag = 1;
	print STDERR "\n\t*****************\n\t$postfix Mode enabled\n";
	print STDERR "\t*****************\n\n";

	if ($gtfFile ne '') {
		`parseGTF.pl "$gtfFile" $peakFile $gffFlag $gidFlag > "$gtfTSSFile"`;
		$size = 4000 if ($updateSize == 0 && $peakFile ne 'rna');
		$peakFile = $gtfTSSFile;
		$toDelete{$gtfTSSFile}=1;
	} elsif ($cpromoter ne '') {
		$size = 4000 if ($updateSize == 0 && $peakFile ne 'rna');
		$peakFile = $cpromoter;
	} elsif ($promoter eq 'default') {
		$size = 4000 if ($updateSize == 0 && $peakFile ne 'rna');
		$peakFile = $genomeDir . "/" . $genome . "." . $postfix;	
		if (-f $peakFile) {
		} else {
			print STDERR "!!! This isn't going to work - can't find $peakFile !!!\n";
			print STDERR "!!! Have you upgraded lately? !!!\n";
			print STDERR "!!! Can't use the \"tss\" option with a custom genome - don't know where the TSS are... !!!\n";
			exit;
		}
	} else {
		$size = 2*$config->{'PROMOTERS'}->{$promoter}->{'end'} if ($updateSize == 0 && $peakFile ne 'rna');
		$peakFile = $config->{'PROMOTERS'}->{$promoter}->{'directory'} . "/" . $promoter . ".pos";
	}


	if ($geneListFile ne '') {
		my $cflag = 0;
		if ($promoterIDtype eq 'null' || $promoterIDtype eq 'custom' || $organism eq 'unknown') {
			`cut -f1 "$geneListFile" > $tmpfile`;
		} else {
			$cflag = 1;
			`convertIDs.pl "$geneListFile" $organism $promoterIDtype no yes yes > "$tmpfile"`;
		}
		`cat "$geneListFile" "$tmpfile" | cut -f1 | sort | uniq > "$tmpfile2"`;
		`cp "$tmpfile2" "$tmpfile3"`;
		`cat "$tmpfile2" "$tmpfile3" | sort | uniq > "$tmpfile"`;
		#print STDERR "`mergeData.pl $tmpfile $peakFile 0 -accVer  | sort | uniq >  $tmpPeakFile`\n";
		`mergeData.pl "$tmpfile" "$peakFile" 0 -accVer  | sort | uniq >  "$tmpPeakFile"`;
		$peakFile = $tmpPeakFile;
		#$toDelete{$tmpPeakFile} = 1;
		`rm "$tmpfile" "$tmpfile2" "$tmpfile3"`;
	}
} else {
	`bed2pos.pl "$peakFile" -check -unique > "$cleanPosFile"`;
	`checkPeakFile.pl "$cleanPosFile"`;
	`cleanUpPeakFile.pl "$cleanPosFile" > "$tmpfile"`;
	`mv "$tmpfile" "$cleanPosFile"`;
	$toDelete{$cleanPosFile}=1;
	$peakFile = $cleanPosFile;
}

if ($updateSize && $size ne 'given') {
	print STDERR "\tResizing peaks...\n";
	`resizePosFile.pl "$peakFile" $size $sizeMove > "$posFile"`;
} else {
	#print STDERR "`cp $peakFile $posFile`;\n";
	open OUT, ">$posFile";
	open IN, $peakFile;
	while (<IN>) {
		chomp;	
		s/\r//g;
		my @line = split /\t/;
		foreach(@line) {
			s/^\s*//g;
			s/\s*$//g;
		}
		print OUT "$line[0]";
		for (my $i=1;$i<@line;$i++) {
			print OUT "\t$line[$i]";
		}
		print OUT "\n";
	}
	close IN;
	close OUT;
}

`rm "$tmpPeakFile"` if ($tmpPeakFile eq $peakFile);


#first extract sequence for positions
if ($seqFlag) {
	$cpgFlag = 1;
	print STDERR "\tExtracting Sequence...\n";
	if ($genome eq 'none') {
		print STDERR "!!! Cannot do anything with genomic sequences (i.e. motifs) if genome is set to: $genome !!!\n";
		print STDERR "!!! Try installing a genome through HOMER or provide a FASTA directory/file\n";
		exit;
	}
	`homerTools extract "$posFile" "$genomeDir" $mflag > "$seqFile"`;
	if ($snpFile ne '') {

	}
}


print STDERR "\tReading Positions...\n";

my $avgSize = 0;
my $avgSizeN = 0;
my $maxSize = 0;
open IN, $posFile;
my @peakOrder = ();
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line = split /\t/;
	my $hit = $line[0];
	my $chr = $line[1];
	my $start = $line[2];
	next if ($start =~ /^[^\d\-\.]/);
	my $end = $line[3];
	my $peakSize = $end-$start+1;
	if ($updateSize == 0) {
		$size = $peakSize;
	}
	my $direction = '+';
	if ($line[4] eq '-' || $line[4] eq '1') {
		$direction = '-';
	}
	my $value = 0;
	my $fdr = 'NA';
	if (@line > 5) {
		$value = $line[5];
	}
	if (@line > 6) {
		$fdr = $line[6];
	}
	my %a = ();
	my %b = ();
	my %c = ();
	$avgSize += $peakSize;
	$avgSizeN++;
	$maxSize = $peakSize if ($peakSize > $maxSize);

	my @compData = ();
	if ($cmpGenomeFlag) {
		foreach(@cmpGenomes) {
			my $compData = {map=>0,score=>0,indel=>0,var=>0,pid=>'NA',paln=>0};
			push(@compData, $compData);
		}
	}


	$peaks{$hit} = {tss=>'NA', tssDist=>'NA',m=>\%a, p=>\%c, cons=>'NA', size=>$peakSize, tssUG=>'NA',
				s=>$start, e=>$end, d=>$direction, v=>$value, c=>$chr, t=>\%b, gc=>'NA',cpg=>'NA',
				centerDist=>1e10,centerDir=>1e10,centerScore=>-10,fdr=>$fdr,ann=>'NA',
				fullann=>'NA',gComp=>\@compData};
	push(@peakOrder, $hit);
}
close IN;
my $halfSize = 0;
if ($size ne 'given') {
	$halfSize = floor($size/2);
	$avgSize = $size;
} else {
	$avgSize /= $avgSizeN if ($avgSizeN > 0);
}

if (-e "$seqFile" && $maxSize > $maxHomer2SeqLength) {
	`tab2fasta.pl "$seqFile" > "$seqFaFile"`;
	$toDelete{$seqFaFile}=1;
}

if ($seqFlag && $snpFile ne '') {
	#open IN, $ARGV[0];	
}


if ($centerMotif ne '') {
	print STDERR "\tLooking for motifs to center regions...\n";

	my $offset = -1*$halfSize;
	#my $offset = -1*$halfSize+$sizeMove;
	#print STDERR "halfSize=$halfSize\nsizeMove=$sizeMove\noffset=$offset\n";
	if ($size eq 'given') {
		$offset = 0;
	}

	my $mfile = $centerMotif;
	my $options = '';

	my $seqFile2Use = " -s \"$seqFile\"";
	if ($maxSize > $maxHomer2SeqLength) {
		$seqFile2Use = " -i \"$seqFaFile\"";
	}

	if ($homer2Flag) {
		$options .= " -strand +" if ($revoppFlag == 0);
		$options .= " -p $maxCPUs";
		if ($mscoreFlag) {
			`homer2 find $seqFile2Use -offset $offset -m "$mfile" -mscore $options > "$tmpfile"`;
		} else {
			`homer2 find $seqFile2Use -offset $offset -m "$mfile" $options > "$tmpfile"`;
			#print STDERR "`homer2 find $seqFile2Use -offset $offset -m $mfile $options > $tmpfile`;\n";
		}
	} else {
		$options .= " -norevopp" if ($revoppFlag == 0);
		if ($mscoreFlag) {
			`homer -s "$seqFile" -a GENESCORE -offset $offset -m "$mfile" $options > "$tmpfile"`;
		} else {
			`homer -s "$seqFile" -a FIND -offset $offset -m "$mfile" $options > "$tmpfile"`;
		}
	}
	my $numMultiMotifs = 0;

	open IN, $tmpfile;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $hit = "";
		my $pos = "";
		my $seq = "";
		my $con = 0;
		my $dir = 0;
		my $score = 0;
		if ($homer2Flag) {
			$hit = $line[0];
			$pos = $line[1];
			$seq = $line[2];
			$con = 0;
			$dir = $line[4];
			$score = $line[5];
		} else {
			$hit = $line[0];
			$pos = $line[1];
			$seq = $line[2];
			$con = $line[3];
			$dir = $line[4];
			$score = $line[6];
		}
		if ($dir eq '+' || $dir eq '0') {
			$dir = 0;
		} else {
			$dir = 1;
			$pos += length($seq)-1 unless ($homer2Flag);
		}
		if ($mirrorFlag==1) {
			$pos *= -1;
		}
		if ($size eq 'given') {
			my $L = $peaks{$hit}->{'e'}-$peaks{$hit}->{'s'};
			$pos = $pos-floor($L/2);
		}
				
		#provide information about centering
		if (abs($pos) < abs($peaks{$hit}->{'centerDist'})) {
			$peaks{$hit}->{'centerDist'} = $pos;
			$peaks{$hit}->{'centerDir'} = $dir;
			$peaks{$hit}->{'centerScore'} = $score;
		}
		if ($multiFlag) {
			my $chr = $peaks{$hit}->{'c'};
			my $start = $peaks{$hit}->{'s'};
			my $end = $peaks{$hit}->{'e'};
			my $v = $peaks{$hit}->{'v'};
			my $dd = $peaks{$hit}->{'d'};
			$dd = 1 if ($dd eq '-');
			$dd = 0 if ($dd eq '+');
			if ($dd == 0) {
				$start += $pos;
				$end += $pos;
				$dd = 1 if ($dir == 1);
			} else {
				$start -= $pos;
				$end -= $pos;
				$dd = 0 if ($dir == 1);
			}
			print "$hit\t$chr\t$start\t$end\t$dd\t$v\t$score\n";
			$numMultiMotifs++;
		}
	}
	close IN;
	`rm "$tmpfile" "$posFile" "$seqFile"`;
	if ($multiFlag) {
		deleteFiles();
		print STDERR "\tCentered on $numMultiMotifs motifs total\n\n";
		exit;
	}

	my $total = 0;
	my $goodPeaks = 0;
	my $totalChange = 0;
	foreach(keys %peaks) {
		$total++;
		my $hit = $_;
		next if ($peaks{$hit}->{'centerDist'} > 1e8);
		if ($size ne 'given') {
			next if (abs($peaks{$hit}->{'centerDist'}) > $size);
		}
		$goodPeaks++;
		my $start = $peaks{$hit}->{'s'};
		my $end = $peaks{$hit}->{'e'};
		my $dir = $peaks{$hit}->{'d'};
		if ($dir eq '+' || $dir eq '0') {
			$dir = 0;
		} else {
			$dir = 1;
		}
		my $chr = $peaks{$hit}->{'c'};
		my $v = $peaks{$hit}->{'v'};

		my $posOffset = $peaks{$hit}->{'centerDist'};
		my $score = $peaks{$hit}->{'centerScore'};
		$totalChange+=abs($posOffset);

		my $changeDir = $peaks{$hit}->{'centerDir'};
		if ($dir == 0) {
			$start += $posOffset;
			$end += $posOffset;
			$dir = 1 if ($changeDir == 1);
		} else {
			$start -= $posOffset;
			$end -= $posOffset;
			$dir = 0 if ($changeDir == 1);
		}
		print "$hit\t$chr\t$start\t$end\t$dir\t$v\t$score\n";
	}
	my $avgChange = 'N/A';
	my $status = 'failed';
	if ($goodPeaks > 0) {
		$avgChange = $totalChange/$goodPeaks;
		$status = "successful";
	}
	print STDERR "\nTotal Peaks/Regions:     $total\n";
	print STDERR "Total Peaks re-centered: $goodPeaks\n";
	print STDERR "Avg. Adjustement size:   $avgChange\n";
	print STDERR "\nPeak/Region centering was $status\n"; 
	deleteFiles();
	exit;
}


if ($consFlag == 1) {
	if (open IN, "$consDir/chr1.fa") {
		close IN;
		print STDERR "\tExtracting Conservation...\n";
		`homerTools extract "$posFile" "$consDir" > "$consFile"`;
	} else {
		print STDERR "Conservation information not present in $consDir (refer to documentation)\n";	
		print STDERR "Skipping Conservation\n";	
		$consFlag = 0;
	}

	print STDERR "\tCalculating average conservation...\n";
	`conservationPerLocus.pl "$consFile" > "$tmpfile"`;
	open IN, $tmpfile;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		next if (!exists($peaks{$line[0]}));
		$peaks{$line[0]}->{'cons'} = $line[1];
	}
	close IN;
	`rm "$tmpfile"`;
}

if ($seqFlag==1 && $cpgFlag == 1) {
	print STDERR "\tCalculating CpG/GC Content of Peaks\n";
	`homerTools freq "$seqFile" -gc "$tmpfile" > /dev/null`;
	open IN, $tmpfile;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		next if (!exists($peaks{$line[0]}));
		$peaks{$line[0]}->{'cpg'} = $line[1];
		$peaks{$line[0]}->{'gc'} = $line[2];
	}
	close IN;
	`rm "$tmpfile"`;
}

#find motifs in fragments
my %motifNames = ();
my @motifOrder = ();

print STDERR "\t-----------------------\n";

my $offset = -1*$halfSize+$sizeMove;
if ($size eq 'given') {
	$offset = 0;
}
my @motifNames = ();

if ($mbedFile ne '') {
	my $crapass = '';
	my $shitface = '';
	open MBED, ">$mbedFile" or print STDERR "Could not open file $mbedFile for writing!\n";
	if ($mpeak==0) {
		print MBED "track name=\"$ogPeakFile motifs\" description=\"$ogPeakFile motifs\""
					. " visibility=\"pack\" useScore=\"1\"\n";
	}
}

@motifFiles = sort {$filterMotifFiles{$a} cmp $filterMotifFiles{$b}} @motifFiles;
my @copyMotifFiles = @motifFiles;
my %motifMask = ();

if ($mfastaFile ne '') {
	open MFASTA, ">$mfastaFile";
}

for (my $i=0;$i<@motifFiles;$i++) {
	if ($i==0) {
		print STDERR "\tLooking for Motifs...\n";
	}
	my $mfile = $motifFiles[$i];


	open IN, $mfile;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if (@line < 2);
		if ($line[0] =~ /^>/) {
			my $mname = $line[1];

			if (!exists($motifNames{$mname})) {
				$motifNames{$mname} = 1;
				push(@motifNames, $mname);
			}
		}
	}
	close IN;



	my $options = '';
	my $code = $filterMotifFiles{$mfile};

	if ($homer2Flag) {
		$options .= " -strand +" if ($revoppFlag == 0);
		$options .= " -p $maxCPUs";
		my $seqFile2Use = " -s \"$seqFile\"";
		if ($maxSize > $maxHomer2SeqLength) {
			print STDERR "!!! Unfortunately, HOMER cannot currently scan sequences >1e7 in length)\n";
			print STDERR "!!! Try using scanMotifGenomeWide.pl to scan for motifs in large FASTA files\n";
			exit;
			$seqFile2Use = " -i \"$seqFaFile\"";
		}
		if ($mscoreFlag) {
			`homer2 find $seqFile2Use -offset $offset -m "$mfile" -mscore $options > "$tmpfile"`;
		} else {
			#print STDERR "Here(maxSize=$maxSize, $seqFile2Use): `homer2 find $seqFile2Use -offset $offset -m $mfile $options > $tmpfile`;\n";
			`homer2 find $seqFile2Use -offset $offset -m "$mfile" $options > "$tmpfile"`;
		}
	} else {
		$options .= " -norevopp" if ($revoppFlag == 0);
		if ($mscoreFlag) {
			`homer -s "$seqFile" -a GENESCORE -offset $offset -m "$mfile" $options > "$tmpfile"`;
		} else {
			`homer -s "$seqFile" -a FIND -offset $offset -m "$mfile" $options > "$tmpfile"`;
		}
	}

	if ($consFlag==1) {
		`mv "$tmpfile" "$tmpfile2"`;
		`getSiteConservation.pl "$consFile" "$tmpfile2" $offset > "$tmpfile"`;
		`rm "$tmpfile2"`;
	}

	open IN, $tmpfile;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $hit = $line[0];
		my $pos = $line[1];
		my $seq = $line[2];
		my $con = $line[3];
		my $dir = $line[4];
		my $mname = $line[5];
		my $score = 0;
		if ($homer2Flag) {
			$mname = $line[3];
			$score = $line[5];
			$con = 0 ;
			if ($dir eq '-' || $dir eq '1') {
				$pos -= length($seq)-1;
			}
		} else {
			$score = $line[6];
		}
		my $midpoint = $pos + floor(length($seq)/2);

		if ($mfastaFile ne '') {
			my $mfastaSeq = $seq;
			if ($dir eq '-' || $dir eq '1') {
				$mfastaSeq = HomerConfig::revopp($mfastaSeq);
			}
			print MFASTA ">$mname $score $hit $pos $dir\n$mfastaSeq\n";

		}


		if ($code eq '-fm') {
			if (!exists($mask{$hit})) {
				my %a = ();
				$mask{$hit} = \%a;
			}
			for (my $j=$pos;$j<=$pos+length($seq);$j++) {
				$mask{$hit}->{$j}=1;
			}
			next;
		} elsif ($code eq '-m') {
			my $bad = 0;
			if (exists($mask{$hit})) {
				for (my $j=$pos;$j<=$pos+length($seq);$j++) {
					if (exists($mask{$hit}->{$j})) {
						$bad = 1;
						last;
					}
				}
			}
			if ($bad) {
				next;
			}
		}


		my $pChr = $peaks{$hit}->{'c'};
		my $pStart = $peaks{$hit}->{'s'};
		my $pDir = $peaks{$hit}->{'d'};
		my $bedStart = $pStart + $pos - $offset-1;
		my $bedEnd = $bedStart + length($seq);
		if ($pDir eq '1' || $pDir eq '-') {
			$bedEnd = $peaks{$hit}->{'e'} - $pos + $offset;
			$bedStart = $bedEnd - length($seq);
		}
		my $bedDir = '+';
		if ($dir eq '-' || $dir eq '1') {
			$bedDir = '-';
			if ($pDir eq '-' || $pDir eq '1') {
				$bedDir = '+';
			}
		} else {
			if ($pDir eq '-' || $pDir eq '1') {
				$bedDir = '-';
			}
		}

		if ($mbedFile ne '') {
			if ($mpeak==0) {
				my $ss=$score;
				print MBED "$pChr\t$bedStart\t$bedEnd\t$mname\t$ss\t$bedDir\n";
			} else {
				$bedStart += 1;
				print MBED "$mname\t$pChr\t$bedStart\t$bedEnd\t$bedDir\t$score\t$seq\n";
			}
		}
		if (!exists($peaks{$hit}->{'m'}->{$mname})) {
			my @a = ();
			$peaks{$hit}->{'m'}->{$mname} = \@a;
		}

		my @compData = ();
		if ($cmpGenomeFlag) {
			foreach(@cmpGenomes) {
				my $compData = {map=>0,score=>0,indel=>0,var=>0,m=>0};
				push(@compData, $compData);
			}
		}


		my $m = {p=>$pos, s=>$seq, c=>$con, d=>$dir, valid=>1,score=>$score,mp=>$midpoint,
					gc=>$pChr,gs=>$bedStart,ge=>$bedEnd,gd=>$bedDir,gComp=>\@compData};
		push(@{$peaks{$hit}->{'m'}->{$mname}}, $m);

		if ($removeCloseMotifs != 0) {
			my $num = scalar(@{$peaks{$hit}->{'m'}->{$mname}});
			for (my $i=0;$i<$num;$i++) {
				next if ($peaks{$hit}->{'m'}->{$mname}->[$i]->{'valid'}==0);
				my $p1 = $peaks{$hit}->{'m'}->{$mname}->[$i]->{'p'};
				my $d1 = $peaks{$hit}->{'m'}->{$mname}->[$i]->{'d'};
				for (my $j=0;$j<$num;$j++) {
					next if ($j==$i);
					my $p2 = $peaks{$hit}->{'m'}->{$mname}->[$j]->{'p'};
					my $d2 = $peaks{$hit}->{'m'}->{$mname}->[$j]->{'d'};
					if (abs($p1-$p2) < $rmMotifThresh) {
						if (($removeCloseMotifs == 1 && $d2 eq $d1) 
								|| ($removeCloseMotifs == -1 && $d2 ne $d1)) {
							$peaks{$hit}->{'m'}->{$mname}->[$i]->{'valid'}=0 if ($removeRevoppMotifs == 0);
							$peaks{$hit}->{'m'}->{$mname}->[$j]->{'valid'}=0;
						}
					}
				}
			}
		}

		if (!exists($motifNames{$mname})) {
			$motifNames{$mname} = 1;
			push(@motifNames, $mname);
		}


		#provide information about centering
	}
	close IN;
	`rm "$tmpfile"`;
}
if ($mbedFile ne '') {
	close MBED;
}
if ($mfastaFile ne '') {
	close MFASTA;
}
if ($mlogicFile ne '') {
	calcMotifLogic(\%peaks,\@motifNames,$strandFlag,$mlogicFile);
}


## cmpGenomes section ....................................................................
my %motifSubMatrix = ();
if ($cmpGenomeFlag) {
	for (my $i=0;$i<@motifFiles;$i++) {
		my $mfile = $motifFiles[$i];
		open IN, $mfile;
		my $cur = '';
		my @matrix1 = ();
		my $matrix = \@matrix1;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			if ($line[0] =~ /^>/) {
				if ($cur ne '') {
					$motifSubMatrix{$cur} = $matrix;
				}
				$cur = $line[1];
				my @matrix = ();
				$matrix = \@matrix;
			} else {
				push(@$matrix, [0, 0, 0, 0]);
			}
		}
		close IN;
		if ($cur ne '') {
			$motifSubMatrix{$cur} = $matrix;
		}
	}
}
my @cpeaks = ();
for (my $z=0;$z<@cmpGenomes;$z++) {
	my $liftOver = $cmpLiftover[$z];
	my $revliftOver = "";
	$revliftOver = $revLiftover[$z] if (@revLiftover > $z);
	my $cgenome = $cmpGenomes[$z];
	my $cgenomeDir = $cmpGenomeInfo{$cgenome}->{'dir'};
	my $cmaskFlag = $cmpGenomeInfo{$cgenome}->{'mask'};
	#$cmpGenomeInfo{$cmpGenomes[$z]} = {genome=>$cgenome,dir=>$cgenomeDir,custom=>$ccustomGenome,mask=>$cmaskFlag};
	

	my $extractOption = "";
	if ($cmaskFlag) {
		$extractOption = "-mask";
	}

	print STDERR "\n\tComparing motifs in genome $cgenome\n\n";

	print STDERR "\tLifting over peak positions for $cgenome\n";
	`convertCoordinates.pl "$liftOver" "$posFile" "$tmpfile" -type peaks -p $maxCPUs`;
	print STDERR "\tExtracting sequence for liftover positions in $cgenome\n";
	`homerTools extract "$tmpfile" "$cgenomeDir" $extractOption > "$tmpfile2"`;


	# compare sequences using blastn.....................
	if ($skipBlastn == 0) {
		print STDERR "\tChecking genome peak alignments\n";
		my %seq1 = ();
		my %seq2 = ();
		open IN, $seqFile;
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			$seq1{$line[0]}=$line[1];
		}
		close IN;
		open IN, $tmpfile2;
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			$seq2{$line[0]}=$line[1];
			}
		close IN;
		my $counter = 0;
		foreach(keys %peaks) {
			$counter++;
			if ($counter % 100 ==0) {
				print STDERR "\t\t$counter\n";
			}
			my $hit = $_;
			my $s1 = '';
			my $s2 = '';
			if (!exists($seq1{$hit}) || !exists($seq2{$hit})) {
				$peaks{$hit}->{'gComp'}->[$z]->{'pid'} = "NA";
					$peaks{$hit}->{'gComp'}->[$z]->{'paln'} = 0;
				next;
			}
			my $length = $peaks{$hit}->{'e'} - $peaks{$hit}->{'s'};
			$length = 1 if ($length < 1);
			$s1 = $seq1{$hit};
			$s2 = $seq2{$hit};
			#print STDERR "s1=$s1\ns2=$s2\n";
			open OUT, ">$tmpfile5";
			print OUT ">$hit\n$s1\n";
			close OUT;
			open OUT, ">$tmpfile6";
			print OUT ">$hit\n$s2\n";
			close OUT;
			`blastn -query "$tmpfile5" -subject "$tmpfile6" -outfmt 6 > $tmpfile7`;
			open IN, $tmpfile7;
			while (<IN>) {
				chomp;
				my @line = split /\t/;
				my $pidentity = $line[2];
				my $aln  = $line[3]/$length;
				$aln = 1 if ($aln > 1);
				$peaks{$hit}->{'gComp'}->[$z]->{'pid'} = $pidentity;
				$peaks{$hit}->{'gComp'}->[$z]->{'paln'} = $aln;
				last;
			}
			close IN;
			`rm $tmpfile5 $tmpfile6 $tmpfile7`;
		}
	}

	my %cpeaks = ();
	open IN, $tmpfile;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		my @line = split /\t/;
		my $hit = $line[0];
		my $chr = $line[1];
		my $start = $line[2];
		next if ($start =~ /^[^\d\-\.]/);
		my $end = $line[3];
		my $peakSize = $end-$start+1;
		if ($updateSize == 0) {
			$size = $peakSize;
		}
		my $direction = '+';
		if ($line[4] eq '-' || $line[4] eq '1') {
			$direction = '-';
		}
		my $value = 0;
		my $fdr = 'NA';
		if (@line > 5) {
			$value = $line[5];
		}
		my %a = ();
		my %b = ();
		my %c = ();

		$cpeaks{$hit} = {tss=>'NA', tssDist=>'NA',m=>\%a, p=>\%c, cons=>'NA', size=>$peakSize, tssUG=>'NA',
				s=>$start, e=>$end, d=>$direction, v=>$value, c=>$chr, t=>\%b, gc=>'NA',cpg=>'NA',
				centerDist=>1e10,centerDir=>1e10,centerScore=>-10,fdr=>$fdr,ann=>'NA',
				fullann=>'NA'};

		if (!exists($peaks{$hit})) {
			print STDERR "!!! Something is wrong - can't find $hit in original peak file!\n";
			exit;
		}
		#let us know that this peak was 'liftoverable' - should add more stats like # mutations;
		$peaks{$hit}->{'gComp'}->[$z]->{'map'} = 1;
	}
	close IN;



	print STDERR "\tChecking motifs across genomes\n";

	my %cmask = ();
	for (my $i=0;$i<@motifFiles;$i++) {
		my $mfile = $motifFiles[$i];
		#print STDERR "$mfile\n";
		my $options = '';
		my $code = $filterMotifFiles{$mfile};

		if ($homer2Flag) {
			$options .= " -strand +" if ($revoppFlag == 0);
			$options .= " -p $maxCPUs";
			if ($mscoreFlag) {
				`homer2 find -s "$tmpfile2" -offset $offset -m "$mfile" -mscore $options > "$tmpfile3"`;
			} else {
				`homer2 find -s "$tmpfile2" -offset $offset -m "$mfile" $options > "$tmpfile3"`;
			}
		} else {
			$options .= " -norevopp" if ($revoppFlag == 0);
			if ($mscoreFlag) {
				`homer -s "$tmpfile2" -a GENESCORE -offset $offset -m "$mfile" $options > "$tmpfile3"`;
			} else {
				`homer -s "$tmpfile2" -a FIND -offset $offset -m "$mfile" $options > "$tmpfile3"`;
			}
		}

		open IN, $tmpfile3;
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			my $hit = $line[0];
			if (!exists($cpeaks{$hit})) {
				print STDERR "!!!!!!! problem can't find $hit in cpeaks\n";
			}
			my $pos = $line[1];
			my $seq = $line[2];
			my $con = $line[3];
			my $dir = $line[4];
			my $mname = $line[5];
			my $score = 0;
				if ($homer2Flag) {
				$mname = $line[3];
				$score = $line[5];
				$con = 0 ;
				if ($dir eq '-' || $dir eq '1') {
					$pos -= length($seq)-1;
				}
			} else {
				$score = $line[6];
			}
			my $midpoint = $pos + floor(length($seq)/2);
	
			if ($code eq '-fm') {
				if (!exists($cmask{$hit})) {
						my %a = ();
					$cmask{$hit} = \%a;
				}
				for (my $j=$pos;$j<=$pos+length($seq);$j++) {
					$cmask{$hit}->{$j}=1;
				}
				next;
			} elsif ($code eq '-m') {
				my $bad = 0;
				if (exists($cmask{$hit})) {
					for (my $j=$pos;$j<=$pos+length($seq);$j++) {
						if (exists($cmask{$hit}->{$j})) {
							$bad = 1;
							last;
						}
					}
				}
				if ($bad) {
					next;
				}
			}

			my $pChr = $cpeaks{$hit}->{'c'};
			my $pStart = $cpeaks{$hit}->{'s'};
			my $pDir = $cpeaks{$hit}->{'d'};
			my $bedStart = $pStart + $pos - $offset-1;
			my $bedEnd = $bedStart + length($seq);
			if ($pDir eq '1' || $pDir eq '-') {
				$bedEnd = $cpeaks{$hit}->{'e'} - $pos + $offset;
				$bedStart = $bedEnd - length($seq);
			}
			my $bedDir = '+';
			if ($dir eq '-' || $dir eq '1') {
				$bedDir = '-';
				if ($pDir eq '-' || $pDir eq '1') {
					$bedDir = '+';
				}
			} else {
				if ($pDir eq '-' || $pDir eq '1') {
					$bedDir = '-';
				}
			}

			if (!exists($cpeaks{$hit}->{'m'}->{$mname})) {
				my @a = ();
				$cpeaks{$hit}->{'m'}->{$mname} = \@a;
			}		
			my @compData=();
			if ($cmpGenomeFlag) {
				my $compData = {map=>0,score=>0,indel=>0,var=>0,m=>0};
				push(@compData, $compData);
			}


			my $m = {p=>$pos, s=>$seq, c=>$con, d=>$dir, valid=>1,score=>$score,mp=>$midpoint,
					gc=>$pChr,gs=>$bedStart,ge=>$bedEnd,gd=>$bedDir,gComp=>\@compData};
			push(@{$cpeaks{$hit}->{'m'}->{$mname}}, $m);

			if ($removeCloseMotifs != 0) {
				my $num = scalar(@{$cpeaks{$hit}->{'m'}->{$mname}});
				for (my $i=0;$i<$num;$i++) {
					next if ($cpeaks{$hit}->{'m'}->{$mname}->[$i]->{'valid'}==0);
					my $p1 = $cpeaks{$hit}->{'m'}->{$mname}->[$i]->{'p'};
					my $d1 = $cpeaks{$hit}->{'m'}->{$mname}->[$i]->{'d'};
					for (my $j=0;$j<$num;$j++) {
						next if ($j==$i);
						my $p2 = $cpeaks{$hit}->{'m'}->{$mname}->[$j]->{'p'};
						my $d2 = $cpeaks{$hit}->{'m'}->{$mname}->[$j]->{'d'};
						if (abs($p1-$p2) < $rmMotifThresh) {
							if (($removeCloseMotifs == 1 && $d2 eq $d1) 
									|| ($removeCloseMotifs == -1 && $d2 ne $d1)) {
								$cpeaks{$hit}->{'m'}->{$mname}->[$i]->{'valid'}=0 if ($removeRevoppMotifs == 0);
								$cpeaks{$hit}->{'m'}->{$mname}->[$j]->{'valid'}=0;
							}
						}
					}
				}
			}
		}
		close IN;
		`rm "$tmpfile3"`;
	}



	open MBED, ">$tmpfile";
	foreach(keys %peaks) {
		my $hit = $_;
		foreach(keys %{$peaks{$hit}->{'m'}}) {
			my $mname= $_;
			next if (!exists($peaks{$hit}->{'m'}->{$mname}));
			my $i=0;
			foreach(@{$peaks{$hit}->{'m'}->{$mname}}) {
				my $m = $_;
				my $seq = $m->{'s'};
				my $c = $m->{'gc'};
				my $s = $m->{'gs'};
				my $e = $m->{'ge'};
				my $d = $m->{'gd'};
				my $ss = floor($m->{'score'});
				my $name = $hit . "||" . $mname . "||" . $i . "||" . $m->{'p'} . "||" . $m->{'d'} . "||" . $seq;
				print MBED "$c\t$s\t$e\t$name\t$ss\t$d\n";
				$i++;
			}
		}
	}
	close MBED;

	print STDERR "\tLifting over motif positions for $cgenome\n";
	print STDERR "`convertCoordinates.pl $liftOver $tmpfile $tmpfile2 -type bed -p $maxCPUs`;\n";
	`convertCoordinates.pl "$liftOver" "$tmpfile" "$tmpfile2" -type bed -p $maxCPUs`;
#exit;
	print STDERR "\tExtracting sequence for liftover positions in $cgenome\n";
	`homerTools extract "$tmpfile2" "$cgenomeDir" $extractOption > "$tmpfile3"`;

	open IN, $tmpfile3;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $sname = $line[0];
		my $seq = $line[1];

		my @info = split /\|\|/, $sname;
		next if (@info < 4);
		my $hit = $info[0];
		my $mname = $info[1];
		my $mindex = $info[2];
		my $p = $info[3];
		my $strand = $info[4];
		my $ogseq = $info[5];

		if (!exists($peaks{$hit})) {
			print STDERR "!! Problem - couldn't find peak name...\n";
			exit;
		}
		if (!exists($peaks{$hit}->{'m'}->{$mname})) {
			print STDERR "!! Couldn't find $mname in peak $hit\n";
			exit;
		}

		$peaks{$hit}->{'m'}->{$mname}->[$mindex]->{'gComp'}->[$z]->{'map'}=1;
		
		my $rv = $peaks{$hit}->{'m'}->{$mname}->[$mindex]->{'d'};
		if ($rv eq '-') {
			$seq = HomerConfig::revopp($seq);
		}
		#print STDERR ">$hit $rv\n\t$ogseq\n\t$seq\n";
		if ($ogseq eq $seq) {
			$peaks{$hit}->{'m'}->{$mname}->[$mindex]->{'gComp'}->[$z]->{'m'}=1;
		} else {
			$peaks{$hit}->{'m'}->{$mname}->[$mindex]->{'gComp'}->[$z]->{'m'}=0;
			if (length($ogseq) != length($seq)) {
				$peaks{$hit}->{'m'}->{$mname}->[$mindex]->{'gComp'}->[$z]->{'indel'}=1;
			} else {
				my @subs = getSubs($ogseq,$seq);
				foreach(@subs) {
					my $p1 = $_->[0];
					my $n1 = $_->[1];
					next if (!exists($alpha2index{$n1}));
					my $p2 = $alpha2index{$n1};
					$motifSubMatrix{$mname}->[$p1]->[$p2]+=1;
				}
			}
		}
		#my $compData = {map=>0,score=>0,indel=>0,var=>0};
	}
	close IN;


	open MBED, ">$tmpfile";
	foreach(keys %cpeaks) {
		my $hit = $_;
		foreach(keys %{$cpeaks{$hit}->{'m'}}) {
			my $mname= $_;
			next if (!exists($cpeaks{$hit}->{'m'}->{$mname}));
			my $i=0;
			foreach(@{$cpeaks{$hit}->{'m'}->{$mname}}) {
				my $m = $_;
				my $seq = $m->{'s'};
				my $c = $m->{'gc'};
				my $s = $m->{'gs'};
				my $e = $m->{'ge'};
				my $d = $m->{'gd'};
				my $ss = floor($m->{'score'});
				my $name = $hit . "||" . $mname . "||" . $i . "||" . $m->{'p'} . "||" . $m->{'d'} . "||" . $seq;
				print MBED "$c\t$s\t$e\t$name\t$ss\t$d\n";
				$i++;
			}
		}
	}
	close MBED;

	print STDERR "\tLifting over motif positions from $cgenome back to $genome\n";
	print STDERR "`convertCoordinates.pl $revliftOver $tmpfile $tmpfile2 -type bed -p $maxCPUs`;\n";
	`convertCoordinates.pl "$revliftOver" "$tmpfile" "$tmpfile2" -type bed -p $maxCPUs`;
#exit;
	print STDERR "\tExtracting sequence for liftover positions in $cgenome\n";
	`homerTools extract "$tmpfile2" "$genomeDir" $extractOption > "$tmpfile3"`;

	open IN, $tmpfile3;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $sname = $line[0];
		my $seq = $line[1];

		my @info = split /\|\|/, $sname;
		next if (@info < 4);
		my $hit = $info[0];
		my $mname = $info[1];
		my $mindex = $info[2];
		my $p = $info[3];
		my $strand = $info[4];
		my $ogseq = $info[5];

		if (!exists($cpeaks{$hit})) {
			print STDERR "!! Problem - couldn't find peak name...\n";
			exit;
		}
		if (!exists($cpeaks{$hit}->{'m'}->{$mname})) {
			print STDERR "!! Couldn't find $mname in peak $hit\n";
			exit;
		}

		$cpeaks{$hit}->{'m'}->{$mname}->[$mindex]->{'gComp'}->[0]->{'map'}=1;
		
		my $rv = $cpeaks{$hit}->{'m'}->{$mname}->[$mindex]->{'d'};
		if ($rv eq '-') {
			$seq = HomerConfig::revopp($seq);
		}
		#print STDERR ">$hit $rv\n\t$ogseq\n\t$seq\n";
		if ($ogseq eq $seq) {
			$cpeaks{$hit}->{'m'}->{$mname}->[$mindex]->{'gComp'}->[0]->{'m'}=1;
		} else {
			$cpeaks{$hit}->{'m'}->{$mname}->[$mindex]->{'gComp'}->[0]->{'m'}=0;
			if (length($ogseq) != length($seq)) {
				$cpeaks{$hit}->{'m'}->{$mname}->[$mindex]->{'gComp'}->[0]->{'indel'}=1;
			} else {
				my @subs = getSubs($ogseq,$seq);
				foreach(@subs) {
					my $p1 = $_->[0];
					my $n1 = $_->[1];
					next if (!exists($alpha2index{$n1}));
					my $p2 = $alpha2index{$n1};
					$motifSubMatrix{$mname}->[$p1]->[$p2]+=1;
				}
			}
		}
		#my $compData = {map=>0,score=>0,indel=>0,var=>0};
	}
	push(@cpeaks, \%cpeaks);
}

if ($cmpGenomeFlag) {
	open OUT, ">submatrix.motif";
	foreach(keys %motifSubMatrix) {
		my $mname = $_;
		print OUT ">$mname\t$mname\t0\n";
		foreach(@{$motifSubMatrix{$mname}}) {
			print OUT "$_->[0]";
			for (my $i=1;$i<4;$i++) {
				print OUT "\t$_->[$i]";
			}
			print OUT "\n";
		}
	}
	close OUT;
}

sub getSubs {
	my ($s1,$s2) = @_;
	my @subs = ();
	for (my $i=0;$i<length($s1);$i++) {
		my $n1 = substr($s1,$i,1);
		my $n2 = substr($s2,$i,1);
		if ($n1 ne $n2) {
			push(@subs, [$i,$n2]);
		}
	}
	return @subs;
}
# End cmpGenomes section.................................


# find nearest Peaks
for (my $i=0;$i<@peakFiles;$i++) {
	my $peakFile = $peakFiles[$i];
	print STDERR "\tFinding nearby peaks in $peakFile\n";

	`bed2pos.pl "$peakFile" -check -unique > "$tmpfile"`;
	`annotateRelativePosition.pl "$posFile", "$tmpfile", 0 > "$tmpfile2"`;

	open IN, $tmpfile2;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $hit = $line[0];
		my $dist = $line[2];
		my $peakID = $line[1];
		my $peakCount = $line[4];
		my $relStrand = $line[5];
		next if (!exists($peaks{$hit}));

		my %peakDist = ();
		if (@line > 6) {
			my @pp = split /\,/,$line[6];
			foreach(@pp) {
				next if ($_ eq '');
				my @a = split /\=/;
				my $pp = $a[0];
				my @b = split /\|/,$a[1];
				my $v = 0;
				foreach(@b) {
					if ($_ ne 'NA') {
						$v += $_;
					}
				}
				$peakDist{$pp} = $v;
			}
		}

		$peaks{$hit}->{'p'}->{$peakFile} = {d=>$dist,id=>$peakID,s=>$peakCount,rs=>$relStrand,pd=>\%peakDist};
	}
	close IN;
	`rm "$tmpfile" "$tmpfile2"`;
}

if ($matrixPrefix ne '') {

	my @peaks = keys %peaks;
	my $numPeaks = scalar(@peaks);
	open OUT, ">$matrixPrefix.stats.txt";
	print OUT "Motif1:Motif2\tTotalPeaks\tPeaks with Motif2\tPeaks with Motif1"
				. "\tPeaks with Both\tExpected Overlap\tObserved/Expected\tlogPvalue\n";

	print STDERR "\tCalculating motif co-occurrence matrix\n";
	my %matrix=();
	my %matrixRatio=();
	my %matrixCount=();
	my %numSites = ();

	for (my $i=0;$i<@motifNames;$i++) {
		for (my $j=$i;$j<@motifNames;$j++) {
			my $n1 = 0;
			my $n2 = 0;
			my $n = 0;
			my $nx = 0;
			my $ny = 0;
			my $nz = 0;
			my $nmotifs1 = 0;
			my $nmotifs2 = 0;

			my $curMinDist = $matrixMinDist;;
			if ($i==$j && $curMinDist < 4) {
				# if checking against self, lets always eliminate double counting
				$curMinDist = 4;
			}
			foreach(@peaks) {
				my $pname = $_;
				my $p1 = 0;
				my $p2 = 0;
				my $nm1 = 0;
				my $nm2 = 0;
				if (exists($peaks{$pname}->{'m'}->{$motifNames[$i]})) {
					$p1=1;
					$n1++;
					$nm1=1;
					# this part counts how many motifs are present in the peak (> than curMinDist apart)
					my $narray = scalar(@{$peaks{$pname}->{'m'}->{$motifNames[$i]}});
					for (my $k=0;$k<$narray-1;$k++) {
						my $pos1 = $peaks{$pname}->{'m'}->{$motifNames[$i]}->[$k]->{'mp'};
						my $z = $k+1;
						while ($z < $narray && abs($pos1 -
									$peaks{$pname}->{'m'}->{$motifNames[$i]}->[$z]->{'mp'}) < $curMinDist) {
							$z++;
						}
						if ($z < $narray) {
							$nm1++;
							$k = $z;
						}
					}
					$nmotifs1 += $nm1;
					#print STDERR "$pname $narray $nm1\n";
				}
				if (exists($peaks{$pname}->{'m'}->{$motifNames[$j]})) {
					$p2=1;
					$n2++;
					$nm2=1;
					my $narray = scalar(@{$peaks{$pname}->{'m'}->{$motifNames[$j]}});
					for (my $k=0;$k<$narray-1;$k++) {
						my $pos1 = $peaks{$pname}->{'m'}->{$motifNames[$j]}->[$k]->{'mp'};
						my $z = $k+1;
						while ($z < $narray && abs($pos1 -
									$peaks{$pname}->{'m'}->{$motifNames[$j]}->[$z]->{'mp'}) < $curMinDist) {
							$z++;
						}
						if ($z < $narray) {
							$nm2++;
							$k = $z;
						}
					}
					$nmotifs2 += $nm2;
				}
				if ($p1==1 && $p2 == 1) {
					my $hit = 0;
					foreach(@{$peaks{$pname}->{'m'}->{$motifNames[$i]}}) {
						my $m1 = $_;
						foreach(@{$peaks{$pname}->{'m'}->{$motifNames[$j]}}) {
							my $m2 = $_;
							my $d = abs($m1->{'mp'} - $m2->{'mp'});
							if ($d >= $curMinDist && $d <= $matrixMaxDist) {
								$hit = 1;
							}
						}
					}
					if ($hit == 1)  {
						$n++;
					}
				}
			}


			my $expectedOverlap = 0.0;
			my $ratio = 1.0;
			my $logp = 0.0;

			# in theory this is only for when there's no matrixMin/MaxDist, but for now it will do
			if (1) { #$matrixMinDist <= 0 && $matrixMaxDist > 1e9) {
				$expectedOverlap = $n1*$n2/$numPeaks;
				if ($i==$j) {
					my $avgMotifPerPeak = $nmotifs1/$numPeaks;
					#use binomial distribution to apprixmate how many peaks have 2 or more motifs present
					#in them by chance given the total number of motifs in the whole set of peaks
					my $ptotal = exp(Statistics::logbinomial($avgSize,2,$avgMotifPerPeak/$avgSize,$avgSize*$numPeaks));
					$expectedOverlap = $ptotal*$numPeaks;
				}
			}
			my $expected = $expectedOverlap;
			if ($expectedOverlap < 0.5) {
				$expectedOverlap = 0.5;
			}
			$ratio = $n/$expectedOverlap;
			my $backRatio = $expectedOverlap/$numPeaks;

			if ($ratio >= 1.0) {
				$logp = Statistics::logbinomial($numPeaks, $n, $backRatio, $numPeaks);
			} else {
				$logp = -1*Statistics::ilogbinomial($numPeaks, $n, $backRatio, $numPeaks);
			}


			if (!exists($matrix{$motifNames[$i]})) {
				my %a = ();	
				my %b = ();	
				my %c = ();	
				$matrix{$motifNames[$i]}=\%a;
				$matrixRatio{$motifNames[$i]}=\%b;
				$matrixCount{$motifNames[$i]}=\%c;
			}
			if (!exists($matrix{$motifNames[$j]})) {
				my %a = ();	
				my %b = ();	
				my %c = ();	
				$matrix{$motifNames[$j]}=\%a;
				$matrixRatio{$motifNames[$j]}=\%b;
				$matrixCount{$motifNames[$j]}=\%c;
			}
			$matrix{$motifNames[$i]}->{$motifNames[$j]} = $logp;
			$matrix{$motifNames[$j]}->{$motifNames[$i]} = $logp;
			$matrixRatio{$motifNames[$i]}->{$motifNames[$j]} = $ratio;
			$matrixRatio{$motifNames[$j]}->{$motifNames[$i]} = $ratio;
			$matrixCount{$motifNames[$i]}->{$motifNames[$j]} = $n;
			$matrixCount{$motifNames[$j]}->{$motifNames[$i]} = $n;
			$numSites{$motifNames[$i]} = $n1;
			$numSites{$motifNames[$j]} = $n2;

			print OUT "$motifNames[$i]:$motifNames[$j]\t$numPeaks\t$n1\t$n2\t$n\t$expected\t$ratio\t$logp\n";
		}
	}
	close OUT;


	open OUT, ">$matrixPrefix.logPvalue.matrix.txt";
	open OUT2, ">$matrixPrefix.ratio.matrix.txt";
	open OUT3, ">$matrixPrefix.count.matrix.txt";
	print OUT "Motif Name (# sites) [values are natural log, + values for divergence]";
	print OUT2 "Motif Name (# sites) [observed/expected]";
	print OUT3 "Motif Name (# sites) [Overlap Counts]";
	foreach(@motifNames) {
		print OUT "\t$_ ($numSites{$_})";
		print OUT2 "\t$_ ($numSites{$_})";
		print OUT3 "\t$_ ($numSites{$_})";
	}
	print OUT "\n";
	print OUT2 "\n";
	print OUT3 "\n";
	for (my $i=0;$i<@motifNames;$i++) {
		print OUT "$motifNames[$i] ($numSites{$motifNames[$i]})";
		print OUT2 "$motifNames[$i] ($numSites{$motifNames[$i]})";
		print OUT3 "$motifNames[$i] ($numSites{$motifNames[$i]})";
		for (my $j=0;$j<@motifNames;$j++) {
			my $logp = 0;
			$logp = $matrix{$motifNames[$i]}->{$motifNames[$j]};
			print OUT "\t$logp";
			my $logRatio = $matrixRatio{$motifNames[$i]}->{$motifNames[$j]};
			print OUT2 "\t$logRatio";
			my $v = $matrixCount{$motifNames[$i]}->{$motifNames[$j]};
			print OUT3 "\t$v";
		}
		print OUT "\n";
		print OUT2 "\n";
		print OUT3 "\n";
	}

	close OUT;
	close OUT2;
	close OUT3;
}



######################################################################
########  Histogram Mode..................  ##########################
######################################################################
if ($histBinSize > 0) {
	if ($size eq 'given') {
		print STDERR "\tCompiling per % Histograms...\n";
	} else {
		print STDERR "\tCompiling per bp Histograms...\n";
	}

	my %histograms = ();
	my @histogramNames = ();

	print STDERR "\tFinding Tags in Peaks from each directory...\n";
	my %fragLengths = ();
	for (my $i=0;$i<@tagDirs;$i++) {
		my $dir = $tagDirs[$i];
		my ($t, $p, $flen, $slen) = HomerConfig::readTagInfo($dir,$init2One);

		$tagTotals{$dir} = $t;
		if ($fragLength eq 'auto') {
			if ($flen eq 'NA') {
				print STDERR "\tDefault fragment length set to 150\n";
				$flen = 150;
			} elsif ($flen < 42) {
				print STDERR "\t!warning, ChIP-Fragment length for $dir seems short ($flen)\n";
			}
			$fragLengths{$dir} = $flen;
		} else {
			$fragLengths{$dir} = $fragLength;
		}
		$normFactors{$dir} = 1;
	}

	if ($adjustFlag == 1) {
		if ($normValue < 1) {
			my $nn = 0;
			foreach(@tagDirs) {
				$normValue += $tagTotals{$_};
				$nn++;
			}
			$normValue /= $nn if ($nn > 0);
		}
		foreach(@tagDirs) {
			my $total = $tagTotals{$_};
			next if ($total < 1);
			my $ratio = $normValue / $total;
			print STDERR "\tRatio for $_ : $ratio\n";
			my $vv= sprintf("%.2f",$ratio);
			$normFactors{$_} = $vv;
		}
	}

	my %ghistData = ();
	my $emptyRow = '';
	my @allDirs = (@peakFiles, @motifNames, @tagDirs, @bedGraphFiles, @wigFiles);
	if ($ghistFlag == 1) {

		# should add motif, wig, and peak files to this at some point
		print STDERR "\n\tOrder of experiments in output file:\n";
		print "Gene";	
		if ($size eq 'given') {
			my $startBin = 0;
			for (my $i=0;$i<@allDirs;$i++) {
				print STDERR "\t\t$allDirs[$i]\n";
				for (my $j=1;$j<=$histBinSize;$j++) {
					my $v = sprintf("%.1f",$j/$histBinSize*100) . '%';
					print "\t$v";
					if ($i==0) {
						$emptyRow .= "\t0";
					}
				}
			}
			print "\n";
		} else {
			my $startBin = -1*floor($halfSize/$histBinSize+0.5)*$histBinSize;
			my $endBin = floor($halfSize/$histBinSize+0.5)*$histBinSize;
			for (my $j=0;$j<@allDirs;$j++) {
				print STDERR "\t\t$allDirs[$j]\n";
				for (my $i=$startBin;$i<=$endBin;$i+=$histBinSize) {
					print "\t$i";
					if ($j==0) {
						$emptyRow .= "\t0";
					}
				}
			}
			print "\n";
		}
	}


	foreach(@peakFiles) {
		my $peakFile = $_;
		my %total = ();
		foreach(@peakOrder) {
		#foreach(keys %peaks) {
			my $peakID = $_;
			my %currentPeak = ();
			if (exists($peaks{$peakID}->{'p'}->{$peakFile})) {
				if ($newPeakHistogramFlag) {
					foreach(keys %{$peaks{$peakID}->{'p'}->{$peakFile}->{'pd'}}) {
						my $d = $_;
						my $v = $peaks{$peakID}->{'p'}->{$peakFile}->{'pd'}->{$d};
						if ($size eq 'given') {
							my $peakWidth = $peaks{$peakID}->{'e'}-$peaks{$peakID}->{'s'};
							$v /= ($peakWidth/$histBinSize) if ($peakWidth > 0);
							$d = $d+($peakWidth/2);
							if ($d >= 0 && $d <= $peakWidth) {
								my $binValue = floor(($d/($peakWidth))*$histBinSize);
								$total{$binValue}+=$v;
								$currentPeak{$binValue}+=$v;
							}
						} else {
							my $binValue = floor($d/$histBinSize+0.5)*$histBinSize;
							$total{$binValue}+=$v;
							$currentPeak{$binValue}+=$v;
						}
					}

				} else {
					my $d = $peaks{$peakID}->{'p'}->{$peakFile}->{'d'};
					if ($size eq 'given') {
						my $peakWidth = $peaks{$peakID}->{'e'}-$peaks{$peakID}->{'s'};
						$d = $d+($peakWidth/2);
						if ($d >= 0 && $d <= $peakWidth) {
							my $binValue = floor(($d/($peakWidth))*$histBinSize);
							$total{$binValue}++;
							$currentPeak{$binValue}++;
						}
					} else {
						my $binValue = floor($d/$histBinSize+0.5)*$histBinSize;
						$total{$binValue}++;
						$currentPeak{$binValue}++;
					}
				}
			}

			if ($ghistFlag == 1) {
				my $startBin = -1*floor($halfSize/$histBinSize+0.5)*$histBinSize;
				my $endBin = floor($halfSize/$histBinSize+0.5)*$histBinSize;
				my $incSize = $histBinSize;
				if ($size eq 'given') {
					$startBin = 1;
					$endBin = $histBinSize;
					$incSize = 1;
				}
				my $ghistStr = "";
				for (my $b=$startBin;$b<=$endBin;$b+=$incSize) {
					my $v = 0;
					if (exists($currentPeak{$b})) {
						$v = $currentPeak{$b};
					}
					$ghistStr .= "\t$v";
				}
				if (!exists($ghistData{$peakID})) {
					my %c = ();
					$ghistData{$peakID} = \%c;
				}
				$ghistData{$peakID}->{$peakFile} = $ghistStr;
			}
		}
		push(@histogramNames, $peakFile);
		$histograms{$peakFile}=\%total;
	}

	foreach(@motifNames) {
		my $mname = $_;
		my %p5motifs = ();
		my %p3motifs = ();
		my %total = ();
		foreach(@peakOrder) {
		#foreach(keys %peaks) {
			my $peakID = $_;
			my $curPeak = $peaks{$peakID};
			my $peakWidth = $curPeak->{'e'}-$curPeak->{'s'};
			my %peakTotal = ();
			if (exists($curPeak->{'m'}->{$mname})) {
				foreach(@{$curPeak->{'m'}->{$mname}}) {
					my $pos = $_->{'p'};
					my $valid = $_->{'valid'};
					next if ($valid == 0);
					my $dir = $_->{'d'};
					my $seqLen = length($_->{'s'});
					if ($dir eq '+' || $dir eq "0") {
						$dir = 0;
					} else {
						$dir = 1;
						$pos += $seqLen-1;
					}
					my $midPoint = 0;
					if ($size eq 'given') {
						my $v = 1;
						$v /= ($peakWidth/$histBinSize) if ($peakWidth > 0);
						if ($pos >= 0 && $pos <= $peakWidth) {
							my $binValue = floor($pos/$peakWidth*$histBinSize);
							if ($dir == 0) {
								$p5motifs{$binValue}+=$v;
							} else {
								$p3motifs{$binValue}+=$v;
							}
						}
						if ($dir == 0) {
							$midPoint = $pos + floor($seqLen/2);
						} else {
							$midPoint = $pos - floor($seqLen/2);
						}
						if ($midPoint >= 0 && $midPoint <= $peakWidth) {
							my $binValue = floor($midPoint/$peakWidth*$histBinSize);
							$total{$binValue}+=$v;
							$peakTotal{$binValue}+=$v;
						}
					} else {
						my $binValue = floor($pos/$histBinSize+0.5)*$histBinSize;
						if ($dir == 0) {
							$p5motifs{$binValue}++;
							$midPoint = $pos + floor($seqLen/2);
						} else {
							$p3motifs{$binValue}++;
							$midPoint = $pos - floor($seqLen/2);
						}
						my $midBinValue = floor($midPoint/$histBinSize+0.5)*$histBinSize;
						$total{$midBinValue}++;
						$peakTotal{$midBinValue}++;
					}
				}
			}
			if ($ghistFlag == 1) {
				my $startBin = -1*floor($halfSize/$histBinSize+0.5)*$histBinSize;
				my $endBin = floor($halfSize/$histBinSize+0.5)*$histBinSize;
				my $incSize = $histBinSize;
				if ($size eq 'given') {
					$startBin = 1;
					$endBin = $histBinSize;
					$incSize = 1;
				}
				my $ghistStr = "";
				for (my $b=$startBin;$b<=$endBin;$b+=$incSize) {
					my $v = 0;
					if (exists($peakTotal{$b})) {
						$v = $peakTotal{$b};
					}
					$ghistStr .= "\t$v";
				}
				if (!exists($ghistData{$peakID})) {
					my %c = ();
					$ghistData{$peakID} = \%c;
				}
				$ghistData{$peakID}->{$mname} = $ghistStr;
			}
		}
		my $p5name = $mname . " + sites";
		my $p3name = $mname . " - sites";
		my $totalname = $mname . " total sites";
		push(@histogramNames, $totalname);
		push(@histogramNames, $p5name);
		push(@histogramNames, $p3name);
		$histograms{$totalname}=\%total;
		$histograms{$p5name}=\%p5motifs;
		$histograms{$p3name}=\%p3motifs;

	}

	for (my $i=0;$i<@tagDirs;$i++) {
		my %p5tags = ();
		my %p3tags = ();
		my %coverage = ();
		my %diffMap = ();

		#keeps track of number of tags added for ratio mode
		my %p5tagsN = ();
		my %p3tagsN = ();
		my %coverageN = ();
		my %diffMapN = ();

		my $dir = $tagDirs[$i];
		my $offset = $halfSize*-1+$sizeMove;
		my $sizeRegion = $size;
		my $optStr = " -strand $strandFlag";
		if ($size eq 'given') {
			$optStr .= ' -fixed';
			$sizeRegion = '';
		} else {
			#my $halfRegionSize = floor($sizeRegion/2);
			my $halfRegionSize = floor($sizeRegion/2)+$histBinSize*2+abs($fragLengths{$dir});
			$optStr .= " -offset $offset -start -$halfRegionSize -end $halfRegionSize ";
		}
		
		#print STDERR "`getRelativeTagPositions.pl $posFile,$offset $dir $sizeRegion > $tmpfile`\n";
		#print STDERR "`getPeakTags $posFile $dir -peaktags $optStr > $tmpfile`;\n";
		#`getRelativeTagPositions.pl "$posFile",$offset "$dir" $sizeRegion > "$tmpfile"`;
		#print STDERR "`getPeakTags $posFile $dir -peaktags $optStr > $tmpfile`;\n";
		`getPeakTags "$posFile" "$dir" -peaktags $optStr > "$tmpfile"`;
		#print STDERR "$tmpfile\n";
		#exit;

		my $p5StrandFlag = 1;
		my $p3StrandFlag = 1;
		if ($strandFlag eq '+') {
			$p3StrandFlag = 0;
		} elsif ($strandFlag eq '-') {
			$p5StrandFlag = 0;
		}

		my $normLengthFactor = 1.0;
		if ($normLength > 1e-10) {
			$normLengthFactor = $normLength/$fragLengths{$dir};
		}
		if ($ratioFlag) {
			$normLengthFactor = 1.0;
		}

		print STDERR "\tProcessing tags from $dir\n";
		my $count = 0;
		open IN, $tmpfile;
		while (<IN>) {
			$count++;
			if ($count % 10000 == 0) {
				print STDERR "\t$count\n";
			}
			chomp;
			my @line = split /\t/;
			next if (@line < 2);
			my $peakID = $line[0];

			my $peakWidth = $peaks{$peakID}->{'e'}-$peaks{$peakID}->{'s'};
			my %peakMap = ();
			my %peakMapN = ();
			my @pos = split /\,/,$line[1];
			my $p3total = 0;
			my $p5total = 0;
			my $alltotal = 0;
			if ($histNorm > 0) {
				foreach(@pos) {
					my @pair = split /\=/;
					my $pos = $pair[0];
					my @values = split /\|/, $pair[1];
					if ($init2One > 0) {
						foreach(@values) {
							next if ($_ eq 'NA');
							$_ = $init2One if ($_ > $init2One);
						}
					}
					if ($adjustFlag ==1 && exists($normFactors{$dir})) {
						foreach(@values) {
							next if ($_ eq 'NA');
							$_ *= $normFactors{$dir};
						}
					}
					if ($revoppFlag == 0) {
						$p3total=0;
					}
					if ($values[0] ne 'NA' && $p5StrandFlag) {
						$p5total += $values[0];
						$alltotal += $values[0];
					}
					if ($values[1] ne 'NA' && $p3StrandFlag) {
						$p3total += $values[1];
						$alltotal += $values[1];
					}
				}
				$p5total = $histNorm if ($p5total < $histNorm);
				$p3total = $histNorm if ($p3total < $histNorm);
				$alltotal = $histNorm if ($alltotal < $histNorm);
			} else {
				$p3total = 1;
				$p5total = 1;
				$alltotal = 1;
				if ($size eq 'given') {
					$p3total = $peakWidth/$histBinSize;
					$p5total = $peakWidth/$histBinSize;
					$alltotal = $peakWidth/$histBinSize;
				}
			}
			#print STDERR "$p3total\t$p5total\t$alltotal\n";	

			foreach(@pos) {
				my @pair = split /\=/;
				my $pos = $pair[0];
				my @values = split /\|/, $pair[1];
				if ($init2One > 0) {
					foreach(@values) {
						next if ($_ eq 'NA');
						$_ = $init2One if ($_ > $init2One);
					}
				}
				if ($adjustFlag ==1 && exists($normFactors{$dir})) {
					foreach(@values) {
						next if ($_ eq 'NA');
						$_ *= $normFactors{$dir} if ($ratioFlag == 0);
					}
				}
				if ($revoppFlag == 0) {
					#$p3total=0; ### Why???
				}
				my $binValue = 0;
				if ($size eq 'given') {
					$binValue = floor($pos/$peakWidth*$histBinSize);
				} else {
					$binValue = floor($pos/$histBinSize+0.5)*$histBinSize;
				}

				if ($values[0] ne 'NA' && $p5StrandFlag) {
					if ($ratioFlag == 1) {
						$p5tags{$binValue} += $values[0];
						$p5tagsN{$binValue} ++;
					} else {
						$p5tags{$binValue} += $values[0]/$p5total;
						$p5tagsN{$binValue} += $p5total;
					}

					my $endPos = $pos+$fragLengths{$dir};
					my $startBin = $binValue;
					my $endBin = $binValue;
					if ($size eq 'given') {
						$endBin = floor($endPos/$peakWidth*$histBinSize)+1;
					} else {
						$endBin = floor($endPos/$histBinSize+0.5)*$histBinSize+$histBinSize;
					}
					my $v = $values[0]/$alltotal*$normLengthFactor;
					if ($ratioFlag == 1) {
						$v = $values[0]/$normLengthFactor;
					}
					$diffMap{$startBin} += $v;
					$diffMap{$endBin} -= $v;
					$peakMap{$startBin} += $v;
					$peakMap{$endBin} -= $v;
					$diffMapN{$startBin}++;
					$diffMapN{$endBin}--;
					$peakMapN{$startBin}++;
					$peakMapN{$endBin}--;
				}
				if ($values[1] ne 'NA' && $p3StrandFlag) {
					if ($ratioFlag == 1) {
						$p3tags{$binValue} += $values[1];
						$p3tagsN{$binValue}++;
					} else {
						$p3tags{$binValue} += $values[1]/$p3total;
						$p3tagsN{$binValue}+= $p3total;
					}

					my $endPos = $pos-$fragLengths{$dir};
					my $startBin = $binValue;
					my $endBin = $binValue;
					if ($size eq 'given') {
						$startBin = floor($endPos/$peakWidth*$histBinSize);
						$endBin ++;
					} else {
						$startBin = floor($endPos/$histBinSize+0.5)*$histBinSize;
						$endBin += $histBinSize;
					}
					my $v = $values[1]/$alltotal*$normLengthFactor;
					if ($ratioFlag == 1) {
						$v = $values[1]/$normLengthFactor;
					}
					$diffMap{$startBin} += $v;
					$diffMap{$endBin} -= $v;
					$peakMap{$startBin} += $v;
					$peakMap{$endBin} -= $v;
					$diffMapN{$startBin}++;
					$diffMapN{$endBin}--;
					$peakMapN{$startBin}++;
					$peakMapN{$endBin}--;
				}
			}
			if ($ghistFlag == 1) {
				my %peakCoverage=();
				my @peakSortBins = sort {$a <=> $b} keys %peakMap;
				my $startBin = $sortBins[0];
				my $endBin = $sortBins[@sortBins-1];
				my $value = 0;
				my $N = 0;
				foreach(@peakSortBins) {
					if (exists($peakMap{$_})) {
						$value += $peakMap{$_};
						$N += $peakMapN{$_};
					}
					my $v = $value;
					if ($ratioFlag) {
						$v /= $N if ($N>0);
					}
					$peakCoverage{$_}=$v;
				}
				$startBin = -1*floor($halfSize/$histBinSize+0.5)*$histBinSize;
				$endBin = floor($halfSize/$histBinSize+0.5)*$histBinSize;
				my $incSize = $histBinSize;
				if ($size eq 'given') {
					$startBin = 1;
					$endBin = $histBinSize;
					$incSize = 1;
				}

				my $last = 0;
				my $ghistStr = "";
				for (my $b=$startBin;$b<=$endBin;$b+=$incSize) {
					my $v = $last;
					if (exists($peakCoverage{$b})) {
						$v = $peakCoverage{$b};
						$last = $v;
					}
					$ghistStr .= "\t$v";
				}
				if (!exists($ghistData{$peakID})) {
					my %c = ();
					$ghistData{$peakID} = \%c;
				}
				$ghistData{$peakID}->{$dir} = $ghistStr;
			}
		}
		close IN;
		`rm "$tmpfile"`;
		#print STDERR "\tDone processing Tags\n";

		if ($ghistFlag == 1) {
			next;
		}

		#build coverage map
		my @sortBins = sort {$a <=> $b} keys %diffMap;
		if (@sortBins > 1) {

			my %tmpCoverage = ();
			my $value = 0;
			my $N = 0;
			foreach(@sortBins) {
				my $b = $_;
				$value += $diffMap{$b};
				$N += $diffMapN{$b};
				my $v = $value;
				if ($ratioFlag == 1) {
					$v /= $N if ($N > 0);
				}
				$tmpCoverage{$b} = $v;
			}
			

			my $startBin = $sortBins[0];
			my $endBin = $sortBins[@sortBins-1];
			my $incSize = $histBinSize;
			my $lastValue = 0;

			if ($size eq 'given') {
				$startBin = 0;
				$endBin = $histBinSize;
				$incSize = 1;
			}
			for (my $b=$startBin;$b<=$endBin;$b+=$incSize) {
				if (exists($tmpCoverage{$b})) {
					$lastValue = $tmpCoverage{$b};
				}
				$coverage{$b} = $lastValue;
				if ($ratioFlag != 1) {
					$coverage{$b} *= $histBinSize;
				}
			}
		}
		if ($ratioFlag==1) {
			# && $size ne 'given') {
			foreach(keys %p5tags) {
				$p5tags{$_} /= $p5tagsN{$_};
			}
			foreach(keys %p3tags) {
				$p3tags{$_} /= $p3tagsN{$_};
			}
		}

		my $p5name = "$dir" . " + Tags";
		my $p3name = "$dir" . " - Tags";
		my $coverageName = "$dir" . " Coverage";
		push(@histogramNames, $coverageName);
		push(@histogramNames, $p5name);
		push(@histogramNames, $p3name);
		$histograms{$p5name} = \%p5tags;
		$histograms{$p3name} = \%p3tags;
		$histograms{$coverageName} = \%coverage;
	}

	my @covFiles = ();
	my @typeCovFiles = ();
	foreach(@bedGraphFiles) {
		push(@covFiles, $_);
		push(@typeCovFiles, "-bedGraph");
	}
	foreach(@wigFiles) {
		push(@covFiles, $_);
		push(@typeCovFiles, "-wig");
	}

	for (my $i=0;$i<@covFiles;$i++) {
		my $bedGraphFile = $covFiles[$i];
		my $type = $typeCovFiles[$i];
		my %p5tags = ();
		my %p3tags = ();
		my %coverage = ();
		my %diffMap = ();

		#keeps track of number of tags added for ratio mode
		my %p5tagsN = ();
		my %p3tagsN = ();
		my %coverageN = ();
		my %diffMapN = ();

		my $dir = $tagDirs[$i];
		my $offset = $halfSize*-1+$sizeMove;
		my $sizeRegion = $size;
		my $optStr = " -strand $strandFlag";
		if ($size eq 'given') {
			$optStr .= ' -fixed';
			$sizeRegion = '';
		} else {
			my $halfRegionSize = floor($sizeRegion/2)+$histBinSize*2;
			$optStr .= " -offset $offset -start -$halfRegionSize -end $halfRegionSize ";
		}
		
		`getPeakTags "$posFile" $type "$bedGraphFile" -peaktags $optStr > "$tmpfile"`;
		#print STDERR "`getPeakTags $posFile -bedGraph $bedGraphFile -peaktags $optStr > $tmpfile`;\n";

		my $count = 0;
		open IN, $tmpfile;
		while (<IN>) {
			$count++;
			if ($count % 10000 == 0) {
				print STDERR "\t$count\n";
			}
			chomp;
			s/\r//g;
			my @line = split /\t/;
			my $peakID = $line[0];
			next if (!exists($peaks{$peakID}));
			my $peakWidth = $peaks{$peakID}->{'e'}-$peaks{$peakID}->{'s'};
			my %peakMap = ();
			my @pos = split /\,/,$line[1];
			my $alltotal = 0;
			my $wigRatio = 1;
			if ($size ne 'given') {
				$wigRatio *= $histBinSize;
			}
			my $p5total = 0;
			my $p3total = 0;
			if ($histNorm > 0) {
				foreach(@pos) {
					my @pair = split /\=/;
					my @p = split /\|/, $pair[0];
					my $value = $pair[1];
					$alltotal += $value*($p[1]-$p[0]);
					$p5total += $value;
					$p3total += $value;
				}
				$alltotal = $histNorm if ($alltotal < $histNorm);
				$p5total = $histNorm if ($p5total < $histNorm);
				$p3total = $histNorm if ($p3total < $histNorm);
			} else {
				$alltotal = 1;
				$p5total = 1;
				$p3total = 1;
				if ($size eq 'given') {
					$p3total = $peakWidth/$histBinSize;
					$p5total = $peakWidth/$histBinSize;
					$alltotal = $peakWidth/$histBinSize;
				}
			}
			foreach(@pos) {
				my @pair = split /\=/;
				my @p = split /\|/, $pair[0];
				my $value = $pair[1];

				my $binValueStart = 0;
				my $binValueEnd = 0;
				if ($size eq 'given') {
					$binValueStart = floor($p[0]/$peakWidth*$histBinSize);
					$binValueEnd = floor(($p[1])/$peakWidth*$histBinSize);
				} else {
					$binValueStart = floor($p[0]/$histBinSize+0.5)*$histBinSize;
					$binValueEnd = floor(($p[1])/$histBinSize+0.5)*$histBinSize;
				}

				$p5tags{$binValueStart} += $value/$p5total;
				$p5tagsN{$binValueStart} += $p5total;
				$p3tags{$binValueEnd} += $value/$p3total;
				$p3tagsN{$binValueEnd} += $p3total;
			
				my $startBin = $binValueStart;
				my $endBin = $binValueEnd;

				$diffMap{$startBin} += $value/$alltotal*$wigRatio;
				$diffMap{$endBin} -= $value/$alltotal*$wigRatio;
				$peakMap{$startBin} += $value/$alltotal*$wigRatio;
				$peakMap{$endBin} -= $value/$alltotal*$wigRatio;
				$diffMapN{$startBin}++;
				$diffMapN{$endBin}--;
				$peakMapN{$startBin}++;
				$peakMapN{$endBin}--;
			}
			if ($ghistFlag == 1) {
				my %peakCoverage=();
				my @peakSortBins = sort {$a <=> $b} keys %peakMap;
				my $startBin = $sortBins[0];
				my $endBin = $sortBins[@sortBins-1];
				my $value = 0;
				my $N = 0;
				foreach(@peakSortBins) {
					if (exists($peakMap{$_})) {
						$value += $peakMap{$_};
						$N += $peakMapN{$_};
					}
					my $v = $value;
					if ($ratioFlag) {
						$v /= $N if ($N>0);
					}
					$peakCoverage{$_}=$v/$histBinSize;
				}
				$startBin = -1*floor($halfSize/$histBinSize+0.5)*$histBinSize;
				$endBin = floor($halfSize/$histBinSize+0.5)*$histBinSize;
				my $incSize = $histBinSize;
				if ($size eq 'given') {
					$startBin = 1;
					$endBin = $histBinSize;
					$incSize = 1;
				}

				my $last = 0;
				my $ghistStr = "";
				for (my $b=$startBin;$b<=$endBin;$b+=$incSize) {
					my $v = $last;
					if (exists($peakCoverage{$b})) {
						$v = $peakCoverage{$b};
						$last = $v;
					}
					$ghistStr .= "\t$v";
				}
				if (!exists($ghistData{$peakID})) {
					my %c = ();
					$ghistData{$peakID} = \%c;
				}
				$ghistData{$peakID}->{$bedGraphFile} = $ghistStr;
			}
		}
		close IN;
		`rm "$tmpfile"`;
		#print STDERR "\tDone processing Tags\n";

		if ($ghistFlag == 1) {
			next;
		}

		#build coverage map
		my @sortBins = sort {$a <=> $b} keys %diffMap;
		if (@sortBins > 1) {

			my %tmpCoverage = ();
			my $value = 0;
			my $N = 0;
			foreach(@sortBins) {
				my $b = $_;
				$value += $diffMap{$b};
				$N += $diffMapN{$b};
				my $v = $value;
				if ($ratioFlag == 1) {
					$v /= $N if ($N > 0);
				}
				$tmpCoverage{$b} = $v;
			}
			

			my $startBin = $sortBins[0];
			my $endBin = $sortBins[@sortBins-1];
			my $incSize = $histBinSize;
			my $lastValue = 0;

			if ($size eq 'given') {
				$startBin = 0;
				$endBin = $histBinSize;
				$incSize = 1;
			}
			for (my $b=$startBin;$b<=$endBin;$b+=$incSize) {
				if (exists($tmpCoverage{$b})) {
					$lastValue = $tmpCoverage{$b};
				}
				$coverage{$b} = $lastValue;
			}
		}
		if ($ratioFlag==1) {
			foreach(keys %p5tags) {
				$p5tags{$_} /= $p5tagsN{$_};
			}
			foreach(keys %p3tags) {
				$p3tags{$_} /= $p3tagsN{$_};
			}
		}

		my $p5name = "$bedGraphFile" . " + Tags";
		my $p3name = "$bedGraphFile" . " - Tags";
		my $coverageName = "$bedGraphFile" . " Coverage";
		push(@histogramNames, $coverageName);
		push(@histogramNames, $p5name);
		push(@histogramNames, $p3name);
		$histograms{$p5name} = \%p5tags;
		$histograms{$p3name} = \%p3tags;
		$histograms{$coverageName} = \%coverage;
	}


	if ($ghistFlag == 1) {
		foreach(@peakOrder) {
			my $peakID = $_;
			next if (!exists($ghistData{$peakID}));
			print "$peakID";
			for (my $i=0;$i<@allDirs;$i++) {
				my $dir = $allDirs[$i];
				if (!exists($ghistData{$peakID}->{$dir})) {
					print $emptyRow;
				} else {
					print $ghistData{$peakID}->{$dir};
				}
			}
			print "\n";
		}
		`rm "$posFile"`;
		deleteFiles();
		exit;
	} else {

		#normalize histograms to per bp numbers
		my $numPeaks = scalar(keys %peaks);
		my $totalBpInBin = $histBinSize * $numPeaks;
		if ($size eq 'given') {
			$totalBpInBin = $numPeaks;
		}
		if ($ratioFlag==1) {
			$totalBpInBin = $histBinSize;
			$totalBpInBin = 1;
		}
		foreach(values %histograms) {
			foreach(values %$_) {
				$_ /= $totalBpInBin;
			}
		}
	}

	if ($seqFlag) {
		my %freqHists = ();
		my %freqCounts = ();
		my @freqNames = ();
		my $limit = 4;
		if ($diFlag == 1) {
			$limit = 24;
		}

		my $option = '';
		if ($size ne 'given') {
			$option = " -maxlen " . ($size+1);
		}
		`homerTools freq "$seqFile" -offset -$halfSize $option > "$tmpfile"`;
		open IN, $tmpfile;
		my $rowCount = 0;
		while (<IN>) {
			$rowCount++;
			chomp;
			my @line = split /\t/;
			if ($rowCount == 1) {
				for (my $i=1;$i<=$limit;$i++) {
					push(@freqNames, $line[$i]);
					my %a = ();
					my %b = ();
					$freqHists{$line[$i]} = \%a;
					$freqCounts{$line[$i]} = \%b;
				}
				next;
			}
			my $pos = $line[0];
			my $binValue = floor($pos/$histBinSize+0.5)*$histBinSize;
			for (my $i=1;$i<=$limit;$i++) {
				$freqHists{$freqNames[$i-1]}->{$binValue}+=$line[$i];
				$freqCounts{$freqNames[$i-1]}->{$binValue}++;
			}
		}
		close IN;
		`rm "$tmpfile"`;
		foreach(@freqNames) {
			my $name = $_;
			foreach(keys %{$freqHists{$name}}) {
				$freqHists{$name}->{$_} /= $freqCounts{$name}->{$_};
			}
			my $endName = "$name frequency";
			push(@histogramNames, $endName);
			$histograms{$endName} = $freqHists{$name};
		}
		`rm "$seqFile"`;
	}

	if ($consFlag) {
		`conservationAverage.pl "$consFile"  -$halfSize > "$tmpfile"`;
		my %cons = ();
		my %consCounts = ();
		open IN, $tmpfile;
		my $rowCount = 0;
		while (<IN>) {
			$rowCount++;
			next if ($rowCount < 2);
			chomp;
			my @line = split /\t/;
			my $pos = $line[0];
			my $binValue = floor($pos/$histBinSize+0.5)*$histBinSize;
			$cons{$binValue} += $line[1];
			$consCounts{$binValue}++;
		}
		close IN;
		`rm "$tmpfile"`;
		foreach(keys %cons) {
			$cons{$_} /= $consCounts{$_};
		}
		push(@histogramNames, "Conservation");
		$histograms{"Conservation"} = \%cons;
		`rm "$consFile"`;
	}
	`rm "$posFile"`;

	print "Distance from Center (cmd=$cmd)";
	foreach(@histogramNames) {
		print "\t$_";
	}
	print "\n";

	my $startBin = -1*floor($halfSize/$histBinSize+0.5)*$histBinSize;
	my $endBin = floor($halfSize/$histBinSize+0.5)*$histBinSize;
	my $incSize = $histBinSize;
	if ($size eq 'given') {
		$startBin = 0;
		$endBin = $histBinSize-1;
		$incSize = 1;
	}
	for (my $i=$startBin;$i<=$endBin;$i+=$incSize) {
		my $v = $i;
		if ($size eq 'given') {
			$v = ($i+1)/$histBinSize;
		}
		print "$v";
		foreach(@histogramNames) {
			if (exists($histograms{$_}->{$i})) {
				print "\t$histograms{$_}->{$i}";
			} else {
				print "\t0";
			}
		}
		print "\n";
	}
	print STDERR "\n";
	deleteFiles();
	exit;
}

######################################################################
########  Annotation mode related code....  ##########################
######################################################################

my %nativeTSSid = ();
# find nearest TSS
if ($tssFlag == 0 && $noGeneFlag == 0) {
	print STDERR "\tFinding Closest TSS...\n";
	my $promoterFile = $genomeDir . "/" . $genome . ".tss";	

	if ($gtfFile ne '') {
		print STDERR "\tProcessing custom annotation file...\n";
		`parseGTF.pl "$gtfFile" tss $gffFlag $gidFlag > "$gtfTSSFile"`;
		$promoterFile = $gtfTSSFile;
		$toDelete{$gtfTSSFile}=1;
	}

	my $promoterOffset = -2000;
	if ($promoter ne 'default') {
		$promoterFile = $config->{'PROMOTERS'}->{$promoter}->{'directory'} . "/$promoter.pos";
		$promoterOffset = $config->{'PROMOTERS'}->{$promoter}->{'start'};
	}
	if ($cpromoter ne '') {
		$promoterFile = $cpromoter;
		#$promoterOffset = -2000;
		$promoterOffset = ""
	}

	my $skipPromoterFile = 0;
	if ($mapFile ne '') {
		open IN, $mapFile;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line= split /\t/;
			if (exists($peaks{$line[0]})) {
				my $dist = 0;
				my $gene = $line[1];
				if (@line > 2) {
					if ($line[2] eq 'interchromosomal') {
						$line[2] = 1e9;
					}
					$dist = $line[2];
				}
				$gene =~ s/\-HOMER.*$//;
				if (exists($peaks{$line[0]}->{'tss'}) && $peaks{$line[0]}->{'tss'} ne 'NA') {
					if ($dist < $peaks{$line[0]}->{'tssDist'}) {
						$peaks{$line[0]}->{'tss'} = $gene;
						$peaks{$line[0]}->{'tssDist'} = $dist;
					}
				} else {
					$peaks{$line[0]}->{'tss'} = $gene;
					$peaks{$line[0]}->{'tssDist'} = $dist;
				}
			}
		}
		close IN;

	} elsif (-f $promoterFile) {
		`annotateRelativePosition.pl "$posFile", "$promoterFile",$promoterOffset  0 > "$tmpfile"`;
		#`cp "$tmpfile" check.tsv`;
		open IN, $tmpfile;
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			my $hit = $line[0];
			my $dist = $line[2];
			my $gene = $line[1];
			my $relStrand = $line[5];
			$gene =~ s/\-HOMER.*$//;
			my $geneDirection = $line[3];
			if ($geneDirection == 0) {
				$dist *= -1;
				#print STDERR "$dist\n";
			}
			if ($peaks{$hit}->{'d'} eq '-') {
				$dist *= -1;
			}
			$peaks{$hit}->{'tss'} = $gene;
			$peaks{$hit}->{'tssDist'} = $dist;
			$peaks{$hit}->{'rstrand'} = $relStrand;
			$nativeTSSid{$gene}='NA';
		}
		close IN;
		`rm "$tmpfile"`;
	} else {
		my $skipPromoterFile = 1;
		print STDERR "\tSkipping TSS assignment (can't find file for genome $genome)\n";
		print STDERR "\t\tCan't find promoterFile $promoterFile\n";
		#print STDERR "!!!! Can't find TSS file for $genome ($promoterFile) !!!\n";
	}

} else {
	foreach(keys %peaks) {
		$peaks{$_}->{'tss'} = $_;
		$peaks{$_}->{'tss'} =~ s/\-HOMER.*$//;
		$peaks{$_}->{'tssDist'} = 0;
		$nativeTSSid{$_} = 'NA';
	}
}


if ($noAnnFlag == 0 && $noGeneFlag == 0) {
	my $annotationFile = $genomeDir . "/" . $genome . ".basic.annotation";

	if ($customAnnotationFile ne '') {
		$annotationFile = $customAnnotationFile;
	} elsif ($gtfFile ne '' && $customAnnotationFile eq '') {
		`parseGTF.pl "$gtfFile" ann $gffFlag $gidFlag > "$tmpfile"`;
		`assignGenomeAnnotation "$tmpfile" "$tmpfile" -prioritize "$tmpfile2"`;
		$annotationFile = $tmpfile2;
	} else {
		if (-f $annotationFile) {
		} else {
			$annotationFile = $genomeDir . "/" . $genome . ".annotation";
		}
	}

	if (-f $annotationFile) {
		my $opt = '';
		if ($annStatFile ne '') {
			$opt = " -stats \"$tmpfile3\"";
		}
		`assignGenomeAnnotation "$posFile" "$annotationFile" -ann "$tmpfile" $opt`;
		#print STDERR "`assignGenomeAnnotation $posFile $annotationFile -ann $tmpfile`\n";
		open IN, $tmpfile;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			next if (!exists($peaks{$line[0]}));
			$line[2] =~ s/\-+\d+$//;
			$line[2] =~ s/\-HOMER\d+$//;
			$peaks{$line[0]}->{'ann'} = $line[2];
		}
		`rm "$tmpfile"`;
	} else {
		print STDERR "\tCould not find basic annotation file ($annotationFile)\n";
	}

	$annotationFile = $genomeDir . "/" . $genome . ".full.annotation";
	if (-f $annotationFile) {
	} else {
		$annotationFile = $genomeDir . "/" . $genome . ".annotation";
	}

	if (-f $annotationFile) {
		print STDERR "\tNOTE: If this part takes more than 2 minutes, there is a good chance\n";
		print STDERR "\t\tyour machine ran out of memory: consider hitting ctrl+C and rerunning\n";
		print STDERR "\t\tthe command with \"-noann\"\n";
		if ($annStatFile eq '') {
			print STDERR "\tTo capture annotation stats in a file, use \"-annStats <filename>\" next time\n";
		}
		my $opt = '';
		if ($annStatFile ne '') {
			$opt = " -stats \"$tmpfile4\"";
		}
		`assignGenomeAnnotation "$posFile" "$annotationFile" -ann "$tmpfile" $opt`;
		open IN, $tmpfile;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			next if (!exists($peaks{$line[0]}));
			$line[2] =~ s/\-+\d+$//;
			$line[2] =~ s/\-HOMER\d+$//;
			$peaks{$line[0]}->{'fullann'} = $line[2];
		}
	} else {
		if ($gtfFile eq '') {
			print STDERR "\tCould not find full/detailed annotation file ($annotationFile)\n";
		}
	}
	if ($annStatFile ne '') {
		`cat "$tmpfile3" "$tmpfile4" > "$annStatFile"`;
	}
	
	`rm -f "$tmpfile" "$tmpfile2" "$tmpfile3" "$tmpfile4"`;
}



if ($goDir ne '' && $noGeneFlag == 0 && $organism ne 'unknown') {
	print STDERR "\tPerforming Gene Ontology Analysis...\n";
	open OUT, ">$tmpfile";
	foreach(keys %nativeTSSid) {
		print OUT "$_\n";
	}
	close OUT;
	`findGO.pl "$tmpfile" $organism "$goDir"`;
	`rm "$tmpfile"`;
}

if ($genomeOntologyDir ne '') {
	print STDERR "\t--------------------------------------\n";
	print STDERR "\tPerforming Genome Ontology Analysis...\n";
	`GenomeOntology.pl "$posFile" $genome "$genomeOntologyDir"`;
}



my @allIndividuals = ();
my %snps = ();
if ($vcfFile ne "") {
	print STDERR "\tExtracting genetic variation information from VCF file...\n";

	my $cppOpt = "";
	if ($size eq 'given') {
		$cppOpt .= " -fixed";
	} else {
		$cppOpt .= " -start -$halfSize -end $halfSize";
	}
	if (scalar(@individuals) > 1) {
		$cppOpt .= " -individuals";
		foreach(@individuals) {
			$cppOpt .= " $_";
		}
	}

	if ($editDistanceFlag) {
	} else {
		$cppOpt .= " -peaksnps";
	}

	#print STDERR "`getPeakTags $posFile $dir $cppOpt > $tmpfile`\n";
	`getPeakTags "$posFile" -vcf "$vcfFile" $cppOpt > "$tmpfile"`;
	
	open IN, $tmpfile;
	my $lineCount = 0;
	while (<IN>) {
		$lineCount++;
		chomp;
		my @line = split /\t/;
		if ($lineCount == 1) {
			last if ($line[1] eq '');
			@allIndividuals = split /\,/,$line[1];
			next;
		}
		my $hit = $line[0];
		next if (!exists($peaks{$hit}));
		next if (@line < 2);
		next if ($line[1] eq '');
		if ($editDistanceFlag) {
			my @pos = split /\,/,$line[1];
			$peaks{$hit}->{'snps'} = \@pos;
		} else {
			$peaks{$hit}->{'snps'} = $line[1];
		}
	}
	close IN;
	`rm "$tmpfile"`;
}
	

print STDERR "\tCounting Tags in Peaks from each directory...\n";
my @newDirs = ();
my %tagTotals = ();
my %normFactors = ();
my @allDirs = ();
my %dirTypes = ();
for (my $i=0;$i<@tagDirs;$i++) {
	push(@allDirs, $tagDirs[$i]);
	$dirTypes{$tagDirs[$i]} = "tagDir";
}
for (my $i=0;$i<@bedGraphFiles;$i++) {
	push(@allDirs, $bedGraphFiles[$i]);
	$dirTypes{$bedGraphFiles[$i]} = "bedGraph";
}
for (my $i=0;$i<@wigFiles;$i++) {
	push(@allDirs, $wigFiles[$i]);
	$dirTypes{$wigFiles[$i]}="wiggle";
}

for (my $i=0;$i<@allDirs;$i++) {
	my $dir = $allDirs[$i];
	my $type = $dirTypes{$dir};

	my $tagTotal = 0;	
	my $tagPosTotal = 0;	
	my $dirFragLength = 0;
	my $dirHalfFragLength = 0;
	my $dirPeakLength = 0;

	if ($type eq 'tagDir') {
		($tagTotal, $tagPosTotal,$dirFragLength,$dirPeakLength) = HomerConfig::readTagInfo($dir,$init2One);
	}

	$tagTotals{$dir} = $tagTotal;
	$normFactors{$dir} = 1;

	my $cppOpt = "";
	$cppOpt .= " -strand $strandFlag ";

	if ($ratioFlag==1) {
		$cppOpt .= " -ratio";
	}
	if ($fragLength eq 'auto') {
		$cppOpt .= " -tagAdjust auto";
	} else {
		$dirHalfFragLength = floor($fragLength/2);
		$cppOpt .= " -tagAdjust $dirHalfFragLength";
	}
	$cppOpt .= " -tbp $init2One";
	my $preOpt = $cppOpt;

	if ($nfrFlag) {
		$cppOpt .= " -nfr -nfrSize $nfrSize";
	}

	if ($size eq 'given') {
		$cppOpt .= " -fixed";
	} else {
		$cppOpt .= " -start -$halfSize -end $halfSize";
	}

	my $inputFile = "\"$dir\"";
	if ($type eq 'bedGraph') {
		$inputFile = "-bedGraph $inputFile";
	} elsif ($type eq 'wiggle') {
		$inputFile = "-wig $inputFile";
	}
	#print STDERR "`getPeakTags $posFile $inputFile $cppOpt > $tmpfile`;\n";
	`getPeakTags "$posFile" $inputFile $cppOpt > "$tmpfile"`;
	
	open IN, $tmpfile;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $hit = $line[0];
		my $v = $line[1];
		$peaks{$hit}->{'t'}->{$dir} = $v;
	}
	close IN;
	`rm "$tmpfile"`;
	push(@newDirs, $dir);

	if ($local > 0 && $size ne 'given') {
		my $factor = sprintf("%.1f", $local/$size);
		my $dirBack = $dir . "-x$factor";


		$preOpt .= " -start -$halfLocal -end $halfLocal";
		`getPeakTags "$posFile" $inputFile $preOpt > "$tmpfile"`;

		open IN, $tmpfile;
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			my $hit = $line[0];
			my $v = $line[1];
			$peaks{$hit}->{'t'}->{$dirBack} = $v;
		}
		close IN;
		`rm "$tmpfile"`;
		$tagTotals{$dirBack} = $tagTotal;
		push(@newDirs, $dirBack);
		$dirTypes{$dirBack} = $type;
	}
}

@allDirs = @newDirs;
if ($ratioFlag == 0 && $adjustFlag == 1 && $nfrFlag == 0) {
	if ($normValue < 1) {
		my $nn = 0;
		foreach(@allDirs) {
			next if ($dirTypes{$_} ne 'tagDir');
			$normValue += $tagTotals{$_};
			$nn++;
		}
		$normValue /= $nn if ($nn > 0);
	}
	my @peaks = keys %peaks;
	foreach(@allDirs) {
		next if ($dirTypes{$_} ne 'tagDir');
		my $total = $tagTotals{$_};
		next if ($total < 1);
		my $ratio = $normValue / $total;
		print STDERR "\tRatio for $_ : $ratio\n";
		my $vv= sprintf("%.2f",$ratio);
		$normFactors{$_} = $vv;
		my $dir = $_;
		foreach(@peaks) {
			if (exists($peaks{$_}->{'t'}->{$dir})) {
				my $v = $peaks{$_}->{'t'}->{$dir} * $ratio;

				if ($logFlag) {
					$v = log(1+$v+rand()*$ratio)/log(2);
				} 
				if ($sqrtFlag) {
					$v = sqrt($v+rand()*$ratio);
				}
				$peaks{$_}->{'t'}->{$dir} = sprintf("%.2f",$v);
			}
		}
	}
}


my %ug2gene = ();
my $convFile = $homeDir . "/data/accession/$organism" . "2gene.tsv";
my @geneData = ();
my @geneDataHeader = (); 
my %acc2gene = ();


print STDERR "\tOrganism: $organism\n";
if ($noGeneFlag == 0 && $organism ne 'unknown') {

	print STDERR "\tLoading Gene Informaiton...\n";
	open IN, $convFile;
	while (<IN>) {
		chomp;
		my @line= split /\t/;
		next if (exists($acc2gene{$line[0]}));
		$acc2gene{$line[0]} = $line[1];
	}
	close IN;

	my $descriptionFile = $homeDir . "/data/accession/$organism.description";
	open IN, $descriptionFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line= split /\t/;
		my $gid = $line[0];
		if (!exists($gene{$gid})) {
			my %a = ();
			$gene{$gid} = \%a;
		}
		$gene{$gid}->{'gid'} = $line[0];
		$gene{$gid}->{'ug'} = $line[1];
		$gene{$gid}->{'refseq'} = $line[2];
		$gene{$gid}->{'ensembl'} = $line[3];
		$gene{$gid}->{'name'} = $line[4];
		$gene{$gid}->{'alias'} = $line[5];
		$gene{$gid}->{'chr'} = $line[7];
		$gene{$gid}->{'desc'} = $line[8];
		$gene{$gid}->{'ttype'} = "NA";
		if (@line > 9) {
			$gene{$gid}->{'ttype'} = $line[9];
		}
	}
	close IN;
}

if (@geneDataFiles > 0) {
	for (my $j=0;$j<@geneDataFiles;$j++) {
		last if ($noGeneFlag);
		if ($organism ne 'unknown') {
			#print STDERR "`convertIDs.pl $geneDataFiles[$j] $organism $promoterIDtype yes yes yes > $tmpfile`;\n";
			`convertIDs.pl "$geneDataFiles[$j]" $organism $promoterIDtype yes yes yes > "$tmpfile"`;
		} else {
			`cp "$geneDataFiles[$j]" "$tmpfile"`;
		}
	
		my @currentHeader = ();	
		my $linecount = 0;
		my %geneData = ();

		open IN, $tmpfile;
		while (<IN>) {
			$linecount++;
			chomp;
			s/\r//g;
			my @line = split /\t/;
			if ($linecount == 1) {
				for (my $i=1;$i<@line;$i++){ 
					push(@currentHeader , $line[$i]);
				}
				next;
			}
			my $id = shift @line;
			my $oid = $line[0];
			while (scalar(@line) < scalar(@currentHeader)) {
				push(@line, "");
			}
			if (!exists($geneData{$oid})) {
				$geneData{$oid} = \@line;
			}
			if (!exists($geneData{$id})) {
				if ($id ne '') {
					$geneData{$id} = \@line;
				}
			}
		}
		close IN;

		push(@geneDataHeader, \@currentHeader);
		push(@geneData, \%geneData);
		`rm "$tmpfile"`;
			
	}
}

if ($gwasCatalog ne '') {
	print STDERR "\tChecking for GWAS risk SNP overlap\n";
	open IN, $gwasCatalog;
	open OUT, ">$tmpfile";
	my %gwas = ();
	my $idCount = 1;
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
	#print STDERR "`mergePeaks -d given $posFile $tmpfile > $tmpfile2 2> /dev/null`;\n";
	`mergePeaks -d given "$posFile" "$tmpfile" > "$tmpfile2" 2> /dev/null`;

	open IN, $tmpfile2;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if (@line < 10);
		my @peaks = split /\,/, $line[8];
		my @snps = split /\,/, $line[9];
		
		for (my $i=0;$i<@peaks;$i++) {
			my $pid = $peaks[$i];
			next if (!exists($peaks{$pid}));

			for (my $j=0;$j<@snps;$j++) {
				my $snpID = $snps[$j];
				next if (!exists($gwas{$snpID}));
			
				my $str = $gwas{$snpID}->{'snp'};
				if (!exists($peaks{$pid}->{'gwas'})) {
					$peaks{$pid}->{'gwas'}= $str;
				} else {
					$peaks{$pid}->{'gwas'} .= "," . $str;
				}
			}
		}
	}
	close IN;

	`rm -f "$tmpfile" "$tmpfile2"`;

}


#print out annotation file...

print STDERR "\tOutputing Annotation File...\n";
print "PeakID (cmd=$cmd)\tChr\tStart\tEnd\tStrand";
if ($tssFlag == 0) {
	print "\tPeak Score";
} else {
	#print "\tCAGE/EST Tag Count";
	print "\tNot Used";
}
print "\tFocus Ratio/Region Size";

if ($noGeneFlag == 0) {
	print "\tAnnotation\tDetailed Annotation"
			. "\tDistance to TSS\tNearest PromoterID\tEntrez ID\tNearest Unigene\tNearest Refseq"
			. "\tNearest Ensembl\tGene Name\tGene Alias\tGene Description\tGene Type";
}

if ($gwasCatalog ne '') {
	print "\tOverlapping GWAS Risk SNPs";
}

if ($cpgFlag==1) {
	print "\tCpG%\tGC%";
}
if ($consFlag) {
	print "\tConservation Average (0=low to 1=high conservation)";
}

if ($cmpGenomeFlag) {
	for (my $i=0;$i<@cmpGenomes;$i++) {
		print "\t% of peak aligned in $cmpGenomes[$i]";
		print "\t% alignment in $cmpGenomes[$i]";
	}
}

if (@allIndividuals > 0) {
	if ($editDistanceFlag) {
		foreach(@allIndividuals) {
			print "\t$_";
			print " Edit Distance from ref";
		}
	} else {
		print "\tSNPs (";
		my $c = 0;
		foreach(@allIndividuals) {
			print "," if ($c>0);
			$c++;
			print "$_";
		}
	}
}

foreach(@allDirs) { 
	if ($dirTypes{$_} eq 'tagDir') {
		if ($fpkmFlag) {
			print "\t$_ FPKM";
		} else {
			print "\t$_ Tag Count in $size bp ($tagTotals{$_} Total, normalization factor = $normFactors{$_}, effective total = $normValue)";
		}
	} elsif ($dirTypes{$_} eq 'bedGraph') {
		print "\t$_ bedGraph avg over $size bp";
	} elsif ($dirTypes{$_} eq 'wiggle') {
		print "\t$_ wig avg over $size bp";
	}
}
foreach(@peakFiles) {
	print "\t$_ Distance to nearest Peak, Peak ID";
}
foreach(@motifNames) {
	if ($mscoreFlag==1) {
		print "\t$_ Best Motif log-odds Score";
	} else {
		print "\t$_ Distance From Peak(sequence,strand,conservation)";
	}
	if ($cmpGenomeFlag) {
		for (my $i=0;$i<@cmpGenomes;$i++) {
			print "\tFound in $cmpGenomes[$i]";
		}
	}
	if (@revLiftover > 0) {
		for (my $i=0;$i<@cmpGenomes;$i++) {
			if ($mscoreFlag==1) {
				print "\t$_ Best Motif log-odds Score in $cmpGenomes[$i]";
			} else {
				print "\t$_ Distance From Peak(sequence,strand,conservation) in $cmpGenomes[$i]";
			}
			print "\tFound in $genome";
		}
	}
}
foreach(@geneDataHeader) {
	foreach(@$_) {
		print "\t$_";
	}
}
print "\n";

my $sortFlag = 'num';
foreach(keys %peaks) {
	#print STDERR "$_ $peaks{$_}->{'v'}\n";
	if ($peaks{$_}->{'v'} !~ /^[\-\d\.\+e]+$/ || $peaks{$_}->{'v'} eq '.') {
		$sortFlag = 'str';
	}
}
my @hits = ();
if ($sortFlag eq 'num') {
	@hits = sort {$peaks{$b}->{'v'} <=> $peaks{$a}->{'v'}} keys %peaks;
} else {
	@hits = sort {$peaks{$b}->{'v'} cmp $peaks{$a}->{'v'}} keys %peaks;
}


my $totalRows = 0;
my $numBadRows = 0;


foreach(@hits) {
	my $numBlanks =0;
	my $rowStr = '';

	my $hit = $_;
	my $chr = $peaks{$hit}->{'c'};
	my $start = $peaks{$hit}->{'s'};
	my $end = $peaks{$hit}->{'e'};
	my $dir = $peaks{$hit}->{'d'};
	my $v = $peaks{$hit}->{'v'};
	my $fdr = $peaks{$hit}->{'fdr'};
	my $ann = $peaks{$hit}->{'ann'};
	my $fullann = $peaks{$hit}->{'fullann'};
	my $cons = 'NA';
	my $peakLength = $end - $start;
	$peakLength = 1 if ($peakLength < 1);
	if ($peaks{$hit}->{'cons'} ne 'NA') {
		$cons = sprintf("%.2f",$peaks{$hit}->{'cons'}/9.0);
	}
	my $tssID = $peaks{$hit}->{'tss'};
	my $dist = $peaks{$hit}->{'tssDist'};


	$rowStr .= "$hit\t$chr\t$start\t$end\t$dir\t$v\t$fdr";
	if ($noGeneFlag == 0) {
		$rowStr .= "\t$ann\t$fullann\t$dist";
		my $ogID = $tssID;
		my $tssAltID = $tssID;
		$tssAltID =~ s/\.\d+?$//;
		if (!exists($acc2gene{$tssID}) && exists($acc2gene{$tssAltID})) {
			$tssID = $tssAltID;
		}
		if (exists($acc2gene{$tssID}) && $acc2gene{$tssID} ne 'NA' && exists($gene{$acc2gene{$tssID}})) {
			my $gid = $acc2gene{$tssID};
			my $ugid = $gene{$gid}->{'ug'};
			my $refseq = $gene{$gid}->{'refseq'};
			my $embl = $gene{$gid}->{'ensembl'};
			my $gname = $gene{$gid}->{'name'};
			my $alias = $gene{$gid}->{'alias'};
			my $desc = $gene{$gid}->{'desc'};
			my $ttype = $gene{$gid}->{'ttype'};
			$rowStr .= "\t$ogID\t$gid\t$ugid\t$refseq\t$embl\t$gname\t$alias\t$desc\t$ttype";
		} else {
			$rowStr .= "\t$ogID\t\t\t\t\t\t\t\t";
		}
	}
	if ($gwasCatalog ne '') {
		if (exists($peaks{$hit}->{'gwas'})) {
			$rowStr .= "\t$peaks{$hit}->{'gwas'}";
		} else {
			$rowStr .= "\t";
		}
	}
	if ($cpgFlag==1) {
		$rowStr .= "\t$peaks{$hit}->{'cpg'}\t$peaks{$hit}->{'gc'}";
		$numBlanks++ if ($peaks{$hit}->{'cpg'} eq 'NA' || $peaks{$hit}->{'cpg'} eq '');
		$numBlanks++ if ($peaks{$hit}->{'gc'} eq 'NA' || $peaks{$hit}->{'gc'} eq '');
	}
	if ($consFlag==1) {
		$rowStr .= "\t$cons";
		$numBlanks++ if ($cons eq 'NA' || $cons eq '');
	}
	if ($cmpGenomeFlag) {
		for (my $i=0;$i<@cmpGenomes;$i++) {
			my $v1 = "NA";
			my $v2 = $peaks{$hit}->{'gComp'}->[$i]->{'map'};
			if ($skipBlastn == 0) {
				$v1 = $peaks{$hit}->{'gComp'}->[$i]->{'pid'};
				$v2 = $peaks{$hit}->{'gComp'}->[$i]->{'paln'};
			}
			$rowStr .= "\t$v1\t$v2";
		}
	}
	if (@allIndividuals > 0) {
		if ($editDistanceFlag) {
			foreach(my $i=0;$i<@allIndividuals;$i++) {
				my $v = 0;
				if (exists($peaks{$hit}->{'snps'})) {
					if (scalar(@{$peaks{$hit}->{'snps'}}) > $i) {
						$v = $peaks{$hit}->{'snps'}->[$i];
					}
				}
				$rowStr .= "\t$v";
			}
		} else {
			$rowStr .= "\t";
			if (exists($peaks{$hit}->{'snps'})) {
				$rowStr .= $peaks{$hit}->{'snps'};
			}
		}
	}
	foreach(@allDirs){ 
		my $D = $_;
		if (exists($peaks{$hit}->{'t'}->{$D})) {
			my $v = $peaks{$hit}->{'t'}->{$D};
			if ($fpkmFlag) {
				$v *= 1000.0/$peakLength;
			}
			$rowStr .= "\t$v";
		} else {
			$rowStr .= "\t";
			$numBlanks++;
		}
	}
	foreach(@peakFiles) {
		my $str = '';
		if (exists($peaks{$hit}->{'p'}->{$_})) {
			if ($pCountFlag==1) {
				$str = $peaks{$hit}->{'p'}->{$_}->{'s'};
			} elsif ($pDistFlag==1) {
				$str = abs($peaks{$hit}->{'p'}->{$_}->{'d'});
			} elsif ($pDistFlag==2) {
				$str = $peaks{$hit}->{'p'}->{$_}->{'d'};
			} else {
				$str = $peaks{$hit}->{'p'}->{$_}->{'d'} . "," . $peaks{$hit}->{'p'}->{$_}->{'id'};
			}
		} else {
			$numBlanks++;
		}
		$rowStr .= "\t$str";
	}
	foreach(@motifNames) {
		my $str = '';
		my $motif = $_;
		my $nmotifs = 0;
		my $minDist2Motif = '';
		if (exists($peaks{$hit}->{'m'})) {
			if (exists($peaks{$hit}->{'m'}->{$motif})) {
				my @pos = sort {$a->{'p'} <=> $b->{'p'}} @{$peaks{$hit}->{'m'}->{$motif}};	
				my $c = 0;
				foreach(@pos) {
					if ($mscoreFlag) {
						my $score = $_->{'score'};
						$str .= $score;
						last;
					} 
					next if ($removeCloseMotifs && $_->{'valid'}==0);
					$nmotifs++;
					$str .= ',' if ($c > 0);
					$c++;
					my $con = sprintf("%.2f",$_->{'c'}/9.0);
					$str .= $_->{'p'} . "(" . $_->{'s'} . "," . $_->{'d'} . "," . $con . ")";
					if ($minDist2Motif eq '') {
						$minDist2Motif = $_->{'p'};
					} elsif (abs($_->{'p'}) < abs($minDist2Motif)) {
						$minDist2Motif = $_->{'p'};
					}
				}
			}
		}

		if ($nscoreFlag == 1) {
			$str = $nmotifs;
		} 
		if ($mdistFlag == 1) {
			$str = $minDist2Motif;
		}
		$rowStr .= "\t$str";
		if ($cmpGenomeFlag) {
			for (my $i=0;$i<@cmpGenomes;$i++) {
				my $str = '';
				if (exists($peaks{$hit}->{'m'})) {
					if (exists($peaks{$hit}->{'m'}->{$motif})) {
						my $total = 0;
						my $match = 0;

						foreach(@{$peaks{$hit}->{'m'}->{$motif}}) {
							$match += $_->{'gComp'}->[$i]->{'m'};
							$total++;
						}
						if ($total < 1) {
							$match = -1;
						} else {
							$match /= $total;
						}
						$str = $match;
					}
				}
				$rowStr .= "\t$str";
			}

			if (@revLiftover > 0) {
				for (my $i=0;$i<@cmpGenomes;$i++) {
					my $str = '';
					my $nmotifs = 0;
					my $minDist2Motif = '';
					if (exists($cpeaks[$i]->{$hit}->{'m'})) {
						if (exists($cpeaks[$i]->{$hit}->{'m'}->{$motif})) {
							my @pos = sort {$a->{'p'} <=> $b->{'p'}} @{$cpeaks[$i]->{$hit}->{'m'}->{$motif}};	
							my $c = 0;
							foreach(@pos) {
								if ($mscoreFlag) {
									my $score = $_->{'score'};
									$str .= $score;
									last;
								} 
								next if ($removeCloseMotifs && $_->{'valid'}==0);
								$nmotifs++;
								$str .= ',' if ($c > 0);
								$c++;
								#my $con = sprintf("%.2f",$_->{'c'}/9.0);
								my $con = 0;
								$str .= $_->{'p'} . "(" . $_->{'s'} . "," . $_->{'d'} . "," . $con . ")";
								if ($minDist2Motif eq '') {
									$minDist2Motif = $_->{'p'};
								} elsif (abs($_->{'p'}) < abs($minDist2Motif)) {
									$minDist2Motif = $_->{'p'};
								}
							}
						}
					}
		
					if ($nscoreFlag == 1) {
						$str = $nmotifs;
					} 
					if ($mdistFlag == 1) {
						$str = $minDist2Motif;
					}
					$rowStr .= "\t$str";
					$str = '';
					if (exists($cpeaks[$i]->{$hit}->{'m'})) {
						if (exists($cpeaks[$i]->{$hit}->{'m'}->{$motif})) {
							my $total = 0;
							my $match = 0;
	
							foreach(@{$cpeaks[$i]->{$hit}->{'m'}->{$motif}}) {
								$match += $_->{'gComp'}->[0]->{'m'};
								$total++;
							}
							if ($total < 1) {
								$match = -1;
							} else {
								$match /= $total;
							}
							$str = $match;
						}
					}
					$rowStr .= "\t$str";
				}
			}
		}
	}

	for (my $i=0;$i<@geneData;$i++) {
		if (exists($geneData[$i]->{$tssID})) {
			foreach(@{$geneData[$i]->{$tssID}}) {
				$rowStr .= "\t$_";
			}
		} else {
			$numBlanks++;
			foreach(@{$geneDataHeader[$i]}) {
				$rowStr .= "\t";
			}
		}
	}
			
	$rowStr .= "\n";


	$totalRows++;
	if ($noblanksFlag == 1) {
		if ($numBlanks > 0) {
			next;
			$numBadRows++;
		}
	}

	print $rowStr;
		
}
if ($noblanksFlag == 1) {
	print STDERR "\tRemoved $numBadRows with missing data (out of $totalRows)\n";
}

`rm "$posFile"`; 
`rm "$seqFile"` if ($seqFlag);
`rm "$consFile"` if ($consFlag);
print STDERR "\tDone annotating peaks file\n\n";
deleteFiles();

sub deleteFiles {
	foreach(keys %toDelete) {
		`rm "$_"`;
	}
}

sub calcMotifLogic {

	my ($peaks,$motifNames,$strand,$mlogicFile) = @_;

	my @codes = ();
	my @rcodes = ();
	my %index = ();
	for (my $i=0;$i<@$motifNames;$i++) {
		my $c = chr(65+$i);
		my $rc = chr(97+$i);
		print STDERR "\tCodes: $motifNames->[$i] $i $c $rc\n";
		push(@codes,$c);
		push(@rcodes,$rc);
		$index{$c} = $i;
		$index{$rc} = $i;
	}

	my %codeCounts = ();
	my $assStr = "";
	my %pCodes = ();
	my $peakTotal = 0;

	foreach(keys %$peaks) {
		my $peakID = $_;
		$peakTotal++;
		my @hits = ();
		my $bad = 0;
		for (my $i=0;$i<@$motifNames;$i++) {
			my $mname = $motifNames[$i];
			if (!exists($peaks->{$peakID}->{'m'}->{$mname})) {
				$bad = 1;
				next;
			}
			foreach(@{$peaks->{$peakID}->{'m'}->{$mname}}) {
				next if ($_->{'valid'}==0);
				my $p = $_->{'p'};
				my $d = $_->{'d'};
				my $c = $codes[$i];
				if ($d eq '-' || $d eq '1') {
					$c = $rcodes[$i];
				}
				my $info = {p=>$p,d=>$d,c=>$c};
				push(@hits, $info);
			}
		}
		#next if ($bad==1);
		@hits = sort {$a->{'p'} <=> $b->{'p'} } @hits;
		my $codeStr = "";
		foreach(@hits) {
			$codeStr .= $_->{'c'};
		}
		my $usedStr = $codeStr;
		my $rv = revOppCode($codeStr);
		if ($strand eq 'both') {
			my $cc = $codeStr cmp $rv;
			if ($cc > 0) {
				$usedStr = $rv;
			}
		}

		if (exists($codeCounts{$usedStr})) {
			$codeCounts{$usedStr}++;
		} else {
			$codeCounts{$usedStr}=1;
		}
		$pCodes{$peakID}= $usedStr;
		$assStr .= "$peakID\t$peaks->{$peakID}->{'c'}\t$peaks->{$peakID}->{'s'}\t$peaks->{$peakID}->{'e'}\t$peaks->{$peakID}->{'d'}\t$usedStr\n";
	}
	my %codeSubs = ();
	foreach(keys %codeCounts) {
		my $code = $_;
		my $rcode = revOppCode($code);	
		my $total = 0;
		foreach(values %pCodes) {
			my $c = $_;
			my $f = 0;
			if (index($c,$code) > -1) {
				$total++;
				next;
			} elsif (index($c,$rcode) > -1) {
				$total++;
				next;
			}
		}
		$codeSubs{$code}=$total;
	}


	open OUT, ">$mlogicFile";
	
	print OUT "Motif Name\t+ strand code\t- strand code\n";
	for (my $i=0;$i<@codes;$i++) {
		print OUT "$motifNames->[$i]\t$codes[$i]\t$rcodes[$i]\n";
	}
	print OUT "\n\n";

	print OUT "CRM codes\t# of Peaks\t% of Peaks\t# of Peaks (exact)\t% of Peaks (exact)\n";
	my @crm = sort {$codeSubs{$b} <=> $codeSubs{$a}} keys %codeCounts;
	foreach(@crm) {
		my $v= $_;
		$v = 'None' if ($v eq '');
		print OUT "$v";
		$v = $codeSubs{$_};
		print OUT "\t$v";
		$v = sprintf("%.3lf",$codeSubs{$_}/$peakTotal);
		print OUT "\t$v";
		$v = $codeCounts{$_};
		print OUT "\t$v";
		$v = sprintf("%.3lf",$codeCounts{$_}/$peakTotal);
		print OUT "\t$v";
		print OUT "\n";
	}


	print OUT "\n\n";
	print OUT "PeakID\tchr\tstart\tend\tstrand\tCRM code\n";
	print OUT $assStr;
	close OUT;
	

}
sub revOppCode {
	my ($code) = @_;
	$code = reverse($code);
	my $len = length($code);
	my $ncode = "";
	for (my $i=0;$i<$len;$i++) {
		my $c = substr($code,$i,1);
		my $n = ord($c);
		if ($n < 95) {
			$n += 32;
		} else {
			$n -= 32;
		}
		$ncode .= chr($n);
	}
	return $ncode;
}

