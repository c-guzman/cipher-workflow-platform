#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";


# Copyright 2009-2014 Christopher Benner <cbenner@ucsd.edu>
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

my $pthresh = 0.00001;

my $config = HomerConfig::loadConfigFile();

if (@ARGV < 2) {
	cmdLineOptions();
}
sub cmdLineOptions {
	print STDERR "\n\tThis program will assess your input peak or tag directory for enrichment in a variety of\n";
	print STDERR "\tgenome annotations and previous experimental results. (Written with Dr. Kasey Hutt)\n";
	print STDERR "\n\tUsage: GenomeOntology.pl <peak file/Tag Directory> <genome> <Output Directory> [additional options]\n";
	print STDERR "\n\tThis produces an HTML page and additional data files in the <Output Directory>.  \n";
	print STDERR "\tTo customize this analysis, place additional peak/annotation files in the directory:\n";
	print STDERR "\t\t$homeDir/data/genomes/\"GENOME\"/annotations/custom/\n";

	print STDERR "\n\tPossible Genomes:\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}\n";
	}

	print STDERR "\n\tBasic options:\n";
	print STDERR "\t\t-gsize <#> (Genome size used for significance calculations)\n";
	print STDERR "\t\t-bg <peakfile/Tag Directory> (Performs additional significance calculations\n";
	print STDERR "\t\t\trelative to Control Peaks/Tags)\n";
	print STDERR "\n";
	exit;
}

my $peakFile = $ARGV[0];
my $genome = $ARGV[1];
my $outputDir = $ARGV[2];
my $controlPeak = "";
if (!exists($config->{'GENOMES'}->{$genome})) {
	print STDERR "!!!!Genome $genome not found in $homeDir/config.txt\n\n";
    print STDERR "\tTo check if is available, run \"perl $homeDir/configureHomer.pl -list\"\n";
    print STDERR "\tIf so, add it by typing \"perl $homeDir/configureHomer.pl -install $genome\"\n";
    print STDERR "\n";
	exit;
}
my $genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'} . "/";
my $dirFlag = 0;

if (-f $peakFile ) {
	print STDERR "\tInput is a Peak file\n";
} elsif (-d $peakFile ) {
	print STDERR "\tInput is a Tag Directory\n";
	$dirFlag = 1;
} else {
	print STDERR "!!! Could not recognize input: $peakFile !!!\n";
	exit;
}

my $gsize = 2e9;


for (my $i=3;$i<@ARGV; $i++) {
	if ($ARGV[$i] eq '-gsize') {
		$gsize = $ARGV[++$i];
		print STDERR "\tEffective Genome size set to $gsize\n";
	} elsif ($ARGV[$i] eq '-bg') {
		$controlPeak = $ARGV[++$i];
		print STDERR "\tUsing $controlPeak as a control\n";
	} elsif ($ARGV[$i] eq '-prefix') {
	} else {
		print STDERR "$ARGV[$i] not recognized!\n\n";
		cmdLineOptions();
		exit;
	}
}


my $annDir = $genomeDir . "annotations/";

my $rand = rand();
my $tmpfile = $rand . ".tmp";
my $tmpfile2 = $rand . ".2.tmp";

my %dirs = ();
`ls $annDir > $tmpfile`;
open IN, $tmpfile;
while (<IN>) {
	chomp;
	#print STDERR "$_\n";
	$dirs{$_} = $annDir . $_;
}
close IN;
`rm $tmpfile`;

my @fileCounts = ();

my @dirs = keys %dirs;
my $numDirs = scalar(@dirs);

if ($numDirs < 1) {
	print STDERR "!!! Something could be wrong - didn't see any directories in the directory:\n";
	print STDERR "\t$annDir\n";
	exit;
}

`mkdir -p $outputDir`;
open HTML, ">$outputDir/GenomeOntology.html";
print HTML "<HTML><HEAD><TITLE>Genome Ontology Results for $peakFile</TITLE></HEAD><BODY>\n";
print HTML "<H1>Genome Ontology Enrichment Results</H1>\n";
print HTML "<p>Input file: $peakFile";
if ($dirFlag == 0) {
	print HTML " (peak file)\n";
} else {
	print HTML " (Tag Directory)\n";
}
print HTML "</p>\n";
print HTML "<p>Links to enrichment analysis for each set of annotations:</p>\n";
print HTML "<UL>\n";


my @header = ();

my %all = ();
for (my $i=0;$i<@dirs;$i++) {
	my $dname = $dirs[$i];
	my $dir = $dirs{$dname};
	$fileCounts[$i]=0;
	`ls $dir > $tmpfile`;
	open IN, $tmpfile;
	while (<IN>) {
		chomp;
		#print "$_\n";
		$fileCounts[$i]++;
	}
	close IN;
	`rm $tmpfile`;

	my $z = $i+1;
	if ($fileCounts[$i] < 1) {
		print STDERR "\tAnalyzing $z of $numDirs annotation sets ($dname - $fileCounts[$i] total) Skipping...\n";
		next;
	}
	print STDERR "\tAnalyzing $z of $numDirs annotation sets ($dname - $fileCounts[$i] total)\n";
	my $bgOption = "";
	if ($controlPeak ne '') {
		$bgOption = " -bg \"$controlPeak\" ";
	}
	#print STDERR "`genomeOntology -d given -gsize $gsize $bgOption $peakFile $dir/* > $tmpfile`;\n";
	`ls "$dir"/* > "$tmpfile2"`;
	#`genomeOntology -d given -gsize $gsize $bgOption "$peakFile" "$dir"/* > $tmpfile`;
	`genomeOntology -d given -gsize $gsize $bgOption "$peakFile" -file "$tmpfile2" > "$tmpfile"`;
	`rm -f "$tmpfile2"`;

	open OUT, ">$outputDir/$dname.genomeOntology.txt";
	open IN, $tmpfile;
	my $cc = 0;
	my %results = ();
	while (<IN>) {
		$cc++;
		if ($cc == 1) {
			print OUT "Name\t$_";
			my @line = split /\t/;
			if ($dirFlag == 0) {
				$line[1]=~ s/\[ref=(\d*)\]//;
				$totalInput = $1;
				$line[2]=~ s/\[ref=(\d*)\]//;
				$totalInputBp = $1;
				if ($controlPeak ne '') {
					$line[10]=~/,\[total=(\d*)\]/;
					$totalControl = $1;
					$line[11]=~/,\[total=(\d*)\]/;
					$totalControlBp = $1;
				}
			} else {
				$line[3]=~ /\[Total=(.*)\]/;
				$totalInput = $1;
				if ($controlPeak ne '') {
					$line[8]=~/,\[Total=(.*)\]/;
					$totalControl = $1;
				}
			}
			if (@header < 1) {
				@header = @line;
			}
			next;
		}
		chomp;
		my $og = $_;
		my @line = split /\t/;
		my $id = $line[0];
		$id =~ s/^.*\///;
		$id =~ s/\.ann\.txt$//;
		my $score = 0;
		if ($dirFlag==1) {
			$score = $line[6];
			if ($controlPeak ne '') {
				$score = $line[14];
			}
		} else {
			$score = $line[8];
			if ($controlPeak ne '') {
				$score = $line[17];
			}
		}
		$results{$id} = {og=>$og,s=>$score};
		$allID = $dname . "-" . $id;
		$all{$allID} = {id=>$id, s=>$score, d=>\@line,n=>$dname};
	}
	close IN;
	`rm $tmpfile`;

	my @results = sort {$results{$a}->{'s'} <=> $results{$b}->{'s'}} keys %results;	
	foreach(@results) {
		print OUT "$_\t$results{$_}->{'og'}\n";
	}
	close OUT;

	print HTML "\t<LI>$dname - <A HREF=\"$dname.genomeOntology.txt\">$dname.genomeOntology.txt</A> ($fileCounts[$i] total annotations)</LI>\n";
}
print HTML "</UL>\n";
print HTML "<BR/><H4>Ranked list of results: (Positive log P-values signify divergence - Stuff at the bottom may be significantly NOT associated)</H4>\n";

if ($dirFlag == 0) {
	print HTML "<P>Total Input Regions ($peakFile): $totalInput, $totalInputBp bp</P>\n";
	if ($controlPeak ne '') {
		print HTML "<P>Total Control Regions ($controlPeak): $totalControl, $totalControlBp bp</P>\n";
	}
} else {
	print HTML "<P>Total Tag Count ($peakFile): $totalInput</P>\n";
	if ($controlPeak ne '') {
		print HTML "<P>Total Control Tag Count ($controlPeak): $totalControl</P>\n";
	}
}
print HTML "<TABLE border=\"1\" cellpading=\"0\" cellspacing=\"0\"><TR>\n";
print HTML "   <TD>P-value</TD><TD>Log P-value</TD>\n";
if ($controlPeak ne '') {
	print HTML "   <TD>P-value vs. Control</TD><TD>Log P-value vs. Control</TD>\n";
}
print HTML "<TD>Annotation</TD><TD>Ann Group</TD>";
foreach(my $i=1;$i<@header;$i++) {
	print HTML "<TD>$header[$i]</TD>";
}
print HTML "</TR>\n";

my @results = sort {$all{$a}->{'s'} <=> $all{$b}->{'s'}} keys %all;
for (my $i=0;$i<@results;$i++) {
	my $allid = $results[$i];
	my $id = $all{$allid}->{'id'};
	my $gid = $all{$allid}->{'n'};
	my $logP = 0;
	my $p = 1.0;
	if ($dirFlag == 0) {
		$logP = $all{$allid}->{'d'}->[8];
		$p =  $all{$allid}->{'d'}->[9];
	} else {
		$logP = $all{$allid}->{'d'}->[6];
		$p =  $all{$allid}->{'d'}->[7];
	}
	next if ($p > $pthresh);

	if ($p < 1e-50) {
		$p = -1*abs($logP);
		$p /= 2.303;
		$p = "1e" . floor($p);
	}

	if (length($id) > 100) {
		$id = substr($id,0,100) . "...";
	}

	print HTML "\t<TR><TD>$p</TD><TD>$logP</TD>\n";
	if ($controlPeak ne '') {
		my $logP = 0.0;
		my $p = 1.0;
		if ($dirFlag == 0) {
			$logP = $all{$allid}->{'d'}->[17];
			$p = $all{$allid}->{'d'}->[18];
		} else {
			$logP = $all{$allid}->{'d'}->[14];
			$p = $all{$allid}->{'d'}->[15];
		}
		if ($p < 1e-50) {
			$p = -1*abs($logP);
			$p /= 2.303;
			$p = "1e" . floor($p);
		}
		print HTML "\t<TD>$p</TD><TD>$logP</TD>\n";
	}


	print HTML "<TD>$id</TD><TD>$gid</TD>\n";
	for (my $i=1;$i<@{$all{$allid}->{'d'}};$i++) {
		print HTML "<TD>$all{$allid}->{'d'}->[$i]</TD>";
	}
	print HTML "</TR>\n";
	print HTML "\t</TR>";

}
print HTML "</TABLE>\n";
print HTML "</BODY></HTML>\n";

exit;

