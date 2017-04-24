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

my $maxDistance = 1000;

sub printCMD {

	print STDERR "\n\tUsage: assignTSStoGene.pl <tss peak file> [options]\n";
	print STDERR "\n\tThis program takes TSS defined by 5' RNA sequencing and assigns them to genes\n";
	print STDERR "\n\tOne of the following options are required to define genes:\n";
	print STDERR "\t\t-genome <genomeVersion> (use default homer gene annotation/RefSeq)\n";
	print STDERR "\t\t-gtf <GTF file> (use custom gene annotation, can also use -gff or -gff3)\n";
	print STDERR "\t\t\t-gid (use gene_id with GTF file)\n";
	print STDERR "\t\t-bed <peak/BED file> (use custom gene annotation in bed/peak file format)\n";
	print STDERR "\t\t-refTSS <tss peak file> (only supply reference TSS positions)\n";
	print STDERR "\n\tOther Options:\n";
	print STDERR "\t\t-d <#> (max dist from tss to gene allowed, default: $maxDistance)\n";
	print STDERR "\t\t-nokeepRef (don't keep reference promoters not found in the tss peak file, default: keep)\n";
	print STDERR "\t\t-keepOrphans (keep TSS without reference annotation, default: remove)\n";
	print STDERR "\t\t-bedOut <file name> (output genes with new 5'end)\n";
	print STDERR "\n\t\t-3p (do 3' end assignment instead of TSS, assumes peaks are on - strand rel to gene)\n";
	print STDERR "\n";
	exit;

}

if (@ARGV < 3) { 
	printCMD();
}

my $cmd = $ARGV[0];
for (my $i=1;$i<@ARGV;$i++) {
	$cmd .= " " . $ARGV[$i];
}

print STDERR "\n";
my %toDelete = ();

my $gtfFile = "";
my $bedFile = "";
my $bedOut = "";
my $gidFlag = "";
my $gtfFormat = "";
my $oldTSSfile = "";
my $genome = "";
my $polyAflag = 0;

my $organism = "unknown";
my $promoter = "default";
my $consDir = "";
my $genomeDir = "";
my $genomeParseDir = "";
my $customGenome = 0;
my $keepRef = 1;
my $keepOrphans = 0;

my $tssFile = $ARGV[0];

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-genome') {
		$genome  = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-nokeepRef') {
		$keepRef = 0;
	} elsif ($ARGV[$i] eq '-bed') {
		$bedFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bedOut') {
		$bedOut = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-keepOrphans') {
		$keepOrphans = 1;
	} elsif ($ARGV[$i] eq '-gtf') {
		$gtfFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-gid') {
		$gidFlag = " -gid";
	} elsif ($ARGV[$i] eq '-gff') {
		$gtfFile = $ARGV[++$i];
		$gtfFormat = "-gff";
	} elsif ($ARGV[$i] eq '-gff3') {
		$gtfFile = $ARGV[++$i];
		$gtfFormat = "-gff3";
	} elsif ($ARGV[$i] eq '-refTSS') {
		$oldTSSfile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-d') {
		$maxDistance = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-3p') {
		$polyAflag = 1;
		print STDERR "\tTreating tss file as 3' end file instead of 5' end file (- strand from gene)\n";
	} else {
		print STDERR "Not recognized: $ARGV[$i]\n";
		printCMD();
	}
}


if ($genome eq "" && $gtfFile eq "" && $oldTSSfile eq '' && $bedFile eq '') {
	print STDERR "!!! Missing reference genes (-genome/-gtf/-refTSS/-bed) !!!\n";
	printCMD();
}


if ($genome ne '') {
	if ($genome eq 'none') {
		print STDERR "!!! Can't specify -genome none for this command !!!\n";
		printCMD();
	} elsif (!exists($config->{'GENOMES'}->{$genome})) {
		print STDERR "!!! Genome not recognized - that's a problem for this command !!!\n";
		printCMD();
	} else {
		$genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'};	
		$organism = $config->{'GENOMES'}->{$genome}->{'org'};	
		$promoter = $config->{'GENOMES'}->{$genome}->{'promoters'};
		$consDir = $config->{'GENOMES'}->{$genome}->{'directory'} . "/conservation/";
		if ($oldTSSfile eq '') {
			$oldTSSfile = $genomeDir . "/" . $genome . ".tss";
			$oldTSSfile = $genomeDir . "/" . $genome . ".tts" if ($polyAflag);
	
		}
	}
	print STDERR "\tGenome = $genome\n";
	print STDERR "\tOrganism = $organism\n";
}

my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";

if ($gtfFile ne '') {
	if ($oldTSSfile eq '') {
		if ($polyAflag) {
			`parseGTF.pl "$gtfFile" tts $gtfFormat $gidFlag > "$tmpFile2"`;
			`adjustPeakFile.pl "$tmpFile2" -flip > "$tmpFile"`;
		} else {
			`parseGTF.pl "$gtfFile" tss $gtfFormat $gidFlag > "$tmpFile"`;
		}
	} else {
		`cp "$oldTSSfile" "$tmpFile"`;
	}
} elsif ($bedFile ne '') {
	if ($oldTSSfile eq '') {
		`adjustPeakFile.pl "$bedFile" -5p > "$tmpFile"`;
		`adjustPeakFile.pl "$bedFile" -3p > "$tmpFile"` if ($polyAflag);
	} else {
		`cp "$oldTSSfile" "$tmpFile"`;
	}
} elsif ($oldTSSfile eq '') {
	my $refSeqTSS = "$genomeDir/$genome.tss";
	`cp "$refSeqTSS" "$tmpFile"`;
} else {
	`cp "$oldTSSfile" "$tmpFile"`;
}		
$toDelete{$tmpFile} = 1;
$toDelete{$tmpFile2} = 1;

`bed2pos.pl -check "$tssFile" > "$tmpFile3"`;

my %tss = ();
open IN, "$tmpFile3";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $id = $line[0];
	my $chr=  $line[1];
	my $start = $line[2];
	my $end = $line[3];
	my $strand = $line[4];
	my $v = $line[5];
	my $focus = 1;
	$focus = $line[6] if (@line > 6);
	if ($polyAflag) {
		if ($strand eq '0' || $strand eq '+') {
			$strand = 1;
		} else {
			$strand = 0;
		}
	}
	$tss{$id} = {c=>$chr,s=>$start,e=>$end,d=>$strand,v=>$v,f=>$focus};
}
close IN;

`annotateRelativePosition.pl "$tmpFile3", "$tmpFile", 1 > "$tmpFile2"`;
open IN, "$tmpFile2";
my %genes = ();
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $tssID =$line[0];
	my $gid =$line[1];
	my $dist = $line[2];
	if (abs($dist) > $maxDistance) {
		next;
	}
	my $d ={id=>$tssID,dist=>$dist};
	if (!exists($genes{$gid})) {
		$genes{$gid} = $d;
	} else {
		if ($tss{$tssID}->{'v'} > $tss{$genes{$gid}->{'id'}}->{'v'}) {
			$genes{$gid} = $d;
		} elsif ($tss{$tssID}->{'v'} == $tss{$genes{$gid}->{'id'}}->{'v'}
				&& abs($dist) < abs($genes{$gid}->{'dist'})) {
			$genes{$gid} = $d;
		}
	}
}
close IN;

my %printedTSS = ();
open IN, $tmpFile;
while (<IN>){
	chomp;
	s/\r//g;
	next if (/^#/);
	my $og = $_;
	my @line = split /\t/;

	if (exists($genes{$line[0]})) {
		my $tssID = $genes{$line[0]}->{'id'};
		$printedTSS{$tssID} = 1;
		print "$line[0]\t$tss{$tssID}->{'c'}\t$tss{$tssID}->{'s'}\t$tss{$tssID}->{'e'}\t$tss{$tssID}->{'d'}\t$tss{$tssID}->{'v'}\t$tss{$tssID}->{'f'}\n";
	} elsif ($keepRef) {
		print "$og\n";
	}
}
close IN;
if ($keepOrphans) {
	foreach(keys %tss) {
		next if (exists($printedTSS{$_}));
		my $tssID = $_;
		print "$tssID\t$tss{$tssID}->{'c'}\t$tss{$tssID}->{'s'}\t$tss{$tssID}->{'e'}\t$tss{$tssID}->{'d'}$tss{$tssID}->{'v'}\t$tss{$tssID}->{'f'}\n";
	}
}

if ($bedOut ne '' && $bedFile ne '') {
	`bed2pos.pl "$bedFile" -check > "$tmpFile"`;
	open IN, $tmpFile;
	open OUT, ">$bedOut" or die "Could not open $bedOut for writing!\n";
	my $c = 0;
	my $cc =0;
	while (<IN>) {
		my $og = $_;
		chomp;
		s/\r//g;
		$c++;
		my @line = split /\t/;
		if (exists($genes{$line[0]})) {
			my $tssID = $genes{$line[0]}->{'id'};
			my $s = $line[2];
			my $e = $line[3];
			my $p = ($tss{$tssID}->{'s'}+$tss{$tssID}->{'e'})/2;
			if ($polyAflag) {
				if ($tss{$tssID}->{'d'} eq '+' || $tss{$tssID}->{'d'} eq '0') {
					$e = $p;
					$s = $e-1 if ($e-$s < 1);
				} else {
					$s = $p;
					$e = $s+1 if ($e-$s < 1);
				}
			} else {
				if ($tss{$tssID}->{'d'} eq '+' || $tss{$tssID}->{'d'} eq '0') {
					$s = $p;
					$e = $s+1 if ($e-$s < 1);
				} else {
					$e = $p;
					$s = $e-1 if ($e-$s < 1);
				}
			}
			print OUT "$line[0]\t$tss{$tssID}->{'c'}\t$s\t$e\t$line[4]";
			for (my $i=5;$i<@line;$i++) {
				print OUT "\t$line[$i]";
			}
			print OUT "\n";
			$cc++;
		} elsif ($keepRef) {
			print OUT $og;
		}
	}
	close IN;
	my $r = $cc/$c;
	print STDERR "\t$cc of $c genes had TSS assigned ($r)\n";
}

deleteFiles();

sub deleteFiles {
	foreach(keys %toDelete) {
		`rm "$_"`;
	}
}
