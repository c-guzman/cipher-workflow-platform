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

	print STDERR "\n\tUsage: getConservedRegions.pl <peak file> <genome version>  [additional options...]\n";
	print STDERR "\n\tAnalyzes regions near peaks and return conserved islands as a peak file\n";
	print STDERR "\n\tAvailable Genomes (required argument): (name,org,directory,default promoter set)\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\t$config->{'GENOMES'}->{$_}->{'directory'}"
					. "\t$config->{'GENOMES'}->{$_}->{'promoters'}\n";
	}
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-size <#|given> (size of region centered on peaks to look for conserved islands)\n";
	print STDERR "\t\t\tdefault: given\n";
	print STDERR "\t\t-bufSize <#> (size of area around conserved islands to include)\n";
	print STDERR "\t\t\tdefault: 25\n";
	print STDERR "\t\t-outSize <#> (size of conserved island segments to output [larger islands will be split])\n";
	print STDERR "\t\t\tdefault: 100\n";
	print STDERR "\t\t-cons <0.0-1.0> (phastCons score needed to define conservation island, 0=all used)\n";
	print STDERR "\t\t\tdefault: 0.8\n";
	print STDERR "\t\t-p <peakfile> [peakfile2]... (peak regions to exclude)\n";
	print STDERR "\t\t-keepExons (By default, exons are excluded)\n";
	print STDERR "\n";
	exit;
}

if (@ARGV < 2) { 
	printCMD();
}


my $peakFile = $ARGV[0];
my $genome = $ARGV[1];

my $size = 'given';
my $sizeMove = 0;
my $updateSize = 0;
my @excludeFiles = ();
my $keepExons = 0;

my $outSize = 100;
my $bufSize = 25;

my $consThresh = 0.8;

for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-size') {
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
	} elsif ($ARGV[$i] eq '-p') {
		print STDERR "\tPeak Files:\n";
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@excludeFiles, $ARGV[$i]);
			print STDERR "\t\t$ARGV[$i]\n";
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} elsif ($ARGV[$i] eq '-keepExons') {
		$keepExons = 1;
	} elsif ($ARGV[$i] eq '-outSize') {
		$outSize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-bufSize') {
		$bufSize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cons') {
		$consThresh = $ARGV[++$i];
	} else {
		print STDERR "!!! Don't recognize: $ARGV[$i]\n";
		printCMD();
	}
}

# since conservation is actually reported from 0-9
$consThresh = 0.9 if ($consThresh > 0.9);
my $minPeakSize = $bufSize+10;

my $genomeDir = "";
my $organism = "unknown";
if (!exists($config->{'GENOMES'}->{$genome})) {
	print STDERR "!!! This program will only work with preconfigured genomes...\n";
	exit;
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

print STDERR "\n\tPeak File = $peakFile\n";
print STDERR "\tGenome = $genome\n";
print STDERR "\tOrganism = $organism\n";
my $conservationDir = $genomeDir . "/conservation/";

if ($keepExons == 0) {
	push(@excludeFiles, $genomeDir . "/annotations/basic/exons.ann.txt");
}

#tmp files
my $rand = rand();
my $tmpfile = $rand . ".tmp";
my $tmpfile2 = $rand . ".2.tmp";
my $tmpfile3 = $rand . ".3.tmp";
my %toDelete = ();

`bed2pos.pl "$peakFile" -check > "$tmpfile"`;
if ($size eq 'given') {
	`cp "$tmpfile" "$tmpfile2"`;
} else {
	`resizePosFile.pl "$tmpfile" $size $sizeMove > "$tmpfile2"`;
}

`homerTools extract "$tmpfile2" "$conservationDir" > "$tmpfile3"`;



my %peaks = ();
open IN, $tmpfile2;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line = split /\t/;
	my $id =$line[0];
	my $chr = $line[1];
	my $start = $line[2];
	my $end = $line[3];
	my $dir = $line[4];
	if ($dir eq '-' || $dir eq '1') {
		$dir = '-';
	} else {
		$dir = '+';
	}
	my $p = {c=>$chr,s=>$start,e=>$end,d=>$dir};
	$peaks{$id} = $p;
}
close IN;

my %islands = ();
my %chrs = ();

open IN, $tmpfile3;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $id = $line[0];
	my @cons = split //,$line[1];
	next if (!exists($peaks{$id}));
	my $chr = $peaks{$id}->{'c'};
	my $start = $peaks{$id}->{'s'};
	my $end = $peaks{$id}->{'e'};
	my $dir = $peaks{$id}->{'d'};
	
	my $curStart = -1e10;
	my $curEnd = -1e10;

	my @sets = ();
	for (my $i=0;$i<@cons;$i++) {
		$cons[$i] = 0 if ($cons[$i] eq 'N');
		next if ($cons[$i]/9.99 < $consThresh);
		if ($curEnd < $i-1) {
			if ($curStart > -1e8) {
				my @a = ($curStart, $curEnd+$bufSize);
				push(@sets, \@a);
			}
			$curStart = $i-$bufSize;
		}
		$curEnd = $i;
	}
	if ($curStart > -1e8) {
		my @a = ($curStart, $curEnd+$bufSize);
		push(@sets, \@a);
	}

	my @sets2 = ();
	my $index = 0;
	foreach(@sets) {
		if ($index == 0) {
			push(@sets2, $_);
			$index++;
		} else {
			if ($_->[0] - $sets2[$index-1]->[1] < 0) {
				$sets2[$index-1]->[1] = $_->[1];
			} else {
				push(@sets2, $_);
				$index++;
			}
		}
	}

	my $subID = 1;
	foreach(@sets2) {
		my $s = 0;
		my $e = 0;
		if ($dir eq '+') {
			$s = $start + $_->[0];
			$e = $start + $_->[1];
		} else {
			$s = $end - $_->[1];
			$e = $end - $_->[0];
		}
		my $p = {c=>$chr,s=>$s,e=>$e,d=>$dir};
		my $sid = $id . "-" . $subID++;
		#print "$sid\t$chr\t$s\t$e\t$dir\n";
		$islands{$sid} = $p;
		if (!exists($chrs{$chr})) {
			my %a = ();
			$chrs{$chr}= \%a;
		}
		$chrs{$chr}->{$sid} = $p;
	}
}
close IN;

my %exons = ();
my $bid = 1;
for (my $i=0;$i<@excludeFiles;$i++) { 
	print STDERR "\tProcessing $excludeFiles[$i]\n";
	`bed2pos.pl "$excludeFiles[$i]" -check > "$tmpfile"`;
	open IN, "$tmpfile";
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $bid++;
		my $chr = $line[1];
		my $p = {s=>$line[2],e=>$line[3]};
		if (!exists($exons{$chr})) {
			my %a = ();
			$exons{$chr}=\%a;
		}
		$exons{$chr}->{$id} = $p;
	}
	close IN;
}

my %report = ();
my @chrs = sort {$a cmp $b} keys %chrs;
foreach(@chrs) {
	my $chr = $_;
	print STDERR "\tFiltering conserved islands on $chr\n";
	
	my @pids = sort {$chrs{$chr}->{$a}->{'s'} <=> $chrs{$chr}->{$b}->{'s'}} keys %{$chrs{$chr}};
	my @eids = ();
	if (exists($exons{$chr})) {
		@eids = sort {$exons{$chr}->{$a}->{'s'} <=> $exons{$chr}->{$b}->{'s'}} keys %{$exons{$chr}};
	}

	my $index = 0;

	#remove redundant
	my $lastStart = -1;
	my $lastEnd = -1;
	foreach(@pids) {
		my $id = $_;
		while ($index < scalar(@eids) && $exons{$chr}->{$eids[$index]}->{'e'} < $chrs{$chr}->{$id}->{'s'}) {
			$index++;
		}

		my $start = $chrs{$chr}->{$id}->{'s'};
		my $end = $chrs{$chr}->{$id}->{'e'};
		next if ($start == $lastStart && $end == $lastEnd);
		$lastStart = $start;
		$lastEnd = $end;
		my $dir = $chrs{$chr}->{$id}->{'d'};
		my $len = $end-$start;
		my @curpeaks = ();
		push(@curpeaks,[$start,$end]);
		my @toSkip = ();
		for (my $i=$index;$i<@eids;$i++) {
			my $fs = $exons{$chr}->{$eids[$i]}->{'s'};
			my $fe = $exons{$chr}->{$eids[$i]}->{'e'};

			last if ($fs > $end);
			next if ($fe < $start);

			foreach(@curpeaks) {
#print STDERR "$id vs. $i: $_->[0],$_->[1] vs. $fs,$fe\n";
				next if ($fe < $_->[0]);	
				next if ($fs > $_->[1]);	
				if ($fs > $_->[0]) {
					if ($fe < $_->[1]) {
						my $e = $_->[1];
#print STDERR "\tnew!: $fe,$e\n";
						push(@curpeaks,[$fe,$e]);
					}
					$_->[1] = $fs;
				} elsif ($fe < $_->[1]) {
					$_->[0] = $fe;
				} else {
					$_->[0] = $fs;
					$_->[1] = $fs;
					#remove
				}
#print STDERR "\tnow: $_->[0],$_->[1]\n";
			}
		}
		my $ssid = 1;
		foreach(@curpeaks) {
			next if ($_->[1]-$_->[0] < $minPeakSize);
			my $p = {c=>$chr,s=>$_->[0],e=>$_->[1],d=>$dir};
			my $sid = $id . "-" . $ssid++;
			$report{$sid} = $p;
		}
	}

}

my $total = 0;
my @ids = sort {$report{$a}->{'c'} cmp $report{$b}->{'c'} || $report{$a}->{'s'} <=> $report{$b}->{'s'}} keys %report;
foreach(@ids) {
	my $id = $_;
	my $chr = $report{$_}->{'c'};
	my $s = $report{$_}->{'s'};
	my $e = $report{$_}->{'e'};
	my $d = $report{$_}->{'d'};
	my @pos = ();
	my $ssid = 1;
	if ($e-$s < $outSize*1.5) {
		push(@pos,[$s,$e]);
	} else {
		my $n = floor(($e-$s)/$outSize+0.5);
		my $interval = ($e-$s)/$n;
		for (my $i=0;$i<$n;$i++) {
			push(@pos, [floor($s+$i*$interval),floor($s+($i+1)*$interval)]);
		}
	}
	foreach(@pos) {
		my $iid = $id . "-" . $ssid++;
		print "$iid\t$chr\t$_->[0]\t$_->[1]\t$d\n";
		$total++;
	}
}
print STDERR "\n\t$total Conserved Islands found\n\n";

`rm -f "$tmpfile" "$tmpfile2" "$tmpfile3"`;
