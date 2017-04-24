#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";


# Copyright 2009-2016 Christopher Benner <cbenner@salk.edu>
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


if (@ARGV < 2) {
	printCMD();
}

my $denovoFile = $ARGV[0];
my $genome = $ARGV[1];
my $refRNAfile = "";
my $repeatFile = "";

my $tssDistance = 500;
my $repeatCheckRegion = 100;
my $noRepeats = 0;
my $minOverlap = 0.1;
my $tagDir = "";
my $maxFold = 0;
$strand = "+";
$minRPKM = 0.01;
$minP = 10;
my $tssFile = "";
my $tssRange = 1000;


sub printCMD {
	print STDERR "\n\tProgram to annotate de novo transcripts from GRO-seq data\n";
	print STDERR "\n\tUsage: annotateTranscripts.pl <transcript peak file> <genome> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-noRepeats (skip repeat annotation)\n";
	print STDERR "\t\t-minOverlap <#> (Min fraction of overlap to assign transcript as genic, 0.1)\n";
	print STDERR "\t\t-promoterSize <#> (to identify TSS antisense transcripts vs. enhancers, default: 500)\n";
	print STDERR "\t\t-repeatSize <#> (size from beginning of transcript to use for repeat ann, default: 500)\n";
	print STDERR "\t\t-d <tag directory> (to asses rpkm, merge fragments of same gene)\n";
	print STDERR "\t\t\t-min <#> (minimum rpkm, default: 0.01)\n";
	print STDERR "\t\t\t-minp <#> (minimum number of unique reads, default: 10)\n";
	print STDERR "\t\t\t-strand <+|-|both> (strand to search for reads, default: +)\n";
	print STDERR "\t\t\t-merge <#> (maximum fold difference for adjacent transcripts to be merged, e.g. 2)\n";
	#print STDERR "\t\t-tss <tss peak file> (adjust TSS to these called TSS locations)\n";
	#print STDERR "\t\t\t-tssRange <#> (TSS reassignment range, default: 1000)\n";

	print STDERR "\n";
	exit;
}


for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-noRepeats') {
		$noRepeats = 1;
	} elsif ($ARGV[$i] eq '-promoterSize') {
		$tssDistance = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-repeatSize') {
		$repeatCheckRegion = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minOverlap') {
		$minOverlap = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-d') {
		$tagDir = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-strand') {
		$strand = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$minRPKM = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minp') {
		$minP = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-tss') {
		$tssFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-merge') {
		$maxFold = $ARGV[++$i];
print STDERR "\tmergeFold = $maxFold\n";
	} else {
		printCMD();
	}
}
if ($tagDir eq '') {
	$minRPKM = -1e10;
	$minP = -1e10;
}

my $genomeDir = "";
my $customGenome = "";
my $organism = 'null';
my $preparsedDir = "";
if (!exists($config->{'GENOMES'}->{$genome})) {
    $customGenome = $genome;
    ($genome,$genomeDir,$preparsedDir) = HomerConfig::parseCustomGenome($genome);
} else {
    $genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'} . "/";
    $organism =  $config->{'GENOMES'}->{$genome}->{'org'};
}





my $denovo = openPeakFile($denovoFile);


# step 0- reannotate TSS using tss File
if ($tssFile ne '') {
	print STDERR "\tReassigning TSS within range of $tssRange...\n";
	my $tssPeaks = openPeakFile($tssFile);	
}

# first step, annotate transcripts relative to RNA file
print STDERR "\tAnnotating transcripts to reference genes...\n";
if ($refRNAfile eq '') {
	$refRNAfile = $genomeDir . '/' . $genome . ".rna";
}
unless (-e $refRNAfile) {
	print STDERR "!!!! Could not open refence RNA file \"$refRNAfile\"\n";
	exit;
}
my $repeatFlag = 0;
my $ref = openPeakFile($refRNAfile);
annotateTranscripts($denovo, $ref,$repeatFlag);

if ($tagDir ne '') {
	if ($maxFold > 0) {
		getRPKM($denovo, $tagDir,0);
		mergeAdjacentTranscripts($denovo,$maxFold);
	}
	getRPKM($denovo, $tagDir,0);
	if ($minP > 0) {
		getRPKM($denovo, $tagDir,1);
	}
}


# second step, annotate transcripts relative to 
if ($noRepeats == 0) {
	print STDERR "\tAnnotating non-gene transcripts to possibly repeats...\n";
	if ($repeatFile eq '') {
		$repeatFile = $genomeDir . '/' . $genome . ".repeats";
	}
	unless (-e $repeatFile) {
		print STDERR "!!!! Could not open repeat annotation file \"$repeatFile\"\n";
		exit;
	}
	$repeatFlag =1;
	$ref = openPeakFile($repeatFile);
	annotateTranscripts($denovo, $ref,$repeatFlag);
}

printTranscripts($denovo,$organism,"",1);
exit;

sub mergeAdjacentTranscripts {
	my ($denovo, $maxFold)  = @_;

	my $eps = 0.001;
	my $totalTranscripts = 0;
	my $totalMerged = 0;;
	foreach(keys %$denovo) {
		my $c = $_;
		my %gids = ();
		foreach(keys %{$denovo->{$c}}) {
			my $id = $_;
			$totalTranscripts++;
			my $gid = $denovo->{$c}->{$id}->{'gid'};
			my $type = $denovo->{$c}->{$id}->{'type'};
			#my $rpkm = $denovo->{$c}->{$id}->{'rpkm'};
			if ($type eq 'Gene') {
				if (!exists($gids{$gid})) {
					my @a = ();
					$gids{$gid} = \@a;
				}
				push(@{$gids{$gid}},$id);
			}
		}
		foreach(keys %gids) {
			my $gid = $_;
			my @ids = sort {$denovo->{$c}->{$a}->{'d'} cmp $denovo->{$c}->{$b}->{'d'} ||
								$denovo->{$c}->{$a}->{'s'} cmp $denovo->{$c}->{$b}->{'s'}}
											@{$gids{$gid}};
			my $lastGood = 0;
			for (my $i=1;$i<@ids;$i++) {
				my $id = $_;
				if ($denovo->{$c}->{$ids[$i]}->{'d'} ne $denovo->{$c}->{$ids[$lastGood]}->{'d'}) {
					$lastGood = $i;
					next;
				}
				my $rpkm = $denovo->{$c}->{$ids[$i]}->{'rpkm'};
				my $rpkmLast = $denovo->{$c}->{$ids[$lastGood]}->{'rpkm'};
				my $fold  = ($rpkm + $eps)/($rpkmLast+$eps);
				if ($fold < $maxFold && 1/$fold < $maxFold) {
					$denovo->{$c}->{$ids[$lastGood]}->{'e'} = $denovo->{$c}->{$ids[$i]}->{'e'};
					#print STDERR "merged $c:$denovo->{$c}->{$ids[$i]}->{'s'}" ."-$denovo->{$c}->{$ids[$i]}->{'e'}\n";
					delete $denovo->{$c}->{$ids[$i]};
					$totalMerged++;
				} else {
					$lastGood = $i;
					next;
				}
			}

		}
	}
	print STDERR "\tMerged $totalMerged out of $totalTranscripts transcripts\n";
}

sub assignTSS {
	my ($transcripts, $denovoFile, $tss, $tssFile, $tssRange) = @_;
	my $rand = rand();
	my $tmpFile = $rand . ".tmp";
	my $tmpFile2 = $rand . ".2.tmp";
	`adjustPeakFile.pl "$denovoFile" -5p > "$tmpFile"`;
	`mergePeaks "$tmpFile" "$tssFile"`;
}

sub getRPKM {
	my ($denovo, $tagDir,$pc) = @_;
	print STDERR "\tQuantifying Transcript levels with data in \"$tagDir\"\n";
	my $rand = rand();
	my $tmpFile = $rand . ".tmp";
	my $tmpFile2 = $rand . ".2.tmp";
	my $tmpRPKM = $minRPKM;
	$minRPKM = -1e10;
	my $tmpMINP = $minP;
	$minRPKM = -1e10;
	$minP = -1e10;
	printTranscripts($denovo,'null',$tmpFile,0);
	$minRPKM = $tmpRPKM;
	$minP = $tmpMINP;
	my $option = "";
	if ($pc > 0) {
		$option = " -pc $pc -noadj ";
	} else {
		$option = " -rpkm ";
	}
	`analyzeRepeats.pl "$tmpFile" $genome -noCondensing -strand $strand -d "$tagDir" $option > "$tmpFile2"`;
	open IN, $tmpFile2;
	while (<IN>) { 
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $line[0];
		my $chr = $line[1];
		my $rpkm = $line[8];
		next if (!exists($denovo->{$chr}));
		next if (!exists($denovo->{$chr}->{$id}));
		if ($pc > 0) {
			$denovo->{$chr}->{$id}->{'uniq'} = $rpkm;
		} else {
			$denovo->{$chr}->{$id}->{'rpkm'} = $rpkm;
		}
		
	}
	close IN;
	`rm "$tmpFile" "$tmpFile2"`;
}

sub printTranscripts {
	my ($denovo, $org,$outputFileName, $ucscFlag) = @_;

	my $emptyAnn = "";
	my %ann = ();
	if ($org ne '' && $org ne 'null') {
		$emptyAnn = "\t\t\t\t";
		my $rand = rand();
		my $tmpFile = $rand . ".tmp";
		my $tmpFile2 = $rand . ".2.tmp";
		open OUT, ">$tmpFile";
		foreach(keys %$denovo) {
			my $c = $_;
			foreach(values %{$denovo->{$c}}) {
				my $n = $_->{'gid'};
				print OUT "$n\n";
			}
		}
		`addGeneAnnotation.pl "$tmpFile" $org no > "$tmpFile2"`;
		open IN, $tmpFile2;
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			my $annStr = "\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[10]";
			$ann{$line[1]} = $annStr;
		}
		close IN;
		`rm "$tmpFile" "$tmpFile2"`;
	}

	if ($ucscFlag) {
		open UCSC, ">ucsc.txt";	
	}

	my $file = '';
	if ($outputFileName ne '') {
		open OUT, ">$outputFileName";
		$file = *OUT;
	} else {
		$file = *STDOUT;
	}

	print $file "TranscriptID\tchr\tstart\tend\tstrand\tOriginalID\tAnnotatedID\tType";
	print $file "\tName\tAlias\tOrf\tDescription\tRefSeqType";
	print $file "\tOverlap\tAntisense\tClosest TSS\tRPKM\tUniq Reads\n";

	my $totalPossible = 0;
	my $totalPassed = 0;

	my %pieChart = ();
	my %okTypes = ();
	$okTypes{'Gene'}=1;
	$okTypes{'Gene-antisense'}=1;
	$okTypes{'Promoter-antisense'}=1;
	$okTypes{'Enhancer'}=1;


	my %unique = ();
	my @chrs = sort {$a cmp $b} keys %$denovo;
	foreach(@chrs) {
		my $c = $_;
		my @ids = sort {$denovo->{$c}->{$a}->{'s'} <=> $denovo->{$c}->{$b}->{'s'}} keys %{$denovo->{$c}};

		for (my $i=0;$i<@ids;$i++) {
			$totalPossible++;

			my $id = $ids[$i];
			my $start = $denovo->{$c}->{$id}->{'s'};
			my $end = $denovo->{$c}->{$id}->{'e'};
			my $strand = $denovo->{$c}->{$id}->{'d'};
			my $nearest = $denovo->{$c}->{$id}->{'n'};
			my $nearestDist = $denovo->{$c}->{$id}->{'ndist'};
			my $gid = $denovo->{$c}->{$id}->{'gid'};
			my $type = $denovo->{$c}->{$id}->{'type'};
			my $rpkm = $denovo->{$c}->{$id}->{'rpkm'};
			my $uniq = $denovo->{$c}->{$id}->{'uniq'};
			next if ($rpkm < $minRPKM);
			next if ($uniq < $minP);
			$totalPassed++;

			my $nid ='';
			if ($type eq 'Gene') {
				$nid = $gid;
			} elsif ($type eq 'Gene-antisense') {
				$nid = $gid . "-gene-antisense";
			} elsif ($type eq 'Promoter-antisense') {
				$nid = $gid . '-promoter-antisense'; 
			} elsif ($type eq 'Enhancer') {
				$nid = $gid . '-enhancer'; 
			} elsif ($type eq 'other') {
				$nid = $gid . '-other'; 
			} elsif ($type eq 'other') {
				$nid = $gid . '-other'; 
			} elsif ($type =~ /^Repeat-(.*)/) {
				$nid = $gid . '-repeat-' . $1; 
			} else {
				$nid = $id . '-other';
			}
			if (exists($unique{$nid})) {
				$nid .= "-HOMER" . ++$unique{$nid};
			} else {
				$unique{$nid}=1;
			}
			
			print $file "$nid\t$c\t$start\t$end\t$strand\t$id\t$gid\t$type";

			my $pieType = $type;
			if (!exists($okTypes{$pieType})) {
				$pieType = 'other';
			}
			$pieChart{$pieType}++;

			my $ann = $emptyAnn;
			if (exists($ann{$gid})) {
				$ann = $ann{$gid};
			}
			print $file "$ann";


			my $name = "$id-";
			my $sense = '';
			my $antisense = '';	
			my @overlap = sort {$denovo->{$c}->{$id}->{'o'}->{$b} <=> $denovo->{$c}->{$id}->{'o'}->{$a}} 
											keys %{$denovo->{$c}->{$id}->{'o'}};
			foreach(@overlap) {
				my $rid = $_;
				my $overlap = $denovo->{$c}->{$id}->{'o'}->{$rid};
				$name .= $rid . "($overlap),";
				if ($overlap > 0) {
					$sense .= $rid . "($overlap),";
				} else {
					$antisense .= $rid . "($overlap),";
				}
			}
			print $file "\t$sense\t$antisense\t$nearest($nearestDist)\t$rpkm\t$uniq\n";
			$name .= "TSS:$nearest($nearestDist)";

			if ($ucscFlag) {
				print UCSC "$c\t$start\t$end\t$gid-$type\t1\t$strand\n";
			}


		}
	}
	if ($minRPKM > 0) {
		print STDERR "\t$totalPassed of $totalPossible passed RPKM threshold of $minRPKM\n";
	}
	if ($ucscFlag) {
		close UCSC;
	}
	if ($outputFileName ne '') {
		close OUT;
	}
	print STDERR "\n\tAnnotation Summary:\n";
	foreach(keys %pieChart) {
		print STDERR "\t$_\t$pieChart{$_}\n";
	}
	print STDERR "\n";
}



sub annotateTranscripts {
	my ($denovo, $ref, $repeatFlag) = @_;
	my @chrs = sort {$a cmp $b} keys %$denovo;
	foreach(@chrs) {
		my $c = $_;
		print STDERR "\t\t$c\n";
		my @ids = sort {$denovo->{$c}->{$a}->{'s'} <=> $denovo->{$c}->{$b}->{'s'}} keys %{$denovo->{$c}};
		my @rids = sort {$ref->{$c}->{$a}->{'s'} <=> $ref->{$c}->{$b}->{'s'}} keys %{$ref->{$c}};
		my $refIndex = 0;
	
		for (my $i=0;$i<@ids;$i++) {

			my $id = $ids[$i];
			my $start = $denovo->{$c}->{$id}->{'s'};
			my $end = $denovo->{$c}->{$id}->{'e'};
			my $strand = $denovo->{$c}->{$id}->{'d'};
			my $tss = $denovo->{$c}->{$id}->{'tss'};
			my $len = $denovo->{$c}->{$id}->{'len'};
			if ($repeatFlag) {
				if ($strand eq '+') {
					#$start -= $tssDistance;
					if ($end-$start>$repeatCheckRegion) {
						$end = $start+$repeatCheckRegion;
					}
				} else {
					#$end += $tssDistance;
					if ($end-$start>$repeatCheckRegion) {
						$start = $end-$repeatCheckRegion;
					}
				}
			}

			my $bestOverlap = 0;
			my $bestOverlapID = '';

			my $overlapFlag = 0;


			my $tssDist = 0;
			for (my $j=$refIndex;$j<@rids;$j++) {
				my $rid = $rids[$j];
				my $s = $ref->{$c}->{$rid}->{'s'};
				my $e = $ref->{$c}->{$rid}->{'e'};
				my $d = $ref->{$c}->{$rid}->{'d'};
				my $t = $ref->{$c}->{$rid}->{'tss'};
				my $l = $ref->{$c}->{$rid}->{'len'};

				my $dist = abs($tss-$t);
				if ($repeatFlag == 0 && ($dist < $denovo->{$c}->{$id}->{'ndist'})) {
					$denovo->{$c}->{$id}->{'ndist'} = $dist;
					$denovo->{$c}->{$id}->{'n'} = $rid;
				}

				my $overlap = 0;
				if ($s <= $end && $e >= $start) {
					#overlap
					my $ss = $start;
					$ss = $s if ($s > $ss);
					my $ee = $end;
					$ee = $e if ($e < $ee);
					$overlap = $ee-$ss;
#print STDERR "$id\t$rid\t$overlap\t$l\t$len\n" if ($id eq 'chr2-883-0');
					if ($overlap/$l > $minOverlap || $overlap/$len > $minOverlap) {
						if ($d ne $strand) {
							$overlap *= -1;
						}
						if ((($overlap > 0 || $bestOverlap > 0) && $overlap > $bestOverlap)
									|| ($overlap <= 0 && $bestOverlap <= 0 && $overlap < $bestOverlap)) {
							if ($repeatFlag == 0 || $overlap > 0) {
								$bestOverlap = $overlap;
								$bestOverlapID =$rid;
							}
						}
						if ($repeatFlag == 0) {
							$denovo->{$c}->{$id}->{'o'}->{$rid}=$overlap;
						}
					}
					$overlapFlag = 1;
				} else {
					if ($overlapFlag == 0 && $tss > $t && $j == $refIndex+1) {
						$refIndex = $j;
					}
					if ($denovo->{$c}->{$id}->{'ndist'}+$tss+100 < $s) {
						last;
					}
				}
			}			
			if ($bestOverlapID ne "") {
				if ($repeatFlag) {
					if ($denovo->{$c}->{$id}->{'type'} ne "Gene") {
						$repeatName = $bestOverlapID;
						$repeatName =~ s/\-HOMER\d+//;
						
						$denovo->{$c}->{$id}->{'type'} = "Repeat-" . $repeatName;
						if ($bestOverlap < 0) {
							$denovo->{$c}->{$id}->{'type'} = "Repeat-" . $repeatName . "-antisense";
						}
					}
				} else {
					$denovo->{$c}->{$id}->{'gid'} = $bestOverlapID;
					if ($denovo->{$c}->{$id}->{'o'}->{$bestOverlapID} < 0) {
						$denovo->{$c}->{$id}->{'type'} = "Gene-antisense";
					} else {
						$denovo->{$c}->{$id}->{'type'} = "Gene";
					}
				}
			} elsif ($repeatFlag == 0) {
				$denovo->{$c}->{$id}->{'gid'} = $denovo->{$c}->{$id}->{'n'};
				if ($denovo->{$c}->{$id}->{'ndist'} < $tssDistance) {
					if ($ref->{$c}->{$denovo->{$c}->{$id}->{'n'}}->{'d'} ne $denovo->{$c}->{$id}->{'d'}) {
						$denovo->{$c}->{$id}->{'type'} = "Promoter-antisense";
					} else {
						$denovo->{$c}->{$id}->{'type'} = "Promoter-sense";
					}
				} else {
					$denovo->{$c}->{$id}->{'type'} = "Enhancer";
				}
			}
		}
	}
}



sub openPeakFile {
	my ($file) = @_;
	my %peaks = ();
	open IN, $file;
	while (<IN>) {
		chomp;	
		s/\r//g;
		next if (/^\s*\#/);
		my @line = split /\t/;
		next if (@line < 5);
		next unless ($line[2] =~ /^[\d\.\-\e\+]+$/);
        next unless ($line[3] =~ /^[\d\.\-\e\+]+$/);
		my $id = $line[0];
		my $chr = $line[1];
		my $s = $line[2];	
		my $e = $line[3];
		my $strand = $line[4];
		$strand = "+" if ($strand eq 0);
		$strand = "-" if ($strand eq 1);
		if (!exists($peaks{$chr})) {
			my %a = ();
		}
		my $tss = $s;
		if ($strand eq '-') {
			$tss = $e
		}
		my $len = $e-$s;
		$len = 1 if ($len < 1);
		my $tid = $id;
		my $count = 1;
		while (exists($peaks{$chr}->{$tid})) {
			$tid .= "-$count";
			$count++;
		}
		my %overlap = ();
		my %near = ();
		my $nearID = '';
		$peaks{$chr}->{$tid} = {id=>$id,s=>$s,e=>$e,d=>$strand,c=>$chr,tss=>$tss,o=>\%overlap,n=>$nearID,ndist=>1e10,
								gid=>'',code=>'',cs=>$s,ce=>$e,len=>$len,rpkm=>0,uniq=>0};
	}
	close IN;
	return \%peaks;	
}

