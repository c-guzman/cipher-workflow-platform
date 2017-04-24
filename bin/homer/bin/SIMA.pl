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
use Statistics;
use HomerConfig;
use Fcntl qw(:flock);

sub printCMD {
	print STDERR "\n\tNormal usage: SIMA.pl <HIC directory> [options]\n";
	print STDERR "\n\tOutput table is sent to stdout\n";
	print STDERR "\t\t(See below for output visualization formatting)\n";
	print STDERR "\n\tRequired Options:\n";
	print STDERR "\t\t-d <domain peak file> (Domains to perform analysis on)\n";
	print STDERR "\t\t-p <peak file1> [peak file2] ... (features to check for enrichment)\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-res <#> (resolution, default=2500)\n";
	print STDERR "\t\t-superRes <#> (super resolution/window size, default=10000)\n";
	print STDERR "\t\t-minDist <#> (minimum interaction read distance, default: 2x superRes)\n";
	print STDERR "\t\t-minDsize <#> (minimum domain size, default: 500000)\n";
	print STDERR "\t\t-min <#> (minimum distance between domains to test significance, default=-1)\n";
	print STDERR "\t\t-max <#> (maximum distance, set to -1 to allow inter-chr, default=1e9)\n";
	print STDERR "\t\t-chr <chromosome> (only analyze this chromosome, default: all)\n";
	print STDERR "\t\t-p2 <peak file1> [peak file2] ... (heterogenous peak enrichments)\n";
	print STDERR "\t\t-AvsA (All versus All, compare -p peaks against one another)\n";
	print STDERR "\t\t-N <#> (Number of randomizations per domain, default: 1000)\n";
	#print STDERR "\t\t-distribution <filename> (Output random iteration results)\n";
	#print STDERR "\t\t-o <outputFile> (default: output sent to stdout)\n";
	print STDERR "\t\t-rdepth (adjust randomization based on HiC read depth)\n";
	print STDERR "\t\t-cpu <#> (number of CPUs to use for analysis, default: 1)\n";
	print STDERR "\t\t-merge (merge results for all domains)\n";
	print STDERR "\n";

	print STDERR "\n\n\tOutput Visualization Formatting (Run SIMA first, then format the output)\n";
	print STDERR "\n\t\tMatrix Mode: Takes output and prints out a matrix for visualization\n";
	print STDERR "\t\t\tSIMA.pl -matrix <SIMA output from analysis> [options]\n";
	print STDERR "\n\t\tMatrix Mode Options: (defaults to resolution of 200000, output to stdout)\n";
	print STDERR "\t\t\t-stat <pvalue|ratio> (output stat for matrix mode, default: ratio)\n";
	print STDERR "\t\t\t-pvalue <#> (p-value cutoff to report the ratio, default: 0.01)\n";
	print STDERR "\t\t\t-minPeaks <#> (minimum number of peaks, default: none)\n";
	print STDERR "\t\t\t-res <#> (resolution of matrix, default=200000)\n";
	print STDERR "\t\t\t-p <peak file1> (features from initial analysis to show)\n";
	print STDERR "\t\t\t-p2 <peak file2> (features from initial analysis to show, if used/different)\n";
	print STDERR "\n\t\tCytoscape Mode: Takes output from single domain and prints files\n";
	print STDERR "\t\t\tSIMA.pl -cytoscape <SIMA output from analysis>\n";
	print STDERR "\t\t\t\t(output to \"cytoscape.filename.*\" files)\n";
	print STDERR "\t\t\t-dname <name> (domain name to show)\n";
	print STDERR "\t\t\t-dname2 <name2> (domain name to show, if different)\n";

	#print STDERR "\n\t\tRelative Mode: (treat whole genome as a single domain)\n";
	#print STDERR "\t\t\t-d \"relative\" (required to activate relative mode)\n";
	#print STDERR "\t\t\t-cytoscape <filename prefix> (cytoscape output files)\n";
	#print STDERR "\t\t\t-maxDist <#> (maximum distance for relative analysis, default: 1e6)\n";
	print STDERR "\n";
	exit;

}
my $localFlag =0;
my @peakFiles = ();
my @peakFiles2 = ();
my $res = 2500;
my $superRes = 10000;
my @domains = ();
my $minDist = -1;
my $maxDist = 1e9;
$rdepthFlag = 0;
my $onlyChr = "";
$allVsAllFlag = 0;
my $minDomainSize = 5e5;
my $matrixRes = 0;
my $Nrandom = 1000;
my $matrixFile = '';
my $pvalueCutOff = 0.01;
my $minPeaks = -1;
my $minInterDistance = -123454321;
$distributionFile = "";
my $outputStat = 'ratio';
my $maxCPUs = 1;
my $dname1 = '';
my $dname2 = '';
my $outputFile = "";
my $maxRelDist = 1e6;
my $relativeFlag = 0;
my $mergeFlag = 0;

if (@ARGV < 2) {
	printCMD();
}
my $ped = $ARGV[0];
my $startIndex = 1;
my $cytoscapeFlag= 0;
if ($ped eq '-matrix') {
	$matrixFile = $ARGV[1];
	$startIndex = 2;
	$res = 200000;
} elsif ($ped eq '-cytoscape') {
	$matrixFile = $ARGV[1];
	$startIndex = 2;
	$cytoscapeFlag=1;
}


my $setRes = 0;
my $setSuperRes = 0;
for (my $i=$startIndex;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-res') {
		$res = $ARGV[++$i];
		$setRes = 1;
	} elsif ($ARGV[$i] eq '-distribution') {
		$distributionFile = $ARGV[++$i];
		open DIST, ">$distributionFile";
	} elsif ($ARGV[$i] eq '-merge') {
		$mergeFlag = 1;
	} elsif ($ARGV[$i] eq '-AvsA') {
		$allVsAllFlag = 1;
	} elsif ($ARGV[$i] eq '-pvalue') {
		$pvalueCutOff = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minDist') {
		$minInterDistance = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-stat') {
		$outputStat = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-rdepth') {
		$rdepthFlag = 1;
	} elsif ($ARGV[$i] eq '-superRes') {
		$superRes = $ARGV[++$i];
		$setSuperRes = 1;
	} elsif ($ARGV[$i] eq '-min') {
		$minDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minPeaks') {
		$minPeaks = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minDsize') {
		$minDomainSize = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-N') {
		$Nrandom = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-o') {
		$outputFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-chr') {
		$onlyChr = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxDist') {
		$maxRelDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-max') {
		$maxDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-matrix') {
		$matrixRes = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-d') {
		$domainFile = $ARGV[++$i];
		if ($domainFile eq 'relative') {
			$relativeFlag = 1;
		}
	} elsif ($ARGV[$i] eq '-dname') {
		$dname1 = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-dname2') {
		$dname2 = $ARGV[++$i];
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
	} elsif ($ARGV[$i] eq '-p2') {
		my $bail = 0;
		while ($ARGV[++$i] !~ /^\-/) {
			push(@peakFiles2, $ARGV[$i]);
			if ($i>=@ARGV-1) {
				$bail = 1;
				last;
			}
		}
		last if ($bail == 1);
		$i--;
	} else {
		print STDERR "!!! What is \"$ARGV[$i]\"??\n";
		printCMD();
	}
}
if ($minInterDistance == -123454321) {
	$minInterDistance = 2*$superRes;
}

foreach(@peakFiles) {
	print STDERR "\t\t$_\n";
}

if ($allVsAllFlag) {
	my @newPeaks1 = ();
	my @newPeaks2 = ();
	for (my $i=0;$i<@peakFiles;$i++) {
		for (my $j=$i;$j<@peakFiles;$j++) {
			push(@newPeaks1, $peakFiles[$i]);
			push(@newPeaks2, $peakFiles[$j]);
		}
	}
	@peakFiles = @newPeaks1;
	@peakFiles2 = @newPeaks2;
}
if ($cytoscapeFlag) {
	if ($mergeFlag) {
		$dname1 = "merge";
		$dname2 = "merge";
	}
	printCytoscape($matrixFile,$dname1,$dname2);
	exit;
}
if ($matrixFile ne '') {
	my $stats = readSIMAstats($matrixFile);
	printMatrix($stats, $res, $onlyChr, \@peakFiles, \@peakFiles2,$outputStat,$pvalueCutOff,$minPeaks);
	exit;
}

if ($domainFile eq '' || @peakFiles < 1) {
	print STDERR "\t!!! Need a domain file and peak files!!!\n";
	exit;
}

if ($resSet && !$superResSet) {
    $superRes = $res;
}
if (!$resSet && $superResSet) {
    $res = $superRes;
}
print STDERR "\n\tres set to $res\n";
print STDERR "\tsuperRes set to $superRes\n\n";
my $possibleRes = HomerConfig::getHiCBgRes($ped,$superRes,$maxCPUs);

my $domains = readPeakFile($domainFile, $minDomainSize, $onlyChr);

my $cpus = 0;
my @pids = ();

my $rand = rand();


my $mergeStatsFile = $rand . ".mergeStats";
#open STATS, ">$mergeStatsFile";


my $tmpMatrixFile = $rand . ".matrix";
open MATRIX, ">$tmpMatrixFile";

		
flock(MATRIX,LOCK_EX);	
seek(MATRIX,0,2);

print MATRIX "#DomainName(1)\tchr(1)\tstart(1)\tend(1)\tDomainName(2)\tchr(2)\tstart(2)\tend(2)"
			. "\tPeakFile1\tPeakFile2\tNpx\tNpy\tp-value\tRatio\tPeakEnrichment\tRandEnrichment\n";
flock(MATRIX,LOCK_UN);
close MATRIX;
		

my @results = ();
for (my $i=0;$i<@$domains;$i++) {
	my $d1 = $domains->[$i];
	my @result = ();
	for (my $j=0;$j<@$domains;$j++) {
		push(@result, 0);
	}
	push(@results, \@result);

	for (my $j=$i;$j<@$domains;$j++) { 
		my $d2 = $domains->[$j];
		if ($maxDist > -1) {
			next if ($d1->{'c'} ne $d2->{'c'});
			next if ($d2->{'p'} - $d1->{'p'} >= $maxDist);
		}
		next if ($d1->{'c'} eq $d2->{'c'} && $d2->{'p'} - $d1->{'p'} < $minDist);

		my $rr = rand();
		my $pid = fork();
		$cpus++;
		if ($pid == 0) {
			compareDomains($ped,$d1,$d2,\@peakFiles,\@peakFiles2,$res, $superRes,$Nrandom,
									$minInterDistance,$rand . "$i" . "x" . "$j");
			exit(0);
		}
		push(@pids, $pid);
		if ($cpus >= $maxCPUs) {
			my $id = wait();
			$cpus--;
		}
	}
}
my $id = 0;
while ($id >= 0) {
    $id = wait();
    if ($id == -1) {
    } else {
    }
}
print STDERR "\n";


open IN, "$tmpMatrixFile";
while (<IN>) {
	print $_;
}
close IN;

my @mergeStatsArray = ();
for (my $i=0;$i<@peakFiles;$i++) {
	my @scores = ();
	my @n = ();
	my @randScores = ();
	my $s = {s=>\@scores,r=>\@randScores,n=>\@n};
	push(@mergeStatsArray, $s);
}
my $mergeStats = \@mergeStatsArray;

open IN, $mergeStatsFile;
while (<IN>) {
	chomp;
	my @line = split /\t/;
	my $i = $line[0];
	push(@{$mergeStats->[$i]->{'s'}}, [$line[1],$line[2]]);
	push(@{$mergeStats->[$i]->{'n'}}, [$line[3],$line[4]]);
	my @samples = ();
	for (my $j=5;$j<@line;$j+=2) {
		push(@samples, [$line[$j],$line[$j+1]]);
	}
	push(@{$mergeStats->[$i]->{'r'}}, \@samples);
}
close IN;

if ($mergeFlag) {
	my $oneFlag = 0;
	if (scalar(@peakFiles2) < 1) {
		$oneFlag = 1;
	}
	
	for (my $i=0;$i<@peakFiles;$i++) {
		my $file1= $peakFiles[$i];	
		my $file2= $peakFiles[$i];	
		if ($oneFlag==0) {
			$file2 = $peakFiles2[$i];
		}
		my @n = @{$mergeStats->[$i]->{'n'}};
		my $N = scalar(@n);
		my $obs = 0;
		my $exp = 0;
		my $robs = 0;
		my $rexp = 0;
		my @rand = ();
		for (my $k=0;$k<$Nrandom;$k++) {
			push(@rand,[0,0]);
		}
		my $Nx = 0;
		my $Ny = 0;
		for (my $j=0;$j<@n;$j++) {
			next if (!defined($mergeStats->[$i]->{'s'}->[$j]->[0]));
			$Nx += $mergeStats->[$i]->{'n'}->[$j]->[0];
			$Ny += $mergeStats->[$i]->{'n'}->[$j]->[1];
			$obs += $mergeStats->[$i]->{'s'}->[$j]->[0];
			$exp += $mergeStats->[$i]->{'s'}->[$j]->[1];
			for (my $k=0;$k<$Nrandom;$k++) {
				next if (!defined($mergeStats->[$i]->{'r'}->[$j]->[$k]->[0]));
				$rand[$k]->[0] += $mergeStats->[$i]->{'r'}->[$j]->[$k]->[0];
				$rand[$k]->[1] += $mergeStats->[$i]->{'r'}->[$j]->[$k]->[1];
				$robs += $mergeStats->[$i]->{'r'}->[$j]->[$k]->[0];
				$rexp += $mergeStats->[$i]->{'r'}->[$j]->[$k]->[1];
			}
		}
		my $oEnrich = ($obs+1)/($exp+1);
		my $eEnrich = ($robs+1)/($rexp+1);
		my $ratio = $oEnrich / $eEnrich;
		my $pvalue = 0;
		for (my $k=0;$k<$Nrandom;$k++) {
			my $s = ($rand[$k][0]+1)/($rand[$k][1]+1);
			$pvalue++ if ($s > $oEnrich);
		}
		$pvalue /= $Nrandom;
		
		print "merge\t\t\t\tmerge\t\t\t\t$file1\t$file2\t$Nx\t$Ny\t$pvalue\t$ratio\t$oEnrich\t$eEnrich\n";

	}
}
`rm -f "$mergeStatsFile"`;



exit;



#Funtions...................................................................


sub compareDomains {

	my ($ped, $d1, $d2, $peakFiles, $peakFiles2, $res, $superRes, $Nrandom, $minInterDistance,$rand) = @_;

	my $posStr1 = $d1->{'c'} . ":" . $d1->{'s'} . "-" . $d1->{'e'};
	my $posStr2 = $d2->{'c'} . ":" . $d2->{'s'} . "-" . $d2->{'e'};
	print STDERR "\t\t$posStr1 vs. $posStr2\n";

	open STATS, ">>$mergeStatsFile";
	open MATRIX, ">>$tmpMatrixFile";

	my $minIndexSeparation = -1;
	if ($posStr1 eq $posStr2) {
		$minIndexSeparation = $minInterDistance/$res;
	}
	
	my $tmpFile = $rand . ".tmp";
	my $tmpFile2 = $rand . ".2.tmp";
	my $tmpFile3 = $rand . ".3.tmp";
	my $tmpFile4 = $rand . ".4.tmp";
	my $tmpFile5 = $rand . ".5.tmp";

	`analyzeHiC "$ped" -pos $posStr1 -pos2 $posStr2 -res $res -superRes $superRes -rawAndExpected "$tmpFile2" -min -1 -std 10000 -pout "$tmpFile4" -pout2 "$tmpFile5" > "$tmpFile" 2> /dev/null`;
	#`analyzeHiC "$ped" -pos $posStr1 -pos2 $posStr2 -res $res -superRes $superRes -rawAndExpected "$tmpFile2" -min -1 -std 10000 -pout "$tmpFile4" -pout2 "$tmpFile5" > "$tmpFile"`;
	#print STDERR "`analyzeHiC $ped -pos $posStr1 -pos2 $posStr2 -res $res -superRes $superRes -rawAndExpected $tmpFile2 -min -1 -std 10000 -pout $tmpFile4 -pout2 $tmpFile5 > $tmpFile 2> /dev/null`;\n";

	my $matrix = {x=>'',y=>''};

	readMatrix($matrix, $tmpFile, 'raw',$minIndexSeparation,"","");
	readMatrix($matrix, $tmpFile2, 'expected',$minIndexSeparation,$tmpFile4,$tmpFile5);
	#readMatrix($matrix, $tmpFile, 'raw',$minIndexSeparation,$tmpFile4,$tmpFile5);
	#readMatrix($matrix, $tmpFile2, 'expected',$minIndexSeparation,$tmpFile4,$tmpFile5);
	`rm "$tmpFile" "$tmpFile2" "$tmpFile4" "$tmpFile5"`;

	my $Nx = @{$matrix->{'x'}};
	my $Ny = @{$matrix->{'y'}};
	my $avg = $matrix->{'raw-Total'}/$matrix->{'expected-Total'};

	my $oneFlag = 0;
	if (scalar(@$peakFiles2) < 1) {
		$oneFlag = 1;
	}


	for (my $i=0;$i<@$peakFiles;$i++) {
	

		my $file1= $peakFiles->[$i];	
		my $file2= $peakFiles->[$i];	
		if ($oneFlag==0) {
			$file2 = $peakFiles2->[$i];
		}
		my @X = readFeaturePeakFile($matrix,$file1,'X');
		my @Y = readFeaturePeakFile($matrix,$file2,'Y');

		my $score = scoreMatrix($matrix,$X[0],$Y[1],$minIndexSeparation);
		my $Npx = @{$X[0]};
		my $Npy = @{$Y[1]};

		#push(@{$mergeStats->[$i]->{'s'}}, $score);
		#push(@{$mergeStats->[$i]->{'n'}}, [$Npx,$Npy]);
		my @randScores = ();	

		#print STDERR "\t$Matrix = $Nx by $Ny  (peaks: $Npx by $Npy)\n";
		#print STDERR "\tMatrix = $Nx by $Ny  (peaks: $Npx by $Npy) avg = $avg ($score->[0] $score->[1] $score->[2])\n";

		my $VV = (($score->[0]+1)/($score->[1]+1));
		my $numBetter = 0;
		my $avgRnd = 0;
		#print STDERR "\t\tRandomizing:";
		for (my $j=0;$j<$Nrandom;$j++) {
			if ($j % 1000 == 0) {
		#		print STDERR "\t$j";
			}
			my $xx = '';
			my $yy = '';
			if ($rdepthFlag) {
				$xx = randomArrayDepth($Npx,$Nx, $matrix->{'pX'});
				$yy = randomArrayDepth($Npy,$Ny, $matrix->{'pY'});
			} else {
				$xx = randomArray($Npx,$Nx);
				$yy = randomArray($Npy,$Ny);
			}
			my $scores= scoreMatrix($matrix,$xx,$yy,$minIndexSeparation);
			push(@randScores, $scores);
			my $v = (($scores->[0]+1)/($scores->[1]+1));
			if ($distributionFile ne '') {
				print DIST "$v\t$scores->[0]\t$scores->[1]\n";
			}
			if ($v >= $VV) {
				$numBetter++;
			}
			$avgRnd += $v;
		}

		my $numBetterRatio = $numBetter/$Nrandom;
		$avgRnd  /= $Nrandom;
		$ratioScore = 1;
		if ($avgRnd != 0) {
			$ratioScore = $VV/$avgRnd;
		}
		#print STDERR "$ratioScore (p=$numBetterRatio)\n";
	
		flock(STDOUT,LOCK_EX);	
		seek(STDOUT,0,2);
		print STDOUT "$d1->{'id'}\t$d1->{'c'}\t$d1->{'s'}\t$d1->{'e'}\t$d2->{'id'}\t$d2->{'c'}\t$d2->{'s'}"
					. "\t$d2->{'e'}\t$file1\t$file2\t$Npx\t$Npy\t$numBetterRatio\t$ratioScore\t$VV\t$avgRnd\n";
		flock(STDOUT,LOCK_UN);


		flock(STATS,LOCK_EX);	
		seek(STATS,0,2);
		print STATS "$i\t$score->[0]\t$score->[1]\t$Npx\t$Npy";
		foreach(@randScores) {
			print STATS "\t$_->[0]\t$_->[1]";
		}
		print STATS "\n";
		flock(STATS,LOCK_UN);

		
	}
	close STATS;
	close MATRIX;
}




sub printMatrix {
	my ($results, $matrixRes, $onlyChr, $peakFiles1, $peakFiles2, $outputStat,$pvalueCutOff,$minPeaks) = @_;

	my @domains = values %$results;
	@domains = sort peakcmp @domains;

	print STDERR "\tMatrix resolution = $matrixRes\n";

	my @matrix = ();
	my $curChr = "";
	my $curPos = 0;
	print "Domains\tDomains";
	my @domainMap = ();
	my @pos = ();
	for (my $i=0;$i<@domains;$i++) {
		#print STDERR "$domains[$i]->{'id'} $domains[$i]->{'c'}\t$domains[$i]->{'p'}\n";
		next if ($onlyChr ne '' && $domains[$i]->{'c'} ne $onlyChr);
		if ($domains[$i]->{'c'} ne $curChr) {
			$curChr = $domains[$i]->{'c'};
			$curPos = 0;
		}
		for (;$curPos<$domains[$i]->{'s'};$curPos+=$matrixRes) {
			print "\t$curChr-$curPos";
			push(@pos, "$curChr-$curPos");
			push(@domainMap, "NA");
		}
		for (;$curPos<$domains[$i]->{'e'};$curPos+=$matrixRes) {
			print "\t$curChr-$curPos";
			push(@pos, "$curChr-$curPos");
			push(@domainMap, $i);
		}
	}
	print "\n";
	for (my $i=0;$i<@pos;$i++) {
		print "$pos[$i]\t$pos[$i]";
		my $dindex1 = $domainMap[$i];
		for (my $j=0;$j<@pos;$j++) {
			my $dindex2 = $domainMap[$j];
			if ($dindex1 eq 'NA' || $dindex2 eq 'NA') {
				print "\tNA";
				next;
			}
			my $dname2 = $domains[$dindex2]->{'id'};
			my $N=0;
			my $v = 0;
			for (my $k=0;$k<@$peakFiles1;$k++) {
				my $file1 = $peakFiles1->[$k];
				my $file2 = $file1;
				if (scalar(@$peakFiles2) > 0) {
					$file2 = $peakFiles2->[$k];
				}
				my $pp = $file1 . "x" . $file2;
		
				next if (!exists($domains[$dindex1]->{'r'}->{$dname2}));
				next if (!exists($domains[$dindex1]->{'r'}->{$dname2}->{$pp}));
				my $r = $domains[$dindex1]->{'r'}->{$dname2}->{$pp};
				next if ($r->{'np1'} < $minPeaks || $r->{'np2'} < $minPeaks);
				if ($outputStat eq 'pvalue') {
					$v += log($r->{'pvalue'}+0.00005)/log(2);
					$N++;
				} elsif ($outputStat eq 'ratio') {
					my $vvv = $r->{'ratio'};
					my $bad = 0;
					if ($vvv > 1) {
						if ($r->{'pvalue'} > $pvalueCutOff) {
							$bad = 1;
						}
					} else {
						if ($r->{'pvalue'} < 1-$pvalueCutOff) {
							$bad = 1;
						}
					}
					if ($bad == 0) {
						$v += $vvv;
						$N++;
					} else {
						$v = 1;
						$N++;
					}
				}
			}
			if (@$peakFiles1 < 1) { #$N == 0) {
				foreach(values %{$domains[$dindex1]->{'r'}->{$dname2}}) {
					if ($outputStat eq 'pvalue') {
						$v += log($_->{'pvalue'}+0.00005)/log(2);
						$N++;
					} elsif ($outputStat eq 'ratio') {
						my $vvv = $_->{'ratio'};
						my $bad = 0;
						if ($vvv > 1) {
							if ($_->{'pvalue'} > $pvalueCutOff) {
								$bad = 1;
							}
						} else {
							if ($_->{'pvalue'} < 1-$pvalueCutOff) {
								$bad = 1;
							}
						}
						if ($bad == 0) {
							$v += $vvv;
							$N++;
						} else {
							$v = 1;
							$N++;
						}
					}
				}
			}
			if ($N < 1) {
				#print STDERR "Might have a problem....\n";
				$v = 'NA';
			} else {
				$v /= $N;
			}
			print "\t$v"
			
		}
		print "\n";
	}
	
	

}

sub readSIMAstats {
	my ($file) = @_;

	my %domains = ();
	open IN, $file or die "Could not open file \"$file\"\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		my @line = split /\t/;
		my $id1 = $line[0];
		my $id2 = $line[4];
		my $index1=0;
		my $index2=0;
		if (!exists($domains{$id1})) {
			my %a = ();
			my $p = {id=>$id1,c=>$line[1],s=>$line[2],e=>$line[3],r=>\%a,p=>floor(($line[2]+$line[3])/2)};
			$domains{$id1}=$p;
		}
		if (!exists($domains{$id2})) {
			my %a = ();
			my $p = {id=>$id2,c=>$line[5],s=>$line[6],e=>$line[7],r=>\%a,p=>floor(($line[6]+$line[7])/2)};
			$domains{$id2}=$p;
		}
		my $pp = $line[8] . "x" . $line[9];
		my $values = {p1=>$line[8],p2=>$line[9],np1=>$line[10],np2=>$line[11],pvalue=>$line[12],
										ratio=>$line[13],vv=>$line[14],avgRnd=>$line[15]};

		if (!exists($domains{$id1}->{'r'}->{$id2})) {
			my %a = ();
			$domains{$id1}->{'r'}->{$id2} = \%a;
		}
		$domains{$id1}->{'r'}->{$id2}->{$pp} = $values;
		if (!exists($domains{$id2}->{'r'}->{$id1})) {
			my %a = ();
			$domains{$id2}->{'r'}->{$id1} = \%a;
		}
		$domains{$id2}->{'r'}->{$id1}->{$pp} = $values;
	}
	close IN;
	return \%domains;
}

sub randomArray {
	my ($n, $N) = @_;
	my %opt = ();
	my @array = ();
	for (my $i=0;$i<$n;$i++) {
		my $x = floor(rand()*$N);
		if (exists($opt{$x})) {
			$i--;
			next;
		}
		$opt{$x} = 1;
		push(@array, $x);
	}
	return \@array;
}
sub randomArrayDepth {
	my ($n, $N, $px) = @_;
	my @peakIndex = ();
	my $curr = 0;
	my %chosen = ();
	while ($curr < $n) {
		my @rand = ();
		for (my $i=$curr;$i<$n;$i++){ 
			push(@rand, rand());
		}
		@rand = sort {$a <=> $b} @rand;
		my $fsum = 0;
		my $index = 0;
		for (my $i=0;$i<@rand;$i++) {
			my $r = $rand[$i];
			$index++ while ($index < scalar(@$px) && $r>$px->[$index]);
			last if ($index >= @$px);
			if (!exists($chosen{$index})) {
				push(@peakIndex,$index);
				$chosen{$index}=1;
#print STDERR "$n $N $index\n";
			}
		}
		$curr = scalar(@peakIndex);
		last if ($curr>=$n);
	}
	return \@peakIndex;
}


sub scoreMatrix {
	my ($matrix, $arrayX, $arrayY, $minIndexSeparation) = @_;


	my @totals = (0,0,0);

	my @types = ("raw", "expected", "norm");
	for (my $z = 0;$z<@types;$z++) {
		my $type = $types[$z];
		next if (!exists($matrix->{$type}));
		my $m = $matrix->{$type};
		for (my $i=0;$i<@$arrayX;$i++) {
			my $x = $arrayX->[$i];
			for (my $j=0;$j<@$arrayY;$j++) {
				my $y = $arrayY->[$j];
				next if (abs($x-$y) < $minIndexSeparation);
				if (!defined($m->[$y][$x])) {
					print STDERR "!! $type\t$y\t$x\n";
				}
				$totals[$z] += $m->[$y][$x];
			}
		}
	}
	return \@totals;
}

sub readPeakFile {
	my ($file, $msize, $onlyChr) =@_;
	my $tmpFile = rand() . ".tmp." . $file;
	$tmpFile =~ s/\///g;
	`bed2pos.pl "$file" -check > "$tmpFile" 2> /dev/null`;
	my @data = ();
	open IN, $tmpFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		my @line = split /\t/;
		next if (@line < 5);
		next if ($line[1] =~ /random/);
		next if ($onlyChr ne '' && $onlyChr ne $line[1]);
		next if ($line[3]-$line[2] < $msize);
		my $m = floor(($line[2]+$line[3])/2);
		my $v = 0;
		if (@line > 5) {
			$v = $line[5];
		}
		my $p = {id=>$line[0],c=>$line[1],p=>$m,s=>$line[2],e=>$line[3],d=>$line[4],v=>$v};
		push(@data, $p);
	}
	close IN;
	`rm "$tmpFile"`;
	@data = sort peakcmp @data;
	return \@data;
}



sub readFeaturePeakFile {
	my ($matrix, $file, $mode) = @_;

	my %diff = ();
	
	my @x = ();
	my @y = ();
	my @xpos = ();
	my @ypos = ();
	my $last = 0;
	for (my $i=0;$i<@{$matrix->{'x'}}; $i++) {
		my $name = $matrix->{'x'}->[$i];
		$name=~ /(.+)\-(\d+)/;
		my $p = {c=>$1,p=>$2};
		if ($i>0) {
			$diff{$2-$last}++;
		}
		$last =$2;
		push(@xpos, $p);
		push(@x,0);
	}

	my @ress = sort {$diff{$b} <=> $diff{$a}} keys %diff;
	my $res = $ress[0];
	#print STDERR "\tEstimated Resolution = $res\n";

	for (my $i=0;$i<@{$matrix->{'y'}}; $i++) {
		my $name = $matrix->{'y'}->[$i];
		$name=~ /(.+)\-(\d+)/;
		my $p = {c=>$1,p=>$2};
		push(@ypos, $p);
		push(@y,0);
	}

	my @peaks = ();
	open IN, $file or die "Could not open $file\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		my @line = split /\t/;
		next if ($line[1] =~ /random/);
		my $p = {c=>$line[1],p=>floor(($line[2]+$line[3])/2)};
		push(@peaks,$p);
	}
	close IN;

	@peaks = sort peakcmp @peaks;
	#@ypos = sort peakcmp @ypos;
	#print "Assigning peaks\n";
	if ($mode eq 'X' || $mode eq 'XY') {
		my $index = 0;
		for (my $i=0;$i<@peaks;$i++) {
			while ($index < @xpos && peakcmp2($peaks[$i],$xpos[$index]) > 0) {
				$index++;
			}
			last if ($index >= @xpos);
			for (my $j=$index;$j<@xpos;$j++) {
				last if ($peaks[$i]->{'c'} ne $xpos[$j]->{'c'});
				last if ($peaks[$i]->{'p'} < $xpos[$j]->{'p'});
				if ($peaks[$i]->{'p'} > $xpos[$j]->{'p'}+$res) {
					$index++;
					next;
				}
				if ($peaks[$i]->{'p'} >= $xpos[$j]->{'p'} && $peaks[$i]->{'p'} < $xpos[$j]->{'p'}+$res) {
					#print STDERR "X: $j\t$xpos[$j]->{'p'}\t$peaks[$i]->{'p'}\n";
					$x[$j]=1;
					last;
				}
			}
		}
	}

	if ($mode eq 'Y' || $mode eq 'XY') {
		my $index = 0;
		for (my $i=0;$i<@peaks;$i++) {
			while ($index < @ypos && peakcmp2($peaks[$i],$ypos[$index]) > 0) {
				$index++;
			}
			last if ($index >= @ypos);
			for (my $j=$index;$j<@ypos;$j++) {
				last if ($peaks[$i]->{'c'} ne $ypos[$j]->{'c'});
				last if ($peaks[$i]->{'p'} < $ypos[$j]->{'p'});
				if ($peaks[$i]->{'p'} > $ypos[$j]->{'p'}+$res) {
					$index++;
					next;
				}
				if ($peaks[$i]->{'p'} >= $ypos[$j]->{'p'} && $peaks[$i]->{'p'} < $ypos[$j]->{'p'}+$res) {
					#print STDERR "Y: $j\t$ypos[$j]->{'p'}\t$peaks[$i]->{'p'}\n";
					$y[$j]=1;
					last;
				}
			}
		}
	}

	my @xx = ();
	my @yy = ();
	for (my $i=0;$i<@x;$i++) {
		if ($x[$i] > 0) {
			push(@xx, $i);
		}
	}
	for (my $i=0;$i<@y;$i++) {
		if ($y[$i] > 0) {
			push(@yy, $i);
		}
	}

	return (\@xx,\@yy);

}

sub peakcmp {
	my $c1 = $a->{'c'};
	my $c2 = $b->{'c'};
	my $p1 = $a->{'p'};
	my $p2 = $b->{'p'};
	$c1 =~ s/^chr//;
	$c2 =~ s/^chr//;
	my $rv = 0;
	if ($c1 =~ /\d/ && $c2 =~ /\d/) {
		$rv = $c1 <=> $c2;
	} elsif ($c1 =~ /\d/) {
		$rv = -1;
	} elsif ($c2 =~ /\d/) {
		$rv =  1;
	} else {
		$rv = $c1 cmp $c2;
	}
	if ($rv == 0) {
		$rv = $p1 <=> $p2;
	}
	return $rv;
}
sub peakcmp2 {
	my ($a,$b) = @_;
	my $c1 = $a->{'c'};
	my $c2 = $b->{'c'};
	my $p1 = $a->{'p'};
	my $p2 = $b->{'p'};
	$c1 =~ s/^chr//;
	$c2 =~ s/^chr//;
	my $rv = 0;
	if ($c1 =~ /\d/ && $c2 =~ /\d/) {
		$rv = $c1 <=> $c2;
	} elsif ($c1 =~ /\d/) {
		$rv = -1;
	} elsif ($c2 =~ /\d/) {
		$rv = 1;
	} else {
		$rv = $c1 cmp $c2;
	}
	if ($rv == 0) {
		#$rv = $p1 <=> $p2;
	}
	return $rv;
}

sub readMatrix {
	my ($matrix,$file,$type,$minSep,$poutFile,$pout2File) = @_;
	open IN, $file or print STDERR "Could not open $file\n";
	my $count = 0;
	my @x = ();
	my @y = ();
	my $total = 0;
	my @data = ();
	my $X = 0;
	my $Y = 0;
	while (<IN>) {
		$count++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($count == 1) {
			shift @line;
			shift @line;
			@x = @line;
			$X = scalar(@x);
			next;
		}
		push(@y, shift @line);
		shift @line;
		for (my $i=0;$i<@line;$i++) {
			next if (abs($Y-$i) < $minSep);
			$total += $line[$i];
		}
		push(@data, \@line);
		$Y++;
	}
	close IN;
	my @pX = ();
	my @pY = ();
	if (-e $pout2File) {
		my $pp =  readPeakFile($pout2File,0,"");
		my $totalCoverage = 0;
		foreach(@$pp) {
			$totalCoverage += $_->{'v'};
		}
		my $sum = 0;
		foreach(@$pp) {
			$sum+=$_->{'v'}/$totalCoverage;
			push(@pX, $sum);
		}
	} else {
		for (my $i=0;$i<$X;$i++){ 
			push(@pX, $i/$X);
		}
	}
	if (-e $poutFile) {
		my $pp =  readPeakFile($poutFile,0,"");
		my $totalCoverage = 0;
		foreach(@$pp) {
			$totalCoverage += $_->{'v'};
		}
		my $sum = 0;
		foreach(@$pp) {
			$sum+=$_->{'v'}/$totalCoverage;
			push(@pY, $sum);
		}
	} else {
		for (my $i=0;$i<$Y;$i++){ 
			push(@pY, $i/$Y);
		}
	}

	$matrix->{'pX'} = \@pX;
	$matrix->{'pY'} = \@pY;
	$matrix->{'x'} = \@x;
	$matrix->{'y'} = \@y;
	$matrix->{$type} = \@data;
	$matrix->{$type."-Total"} = $total;
}



sub printCytoscape {


	my ($inputFile, $dname, $dname2) = @_;


	my $iname =$inputFile;
	$iname =~ s/\//\_/g;
	my $prefix = 'cytoscape.' . $iname;
	open NETWORK, ">$prefix.network.txt";
	open NODESIZE, ">$prefix.nodes.size.txt";
	open EDGERATIO, ">$prefix.edges.ratio.txt";
	open EDGEPVALUE, ">$prefix.edges.pvalue.txt";
	print NODESIZE "NodeSize\n";
	print EDGERATIO "EdgeRatio\n";
	print EDGEPVALUE "EdgePValue\n";
	my %done = ();
	open IN, $inputFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		my @line = split /\t/;
		if ($dname2 ne '') {
			next if ($line[4] ne $dname2);
		}
		if ($dname ne '') {
			next if ($line[0] ne $dname);
			next if ($dname2 eq '' && $line[0] ne $dname);
		}
		my $peak1 = cleanPeaks($line[8]);
		my $peak2 = cleanPeaks($line[9]);
		if (!exists($done{$peak1})) {
			print NODESIZE "$peak1 = $line[10]\n";
			$done{$peak1} = 1;
		}
		if (!exists($done{$peak2})) {
			print NODESIZE "$peak2 = $line[11]\n";
			$done{$peak2} = 1;
		}
		print NETWORK "$peak1\tpp\t$peak2\n";
		my $logR = $line[13];
		if ($logR > 0) {
			$logR = log($logR)/log(2);
		}
		my $logp = log($line[12]+0.0001)/log(2);
		if ($line[12] > 0.5) {
			$logp = -1*log((1-$line[12])+0.0001)/log(2);
		}
		if ($logR == 0) {
			$logR = '0.0';
		}
		if ($logp == 0) {
			$logp = '0.0';
		}
		print EDGERATIO "$peak1 (pp) $peak2 = $logR\n";
		print EDGEPVALUE "$peak1 (pp) $peak2 = $logp\n";
	}
	close IN;
	close NETWORK;
	close EDGERATIO;
	close NODESIZE;
	return;
}
sub cleanPeaks {
	my ($name) = @_;
	$name =~ s/\.\.\///g;
	$name =~ s/\/peaks\.txt//;
	$name =~ s/\/regions\.txt//;
	$name =~ s/\.txt//;
	$name =~ s/Bcell.*ko-//;
	$name =~ s/Ig\///;
	$name =~ s/-\d+//;
	return $name;
}
