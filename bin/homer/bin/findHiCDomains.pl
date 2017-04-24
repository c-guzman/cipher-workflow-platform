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

my $stdFilter = 4;
my $minFilter = 0.15;
my $res = 5000;
my $superRes = 25000;
my $corrDepth = 3;
my $activeZones = "";
my $inActiveZones = "";
my $maxCPUs = 1;
my $domainDist = 1000000;
my $minDist = 0;
my $nolog = '-nolog';
my $normType = '-simpleNorm';
$window = 25000;
$minDomainLength = 25000;
$maxError = 0.25;
$minBias = 0.50;
$minExpBias = 0.00001;
$minDelta = 1;

sub printCMD() {
	print STDERR "\n\tUsage findHiCDomains.pl <output prefix> <HiC directory> [options]\n";
	print STDERR "\t -or- findHiCDomains.pl <output prefix> <directionality index bedgraph> [options]\n";
	print STDERR "\t\t(Use the 2nd usage to change parameters for domain calls after running)\n";
	print STDERR "\n\tDirectionality index calculation options (used to call domains):\n";
	print STDERR "\t\t-res <#> (resolution in bp, default: $res)\n";
	print STDERR "\t\t-superRes <#> (super resolution in bp, i.e. window size, default: $superRes)\n";
	print STDERR "\t\t-maxDist <#> (max distance from loci to check direction index, default: $domainDist)\n";
	print STDERR "\t\t-minDist <#> (minimum distance from loci to check direction index, default: $minDist)\n";
	print STDERR "\t\t-log (score direction index with log ratios, default: use linear ratios)\n";
	print STDERR "\t\t-std <#> (exclude regions with sequencing depth exceeding # std deviations, default: $stdFilter)\n";
	print STDERR "\t\t-min <#> (exclude regions with sequencing depth less than this fraction of mean, default: $minFilter)\n";
	print STDERR "\t\t-cpu <#> (number of CPUs to use, default: 1)\n";


	print STDERR "\n\tDomain calling options:\n";
	print STDERR "\t\t-window <#> (smoothing window, removes noise in directionality index, default: $window)\n";
	print STDERR "\t\t-minIndex <#> (minimum index score to be considered for boundary (as zscore), default: $minBias)\n";
	print STDERR "\t\t-minExpIndex <#> (minimum index score to be considered for boundary (as zscore), default: $minExpBias)\n";
	print STDERR "\t\t-minDelta <#> (minimum change in dir-index between boundaries (as zscore), default: $minDelta)\n";
	print STDERR "\t\t-minLength <#> (minimum domain length, default: $minDomainLength)\n";
	print STDERR "\t\t-maxError <#> (maximum variation in directionality index within domain, default: $maxError)\n";
	print STDERR "\n\tOutput files:\n";
	print STDERR "\t\t<prefix>.directionIndex.bedGraph - UCSC upload file showing directionality index\n";
	print STDERR "\t\t<prefix>.smoothedIndex.bedGraph - UCSC upload file showing smoothed directionality index\n";
	print STDERR "\t\t<prefix>.domains.bed - BED file showing domain calls\n";
	print STDERR "\t\t<prefix>.innerDomains.bed - BED file showing strong inner domain calls\n";
	print STDERR "\t\t<prefix>.boundaries.txt - peak file centered on boundaries facing domain interiors\n";
	print STDERR "\n";
	exit;
}


if (@ARGV < 2) {
	printCMD();
}

my $prefix = $ARGV[0];
my $directory = $ARGV[1];
my $rpath = "R";
my $numPC=1;
for (my $i=2;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-res') {
		$res = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-window') {
		$window = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minDelta') {
		$minDelta = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minIndex') {
		$minBias = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minExpIndex') {
		$minExpBias = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxError') {
		$maxError = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-superRes') {
		$superRes = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minLength') {
		$minDomainLength = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-log') {
		$nolog = '';
	} elsif ($ARGV[$i] eq '-minDist') {
		$minDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-domainDist' || $ARGV[$i] eq '-maxDist') {
		$domainDist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-std') {
		$stdFilter = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$minFilter = $ARGV[++$i];
	} else {
		print STDERR "!!! Didn't recognize option \"$ARGV[$i]\"!!!!\n";
		printCMD();
	}
}

if ($res > $superRes) {
	print STDERR "!!! Warning -superRes ($superRes) is smaller than -res ($res)...\n\t\tPausing for 5s...";
	`sleep 5`;
}

print STDERR "\tOutput files will start with: $prefix\n";
my $outFile = $prefix . ".directionIndex.bedGraph";
my $fileFlag = 0;
if (-d $directory) {
	print STDERR "\tAnalyzing HiC directory: $directory\n\n";
} elsif (-f $directory) {
	print STDERR "\tTreating $directory as a Directionality Index bedGraph file...\n";
	$outFile = $directory;
	$fileFlag = 1;
} else {
	print STDERR "\tNot sure what the input is, giving up...\n";
	exit;
}


if ($fileFlag==0) {

	my $force = 1;
	my $possibleRes = HomerConfig::getHiCBgRes($directory,$superRes,$maxCPUs);


	my $rand = rand();
	my $tmpFile = $rand . ".tmp";

	my @chrs = ();
	`ls -1 "$directory"/*.tags.tsv > "$tmpFile"`;
	open IN, $tmpFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		s/\.tags\.tsv//;
		s/^.*\///;
		push(@chrs,$_);
	}
	close IN;
	`rm "$tmpFile"`;

	print STDERR "\tWill analyze chrs: @chrs\n";
	my @files = ();


	my $cpus = 0;
	foreach(@chrs) {
		my $chr = $_;

		my $chrTmpFile = $rand . ".directionBias.$chr.tmp";
		push(@files, $chrTmpFile);

		my $pid = fork();
		$cpus++;
		if ($pid == 0) {
			#child process
			print STDERR "\tAnalyzing $chr for directionality index...\n";
			`analyzeHiC "$directory" -res $res -superRes $superRes -boundary $domainDist -maxDist $domainDist -minDist $minDist $nolog $normType -nomatrix -chr $chr -std $stdFilter -min $minFilter > "$chrTmpFile" 2> /dev/null`;

			exit(0);
	
		}

		if ($cpus >= $maxCPUs) {
			my $id = wait();
			$cpus--;
		}
	}
	my $id = 0;
	while ($id >= 0) {
		$id = wait();
	}

	my %x = ();
	my $data = \%x;
	
	foreach(@files) {
		my $file = $_;
		readBedGraph($file,$data,'v');
		`rm "$file"`;
	}
	open OUT, ">$outFile";
	print OUT "track name=\"$directory Hi-C TAD Directionality Index\" type=bedGraph visibility=full\n";
	my $fh = *OUT;
	printBedGraph($fh,$data,'v');
	close $fh;
}


my %z = ();
my $data = \%z;
readBedGraph($outFile, $data, 'v');
my ($avg,$std) = getVariation($data,'v');
print STDERR "\tBasic Directionality Index Stats:\n\t\tAvg:$avg\n\t\tSTD:$std\n";
$minBias = $std*$minBias;
$minExpBias = $std*$minExpBias;
$minDelta = $std*$minDelta;
print STDERR "\t\tMin directionality index for boundary set to $minBias (expanded set to $minExpBias)\n";

my $outBed = $prefix . ".innerDomains.bed";
open OUT, ">$outBed";
my $fh = *OUT;
print $fh "track name=\"$directory Hi-C Topological Inner Domains\" type=bed visibility=dense\n";

my $outBed2 = $prefix . ".domains.bed";
open OUTEXP, ">$outBed2";
my $fhexp = *OUTEXP;
print $fhexp "track name=\"$directory Hi-C Topological Domains\" type=bed visibility=dense\n";


getDomains($fh,$fhexp,$data,$window,$std,'v');
exit;




sub readBedGraph {
	my ($file,$data,$hash) = @_;
	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		next if (/^track/);
		my @line = split /\t/;
	
		my $chr = $line[0];
		my $pos = floor(($line[1]+$line[2])/2);
		my $v = $line[3];
		if (!exists($data->{$chr})) {
			my %a = ();
			$data->{$chr} = \%a;
		}
		$data->{$chr}->{$pos}->{$hash} = $v;
	}
	close IN;
}


sub printBedGraph {
	my ($fh,$data, $hash) = @_;
	foreach(keys %$data) {
		my $c = $_;
		my @pos = sort {$a <=> $b} keys %{$data->{$c}};
		for (my $i=1;$i<@pos-1;$i++) {
			my $v = $data->{$c}->{$pos[$i]}->{$hash};
			my $s = floor(($pos[$i]+$pos[$i-1])/2);
			my $e = floor(($pos[$i+1]+$pos[$i])/2);
			print $fh "$c\t$s\t$e\t$v\n";
		}
	}
}


sub getDomains {
	my ($fh,$fhexp,$data, $window,$std, $hash) = @_;
	print STDERR "\tDomain Window: $window\n";


	#first, create moving average and store to hash 'w' (get rid of noise)
	#also store the max and min positions within the average to define boundries later
	foreach(keys %$data) {
		my $c = $_;
		my @pos = sort {$a <=> $b} keys %{$data->{$c}};
		for (my $i=0;$i<@pos;$i++){ 
			my $p = $pos[$i];
			my $v = $data->{$c}->{$pos[$i]}->{$hash};
			my $avg= $v;
			my $max= $v;
			my $min= $v;
			my $maxIndex =$i;
			my $minIndex =$i;
			my $N = 1;
			for (my $j=$i-1;$j>=0;$j--) {
				last if ($p-$pos[$j]>$window);
				$N++;
				my $v = $data->{$c}->{$pos[$j]}->{$hash};
				$avg += $v;
				if ($v > $max) {
					$max = $v;
					$maxIndex = $j;
				}
				if ($v < $min) {
					$min = $v;
					$minIndex = $j;
				}
			}
			for (my $j=$i+1;$j<@pos;$j++) {
				last if ($pos[$j]-$p>$window);
				$N++;
				my $v = $data->{$c}->{$pos[$j]}->{$hash};
				$avg += $v;
				if ($v > $max) {
					$max = $v;
					$maxIndex = $j;
				}
				if ($v < $min) {
					$min = $v;
					$minIndex = $j;
				}
			}
			my $pexp = $p;
			$v = $data->{$c}->{$pos[$i]}->{$hash};
			if ($v > 0) {
				for (my $j=$i;$j>=0;$j--) {
					if ($data->{$c}->{$pos[$j]}->{$hash} < $minExpBias) {
						last;
					}
					$pexp = $pos[$j];
				}
			} else {
				for (my $j=$i;$j<@pos;$j++) {
					if ($data->{$c}->{$pos[$j]}->{$hash} > -1*$minExpBias) {
						last;
					}
					$pexp = $pos[$j];
				}
			}
			$avg /= $N if ($N>0);
			$data->{$c}->{$p}->{'w'} = $avg;
			$data->{$c}->{$p}->{'wmax'} = $pos[$maxIndex];
			$data->{$c}->{$p}->{'wmin'} = $pos[$minIndex];
			$data->{$c}->{$p}->{'pexp'} = $pexp;
			
		}

	}

	open OUT4, ">$prefix.smoothedIndex.bedGraph";
	print OUT4 "track name=\"$directory Hi-C Directionality Index ($window smoothed)\" type=bedGraph visibility=full\n";
	my $fh4 = *OUT4;
	printBedGraph($fh4,$data,'w');
	close $fh4;


	open OUT3, ">$prefix.boundaries.txt";
	my $zz = 1;

	#Now we identify likely maxima and minima in the smoothed directional bias data
	$hash = 'w';
	my $index = 1;
	foreach(keys %$data) {
		my $c = $_;
		my @pos = sort {$a <=> $b} keys %{$data->{$c}};
		my $curSum = $data->{$c}->{$pos[0]}->{$hash};
		my $curN = 0;
		my $lowIndex =0;
		my $highIndex =-1;
		for (my $i=0;$i<@pos;$i++) {
			my $p = $pos[$i];

			#get sum across the window (not used right now...)
			#while ($lowIndex < @pos && $pos[$lowIndex] < $p-$window) {
			#	$curSum -= $data->{$c}->{$pos[$lowIndex]}->{$hash};
			#	$curN--;
			#	$lowIndex++;
			#}
			#while ($highIndex < @pos-1 &&  $pos[$highIndex+1] < $p+$window) {
			#	$highIndex++;
			#	$curSum += $data->{$c}->{$pos[$highIndex]}->{$hash};
			#	$curN++;
			#}
			#check if average directional bias is high enough to be considered
			#my $avg = $curSum / $curN;
			#if (abs($avg) > $minBias) {
			
			
			#check if directional bias is high enough to be considered
			my $v = $data->{$c}->{$pos[$i]}->{$hash};
			if (abs($v) > $minBias) {
			} else {
				next;
			}
	

			#check if we are located at a local min or max
			my $max=1;	
			my $min=1;
			my $bail = 0;
			for (my $j=$i-1;$j>=0;$j--) {
				if ($pos[$i]-$pos[$j] > $window) {
					$bail = 1;
				}
				last if ($bail && $data->{$c}->{$pos[$j]}->{$hash}*$v <= 0);
				$max = 0 if ($data->{$c}->{$pos[$j]}->{$hash} > $v);
				$min = 0 if ($data->{$c}->{$pos[$j]}->{$hash} < $v);
				last if ($max==0 && $min==0);
			}
			$bail = 0;
			for (my $j=$i+1;$j<@pos;$j++) {
				if ($pos[$j]-$pos[$i] > $window) {
					$bail = 1;
				}
				last if ($bail && $data->{$c}->{$pos[$j]}->{$hash}*$v <= 0);
				$max = 0 if ($data->{$c}->{$pos[$j]}->{$hash} > $v);
				$min = 0 if ($data->{$c}->{$pos[$j]}->{$hash} < $v);
				last if ($max==0 && $min==0);
			}
			next if ($max==0 && $min==0);
			next if ($i==0 || $i == @pos-1);

			my $hres = floor($res/2);
			if ($max) {
				$data->{$c}->{$pos[$i]}->{'max'} = $v;
				my $p = $data->{$c}->{$pos[$i]}->{'wmax'};
				my $s = $p-$hres;
				my $e = $p+$hres;
				print OUT3 "Boundary$zz\t$c\t$s\t$e\t+\t$v\n";
				$zz++;
			}
			if ($min) {
				$data->{$c}->{$pos[$i]}->{'min'} = $v;
				my $p = $data->{$c}->{$pos[$i]}->{'wmin'};
				my $s = $p-$hres;
				my $e = $p+$hres;
				print OUT3 "Boundary$zz\t$c\t$s\t$e\t-\t$v\n";
				$zz++;
			}

		}



		# Link identifed boundaries into domains.  If two 'max'(i.e. starts) boundaries occur
		# before a min(i.e. end), then stop the first domain at it's lowest point before
		# the second start.
		#
		# Also, domains are 'checked' to make sure they are reasonable.
		#
		my $curStart=-1;
		my $curStartIndex=0;
		my $curMin=0;
		my $curMinPos = -1;
		my $curMinIndex=0;
		my $curMax = 0;
		my $curMaxPos = -1;
		my $curMaxIndex=0;
		my $last='';
		for (my $i=0;$i<@pos-1;$i++) {
			my $v = $data->{$c}->{$pos[$i]}->{$hash};
			if (exists($data->{$c}->{$pos[$i]}->{'max'})) {
				if ($last ne 'min') {
					if ($curMinPos > 0 && $curStart > 0) {
						my ($ok,$err) = checkDomain($data,$c,$curStartIndex,$curMinIndex, \@pos,$hash);
						if ($ok) {
							printDomainBED($data,$fh,$c,$curStart,$curMinPos,"Domain-$index-$err","");
							printDomainBED($data,$fhexp,$c,$curStart,$curMinPos,"Domain-$index-$err","exp");
							$index++;
						}
					}
				}
				$curStart = $pos[$i];
				$curStartIndex = $i;
				$curMin = $v;
				$curMinPos = $pos[$i];
				$curMinIndex = $i;
				$curMax = $v;
				$curMaxPos = $pos[$i];
				$curMaxIndex = $i;
				$last = 'max';
			} elsif (exists($data->{$c}->{$pos[$i]}->{'min'})) {
				if ($last ne 'max') {
					$curStart = $curMaxPos;
					$curStartIndex = $curMaxIndex;
				}
				if ($curStart > 0) {
					my ($ok,$err) = checkDomain($data,$c,$curStartIndex,$i, \@pos,$hash);
					if ($ok) {
						printDomainBED($data,$fh,$c,$curStart,$pos[$i],"Domain-$index-$err","");
						printDomainBED($data,$fhexp,$c,$curStart,$pos[$i],"Domain-$index-$err","exp");
						$index++;
					}
				}
				$curStart = $pos[$i];
				$curStartIndex = $i;
				$curMax = $v;
				$curMaxPos = $pos[$i];
				$curMaxIndex = $i;
				$curMin = $v;
				$curMinPos = $pos[$i];
				$curMinIndex = $i;
				$last = 'min';
			} else {
				if ($v < $curMin) {
					$curMin = $v;
					$curMinPos = $pos[$i];
					$curMinIndex = $i;
				}
				if ($v > $curMax) {
					$curMax = $v;
					$curMaxPos = $pos[$i];
					$curMaxIndex = $i;
				}
			}
		}
	}
	print STDERR "\tTotal Domains: $index\n";
	close OUT3;
}
sub checkDomain {
	my ($data, $c, $sIndex, $eIndex, $pos,$hash) = @_;
	#return (1,0);
	if ($sIndex == $eIndex) {
		#print STDERR "Indexes are the same!!\n";
		return (0,0);
	}
	my $StartV = $data->{$c}->{$pos->[$sIndex]}->{$hash};
	my $EndV = $data->{$c}->{$pos->[$eIndex]}->{$hash};
	my $delta = $StartV - $EndV;
	my $length = $pos->[$eIndex]-$pos->[$sIndex];
	if ($length < $minDomainLength) {
		#print STDERR "Window too small ($length)!!\n";
		return (0,0);
	}
	my $error = 0;
	my $N = 0;
	for (my $i=$sIndex+1;$i<$eIndex;$i++) {
		my $v = $data->{$c}->{$pos->[$i]}->{$hash};
		my $dist = $pos->[$i]-$pos->[$sIndex];
		my $expected = $StartV-$delta*($dist/$length);
		#print STDERR "$c\t$pos->[$sIndex]\t$pos->[$eIndex]\t$i\t$dist\t$length\t$expected\t$v\t$delta\n";
		$error += abs($expected-$v);
		$N++;
	}
	
	$error /= $N if ($N>0);
	$error /= $delta if ($delta > 0);
	my $ok = 1;
	if ($delta < $minDelta) {
		$ok = 0;
	}
	if ($error > $maxError) {
		$ok = 0;
	}
	return ($ok,$error);
}

sub printDomainBED {
	my ($data,$fh, $c, $start, $end, $name, $type) = @_;
	my $s = $data->{$c}->{$start}->{'wmax'};
	my $e = $data->{$c}->{$end}->{'wmin'};
	if ($type eq 'exp') {
		$s = $data->{$c}->{$start}->{'pexp'};
		$e = $data->{$c}->{$end}->{'pexp'};
	}
	return if ($e-$s < 100);
	my $ss = $e-$s-50;
	print $fh "$c\t$s\t$e\t$name\t1\t+\t$s\t$e\t150,150,150\t2\t50,50\t0,$ss\n";
}
sub getVariation {
	my ($data, $hash) = @_;
	my $avg = 0;
	my $N = 0;
	foreach(keys %$data) {
		my $c = $_;
		my @pos = sort {$a <=> $b} keys %{$data->{$c}};
		for (my $i=1;$i<@pos-1;$i++) {
			my $v = $data->{$c}->{$pos[$i]}->{$hash};
			$avg += $v;
			$N++;
		}
	}
	$avg /= $N if ($N > 0);	
	my $std=0;

	foreach(keys %$data) {
		my $c = $_;
		my @pos = sort {$a <=> $b} keys %{$data->{$c}};
		for (my $i=1;$i<@pos-1;$i++) {
			my $v = $data->{$c}->{$pos[$i]}->{$hash};
			$std += ($avg-$v)*($avg-$v);
		}
	}
	$std /= $N if ($N > 0);
	$std = sqrt($std);
	return ($avg, $std);
}

sub calcDerivative {
	my ($data,$window,$hash,$hash2) = @_;
	foreach(keys %$data) {
		my $c = $_;
		my @pos = sort {$a <=> $b} keys %{$data->{$c}};
		for (my $i=0;$i<@pos;$i++) {
			my $sumUp=0;
			my $sumDown=0;
			my $nUp=0;
			my $nDown=0;
			for (my $j=$i-1;$j>=0;$j--) {
				last if ($pos[$i]-$pos[$j] > $window);
				$sumUp += $data->{$c}->{$pos[$j]}->{$hash};
				$nUp++;
			}
			for (my $j=$i+1;$j<@pos;$j++) {
				last if ($pos[$j]-$pos[$i] > $window);
				$sumDown += $data->{$c}->{$pos[$j]}->{$hash};
				$nDown++;
			}
			$sumUp /= $nUp if ($nUp>0);
			$sumDown /= $nDown if ($nDown>0);
	
			$data->{$c}->{$pos[$i]}->{$hash2} = $sumDown-$sumUp;
		}
	}
}

