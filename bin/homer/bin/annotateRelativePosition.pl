#!/usr/bin/env perl
use warnings;


use POSIX;
if (@ARGV < 2) {
	print STDERR "\n\t<tag file 1> <tag file 2> [same direction=1 | 0]\n";
	print STDERR "\tIf position file, add \",\" and offset to center (i.e. pos.file,-2000)\n";
	print STDERR "\tIf no offset is given (i.e. \"pos.file, \"), then the midpoint is used\n\n";
	exit;
}

my $dirFlag = 0;
if (@ARGV > 2) {
	$dirFlag = $ARGV[2];
	#print STDERR "\tDirection Flag is set to $dirFlag (0 = strand unimportant, 1 = strand important)\n";
}

my %query = ();
my $posFlag = 0;
my $offset = 0;
my $midFlag = 0;
my $filename = $ARGV[0];
if ($filename =~ /^(.*?)\,(.*?)$/) {
	$filename = $1;
	$offset = $2;
	$posFlag = 1;
	if ($offset eq '') {
		$midFlag = 1;
		#print STDERR "\t$filename is position file - using midpoint as reference\n";
	} else {
		#print STDERR "\t$filename is position file with offset $offset\n";
	}
}
open IN, $filename;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line= split /\t/;
	my $p = $line[2];
	if ($p =~ /[^\d\-]/) {
		next;
	}
	my $c = $line[1];
	my $id = $line[0];
	my $d = $line[3];
	my $startOffset = 0;
	my $endOffset = 0;
	if ($posFlag) {
		$d = $line[4];
		if ($d eq '-' || $d eq '1') {
			$d = 1;
		} else {
			$d = 0;
		}
		if ($midFlag) {
			$p = floor(($line[2]+$line[3])/2);
			$endOffset = $line[3]-$p;
			$startOffset = -1*$endOffset;
		} else {
			$p = $line[2] - $offset;
			$startOffset = $line[2]-$p;
			$endOffset = $line[3]-$p;
			if ($d == 1) {
				my $tmp = $startOffset *-1;
				$startOffset = -1*$endOffset;
				$endOffset = $tmp;
				$p = $line[3] + $offset;
			}
		}
	} else {
		if ($d eq '-' || $d eq '1') {
			$d = 1;
		} else {
			$d = 0;
		}
	}
	if (!exists($query{$c})) {
		my %a = ();
		my %b = ();
		my @dir = (\%a, \%b);
		$query{$c} = \@dir;
	}
	my $entry = {p=>$p,d=>$d,s=>$startOffset,e=>$endOffset};
	if ($dirFlag == 1) {
		$query{$c}->[$d]->{$id} = $entry;
	} else {
		$query{$c}->[0]->{$id} = $entry;
	}
}
close IN;

my %target = ();
$posFlag = 0;
$offset = 0;
$midFlag = 0;
$filename = $ARGV[1];
if ($filename =~ /^(.*?)\,(.*?)$/) {
	$filename = $1;
	$offset = $2;
	$posFlag = 1;
	if ($offset eq '') {
		$midFlag = 1;
		#print STDERR "\t$filename is position file - using midpoint as reference\n";
	} else {
		#print STDERR "\t$filename is position file with offset $offset\n";
	}
}
open IN, $filename;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line= split /\t/;
	my $p = $line[2];
	if ($p =~ /[^\d\-]/) {
		next;
	}
	my $c = $line[1];
	my $id = $line[0];
	my $d = $line[3];
	if ($posFlag) {
		$d = $line[4];
		if ($d eq '-' || $d eq '1') {
			$d = 1;
		} else {
			$d = 0;
		}
		if ($midFlag) {
			$p = floor(($line[2]+$line[3])/2);
		} else {
			$p = $line[2] - $offset;
			if ($d == 1) {
				$p = $line[3] + $offset;
			}
		}
	} else {
		if ($d eq '-' || $d eq '1') {
			$d = 1;
		} else {
			$d = 0;
		}
	}
	if (!exists($target{$c})) {
		my %a = ();
		my %b = ();
		my @dir = (\%a, \%b);
		$target{$c} = \@dir;
	}
	my $entry = {p=>$p,d=>$d};
	if ($dirFlag == 1) {
		$target{$c}->[$d]->{$id} = $entry;
	} else {
		$target{$c}->[0]->{$id} = $entry;
	}
}
close IN;


foreach(keys %query) {
	my $c = $_;
	next if (!exists($target{$c}));

	for (my $d = 0;$d <2;$d ++) {
		my @query = sort {$query{$c}->[$d]->{$a}->{'p'} <=> $query{$c}->[$d]->{$b}->{'p'}} keys %{$query{$c}->[$d]};
		next if (@query < 1);
		my @target = sort {$target{$c}->[$d]->{$a}->{'p'} <=> $target{$c}->[$d]->{$b}->{'p'}} keys %{$target{$c}->[$d]};
		next if (@target < 1);
		my $NT = @target;
		my $tIndex = 0;
		for (my $qIndex = 0;$qIndex < @query;$qIndex++) {
			my $qpos = $query{$c}->[$d]->{$query[$qIndex]}->{'p'};
			my $qdir = $query{$c}->[$d]->{$query[$qIndex]}->{'d'};
#print STDERR "qdir $qdir $qIndex $qpos $c $d $query[$qIndex]\n";
			my $lowLimit = $qpos + $query{$c}->[$d]->{$query[$qIndex]}->{'s'};
			my $highLimit = $qpos + $query{$c}->[$d]->{$query[$qIndex]}->{'e'};
			if ($qdir == 1) {
				$lowLimit = $qpos - $query{$c}->[$d]->{$query[$qIndex]}->{'e'};
				$highLimit = $qpos - $query{$c}->[$d]->{$query[$qIndex]}->{'s'};
			}
			my $numPeaks = 0;
			my %peakPos = ();
			my $resetIndex = $tIndex;
			
			my $bestT = '';
			my $tpos = 0;
			my $tdir = 0;
			my $bestDist = 1e100;
			if ($tIndex >= @target) {
				$bestT = scalar(@target-1);
				$tpos= $target{$c}->[$d]->{$target[$tIndex]}->{'p'};
				$tdir= $target{$c}->[$d]->{$target[$tIndex]}->{'d'};
				$bestDist = abs($qpos-$tpos);
			} else {
				$bestT = $tIndex;
				$tpos= $target{$c}->[$d]->{$target[$tIndex]}->{'p'};
				$tdir= $target{$c}->[$d]->{$target[$tIndex]}->{'d'};
				$bestDist = abs($qpos-$tpos);
				while ($tpos < $qpos) {
					$bestDist = abs($qpos-$tpos);
					$bestT = $tIndex;
					$tIndex++;
					if ($tpos > $lowLimit) {
						$numPeaks++;
						my $poffset = $tpos-$qpos;
						$poffset*=-1 if($qdir eq '1' || $qdir eq '-');
						if (!exists($peakPos{$poffset})) {
							my @zz = (0,0);
							$peakPos{$poffset} = \@zz;
						}
						my $ddindex = 0;
						$ddindex = 1 if ($tdir != $qdir);
						$peakPos{$poffset}->[$ddindex]++;
					} else {
						$resetIndex++;
					}
					last if ($tIndex >= @target);
					$tpos= $target{$c}->[$d]->{$target[$tIndex]}->{'p'};
					$tdir= $target{$c}->[$d]->{$target[$tIndex]}->{'d'};
				}
				if ($tIndex < @target) {
					my $dist = abs($qpos-$tpos);
					if ($dist < $bestDist) {
						$bestDist = $dist;
						$bestT = $tIndex;
					}
					my $tmpIndex = $tIndex;
					while ($tmpIndex < @target) {
						$tpos = $target{$c}->[$d]->{$target[$tmpIndex]}->{'p'};
						$tdir = $target{$c}->[$d]->{$target[$tmpIndex]}->{'d'};
						if ($tpos < $highLimit) {
							$numPeaks++;
							my $poffset = $tpos-$qpos;
							$poffset*=-1 if($qdir eq '1' || $qdir eq '-');
							if (!exists($peakPos{$poffset})) {
								my @zz = (0,0);
								$peakPos{$poffset} = \@zz;
							}
							my $ddindex = 0;
							$ddindex = 1 if ($tdir != $qdir);
							$peakPos{$poffset}->[$ddindex]++;
						} else {
							last;
						}
						$tmpIndex++;
					}
				}
				$tIndex = $resetIndex;
				if ($tIndex > 0) {
					$tIndex--;
				}
			}
			$tpos = $target{$c}->[$d]->{$target[$bestT]}->{'p'};
			$tdir = $target{$c}->[$d]->{$target[$bestT]}->{'d'};

			my $dist = $tpos - $qpos;
			if ($qdir == 1) {
				$dist *= -1;
			}
			$relativeStrand = '-';
			$relativeStrand = '+' if ($qdir eq $tdir);
			my $tid = $target[$bestT];
			my $qid = $query[$qIndex];
			print "$qid\t$tid\t$dist\t$tdir\t$numPeaks\t$relativeStrand\t";
			my @pp = sort {$a <=> $b} keys %peakPos;
			my $cc = 0;
			foreach(@pp) {
				print "," if ($cc > 0);
				$cc++;
				my $v1 = $peakPos{$_}->[0];
				$v1 = 'NA' if ($v1 < 1);
				my $v2 = $peakPos{$_}->[1];
				$v2 = 'NA' if ($v2 < 1);
				print "$_=$v1|$v2";
			}

			print "\n";
		}
	}
}



