#!/usr/bin/env perl
use warnings;


sub printCMD() {
	print STDERR "\n\tUsage: findHiCCompartments.pl <PC1.txt file> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-opp (return inactive, not active regions)\n";
	print STDERR "\t\t-thresh <#> (threshold for active regions, default: 0)\n";
	print STDERR "\t\t-bg <2nd PC1.txt file> (for differential domains)\n";
	print STDERR "\t\t\t-diff <#> (difference threshold, default: 50)\n";
	print STDERR "\t\t-corr <corrDiff.txt file> (for differential domains)\n";
	print STDERR "\t\t\t-corrDiff <#> (correlation threshold, default: 0.4)\n";
	print STDERR "\t\t-peaks (output as peaks/regions, not continuous domains)\n";
	print STDERR "\n";
	exit;
}
if (@ARGV < 1) {
	printCMD();
}
my $oppFlag = 0;
my $threshold = 0;
my $eigenFile2 = "";
my $bgDiffThresh = 50;
my $corrDiffFile = "";
my $corrDiffThresh = 0.4;
my $peakFlag = 0;

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-opp') {
		$oppFlag = 1;
		print STDERR "\tFinding opposite (inactive regions)\n";
	} elsif ($ARGV[$i] eq '-thresh') {
		$threshold = $ARGV[++$i];
		print STDERR "\tActive threshold set at $threshold\n";
	} elsif ($ARGV[$i] eq '-bg') {
		$eigenFile2 = $ARGV[++$i];
		print STDERR "\tBackground regions defined by $eigenFile2\n";
	} elsif ($ARGV[$i] eq '-corr') {
		$corrDiffFile = $ARGV[++$i];
		print STDERR "\tBackground region correlation defined by $corrDiffFile\n";
	} elsif ($ARGV[$i] eq '-peaks') {
		$peakFlag =1;
	} elsif ($ARGV[$i] eq '-diff') {
		$bgDiffThresh = $ARGV[++$i];
		print STDERR "\tBackground differential set at $bgDiffThresh\n";
	} elsif ($ARGV[$i] eq '-corrDiff') {
		$corrDiffThresh = $ARGV[++$i];
		print STDERR "\tBackground correlation differential set at $corrDiffThresh\n";
	} else {
		printCMD();
	}
}

my $fThresh = 1;
my %data = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^#/);
	my @line = split /\t/;
	if (@line < 5) {
		print STDERR "\t!!! Is this file a PC1.txt formatted peak file? (bedGraph won't work...)\n";
		next;
	}
	my $name =$line[0];
	my $chr =$line[1];
	my $start = $line[2];
	my $end = $line[3];
	my $dir = $line[4];
	my $v = $line[5];
	my $good = 0;
	if ($v > $threshold && $oppFlag == 0) {
		$good = 1;
	}
	if ($v < $threshold && $oppFlag == 1) {
		$good = 1;
	}
	my $p = {name=>$name,c=>$chr,s=>$start,e=>$end,d=>$dir,v=>$v,ok=>$good,f=>1};
	if (!exists($data{$chr})) {
		my %a = ();
		$data{$chr} = \%a;
	}
	$data{$chr}->{$name} = $p;
}
close IN;

if ($eigenFile2 ne '') {

	$fThresh++;
	open IN, $eigenFile2;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		my @line = split /\t/;
		if (@line < 5) {
			print STDERR "\t!!! Is this file a PC1.txt formatted peak file? (bedGraph won't work...)\n";
			next;
		}
		my $name =$line[0];
		my $chr =$line[1];
		my $start = $line[2];
		my $end = $line[3];
		my $dir = $line[4];
		my $v = $line[5];
		next if (!exists($data{$chr}));
		next if (!exists($data{$chr}->{$name}));
		$data{$chr}->{$name}->{'v2'} = $v;
		if (($data{$chr}->{$name}->{'v'} - $v > $bgDiffThresh && $oppFlag == 0) 
		 			|| ($v - $data{$chr}->{$name}->{'v'}  > $bgDiffThresh && $oppFlag == 1) ) {
			if ($data{$chr}->{$name}->{'ok'} > 0) {
				$data{$chr}->{$name}->{'ok'} = 1;
			}
		} else {
			$data{$chr}->{$name}->{'ok'} = 0;
		}
		$data{$chr}->{$name}->{'f'}++;
	}
	close IN;
}
if ($corrDiffFile ne '') {
	$fThresh++;
	open IN, $corrDiffFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		next if (/^#/);
		my @line = split /\t/;
		my $name =$line[0];
		my $chr =$line[1];
		my $start = $line[2];
		my $end = $line[3];
		my $dir = $line[4];
		my $v = $line[5];
		next if (!exists($data{$chr}));
		next if (!exists($data{$chr}->{$name}));
		$data{$chr}->{$name}->{'cd'} = $v;
		if ($v < $corrDiffThresh) {
			if ($data{$chr}->{$name}->{'ok'} > 0) {
				$data{$chr}->{$name}->{'ok'} = 1;
			}
		} else {
			$data{$chr}->{$name}->{'ok'} = 0;
		}
		$data{$chr}->{$name}->{'f'}++;
	}
	close IN;
}
if ($peakFlag) {
	print STDERR "\tPrinting peaks (not domains)...\n";
	foreach(keys %data) {
		my $chr = $_;
		my @names = sort {$data{$chr}->{$a}->{'s'} <=> $data{$chr}->{$b}->{'s'}} keys %{$data{$chr}};
		for (my $i=0;$i<@names;$i++) {
			if ($data{$chr}->{$names[$i]}->{'ok'}==1 && $data{$chr}->{$names[$i]}->{'f'}>= $fThresh) {
				my $curStart = $data{$chr}->{$names[$i]}->{'s'};
				my $curEnd = $data{$chr}->{$names[$i]}->{'e'};
				my $avg = $data{$chr}->{$names[$i]}->{'v'};
				my $len = $curEnd-$curStart;
				print "$names[$i]\t$chr\t$curStart\t$curEnd\t+\t$avg\t$len\n";
			}
		}
	}
	exit;
}

my $id = 0;
foreach(keys %data) {
	my $chr = $_;
	my @names = sort {$data{$chr}->{$a}->{'s'} <=> $data{$chr}->{$b}->{'s'}} keys %{$data{$chr}};
	
	my $curStart = -1;
	my $curEnd = -1;
	my $curSig = 0;
	my $curN = 0;
	for (my $i=0;$i<@names;$i++) {
		if ($data{$chr}->{$names[$i]}->{'ok'}==0 && $data{$chr}->{$names[$i]}->{'f'}>= $fThresh) {
			if ($curStart >= 0) {
				$id++;
				my $avg = $curSig/$curN;
				my $len= $curEnd-$curStart;
				print "Domain-$id\t$chr\t$curStart\t$curEnd\t+\t$avg\t$len\n";
			}
			$curStart = -1;
			$curEnd = -1;
			$curSig = 0;
			$curN = 0;
		} elsif ($data{$chr}->{$names[$i]}->{'ok'}==1 && $data{$chr}->{$names[$i]}->{'f'}>= $fThresh) {
			if ($curStart < 0) {
				$curStart = $data{$chr}->{$names[$i]}->{'s'};
			}
			$curEnd = $data{$chr}->{$names[$i]}->{'e'};
			$curSig += $data{$chr}->{$names[$i]}->{'v'};
			$curN++;
		}

	}
	if ($curStart >= 0) {
		$id++;
		my $avg = $curSig/$curN;
		my $len= $curEnd-$curStart;
		print "Domain-$id\t$chr\t$curStart\t$curEnd\t+\t$avg\t$len\n";
	}
}

