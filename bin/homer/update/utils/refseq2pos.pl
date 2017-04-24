#!/usr/bin/perl -w


my $most5pFlag = 1;
if (@ARGV < 1) {
	print STDERR "<refseqAli.txt> [1 for TTS]\n";
	print STDERR "<refseqAli.txt> [2 for for full gene]\n";
	print STDERR "most 5' flag set to $most5pFlag\n";
	exit;
}
my $ttsFlag = 0;
if (@ARGV > 1) {
	$ttsFlag = $ARGV[1];
}
my %hits = ();
open IN, $ARGV[0];
while (<IN>){ 
	chomp;
	my @line = split /\t/;
	my $id = $line[10];
	my $chr= $line[14];
	my $tss = $line[16]+1;
	my $dir = 0;
	if ($ttsFlag == 0) {
		if ($line[9] eq '+') {
			$dir = 0;
			$tss = $line[16]+1;
		} else {
			$dir = 1;
			$tss = $line[17];
		}
	} else {
		if ($line[9] eq '+') {
			$dir = 0;
			$tss = $line[17];
		} else {
			$dir = 1;
			$tss = $line[16]+1;
		}
	}
	my $start = $tss-2000;
	my $end = $tss+2000;
	if ($ttsFlag == 2) {
		$start = $line[16]+1;
		$end = $line[17];
	}
	my $str =  "$id\t$chr\t$start\t$end\t$dir\n";

	my $score = $line[1];
	if ($most5pFlag) {
		$score = -1*$tss;
		if ($dir == 1) {
			$score = $tss;
		}
	}

	if (!exists($hits{$id})) {
		$hits{$id}={d=>$str,s=>$score,c=>$chr};
	} else {
		if ($score > $hits{$id}->{'s'}) {
			$hits{$id}={d=>$str,s=>$line[1],c=>$chr};
		} elsif ($score == $hits{$id}->{'s'}) {
			if ($hits{$id}->{'c'} =~ /_/) {
				$hits{$id}={d=>$str,s=>$line[1],c=>$chr};
			}
		}
	}
}
close IN;

foreach(values %hits) {
	print $_->{'d'};
}
if ($ttsFlag == 0) {
	print STDERR "Returned TSS\n";
} else {
	print STDERR "Returned TTS - i.e. END of GENE\n";
}
