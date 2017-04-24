#!/usr/bin/perl -w

$atgFlag = 1;

my %hits = ();
open IN, $ARGV[0];
while (<IN>){ 
	chomp;
	my @line = split /\t/;
	my $id = $line[1];
	my $chr= $line[2];
	my $tss = $line[4]+1;
	my $atg = $line[6]+1;
	my $dir = 0;
	my $score = $line[0];

	if ($line[3] eq '+') {
	} else {
		$dir = 1;
		$tss = $line[5];
		$atg = $line[7];
	}
	my $start = $tss-500;
	my $end = $tss+500;
	if ($atgFlag==1) { 
		$start = $atg-500;
		$end = $atg+500;
	}

	my $str =  "$id\t$chr\t$start\t$end\t$dir\n";
	if (!exists($hits{$id})) {
		$hits{$id}={d=>$str,s=>$score,c=>$chr};
	} else {
		if ($score > $hits{$id}->{'s'}) {
			$hits{$id}={d=>$str,s=>$score,c=>$chr};
		} elsif ($score == $hits{$id}->{'s'}) {
			if ($hits{$id}->{'c'} =~ /_/) {
				$hits{$id}={d=>$str,s=>$score,c=>$chr};
			}
		}
	}
}
close IN;

foreach(values %hits) {
	print $_->{'d'};
}
