#!/usr/bin/env perl
#
#
if (@ARGV < 3) {
	print STDERR "<conservation file (seq-like 0-9)> <site file> <offset>\n";
	exit;
}
my $offset = $ARGV[2];
my $maxRatio = 10;
my $window = 1;

my %sites = ();	
open IN, $ARGV[1];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$site = {offset=>$line[1], seq=>$line[2],cons=>$line[3], dir=>$line[4], name=>$line[5],score=>$line[6]};
	if (!exists($sites{$line[0]})) {
		my @a = ();
		$sites{$line[0]} = \@a;
	}
	push(@{$sites{$line[0]}}, $site);
}
close IN;

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (!exists($sites{$line[0]}));
	my @cons = split undef, $line[1];
	foreach(@{$sites{$line[0]}}) {
		my $L = length($_->{'seq'});		
		my $motifCons = 0;
		my $start = $_->{'offset'} - $offset;
		for (my $i=$start;$i<$start+$L;$i++) {
			$motifCons += $cons[$i];
		}
		$motifCons /= $L;

		my $winCons = 0;
		my $n=0;
		for (my $i=$start-$L*$window;$i<$start;$i++) {
			next if ($i<0);
			$winCons += $cons[$i];
			$n++;
		}
		for (my $i=$start+$L;$i<$start+$L+$L*$window;$i++) {
			next if ($i>=@cons);
			$winCons += $cons[$i];
			$n++;
		}
		$winCons /= $n if ($n != 0);
		print "$line[0]\t$_->{'offset'}\t$_->{'seq'}\t$motifCons\t$_->{'dir'}\t$_->{'name'}\t$_->{'score'}\n";
	}
}
close IN;
