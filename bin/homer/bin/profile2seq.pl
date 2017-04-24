#!/usr/bin/env perl
use warnings;

if (@ARGV < 2) {
	print STDERR "<profile> <# seq> [-rna]\n";
	exit;
}

my $rnaFlag = 0;
if (@ARGV > 2 && $ARGV[2] eq '-rna') {
	$rnaFlag=1;
}
my $type = 0;
my $header = "";
my @profile = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	if (/^>/) {
		$type = 0;
		$header = $_;
		next;
	}
	if (/^</) {
		$type = 1;
		$header = $_;
		next;
	}
	my @line =split /\t/;
	next if (@line < 4);
	my $sum =0;
	foreach(@line) {
		$sum += $_;
	}
	foreach(@line) {
		$_ /= $sum;
	}
	push(@profile, \@line);
}
close IN;

my @dna =('A','C','G','T');
if ($rnaFlag) {
	@dna =('A','C','G','U');
}
my @prot = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');


my $N = $ARGV[1];
for (my $i=0;$i<$N;$i++) {
	print ">$i\n";
	for (my $j=0;$j<@profile;$j++) {

		my $r = $i/$N;
		my $sum = 0;
		if ($type == 0) {
			for (my $k=0;$k<@dna;$k++) {
				if ($k<@{$profile[$j]}) {
					$sum += $profile[$j][$k];
					if ($r < $sum) {
						print "$dna[$k]";
						last;
					}
				} else {
					print STDERR "Formatting problem with $ARGV[0]\n$header\n";
					sleep(5);
				}
			}
		} elsif ($type==1) {
			for (my $k=0;$k<@prot;$k++) {
				$sum += $profile[$j][$k];
				if ($r < $sum) {
					print "$prot[$k]";
					last;
				}
			}
		}
	}
	print "\n";
}

	
