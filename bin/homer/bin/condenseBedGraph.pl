#!/usr/bin/env perl
use warnings;




if (@ARGV < 1) {
	print STDERR "\n\tcondenseBedGraph.pl <bedGraph file> [2nd bedGraph file] ...\n";
	print STDERR "\t\t-s <bedGraph file to subtract>\n";
	print STDERR "\t\t-log (output log2)\n";
	print STDERR "\n";
	exit;
}

my $mode = 'sum';
my $trackLine = 0;
my %data = ();
my $logFlag = 0;
my $coeff = 1;

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-s') {
		$coeff = -1;
		next;
	}
	if ($ARGV[$i] eq '-log') {
		$logFlag = 1;
		print STDERR "\tWill output log2 of signal\n";
		next;
	}
	open IN, $ARGV[$i];
	while (<IN>) {
		chomp;
		s/\r//g;
		if (/^track/) {
			$trackLine = $_;
			next;
		}
		my @line = split /\t/;
		next if (@line < 4);
		my $chr = $line[0];
		my $start = $line[1];
		my $end = $line[2];
		my $v = $line[3]*$coeff;
		if (!exists($data{$chr})) {
			my %a = ();
			$data{$chr} = \%a;
		}
		$data{$chr}->{$start}=0  if (!exists($data{$chr}->{$start}));
		$data{$chr}->{$end}=0  if (!exists($data{$chr}->{$end}));
		
		if ($mode eq 'sum') {
			$data{$chr}->{$start}+=$v;
			$data{$chr}->{$end}-=$v;
		}

	}
	close IN;
}


print "$trackLine\n";
foreach(keys %data) {
	my $c = $_;
	my $v = 0;
	my @pos = sort {$a <=> $b} keys %{$data{$c}};
	for (my $i=0;$i<@pos-1;$i++) {
		$v += $data{$c}->{$pos[$i]};
		print "$c\t$pos[$i]\t$pos[$i+1]\t$v\n";
	}
}
