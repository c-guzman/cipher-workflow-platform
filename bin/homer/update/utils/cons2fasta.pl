#!/usr/bin/perl -w
use POSIX;
if (@ARGV < 1) {
	print STDERR "<conservation.data file> - makes a FASTA like file\n";
	exit;
}
$maxlen = 100;
$currentlen = 0;

my $index = 0;

print ">$ARGV[0]-conservation\n";

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	if (/^f/){ 
		/fixedStep .* start=(\d+?) step/;
		my $start = $1;
		my $str = '';
		my $L = $start - $index-1;
#print STDERR "$L - adding start=$start, index=$index\n";
		for (my $i=0;$i<$L;$i++) {
			$str .= '0';
		}
		$index += $L;
		printdata($str);
	} else {
		my $v = floor(10*$_);
		$v = 9 if ($v > 9);
		$v = 0 if ($v < 0);
		printdata($v);
		$index++;
	}
}
close IN;

sub printdata {
	my ($data) = @_;

	my $total = length($data);
	my $spaceleft = $maxlen - $currentlen;
	if ($total <= $spaceleft) {
		print "$data";
		$currentlen += $total;
		if ($currentlen >= $maxlen) {
			print "\n";
			$currentlen = 0;
		}
	} else {
		my $s = substr($data, 0, $spaceleft);
		print "$s\n";
		my $dataindex += $spaceleft;
		my $lines = floor(($total-$dataindex)/$maxlen);
		for (my $i=0;$i<$lines;$i++) {
			$s = substr($data, $dataindex, $maxlen);
			print "$s\n";
			$dataindex += $maxlen;
		}
		$currentlen = 0;
		if ($dataindex < $total) {
			$s = substr($data, $dataindex);
			$currentlen += $total-$dataindex;
			print "$s";
			if ($currentlen >= $maxlen) {
				print "\n";
				$currentlen = 0;
			}
		}
	}
}
