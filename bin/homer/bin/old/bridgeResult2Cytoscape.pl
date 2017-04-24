#!/usr/bin/perl -w

if (@ARGV < 2) {
	print STDERR "<result file> <prefix name>\n";
	exit;
}

my $prefix = $ARGV[1];
open NETWORK, ">$prefix.network.txt";
open NODESIZE, ">$prefix.nodes.size.txt";
open EDGERATIO, ">$prefix.edgeRatio.ratio.txt";
open EDGEPVALUE, ">$prefix.edgeRatio.pvalue.txt";
print NODESIZE "NodeSize\n";
print EDGERATIO "EdgeRatio\n";
print EDGEPVALUE "EdgePValue\n";
my %done = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
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
