#!/usr/bin/env perl
use warnings;


my $pseudoCount = 3;
if (@ARGV < 2) {
	print STDERR "<bedgraph file1> <bedgraph file2> [method]\n";
	print STDERR "method = diff or ratio\n";
	print STDERR "Ratio are calculated with a pseudo count of $pseudoCount\n";
	print STDERR "Can both be gzipped\n";
	exit;
}
my $method = "ratio";
if (@ARGV > 2) {
	$method = $ARGV[2];
	print STDERR "Method set to $method\n";
}
my $log2 = log(2);
my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
my $file1Zip = 0;
my $file2Zip = 0;
if ($file1 =~ /\.gz$/) {
	$file1Zip = 1;
	`gunzip "$file1"`;
	$file1 =~ s/\.gz//;
}
if ($file2 =~ /\.gz$/) {
	$file2Zip = 1;
	`gunzip "$file2"`;
	$file2 =~ s/\.gz//;
}
my $d1 = readBedGraph($file1);
my $d2 = readBedGraph($file2);

print "track type=bedGraph name=\"Difference ($method) between $file1 and $file2\" visibility=full yLineOnOff=on autoScale=on yLineMark=\"0.0\" color=200,0,0 maxHeightPixels=128:100:11\n";

foreach(keys %$d1) {
	my $chr = $_;
	my $i1=0;
	my $i2=0;
	my @pos1 = sort {$a <=> $b} keys %{$d1->{$chr}};
	my @pos2 = sort {$a <=> $b} keys %{$d2->{$chr}};
	my %pos = ();
	foreach(@pos1) {
		$pos{$_}=1;
	}
	foreach(@pos2) {
		$pos{$_}=1;
	}	
	my @pos = sort {$a <=> $b} keys %pos;
	my $v1 = 0;
	my $v2 = 0;
	my $V = 0;
	my $P1 = $pos[0];
	my $lastp = -1;
	foreach(@pos) {
		my $p = $_;
		if ($lastp != -1) {
			my $V=0;
			if ($method eq 'diff') {
				$V = $v2 - $v1;
			} elsif ($method eq 'ratio') {
				$V = log(($v2+$pseudoCount)/($v1+$pseudoCount))/$log2;
			}
			$V = sprintf("%.2f", $V);
			print "$chr\t$lastp\t$p\t$V\n";
		}

		if (exists($d1->{$chr}->{$p})) {
			$v1 = $d1->{$chr}->{$p};
		}
		if (exists($d2->{$chr}->{$p})) {
			$v2 = $d2->{$chr}->{$p};
		}
		$lastp = $p;

	}
}
	
	



if ($file1Zip) {
	`gzip "$file1"`;
}
if ($file2Zip) {
	`gzip "$file2"`;
}


sub readBedGraph {
	my ($file) = @_;
	my %data = ();
	open IN, $file;
	while (<IN>) {
		chomp;
		next if (/^track/);
		my @line=split /\s+/;
		my $chr = $line[0];
		if (!exists($data{$chr})) {
			my %a = ();
			$data{$chr} = \%a;
		}
		my $p1 = $line[1];
		my $p2 = $line[2];
		my $v = $line[3];
		$data{$chr}->{$p1}=$v;
	}
	close IN;
	return \%data;
}
