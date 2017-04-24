#!/usr/bin/perl -w

if (@ARGV < 1) {
	print STDERR "<sam file>\n";
	exit;
}

my $minSpliceLength = 20;

my $totalReads = 0;
my $spliceReads = 0;

my %data = ();
for (my $z=0;$z<@ARGV;$z++) {

open IN, $ARGV[$z];
while (<IN>) {
	chomp;
	next if (/^\@/);
	my @line = split /\t/;

	my $strand = 0;
	if ($line[1] & 16) {
		$strand = 1;
	}
	my $chr = $line[2];
	my $pos = $line[3];
	my $mapq = $line[4];
	my $cigar = $line[5];

	$totalReads++;
	my $sflag = 0;
	my $curPos = $pos;
	while ($cigar =~ s/(\d+)([A-Z])//) {
		my $len = $1;
		my $code = $2;
		if ($code eq 'M') {
			$curPos += $len;
		}
		if ($code eq 'N') {
			my $s = $curPos;
			$curPos += $len;
			my $e = $curPos-1;

			if ($len > $minSpliceLength) {
				$sflag=1;
				
				my $id = "$chr:$s-$e:$strand";
				if (!exists($data{$id})) {
					$data{$id} = {c=>$chr,s=>$s,e=>$e,d=>$strand,v=>1};
				} else {
					$data{$id}->{'v'}++;
				}
			}
		}
	}
	$spliceReads += $sflag;
}
close IN;
}

print STDERR "\tTotalReads= $totalReads\n";
print STDERR "\tSpliceReads=$spliceReads\n";
my $p = $spliceReads/$totalReads;
print STDERR "\t$p\n";

print "track name=$ARGV[0] description=\"$ARGV[0]\" useScore=1\n";
my $ID = 1;
foreach(values %data) {
	my $c = $_->{'c'};
	my $s = $_->{'s'}-1;
	my $e = $_->{'e'};
	my $d = $_->{'d'};
	my $v = $_->{'v'};
	my $ss = $s-1;
	my $ee = $e+1;
	my $strand = "+";
	if ($d eq '1') {
		$strand = "-";
	}
	my $l = $e-$s+1;
	print "$c\t$ss\t$ee\tSJ$ID-$v\t$v\t$strand\t$ss\t$ee\t0\t2\t1,1,\t0,$l,\n";
	$ID++;
}
