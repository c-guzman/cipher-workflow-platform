#!/usr/bin/perl -w

my %reads = ();

my $percent = 0.75;
if (@ARGV < 1) {
	print STDERR "<blat *.psl file>\n";
	exit;
}
my $program = "blat";
my %best = ();

open IN, $ARGV[0];
my $count = 0;
while (<IN>) {
	$count++;
	next if ($count < 4);
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (@line < 20);
	next if ($line[0] eq 'match');
	next if ($line[0] eq '');
	my $id = $line[9];
	my $score = $line[0];

	if (!exists($best{$id})) {
		$best{$id} = $count;
	} else {
		if ($score > $best{$id}) {
			$best{$id} = $count;
		}
	}
}
close IN;


open IN, $ARGV[0];
$count = 0;
while (<IN>) {
	$count++;
	next if ($count < 4);
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (@line < 20);
	next if ($line[0] eq 'match');
	next if ($line[0] eq '');
	my $chr = $line[13];
	my $d = $line[8];
	my $id = $line[9];
	my $score = $line[0];
	next if ($best{$id} != $count);

	if (!exists($reads{$id})) {
		$reads{$id} = {s=>$score,str=>""};
	} else {
		if ($score < $reads{$id}->{'s'}) {
			next;
		} else {
			$reads{$id}->{'s'} = $score;
			$reads{$id}->{'str'} = "";
		}
	}
	

	my @lens = split /\,/, $line[18];
	my @gstarts = split /\,/, $line[20];
	my $n = scalar(@lens);
	my $str = "";
	for (my $i=0;$i<@lens;$i++) {
		my $s = $gstarts[$i]+1;
		my $e = $gstarts[$i] + $lens[$i];
		$str .=  "$chr\t$program\texon\t$s\t$e\t0\t$d\t.\tgene_id \"$id\"; transcript_id \"$id\";\n";
	}
	$reads{$id}->{'str'} .= $str;
}
close IN;

foreach(keys %reads) {	
	print $reads{$_}->{'str'};
}
