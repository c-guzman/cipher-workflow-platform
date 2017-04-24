#!/usr/bin/env perl
use warnings;


use POSIX;
if (@ARGV < 2) {
	print STDERR "<suffix i.e. \".cgfreq\"> <seq file1> <seq file2>...\n";
	print STDERR "outputs to \"group[suffix].[group value].seq\"\n";
	exit;
}


my $colIndex = 2;
my $min = 0;
my $max = 1;
my $step = 0.1;

my $suffix = $ARGV[0];

my %group= ();
for (my $i=$min;$i<=$max;$i+=$step) {
	my $file = ">group" . $suffix . ".$i.seq";
	my $handle;
	local $blah;
	open $blah, $file;
	print STDERR " $file $i\n";
	$group{$i} = $blah;
}

for (my $z=1;$z<@ARGV;$z++) {
print STDERR "\t$ARGV[$z]\n";
	my $valFile = $ARGV[$z] . $suffix;
	my %values = ();
	open IN, $valFile;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		$values{$line[0]} = $line[$colIndex];
	}
	close IN;
	
	open IN, $ARGV[$z];
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		next if (!exists($values{$line[0]}));
		my $v = $values{$line[0]};
		my $index = floor(($v-$min)/$step)*$step;
		if (!exists($group{$index})) {
			print STDERR "OUT of RANGE!!!\n";
		} else {
			my $f = $group{$index};
			print $f "$_\n";
		}
	}
	close IN;

}

for (my $i=$min;$i<=$max;$i+=$step) {
	my $file = "group" . $suffix . ".$i.seq";
	close $group{$i};
}
