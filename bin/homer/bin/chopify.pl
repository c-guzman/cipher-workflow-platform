#!/usr/bin/env perl
use warnings;


if (@ARGV < 2) {
	print STDERR "\n\tchopify.pl - chops peak file into several smaller peaks that tile the original regions\n";
	print STDERR "\n\tUsage: chopify.pl <peak/BED file> <size # in bp>\n\n\tOutput peak file sent to stdout\n\n";
	exit;
}

my $res = $ARGV[1];

my $rand = rand();
my $tmpFile = $rand . ".tmp";

`bed2pos.pl "$ARGV[0]" -check > "$tmpFile"`;
open IN, $tmpFile;
while (<IN>) {
	chomp;
	s/\r//g;
	my $og = $_;
	my @line = split /\t/;
	my $s = $line[2];
	my $e = $line[3];
	if ($e-$s > $res) {
		my $z = 1;
		for (my $i=$s;$i<$e;$i+=$res) {
			my $s1 = $i;
			my $e1 = $s1+$res;
			print "$line[0]" . "-chop" . "$z\t$line[1]\t$s1\t$e1\t$line[4]";
			for (my $i=5;$i<@line;$i++) {
				print "\t$line[$i]";
			}
			print "\n";
			$z++;
		}
	} else {
		print "$og\n";
	}
}
close IN;

`rm -f $tmpFile`;
