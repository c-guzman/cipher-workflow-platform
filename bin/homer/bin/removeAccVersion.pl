#!/usr/bin/env perl
use warnings;
#
if (@ARGV < 1) {
	print STDERR "<intput file>\n";
	print STDERR "Removes the NM_012123.1 to be NM_012123 in the first column\n";
	exit;
}

my $ensemblFlag = 0;
if (@ARGV > 1) {
	$ensemblFlag = 1;
}

open IN, $ARGV[0] or die "Couldn't open file: \"$ARGV[0]\"\n";
while (<IN> ){	
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if ($_ eq '' || @line < 1);

	#if ($ensemblFlag) {
	if ($line[0] =~ /^transcript\:.*_T\d\d/) {
		$line[0] =~ s/^transcript\://;
		$line[0] =~ s/^gene\://;
		$line[0] =~ s/_T\d\d$//;
	} else {
		#need top protect unigene ids!
		if ($line[0] =~ /^..\.\d+/) {  #unigene
			if ($line[0] =~ /^..\.\d+\.\d+$/) {
				$line[0] =~ s/\.\d+$//;
			}
		} elsif ($line[0] =~ /^...\.\d+/) { #unignee with 3 letter ids
			if ($line[0] =~ /^...\.\d+\.\d+$/) {
				$line[0] =~ s/\.\d+$//;
			}
		} else  {
			$line[0] =~ s/\.\d+$//;
		}
	}
	my $c = 0;
	foreach(@line) {
		print "\t" if ($c>0);
		$c++;
		print "$_";
	}
	print "\n";
}
close IN;
