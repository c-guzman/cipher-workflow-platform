#!/usr/bin/perl -w
#
if (@ARGV < 1) {
	printCMD();
}

sub printCMD {
	print STDERR "\n\tUsage: getLikelyAdapters.pl <fastq file> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-len <#,#,...,#> (lengths of oligos to look for, default: 10,20,25,30)\n";
	print STDERR "\t\t-n <#> (Number of reads to check for overrepresented oligos, default: 1e5)\n";
	print STDERR "\t\t-S <#> (number of overrepresented oligos to report, default: 10)\n";
	print STDERR "\n";
	exit;
}

my $file = $ARGV[0];
my $numSeqsToCheck=1e5;
my $num2report = 10;
my @lens = (10,20,25,30);


for (my $i=1;$i<@ARGV;$i++){ 
	if ($ARGV[$i] eq '-n') {
		$numSeqsToCheck = $ARGV[++$i];
		next;
	} elsif ($ARGV[$i] eq '-S') {
		$num2report = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-len') {
		@lens = split /\,/, $ARGV[++$i];
	} else {
		printCMD();
	}
}


my %seqs = ();
foreach(@lens) {
	my %a = ();
	$seqs{$_}=\%a;
}
my $numSeqs = 0;

my $count = 0;
if ($file =~ /\.gz$/) {
    open IN, "gunzip -c $file |" or die "Could not open $file with gunzip\n";
} elsif ($file =~ /\.bz2$/) {
    open IN, "bunzip2 -c $file |" or die "Could not open $file with bunzip2\n";
} else {
    open IN, $file or die "Could not open $file\n";
}


while (<IN>) {
	$count++;
	if ($count % 4 == 2) {
		$numSeqs++;
		last if ($numSeqs > $numSeqsToCheck);
		print STDERR "\t$numSeqs\n" if ($numSeqs % 1e5 == 0);
			
		my $seq = $_;
		my $len = length($seq);
		foreach(@lens) {
			my $L = $_;
			for (my $i=0;$i<$len-$L;$i++) {
				my $m = substr($seq,$i,$L);
				$seqs{$L}->{$m}++;
			}
		}
	}
}
close IN;
	
foreach(@lens) {
	my $L = $_;
	my @a = sort {$seqs{$L}->{$b} <=> $seqs{$L}->{$a}} keys %{$seqs{$L}};	
	print "\nLength = $L (total=$numSeqs)\n";
	for (my $i=0;$i<$num2report;$i++) {
		print "$a[$i]\t$seqs{$L}->{$a[$i]}\n";
	}
}

