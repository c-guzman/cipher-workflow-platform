#!/usr/bin/env perl
use warnings;

if (@ARGV < 1) {
	print STDERR " [-len #] <qseq.txt file 1> [qseq.txt file 2] ...\n";
	print STDERR " or qseq2fastq.pl -all : this will do all lanes\n";
	exit;
}
if ($ARGV[0] eq '-all') {
	print STDERR "\tCreating fq files from all qseq files.\n";
	for (my $i=1;$i<=8;$i++){ 
		print STDERR "\tLane $i:\n";
		`ls -1 s_$i*_qseq.txt > .ls`;
		open IN, ".ls";
		my @files = ();
		while (<IN>) {
			chomp;
			push(@files, $_);
		}
		close IN;
		if (@files < 1) {
			print STDERR "!! No files found for lane $i\n";
			next;
		}
		open OUT, ">s_$i" . "_sequence.qseq.txt";
		for (my $j=0;$j<@files;$j++) {
			print STDERR "\t\tparsing $files[$j]\n";
			open IN, $files[$j];
			while (<IN>) {
				chomp;
				my @line = split /\t/;
				next if (@line < 10);
				my $id = $line[0] . "_" . $line[1] . "_" . $line[2] . "_" . $line[3] . "_" . $line[4] . "_" . $line[5];
				my $seq = $line[8];
				my $qual = $line[9];
				print OUT "\@$id\n";
				print OUT "$seq\n";
				print OUT "+$id\n";
				print OUT "$qual\n";
			}
			close IN;
		}
		close OUT;
		`rm .ls`;
	}
	exit;
}

my $len = 0;
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-len') {
		$len = $ARGV[++$i];
		next;
	}
		
	open IN, $ARGV[$i];
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $line[0] . "_" . $line[1] . "_" . $line[2] . "_" . $line[3] . "_" . $line[4] . "_" . $line[5];
		my $seq = $line[8];
		my $qual = $line[9];
		if ($len > 0) {
			next if (length($seq) != $len);
		}
		print "\@$id\n";
		print "$seq\n";
		print "+$id\n";
		print "$qual\n";
	
	}
	close IN;

}
