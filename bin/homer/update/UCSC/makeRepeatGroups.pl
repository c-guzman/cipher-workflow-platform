#!/usr/bin/perl -w
my $col = 11;
my $rnaDefFlag = 1;

if (@ARGV < 2) {
	print STDERR "\n\tUsage: makeRepeatGroups.pl <directory with *_rmsk.txt files> <output directory>\n";
	print STDERR "\tEach repeat family will get a *.ann.txt file in the output directory\n";
	print STDERR "\tAlso, a repeats.ann.txt file will be created containing all repeats,\n";
	print STDERR "\tand a repeats.rna file will be created for calculating repeat expression\n";
	print STDERR "\tThese files will be in current directory\n";
	print STDERR "\n";
	exit;
} 
my $outputDir = $ARGV[1];

my $directory = $ARGV[0];
my @files = ();

my $tmp = rand();
my $tmpfile = $tmp . ".tmp";
`ls "$directory"/*rmsk.txt > $tmpfile`;
open IN, $tmpfile;
while (<IN>) {
	chomp;
	push(@files, $_);
}
close IN;
`rm $tmpfile`;


my %repeats = ();
my $startFlag = 0;
my $ID = 0;

open OUT2, ">repeats.ann.txt";
open OUT3, ">repeats.rna";
for (my $i=0;$i<@files;$i++) {

	my $file =  $files[$i];
	print STDERR "\tProcessing $file\n";

	open IN, $file;
	while (<IN>) {
		chomp;
		my @line = split /\t/;
		my $chr = $line[5];
		my $start = $line[6]+1;
		my $end = $line[7];

		my $repStart = $line[15];

		my $dir = 0;
		if ($line[9] eq '-') {
			$dir = 1;
		}
		if ($startFlag==1) {
			if ($dir == 0) {
				$tss = $start + abs($line[13]);
			} else {
				$tss = $end - abs($line[15])+1;
			}
			$start = $tss-2000;
			$end = $tss+2000;
				
		}
	
		$ID++;

		my $div = $line[2]/1000;	

		my $ogStart = $start-$line[13]+1;
		my $ogEnd = $end-$line[15];
		if ($dir == 1) {
			$ogStart = $start+$line[13]-1;
			$ogEnd = $end+$line[15];
		}

		my $name1 = $line[10] . "|" . $line[11] . "|" . $line[12];
		my $name2 = $line[11] . "|" . $line[12];
		my $name3 = $line[11];
		my $name = $name1 . "-HOMER$ID";
		my $output = "$name\t$chr\t$start\t$end\t$dir\t$line[11]\n";
		my $output2 = "$name\t$chr\t$start\t$end\t$dir\tE1:$start\t$div\t$ogStart\t$ogEnd\n";
		#print "$name\t$chr\t$start\t$end\t$dir\n";
		print OUT2 $output;
		print OUT3 $output2;

		my @names = ($name1, $name2, $name3);
		foreach(@names) {
			if (!exists($repeats{$_})) {
				my @a = ();
				$repeats{$_} = \@a;
			}
			push(@{$repeats{$_}}, $output);
		}
	}
	close IN;
}	
close OUT2;
close OUT3;


foreach(keys %repeats) {
	my $rep = $_;
	my $filename= $_;
	$filename =~ s/\(//g;
	$filename =~ s/\)//g;
	$filename =~ s/\//\_/g;
	my $fullFile = $outputDir . "/" . $filename . ".ann.txt";
	open OUT, ">$fullFile" or die "Couldn't open $fullFile\n";
	foreach(@{$repeats{$rep}}) {
		print OUT "$_";
	}
	close OUT
}



