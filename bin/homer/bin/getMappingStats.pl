#!/usr/bin/perl -w


if (@ARGV < 1) {
	print STDERR "\n\tUsage: getMappingStats.pl [options] <sample directory> [sample directory2] ...\n";
	print STDERR "\tWill print stats to stdout\n";
	print STDERR "\n\tProgram looks for:\n";
	print STDERR "\t\t*.lengths file containing trimming stats\n";
	print STDERR "\t\t*.bowtie2.log containing bowtie2 mapping stats\n";
	print STDERR "\t\t*.final.out containing STAR mapping stats\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-min <#> (minimum length to consider a sequence an adapter-dimer, def. 15)\n";
	print STDERR "\t\t-genome <genome version> (limit analysis to alignment for this genome)\n";
	print STDERR "\n";
	exit;
}

$minlen = 15;
my $GENOME = "";

$tmpFile = rand() . ".tmp";

my @files = ();
for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-min') {
		$minlen = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-genome') {
		$GENOME = $ARGV[++$i];
	} else {
		push(@files, $ARGV[$i]);
	}
}



print "Experiment File\tGenome\tTotal reads\t\% Adapter Dimers";
print "\t\% Aligned\t\% Unique Match\t\% Multimappers\t\% unmapped";
print "\tAligner\n";



foreach(@files) {
	my $file = $_;
	print "$file";
	#print STDERR "FILE: $file";

	my $trimFile = '';
	my $bowtie2File = '';
	my $starFile = '';
	my $genome = 'NA';
	my $aligner = 'NA';

	my $m = $file;
	my $DIR = ".";
	my $file2 = $file;
	if ($m =~ /^(.+)\/([^\/]+)$/) {
		$DIR = $1;
		$file2 = $2;
	}
	#print STDERR "||$file||\n";
	#print STDERR "||$DIR||\t||$file2||\n";
	#`sleep 100`;
	#exit;
	#
	$file2 =~ s/\.gz$//;
	$file2 =~ s/fastq.*$/fastq/;
	$file2 =~ s/fq.*$/fq/;
	$file2 =~ s/\.bam$//;
	$file2 =~ s/\.sam$//;
	`ls $DIR/$file2* > $tmpFile`;
	my $error = "ls $m*\n";
	open IN, $tmpFile;
	while (<IN>) {
		$error .= $_;
		chomp;
		s/\r//g;
		if ($_ =~ /\.lengths$/) {
			$trimFile = $_;
		} elsif ($_ =~ /\.(.*?)\.bowtie2\.log/) {
			my $cgenome = $1;
			$cgenome =~ s/^.*\.//;
			if ($GENOME ne '') {
				if ($cgenome ne $GENOME) {
					next;
				}
			}
			$genome = $cgenome;
			$bowtie2File = $_;
		} elsif ($_ =~ /\.(.*?)\.Log\.final\.out/ || $_ =~ /\.(.*?)\.STAR\.log/) {
			my $cgenome = $1;
			$cgenome =~ s/^.*\.//;
			if ($GENOME ne '') {
				if ($cgenome ne $GENOME) {
					next;
				}
			}
			$genome = $cgenome;
			$starFile = $_;
		}
	}
	close IN;
	`rm $tmpFile`;

	my $adapterDimers = 'NA';
	my $total = 0;
	if ($trimFile ne '') {
		($adapterDimers,$goodReads) = readTrimLengthFile($trimFile);
		if ($goodReads > 0) {
			$total = $goodReads;
		}
		if ($adapterDimers ne 'NA') {
			$total += $adapterDimers;
			$adapterDimers /= $total if ($total > 0);
			$adapterDimers = sprintf("%.2f",$adapterDimers*100) . '%';
		}
	}
	my ($mt,$um,$mm,$un) = (0,0,0,0);
	if ($starFile ne '') {
		($mt,$um,$mm,$un) = readSTARFile($starFile);
		$aligner = 'STAR';
	} elsif ($bowtie2File ne '') {
		($mt,$um,$mm,$un) = readBowtie2File($bowtie2File);
		$aligner = 'bowtie2';
	} else {
		print STDERR "!!! No alignment log found for $file: !!!\n";
		print STDERR "$error";
		print STDERR "\n";
	}

	$total = 1 if ($total < 1);
	$total = $mt if ($total < $mt);

	print "\t$genome";
	print "\t$total";
	print "\t$adapterDimers";

	$total = 1 if ($total < 0);
	$p = sprintf("%.2f",($um+$mm)/$total*100) . '%';
	print "\t$p";
	$p = sprintf("%.2f",($um)/$total*100) . '%';
	print "\t$p";
	$p = sprintf("%.2f",$mm/$total*100) . '%';
	print "\t$p";
	$p = sprintf("%.2f",$un/$total*100) . '%';
	print "\t$p";


	print "\t$aligner\n";
}


sub readBowtie2File {
	my ($file) = @_;
	open IN, $file or print STDERR "Can't open: $file\n";
	my $totalReads = 0;
	my $uniqReads = 0;
	my $multiReads = 0;
	while (<IN>) {
		chomp;
		s/\r//g;
		if (/(\d+) reads; of these:/) {
			$totalReads = $1;
		} elsif (/(\d+) \(.*\) aligned exactly 1 time/) {
			$uniqReads = $1;
		} elsif (/(\d+) \(.*\) aligned >1 times/) {
			$multiReads = $1;
		}
	}
	close IN;
	my $unmappedReads = $totalReads - $uniqReads - $multiReads;
	return ($totalReads, $uniqReads, $multiReads, $unmappedReads);
}

sub readSTARFile {
	my ($file) = @_;
	open IN, $file;
	#print STDERR "\t$file\n";
	my $totalReads = 0;
	my $uniqReads = 0;
	my $multiReads = 0;
	while (<IN>) {
		chomp;
		s/\r//g;
		if (/Number of input reads \|\s(\d+)/) {
			$totalReads = $1;
		} elsif (/Uniquely mapped reads number \|\s(\d+)/) {
			$uniqReads = $1;
		} elsif (/Number of reads mapped to multiple loci \|\s(\d+)/) {
			$multiReads = $1;
		}
	}
	close IN;
	my $unmappedReads = $totalReads - $uniqReads - $multiReads;
	#print STDERR "\t$unmappedReads = $totalReads - $uniqReads - $multiReads\n";
	return ($totalReads, $uniqReads, $multiReads, $unmappedReads);
}	
	


sub readTrimLengthFile {

	my ($file) = @_;

	my $adapterDimers = 0;
	my $goodReads = 0;
	
	if (open IN, $file) {
		while (<IN>) {
			chomp;
			s/\r//g;
			my @line = split /\t/;
			next if ($line[0] =~ /Length/);
			if ($line[0] < $minlen) {
				$adapterDimers += $line[1];
			} else {
				$goodReads += $line[1];
			}
		}
		close IN;
	} else {
		$adapterDimers = 'NA';
	}
	return ($adapterDimers, $goodReads);
}
		
