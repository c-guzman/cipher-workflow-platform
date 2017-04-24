#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";

use POSIX;
use Statistics;

my $dist = 1000;
my $minSNPs =10;
my $maxCPUs = 1;

sub printCMD() {
	print STDERR "\n\tgetGWASoverlap.pl <gwas catolog file> -p <peak file1> [peak file2] ... [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-d <#> (Overlap distance, default: $dist)\n";
	print STDERR "\t\t-min <#> (minimum number of significant SNPs to consider study, default: $minSNPs)\n";
	print STDERR "\t\t-cpu <#> (number of threads to use, default: $maxCPUs)\n";
	print STDERR "\t\t-empirical <#> (calculate FDR q-values empirically, <#> number of randomizations)\n";
	print STDERR "\n\tThe gwas catalog file can be downloaded from UCSC annotation database:\n";
	print STDERR "\t\ti.e.: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/gwasCatalog.txt.gz\n";
	print STDERR "\n";
	exit;
}



if (@ARGV < 1) {
	printCMD();
}


my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";
my $tmpFile3 = $rand . ".3.tmp";
my $tmpEmpirical = $rand . ".empirical";

my $gwasFile = $ARGV[0];
my @peakFiles = ();
my $numEmpirical = 0;

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-d') {
		$dist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$minSNPs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-empirical') {
		$numEmpirical = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-p') {
		$i++;
		while ($i<@ARGV && $ARGV[$i] !~ /^-/) {
			push(@peakFiles, $ARGV[$i++]);
		}
		$i--;
	} else {
		print STDERR "Couldn't recognize $ARGV[$i]\n";
		printCMD();
	}
}

my %limits = ();


print STDERR "\n\tMinimum Risk SNPs per study: $minSNPs\n";
print STDERR "\tMaximum Distance between Peaks and risk SNPs: $dist\n\n";

print STDERR "\n\tAnalyzing peak files:\n";
my @fdrPvalues = ();

for (my $i=0;$i<@peakFiles;$i++) {
	print STDERR "\t\t$peakFiles[$i]\n";
	if ($numEmpirical > 0) {
		`bed2pos.pl "$peakFiles[$i]" -check > "$tmpFile"`;
		open IN, $tmpFile;
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			next if ($line[0] =~ /^#/);
			next if (@line < 5);
			if (!exists($limits{$line[1]})) {
				$limits{$line[1]} =0;
			}
			$limits{$line[1]} = $line[3] if ($limits{$line[1]} < $line[3]);
		}
		close IN;
		`rm $tmpFile`;
		my @a = ();
		push(@fdrPvalues, \@a);
	}
}
print STDERR "\n";

my $limitTotal = 0;
if ($numEmpirical > 0) {
	#print STDERR "\tChromosome Limits:\n";
	foreach(keys %limits) {
		#print STDERR "\t\t$_\t$limits{$_}\n";
		$limitTotal += $limits{$_};
	}
	#print STDERR "\tApprox Size: $limitTotal\n\n";
}


my $gwas = parseGWAScatalog($gwasFile);

my @results = ();
my @studies = ();

my $numStudy = scalar(keys %$gwas);
my $count = 0;
my $count2 = 0;
foreach(keys %$gwas) {
	$count2++;
	my $study = $_;
	my @snps = keys %{$gwas->{$study}};
	my $N = scalar(@snps);
	next if (@snps < $minSNPs);
	print STDERR "\tAnalyzing Study $count2 of $numStudy ($study, $N total risk SNPs)\n";

	my @peakResults = ();

	open OUT, ">$tmpFile";
	foreach(@snps) {
		my $allele = $_;
		my $p = $gwas->{$study}->{$allele};
		print OUT "$p\n";
	}
	close OUT;

	$count++;
	if ($count == 1) {
		print "Study";
		print "\tTotal risk SNPs";
		for (my $i=0;$i<@peakFiles;$i++) {
			print "\t$peakFiles[$i] Overlap";
			print "\t$peakFiles[$i] p-value";
			if ($numEmpirical > 0) {
				print "\t$peakFiles[$i] q-value FDR/Empirical(n=$numEmpirical)";
			} else {
				print "\t$peakFiles[$i] q-value FDR/Benjamini";
			}
			print "\t$peakFiles[$i] SNPs";
		}
		print "\n";
	}
	my $numDsnps = $N;
	push(@studies, "$study\t$N");
	

	for (my $i=0;$i<@peakFiles;$i++) {
		`cp "$peakFiles[$i]" "$tmpFile2"`;

		`mergePeaks "$tmpFile" "$tmpFile2" -matrix $rand  -d $dist > /dev/null 2> /dev/null`;

		my $logp = 0;

		my $c = 0;
		open IN, "$rand.logPvalue.matrix.txt";
		while (<IN>) {
			$c++;
			next if ($c < 2);
			s/\r//g;
			my @line = split /\t/;
			$logp = $line[2];
			last;
		}
		close IN;
		my $pvalue = 1;
		if ($logp < 0.1) {
			$pvalue = exp($logp);
		} else {
			$pvalue = 1-exp(-1*$logp);
		}
		`rm $rand.*matrix.txt`;

		`mergePeaks "$tmpFile" "$tmpFile2" -cobound 1 -prefix $rand -d $dist > /dev/null 2> /dev/null`;

		my @bound = ();
		my $snpStr = "";
		open IN, "$rand.coBoundBy1.txt";	
		while (<IN>) {
			chomp;
			next if (/^#/);
			my @line = split /\t/;
			push(@bound, $line[0]);
			$snpStr .= "$line[0]" . ",";
		}
		close IN;
		my $N = scalar(@bound);	
		`rm $rand.coBoundBy*`;

		my $fdr = 1;

		if ($numEmpirical > 0) {
			my $fileStr = "";
			my @chrs = keys %limits;
			for (my $j=0;$j<$numEmpirical;$j++) {
				my @rands = ();
				for (my $k=0;$k<$numDsnps;$k++) {
					push(@rands, floor($limitTotal*rand()));
					#print STDERR "$j $k $rands[$k]\n";
				}
				@rands = sort {$a <=> $b} @rands;
				my $f = "$tmpEmpirical.$j.tmp";
				$fileStr .= " $f";
				open OUT, ">$f";
				my $chrIndex = 0;
				my $totalIndex = 0;
				foreach(@rands) {
					my $r = $_;	
					while ($totalIndex+$limits{$chrs[$chrIndex]}< $r) {
						$totalIndex += $limits{$chrs[$chrIndex++]};
					}
					my $p = $r-$totalIndex;
					my $p2 = $p+1;
					print OUT "$r\t$chrs[$chrIndex]\t$p\t$p2\t+\n";
				}
				close OUT;
			}
			#print STDERR "`mergePeaks $tmpFile2 -cobound 1 -matrix $rand -d $dist $fileStr > /dev/null`;\n";
			`mergePeaks "$tmpFile2" -cobound 1 -matrix $rand -d $dist $fileStr -gsize $limitTotal > /dev/null 2> /dev/null`;

			my @pvalues = ();
			open IN, "$rand.logPvalue.txt";
			my $c = 0;
			while (<IN>) {
				$c++;
				next if ($c < 2);
				chomp;
				my @line = split /\t/;
				my $p = 1;
				if ($line[1] <= 0.1) {
					$p = exp($line[1]);
				} else {
					$p = 1-exp(-1*$line[1]);
				}
				push(@pvalues, $p);
				#print STDERR "\t$pvalue\t$p\n";
			}
			close IN;
			`rm -f $fileStr`;
			push(@{$fdrPvalues[$i]},\@pvalues);

		}

		my $d = {n=>$N,p=>$pvalue,snps=>$snpStr,f=>$fdr};	
		push(@peakResults, $d);
	}
	push(@results, \@peakResults);
	`rm $tmpFile $tmpFile2`;

}


for (my $i=0;$i<@peakFiles;$i++) {
	my @pvalues = ();
	foreach(@results) {
		push(@pvalues, $_->[$i]->{'p'});
	}
	my $fdrs = "";
	if ($numEmpirical < 1) {
		$fdrs = Statistics::benjaminiFDR(\@pvalues);
	} else {
		my $aa = scalar(@pvalues);
		my $bb = scalar(@{$fdrPvalues[$i]});
		$fdrs = Statistics::empiricalFDR(\@pvalues, $fdrPvalues[$i]);
	}
	for (my $j=0;$j<@results;$j++) {
		$results[$j]->[$i]->{'f'} = $fdrs->[$j];
	}
}

for (my $i=0;$i<@results;$i++) {
	print $studies[$i];
	foreach(@{$results[$i]}) {
		print "\t$_->{'n'}\t$_->{'p'}\t$_->{'f'}\t$_->{'snps'}";
	}
	print "\n";
}	

`rm -f $rand*`;

sub parseGWAScatalog {
	my ($file) = @_;
	my %gwas = ();
	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $chr = $line[1];
		my $start = $line[2];
		my $end = $line[3];
		my $snp = $line[4];
		my $pubmed = $line[5];
		my $study = $line[10] . "($pubmed)";
		my $allele = $line[15];
		my $pvalue = $line[17];
		my $position = "$allele\t$chr\t$start\t$end\t0";
	
		if (!exists($gwas{$study})) {
			my %snps = ();
			$gwas{$study} = \%snps;
		}
		$gwas{$study}->{$allele} = $position;
	}
	close IN;
	return \%gwas;
}
