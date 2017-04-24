#!/usr/bin/perl -w


print "Filename\tName\tConsensus\tFactor Name\tDBD\tSubmotif\tCell Type\tIP\tAssay\tGEO ID\tSource\n";
foreach(@ARGV) {
	my $file = $_;
	my $name = '';
	open IN, $file;
	while (<IN>) {
		chomp;
		s/\r//g;
		if (/^>/) {
			my @line = split /\t/;
			my $cons = $line[0];
			my $name = $line[1];
			$cons =~ s/^\>//;
			my @cat = split /\//,$name;
			my $subsite = 'canonical';
			my $factorname = $cat[0];

			my $DBD = "NA";
			if ($cat[0] =~ s/(.*?)\((.*?)\)//) {
				$factorname = $1;
				$DBD = $2;
			}
			if ($cat[0] =~ /,(.*?)$/) {
				$subsite = $1;
			}
			my $cells = "NA";
			my $ip = "NA";
			my $assay = "NA";
			my $geo = "NA";
			my $source = $cat[2];
			if ($cat[1] eq 'Promoter') {
				$assay = "Promoter Enrichment";
			} else {
				if ($cat[1] =~ s/\((.*?)\)$//) {
					$geo = $1;
				}
				if ($cat[1] =~ /^(.*?)-(.*?)-(.*)$/) {
					$cells = $1;
					$ip = $2;
					$assay = $3;
				}
			}


			print "$file\t$name\t$cons\t$factorname\t$DBD\t$subsite\t$cells\t$ip\t$assay\t$geo\t$source\n";
		}
	}
	close IN;
}
	

