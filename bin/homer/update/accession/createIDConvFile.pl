#!/usr/bin/perl -w

if (@ARGV < 2) {
	print STDERR "\n\tusage: createIDConvFile.pl <taxid> <org> [extra file] ...\n";
	print STDERR "\n\tFormat of extra file: feature <tab> locuslink \n";
	print STDERR "\n\tThis program assumes the following files are in the current directory:\n";
	print STDERR "\t(most are from the NCBI Gene Database FTP site: ftp://ftp.ncbi.nih.gov/gene/DATA/)\n";
	print STDERR "\t\tgene2accession\n";
	print STDERR "\t\tgene2unigene\n";
	print STDERR "\t\tgene2refseq\n";
	print STDERR "\t\tgene2ensembl\n";
	print STDERR "\t\tgene_info\n";
	print STDERR "\t\tgene_history\n";
	print STDERR "\t\t\n";
	exit;
}
my $org = $ARGV[1];
my $dir = "data/";

my $taxid = $ARGV[0];

my %unigene = ();
my %gene = ();
my %ensembl = ();
my %refseq = ();

my %geneIDs = ();
my %unigeneIDs = ();
my %refseqIDs = ();
my %ensemblIDs = ();

print STDERR "\tparsing gene2accession\n";
open IN, "gene2accession" or die "!!! Could not find gene2accession\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $oid = $line[0];
	next if ($oid ne $taxid);
	my $locus = $line[1];
	$geneIDs{$locus} = 1;
	$gene{$locus} = $locus;

	my @ids = ();
	push(@ids, $line[3]);
	push(@ids, $line[5]);
	push(@ids, $line[7]);
	for (my $i=0;$i<@ids;$i++) {	
		$ids[$i] =~ s/\.(.*)$//;
		next if ($ids[$i] eq '-');
		next if ($ids[$i] eq '');
		$gene{$ids[$i]} = $locus;
	}
}
close IN;

print STDERR "\tparsing gene2unigene\n";
open IN, "gene2unigene" or die "Could not find gene2unigene\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (!exists($geneIDs{$line[0]}));
	$gene{$line[1]} = $line[0];
	if (!exists($unigene{$line[0]})) {
		$unigene{$line[0]} = $line[1];
	} else {
		my $alt = $unigene{$line[0]};
		my $cur = $line[1];
		my $n1 = $alt;
		my $n2 = $cur;
		$n1 =~ s/^.*?\.//;
		$n2 =~ s/^.*?\.//;
		if ($n2 < $n1) {
			$unigene{$line[0]} = $line[1];
		}
	}
}
close IN;

print STDERR "\tparsing gene2refseq\n";
open IN, "gene2refseq" or die "Could not find gene2refseq\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if ($line[0] ne $taxid);
	my $locus = $line[1];
	$geneIDs{$locus} = 1;
	$gene{$locus} = $locus;
	
	my $mrna = $line[3];
	my $prot = $line[5];
	$mrna =~ s/\.(.*)$//;
	next if ($mrna eq '-');
	next if ($mrna eq '');
	$gene{$mrna} = $locus;
	if (!exists($refseq{$locus})) {
		$refseq{$locus} = $mrna;
	} else {
		my $altRefSeq = $refseq{$locus};
		if ($altRefSeq =~ /^X/ && $mrna !~ /^X/) {
			$refseq{$locus} = $mrna;
		} elsif ($altRefSeq !~ /^X/ && $mrna =~ /^X/) {
			next;
		} elsif ($altRefSeq =~ /^NR/ && $mrna =~ /^NM/) {
			$refseq{$locus} = $mrna;
		} elsif ($altRefSeq =~ /^NM/ && $mrna =~ /^NR/) {
			next;
		} else {
			my $n1 = $mrna;
			my $n2 = $altRefSeq;
			$n1 =~ s/^.*?_//;
			$n2 =~ s/^.*?_//;
			if ($n1 < $n2) {
				$refseq{$locus} = $mrna;
			} else {
				next;
			}	
		}
	}
}
close IN;

print STDERR "\tparsing gene2ensembl\n";
open IN, "gene2ensembl" or die "Could not find gene2ensembl\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my $oid = $line[0];
	next if ($oid ne $taxid);
	my $locus = $line[1];
	$geneIDs{$locus} = 1;
	$gene{$locus} = $locus;

	my @ids = ();
	push(@ids, $line[2]);
	push(@ids, $line[4]);
	push(@ids, $line[6]);
	for (my $i=0;$i<@ids;$i++) {	
		$ids[$i] =~ s/\.(.*)$//;
		next if ($ids[$i] eq '-');
		next if ($ids[$i] eq '');
		$gene{$ids[$i]} = $locus;
	}
}
close IN;



my %orfs = ();
my %names = ();

print STDERR "\tparsing gene_info\n";
open IN, "gene_info" or die "Could not find gene_info\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if ($line[0] ne $taxid);
	my $locus = $line[1];
	$geneIDs{$locus} = 1;
	$gene{$locus} = $locus;
	my @dbrefs = split /\|/, $line[5];
	foreach(@dbrefs) {
		if (/Ensembl\:(.*)/) {
			my $ensb = $1;
			#my @ids = split /\|/, $ensb;
			#foreach(@ids) {
			#s/Ensembl\://;
			#s/^\s*//;
			#s/\s*$//;
			#if ($_ !~ /^Vega/i) {
			$gene{$ensb} = $locus;
			if (!exists($ensembl{$locus})) {
				$ensembl{$locus} = $ensb;
			} else {
				my $cid = $ensb;
				my $alt = $ensembl{$locus};
				my $n1 = $cid;
				my $n2 = $alt;
				$n1 =~ s/^.*G//;
				$n2 =~ s/^.*G//;
				$n1 =~ s/^.*I//;
				$n2 =~ s/^.*I//;
				$n1 =~ s/^.*T//;
				$n2 =~ s/^.*T//;
				if ($n1 < $n2) {
					$ensembl{$locus} = $cid;
				}
			}
		}
	}
	my $name = $line[2];
	my $orf = $line[3];
	my $alias = $line[4];
	if ($orf ne '' && $orf ne '-') {
		$gene{$orf} = $locus;
		$orfs{$locus} = $orf;
	}
	if ($name ne '' && $name ne '-') {
		$gene{$name} = $locus;
		$names{$locus} = $name;
	}
	
	my $desc = "\t$name\t$alias\t$orf\t$line[7]\t$line[8]\t$line[9]\n";
	$description{$locus} = $desc;
}
close IN;

print STDERR "\tLooking through gene_history\n";
open IN, "gene_history" or die "Could not open gene_history\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if ($line[0] ne $taxid);
	my $locus = $line[1];
	my $old = $line[2];
	$gene{$old} = $locus;
	if (exists($geneIDs{$old})) {
		print STDERR "geneIDs contains old gene ids!!\n";
	}
	#$geneIDs{$old} = $locus;
}
close IN;
	



my $descHeader = "\tname\talias\torf\tchromosome\tdescription\ttype\n";

for(my $i=2;$i<@ARGV;$i++) {
	print STDERR "\tParsing $ARGV[$i]\n";
	open IN, $ARGV[$i];
	while (<IN>) {
		chomp;
		s/\r//g;
		s/\"//g;
		my @line = split /\t/;
		if ($line[1] =~ /^(..\.\d+)/) {
			$line[1] = $1;
		} else {
			$line[1] =~ s/\.(.*?)$//;
		}
		if (!exists($gene{$line[1]})) {
			#print STDERR "$line[0]\t$line[1]\n";
			next;
		}
		$gene{$line[0]} = $gene{$line[1]};
	}
	close IN;
}

print STDERR "\tOutputing data\n";
open OUT, ">$org" . "2gene.tsv";
foreach (keys %gene) {
	my $acc = $_;
	my $gid = $gene{$acc};
	print OUT "$acc\t$gid";
	my $ug = '';
	if (exists($unigene{$gid})) {
		$ug = $unigene{$gid};
	}
	if ($org eq 'yeast') {
		if (exists($orfs{$gid})) {
			#$ug = $orfs{$gid};
		}
	}
	print OUT "\t$ug";

	my $refseq = '';
	if (exists($refseq{$gid})) {
		$refseq = $refseq{$gid};
	}
	print OUT "\t$refseq";

	my $ensembl = '';
	if (exists($ensembl{$gid})) {
		$ensembl = $ensembl{$gid};
	}
	print OUT "\t$ensembl";

	my $orf = "";
	if (exists($orfs{$gid})) {
		$orf = $orfs{$gid};
	}
	print OUT "\t$orf";

	my $name = "";
	if (exists($names{$gid})) {
		$name = $names{$gid};
	}
	print OUT "\t$name";
	print OUT "\n";
}
close OUT;

open OUT, ">$org.description";
print OUT "GeneID\tUnigene\tRefSeq\tEnsembl" . $descHeader;
foreach(keys %description) {
	my $gid = $_;
	print OUT "$gid";
	my $ug = '';
	if (exists($unigene{$gid})) {
		$ug = $unigene{$gid};
	}
	print OUT "\t$ug";
	my $refseq = '';
	if (exists($refseq{$gid})) {
		$refseq = $refseq{$gid};
	}
	print OUT "\t$refseq";
	my $ensembl = '';
	if (exists($ensembl{$gid})) {
		$ensembl = $ensembl{$gid};
	}
	print OUT "\t$ensembl";
	print OUT "$description{$gid}";
}
close OUT;

