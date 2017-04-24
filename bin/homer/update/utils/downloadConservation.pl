#!/usr/bin/perl -w

my $baseURL = "http://hgdownload.cse.ucsc.edu/goldenPath/hg17/phastCons17way/";

my @chr = ();
for (my $i=1;$i<=22;$i++) {
	push(@chr, "chr" . $i);
}
push(@chr, "chrM");
push(@chr, "chrUn");
push(@chr, "chrX");
push(@chr, "chrY");
my @newChr = ();
foreach(@chr) {
	push(@newChr, $_, "$_" . "_random");
}

foreach(@newChr) {
	my $chr = $_;
	#my $fileName = $chr . ".placental.pp.data.gz";
	my $fileName = $chr . ".gz";
	my $url = $baseURL . $fileName;
	`wget -O $fileName $url`;
}
