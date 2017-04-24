#!/usr/bin/perl -w
if (@ARGV < 1) {
	print STDERR "<uniprot_sprot.dat>\n";
	exit;
}
my %interpro = ();
my %gene3d = ();
my %pfam = ();
my %prints = ();
my %smart = ();
my %prosite = ();

my $z = 0;
for (my $i=0;$i<@ARGV;$i++){ 

my %prot = ();
my $id = 'NA';
open IN, $ARGV[$i];
while (<IN>) {
	chomp;
	s/\r//g;
	if (/^\/\//) {
		$id = 'NA';
	} elsif (/^ID   (.*?)\s/) {
		$id = $1;
		my %a = ();
		my %acc= ();
		$prot{$id} = \%a;
		$prot{$id}->{'acc'} = \%acc;
	} elsif (/^AC   (.*)\;/) {
		my $all = $1;
		my @acc = split /\;/, $all;
		foreach(@acc) {
			$_ =~ s/\s//g;
			$prot{$id}->{'acc'}->{$_}=1;
		}
	} elsif (/^DR   GeneID\; (.*?)\;/) {
		$a = $1;
		$a =~ s/\.\d+$//;
		$prot{$id}->{'gene'} = $a; 
		print "$id\t$a\n";
		foreach(keys %{$prot{$id}->{'acc'}}) {
			print "$_\t$a\n";
		}
	} elsif (/^DR   PROSITE\; (.*?)\; (.*?)\;/) {
		my $cid = $1;
		my $cterm = $2;
		next if (!exists($prot{$id}->{'gene'}));
		if (!exists($prosite{$cid})) {
			my %a = ();
			$prosite{$cid} = {term=>$cterm,genes=>\%a};
		}
		$prosite{$cid}->{'genes'}->{$prot{$id}->{'gene'}} = 1;
	} elsif (/^DR   SMART\; (.*?)\; (.*?)\;/) {
		my $cid = $1;
		my $cterm = $2;
		next if (!exists($prot{$id}->{'gene'}));
		if (!exists($smart{$cid})) {
			my %a = ();
			$smart{$cid} = {term=>$cterm,genes=>\%a};
		}
		$smart{$cid}->{'genes'}->{$prot{$id}->{'gene'}} = 1;
	} elsif (/^DR   PRINTS\; (.*?)\; (.*?)\./) {
		my $cid = $1;
		my $cterm = $2;
		next if (!exists($prot{$id}->{'gene'}));
		if (!exists($prints{$cid})) {
			my %a = ();
			$prints{$cid} = {term=>$cterm,genes=>\%a};
		}
		$prints{$cid}->{'genes'}->{$prot{$id}->{'gene'}} = 1;
	} elsif (/^DR   Pfam\; (.*?)\; (.*?)\;/) {
		my $cid = $1;
		my $cterm = $2;
		next if (!exists($prot{$id}->{'gene'}));
		if (!exists($pfam{$cid})) {
			my %a = ();
			$pfam{$cid} = {term=>$cterm,genes=>\%a};
		}
		$pfam{$cid}->{'genes'}->{$prot{$id}->{'gene'}} = 1;
	} elsif (/^DR   Gene3D\; (.*?)\; (.*?)\;/) {
		my $cid = $1;
		my $cterm = $2;
		next if (!exists($prot{$id}->{'gene'}));
		if (!exists($gene3d{$cid})) {
			my %a = ();
			$gene3d{$cid} = {term=>$cterm,genes=>\%a};
		}
		$gene3d{$cid}->{'genes'}->{$prot{$id}->{'gene'}} = 1;
	} elsif (/^DR   InterPro\; (.*?)\; (.*?)\./) {
		my $cid = $1;
		my $cterm = $2;
		next if (!exists($prot{$id}->{'gene'}));
		if (!exists($interpro{$cid})) {
			my %a = ();
			$interpro{$cid} = {term=>$cterm,genes=>\%a};
		}
		$interpro{$cid}->{'genes'}->{$prot{$id}->{'gene'}} = 1;
		$z++;
	}
}
close IN;

}
print STDERR "$z genes added\n";
outputDatabase(\%interpro,"interpro.genes");
outputDatabase(\%gene3d,"gene3d.genes");
outputDatabase(\%pfam,"pfam.genes");
outputDatabase(\%prints,"prints.genes");
outputDatabase(\%smart,"smart.genes");
outputDatabase(\%prosite,"prosite.genes");


sub outputDatabase {
	my ($data, $file) = @_;
	open OUT, ">$file";
	foreach(keys %$data) {
		my $id = $_;
		print OUT "$id\t$data->{$id}->{'term'}\t";
		my $c = 0;
		foreach(keys %{$data->{$id}->{'genes'}}) {
			print OUT "," if ($c > 0);
			$c++;
			print OUT "$_";
		}
		print OUT "\n";
	}
	close OUT;
}


