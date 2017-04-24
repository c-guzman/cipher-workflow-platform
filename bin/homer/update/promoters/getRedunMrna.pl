#!/usr/bin/perl -w
my %conv = ();
my %ids = ();
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split  /\t/;
	my $id = $line[1];
	my $name = $line[5];
	if ($name ne '') {
		if (!exists($conv{$name})) {
			my %a = ();
			$conv{$name} = \%a;
		}
		$conv{$name}->{$id} = 1;
#print STDERR "mapping |$name| to |$id|\n";
		$ids{$id} = $name;
	}
}
close IN;

foreach(keys %ids) {
	my @a = ();
	my $id = $_;
	my $name = $ids{$id};
#print STDERR "looking up |$name| to |$id|\n";
	next if ($name eq '');
	foreach(keys %{$conv{$name}}) {
		if ($_ ne $id) {
			push(@a, $_);
		}
	}
	if (@a > 0) {
		print "$id\t$a[0]";
		for (my $i=1;$i<@a;$i++) {
			print ",$a[$i]";
		}
		print "\n";
	}
}
