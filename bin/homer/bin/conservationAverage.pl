#!/usr/bin/env perl




if (@ARGV < 2) {
	print STDERR "<cons seq file> <offset>\n";
	exit;
}
my $offset = $ARGV[1];
my $N = ();
my @avg= ();
open IN, $ARGV[0];
my $count = 0;
print STDERR "\tCalculating Average Conservation...\n";
while (<IN>) {
	$count++;
print STDERR "\t$count\n" if ($count % 1000 == 0);
	chomp;
	s/\r//g;
	my @line = split /\t/;
	my @c = split undef, $line[1];
	foreach(@c) {
		$_ /= 9;
	}
	$N++;
	if ($count == 1) {
		@avg = @c;
	} else {
		for (my $i=0;$i<@c;$i++) {
			$avg[$i] += $c[$i];
		}
	}
}

print "Distance\tAvg\n";
for (my $i=0;$i<@avg;$i++) {
	$avg[$i] /= $N;
	my $p = $offset + $i;
	print "$p\t$avg[$i]\n";
}
