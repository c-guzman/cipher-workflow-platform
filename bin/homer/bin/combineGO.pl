#!/usr/bin/perl -w


sub printCMD {
	print STDERR "\n\tcombineGO.pl [options] -f <filename i.e. biological_process.txt> -d <directory1> [directory2] ...\n";
	print STDERR "\n\toptions:\n";
	print STDERR "\t\t-f <filename> (Filename to join enrichment values from)\n";
	print STDERR "\t\t-d <directory1> [directory2] ... (findMotifs.pl output directories to join)\n";
	print STDERR "\t\t-top <#> (Only keep top # terms per directory, default: keep all)\n";
	print STDERR "\t\t-column <#> (Column in files to join, default: 4)\n";
	print STDERR "\n";
	exit;
}
if (@ARGV < 3) {
	printCMD();
}

#my $col = 10;
my $col = 4;
my $top = 0;
my $filename= $ARGV[0];
my @tagDirs = ();

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-f') {
		$filename=$ARGV[++$i];
	} elsif ($ARGV[$i] eq '-d') {
        my $bail = 0;
        while ($ARGV[++$i] !~ /^\-/) {
            push(@tagDirs, $ARGV[$i]);
            print STDERR "\t\t$ARGV[$i]\n";
            if ($i>=@ARGV-1) {
                $bail = 1;
                last;
            }
        }
        last if ($bail == 1);
        $i--;
	} elsif ($ARGV[$i] eq '-top') {
		$top = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-column') {
		$col = $ARGV[++$i];
	} else {
		print STDERR "!! Can't recognize $ARGV[$i]\n";
		printCMD();
	}
}


my %top = ();


print "TermID\tTermDesc";
for (my $i=0;$i<@tagDirs;$i++) {
	my $dir = $tagDirs[$i];
	open IN, "$dir/$filename" or die "Could not open $dir/$filename!!!\n";

	print "\t$dir";

	my $c=0;
	my %d = ();
	while (<IN>) {
		$c++;
		next if ($c==1);
		chomp;
		s/\r//g;
		my @line = split /\t/;
		my $id = $line[0] . "\t" . $line[1];
		$d{$id} = $line[$col-1];
	}
	close IN;

	my @a = sort {$d{$a} <=> $d{$b}} keys %d;
	if ($top > 0) {
		for (my $j=0;$j<$top;$j++) {
			$top{$a[$j]}=1;
		}
	} else {
		foreach(keys %d) {
			$top{$_}=1;
		}
	}
	push(@data, \%d);
}
print "\n";
foreach(keys %top) {
	my $id = $_;
	print "$_";
	foreach(@data) {
		my $v = 0;
		if (exists($_->{$id})) {
			$v = $_->{$id};
		}
		print "\t$v";
	}
	print "\n";
}

