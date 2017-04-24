#!/usr/bin/env perl
use warnings;
use lib "/home/carlos/tools/homer/.//bin";
my $homeDir = "/home/carlos/tools/homer/./";

package HomerConfig;


sub loadConfigFile {

	my ($file) = @_;
	if (!defined($file)) {
		$file = $homeDir . "/config.txt";
	}
	my %a = ();
	my %b = ();
	my %c = ();
	my %d = ();
	my %e = ();
	my $config = {PROMOTERS=>\%a, GENOMES=>\%b, ORGANISMS=>\%c,SETTINGS=>\%d,SOFTWARE=>\%e};

	parseConfigFile($file,$config,$homeDir);

	my $localConfig = $ENV{"HOME"} . "/.homerConfig.txt";
	if (-e $localConfig) {
		#parseConfigFile($localConfig,$config,"");
	}
	return $config;
}
sub printConfigFile {
	my ($config,$outFile, $updateFlag) = @_;

	my $file = $homeDir . "/config.txt";
	if (defined($outFile)) {
		$file = $outFile;
		#print STDERR "\tPrinting Configuration to $outFile (may not be standard...)\n";
	}

	open CONFIG, ">$file";

	if (defined($updateFlag)) {
		print CONFIG "# Homer Update File\n";
	} else {
		print CONFIG "# Homer Configuration File (automatically generated)\n";
	}
	print CONFIG "#\n";
	print CONFIG "# This file is updated from the Homer website and contains information about data available for\n";
	print CONFIG "# use with the program.\n";
	print CONFIG "#\n";
	print CONFIG "# This file has been further modified by the local update script\n";
	print CONFIG "#\n";
	print CONFIG "# Each section has the same format, which is <tab> separated values specifying:\n";
	print CONFIG "# package name <tab> version <tab> description <tab> url <tab> optional parameters (, separated)\n";
	print CONFIG "#\n";

	my @groups = ("SOFTWARE","ORGANISMS","PROMOTERS","GENOMES");
	foreach(@groups) {
		my $group = $_;
		print CONFIG "\n$group\n";
		foreach(keys %{$config->{$group}}) {
			my $package = $_;
			my $version = $config->{$group}->{$package}->{'version'};
			my $desc = $config->{$group}->{$package}->{'desc'};
			my $url = $config->{$group}->{$package}->{'url'};
			my $location = $config->{$group}->{$package}->{'location'};
			my $params = "";
			my $c = 0;
			foreach(@{$config->{$group}->{$package}->{'parameters'}}) {
				$params .= "," if ($c > 0);
				$c++;
				$params .= $_;
			}
			print CONFIG "$package\t$version\t$desc\t$url\t$location\t$params\n";
		}
	}
	if (defined($updateFlag)) {
	} else {
		print CONFIG "SETTINGS\n";
		if (exists($config->{'SETTINGS'})) {
			foreach(keys %{$config->{'SETTINGS'}}) {
				my $var = $_;
				my $val = $config->{'SETTINGS'}->{$var};
				print CONFIG "$var=$val\n";
			}
		}
	}
	close CONFIG;
}
sub parseConfigFile {
	my ($file, $config,$homeDir) = @_;
	my $mode = '';	
	open IN, $file or die "Could not open configuration file ($file)\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		s/#.*$//;
		s/^\s+//;
		next if ($_ eq '');

		if (/^SOFTWARE/) {
			$mode = 'SOFTWARE';
			next;
		}
		if (/^ORGANISMS/) {
			$mode = 'ORGANISMS';
			next;
		}
		if (/^PROMOTERS/) {
			$mode = 'PROMOTERS';
			next;
		}
		if (/^GENOMES/) {
			$mode = 'GENOMES';
			next;
		}
		if (/^SETTINGS/) {
			$mode = 'SETTINGS';
			next;
		}
		#if (/^PREPARSED/) {
		#	$mode = 'PREPARSED';
		#	next;
		#}
		my @line = split /\t/;
		next if ($mode eq '');

		if (@line > 4) {
			my $package = $line[0];
			my $version = $line[1];
			my $description = $line[2];
			my $url = $line[3];
			my $location = $line[4];
			my @params = ();
			if (@line > 5) {
				@params = split /\,/, $line[5];
			}
			my $d = $homeDir . "/" . $location;
			if ($mode eq 'GENOMES') {
				$config->{'GENOMES'}->{$package} = {org=>$params[0],promoters=>$params[1], directory=>$d,
							location=>$location, version=>$version,url=>$url,desc=>$description, parameters=>\@params};
			} elsif ($mode eq 'PROMOTERS') {
				$config->{'PROMOTERS'}->{$package} = {org=>$params[0], directory=>$d, genome=>$params[1], 
							idtype=>$params[2], start=>$params[3], end=>$params[4],
							location=>$location, version=>$version,url=>$url,desc=>$description, parameters=>\@params};
			} elsif ($mode eq 'ORGANISMS') {
				my $taxid = 'NA';
				my $source = 'NA';
				$taxid = $params[0] if (@params > 0);
				$source = $params[1] if (@params > 1);
				$config->{'ORGANISMS'}->{$package} = {taxid=>$taxid,source=>$source, directory=>$d,
							location=>$location, version=>$version,url=>$url,desc=>$description, parameters=>\@params};
			} elsif ($mode eq 'SOFTWARE') {
				$config->{'SOFTWARE'}->{$package} = {
							location=>$location, version=>$version,url=>$url,desc=>$description, parameters=>\@params};
			#} elsif ($mode eq 'PREPARSED') {
			#	$config->{'PREPARSED'}->{$package} = {genome=>$params[0],size=>$params[1], directory=>$d,
			#				location=>$location, version=>$version,url=>$url,desc=>$description, parameters=>\@params};
			}
		} elsif ($mode eq 'SETTINGS') {
			my @var = split /\=/, $line[0];
			next if (@var < 2);
			$config->{'SETTINGS'}->{$var[0]} = $var[1];
		}
	}
	return $config;
}

sub readMotifFile {
	my ($file) = @_;
	my @motifs = ();
	my $index = -1;
	open MOTIFFILE, $file;
	while (<MOTIFFILE>) {
		my $og = $_;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($line[0] =~ /^>/) {
			my @a = ();
			my $m = {name=>$line[1],consensus=>$line[0],len=>0,matrix=>\@a,text=>$og,threshold=>$line[2]};
			push(@motifs, $m);
			$index++;
		} else {
			$motifs[$index]->{'len'}++;
			$motifs[$index]->{'text'} .= $og;
			push(@{$motifs[$index]->{'matrix'}}, \@line);
		}
	}
	close MOTIFFILE;
	return \@motifs;
}

sub readTagInfo {
    my ($dir,$pc) = @_;

	my $file = "$dir/" . "/tagInfo.txt";
    my $p=0;
    my $t = 0;
    my $fragEst = "NA";
    my $peakEst = "NA";
	if (!defined($pc)) {
		$pc = 0;
	}
    open IN, $file or die "Could not open $file\n";
    my $count = 0;
    while (<IN>) {
        $count++;
        if ($count == 1) {
            if (/Estimated Fragment Length=(\d+?)[,\t]/) {
                $fragEst = $1;
            }
            if (/Estimated Peak Width=(\d+?)[,\t]/) {
                $peakEst = $1;
            }

        }
        next if ($count < 2);
        chomp;
        s/\r//g;
        my @line= split /\t/;
		if ($count == 2 || $line[0] eq 'genome') {
        	$p = $line[1];
        	$t = $line[2];
		} else {
			if ($line[0] =~ /fragmentLengthEstimate=(.+)$/) {
				$fragEst = $1;
			}
			if ($line[0] =~ /peakSizeEstimate=(.+)$/) {
				$peakEst = $1;
			}
		}
    }
    close IN;
	if ($pc == 1) {	
		$t = $p;
	} elsif ($pc > 1) {
		$file = "$dir/tagCountDistribution.txt";
		if (-e $file) {
			open IN, $file;	
			my $total = 0;
			my $count = 0;
			while (<IN>) {
				$count++;
				chomp;
				s/\r//g;
				my @line = split /\t/;
				next if ($count==1);
				$total += $line[1];
				if ($line[0] eq $pc) {
					last;
				}
			}
			close IN;
			#print STDERR "total=$total\n";
			$t *= $total;
		}
	}

    #print STDERR "$t, $p, $fragEst, $peakEst, $pc\n";
    return ($t, $p, $fragEst, $peakEst);
}

sub getHiCBgRes {
	my ($dir,$currentRes,$cpus) = @_;
	my $tmpFile = rand() . ".tmp";
	`ls -1 "$dir"/HiCbackground_* > "$tmpFile"`;
	open TMPIN, $tmpFile;
	my %res = ();
	while (<TMPIN>) {
		chomp;
		if (/HiCbackground_(\d+)_bp/) {
			$res{$1}=1;
		}
	}
	close TMPIN;
	`rm "$tmpFile"`;
	if (!exists($res{$currentRes})) {
		print STDERR "\t! Background for resolution $currentRes not found for Directory $dir\n";
		print STDERR "\t\tExisting Bg models:\n";
		foreach(keys %res) {
			print STDERR "\t\t\t$_\n";
		}
		print STDERR "\tWill create model for $currentRes in 5 seconds...\n";
		`sleep 5`;
		`analyzeHiC "$dir" -bgonly -res $currentRes -cpu $cpus`;
	}
	return \%res;

}

sub revopp {
	my ($seq) = @_;
	my $rv = reverse($seq);
	$rv =~ s/A/X/g;
	$rv =~ s/T/A/g;
	$rv =~ s/X/T/g;
	$rv =~ s/C/X/g;
	$rv =~ s/G/C/g;
	$rv =~ s/X/G/g;
	return $rv;
}


sub parseCustomGenome {
	my ($genome) = @_;

	my $genomeDir = "";
	my $genomeName = "";
	my $genomeParseDir = "";
	if (-e $genome) {
		$customGenome=1;
		$genomeDir = $genome;
		$genomeParseDir = $genome;
		if (-f $genomeParseDir) {
			if ($genomeParseDir =~ s/\/([^\/]*)$//) {
				$genomeName = $1;
			} else {
				$genomeName = $genome;
				$genomeParseDir = ".";
			}
		} else {
			my $a = $genomeParseDir;
			$a  =~ s/\/+$//;
			$a  =~ s/\/([^\/]*)$//;
			$genomeName = $1;
		}
		#print STDERR "$genomeName\n";
		$genomeParseDir .= "/preparsed/";
		print STDERR "\tUsing Custom Genome\n";
	} else {
		print STDERR "\n!!!!Genome $genome not found in $homeDir/config.txt\n\n";
		print STDERR "\tTo check if is available, run \"perl $homeDir/configureHomer.pl -list\"\n";
		print STDERR "\tIf so, add it by typing \"perl $homeDir/configureHomer.pl -install $genome\"\n";
		print STDERR "\n";
		exit;
	}
	return ($genomeName,$genomeDir,$genomeParseDir);
}

sub parseUCSCStr {
	my ($str) = @_;
	my $chr = "";
	my $start = 0;
	my $end = 1e10;
	$str =~ s/\,//g;
	$str =~ /(.*?)\:(\d+)\-(\d+)/;
	$chr = $1;
	$start = $2;
	$end = $3;
	return ($chr,$start,$end);

}

sub checkMSet {
	my ($mset, $org) = @_;
	if ($mset eq 'auto' && $org ne '' && $org ne 'null') {
		my %table = ();
		open ORGTABLE, "$homeDir/data/knownTFs/organism.table.txt" or die "Where is the file?\n";
		while (<ORGTABLE>) {
			next if (/^#/);
			chomp;
			s/\r//g;
			my @line = split /\t/;
			$table{$line[0]} = $line[1];
		}
		close ORGTABLE;
		if (exists($table{$org})) {
			$mset = $table{$org};
			print STDERR "\tFound mset for \"$org\", will check against $mset motifs\n";
		} else {
			print STDERR "\tCould not find mset for \"$org\", will check against all motifs\n";
			$mset = 'all';
		}
	} elsif ($mset eq 'auto') {
		$mset='all';
	}
	my $msetDir = $homeDir . "/data/knownTFs/$mset/";
	my $allMotifsFile = $msetDir . "all.motifs";
	my $knownMotifsFile = $msetDir . "known.motifs";
	unless (-e $allMotifsFile) {
		print STDERR "\tWarning, couldn't find $allMotifsFile for mset: \"$mset\"\n";
		print STDERR "\t\tUsing all motifs...\n";
		$mset = 'all';
		$msetDir = $homeDir . "/data/knownTFs/$mset/";
		$allMotifsFile = $msetDir . "all.motifs";
		$knownMotifsFile = $msetDir . "known.motifs";
	}
	#print STDERR "$allMotifsFile, $knownMotifsFile\n";
	return ($allMotifsFile, $knownMotifsFile);

}


1;
