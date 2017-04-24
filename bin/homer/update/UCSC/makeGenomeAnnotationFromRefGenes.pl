#!/usr/bin/perl -w


my %legit = ();
$legit{'ncRNA'} = 1;
$legit{'miscRNA'} = 1;
$legit{'pseudo'} = 1;
$legit{'other'} = 1;
$legit{'rRNA'} = 1;
$legit{'snoRNA'} = 1;
$legit{'scRNA'} = 1;


sub printCMD {
	print STDERR "\n\tUsage: makeGenomeAnnotationFromRefGenes.pl <refGenes.txt> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t(1st argument is the refGenes.txt file location! [required])\n";
	print STDERR "\t\t-prefix <name> (name of genome dataset, i.e. hg19 [required])\n";
	print STDERR "\t\t-org <name> (name of organism, i.e. human [required])\n";
	print STDERR "\t\t-tssStart <#> (start of promoter definition, def: -1000)\n";
	print STDERR "\t\t-tssEnd <#> (end of promoter definition, def: 100)\n";
	print STDERR "\t\t-ttsStart <#> (start of transcription termination definition, def: -100)\n";
	print STDERR "\t\t-ttsEnd <#> (end of transcription termination definition, def: 1000)\n";
	print STDERR "\t\t-knownGenes (flag specifying 1st argument is knownGenes.txt file instead of refGene.txt)\n";
	print STDERR "\n\t\tProduces files:\n";
	print STDERR "\t\t[base directory files]:\n";
	print STDERR "\t\t\t<prefix>.tss\n";
	print STDERR "\t\t\t<prefix>.tts\n";
	print STDERR "\t\t\t<prefix>.rna\n";
	print STDERR "\t\t\t<prefix>.splice5p\n";
	print STDERR "\t\t\t<prefix>.splice3p\n";
	print STDERR "\t\t\t<prefix>.aug\n";
	print STDERR "\t\t\t<prefix>.stop\n";
	print STDERR "\t\t\t<prefix>.miRNA\n";
	print STDERR "\t\t\t<prefix>.raw.annotation\n";
	print STDERR "\t\t\t<prefix>.intron.annotation\n";
	print STDERR "\t\t[annotation directory files]\n";
	print STDERR "\t\t\tpromoters.ann.txt\n";
	print STDERR "\t\t\ttts.ann.txt\n";
	print STDERR "\t\t\texons.ann.txt\n";
	print STDERR "\t\t\tcoding.ann.txt\n";
	print STDERR "\t\t\tutr5.ann.txt\n";
	print STDERR "\t\t\tutr3.ann.txt\n";
	print STDERR "\t\t\tintrons.ann.txt\n";
	print STDERR "\t\t\tncRNA.ann.txt\n";
	print STDERR "\t\t\tmiscRNA.ann.txt\n";
	print STDERR "\t\t\tmiRNA.ann.txt\n";
	print STDERR "\t\t\totherRNA.ann.txt\n";
	print STDERR "\t\t\tpseudo.ann.txt\n";
	print STDERR "\t\t\tscRNA.ann.txt\n";
	print STDERR "\t\t\tsnoRNA.ann.txt\n";
	print STDERR "\t\t\tsnRNA.ann.txt\n";
	print STDERR "\n";
	exit;
}
if (@ARGV < 1) {
	print STDERR "!!! Not enough command line arguments\n";
	printCMD();
}

my $refGeneFile = $ARGV[0];
my $refGeneFlag = 1;
my $organism = '';
my $prefix = '';
my $promoterStart = -1000;
my $promoterEnd = 100;
my $ttsStart = -100;
my $ttsEnd = 1000;
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-prefix') {
		$prefix = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-org') {
		$organism = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-knownGenes') {
		$refGeneFlag = 0;
	} elsif ($ARGV[$i] eq '-tssStart') {
		$promoterStart = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-tssEnd') {
		$promoterEnd = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-ttsStart') {
		$ttsStart = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-ttsEnd') {
		$ttsEnd = $ARGV[++$i];
	} else {
		printCMD();
	}
}

if ($organism eq '') {
	print STDERR "!!! Need to specify organism (-org)!\n";
	printCMD();
}
if ($prefix eq '') {
	print STDERR "!!! Need to specify dataset prefix (-prefix)\n";
	printCMD();
}

my $rand = rand();
my $tmpFile = $rand . ".tmp";
my $tmpFile2 = $rand . ".2.tmp";

if ($refGeneFlag == 0) {
	`cut -f1 "$refGeneFile" > "$tmpFile"`;
} else {
	`cut -f2 "$refGeneFile" > "$tmpFile"`;
}
`addGeneAnnotation.pl "$tmpFile" $organism > "$tmpFile2"`;
open IN, $tmpFile2;
my %type  = ();
my %typeArrays = ();
while (<IN>) {
	chomp;
	my @line = split /\t/;
	my $gtype = 'other';
	if (@line > 10) {
		$gtype = $line[10];
		if ($line[5] =~ /^mir/i || $line[6] =~ /mir-/i) {
			$gtype = 'miRNA';
		}
	}
	$type{$line[1]} = $gtype;
	if (!exists($typeArrays{$gtype})) {
		my @a = ();
		$typeArrays{$gtype}=\@a;
	}
}
close IN;
`rm $tmpFile $tmpFile2`;

my $upstreamLength = 5000;

my $rnaFlag = 0;

my @promoters = ();
my @tts = ();
my @exons = ();
my @introns = ();
my @upstream = ();
my @downstream = ();
my @UTR3 = ();
my @UTR5 = ();
my %uniqRNA = ();

open OUT, ">$prefix.rna";
open MIRNA, ">$prefix.miRNA";
open TSS, ">$prefix.tss";
open TTS, ">$prefix.tts";
open SPLICE5, ">$prefix.splice5p";
open SPLICE3, ">$prefix.splice3p";
open AUG, ">$prefix.aug";
open STOP, ">$prefix.stop";
open IN, $refGeneFile;
while (<IN>) {
	chomp;
	my @line= split /\t/;

	my $id = "";
	my $chr = "";
	my $dir = 0;
	my $gstart = 0;
	my $gend = 0;
	my $tstart = 0;
	my $tend = 0;
	my $numExons = 0;
	my $numIntrons = 0;
	my @starts = ();
	my @ends = ();
	my $gtype = "";

	if ($refGeneFlag == 0) {
		$id = $line[0];
		$chr = $line[1];
		$dir = $line[2];
		$gstart = $line[3]+1;
		$gend = $line[4];
		$tstart = $line[5]+1;
		$tend = $line[6];
		$numExons = $line[7];
		$numIntrons = $numExons-1;
		@starts = split /\,/,$line[8];
		foreach(@starts) {
			$_+=1;
		}
		@ends = split /\,/,$line[9];
	} else {
		$id = $line[1];
		$chr = $line[2];
		$dir = $line[3];
		$gstart = $line[4]+1;
		$gend = $line[5];
		$tstart = $line[6]+1;
		$tend = $line[7];
		$numExons = $line[8];
		$numIntrons = $numExons-1;
		@starts = split /\,/,$line[9];
		foreach(@starts) {
			$_+=1;
		}
		@ends = split /\,/,$line[10];
	}

	my $ncFlag = 0;
	if ($tstart > $tend) {
		$ncFlag = 1;
	}

	my $rtype = "ncRNA";
	if (exists($type{$id})) {
		$gtype = $type{$id};
		$rtype = $gtype;
		if ($rtype eq 'protein-coding') {
			$rtype = 'E';
		}
	}

	my $oid = $id;
	if (exists($uniqRNA{$id})) {
		$oid .= "-HOMER" . $uniqRNA{$id}++;
	} else {
		$uniqRNA{$id} = 2;
	}


	my @lens = ();
	for (my $i=0;$i<@starts;$i++) {
		push(@lens, $ends[$i]-$starts[$i]);
	}

	my @data = ();

	my $exonCount=1;

	if ($dir eq '+' || $dir eq '0') {
		my $tss = $gstart;
		my $tts = $gend;

		my $pstart = $tss+$promoterStart;
		my $pend = $tss+$promoterEnd;
		my $r = {c=>$chr,s=>$pstart,e=>$pend,d=>0, t=>'P', n=>"promoter-TSS ($id)"};
		push(@promoters, $r);

		$pstart = $tss-2000;
		$pend = $tss+2000;
		print TSS "$oid\t$chr\t$pstart\t$pend\t0\n";
		$pstart = $tts-2000;
		$pend = $tts+2000;
		print TTS "$oid\t$chr\t$pstart\t$pend\t0\n";

		$pstart = $tstart-100;
		$pend = $tstart+100;
		print AUG "$oid\t$chr\t$pstart\t$pend\t0\n";
		$pstart = $tend-100;
		$pend = $tend+100;
		print STOP "$oid\t$chr\t$pstart\t$pend\t0\n";


		$pstart = $tts+$ttsStart;
		$pend = $tts+$ttsEnd;
		$r = {c=>$chr,s=>$pstart,e=>$pend,d=>0, t=>'TTS', n=>"TTS ($id)"};
		push(@tts, $r);

		$pstart = $tss-$upstreamLength;
		$pend = $tss;
		$r = {c=>$chr,s=>$pstart,e=>$pend,d=>0, t=>'P', n=>"upstream ($id)"};
		push(@upstream, $r);

		$pstart = $tts;
		$pend = $tts+$upstreamLength;
		$r = {c=>$chr,s=>$pstart,e=>$pend,d=>0, t=>'P', n=>"downstream ($id)"};
		push(@downstream, $r);

		for (my $i=0;$i<$numExons;$i++) {
			my $s = $starts[$i];
			my $e = $ends[$i];
			my $IDnum = $i+1;
			if ($i < $numExons-1) {
				my $oStart = $e-100;
				my $oEnd = $e+100;
				print SPLICE5 "$id-$IDnum\t$chr\t$oStart\t$oEnd\t0\n";
			}
			if ($i>0) {
				my $oStart = $s-100;
				my $oEnd = $s+100;
				print SPLICE3 "$id-$IDnum\t$chr\t$oStart\t$oEnd\t0\n";
			}

	
			if ($tstart >= $s && $tstart <= $e) {
				my $s1 = $s;
				my $e1 = $tstart-1;
				my $r = {c=>$chr,s=>$s1,e=>$e1,d=>0, t=>'5UTR', n=>"5' UTR ($id, exon $exonCount of $numExons)"};
				if ($ncFlag) {
					push(@data, "E$exonCount" . ":$s1");
					$r = {c=>$chr,s=>$s1,e=>$e1,d=>0, t=>$rtype, n=>"non-coding ($id, exon $exonCount of $numExons)"};
					push(@exons, $r);
				} else {
					push(@data, "E$exonCount" . "_5UTR:$s1");
					push(@UTR5, $r);
				}
				push(@{$typeArrays{$gtype}}, $r) if ($gtype ne '');
				$s = $e1+1;
			}
			if ($tend >= $s && $tend+1 < $e) {
				my $s1 = $s;
				my $e1 = $tend;
				my $r = {c=>$chr,s=>$s1,e=>$e1,d=>0, t=>'E', n=>"exon ($id, exon $exonCount of $numExons)"};
				push(@exons, $r);
				push(@{$typeArrays{$gtype}}, $r) if ($gtype ne '');
				push(@data, "E$exonCount" . ":$s1");
				$s = $e1+1;
				$s1 = $s;
				$e1= $e;
				$r = {c=>$chr,s=>$s1,e=>$e1,d=>0, t=>'3UTR', n=>"3' UTR ($id, exon $exonCount of $numExons)"};
				push(@UTR3, $r);
				push(@{$typeArrays{$gtype}}, $r) if ($gtype ne '');
				push(@data, "E$exonCount" . "_3UTR:$s1");
				$s = $e1;
			} else {
				$r = {c=>$chr,s=>$s,e=>$e,d=>0, t=>'E', n=>"exon ($id, exon $exonCount of $numExons)"};
				if ($ncFlag) {
					$r = {c=>$chr,s=>$s,e=>$e,d=>0, t=>$rtype, n=>"non-coding ($id, exon $exonCount of $numExons)"};
					push(@exons, $r);
					push(@data, "E$exonCount" . ":$s");
				} else {
					if ($e < $tstart) {
						$r = {c=>$chr,s=>$s,e=>$e,d=>0, t=>'5UTR', n=>"5' UTR ($id, exon $exonCount of $numExons)"};
						push(@data, "E$exonCount" . "_5UTR:$s");
						push(@UTR5, $r);
					} elsif ($s > $tend) {
						$r = {c=>$chr,s=>$s,e=>$e,d=>0, t=>'3UTR', n=>"3' UTR ($id, exon $exonCount of $numExons)"};
						push(@UTR3, $r);
						push(@data, "E$exonCount" . "_3UTR:$s");
					} else {
						push(@exons, $r);
						push(@data, "E$exonCount" . ":$s");
					}
				}
				push(@{$typeArrays{$gtype}}, $r) if ($gtype ne '');
			}

			if ($i<$numExons-1) {
				my $s = $ends[$i]+1;
				my $e = $starts[$i+1]-1;
				my $r = {c=>$chr,s=>$s,e=>$e,d=>0,t=>'I',n=>"intron ($id, intron $exonCount of $numIntrons)"};
				push(@introns, $r);
				push(@data, "I$exonCount" . ":$s");
			}
			$exonCount++;
		}


	} else {

		my $tss = $gend;
		my $tts = $gstart;

		my $pstart = $tss-$promoterStart;
		my $pend = $tss-$promoterEnd;
		my $r = {c=>$chr,s=>$pend,e=>$pstart,d=>1,t=>'P',n=>"promoter-TSS ($id)"};
		push(@promoters, $r);

		$pstart = $tss-2000;
		$pend = $tss+2000;
		print TSS "$oid\t$chr\t$pstart\t$pend\t1\n";
		$pstart = $tts-2000;
		$pend = $tts+2000;
		print TTS "$oid\t$chr\t$pstart\t$pend\t1\n";


		$pstart = $tend-100;
		$pend = $tend+100;
		print AUG "$oid\t$chr\t$pstart\t$pend\t1\n";
		$pstart = $tstart-100;
		$pend = $tstart+100;
		print STOP "$oid\t$chr\t$pstart\t$pend\t1\n";

		$pstart = $tts-$ttsStart;
		$pend = $tts-$ttsEnd;
		$r = {c=>$chr,s=>$pend,e=>$pstart,d=>1, t=>'TTS', n=>"TTS ($id)"};
		push(@tts, $r);

		$pstart = $tss;
		$pend = $tss+$upstreamLength;
		$r = {c=>$chr,s=>$pstart,e=>$pend,d=>1, t=>'P', n=>"upstream ($id)"};
		push(@upstream, $r);

		$pstart = $tts-$upstreamLength;
		$pend = $tts;
		$r = {c=>$chr,s=>$pstart,e=>$pend,d=>1, t=>'P', n=>"downstream ($id)"};
		push(@downstream, $r);

		for (my $i=$numExons-1;$i>=0;$i--) {
			my $s = $starts[$i];
			my $e = $ends[$i];

			my $IDnum = $numExons-$i;
			if ($i < $numExons-1) {
				my $oStart = $e-100;
				my $oEnd = $e+100;
				print SPLICE3 "$id-$IDnum\t$chr\t$oStart\t$oEnd\t1\n";
			}
			if ($i>0) {
				my $oStart = $s-100;
				my $oEnd = $s+100;
				print SPLICE5 "$id-$IDnum\t$chr\t$oStart\t$oEnd\t1\n";
			}

			if ($tend >= $s && $tend < $e) {
				my $s1 = $tend+1;
				my $e1 = $e;
				my $r = {c=>$chr,s=>$s1,e=>$e1,d=>1, t=>'5UTR', n=>"5' UTR ($id, exon $exonCount of $numExons)"};
				if ($ncFlag) {
					$r = {c=>$chr,s=>$s1,e=>$e1,d=>1, t=>$rtype, n=>"non-coding ($id, exon $exonCount of $numExons)"};
					push(@data, "E$exonCount" . ":$s1");
					push(@exons, $r);
				} else {
					push(@data, "E$exonCount" . "_5UTR:$s1");
					push(@UTR5, $r);
				}
				push(@{$typeArrays{$gtype}}, $r) if ($gtype ne '');
				$e = $s1-1;
			}

			if ($tstart >= $s && $tstart <= $e) {
				my $s1 = $tstart;
				my $e1 = $e;
				my $r = {c=>$chr,s=>$s1,e=>$e1,d=>1, t=>'E', n=>"exon ($id, exon $exonCount of $numExons)"};
				push(@exons, $r);
				push(@{$typeArrays{$gtype}}, $r) if ($gtype ne '');
				push(@data, "E$exonCount" . ":$s1");
				$e = $s1-1;
				$s1 = $s;
				$e1= $e;
				$r = {c=>$chr,s=>$s1,e=>$e1,d=>1, t=>'3UTR', n=>"3' UTR ($id, exon $exonCount of $numExons)"};
				push(@UTR3, $r);
				push(@{$typeArrays{$gtype}}, $r) if ($gtype ne '');
				push(@data, "E$exonCount" . "_3UTR:$s1");
				$e = $s1-1;
			} else {
				$r = {c=>$chr,s=>$s,e=>$e,d=>1, t=>'E', n=>"exon ($id, exon $exonCount of $numExons)"};
				if ($ncFlag) {
					$r = {c=>$chr,s=>$s,e=>$e,d=>1, t=>$rtype, n=>"non-coding ($id, exon $exonCount of $numExons)"};
					push(@exons, $r);
					push(@data, "E$exonCount" . ":$s");
				} else {
					if ($e < $tstart) {
						$r = {c=>$chr,s=>$s,e=>$e,d=>1, t=>'3UTR', n=>"3' UTR ($id, exon $exonCount of $numExons)"};
						push(@UTR3, $r);
						push(@data, "E$exonCount" . "_3UTR:$s");
					} elsif ($s > $tend) {
						$r = {c=>$chr,s=>$s,e=>$e,d=>1, t=>'5UTR', n=>"5' UTR ($id, exon $exonCount of $numExons)"};
						push(@UTR5, $r);
						push(@data, "E$exonCount" . "_5UTR:$s");
					} else {
						push(@exons, $r);
						push(@data, "E$exonCount" . ":$s");
					}
				}
				push(@{$typeArrays{$gtype}}, $r) if ($gtype ne '');
			}


			if ($i>0) {
				my $s = $ends[$i-1]+1;
				my $e = $starts[$i]-1;
				my $r = {c=>$chr,s=>$s,e=>$e,d=>1,t=>'I',n=>"intron ($id, intron $exonCount of $numIntrons)"};
				push(@introns, $r);
				push(@data, "I$exonCount" . ":$s");
			}
			$exonCount++;
		}
	}

	print OUT "$oid\t$chr\t$gstart\t$gend\t$dir\t";
	print MIRNA "$oid\t$chr\t$gstart\t$gend\t$dir\t" if ($gtype eq 'miRNA');
	my $c = 0;
	foreach(@data) {
		$c++;
		print OUT "," if ($c > 1);
		print OUT $_;
		if ($gtype eq 'miRNA') {
			print MIRNA "," if ($c > 1);
			print MIRNA $_;
		}
			
	}
	print OUT "\n";
	print MIRNA "\n" if ($gtype eq 'miRNA');

}
close IN;



my %uniqueIDs = ();
my @all = ();
push(@all, @promoters, @tts, @exons, @UTR5, @UTR3);
open ANN, ">$prefix.raw.annotation";
foreach(@all) {
	my $name = $_->{'n'};
	if (exists($uniqueIDs{$name})) {
		$name .= "." . $uniqueIDs{$name}++;
	} else {
		$uniqueIDs{$name} = 2;
	}
	print ANN "$name\t$_->{'c'}\t$_->{'s'}\t$_->{'e'}\t$_->{'d'}\t$_->{'t'}\n";
}
close ANN;

open ANN, ">$prefix.intron.annotation";
foreach(@introns) {
	my $name = $_->{'n'};
	if (exists($uniqueIDs{$name})) {
		$name .= "." . $uniqueIDs{$name}++;
	} else {
		$uniqueIDs{$name} = 2;
	}
	print ANN "$name\t$_->{'c'}\t$_->{'s'}\t$_->{'e'}\t$_->{'d'}\t$_->{'t'}\n";
}
close ANN;

printAnnFile(\@promoters, "promoters.ann.txt");
printAnnFile(\@tts, "tts.ann.txt");
my @EXON = ();
push(@EXON, @exons, @UTR5, @UTR3);
printAnnFile(\@EXON, "exons.ann.txt");
printAnnFile(\@exons, "coding.ann.txt");
printAnnFile(\@UTR5, "utr5.ann.txt");
printAnnFile(\@UTR3, "utr3.ann.txt");
printAnnFile(\@introns, "introns.ann.txt");
foreach(keys %typeArrays) {
	my $t = $_;
	printAnnFile($typeArrays{$t},"$t.ann.txt");
}

exit;

sub printAnnFile {
	my ($anns, $file) = @_;
	open OUT, ">$file";
	my %uniq = ();
	foreach(@$anns) {
		my $name = $_->{'n'};
		if (exists($uniq{$name})) {
			$name .= "." . $uniq{$name}++;
		} else {
			$uniq{$name} = 2;
		}
		print OUT "$name\t$_->{'c'}\t$_->{'s'}\t$_->{'e'}\t$_->{'d'}\t$_->{'t'}\n";
	}
	close OUT;
}
