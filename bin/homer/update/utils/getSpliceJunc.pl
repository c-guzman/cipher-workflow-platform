#!/usr/bin/perl -w

if (@ARGV < 2) {
	print STDERR "\n\t<refGenes.txt> <mrna.fa>\n";
	print STDERR "\n";
	exit;
}

my $refGeneFlag = 1;
if (@ARGV > 2) {
	$refGeneFlag = 0;
}

my $prefix = $ARGV[1];

my $upstreamLength = 5000;

my $promoterStart = -1000;
my $promoterEnd = 100;
my $ttsStart = -100;
my $ttsEnd = 1000;
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
open TSS, ">$prefix.tss";
open TTS, ">$prefix.tts";
open SPLICE5, ">$prefix.splice5p";
open SPLICE3, ">$prefix.splice3p";
open AUG, ">$prefix.aug";
open STOP, ">$prefix.stop";
open ANN, ">$prefix.raw.annotation";
open IN, $ARGV[0];
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

	if ($refGeneFlag == 0) {
		$id = $line[0];
		$chr = $line[1];
		$dir = $line[2];
		$gstart = $line[3];
		$gend = $line[4];
		$tstart = $line[5];
		$tend = $line[6];
		$numExons = $line[7];
		$numIntrons = $numExons-1;
		@starts = split /\,/,$line[8];
		@ends = split /\,/,$line[9];
	} else {
		$id = $line[1];
		$chr = $line[2];
		$dir = $line[3];
		$gstart = $line[4];
		$gend = $line[5];
		$tstart = $line[6];
		$tend = $line[7];
		$numExons = $line[8];
		$numIntrons = $numExons-1;
		@starts = split /\,/,$line[9];
		@ends = split /\,/,$line[10];
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
		my $tss = $gstart+1;
		my $tts = $gend+1;

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

		$pstart = $tstart+1-100;
		$pend = $tstart+1+100;
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
				my $oStart = $e+1-100;
				my $oEnd = $e+1+100;
				print SPLICE5 "$id-$IDnum\t$chr\t$oStart\t$oEnd\t0\n";
			}
			if ($i>0) {
				my $oStart = $s+1-100;
				my $oEnd = $s+1+100;
				print SPLICE3 "$id-$IDnum\t$chr\t$oStart\t$oEnd\t0\n";
			}

	
			if ($tstart >= $s && $tstart <= $e) {
				my $s1 = $s;
				my $e1 = $tstart;
				my $r = {c=>$chr,s=>$s1,e=>$e1,d=>0, t=>'E', n=>"5' UTR ($id, exon $exonCount of $numExons)"};
				push(@UTR5, $r);
				push(@data, "E$exonCount" . "_5UTR:$s1");
				$s = $e1;
			}
			if ($tend >= $s && $tend <= $e) {
				my $s1 = $s;
				my $e1 = $tend;
				my $r = {c=>$chr,s=>$s1,e=>$e1,d=>0, t=>'E', n=>"exon ($id, exon $exonCount of $numExons)"};
				push(@exons, $r);
				push(@data, "E$exonCount" . ":$s1");
				$s = $e1;
				$s1 = $s;
				$e1= $e;
				$r = {c=>$chr,s=>$s1,e=>$e1,d=>0, t=>'E', n=>"3' UTR ($id, exon $exonCount of $numExons)"};
				push(@UTR3, $r);
				push(@data, "E$exonCount" . "_3UTR:$s1");
				$s = $e1;
			} else {
				$r = {c=>$chr,s=>$s,e=>$e,d=>0, t=>'E', n=>"exon ($id, exon $exonCount of $numExons)"};
				if ($e < $tstart) {
					$r = {c=>$chr,s=>$s,e=>$e,d=>0, t=>'E', n=>"5' UTR ($id, exon $exonCount of $numExons)"};
					push(@UTR5, $r);
					push(@data, "E$exonCount" . "_5UTR:$s");
				} elsif ($s > $tend) {
					$r = {c=>$chr,s=>$s,e=>$e,d=>0, t=>'E', n=>"3' UTR ($id, exon $exonCount of $numExons)"};
					push(@UTR3, $r);
					push(@data, "E$exonCount" . "_3UTR:$s");
				} else {
					push(@exons, $r);
					push(@data, "E$exonCount" . ":$s");
				}
			}

			if ($i<$numExons-1) {
				my $s = $ends[$i];
				my $e = $starts[$i+1];
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
		$pstart = $tstart+1-100;
		$pend = $tstart+1+100;
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

			if ($tend >= $s && $tend <= $e) {
				my $s1 = $tend;
				my $e1 = $e;
				my $r = {c=>$chr,s=>$s1,e=>$e1,d=>1, t=>'5UTR', n=>"5' UTR ($id, exon $exonCount of $numExons)"};
				push(@UTR5, $r);
				push(@data, "E$exonCount" . "_5UTR:$s1");
				$e = $s1;
			}

			if ($tstart >= $s && $tstart <= $e) {
				my $s1 = $tstart;
				my $e1 = $e;
				my $r = {c=>$chr,s=>$s1,e=>$e1,d=>1, t=>'E', n=>"exon ($id, exon $exonCount of $numExons)"};
				push(@exons, $r);
				push(@data, "E$exonCount" . ":$s1");
				$e = $s1;
				$s1 = $s;
				$e1= $e;
				$r = {c=>$chr,s=>$s1,e=>$e1,d=>1, t=>'3UTR', n=>"3' UTR ($id, exon $exonCount of $numExons)"};
				push(@UTR3, $r);
				push(@data, "E$exonCount" . "_3UTR:$s1");
				$e = $s1;
			} else {
				$r = {c=>$chr,s=>$s,e=>$e,d=>1, t=>'E', n=>"exon ($id, exon $exonCount of $numExons)"};
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


			if ($i>0) {
				my $s = $ends[$i-1];
				my $e = $starts[$i];
				my $r = {c=>$chr,s=>$s,e=>$e,d=>1,t=>'I',n=>"intron ($id, intron $exonCount of $numIntrons)"};
				push(@introns, $r);
				push(@data, "I$exonCount" . ":$s");
			}
			$exonCount++;
		}
	}

	print OUT "$oid\t$chr\t$gstart\t$gend\t$dir\t";
	my $c = 0;
	foreach(@data) {
		$c++;
		print OUT "," if ($c > 1);
		print OUT $_;
	}
	print OUT "\n";

}
close IN;



my %uniqueIDs = ();
my @all = ();
push(@all, @promoters, @tts, @exons, @UTR5, @UTR3, @introns);
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

printAnnFile(\@promoters, "promoters.ann.txt");
printAnnFile(\@tts, "tts.ann.txt");
my @EXON = ();
push(@EXON, @exons, @UTR5, @UTR3);
printAnnFile(\@EXON, "exons.ann.txt");
printAnnFile(\@exons, "coding.ann.txt");
printAnnFile(\@UTR5, "utr5.ann.txt");
printAnnFile(\@UTR3, "utr3.ann.txt");
printAnnFile(\@introns, "introns.ann.txt");

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
