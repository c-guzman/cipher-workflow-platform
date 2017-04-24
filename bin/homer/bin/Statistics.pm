package Statistics;
use POSIX;

my $pairedData = 0;

sub logPoisson {
	my ($lambda, $k) = @_;
	my $v = $k*log($lambda) - $lambda - factln($k);
	return $v;
}

sub cumulativePoisson {
	my ($lambda, $k) = @_;
	my $a = $lambda;
	my $b = $k;
	if ($b > $a) {
		my $c = $b;
		$b = $a;
		$a = $c;
	}
	if ($a == 0) {
		return 1;
	}

	my $P = 0;
	for ($i=0;$i<=$b;$i++) {
		my $p = $i*log($a) - $a - factln($i);
		$P += exp($p);
		next;
		if ($i==0) {
			$P = $p;
		} else {
        	$P += log(1+exp($p-$P));
		}
	}
	return $P;
}

# Kolmogorov-Smirnov - more or less from nr.com
sub ks {
	my ($data1, $data2) = @_;
	
	my $j1=0;
	my $j2=0;
	my $fn1=0;
	my $fn2=0;
	my @data1 = sort {$a <=> $b} @$data1;
	my @data2 = sort {$a <=> $b} @$data2;
	my $en1 = scalar(@data1);
	my $en2 = scalar(@data2);
	my $n1 = $en1;
	my $n2 = $en2;
	my $d = 0;
	while ($j1 < $n1 && $j2 < $n2) {
		$d1 = $data1[$j1];
		$d2 = $data2[$j2];
		if ($d1 <= $d2) {
			$fn1 = $j1/$en1;
			$j1++;
		}
		if ($d2 <= $d1) {
			$fn2=$j2/$en2;
			$j2++;
		}
		$dt = abs($fn2-$fn1);
		if ($dt > $d) {
			$d = $dt;
		}
	}
	my $en = sqrt(($en1*$en2)/($en1+$en2));
	return probks(($en+0.12+(0.11/$en))*$d);
}

sub probks {
	my ($alam) = @_;
	my $a2 = -2*$alam*$alam;
	my $term = 0;
	my $fac = 2;
	my $sum = 0;
	my $EPS1 = 0.001;
	my $EPS2 = 1e-8;
	my $termbf = 0;
	for (my $j=1;$j<=100;$j++) {
		$term = $fac*exp($a2*$j*$j);
		$sum+=$term;
		if (abs($term) <= $EPS1*$termbf || abs($term) <= $EPS2*$sum) {
			return $sum;
		}
		$fac = -1*$fac;
		$termbf=abs($term);
	}
	return 1;
}

sub friedman {
	my ($data) = @_;

	my $numA = scalar(@$data);
	my $numB = scalar(@{$data->[0]});
	my $numR = scalar(@{$data->[0]->[0]});
	if ($numA < 2 || $numB < 2 || $numR < 1) {
		print STDERR "Data format error: $numA $numB $numR\n";
	}

	my @R = ();
	for (my $i=0;$i<$numA;$i++) {
		push(@R, 0);
	}
	for (my $i=0;$i<$numB;$i++) {
		my @sub = ();
		for (my $j=0;$j<$numA;$j++) {
			for (my $k=0;$k<$numR;$k++){ 
				push(@sub, $data->[$j]->[$i]->[$k]);
			}
		}
		my $order = getOrder(\@sub);
		for (my $j=0;$j<$numA;$j++) {
			my $R = 0;
			for (my $k=0;$k<$numR;$k++){ 
				$R += $order->[$j*$numR+$k];
			}
			$R[$j] += $R;
		}
	}
	my $chi2 = 0;
	my $expect = $numB*$numR*($numA*$numR+1)/2;
	#print STDERR "$expect\n";
	
	for (my $i=0;$i<$numA;$i++) {
		#print STDERR "$R[$i]\n";
		my $v = $R[$i]-$expect;
		$chi2 += $v*$v;
	}
	#print STDERR "$chi2\t$numB\t$expect\n";
	#print STDERR "$chi2\n";
	my $aa = ($numB*$numA*$numR*$numR*($numR*$numA+1));
	my $coeff = 12/$aa;
	#print STDERR "$aa\t$coeff";
	$chi2 *= $coeff;
	#print STDERR "$chi2\n";
	my $df = $numA-1;

	my $pvalue = chi2pvalue($chi2, $df);
	return ($chi2, $pvalue);
}


# $pvalues is an array of p-values
# $results is an array of array of p-values (2nd index is number of randomizations)
sub empiricalFDR {
	my ($pvalues, $results) = @_;

	my $numRandomizations = scalar(@{$results->[0]});
	my @p = ();
	foreach(@$results) {
		foreach(@$_) {
			push(@p, $_);
		}
	}
	return empiricalFDR2($pvalues, \@p, $numRandomizations);
}

# $pvalues is an array of p-values
# $results is an array of all p-values from randomizations
sub empiricalFDR2 {
	my ($pvalues, $results, $numRandomizations) = @_;

	my @order = ();
	for (my $i=0;$i<@$pvalues;$i++) {
		my $d = {p=>$pvalues->[$i],i=>$i,f=>1,c=>0};
		push(@order, $d);
	}
	@order = sort {$a->{'p'} <=> $b->{'p'}} @order;

	my @p = sort {$a <=> $b} @$results;

	my $pIndex = 0;
	my $lastFDR = -1e10;
	for (my $i=0;$i<@order;$i++) {
		while ($pIndex < @p && $p[$pIndex] < $order[$i]->{'p'}) {
			$pIndex++;
		}
		my $FDR = (($pIndex+1)/$numRandomizations)/($i+1);
		$FDR = $lastFDR if ($FDR < $lastFDR);
		if ($FDR > 1) {
			$FDR = 1;
		}
		$lastFDR = $FDR;
		$order[$i]->{'f'} = $FDR;
	}
	my @fdr = ();
	@order = sort {$a->{'i'} <=> $b->{'i'}} @order;
	foreach(@order) {
		push(@fdr, $_->{'f'});
	}
	return \@fdr;
}
	

#benjamini-hochberg FDR
sub benjaminiFDR {
	my ($pvalues, $numTests) = @_;

	my @order = ();
	for (my $i=0;$i<@$pvalues;$i++) {
		my $d = {p=>$pvalues->[$i],i=>$i,f=>1};
		push(@order, $d);
	}
	@order = sort {$a->{'p'} <=> $b->{'p'}} @order;

	if (!defined($numTests)) {
		$numTests = scalar(@$pvalues);
	}

	my $lastFDR = -1e10;
	my $lastIndex = -1;
	for (my $i=0;$i<@order;$i++) {
		if ($i<@order-1) {
			# if it's a tie, need to count all of them in the calculation
			next if ($order[$i+1]->{'p'} == $order[$i]->{'p'});
		}
		my $fdr = $order[$i]->{'p'}*$numTests/($i+1);
		if ($fdr < $lastFDR) {
			$fdr = $lastFDR;
		}
		if ($fdr > 1) {
			$fdr = 1;
		}
		for (my $j=$i;$j>$lastIndex;$j--) {
			$order[$j]->{'f'} = $fdr;
		}
		$lastFDR = $fdr;
		$lastIndex = $i;
	}

	my @fdr = ();
	@order = sort {$a->{'i'} <=> $b->{'i'}} @order;
	foreach(@order) {
		push(@fdr, $_->{'f'});
	}
	return \@fdr;
}


sub anova2 {
	my ($data) = @_;

	my $numA = scalar(@$data);
	my $numB = scalar(@{$data->[0]});
	my $numR = scalar(@{$data->[0]->[0]});
	if ($numA < 2 || $numB < 2 || $numR < 1) {
		print STDERR "Data format error: $numA $numB $numR\n";
	}

	#calculate Means
	my @MB = ();
	my @MAB = ();
	my @MA = ();
	my $MT = 0;
	for (my $i=0;$i<$numA;$i++) {
		my @a = ();
		for (my $j=0;$j<$numB;$j++) {
			push(@a, 0);
		}
		push(@MAB, \@a);
		push(@MA, 0);
	}
	for (my $i=0;$i<$numB;$i++) {
		push(@MB, 0);
	}
	for (my $i=0;$i<$numA;$i++) {
		for (my $j=0;$j<$numB;$j++) {
			for (my $k=0;$k<$numR;$k++) {
			#if (!defined($numR)) {
			#	print STDERR "$bad\n";
			#}
				my $v = $data->[$i]->[$j]->[$k];
				$MA[$i] += $v;
				$MB[$j] += $v;
				$MAB[$i][$j] += $v;
				$MT += $v;
			}
		}

	}
	for (my $i=0;$i<$numA;$i++) {
		for (my $j=0;$j<$numB;$j++) {
			$MAB[$i][$j] /= $numR;
		}
		$MA[$i] /= $numB*$numR;
	}
	for (my $i=0;$i<$numB;$i++) {
		$MB[$i] /= $numA*$numR;
	}
	$MT /= $numA*$numB*$numR;


	#calculate sum of squares
	my $SST = 0;
	my $SSA = 0;
	my $SSB = 0;
	my $SSAB = 0;
	my $SSE = 0;
	for (my $i=0;$i<$numA;$i++) {
		for (my $j=0;$j<$numB;$j++) {
			for (my $k=0;$k<$numR;$k++) {
				my $v = $data->[$i]->[$j]->[$k];
				$SST += ($v-$MT)*($v-$MT);
				$SSE += ($v - $MAB[$i][$j])*($v-$MAB[$i][$j]);
			}
			my $vv= ($MAB[$i][$j] - $MA[$i] -$MB[$j] +$MT);
			$SSAB += $vv*$vv;
		}
		$SSA += ($MA[$i]-$MT)*($MA[$i]-$MT);
	}
	for (my $i=0;$i<$numB;$i++) {
		$SSB += ($MB[$i]-$MT)*($MB[$i]-$MT);
	}
	$SSA *= $numB*$numR;
	$SSB *= $numA*$numR;
	$SSAB *= $numR;

	my $DFT = $numA*$numB*$numR -1;
	my $DFA = $numA - 1;
	my $DFB = $numB - 1;
	my $DFAB = ($numA-1)*($numB-1);
	my $DFE = $numA*$numB*($numR - 1);

	my $MSA = $SSA/$DFA;
	my $MSB = $SSB/$DFB;
	my $MSAB = $SSAB/$DFAB;
	my $MSE = $SSE/$DFE;

	my $FA = $MSA/$MSE;
	my $FB = $MSB/$MSE;
	my $FAB = $MSAB/$MSE;

	my $pA = FPvalue($FA,$DFA, $DFE);
	my $pB = FPvalue($FB, $DFB, $DFE);
	my $pAB = FPvalue($FAB, $DFAB, $DFE);
	#print STDERR "$MSAB\t$MSE\t$FAB\t$pAB\n";
	
	return ($FA, $pA);

}

sub FPvalue {
	my ($F, $v1, $v2) = @_;
	my $x = $v2/($v2+$v1*$F);
	my $a = $v2/2;
	my $b = $v1/2;
	my $p = betai($a, $b,$x);
	return $p;
}

#input is $data, which is an array of arrays holding the data.
sub kruskalWallace {
	my ($data) = @_;
	my @N = ();
	my @d = ();
	my $N = 0;
	foreach(@$data) {
		my $n =0;
		foreach(@$_) {
			push(@d, $_);
			$n++;
		}
		push(@N, $n);
		$N+= $n;
	}
	my $rank = getOrder(\@d);
	my $i=0;
	my $KW = 0;
	my @avgRank=();
	for (my $j=0;$j<@N;$j++) {
		my $R = 0;
		for (my $k=0;$k<$N[$j];$k++) {
			$R+=$rank->[$i];
			$i++;
		}
		$avgRank[$j] = ($R / $N[$j])/$N;
			
		my $v = ($R-$N[$j]*($N+1)/2);
		$KW += 1/$N[$j]*$v*$v;
	}

	$KW *= 12/($N*($N+1));
	my $df = scalar(@N)-1;

	my $pvalue = chi2pvalue($KW, $df);
	return ($KW, $pvalue, \@avgRank);
}

sub kruskalWallace2 {
	my ($R,$n,$N) = @_;
	my $R2 = $N*($N+1)/2 - $R;

	my $v = $R-$n*($N+1)/2;
	my $KW = 1/$n*$v*$v;
	$v= $R2-($N-$n)*($N+1)/2;
	$KW += 1/($N-$n)*$v*$v;
	$KW *= 12/($N*($N+1));
	my $df = 1;
	my $pvalue = chi2pvalue($KW,$df);
	return ($KW,$pvalue);
}

sub randn {
	my ($v1,$v2,$r) = (0,0,0);
	do {
		$v1 = rand(2)-1;
		$v2 = rand(2)-1;
		$r = $v1*$v1+$v2*$v2;
	} while ($r >=1 || $r==0);
	return  $r*sqrt(-2*(log($r))/$r);
}

sub getOrder {
	my ($data) = @_;
	my @a = sort {$a <=> $b} @$data;
	my %value = ();
	for (my $i=0;$i<@a;$i++) {
		my $n = 0;
		my $j= $i+1;
		while ($j < @a && $a[$j] == $a[$i]) {
			$j++;
		}
		$j--;
		my $r = ($i+$j)/2+1;
		$value{$a[$i]} = $r;
		$i=$j;
	}
	my @rank = ();
	foreach(@$data) {
		push(@rank, $value{$_});
	}
	return \@rank;
}

sub foldChange {
	my ($data1, $data2, $LFlag) = @_;

	my $m1 = @$data1;
	my $m2 = @$data2;
	#print STDERR "== $m1 $m2\n";

	my @fc = ();
	my ($avg1, $var1) = avevar($data1);
	my ($avg2, $var2) = avevar($data2);

	my ($min1,$max1,$min2,$max2) = (1e30, -1e30,1e30,-1e30);


	my $cc = 0;
	foreach(@$data1) {
		$max1 = $_ if ($_>$max1);
		$min1 = $_ if ($_<$min1);
		my $d = $_;
		next if ($d == 0);
		if ($pairedData) {
			if ($LFlag) {
				push(@fc, exp(($data2->[$cc]-$d)*log(2)));
			} else {
				push(@fc, $data2->[$cc]/$d);
			}
		} else {
			foreach(@$data2) {
				if ($LFlag) {
					push(@fc, exp(($_-$d)*log(2)));
				} else {
					push(@fc, $_/$d);
				}
			}
		} 
		$cc++;
	}

	foreach(@$data2) {
		$max2 = $_ if ($_>$max2);
		$min2 = $_ if ($_<$min2);
	}
	my @sfc = sort {$a <=> $b} @fc;
	

	my $FC =0;
	if ($LFlag) {
		$FC = exp(($avg2-$avg1)*log(2));
	} else {
		if ($avg1 < 0.00001) {
			print STDERR "AVG1 is 0\n";
			$FC = 1;
		} else {
			$FC = $avg2 / $avg1;
		}
	}

	my $extremeFC = 1;
	if ($FC > 1) {
		$max1 += .00001;
		if ($LFlag) {
			$extremeFC = exp(($min2-$max1)*log(2));
		} else {
			$extremeFC = $min2/$max1;
		}
		if ($pairedData) {
			$extremeFC = $sfc[0];
		}
	} else {
		$min1 += .00001;
		if ($LFlag) {
			$extremeFC = exp(($max2-$min1)*log(2));
		} else {
			$extremeFC = $max2/$min1;
		}
		if ($pairedData) {
			$extremeFC = $sfc[@sfc-1];
		}
	}

	my ($fcavg, $fcvar) = avevar(\@fc);
	return ($FC, $fcvar,$extremeFC);
}

sub permute {
	my ($total, $sample) = @_;
	my %rand = ();
	my @b = ();
	for (my $i=0;$i<$sample;$i++) {
		my $r= floor(rand()*$total);
		while(exists($rand{$r})) {
			$r = floor(rand()*$total);
		}
		$rand{$r} = 1;
		push(@b, $r);
	}
	return \@b;
}

# N = number of items per group
sub ttestVar {
	my ($avg1, $var1, $avg2, $var2, $N) = @_;

	my $df = $N*2-2;
	$svar = ($var1+$var2)/2;
	$t = ($avg1-$avg2) / sqrt($svar*2/$N);
	$p = betai(0.5*$df,0.5,$df/($df+$t*$t));
	return ($t,$p);
}


sub ttest {
	my ($data1,$data2) = @_;

	my ($ave1,$var1) = avevar($data1);
	my ($ave2,$var2) = avevar($data2);

	my $n1 = @$data1;
	my $n2 = @$data2;
	my $df = $n1+$n2-2;
	my $svar = (($n1-1)*$var1+($n2-1)*$var2)/$df;
	my $t =0;
	my $p = 1.0;
	if ($svar == 0) {
		$p = 0;
		$t = 1e10;
	} else {
		$t = ($ave1-$ave2)/sqrt($svar*(1.0/$n1+1.0/$n2));
		$p = betai(0.5*$df,0.5,$df/($df+$t*$t));
	}
	return ($t,$p);
}

sub avevar {
	my ($data) = @_;
	my $n = @$data;
	if ($n == 0) {
		print STDERR "N equal to zero\n";
	}
	my $ave = 0.0;
	for (my $j=0;$j<@$data;$j++) {
		$ave+= $data->[$j];
	}
	$ave /= @$data;
	if (@$data < 2) {
		return ($ave,0);
	}
	if (@$data < 3) {
		return ($ave, abs(($data->[0]-$data->[1])/2));
	}
	my $ep = 0.0;
	my $var = 0.0;
	for (my $j=0;$j<@$data;$j++) {
		my $s = $data->[$j]-$ave;
		$ep += $s;
		$var += $s*$s;
	}
	$var = ($var-$ep*$ep/$n)/($n-1);
	return ($ave,$var);
}

sub betai {
	my ($a,$b,$x) = @_;
	my ($bt) = (0.0);
	my $eeps = 1e-10;
	if ($x < 0.0-$eeps || $x > 1.0+$eeps) {
		print STDERR "X is not right value in betai\n";
		exit;
	}
	$x = 0 if ($x < 0);
	$x = 1.0 if ($x > 1);
	if ($x==0.0 || $x==1.0) {
		$bt = 0.0;
	}
	else {
		$bt = exp(gammln($a+$b)-gammln($a)-gammln($b)+$a*log($x)+$b*log(1.0-$x));
	}
	if ($x < ($a+1.0)/($a+$b+2.0)) {
		return $bt*betacf($a,$b,$x)/$a;
	}
	else {
		return 1.0-$bt*betacf($b,$a,1.0-$x)/$b;
	}
}

sub betacf {
	my ($a,$b,$x) = @_;
	my ($m,$m2,$aa,$c,$d,$del,$h,$qab,$qam,$qap) = (0,0,0,0,0,0,0,0,0,0);
	my ($MAXIT,$EPS,$FPMIN) = (100,3.0e-7,1.0e-30);

	$qab=$a+$b;
	$qap=$a+1.0;
	$qam=$a-1.0;
	$c=1.0;
	$d=1.0-$qab*$x/$qap;
	$d=$FPMIN if (abs($d) < $FPMIN);
	$d=1.0/$d;
	$h=$d;

	for (my $m=1;$m<=$MAXIT;$m++) {
		$m2=2*$m;
		$aa=$m*($b-$m)*$x/(($qam+$m2)*($a+$m2));
		$d=1.0+$aa*$d;
		$d=$FPMIN if (abs($d) < $FPMIN);
		$c=1.0+$aa/$c;
		$c=$FPMIN if (abs($c) < $FPMIN);
		$d=1.0/$d;
		$h *= $d*$c;
		$aa = -($a+$m)*($qab+$m)*$x/(($a+$m2)*($qap+$m2));
		$d=1.0+$aa*$d;
		$d=$FPMIN if (abs($d) < $FPMIN);
		$c=1.0+$aa/$c;
		$c=$FPMIN if (abs($c) < $FPMIN);
		$d=1.0/$d;
		$del=$d*$c;
		$h*=$del;
		last if (abs($del-1.0)<$EPS);
	}
    if ($m > $MAXIT) {
        print STDERR "Error calculating betacf\n";
        exit;
    }
    return $h;
}

sub gammln {
	my ($xx) = @_;
	
    my @cof=(76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,
            0.1208650973866179e-2,-0.5396239384953e-5);
	my ($x,$y,$tmp,$ser)=(0,0,0,0);
	$y=$xx;
	$x=$xx;
	$tmp=$x+5.5;
	$tmp -= ($x+0.5)*log($tmp);
	$ser=1.000000000190015;
	for (my $j=0;$j<=5;$j++) {
		$ser += $cof[$j]/(++$y);
	}
	return -1*$tmp+log(2.5066282746310005*$ser/$x);
}




sub chi2pvalue {
	my ($chi2, $df) = @_;
	return gammq($df/2, $chi2/2);
}

sub gammq {
	my ($a, $x) = @_;
	my $gamser = 0;
	my $gln = 0;
	if ($x < 0.0 || $a <= 0.0) {
		print STDERR "Bad data\n";
	}
	if ($x < ($a+1.0)) {
		my ($gamser, $gln) = gser($a, $x);
		return 1.0 - $gamser;
	} else {
		my ($gamser, $gln) = gcf($a,$x);
		return $gamser;
	}
}
		

sub gammp {
	my ($a, $x) = @_;

	if ($x < 0.0 || $a <= 0.0) {
		print STDERR "X:$x and A:$a are invalid\n";
	}
	if ($x < ($a+1.0)) {
		my ($gamser,$gln) = gser($a,$x);
		return $gamser;
	} else {
		my ($gamser,$gln) = gcf($a,$x);
		return 1.0-$gamser;
	}
}

sub gser {
	my ($a, $x) = @_;
	my $gln = gammln($a);
	my $gamser = 0;
	my $ITMAX = 100;
	my $EPS = 3.0e-7;
	if ($x <= 0.0) {
		$gamser = 0;
		return ($gamser, $gln);
	} else {
		my $ap = $a;
		my $del = 1.0/$a;
		my $sum = $del;
		for (my $i=1;$i<=$ITMAX;$i++) {
			$ap++;
			$del *= $x/$ap;
			$sum += $del;
			if (abs($del) < abs($sum)*$EPS) {
				$gamser = $sum*exp(-1*$x+$a*log($x)-$gln);
				return ($gamser, $gln);
			}
		}
		print STDERR "Reached ITMAX limit\n";
		return ($gamser, $gln);
	}
}

sub gcf {
	my ($a, $x) = @_;
	my $ITMAX=100;
	my $EPS=3.0e-7;
	my $FPMIN = 1.0e-30;
	my $gammcf = 0;
	my $gln = gammln($a);
	my $b = $x+1.0-$a;
	my $c = 1.0/$FPMIN;
	my $d = 1.0/$b;
	my $h = $d;
	my $i=1;
	for ($i=1;$i<=$ITMAX;$i++) {
		my $an = $i*($i-$a);
		$b += 2.0;
		$d=$an*$d+$b;
		if (abs($d) < $FPMIN) {
			$d=$FPMIN;
		}
		$c=$b+$an/$c;
		if (abs($c) < $FPMIN) {
			$c=$FPMIN;
		}
		$d=1.0/$d;
		$del=$d*$c;
		$h*=$del;
		last if (abs($del-1.0) < $EPS);
	}
	if ($i>$ITMAX) {
		print STDERR "Reached ITMAX limit\n";
	}
	$gammcf = exp(-1*$x+$a*log($x)-$gln)*$h;
	return ($gammcf, $gln);
}

sub randomizeArray {
	my ($array) = @_;
	my @index = ();
	for (my $i=0;$i<@$array;$i++) {
		push(@index, $i);
	}
	my @array = ();
	for (my $i=0;$i<@$array;$i++) {
		my $r = floor(rand()*scalar(@index));
		my $j = splice(@index, $r, 1);
		push(@array, $array->[$j]);
		
	}
	return \@array;
}


# n = sample size
# k = number of sucesses
# r = rate of success in null distribution
# N = number of examples in Null Distribtion
sub logbinomial {
	my ($n, $k, $r, $N) = @_;
	my $Llimit = 1/$N;
	my $Hlimit = ($N-1)/$N;
	if ($r < $Llimit) {
		$r = $Llimit;
	}
	if ($r > $Hlimit) {
		$r = $Hlimit;
	}
	
	if ($k==0) {
		return 0;
	}
	return logbetai($k,$n-$k+1,$r);
}
sub ilogbinomial {
	my ($n, $k, $r, $N) = @_;
	my $Llimit = 1/$N;
	my $Hlimit = ($N-1)/$N;
	if ($r < $Llimit) {
		$r = $Llimit;
	}
	if ($r > $Hlimit) {
		$r = $Hlimit;
	}
	$r = 1.0-$r;
	$k = $n-$k;
	
	if ($k==0) {
		return 0;
	}
	return logbetai($k,$n-$k+1,$r);
}

sub binomial {
	my ($n, $k, $r, $N) = @_;
	my $Llimit = 1/$N;
	my $Hlimit = ($N-1)/$N;
	if ($r < $Llimit) {
		$r = $Llimit;
	}
	if ($r > $Hlimit) {
		$r = $Hlimit;
	}
	
	if ($k==0) {
		return 0;
	}
	return betai($k,$n-$k+1,$r);
}

# N=total population
# n1=population size sample 1
# n2=population size sample 2
# n = size of overlap
sub loghypergeo {
    my ($N, $n1, $n2, $n) = @_;
    my $c3 = bicoln($N,$n2);
    my $P = 0;
    my $AA = ($n1>$n2)?$n2:$n1;
    for (my $i=$n;$i<=$AA;$i++) {
        my $c1 = bicoln($n1, $i);
        my $c2 = bicoln($N-$n1,$n2-$i);
		my $p = (($c1+$c2)-$c3);
		if ($i==$n) {
			$P = $p;
		} else {
        	$P += log(1+exp($p-$P));
		}
    }
    return $P;
}

sub convertLogPvalue {
	my ($logP) = @_;
	if ($logP > -500) {
		return exp($logP);
	} else {
		my $exponent = $logP/2.3026;
		return "1e$exponent";
	}
}


sub hypergeo {
    my ($N, $n1, $n2, $n) = @_;
    my $c3 = bicoln($N,$n2);
    my $P = 0;
    my $AA = ($n1>$n2)?$n2:$n1;
    for (my $i=$n;$i<=$AA;$i++) {
        my $c1 = bicoln($n1, $i);
        my $c2 = bicoln($N-$n1,$n2-$i);
        $P += exp(($c1+$c2)-$c3);
    }
    return $P;
}

# g is the fraction of spectral density in maximum
# N is number of bp/2 or total number of densities.
sub logGPvalue {
	my ($g, $N) = @_;
	my $P = 0;

	for (my $p=1;$p<=1/$g;$p++) {
		my $sign = ((-1)**$p); 
		my $s = factln($p);
	}


}

# - add up the other tail
sub iloghypergeo {
    my ($N, $n1, $n2, $n) = @_;
    my $c3 = bicoln($N,$n2);
    my $P = 0;
	my $lowerLimit = $n1+$n2-$N;
	$lowerLimit = 0 if ($lowerLimit < 0);
    for (my $i=$n;$i>=$lowerLimit;$i--) {
        my $c1 = bicoln($n1, $i);
        my $c2 = bicoln($N-$n1,$n2-$i);
		my $p = (($c1+$c2)-$c3);
		if ($i==$n) {
			$P = $p;
		} else {
        	$P += log(1+exp($p-$P));
		}
    }
    return $P;
}

sub logbetai {
	my ($a, $b, $x) = @_;
	if ($x <= 0.0 || $x >= 1.0) {
		print STDERR "variable x has invalid value in function logbetai (x=$x)\n";
	}
	my $bt = gammln($a+$b)-gammln($a)-gammln($b)+$a*log($x)+$b*log(1.0-$x);
	if ($x < ($a+1.0)/($a+$b+2.0)) {
		return $bt+log(betacf($a,$b,$x))-log($a);
	} else {
		return log(1.0-exp($bt)*betacf($b,$a,1.0-$x)/$b);
	}

}

sub factln {
	my @a = ();
	my $n = shift;
	if ($n < 0) {
		print "Negative factorial\n";
	}
	if ($n <= 1) {
		return 0;
	}
	return gammln($n+1);

	#if ($n <= 100) {
	#	return $a[$n] ? $a[$n] : ($a[$n] = gammln($n+1));
		#return $a[$n] ? $a[$n] : ($a[$n] = gammln($n+1));
	#} else {
	#	return gammln($n+1);
	#}
}

sub bico {
	my $n = shift;
	my $k = shift;
	return floor(0.5 + exp(factln($n) - factln($k) - factln($n-$k)));
}

sub bicoln {
    my ($nn, $k) = @_;
    #return log(floor(0.5+exp(factln($nn)-factln($k)-factln($nn-$k))));
    return factln($nn)-factln($k)-factln($nn-$k);
}

sub correlation {
	my ($x,$y) = @_;

	my $n = scalar(@$x);
	return 0 if ($n==0 || scalar(@$y) == 0);
	my $xysum = 0;
	my $xsum = 0;
	my $ysum = 0;
	my $x2sum = 0;
	my $y2sum = 0;
	for (my $i=0;$i<$n;$i++) {
		$xysum += $x->[$i]*$y->[$i];
		$xsum += $x->[$i];
		$ysum += $y->[$i];
		$x2sum += $x->[$i]*$x->[$i];
		$y2sum += $y->[$i]*$y->[$i];
	}
	my $r = "NA";
	my $numerator = $xysum - $xsum*$ysum/$n;
	my $denomerator = ($x2sum - $xsum*$xsum/$n)*($y2sum-$ysum*$ysum/$n);
	return ("NA","NA") if ($denomerator <= 0);
	$denomerator = sqrt($denomerator);
	if ($denomerator > 0) {
		$r = $numerator / $denomerator;
	}
	my $df = $n-2;
	my $den = ((1-$r)*(1+$r));
	if ($den < 1e-10) {
		$den = 1e-10;
	}
	my $t = $r*sqrt($df/($den));
	my $logp = 1;
	if (abs($t) > 1e-5) {
		$logp = logbetai(0.5*$df,0.5,$df/($df+$t*$t));
	}

	return ($r,$logp);

}

##############################
#
sub quantileNorm {
	my ($matrix) = @_;

	my @data = ();
	my @exp = ();
	my $numGenes = 0;
	my @names = ();

	my $nsamples = scalar(@{$matrix->[0]});
	for (my $i=0;$i<@$matrix;$i++) {
		if ($i ==0) {
			for (my $j=0;$j<@{$matrix->[$i]};$j++) {
				my @a = ();
				push(@data, \@a);
				my %x = ();
				push(@exp, \%x);
			}
		}
		push(@names, $i);
		for (my $j=0;$j<@{$matrix->[$i]};$j++) {
			$exp[$j]->{$i} = $matrix->[$i][$j];
			push(@{$data[$j]}, $matrix->[$i][$j]);
		}
		$numGenes++;
	}


	#sort each array
	foreach(@data) {
		my @rdata = sort {$b <=> $a} @$_;
		$_ = \@rdata;
	}

	my @model = ();
	for (my $i=0;$i<@{$data[0]};$i++) {
		my $avg = 0;
		my $n = 0;
		for (my $j=0;$j<@data;$j++) {
			$avg += $data[$j][$i];
			$n++;
		}
		push(@model, $avg/$n);
	}

	#map values back to the genes
	for (my $i=0;$i<@data;$i++) {
		my @genes= keys %{$exp[$i]};
		@genes = sort {$exp[$i]->{$b} <=> $exp[$i]->{$a}} @genes;
		for (my $j=0;$j<@genes;$j++) {
			$exp[$i]->{$genes[$j]} = $model[$j];
		}
	}

	my @nmatrix = ();

	@names = sort {$a <=> $b} @names;	
	foreach(@names) {
		my $name = $_;
		my @a = ();
		foreach(@exp) {
			my $v = $_->{$name};
			push(@a, $v);
		}
		push(@nmatrix, \@a);
	}
	return \@nmatrix;
}


1;
