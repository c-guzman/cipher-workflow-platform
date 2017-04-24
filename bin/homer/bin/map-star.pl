#!/usr/bin/env perl
use warnings;
#


my $exec = "STAR";

my $options = "";
#my $options = " --alignSplicedMateMapLminOverLmate 0.05";

if (@ARGV < 3) {
	print STDERR "\n\tUsage: map-star.pl [optoins] -p <cpu#> -x <genome index dir> <file1> <file2> ...\n";
	print STDERR "\tFor paired-end, enter files as  <file1.r1>,<file1.r2> <file2.r1>,<file2.r2>...\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-p <#> (number of CPU/cores to use)\n";
	print STDERR "\t\t-pe (Look for paired end file if only one is given)\n";
	print STDERR "\t\t-x <index directory> (Directory with STAR index for mapping)\n";
	print STDERR "\t\t-mis <#> (max mismatches allowed, default: STAR default, 10)\n";
	print STDERR "\t\t-cufflinksopts (use cufflinks friendly options when mapping)\n";
	print STDERR "\t\t-leave (Leave genome index in shared memory, by default it is removed)\n";
	print STDERR "\t\t-path <path-to-star> (Default: expects 'STAR' in your PATH)\n";
	print STDERR "\t\t-NoSharedMemory (Do not used shared memory for genome [slower])\n";
	print STDERR "\t\t-pass \"...\" (pass cmd options to STAR, include the quotes!)\n";
	print STDERR "\n";
	exit;
}


my @files = ();
my $numCPUs = 1;
my $gdir = "";
my $removeSharedMemory = 1;
my $peFlag = 0;
my $unmappedFlag = ' --outReadsUnmapped Fastx ';
#my $unmappedFlag = '
my $genomeLoad = "LoadAndKeep";
my $pass = "";

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-p' || $ARGV[$i] eq '-cpu') {
		$numCPUs = $ARGV[++$i];
		print STDERR "\tWill use $numCPUs CPUs with STAR...\n";
	} elsif ($ARGV[$i] eq '-x' || $ARGV[$i] eq '-index') {
		$gdir = $ARGV[++$i];
		print STDERR "\tGenome Index: $gdir\n";
	} elsif ($ARGV[$i] eq '-pe') {
		$peFlag = 1;
		print STDERR "\tWill look for a paired-end file\n";
	} elsif ($ARGV[$i] eq '-NoSharedMemory') {
		$genomeLoad = "NoSharedMemory";
	} elsif ($ARGV[$i] eq '-leave') {
		$removeSharedMemory = 0;
	} elsif ($ARGV[$i] eq '-path') {
		$exec = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pass') {
		$pass .= " " . $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-mis') {
		$options .= " --outFilterMismatchNmax $ARGV[++$i]";
	} elsif ($ARGV[$i] eq '-cufflinksopts') {
		print STDERR "\tTo help cufflinks work better, the following options will be added:\n";
		print STDERR "\t\t--outSAMstrandField intronMotif\n";
		#$options .= " --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical";
		$options .= " --outSAMstrandField intronMotif ";
	} else {
		push(@files, $ARGV[$i]);
	}
}
	

if (`which "$exec"`) {
} else {
	print STDERR "\n!!! Could not find STAR installed !!!\n";
	exit;
}

my $indexName = $gdir;
$indexName =~ s/\/+$//;
$indexName =~ s/^.*\///;
$indexName =~ s/Genome_//;
print STDERR "\tIndex (short) Name: $indexName\n";

print STDERR "\tFiles to be analyzed:\n";
foreach(@files) {
	print STDERR "\t\t$_\n";
}
print "\n";


if ($removeSharedMemory && $genomeLoad ne "NoSharedMemory") {
	`$exec --genomeLoad Remove --genomeDir $gdir`;
	`rm -f Log.progress.out Aligned.out.sam Log.out`;
	`rm -fr _tmp/`;
}


my %files = ();
foreach(@files) {
	my $file = $_;
	my %delete = ();

	my $file1 = $file;
	my $file2 = "";
	if ($file =~ /\,/) {
		# paired end
		my @a = split /\,/, $file;
		$file1 = $a[0];
		$file2 = $a[1];
	}
	my $f1 = $file1;
	my $f2 = $file2;

	if ($peFlag) {
		if ($f2 ne '') {
			#no problem
		} else {


			my $testFile = $f1;

			if ($testFile =~ s/_R1_/_R2_/ || $testFile =~ s/\.R1\./\.R2\./ || $testFile =~ s/_R1\./_R2\./) {
				if (-e $testFile) {
					print STDERR "\tUsing PE file: $f1 & $testFile\n";
					$file2 = $testFile;
					$f2 = $testFile;
				}
			} elsif ($testFile =~ /SRR.*_1\.fastq/) {
				$testFile =~ s/_1\.fastq/_2\.fastq/;
				if (-e $testFile) {
					print STDERR "\tUsing PE file: $f1 & $testFile\n";
					$file2 = $testFile;
					$f2 = $testFile;
				}
			} else {
				print STDERR "!!! Warning !!! Could not find paired end file for $f1\n";
			}
		}
	}

	if ($f1 =~ s/\.gz$//) {
		print STDERR "\tUnzipping $file1 -> $f1\n";
		`gunzip -c "$file1" > "$f1"`;
		$delete{$f1} = 1;
	}
	if ($f2 =~ s/\.gz$//) {
		print STDERR "\tUnzipping $file2 -> $f2\n";
		`gunzip -c "$file2" > "$f2"`;
		$delete{$f2} = 1;
	}
	my $outputPrefix = $f1 . "." . $indexName;
	print STDERR "`$exec --genomeDir $gdir --runThreadN $numCPUs --readFilesIn $f1 $f2 --outFileNamePrefix $outputPrefix. $options $pass`;\n";
	`$exec --genomeLoad $genomeLoad $unmappedFlag --genomeDir "$gdir" --runThreadN $numCPUs --readFilesIn "$f1" "$f2" --outFileNamePrefix "$outputPrefix." $options $pass`;

	# clean up unwanted output files
	`rm -f "$outputPrefix.Log.out" "$outputPrefix.Log.progress.out" "$outputPrefix.SJ.out.tab"`;
	my $afile = $outputPrefix . ".Aligned.out.sam";
	$files{$afile}=$outputPrefix;

	foreach(keys %delete) {
		my $f = $_;
		`rm -f "$f"`;
	}

}
if ($removeSharedMemory && $genomeLoad ne "NoSharedMemory") {
	`$exec --genomeLoad Remove --genomeDir $gdir`;
	`rm -f Log.progress.out Aligned.out.sam Log.out`;
	`rm -rf _tmp/`;
}


my $numCPUs2 = $numCPUs;
$numCPUs = 1;

my $cpus = 0;
my @pids = ();
foreach(keys %files) {
	my $file = $_;
	my $ogFile = $files{$file};
	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		`samtools view -b -S -@ $numCPUs2 $file > $file.unsorted.bam`;
		`samtools sort -@ $numCPUs2 $file.unsorted.bam $ogFile`;
		`samtools index $ogFile.bam`;
		`rm -f $file.unsorted.bam`;
		exit(0);
	}
	push(@pids, $pid);
	if ($cpus >= $numCPUs) {
		my $id = wait();
		$cpus--;
	}
}
my $id = 0;
while ($id >= 0) {
	$id = wait();
}
