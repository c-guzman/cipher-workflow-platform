# C I P H E R
Version 1.0.0 | Updated April 2017

Author: Carlos Guzman
E-mail: cag104@ucsd.edu

CIPHER is a data processing pipeline suite for next generation sequencing data including ChIP-seq, RNA-seq, DNase-seq, MNase-seq, and GRO-seq. By taking advantage of the Nextflow language, and Docker containers, CIPHER is an extremely easy to use, and reproducible pre-processing workflow toolkit.

CIPHER also includes a simple random forest-based ML model for the de-novo identification of enhancer elements based on DNase, H3K4me1, H3K27Ac, and H3K4me3 bedGraph files.

## HELP

CIPHER has a built in help command. For more information regarding possible parameters and their meanings, open up the command line terminal and type:

```
nextflow run main.nf --help
```

This will bring up the following help information:

```
~ C I P H E R ~ Version 1.0.0
**************************************************************

REQUIRED PARAMETERS:
--mode				Choose from available: chip, rna, gro, mnase, dnase, analysis.
--config			Configuration file with sample information. Check README for more information.
--fasta				Reference genome in FASTA format.
--gtf				Reference genome in GTF format.
--lib				Library information. "s" for single-stranded data, "p" for pair-ended data.
--readLen			The length of your reads.

RNA-seq ONLY:
--strandInfo			Strandedness information. Choose from "unstranded", "frFirstStrand", or "frSecondStrand".
--expInfo			Experiment config file for RNA-seq data DGE analysis. Check README for more information.

ANALYSIS MODE ONLY:
--analysis			Choose from available: "predictEnhancers".

OPTIONAL PARAMETERS:
--threads			Number of threads. (Default: 1)
--minid				Minimum alignment identity to look for during mapping. Higher is faster and less sensitive. (Default: 0.76)
--qvalue			Minimum FDR cutoff for peak detection in MACS2 and EPIC. (Default: 0.05)
--epic_w			Size of the windows used to scan the genome for peak detection in EPIC. (Default: 200)
--epic_g			A multiple of epic_w used to determine the gap size in EPIC. (Default: 3)
--maxindel			Maximum indel length searched during mapping. 200k recommended for vertebrate genomes. (Default: 200k)
--intronlen			Maximum intron length during mapping. 20 recommended for vertebrate genomes. (Default: 20)
--outdir			Name of output directory. (Default: results)

--subsample			Set this flag to subsample reads for testing.

**************************************************************
```

**IMPORTANTLY** if your analysis is running out of memory / RAM you will need to make use of the **-qs** flag. This specifies the number of processes that can be run at one time in parallel. For smaller RAM computers I suggest to decrease the number to something more appopriate. On a 128gb of RAM local desktop, I typically use a -qs flag of 10.

## INSTALLATION

The CIPHER program requires no manual installation. Just download and unzip the tar file from GitHub or the link at the top of the page. The download may be large.

By default CIPHER will run off a Docker container.

The only required installs on your local computer / cluster are **Docker** and **Nextflow**.

## INPUT

The input and metadata should be specified in a config file passed in via the `--config` flag. The pipelines require slightly differently formatted config files.

**Single-End ChIP-seq, GRO-seq, DNase-seq, and MNase-seq** config file format:

```
sample1 sample1_run1  /path/to/sample1_run1.fastq.gz   control1    Pol-II
sample1 sample1_run2  /path/to/sample1_run2.fastq.gz   control1    Pol-II
sample3 sample2_run1  /path/to/sample2_run1.fastq.gz   control2    H3K4me3
control1  control1_run1 /path/to/control1_run1.fastq.gz  control1    input
control2  control2_run1 /path/to/control2_run1.fastq.gz  control2    input
```

1. **SampleID**: The prefix of your file. Typically the name of your mark. This ID will be used to merge technical replicates.
2. **ID**: The prefix of your fastq files. Must be the same as the prefix of your .fastq.gz file. This ID is used to differentiate between different technical replicates.
3. **Path**: The path to the fastq file to be processed.
4. **ControlID**: The prefix of the Input file to be used for peak calling. Must be the same as the control sampleID. `-` if there is no control.
5. **Mark**: The mark/histone being processed or `input` if the line refers to a control. If this is not a control or input file, then this should be the same as your SampleID.

**Pair-End ChIP-seq, GRO-seq, DNase-seq, and MNase-seq** config file format:

```
sample1 sample1_run1  /path/to/sample1_run1_R1.fastq.gz /path/to/sample1_run1_R2.fastq.gz   control1    Pol-II
sample1 sample1_run2  /path/to/sample1_run2_R1.fastq.gz /path/to/sample1_run2_R2.fastq.gz   control1    Pol-II
sample3 sample2_run1  /path/to/sample2_run1_R1.fastq.gz /path/to/sample2_run1_R2.fastq.gz   control2    H3K4me3
control1  control1_run1 /path/to/control1_run1_R1.fastq.gz  /path/to/control1_run1_R2.fastq.gz  control1    input
control2  control2_run1 /path/to/control2_run1_R1.fastq.gz  /path/to/control2_run1_R2.fastq.gz  control2    input
```

1. **SampleID**: The prefix of your file. Typically the name of your mark. This ID will be used to merge technical replicates.
2. **ID**: The prefix of your fastq files. Must be the same as the prefix of your .fastq.gz file. This ID is used to differentiate between different technical replicates.
3. **Path1**: The path to the first fastq file to be processed. Must include _R1 at the end of file prefixes in order to differentiate between pair ended files.
4. **Path2**: The path to the second fastq file to be processed. Must include _R2 at the end of the file prefixes in order to differentiate between pair ended files.
5. **ControlID**: The prefix of the Input file to be used for peak calling. Must be the same as the control sampleID. `-` if there is no control.
6. **Mark**: The mark/histone being processed or `input` if the line refers to a control. If this is not a control or input file, then this should be the same as your SampleID.

**Single-End RNA-seq** config file format:

```
sample1 sample1_run1  /path/to/sample1_run1.fastq.gz
sample1 sample1_run2  /path/to/sample1_run2.fastq.gz
```

1. **SampleID**: The prefix of your file. Typically the name of your mark. This ID will be used to merge technical replicates for future splicing analysis. DOES NOT MERGE REPLICATES FOR DGE.
2. **ID**: The prefix of your fastq files. Must be the same as the prefix of your .fastq.gz file.
3. **Path1**: The path to the first fastq file to be processed.

**Pair-End RNA-seq** config file format:

```
sample1 sample1_run1  /path/to/sample1_run1_R1.fastq.gz /path/to/sample1_run1_R2.fastq.gz
sample1 sample1_run2  /path/to/sample1_run2_R1.fastq.gz /path/to/sample1_run2_R2.fastq.gz
```

1. **SampleID**: The prefix of your file. Typically the name of your mark. This ID will be used to merge technical replicates for future splicing analysis. DOES NOT MERGE REPLICATES FOR DGE.
2. **ID**: The prefix of your fastq files. Must be the same as the prefix of your .fastq.gz file. This ID is used to differentiate between different technical replicates and conditions.
3. **Path1**: The path to the first fastq file to be processed. Must include _R1 at the end of file prefixes in order to differentiate between pair ended files.
4. **Path2**: The path to the second fastq file to be processed. Must include _R2 at the end of the file prefixes in order to differentiate between pair ended files.

**RNA-Seq Experiment Information File**

RNA-seq experiments require a text file with metadata regarding your samples and conditions. Currently, CIPHER only really supports 2 condition DGE (such as WT/KO). The experiment file should look like this:

```
sample	condition
Ctrl1   WT
Ctrl2	WT
Ctrl3	WT
KO1 	KO
KO2		KO
KO3		KO
```

The headers sample and condition must be included.

1. **Sample**: Refers to the SampleID in your config file. Used to differentiate between conditions and replicates.
2. **Condition**: Refers to the condition of your sample. CIPHER currently only supports two condition DGE analysis.

## Getting Started Guide

In this section we'll go over the basics of how to install CIPHER's dependencies and then run a simple ChIP-seq analysis. 

1) Download or clone CIPHER from either GitHub or the website.

```
git clone https://github.com/c-guzman/cipher-workflow-platform.git
```

2) Create a **config** file by following the format specified in the **Input** section.

3) Download and save a **fasta** file, and a **gtf** file of your reference genome. If you are unsure about where to look for these files, then take a look at the **Annotation Files** section. CIPHER does not currently have a de-novo assembly workflow. Reference files can be downloaded for a select few genomes using the `--download` flag (e.g. `--download hg19` for hg19 fasta and gtf files).

4) Run the CIPHER ChIP-seq workflow for single-stranded human sample data.

```
nextflow run /path/to/main.nf --mode chip --config config.txt --fasta genome.fa --gtf genome.gtf --lib s --readLeng 75 --threads 10
```

This command will run the ChIP-seq workflow for the samples specified in your config.txt using the reference fasta and gtf files you provided for single-ended data in the hg19 human genome. Other genomes are available (hg38, mm9, mm10), and custom genomes can be set using the **--macs_size** and **--epic_size** flags.

## Annotation Files

Annotation files are typically FASTA and GTF files for various genomes / species.

**Gencode**: http://www.gencodegenes.org/ (human and mouse genomes)

**UCSC**: https://genome.ucsc.edu/cgi-bin/hgTables (very comphrensive list of genomes)

**Ensembl**: http://useast.ensembl.org/info/data/ftp/index.html (very comphrensive list of genomes .. ChIP-seq requires you to format chromosome names into UCSC based chromosomes (chr1, chr2, chr3, chr ... etc))

## Analysis Mode

The analysis mode in CIPHER is meant to streamline more sophisticated and complex downstream analysis. 

**De-Novo Enhancer Prediction**

CIPHER is able to predict potential enhancer elements via a machine learning model using bedGraph files for DNase, H3K4me1, H3K4me3, and H3K27Ac markers. The bedGraph files must be RPKM normalized, but running CIPHER automatically outputs RPKM normalized bedGraphs.

To identify enhancers you will need a **config** file with marker information, and a reference genome **FASTA** file.

The tab separated **config** file looks like this:

```
DNase	  /path/to/DNase.bedGraph
H3K27Ac	  /path/to/H3K27Ac.bdg
H3K4me1	  /path/to/H3K4me1.bdg
H3K4me3	  /path/to/H3K4me3.bdg
```

Where column one indicates the marker ID **(THESE MUST ALWAYS BE THE SAME, THE PIPELINE WILL FAIL IF YOU DO NOT INCLUDE THE APPROPRIATE IDs)** and column 2 indicates the path to the individual bedGraph file.

Once you have the necessary files, you can identify enhancers by running:

```
nextflow run /path/to/main.nf --mode analysis --config config.txt --fasta genome.fa --predictEnhancers

```
## Cluster support

CIPHER is possible to execute it on your computer or any cluster resource
manager without modifying it.

Currently the following platforms are supported:

  + Oracle/Univa/Open Grid Engine (SGE)
  + Platform LSF
  + SLURM
  + PBS/Torque


By default the pipeline is parallelized by spanning multiple threads in the machine where the script is launched.

To submit the execution to a SGE cluster create a file named `nextflow.config`, in the directory
where the pipeline is going to be launched, with the following content:

    process {
      executor='sge'
      queue='<your queue name>'
    }

In doing that, tasks will be executed through the `qsub` SGE command, and so your pipeline will behave like any
other SGE job script, with the benefit that *Nextflow* will automatically and transparently manage the tasks
synchronisation, file(s) staging/un-staging, etc.
