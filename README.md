# C I P H E R
Version 1.0.0 | Updated August 2017

Author: Carlos Guzman
E-mail: cag104@ucsd.edu

CIPHER is a data processing workflow platform for next generation sequencing data including ChIP-seq, RNA-seq, DNase-seq, MNase-seq, ATAC-seq and GRO-seq. By taking advantage of the Nextflow language, and Singularity containers, CIPHER is an extremely easy to use, and reproducible pre-processing workflow toolkit.

## HELP

CIPHER has a built in help command. For more information regarding possible parameters and their meanings, open up the command line terminal and type:

```
nextflow run cipher.nf --help
```

## Installation

Download or `git clone` this repository and install dependencies.

The only required dependencies to run CIPHER is:

  + Nextflow (https://www.nextflow.io/)
  + Singularity (http://singularity.lbl.gov/index.html)

## CONFIG Files

Config files are tab separated text files with 5 columns for single-ended data and 6 columns for pair ended data.

Single-ended CONFIG:

```
sample1		sample1_rep1	/path/to/fastq.gz 	control1	sample1
sample2		sample2_rep1	/path/to/fastq.gz 	control1	input
```

Pair-ended CONFIG:

```
sample1		sample1_rep1	/path/to/fastq_R1.gz 	/path/to/fastq_R2.gz	control1	sample1
sample2		sample2_rep1	/path/to/fastq_R1.gz  	/path/to/fastq_R2.gz	control1	input
```

**DO NOT MIX AND MATCH SINGLE AND PAIR ENDED DATA INTO THE SAME CONFIG FILE. CIPHER DOES NOT HANDLE THIS USE-CASE YET.**

Where columns refer to:

* 1. MergeID - Prefix used for naming files that are merged together.
* 2. SampleID - Prefix used for naming files that are not merged together. Typically includes replicate information.
* 3. Path 1 - The file path to first FASTQ file. Typically the R1 file in pair-ended data.
* 4. Path 2 - The file path to second FASTQ file. Only required for pair-ended data. Typically the R2 file in pair-ended data.
* 5. InputID - Used to pair sample and input files for various types of sequencing data. Use `-` if no input file is available or needed (as is the case in RNA-seq/GRO-seq/MNase-seq/etc.
* 6. Mark - Used to differentiate sample files from input files. Use the keyword `input` if that sample corresponds to an input file. Otherwise use `MergeID`.

## Running CIPHER

1) Install required dependencies

2) Create Singularity container (will require `sudo` access, so a container can be created on a local laptop/desktop and then transferred to the appopriate location/machine/cluster)

	```
	sudo singularity create -s 8000 cipher.img
	```

	```
	sudo singularity bootstrap cipher.img Singularity
	```

3) Run your workflow

	```
	nextflow run cipher.nf -with-singularity <cipher.img> --mode <MODE> --config <CONFIG> --fa <FASTA> --gtf <GTF> --lib <LIB> --readLen <LENGTH> [options]
	```

**NOTE:** If not running on a cluster please set the `-qs <INT>` flag in order to control the number of processes that CIPHER parallelizes. Too many and the workflow will abruptly end because it runs out of memory. `nextflow run -qs <INT> cipher.nf ...`

## Example Data

Some example data to test CIPHER's workflows can be found in the `example_data` folder. The user should alter the config file fastq paths before running the workflow otherwise the run will fail.

## Running CIPHER on a Cluster

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
