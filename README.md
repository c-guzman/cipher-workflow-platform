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
