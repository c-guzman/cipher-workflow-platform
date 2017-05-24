RNA-seq Workflow
================

RNA-seq workflows require the following files:

1. A config file (described below)
2. Reference genome FASTA file
3. Reference genome GTF file
4. A experimental config file (described below)

*Config File*
A config file is a tab separated text file that includes information regarding the name, location, and input of your experiment.

**Single-End Config**
This is the config file format for single-ended data.

    sample1		sample1_rep1	/path/to/sample1_rep1.fastq.gz
    sample1		sample1_rep2	/path/to/sample1_rep2.fastq.gz
    sample1     sample1_rep3    /path/to/sample1_rep2.fastq.gz
    sample2		sample2_rep1	/path/to/sample2_rep1.fastq.gz
    sample2		sample2_rep2	/path/to/sample2_rep2.fastq.gz
    sample2     sample2_rep3    /path/to/sample1_rep2.fastq.gz

    ::

The columns represent:

1. MergeID: The merge ID that will be used should your files be merged together. Should be the same for all replicates.
2. ID: The ID that will be used to name the majority of your files that are not merged. Recommended to be used to differentiate between different technical replicates.
3. Path: The path to the fastq file to be processed. Can be gzipped or not.

**Pair-End Config**
This is the config file format for pair-ended data.

    sample1		sample1_rep1	/path/to/sample1_rep1_R1.fastq.gz /path/to/sample1_rep1_R2.fastq.gz
    sample1		sample1_rep2	/path/to/sample1_rep2_R1.fastq.gz /path/to/sample1_rep2_R2.fastq.gz
    sample1     sample1_rep3    /path/to/sample1_rep3_R1.fastq.gz /path/to/sample1_rep3_R2.fastq.gz
    sample2		sample2_rep1	/path/to/sample2_rep1_R1.fastq.gz /path/to/sample2_rep1_R2.fastq.gz
    sample2		sample2_rep2	/path/to/sample2_rep2_R1.fastq.gz /path/to/sample2_rep2_R2.fastq.gz
    sample2     sample2_rep3    /path/to/sample2_rep1_R1.fastq.gz /path/to/sample2_rep3_R2.fastq.gz

    ::

The columns represent:

1. MergeID: The merge ID that will be used should your files be merged together. Should be the same for all replicates.
2. ID: The ID that will be used to name the majority of your files that are not merged. Recommended to be used to differentiate between different technical replicates.
3. Path1: The path to the first fastq file to be processed. Can be gzipped or not.
4. Path2: The path to the second fastq file to be processed. Can be gzipped or not.

*Experimental Config File*
A experimental config file is a tab separated text file that has sample and condition information. Your experimental config file must be ordered in the same way that your config file is. (sample1 in config file == sample1 in expInfo file). Additionally you must include the headers in the expInfo config file.

   sample  condition
   ctrl1    WT
   ctrl2    WT
   ctrl3    WT
   ko1      KO
   ko2      KO
   ko3      KO

The headers "sample" and "condition" must be included.

The columns represent:

1. SampleID: Refers to the ID in your config file. Used to differentiate between different replicates.
2. Condition: Refers to the condition or experimental variable of your sample. CIPHER currently only supports two condition DGE analysis.

**Simple RNA-seq Tutorial (single-end, 75 length reads, frFirstStrand)**

    nextflow run /path/to/main.nf --mode mnase --config /path/to/config.txt --fasta /path/to/fasta.fa --gtf /path/to/gtf.gtf --lib s --readLen 75 --strandInfo frFirstStrand --expInfo exp_config.txt

    ::

**Simple RNA-seq Tutorial (pair-end, 75 length reads)**

    nextflow run /path/to/main.nf --mode mnase --config /path/to/config.txt --fasta /path/to/fasta.fa --gtf /path/to/gtf.gtf --lib p --readLen 75 --strandInfo frFirstStrand --expInfo exp_config.txt

    ::

**Simple RNA-seq Tutorial (single-end, 75 length reads, use star aligner instead of default bbmap, use 5 threads)**

    nextflow run /path/to/main.nf --mode mnase --config /path/to/config.txt --fasta /path/to/fasta.fa --gtf /path/to/gtf.gtf --lib s --readLen 75 --strandInfo frFirstStrand --expInfo exp_config.txt --aligner star --threads 5

    ::