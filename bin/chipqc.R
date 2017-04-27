#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
# args[1] = bam
# args[2] = narrow peaks
# args[3] = id

# Load ChIPQC package
library(ChIPQC)

bam = args[1]
peak_file = args[2]
id = args[3]

sample = ChIPQCsample(bam, peaks = peak_file)
ChIPQCreport(sample, reportName=id, reportFolder=id)