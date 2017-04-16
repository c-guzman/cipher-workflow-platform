#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
# args[1] = bam
# args[2] = narrow peaks
# args[3] = id

# Load ChIPQC package
if (!require("ChIPQC")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("ChIPQC", suppressUpdates=TRUE)
  library("ChIPQC")
    }

sample = ChIPQCsample(args[1], peaks = args[2])
ChIPQCreport(sample, reportName=args[3], reportFolder=args[3])