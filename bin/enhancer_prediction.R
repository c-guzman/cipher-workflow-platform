#!/usr/bin/env Rscript

## Carlos Guzman
## Conducts enhancer-identification for CIPHER
## cguzman.roma@gmail.com
## Updated 11-3-16

## Load packages
library(plyr)
library(randomForest)
library(data.table)

## Set up variable to control command line arguments
## args[1] is the matrix file output by CIPHER
## args[2] is the model rds file
args <- commandArgs(TRUE)
model <- args[2]

genome <- fread(args[1])

genome <- rename(genome, c("V1"="Chr", "V2"="Start", "V3"="End", "V4"="DNase",
                           "V5"="H3K27Ac", "V6"="H3K4me1", "V7"="H3K4me3"))

genome <- genome[complete.cases(genome),]

genome <- subset(genome, subset = genome$DNase != 0 & genome$H3K27Ac != 0 & genome$H3K4me1 != 0 & genome$H3K4me3 != 0)

data.rf <- readRDS(model)
genome$predicted.response <- predict(data.rf, genome)

enhancers <- subset(x = genome, subset = genome$predicted.response == "Enhancers")
write.table(enhancers, file = "predicted.enhancers.txt", quote = F, row.names = F, col.names = F,sep = '\t')