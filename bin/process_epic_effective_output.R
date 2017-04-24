#!/usr/bin/env Rscript

# Author: Carlos Guzman
# Process epic-effective output

args <- commandArgs(trailingOnly=TRUE)
egs_file <- args[-1]

data = read.table(egs_file,
                  header = F, sep = ":", stringsAsFactors = F)

data = trimws(data$V2)
egs_ratio = as.numeric(data[4])
data = as.matrix(as.numeric(data[c(2,4)]))

egs_size = round(data[1] * data[2])

write.table(egs_size, "egs_size.txt", quote = F, sep = '\t',
            row.names = F, col.names = F)

write.table(egs_ratio, "egs_ratio.txt", quote = F, sep ='\t',
            row.names = F, col.names = F)