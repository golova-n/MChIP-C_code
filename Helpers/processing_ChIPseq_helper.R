#! /usr/bin/env Rscript
library('dplyr')
library('tidyr')
DNase <- read.csv("Auxiliary_data/ENCFF621ZJY.bed.gz", header = F, sep = "\t")
DNase <- filter(DNase, V7 > 200)
write.table(DNase[,1:3], "tmp/DNase.filtered.bed", quote = F, col.names = F, row.names = F, sep = "\t")