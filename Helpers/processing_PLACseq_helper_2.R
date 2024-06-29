#! /usr/bin/env Rscript
library("dplyr")
library("tidyr")
anchors <- read.csv("tmp/PLACseq_anchors_hg19.bed", header = F, sep = "\t")
OE <- read.csv("tmp/PLACseq_OE_hg19.bed", header = F, sep = "\t")
all_loops <- cbind(anchors,OE)
write.table(all_loops, "PLACseq_output/PLACseq_interactions.bedpe", col.names = F, row.names = F, sep = "\t", quote = F)