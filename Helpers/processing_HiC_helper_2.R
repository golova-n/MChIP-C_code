#! /usr/bin/env Rscript
library("dplyr")
library("tidyr")
library("GenomicRanges")
mergeLoops <- function(loops_a, loops_b){
  loops_a_ranges1 <- GRanges(seqnames = loops_a$BIN1_CHR, ranges=IRanges(start=loops_a$BIN1_START, end=loops_a$BIN1_END, id = seq(1,nrow(loops_a))))
  loops_a_ranges2 <- GRanges(seqnames = loops_a$BIN2_CHR, ranges=IRanges(start=loops_a$BIN2_START, end=loops_a$BIN2_END, id = seq(1,nrow(loops_a))))
  loops_b_ranges1 <- GRanges(seqnames = loops_b$BIN1_CHR, ranges=IRanges(start=loops_b$BIN1_START, end=loops_b$BIN1_END, id = seq(1,nrow(loops_b))))
  loops_b_ranges2 <- GRanges(seqnames = loops_b$BIN2_CHR, ranges=IRanges(start=loops_b$BIN2_START, end=loops_b$BIN2_END, id = seq(1,nrow(loops_b))))
  rbind(loops_a[-(as.data.frame(findOverlaps(loops_a_ranges1, loops_b_ranges1, minoverlap = 2))[,1]
                  [split(as.data.frame(findOverlaps(loops_a_ranges1, loops_b_ranges1, minoverlap = 2)),seq(nrow(as.data.frame(findOverlaps(loops_a_ranges1, loops_b_ranges1, minoverlap = 2))))) %in% 
                      split(as.data.frame(findOverlaps(loops_a_ranges2, loops_b_ranges2, minoverlap = 2)),seq(nrow(as.data.frame(findOverlaps(loops_a_ranges2, loops_b_ranges2, minoverlap = 2)))))]),], loops_b)
}
loops_20kb <- read.csv("tmp/HiC_loops_20kb.tsv", header = T, sep = "\t")
loops_10kb <- read.csv("tmp/HiC_loops_10kb.tsv", header = T, sep = "\t")
loops_5kb <- read.csv("tmp/HiC_loops_5kb.tsv", header = T, sep = "\t")
loops_2kb <- read.csv("tmp/HiC_loops_2kb.tsv", header = T, sep = "\t")
loops_1kb <- read.csv("tmp/HiC_loops_1kb.tsv", header = T, sep = "\t")
loops <- mergeLoops(loops_1kb,loops_2kb)
loops <- mergeLoops(loops,loops_5kb)
loops <- mergeLoops(loops,loops_10kb)
loops <- mergeLoops(loops,loops_20kb)


write.table(loops, file="HiC_output/HiC_interactions.bedpe", col.names = F, row.names = F, sep = "\t", quote = F)