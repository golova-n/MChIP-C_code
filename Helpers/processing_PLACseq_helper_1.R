#! /usr/bin/env Rscript
library("dplyr")
library("tidyr")
rep1_loops <- read.csv("tmp/GSE161873_format_K562wtr1.10k.2.peaks.bedpe.txt.gz", header = F, sep = "\t")
rep2_loops <- read.csv("tmp/GSE161873_format_K562wtr2.10k.2.peaks.bedpe.txt.gz", header = F, sep = "\t")
rep1_loops <- (separate(rep1_loops,V4, into=c("V4","V4d"),sep = ",") %>% separate(col=V4, into=c("V4","V4c"),sep = "-") %>% separate(col=V4, into=c("V4a","V4b"),sep = ":"))[,1:6]
rep2_loops <- (separate(rep2_loops,V4, into=c("V4","V4d"),sep = ",") %>% separate(col=V4, into=c("V4","V4c"),sep = "-") %>% separate(col=V4, into=c("V4a","V4b"),sep = ":"))[,1:6]
all_loops <- distinct(rbind(rep1_loops,rep2_loops)) %>% arrange(V1,V2,V4a,V4b)
# filtering out raws with un-liftable entries
all_loops <-filter(all_loops, !((V1=="chr1" & V2==205980000 & V3==205990000)|(V1=="chr1" & V2==235060000 & V3==235070000)|(V1=="chr17" & V2==64420000 & V3==64430000)|(V1=="chr20" & V2==29870000 & V3==29880000)|(V1=="chr20" & V2==30810000 & V3==30820000)|(V1=="chr6" & V2==106950000 & V3==106960000)|(V1=="chr8" & V2==144020000 & V3==144030000)))
all_loops <-filter(all_loops, !((V4a=="chr1" & V4b==205980000 & V4c==205990000)|(V4a=="chr1" & V4b==235060000 & V4c==235070000)|(V4a=="chr17" & V4b==64420000 & V4c==64430000)|(V4a=="chr6" & V4b==106950000 & V4c==106960000)|(V4a=="chr8" & V4b==144020000 & V4c==144030000)))
write.table(all_loops[,1:3], "tmp/PLACseq_anchors_hg38.bed", col.names = F, row.names = F, sep = "\t", quote = F)
write.table(all_loops[,4:6], "tmp/PLACseq_OE_hg38.bed", col.names = F, row.names = F, sep = "\t", quote = F)