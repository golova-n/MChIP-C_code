#! /usr/bin/env Rscript
library("dplyr")
library("tidyr")
bedpe <- read.csv("tmp/MicroC_filtered.bedpe", header = F, sep = "\t")
bedpe <- filter(bedpe, V8>0 & (V5-V2)>5000 & (V5-V2)<1100000)
bedpe <- filter(bedpe,!(duplicated(bedpe) | duplicated(bedpe, fromLast = TRUE)))
bedgraph <- rbind(bedpe[,c(1:3,8)],setNames(bedpe[,c(4:6,8)],names(bedpe[,c(1:3,8)])))
bedgraph <- bedgraph %>% group_by(V1,V2,V3) %>% summarise(sum=sum(V8))
write.table(bedgraph, "tmp/MicroC.ICE.bedgraph", col.names = F, row.names = F, sep = "\t", quote = F)