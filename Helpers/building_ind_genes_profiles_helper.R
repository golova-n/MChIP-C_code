#! /usr/bin/env Rscript
library('dplyr')
library('tidyr')
library("GenomicRanges")

TSS <- read.csv("auxiliary_data/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed", header = F, col.names = c("chr","start","end","gene","V5","strand"), sep = "\t")[,c(1:4,6)]
vp <- read.csv("MChIPC_output/MChIPC_viewpoints.bed", header = F, sep = "\t")[,c(1:4)]
TSS_ranges <- GRanges(seqnames = TSS$chr, ranges=IRanges(start=TSS$start, end=TSS$end, enh_id = seq(1,nrow(TSS))))
vp_ranges <- GRanges(seqnames = vp$V1, ranges=IRanges(start=vp$V2, end=vp$V3, enh_id = seq(1,nrow(vp))))
overlaps <- as.data.frame(findOverlaps(vp_ranges,TSS_ranges))
vp$queryHits <- seq(nrow(vp))
TSS$subjectHits <- seq(nrow(TSS))
vp <- left_join(vp, overlaps) %>% left_join(TSS)
vp <- vp[,c(1:5,10)]%>% group_by(queryHits) %>% mutate(genes=paste0(gene, collapse=","))
vp <- vp[,c(1:3,7,4)] %>% distinct()
for (i in 1:nrow(vp)){
  if (vp[i,4]=="NA"){
    vp[i,4]<-paste0(vp[i,1],":",vp[i,2],"-",vp[i,3])}}
write.table(vp, file="MChIPC_output/MChIPC_viewpoints.bed", col.names = F, row.names = F, sep = "\t", quote = F)