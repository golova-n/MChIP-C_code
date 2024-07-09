#! /usr/bin/env Rscript
system("conda activate mchip-c", intern=T) 
options(warn=-1)
suppressMessages({
library('dplyr')
library('tidyr')
library('fitdistrplus')
library("GenomicRanges")})
message("processing combined dataset")
message("\tloading the data")
suppressMessages({
# loading contact sum files and coverage files ('vp' stands for 'viewpoint', 'OE' - for 'Other End')
interactions_rep1 <- read.csv("tmp/interactions.sum.rep1.txt", header = F, sep = "\t")
interactions_rep2 <- read.csv("tmp/interactions.sum.rep2.txt", header = F, sep = "\t")
interactions_rep3 <- read.csv("tmp/interactions.sum.rep3.txt", header = F, sep = "\t")
interactions_rep4 <- read.csv("tmp/interactions.sum.rep4.txt", header = F, sep = "\t")

coverage_OE_rep1 <- read.csv("tmp/genomic_bins_coverage_rep1.bed", header = F, sep = "\t")
coverage_vp_rep1 <- read.csv("tmp/viewpoints_coverage_rep1.bed", header = F, sep = "\t",  dec=",")
coverage_vp_rep1$V5 <- as.numeric(coverage_vp_rep1$V5)
coverage_vp_rep1 <- coverage_vp_rep1[,c(1,2,3,5)]
coverage_OE_rep2 <- read.csv("tmp/genomic_bins_coverage_rep2.bed", header = F, sep = "\t")
coverage_vp_rep2 <- read.csv("tmp/viewpoints_coverage_rep2.bed", header = F, sep = "\t",  dec=",")
coverage_vp_rep2$V5 <- as.numeric(coverage_vp_rep2$V5)
coverage_vp_rep2 <- coverage_vp_rep2[,c(1,2,3,5)]
coverage_OE_rep3 <- read.csv("tmp/genomic_bins_coverage_rep3.bed", header = F, sep = "\t")
coverage_vp_rep3 <- read.csv("tmp/viewpoints_coverage_rep3.bed", header = F, sep = "\t",  dec=",")
coverage_vp_rep3$V5 <- as.numeric(coverage_vp_rep3$V5)
coverage_vp_rep3 <- coverage_vp_rep3[,c(1,2,3,5)]
coverage_OE_rep4 <- read.csv("tmp/genomic_bins_coverage_rep4.bed", header = F, sep = "\t")
coverage_vp_rep4 <- read.csv("tmp/viewpoints_coverage_rep4.bed", header = F, sep = "\t",  dec=",")
coverage_vp_rep4$V5 <- as.numeric(coverage_vp_rep4$V5)
coverage_vp_rep4 <- coverage_vp_rep4[,c(1,2,3,5)]

colnames(interactions_rep1) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep1")
colnames(interactions_rep2) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep2")
colnames(interactions_rep3) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep3")
colnames(interactions_rep4) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep4")

colnames(coverage_OE_rep1) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep1")
colnames(coverage_vp_rep1) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep1")
colnames(coverage_OE_rep2) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep2")
colnames(coverage_vp_rep2) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep2")
colnames(coverage_OE_rep3) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep3")
colnames(coverage_vp_rep3) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep3")
colnames(coverage_OE_rep4) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep4")
colnames(coverage_vp_rep4) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep4")

# combining all data in one df
interactions <- full_join(interactions_rep1, interactions_rep2)
interactions <- full_join(interactions, interactions_rep3)
interactions <- full_join(interactions, interactions_rep4)
interactions[is.na(interactions)] <- 0
interactions <- interactions %>% arrange(vp_chr,vp_start,OE_start)
interactions$N_sum <- interactions$N_rep1 + interactions$N_rep2 + interactions$N_rep3 + interactions$N_rep4

interactions <- left_join(interactions,coverage_OE_rep1)
interactions <- left_join(interactions,coverage_OE_rep2)
interactions <- left_join(interactions,coverage_OE_rep3)
interactions <- left_join(interactions,coverage_OE_rep4)
interactions <- left_join(interactions,coverage_vp_rep1)
interactions <- left_join(interactions,coverage_vp_rep2)
interactions <- left_join(interactions,coverage_vp_rep3)
interactions <- left_join(interactions,coverage_vp_rep4)

rm(list=ls()[grep("rep", ls())])})

message("\tfiltering the data")
# assessing veiwpoint coverage and filtering poorly covered veiwpoints - in future this should be done after normalization, but now we have strange poorly ChIP covered veiwpoints
suppressMessages({
vp_coverage <- group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage = sum(N_sum))
vp <- read.csv("MChIPC_output/MChIPC_viewpoints.bed", header = F, sep = "\t")
colnames(vp)[1:3] <- colnames(vp_coverage)[1:3]
vp <- left_join(vp, vp_coverage)
write.table(vp, file="MChIPC_output/MChIPC_viewpoints.bed", col.names = F, row.names = F, sep = "\t", quote = F)
rm(vp)
interactions <- left_join(interactions, group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage = sum(N_sum)))
interactions <- filter(interactions, vp_coverage>=1000)

# creating dist_rank column and filtering by distance
interactions$dist_rank <- 0
interactions$dist_rank[interactions$vp_start>interactions$OE_start] <- (interactions$vp_start[interactions$vp_start>interactions$OE_start] - interactions$OE_start[interactions$vp_start>interactions$OE_start])/250
interactions$dist_rank[interactions$vp_end<interactions$OE_start] <- (interactions$OE_end[interactions$vp_end<interactions$OE_start] - interactions$vp_end[interactions$vp_end<interactions$OE_start])/250
interactions <- filter(interactions, dist_rank>=20 & dist_rank<=10000)})

message("\tnormalizing the data")
# ChIP-signal normalization (want the sum to remain the same)
suppressMessages({
interactions$N_rep1_norm <- interactions$N_rep1 / (interactions$vp_cov_rep1 + interactions$OE_cov_rep1)
interactions$N_rep1_norm <- interactions$N_rep1_norm*sum(interactions$N_rep1)/sum(interactions$N_rep1_norm)
interactions$N_rep2_norm <- interactions$N_rep2 / (interactions$vp_cov_rep2 + interactions$OE_cov_rep2)
interactions$N_rep2_norm <- interactions$N_rep2_norm*sum(interactions$N_rep2)/sum(interactions$N_rep2_norm)
interactions$N_rep3_norm <- interactions$N_rep3 / (interactions$vp_cov_rep3 + interactions$OE_cov_rep3)
interactions$N_rep3_norm <- interactions$N_rep3_norm*sum(interactions$N_rep3)/sum(interactions$N_rep3_norm)
interactions$N_rep4_norm <- interactions$N_rep4 / (interactions$vp_cov_rep4 + interactions$OE_cov_rep4)
interactions$N_rep4_norm <- interactions$N_rep4_norm*sum(interactions$N_rep4)/sum(interactions$N_rep4_norm)

# summarising replicates
interactions$N_sum_norm <- interactions$N_rep1_norm + interactions$N_rep2_norm + interactions$N_rep3_norm + interactions$N_rep4_norm

#veiwpoint coverage after normalization
vp_coverage <- group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage = sum(N_sum))
vp_coverage <- left_join(vp_coverage, group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage_norm= sum(N_sum_norm)))
interactions <- left_join(interactions, vp_coverage[,c(1,2,3,5)])
interactions$N_sum_double_norm <- interactions$N_sum_norm*median(interactions$vp_coverage_norm) / interactions$vp_coverage_norm
vp_coverage <- left_join(vp_coverage, group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage_double_norm= sum(N_sum_double_norm)))})

message("\tbuilding background models, calculating p-values and calling interactions")
#building models, calculating p-values
interactions$p_val <- NA
interactions$p_val <- as.numeric(interactions$p_val)
for (i in seq(20,4000)){
  interactions_rank <- filter(interactions, dist_rank==i)
  x <- fitdistr(interactions_rank$N_sum_double_norm[interactions_rank$N_sum_double_norm < quantile(interactions_rank$N_sum_double_norm, 0.95)], "weibull")
  interactions_rank$p_val <- 1-pweibull(interactions_rank$N_sum_double_norm, shape = x$estimate["shape"], scale=x$estimate["scale"])
  interactions <- interactions %>% rows_update(interactions_rank, by=c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end"))}

message("\twriting output files")
# writing full interaction statistics file
write.table(interactions, file="MChIPC_output/MChIPC_interaction_data.txt", col.names = T, row.names = F, sep = "\t", quote = F)
# true loops p_val<0.01 & 6+ ligation products | interactions with 6+ total ligation products (absolutely arbitrary) further than 1Mbp from viewpoint
loops <- filter(interactions, (p_val<=0.01 & N_sum > 5) | (dist_rank > 4000 & N_sum > 5))
write.table(loops[,1:6], file="MChIPC_output/MChIPC_interactions.bedpe", col.names = F, row.names = F, sep = "\t", quote = F)

message("\tcalling interactions in individual replicates and writing output files")
# also  calliing loops in individual replicates:
suppressMessages({
vp_coverage <- left_join(vp_coverage, group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage_1_norm= sum(N_rep1_norm),
                      vp_coverage_2_norm= sum(N_rep2_norm),vp_coverage_3_norm= sum(N_rep3_norm),vp_coverage_4_norm= sum(N_rep4_norm)))
interactions <- left_join(interactions, vp_coverage[,c(1,2,3,7:10)])
interactions$N_sum_double_norm_1 <- interactions$N_rep1_norm*median(interactions$vp_coverage_1_norm) / interactions$vp_coverage_1_norm
interactions$N_sum_double_norm_2 <- interactions$N_rep2_norm*median(interactions$vp_coverage_2_norm) / interactions$vp_coverage_2_norm
interactions$N_sum_double_norm_3 <- interactions$N_rep3_norm*median(interactions$vp_coverage_3_norm) / interactions$vp_coverage_3_norm
interactions$N_sum_double_norm_4 <- interactions$N_rep4_norm*median(interactions$vp_coverage_4_norm) / interactions$vp_coverage_4_norm
interactions$p_val <- NA
interactions$p_val <- as.numeric(interactions$p_val)})

message("\t\tanalysing replicate 1")
interactions_sub <- filter(interactions, N_sum_double_norm_1>0)[,c(1:6,7,21,29,34)]
for (i in seq(20,4000)){
  interactions_rank <- filter(interactions_sub, dist_rank==i)
  x <- fitdistr(interactions_rank$N_sum_double_norm_1[interactions_rank$N_sum_double_norm_1 < quantile(interactions_rank$N_sum_double_norm_1, 0.95)], "weibull")
  interactions_rank$p_val <- 1-pweibull(interactions_rank$N_sum_double_norm_1, shape = x$estimate["shape"], scale=x$estimate["scale"])
  interactions_sub <- interactions_sub %>% rows_update(interactions_rank, by=c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end"))}
loops <- filter(interactions_sub, (p_val<=0.01 & N_rep1 > 3) | (dist_rank > 4000 & N_rep1 > 3))
write.table(loops[,1:6], file="tmp/MChIPC_interactions_rep1.bedpe", col.names = F, row.names = F, sep = "\t", quote = F)

message("\t\tanalysing replicate 2")
interactions_sub <- filter(interactions, N_sum_double_norm_2>0)[,c(1:6,8,21,29,35)]
for (i in seq(20,4000)){
  interactions_rank <- filter(interactions_sub, dist_rank==i)
  x <- fitdistr(interactions_rank$N_sum_double_norm_2[interactions_rank$N_sum_double_norm_2 < quantile(interactions_rank$N_sum_double_norm_2, 0.95)], "weibull")
  interactions_rank$p_val <- 1-pweibull(interactions_rank$N_sum_double_norm_2, shape = x$estimate["shape"], scale=x$estimate["scale"])
  interactions_sub <- interactions_sub %>% rows_update(interactions_rank, by=c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end"))}
loops <- filter(interactions_sub, (p_val<=0.01 & N_rep2 > 3) | (dist_rank > 4000 & N_rep2 > 3))
write.table(loops[,1:6], file="tmp/MChIPC_interactions_rep2.bedpe", col.names = F, row.names = F, sep = "\t", quote = F)

message("\t\tanalysing replicate 3")
interactions_sub <- filter(interactions, N_sum_double_norm_3>0)[,c(1:6,9,21,29,36)]
for (i in seq(20,4000)){
  interactions_rank <- filter(interactions_sub, dist_rank==i)
  x <- fitdistr(interactions_rank$N_sum_double_norm_3[interactions_rank$N_sum_double_norm_3 < quantile(interactions_rank$N_sum_double_norm_3, 0.95)], "weibull")
  interactions_rank$p_val <- 1-pweibull(interactions_rank$N_sum_double_norm_3, shape = x$estimate["shape"], scale=x$estimate["scale"])
  interactions_sub <- interactions_sub %>% rows_update(interactions_rank, by=c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end"))}
loops <- filter(interactions_sub, (p_val<=0.01 & N_rep3 > 3) | (dist_rank > 4000 & N_rep3 > 3))
write.table(loops[,1:6], file="tmp/MChIPC_interactions_rep3.bedpe", col.names = F, row.names = F, sep = "\t", quote = F)

message("\t\tanalysing replicate 4")
interactions_sub <- filter(interactions, N_sum_double_norm_4>0)[,c(1:6,10,21,29,37)]
for (i in seq(20,4000)){
  interactions_rank <- filter(interactions_sub, dist_rank==i)
  x <- fitdistr(interactions_rank$N_sum_double_norm_4[interactions_rank$N_sum_double_norm_4 < quantile(interactions_rank$N_sum_double_norm_4, 0.95)], "weibull")
  interactions_rank$p_val <- 1-pweibull(interactions_rank$N_sum_double_norm_4, shape = x$estimate["shape"], scale=x$estimate["scale"])
  interactions_sub <- interactions_sub %>% rows_update(interactions_rank, by=c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end"))}
loops <- filter(interactions_sub, (p_val<=0.01 & N_rep4 > 3) | (dist_rank > 4000 & N_rep4 > 3))
write.table(loops[,1:6], file="tmp/MChIPC_interactions_rep4.bedpe", col.names = F, row.names = F, sep = "\t", quote = F)
rm(list=ls())

message("\t\tcalculating proportions of replicate-specific interactions")
# calculating proportion of loops specific for each replicate
suppressMessages({
loops_1 <- read.csv("tmp/MChIPC_interactions_rep1.bedpe", header = F, sep = "\t")
loops_2 <- read.csv("tmp/MChIPC_interactions_rep2.bedpe", header = F, sep = "\t")
loops_3 <- read.csv("tmp/MChIPC_interactions_rep3.bedpe", header = F, sep = "\t")
loops_4 <- read.csv("tmp/MChIPC_interactions_rep4.bedpe", header = F, sep = "\t")
loop_1_ranges <- GRanges(seqnames = loops_1$V4, ranges=IRanges(start=loops_1$V5, end=loops_1$V6, enh_id = seq(1,nrow(loops_1))))
loop_2_ranges <- GRanges(seqnames = loops_2$V4, ranges=IRanges(start=loops_2$V5, end=loops_2$V6, enh_id = seq(1,nrow(loops_2))))
loop_3_ranges <- GRanges(seqnames = loops_3$V4, ranges=IRanges(start=loops_3$V5, end=loops_3$V6, enh_id = seq(1,nrow(loops_3))))
loop_4_ranges <- GRanges(seqnames = loops_4$V4, ranges=IRanges(start=loops_4$V5, end=loops_4$V6, enh_id = seq(1,nrow(loops_4))))
loops_merged <- full_join(loops_1,loops_2)
loops_merged <- full_join(loops_merged,loops_3)
loops_merged <- full_join(loops_merged,loops_4)
loop_merged_ranges <- GRanges(seqnames = loops_merged$V4, ranges=IRanges(start=loops_merged$V5-999, end=loops_merged$V6+999, enh_id = seq(1,nrow(loops_merged))))
loops_merged$rep1 <-FALSE
loops_merged$rep2 <-FALSE
loops_merged$rep3 <-FALSE
loops_merged$rep4 <-FALSE
loops_merged$rep1[unique(as.data.frame(findOverlaps(loop_merged_ranges, loop_1_ranges))[paste(loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_1_ranges))[,1],1], 
                 loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_1_ranges))[,1],2], loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_1_ranges))[,1],3]) == paste(loops_1[as.data.frame(findOverlaps(loop_merged_ranges, loop_1_ranges))[,2],1], 
                 loops_1[as.data.frame(findOverlaps(loop_merged_ranges, loop_1_ranges))[,2],2], loops_1[as.data.frame(findOverlaps(loop_merged_ranges, loop_1_ranges))[,2],3]),1])] <- TRUE
loops_merged$rep2[unique(as.data.frame(findOverlaps(loop_merged_ranges, loop_2_ranges))[paste(loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_2_ranges))[,1],1],
                 loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_2_ranges))[,1],2], loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_2_ranges))[,1],3]) == paste(loops_2[as.data.frame(findOverlaps(loop_merged_ranges, loop_2_ranges))[,2],1],
                 loops_2[as.data.frame(findOverlaps(loop_merged_ranges, loop_2_ranges))[,2],2], loops_2[as.data.frame(findOverlaps(loop_merged_ranges, loop_2_ranges))[,2],3]),1])] <- TRUE
loops_merged$rep3[unique(as.data.frame(findOverlaps(loop_merged_ranges, loop_3_ranges))[paste(loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_3_ranges))[,1],1], 
                 loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_3_ranges))[,1],2], loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_3_ranges))[,1],3]) == paste(loops_3[as.data.frame(findOverlaps(loop_merged_ranges, loop_3_ranges))[,2],1], 
                 loops_3[as.data.frame(findOverlaps(loop_merged_ranges, loop_3_ranges))[,2],2], loops_3[as.data.frame(findOverlaps(loop_merged_ranges, loop_3_ranges))[,2],3]),1])] <- TRUE
loops_merged$rep4[unique(as.data.frame(findOverlaps(loop_merged_ranges, loop_4_ranges))[paste(loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_4_ranges))[,1],1], 
                 loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_4_ranges))[,1],2], loops_merged[as.data.frame(findOverlaps(loop_merged_ranges, loop_4_ranges))[,1],3]) == paste(loops_4[as.data.frame(findOverlaps(loop_merged_ranges, loop_4_ranges))[,2],1],
                 loops_4[as.data.frame(findOverlaps(loop_merged_ranges, loop_4_ranges))[,2],2], loops_4[as.data.frame(findOverlaps(loop_merged_ranges, loop_4_ranges))[,2],3]),1])] <- TRUE

loops_1$rep1_perf <- TRUE
loops_2$rep2_perf <- TRUE
loops_3$rep3_perf <- TRUE
loops_4$rep4_perf <- TRUE
loops_merged <- left_join(loops_merged, loops_1)
loops_merged <- left_join(loops_merged, loops_2)
loops_merged <- left_join(loops_merged, loops_3)
loops_merged <- left_join(loops_merged, loops_4)
loops_merged[is.na(loops_merged)] <- FALSE

N1 <- nrow(filter(loops_merged, rep1_perf&(rep2|rep3|rep4)))*100/sum(loops_merged$rep1_perf)
N2 <- nrow(filter(loops_merged, rep2_perf&(rep1|rep3|rep4)))*100/sum(loops_merged$rep2_perf)
N3 <- nrow(filter(loops_merged, rep3_perf&(rep1|rep2|rep4)))*100/sum(loops_merged$rep3_perf)
N4 <- nrow(filter(loops_merged, rep4_perf&(rep1|rep2|rep3)))*100/sum(loops_merged$rep4_perf)})

message("\t\t",paste0(round(N1,2)),"% of loops identified in rep1 can be detected in at least one other replicate")
message("\t\t",paste0(round(N2,2)),"% of loops identified in rep2 can be detected in at least one other replicate")
message("\t\t",paste0(round(N3,2)),"% of loops identified in rep3 can be detected in at least one other replicate")
message("\t\t",paste0(round(N4,2)),"% of loops identified in rep4 can be detected in at least one other replicate")

rm(list=ls())

message("processing downsampled datasets")
message("processing 10% subsampled dataset")
message("\tloading the data")
# loading contact sum files and coverage files ('vp' stands for 'viewpoint', 'OE' - for 'Other End')
suppressMessages({
interactions_rep1 <- read.csv("tmp/interactions.sum.rep1.sub_0.1.txt", header = F, sep = "\t")
interactions_rep2 <- read.csv("tmp/interactions.sum.rep2.sub_0.1.txt", header = F, sep = "\t")
interactions_rep3 <- read.csv("tmp/interactions.sum.rep3.sub_0.1.txt", header = F, sep = "\t")
interactions_rep4 <- read.csv("tmp/interactions.sum.rep4.sub_0.1.txt", header = F, sep = "\t")

coverage_OE_rep1 <- read.csv("tmp/genomic_bins_coverage_rep1.sub_0.1.bed", header = F, sep = "\t")
coverage_vp_rep1 <- read.csv("tmp/viewpoints_coverage_rep1.sub_0.1.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep1 <- coverage_vp_rep1[,c(1,2,3,5)]
coverage_OE_rep2 <- read.csv("tmp/genomic_bins_coverage_rep2.sub_0.1.bed", header = F, sep = "\t")
coverage_vp_rep2 <- read.csv("tmp/viewpoints_coverage_rep2.sub_0.1.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep2 <- coverage_vp_rep2[,c(1,2,3,5)]
coverage_OE_rep3 <- read.csv("tmp/genomic_bins_coverage_rep3.sub_0.1.bed", header = F, sep = "\t")
coverage_vp_rep3 <- read.csv("tmp/viewpoints_coverage_rep3.sub_0.1.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep3 <- coverage_vp_rep3[,c(1,2,3,5)]
coverage_OE_rep4 <- read.csv("tmp/genomic_bins_coverage_rep4.sub_0.1.bed", header = F, sep = "\t")
coverage_vp_rep4 <- read.csv("tmp/viewpoints_coverage_rep4.sub_0.1.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep4 <- coverage_vp_rep4[,c(1,2,3,5)]

colnames(interactions_rep1) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep1")
colnames(interactions_rep2) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep2")
colnames(interactions_rep3) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep3")
colnames(interactions_rep4) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep4")

colnames(coverage_OE_rep1) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep1")
colnames(coverage_vp_rep1) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep1")
colnames(coverage_OE_rep2) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep2")
colnames(coverage_vp_rep2) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep2")
colnames(coverage_OE_rep3) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep3")
colnames(coverage_vp_rep3) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep3")
colnames(coverage_OE_rep4) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep4")
colnames(coverage_vp_rep4) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep4")

# combining all data in one df
interactions <- full_join(interactions_rep1, interactions_rep2)
interactions <- full_join(interactions, interactions_rep3)
interactions <- full_join(interactions, interactions_rep4)
interactions[is.na(interactions)] <- 0
interactions <- interactions %>% arrange(vp_chr,vp_start,OE_start)
interactions$N_sum <- interactions$N_rep1 + interactions$N_rep2 + interactions$N_rep3 + interactions$N_rep4

interactions <- left_join(interactions,coverage_OE_rep1)
interactions <- left_join(interactions,coverage_OE_rep2)
interactions <- left_join(interactions,coverage_OE_rep3)
interactions <- left_join(interactions,coverage_OE_rep4)
interactions <- left_join(interactions,coverage_vp_rep1)
interactions <- left_join(interactions,coverage_vp_rep2)
interactions <- left_join(interactions,coverage_vp_rep3)
interactions <- left_join(interactions,coverage_vp_rep4)

rm(list=ls()[grep("rep", ls())])})

message("\tfiltering the data")
suppressMessages({
vp <- read.csv("MChIPC_output/MChIPC_viewpoints.bed", header=F, sep = "\t")
#should replace these two lines later
vp <- vp[,c(1,2,3,5)]
colnames(vp) <- c("vp_chr","vp_start","vp_end","vp_coverage")
interactions <- left_join(interactions, vp)
interactions <- filter(interactions, vp_coverage>=1000)

# creating dist_rank column and filtering by distance
interactions$dist_rank <- 0
interactions$dist_rank[interactions$vp_start>interactions$OE_start] <- (interactions$vp_start[interactions$vp_start>interactions$OE_start] - interactions$OE_start[interactions$vp_start>interactions$OE_start])/250
interactions$dist_rank[interactions$vp_end<interactions$OE_start] <- (interactions$OE_end[interactions$vp_end<interactions$OE_start] - interactions$vp_end[interactions$vp_end<interactions$OE_start])/250
interactions <- filter(interactions, dist_rank>=20 & dist_rank<=10000)})

message("\tnormalizing the data")
suppressMessages({
# ChIP-signal normalization (want the sum to remain the same)
interactions$N_rep1_norm <- interactions$N_rep1 / (interactions$vp_cov_rep1 + interactions$OE_cov_rep1)
interactions$N_rep1_norm <- interactions$N_rep1_norm*sum(interactions$N_rep1)/sum(interactions$N_rep1_norm)
interactions$N_rep2_norm <- interactions$N_rep2 / (interactions$vp_cov_rep2 + interactions$OE_cov_rep2)
interactions$N_rep2_norm <- interactions$N_rep2_norm*sum(interactions$N_rep2)/sum(interactions$N_rep2_norm)
interactions$N_rep3_norm <- interactions$N_rep3 / (interactions$vp_cov_rep3 + interactions$OE_cov_rep3)
interactions$N_rep3_norm <- interactions$N_rep3_norm*sum(interactions$N_rep3)/sum(interactions$N_rep3_norm)
interactions$N_rep4_norm <- interactions$N_rep4 / (interactions$vp_cov_rep4 + interactions$OE_cov_rep4)
interactions$N_rep4_norm <- interactions$N_rep4_norm*sum(interactions$N_rep4)/sum(interactions$N_rep4_norm)

# summarising replicates
interactions$N_sum_norm <- interactions$N_rep1_norm + interactions$N_rep2_norm + interactions$N_rep3_norm + interactions$N_rep4_norm

#veiwpoint coverage after normalization
vp_coverage <- group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage = sum(N_sum))
vp_coverage <- left_join(vp_coverage, group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage_norm= sum(N_sum_norm)))
interactions <- left_join(interactions, vp_coverage[,c(1,2,3,5)])
interactions$N_sum_double_norm <- interactions$N_sum_norm*median(interactions$vp_coverage_norm) / interactions$vp_coverage_norm
vp_coverage <- left_join(vp_coverage, group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage_double_norm= sum(N_sum_double_norm)))})

message("\tbuilding background models, calculating p-values and calling interactions")
#building models, calculating p-values
interactions$p_val <- NA
interactions$p_val <- as.numeric(interactions$p_val)
for (i in seq(20,4000)){
  interactions_rank <- filter(interactions, dist_rank==i)
  x <- fitdistr(interactions_rank$N_sum_double_norm[interactions_rank$N_sum_double_norm < quantile(interactions_rank$N_sum_double_norm, 0.95)], "weibull")
  interactions_rank$p_val <- 1-pweibull(interactions_rank$N_sum_double_norm, shape = x$estimate["shape"], scale=x$estimate["scale"])
  interactions <- interactions %>% rows_update(interactions_rank, by=c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end"))}

message("\twriting output files")

# true loops p_val<0.01 & 3+ ligation products | interactions with 6+ total ligation products (absolutely arbitrary) further than 1Mbp from viewpoint
loops <- filter(interactions, (p_val<=0.01 & N_sum > 3) | (dist_rank > 4000 & N_sum > 3))
write.table(loops[,1:6], file="MChIPC_output/MChIPC_interactions.sub_0.1.bedpe", col.names = F, row.names = F, sep = "\t", quote = F)

rm(list=ls())
message("processing 25% subsampled dataset")
message("\tloading the data")
# loading contact sum files and coverage files ('vp' stands for 'viewpoint', 'OE' - for 'Other End')
suppressMessages({
interactions_rep1 <- read.csv("tmp/interactions.sum.rep1.sub_0.25.txt", header = F, sep = "\t")
interactions_rep2 <- read.csv("tmp/interactions.sum.rep2.sub_0.25.txt", header = F, sep = "\t")
interactions_rep3 <- read.csv("tmp/interactions.sum.rep3.sub_0.25.txt", header = F, sep = "\t")
interactions_rep4 <- read.csv("tmp/interactions.sum.rep4.sub_0.25.txt", header = F, sep = "\t")

coverage_OE_rep1 <- read.csv("tmp/genomic_bins_coverage_rep1.sub_0.25.bed", header = F, sep = "\t")
coverage_vp_rep1 <- read.csv("tmp/viewpoints_coverage_rep1.sub_0.25.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep1 <- coverage_vp_rep1[,c(1,2,3,5)]
coverage_OE_rep2 <- read.csv("tmp/genomic_bins_coverage_rep2.sub_0.25.bed", header = F, sep = "\t")
coverage_vp_rep2 <- read.csv("tmp/viewpoints_coverage_rep2.sub_0.25.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep2 <- coverage_vp_rep2[,c(1,2,3,5)]
coverage_OE_rep3 <- read.csv("tmp/genomic_bins_coverage_rep3.sub_0.25.bed", header = F, sep = "\t")
coverage_vp_rep3 <- read.csv("tmp/viewpoints_coverage_rep3.sub_0.25.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep3 <- coverage_vp_rep3[,c(1,2,3,5)]
coverage_OE_rep4 <- read.csv("tmp/genomic_bins_coverage_rep4.sub_0.25.bed", header = F, sep = "\t")
coverage_vp_rep4 <- read.csv("tmp/viewpoints_coverage_rep4.sub_0.25.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep4 <- coverage_vp_rep4[,c(1,2,3,5)]

colnames(interactions_rep1) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep1")
colnames(interactions_rep2) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep2")
colnames(interactions_rep3) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep3")
colnames(interactions_rep4) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep4")

colnames(coverage_OE_rep1) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep1")
colnames(coverage_vp_rep1) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep1")
colnames(coverage_OE_rep2) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep2")
colnames(coverage_vp_rep2) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep2")
colnames(coverage_OE_rep3) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep3")
colnames(coverage_vp_rep3) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep3")
colnames(coverage_OE_rep4) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep4")
colnames(coverage_vp_rep4) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep4")

# combining all data in one df
interactions <- full_join(interactions_rep1, interactions_rep2)
interactions <- full_join(interactions, interactions_rep3)
interactions <- full_join(interactions, interactions_rep4)
interactions[is.na(interactions)] <- 0
interactions <- interactions %>% arrange(vp_chr,vp_start,OE_start)
interactions$N_sum <- interactions$N_rep1 + interactions$N_rep2 + interactions$N_rep3 + interactions$N_rep4

interactions <- left_join(interactions,coverage_OE_rep1)
interactions <- left_join(interactions,coverage_OE_rep2)
interactions <- left_join(interactions,coverage_OE_rep3)
interactions <- left_join(interactions,coverage_OE_rep4)
interactions <- left_join(interactions,coverage_vp_rep1)
interactions <- left_join(interactions,coverage_vp_rep2)
interactions <- left_join(interactions,coverage_vp_rep3)
interactions <- left_join(interactions,coverage_vp_rep4)

rm(list=ls()[grep("rep", ls())])})

message("\tfiltering the data")
suppressMessages({
vp <- read.csv("MChIPC_output/MChIPC_viewpoints.bed", header=F, sep = "\t")
#should replace these two lines later
vp <- vp[,c(1,2,3,5)]
colnames(vp) <- c("vp_chr","vp_start","vp_end","vp_coverage")
interactions <- left_join(interactions, vp)
interactions <- filter(interactions, vp_coverage>=1000)

# creating dist_rank column and filtering by distance
interactions$dist_rank <- 0
interactions$dist_rank[interactions$vp_start>interactions$OE_start] <- (interactions$vp_start[interactions$vp_start>interactions$OE_start] - interactions$OE_start[interactions$vp_start>interactions$OE_start])/250
interactions$dist_rank[interactions$vp_end<interactions$OE_start] <- (interactions$OE_end[interactions$vp_end<interactions$OE_start] - interactions$vp_end[interactions$vp_end<interactions$OE_start])/250
interactions <- filter(interactions, dist_rank>=20 & dist_rank<=10000)})

message("\tnormalizing the data")
# ChIP-signal normalization (want the sum to remain the same)
suppressMessages({
interactions$N_rep1_norm <- interactions$N_rep1 / (interactions$vp_cov_rep1 + interactions$OE_cov_rep1)
interactions$N_rep1_norm <- interactions$N_rep1_norm*sum(interactions$N_rep1)/sum(interactions$N_rep1_norm)
interactions$N_rep2_norm <- interactions$N_rep2 / (interactions$vp_cov_rep2 + interactions$OE_cov_rep2)
interactions$N_rep2_norm <- interactions$N_rep2_norm*sum(interactions$N_rep2)/sum(interactions$N_rep2_norm)
interactions$N_rep3_norm <- interactions$N_rep3 / (interactions$vp_cov_rep3 + interactions$OE_cov_rep3)
interactions$N_rep3_norm <- interactions$N_rep3_norm*sum(interactions$N_rep3)/sum(interactions$N_rep3_norm)
interactions$N_rep4_norm <- interactions$N_rep4 / (interactions$vp_cov_rep4 + interactions$OE_cov_rep4)
interactions$N_rep4_norm <- interactions$N_rep4_norm*sum(interactions$N_rep4)/sum(interactions$N_rep4_norm)

# summarising replicates
interactions$N_sum_norm <- interactions$N_rep1_norm + interactions$N_rep2_norm + interactions$N_rep3_norm + interactions$N_rep4_norm

#veiwpoint coverage after normalization
vp_coverage <- group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage = sum(N_sum))
vp_coverage <- left_join(vp_coverage, group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage_norm= sum(N_sum_norm)))
interactions <- left_join(interactions, vp_coverage[,c(1,2,3,5)])
interactions$N_sum_double_norm <- interactions$N_sum_norm*median(interactions$vp_coverage_norm) / interactions$vp_coverage_norm
vp_coverage <- left_join(vp_coverage, group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage_double_norm= sum(N_sum_double_norm)))})

message("\tbuilding background models, calculating p-values and calling interactions")
#building models, calculating p-values
interactions$p_val <- NA
interactions$p_val <- as.numeric(interactions$p_val)
for (i in seq(20,4000)){
  interactions_rank <- filter(interactions, dist_rank==i)
  x <- fitdistr(interactions_rank$N_sum_double_norm[interactions_rank$N_sum_double_norm < quantile(interactions_rank$N_sum_double_norm, 0.95)], "weibull")
  interactions_rank$p_val <- 1-pweibull(interactions_rank$N_sum_double_norm, shape = x$estimate["shape"], scale=x$estimate["scale"])
  interactions <- interactions %>% rows_update(interactions_rank, by=c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end"))}

message("\twriting output files")

# true loops p_val<0.01 & 6+ ligation products | interactions with 6+ total ligation products (absolutely arbitrary) further than 1Mbp from viewpoint
loops <- filter(interactions, (p_val<=0.01 & N_sum > 3) | (dist_rank > 4000 & N_sum > 3))
write.table(loops[,1:6], file="MChIPC_output/MChIPC_interactions.sub_0.25.bedpe", col.names = F, row.names = F, sep = "\t", quote = F)


rm(list=ls())

message("processing 50% subsampled dataset")
message("\tloading the data")
# loading contact sum files and coverage files ('vp' stands for 'viewpoint', 'OE' - for 'Other End')
suppressMessages({
interactions_rep1 <- read.csv("tmp/interactions.sum.rep1.sub_0.5.txt", header = F, sep = "\t")
interactions_rep2 <- read.csv("tmp/interactions.sum.rep2.sub_0.5.txt", header = F, sep = "\t")
interactions_rep3 <- read.csv("tmp/interactions.sum.rep3.sub_0.5.txt", header = F, sep = "\t")
interactions_rep4 <- read.csv("tmp/interactions.sum.rep4.sub_0.5.txt", header = F, sep = "\t")

coverage_OE_rep1 <- read.csv("tmp/genomic_bins_coverage_rep1.sub_0.5.bed", header = F, sep = "\t")
coverage_vp_rep1 <- read.csv("tmp/viewpoints_coverage_rep1.sub_0.5.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep1 <- coverage_vp_rep1[,c(1,2,3,5)]
coverage_OE_rep2 <- read.csv("tmp/genomic_bins_coverage_rep2.sub_0.5.bed", header = F, sep = "\t")
coverage_vp_rep2 <- read.csv("tmp/viewpoints_coverage_rep2.sub_0.5.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep2 <- coverage_vp_rep2[,c(1,2,3,5)]
coverage_OE_rep3 <- read.csv("tmp/genomic_bins_coverage_rep3.sub_0.5.bed", header = F, sep = "\t")
coverage_vp_rep3 <- read.csv("tmp/viewpoints_coverage_rep3.sub_0.5.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep3 <- coverage_vp_rep3[,c(1,2,3,5)]
coverage_OE_rep4 <- read.csv("tmp/genomic_bins_coverage_rep4.sub_0.5.bed", header = F, sep = "\t")
coverage_vp_rep4 <- read.csv("tmp/viewpoints_coverage_rep4.sub_0.5.bed", header = F, sep = "\t",  dec=".")
coverage_vp_rep4 <- coverage_vp_rep4[,c(1,2,3,5)]

colnames(interactions_rep1) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep1")
colnames(interactions_rep2) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep2")
colnames(interactions_rep3) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep3")
colnames(interactions_rep4) <- c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end", "N_rep4")

colnames(coverage_OE_rep1) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep1")
colnames(coverage_vp_rep1) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep1")
colnames(coverage_OE_rep2) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep2")
colnames(coverage_vp_rep2) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep2")
colnames(coverage_OE_rep3) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep3")
colnames(coverage_vp_rep3) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep3")
colnames(coverage_OE_rep4) <- c("OE_chr","OE_start","OE_end", "OE_cov_rep4")
colnames(coverage_vp_rep4) <- c("vp_chr","vp_start","vp_end", "vp_cov_rep4")

# combining all data in one df
interactions <- full_join(interactions_rep1, interactions_rep2)
interactions <- full_join(interactions, interactions_rep3)
interactions <- full_join(interactions, interactions_rep4)
interactions[is.na(interactions)] <- 0
interactions <- interactions %>% arrange(vp_chr,vp_start,OE_start)
interactions$N_sum <- interactions$N_rep1 + interactions$N_rep2 + interactions$N_rep3 + interactions$N_rep4

interactions <- left_join(interactions,coverage_OE_rep1)
interactions <- left_join(interactions,coverage_OE_rep2)
interactions <- left_join(interactions,coverage_OE_rep3)
interactions <- left_join(interactions,coverage_OE_rep4)
interactions <- left_join(interactions,coverage_vp_rep1)
interactions <- left_join(interactions,coverage_vp_rep2)
interactions <- left_join(interactions,coverage_vp_rep3)
interactions <- left_join(interactions,coverage_vp_rep4)

rm(list=ls()[grep("rep", ls())])})

message("\tfiltering the data")
suppressMessages({
vp <- read.csv("MChIPC_output/MChIPC_viewpoints.bed", header=F, sep = "\t")
#should replace these two lines later
vp <- vp[,c(1,2,3,5)]
colnames(vp) <- c("vp_chr","vp_start","vp_end","vp_coverage")
interactions <- left_join(interactions, vp)
interactions <- filter(interactions, vp_coverage>=1000)

# creating dist_rank column and filtering by distance
interactions$dist_rank <- 0
interactions$dist_rank[interactions$vp_start>interactions$OE_start] <- (interactions$vp_start[interactions$vp_start>interactions$OE_start] - interactions$OE_start[interactions$vp_start>interactions$OE_start])/250
interactions$dist_rank[interactions$vp_end<interactions$OE_start] <- (interactions$OE_end[interactions$vp_end<interactions$OE_start] - interactions$vp_end[interactions$vp_end<interactions$OE_start])/250
interactions <- filter(interactions, dist_rank>=20 & dist_rank<=10000)})

message("\tnormalizing the data")
# ChIP-signal normalization (want the sum to remain the same)
suppressMessages({
interactions$N_rep1_norm <- interactions$N_rep1 / (interactions$vp_cov_rep1 + interactions$OE_cov_rep1)
interactions$N_rep1_norm <- interactions$N_rep1_norm*sum(interactions$N_rep1)/sum(interactions$N_rep1_norm)
interactions$N_rep2_norm <- interactions$N_rep2 / (interactions$vp_cov_rep2 + interactions$OE_cov_rep2)
interactions$N_rep2_norm <- interactions$N_rep2_norm*sum(interactions$N_rep2)/sum(interactions$N_rep2_norm)
interactions$N_rep3_norm <- interactions$N_rep3 / (interactions$vp_cov_rep3 + interactions$OE_cov_rep3)
interactions$N_rep3_norm <- interactions$N_rep3_norm*sum(interactions$N_rep3)/sum(interactions$N_rep3_norm)
interactions$N_rep4_norm <- interactions$N_rep4 / (interactions$vp_cov_rep4 + interactions$OE_cov_rep4)
interactions$N_rep4_norm <- interactions$N_rep4_norm*sum(interactions$N_rep4)/sum(interactions$N_rep4_norm)

# summarising replicates
interactions$N_sum_norm <- interactions$N_rep1_norm + interactions$N_rep2_norm + interactions$N_rep3_norm + interactions$N_rep4_norm

#veiwpoint coverage after normalization
vp_coverage <- group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage = sum(N_sum))
vp_coverage <- left_join(vp_coverage, group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage_norm= sum(N_sum_norm)))
interactions <- left_join(interactions, vp_coverage[,c(1,2,3,5)])
interactions$N_sum_double_norm <- interactions$N_sum_norm*median(interactions$vp_coverage_norm) / interactions$vp_coverage_norm
vp_coverage <- left_join(vp_coverage, group_by(interactions,vp_chr,vp_start,vp_end) %>% summarise(vp_coverage_double_norm= sum(N_sum_double_norm)))})

message("\tbuilding background models, calculating p-values and calling interactions")
#building models, calculating p-values
interactions$p_val <- NA
interactions$p_val <- as.numeric(interactions$p_val)
for (i in seq(20,4000)){
  interactions_rank <- filter(interactions, dist_rank==i)
  x <- fitdistr(interactions_rank$N_sum_double_norm[interactions_rank$N_sum_double_norm < quantile(interactions_rank$N_sum_double_norm, 0.95)], "weibull")
  interactions_rank$p_val <- 1-pweibull(interactions_rank$N_sum_double_norm, shape = x$estimate["shape"], scale=x$estimate["scale"])
  interactions <- interactions %>% rows_update(interactions_rank, by=c("vp_chr","vp_start","vp_end","OE_chr","OE_start","OE_end"))}

message("\twriting output files")

# true loops p_val<0.01 & 6+ ligation products | interactions with 6+ total ligation products (absolutely arbitrary) further than 1Mbp from viewpoint
loops <- filter(interactions, (p_val<=0.01 & N_sum > 3) | (dist_rank > 4000 & N_sum > 3))
write.table(loops[,1:6], file="MChIPC_output/MChIPC_interactions.sub_0.5.bedpe", col.names = F, row.names = F, sep = "\t", quote = F)
