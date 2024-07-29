#!/usr/bin/env bash

exec 3>&1 &>/dev/null
source $(conda info --base)/etc/profile.d/conda.sh
conda activate mchip-c

mkdir -p ChIPseq_output/mapping_stats
mkdir -p ChIPseq_output/sam
mkdir -p ChIPseq_output/profiles

echo "downloading H3K4me3 ChIP-seq fastq files" >&3
wget -q -P tmp/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23410270/SRR23410270
fasterq-dump -e 16 -O tmp/ tmp/SRR23410270
rm tmp/SRR23410270
pigz -p 16 tmp/*.fastq

echo "mapping and parsing ChIP-seq data" >&3
bwa mem -SP5M -t 16 Auxiliary_data/hg19/male.hg19.fa tmp/SRR23410270_1.fastq.gz tmp/SRR23410270_2.fastq.gz | pigz -p 16 > tmp/ChIP_H3K4me3.sam.gz
rm tmp/*.fastq.gz
pairtools parse --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -c Auxiliary_data/hg19/hg19.chrom.sizes --output-stats ChIPseq_output/mapping_stats/stats_map_ChIPseq.txt -o tmp/ChIP_H3K4me3.pairsam.gz tmp/ChIP_H3K4me3.sam.gz
pairtools sort --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --nproc 16 -o tmp/ChIP_H3K4me3.sorted.pairsam.gz tmp/ChIP_H3K4me3.pairsam.gz
pairtools select --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 1000) and (strand1=="+") and (strand2=="-"))' -o tmp/ChIP_H3K4me3.pairsam.gz tmp/ChIP_H3K4me3.sorted.pairsam.gz
pairtools dedup --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-stats ChIPseq_output/mapping_stats/stats_dedup_ChIPseq.txt -o tmp/ChIP_H3K4me3.dedup.pairsam.gz tmp/ChIP_H3K4me3.pairsam.gz 
pairtools split --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-sam tmp/ChIP_H3K4me3.sam.gz tmp/ChIP_H3K4me3.dedup.pairsam.gz
samtools view -@ 16 -bh tmp/ChIP_H3K4me3.sam.gz > tmp/ChIP_H3K4me3.bam
samtools sort -@ 16 tmp/ChIP_H3K4me3.bam > ChIPseq_output/sam/ChIP_H3K4me3.sorted.bam
rm tmp/*H3K4me3*

echo "generating H3K4me3 ChIP-seq profiles">&3
samtools index -@ 16 ChIPseq_output/sam/ChIP_H3K4me3.sorted.bam
bamCoverage -b ChIPseq_output/sam/ChIP_H3K4me3.sorted.bam -bs 50 -e -p 16 -o ChIPseq_output/profiles/ChIP_H3K4me3.bw

echo "calculating read coverage in DHSs for MChIP-C mononucleosomal data and ChIP-seq data" >&3
wget -P Auxiliary_data/ https://www.encodeproject.org/files/ENCFF621ZJY/@@download/ENCFF621ZJY.bed.gz
./Helpers/processing_ChIPseq_helper.R
multiBigwigSummary BED-file -p 16 -b MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep1.bw MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep2.bw MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep3.bw MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep4.bw ChIPseq_output/profiles/ChIP_H3K4me3.bw --outRawCounts MChIPC_output/mononucleosomal_and_ChIP/H3K4me3_coverage.txt --BED tmp/DNase.filtered.bed -o tmp/results.npz
rm tmp/results.npz
