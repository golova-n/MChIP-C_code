#!/usr/bin/env bash

exec 3>&1 &>/dev/null
source $(conda info --base)/etc/profile.d/conda.sh
conda activate mchip-c

echo "generating mononucleosomal profiles" >&3
for rep in rep1 rep2 rep3 rep4;
do
	samtools index -@ 16 MChIPC_output/sam/Mononucleosomal_$rep.sorted.bam
	bamCoverage -p 16 -b MChIPC_output/sam/Mononucleosomal_$rep.sorted.bam -bs 50 -e -o MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_$rep.bw
done
bigwigCompare -b1 MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep1.bw -b2 MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep2.bw --operation mean -p 16 -o tmp/Mononucleosomal_1_2.bw
bigwigCompare -b1 MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep3.bw -b2 MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep4.bw --operation mean -p 16 -o tmp/Mononucleosomal_3_4.bw
bigwigCompare -b1 tmp/Mononucleosomal_1_2.bw -b2 tmp/Mononucleosomal_3_4.bw --operation mean -p 16 -o MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_mean.bw
rm tmp/*.bw

echo "calling peaks in mononuclesomal profiles and identifying MChIP-C viewpoints" >&3
mkdir -p tmp/macs/
for rep in rep1 rep2 rep3 rep4;
do
	macs2 callpeak -t MChIPC_output/sam/Mononucleosomal_$rep.sorted.bam --outdir tmp/macs/ -n Mononucleosomal_$rep -f BAMPE -q 0.0001 --max-gap 1000
done
mv tmp/macs/*.narrowPeak MChIPC_output/mononucleosomal_and_ChIP
bedtools makewindows -g Auxiliary_data/hg19/hg19.chrom.sizes -w 250 > Auxiliary_data/hg19/genomic_bins.bed
# finding concensus peaks
bedtools multiinter -i MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep1_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep2_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep3_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep4_peaks.narrowPeak | awk '$4 >= 3 {print $1"\t"$2"\t"$3}' | bedtools merge -d 1000 -i stdin | bedtools slop -b 750 -i stdin -g Auxiliary_data/hg19/hg19.chrom.sizes > tmp/concensus_peaks.bed
# creating viewpoint file (binned_peaks)
bedtools intersect -wa -a Auxiliary_data/hg19/genomic_bins.bed -b tmp/concensus_peaks.bed | bedtools merge > MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed
# creating H3K4me3 mask
bedtools multiinter -i MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep1_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep2_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep3_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep4_peaks.narrowPeak | awk '{print $1"\t"$2"\t"$3}' | bedtools merge -i stdin > MChIPC_output/mononucleosomal_and_ChIP/H3K4me3.mask.bed
./Helpers/processing_MChIPC_helper.R

echo "generating raw and merged aggregate MChIP-C profiles" >&3
for rep in rep1 rep2 rep3 rep4;
do
	pairtools split --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-sam tmp/MChIPC_$rep.sam.gz --output-pairs MChIPC_output/pairs/MChIPC_$rep.pairs.gz MChIPC_output/sam/MChIPC_$rep.dedup.pairsam.gz
	pigz -p 16 -dc MChIPC_output/pairs/MChIPC_$rep.pairs.gz | grep -v '^#'| awk '{OFS="\t"}{print $2,$3,$3+1,$4,$5,$5+1}' > tmp/contacts.$rep.bedpe
	samtools view -@ 16 -bh tmp/MChIPC_$rep.sam.gz | samtools sort -@ 16 > MChIPC_output/sam/MChIPC_$rep.raw.bam
	rm tmp/MChIPC_$rep.sam.gz
	samtools index -@ 16 MChIPC_output/sam/MChIPC_$rep.raw.bam
	samtools view -b -P -h -@ 16 --region-file MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed MChIPC_output/sam/MChIPC_$rep.raw.bam > MChIPC_output/sam/MChIPC_$rep.merged.bam
	samtools index -@ 16 MChIPC_output/sam/MChIPC_$rep.merged.bam
	bamCoverage -p 16 -b MChIPC_output/sam/MChIPC_$rep.merged.bam -bs 250 -o MChIPC_output/aggregate_profiles/MChIPC_$rep.merged.250bp.bw
done
multiBigwigSummary bins -bs 250 -p 16 -b MChIPC_output/aggregate_profiles/MChIPC_rep1.merged.250bp.bw MChIPC_output/aggregate_profiles/MChIPC_rep2.merged.250bp.bw MChIPC_output/aggregate_profiles/MChIPC_rep3.merged.250bp.bw MChIPC_output/aggregate_profiles/MChIPC_rep4.merged.250bp.bw --outRawCounts MChIPC_output/aggregate_profiles/replicates_coverage.txt -o tmp/results.npz
rm tmp/results.npz
samtools merge -@ 16 -o MChIPC_output/sam/MChIPC_raw.bam MChIPC_output/sam/MChIPC_rep1.raw.bam MChIPC_output/sam/MChIPC_rep2.raw.bam MChIPC_output/sam/MChIPC_rep3.raw.bam MChIPC_output/sam/MChIPC_rep4.raw.bam
samtools index -@ 16 MChIPC_output/sam/MChIPC_raw.bam
bamCoverage -p 16 -b MChIPC_output/sam/MChIPC_raw.bam -bs 250 -o MChIPC_output/aggregate_profiles/MChIPC_raw.250bp.bw
bamCoverage -p 16 -b MChIPC_output/sam/MChIPC_raw.bam -bs 100 -o MChIPC_output/aggregate_profiles/MChIPC_raw.100bp.bw
samtools merge -@ 16 -o MChIPC_output/sam/MChIPC_merged.bam MChIPC_output/sam/MChIPC_rep1.merged.bam MChIPC_output/sam/MChIPC_rep2.merged.bam MChIPC_output/sam/MChIPC_rep3.merged.bam MChIPC_output/sam/MChIPC_rep4.merged.bam
samtools index -@ 16 MChIPC_output/sam/MChIPC_merged.bam
bamCoverage -p 16 -b MChIPC_output/sam/MChIPC_merged.bam -bs 250 -o MChIPC_output/aggregate_profiles/MChIPC_merged.250bp.bw
bamCoverage -p 16 -b MChIPC_output/sam/MChIPC_merged.bam -bs 100 -o MChIPC_output/aggregate_profiles/MChIPC_merged.100bp.bw

echo "creating cooler MChIP-C files" >&3
for rep in rep1 rep2 rep3 rep4;
do
	cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 Auxiliary_data/hg19/hg19.chrom.sizes:250 MChIPC_output/pairs/MChIPC_$rep.pairs.gz tmp/MChIPC_$rep.250bp.cool
	cooler zoomify --nproc 16 --out MChIPC_output/cools/MChIPC_$rep.mcool --resolutions 250,500,1000,2000,5000,10000,50000,100000,500000,1000000,5000000,10000000 tmp/MChIPC_$rep.250bp.cool
done
cooler merge tmp/MChIPC.250bp.cool tmp/MChIPC_rep1.250bp.cool tmp/MChIPC_rep2.250bp.cool tmp/MChIPC_rep3.250bp.cool tmp/MChIPC_rep4.250bp.cool
cooler zoomify --nproc 16 --out MChIPC_output/cools/MChIPC.mcool --resolutions 250,500,1000,2000,5000,10000,50000,100000,500000,1000000,5000000,10000000 tmp/MChIPC.250bp.cool
rm tmp/MChIPC*.cool

echo "preparing MChIP-C input files for loop calling" >&3
for rep in rep1 rep2 rep3 rep4;
do	
	for sample in 0.1 0.25 0.5;
	do
		pairtools split --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-sam tmp/MChIPC_$rep.sub_$sample.sam.gz --output-pairs tmp/MChIPC_$rep.sub_$sample.pairs.gz tmp/MChIPC_$rep.sub_$sample.dedup.pairsam.gz;
		pigz -p 16 -dc tmp/MChIPC_$rep.sub_$sample.pairs.gz | grep -v '^#'| awk '{OFS="\t"}{print $2,$3,$3+1,$4,$5,$5+1}' > tmp/contacts.$rep.sub_$sample.bedpe;
	done
	samples=("" .sub_0.1 .sub_0.25 .sub_0.5)
	for sample in "${samples[@]}";
	do			
		# overlapping contacts with binned peaks
		bedtools pairtobed -b MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed -a tmp/contacts.$rep$sample.bedpe -bedpe > tmp/contacts.peakoverlap.$rep$sample.bedpe
		#flipping contacts with baits last and removing coordinates of bait read
		awk '{OFS="\t"} {if ($2<=$8) {print $1,$2,$3,$7,$8,$9} else {print $4,$5,$6,$7,$8,$9}}' tmp/contacts.peakoverlap.$rep$sample.bedpe > tmp/contacts.peakoverlap.flipped.$rep$sample.bedpe
		# overlapping OE with bins (do not understand sort completely)
		bedtools intersect -wo -a tmp/contacts.peakoverlap.flipped.$rep$sample.bedpe -b Auxiliary_data/hg19/genomic_bins.bed | sort -k4,4 -k5,8n -k8,8n | cut -f 4,5,6,7,8,9 | uniq -c | awk {'OFS="\t"}{print $2,$3,$4,$5,$6,$7,$1}' > tmp/interactions.sum.$rep$sample.txt
		rm tmp/contacts.$rep$sample.bedpe tmp/contacts.peakoverlap.$rep$sample.bedpe tmp/contacts.peakoverlap.flipped.$rep$sample.bedpe
		# estimating coverage for viewpoints (average coverage per 250bp) and all genomic bins (250 bp each) in all 4 replicates
		cp MChIPC_output/sam/Mononucleosomal_$rep* tmp/
		samtools bedcov MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed tmp/Mononucleosomal_$rep$sample.sorted.bam > tmp/viewpoints_coverage_$rep$sample.bed
		awk '{OFS="\t"}{$5=$4*250/($3-$2)}{print}' tmp/viewpoints_coverage_$rep$sample.bed > tmp/tmp.bed && mv -f tmp/tmp.bed tmp/viewpoints_coverage_$rep$sample.bed
		samtools bedcov Auxiliary_data/hg19/genomic_bins.bed tmp/Mononucleosomal_$rep$sample.sorted.bam > tmp/genomic_bins_coverage_$rep$sample.bed
		#rm tmp/MChIPC_$rep.sub_$sample.dedup.pairsam.gz tmp/MChIPC_$rep.sub_$sample.pairs.gz tmp/Mononucleosomal_$rep.sub_$sample.sorted.bam tmp/stats_dedup_MChIPC_$rep.sub_$sample.txt
	done
done

echo "building profiles for individual viewpoints" >&3
cat MChIPC_output/MChIPC_viewpoints.bed | while read line; do  promoter=$(echo $line | cut -d " " -f 4); CHR=$(echo $line | cut -d " " -f 1); START=$(echo $line | cut -d " " -f 2); END=$(echo $line | cut -d " " -f 3); samtools view -bh -@ 16 -P MChIPC_output/sam/MChIPC_merged.bam $CHR:$START-$END |samtools sort -@ 16 > tmp/MChIPC.$promoter.bam; samtools index -@ 16 tmp/MChIPC.$promoter.bam; bamCoverage -p 16 -b tmp/MChIPC.$promoter.bam -bs 250 -o MChIPC_output/individual_viewpoint_profiles/MChIPC.$promoter.250bp.bw; bamCoverage -p 16 -b tmp/MChIPC.$promoter.bam -bs 100 -o MChIPC_output/individual_viewpoint_profiles/MChIPC.$promoter.100bp.bw; rm tmp/MChIPC.$promoter.*;  done
