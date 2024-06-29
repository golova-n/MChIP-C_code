echo "generating mononucleosomal profiles"
for rep in rep1 rep2 rep3 rep4;
do
	samtools index -@ 8 MChIPC_output/sam/Mononucleosomal_$rep.sorted.bam
	bamCoverage -p 8 -b MChIPC_output/sam/Mononucleosomal_$rep.sorted.bam -bs 50 -e -p 8 -o MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_$rep.bw
done

echo "calling peaks in mononuclesomal profiles and identifying MChIP-C viewpoints"
mkdir -p tmp/macs/
for rep in rep1 rep2 rep3 rep4;
do
	macs2 callpeak -t MChIPC_output/sam/Mononucleosomal_$rep.sorted.bam --outdir tmp/macs/ -n Mononucleosomal_$rep -f BAMPE -q 0.0001 --max-gap 1000 --verbose 0
done
mv tmp/macs/*.narrowPeak MChIPC_output/mononucleosomal_and_ChIP
bedtools makewindows -g tmp/hg19/hg19.chrom.sizes -w 250 > tmp/hg19/genomic_bins.bed
# finding concensus peaks
bedtools multiinter -i MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep1_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep2_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep3_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep4_peaks.narrowPeak | awk '$4 >= 3 {print $1"\t"$2"\t"$3}' | bedtools merge -d 1000 -i stdin | bedtools slop -b 750 -i stdin -g tmp/hg19/hg19.chrom.sizes > tmp/concensus_peaks.bed
# creating viewpoint file (binned_peaks)
bedtools intersect -wa -a tmp/hg19/genomic_bins.bed -b tmp/concensus_peaks.bed | bedtools merge > MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed
# creating H3K4me3 mask
bedtools multiinter -i MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep1_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep2_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep3_peaks.narrowPeak MChIPC_output/mononucleosomal_and_ChIP/Mononucleosomal_rep4_peaks.narrowPeak | awk '{print $1"\t"$2"\t"$3}' | bedtools merge -i stdin > MChIPC_output/mononucleosomal_and_ChIP/H3K4me3.mask.bed

echo "generating raw and merged aggregate MChIP-C profiles"
for rep in rep1 rep2 rep3 rep4;
do
	pairtools split --cmd-in 'pigz -d -p 8' --cmd-out 'pigz -p 8' --output-sam tmp/MChIPC_$rep.sam.gz --output-pairs MChIPC_output/pairs/MChIPC_$rep.pairs.gz MChIPC_output/sam/MChIPC_$rep.dedup.pairsam.gz
	samtools view -@ 8 -bh tmp/MChIPC_$rep.sam.gz > tmp/MChIPC_$rep.bam
	samtools sort -@ 8 tmp/MChIPC_$rep.bam > MChIPC_output/sam/MChIPC_$rep.raw.bam
	samtools index -@ 8 MChIPC_output/sam/MChIPC_$rep.raw.bam
	samtools view -b -P -h -@ 8 --region-file MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed MChIPC_output/sam/MChIPC_$rep.raw.bam > MChIPC_output/sam/MChIPC_$rep.merged.bam
	samtools index -@ 8 MChIPC_output/sam/MChIPC_$rep.merged.bam
	bamCoverage -p 8 -b MChIPC_output/sam/MChIPC_$rep.merged.bam -bs 250 -o MChIPC_output/aggregate_profiles/MChIPC_$rep.merged.250bp.bw
done
samtools merge -@ 8 -o MChIPC_output/sam/MChIPC_raw.bam MChIPC_output/sam/MChIPC_rep1.raw.bam MChIPC_output/sam/MChIPC_rep2.raw.bam MChIPC_output/sam/MChIPC_rep3.raw.bam MChIPC_output/sam/MChIPC_rep4.raw.bam
samtools index -@ 8 MChIPC_output/sam/MChIPC_raw.bam
bamCoverage -p 8 -b MChIPC_raw.bam -bs 250 -o MChIPC_output/aggregate_profiles/MChIPC_raw.250bp.bw
bamCoverage -p 8 -b MChIPC_raw.bam -bs 100 -o MChIPC_output/aggregate_profiles/MChIPC_raw.100bp.bw
samtools merge -@ 8 -o MChIPC_output/sam/MChIPC_merged.bam MChIPC_output/sam/MChIPC_rep1.merged.bam MChIPC_output/sam/MChIPC_rep2.merged.bam MChIPC_output/sam/MChIPC_rep3.merged.bam MChIPC_output/sam/MChIPC_rep4.merged.bam
samtools index -@ 8 MChIPC_output/sam/MChIPC_merged.bam
bamCoverage -p 8 -b MChIPC_output/sam/MChIPC_merged.bam -bs 250 -o MChIPC_output/aggregate_profiles/MChIPC_merged.250bp.bw
bamCoverage -p 8 -b MChIPC_output/sam/MChIPC_merged.bam -bs 100 -o MChIPC_output/aggregate_profiles/MChIPC_merged.100bp.bw


echo "creating cooler MChIP-C files"
for rep in rep1 rep2 rep3 rep4;
do
	cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 tmp/hg19/hg19.chrom.sizes:250 MChIPC_output/pairs/MChIPC_$rep.pairs.gz tmp/MChIPC_$rep.250bp.cool
	cooler zoomify --nproc 8 --out MChIPC_output/cools/MChIPC_$rep.mcool --resolutions 250,500,1000,2000,5000,10000,50000,100000,500000,1000000,5000000,10000000 tmp/MChIPC_$rep.250bp.cool
done
cooler merge tmp/MChIPC.250bp.cool tmp/MChIPC_rep1.250bp.cool tmp/MChIPC_rep2.250bp.cool tmp/MChIPC_rep3.250bp.cool tmp/MChIPC_rep4.250bp.cool
cooler zoomify --nproc 8 --out MChIPC_output/cools/MChIPC.mcool --resolutions 250,500,1000,2000,5000,10000,50000,100000,500000,1000000,5000000,10000000 tmp/MChIPC.250bp.cool

echo "preparing MChIP-C input files for loop calling"
for rep in rep1 rep2 rep3 rep4;
do
	pigz -p 8 -dc MChIPC_output/pairs/MChIPC_$rep.pairs.gz | grep -v '^#'| awk '{OFS="\t"}{print $2,$3,$3+1,$4,$5,$5+1}' > tmp/contacts.$rep.bedpe
	# overlapping contacts with binned peaks
	bedtools pairtobed -b MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed -a tmp/contacts.$rep.bedpe -bedpe > tmp/contacts.peakoverlap.$rep.bedpe
	#flipping contacts with baits last and removing coordinates of bait read
	awk '{OFS="\t"} {if ($2<=$8) {print $1,$2,$3,$7,$8,$9} else {print $4,$5,$6,$7,$8,$9}}' tmp/contacts.peakoverlap.$rep.bedpe > tmp/contacts.peakoverlap.flipped.$rep.bedpe
	# overlapping OE with bins (do not understand sort completely)
	bedtools intersect -wo -a tmp/contacts.peakoverlap.flipped.$rep.bedpe -b tmp/hg19/genomic_bins.bed | sort -k4,4 -k5,8n -k8,8n | cut -f 4,5,6,7,8,9 | uniq -c | awk {'OFS="\t"}{print $2,$3,$4,$5,$6,$7,$1}' > tmp/interactions.sum.$rep.txt
	rm tmp/contacts.$rep.bedpe tmp/contacts.peakoverlap.$rep.bedpe tmp/contacts.peakoverlap.flipped.$rep.bedpe
	# estimating coverage for viewpoints (average coverage per 250bp) and all genomic bins (250 bp each) in all 4 replicates
	samtools bedcov MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed MChIPC_output/sam/Mononucleosomal_$rep.sorted.bam > tmp/viewpoints_coverage_$rep.bed
	awk '{OFS="\t"}{$5=$4*250/($3-$2)}{print}' tmp/viewpoints_coverage_$rep.bed > tmp/tmp.bed && mv -f tmp/tmp.bed tmp/viewpoints_coverage_$rep.bed
	samtools bedcov tmp/hg19/genomic_bins.bed MChIPC_output/sam/Mononucleosomal_$rep.sorted.bam > tmp/genomic_bins_coverage_$rep.bed
done
