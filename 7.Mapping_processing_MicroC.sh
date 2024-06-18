repA=(SRR20340567 SRR20340579)
repB=(SRR20340569 SRR20340581)
repC=(SRR20340571 SRR20340583)
repD=(SRR20340573 SRR20340585)
repE=(SRR20340575 SRR20340586)
repF=(SRR20340577 SRR20340588)
reps=(repA repB repC repD repE repF)
declare -n rep
for rep in ${reps[@]}
do
	echo "downloading Micro-C ${!rep} fastq files"
	fasterq-dump -e 16 -O tmp/ ${rep[0]}
	fasterq-dump -e 16 -O tmp/ ${rep[1]}
	cat tmp/${rep[0]}_1.fastq tmp/${rep[1]}_1.fastq | pigz -p 16 > tmp/MicroC_${!rep}_1.fastq.gz
	cat tmp/${rep[0]}_2.fastq tmp/${rep[1]}_2.fastq | pigz -p 16 > tmp/MicroC_${!rep}_2.fastq.gz
	echo "mapping and parsing Micro-C ${!rep}"
	bwa mem -SP5M -t 16 tmp/hg19/male.hg19.fa tmp/MicroC_${!rep}_1.fastq.gz tmp/MicroC_${!rep}_2.fastq.gz 2>/dev/null | pigz -p 16 > tmp/MicroC_${!rep}.sam.gz
	rm tmp/*.fastq.gz
	pairtools parse --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -c tmp/hg19/hg19.chrom.sizes --output-stats MicroC_output/mapping_stats/stats_map_MicroC_${!rep}.txt -o tmp/MicroC_${!rep}.pairsam.gz tmp/MicroC_${!rep}.sam.gz
	rm tmp/*.sam.gz
	pairtools sort --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --nproc 16 -o tmp/MicroC_${!rep}.sorted.pairsam.gz tmp/MicroC_${!rep}.pairsam.gz
	pairtools select --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16'  '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o tmp/MicroC_${!rep}.cis.pairsam.gz tmp/MicroC_${!rep}.sorted.pairsam.gz
	rm tmp/*.sorted.pairsam.gz
	echo "deduplicating and selecting viepoint reads in Micro-C ${!rep}"
	pairtools dedup --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-stats MicroC_output/mapping_stats/stats_dedup_MicroC_${!rep}.txt -o tmp/MicroC_${!rep}.cis.dedup.pairsam.gz tmp/MicroC_${!rep}.cis.pairsam.gz
	rm tmp/*.cis.pairsam.gz
	pairtools split --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-sam tmp/MicroC_${!rep}.cis.sam.gz tmp/MicroC_${!rep}.cis.dedup.pairsam.gz
	rm tmp/*.dedup.pairsam.gz
	samtools sort -@ 16 tmp/MicroC_${!rep}.cis.sam.gz > tmp/MicroC_${!rep}.cis.sorted.bam
	rm tmp/*.sam.gz
	samtools index -@ 16 tmp/MicroC_${!rep}.cis.sorted.bam
	samtools view -b -P -h -@ 16 --region-file MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed tmp/MicroC_${!rep}.cis.sorted.bam > MicroC_output/sam/MicroC_${!rep}.cis.filtered.bam
	rm tmp/*.bam
done
echo "merging Micro-C replicates"
samtools merge -@ 16 -o MicroC_output/sam/MicroC_merged.bam MicroC_output/sam/MicroC_*.cis.filtered.bam

exec 2>/dev/null

echo "generating aggregate Micro-C profiles"
samtools index -@ 16 MicroC_output/sam/MicroC_merged.bam
bamCoverage -p 16 -b MicroC_output/sam/MicroC_merged.bam -bs 250 -o MicroC_output/aggregate_profiles/MicroC_merged.250bp.bw
echo "loading mcool Micro-C file"
echo "extracting promoter-relevant data from the genome-wide matrix"
cooler dump --join -b MicroC_output/cools/MicroC.mcool::/resolutions/250 | awk '($1 == $4)' | awk '$8=="" {$8="0"} 1' OFS="\t" | pairToBed -a stdin -b MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed > tmp/MicroC_filtered.bedpe
./Helpers/processing_MicroC_helper_1.R
echo "generating aggregate ICE Micro-C profiles"
bedGraphToBigWig tmp/MicroC.ICE.bedgraph tmp/hg19/hg19.chrom.sizes tmp/MicroC.ICE.bw
bigwigCompare -p 16 -b1 tmp/MicroC.ICE.bw -b2 tmp/MicroC.ICE.bw -bs 250 --operation mean -o MicroC_output/aggregate_profiles/MicroC_merged.ICE.250bp.bw
echo "calling Micro-C point interactions"
:"
for res in 250 500 1kb 2kb 5kb 10kb
do
	mustache -f MicroC_output/cools/MicroC.mcool -r $res -pt 0.1 -st 0.88 -oc 2 -p 16 -v 0 -o tmp/MicroC_loops_$res.tsv
done
./Helpers/processing_MicroC_helper_2.R
"
echo "generating Micro-C signal file for plotting E(DHS)-P heatmaps"
awk '{OFS="\t"} {if ($10<=$3) {print $4,$5,$6,$9,$10,$11,$7,$8} else {print $1,$2,$3,$9,$10,$11,$7,$8}}' tmp/MicroC_filtered.bedpe | bedtools intersect -wo -a stdin -b tmp/hg19/genomic_bins.bed -f 0.05 | sort -k4,4 -k5,5n -k10,10n |awk {'OFS="\t"}{if (sqrt(($5-$10)^2) < 1010000) print $4,$5,$6,$9,$10,$11,$7,$8}' > tmp/MicroC.interactions.txt
for i in {1..22} X Y M;
do
LC_NUMERIC=en_US.UTF-8 awk -v i=chr$i '{OFS="\t"}{if ($1==i) key=$1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6}{raw[key]+=$7; ICE[key]+=$8}END{ for(key in raw) print key,raw[key],ICE[key]}' tmp/MicroC.interactions.txt | sort -k2,2n -k5,5n >> MicroC_output/MicroC.signal.df.txt;
done
echo "generating Micro-C profiles for individual viewpoints"
cat Auxiliary_data/genes_of_interest.txt | while read line; 
do
promoter=$(echo $line | cut -d " " -f 4); CHR=$(echo $line | cut -d " " -f 1); START=$(echo $line | cut -d " " -f 2); END=$(echo $line | cut -d " " -f 3); samtools view -@ 16 -bh -P MicroC_output/sam/MicroC_merged.bam $CHR:$START-$END | samtools sort -@ 16 > tmp/MicroC.$promoter.bam; samtools index -@ 16 tmp/MicroC.$promoter.bam; bamCoverage -p 16 -b tmp/MicroC.$promoter.bam -bs 250 -o MicroC_output/individual_viewpoint_profiles/MicroC.$promoter.250bp.bw; bamCoverage -p 16 -b tmp/MicroC.$promoter.bam -bs 100 -o MicroC_output/individual_viewpoint_profiles/MicroC.$promoter.100bp.bw; done
rm tmp/MicroC*
