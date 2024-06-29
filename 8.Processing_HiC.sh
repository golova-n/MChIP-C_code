exec 2>/dev/null

#mkdir -p HiC_output/cools
#mkdir -p HiC_output/aggregate_profiles

echo "downloading mcool Hi-C file"
# find a source for our mcool file
echo "extracting promoter-relevant data from the genome-wide matrix"
cooler dump --join -b HiC_output/cools/HiC.mcool::/resolutions/1000 | awk '($1 == $4)' | awk '$8=="" {$8="0"} 1' OFS="\t" | pairToBed -a stdin -b MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed > tmp/HiC_filtered.bedpe
./Helpers/processing_HiC_helper_1.R
echo "generating aggregate Hi-C profiles"
bedGraphToBigWig tmp/HiC.bedgraph Auxiliary_data/hg19/hg19.chrom.sizes tmp/HiC.bw
bigwigCompare -p 16 -b1 tmp/HiC.bw -b2 tmp/HiC.bw -bs 250 --operation mean -o HiC_output/aggregate_profiles/HiC_merged.250bp.bw
echo "generating aggregate ICE Hi-C profiles"
bedGraphToBigWig tmp/HiC.ICE.bedgraph Auxiliary_data/hg19/hg19.chrom.sizes tmp/HiC.ICE.bw
bigwigCompare -p 16 -b1 tmp/HiC.ICE.bw -b2 tmp/HiC.ICE.bw -bs 250 --operation mean -o HiC_output/aggregate_profiles/HiC_merged.ICE.250bp.bw
echo "calling Hi-C point interactions"
for res in 1kb 2kb 5kb 10kb 20kb
do
	mustache -f HiC_output/cools/HiC.mcool -r $res -pt 0.1 -st 0.88 -oc 2 -p 16 -v 0 -o tmp/HiC_loops_$res.tsv
done
./Helpers/processing_HiC_helper_2.R
echo "generating Hi-C signal file for plotting E(DHS)-P heatmaps"
awk '{OFS="\t"} {if ($10<=$3) {print $4,$5,$6,$9,$10,$11,$7,$8} else {print $1,$2,$3,$9,$10,$11,$7,$8}}' tmp/HiC_filtered.bedpe | bedtools intersect -wo -a stdin -b Auxiliary_data/hg19/genomic_bins.bed -f 0.05 | sort -k4,4 -k5,6n -k10,10n |awk {'OFS="\t"}{if (sqrt(($5-$10)^2) < 1010000) print $4,$5,$6,$9,$10,$11,$7,$8}' > tmp/HiC.interactions.txt
rm -f HiC_output/HiC.signal.df.txt;
for i in {1..22} X Y M;
do
LC_NUMERIC=en_US.UTF-8 awk -v i=chr$i '{OFS="\t"}{if ($1==i) key=$1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6}{raw[key]+=$7; ICE[key]+=$8}END{ for(key in raw) print key,raw[key],ICE[key]}' tmp/HiC.interactions.txt | sort -k2,2n -k5,5n >> HiC_output/HiC.signal.df.txt;
done
rm tmp/HiC*
