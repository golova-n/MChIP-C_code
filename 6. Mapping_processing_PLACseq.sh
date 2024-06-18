
exec 2>/dev/null
: '
echo "downloading PLAC-seq fastq files"
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/94970dd9-f895-454c-86b9-1cb8190c2fad/4DNFI2I5A8ZU.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/fbcf7690-b6d2-4dfa-93c0-7bf4aef6960c/4DNFI9TNKOBD.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/cb7319ee-22ce-4bfd-b1f1-6b30db4eb800/4DNFIM5BWS5H.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/62ba96e5-0bdd-4c3c-81fe-c6e45942253a/4DNFIMWVLMTQ.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/d49823ca-91c1-4c20-9ef2-7454026e12bd/4DNFITMUHEBY.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/04739d3d-ff7b-4ef7-a6c9-d1008d37b66f/4DNFIZEP644W.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/684443f1-35ec-4170-a49f-1959293730ab/4DNFI69HRUMH.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/7c9a6106-0934-45d3-8961-2fe4710be3f0/4DNFICGRXFU6.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/755771dc-85bc-4723-90a5-a05043abf737/4DNFIJYF2QVF.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/ec3bff0b-3b5a-4ed6-a3d2-ca2ad72752c3/4DNFIMTJX2X7.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/6117857e-d531-4f83-a7c9-9374e9eafc1f/4DNFIP2TS3O8.fastq.gz
wget -q -P tmp/ https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/8efdf07a-e192-4989-a244-877a164c565f/4DNFIU98V4FB.fastq.gz
'
echo "mapping and parsing PLAC-seq rep1 data"
bwa mem -SP5M -t 8 Auxiliary_data/hg19/male.hg19.fa tmp/4DNFIM5BWS5H.fastq.gz tmp/4DNFITMUHEBY.fastq.gz | pigz -p 16 > tmp/PLACseq_rep1_1.sam.gz
pairtools parse --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -c Auxiliary_data/hg19/hg19.chrom.sizes --output-stats PLACseq_output/mapping_stats/stats_map_PLACseq_rep1_1.txt -o tmp/PLACseq_rep1_1.pairsam.gz tmp/PLACseq_rep1_1.sam.gz
rm tmp/*.sam.gz
pairtools sort --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --nproc 16 -o tmp/PLACseq_rep1_1.sorted.pairsam.gz tmp/PLACseq_rep1_1.pairsam.gz
rm tmp/PLACseq_rep1_1.pairsam.gz
bwa mem -SP5M -t 8 Auxiliary_data/hg19/male.hg19.fa tmp/4DNFI9TNKOBD.fastq.gz tmp/4DNFI2I5A8ZU.fastq.gz | pigz -p 16 > tmp/PLACseq_rep1_2.sam.gz
pairtools parse --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -c Auxiliary_data/hg19/hg19.chrom.sizes --output-stats PLACseq_output/mapping_stats/stats_map_PLACseq_rep1_2.txt -o tmp/PLACseq_rep1_2.pairsam.gz tmp/PLACseq_rep1_2.sam.gz
rm tmp/*.sam.gz
pairtools sort --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --nproc 16 -o tmp/PLACseq_rep1_2.sorted.pairsam.gz tmp/PLACseq_rep1_2.pairsam.gz
rm tmp/PLACseq_rep1_2.pairsam.gz
bwa mem -SP5M -t 8 Auxiliary_data/hg19/male.hg19.fa tmp/4DNFIZEP644W.fastq.gz tmp/4DNFIMWVLMTQ.fastq.gz | pigz -p 16 > tmp/PLACseq_rep1_3.sam.gz
pairtools parse --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -c Auxiliary_data/hg19/hg19.chrom.sizes --output-stats PLACseq_output/mapping_stats/stats_map_PLACseq_rep1_3.txt -o tmp/PLACseq_rep1_3.pairsam.gz tmp/PLACseq_rep1_3.sam.gz
rm tmp/*.sam.gz
pairtools sort --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --nproc 16 -o tmp/PLACseq_rep1_3.sorted.pairsam.gz tmp/PLACseq_rep1_3.pairsam.gz
rm tmp/PLACseq_rep1_3.pairsam.gz

echo "mapping and parsing PLAC-seq rep2 data"
bwa mem -SP5M -t 8 Auxiliary_data/hg19/male.hg19.fa tmp/4DNFIP2TS3O8.fastq.gz tmp/4DNFIU98V4FB.fastq.gz | pigz -p 16 > tmp/PLACseq_rep2_1.sam.gz
pairtools parse --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -c Auxiliary_data/hg19/hg19.chrom.sizes --output-stats PLACseq_output/mapping_stats/stats_map_PLACseq_rep2_1.txt -o tmp/PLACseq_rep2_1.pairsam.gz tmp/PLACseq_rep2_1.sam.gz
rm tmp/*.sam.gz
pairtools sort --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --nproc 16 -o tmp/PLACseq_rep2_1.sorted.pairsam.gz tmp/PLACseq_rep2_1.pairsam.gz
rm tmp/PLACseq_rep2_1.pairsam.gz
bwa mem -SP5M -t 8 Auxiliary_data/hg19/male.hg19.fa tmp/4DNFI69HRUMH.fastq.gz tmp/4DNFIMTJX2X7.fastq.gz | pigz -p 16 > tmp/PLACseq_rep2_2.sam.gz
pairtools parse --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -c Auxiliary_data/hg19/hg19.chrom.sizes --output-stats PLACseq_output/mapping_stats/stats_map_PLACseq_rep2_2.txt -o tmp/PLACseq_rep2_2.pairsam.gz tmp/PLACseq_rep2_2.sam.gz
rm tmp/*.sam.gz
pairtools sort --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --nproc 16 -o tmp/PLACseq_rep2_2.sorted.pairsam.gz tmp/PLACseq_rep2_2.pairsam.gz
rm tmp/PLACseq_rep2_2.pairsam.gz
bwa mem -SP5M -t 8 Auxiliary_data/hg19/male.hg19.fa tmp/4DNFICGRXFU6.fastq.gz tmp/4DNFIJYF2QVF.fastq.gz | pigz -p 16 > tmp/PLACseq_rep2_3.sam.gz
pairtools parse --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -c Auxiliary_data/hg19/hg19.chrom.sizes --output-stats PLACseq_output/mapping_stats/stats_map_PLACseq_rep2_3.txt -o tmp/PLACseq_rep2_3.pairsam.gz tmp/PLACseq_rep2_3.sam.gz
rm tmp/*.sam.gz
pairtools sort --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --nproc 16 -o tmp/PLACseq_rep2_3.sorted.pairsam.gz tmp/PLACseq_rep2_3.pairsam.gz
rm tmp/PLACseq_rep2_3.pairsam.gz

for rep in rep1 rep2;
do
	echo "merging and deduplicating and selecting viewpoint signal in $rep PLAC-seq"
	pairtools merge --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/PLACseq_$rep.sorted.pairsam.gz tmp/PLACseq_$rep*.sorted.pairsam.gz
	pairtools dedup --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-stats PLACseq_output/mapping_stats/stats_dedup_PLACseq_$rep.txt -o tmp/PLACseq_$rep.dedup.pairsam.gz tmp/PLACseq_$rep.sorted.pairsam.gz
	# rm
	pairtools select --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o tmp/PLACseq_$rep.cis.pairsam.gz tmp/PLACseq_$rep.dedup.pairsam.gz
	pairtools split --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-sam tmp/PLACseq_$rep.cis.sam.gz --output-pairs tmp/PLACseq_$rep.pairs.gz tmp/PLACseq_$rep.cis.pairsam.gz
	samtools view -@ 16 -bh tmp/PLACseq_$rep.cis.sam.gz > tmp/PLACseq_$rep.cis.bam
	samtools sort -@ 16 tmp/PLACseq_$rep.cis.bam > tmp/PLACseq_$rep.cis.sorted.bam
	samtools index -@ 16 tmp/PLACseq_$rep.cis.sorted.bam
	samtools view -b -P -h -@ 16 --region-file MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed tmp/PLACseq_$rep.cis.sorted.bam > PLACseq_output/sam/PLACseq_$rep.cis.filtered.bam
done

echo "merging viewpoint signals of PLAC-seq replicates"
pairtools merge --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o PLACseq_output/pairs/PLACseq.pairs.gz tmp/PLACseq_rep*.pairs.gz
samtools merge -@ 16 -o PLACseq_output/sam/PLACseq_merged.bam PLACseq_output/sam/PLACseq_rep1.cis.filtered.bam PLACseq_output/sam/PLACseq_rep2.cis.filtered.bam
samtools index -@ 16 PLACseq_output/sam/PLACseq_merged.bam
echo "downloading hg19 DpnII digest-file"
#
echo "generating aggregate PLAC-seq profiles"
bedtools coverage -b PLACseq_output/sam/PLACseq_merged.bam -a Auxiliary_data/hg19/hg19.DpnII_digest.bed -sorted -counts -g Auxiliary_data/hg19/hg19.chrom.sizes > tmp/PLACseq_merged.bedgraph
sort -k1,1 -k2,2n tmp/PLACseq_merged.bedgraph > tmp/PLACseq_merged.sorted.bedgraph
bedGraphToBigWig tmp/PLACseq_merged.sorted.bedgraph Auxiliary_data/hg19/hg19.chrom.sizes tmp/PLACseq_merged.bw
bigwigCompare -b1 tmp/PLACseq_merged.bw -b2 tmp/PLACseq_merged.bw -bs 250 --operation mean --numberOfProcessors 16 -o PLACseq_output/aggregate_profiles/PLACseq_merged.250bp.bw
rm tmp/PLACseq*
echo "loading PLAC-seq interaction files and merging them"
wget -q -P tmp/ https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161873/suppl/GSE161873%5Fformat%5FK562wtr1%2E10k%2E2%2Epeaks%2Ebedpe%2Etxt%2Egz
wget -q -P tmp/ https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161873/suppl/GSE161873%5Fformat%5FK562wtr2%2E10k%2E2%2Epeaks%2Ebedpe%2Etxt%2Egz
./Helpers/processing_PLACseq_helper_1.R
liftOver tmp/PLACseq_anchors_hg38.bed Auxiliary_data/hg19/hg38ToHg19.over.chain.gz tmp/PLACseq_anchors_hg19.bed /dev/null -minMatch=0.1
liftOver tmp/PLACseq_OE_hg38.bed Auxiliary_data/hg19/hg38ToHg19.over.chain.gz tmp/PLACseq_OE_hg19.bed /dev/null -minMatch=0.1
./Helpers/processing_PLACseq_helper_2.R
echo "generating PLAC-seq signal file for plotting E(DHS)-P heatmaps"
pigz -dc -p 16 PLACseq_output/pairs/PLACseq.pairs.gz | grep -v '^#'| awk '{OFS="\t"}{print $2,$3,$3+1,$4,$5,$5+1}' | bedtools pairtobed -b MChIPC_output/mononucleosomal_and_ChIP/binned_peaks.bed -a stdin -bedpe | awk '{OFS="\t"} {if ($2<=$8) {print $1,$2,$3,$7,$8,$9} else {print $4,$5,$6,$7,$8,$9}}'|bedtools intersect -wo -a stdin -b Auxiliary_data/hg19/hg19.DpnII_digest.bed | sort -k4,4 -k5,5n -k8,8n | cut -f 4,5,6,7,8,9 | uniq -c | awk '{OFS="\t"}{print $5,$6,$7,$2,$3,$4,$1}' | bedtools intersect -wo -a stdin -b Auxiliary_data/hg19/genomic_bins.bed | awk '{OFS="\t"}{if (sqrt(($9-$5)^2)<1010000) {print $4,$5,$6,$8,$9,$10,$7,$11}}' > tmp/PLACseq.contacts.binned.txt
rm -f PLACseq_output/PLACseq.signal.df.txt
for i in {1..22} X Y M;
do 
awk -v i=chr$i '{OFS="\t"}{if ($1==i) arr[$1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6]+=$7*$8/250}END{for (key in arr) print key,arr[key]}' tmp/PLACseq.contacts.binned.txt | sort -k2,2n -k5,5n >> PLACseq_output/PLACseq.signal.df.txt;
done
echo "generating PLAC-seq profiles for individual viewpoints"
cat Auxiliary_data/genes_of_interest.txt | while read line; do  promoter=$(echo $line | cut -d " " -f 4); CHR=$(echo $line | cut -d " " -f 1); START=$(echo $line | cut -d " " -f 2); END=$(echo $line | cut -d " " -f 3); samtools view -@ 16 -bh -P PLACseq_output/sam/PLACseq_merged.bam $CHR:$START-$END | samtools sort > tmp/PLACseq.$promoter.bam; samtools index -@ 16 tmp/PLACseq.$promoter.bam; bedtools coverage -b tmp/PLACseq.$promoter.bam -a Auxiliary_data/hg19/hg19.DpnII_digest.bed -sorted -counts -g Auxiliary_data/hg19/hg19.chrom.sizes | awk '$4>0'> PLACseq_output/individual_viewpoint_profile/PLACseq.$promoter.bedgraph; done
rm tmp/PLACseq*
