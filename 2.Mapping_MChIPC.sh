exec 3>&1 &>/dev/null

mkdir -p MChIPC_output/mapping_stats
mkdir -p MChIPC_output/pairs
mkdir -p MChIPC_output/cools
mkdir -p MChIPC_output/aggregate_profiles
mkdir -p MChIPC_output/mononucleosomal_and_ChIP
mkdir -p MChIPC_output/individual_viewpoint_profiles
mkdir -p MChIPC_output/sam

echo "downloading MChIPC_rep1_1 fastq files" >&3
wget -q -P tmp/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23410275/SRR23410275
fasterq-dump -e 16 -O tmp/ tmp/SRR23410275
rm tmp/SRR23410275
pigz -p 16 tmp/*.fastq

echo "mapping MChIPC_rep1_1" >&3
bwa mem -SP5M -t 16 Auxiliary_data/hg19/male.hg19.fa tmp/SRR23410275_1.fastq.gz tmp/SRR23410275_2.fastq.gz 2>/dev/null | pigz -p 16 > tmp/MChIPC_rep1_1.sam.gz
rm tmp/*.fastq.gz

echo "downloading MChIPC_rep1_2 fastq files" >&3
wget -q -P tmp/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23410276/SRR23410276
fasterq-dump -e 16 -O tmp/ tmp/SRR23410276 
rm tmp/SRR23410276
pigz -p 16 tmp/*.fastq

echo "mapping MChIPC_rep1_2" >&3
bwa mem -SP5M -t 16 Auxiliary_data/hg19/male.hg19.fa tmp/SRR23410276_1.fastq.gz tmp/SRR23410276_2.fastq.gz 2>/dev/null | pigz -p 16 > tmp/MChIPC_rep1_2.sam.gz
rm tmp/*.fastq.gz

echo "downloading MChIPC_rep1_3 fastq files" >&3
wget -q -P tmp/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23410274/SRR23410274
fasterq-dump -e 16 -O tmp/ tmp/SRR23410274 
rm tmp/SRR23410274
pigz -p 16 tmp/*.fastq

echo "mapping MChIPC_rep1_3" >&3
bwa mem -SP5M -t 16 Auxiliary_data/hg19/male.hg19.fa tmp/SRR23410274_1.fastq.gz tmp/SRR23410274_2.fastq.gz 2>/dev/null | pigz -p 16 > tmp/MChIPC_rep1_3.sam.gz
rm tmp/*.fastq.gz

echo "downloading MChIPC_rep2 fastq files" >&3
wget -q -P tmp/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23410273/SRR23410273
fasterq-dump -e 16 -O tmp/ tmp/SRR23410273 
rm tmp/SRR23410273
pigz -p 16 tmp/*.fastq

echo "mapping MChIPC_rep2" >&3
bwa mem -SP5M -t 16 Auxiliary_data/hg19/male.hg19.fa tmp/SRR23410273_1.fastq.gz tmp/SRR23410273_2.fastq.gz 2>/dev/null | pigz -p 16 > tmp/MChIPC_rep2.sam.gz
rm tmp/*.fastq.gz

echo "downloading MChIPC_rep3 fastq files"
wget -q -P tmp/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23410272/SRR23410272
fasterq-dump -e 16 -O tmp/ tmp/SRR23410272
rm tmp/SRR23410272
pigz -p 16 tmp/*.fastq

echo "mapping MChIPC_rep3"
bwa mem -SP5M -t 16 Auxiliary_data/hg19/male.hg19.fa tmp/SRR23410272_1.fastq.gz tmp/SRR23410272_2.fastq.gz 2>/dev/null | pigz -p 16 > tmp/MChIPC_rep3.sam.gz
rm tmp/*.fastq.gz

echo "downloading MChIPC_rep4 fastq files" >&3
wget -q -P tmp/ https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR23410271/SRR23410271
fasterq-dump -e 16 -O tmp/ tmp/SRR23410271
rm tmp/SRR23410271
pigz -p 16 tmp/*.fastq

echo "mapping MChIPC_rep4" >&3
bwa mem -SP5M -t 16 Auxiliary_data/hg19/male.hg19.fa tmp/SRR23410271_1.fastq.gz tmp/SRR23410271_2.fastq.gz 2>/dev/null | pigz -p 16 > tmp/MChIPC_rep4.sam.gz
rm tmp/*.fastq.gz

for rep in rep1_1 rep1_2 rep1_3 rep2 rep3 rep4;
do
	echo "parsing MChIPC_$rep data" >&3
	pairtools parse --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -c Auxiliary_data/hg19/hg19.chrom.sizes --output-stats MChIPC_output/mapping_stats/stats_map_MChIPC_$rep.txt -o tmp/MChIPC_$rep.pairsam.gz tmp/MChIPC_$rep.sam.gz;
	rm tmp/MChIPC_$rep.sam.gz;
	pairtools sort --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --nproc 16 -o tmp/MChIPC_$rep.sorted.pairsam.gz tmp/MChIPC_$rep.pairsam.gz;
	pairtools select --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16'  '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o tmp/MChIPC_$rep.cis.pairsam.gz tmp/MChIPC_$rep.sorted.pairsam.gz;
	pairtools select --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 200) and (strand1=="+") and (strand2=="-"))' -o tmp/MN_$rep.pairsam.gz tmp/MChIPC_$rep.sorted.pairsam.gz;
	echo "sampling from MChIPC_$rep data" >&3
	for sample in 0.1 0.25 0.5;
	do 
		echo -e "\tsampling $sample proportion of reads" >&3
		pairtools sample --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MChIPC_$rep.sub_$sample.pairsam.gz $sample tmp/MChIPC_$rep.pairsam.gz;
		pairtools sort --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --nproc 16 -o tmp/MChIPC_$rep.sub_$sample.sorted.pairsam.gz tmp/MChIPC_$rep.sub_$sample.pairsam.gz;
		pairtools select --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16'  '((chrom1==chrom2) and (abs(pos1 - pos2) > 5000))' -o tmp/MChIPC_$rep.sub_$sample.cis.pairsam.gz tmp/MChIPC_$rep.sub_$sample.sorted.pairsam.gz;
		pairtools select --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' '((chrom1==chrom2) and (abs(pos1 - pos2) > 100) and (abs(pos1-pos2) < 200) and (strand1=="+") and (strand2=="-"))' -o tmp/MN_$rep.sub_$sample.pairsam.gz tmp/MChIPC_$rep.sub_$sample.sorted.pairsam.gz;
		rm tmp/MChIPC_$rep.sub_$sample.pairsam.gz
	done
	rm tmp/MChIPC_$rep.pairsam.gz;
	rm tmp/*.sorted.pairsam.gz;
done

echo "merging technical replicates for MChIPC_rep1 mononucleosomal signal" >&3
pairtools merge --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MN_rep1.pairsam.gz tmp/MN_rep1_1.pairsam.gz tmp/MN_rep1_2.pairsam.gz tmp/MN_rep1_3.pairsam.gz
pairtools merge --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MN_rep1.sub_0.1.pairsam.gz tmp/MN_rep1_1.sub_0.1.pairsam.gz tmp/MN_rep1_2.sub_0.1.pairsam.gz tmp/MN_rep1_3.sub_0.1.pairsam.gz
pairtools merge --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MN_rep1.sub_0.25.pairsam.gz tmp/MN_rep1_1.sub_0.25.pairsam.gz tmp/MN_rep1_2.sub_0.25.pairsam.gz tmp/MN_rep1_3.sub_0.25.pairsam.gz
pairtools merge --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MN_rep1.sub_0.5.pairsam.gz tmp/MN_rep1_1.sub_0.5.pairsam.gz tmp/MN_rep1_2.sub_0.5.pairsam.gz tmp/MN_rep1_3.sub_0.5.pairsam.gz
echo "merging technical replicates for MChIPC_rep1 MChIP-C signal" >&3
pairtools merge --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MChIPC_rep1.cis.pairsam.gz tmp/MChIPC_rep1_1.cis.pairsam.gz tmp/MChIPC_rep1_2.cis.pairsam.gz tmp/MChIPC_rep1_3.cis.pairsam.gz
pairtools merge --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MChIPC_rep1.sub_0.1.cis.pairsam.gz tmp/MChIPC_rep1_1.cis.sub_0.1.pairsam.gz tmp/MChIPC_rep1_2.sub_0.1.cis.pairsam.gz tmp/MChIPC_rep1_3.sub_0.1.cis.pairsam.gz
pairtools merge --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MChIPC_rep1.sub_0.25.cis.pairsam.gz tmp/MChIPC_rep1_1.cis.sub_0.25.pairsam.gz tmp/MChIPC_rep1_2.sub_0.25.cis.pairsam.gz tmp/MChIPC_rep1_3.sub_0.25.cis.pairsam.gz
pairtools merge --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MChIPC_rep1.sub_0.5.cis.pairsam.gz tmp/MChIPC_rep1_1.cis.sub_0.5.pairsam.gz tmp/MChIPC_rep1_2.sub_0.5.cis.pairsam.gz tmp/MChIPC_rep1_3.sub_0.5.cis.pairsam.gz
#rm tmp/*rep1_*

for rep in rep1 rep2 rep3 rep4;
do
	echo "primary processing of $rep data" >&3
	echo -e "\tdeduplicating mononucleosomal signal from full dataset" >&3
	nreads=$(gunzip -c tmp/MN_$rep.pairsam.gz | wc -l);
	prop=$(bc <<< "scale=2; 20000000/$nreads");
	pairtools sample --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MN_$rep.subsampled.pairsam.gz $prop tmp/MN_$rep.pairsam.gz;
	pairtools dedup --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-stats MChIPC_output/mapping_stats/stats_dedup_Mononucleosomal_$rep.txt -o tmp/MN_$rep.dedup.pairsam.gz tmp/MN_$rep.subsampled.pairsam.gz;
	pairtools split --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-sam tmp/MN_$rep.sam.gz tmp/MN_$rep.dedup.pairsam.gz;
	# I have to add flag 2 (properly paired) for MACS2 to work normally
	samtools view -@ 16 -bh tmp/MN_$rep.sam.gz --add-flag 2 > tmp/MN_$rep.bam;
	samtools sort -@ 16 tmp/MN_$rep.bam > MChIPC_output/sam/Mononucleosomal_$rep.sorted.bam;
	echo -e "\tdeduplicating MChIP-C signal" >&3
	pairtools dedup --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-stats MChIPC_output/mapping_stats/stats_dedup_MChIPC_$rep.txt -o MChIPC_output/sam/MChIPC_$rep.dedup.pairsam.gz tmp/MChIPC_$rep.cis.pairsam.gz;
	rm tmp/MChIPC_$rep.cis.pairsam.gz;
done

echo "primary processing of downsampled datasets" >&3
for sample in 0.1 0.25 0.5;
do
	echo -e "\tprocessing $sample (downsampled) datasets" >&3
	for rep in rep1 rep2 rep3 rep4;
	do
		eval "nreads_$rep=$(gunzip -c tmp/MN_$rep.sub_$sample.pairsam.gz | wc -l)";
	done
	max=$(printf "%s\n" $nreads_rep1 $nreads_rep2 $nreads_rep3 $nreads_rep4 | sort -n | head -1)
	if [ $max -gt 20000000 ]; 
	then
		max=20000000; 
	fi
	for rep in rep1 rep2 rep3 rep4;
	do
		eval "nreads=\$nreads_$rep";
		prop=$(bc <<< "scale=2; $max/$nreads");
		pairtools sample --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MN_$rep.sub_$sample.subsampled.pairsam.gz $prop tmp/MN_$rep.sub_$sample.pairsam.gz;
		pairtools dedup --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' -o tmp/MN_$rep.sub_$sample.dedup.pairsam.gz tmp/MN_$rep.sub_$sample.subsampled.pairsam.gz;
		pairtools split --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-sam tmp/MN_$rep.sub_$sample.sam.gz tmp/MN_$rep.sub_$sample.dedup.pairsam.gz;
		# I have to add flag 2 (properly paired) for MACS2 to work normally
		samtools view -@ 16 -bh tmp/MN_$rep.sub_$sample.sam.gz --add-flag 2 > tmp/MN_$rep.sub_$sample.bam;
		samtools sort -@ 16 tmp/MN_$rep.sub_$sample.bam > tmp/Mononucleosomal_$rep.sub_$sample.sorted.bam;
		samtools index -@ 16 tmp/Mononucleosomal_$rep.sub_$sample.sorted.bam;
		pairtools dedup --cmd-in 'pigz -d -p 16' --cmd-out 'pigz -p 16' --output-stats tmp/stats_dedup_MChIPC_$rep.sub_$sample.txt -o tmp/MChIPC_$rep.sub_$sample.dedup.pairsam.gz tmp/MChIPC_$rep.sub_$sample.cis.pairsam.gz;
		rm tmp/MChIPC_$rep.sub_$sample.cis.pairsam.gz;
	done
done
rm tmp/MN_*
echo "done" >&3
