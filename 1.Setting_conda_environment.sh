#!/usr/bin/env bash

exec 3>&1 &>/dev/null

mkdir -p Figures/Fig.2
mkdir -p Figures/Fig.3
mkdir -p Figures/Fig.4
mkdir -p Figures/Fig.5
mkdir -p Figures/Fig.S1
mkdir -p Figures/Fig.S2
mkdir -p Figures/Fig.S3
mkdir -p Figures/Fig.S4
mkdir -p Figures/Fig.S5
mkdir -p Summary_output_datasets

echo "setting conda environment" >&3
conda install -y -n base conda-libmamba-solver
conda config --set solver libmamba
conda env create -f mchipc_environment.yml
conda activate mchip-c

echo "installing R packages" >&3
./Helpers/setting_conda_helper.R

echo "downloading auxiliary data" >&3
wget https://storage.googleapis.com/encode-pipeline-genome-data/hg19/bwa_index/male.hg19.fa.tar -q -P Auxiliary_data/hg19
tar -xf Auxiliary_data/hg19/male.hg19.fa.tar -C Auxiliary_data/hg19/
rm Auxiliary_data/hg19/male.hg19.fa.tar
wget https://storage.googleapis.com/encode-pipeline-genome-data/hg19/hg19.chrom.sizes -q -P Auxiliary_data/hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -q -P Auxiliary_data/hg19
gunzip -d Auxiliary_data/hg19/hg19.DpnII_digest.bed.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FK562%5FHiCCUPS%5Flooplist%5Fwith%5Fmotifs%2Etxt%2Egz -q -P Auxiliary_data/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861%5Fall%5Fdeg%5Fresults%2Eat%5Fscale%2Etxt%2Egz -q -P Auxiliary_data/CRISPRi_data
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmK562HMM.bed.gz -q -P Auxiliary_data/
# load digest file
echo "loading epigenomic profile and peak files" >&3

for item in `cat Auxiliary_data/bed_list.txt`
do
	wget -P Auxiliary_data/peaks $item
done


for item in `cat Auxiliary_data/bigwig_list.txt`
do
	wget -P Auxiliary_data/bigwigs $item
done

wget -O tmp/GSM2635249_ChIP-seq_K562_DMSO_BRD4.wig.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2635249&format=file&file=GSM2635249%5FChIP%2Dseq%5FK562%5FDMSO%5FBRD4%2Ewig%2Egz"
wigToBigWig tmp/GSM2635249_ChIP-seq_K562_DMSO_BRD4.wig.gz Auxiliary_data/hg19/hg19.chrom.sizes Auxiliary_data/bigwigs/BRD4.hg19.bigWig -clip

wget -O Auxiliary_data/peaks/BRD4.hg19.peaks.bed.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2635249&format=file&file=GSM2635249%5FChIP%2Dseq%5FK562%5FDMSO%5FBRD4%5Fpeaks%2Ebed%2Egz"

wget -P tmp/ https://www.encodeproject.org/files/ENCFF032HVZ/@@download/ENCFF032HVZ.bigWig
CrossMap.py bigwig Auxiliary_data/hg19/hg38ToHg19.over.chain.gz tmp/ENCFF032HVZ.bigWig Auxiliary_data/bigwigs/CHD2.hg19.bigWig
mv Auxiliary_data/bigwigs/CHD2.hg19.bigWig.bw Auxiliary_data/bigwigs/CHD2.hg19.bigWig
wget -P tmp/ https://www.encodeproject.org/files/ENCFF947AEO/@@download/ENCFF947AEO.bed.gz
liftOver tmp/ENCFF947AEO.bed.gz -bedPlus=3 Auxiliary_data/hg19/hg38ToHg19.over.chain.gz Auxiliary_data/peaks/CHD2.hg19.peaks.bed.gz tmp/unmapped

wget https://zenodo.org/records/12805899/files/56379.bw.bigwig -c -O 'tmp/56379.bw.bigwig'
CrossMap.py bigwig Auxiliary_data/hg19/hg38ToHg19.over.chain.gz tmp/56379.bw.bigwig Auxiliary_data/bigwigs/CDK8.hg19.bigWig
mv Auxiliary_data/bigwigs/CDK8.hg19.bigWig.bw Auxiliary_data/bigwigs/CDK8.hg19.bigWig

wget https://zenodo.org/records/12805899/files/56379_peaks.bed -c -O 'tmp/56379_peaks.bed'
liftOver tmp/56379_peaks.bed -bedPlus=3  Auxiliary_data/hg19/hg38ToHg19.over.chain.gz Auxiliary_data/peaks/CDK8.hg19.peaks.bed.gz tmp/unmapped

wget https://zenodo.org/records/12805899/files/74667.bw.bigwig -c -O 'tmp/74667.bw.bigwig'
CrossMap.py bigwig Auxiliary_data/hg19/hg38ToHg19.over.chain.gz tmp/74667.bw.bigwig Auxiliary_data/bigwigs/MED1.hg19.bigWig
mv Auxiliary_data/bigwigs/MED1.hg19.bigWig.bw Auxiliary_data/bigwigs/MED1.hg19.bigWig

wget https://zenodo.org/records/12805899/files/74667_peaks.bed -c -O 'tmp/74667_peaks.bed'
liftOver tmp/74667_peaks.bed -bedPlus=3  Auxiliary_data/hg19/hg38ToHg19.over.chain.gz Auxiliary_data/peaks/MED1.hg19.peaks.bed.gz tmp/unmapped

rm tmp/*
rm Auxiliary_data/bigwigs/*.bgr
echo "done" >&3
