exec 3>&1 &>/dev/null

mkdir -p Auxiliary_data/hg19
mkdir -p Figures/Fig.2
mkdir -p Figures/Fig.3
mkdir -p Figures/Fig.4
mkdir -p Figures/Fig.5
mkdir -p Figures/Fig.S1
mkdir -p Figures/Fig.S2
mkdir -p Figures/Fig.S3
mkdir -p Figures/Fig.S4
mkdir -p Figures/Fig.S5

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
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FK562%5FHiCCUPS%5Flooplist%5Fwith%5Fmotifs%2Etxt%2Egz -q -P Auxiliary_data/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120861/suppl/GSE120861%5Fall%5Fdeg%5Fresults%2Eat%5Fscale%2Etxt%2Egz -q -p Auxiliary_data/CRISPRi_data
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
# change to google-disk (problems with 
wget -O Auxiliary_data/peaks/BRD4.hg19.peaks.bed.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2635249&format=file&file=GSM2635249%5FChIP%2Dseq%5FK562%5FDMSO%5FBRD4%5Fpeaks%2Ebed%2Egz"

wget -P tmp/ https://www.encodeproject.org/files/ENCFF032HVZ/@@download/ENCFF032HVZ.bigWig 
CrossMap.py bigwig Auxiliary_data/hg19/hg38ToHg19.over.chain.gz tmp/ENCFF032HVZ.bigWig Auxiliary_data/bigwigs/CHD2.hg19.bigWig
mv Auxiliary_data/bigwigs/CHD2.hg19.bigWig.bw Auxiliary_data/bigwigs/CHD2.hg19.bigWig
wget -P tmp/ https://www.encodeproject.org/files/ENCFF947AEO/@@download/ENCFF947AEO.bed.gz
liftOver tmp/ENCFF947AEO.bed.gz -bedPlus=3 Auxiliary_data/hg19/hg38ToHg19.over.chain.gz Auxiliary_data/peaks/CHD2.hg19.peaks.bed.gz tmp/unmapped

# at this point I was not able to figure out how to do automatic download from cistrome DB (I'm putting those two bigwigs CDK8 and MED1 in google-disk for now)
wget --header='Host: drive.usercontent.google.com' --header='User-Agent: Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/123.0.0.0 Safari/537.36' --header='Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7' --header='Accept-Language: en-US,en;q=0.9,fr;q=0.8,he;q=0.7,uk;q=0.6,ru;q=0.5' --header='Cookie: HSID=AreFZnLEAYNCZd8lD; SSID=AB0KHDJ_k2CmlfeN2; APISID=w0OR5kHsdu-ViFOc/ADZtMAgS6b409t16Q; SAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; __Secure-1PAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; __Secure-3PAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; SEARCH_SAMESITE=CgQIiZsB; AEC=AQTF6HyN-7_1hZyuKLji2rpC5aGjsaO5E70qWDrQ84fwMKw_w7R-JzTT6WE; __Secure-ENID=19.SE=CKt4IM96j1JmgStIOmsvbq1eHWW-yI7A_A-opmoxGtw1RYvGsocdl4BaC4yGquZ_WIiJUoP7vudMHSRG0R0oCkfDAjn1mZRlHsrPnNPjlJjl3bVynKYoi0mr1coOrv6UtdmvtH6b4sOhKGHmnM83pviMnbBAhMBam5XgFLyeIiChrt3lG-JejR9WK9Zd4_4NT3aQcUE1VIGpNfjBdSLYzSHhEa7qg-ew5zcYgaYVDSTq4guxMqjg1hDR4KUc9U4crWSJGb03tDSntMmYwcSd2W7ASoSqUluAcyUlNQ; SID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSBSnOHQwBQF192QvY95oVJhwACgYKATISAQASFQHGX2Mi_3VhqO3HX710vnTopsYTUBoVAUF8yKocFRIn3AefJYol0psthxtF0076; __Secure-1PSID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSB0Yb863nI---ZVBChR265tQACgYKAUoSAQASFQHGX2MiIVcVPKGRidM0D6GDIA4ulRoVAUF8yKoAfd6uG1D-4hxygbgQjAf30076; __Secure-3PSID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSBcx3LThkwemjJg_so7918AgACgYKAecSAQASFQHGX2Mi0yLhhbEuEoG0UwIbcuenyRoVAUF8yKrEx4mWAY8_FAK8oT-XmIUd0076; NID=514=o9tNZU5jtMG3Dj4YIejY3AsbT-KAUQL4L9Cj14eGA-RVqo08pPrCNgBfzFq-1M1OKYIhKXI9FR6Wu9DFlsy1Q5GwUzUOaIiH8FyujiHfB-8mhAbK3BYd_9hNshOAiDbTVZ1kwoAUlm08bAWD_W1XaTlPc0fP7dJErB8rZ1RBuhmoGCgmhZL57rOvWv9PD1WWkQdU3sD9ybOEVgQ_mAaAxn7Tnjk_8GUBW1IIEoCu8us0; __Secure-1PSIDTS=sidts-CjIBLwcBXJNJQWKdBrhrrlo1G7BlDA-YYC7aI_6pxpUXHQzXMSmE2FVg5EmbPPESr4toIBAA; __Secure-3PSIDTS=sidts-CjIBLwcBXJNJQWKdBrhrrlo1G7BlDA-YYC7aI_6pxpUXHQzXMSmE2FVg5EmbPPESr4toIBAA; SIDCC=AKEyXzUFCmpD9HTkgZwWLnoRxY3IdB-R7SLEhIVgOKZ1JPXTua7TQOqIsO5i52ipJ-eEy7KJn70; __Secure-1PSIDCC=AKEyXzU6Wc9nKK8TSviJXupDnw8HpJwmaplu0DNXYU3AOmdOH3uC7nvQor5oSYRzNxd_HL7m6Gs; __Secure-3PSIDCC=AKEyXzWQZ1ngC-mqzc6rBG8nUVq4FsiIvU4i1SjlI-DaX1gk65d9-j-zmDldBEx60GlznPUmfg' --header='Connection: keep-alive' 'https://drive.usercontent.google.com/download?id=10cX-ebBkF3DtDUOvOL8NSvzvm8W4in9s&export=download&authuser=0&confirm=t&uuid=2d3a06e3-6207-480a-8cb1-30a489906f7d&at=APZUnTV53BB2cHtMx8TIZmy5y_uu:1715794223070' -c -O 'tmp/56379.bw.bigwig'
CrossMap.py bigwig Auxiliary_data/hg19/hg38ToHg19.over.chain.gz tmp/56379.bw.bigwig Auxiliary_data/bigwigs/CDK8.hg19.bigWig
mv Auxiliary_data/bigwigs/CDK8.hg19.bigWig.bw Auxiliary_data/bigwigs/CDK8.hg19.bigWig
wget --header='Host: drive.usercontent.google.com' --header='User-Agent: Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/123.0.0.0 Safari/537.36' --header='Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7' --header='Accept-Language: en-US,en;q=0.9,fr;q=0.8,he;q=0.7,uk;q=0.6,ru;q=0.5' --header='Cookie: HSID=AreFZnLEAYNCZd8lD; SSID=AB0KHDJ_k2CmlfeN2; APISID=w0OR5kHsdu-ViFOc/ADZtMAgS6b409t16Q; SAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; __Secure-1PAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; __Secure-3PAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; SEARCH_SAMESITE=CgQIiZsB; AEC=AQTF6HyN-7_1hZyuKLji2rpC5aGjsaO5E70qWDrQ84fwMKw_w7R-JzTT6WE; __Secure-ENID=19.SE=CKt4IM96j1JmgStIOmsvbq1eHWW-yI7A_A-opmoxGtw1RYvGsocdl4BaC4yGquZ_WIiJUoP7vudMHSRG0R0oCkfDAjn1mZRlHsrPnNPjlJjl3bVynKYoi0mr1coOrv6UtdmvtH6b4sOhKGHmnM83pviMnbBAhMBam5XgFLyeIiChrt3lG-JejR9WK9Zd4_4NT3aQcUE1VIGpNfjBdSLYzSHhEa7qg-ew5zcYgaYVDSTq4guxMqjg1hDR4KUc9U4crWSJGb03tDSntMmYwcSd2W7ASoSqUluAcyUlNQ; SID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSBSnOHQwBQF192QvY95oVJhwACgYKATISAQASFQHGX2Mi_3VhqO3HX710vnTopsYTUBoVAUF8yKocFRIn3AefJYol0psthxtF0076; __Secure-1PSID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSB0Yb863nI---ZVBChR265tQACgYKAUoSAQASFQHGX2MiIVcVPKGRidM0D6GDIA4ulRoVAUF8yKoAfd6uG1D-4hxygbgQjAf30076; __Secure-3PSID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSBcx3LThkwemjJg_so7918AgACgYKAecSAQASFQHGX2Mi0yLhhbEuEoG0UwIbcuenyRoVAUF8yKrEx4mWAY8_FAK8oT-XmIUd0076; NID=514=o9tNZU5jtMG3Dj4YIejY3AsbT-KAUQL4L9Cj14eGA-RVqo08pPrCNgBfzFq-1M1OKYIhKXI9FR6Wu9DFlsy1Q5GwUzUOaIiH8FyujiHfB-8mhAbK3BYd_9hNshOAiDbTVZ1kwoAUlm08bAWD_W1XaTlPc0fP7dJErB8rZ1RBuhmoGCgmhZL57rOvWv9PD1WWkQdU3sD9ybOEVgQ_mAaAxn7Tnjk_8GUBW1IIEoCu8us0; __Secure-1PSIDTS=sidts-CjIBLwcBXJNJQWKdBrhrrlo1G7BlDA-YYC7aI_6pxpUXHQzXMSmE2FVg5EmbPPESr4toIBAA; __Secure-3PSIDTS=sidts-CjIBLwcBXJNJQWKdBrhrrlo1G7BlDA-YYC7aI_6pxpUXHQzXMSmE2FVg5EmbPPESr4toIBAA; SIDCC=AKEyXzUXVAWh7tyrncxwwq0i8HSan86eGf2dlSIHHdCeacv1Z5xcUtdwVUOqxA6RdmpQsmg08SI; __Secure-1PSIDCC=AKEyXzVHFWTWhCF9Cv7pfwVsuH8f2rFjX-MhrB6gKe3YkIpsnOgCZBFVmK2FrbPdkmbSsdmm9Qg; __Secure-3PSIDCC=AKEyXzU3CSFmVhx_mlGqsKvrX6n6Pxjh5VL-HtA7XowRPxHoojvMA-tUrm_dVUjU95rtYcHrXQ' --header='Connection: keep-alive' 'https://drive.usercontent.google.com/download?id=15cHoIAR6WpphBWQgSYqXG4JVu6HHmutP&export=download&authuser=0&confirm=t&uuid=8913c643-c925-47d5-8f4e-ba033d34e6a5&at=APZUnTXV0hCfcrB3ui7poqoS9M4N:1715794331811' -c -O 'tmp/74667_peaks.bed'
liftOver tmp/74667_peaks.bed -bedPlus=3  Auxiliary_data/hg19/hg38ToHg19.over.chain.gz Auxiliary_data/peaks/CDK8.hg19.peaks.bed.gz tmp/unmapped
wget --header='Host: drive.usercontent.google.com' --header='User-Agent: Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/123.0.0.0 Safari/537.36' --header='Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7' --header='Accept-Language: en-US,en;q=0.9,fr;q=0.8,he;q=0.7,uk;q=0.6,ru;q=0.5' --header='Cookie: HSID=AreFZnLEAYNCZd8lD; SSID=AB0KHDJ_k2CmlfeN2; APISID=w0OR5kHsdu-ViFOc/ADZtMAgS6b409t16Q; SAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; __Secure-1PAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; __Secure-3PAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; SEARCH_SAMESITE=CgQIiZsB; AEC=AQTF6HyN-7_1hZyuKLji2rpC5aGjsaO5E70qWDrQ84fwMKw_w7R-JzTT6WE; __Secure-ENID=19.SE=CKt4IM96j1JmgStIOmsvbq1eHWW-yI7A_A-opmoxGtw1RYvGsocdl4BaC4yGquZ_WIiJUoP7vudMHSRG0R0oCkfDAjn1mZRlHsrPnNPjlJjl3bVynKYoi0mr1coOrv6UtdmvtH6b4sOhKGHmnM83pviMnbBAhMBam5XgFLyeIiChrt3lG-JejR9WK9Zd4_4NT3aQcUE1VIGpNfjBdSLYzSHhEa7qg-ew5zcYgaYVDSTq4guxMqjg1hDR4KUc9U4crWSJGb03tDSntMmYwcSd2W7ASoSqUluAcyUlNQ; SID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSBSnOHQwBQF192QvY95oVJhwACgYKATISAQASFQHGX2Mi_3VhqO3HX710vnTopsYTUBoVAUF8yKocFRIn3AefJYol0psthxtF0076; __Secure-1PSID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSB0Yb863nI---ZVBChR265tQACgYKAUoSAQASFQHGX2MiIVcVPKGRidM0D6GDIA4ulRoVAUF8yKoAfd6uG1D-4hxygbgQjAf30076; __Secure-3PSID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSBcx3LThkwemjJg_so7918AgACgYKAecSAQASFQHGX2Mi0yLhhbEuEoG0UwIbcuenyRoVAUF8yKrEx4mWAY8_FAK8oT-XmIUd0076; NID=514=QItrnCGgAVukXbED86MZYq3wjODqiZDb-0OE_v9hr0QceIawMrZzshZCPN_UpdM_58lTMecfKNusSr67JrpJUau0GoEvTnFKd1BMuGhxd85yuuQAkOnac9RmTfliS_S8Aa72ePK31xkCbDwRqcYXwzsIb0jkgqKXBU-mW5y6z1gEI5Tp5AGbYcz02U6PY--1juDnF_Zx_9CqHo0IYB41CETHT4x6ft6BvEPZS8oeXWqx; __Secure-1PSIDTS=sidts-CjIBLwcBXBGLTGU1VCR9D01W7nrLYZzodQhOQ8qKWGOxFBtiNbKt5MkueqV_t9_uofQfzBAA; __Secure-3PSIDTS=sidts-CjIBLwcBXBGLTGU1VCR9D01W7nrLYZzodQhOQ8qKWGOxFBtiNbKt5MkueqV_t9_uofQfzBAA; SIDCC=AKEyXzXX-46UHtjUpC6yljkzOBYWDakKSZch3aHVT90xjCVThCXXSqWqFcAuwhPwNA1xAoQ_p7c; __Secure-1PSIDCC=AKEyXzVLAfSulOnGZf0EgkLL7HD7VlRPL_I6_aDExVtucn3hJiuqKmJKBY44isrqAIJFH9Qg5xE; __Secure-3PSIDCC=AKEyXzUI7qSGeeOCQWMcLtqhxkcxhBj1s1-N03wGcSzKcm_d1dPKz8Qa8OgGze9bE_rAcIH1lg' --header='Connection: keep-alive' 'https://drive.usercontent.google.com/download?id=1p6H6TaEzHTRsqwmfqR0BYpIzvqwwQtRp&export=download&authuser=0&confirm=t&uuid=b10f654b-ab48-43a9-82ad-0a5878e299c6&at=APZUnTXWzwoYN-oK4w5OPNCs4doF:1715801407702' -c -O 'tmp/74667.bw.bigwig'
CrossMap.py bigwig Auxiliary_data/hg19/hg38ToHg19.over.chain.gz tmp/74667.bw.bigwig Auxiliary_data/bigwigs/MED1.hg19.bigWig
mv Auxiliary_data/bigwigs/MED1.hg19.bigWig.bw Auxiliary_data/bigwigs/MED1.hg19.bigWig
wget --header='Host: drive.usercontent.google.com' --header='User-Agent: Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/123.0.0.0 Safari/537.36' --header='Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7' --header='Accept-Language: en-US,en;q=0.9,fr;q=0.8,he;q=0.7,uk;q=0.6,ru;q=0.5' --header='Cookie: HSID=AreFZnLEAYNCZd8lD; SSID=AB0KHDJ_k2CmlfeN2; APISID=w0OR5kHsdu-ViFOc/ADZtMAgS6b409t16Q; SAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; __Secure-1PAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; __Secure-3PAPISID=4LAymAWcB5d1lXsp/A4SLPLDPqA-mPT_rR; SEARCH_SAMESITE=CgQIiZsB; AEC=AQTF6HyN-7_1hZyuKLji2rpC5aGjsaO5E70qWDrQ84fwMKw_w7R-JzTT6WE; __Secure-ENID=19.SE=CKt4IM96j1JmgStIOmsvbq1eHWW-yI7A_A-opmoxGtw1RYvGsocdl4BaC4yGquZ_WIiJUoP7vudMHSRG0R0oCkfDAjn1mZRlHsrPnNPjlJjl3bVynKYoi0mr1coOrv6UtdmvtH6b4sOhKGHmnM83pviMnbBAhMBam5XgFLyeIiChrt3lG-JejR9WK9Zd4_4NT3aQcUE1VIGpNfjBdSLYzSHhEa7qg-ew5zcYgaYVDSTq4guxMqjg1hDR4KUc9U4crWSJGb03tDSntMmYwcSd2W7ASoSqUluAcyUlNQ; SID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSBSnOHQwBQF192QvY95oVJhwACgYKATISAQASFQHGX2Mi_3VhqO3HX710vnTopsYTUBoVAUF8yKocFRIn3AefJYol0psthxtF0076; __Secure-1PSID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSB0Yb863nI---ZVBChR265tQACgYKAUoSAQASFQHGX2MiIVcVPKGRidM0D6GDIA4ulRoVAUF8yKoAfd6uG1D-4hxygbgQjAf30076; __Secure-3PSID=g.a000jgjZk-sEw8WZQlR7K8xUpALupac4dnRqLHxUoupayhzlgDSBcx3LThkwemjJg_so7918AgACgYKAecSAQASFQHGX2Mi0yLhhbEuEoG0UwIbcuenyRoVAUF8yKrEx4mWAY8_FAK8oT-XmIUd0076; NID=514=o9tNZU5jtMG3Dj4YIejY3AsbT-KAUQL4L9Cj14eGA-RVqo08pPrCNgBfzFq-1M1OKYIhKXI9FR6Wu9DFlsy1Q5GwUzUOaIiH8FyujiHfB-8mhAbK3BYd_9hNshOAiDbTVZ1kwoAUlm08bAWD_W1XaTlPc0fP7dJErB8rZ1RBuhmoGCgmhZL57rOvWv9PD1WWkQdU3sD9ybOEVgQ_mAaAxn7Tnjk_8GUBW1IIEoCu8us0; __Secure-1PSIDTS=sidts-CjIBLwcBXJNJQWKdBrhrrlo1G7BlDA-YYC7aI_6pxpUXHQzXMSmE2FVg5EmbPPESr4toIBAA; __Secure-3PSIDTS=sidts-CjIBLwcBXJNJQWKdBrhrrlo1G7BlDA-YYC7aI_6pxpUXHQzXMSmE2FVg5EmbPPESr4toIBAA; SIDCC=AKEyXzUDLxRdlsET6OEOM08fjJAg6TwgqMtVGBCVtGT0bcDCLqRY5l_1MqaPTF4jrF39i-GVb7w; __Secure-1PSIDCC=AKEyXzVFOkCcxo2_Ta5BsN7Od8otsC1Zal7KwMU71D8G8MHBJGthDKD2kLvKs0ik-khsWO-FwF0; __Secure-3PSIDCC=AKEyXzXncQyn3bwSkfKgPjpS6zOPCrwz3UgLWKfwCQJdNpD3Ob9Hdw60WLKpOUIL5aMrFPy-Fw' --header='Connection: keep-alive' 'https://drive.usercontent.google.com/download?id=1uS1bTuViHk2MtSZPwLlmapOG4TZ0xKre&export=download&authuser=0&confirm=t&uuid=61ed3e01-50c0-4eb8-93bd-3901d904ccba&at=APZUnTXH7VtsXjX7Qa9QyPsb2tXe:1715794302962' -c -O 'tmp/56379_peaks.bed'
liftOver tmp/56379_peaks.bed -bedPlus=3  Auxiliary_data/hg19/hg38ToHg19.over.chain.gz Auxiliary_data/peaks/MED1.hg19.peaks.bed.gz tmp/unmapped
echo "done" >&3
