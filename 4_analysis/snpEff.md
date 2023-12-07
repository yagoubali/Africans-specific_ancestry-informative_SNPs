##### Functional annotation using snpEff

###### 1. Download and install SnpEff and its databases

```bash
# Download and install SnpEff
curl -v -L 'https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip' > snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
## download databases
java -jar snpEff.jar download -v GRCh37.87
java -jar snpEff.jar download -v GRCh37.p13
#dbNSFP_4.1a
wget -c https://snpeff.blob.core.windows.net/databases/dbs/GRCh37/dbNSFP_4.1a/dbNSFP4.1a.txt.gz
wget -c https://snpeff.blob.core.windows.net/databases/dbs/GRCh37/dbNSFP_4.1a/dbNSFP4.1a.txt.gz.tbi
# GWAS catlog
# link to download is here --> http://pcingola.github.io/SnpEff/snpsift/gwascatalog/
```

###### 2. snpEff annotation based on GRCh37.87

```bash
snpEff_bin="analysis/snpEff"
base_outdir="analysis/snpEff"
mkdir -p  "${base_outdir}/overlapped_snps"
mkdir -p  "${base_outdir}/Africans_snps"
vcf_file_africans="analysis/AFR/africans_snps_specific.vcf"
vcf_file_overlapped="analysis/AFR/africans_snps_overlapped.vcf"
## annotations

java -Xmx8g -jar  ${snpEff_bin}/snpEff.jar  GRCh37.87 ${vcf_file_overlapped} > \
 ${base_outdir}/overlapped_snps/snps_GRCh37.87_ann.vcf

 rename 's/snpEff/GRCh37.87/g' snpEff*
 mv GRCh37* ${base_outdir}/overlapped_snps/

java -Xmx8g -jar ${snpEff_bin}/snpEff.jar  GRCh37.87 ${vcf_file_africans} > \
 ${base_outdir}/Africans_snps/snps_GRCh37.87_ann.vcf

 rename 's/snpEff/GRCh37.87/g' snpEff*
 mv GRCh37* ${base_outdir}/Africans_snps/
```

###### 2. snpEff annotation based on GRCh37.p13

````bash
snpEff_bin="analysis/snpEff"
base_outdir="analysis/snpEff"
mkdir -p  "${base_outdir}/overlapped_snps"
mkdir -p  "${base_outdir}/Africans_snps"
vcf_file_africans="analysis/AFR/africans_snps_specific.vcf"
vcf_file_overlapped="analysis/AFR/africans_snps_overlapped.vcf"

#annotation
java -Xmx8g -jar  ${snpEff_bin}/snpEff.jar  GRCh37.p13 ${vcf_file_overlapped} > \
 ${base_outdir}/overlapped_snps/snps_GRCh37.p13_ann.vcf

 rename 's/snpEff/GRCh37.p13/g' snpEff*
 mv GRCh37* ${base_outdir}/overlapped_snps/

java -Xmx8g -jar  ${snpEff_bin}/snpEff.jar  GRCh37.p13 ${vcf_file_africans} >  \
${base_outdir}/Africans_snps/snps_GRCh37.p13_ann.vcf

rename 's/snpEff/GRCh37.p13/g' snpEff*
 mv GRCh37* ${base_outdir}/Africans_snps/
```
######  3. snpEff annotation based on  GWAS catalog
```bash
snpEff_bin="analysis/snpEff"
base_outdir="analysis/snpEff"
mkdir -p  "${base_outdir}/overlapped_snps"
mkdir -p  "${base_outdir}/Africans_snps"
vcf_file_africans="analysis/AFR/africans_snps_specific.vcf"
vcf_file_overlapped="analysis/AFR/africans_snps_overlapped.vcf"

#annotation
mkdir -p  "${base_outdir}/catalog"

vcf_file_africans="analysis/AFR/africans_snps_specific.vcf"
vcf_file_overlapped="analysis/AFR/africans_snps_overlapped.vcf"
cp ${vcf_file_africans}  ${base_outdir}/catalog/
cp ${vcf_file_overlapped} ${base_outdir}/catalog/
gwas_catalog="gwas_catalog_v1.0.2-associations_e110_r2023-11-08.tsv"

java -Xmx8g -jar  ${snpEff_bin}/SnpSift.jar gwasCat \
-db ${snpEff_bin}/data/${gwas_catalog}  \
${base_outdir}/catalog/africans_snps_specific.vcf | \
tee ${base_outdir}/catalog/africans.gwas.vcf

java -Xmx8g -jar  ${snpEff_bin}/SnpSift.jar gwasCat \
 -db ${snpEff_bin}/data/${gwas_catalog}  \
${base_outdir}/catalog/africans_snps_overlapped.vcf | \
tee ${base_outdir}/catalog/overlapped.gwas.vcf

grep -v '#' ${base_outdir}/catalog/africans.gwas.vcf | grep 'GWASCAT' | cut -f 1,2,3,4,5,8 > ${base_outdir}/catalog/africans.gwas

grep -v '#' ${base_outdir}/catalog/overlapped.gwas.vcf | grep 'GWASCAT' | cut -f 1,2,3,4,5,8 > ${base_outdir}/catalog/overlapped.gwas
````
