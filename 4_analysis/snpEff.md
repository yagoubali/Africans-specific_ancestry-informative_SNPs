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
cd analysis/
snpEff_bin="../snpEff"
base_outdir="../analysis/out_snpEff"
mkdir -p  "${base_outdir}/pairwise_snps"
mkdir -p  "${base_outdir}/global_snps"
vcf_file_global="../analysis/AFR/afr_global_snps.vcf"
vcf_file_pairwise="../analysis/AFR/afr_pairwise_snp.vcf"
## annotations

java -Xmx8g -jar  ${snpEff_bin}/snpEff.jar  GRCh37.87 ${vcf_file_pairwise} > \
 ${base_outdir}/pairwise_snps/pairwise_snps_GRCh37.87_ann.vcf

 rename 's/snpEff/GRCh37.87/g' snpEff*
 mv GRCh37* ${base_outdir}/pairwise_snps/

java -Xmx8g -jar ${snpEff_bin}/snpEff.jar  GRCh37.87 ${vcf_file_global} > \
 ${base_outdir}/global_snps/global_snps_GRCh37.87_ann.vcf

 rename 's/snpEff/GRCh37.87/g' snpEff*
 mv GRCh37* ${base_outdir}/global_snps/
```

###### 3. snpEff annotation based on GRCh37.p13

```bash
snpEff_bin="../snpEff"
base_outdir="../analysis/out_snpEff"
vcf_file_global="../analysis/AFR/afr_global_snps.vcf"
vcf_file_pairwise="../analysis/AFR/afr_pairwise_snp.vcf"

#annotation
java -Xmx8g -jar  ${snpEff_bin}/snpEff.jar  GRCh37.p13 ${vcf_file_global} > \
 ${base_outdir}/global_snps/global_snps_GRCh37.p13_ann.vcf

 rename 's/snpEff/GRCh37.p13/g' snpEff*
 mv GRCh37* ${base_outdir}/global_snps/

java -Xmx8g -jar  ${snpEff_bin}/snpEff.jar  GRCh37.p13 ${vcf_file_pairwise} >  \
${base_outdir}/pairwise_snps/pairwise_snps_GRCh37.p13_ann.vcf

rename 's/snpEff/GRCh37.p13/g' snpEff*
 mv GRCh37* ${base_outdir}/pairwise_snps/
```
######  4. snpEff annotation based on  GWAS catalog
```bash
snpEff_bin="../snpEff"
base_outdir="../analysis/out_snpEff"
vcf_file_global="../analysis/AFR/afr_global_snps.vcf"
vcf_file_pairwise="../analysis/AFR/afr_pairwise_snp.vcf"

#annotation
mkdir -p  "${base_outdir}/catalog"



cp ${vcf_file_global}  ${base_outdir}/catalog/
cp ${vcf_file_pairwise} ${base_outdir}/catalog/
gwas_catalog="gwas_catalog_v1.0.2-associations_e112_r2024-08-12.tsv"

java -Xmx8g -jar  ${snpEff_bin}/SnpSift.jar gwasCat \
-db ${snpEff_bin}/data/${gwas_catalog}  \
${base_outdir}/catalog/afr_global_snps.vcf | \
tee ${base_outdir}/catalog/global.gwas.vcf

java -Xmx8g -jar  ${snpEff_bin}/SnpSift.jar gwasCat \
 -db ${snpEff_bin}/data/${gwas_catalog}  \
${base_outdir}/catalog/afr_pairwise_snp.vcf | \
tee ${base_outdir}/catalog/pairwise.gwas.vcf

grep -v '#' ${base_outdir}/catalog/global.gwas.vcf | grep 'GWASCAT' | cut -f 1,2,3,4,5,8 > ${base_outdir}/catalog/global.gwas

grep -v '#' ${base_outdir}/catalog/pairwise.gwas.vcf | grep 'GWASCAT' | cut -f 1,2,3,4,5,8 > ${base_outdir}/catalog/pairwise.gwas
````
