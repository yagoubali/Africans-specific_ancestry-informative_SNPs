##### 1. SNP annotation for hg19 in vcf format

```bash

download_dir="miscellaneous/"  ## a  directory to store Miscellaneous, ie. miscellaneous folder
#wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz -P ${download_dir}
wget -O - "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz" | gunzip -t && echo SUCCESS

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz.tbi -P ${download_dir}
## construct annotations
##cat common_all_20180423.vcf.gz | gzip -d | bgzip -c > snps_common_all_20180423.vcf.gz && tabix snps_common_all_20180423.vcf.gz
```

##### 2. SNP annotation for hg19 in bed format

```bash
#cd project_dir ### the main project folder
mkdir -p miscellaneous/bed
download_dir="miscellaneous/bed"  ## a  directory to store miscellaneous/bed, ie. miscellaneous folder
url="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/BED/"
prefix="bed_chr_";
suffix="bed.gz";
suffix2="bed"

for chrom in {1..22}; do
  file=${url}${prefix}${chrom}.${suffix};
  echo " downloading ... $file ...";
  wget -c ${file} -P ${download_dir}
done
# Prepare downloaded bed files by:
  #.1. remove header line as it is just an info line
  # 2. Rename chromosomes to numbers by removing chr prefix

for chrom in {1..22}; do
  file=${download_dir}/${prefix}${chrom}.${suffix2};
  gzip -d  ${file}.gz
  sed -i '1d' ${file}
  sed -i 's/chr//g' ${file}
  head ${file}
done
```

##### 3. UCSC LiftOver files

```bash
download_dir="miscellaneous/"  ## a  directory to store Miscellaneous, ie. miscellaneous folder
wget -c https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver -P ${download_dir}
chmod +x ${download_dir}liftOver
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz -P ${download_dir}
```
