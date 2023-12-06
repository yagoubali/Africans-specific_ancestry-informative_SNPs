## 1. 1000 genomes (1KG)

```bash
# Download 1000 genomes vcf files

set -x;
# cd projectDir
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr";
suffix=".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz";

download_dir="1KG"

for chr in {1..22}; do
    wget -c $prefix$chr$suffix  -P ${download_dir};
    wget -c $prefix$chr$suffix.tbi -P ${download_dir};
    gzip -t ${download_dir}/ALL.chr$chr$suffix
done

wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_male_samples_v3.20130502.ALL.panel
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped
```

# Check any downloading error in gz files

```bash
for chrom in {1..22}; do
  file=${prefix}${chrom}.${suffix};
  echo " Testing ... $file ...";

  gzip -t  ${file}
done
```

## 2. HapMap3

```bash

set -x;
download_dir="HapMap3/"  ## a  directory to store hapmap 3 files

wget -c https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/00README.txt   -P ${download_dir}
wget -c https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/hapmap3_r1_b36_fwd.qc.poly.tar.bz2  -P ${download_dir}
wget -c https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/phase_3_samples.doc  -P ${download_dir}
wget -c https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/relationships_w_pops_051208.txt  -P ${download_dir}
wget -c https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/relationships_w_pops_121708.txt  -P ${download_dir}
cd ${download_dir}
tar -xvf hapmap3_r1_b36_fwd.qc.poly.tar.bz2
```

## 3. Simon Genome Diversity Project

```bash

set -x;

download_dir="." # a  directory to store Simon Genome Diversity Project files


wget -c https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.bim.zip   -P ${download_dir}
wget -c https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.bed  -P ${download_dir}
wget -c https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.fam  -P ${download_dir}

wget -c https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt  -P ${download_dir}
wget -c https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.samples.txt -P ${download_dir}

## unzip
unzip cteam_extended.v4.maf0.1perc.bim.zip

```

## 4. Miscellaneous

### 1.1 SNP annotation for hg19 in vcf format

```bash

download_dir="miscellaneous/"  ## a  directory to store Miscellaneous, ie. miscellaneous folder
#wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz -P ${download_dir}
wget -O - "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz" | gunzip -t && echo SUCCESS

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz.tbi -P ${download_dir}
## construct annotations
##cat common_all_20180423.vcf.gz | gzip -d | bgzip -c > snps_common_all_20180423.vcf.gz && tabix snps_common_all_20180423.vcf.gz
```

### 1.1 SNP annotation for hg19 in bed format

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

### 1. UCSC LiftOver files

```bash
download_dir="miscellaneous/"  ## a  directory to store Miscellaneous, ie. miscellaneous folder
wget -c https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver -P ${download_dir}
chmod +x ${download_dir}liftOver
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz -P ${download_dir}
```
