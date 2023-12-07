##### 1. 1000 genomes (1KG)

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

##### 2. Check any downloading error in gz files

```bash
set -x;
# cd projectDir/

prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr";
suffix=".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz";

download_dir="1KG"

for chrom in {1..22}; do
  file="${download_dir}/${prefix}${chrom}.${suffix}";
  echo "Testing ... $file ...";
  gzip -t  ${file}
done
```
