#### The preprocessing steps for 1KGP data include:

1. Normalize SNPs by splitting multiallelic sites into biallelic records
2. Extract SNPs from the raw VCF files.
3. Assign and update input VCF files with the correct SNPs rsIDs.
4. Convert VCF files format to PLINK format
5. Merge per chromosome PLINK files
6. Remove ambiguous SNPs.
7. Remove duplicate SNPs

##### 1. Normalize SNPs

```bash

mkdir -p preprocess_raw_data/1KG
outdir="preprocess_raw_data/1KG";
KG_dir="1KG";
prefix="ALL.chr";
suffix="phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf";
for chrom in {1..22}; do
  file=${prefix}${chrom}.${suffix};
  echo " Processing ... $file ...";
  echo "Step 1:  unzip ${file}.gz ...";
  gzip -d -k ${KG_dir}/${file}.gz
  echo "Step 2:  Normalize SNPs in ${file}.gz ...";
  bcftools norm -m-any ${KG_dir}/${file} -O z  -o ${outdir}/${prefix}${chrom}_normalized.vcf.gz

  echo "Step 3:  cleaning the intermediate extracted file"
  rm ${KG_dir}/${file}
  echo "............................."
done

```

##### 2. Extract SNPs from VCF files

```bash
#cd project_dir ### the main project folder

outdir="preprocess_raw_data/1KG";
prefix="ALL.chr";
suffix="_normalized.vcf.gz";

for chrom in {1..22}; do
  file=${prefix}${chrom}${suffix};
  echo "Step 1:  extract SNPs from ${file} ...";
  bcftools view --output-type v --output ${outdir}/${prefix}${chrom}_snps.vcf   \
   --threads 4 --min-alleles 2  --max-alleles 2 --types snps  ${outdir}/${file}

  echo "Step 2:  cleaning the intermediate extracted file"
  rm  ${outdir}/${file}
  echo "............................."
done
```

##### 3. Annotate SNPs: Assign and update input VCF files with the correct SNPs rsIDs

```bash
#cd project_dir ### the main project folder
#mkdir -p preprocess_raw_data/1KG
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz

outdir="preprocess_raw_data/1KG";
#bed_dir="miscellaneous/bed";
prefix="ALL.chr";
suffix="_snps.vcf";
suffix2="_annotated.vcf";

for chrom in {1..22}; do
  file=${prefix}${chrom}${suffix};
  file_out=${prefix}${chrom}${suffix2};
  echo " Annotate ... $file ...";
  echo "Step 1:  bgzip ${file} ...";
  bgzip -c -f  ${outdir}/${file} >  ${outdir}/${file}.gz

  echo "Step 2:  index ${file}.gz ...";
  tabix   ${outdir}/${file}.gz

  echo "Step 3:  Assign and update input VCF files with the correct SNPs rsIDs ...";

  #bcftools annotate -c CHROM,FROM,TO,ID -a ${bed_dir}/bed_chr_${chrom}.bed -o ${outdir}/${file_out} ${outdir}/${file}.gz

   bcftools annotate -c CHROM,FROM,TO,ID -a All_20180423.vcf.gz -o ${outdir}/${file_out} \
  ${outdir}/${file}.gz


  echo "cleaning the intermediate extracted file"
  #rm   ${outdir}/${file}*
  echo "............................."
done
  rm ${outdir}/*snps.vcf*
```

#### 4. Convert VCF files to PLINK format

```bash
#cd project_dir ### the main project folder
#mkdir -p preprocess_raw_data/1KG


outdir="preprocess_raw_data/1KG";
prefix="ALL.chr";
suffix="_annotated.vcf";

for chrom in {1..22}; do
  file=${prefix}${chrom}${suffix};
  file_out=${prefix}${chrom}
  echo "Convert VCF to PLINK $file  ... ";
  plink2 --vcf ${outdir}/${file}    --max-alleles 2 --make-bed --rm-dup exclude-all --out  ${outdir}/${file_out}

  echo "cleaning the intermediate extracted file"
  #rm   ${outdir}/${file}
  echo "............................."
done
rm  ${outdir}/*log
rm  ${outdir}/*vcf
```

##### 5. Merge per chromosome PLINK files
#Warning: Variants 'rs3883815' and 'rs1381052652' have the same position.
#Warning: Variants 'rs587637745' and 'rs1273466016' have the same position.
#Warning: Variants 'rs1334470708' and 'rs12082723' have the same position.

```bash
outdir="preprocess_raw_data/1KG";
prefix="ALL.chr";
echo "step 1: remove SNPs without rs ids"
touch ${outdir}/snp_list
echo "." > ${outdir}/snp_list
for chrom in {1..22}; do
mydata=${prefix}${chrom}
filter_data=${prefix}${chrom}_cleaned
plink --bfile ${outdir}/${mydata} --exclude ${outdir}/snp_list --allow-no-sex  --make-bed --out ${outdir}/${filter_data}
#rm ${mydata}.bim ${mydata}.bed ${mydata}.fam
done

ls ${outdir}/*_cleaned.fam | sed -e 's/.fam//' > ${outdir}/merge1KG.list
sed -i '1d' ${outdir}/merge1KG.list # To delete first line i.e ALL.chr10

plink   --bfile ${outdir}/ALL.chr10_cleaned \
  --make-bed \
  --merge-list ${outdir}/merge1KG.list \
  --out ${outdir}/1KG

  echo "cleaning the intermediate extracted file"
  rm   ${outdir}/${prefix}*
  echo "............................."

rm  preprocess_raw_data/1KG/ALL.chr*

```

##### 6. Remove ambiguous SNPs

```bash
#### 11695664 preprocess_raw_data/1KG/1KG.ac_gt_snps
bash
outdir="preprocess_raw_data/1KG";
prefix_plink="1KG"

awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" \
  || $5$6 == "AT" || $5$6 == "TA") {print $2}' \
   ${outdir}/${prefix_plink}.bim > \
   ${outdir}/${prefix_plink}.ac_gt_snps

  # old -->  Exclude ambiguous SNPs  plink files ---> 11918284 SNPs, i.e, 64604251 variants remaining
  # Sep 6,24 --exclude: 64554313 variants remaining. 64554313 variants remaining after main filters
plink2 --bfile  ${outdir}/${prefix_plink} \
  --exclude  ${outdir}/${prefix_plink}.ac_gt_snps \
  --make-bed \
  --out ${outdir}/${prefix_plink}_no_ambiguous_snps
```

##### 6. Remove duplicate SNPs

```bash
outdir="preprocess_raw_data/1KG";
prefix_plink="1KG_no_ambiguous_snps"
cut -f 2 ${outdir}/${prefix_plink}.bim | sort | uniq -d > ${outdir}/duplicates.ids

### How to remove duplicate SNPs, No need to run the code below as 0 duplicates
plink2 --bfile ${outdir}/${prefix_plink} \
--exclude  ${outdir}/duplicates.ids  \
--make-bed \
--out ${outdir}/${prefix_plink}_noDup

rm ${outdir}/${prefix_plink}.{bed,bim,fam}
```

##### 7. Update FID in fam file

```bash
outdir="preprocess_raw_data/1KG";
prefix_plink="1KG_no_ambiguous_snps_noDup"
awk 'BEGIN{OFS="\t"; } {print $1,$2,$2,$2}' \
${outdir}/${prefix_plink}.fam > ${outdir}/update.fam

plink2  --bfile ${outdir}/${prefix_plink} \
  --make-bed \
  --out  ${outdir}/1KG_updatedFID \
  --update-ids  ${outdir}/update.fam \
  --allow-no-sex

rm ${outdir}/${prefix_plink}*

```
