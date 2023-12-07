#### HapMap3 is preprocessed by performing the following steps:

1. Convert map and ped files into PLINK binary filesets
2. Merge HapMap files.
3. Removing duplicates and ambiguous SNPs
4. Performing liftOver to convert the SNPs genomic positions from NCBI36 (hg18) into NCBI37 (hg19)
5. Removing individuals that are existed in 1KGP and SGDP
6. Removing SNPs that are not existed in both 1KGP and SGDP,
7. Removing SNPs that their chromosome doses do not match both 1KGP and SGDP.
8. Correcting SNPs genomic positions based on 1KGP and SGDP.
9. Flipping alleles based on 1KGP and SGDP.
10. Remove SNPs that are inconsistent with 1KGP and SGDP after performing the steps (2-9).

##### 1. HapMap3 map & ped files into PLINK binary files: fam, bed, and bim

```bash
#cd PathTO_hapmap3_directory
mkdir -p preprocess_raw_data/HapMap3
outdir="preprocess_raw_data/HapMap3";
hapmap_dir="HapMap3"


files=($(ls ${hapmap_dir}/*.recode.map| sed -e 's/recode.map/recode/g' ))

for ped in "${files[@]}"; do
ped2="$(basename ${ped})"
echo ${ped2};
plink --file ${hapmap_dir}/${ped2} --make-bed --chr 1-22 --out ${outdir}/${ped2} --allow-no-sex
done
```

##### 2. HapMap3 merge PLINK binary files

```bash
outdir="preprocess_raw_data/HapMap3";
ls ${outdir}/*.fam | sed -e 's/.fam//g' | sed -e '1d' > ${outdir}/merge.list

plink --bfile ${outdir}/hapmap3_r1_b36_fwd.ASW.qc.poly.recode \
  --merge-list ${outdir}/merge.list \
  --out ${outdir}/hapmap3  --make-bed --allow-no-sex
```

##### 3. HapMap3 Remove duplicate SNPs and Ambiguous SNPs

```bash
outdir="preprocess_raw_data/HapMap3";
plink_hapmap3="${outdir}/hapmap3" # Base name of PLINK  binary files from step 2

cut -f 2 ${plink_hapmap3}.bim | sort | uniq -d  > ${outdir}/duplicates.ids

### How to remove duplicate SNPs 0, duplicates
##
plink2 --bfile ${plink_hapmap3} \
--exclude  ${outdir}/duplicates.ids  \
--make-bed \
--out ${outdir}/hapmap3_noDup

rm ${plink_hapmap3}.bim ${plink_hapmap3}.bed ${plink_hapmap3}.fam

awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" \
  || $5$6 == "AT" || $5$6 == "TA") {print $2}' \
  ${outdir}/hapmap3_noDup.bim > \
  ${outdir}/hapmap3.ac_gt_snps

# 2. Exclude ambiguous SNPs  plink files, 126385 snps
  plink2 --bfile ${outdir}/hapmap3_noDup \
  --exclude ${outdir}/hapmap3.ac_gt_snps \
  --make-bed \
  --out ${outdir}/hapmap3_wo_ambiguous_snps

  rm ${outdir}/hapmap3_noDup*
```

##### 4. HapMap3 LiftOver

```bash
outdir="preprocess_raw_data/HapMap3";
plink_hapmap3="${outdir}/hapmap3_wo_ambiguous_snps" ## Base file name of hapmap3_no_ambiguous_snps
awk '{print "chr" $1, $4-1, $4, $2 }' ${plink_hapmap3}.bim > ${outdir}/hapmap3.liftover
chain='miscellaneous/hg18ToHg19.over.chain.gz'
miscellaneous/liftOver ${outdir}/hapmap3.liftover $chain ${outdir}/hapmap.hg19 ${outdir}/hapmap.unmapped

# extract mapped variants
awk '{print $4}' ${outdir}/hapmap.hg19 \
> ${outdir}/hapmap.hg19_snps

# extract updated positions
awk '{print $4, $3}'  ${outdir}/hapmap.hg19 \
 > ${outdir}/hapmap.hg19_pos

#﻿Update the reference data

plink2 --bfile ${plink_hapmap3}  \
--extract ${outdir}/hapmap.hg19_snps \
--update-map ${outdir}/hapmap.hg19_pos  \
--make-bed \
--out ${outdir}/hapmap3_hg19

rm ${plink_hapmap3}*
```

##### 5. Update HapMap3 IDs and remove overlapped samples with samples from merged 1KG and SGDP

```bash
###Update FID in hapmap
#If there's no recognized header line, the file body must contain exactly two or four columns.
# If it has two, it's interpreted as <old IID>, <new IID>.
# If it has four, it's interpreted as <old FID>, <old IID>, <new FID>, <new IID>.
outdir="preprocess_raw_data/HapMap3";
plink_hapmap3="${outdir}/hapmap3_hg19"
SGDP_1KG="merged_data/1KG_SGDP"
cut -d " " -f1 ${SGDP_1KG}.fam > ${outdir}/SGDP_1KG.ids

 # Base name of PLINK  binary files from step 2

awk '{print $1"\t"$2"\t"$2"\t"$2}' ${plink_hapmap3}.fam > ${outdir}/hapmap3_hg19.fid

plink2 --bfile ${plink_hapmap3} --update-ids ${outdir}/hapmap3_hg19.fid \
  --make-bed --out ${outdir}/hapmap3_hg19_updatedFID



#plink2 --bfile ${plink}_FID_updated \
#  --keep-fam   ${outdir}/SGDP_1KG.ids \
#  --make-bed \
#  --out  ${outdir}/hapmap3_with_1KG

plink2 --bfile ${outdir}/hapmap3_hg19_updatedFID \
 --remove-fam  ${outdir}/SGDP_1KG.ids \
  --make-bed \
  --out ${outdir}/hapmap3_cleanIds

  rm ${outdir}/hapmap3_hg19_updatedFID*
  rm ${plink_hapmap3}*{bim,bed,fam,log}

```

##### 6. Removing SNPs that are not existed in both 1KGP and SGDP,

```bash
outdir="preprocess_raw_data/HapMap3";
plink_hapmap3="${outdir}/hapmap3_cleanIds"
SGDP_1KG="merged_data/1KG_SGDP"
cut -f 2 ${SGDP_1KG}.bim > ${outdir}/SGDP_1KG.rsids

plink2 --bfile ${plink_hapmap3} \
--extract ${outdir}/SGDP_1KG.rsids \
--make-bed  \
--out ${outdir}/hapmap3_cleaned_rsids

rm ${plink_hapmap3}*{bim,bed,fam,log}
```

##### 7. Removing SNPs that their chromosome doses do not match both 1KGP and SGDP.

```bash
outdir="preprocess_raw_data/HapMap3";
plink_hapmap3="${outdir}/hapmap3_cleaned_rsids"
SGDP_1KG="merged_data/1KG_SGDP"
awk 'BEGIN {OFS="\t"; print "rsid","chr_1KG","chr_hapmap"} FNR==NR {a[$2]=$1; next} \
($2 in a && a[$2] != $1) {print $2,a[$2],$1 }' \
 ${SGDP_1KG}.bim ${plink_hapmap3}.bim \
> \
${outdir}/hapmap3_chr.errors ## 1 error
#rsid	chr_1KG	chr_hapmap
#rs1640558	7	10
cut -f1 ${outdir}/hapmap3_chr.errors | sed -e '1d' > ${outdir}/chr.mismatch

plink2 --bfile ${plink_hapmap3} \
--exclude ${outdir}/chr.mismatch \
--make-bed  \
--out ${outdir}/hapmap3_cleaned_chr

rm ${plink_hapmap3}*{bim,bed,fam,log}
```

##### 8. Correcting SNPs genomic positions based on 1KGP and SGDP.

```bash
outdir="preprocess_raw_data/HapMap3";
plink_hapmap3="${outdir}/hapmap3_cleaned_chr"
SGDP_1KG="merged_data/1KG_SGDP"

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
($2 in a && a[$2] != $4) {print a[$2],$2}' \
  ${SGDP_1KG}.bim ${plink_hapmap3}.bim \
> \
${outdir}/hapmap3_pos.update

## 15 SNPs
plink --bfile ${plink_hapmap3} \
--update-map ${outdir}/hapmap3_pos.update 1 2 \
--make-bed \
--out ${outdir}/hapmap3_cleaned_pos

rm ${plink_hapmap3}*{bim,bed,fam,log}
```

##### 9. Flipping alleles based on 1KGP and SGDP.

```bash
outdir="preprocess_raw_data/HapMap3";
plink_hapmap3="${outdir}/hapmap3_cleaned_pos"
SGDP_1KG="merged_data/1KG_SGDP"

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
  ${SGDP_1KG}.bim ${plink_hapmap3}.bim \
> ${outdir}/hapmap3_flip.rsids

## 668 SNPs

plink --bfile ${plink_hapmap3} \
--flip ${outdir}/hapmap3_flip.rsids \
--make-bed --allow-no-sex \
--out ${outdir}/hapmap3_flipped

rm ${plink_hapmap3}*{bim,bed,fam,log}
```

##### 10. Remove SNPs that are inconsistent with 1KGP and SGDP after performing the steps (2-9).

```bash
outdir="preprocess_raw_data/HapMap3";
plink_hapmap3="${outdir}/hapmap3_flipped"
SGDP_1KG="merged_data/1KG_SGDP"

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
 ${SGDP_1KG}.bim ${plink_hapmap3}.bim  > \
${outdir}/hapmap3_mismatch.still

## 64 still mismatch
plink2 --bfile ${plink_hapmap3} \
--exclude ${outdir}/hapmap3_mismatch.still \
--make-bed  \
--out ${outdir}/hapmap3_final

rm ${plink_hapmap3}*{bim,bed,fam,log}
```
