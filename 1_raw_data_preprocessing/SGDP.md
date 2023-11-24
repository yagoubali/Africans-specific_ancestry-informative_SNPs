SGDP data is preprocessed by performing the following steps:

1. Updating IDs
2. Removing individuals that are existed in 1KGP
3. Removing duplicate SNPs,
4. Removing ambiguous SNPs
5. Assigning SNPs names based on the correct SNPs rsIDs,
6. Removing SNPs that their chromosome doses do not match 1KGP. Also Correcting SNPs genomic positions based on 1KGP
7. Filipping alleles based on 1KGP
8. Remove SNPs that are inconsistent with 1KGP after performing the steps (2-7).

## 1. Update samples ID in PLINK fam file

```bash
mkdir -p preprocess_raw_data/SGDP
outdir="preprocess_raw_data/SGDP";
SGDP='cteam_extended.v4.maf0.1perc'
##<old FID>, <old IID>, <new FID>, <new IID>.
sed -i 's/  */\t/g' SGDP/${SGDP}.fam
awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2,$2,$2}' SGDP/${SGDP}.fam > ${outdir}/update.ids
# 345 samples updated.
plink2  --bfile SGDP/${SGDP}  --make-bed --out ${outdir}/SGDP_updated_ids  --update-ids ${outdir}/update.ids
```

## 2. Removing SGDP samples that are overlapped with 1KG samples

# ---> 23 individuals have been removed, remaining 322 samples

```bash
#extract 1KG samples ids
outdir="preprocess_raw_data/SGDP";
path_1KG="preprocess_raw_data/1KG"
cut -d " " -f2 ${path_1KG}/1KG.fam > ${outdir}/1KG.ids
ids_1kg="${outdir}/1KG.ids"
plink2 --bfile ${outdir}/SGDP_updated_ids --remove-fam ${ids_1kg} --make-bed --out ${outdir}/SGDP_wo1KG_ids

## clean old files
rm ${outdir}/SGDP_updated_ids*
```

## 3. checking duplicate variants

```bash
outdir="preprocess_raw_data/SGDP";
plink="${outdir}/SGDP_wo1KG_ids" ## Base name of PLINK binary files from step2

cut -f 2 ${plink}.bim | sort | uniq -d  > ${outdir}/duplicates.ids
# 0 duplicates
### remove duplicate SNPs --> no need to run th code below
plink2 --bfile ${plink} \
--exclude  ${outdir}/duplicates.ids  \
--make-bed \
--out ${outdir}/SGDP_noDup
```

## 4. Remove ambiguous variants

```bash
outdir="preprocess_raw_data/SGDP";
plink="${outdir}/SGDP_wo1KG_ids"

awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" \
  || $5$6 == "AT" || $5$6 == "TA") {print $2}' \
  ${plink}.bim > \
  ${outdir}/SGDP.ac_gt_snps

# Exclude ambiguous SNPs  plink files, i.e, 5291339 ambiguous SNPs
  plink2 --bfile ${plink} \
  --exclude ${outdir}/SGDP.ac_gt_snps \
  --make-bed \
  --out ${outdir}/SGDP_wo_ambiguous_snp
```

## 5. Update RsIds to match with 1KG

```bash
outdir="preprocess_raw_data/SGDP";
bim_1KG="preprocess_raw_data/1KG/1KG_updatedFID.bim"
plink="${outdir}/SGDP_wo_ambiguous_snp";  ## Base name of PLINK binary files from step3

awk 'BEGIN{OFS="\t"}{print $1"_"$4,$2}' ${bim_1KG} > ${outdir}/SGDP_update.rsids

plink2  --bfile ${plink} \
--update-name ${outdir}/SGDP_update.rsids 2 1 \
--make-bed \
--out ${outdir}/SGDP_rsids_updated
```

## 6. SNPs filliping and error corrections, i.e, correction of SNPs position and chromosome errors

```bash
outdir="preprocess_raw_data/SGDP";
plink_1KG="preprocess_raw_data/1KG/1KG_updatedFID"
plink_SGDP="${outdir}/SGDP_rsids_updated"

##remove chromosomes errors ----> 0 errors no need to run the code below
awk 'BEGIN {OFS="\t"; print "rsid","chr_1","chr_2"} FNR==NR {a[$2]=$1; next} \
($2 in a && a[$2] != $1) {print $2,a[$2],$1 }' \
 ${plink_1KG}.bim ${plink_SGDP}.bim \
> \
${outdir}/chr.errors

cut -f1 ${outdir}/chr.errors | sed -e '1d' > ${outdir}/chr.mismatch
plink2 --bfile ${plink_SGDP} \
--exclude ${outdir}/chr.mismatch \
--make-bed  \
--out ${outdir}/SGDP_cleaned_chr

## pos errors ----> 0 errors

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
($2 in a && a[$2] != $4) {print a[$2],$2}' \
 ${plink_1KG}.bim ${plink_SGDP}.bim\
> \
${outdir}/pos.errors


##flip  ----> 97897 rsids to be flipped
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
 ${plink_1KG}.bim ${plink_SGDP}.bim \
 > ${outdir}/flip.rsids

plink --bfile ${plink_SGDP} \
--flip ${outdir}/flip.rsids \
--make-bed \
--out ${outdir}/SGDP_cleaned_flipped \
--allow-no-sex

rm ${plink_SGDP}*
```

## 7. Check for SNPs errors after corrections

```bash
outdir="preprocess_raw_data/SGDP";
plink_1KG="preprocess_raw_data/1KG/1KG_updatedFID"
plink_SGDP="${outdir}/SGDP_cleaned_flipped"


### Check for SNPs mismatch after correction ----> 97897
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
 ${plink_1KG}.bim ${plink_SGDP}.bim  > \
${outdir}/still.mismatch

plink2 --bfile ${plink_SGDP} \
--exclude ${outdir}/still.mismatch \
--make-bed  \
--out ${outdir}/SGDP_clean

rm  ${plink_SGDP}*
```

## 8. Extract final list of SNPs to be included from SGDP

```bash
### extract final SNPs ----> 18487774 variants remaining after main filters.
outdir="preprocess_raw_data/SGDP";
plink_SGDP="${outdir}/SGDP_clean"
cut -f2 ${plink_SGDP}.bim | grep 'rs' > ${outdir}/SGDP_final.rsids
plink2 --bfile ${plink_SGDP} \
--extract ${outdir}/SGDP_final.rsids \
--make-bed  \
--out ${outdir}/SGDP_final

rm  ${plink_SGDP}*

```
