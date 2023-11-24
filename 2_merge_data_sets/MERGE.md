## 1. Merge 1KG updated FID with SGDP cleaned data. sets

```bash
outdir="merged_data"
dir_1KG="preprocess_raw_data/1KG"
dir_SGDP="preprocess_raw_data/SGDP"

#extract rs ids from cleaned SGPD
cut -f 2 ${dir_SGDP}/SGDP_final.bim > ${outdir}/SGDP_final.rsids

 #  Base name of PLINK binary files of final step of 1KG preprocssing
plink2 --bfile ${dir_1KG}/1KG_updatedFID \
--extract ${outdir}/SGDP_final.rsids \
--make-bed  \
--out ${outdir}/1KG_2merge

plink2 --bfile ${dir_SGDP}/SGDP_final \
--extract ${outdir}/SGDP_final.rsids \
--make-bed  \
--out ${outdir}/SGDP_2merge

 # Merge cleaned SGDP with 1KG

plink --bfile ${outdir}/1KG_2merge  --bmerge ${outdir}/SGDP_2merge --make-bed --out ${outdir}/1KG_SGDP --allow-no-sex
```

## 2. Merge cleaned HapMap3 with 1KG and SGDP.

```bash
outdir="merged_data"
SGDP_1KG="${outdir}/1KG_SGDP"
HapMap="preprocess_raw_data/HapMap3/hapmap3_final"

#extract rs ids from final HapMap3 --> 1247678
cut -f 2 ${HapMap}.bim > ${outdir}/HapMap_final.rsids

 #
plink2 --bfile ${SGDP_1KG} \
--extract ${outdir}/HapMap_final.rsids \
--make-bed  \
--out ${outdir}/1KG_SGDP_2merge


 # Merge

plink --bfile ${outdir}/1KG_SGDP_2merge --bmerge ${HapMap} --make-bed \
--out ${outdir}/merged --allow-no-sex

```

## 3. Prepare samples panels

```bash
outdir="merged_data"
panel_1KG="1KG/integrated_call_samples_v3.20130502.ALL.panel"
panel_SGDP="SGDP/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt"
panel_HapMap="HapMap3/relationships_w_pops_121708.txt"
fam_SGDP_final="preprocess_raw_data/SGDP/SGDP_final.fam"
fam_HapMap_final="preprocess_raw_data/HapMap3/hapmap3_final.fam"

cp ${panel_1KG} ${outdir}/
cut -f4,6,7,8,11  ${panel_SGDP} > ${outdir}/SGDP.panel1
sed -i 's/(//g' ${outdir}/SGDP.panel1
sed -i 's/)//g' ${outdir}/SGDP.panel1
cat  ${outdir}/SGDP.panel1 | sed -e '1d' |  \
awk -f  miscellaneous/SGDP_panel.awk - > ${outdir}/SGDP_ALL.panel

sed -i 's/-/_/g' ${outdir}/SGDP_ALL.panel
cp ${fam_SGDP_final} ${outdir}/
sed -i 's/-/_/g' ${outdir}/SGDP_final.fam

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1]=$1; next} \
($1 in a ) {print $0}' \
${outdir}/SGDP_final.fam ${outdir}/SGDP_ALL.panel > ${outdir}/SGDP_panel.merge

## HapMap3

awk -f  miscellaneous/HapMap3_panel.awk  ${panel_HapMap} > ${outdir}/HapMap_ALL.panel
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1]=$1; next} \
($1 in a ) {print $0}' \
 ${fam_HapMap_final} ${outdir}/HapMap_ALL.panel > ${outdir}/HapMap_panel.merge

##
 cat ${panel_1KG} >  ${outdir}/Final.panel
 cat  ${outdir}/SGDP_panel.merge >>  ${outdir}/Final.panel
 cat ${outdir}/HapMap_panel.merge >>  ${outdir}/Final.panel
```

## 3. update fam file of the merged data set.

```bash
outdir="merged_data"
plink_merged_data="${outdir}/merged"
panel_final="${outdir}/Final.panel"
head ${plink_merged_data}.fam
sed -i 's/-/_/g' ${plink_merged_data}.fam
mv ${plink_merged_data}.fam ${plink_merged_data}.fam_old
 awk 'BEGIN{OFS="\t"; } {print $1,$2,$3,$4,$5,-9}' ${plink_merged_data}.fam_old \
 > ${plink_merged_data}.fam
```

## 3. Remove 2 samples of SGDP that are missed in the sample panel file.

```bash
outdir="merged_data"
plink_merged_data="${outdir}/merged"
touch ${outdir}/remove_2SGDP.samples
echo "HGDP01344" > ${outdir}/remove_2SGDP.samples
echo "IHW9195" >> ${outdir}/remove_2SGDP.samples

plink2 --bfile ${plink_merged_data}  \
--remove-fam  ${outdir}/remove_2SGDP.samples \
--make-bed \
--out ${outdir}/merged_clean

```

## 4. ### Update FID in the cleaned fam file by using super populations ids.

```bash
outdir="merged_data"
plink_merged_data="${outdir}/merged_clean"
panel_final="${outdir}/Final.panel"

## Create a file with the follows
## Old family ID Old within-family ID New family ID New within-family ID
sed -e '1d'  ${panel_final} | awk 'BEGIN{OFS="\t"}  {print $1,$1,$3,$1}' \
> ${outdir}/FID.superpop

 plink2  --bfile ${plink_merged_data} \
  --make-bed \
  --out ${outdir}/merged_final \
  --update-ids ${outdir}/FID.superpop \
  --allow-no-sex

  awk 'BEGIN{OFS="\t"; } {print $1,$2,$2,$2}'  ${outdir}/merged_final.fam \
  > ${outdir}/merged_final_update.fam
```

## 5. Remove non 1KG super populations

```bash
outdir="merged_data"
plink_merged_data="${outdir}/merged_final"
touch ${outdir}/FID.1KG
echo "EUR" > ${outdir}/FID.1KG
echo "EAS" >> ${outdir}/FID.1KG
echo "SAS" >> ${outdir}/FID.1KG
echo "AMR" >> ${outdir}/FID.1KG
echo "AFR" >> ${outdir}/FID.1KG
plink2  --bfile ${plink_merged_data} \
  --make-bed \
  --keep-fam ${outdir}/FID.1KG \
    --out ${outdir}/merged_final_5FID

# 3251 samples remaining.

```
