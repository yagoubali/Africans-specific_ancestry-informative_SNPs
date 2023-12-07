##### Merge cleaned HapMap3 with 1KG and SGDP.

```bash
outdir="merged_data"
SGDP_1KG="${outdir}/1KG_SGDP"
HapMap="preprocess_raw_data/HapMap3/hapmap3_final"

#extract rs ids from final HapMap3 --> 1247678
cut -f 2 ${HapMap}.bim > ${outdir}/HapMap_final.rsids

plink2 --bfile ${SGDP_1KG} \
      --extract ${outdir}/HapMap_final.rsids \
      --make-bed  \
      --out ${outdir}/1KG_SGDP_2merge

 # Merge

plink --bfile ${outdir}/1KG_SGDP_2merge --bmerge ${HapMap} --make-bed \
    --out ${outdir}/merged --allow-no-sex

```