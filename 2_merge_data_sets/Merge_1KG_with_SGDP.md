##### Merge 1KG updated FID with SGDP cleaned data sets

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

plink --bfile ${outdir}/1KG_2merge  --bmerge ${outdir}/SGDP_2merge \
      --make-bed \
      --out ${outdir}/1KG_SGDP --allow-no-sex
```