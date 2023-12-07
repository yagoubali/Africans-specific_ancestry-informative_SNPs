##### MAF distribution

```bash
outdir="analysis/qc_merged_data"
plink_file="merged_data/merged_final_5FID"
plink --bfile ${plink_file} --freq --out ${outdir}/MAF_check \
      --allow-no-sex --make-founders
```
