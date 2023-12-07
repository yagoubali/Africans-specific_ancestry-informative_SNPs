##### Remove non 1KG super populations

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
