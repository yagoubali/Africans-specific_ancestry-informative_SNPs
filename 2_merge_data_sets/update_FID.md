##### Update FID in the cleaned fam file by using super populations ids.

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
