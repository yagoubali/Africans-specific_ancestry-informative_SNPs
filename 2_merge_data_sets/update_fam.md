##### update fam file of the merged data set.

```bash
outdir="merged_data"
plink_merged_data="${outdir}/merged"
panel_final="${outdir}/Final.panel"
sed -i 's/-/_/g' ${plink_merged_data}.fam
mv ${plink_merged_data}.fam ${plink_merged_data}.fam_old
awk 'BEGIN{OFS="\t"; } {print $1,$2,$3,$4,$5,-9}' ${plink_merged_data}.fam_old \
     > ${plink_merged_data}.fam
```