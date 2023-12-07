##### Remove 2 samples of SGDP that are missed in the sample panel file.

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