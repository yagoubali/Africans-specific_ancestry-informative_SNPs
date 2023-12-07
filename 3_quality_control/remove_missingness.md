##### Remove misssinig

```bash
outdir="analysis/qc_merged_data"
plink_file="merged_data/merged_final_5FID"

## We recommend to first filter SNPs and individuals based on a relaxed threshold (0.2; >20%)
## Then a filter with a more stringent threshold can be applied (0.02).  -----> this

plink2 --bfile ${plink_file} --geno 0.02 --mind 0.02 --make-bed  \
      --out ${outdir}/merged_final_5FID_Fst

#478 samples removed due to missing genotype data (--mind).
#--geno: 53 variants removed due to missing genotype data.
```
