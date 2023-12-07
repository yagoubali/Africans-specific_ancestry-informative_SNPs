##### Extract Ancestry Informative SNPs

```bash
mkdir -p analysis/AFR
outdir="analysis/AFR"
plink_file="analysis/qc_merged_data/merged_final_5FID_Fst"
afr_snps_overlapped="analysis/Fst/SNPs_list/AFR_0.4.txt";
afr_snps_specific="analysis/Fst/SNPs_list/AFR_specific_0.4.txt";


touch ${outdir}/FID.Afr
echo "AFR" > ${outdir}/FID.Afr

plink --bfile ${plink_file}  \
      --extract ${afr_snps_overlapped}  \
       --keep-fam ${outdir}/FID.Afr  \
      --make-bed   \
       --recode vcf \
      --out ${outdir}/africans_snps_overlapped \
      --allow-no-sex
#10753 variants and 747 people pass filters and QC

plink --bfile ${plink_file}  \
      --extract ${afr_snps_specific}  \
       --keep-fam ${outdir}/FID.Afr  \
      --make-bed   \
       --recode vcf \
      --out ${outdir}/africans_snps_specific \
      --allow-no-sex
#5784 variants and 747 people pass filters and QC.
```
