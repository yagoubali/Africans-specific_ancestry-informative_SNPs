##### Extract Ancestry Informative SNPs

```bash
mkdir -p analysis/AFR
outdir="analysis/AFR"
plink_file="analysis/qc_merged_data/merged_final_5FID_Fst"
afr_global_snps="analysis/Fst/SNPs_list/global_AFR_0.4.txt";
afr_pairwise_snps="analysis/Fst/SNPs_list/pairwise_AFR_SNPs_0.4.txt";

touch ${outdir}/FID.Afr
echo "AFR" > ${outdir}/FID.Afr

plink --bfile ${plink_file}  \
      --extract ${afr_global_snps}  \
       --keep-fam ${outdir}/FID.Afr  \
      --make-bed   \
       --recode vcf \
      --out ${outdir}/afr_global_snps\
      --allow-no-sex
#10753 variants and 747 people pass filters and QC

## Sep 5, 25 --->
# --extract: 11372 variants remaining.
# --keep-fam: 747 people remaining.


plink --bfile ${plink_file}  \
      --extract ${afr_pairwise_snps}  \
       --keep-fam ${outdir}/FID.Afr  \
      --make-bed   \
       --recode vcf \
      --out ${outdir}/afr_pairwise_snp \
      --allow-no-sex
#5784 variants and 747 people pass filters and QC.

## Sep 5, 24 ---> 

# --extract: 68126 variants remaining.
# --keep-fam: 747 people remaining.


## Fst 0.6

mkdir -p analysis/AFR_0.6
outdir="analysis/AFR_0.6"
plink_file="analysis/qc_merged_data/merged_final_5FID_Fst"
afr_global_snps="analysis/Fst_0.6/SNPs_list/global_AFR_0.6.txt";
afr_pairwise_snps="analysis/Fst_0.6/SNPs_list/pairwise_AFR_SNPs_0.6.txt";

touch ${outdir}/FID.Afr
echo "AFR" > ${outdir}/FID.Afr

plink --bfile ${plink_file}  \
      --extract ${afr_global_snps}  \
       --keep-fam ${outdir}/FID.Afr  \
      --make-bed   \
       --recode vcf \
      --out ${outdir}/afr_global_snps_0.6\
      --allow-no-sex
#10753 variants and 747 people pass filters and QC

## Sep 5, 25 --->
# --extract: 11372 variants remaining.
# --keep-fam: 747 people remaining.


plink --bfile ${plink_file}  \
      --extract ${afr_pairwise_snps}  \
       --keep-fam ${outdir}/FID.Afr  \
      --make-bed   \
       --recode vcf \
      --out ${outdir}/afr_pairwise_snps_0.6 \
      --allow-no-sex
#5784 variants and 747 people pass filters and QC.

## Sep 5, 24 ---> 

# --extract: 68126 variants remaining.
# --keep-fam: 747 people remaining.
```
