##### 1. Check individuals' missingness

```bash
mkdir -p analysis/qc_merged_data
outdir="analysis/qc_merged_data"
plink_file="merged_data/merged_final_5FID"
plink --bfile ${plink_file} --missing  --out ${outdir}/missing_all --allow-no-sex
```

##### 2. Plot results

```R
lmiss= read.table("analysis/qc_merged_data/missing_all.lmiss", header=T)
imiss = read.table("analysis/qc_merged_data/missing_all.imiss", header=T)
ggbg2 <- function() {
  points(0,0,pch=16, cex=1e6, col="lightgray")
  grid(col="white", lty=1)
}
png("analysis/qc_merged_data/missingness.png", width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(2,1))
plot(density(imiss$F_MISS), main=" Missingness of individuals" , panel.last=ggbg2(), ylab="Density")
par(new=TRUE)
polygon(density(imiss$F_MISS), col = "slateblue1")

plot(density(lmiss$F_MISS), main=" Missingness of SNPs" , panel.last=ggbg2(), ylab="Density")
par(new=TRUE)
polygon(density(lmiss$F_MISS), col = "slateblue1")
dev.off()
```
