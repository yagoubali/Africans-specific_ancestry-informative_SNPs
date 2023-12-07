##### HWE distribution

```bash
outdir="analysis/qc_merged_data"
plink_file="merged_data/merged_final_5FID"
plink --bfile ${plink_file}  --hardy \
      --out ${outdir}/HWE_check \
      --allow-no-sex --make-founders
```

##### Plot Results
Plot results

```R
ggbg2 <- function() {
  points(0,0,pch=16, cex=1e6, col="lightgray")
  grid(col="white", lty=1)
}
output="analysis/qc_merged_data/"
maf_file=paste0(output, "MAF_check.frq")
hwf_file=paste0(output, "HWE_check.hwe")
maf=read.table(maf_file, header=T)
d=density(maf$MAF)
hwe=read.table( hwf_file, header=T)
dd=(hwe$O.HET./hwe$E.HET.)
summary(maf$MAF)
summary(dd)
png(paste0(output, "maf_hwe.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(2,2))
plot(d, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(d, col = "slateblue1")
boxplot(maf$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
points(0,0,pch=16, cex=1e6, col="lightgray")
par(new=TRUE)
boxplot(maf$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
hist(dd, panel.last=ggbg2(), ylab="Density",main="", xlab="Observed hwe/Expected hwe",col = "slateblue1")
hist(hwe$P, xlab="Hardy-Weinberg Equilibrium (HWE) p-value", ylab="Frequency", col = "slateblue1",main="", panel.last=ggbg2())
dev.off()

```