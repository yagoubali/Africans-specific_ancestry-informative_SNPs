##### Characterizing ancestry informative SNPs

###### 1. MAF distribution

```bash
outdir="analysis/qc_ancestry_snps"
plink_file_africans="analysis/AFR/africans_snps_specific"
plink_file_overlapped="analysis/AFR/africans_snps_overlapped"


#MAF
plink --bfile ${plink_file_africans} --freq --out ${outdir}/MAF_check_africans --allow-no-sex --make-founders

plink --bfile ${plink_file_overlapped} --freq --out ${outdir}/MAF_check_overlapped --allow-no-sex --make-founders
```

###### 2. HWE distribution

```bash
outdir="analysis/qc_ancestry_snps"
plink_file_africans="analysis/AFR/africans_snps_specific"
plink_file_overlapped="analysis/AFR/africans_snps_overlapped"
plink --bfile ${plink_file_africans}  --hardy \
      --out ${outdir}/HWE_check_africans --allow-no-sex --make-founders

plink --bfile ${plink_file_overlapped}  --hardy \
      --out ${outdir}/HWE_check_overlapped --allow-no-sex --make-founders
```

###### 3. LD distribution

```bash
outdir="analysis/qc_ancestry_snps"
plink_file_africans="analysis/AFR/africans_snps_specific"
plink_file_overlapped="analysis/AFR/africans_snps_overlapped"
plink --bfile ${plink_file_africans} --r2  --out ${outdir}/LD_africans

plink --bfile ${plink_file_overlapped} --r2  --out ${outdir}/LD_overlapped

```

###### 4. Plot results

```R
rm(list = ls())
ggbg2 <- function() {
  points(0,0,pch=16, cex=1e6, col="lightgray")
  grid(col="white", lty=1)
}
output="analysis/qc_ancestry_snps/"
maf_file_africans=paste0(output, "MAF_check_africans.frq")
hwe_file_africans=paste0(output, "HWE_check_africans.hwe")
ld_file_africans=paste0(output, "LD_africans.ld")

maf_file_overlapped=paste0(output, "MAF_check_overlapped.frq")
hwe_file_overlapped=paste0(output, "HWE_check_overlapped.hwe")
ld_file_overlapped=paste0(output, "LD_overlapped.ld")

maf_africans=read.table(maf_file_africans, header=T)
d_africans=density(maf_africans$MAF)
hwe_africans=read.table(hwe_file_africans, header=T)
dd_africans=(hwe_africans$O.HET./hwe_africans$E.HET.)


png(paste0(output, "maf_hwe_africans.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(2,2))
plot(d_africans, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(d_africans, col = "slateblue1")
boxplot(maf_africans$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
points(0,0,pch=16, cex=1e6, col="lightgray")
par(new=TRUE)
boxplot(maf_africans$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
hist(dd_africans, panel.last=ggbg2(), ylab="Density",main="", xlab="Observed hwe/Expected hwe",col = "slateblue1")
hist(hwe_africans$P, xlab="Hardy-Weinberg Equilibrium (HWE) p-value", ylab="Frequency", col = "slateblue1",main="", panel.last=ggbg2())
dev.off()

maf_overlapped=read.table(maf_file_overlapped, header=T)
d_overlapped=density(maf_overlapped$MAF)
hwe_overlapped=read.table(hwe_file_overlapped, header=T)
dd_overlapped=(hwe_overlapped$O.HET./hwe_overlapped$E.HET.)

png(paste0(output, "maf_hwe_overlapped.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(2,2))
plot(d_overlapped, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(d_overlapped, col = "slateblue1")
boxplot(maf_overlapped$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
points(0,0,pch=16, cex=1e6, col="lightgray")
par(new=TRUE)
boxplot(maf_overlapped$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
hist(dd_overlapped, panel.last=ggbg2(), ylab="Density",main="", xlab="Observed hwe/Expected hwe",col = "slateblue1")
hist(hwe_overlapped$P, xlab="Hardy-Weinberg Equilibrium (HWE) p-value", ylab="Frequency", col = "slateblue1",main="", panel.last=ggbg2())
dev.off()

## ld

ld_africans=read.table(ld_file_africans, header=T)
x_africans=ld_africans[,7]
density_x_africans= density(x_africans)
png(paste0(output, "ld_africans.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(1,2))
plot (density_x_africans, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(density_x_africans, col = "slateblue1")
qqnorm(x_africans, panel.last=ggbg2(), ylab="",main="", xlab="")
par(new=TRUE)
qqnorm(x_africans,main="")
par(new=TRUE)
qqline(x_africans)
dev.off()

ld_overlapped=read.table(ld_file_overlapped, header=T)
x_overlapped=ld_overlapped[,7]
density_x_overlapped= density(x_overlapped)
png(paste0(output, "ld_overlapped.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(1,2))
plot (density_x_overlapped, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(density_x_overlapped, col = "slateblue1")
qqnorm(x_overlapped, panel.last=ggbg2(), ylab="",main="", xlab="")
par(new=TRUE)
qqnorm(x_overlapped,main="")
par(new=TRUE)
qqline(x_overlapped)
dev.off()
#summary
summary(maf_africans$MAF)
summary(dd_africans)
summary(maf_overlapped$MAF)
summary(dd_overlapped)
summary(x_africans)
summary(x_overlapped)
```
