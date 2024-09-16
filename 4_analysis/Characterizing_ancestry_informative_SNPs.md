##### Characterizing ancestry informative SNPs

###### 1. MAF distribution

```bash
mkdir -p analysis/qc_ancestry_snps
outdir="analysis/qc_ancestry_snps"
plink_file_africans_global="analysis/AFR/afr_global_snps"
plink_file_africans_pairwise="analysis/AFR/afr_pairwise_snp"


#MAF
plink --bfile ${plink_file_africans_global} --freq --out ${outdir}/MAF_check_africans_global --allow-no-sex --make-founders

plink --bfile ${plink_file_africans_pairwise} --freq --out ${outdir}/MAF_check_pairwise --allow-no-sex --make-founders
```

###### 2. HWE distribution

```bash
outdir="analysis/qc_ancestry_snps"
plink_file_africans_global="analysis/AFR/afr_global_snps"
plink_file_africans_pairwise="analysis/AFR/afr_pairwise_snp"
plink --bfile ${plink_file_africans_global}  --hardy \
      --out ${outdir}/HWE_check_africans_global --allow-no-sex --make-founders

plink --bfile ${plink_file_africans_pairwise}  --hardy \
      --out ${outdir}/HWE_check_pairwise --allow-no-sex --make-founders
```

###### 3. LD distribution

```bash
outdir="analysis/qc_ancestry_snps"
plink_file_africans_global="analysis/AFR/afr_global_snps"
plink_file_africans_pairwise="analysis/AFR/afr_pairwise_snp"
plink --bfile ${plink_file_africans_global} --r2  --out ${outdir}/LD_africans_global

plink --bfile ${plink_file_africans_pairwise} --r2  --out ${outdir}/LD_africans_pairwise

```

###### 4. Plot results

```R
rm(list = ls())
ggbg2 <- function() {
  points(0,0,pch=16, cex=1e6, col="lightgray")
  grid(col="white", lty=1)
}
output="analysis/qc_ancestry_snps/"
maf_file_global=paste0(output, "MAF_check_africans_global.frq")
hwe_file_global=paste0(output, "HWE_check_africans_global.hwe")
ld_file_global=paste0(output, "LD_africans_global.ld")

maf_file_pairwise=paste0(output, "MAF_check_pairwise.frq")
hwe_file_pairwise=paste0(output, "HWE_check_pairwise.hwe")
ld_file_pairwise=paste0(output, "LD_africans_pairwise.ld")

maf_africans_global=read.table(maf_file_global, header=T)
d_africans_global=density(maf_africans_global$MAF)
hwe_africans_global=read.table(hwe_file_global, header=T)
dd_africans_global=(hwe_africans_global$O.HET./hwe_africans_global$E.HET.)


png(paste0(output, "maf_hwe_africans_global.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(2,2))
plot(d_africans_global, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(d_africans_global, col = "slateblue1")
boxplot(maf_africans_global$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
points(0,0,pch=16, cex=1e6, col="lightgray")
par(new=TRUE)
boxplot(maf_africans_global$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
hist(dd_africans_global, panel.last=ggbg2(), ylab="Density",main="", xlab="Observed hwe/Expected hwe",col = "slateblue1")
hist(hwe_africans_global$P, xlab="Hardy-Weinberg Equilibrium (HWE) p-value", ylab="Frequency", col = "slateblue1",main="", panel.last=ggbg2())
dev.off()

maf_pairwise=read.table(maf_file_pairwise, header=T)
d_pairwise=density(maf_pairwise$MAF)
hwe_pairwise=read.table(hwe_file_pairwise, header=T)
dd_pairwise=(hwe_pairwise$O.HET./hwe_pairwise$E.HET.)

png(paste0(output, "maf_hwe_pairwise.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(2,2))
plot(d_pairwise, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(d_pairwise, col = "slateblue1")
boxplot(maf_pairwise$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
points(0,0,pch=16, cex=1e6, col="lightgray")
par(new=TRUE)
boxplot(maf_pairwise$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
hist(dd_pairwise, panel.last=ggbg2(), ylab="Density",main="", xlab="Observed hwe/Expected hwe",col = "slateblue1")
hist(hwe_pairwise$P, xlab="Hardy-Weinberg Equilibrium (HWE) p-value", ylab="Frequency", col = "slateblue1",main="", panel.last=ggbg2())
dev.off()

## ld

ld_africans_global=read.table(ld_file_global, header=T)
x_africans_global=ld_africans_global[,7]
density_x_africans_global= density(x_africans_global)
png(paste0(output, "ld_africans_global.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(1,2))
plot (density_x_africans_global, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(density_x_africans_global, col = "slateblue1")
qqnorm(x_africans_global, panel.last=ggbg2(), ylab="",main="", xlab="")
par(new=TRUE)
qqnorm(x_africans_global,main="")
par(new=TRUE)
qqline(x_africans_global)
dev.off()

ld_pairwise=read.table(ld_file_pairwise, header=T)
x_pairwise=ld_pairwise[,7]
density_x_pairwise= density(x_pairwise)
png(paste0(output, "ld_pairwise.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(1,2))
plot (density_x_pairwise, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(density_x_pairwise, col = "slateblue1")
qqnorm(x_pairwise, panel.last=ggbg2(), ylab="",main="", xlab="")
par(new=TRUE)
qqnorm(x_pairwise,main="")
par(new=TRUE)
qqline(x_pairwise)
dev.off()
#summary
summary(maf_africans_global$MAF)
summary(dd_africans_global)
summary(maf_pairwise$MAF)
summary(dd_pairwise)
summary(x_africans_global)
summary(x_pairwise)
```

### Fst 0.6
```bash
mkdir -p analysis/qc_ancestry_snps_0.6
outdir="analysis/qc_ancestry_snps_0.6"
plink_file_africans_global="analysis/AFR_0.6/afr_global_snps_0.6"
plink_file_africans_pairwise="analysis/AFR_0.6/afr_pairwise_snps_0.6"


#MAF

plink --bfile ${plink_file_africans_global} --freq --out ${outdir}/MAF_check_africans_global_0.6 --allow-no-sex --make-founders

plink --bfile ${plink_file_africans_pairwise} --freq --out ${outdir}/MAF_check_pairwise_0.6 --allow-no-sex --make-founders
```
###### 2. HWE distribution

```bash
plink --bfile ${plink_file_africans_global}  --hardy \
      --out ${outdir}/HWE_check_africans_global_0.6 --allow-no-sex --make-founders

plink --bfile ${plink_file_africans_pairwise}  --hardy \
      --out ${outdir}/HWE_check_pairwise_0.6 --allow-no-sex --make-founders
```
###### 3. LD distribution

```bash
plink --bfile ${plink_file_africans_global} --r2  --out ${outdir}/LD_africans_global_0.6

plink --bfile ${plink_file_africans_pairwise} --r2  --out ${outdir}/LD_africans_pairwise_0.6

```

###### 4. Plot results

```R
rm(list = ls())
ggbg2 <- function() {
  points(0,0,pch=16, cex=1e6, col="lightgray")
  grid(col="white", lty=1)
}
output="../qc_ancestry_snps_0.6/"
maf_file_global=paste0(output, "MAF_check_africans_global_0.6.frq")
hwe_file_global=paste0(output, "HWE_check_africans_global_0.6.hwe")
ld_file_global=paste0(output, "LD_africans_global_0.6.ld")

maf_file_pairwise=paste0(output, "MAF_check_pairwise_0.6.frq")
hwe_file_pairwise=paste0(output, "HWE_check_pairwise_0.6.hwe")
ld_file_pairwise=paste0(output, "LD_africans_pairwise_0.6.ld")

maf_africans_global=read.table(maf_file_global, header=T)
d_africans_global=density(maf_africans_global$MAF)
hwe_africans_global=read.table(hwe_file_global, header=T)
dd_africans_global=(hwe_africans_global$O.HET./hwe_africans_global$E.HET.)


png(paste0(output, "maf_hwe_africans_global_0.6.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(2,2))
plot(d_africans_global, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(d_africans_global, col = "slateblue1")
boxplot(maf_africans_global$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
points(0,0,pch=16, cex=1e6, col="lightgray")
par(new=TRUE)
boxplot(maf_africans_global$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
hist(dd_africans_global, panel.last=ggbg2(), ylab="Density",main="", xlab="Observed hwe/Expected hwe",col = "slateblue1")
hist(hwe_africans_global$P, xlab="Hardy-Weinberg Equilibrium (HWE) p-value", ylab="Frequency", col = "slateblue1",main="", panel.last=ggbg2())
dev.off()

maf_pairwise=read.table(maf_file_pairwise, header=T)
d_pairwise=density(maf_pairwise$MAF)
hwe_pairwise=read.table(hwe_file_pairwise, header=T)
dd_pairwise=(hwe_pairwise$O.HET./hwe_pairwise$E.HET.)

png(paste0(output, "maf_hwe_pairwise_0.6.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(2,2))
plot(d_pairwise, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(d_pairwise, col = "slateblue1")
boxplot(maf_pairwise$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
points(0,0,pch=16, cex=1e6, col="lightgray")
par(new=TRUE)
boxplot(maf_pairwise$MAF, horizontal=TRUE,boxwex = 0.25,col = "slateblue1",xlab= "Minor allele frequency (MAF)")
hist(dd_pairwise, panel.last=ggbg2(), ylab="Density",main="", xlab="Observed hwe/Expected hwe",col = "slateblue1")
hist(hwe_pairwise$P, xlab="Hardy-Weinberg Equilibrium (HWE) p-value", ylab="Frequency", col = "slateblue1",main="", panel.last=ggbg2())
dev.off()

## ld

ld_africans_global=read.table(ld_file_global, header=T)
x_africans_global=ld_africans_global[,7]
density_x_africans_global= density(x_africans_global)
png(paste0(output, "ld_africans_global_0.6.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(1,2))
plot (density_x_africans_global, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(density_x_africans_global, col = "slateblue1")
qqnorm(x_africans_global, panel.last=ggbg2(), ylab="",main="", xlab="")
par(new=TRUE)
qqnorm(x_africans_global,main="")
par(new=TRUE)
qqline(x_africans_global)
dev.off()

ld_pairwise=read.table(ld_file_pairwise, header=T)
x_pairwise=ld_pairwise[,7]
density_x_pairwise= density(x_pairwise)
png(paste0(output, "ld_pairwise_0.6.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(1,2))
plot (density_x_pairwise, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(density_x_pairwise, col = "slateblue1")
qqnorm(x_pairwise, panel.last=ggbg2(), ylab="",main="", xlab="")
par(new=TRUE)
qqnorm(x_pairwise,main="")
par(new=TRUE)
qqline(x_pairwise)
dev.off()
summary(maf_africans_global$MAF)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03414 0.08852 0.13220 0.13812 0.18270 0.27180 

summary(dd_africans_global)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6377  0.8941  0.9203  0.9220  0.9529  1.0392 
summary(maf_pairwise$MAF)
#  Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0006693 0.0642600 0.1305000 0.1437264 0.2129000 0.3661000 
summary(dd_pairwise)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4054  0.9332  0.9666  0.9618  0.9982  1.1008 
summary(x_africans_global)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2005  0.5200  0.7154  0.6919  0.8981  1.0000 
summary(x_pairwise)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2000  0.4068  0.6825  0.6584  0.9275  1.0000
```