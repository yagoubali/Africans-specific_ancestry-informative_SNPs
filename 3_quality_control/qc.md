## 1. Check individuals' missingness

```bash
mkdir -p analysis/qc_merged_data
outdir="analysis/qc_merged_data"
plink_file="merged_data/merged_final_5FID"
plink --bfile ${plink_file} --missing  --out ${outdir}/missing_all --allow-no-sex
```

Plot results

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

## 2. remove misssinig

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

## 3. MAF distribution

```bash
outdir="analysis/qc_merged_data"
plink_file="merged_data/merged_final_5FID"
plink --bfile ${plink_file} --freq --out ${outdir}/MAF_check --allow-no-sex --make-founders
```

## 3. HWE distribution

```bash
outdir="analysis/qc_merged_data"
plink_file="merged_data/merged_final_5FID"
plink --bfile ${plink_file}  --hardy --out ${outdir}/HWE_check --allow-no-sex --make-founders

```

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

## 4. LD distribution

```bash
outdir="analysis/qc_merged_data"
plink_file="merged_data/merged_final_5FID"
plink --bfile ${plink_file} --r2  --out ${outdir}/LD
```

Plot results

```R
ggbg2 <- function() {
  points(0,0,pch=16, cex=1e6, col="lightgray")
  grid(col="white", lty=1)
}
output="analysis/qc_merged_data/"
ld_file=paste0(output, "LD.ld")

ld=read.table(ld_file, header=T)
x=ld[,7]
density_x= density(x)
png(paste0(output, "ld.png"), width = 465, height = 225, units='mm', res = 300)
par(mfrow=c(1,2))
plot (density_x, panel.last=ggbg2(), ylab="Density",main="")
par(new=TRUE)
polygon(density_x, col = "slateblue1")
qqnorm(x, panel.last=ggbg2(), ylab="",main="", xlab="")
par(new=TRUE)
qqnorm(x,main="")
par(new=TRUE)
qqline(x)
dev.off()

```
