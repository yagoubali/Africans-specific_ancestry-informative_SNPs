##### LD distribution

```bash
outdir="analysis/qc_merged_data"
plink_file="merged_data/merged_final_5FID"
plink --bfile ${plink_file} --r2  --out ${outdir}/LD
```

##### Plot results

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