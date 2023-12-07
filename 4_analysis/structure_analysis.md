##### Population structure analaysis

###### 1. PCA analysis using PLINK

```bash
outdir="analysis/structure"
mkdir -p ${outdir}
plink_merged="analysis/qc_merged_data/merged_final_5FID_Fst"
bim_africans="analysis/AFR/africans_snps_specific.bim"
bim_overlapped="analysis/AFR/africans_snps_overlapped.bim"
cut -f2  ${bim_africans} > ${outdir}/snps.africans
cut -f2  ${bim_overlapped} > ${outdir}/snps.overlapped

#extract
# from merged that passed QC
#2773 samples (102 females, 136 males, 2535 ambiguous; 2773 founders) loaded

plink2 --bfile ${plink_merged} \
--extract ${outdir}/snps.africans \
--make-bed  \
--out ${outdir}/snps.africans


plink2 --bfile ${plink_merged} \
--extract ${outdir}/snps.overlapped \
--make-bed  \
--out ${outdir}/snps.overlapped

#PCA
plink2 --bfile ${outdir}/snps.africans --pca --out ${outdir}/pca_africans
plink2 --bfile ${outdir}/snps.overlapped --pca --out ${outdir}/pca_overlapped
```

###### 2. Plot PCA results

```R
rm(list = ls())
outdir="analysis/structure/"
eigenvec_africans <-read.table(
  paste0(outdir,"pca_africans.eigenvec"),header=F)

eigenvec_overlapped <-read.table(
  paste0(outdir,"pca_overlapped.eigenvec"),header=F)

eigenval_africans <-read.table(
    paste0(outdir,"pca_africans.eigenval"),header=F)
eigenval_overlapped <-read.table(
    paste0(outdir,"pca_overlapped.eigenval"),header=F)


ggbg2 <- function() {
  points(0,0,pch=16, cex=1e6, col="lightgray")
  grid(col="white", lty=1)
}

plotting_pca=function(eigenvec, eigenval, outputfile){
    out_png=paste0(outdir,outputfile)
    id_orders=order(eigenvec$V1)
    png(paste0(outdir,outputfile,"_eigenvec.png"), width = 465, height = 225,
         units='mm', res = 300)
     plot(eigenvec$V3[id_orders],eigenvec$V4[id_orders],
          col=factor(eigenvec$V1[id_orders]),
          panel.first=ggbg2(),xlab="PC1", ylab="PC2",
          main= "Global ancestry estimation"
          )
     legend(x = "topright",
          legend=levels(factor(eigenvec$V1[id_orders])),
          fill=unique(factor(eigenvec$V1[id_orders])),
          lwd=0, cex = 1, xpd = TRUE, ncol = 2)
        dev.off()
    png(paste0(outdir,outputfile,"_eigenval.png"), width = 465, height = 225,
         units='mm', res = 300)
    plot(eigenval$V1,
         panel.first=ggbg2(),
         xlab="Eigen value index",
         ylab="Eigen value", type='l')
    dev.off()
}

plotting_pca(eigenvec_africans, eigenval_africans,
         "pca_africans")
plotting_pca(eigenvec_overlapped, eigenval_overlapped,
      "pca_overlapped")
```

###### 3. Population stratifications using Structure

```bash
outdir="analysis/structure/"
plink_africans="${outdir}/snps.africans.bed"
plink_overlapped="${outdir}/snps.overlapped.bed"
admixture -s 100 ${plink_africans}  5
mv snps.africans* ${outdir}/
admixture -s 100 ${plink_overlapped}  5
mv snps.overlapped* ${outdir}/
```

###### 4. plot population stratifications results

```R
rm(list = ls())
outdir="analysis/structure/"
#plink_africans="${outdir}/snps.africans.bed"
#plink_overlapped="${outdir}/snps.overlapped.bed"
Q_african=read.table(paste0(outdir,"snps.africans.5.Q")
     ,header=F)
fam_african=read.table(paste0(outdir,"snps.africans.fam")
    ,header=F)

Q_overlapped=read.table(paste0(outdir,"snps.overlapped.5.Q")
     ,header=F)
fam_overlapped=read.table(paste0(outdir,"snps.overlapped.fam")
    ,header=F)
plot_structure=function(Q, fam, outfile){
     out_png=paste0(outdir,outfile)
     png(out_png, width = 465, height = 225, units='mm', res = 300)
     ordered_index_groups <- order(fam$V1)
     populations = as.character(fam$V1[ordered_index_groups])
     populations_names=unique(populations)
     populations_labels=rep("",length(populations))
     j=0;
     pos=c()
     for (i in 1:length(populations_names)){
          #cat(populations_names[i])
          index_label=length(which(populations==populations_names[i]))
          populations_labels[j+round(index_label/2, digits = 0)] <- populations_names[i]
          pos=c(pos,(j+index_label/2))
            j=j+index_label
           #cat(j)
     }
    par(mar = c(4, 4, 4, 1))
    barplot(t(as.matrix(Q[ordered_index_groups,])),col=rainbow(5),
        space=0,border=NA,  ylab="Ancestry proportions",xaxt="n")
    axis(1,at = pos,labels =populations_names, lwd = 0,las=1)
    par(xpd=TRUE)
    legend(1, 1.09,
        legend=populations_names,
        fill=rainbow(5),lwd=0, cex = 1, xpd = TRUE, ncol = 5)
    dev.off()
}

plot_structure(Q_african,fam_african,"structure_african.png"  )
plot_structure(Q_overlapped, fam_overlapped, "structure_overlapped.png")
```
