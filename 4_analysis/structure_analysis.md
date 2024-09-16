##### Population structure analaysis

###### 1. PCA analysis using PLINK

```bash
outdir="analysis/structure"
mkdir -p ${outdir}
plink_merged="analysis/qc_merged_data/merged_final_5FID_Fst"
bim_global="analysis/AFR/afr_global_snps.bim"
bim_pairwise="analysis/AFR/afr_pairwise_snp.bim"
cut -f2  ${bim_global} > ${outdir}/snps.global
cut -f2  ${bim_pairwise} > ${outdir}/snps.pairwise

#extract
# from merged that passed QC
# 2773 samples (102 females, 136 males, 2535 ambiguous; 2773 founders) loaded

plink2 --bfile ${plink_merged} \
--extract ${outdir}/snps.global \
--make-bed  \
--out ${outdir}/snps.global


plink2 --bfile ${plink_merged} \
--extract ${outdir}/snps.pairwise \
--make-bed  \
--out ${outdir}/snps.pairwise

#PCA
plink2 --bfile ${outdir}/snps.global --pca --out ${outdir}/pca_global
plink2 --bfile ${outdir}/snps.pairwise --pca --out ${outdir}/pca_pairwise
```

###### 2. Plot PCA results

```R
rm(list = ls())
outdir="analysis/structure/"
eigenvec_global <-read.table(
  paste0(outdir,"pca_global.eigenvec"),header=F)

eigenvec_pairwise <-read.table(
  paste0(outdir,"pca_pairwise.eigenvec"),header=F)

eigenval_global <-read.table(
    paste0(outdir,"pca_global.eigenval"),header=F)
eigenval_pairwise <-read.table(
    paste0(outdir,"pca_pairwise.eigenval"),header=F)


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

plotting_pca(eigenvec_global, eigenval_global,
         "pca_global")
plotting_pca(eigenvec_pairwise, eigenval_pairwise,
      "pca_pairwise")
```

###### 3. Population stratifications using Structure

```bash
outdir="analysis/structure/"
plink_global="${outdir}/snps.global.bed"
plink_pairwise="${outdir}/snps.pairwise.bed"
admixture -s 100 ${plink_global}  5
mv snps.global* ${outdir}/
admixture -s 100 ${plink_pairwise}  5
mv snps.pairwise* ${outdir}/
```

###### 4. plot population stratifications results

```R
rm(list = ls())
outdir="analysis/structure/"
#plink_africans="${outdir}/snps.africans.bed"
#plink_overlapped="${outdir}/snps.overlapped.bed"
Q_global=read.table(paste0(outdir,"snps.global.5.Q")
     ,header=F)
P_global=read.table(paste0(outdir,"snps.global.5.P")
     ,header=F)
fam_global=read.table(paste0(outdir,"snps.global.fam")
    ,header=F)

Q_pairwise=read.table(paste0(outdir,"snps.pairwise.5.Q")
     ,header=F)
P_pairwise=read.table(paste0(outdir,"snps.pairwise.5.P")
     ,header=F)
fam_pairwise=read.table(paste0(outdir,"snps.pairwise.fam")
    ,header=F)
plot_structure=function(Q, fam, outfile){
     out_png=paste0(outdir,outfile)
     png(out_png, width = 465, height = 225, units='mm', res = 300)
     ordered_index_groups <- order(fam$V1)
     populations = as.character(fam$V1[ordered_index_groups])
     populations_names= unique(populations)
     population_names_in_matrix= unique(fam[,1])  #
     populations_i=c()
     for (i in 1:length(populations_names)){
          populations_i[i]= which(population_names_in_matrix==populations_names[i])
     }
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
    barplot(t(as.matrix(Q[ordered_index_groups,populations_i])),col=rainbow(5),
        space=0,border=NA,  ylab="Ancestry proportions",xaxt="n")
    axis(1,at = pos,labels =populations_names, lwd = 0,las=1)
    #par(xpd=TRUE)
    #legend(1, 1.09,
     #   legend=populations_names,
     #   fill=rainbow(5),lwd=0, cex = 1, xpd = TRUE, ncol = 5)
    dev.off()
}

plot_structure(Q_global,fam_global,"structure_global_Q.png"  )
plot_structure(P_global,fam_global,"structure_global_P.png"  )
plot_structure(Q_pairwise, fam_pairwise, "structure_pairwise_Q.png")
plot_structure(P_pairwise, fam_pairwise, "structure_pairwise_P.png")
```
## 0.6
```bash
outdir="analysis/structure_0.6"
mkdir -p ${outdir}
plink_merged="analysis/qc_merged_data/merged_final_5FID_Fst"
bim_global="analysis/AFR/afr_global_snps.bim"
bim_pairwise="analysis/AFR/afr_pairwise_snp.bim"
cut -f2  ${bim_global} > ${outdir}/snps.global
cut -f2  ${bim_pairwise} > ${outdir}/snps.pairwise

#extract
# from merged that passed QC
# 2773 samples (102 females, 136 males, 2535 ambiguous; 2773 founders) loaded

plink2 --bfile ${plink_merged} \
--extract ${outdir}/snps.global \
--make-bed  \
--out ${outdir}/snps.global


plink2 --bfile ${plink_merged} \
--extract ${outdir}/snps.pairwise \
--make-bed  \
--out ${outdir}/snps.pairwise

#PCA
plink2 --bfile ${outdir}/snps.global --pca --out ${outdir}/pca_global
plink2 --bfile ${outdir}/snps.pairwise --pca --out ${outdir}/pca_pairwise
```
#Fst 0.6

```bash
outdir="analysis/structure_0.6"
mkdir -p ${outdir}
plink_merged="analysis/qc_merged_data/merged_final_5FID_Fst"
bim_global="analysis/AFR_0.6/afr_global_snps_0.6.bim"
bim_pairwise="analysis/AFR_0.6/afr_pairwise_snps_0.6.bim"
cut -f2  ${bim_global} > ${outdir}/snps.global
cut -f2  ${bim_pairwise} > ${outdir}/snps.pairwise

#extract
# from merged that passed QC
# 2773 samples (102 females, 136 males, 2535 ambiguous; 2773 founders) loaded

plink2 --bfile ${plink_merged} \
--extract ${outdir}/snps.global \
--make-bed  \
--out ${outdir}/snps.global_0.6 


plink2 --bfile ${plink_merged} \
--extract ${outdir}/snps.pairwise \
--make-bed  \
--out ${outdir}/snps.pairwise_0.6 

#PCA
plink2 --bfile ${outdir}/snps.global_0.6  --pca --out ${outdir}/pca_global
plink2 --bfile ${outdir}/snps.pairwise_0.6  --pca --out ${outdir}/pca_pairwise
```

2. Plot PCA results

```R
rm(list = ls())
outdir="../structure_0.6/"
eigenvec_global <-read.table(
  paste0(outdir,"pca_global.eigenvec"),header=F)

eigenvec_pairwise <-read.table(
  paste0(outdir,"pca_pairwise.eigenvec"),header=F)

eigenval_global <-read.table(
    paste0(outdir,"pca_global.eigenval"),header=F)
eigenval_pairwise <-read.table(
    paste0(outdir,"pca_pairwise.eigenval"),header=F)


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

plotting_pca(eigenvec_global, eigenval_global,
         "pca_global_0.6")
plotting_pca(eigenvec_pairwise, eigenval_pairwise,
      "pca_pairwise_0.6")
```

###### 3. Population stratifications using Structure

```bash
outdir="analysis/structure_0.6/"
plink_global="${outdir}/snps.global_0.6.bed"
cut -f1 ${outdir}/snps.global_0.6.fam >  ${outdir}/snps.global_0.6.pop
plink_pairwise="${outdir}/snps.pairwise_0.6.bed"
cut -f1 ${outdir}/snps.pairwise_0.6.fam >  ${outdir}/snps.pairwise_0.6.pop
admixture -s 100 ${plink_global}  5 
mv snps.global* ${outdir}/
admixture -s 100 ${plink_pairwise}  5 
mv snps.pairwise* ${outdir}/
```

4. plot population stratifications results

```R
rm(list = ls())
outdir="analysis/structure_0.6/"
#plink_africans="${outdir}/snps.africans.bed"
#plink_overlapped="${outdir}/snps.overlapped.bed"
Q_global=read.table(paste0(outdir,"snps.global_0.6.5.Q")
     ,header=F)
P_global=read.table(paste0(outdir,"snps.global_0.6.5.P")
     ,header=F)
fam_global=read.table(paste0(outdir,"snps.global_0.6.fam")
    ,header=F)

Q_pairwise=read.table(paste0(outdir,"snps.pairwise_0.6.5.Q")
     ,header=F)
P_pairwise=read.table(paste0(outdir,"snps.pairwise_0.6.5.P")
     ,header=F)
fam_pairwise=read.table(paste0(outdir,"snps.pairwise_0.6.fam")
    ,header=F)


plot_structure(Q_global,fam_global,"structure_global_0.6.png"  )
plot_structure(P_global,fam_global,"structure_global_P_0.6.png"  )
plot_structure(Q_pairwise, fam_pairwise, "structure_pairwise_Q_0.6.png")
plot_structure(P_pairwise, fam_pairwise, "structure_pairwise_P_0.6.png")
```
