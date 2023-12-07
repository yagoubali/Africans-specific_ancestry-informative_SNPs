##### 1. Fst test

```bash
mkdir -p analysis/Fst
outdir="analysis/Fst"
plink_file="analysis/qc_merged_data/merged_final_5FID_Fst"

## prepare phenotype file
awk 'BEGIN{OFS="\t";  print "#FID","IID","pop"} {print $1,$2,$1}' \
${plink_file}.fam > ${outdir}/pheno.all

cut -f1 ${plink_file}.fam  > ${outdir}/populations

plink2 --bfile ${plink_file} \
  --fst pop report-variants method=wc \
  --out ${outdir}/Fst_wc \
  --pheno ${outdir}/pheno.all
```

##### 2. Filter Fst results

```R
setwd("analysis/Fst")
ggbg2 <- function() {
  points(0,0,pch=16, cex=1e6, col="lightgray")
  grid(col="white", lty=1)
}

fitter=function(Fst, cutoff){
  data=read.table(Fst, header=F)
  cat(dim(data), "\n")
  headerName=c("CHROM",	"POS"	,"ID"	,"OBS_CT",	"HUDSON_FST")
  colnames(data)=headerName
  Output=substitute(Fst)
  File_output <- gsub('var', 'txt', Output)
  indexes=which(data$HUDSON_FST>=cutoff)
  write.table(data[indexes,],file=substitute(File_output),row.names=FALSE,quote = FALSE)

  Image_output <- gsub('var', 'jpeg', Output)
  head(data)

   jpeg(Image_output)
   par(mar=c(4,4,2,0),oma=c(0.1,1,1,1),mfrow=c(2,1))
   plot(rev(sort(data$HUDSON_FST)),col='#000080', panel.first=ggbg2(),xlab="Index",
    ylab="Fst\  value", pch='.',cex = 0.5) #  main= "Weir and Cockerham Fst Estimator",
    hist(data$HUDSON_FST,col='blue', panel.first=ggbg2(),xlab="Fst\  value",
     ylab="Frequency",main="",cex = 0.5,cex.axis=0.7,las=2)
  dev.off()
}
AFR.AMR="Fst_wc.AFR.AMR.fst.var"
AFR.EAS="Fst_wc.AFR.EAS.fst.var"
AFR.EUR="Fst_wc.AFR.EUR.fst.var"
AFR.SAS="Fst_wc.AFR.SAS.fst.var"
AMR.EAS="Fst_wc.AMR.EAS.fst.var"
AMR.EUR="Fst_wc.AMR.EUR.fst.var"
AMR.SAS="Fst_wc.AMR.SAS.fst.var"
EAS.EUR="Fst_wc.EAS.EUR.fst.var"
EAS.SAS="Fst_wc.EAS.SAS.fst.var"
EUR.SAS="Fst_wc.EUR.SAS.fst.var"

fitter("Fst_wc.AFR.AMR.fst.var",0.4)
fitter("Fst_wc.AFR.EAS.fst.var",0.4)
fitter("Fst_wc.AFR.EUR.fst.var",0.4)
fitter("Fst_wc.AFR.SAS.fst.var",0.4)
fitter("Fst_wc.AMR.EAS.fst.var",0.4)
fitter("Fst_wc.AMR.EUR.fst.var",0.4)
fitter("Fst_wc.AMR.SAS.fst.var",0.4)
fitter("Fst_wc.EAS.EUR.fst.var",0.4)
fitter("Fst_wc.EAS.SAS.fst.var",0.4)
fitter("Fst_wc.EUR.SAS.fst.var",0.4)

AFR.AMR_t=read.table("Fst_wc.AFR.AMR.fst.txt", header=T)
AFR.EAS_t=read.table("Fst_wc.AFR.EAS.fst.txt", header=T)
AFR.EUR_t=read.table("Fst_wc.AFR.EUR.fst.txt", header=T)
AFR.SAS_t=read.table("Fst_wc.AFR.SAS.fst.txt", header=T)

AMR.EAS_t=read.table("Fst_wc.AMR.EAS.fst.txt", header=T)
AMR.EUR_t=read.table("Fst_wc.AMR.EUR.fst.txt", header=T)
AMR.SAS_t=read.table("Fst_wc.AMR.SAS.fst.txt", header=T)

EAS.EUR_t=read.table("Fst_wc.EAS.EUR.fst.txt", header=T)
EAS.SAS_t=read.table("Fst_wc.EAS.SAS.fst.txt", header=T)

EUR.SAS_t=read.table("Fst_wc.EUR.SAS.fst.txt", header=T)

###
AFR.AMR_RsId=as.character(AFR.AMR_t$ID)
AFR.EAS_RsId=as.character(AFR.EAS_t$ID)
AFR.EUR_RsId=as.character(AFR.EUR_t$ID)
AFR.SAS_RsId=as.character(AFR.SAS_t$ID)

AMR.EAS_RsId=as.character(AMR.EAS_t$ID)
AMR.EUR_RsId=as.character(AMR.EUR_t$ID)
AMR.SAS_RsId=as.character(AMR.SAS_t$ID)

EAS.SAS_RsId=as.character(EAS.SAS_t$ID)
EAS.EUR_RsId=as.character(EAS.EUR_t$ID)

EUR.SAS_RsId=as.character(EUR.SAS_t$ID)

#AFR, AMR, EAS, EUR, SAS
#### AFR
AFR <-list(AFR.AMR=AFR.AMR_RsId,
           AFR.EAS= AFR.EAS_RsId,
           AFR.EUR= AFR.EUR_RsId,
           AFR.SAS=AFR.SAS_RsId)

#AFR, AMR, EAS, EUR, SAS
#AMR
AMR <-list(AMR.AFR=AFR.AMR_RsId,
           AMR.EAS=AMR.EAS_RsId,
           AMR.EUR=AMR.EUR_RsId,
           AMR.SAS=AMR.SAS_RsId)
#AFR, AMR, EAS, EUR, SAS
# EAS
EAS <- list(EAS.AFR= AFR.EAS_RsId,
           EAS.AMR=AMR.EAS_RsId,
           EAS.EUR=EAS.EUR_RsId,
           EAS.SAS=EAS.SAS_RsId)

#AFR, AMR, EAS, EUR, SAS
# EUR
EUR <- list(EUR.AFR= AFR.EUR_RsId,
           EUR.AMR=AMR.EUR_RsId,
           EUR.EAS=EAS.EUR_RsId,
           EUR.SAS=EUR.SAS_RsId)

##AFR, AMR, EAS, EUR, SAS
# SAS
SAS <-list(SAS.AFR= AFR.SAS_RsId,
           SAS.AMR=AMR.SAS_RsId,
           SAS.EAS=EAS.SAS_RsId,
           SAS.EUR=EUR.SAS_RsId)

install.packages("UpSetR")
library(UpSetR)  ##---->
system("mkdir UpSetR")
png("UpSetR/AFR.png",width=8, height=8, units="in", res=500)
upset(fromList(AFR),point.size=5,shade.color="darkblue", show.numbers=FALSE)
dev.off()
png("UpSetR/AMR.png",width=8, height=8, units="in", res=500)
upset(fromList(AMR),point.size=5,shade.color="darkblue", show.numbers=FALSE)
dev.off()
png("UpSetR/EAS.png",width=8, height=8, units="in", res=500)
upset(fromList(EAS),point.size=5,shade.color="darkblue", show.numbers=FALSE)
dev.off()
png("UpSetR/EUR.png",width=8, height=8, units="in", res=500)
upset(fromList(EUR),point.size=5,shade.color="darkblue", show.numbers=FALSE)
dev.off()
png("UpSetR/SAS.png",width=8, height=8, units="in", res=500)
upset(fromList(SAS),point.size=5, shade.color="darkblue", show.numbers=FALSE)
dev.off()

#AFR, AMR, EAS, EUR, OCE, SAS, WAS
AFR_set<-Reduce(intersect,AFR)
AFR_specific=Reduce(setdiff, AFR)
AMR_set<-Reduce(intersect,AMR)
EAS_set<-Reduce(intersect,EAS)
EUR_set<-Reduce(intersect,EUR)
SAS_set<-Reduce(intersect,SAS)

AFR_df=as.data.frame(AFR_set)
AFR_specific_df=as.data.frame(AFR_specific)
AMR_df=as.data.frame(AMR_set)
EAS_df=as.data.frame(EAS_set)
EUR_df=as.data.frame(EUR_set)
SAS_df=as.data.frame(SAS_set)

system("mkdir SNPs_list")
write.table(AFR_df, file="SNPs_list/AFR_0.4.txt",row.names=F, col.names=F,quote = FALSE)
write.table(AFR_specific_df, file="SNPs_list/AFR_specific_0.4.txt",row.names=F, col.names=F,quote = FALSE)
write.table(AMR_df, file="SNPs_list/AMR_0.4.txt",row.names=F, col.names=F,quote = FALSE)
write.table(EAS_df, file="SNPs_list/EAS_0.4.txt",row.names=F, col.names=F,quote = FALSE)
write.table(EUR_df, file="SNPs_list/EUR_0.4.txt",row.names=F, col.names=F,quote = FALSE)
write.table(SAS_df, file="SNPs_list/SAS_0.4.txt",row.names=F, col.names=F,quote = FALSE)

 install.packages("venn")
library("venn")
png("AFR_venn.png",width=8, height=8, units="in", res=500)
venn(AFR,zcolor="style",box=F, ilcs = 0.8)
dev.off()

png("AMR_venn.png",width=8, height=8, units="in", res=500)
 venn(AMR,zcolor="style",box=F, ilcs = 0.8)
dev.off()

png("EAS_venn.png",width=8, height=8, units="in", res=500)
 venn(EAS,zcolor="style",box=F, ilcs = 0.8)
dev.off()

png("EUR_venn.png",width=8, height=8, units="in", res=500)
 venn(EUR,zcolor="style",box=F, ilcs = 0.8)
dev.off()

png("SAS_venn.png",width=8, height=8, units="in", res=500)
venn(SAS,zcolor="style",box=F, ilcs = 0.8)
dev.off()
```
