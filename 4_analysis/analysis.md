## 1. Fst test

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

Filter Fst results

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

## 2. Extract African SNPs

```bash
mkdir -p analysis/AFR
outdir="analysis/AFR"
plink_file="analysis/qc_merged_data/merged_final_5FID_Fst"
afr_snps_overlapped="analysis/Fst/SNPs_list/AFR_0.4.txt";
afr_snps_specific="analysis/Fst/SNPs_list/AFR_specific_0.4.txt";


touch ${outdir}/FID.Afr
echo "AFR" > ${outdir}/FID.Afr

plink --bfile ${plink_file}  \
      --extract ${afr_snps_overlapped}  \
       --keep-fam ${outdir}/FID.Afr  \
      --make-bed   \
       --recode vcf \
      --out ${outdir}/africans_snps_overlapped \
      --allow-no-sex
#10753 variants and 747 people pass filters and QC

plink --bfile ${plink_file}  \
      --extract ${afr_snps_specific}  \
       --keep-fam ${outdir}/FID.Afr  \
      --make-bed   \
       --recode vcf \
      --out ${outdir}/africans_snps_specific \
      --allow-no-sex
#5784 variants and 747 people pass filters and QC.
```

## 2. Characterizing ancestry informative SNPs

### 1. MAF, LD, and HWE distribution

```bash
outdir="analysis/qc_ancestry_snps"
plink_file_africans="analysis/AFR/africans_snps_specific"
plink_file_overlapped="analysis/AFR/africans_snps_overlapped"


#MAF
plink --bfile ${plink_file_africans} --freq --out ${outdir}/MAF_check_africans --allow-no-sex --make-founders

plink --bfile ${plink_file_overlapped} --freq --out ${outdir}/MAF_check_overlapped --allow-no-sex --make-founders
#HWE
plink --bfile ${plink_file_africans}  --hardy --out ${outdir}/HWE_check_africans --allow-no-sex --make-founders

plink --bfile ${plink_file_overlapped}  --hardy --out ${outdir}/HWE_check_overlapped --allow-no-sex --make-founders

#LD
plink --bfile ${plink_file_africans} --r2  --out ${outdir}/LD_africans

plink --bfile ${plink_file_overlapped} --r2  --out ${outdir}/LD_overlapped

```

Plot results

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

## 3. Functional annotation using snpEff

```bash

# Download and install SnpEff
curl -v -L 'https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip' > snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
## download databases
java -jar snpEff.jar download -v GRCh37.87
java -jar snpEff.jar download -v GRCh37.p13
#dbNSFP_4.1a
wget -c https://snpeff.blob.core.windows.net/databases/dbs/GRCh37/dbNSFP_4.1a/dbNSFP4.1a.txt.gz
wget -c https://snpeff.blob.core.windows.net/databases/dbs/GRCh37/dbNSFP_4.1a/dbNSFP4.1a.txt.gz.tbi
# GWAS catlog
# link to download is here --> http://pcingola.github.io/SnpEff/snpsift/gwascatalog/
snpEff_bin="analysis/snpEff"
base_outdir="analysis/snpEff"
mkdir -p  "${base_outdir}/overlapped_snps"
mkdir -p  "${base_outdir}/Africans_snps"
vcf_file_africans="analysis/AFR/africans_snps_specific.vcf"
vcf_file_overlapped="analysis/AFR/africans_snps_overlapped.vcf"
## annotations

java -Xmx8g -jar  ${snpEff_bin}/snpEff.jar  GRCh37.87 ${vcf_file_overlapped} > \
 ${base_outdir}/overlapped_snps/snps_GRCh37.87_ann.vcf

 rename 's/snpEff/GRCh37.87/g' snpEff*
 mv GRCh37* ${base_outdir}/overlapped_snps/

java -Xmx8g -jar ${snpEff_bin}/snpEff.jar  GRCh37.87 ${vcf_file_africans} > \
 ${base_outdir}/Africans_snps/snps_GRCh37.87_ann.vcf

 rename 's/snpEff/GRCh37.87/g' snpEff*
 mv GRCh37* ${base_outdir}/Africans_snps/
 # GRCh37.p13

java -Xmx8g -jar  ${snpEff_bin}/snpEff.jar  GRCh37.p13 ${vcf_file_overlapped} > \
 ${base_outdir}/overlapped_snps/snps_GRCh37.p13_ann.vcf

 rename 's/snpEff/GRCh37.p13/g' snpEff*
 mv GRCh37* ${base_outdir}/overlapped_snps/

java -Xmx8g -jar  ${snpEff_bin}/snpEff.jar  GRCh37.p13 ${vcf_file_africans} >  \
${base_outdir}/Africans_snps/snps_GRCh37.p13_ann.vcf

rename 's/snpEff/GRCh37.p13/g' snpEff*
 mv GRCh37* ${base_outdir}/Africans_snps/

## GWAS catalog
mkdir -p  "${base_outdir}/catalog"

vcf_file_africans="analysis/AFR/africans_snps_specific.vcf"
vcf_file_overlapped="analysis/AFR/africans_snps_overlapped.vcf"
cp ${vcf_file_africans}  ${base_outdir}/catalog/
cp ${vcf_file_overlapped} ${base_outdir}/catalog/
gwas_catalog="gwas_catalog_v1.0.2-associations_e110_r2023-11-08.tsv"

java -Xmx8g -jar  ${snpEff_bin}/SnpSift.jar gwasCat \
-db ${snpEff_bin}/data/${gwas_catalog}  \
${base_outdir}/catalog/africans_snps_specific.vcf | \
tee ${base_outdir}/catalog/africans.gwas.vcf

java -Xmx8g -jar  ${snpEff_bin}/SnpSift.jar gwasCat \
 -db ${snpEff_bin}/data/${gwas_catalog}  \
${base_outdir}/catalog/africans_snps_overlapped.vcf | \
tee ${base_outdir}/catalog/overlapped.gwas.vcf

grep -v '#' ${base_outdir}/catalog/africans.gwas.vcf | grep 'GWASCAT' | cut -f 1,2,3,4,5,8 > ${base_outdir}/catalog/africans.gwas

grep -v '#' ${base_outdir}/catalog/overlapped.gwas.vcf | grep 'GWASCAT' | cut -f 1,2,3,4,5,8 > ${base_outdir}/catalog/overlapped.gwas
```

## 5. Functional annotation using Annovar

```bash
## To download annovar https://annovar.openbioinformatics.org/en/latest/#reference
annovar_path="annovar";
outdir="analysis/annovar"
mkdir -p ${outdir}

## Download db

 perl ${annovar_path}/annotate_variation.pl -buildver hg19 -downdb \
  -webfrom annovar refGene \
   ${annovar_path}/humandb/

perl ${annovar_path}/annotate_variation.pl -buildver hg19 -downdb \
  -webfrom annovar knownGene \
   ${annovar_path}/humandb/

 perl ${annovar_path}/annotate_variation.pl -buildver hg19 -downdb \
  -webfrom annovar ensGene \
   ${annovar_path}/humandb/
 perl ${annovar_path}/annotate_variation.pl -buildver hg19 -downdb \
  -webfrom annovar clinvar_20221231 \
  ${annovar_path}/humandb/

 perl ${annovar_path}/annotate_variation.pl -buildver hg19 -downdb \
  -webfrom annovar dbnsfp42c \
  ${annovar_path}/humandb/

  perl ${annovar_path}/annotate_variation.pl  -build hg19 -downdb \
    targetScanS \
   ${annovar_path}/humandb/

 ## annotations

snps_overlapped="analysis/AFR/africans_snps_overlapped.vcf"
snps_specific="analysis/AFR/africans_snps_specific.vcf"

perl ${annovar_path}/table_annovar.pl \
-vcfinput ${snps_overlapped} ${annovar_path}/humandb/ -buildver hg19 \
-out ${outdir}/overlapped.annovar \
-remove \
-protocol refGene,knownGene,ensGene,dbnsfp42c,clinvar_20221231,targetScanS \
-operation g,g,g,f,f,r -nastring .  -polish

perl ${annovar_path}/table_annovar.pl \
-vcfinput ${snps_specific} ${annovar_path}/humandb/ -buildver hg19 \
-out ${outdir}/africans.annovar \
-remove \
-protocol refGene,knownGene,ensGene,dbnsfp42c,clinvar_20221231,targetScanS \
-operation g,g,g,f,f,r -nastring .  -polish
```

Summarize annotation predictions

```R
rm(list = ls())


outdir="analysis/annovar/summary/"
#system("mkdir -p outdir")
snps_overlapped="analysis/annovar/overlapped.annovar.hg19_multianno.txt"
snps_africans="analysis/annovar/africans.annovar.hg19_multianno.txt"
snps_overlapped_db=read.table(snps_overlapped, header=T, sep="\t")
snps_africans_db=read.table(snps_africans, header=T, sep="\t")

annovar_prediction= function (annovar_results, out_name){
       all_db= data.frame()
       all_headers=(as.character(colnames(annovar_results))[21:116])
       pred_headers=c()
       #for( i in 1:length(all_headers)){
       #  cat(all_headers[i], " Unique values ..  ", length(unique(annovar_results[,which(colnames(annovar_results) == all_headers[i])])), "\n")
       #  }

       for( i in 1:length(all_headers)){
           if(length(grep("pred", all_headers[i]))){
                pred_headers=c(pred_headers,  all_headers[i])
            }
      }

      for( i in 1:length(pred_headers)){
        results= as.data.frame(table(annovar_results[,which(colnames(annovar_results) == pred_headers[i])]))
        percent=round((results$Freq/sum(results$Freq))*100,2)
        results= cbind(results, as.data.frame(percent))
        outfile= paste0(outdir, pred_headers[i], out_name, ".txt")
        write.table(results, file=outfile, row.names=F)
        annotation=as.data.frame(rep(pred_headers[i], dim(results)[1]))
        results=cbind(annotation,results)
        all_db=rbind(all_db, results)

  }


  return(all_db)

}
all_overlapped=paste0(outdir, "all_overlapped.txt")
all_africans=paste0(outdir, "all_africans.txt")
overlapped=annovar_prediction(snps_overlapped_db,"_overlapped" )
africans=annovar_prediction(snps_africans_db,"_africans" )
write.table(overlapped, file=all_overlapped, row.names=F)
write.table(africans, file=all_africans, row.names=F)
```

## Population structure analaysis

### 1. PCA

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

### plot the results

````R
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

### 2. Structure
```bash
outdir="analysis/structure/"
plink_africans="${outdir}/snps.africans.bed"
plink_overlapped="${outdir}/snps.overlapped.bed"
admixture -s 100 ${plink_africans}  5
mv snps.africans* ${outdir}/
admixture -s 100 ${plink_overlapped}  5
mv snps.overlapped* ${outdir}/
````

### plot results

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

<!--## 4. maftools
```bash
vcf2maf="miscellaneous/mskcc-vcf2maf-754d68a/vcf2maf.pl"
mkdir -p analysis/maftools
vcf_overlapped="analysis/AFR/africans_snps_overlapped.vcf"
vcf_specific="analysis/AFR/africans_snps_specific.vcf"
vep_cache="miscellaneous/ensembl-vep/cache"
ref="${vep_cache}/homo_sapiens/110_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz"
vep_path="miscellaneous/ensembl-vep/"
outdir="analysis/maftools"
#https://gist.github.com/ckandoth/61c65ba96b011f286220fa4832ad2bc0?permalink_comment_id=3827238
    perl ${vcf2maf} --input-vcf ${vcf_overlapped} \
    --vep-path   ${vep_path} \
    --ref-fasta  ${ref} \
    --vep-data  ${vep_cache} \
    --output-maf ${outdir}/snps_overlapped.vep.maf

perl ${vcf2maf} --input-vcf ${vcf_specific} --output-maf ${outdir}/snps_specific.vep.maf
``` -->
