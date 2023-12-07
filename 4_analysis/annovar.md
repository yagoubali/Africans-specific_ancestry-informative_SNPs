##### Functional annotation using Annovar

###### 1. Download Annovar and its databases
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
```
##### 2. Functional annotation using Annovar tool
```bash
annovar_path="annovar";
outdir="analysis/annovar"

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

###### 3. Summarize annotation predictions

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
