BiocManager::install("karyoploteR")
library(karyoploteR)
#library(karyoploteR)
#library(GenomicRanges)
#library(IRanges)

base_dir="analysis/AFR/"
outdir="analysis/qc_ancestry_snps/"
overlapped_vcf=paste0(base_dir, "africans_snps_overlapped.vcf")
africans_vcf=paste0(base_dir, "africans_snps_specific.vcf")


overlapped_vcf_df=read.table(overlapped_vcf, comment.char = '#')
overlapped_chr=paste0("chr",overlapped_vcf_df[,1] )
overlapped_snps_pos=overlapped_vcf_df[,2]
overlapped_snps_data= paste0(overlapped_chr, ':', overlapped_snps_pos, '-', overlapped_snps_pos) 
overlapped_snps_gr=GRanges(overlapped_snps_data)


africans_vcf_df=read.table(africans_vcf, comment.char = '#')
africans_chr=paste0("chr",africans_vcf_df[,1] )
africans_snps_pos=africans_vcf_df[,2]
africans_snps_data= paste0(africans_chr, ':', africans_snps_pos, '-', africans_snps_pos) 
africans_snps_gr=GRanges(africans_snps_data)

png(paste0(outdir, "overllaped_snps_karyoplot.png"), width = 465, height = 225, units='mm', res = 300)
kp_overlapped <- plotKaryotype(genome = "hg19", chromosomes="autosomal")
kpAddBaseNumbers(kp_overlapped)
kpPoints(kp_overlapped, data=overlapped_snps_gr, y=0)
dev.off()

png(paste0(outdir, "africans_snps_karyoplot.png"), width = 465, height = 225, units='mm', res = 300)
kp_africans <- plotKaryotype(genome = "hg19", chromosomes="autosomal")
kpAddBaseNumbers(kp_africans)
kpPoints(kp_africans, data=africans_snps_gr, y=0)
dev.off()