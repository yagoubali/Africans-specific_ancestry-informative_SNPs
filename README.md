## Characterization and functional annotations of Africans-specific ancestry informative SNPs

#### Abstract

Human genetic diversity is critical in determining various phenotypes and genetic traits, including rare diseases. Thus, addressing the source of genetic variation is essential for any genomic-based evidence and its applications. For instance, an individual’s genetic diversity should be considered when using personalized medicine to estimate polygenic risk scores. Two levels of diversification can be regarded as (i) individual-level genetic diversity and (ii) ethnic/community-level genetic diversity. However, many researchers pay more attention to ethnic/community-level genetic diversity. In this work, we estimated African-specific ancestry single-nucleotide polymorphisms (SNPs) using three public data sets: the 1000 genomes project (1KGP), the Simons Genome Diversity Project (SGDP), and the International HapMap Project (HapMap3). The merged data is grouped into 5 super populations as provided in 1KGP, namely Africans (AFR), Americans (AMR), East Asians (EAS), Europeans (EUR), and South Asians (SAS). A total of 5784 variants were reported as African-specific ancestry informative SNPs, while 10753 variants were reported as overlapped ancestry informative among all five populations. The functional annotation of these SNPs revealed that most of these SNPs were not reported in many databases because these SNPs were not associated with protein-coding genes.

#### Tools

1. wget (v. 1.21.2)
2. vcftools (v. 0.1.16)
3. bgzip (bgzip (htslib) 1.13+ds)
4. bcftools (version 1.13 Using htslib 1.13+ds)
5. plink2, plink, i.e, plink 1.9
6. LiftOver
7. sed
8. awk
9. [vcf2maf](https://github.com/mskcc/vcf2maf)
10. [maftools](http://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html)
11. [SnpEff & SnpSif](http://pcingola.github.io/SnpEff/examples/)
12. [vep](http://www.ensembl.org/info/docs/tools/vep/script/index.html)

#### Analysis steps

1. [Download data sets](0_download_raw_data/README.md)
2. [Preprocessing and data cleaning](1_raw_data_preprocessing/README.md)
   1. [1000 genomes](1_raw_data_preprocessing/1KG.md)
   2. [Simon Genome Diversity Project (SGDP)](1_raw_data_preprocessing/SGDP.md)
   3. [HapMap3](1_raw_data_preprocessing/HapMap3.md)
3. [Data sets merging](2_merge_data_sets/MERGE.md)
4. [Quality control check of merged data set.](3_quality_control/qc.md)
   1. Check individuals' missingness and Check SNPs' missingness
   2. Filter out data by removing missingness.
5. [Analysis](4_analysis/analysis.md)
