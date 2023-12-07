> Directory layout

    projectDir .
               ├── 1KG
               ├── HapMap3
               ├── SGDP
               ├── miscellaneous
               ├── annovar
               ├── preprocess_raw_data      ├── qc_merged_data # Here
               ├── merged_data              ├── Fst
               ├── analysis ────────────────├── AFR
                                            ├── snpEff
                                            ├── Annovar
                                            ├── qc_ancestry_snps
                                            ├── structure

##### Quality control assessments on merged data set
1. [Check individuals' missingness](check_individuals_missingness.md)
2. [Remove missingness](remove_missingness.md)
3. [MAF distribution](MAF_distribution.md)
4. [HWE distribution](HWE_distribution.md)
5. [LD distribution]()