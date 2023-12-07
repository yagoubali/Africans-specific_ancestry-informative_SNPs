> Directory layout

    projectDir .
               ├── 1KG
               ├── HapMap3
               ├── SGDP
               ├── miscellaneous
               ├── annovar               ├ 1KG     # Preprocess 1KG data
               ├── preprocess_raw_data ──├ SGDP    # Preprocess SGDP data
               ├── merged_data           └ HapMap3 # Preprocess HapMap3 data
               ├── analysis

#### [Preprocessing 1000 genomes data](Preprocessing_1KG.md)
> The preprocessing steps for 1KGP data include:
1. Normalize SNPs by splitting multiallelic sites into biallelic records
2. Extract SNPs from the raw VCF files.
3. Assign and update input VCF files with the correct SNPs rsIDs.
4. Convert VCF files format to PLINK format
5. Merge per chromosome PLINK files
6. Remove ambiguous SNPs.
7. Remove duplicate SNPs

#### [Preprocessing Simons Genome Diversity Project (SGDP)](Preprocessing_SGDP.md)
> SGDP data is preprocessed by performing the following steps:
1. Updating IDs
2. Removing individuals that are existed in 1KGP
3. Removing duplicate SNPs,
4. Removing ambiguous SNPs
5. Assigning SNPs names based on the correct SNPs rsIDs,
6. Removing SNPs that their chromosome doses do not match 1KGP. Also Correcting SNPs genomic positions based on 1KGP
7. Filipping alleles based on 1KGP
8. Remove SNPs that are inconsistent with 1KGP after performing the steps (2-7).

#### [Preprocessing International HapMap Project (HapMap3)](Preprocessing_HapMap3.md)
> HapMap3 is preprocessed by performing the following steps:
1. Convert map and ped files into PLINK binary filesets
2. Merge HapMap files.
3. Removing duplicates and ambiguous SNPs
4. Performing liftOver to convert the SNPs genomic positions from NCBI36 (hg18) into NCBI37 (hg19)
5. Removing individuals that are existed in 1KGP and SGDP
6. Removing SNPs that are not existed in both 1KGP and SGDP,
7. Removing SNPs that their chromosome doses do not match both 1KGP and SGDP.
8. Correcting SNPs genomic positions based on 1KGP and SGDP.
9. Flipping alleles based on 1KGP and SGDP.
10. Remove SNPs that are inconsistent with 1KGP and SGDP after performing the steps (2-9).
