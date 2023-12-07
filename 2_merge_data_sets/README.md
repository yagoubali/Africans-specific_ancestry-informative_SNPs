> Directory layout

    projectDir .
               ├── 1KG
               ├── HapMap3
               ├── SGDP
               ├── miscellaneous
               ├── annovar
               ├── preprocess_raw_data
               ├── merged_data   ## Merge data sets here
               ├── analysis

#### Merging cleaned data sets:

1. [Merge 1KG updated FID with SGDP cleaned data set.](Merge_1KG_with_SGDP.md)
2. [Merge cleaned HapMap3 with 1KG and SGDP.](Merge_1KG-SGDP_with_HapMap.md)
3. [Prepare samples panels.](samples_panels.md)
4. [Update fam file of the merged data set.](update_fam.md)
5. [Remove 2 samples of SGDP that are missed in the sample panel file.](clean_2_SGDP_samples.md)
6. [Update FID in the cleaned fam file by using super populations IDs.](update_FID.md)
7. [Remove non 1KG super populations](keep_1KG-superpopulations.md)
