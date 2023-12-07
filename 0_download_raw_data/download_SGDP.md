##### 1.  Simon Genome Diversity Project

```bash

set -x;

download_dir="SGDP" # a  directory to store Simon Genome Diversity Project files


wget -c https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.bim.zip   -P ${download_dir}
wget -c https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.bed  -P ${download_dir}
wget -c https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.fam  -P ${download_dir}

wget -c https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt  -P ${download_dir}
wget -c https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.samples.txt -P ${download_dir}

## unzip
cd ${download_dir}
unzip cteam_extended.v4.maf0.1perc.bim.zip
cd ..
```
