##### 1. HapMap3

```bash

set -x;
download_dir="HapMap3/"  ## a  directory to store hapmap 3 files

wget -c https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/00README.txt   -P ${download_dir}
wget -c https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/hapmap3_r1_b36_fwd.qc.poly.tar.bz2  -P ${download_dir}
wget -c https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/phase_3_samples.doc  -P ${download_dir}
wget -c https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/relationships_w_pops_051208.txt  -P ${download_dir}
wget -c https://ftp.ncbi.nlm.nih.gov/hapmap/phase_3/relationships_w_pops_121708.txt  -P ${download_dir}

cd ${download_dir}
tar -xvf hapmap3_r1_b36_fwd.qc.poly.tar.bz2
cd ..
```
