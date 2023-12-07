##### Prepare samples panels

```bash
outdir="merged_data"
panel_1KG="1KG/integrated_call_samples_v3.20130502.ALL.panel"
panel_SGDP="SGDP/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt"
panel_HapMap="HapMap3/relationships_w_pops_121708.txt"
fam_SGDP_final="preprocess_raw_data/SGDP/SGDP_final.fam"
fam_HapMap_final="preprocess_raw_data/HapMap3/hapmap3_final.fam"

cp ${panel_1KG} ${outdir}/
cut -f4,6,7,8,11  ${panel_SGDP} > ${outdir}/SGDP.panel1
sed -i 's/(//g' ${outdir}/SGDP.panel1
sed -i 's/)//g' ${outdir}/SGDP.panel1
cat  ${outdir}/SGDP.panel1 | sed -e '1d' |  \
awk -f  miscellaneous/SGDP_panel.awk - > ${outdir}/SGDP_ALL.panel

sed -i 's/-/_/g' ${outdir}/SGDP_ALL.panel
cp ${fam_SGDP_final} ${outdir}/
sed -i 's/-/_/g' ${outdir}/SGDP_final.fam

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1]=$1; next} \
($1 in a ) {print $0}' \
${outdir}/SGDP_final.fam ${outdir}/SGDP_ALL.panel > ${outdir}/SGDP_panel.merge

## HapMap3

awk -f  miscellaneous/HapMap3_panel.awk  ${panel_HapMap} > ${outdir}/HapMap_ALL.panel
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1]=$1; next} \
($1 in a ) {print $0}' \
 ${fam_HapMap_final} ${outdir}/HapMap_ALL.panel > ${outdir}/HapMap_panel.merge

##
 cat ${panel_1KG} >  ${outdir}/Final.panel
 cat  ${outdir}/SGDP_panel.merge >>  ${outdir}/Final.panel
 cat ${outdir}/HapMap_panel.merge >>  ${outdir}/Final.panel
```