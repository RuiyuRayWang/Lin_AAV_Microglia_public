#!/bin/bash
# 20210609
# awk command used to extract ensembl id and gene name

awk -F "\t" '$0 {print $9}' Mus_musculus.GRCm38.102.chr.mScarlet_WPRE.gtf | awk -F "; " '{for(i=1;i<=NF;i++){ if($i ~ /gene_id/ ){for(j=i+1;j<=NF;j++){ if($j ~ /gene_name/) {print $i"\t"$j} } } } }' | tr -d "\"" | awk -F "\t" '{sub(/gene_id/,"");sub(/gene_name/,"");print $0}' | tr -d " " |sort | uniq > MM.GRCm38.102.annotation.tab
