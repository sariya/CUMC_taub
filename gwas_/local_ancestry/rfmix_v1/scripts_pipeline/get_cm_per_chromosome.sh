	#!/bin/bash

######################################
#Date 12/23/2018
######################################

for i in {1..22}
do

dir="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/"

cd $dir/CHR"$i"
pwd

set -x

awk '{print $4}' hgwas.hgdp.ceuyri_CHR"$i".bim | xargs -I % grep -w ^%  genetic_map_chr"$i"_combined_b37.txt | awk '{print $3}' > centimorgans
wc -l hgwas.hgdp.ceuyri_CHR"$i".bim 

wc -l centimorgans

cd $dir/
pwd

printf "CHR$i complete\n"

done


