#!/bin/bash

for i in {1..22}
do

dir="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/"
cd $dir/CHR"$i"

haps_file=($(ls *.haps))

printf "$haps_file\n"

cut -d" " -f 6- $haps_file | sed 's/\s//g' > inputhaplotype

cd $dir/
pwd
printf "$i complete\n"

done


