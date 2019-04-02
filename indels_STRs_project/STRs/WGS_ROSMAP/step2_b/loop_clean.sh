#!/bin/bash

#
#Date 04/01/2019
#Sanjeev Sariya

#VCFtools - 0.1.15

####path_str="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/merged_STRs_hipSTR/"

for i in {1..9}
do

#--take in file from filtered VCF - hipstr
str_merge=$path_str/filtered_CHR"$1".vcf
printf "$str_merge\n"

###filename=$path_str/monomorphic_sites_CHR"$i"

qsub -cwd -N clean_CHR"$i" -v bash perform_clean.sh "$i"

##qsub -cwd -N clean_CHR"$i" -b y /mnt/mfs/cluster/bin/bin/vcftools --vcf $str_merge --exclude $filename --recode --recode-INFO-all --out filter.str.no_monomorphic_CHR"$i".vcf


done
