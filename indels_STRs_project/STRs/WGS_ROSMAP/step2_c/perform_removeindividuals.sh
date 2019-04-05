#!/bin/bash

#
#Date 04/02/2019
#Sanjeev Sariya

printf "CHR"$1" begun\n"

file_keep="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/pheno/keep_individualsSTR"

str_vcf="filter.str.no_monomorphic.poor_depth.highmissing.postcleaning_CHR"$1".recode.vcf"

/mnt/mfs/cluster/bin/bin/vcftools \
--vcf $str_vcf \
--keep  $file_keep \
--recode \
--recode-INFO-all \
--out "filter.str.no_monomorphic.poor_depth.highmissing.postcleaning.removeindi_CHR"$1""

printf "Completed CHR"$1"\n"

