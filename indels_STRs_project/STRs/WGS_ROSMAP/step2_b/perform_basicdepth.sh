#!/bin/bash

#
#Date 03-09-2019
#


cd /mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/merged_STRs_hipSTR

file_vcf=merged_ROSmap_hipSTRCHR"$1".vcf.gz
printf "begin $file_vcf\n"

##get locus missingness
vcftools --gzvcf $file_vcf  --missing-site  -c > persite_missing_CHR"$1".lmiss

printf "missingness per locus complete \n"
##generate depth per individual
vcftools --gzvcf $file_vcf   --depth -c > perperson_depth_CHR"$1".idepth
printf "depth per person complete \n"

## get site depth
vcftools --gzvcf $file_vcf --site-depth -c > persitedepth_CHR"$1".ldepth
printf "depth per locus complete \n"

##get per site mean depth
vcftools --gzvcf $file_vcf --site-mean-depth -c > persitemeandepth_CHR"$1".ldepth.mean
printf "mean depth per locus complete \n"

##Generates a file reporting the missingness on a per-individual basis. The file has the suffix ".imiss".
vcftools --gzvcf $file_vcf --missing-indv -c > perindividual_missingness_CHR"$1".imiss
printf "CHR$1 completed\n"
