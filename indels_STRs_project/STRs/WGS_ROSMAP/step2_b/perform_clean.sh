#!/bin/bash

#04/01/2019
#

printf "start CHR"$1"\n"

path_str="/mnt/mfs/hgrcgrid/shared/GT_ADMIX/INDEL_comparisons/ROSMAP_WGS/merged_STRs_hipSTR/"

str_merge=$path_str/filtered_CHR"$1".vcf
printf "$str_merge\n"

filename=$path_str/monomorphic_sites_CHR"$1"

/mnt/mfs/cluster/bin/bin/vcftools --vcf $str_merge --exclude $filename \
--recode \
--recode-INFO-all  --out filter.str.no_monomorphic_CHR"$1"

printf "exluded monomorphic\n"

meandepth=$path_str/"persitemeandepth_CHR"$1".ldepth.mean" ##file containing depth per site

awk '{
if(FNR>1){
if($3<10){
##keep anything than more than or equal 10

print $1"\t"$2
}
}}' $meandepth > depthless_thanten.CHR"$1" ##throw these sites first column is CHR and other position

/mnt/mfs/cluster/bin/bin/vcftools --vcf filter.str.no_monomorphic_CHR"$1".recode.vcf --exclude-positions ./depthless_thanten.CHR"$1" --recode \
--recode-INFO-all  --out filter.str.no_monomorphic.poor_depth_CHR"$1"

printf "$meandepth excluded poor depth\n"

file_missing=$path_str/persite_missing_CHR"$1".lmiss ##this is file that contains missingness per site

awk '{
if(FNR>1){
if($6>=0.2){

#exclude anything that has missingness more than or equal to 20

print $1"\t"$2
}
}}' $file_missing > missingmore_than.twnty.CHR"$1" ## ##throw these sites first column is CHR and other position

/mnt/mfs/cluster/bin/bin/vcftools --vcf filter.str.no_monomorphic.poor_depth_CHR"$1".recode.vcf \
--exclude-positions missingmore_than.twnty.CHR"$1" --recode \
--recode-INFO-all  --out filter.str.no_monomorphic.poor_depth.highmissing_CHR"$1" 

printf "$file_missing excluded missing \n"

/mnt/mfs/cluster/bin/bin/vcftools  --vcf  filter.str.no_monomorphic.poor_depth.highmissing_CHR"$1".recode.vcf --missing-site --out CHR"$1"_postfiltering

printf "Created post filtered missingness file\n"

awk '{
if(FNR>1){
if($6>=0.2){

#exclude anything that has missingness more than or equal to 20
print $1"\t"$2
}
}}' CHR"$1"_postfiltering.lmiss > exclude.missing_CHR"$1"


/mnt/mfs/cluster/bin/bin/vcftools --vcf  filter.str.no_monomorphic.poor_depth.highmissing_CHR"$1".recode.vcf \
--exclude-positions exclude.missing_CHR"$1" --recode --recode-INFO-all  --out filter.str.no_monomorphic.poor_depth.highmissing.postcleaning_CHR"$1" 

printf "Completed post filtered SNPs cleaning\n"
printf "Exiting CHR"$1"\n"
