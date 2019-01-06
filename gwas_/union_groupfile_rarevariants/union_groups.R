#!/bin/Rscript

#Date 01/06/2019
#Sanjeev Sariya
#Get union between group files of two batches
#


library(argparse)
library(dplyr)
#library(plyr)

group_file1<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/rare_variants_allbatches/pedigree_reconstruction/genebased_rarevariants/filter_based_INFO_butnomaf/batch1/find_snp_groups_annotated/batch1_CHR22groupfile.txt"
group_file2<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/rare_variants_allbatches/pedigree_reconstruction/genebased_rarevariants/filter_based_INFO_butnomaf/batch2/find_snp_groups_annotated/batch2_CHR22groupfile.txt"


df.group1<-read.table(group_file1,header=FALSE)
df.group2<-read.table(group_file2,header=FALSE)

print(dim(df.group1))
print(dim(df.group2))

lollol<-bind_rows(df.group2,df.group1)
#https://stackoverflow.com/a/36868001/2740831
#
print("binding complete")


print(paste("Duplicate rows are ",sum(duplicated(lollol))))
#https://stackoverflow.com/questions/15589601/print-string-and-variable-contents-on-the-same-line-in-r
kiki<-lollol[!duplicated(lollol),]

#https://stats.stackexchange.com/questions/6759/removing-duplicated-rows-data-frame-in-r
#
print(dim(kiki))

 