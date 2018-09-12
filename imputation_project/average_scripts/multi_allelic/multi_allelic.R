#!/bin/Rscript

#Date 06/20/2018
#
#get counts of input data. 1000G is what this is used for mainly. Find multiallelic snps .
#Get their counts. Get their counts if INFO >=0.40 and info 0.80
#Get counts at each step
#

setwd("/mnt/mfs/scratch/GT_Admix/ss5505/imputation/batch2/counts_data/info_threshold_counts/")

library("data.table")
library(plyr)
library(ggplot2)
rm(list=ls())

setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Imputation_1KGP_HRC/imputation_1000indiv/impute2/1KGP/merged_files/unharmonized/info/")

filename<-"merged_info_CHR"
ext<-".info"
df_all <- as.list(rep("", 22)) 

##loop through different chromosomes
for(i in seq(1,22)){
file_impute2_info<-paste(filename,i,ext,sep="")
df_all[[i]]<-data.table::fread(file_impute2_info, showProgress = TRUE,header=TRUE)
}

#######################################################################################################################################
#--merge list of data frames  https://stackoverflow.com/questions/2851327/convert-a-list-of-data-frames-into-one-data-frame #
#######################################################################################################################################

chr_info_melt2<-data.table::rbindlist(df_all)

print(paste("Total imputed SNPs ",nrow(chr_info_melt2)) )
a0_length<-unlist(lapply(chr_info_melt2$a0, function(x) nchar(x)))
a1_length<-unlist(lapply(chr_info_melt2$a1, function(x) nchar(x)))

multi_allelic<-chr_info_melt2[which( a1_length==1 & a0_length>1 | a1_length>1 & a0_length==1 | a1_length>1 & a0_length!=1 | a1_length!=1 & a0_length>1 ),]
print(paste("Total Multi-allelic SNPs ",nrow(multi_allelic)) )

#--

multi_allelic<-multi_allelic[which(multi_allelic$info>=0.40),] #for 0.40 Info

setwd("/mnt/mfs/scratch/GT_Admix/ss5505/imputation/batch2/counts_data/info_threshold_counts/")
print(paste("Total Bi-allelic with 40% info SNPs after MAF filter for 1000G Multialleic",nrow(multi_allelic)) )

multi_allelic<-multi_allelic[which(multi_allelic$info>=0.80),] # get count for 80% info multi-allelic too

print(paste("Total Bi-allelic with 80% info SNPs after MAF filter for 1000G Multialleic",nrow(multi_allelic)) )
print("we are done cranking data")