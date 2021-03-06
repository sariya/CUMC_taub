#!/bin/Rscript

#Date 06/20/2018
#

#Get counts at each step
#

setwd("/mnt/mfs/scratch/GT_Admix/ss5505/imputation/batch2/counts_data/info_threshold_counts/")

library("data.table")
library(plyr)
library(ggplot2)
rm(list=ls())

setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/HRC_WellcomeTrust/imputation_1000indiv/impute2/1KGP/unharmonized/merged_info_files/")

filename<-"merged_info_CHR"
ext<-".info"

df_all <- as.list(rep("", 22)) 

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

###throw whatver is structual variants
a1_length<-unlist(lapply(chr_info_melt2$a1, function(x) nchar(x)))

###throw whatver was an structual variants

multi_allelic<-chr_info_melt2[which( a1_length==1 & a0_length>1 | a1_length>1 & a0_length==1 | a1_length>1 & a0_length!=1 | a1_length!=1 & a0_length>1 ),]

#chr_info_melt2<-chr_info_melt2[which(a0_length==1),] #chr_info_melt2<-chr_info_melt2[which(a1_length==1),]

print(paste("Total Multi-allelic SNPs ",nrow(multi_allelic)) )

#--

#throw SNPs with MAF<0.40 and INFO ==0 
multi_allelic<-multi_allelic[which(multi_allelic$info>=0.40),]

setwd("/mnt/mfs/scratch/GT_Admix/ss5505/imputation/batch2/counts_data/info_threshold_counts/")

print(paste("Total Bi-allelic SNPs after MAF filter for 1000G Multialleic",nrow(chr_info_melt2)) )
print("we are done cranking data")