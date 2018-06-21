#!/bin/Rscript

##
#Date 05/18/2018
#
#MAF box plot for KGP and HRC
#
#Read KGP data first. 
#
library(dplyr)
library("data.table")
library("Cairo")
library(plyr)
library(ggplot2)

setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/HRC_WellcomeTrust/imputation_1000indiv/impute2/1KGP/unharmonized/merged_info_files/")

filename<-"merged_info_CHR"
ext<-".info"
df_all <- as.list(rep("", 22)) 

for(i in seq(1,22)){
file_impute2_info<-paste(filename,i,ext,sep="")

newdf<-paste("df",i,sep="")

df_all[[i]]<-data.table::fread(file_impute2_info, showProgress = TRUE)
#https://stackoverflow.com/questions/16566799/change-variable-name-in-for-loop-using-r

###get length of a0 allele 
a0_length<-unlist(lapply(df_all[[i]]$a0, function(x) nchar(x)))

###throw whatver is structual variants
#df_all[[i]]<-df_all[[i]][which(a0_length==1),]

###get length of a1 allele 
a1_length<-unlist(lapply(df_all[[i]]$a1, function(x) nchar(x)))

###throw whatver was an structual variants
#df_all[[i]]<-df_all[[i]][which(a1_length==1),]

df_all[[i]]<-df_all[[i]][which(a1_length==1 &a0_length==1),]

df_all[[i]]$info<-as.numeric(as.character(df_all[[i]]$info))

#--work with polymorphic SNPs
df_all[[i]]$exp_freq_a1<-as.numeric(as.character(df_all[[i]]$exp_freq_a1))
#df_all[[i]]<-df_all[[i]][which(df_all[[i]]$exp_freq_a1!=0),]

df_all[[i]]<-df_all[[i]][,-c(9,10,11,12)] #delete type, certainty and other useless

###Fix alelel frequency
###
indices<-df_all[[i]]$exp_freq_a1>0.5
df_all[[i]]$exp_freq_a1[indices]<-(1-df_all[[i]]$exp_freq_a1)[indices]
}

df_merged2 <- data.table::rbindlist(df_all, idcol = TRUE)
df_merged2$panel<-"1000G"
print("we are done reading info files from 1000G")
print(dim(df_merged2))
#
#Got to HRC data now
#

setwd("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/HRC_WellcomeTrust/imputation_1000indiv/impute2/HRC/merged_data/merged_info_files/")

filename<-"merged_info_CHR"
ext<-".info"
df_all <- as.list(rep("", 22)) 

for(i in seq(1,22)){
file_impute2_info<-paste(filename,i,ext,sep="")

newdf<-paste("df",i,sep="")

df_all[[i]]<-data.table::fread(file_impute2_info, showProgress = TRUE)
#https://stackoverflow.com/questions/16566799/change-variable-name-in-for-loop-using-r

###get length of a0 allele 
a0_length<-unlist(lapply(df_all[[i]]$a0, function(x) nchar(x)))

###throw whatver is structual variants
df_all[[i]]<-df_all[[i]][which(a0_length==1),]

###get length of a1 allele 
a1_length<-unlist(lapply(df_all[[i]]$a1, function(x) nchar(x)))

###throw whatver was an structual variants
df_all[[i]]<-df_all[[i]][which(a1_length==1),]
df_all[[i]]$info<-as.numeric(as.character(df_all[[i]]$info))

#--work with polymorphic SNPs
df_all[[i]]$exp_freq_a1<-as.numeric(as.character(df_all[[i]]$exp_freq_a1))
#df_all[[i]]<-df_all[[i]][which(df_all[[i]]$exp_freq_a1!=0),]

df_all[[i]]<-df_all[[i]][,-c(9,10,11,12)] #delete type, certainty and other useless

###Fix alelel frequency
###
indices<-df_all[[i]]$exp_freq_a1>0.5
df_all[[i]]$exp_freq_a1[indices]<-(1-df_all[[i]]$exp_freq_a1)[indices]
}

df_merged <- data.table::rbindlist(df_all, idcol = TRUE)
df_merged$panel<-"HRC"

print("we are done reading info files from HRC")
print(dim(df_merged))

####################################
#Merge data frames
####################################

setwd("/mnt/mfs/scratch/GT_Admix/ss5505/imputation/batch2/plots/impute2/get_median_range_infoThreshold/")
merged_dfs<-rbind(df_merged2,df_merged)
merged_dfs<-merged_dfs[,-c(1,2,3,4,5,6,9)]
print("done with merging")

info_uncommon<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.05 & merged_dfs$info>=0.40),]
info_rare<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.01  & merged_dfs$info>=0.40),]
info_ultra_rare<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.001  & merged_dfs$info>=0.40),]

write.table(as.data.frame(info_uncommon %>% group_by(panel) %>% summarize_at(vars(info),median)),"infocutoff_uncommon",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_rare %>% group_by(panel) %>% summarize_at(vars(info),median)),"infocutoff_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_ultra_rare %>% group_by(panel) %>% summarize_at(vars(info),median)),"infocutoff_ultra_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

print("Done printing higher level MAF info and values")

print("Filtering based on INfo 0.40")
merged_dfs<-merged_dfs[which(merged_dfs$info>=0.40),]

binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 1, by=0.05)))
binned_snps<-binned_snps[complete.cases(binned_snps),]
median_all<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_all,"boxPlots_all",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
print("done with All MAF")
##################################
binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1, include.lowest = TRUE,seq(0, 0.05, by=0.005)))
binned_snps<-binned_snps[complete.cases(binned_snps),]
median_uncommon<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_uncommon,"boxPlots_uncommon",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
print("done with 5% MAF")
##################################
binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 0.01, by=0.001)))
binned_snps<-binned_snps[complete.cases(binned_snps),]
median_ultra<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_ultra,"boxPlots_ultra",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
print("Done with ultra rare")

print("Leaving script")

##################################
##################################
##################################