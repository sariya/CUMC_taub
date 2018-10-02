#!/bin/Rscript

##
#Date 06/16/2018
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

df_all[[i]]<-df_all[[i]][which(a1_length==1 & a0_length==1),]
df_all[[i]]$info<-as.numeric(as.character(df_all[[i]]$info))

#--work with polymorphic SNPs
df_all[[i]]$exp_freq_a1<-as.numeric(as.character(df_all[[i]]$exp_freq_a1))
#df_all[[i]]<-df_all[[i]][which(df_all[[i]]$exp_freq_a1!=0),]

df_all[[i]]<-df_all[[i]][,-c(9)] #delete type, certainty and other useless
df_all[[i]]$position<-as.numeric(as.character(df_all[[i]]$position))
###Fix alelel frequency
###
indices<-df_all[[i]]$exp_freq_a1>0.5
df_all[[i]]$exp_freq_a1[indices]<-(1-df_all[[i]]$exp_freq_a1)[indices]
}

df_info_impute_1kgp <- data.table::rbindlist(df_all, idcol = TRUE)
df_info_impute_1kgp$panel<-"1000G"
print("we are done reading info files from 1000G")
print(dim(df_info_impute_1kgp))
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
#df_all[[i]]<-df_all[[i]][which(a0_length==1),]

###get length of a1 allele 
a1_length<-unlist(lapply(df_all[[i]]$a1, function(x) nchar(x)))

###throw whatver was an structual variants
#df_all[[i]]<-df_all[[i]][which(a1_length==1 ),]
df_all[[i]]<-df_all[[i]][which(a1_length==1 & a0_length==1 ),]
df_all[[i]]$info<-as.numeric(as.character(df_all[[i]]$info))

#--work with polymorphic SNPs
df_all[[i]]$position<-as.numeric(as.character(df_all[[i]]$position))
df_all[[i]]$exp_freq_a1<-as.numeric(as.character(df_all[[i]]$exp_freq_a1))
#df_all[[i]]<-df_all[[i]][which(df_all[[i]]$exp_freq_a1!=0),]

df_all[[i]]<-df_all[[i]][,-c(9)] #delete type, certainty and other useless

###Fix alelel frequency
###
indices<-df_all[[i]]$exp_freq_a1>0.5
df_all[[i]]$exp_freq_a1[indices]<-(1-df_all[[i]]$exp_freq_a1)[indices]
}

df_info_impute_hrc <- data.table::rbindlist(df_all, idcol = TRUE)
df_info_impute_hrc$panel<-"HRC"

print("we are done reading info files from HRC")
print(dim(df_info_impute_hrc))

setwd("/mnt/mfs/scratch/GT_Admix/ss5505/imputation/batch2/plots/impute2/compare_common_snps/")
###
###Find common positions
###


#Harmonize SNPs based on position, Allele and chromosome number
## .id is the chromosome number 
df_info_impute_1kgp$harmonizedSNP<-paste(df_info_impute_1kgp$.id,df_info_impute_1kgp$position,df_info_impute_1kgp$a0,df_info_impute_1kgp$a1,sep=":")
df_info_impute_hrc$harmonizedSNP<-paste(df_info_impute_hrc$.id,df_info_impute_hrc$position,df_info_impute_hrc$a0,df_info_impute_hrc$a1,sep=":")

intersect_harmonizedSNP<-(Reduce(intersect,list(df_info_impute_1kgp$harmonizedSNP,df_info_impute_hrc$harmonizedSNP)))
print(length(intersect_harmonizedSNP))
index_hrc<-match(intersect_harmonizedSNP,df_info_impute_hrc$harmonizedSNP)
index_1kgp<-match(intersect_harmonizedSNP,df_info_impute_1kgp$harmonizedSNP)
print(length(index_hrc))
print(length(index_1kgp))

##########################
#########################
######################
intersect_position<-(Reduce(intersect,list(df_info_impute_1kgp$position,df_info_impute_hrc$position)))
print(length(intersect_position))

###
###Find index of common positions
###

index_hrc<-match(intersect_position,df_info_impute_hrc$position)
index_1kgp<-match(intersect_position,df_info_impute_1kgp$position)
print(length(index_hrc))
print(length(index_1kgp))

print("Found index of common positions")
###
###Keep only matched positions 
###

df_info_impute_1kgp<-df_info_impute_1kgp[index_1kgp,]
df_info_impute_hrc<-df_info_impute_hrc[index_hrc,]

##
##Check if alleles match or not
##

sum(df_info_impute_1kgp[,4]==df_info_impute_hrc[,4])
sum(df_info_impute_1kgp[,4]!=df_info_impute_hrc[,4])
passed_a0_allele<-which(df_info_impute_1kgp[,4]==df_info_impute_hrc[,4])
failed_a0_allele<-which(df_info_impute_1kgp[,4]!=df_info_impute_hrc[,4])

dim(df_info_impute_1kgp)
dim(df_info_impute_hrc)

##
##Remove if alleles didn't match
##

if(length(failed_a0_allele)>0){
df_info_impute_1kgp<-df_info_impute_1kgp[-c(failed_a0_allele),]
df_info_impute_hrc<-df_info_impute_hrc[-c(failed_a0_allele),]

}

print("Removed SNPs where A0 didn't match")
dim(df_info_impute_1kgp)
dim(df_info_impute_hrc)

print(sum(df_info_impute_1kgp[,5]==df_info_impute_hrc[,5]))
print(sum(df_info_impute_1kgp[,5]!=df_info_impute_hrc[,5]))
passed_a1_allele<-which(df_info_impute_1kgp[,5]==df_info_impute_hrc[,5])
failed_a1_allele<-which(df_info_impute_1kgp[,5]!=df_info_impute_hrc[,5])

dim(df_info_impute_1kgp)
dim(df_info_impute_hrc)

##
##Remove if alleles didn't match
##

if(length(failed_a1_allele)>0){
df_info_impute_1kgp<-df_info_impute_1kgp[-c(failed_a1_allele),]
df_info_impute_hrc<-df_info_impute_hrc[-c(failed_a1_allele),]
}

dim(df_info_impute_1kgp)
dim(df_info_impute_hrc)

print("Removed SNPs where A1 didn't match")
df_info_impute_1kgp$panel<-"1000G"
df_info_impute_hrc$panel<-"HRC"
merged_dfs<-rbind(df_info_impute_1kgp,df_info_impute_hrc)

print("Merged DFs")

uncommon<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.05),]
rare<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.01),]
ultra_rare<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.001),]

write.table(as.data.frame(uncommon %>% group_by(panel) %>% summarize_at(vars(info),median)),"info_uncommon",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(rare %>% group_by(panel) %>% summarize_at(vars(info),median)),"info_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(ultra_rare %>% group_by(panel) %>% summarize_at(vars(info),median)),"info_ultra_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

print("Compared for all commmon SNPs")

info_uncommon<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.05 & merged_dfs$info>=0.40),]
info_rare<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.01  & merged_dfs$info>=0.40),]
info_ultra_rare<-merged_dfs[which(merged_dfs$exp_freq_a1<=0.001  & merged_dfs$info>=0.40),]

write.table(as.data.frame(info_uncommon %>% group_by(panel) %>% summarize_at(vars(info),median)),"0.4info_uncommon",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_rare %>% group_by(panel) %>% summarize_at(vars(info),median)),"0.4info_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
write.table(as.data.frame(info_ultra_rare %>% group_by(panel) %>% summarize_at(vars(info),median)),"0.4info_ultra_rare",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

print("Performed for all common SNPs with Info >0.40")
print("Done printing higher level MAF info and values")
###

binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 1, by=0.05)))
binned_snps<-binned_snps[complete.cases(binned_snps),]
median_all<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_all,"boxPlots_all",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

##################################
binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1, include.lowest = TRUE,seq(0, 0.05, by=0.005)))
binned_snps<-binned_snps[complete.cases(binned_snps),]
median_uncommon<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_uncommon,"boxPlots_uncommon",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
##################################
binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 0.01, by=0.001)))
binned_snps<-binned_snps[complete.cases(binned_snps),]
median_ultra<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_ultra,"boxPlots_ultra",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
print("Done with ultra rare")

print("Done with cutoff printing")
####

merged_dfs<-merged_dfs[which(merged_dfs$info>=0.40),]

binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 1, by=0.05)))
binned_snps<-binned_snps[complete.cases(binned_snps),]
median_all<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_all,"0.40boxPlots_all",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
print("done with All MAF")
##################################
binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1, include.lowest = TRUE,seq(0, 0.05, by=0.005)))
binned_snps<-binned_snps[complete.cases(binned_snps),]
median_uncommon<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_uncommon,"0.40boxPlots_uncommon",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
print("done with 5% MAF")
##################################
binned_snps<-transform(merged_dfs,binned=cut(merged_dfs$exp_freq_a1,include.lowest = TRUE, seq(0, 0.01, by=0.001)))
binned_snps<-binned_snps[complete.cases(binned_snps),]
median_ultra<-as.data.frame(binned_snps %>% group_by(binned,panel) %>% summarize_at(vars(info),median))
write.table(median_ultra,"0.40boxPlots_ultra",quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
print("Done with ultra rare")

print("Leaving script")

##################################
##################################
##################################




