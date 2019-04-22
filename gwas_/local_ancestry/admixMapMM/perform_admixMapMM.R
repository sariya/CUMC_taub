#!/usr/bin/Rscript

#
#Date 04/22/2019
#Sanjeev Sariya
#Local Ancestry
#https://rdrr.io/bioc/GENESIS/man/admixMapMM.html
library(GWASTools)
library(gdsfmt)

print("Begin Script")
library(dplyr)

#
#Set file names
#
file.afr<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/merged_CHR22_YRI.txt"
file.nat<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/merged_CHR22_NAT.txt"
file.ceu<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/merged_CHR22_CEU.txt"
file.snps<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/CHR22_snps_rfmix"

#
#Read in files and get their dimensions
#

df.snps<-read.csv(file.snps,header=FALSE)
df.afr<-read.csv(file.afr,header=FALSE)
df.nat<-read.csv(file.nat,header=FALSE)
df.ceu<-read.csv(file.ceu,header=FALSE)

print(dim(df.afr))
print(dim(df.nat))
print(dim(df.ceu))
print(dim(df.snps))

#
#Check nrow in order to ensure data look good.
#		
if(nrow(df.afr)!=nrow(df.nat) | nrow(df.ceu)!=nrow(df.nat)){
print("issue with nrow nat ceu yri")
}

if(nrow(df.afr)!=nrow(df.nat) | nrow(df.ceu)!=nrow(df.snps)){
stop("issue with SNPs nrow")
}

#
#Check and reading ends
#
print("Ending Script")

