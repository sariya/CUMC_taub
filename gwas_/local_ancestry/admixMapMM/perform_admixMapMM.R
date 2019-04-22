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

chr<-22
df.snps<-read.csv(file.snps,header=FALSE)
df.afr<-read.csv(file.afr,header=FALSE)
df.nat<-read.csv(file.nat,header=FALSE)
df.ceu<-read.csv(file.ceu,header=FALSE)

print(dim(df.afr))
print(dim(df.nat))
print(dim(df.ceu))
print(dim(df.snps))

colnames(df.snps)<-c("rsids","a1","a2")
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

#
#Read file with rs ids and chr:pos
#

file.rsPos<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/rsids_allchrs_positions"

df.rsPos<-read.table(file.rsPos,header=FALSE)
print(dim(df.rsPos))
df.rsPos$V1<-as.character(df.rsPos$V1)
df.rsPos$V2<-as.character(df.rsPos$V2)

temp<-strsplit(df.rsPos$V1,"\\:") #split column for RS
rs_chrmatrix <- do.call(rbind,temp) #do binding
rs_chrmatrix <- as.data.frame(rs_chrmatrix,stringsAsFactors = FALSE)
rs_chrmatrix$rsids<- df.rsPos$V2
rs_chrmatrix$V1<-as.numeric(as.character(rs_chrmatrix$V1))
rs_chrmatrix$V2<-as.numeric(as.character(rs_chrmatrix$V2))

#
#get length of positions and RSids that match given CHRs
#
print(length(which(rs_chrmatrix$V1==chr)))

rs_chrmatrix<-rs_chrmatrix[which(rs_chrmatrix$V1==chr),]

colnames(rs_chrmatrix)<-c("CHR","POS","rsids")
jointed_rschrpos<-left_join(df.snps,rs_chrmatrix,by=c("rsids"))

if(sum(!complete.cases(jointed_rschrpos)) !=0){

stop("we have some issue jhere with Rsids and positions sum not equal to zero!!!!!")
}


if(nrow(jointed_rschrpos[complete.cases(jointed_rschrpos),]) != nrow(df.ceu)){
stop("we have issue when taking complete cases and nrow with CEU")
}

#
#First three columns are useless in NAT/CEU/AFR 
#

df.afr<-df.afr[,-c()]
df.ceu<-df.ceu[,]
df.nat<-df.nat[,]
print("Ending Script")

