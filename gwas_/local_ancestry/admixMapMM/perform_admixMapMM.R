#!/usr/bin/Rscript

#
#Date 04/22/2019
#Sanjeev Sariya
#Local Ancestry
#https://rdrr.io/bioc/GENESIS/man/admixMapMM.html
print("Begin Script")
library(GWASTools)
library(gdsfmt)
library(dplyr)
library(argparse)
print("loaded libraries")
parser <- ArgumentParser(description="perform admix MM") ##
parser$add_argument('-c',"--ceu",help="CEU input file",required=TRUE) ##CEU Rfmix V1 components
parser$add_argument('-t',"--nat",help="NAT input file",required=TRUE) ## NAT Rfmix V1 components##hgdp
parser$add_argument('-y',"--yri",help="YRI/AFR input file",required=TRUE)  ##YRI RFmix components
parser$add_argument('-x',"--pre",help="prefix for output file",required=TRUE) ##prefix for output
parser$add_argument('-s',"--snps",help="input SNPs file",required=TRUE) # store input snps file
parser$add_argument('-n',"--chr",help="input chr number",required=TRUE) #store chr number

args <- parser$parse_args() #make it a data structure

#
#Set file names
#
#file.afr<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/merged_CHR22_YRI.txt"
#file.nat<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/merged_CHR22_NAT.txt"
#file.ceu<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/merged_CHR22_CEU.txt"
#file.snps<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/local_ancestry_analysis/localancestry_analysis/RFMIX_datapreparations/ancestry_components/CHR22/CHR22_snps_rfmix"

file.afr<-normalizePath(args$yri)
file.nat<-normalizePath(args$nat)
file.ceu<-normalizePath(args$ceu)
file.snps<-normalizePath(args$snps)
chr<-args$chr
out_prefix<-args$pre

print(file.afr)
print(file.nat)
print(file.ceu)
print(file.snps)
print(chr)
print(out_prefix)
#
#Read in files and get their dimensions
#

##chr<-22
chr<-as.numeric(as.character(chr)) ##make change to character chromosome number
print(chr)
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
print(chr)

stop("we have some issue jhere with Rsids and positions sum not equal to zero!!!!!")
}

if(nrow(jointed_rschrpos[complete.cases(jointed_rschrpos),]) != nrow(df.ceu)){
print(chr)
stop("we have issue when taking complete cases and nrow with CEU")

}

#
#First three columns are useless in NAT/CEU/AFR 
#
print(dim(jointed_rschrpos))

print("All data look Good")
#df.afr<-df.afr[,-c()]
#df.ceu<-df.ceu[,]
#df.nat<-df.nat[,]

## Creating a GDS file and variable hierarchy
##http://corearray.sourceforge.net/tutorials/gdsfmt/#creating-a-gds-file-and-variable-hierarchy
##

print("Ending Script")

