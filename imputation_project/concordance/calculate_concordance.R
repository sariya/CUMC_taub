#!/bin/Rscript

#
#Date 09/27/2018
#
##Two files 
#1) .raw file from plink genotyped SNPs
#2) .raw for imputed SNPs in genoypted data. Grep SNPs with CHR22 from .dosage files. 
### plink --gen imputed_genotyped_SNPs --sample batch2_1000_CHR22phasedCHR22-duohmm.sample --hard-call-threshold 0.1 --make-bed --out impute2_plink_dosage_plink
###  plink --bfile impute2_plink_dosage_plink --recode A --make-bed --out recoded_imputed_dosage_plink
###

#
#Find concordance. This is to compare tool's accuracy. 1000G reference panel with IMPUTE2 and MaCH-Admix
#
#

######
###Input files for comparison
######

file.recode_imputed<-"recoded_imputed_dosage_plink.raw"
file.recode_genotyped<-"genotyped.raw"

##################
#function to read .raw file - return data frame
read_raw_file<-function(input.raw){
return(read.table(input.raw,header=TRUE))
}
##################

df.imputed_raw<-read_raw_file(file.recode_imputed)
df.genotype_raw<-read_raw_file(file.recode_genotyped)

print(dim(df.imputed_raw))
print(dim(df.genotype_raw))

if( length(which(df.genotype_raw[,2]!=df.imputed_raw[,2])) !=0){

print("IIDs are not in same order")
stop("Exiting due to IIDs mismatch")
}

#
#Ignore columns with pat, sex, and phenotype information
#
df.imputed_raw<-df.imputed_raw[,-c(1,2,3,4,5,6)]
df.genotype_raw<-df.genotype_raw[,-c(1,2,3,4,5,6)]

##ignore Sum with na.rm and divide by total cobination of SNPs
#https://stackoverflow.com/questions/41790445/compare-each-cell-for-equality-in-two-data-frames-of-equal-size-in-r
concordance<-sum(df.genotype_raw==df.imputed_raw,na.rm=TRUE)/(nrow(df.genotype_raw)*ncol(df.genotype_raw))

print(concordance*100)


