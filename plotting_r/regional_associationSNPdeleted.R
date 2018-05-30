#!/bin/R

###Date March 15 2018
###Sanjeev Sariya
###Dr. Tosto!
###Regional Association plot

###
###Load libraries
###

rm(list=ls())
library("Cairo")
source("/home/ss5505/R_LIB/QC_GWAS_v108_CAIRO.R")
###
###Store Chromosome for which you'd like regional assocition plots
###

list_chr<-c(2,3,6,14)

###
###Parse data- split them, give headers!
###

getwd() #get wd ---

dta<-read.table("merged.txt",header = TRUE,as.is = TRUE)

dta<-dta[,c(2,7,11)] #get rid of extra columns.. 
print(paste("We have total SNPs as ",dim(dta),sep=" "))
dta$af<-as.numeric(as.character(dta$af)) #make AF as numeric

dta<-dta[( dta$af>0.1) & ( dta$af<0.9), ] #throw snps of rare AFs
dta<-dta[,-c(2)] #delete AF column

temp<-strsplit(dta$rs,"\\:") #split column for RS
rs_chrmatrix <- do.call(rbind,temp) #do binding

rs_chrmatrix<-rs_chrmatrix[,-c(3,4)] #throw away A1 and A2 alleles

rs_chrmatrix <- as.data.frame(rs_chrmatrix,stringsAsFactors = FALSE)
rs_chrmatrix$P<-dta$p_wald #get Pvalues
rs_chrmatrix$SNP<-dta$rs #get SNP names

colnames(rs_chrmatrix) <- c("CHR","POSITION","PVALUE","MARKER")

###
###Make them numeric 
###
rs_chrmatrix$PVALUE<-as.numeric(as.character(rs_chrmatrix$PVALUE))
rs_chrmatrix$POSITION<-as.numeric(as.character(rs_chrmatrix$POSITION))
rs_chrmatrix$CHR<-as.numeric(rs_chrmatrix$CHR)

head(rs_chrmatrix)
library("dplyr")

for (i in list_chr){

    ###
    ###Plotting file names
    ###
    plotname<-paste("M1_RAPmerged_CHR",i,sep = "")
    
    dta_chr<-rs_chrmatrix[which(rs_chrmatrix$CHR == i),]
    print(paste("We are doing for chromosome",i,sep = " "))

    sort_chr_dta <- arrange(dta_chr,dta_chr$PVALUE)
    max_hitSNP_pos<-sort_chr_dta[1,2]

    print(paste("Max position on CHR ",i,max_hitSNP_pos),sep=" ")
    
    if( (max_hitSNP_pos - 500000) > 0){
        print("SNP position is OK to do +- 500KB")
        rpstart_pos<-max_hitSNP_pos - 500000
        rpend_pos<-max_hitSNP_pos + 500000

        plot_regional(sort_chr_dta,,data_name="Model 1 Merged", chr=i,start_pos=rpstart_pos,end_pos = rpend_pos, save_dir = getwd(),save_name=plotname)
    }
    
    ###
    ##If loop ends
    ###

}

###
###For loop ends for list of chrosmosomes
###
