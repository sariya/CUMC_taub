#!/bin/R

###Date March 07 2018
###Sanjeev Sariya
###Dr. Tosto!
###Regional Association plot
###
## xvfb-run /mnt/mfs/cluster/bin/R-3.4/bin/Rscript regional_association.R
###
###
###Load libraries
###

rm(list=ls())
library("Cairo")
library("QCGWAS",lib.loc = "~/R_LIB")

###
###Store Chromosome for which you'd like regional assocition plots
###


list_chr<-c(11,6,19)
###
###Parse data- split them, give headers!
###

getwd() #get wd ---

dta<-read.table("small_chr.txt",header = TRUE,as.is = TRUE)
dta<-dta[,c(2,11)] #get rid of extra columns..
head(dta)
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
    plotname<-paste("RAP_CHR",i,sep = "")
    CairoJPEG(filename = plotname,quality = 75)
    
    dta_chr<-rs_chrmatrix[which(rs_chrmatrix$CHR == i),]
    print(paste("We are doing for chromosome",i,sep = " "))

    sort_chr_dta <- arrange(dta_chr,dta_chr$PVALUE)
    max_hitSNP_pos<-sort_chr_dta[1,2]

    print(paste("Max position on CHR ",i,max_hitSNP_pos),sep=" ")
    
    if( (max_hitSNP_pos - 500000) > 0){
        print("SNP position is OK to do +- 500KB")
        rpstart_pos<-max_hitSNP_pos - 500000
        rpend_pos<-max_hitSNP_pos + 500000

        plot_regional(sort_chr_dta,,data_name="Model 2 PR", chr=i,start_pos=rpstart_pos,end_pos = rpend_pos, save_dir = getwd(),save_name=plotname)
        dev.off()
    }
    
    ###
    ##If loop ends
    ###

}

###
###For loop ends for list of chrosmosomes
###
