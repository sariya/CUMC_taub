#!/bin/R

###
###Date 03 15 2018
###

###
###Read merged frame from assoc data
###Split them. Get Chromsomsome
###Delete SNPs of uninterested AF
###

library("Cairo")
library("qqman",lib.loc = "~/R_LIB")

getwd() #get wd ---

dta<-read.table("merged.txt",header = TRUE,as.is = TRUE)

dta<-dta[,c(2,7,11)] #get rid of extra columns.. 
print(paste("We have total SNPs as ",dim(dta),sep=" "))
dta$af<-as.numeric(as.character(dta$af)) #make AF as numeric

dta<-dta[( dta$af>0.1) & ( dta$af<0.9), ] #throw snps of rare AFs
dta<-dta[,-c(2)] #delete AF column

print(paste("We have total SNPs after AF cleaning as ",dim(dta),sep=" "))

temp<-strsplit(dta$rs,"\\:") #split column for RS 
rs_chrmatrix <- do.call(rbind,temp) #do binding 

rs_chrmatrix<-rs_chrmatrix[,-c(3,4)] #throw away A1 and A2 alleles
rs_chrmatrix <- as.data.frame(rs_chrmatrix,stringsAsFactors = FALSE)
rs_chrmatrix$P<-dta$p_wald #get Pvalues 
rs_chrmatrix$SNP<-dta$rs #get SNP names

colnames(rs_chrmatrix) <- c("CHR","BP","P","SNP")
rs_chrmatrix$P<-as.numeric(rs_chrmatrix$P)
rs_chrmatrix$BP<-as.numeric(rs_chrmatrix$BP)
rs_chrmatrix$CHR<-as.numeric(rs_chrmatrix$CHR)

####
###Do Q-Q plots
###

CairoJPEG(filename = "_plotdeleted_qqman.jpeg",quality = 75)
qq(rs_chrmatrix$P,main="Q-Q plot of GWAS p-values Model 1 Merged", xlim = c(0, 7), ylim = c(0,25), pch = 18,col = "blue4", cex = 1.5, las = 1)
dev.off()

####
#Do manhttan plots
###
CairoJPEG(filename = "deletedManht_plot_qqman.jpeg",quality = 75)

manhattan(rs_chrmatrix,main="Manhattan plot Model 1 Merged", ylim = c(0,25), cex = 1.5, cex.axis = 0.9, col = c("blue4", "orange3"),chrlabs = c(1:22))
dev.off()
###
###
