#!/bin/R

###- Date  March 05 2018
### - Prepare data for QQman
#SNP P CHR and position headers
####


####
#Load libraries
###

library("Cairo")
library("qqman",lib.loc = "~/R_LIB")

getwd() #get wd ---

dta<-read.table("merged.txt",header = TRUE,as.is = TRUE)
dta<-dta[,c(2,12)] #get rid of extra columns.. 

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

CairoJPEG(filename = "_plot_qqman.jpeg",quality = 75)
qq(rs_chrmatrix$P,main="Q-Q plot of GWAS p-values Model 1", xlim = c(0, 7), ylim = c(0,18), pch = 18,col = "blue4", cex = 1.5, las = 1)
dev.off()

####
#Do manhttan plots
###
CairoJPEG(filename = "Manht_plot_qqman.jpeg",quality = 75)

manhattan(rs_chrmatrix,main="Manhattan plot Model 1", ylim = c(0,18), cex = 1.5, cex.axis = 0.9, col = c("blue4", "orange3"),chrlabs = c(1:22))
dev.off()
###
###
