#!/bin/R

###Sanjeev Sariya
###Date March 28  2018
###Dr. Tosto

###pathway analyses 
###We need to get RS ids from already prepped file from 1KGP
###Take another file with chr:position:A1:A2 and pvalue
####
#### Rscript grep_rs_pathway.R -o ./ -x hgwasM1 -g
#### /mnt/mfs/hgrcgrid/shared/GT_ADMIX/replication_data/pathway/HGWAS/model1/model1_pvalue 
####-k /mnt/mfs/hgrcgrid/shared/GT_ADMIX/replication_data/pathway/HGWAS/rs_allchr_1kg

getwd()

library("data.table")
library("dplyr")
library("argparse")

parser <- ArgumentParser(description="Get RS from 1000GP for HGWAS SNPs associated")
parser$add_argument('-k',"--kgp",help="Location for 1000kgp RS and SNPs",required=TRUE) # RS 1000GP
parser$add_argument('-g',"--hgwas",help="Location for hgwas position and pvalue file",required=TRUE) # input for HGWAS pvalue
parser$add_argument('-o',"--odir",help="Store output directory",required=TRUE) #store output directory
parser$add_argument('-x',"--pre",help="Store output directory",required=TRUE) #store prefix
args <- parser$parse_args() #make it a data structure


file_kgp<-normalizePath(args$kgp) #  store 1KGP file
file_hgwas<-normalizePath(args$hgwas) #store hgwas p value and snps
out_dir<-paste(normalizePath(args$odir),"/",sep="") #make into full path and store it

setwd(out_dir)
dta_kgp<-read.table(file_kgp,header=TRUE) # read disc data
dta_hgwas<-read.table(file_hgwas,header=TRUE) # read metal file

#print(dim(dta_kgp)) print(dim(dta_hgwas))

kgp_hgwas<-left_join(dta_hgwas,dta_kgp,by = c("rs" = "SNP"))

#print(head(kgp_hgwas))

comp_kgp_hgwas<-kgp_hgwas[complete.cases(kgp_hgwas), ]
print(head(comp_kgp_hgwas))


outfile_name<-paste(args$pre,"_merged_RS_Pvalue",sep="")

getwd()
write.table(comp_kgp_hgwas,outfile_name,col.names = FALSE,sep="\t",row.names = FALSE,quote = FALSE)
