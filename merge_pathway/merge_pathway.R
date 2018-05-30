#!/bin/R

###Sanjeev Sariya
###Date March 28  2018
###Dr. Tosto

###pathway analyses 
### merge Pathway analysies
### 
####
####


getwd()

library("data.table")
library("dplyr")
library("argparse")

parser <- ArgumentParser(description="Get RS from 1000GP for HGWAS SNPs associated")
parser$add_argument('-r',"--repl",help="Location for 1000kgp RS and SNPs",required=TRUE) # replication data 
parser$add_argument('-d',"--disc",help="Location for hgwas position and pvalue file",required=TRUE) # input hgwas data
parser$add_argument('-o',"--odir",help="Store output directory",required=TRUE) #store output directory
parser$add_argument('-x',"--pre",help="Store output directory",required=TRUE) #store prefix
args <- parser$parse_args() #make it a data structure


file_rep<-normalizePath(args$repl) #  store 1KGP file
file_disc<-normalizePath(args$disc) #store hgwas p value and snps
out_dir<-paste(normalizePath(args$odir),"/",sep="") #make into full path and store it

setwd(out_dir)
dta_repl<-read.table(file_rep,header=TRUE,sep='\t') # read disc data
dta_disc<-read.table(file_disc,header=TRUE,sep='\t') # read metal file

colnames(dta_repl)<-c("Set","repl_Size","repl_Count","replz.score","repl_Adj..z.score","repl_p.value","replq.value","repl_List.of.genes")

repl_hgwas<-left_join(dta_disc,dta_repl,by = c("Set" = "Set"))

comp_repl_hgwas<-repl_hgwas[complete.cases(repl_hgwas), ]

pvalue_sorted_disc<- comp_repl_hgwas[order(comp_repl_hgwas[,6]),]

outfile_name<-paste(args$pre,"_merged_path_pvalue",sep="")
print(colnames(pvalue_sorted_disc))

pvalue_sorted_disc_order<-pvalue_sorted_disc[c(1,2,9,3,10,4,11,5,12,6,13,7,14,8,15)]
getwd()
write.table(pvalue_sorted_disc_order,outfile_name,col.names = TRUE,sep="\t",row.names = FALSE,quote = FALSE)