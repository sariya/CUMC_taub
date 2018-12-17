#!/bin/Rscript

#
#Date 12/17/2018
#Sanjeev Sariya
#Dr. Giuseppe Tosto, Karen Marder
#
######################################
# Gene level normalized data
######################################

## https://cran.r-project.org/web/packages/myTAI/vignettes/Expression.html
#
#Read pheno data
#Get case and controls calculate mean per gene in case, then calculate mean per gene in control

#
library(data.table)
library(dplyr)

file.expression<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/gene_expr_rin_site_adjusted/PD_RMAnormalized_core.txt"
file.pheno<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/PDexpression/sanjeev_analyses/gene_expr_rin_site_adjusted/pheno_rin_study"

pd<-data.table::fread(file.expression, showProgress = TRUE) 
transposed<-t(pd)

gene_ids<-pd$V1 #get gene ids
colnames(transposed)<-gene_ids

pheno<-read.table( file.pheno ,header=TRUE)
row.names(pheno)<-pheno$IID
PDmerge<-merge(transposed,pheno, by=0,all=TRUE)
#find missing Cel pheno IIDs
missing_celpehno<-which(is.na(PDmerge$IID))
PDmerge<-PDmerge[-c(missing_celpehno),] #remove data for missing cel phenotypes

PDmerge<-PDmerge[,-c(1)]

print(dim(PDmerge))
PDmerge$AGE_SCALE<-(scale(PDmerge$AGE))[,1]
names<-colnames(PDmerge)

#
# https://stackoverflow.com/questions/9723208/aggregate-summarize-multiple-variables-per-group-e-g-sum-mean
#

#calculate mean of each gene for case and contol
mean_case_control<-as.data.frame(PDmerge %>% group_by(PD) %>% summarise_at(vars(-IID, -CATEGORY, -SEX, -AGE, -PD, -LRKK2, -RIN, -FID, -STUDY), mean ))

