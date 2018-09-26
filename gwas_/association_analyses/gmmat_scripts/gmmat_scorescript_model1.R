#!/bin/Rscript

#
#Date: 08/01/2018
#Sanjeev Sariya
#Model 1 - no APOE. Adjust for age, sex, MDS and kinship

library(GMMAT)
library("argparse")

parser <- ArgumentParser(description="Run GMMAT analysis")

parser$add_argument('-d',"--dosage",help="Location for dosage file",required=TRUE) # dosage file 
parser$add_argument('-s',"--snp",help="Location for snp file",required=TRUE) # file with SNP names
parser$add_argument('-w',"--wdir",help="location where .map and .raw files are present",required=TRUE) #directory that has all .raw and .map files
parser$add_argument('-o',"--odir",help="Store output directory",required=TRUE) #store output directory
parser$add_argument('-p',"--pheno",help="Store phenotype file",required=TRUE) #store phenotype file
parser$add_argument('-x',"--pre",help="Store prefix for output file",required=TRUE) #store prefix for output file
parser$add_argument('-m',"--matrix",help="Store prefix for output file",required=TRUE) #store matrix of genetic relatedness

args <- parser$parse_args() #make it a data structure

matrix_file<-normalizePath(args$matrix) #make into full path and store it
dosage_file<-normalizePath(args$dosage) #make into full path and store it
snp_file<-normalizePath(args$snp) #make into full path and store it
work_dir<-paste(normalizePath(args$wdir),"/",sep="") #make into full path and store i
out_dir<-paste(normalizePath(args$odir),"/",sep="") #make into full path and store it
pheno_file<-normalizePath(args$pheno) #make into full path and store it
prefix<-args$pre #store prefix and use it in outputting results

print(matrix_file)
print(dosage_file)
print(snp_file)
print(work_dir)
print(out_dir)
print(pheno_file)
print(prefix)
outfile<-paste(out_dir,prefix,sep="/")
print(outfile)

setwd(work_dir)
pheno_all <- read.table(pheno_file, header = TRUE, na.strings="NA")
print(dim(pheno_all))

#file_grm_all<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/CHGWAS_analyses_data/multi_merged_data/merge_MESA_NOMAS_HGWAS126/GWAS_files/different_kinship_models/giuseppe_gwas_trials/result.cXX.txt"
grm_all <- as.matrix(read.table(matrix_file))
print(dim(grm_all))

#####Arguments: 
# 1)kinship
# 2) Pheno file with MDS, headers: AD, age, sex and mds1, mds2, mds3
# 3) SNP names from dosage file
# 4) dosage input file
#

print("Running model0")
model0_all <- glmmkin(AD ~ age + sex +mds1 +mds2 + mds3, data = pheno_all, kins = grm_all, family = binomial(link = "logit"))
print("model0 complete")

select <- match(1:10758, names(model0_all$Y))
select[is.na(select)] <- 0

print("Running score")
glmm.score(model0_all, infile = dosage_file,  outfile =outfile,select=select, infile.ncol.skip = 3, infile.ncol.print = 1:3, infile.header.print = c("SNP", "Allele1", "Allele2"),infile.sep = ",")

print("Running score complete")
#output<-glmm.wald(fixed =  AD ~ age + sex + mds1 +mds2 + mds3, data = pheno_all, kins = grm_all, snps=as.character((read.table(snp_file,header=FALSE))$V1),  infile = dosage_file,  infile.ncol.skip = 3,  infile.ncol.print = 1:3, infile.header.print = c("SNP", "Allele1", "Allele2"))
########## ########## ########## ########## ########## 
#print(dim(output))

#write.table(output,outfile,quote=FALSE, row.names=FALSE, col.names=TRUE,sep='\t')
############print("we are done with wald test")










