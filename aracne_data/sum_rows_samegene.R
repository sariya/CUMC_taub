#!/bin/Rscript

#
#Take in cleaned data - no sequence, no annotation information, etc.
#Just matrix with first column as gene value and other as individuals and their expression from gene network
#gene column has duplicate genes. So need to add them before we submit jobs. 

#
#Date - 09/11/2018
#Sanjeev Sariya
#

library("dplyr")
library("argparse")

library(plyr)

matrix_file<-"/mnt/mfs/hgrcgrid/shared/GT_ADMIX/EPISTASIS/gene_expression_data/amanda_myers/all_cases/fixed_spaces_leading_spaces"
gene.exp_mat<-read.table(matrix_file,header=TRUE)
print(dim(gene.exp_mat))

#
####https://stackoverflow.com/a/44230563/2740831
################
###ddply(gene.exp_mat,.(gene),nrow)

gene.exp_mat %>%   group_by(gene) %>%   summarise_all(sum) %>%  data.frame() -> newdf_genes_summed # so that newdf can further be used, if needed

print("Summed genes")
write.table(newdf_genes_summed,"duplicate_genes_summed",sep = "\t",quote=FALSE,col.names = TRUE,row.names = FALSE)