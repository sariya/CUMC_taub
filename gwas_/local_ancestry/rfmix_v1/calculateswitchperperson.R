#!/bin/Rscript

#
#Date 02/21/2019
#Sanjeev Sariya 
#Take in coded RFmix output and calculate per person switch

#
#Column are persons and rows are SNPs #take in Viterbi output
library(data.table)
library(argpasre)


rfmix_coded.file<-"ll"

df.rfmixcoded<-data.table::fread(rsisd.file, showProgress = TRUE,header=FALSE)
df.dosage<-data.table::fread("test", showProgress = TRUE,header=FALSE,sep=",")
print(dim(df.dosage))



