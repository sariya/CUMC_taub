#!/bin/Rscript

#
#Date 09 13 2018
#Sanjeev Sariya
#
#Perform WGCNA using data from Dr. Amanda Myers. 
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

#

LIBRARY(WGCNA)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
exprs_data<-read.table("test",header=TRUE)

