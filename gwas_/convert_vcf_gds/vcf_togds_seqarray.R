#!/bin/Rscript

#
#Date 12/10/2018
#Sanjeev Sariya
#convert .vcf file for annoted snps to GDS. This GDS format is based on seqArray 

#https://www.bioconductor.org/packages/devel/bioc/vignettes/SeqArray/inst/doc/SeqArrayTutorial.html

library(SeqArray)
library(Rcpp)
library("argparse")

parser <- ArgumentParser(description="filter snps based on quality info and imputed prob file")

parser$add_argument('-v',"--vcf",help="input vcf file ",required=TRUE) #
parser$add_argument('-c',"--chr",help="output prefix",required=TRUE) # chr
parser$add_argument('-x',"--outpre",help="output prefix",required=TRUE) # prefix for output
args <- parser$parse_args() #make it a data structure

output_prefix<-args$outpre
chr<-args$chr
output_file<-paste(output_prefix,".gds",sep="")
print("entering seqarray conversion")

vcf.fn <- args$vcf

seqVCF2GDS(vcf.fn, output_file, verbose=FALSE)

print("exiting seqarray conversion")