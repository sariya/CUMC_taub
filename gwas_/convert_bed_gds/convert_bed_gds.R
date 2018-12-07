#!/bin/Rscript

#
#Convert .bim/.bed/.fam files to .gds format
#Date 12/07/2018
#Sanjeev Sariya
#Dr. Gius Tosto
#Location PH 19th floor

library(SNPRelate)
library("argparse")

parser <- ArgumentParser(description="filter snps based on quality info and imputed prob file")

parser$add_argument('-f',"--fam",help="input fam file",required=TRUE) #
parser$add_argument('-b',"--bed",help="input bed file",required=TRUE) #
parser$add_argument('-m',"--bim",help="input fam file ",required=TRUE) #
parser$add_argument('-c',"--chr",help="output prefix",required=TRUE) # chr
parser$add_argument('-x',"--outpre",help="output prefix",required=TRUE) # prefix for output

args <- parser$parse_args() #make it a data structure

output_prefix<-args$outpre
chr<-args$chr

output_file<-paste(output_prefix,".gds",sep="")

bed.fn<-normalizePath(args$bed)
fam.fn<-normalizePath(args$fam)
bim.fn <-normalizePath(args$bim)

print("entering snpRelate")
print(output_file)

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, output_file)

print("exiting snpRelate")
