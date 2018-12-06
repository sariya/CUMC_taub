#!/bin/Rscript


#
#Date 12/06/2018
#PI Dr. Tosto
#Sanjeev Sariya

#
#take quality
#take info file
#take .posterior probability

library("argparse")
library(dplyr)
parser <- ArgumentParser(description="filter snps based on quality")

parser$add_argument('-i',"--info",help="input info file",required=TRUE) #file - that has info information 
parser$add_argument('-p',"--prob",help="probability file",required=TRUE) #file output from king
parser$add_argument('-q',"--qual",help="quality value",required=TRUE) #quality threshold
parser$add_argument('-x',"--outpre",help="output prefix",required=TRUE) # prefix for output

args <- parser$parse_args() #make it a data structure

file.info<-normalizePath(args$info)
file.probs<-normalizePath(args$prob)
quality<-args$qual
outputprefix<-args$outpre




