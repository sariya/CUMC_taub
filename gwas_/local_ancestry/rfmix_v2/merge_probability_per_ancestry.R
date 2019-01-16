#!/bin/Rscript


#
#Date 12/13/2018
#Sanjeev Sariya
#
# https://github.com/slowkoni/rfmix/blob/master/MANUAL.md

#
#
#

library("argparse")
parser <- ArgumentParser(description="")
parser$add_argument('-f',"--fb",help="input forward backward file",required=TRUE) #
parser$add_argument('-c',"--chr",help="output prefix",required=TRUE) # chr
parser$add_argument('-x',"--outpre",help="output prefix",required=TRUE) # prefix for output
args <- parser$parse_args() #make it a data structure


