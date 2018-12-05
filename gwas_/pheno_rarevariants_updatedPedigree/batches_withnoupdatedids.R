#!/bin/Rscript

#
#Date 12/05/2018
#Sanjeev Sariya
# __location__ CUMC 19th Floor 
#__PI__ Dr. Giuseppe Tosto
#
#
#Pedigree updated from King (v2.1.2)
#work with each batch individually
#

library("argparse")
library(dplyr)
parser <- ArgumentParser(description="merge pheno after king updated")

parser$add_argument('-o',"--oldids",help="ids from clinical geneticas",required=TRUE) #file - that has FID-iid
parser$add_argument('-k',"--kingids",help="ids updated by king",required=TRUE) #file output from king
parser$add_argument('-p',"--pheno",help="pheno based on old ids",required=TRUE) #file old iid FID pheno
parser$add_argument('-x',"--outpre",help="output prefix",required=TRUE) # prefix for output
