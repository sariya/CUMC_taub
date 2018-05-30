#!/bin/R

###- Date  March 13 2018
### - Prepare data for QQman from Gene based results of GCTA
###

library("argparse")
library("Cairo")
library("qqman",lib.loc = "~/R_LIB")

getwd() #get wd ---

parser <- ArgumentParser(description="Run Q-Q plot for gene based results from GCTA")
parser$add_argument('-o',"--odir",help="Store output directory",required=TRUE) #store output directory
parser$add_argument('-x',"--pre",help="Store prefix for output file",required=TRUE) #store prefix for output file
parser$add_argument('-g',"--gen",help="Location for geenbased  file",required=TRUE) #genotype file - that has names of .raw file

args <- parser$parse_args() #make it a data structure

gen_file<-normalizePath(args$gen) #make into full path and store it
out_dir<-paste(normalizePath(args$odir),"/",sep="") #make into full path and store it
prefix<-args$pre #store prefix and use it in outputting results

setwd(out_dir)
dta<-read.table(gen_file,header = TRUE,as.is = TRUE)
colnames(dta)
dim(dta)
head(dta$Pvalue)
image_name<-paste(prefix,"_qqman.jpeg",sep="")
CairoJPEG(filename = image_name,quality = 75)

main_plot<-paste("Q-Q plot of Genebased p-values",prefix,sep=" ")
qq(dta$Pvalue,main=main_plot, xlim = c(0, 7), ylim = c(0,8), pch = 18,col = "blue4", cex = 1.5, las = 1)
dev.off()

getwd()
