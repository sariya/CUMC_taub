#!/usr/bin/bash

##
##Sanjeev Sariya
##02/28/2020
##  module load R/3.4

## Rscript /home/ss5505/scripts_cumc/automate_analyses/perform_left_join_forfiles/leftjoin_data_with.headers.R -f1 sorted_CHR21 -f2 rsids_CHRPOS -c1 1  -c2 1 -o CHR21_sorted_CHRPOS

suppressMessages(library(argparse))
suppressMessages(library(dplyr))

parser <- ArgumentParser(description="Left joint for files take col number for both")

parser$add_argument('-f1',"--file1",help="File that is base or reference",required=TRUE)
parser$add_argument('-f2',"--file2",help="Second file",required=TRUE)
parser$add_argument('-c1',"--col1",help="column number in first file",required=TRUE)
parser$add_argument('-c2',"--col2",help="column number in second file",required=TRUE)
parser$add_argument('-o',"--out",help="output file name",required=TRUE)

args <- parser$parse_args() #make it a data structure
file.ref<-normalizePath(args$file1) ##keep this as reference
file.merge<-normalizePath(args$file2)
col1<-as.numeric(args$col1) ##column in the reference/base columns
col2<-as.numeric(args$col2) ##column in file to be matched columns
output_filename<-args$out

df.reference<-read.table(file.ref,header=TRUE)
df.merge<-read.table(file.merge,header=TRUE)

colname1<-colnames(df.reference)[as.numeric(col1)] ##set these later file performing inner join
colname2<-colnames(df.merge)[as.numeric(col2)]

#be careful with joining, colnames are flipped here
df.leftjoin <-suppressWarnings(left_join(df.reference,df.merge, by=setNames(colname2, colname1) ))
## https://stackoverflow.com/a/28399122/2740831

write.table(x=df.leftjoin, file=output_filename,sep="\t",quote=FALSE, append=FALSE,row.names=FALSE, col.names=TRUE)
print("Exiting code now")

