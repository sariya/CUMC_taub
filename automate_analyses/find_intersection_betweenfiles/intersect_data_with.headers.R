#!/usr/bin/env Rscript

##
##Sanjeev Sariya
##02/18/2020
##  module load R/3.4

## /mnt/mfs/cluster/bin/R-3.4/bin/Rscript ~/scripts_cumc/merge/intersect_data_with.headers.R  \
## -f1 db_CH_N37/trend/model1/merged.sorted \
## -f2 db_EUR/model1/merged.sorted -c1 1 -c2 1
#

## /mnt/mfs/cluster/bin/R-3.4/bin/Rscript ~/scripts_cumc/merge/intersect_data_with.headers.R  \
## -f1 filename1 f2 filename2 -c1 colnumber1 -c2 colnumber2 -o outputfilename

suppressMessages(library(argparse))
suppressMessages(library(dplyr))

parser <- ArgumentParser(description="merge files with intersecting values only")

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
print(colname1)
print(colname2)

df.intersect <-suppressWarnings(inner_join(df.reference,df.merge, by=setNames(colname1, colname2) ))
## https://stackoverflow.com/a/28399122/2740831

write.table(x=df.intersect, file=output_filename,sep="\t",quote=FALSE, append=FALSE,row.names=FALSE, col.names=TRUE)
print("Exiting code now")

