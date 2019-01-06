#!/bin/Rscript

#Date 01/06/2019
#Sanjeev Sariya
#Get union between group files of two batches
#This is hard coded for group files for 5 batches
#

library(argparse)
library(dplyr)


#
#make parser and hold arguments for later data printing and cleaning
#
parser <- ArgumentParser(description="make group files' union")
parser$add_argument('-f1',"--file1",help="input viterbifile",required=TRUE) #group file one 
parser$add_argument('-f2',"--file2",help="output prefix",required=TRUE) # group file two
parser$add_argument('-f3',"--file3",help="output prefix",required=TRUE) # group file three
parser$add_argument('-f4',"--file4",help="output prefix",required=TRUE) # group file four 
parser$add_argument('-f5',"--file5",help="output prefix",required=TRUE) # group file five 
parser$add_argument('-x',"--outpre",help="output prefix",required=TRUE) # prefix for output
parser$add_argument('-o',"--outdir",help="output prefix",required=TRUE) # output dir

args <- parser$parse_args() #make it a data structure

group_file1<-args$file1
group_file2<-args$file2
group_file3<-args$file3
group_file4<-args$file4
group_file5<-args$file5

out_dir<-args$outdir
prefix<-args$outpre

#
#Read group files
#

df.group1<-read.table(group_file1,header=FALSE)
df.group2<-read.table(group_file2,header=FALSE)
df.group3<-read.table(group_file3,header=FALSE)
df.group4<-read.table(group_file4,header=FALSE)
df.group5<-read.table(group_file5,header=FALSE)

print(dim(df.group1))
print(dim(df.group2))
print(dim(df.group3))
print(dim(df.group4))
print(dim(df.group5))

temp_union<-bind_rows(df.group1,df.group2,df.group3,df.group4,df.group5) #https://stackoverflow.com/a/36868001/2740831

print(dim(temp_union))
print("binding complete")

print(paste("Duplicate rows are ",sum(duplicated(temp_union)))) #https://stackoverflow.com/questions/15589601/print-string-and-variable-contents-on-the-same-line-in-r
union_nodups<-temp_union[!duplicated(temp_union),]

#
#https://stats.stackexchange.com/questions/6759/removing-duplicated-rows-data-frame-in-r
#
print(dim(union_nodups))

outfile=paste(out_dir,prefix,sep="/")

write.table(union_nodups, file = outfile, append = FALSE, quote = FALSE, sep = "\t",
eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"),
fileEncoding = "")

print("Printed output to directory") 
