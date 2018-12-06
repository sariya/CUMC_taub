#!/bin/Rscript


#
#Date 12/06/2018
#PI Dr. Tosto
#Sanjeev Sariya

#take quality #take info file #take .posterior probability

library("argparse")
library(dplyr)
library(data.table)
parser <- ArgumentParser(description="filter snps based on quality info and imputed prob file")

parser$add_argument('-i',"--info",help="input info file",required=TRUE) #file - that has info information 
parser$add_argument('-p',"--prob",help="probability file",required=TRUE) #file probab output
parser$add_argument('-q',"--qual",help="quality value",required=TRUE) #quality threshold
parser$add_argument('-x',"--outpre",help="output prefix",required=TRUE) # prefix for output
parser$add_argument('-c',"--chr",help="chromsome",required=TRUE) # chromosome

args <- parser$parse_args() #make it a data structure

file.info<-normalizePath(args$info)
file.probs<-normalizePath(args$prob)
quality<-args$qual
outputprefix<-args$outpre
chromsome<-args$chr

print(file.info)
print(file.probs)
print(as.numeric(quality))

df.info<-data.table::fread(file.info, showProgress = TRUE,header=TRUE)
print(dim(df.info))

df.info$info<-as.numeric(df.info$info)
df.info<-df.info[which(df.info$info>=quality),]
print(dim(df.info))

print("Done filtering info")

df.imputeprob<-data.table::fread(file.probs, showProgress = TRUE,header=FALSE)
print(dim(df.imputeprob))
index_highquality<-match(df.info$rs_id,df.imputeprob$V2)
df.imputeprob<-df.imputeprob[index_highquality,]
print("Done filtering iomputed probability")

print(dim(df.imputeprob))

outinfo<-paste(paste(outputprefix,paste("CHR",chromsome,sep=""),sep="_"),"_filtered.info",sep="")
outimputed<-paste(paste(outputprefix,paste("CHR",chromsome,sep=""),sep="_"),"_filtered.imputed",sep="")

write.table(df.info, file = outinfo, append = FALSE, quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE )
write.table(df.imputeprob, file = outimputed, append = FALSE, quote = FALSE, sep = " ",row.names = FALSE, col.names = FALSE)


