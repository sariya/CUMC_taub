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

args <- parser$parse_args() #make it a data structure

file.old_ids<-normalizePath(args$oldids)
file.king_updated<-normalizePath(args$kingids)
file.pheno_oldids<-normalizePath(args$pheno)
outputprefix<-args$outpre
#
#https://stackoverflow.com/questions/16819956/warning-message-in-invalid-factor-level-na-generated
#this is important for replacement later wasted hours on this on 12/05/2018

df.old_ids<-read.table(file.old_ids,header=FALSE,stringsAsFactors=FALSE)
df.king_ids<-read.table(file.king_updated,header=FALSE,stringsAsFactors=FALSE)
df.pheno_oldids<-read.table(file.pheno_oldids,header=TRUE,stringsAsFactors=FALSE)
colnames(df.king_ids)<-paste("king",colnames(df.king_ids),sep="_")
print(head(df.king_ids))

print(dim(df.old_ids))
print(dim(df.king_ids))
print(dim(df.pheno_oldids))

df.pheno_oldids<-df.pheno_oldids[,-c(grep("V[0-9]+",colnames(df.pheno_oldids)))]
print(paste("pheno of all data",nrow(df.pheno_oldids)))

#
#First do a join on pheno from all batches and old ids data 
#

joined_oldpheno<-left_join(df.old_ids,df.pheno_oldids,by=c("V4"="New_IID"))

print(paste("dim for batch input iids",nrow(df.old_ids)))
print(paste("dim for joined input iids",nrow(joined_oldpheno)))

print(head(joined_oldpheno))

#
#Do a join on joinedpheno from input batches and king
#

pheno_oldids_kingids<-left_join(joined_oldpheno,df.king_ids,by=c("V4"="king_V2"))
print(dim(pheno_oldids_kingids))
#write.table(pheno_oldids_kingids, file = "print_test", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA",  row.names = FALSE,col.names = TRUE)

#
#Wherever king columns are NA . copy FID and IIDs
#

for(i in 1:(nrow(pheno_oldids_kingids))){
if(is.na(pheno_oldids_kingids$king_V1[i])){
pheno_oldids_kingids$king_V1[i]<-pheno_oldids_kingids$V3[i]

}
if(is.na(pheno_oldids_kingids$king_V3[i])){
pheno_oldids_kingids$king_V3[i]<-pheno_oldids_kingids$V3[i]
}

if(is.na(pheno_oldids_kingids$king_V4[i])){
pheno_oldids_kingids$king_V4[i]<-pheno_oldids_kingids$V4[i]
}

}

##
##https://stackoverflow.com/questions/5234117/how-to-drop-columns-by-name-in-a-data-frame
##
pheno_oldids_kingids<-pheno_oldids_kingids[, -which(names(pheno_oldids_kingids) %in% c("New_FID","V1","V2","OldFID","OldIID","FID.x","ORDER","FID.y"))]
write.table(pheno_oldids_kingids, file = outputprefix, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA",  row.names = FALSE,col.names = TRUE)








