#!/bin/Rscript

#
#Date 26/10/2018
#Sanjeev Sariya

library(plyr)
#
#
#Merge data from merge_pheno_admixture.R and global admixture
#
#
#

file.pheno.mds<-"mergedmds_oldnewIIDs_order_10252018"

#Old IDs
file.globaladmixture<-"global_admixture.txt"

#
#be careful of merging on key ids
#


df.pheno.mds<-read.table(file.pheno.mds,header=TRUE)
df.global_admxiture<-read.table(file.globaladmixture,header=TRUE)
print(dim(df.pehno.mds))
print(dim(df.global_admxiture))

#
#Delete columns
#

df.pheno.mds_cleaned<-df.pheno.mds[ , -which(names(df.pheno.mds) %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20","V21","V22","V23","V24","V25","V26"))]
print(dim(df.pheno.mds_cleaned))

merged.globaladmixture_pheno.mds<-left_join(df.pheno.mds_cleaned,df.global_admxiture,by=c("OldIID"="IID"))
ordered_merged.globaladmixture_pheno.mds<- merged.globaladmixture_pheno.mds[order(merged.globaladmixture_pheno.mds$ORDER),] 

write.table(ordered_merged.globaladmixture_pheno.mds,"ordered.pheno.mds.global_admixture",row.names=FALSE, col.names=TRUE,quote=FALSE, sep="\t")
